!----------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
! This module drives the generic version of tracers TOPAZ and CFC.
! </OVERVIEW>
!----------------------------------------------------------------
#include <fms_platform.h>
module ocean_generic_mod

  use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer

  use ocean_types_mod,    only: ocean_thickness_type, ocean_public_type
  use ocean_types_mod,    only: ocean_time_type
  use ocean_types_mod,    only: ocean_grid_type, ocean_domain_type
  use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type
  use ocean_types_mod,    only: ocean_density_type
  use ocean_types_mod,    only: ocean_velocity_type
  use ocean_parameters_mod,      only: rho0

  use mpp_mod,            only: stdout, mpp_error, FATAL, NOTE, WARNING
  use mpp_mod,            only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE

  use field_manager_mod, only: fm_get_index,fm_string_len, fm_new_value
  use generic_tracer, only: generic_tracer_init, generic_tracer_source, generic_tracer_update_from_bottom
  use generic_tracer, only: generic_tracer_coupler_get, generic_tracer_coupler_set, generic_tracer_register_diag
  use generic_tracer, only: generic_tracer_end, generic_tracer_get_list, do_generic_tracer, generic_tracer_register
  use generic_tracer, only: generic_tracer_coupler_zero, generic_tracer_vertdiff_G, generic_tracer_vertdiff_M
  use generic_tracer, only: generic_tracer_diag

  use g_tracer_utils,   only: g_tracer_get_name,g_tracer_get_alias,g_tracer_set_values,g_tracer_get_common
  use g_tracer_utils,   only: g_tracer_get_next,g_tracer_type,g_tracer_is_prog,g_tracer_flux_init
  use g_tracer_utils,   only: g_tracer_send_diag,g_tracer_get_values,g_tracer_get_pointer
  use g_tracer_utils,   only: g_tracer_set_pointer

  use coupler_types_mod, only: coupler_2d_bc_type

  implicit none ; private

  logical :: module_initialized = .false.
  ! identification numbers for mpp clocks
  integer :: id_clock_gt_vertdiff
  integer :: id_clock_gt_source,id_clock_gt_btm,id_clock_gt_diag
  integer :: id_clock_gt_get_vals,id_clock_gt_set_vals,id_clock_gt_sum_sfc_setval

  public do_generic_tracer
  public ocean_generic_sum_sfc
  public ocean_generic_zero_sfc
  public ocean_generic_sbc
  public ocean_generic_init
  public ocean_generic_flux_init
  public ocean_generic_column_physics
  public ocean_generic_end
  public ocean_generic_get_field
  public ocean_generic_get_field_pointer
  public ocean_generic_set_pointer

  interface ocean_generic_get_field
     module procedure ocean_generic_get_field_3D
     module procedure ocean_generic_get_field_4D
  end interface

  interface ocean_generic_set_pointer
     module procedure ocean_generic_set_pointer_3D
     module procedure ocean_generic_set_pointer_4D
  end interface

contains

!ALL PE subroutine on Ocean!  Due to otpm design the fluxes should be initialized like this on ALL PE's!
  subroutine ocean_generic_flux_init

    integer :: ind
    character(len=fm_string_len)   :: g_tracer_name,longname, package,units,old_package
    real :: const_init_value
    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_flux_init'
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer,g_tracer_next

    if (.not. module_initialized) then 
       call generic_tracer_register
    endif

    if(.not. do_generic_tracer) return 

    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) then
       call mpp_error(WARNING, trim(sub_name)// ": No generic tracer in the list.")
       return
    endif

    g_tracer=>g_tracer_list  
    do

       call g_tracer_flux_init(g_tracer)
    
       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo
    
  end subroutine ocean_generic_flux_init

  ! <SUBROUTINE NAME="ocean_generic_init">
  !  <OVERVIEW>
  !   Initialize: Add the generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine:
  !     Initializes the generic tracer packages and adds their tracers to the list
  !     Adds the tracers in the list of generic tracers to the set of MOM tracers 
  !     (i.e., make them elements of T_prog or T_diag)
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_init(Domain,Grid,Time) 
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine ocean_generic_init(Domain,Grid,Time)
    type(ocean_domain_type),     intent(in)     :: Domain
    type(ocean_grid_type),       intent(in)     :: Grid
    type(ocean_time_type),       intent(in)     :: Time


    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_init'
    character(len=fm_string_len)   :: g_tracer_name,longname, package,units,old_package,restart_file
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer,g_tracer_next
    real :: const_init_value
    integer :: ntau, ind, length


    if (module_initialized) return 

    !Register all the generic tracers used and create the list of them. 
    !This can be called by ALL PE's. No array fields allocated.
    call generic_tracer_register

    ntau=3 !for MOM
    !!nnz: The following requires Grid%tmask to be TARGET which it is not.
    !So I cannot pass a pointer to this array
    !  grid_tmask => Grid%tmask 

    call generic_tracer_init(Domain%isc,Domain%iec,Domain%jsc,Domain%jec,&
         Domain%isd,Domain%ied,Domain%jsd,Domain%jed,Grid%nk,ntau,&
         Grid%tracer_axes,Grid%tmask,Grid%kmt,Time%model_time)

    ! Register generic tracer modules diagnostics

    call generic_tracer_register_diag()

    !
    !Add the tracers in the list to tpm control and mark them as type='generic'
    !
    !Get the tracer list
    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) then
       call mpp_error(NOTE, trim(sub_name)// ": No generic tracer in the list.")
       return
    endif

    !For each tracer name get its T_prog index and set its field
    old_package=''

    g_tracer=>g_tracer_list  
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname',longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'ocean_restart_file',restart_file)
       call g_tracer_get_values(g_tracer,g_tracer_name,'package', package)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',   units)
       call g_tracer_get_values(g_tracer,g_tracer_name,'const_init_value', const_init_value )

       !The following is a otpm hack. 
       !I can't find any other way to shut off otpm from complaining about the new namelists
       !added in field_table. 
       !otpm is severly restricting developers to add new fields to the field_table for ocea_mod

       if(package .ne. old_package) then

          if (fm_new_value('/ocean_mod/GOOD/good_namelists', package, append = .true.) .le. 0) then
             call mpp_error(FATAL, trim(sub_name) //                           &
                  ' Could not add ' // trim(package) // ' to "good_namelists" list')
          endif

          ind=otpm_set_tracer_package(package, restart_file=restart_file)
          old_package=package
       endif

       !Tell otpm about the generic_tracers so that they can be added to T_prog and T_diag.

       if(g_tracer_is_prog(g_tracer)) then
          if(const_init_value .ne. 0.0) then
             ind=otpm_set_prog_tracer(g_tracer_name, package, longname = longname, &
                  restart_file=restart_file,&
                  units =units, type='generic', const_init_tracer=.true.,const_init_value=const_init_value)       
          else
             ind=otpm_set_prog_tracer(g_tracer_name, package, longname = longname, &
                  restart_file=restart_file,&
                  units =units, type='generic')
          endif
       else
          if(const_init_value .ne. 0.0) then
             ind=otpm_set_diag_tracer(g_tracer_name, package, longname = longname, &
                  restart_file=restart_file,&
                  units =units, type='generic', const_init_tracer=.true.,const_init_value=const_init_value)       
          else
             ind=otpm_set_diag_tracer(g_tracer_name, package, longname = longname, &
                  restart_file=restart_file,&
                  units =units, type='generic')
          endif
       endif


       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo

    ! set clock ids
    id_clock_gt_vertdiff         = mpp_clock_id('(Ocean generic tracer: vertdiff) '     ,grain=CLOCK_ROUTINE)
    id_clock_gt_source           = mpp_clock_id('(Ocean generic tracer: source) '       ,grain=CLOCK_ROUTINE)
    id_clock_gt_btm              = mpp_clock_id('(Ocean generic tracer: bottom up) '    ,grain=CLOCK_ROUTINE)
    id_clock_gt_diag             = mpp_clock_id('(Ocean generic tracer: diag) '         ,grain=CLOCK_ROUTINE)
    id_clock_gt_set_vals         = mpp_clock_id('(Ocean generic tracer: set_values) '   ,grain=CLOCK_ROUTINE)
    id_clock_gt_get_vals         = mpp_clock_id('(Ocean generic tracer: get_values) '   ,grain=CLOCK_ROUTINE)
    id_clock_gt_sum_sfc_setval   = mpp_clock_id('(Ocean generic tracer: sumsfcsetv) '   ,grain=CLOCK_ROUTINE)  

    !---tpm remnant ends.
    module_initialized = .true.


  end subroutine ocean_generic_init

  ! <SUBROUTINE NAME="ocean_generic_zero_sfc">
  !  <OVERVIEW>
  !   zero out the coupler values for all generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This is need since MOM coupler values are acumulated and then divided by time ocean steps.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_zero_sfc(IOB_struc)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine ocean_generic_zero_sfc(IOB_struc)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc

    call generic_tracer_coupler_zero(IOB_struc)

  end subroutine ocean_generic_zero_sfc


  ! <SUBROUTINE NAME="ocean_generic_sum_sfc">
  !  <OVERVIEW>
  !   Calculate the surface state and set coupler values 
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine calculates the surface state and set coupler values for 
  !   those generic tracers that havd flux exchange with atmosphere.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_sum_sfc(Disd,Djsd, Ocean, T_prog, Dens, Time ) 
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine ocean_generic_sum_sfc(Disd,Djsd, Ocean, T_prog, Dens, Time )
    integer,                  intent(in)                    :: Disd,Djsd
    type(ocean_public_type), intent(inout)                  :: Ocean
    type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
    type(ocean_time_type), intent(in)                       :: Time
    type(ocean_density_type), intent(in)                    :: Dens

    type(g_tracer_type), pointer    :: g_tracer_list,g_tracer,g_tracer_next
    integer                           :: g_tracer_index
    integer :: indtemp,indsal
    character(len=fm_string_len)      :: g_tracer_name
    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_sum_sfc'

    indtemp=-1
    indsal=-1

    indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
    if (indtemp .le. 0) then
       call mpp_error(FATAL,trim(sub_name) // ' Could not get the temperature index')
    endif

    indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
    if (indsal .le. 0) then
       call mpp_error(FATAL,trim(sub_name) // ' Could not get the salinity index')
    endif

    !Get the tracer list
    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")

    call generic_tracer_coupler_set(Ocean%fields,&
         ST=T_prog(indtemp)%field(:,:,1,Time%taum1),&
         SS=T_prog(indsal)%field(:,:,1,Time%taum1),&
         rho=Dens%rho,& !nnz: required for MOM
         ilb=Disd, jlb=Djsd,&
         tau=Time%taum1)


    !Send diagnostics
    call g_tracer_send_diag(g_tracer_list, Time%model_time, Time%tau)

  end subroutine ocean_generic_sum_sfc

  ! <SUBROUTINE NAME="ocean_generic_sbc">
  !  <OVERVIEW>
  !   Get the coupler values 
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine gets coupler values for 
  !   those generic tracers that have flux exchange with atmosphere.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_sbc(Ice_ocean_boundary_fluxes,Disd,Djsd, T_prog )
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine ocean_generic_sbc(Ice_ocean_boundary_fluxes,Disd,Djsd, T_prog, runoff )
    type(coupler_2d_bc_type), intent(in)                            :: Ice_ocean_boundary_fluxes
    integer,                  intent(in)                            :: Disd,Djsd
    type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
    real, intent(in), dimension(Disd:,Djsd:)                        :: runoff

    type(g_tracer_type), pointer    :: g_tracer_list,g_tracer,g_tracer_next
    integer                           :: g_tracer_index
    character(len=fm_string_len)      :: g_tracer_name
    character(len=fm_string_len), parameter :: sub_name = 'update_generic_tracer_sbc'

    !Extract the tracer surface fields from coupler 
    call generic_tracer_coupler_get(Ice_ocean_boundary_fluxes)

    !Update T_prog fields from generic tracer fields
    !
    !Get the tracer list
    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its T_prog index and get its flux fields
    g_tracer=>g_tracer_list  
    do
       if(g_tracer_is_prog(g_tracer)) then
          call g_tracer_get_alias(g_tracer,g_tracer_name)
          g_tracer_index = fm_get_index(trim('/ocean_mod/prog_tracers/'//g_tracer_name))
          if (g_tracer_index .le. 0) &
               call mpp_error(FATAL,trim(sub_name) // ' Could not get the index for '//g_tracer_name)

          if (_ALLOCATED(g_tracer%stf) )&
               call g_tracer_get_values(g_tracer,g_tracer_name,'stf',   T_prog(g_tracer_index)%stf,   Disd,Djsd)

          if (_ALLOCATED(g_tracer%btf) )&
               call g_tracer_get_values(g_tracer,g_tracer_name,'btf',   T_prog(g_tracer_index)%btf,   Disd,Djsd)

          !If the tracer has runoff fill in the T_prog(n)%trunoff and  T_prog(n)%runoff_tracer_flux
          if(_ALLOCATED(g_tracer%trunoff)) then
             !Fill in T_prog(n)%trunoff

             call g_tracer_get_values(g_tracer,g_tracer_name,'trunoff',T_prog(g_tracer_index)%trunoff,Disd,Djsd)

             !Fill in T_prog(n)%runoff_tracer_flux
             T_prog(g_tracer_index)%runoff_tracer_flux = T_prog(g_tracer_index)%trunoff * runoff

             !Set g_tracer%runoff_tracer_flux
             call g_tracer_set_values(g_tracer,g_tracer_name,'runoff_tracer_flux',T_prog(g_tracer_index)%runoff_tracer_flux,Disd,Djsd)
             !
             !Fill in T_prog(n)%triver in MOM
             !Note: This is done so that MOM can apply the river fluxes through setting either
             !        the runoff and calving fluxes (when discharge_combine_runoff_calve=.false.) 
             !      or
             !        the total river concentration (when discharge_combine_runoff_calve=.true.)
             !
             !Assume zero calving flux for the generic tracers.
             !T_prog(g_tracer_index)%tcalving = 0 !MOM default
             T_prog(g_tracer_index)%triver = T_prog(g_tracer_index)%trunoff 
             
          endif
       
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo
  end subroutine ocean_generic_sbc

  ! <SUBROUTINE NAME="ocean_generic_column_physics">
  !  <OVERVIEW>
  !   Column physics for generic tracers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine:
  !       Update generic tracer concentration fields from sources and sinks.
  !       Vertically diffuse generic tracer concentration fields.
  !       Update generic tracers from bottom and their bottom reservoir.  
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_column_physics(Thickness, hblt_depth, Time, Grid, dtts, Disd,Djsd, T_prog, T_diag,&
  !                                     sw_pen,opacity, diff_cbt, river, Velocity )
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine ocean_generic_column_physics(Thickness, hblt_depth, Time, Grid, dtts, Disd,Djsd, T_prog, T_diag,&
       sw_pen,opacity, diff_cbt, Velocity)
    integer,                  intent(in)                            :: Disd,Djsd
    type(ocean_grid_type), intent(in)                               :: Grid
    type(ocean_time_type), intent(in)                               :: Time
    type(ocean_thickness_type), intent(in)                          :: Thickness
    type(ocean_Velocity_type), intent(in)                           :: Velocity
    real, intent(in), dimension(Disd:,Djsd:)                        :: hblt_depth,sw_pen
    real, intent(in), dimension(Disd:,Djsd:,:)                      :: opacity
    real, intent(in), dimension(Disd:,Djsd:,:,:)                    :: diff_cbt
    real, intent(in)                                                :: dtts
    type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
    type(ocean_diag_tracer_type), dimension(:), intent(inout)       :: T_diag

    type(g_tracer_type), pointer    :: g_tracer_list,g_tracer,g_tracer_next
    integer                           :: g_tracer_index
    integer :: indtemp,indsal
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,nbands,i,j,k
    real, dimension(:),       Allocatable :: max_wavelength_band
    real, dimension(:,:,:),   Allocatable :: sw_pen_band
    real, dimension(:,:,:,:), Allocatable :: opacity_band
    logical, save                         :: initialize_tau_level = .true.

    character(len=fm_string_len)      :: g_tracer_name
    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_column_physics'

    !
    !Update the fields of the generic tracers from T_prog
    !
    !This step is not needed any more since generic tracers %field points to and reuse the corresponding MOM arrays.
    !
    !Update from sources
    indtemp=-1
    indsal=-1

    indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
    if (indtemp .le. 0) then 
       call mpp_error(FATAL,trim(sub_name) // ' Could not get the temperature index')
    endif

    indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
    if (indsal .le. 0) then  
       call mpp_error(FATAL,trim(sub_name) // ' Could not get the salinity index')
    endif

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau)

    !For MOM opacity and sw_pen have 1 band only.
    nbands=1
    allocate(max_wavelength_band(nbands)); max_wavelength_band =0.0
    allocate( sw_pen_band(nbands,isd:ied, jsd:jed));    
    allocate(opacity_band(nbands,isd:ied, jsd:jed,nk)); 

    do j = jsc, jec ; do i = isc, iec  
       sw_pen_band(1,i,j)   = sw_pen(i,j)
       do k = 1, nk 
          opacity_band(1,i,j,k)= opacity(i,j,k)
       enddo
    enddo; enddo
       
    call mpp_clock_begin(id_clock_gt_source)
    
    call generic_tracer_source(T_prog(indtemp)%field(:,:,:,Time%taup1),&
         T_prog(indsal)%field(:,:,:,Time%taup1), Thickness%rho_dzt(:,:,:,Time%taup1), Thickness%dzt,&
         hblt_depth, Disd, Djsd, Time%taup1,dtts,Grid%dat, Time%model_time,&
         nbands, max_wavelength_band, sw_pen_band, opacity_band, Grid%ht, Velocity%current_wave_stress(:,:))

    call mpp_clock_end(id_clock_gt_source)
    deallocate(max_wavelength_band,sw_pen_band,opacity_band)

    ! Explicit Vertical Diffusion

    call mpp_clock_begin(id_clock_gt_vertdiff)

    call generic_tracer_vertdiff_M(Thickness%rho_dzt(:,:,:,Time%taup1), Thickness%dzwt, diff_cbt(:,:,:,1),&
         dtts, rho0, Time%taup1)

    call mpp_clock_end(id_clock_gt_vertdiff)

    ! Update bottom fields after vertical processes
    call mpp_clock_begin(id_clock_gt_btm)
    call generic_tracer_update_from_bottom(dtts, Time%taup1, Time%model_time)
    call mpp_clock_end(id_clock_gt_btm)

    !
    ! Finish up generic tracers
    !
    call mpp_clock_begin(id_clock_gt_diag)
    call generic_tracer_diag(Disd, Djsd, Time%tau, Time%taup1, dtts, Time%model_time, Thickness%dzt,  &
         Thickness%rho_dzt(:,:,:,Time%tau), Thickness%rho_dzt(:,:,:,Time%taup1))
    call mpp_clock_end(id_clock_gt_diag)

    !!nnz: the following is necessary if generic tracers are allocated by MOM
    !
    !Update T_prog fields from generic tracer fields
    !This step is not needed for %fields any more since generic tracers %field points to and reuse the corresponding MOM arrays.
    !
    !Get the tracer list
    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    g_tracer=>g_tracer_list  
    call mpp_clock_begin(id_clock_gt_get_vals)
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)
       if(g_tracer_is_prog(g_tracer)) then
          g_tracer_index = fm_get_index(trim('/ocean_mod/prog_tracers/'//g_tracer_name))
          if (g_tracer_index .le. 0) &
               call mpp_error(FATAL,trim(sub_name) // ' Could not get the index for '//g_tracer_name)

          if (_ALLOCATED(g_tracer%btf) )&
               call g_tracer_get_values(g_tracer,g_tracer_name,'btf',   T_prog(g_tracer_index)%btf,   Disd,Djsd)
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo
    call mpp_clock_end(id_clock_gt_get_vals)

  end subroutine ocean_generic_column_physics

  subroutine ocean_generic_get_field_pointer(name, field)
    character(len=*),         intent(in)    :: name    
    real, dimension(:,:,:), pointer :: field
    type(g_tracer_type), pointer      :: g_tracer_list

    call generic_tracer_get_list(g_tracer_list)
    call g_tracer_get_pointer(g_tracer_list,name,'field', field )
    
  end subroutine ocean_generic_get_field_pointer

  subroutine ocean_generic_get_field_4D(name, field)
    character(len=*),         intent(in)    :: name    
    real, dimension(:,:,:,:), intent(inout) :: field
    type(g_tracer_type), pointer      :: g_tracer_list
    real, dimension(:,:,:,:), pointer :: ptr

    call generic_tracer_get_list(g_tracer_list)
    call g_tracer_get_pointer(g_tracer_list,name,'field', ptr )
    field = ptr
  end subroutine ocean_generic_get_field_4D

  subroutine ocean_generic_get_field_3D(name, field)
    character(len=*),         intent(in)    :: name    
    real, dimension(:,:,:), intent(inout) :: field
    real, dimension(:,:,:), pointer :: ptr
    type(g_tracer_type), pointer    :: g_tracer_list

    call generic_tracer_get_list(g_tracer_list)
    call g_tracer_get_pointer(g_tracer_list,name,'field', ptr )
    field = ptr
    
  end subroutine ocean_generic_get_field_3D

  subroutine ocean_generic_set_pointer_4d(name, member, field, ilb, jlb)
    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_set_pointer_4d'
    character(len=*),                       intent(in)  :: name
    character(len=*),                       intent(in)  :: member
    real, dimension(ilb:,jlb:,:,:), target, intent(in)  :: field
    integer,                                intent(in)  :: ilb
    integer,                                intent(in)  :: jlb
    type(g_tracer_type), pointer      :: g_tracer_list

    call generic_tracer_get_list(g_tracer_list)
    if (associated(g_tracer_list)) then
      call g_tracer_set_pointer(g_tracer_list, name, member, field, ilb, jlb)
    else
      call mpp_error(NOTE, trim(sub_name)// ": No generic tracer in the list. No generic tracers?")
    endif
  end subroutine ocean_generic_set_pointer_4d

  subroutine ocean_generic_set_pointer_3d(name, member, field, ilb, jlb)
    character(len=fm_string_len), parameter :: sub_name = 'ocean_generic_set_pointer_3d'
    character(len=*),                     intent(in)    :: name
    character(len=*),                     intent(in)  :: member
    real, dimension(ilb:,jlb:,:), target, intent(in)    :: field
    integer,                              intent(in)    :: ilb
    integer,                              intent(in)    :: jlb
    type(g_tracer_type), pointer      :: g_tracer_list

    call generic_tracer_get_list(g_tracer_list)
    if (associated(g_tracer_list)) then
      call g_tracer_set_pointer(g_tracer_list, name, member, field, ilb, jlb)
    else
      call mpp_error(NOTE, trim(sub_name)// ": No generic tracer in the list. No generic tracers?")
    endif
  end subroutine ocean_generic_set_pointer_3d

  ! <SUBROUTINE NAME="ocean_generic_end">
  !  <OVERVIEW>
  !   Ends the generic tracer module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !  Call the end for generic tracer module and deallocate all temp arrays 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call ocean_generic_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine ocean_generic_end
    call generic_tracer_end
  end subroutine ocean_generic_end


end module ocean_generic_mod

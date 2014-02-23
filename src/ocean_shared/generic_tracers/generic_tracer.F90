!----------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
! This module provides the main interfaces between Ocean models and 
! generic tracers.
! </OVERVIEW>
!<DESCRIPTION>
! Generic Tracers are designed to be used by both GFDL Ocean models, GOLD and MOM.
! This module provides the main interfaces for using generic tracers.
! Generic Tracers are contained in separate modules according to their
! chemical/physical similarity (currently generic_TOPAZ, generic_COBALT, 
! generic_ERGOM and generic_CFC)
! This module acts as a router for these various tracer modules and 
! routes the subroutine calls to the appropriate tracer module.
! It also maintains a (linked) list of all generic tracers created in 
! the experiment. This list acts as the "state" of generic tracers and
! contains all the information for all such tracers. This module provides
! a subroutine to query its state at any time. 
!</DESCRIPTION>
! <REFERENCE>
! http://cobweb.gfdl.noaa.gov/~nnz/MITeam_GUTS_022708.pdf 
! </REFERENCE>
! <REFERENCE>
! http://cobweb.gfdl.noaa.gov/~nnz/wiki/doku.php?id=genericunifiedtracers 
! </REFERENCE>
! 
!----------------------------------------------------------------

module generic_tracer

  use fms_mod,           only: open_namelist_file, close_file, check_nml_error  
  use field_manager_mod, only: fm_string_len
  use mpp_mod, only : input_nml_file, mpp_error, NOTE, WARNING, FATAL, stdout, stdlog
  use time_manager_mod, only : time_type
  use coupler_types_mod, only : coupler_2d_bc_type

  use g_tracer_utils, only : g_tracer_type, g_tracer_init, g_diag_type
  use g_tracer_utils, only : g_tracer_get_common, g_tracer_set_common, g_tracer_is_prog
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get, g_tracer_register_diag
  use g_tracer_utils, only : g_tracer_vertdiff_M, g_tracer_vertdiff_G, g_tracer_get_next     
  use g_tracer_utils, only : g_tracer_diag

  use generic_CFC,    only : generic_CFC_register
  use generic_CFC,    only : generic_CFC_init, generic_CFC_update_from_source,generic_CFC_update_from_coupler
  use generic_CFC,    only : generic_CFC_set_boundary_values, generic_CFC_end, do_generic_CFC

  use generic_ERGOM, only : generic_ERGOM_register, generic_ERGOM_register_diag
  use generic_ERGOM, only : generic_ERGOM_init, generic_ERGOM_update_from_source,generic_ERGOM_update_from_coupler
  use generic_ERGOM, only : generic_ERGOM_set_boundary_values, generic_ERGOM_end, do_generic_ERGOM
  use generic_ERGOM, only : generic_ERGOM_update_from_bottom

  use generic_TOPAZ,  only : generic_TOPAZ_register
  use generic_TOPAZ,  only : generic_TOPAZ_init, generic_TOPAZ_update_from_source,generic_TOPAZ_register_diag
  use generic_TOPAZ,  only : generic_TOPAZ_update_from_bottom,generic_TOPAZ_update_from_coupler
  use generic_TOPAZ,  only : generic_TOPAZ_set_boundary_values, generic_TOPAZ_end, do_generic_TOPAZ

  use generic_BLING,  only : generic_BLING_register
  use generic_BLING,  only : generic_BLING_init, generic_BLING_update_from_source,generic_BLING_register_diag
  use generic_BLING,  only : generic_BLING_update_from_bottom,generic_BLING_update_from_coupler
  use generic_BLING,  only : generic_BLING_set_boundary_values, generic_BLING_end, do_generic_BLING

  use generic_miniBLING_mod,  only : generic_miniBLING_register
  use generic_miniBLING_mod,  only : generic_miniBLING_init, generic_miniBLING_update_from_source,generic_miniBLING_register_diag
  use generic_miniBLING_mod,  only : generic_miniBLING_update_from_bottom,generic_miniBLING_update_from_coupler
  use generic_miniBLING_mod,  only : generic_miniBLING_set_boundary_values, generic_miniBLING_end, do_generic_miniBLING
  use generic_miniBLING_mod,  only : generic_miniBLING_diag

  use generic_COBALT,  only : generic_COBALT_register
  use generic_COBALT,  only : generic_COBALT_init, generic_COBALT_update_from_source,generic_COBALT_register_diag
  use generic_COBALT,  only : generic_COBALT_update_from_bottom,generic_COBALT_update_from_coupler
  use generic_COBALT,  only : generic_COBALT_set_boundary_values, generic_COBALT_end, do_generic_COBALT

  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_tracer'
  character(len=fm_string_len), parameter :: package_name   = 'generic_tracer'

  public generic_tracer_register
  public generic_tracer_init
  public generic_tracer_register_diag
  public generic_tracer_source
  public generic_tracer_diag
  public generic_tracer_update_from_bottom
  public generic_tracer_coupler_get
  public generic_tracer_coupler_set
  public generic_tracer_coupler_zero
  public generic_tracer_end
  public generic_tracer_get_list
  public do_generic_tracer
  public generic_tracer_vertdiff_G
  public generic_tracer_vertdiff_M
  public generic_tracer_get_diag_list


  !Linked Lists of all prog and diag tracers in this module
  !Ensure these pointers are "save"d between the calls
  type(g_tracer_type), save, pointer :: tracer_list => NULL()

  !Linked Lists of diagnostics fields (specially those that need to be manipulated by the ocean model)
  !Ensure these pointers are "save"d between the calls
  type(g_diag_type), save, pointer :: diag_list => NULL()

  logical, save :: do_generic_tracer = .false.

  !JGJ 2013/05/31  merged COBALT into siena_201303
  namelist /generic_tracer_nml/ do_generic_tracer, do_generic_CFC, do_generic_TOPAZ,    &
       do_generic_ERGOM, do_generic_BLING, do_generic_miniBLING, do_generic_COBALT

contains


  subroutine generic_tracer_register

    integer :: ioun, io_status, ierr
    integer :: stdoutunit,stdlogunit
    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_register'

    stdoutunit=stdout();stdlogunit=stdlog()
    ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=generic_tracer_nml, iostat=io_status)
ierr = check_nml_error(io_status,'generic_tracer_nml')
#else
    ioun = open_namelist_file()
    read  (ioun, generic_tracer_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'generic_tracer_nml')
    call close_file (ioun)
#endif

    write (stdoutunit,'(/)')
    write (stdoutunit, generic_tracer_nml)
    write (stdlogunit, generic_tracer_nml)

    if(do_generic_CFC) &
         call generic_CFC_register(tracer_list)

    if(do_generic_TOPAZ) &
         call generic_TOPAZ_register(tracer_list)

    if(do_generic_ERGOM) &
         call generic_ERGOM_register(tracer_list)

    if(do_generic_BLING) &
         call generic_BLING_register(tracer_list)

    if(do_generic_miniBLING) &
         call generic_miniBLING_register(tracer_list)

    if(do_generic_COBALT) &
         call generic_COBALT_register(tracer_list)

  end subroutine generic_tracer_register


  ! <SUBROUTINE NAME="generic_tracer_init">
  !  <OVERVIEW>
  !   Initialize generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Reads the namelist generic_tracer_nml to find the requested tracer packages
  !   Sets the common properties to be used by ALL generic tracers: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,init_time
  !   Initialize each requested generic tracer package by calling their init routine
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,init_time)
  !  </TEMPLATE>
  !  <IN NAME="isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3)" TYPE="integer">
  !   Domain boundaries,
  !   ntau: number of time steps retained in the concentration field
  !   axes(3): diag axes
  !  </IN>
  !  <IN NAME="init_time" TYPE="">
  !   
  !  </IN>
  !  <IN NAME="" TYPE="type(time_type)">
  !   Initiation time
  !  </IN>
  !  <IN NAME="grid_tmask" TYPE="real, dimension(:,:,:),target">
  !   Grid mask
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_tracer_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                       intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3)
    type(time_type),               intent(in) :: init_time
    real, dimension(:,:,:),target, intent(in) :: grid_tmask
    integer, dimension(:,:)      , intent(in) :: grid_kmt
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next

    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_init'

    call g_tracer_set_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time) 

    !Allocate and initialize all registered generic tracers
    !JGJ 2013/05/31  merged COBALT into siena_201303
    if(do_generic_CFC .or. do_generic_TOPAZ .or. do_generic_ERGOM .or. do_generic_BLING .or. do_generic_miniBLING .or. do_generic_COBALT) then
       g_tracer => tracer_list        
       !Go through the list of tracers 
       do  
          call g_tracer_init(g_tracer)

          !traverse the linked list till hit NULL
          call g_tracer_get_next(g_tracer, g_tracer_next)
          if(.NOT. associated(g_tracer_next)) exit
          g_tracer=>g_tracer_next  

       enddo
    endif    

    !Initilalize specific tracers
    if(do_generic_CFC) &
         call generic_CFC_init(tracer_list)

    if(do_generic_TOPAZ) &
         call generic_TOPAZ_init(tracer_list)

    if(do_generic_ERGOM) &
         call generic_ERGOM_init(tracer_list)

    if(do_generic_BLING) &
         call generic_BLING_init(tracer_list)

    if(do_generic_miniBLING) &
         call generic_miniBLING_init(tracer_list)

    if(do_generic_COBALT) &
         call generic_COBALT_init(tracer_list)

  end subroutine generic_tracer_init

  subroutine generic_tracer_register_diag
    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_register_diag'
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next
    
    !Diagnostics register for the fields common to All generic tracers
    !JGJ 2013/05/31  merged COBALT into siena_201303

    if(do_generic_CFC .or. do_generic_TOPAZ .or. do_generic_ERGOM .or. do_generic_BLING .or. do_generic_miniBLING .or. do_generic_COBALT) then

       g_tracer => tracer_list        
       !Go through the list of tracers 
       do  
          call g_tracer_register_diag(g_tracer)

          !traverse the linked list till hit NULL
          call g_tracer_get_next(g_tracer, g_tracer_next)
          if(.NOT. associated(g_tracer_next)) exit
          g_tracer=>g_tracer_next  

       enddo
    endif    

    !Diagnostics register for fields particular to each tracer module
    
    if(do_generic_TOPAZ)  call generic_TOPAZ_register_diag(diag_list)    

    if(do_generic_ERGOM)  call generic_ERGOM_register_diag(diag_list)    

    if(do_generic_BLING)  call generic_BLING_register_diag()    
    
    if(do_generic_miniBLING)  call generic_miniBLING_register_diag()    

    if(do_generic_COBALT)  call generic_COBALT_register_diag(diag_list)

  end subroutine generic_tracer_register_diag

  ! <SUBROUTINE NAME="generic_tracer_coupler_get">
  !  <OVERVIEW>
  !   Get coupler values
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calls the corresponding generic_X_update_from_coupler routine for each package X.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_coupler_get(IOB_struc) 
  !  </TEMPLATE>
  !  <IN NAME="IOB_struc" TYPE="type(coupler_2d_bc_type)">
  !   Ice Ocean Boundary flux structure
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_tracer_coupler_get(IOB_struc)
    type(coupler_2d_bc_type), intent(in)    :: IOB_struc

    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_coupler_get'
    !All generic tracers
    !Update tracer boundary values (%stf and %triver) from coupler fluxes foreach tracer in the prog_tracer_list
    call g_tracer_coupler_get(tracer_list,IOB_struc)

    !Specific tracers
    !    if(do_generic_CFC)    call generic_CFC_update_from_coupler(tracer_list) !Nothing to do

    if(do_generic_TOPAZ)  call generic_TOPAZ_update_from_coupler(tracer_list)

    if(do_generic_BLING)  call generic_BLING_update_from_coupler(tracer_list)

    if(do_generic_miniBLING)  call generic_miniBLING_update_from_coupler(tracer_list)

    if(do_generic_COBALT)  call generic_COBALT_update_from_coupler(tracer_list)

  end subroutine generic_tracer_coupler_get


  ! <SUBROUTINE NAME="generic_tracer_diag">
  !  <OVERVIEW>
  !   Do things which must be done after all transports and sources have been calculated
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calls the corresponding generic_X_diag routine for each package X.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call  generic_tracer_diag(tau,model_time)
  !  </TEMPLATE>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !  <IN NAME="model_time" TYPE="time_type">
  !   Model time
  !  </IN>
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_tracer_diag(ilb, jlb, tau, taup1, dtts, model_time, dzt, rho_dzt_tau, rho_dzt_taup1)
    integer,                        intent(in) :: ilb
    integer,                        intent(in) :: jlb
    integer,                        intent(in) :: tau
    integer,                        intent(in) :: taup1
    real,                           intent(in) :: dtts
    type(time_type),                intent(in) :: model_time
    real, dimension(ilb:,jlb:,:),   intent(in) :: dzt
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_tau
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_taup1

    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_update_from_diag'

    if(do_generic_miniBLING)  call generic_miniBLING_diag(tracer_list, ilb, jlb, taup1, model_time, dzt, rho_dzt_taup1)

    call g_tracer_diag(tracer_list, ilb, jlb, rho_dzt_tau, rho_dzt_taup1, model_time, tau, taup1, dtts)

    return

  end subroutine generic_tracer_diag


  ! <SUBROUTINE NAME="generic_tracer_source">
  !  <OVERVIEW>
  !   Update the tracers from sources/sinks
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calls the corresponding generic_X_update_from_source routine for each package X.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call  generic_tracer_source(Temp,Salt,rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dtts,&
  !                               grid_dat,sw_pen,opacity,river,neutralrho,grid_yt)
  !  </TEMPLATE>
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="Temp" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean temperature   
  !  </IN>
  !  <IN NAME="Salt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean salinity
  !  </IN>
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  !  <IN NAME="opacity" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean opacity
  !  </IN>
  !  <IN NAME="sw_pen" TYPE="real, dimension(ilb:,jlb:)">
  !   Shortwave peneteration
  !  </IN>
  !  <IN NAME="hblt_depth" TYPE="real, dimension(ilb:,jlb:)">
  !   
  !  </IN>
  !  <IN NAME="grid_dat" TYPE="real, dimension(ilb:,jlb:)">
  !   Grid area
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_tracer_source(Temp,Salt,rho_dzt,dzt,hblt_depth,ilb,jlb,tau,dtts,&
       grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band,      &
       grid_ht, current_wave_stress)
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp,Salt,rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dtts
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time
    integer,                        intent(in) :: nbands
    real, dimension(:),             intent(in) :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band
    real, dimension(ilb:,jlb:),optional, intent(in) :: grid_ht
    real, dimension(ilb:,jlb:),optional ,    intent(in) :: current_wave_stress


    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_update_from_source'

    !    if(do_generic_CFC)    call generic_CFC_update_from_source(tracer_list) !Nothing to do for CFC

    if(do_generic_TOPAZ)  call generic_TOPAZ_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,&
         hblt_depth,ilb,jlb,tau,dtts,grid_dat,model_time,&
         nbands,max_wavelength_band,sw_pen_band,opacity_band)

    if(do_generic_ERGOM)  call generic_ERGOM_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,&
         hblt_depth,ilb,jlb,tau,dtts,grid_dat,model_time,&
         nbands,max_wavelength_band,sw_pen_band,opacity_band,current_wave_stress)

    if(do_generic_BLING)  call generic_BLING_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,&
         hblt_depth,ilb,jlb,tau,dtts,grid_dat,model_time,&
         nbands,max_wavelength_band,sw_pen_band,opacity_band)

    if(do_generic_miniBLING)  call generic_miniBLING_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,&
         hblt_depth,ilb,jlb,tau,dtts,grid_dat,model_time,&
         nbands,max_wavelength_band,sw_pen_band,opacity_band, grid_ht)

    if(do_generic_COBALT)  call generic_COBALT_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,&
         hblt_depth,ilb,jlb,tau,dtts,grid_dat,model_time,&
         nbands,max_wavelength_band,sw_pen_band,opacity_band)

    return

  end subroutine generic_tracer_source

  ! <SUBROUTINE NAME="generic_tracer_update_from_bottom">
  !  <OVERVIEW>
  !   Update the tracers from bottom fluxes
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calls the corresponding generic_X_update_from_bottom routine for each package X.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_update_from_bottom(dt, tau, model_time)
  !  </TEMPLATE>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index used for the concentration field   
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_tracer_update_from_bottom(dt, tau, model_time)
    real,    intent(in) :: dt
    integer, intent(in) ::tau
    type(time_type),                intent(in) :: model_time

    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_update_from_bottom'

    !    if(do_generic_CFC)    call generic_CFC_update_from_bottom(tracer_list)!Nothing to do for CFC 

    if(do_generic_TOPAZ)  call generic_TOPAZ_update_from_bottom(tracer_list,dt, tau, model_time)

    if(do_generic_ERGOM)  call generic_ERGOM_update_from_bottom(tracer_list,dt, tau, model_time)
   
    if(do_generic_BLING)  call generic_BLING_update_from_bottom(tracer_list,dt, tau)

    if(do_generic_miniBLING)  call generic_miniBLING_update_from_bottom(tracer_list,dt, tau)

    if(do_generic_COBALT)  call generic_COBALT_update_from_bottom(tracer_list,dt, tau, model_time)

    return

  end subroutine generic_tracer_update_from_bottom


  ! <SUBROUTINE NAME="generic_tracer_vertdiff_G">
  !  <OVERVIEW>
  !   Vertically diffuse all generic tracers for GOLD ocean
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine use a tridiagonal solver to update the values
  !   of concentration field from vertical diffusion.
  !   The implicit vertdiff for these tracers should be disabled in the Ocean model.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_vertdiff_G(h_old, ea, eb, dt, Rho_0,tau)
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_tracer_vertdiff_G(h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau)
    real, dimension(:,:,:), intent(in) :: h_old, ea, eb
    real,                   intent(in) :: dt, kg_m2_to_H, m_to_H
    integer,                intent(in) :: tau
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next

    !nnz: Should I loop here or inside the sub g_tracer_vertdiff ?    
    !JGJ 2013/05/31  merged COBALT into siena_201303
    if(do_generic_CFC .or. do_generic_TOPAZ .or. do_generic_ERGOM .or. do_generic_BLING .or. do_generic_miniBLING .or. do_generic_COBALT) then

       g_tracer => tracer_list        
       !Go through the list of tracers 
       do  
          if(g_tracer_is_prog(g_tracer)) &
             call g_tracer_vertdiff_G(g_tracer,h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau)

          !traverse the linked list till hit NULL
          call g_tracer_get_next(g_tracer, g_tracer_next)
          if(.NOT. associated(g_tracer_next)) exit
          g_tracer=>g_tracer_next  

       enddo
    endif

  end subroutine generic_tracer_vertdiff_G

  ! <SUBROUTINE NAME="">
  !  <OVERVIEW>
  !   
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call 
  !  </TEMPLATE>
  !  <IN NAME="" TYPE="">
  !   
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_tracer_vertdiff_M(dh, dhw, diff_cbt, dt, Rho_0,tau)
    real, dimension(:,:,:), intent(in) :: dh, dhw, diff_cbt
    real,                   intent(in) :: dt,Rho_0
    integer,                intent(in) :: tau
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next

    !nnz: Should I loop here or inside the sub g_tracer_vertdiff ?    
    !JGJ 2013/05/31  merged COBALT into siena_201303
    if(do_generic_CFC .or. do_generic_TOPAZ .or. do_generic_ERGOM .or. do_generic_BLING .or. do_generic_miniBLING .or. do_generic_COBALT) then

       g_tracer => tracer_list        
       !Go through the list of tracers 
       do  
          if(g_tracer_is_prog(g_tracer)) &
               call g_tracer_vertdiff_M(g_tracer,dh, dhw, diff_cbt, dt, Rho_0,tau) 

          !traverse the linked list till hit NULL
          call g_tracer_get_next(g_tracer, g_tracer_next)
          if(.NOT. associated(g_tracer_next)) exit
          g_tracer=>g_tracer_next  

       enddo
    endif

  end subroutine generic_tracer_vertdiff_M

  ! <SUBROUTINE NAME="generic_tracer_coupler_set">
  !  <OVERVIEW>
  !   Set the coupler values for each generic tracer
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Calls the corresponding generic_X_set_boundary_values routine for each package X.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_coupler_set(IOB_struc, ST,SS,rho,ilb,jlb,tau)
  !  </TEMPLATE>
  !  <IN NAME="IOB_struc" TYPE="type(coupler_2d_bc_type)">
  !   Ice Ocean Boundary flux structure
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="ST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature   
  !  </IN>
  !  <IN NAME="SS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_tracer_coupler_set(IOB_struc, ST,SS,rho,ilb,jlb,tau)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
    integer, intent(in) :: ilb,jlb,tau
    real, dimension(ilb:,jlb:),  intent(in) :: ST,SS
    real, dimension(ilb:,jlb:,:,:), intent(in)              :: rho

    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_coupler_set'

    !Set coupler fluxes from tracer boundary values (%stf and %triver)for each tracer in the prog_tracer_list
    !User must identify these tracers (not all tracers in module need to set coupler)
    !User must provide the calculations for these boundary values.

    if(do_generic_CFC) &
         call generic_CFC_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)

    if(do_generic_TOPAZ) &
         call generic_TOPAZ_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)

    if(do_generic_ERGOM) &
         call generic_ERGOM_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)

    if(do_generic_BLING) &
         call generic_BLING_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)
    !
    if(do_generic_miniBLING) &
         call generic_miniBLING_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)

    if(do_generic_COBALT) &
         call generic_COBALT_set_boundary_values(tracer_list,ST,SS,rho,ilb,jlb,tau)

    !
    !Set coupler fluxes from tracer boundary values (%alpha and %csurf)
    !for each tracer in the tracer_list that has been marked by the user routine above
    !JGJ 2013/05/31  merged COBALT into siena_201303
    !
    if(do_generic_CFC .or. do_generic_TOPAZ .or. do_generic_ERGOM .or. do_generic_BLING .or. do_generic_miniBLING .or. do_generic_COBALT) &
       call g_tracer_coupler_set(tracer_list,IOB_struc)

  end subroutine generic_tracer_coupler_set

  ! <SUBROUTINE NAME="generic_tracer_coupler_zero">
  !  <OVERVIEW>
  !   Zero out the coupler values for each tracer
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Set the coupler arrays for ALL generic tracers to 0
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call 
  !  </TEMPLATE>
  !  <IN NAME="IOB_struc" TYPE="type(coupler_2d_bc_type)">
  !   Ice Ocean Boundary flux structure
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_tracer_coupler_zero(IOB_struc)
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
    !Generic tracer coupler values are not accumulative. No need to set them to zero
!    call g_tracer_coupler_set(tracer_list,IOB_struc,value=0.0)
  end subroutine generic_tracer_coupler_zero

  ! <SUBROUTINE NAME="generic_tracer_end">
  !  <OVERVIEW>
  !   End this module by calling the corresponding generic_X_end for each package X.
  !  </OVERVIEW>
  !  <TEMPLATE>
  !   call generic_tracer_end
  !  </TEMPLATE>
  ! </SUBROUTINE>
  subroutine generic_tracer_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_tracer_end'
    if(do_generic_CFC) call generic_CFC_end
    if(do_generic_TOPAZ)  call generic_TOPAZ_end
    if(do_generic_ERGOM)  call generic_ERGOM_end
    if(do_generic_BLING)  call generic_BLING_end
    if(do_generic_miniBLING)  call generic_miniBLING_end
    if(do_generic_COBALT)  call generic_COBALT_end

  end subroutine generic_tracer_end

  ! <SUBROUTINE NAME="generic_tracer_get_list">
  !  <OVERVIEW>
  !   Get a pointer to the head of the generic tracers list
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   The (linked) list of all generic tracers acts as the "state"
  !   of the generic tracers. This routine provides a way to get
  !   this state at any time.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_tracer_get_list(list)
  !  </TEMPLATE>
  !  <IN NAME="list" TYPE="type(g_tracer_type),    pointer">
  !   Pointer to head of the linked list
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_tracer_get_list(list)
    type(g_tracer_type),    pointer    :: list
    list => tracer_list
  end subroutine generic_tracer_get_list

  subroutine generic_tracer_get_diag_list(list)
    type(g_diag_type),    pointer    :: list
    list => diag_list
  end subroutine generic_tracer_get_diag_list

end module generic_tracer

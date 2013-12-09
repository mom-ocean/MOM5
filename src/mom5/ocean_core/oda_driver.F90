module oda_driver_mod

#ifdef ENABLE_ODA
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Top-level interface between MOM and ocean data assimilation modules
!</OVERVIEW>
!
!<DESCRIPTION>
! Top level module for ocean data assimilation.  Contains routines for 
! initialization, termination and update of ocean data assimilation increments.
!
! 1.) Get Observations for current model time (and domain)
!  based on data window (default=15 days).  For PE subdomains,
!  data are included if they lie within the region bounded by the
!  cell centers on the first x/y halo extent.
!
! 2.)  
!  IF EAKF (!!not yet implemented!!)
!           a.) ensemble state is distributed among NPES as follows:
!           Ens01:0:(NPE1-1), Ens02:NPE1:(NPE1+NPE2-1),...,
!           EnsN:sum(NPE1+...+NPE(N-1)):NPES
!           b.) Redistribute state to ensemble vector using global PEset.
!           (Ens01,...,EnsN):0:NPES
!           c.) Calculate increments using Ensemble information.
!           d.) Re-distribute increments to members and revert to
!               original PEset.
!  ELSE IF Var2d
!          Minimize quadratic CF based on specfied prior correlation maps
!          and optionally time-variant prior error variance.
!  ELSE ???           
!
!  ENDIF
!
! 3.) Apply increments to model state based on prior analysis.
!     Distribute increments uniformly over analysis period if "do_iau=T".
!
! 4.) Check for convective instability after applying increments.
!  
! 5.) Write Profile (observation space) first guess, analysis and mis-fits
!   
!========================================================================  
!</DESCRIPTION>
!
!  <DATA NAME="use_this_module" TYPE="logical">
!   use the oda module
!  </DATA>    
!<NAMELIST NAME="oda_nml">
!  <DATA NAME="assim_method" TYPE="character">
!   Options are: Var2d, EAKF and NO_ASSIM
!  </DATA>
!  <DATA NAME="assim_start_lat" TYPE="real">
!   Southern data mask  boundary
!  </DATA>  
!
!  <DATA NAME="assim_end_lat" TYPE="real">
!   Northern data mask  boundary
!  </DATA>  
!
!  <DATA NAME="nk_asm" TYPE="integer">
!   Bottom model level for data mask
!  </DATA>
!  
!  <DATA NAME="assim_interval" TYPE="integer">
!   Time between calls to oda (hours)
!  </DATA>  
!
!  <DATA NAME="do_iau" TYPE="logical">
!   Incremental analysis update (evenly distribute increments between
!   calls to ODA.)  
!  </DATA>
!
!  <DATA NAME="do_convect_adjust" TYPE="logical">
!   Adjust for gravitational instability in model after applying increments.
!  </DATA>  
!
!  <DATA NAME="max_profiles" TYPE="integer">
!   Size of allocation for profile data
!  </DATA>  
!
!  <DATA NAME="max_sfc_obs" TYPE="integer">
!   Size of allocation for surface observations
!  </DATA>  
!</NAMELIST>  

  use constants_mod,        only : epsln
  use ensemble_manager_mod, only : ensemble_manager_init, get_ensemble_id
  use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_pelist,get_ensemble_filter_pelist
  use mpp_domains_mod,      only : domain2d, mpp_update_domains
  use mpp_domains_mod,      only : mpp_define_domains, CYCLIC_GLOBAL_DOMAIN, BGRID_NE
  use mpp_domains_mod,      only : FOLD_NORTH_EDGE, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod,      only : null_domain2d, mpp_redistribute, mpp_broadcast_domain
  use mpp_mod,              only : input_nml_file, mpp_npes, stdout, stdlog, mpp_error, FATAL
  use mpp_mod,              only : mpp_pe, mpp_declare_pelist,mpp_get_current_pelist
  use mpp_mod,              only : mpp_set_current_pelist, mpp_root_pe,mpp_sum, mpp_sync
  use mpp_io_mod,           only : MPP_MULTI, MPP_SINGLE, mpp_io_exit
  use oda_types_mod
  use oda_core_mod
  use time_manager_mod,     only : time_type, decrement_time, increment_time
  use time_manager_mod,     only : get_time, get_date
  use write_ocean_data_mod

  use ocean_convect_mod,    only : convection
  use ocean_domains_mod,    only : set_ocean_domain, get_local_indices
  use ocean_grids_mod,      only : set_ocean_grid_size, set_ocean_hgrid_arrays
  use ocean_grids_mod,      only : set_ocean_vgrid_arrays
  use ocean_parameters_mod, only : rho_cp
  use ocean_topog_mod,      only : ocean_topog_init
  use ocean_types_mod
  use ocean_workspace_mod, only : wrk1

  
  IMPLICIT NONE
  
  private

  integer, parameter :: NO_ASSIM = 0, Var2d=1, EAKF=2, TEST_ENSEMBLE_ENGINE=3
  integer, parameter :: max_prog_tracers = 2

  integer :: max_profiles = 100000, max_sfc_obs = 100 ! for run-time allocation of observation arrays
  integer :: asm_method = NO_ASSIM
  character(len=8) :: assim_method = 'NO_ASSIM'

  logical :: do_iau = .true. ! incremental analysis update (recommended)
  logical :: do_convect_adjust=.true. ! for doing convective adjustment
                                      ! after applying increments (recommended)

  logical :: save_omf_snaps = .true.
  logical :: save_oma_snaps = .false.
  
  type(grid_type), target, save  :: T_grid
  type(domain2d), pointer, save  :: Asm_domain

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  integer :: ni, nj, nk
  integer :: num_prog_tracers

  real, dimension(:,:,:), allocatable, save :: corr, omf, oma, amf

  integer :: nk_asm
  integer :: xhalo, yhalo
  integer :: index_temp, index_salt

  integer :: id_omf = -1
  integer :: id_oma = -1
  integer :: id_amf = -1
  integer :: id_put_array = -1 
  integer :: id_got_array = -1 
  integer, parameter :: max_ensemble_size = 10
  integer :: pe_start(max_ensemble_size)=0, pe_end(max_ensemble_size)=0
  integer :: ensemble_size = 1
  logical, dimension(max_ensemble_size) :: ensemble_pe = .false.
  integer, allocatable, dimension(:,:) :: ensemble_pelist
  integer :: pe, npes, npes_per_member, ensemble_id
  logical :: sequential_filter = .false.
  
  integer :: isd_filt, ied_filt, jsd_filt, jed_filt
  integer :: isc_filt, iec_filt, jsc_filt, jec_filt
  
  real :: assim_start_lat = -90.0 , assim_end_lat = 90.0
  integer :: assim_interval = 24
  integer :: assim_tsteps=0
  
  character(len=64) :: fname_fg, fname_anal, fname_omf, fname_oma
  integer ::  iunit_fg, iunit_omf, iunit_anal, iunit_oma  
  real, save :: dtimer
  
  type(ocean_domain_type), pointer, save :: Dom
  type(ocean_grid_type), pointer, save :: Grd

  type(ocean_dist_type), target, save    :: First_guess, Analysis
  type(ocean_profile_type), save, allocatable :: Profiles(:)
  type(ocean_surface_type), save, allocatable :: Sfc_obs(:)

  type(ocean_profile_type), save, allocatable, dimension(:) :: model_obs, &
                                                               model_obsd
  
  type(grid_type), save :: Ens_grid
  type(domain2d), save :: Filter_domain, Send_domain(max_ensemble_size)
  type(field_type), save :: Tracer_ens(max_prog_tracers, max_ensemble_size)
  type(field_type), save :: U_ens(max_ensemble_size), V_ens(max_ensemble_size)
  type(field_type), save :: SSH_ens(max_ensemble_size),T_ens(max_ensemble_size)

  integer :: assim_layout(2)
  integer :: filter_halo_x = 1, filter_halo_y = 1
  integer :: ensemble_siz(6)
  integer :: total_npes_pm,ocean_npes_pm,atmos_npes_pm
  integer, allocatable :: ocean_pelist(:, :) , ocean_pelist_filter(:), original_pelist(:)
  real,    allocatable :: put_array(:,:,:), got_array(:,:,:) 

  logical :: use_this_module=.false.
  
  public :: init_oda, oda

  namelist /oda_nml/ use_this_module,assim_method, assim_start_lat,         &
            assim_end_lat, nk_asm, assim_interval, do_iau, do_convect_adjust,&
            max_profiles, max_sfc_obs, save_omf_snaps,                       &
            save_oma_snaps, assim_layout
  
contains  
  

!#######################################################################
! <SUBROUTINE NAME="init_oda">
!
! <DESCRIPTION>
! Initialize ODA core. Grid and Domain association.  
! </DESCRIPTION>
! <IN NAME="Domain" TYPE="ocean_domain_type"></IN>
! <IN NAME="Grid" TYPE="ocean_grid_type"></IN>  
! <IN NAME="Ocn_Time" TYPE="ocean_time_type"></IN>
! <IN NAME="T_prog" TYPE="ocean_tracer_type(:)"></IN>
  
  subroutine init_oda(Domain, Grid, Ocn_Time, T_prog)
    ! initialize First_guess and Analysis grid and domain information
    ! grids are identical to ocean_model grid.  Halo size may differ.

    use fms_mod, only : open_namelist_file,close_file,check_nml_error
    use mpp_mod, only : mpp_set_current_pelist
    use diag_manager_mod, only : diag_axis_init, register_diag_field
    use mpp_domains_mod, only : mpp_global_field, mpp_define_layout
    use ocean_grids_mod, only : set_ocean_grid_size, set_ocean_hgrid_arrays, &
                                set_ocean_vgrid_arrays
    use ocean_topog_mod, only : ocean_topog_init
    use ocean_domains_mod, only : set_ocean_domain, get_local_indices

    type(ocean_domain_type), target :: Domain
    type(ocean_grid_type), target :: Grid ! domain and grid information for ocean model
    type(ocean_time_type), target :: Ocn_Time
    type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

    type(time_type), pointer :: Time

    real :: adum,bdum,cdum
    
    integer :: len, xhalo, yhalo, halo(2), n, m, k, i, j, id_zt
    integer :: layout(2)
    integer :: ioun, io_status, ierr
    integer :: secs, days
    integer :: dt

    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 

    num_prog_tracers = size(T_prog)

    index_temp = -1; index_salt = -1

    do n=1, num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo

    Dom => Domain
    Grd => Grid
    
    Time => Ocn_Time%model_time

    call get_time(Ocn_time%time_step,secs,days)

    dt=secs+days*86400
    
        
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=oda_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'nml=oda_nml')
#else
    ioun = open_namelist_file()
    read(ioun,nml=oda_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'oda_nml')
    call close_file(ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, oda_nml)  
    write (stdlogunit, oda_nml)

    if (.not.use_this_module) return
    
    assim_interval=assim_interval*3600    
    write(stdoutunit,*) 'Assimilation interval= ',assim_interval,' seconds'
    write(stdoutunit,*) 'Model timestep= ',dt,' seconds'
    adum=assim_interval
    adum=real(adum/dt)
    if (adum - floor(adum) .ne. 0.0) then
        call mpp_error(FATAL,'assim_interval must be integer multiple &
             &of model timestep')
    else
        assim_tsteps=assim_interval/dt
    endif
        
    dtimer=1.0/dt
    
    if (nk_asm > Grid%nk) nk_asm = Grid%nk
    if (nk_asm < 2) call mpp_error(FATAL,&
      'need to specify number of vertical levels for assimilation')

    ni = Grid%ni
    nj = Grid%nj
    nk = Grid%nk

    call mpp_get_compute_domain(Domain%domain2d, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain%domain2d, isd, ied, jsd, jed)

    T_grid%ni = ni
    T_grid%nj = nj
    T_grid%nk = nk

    T_grid%Dom => Domain%domain2d
    
    allocate(T_grid%x(0:ni+1,0:nj+1), T_grid%y(0:ni+1,0:nj+1),&
             T_grid%mask(0:ni+1,0:nj+1,nk))
    allocate(T_grid%z(nk))

    T_grid%x=0.0
    T_grid%y=0.0
    T_grid%z=0.0
    T_grid%mask=0.0
    
    allocate(corr(isd:ied,jsd:jed,nk))       
    allocate(omf(isd:ied,jsd:jed,nk))       
    allocate(oma(isd:ied,jsd:jed,nk))
    allocate(amf(isd:ied,jsd:jed,nk))
    corr(:,:,:) = 0.0
    omf(:,:,:) = 0.0
    oma(:,:,:) = 0.0
    amf(:,:,:) = 0.0

! get global grid information from ocean_model
       
    call mpp_global_field(Domain%domain2d, Grid%tmask(:,:,:), &
         T_grid%mask(1:ni,1:nj,:))

    do j=1,nj
       T_Grid%x(1:ni,j) = Grid%grid_x_t(1:ni)
    enddo

    do i=1,ni
       T_Grid%y(i,1:nj) = Grid%grid_y_t(1:nj)
    enddo
    

    
    if (Grid%cyclic_x) then
        T_grid%x(0,:) = T_grid%x(ni,:)-360.
        T_grid%x(ni+1,:) = T_grid%x(1,:)+360.
    else
        T_grid%x(0,:) = T_grid%x(1,:)-1.e-10
        T_grid%x(ni+1,:) = T_grid%x(ni,:)+1.e-10
    endif

    T_grid%y(:,0) = T_grid%y(:,1)-1.e-3
    T_grid%y(:,nj+1) = T_grid%y(:,nj)+1.e-3

    T_grid%z = Grid%zt

    
    do j=0,nj+1
       do i=0,ni+1
          if (T_grid%y(i,j) < assim_start_lat .or. &
               T_grid%y(i,j) > assim_end_lat ) then
              T_grid%mask(i,j,:) = 0.0
          endif
          do k=nk_asm+1,nk
             T_grid%mask(i,j,k) = 0.0
          enddo
       enddo
    enddo

    T_grid%cyclic = Grid%cyclic_x

    if (assim_tsteps < 1) then
        call mpp_error(FATAL,'invalid assim interval')
    else
        write(stdoutunit,*) 'Assimilating data every ', assim_tsteps, ' timesteps'
    endif


    select case (assim_method)
    case ('NO_ASSIM')
       asm_method = NO_ASSIM
    case ('EAKF')
       asm_method = EAKF
    case ('Var2d')
       asm_method = Var2d
    case ('TEST_ENS')
       asm_method = TEST_ENSEMBLE_ENGINE
    case default
       call mpp_error(FATAL,'Wrong assim_method')
    end select

    select case (asm_method)

    case(TEST_ENSEMBLE_ENGINE)
       !This is just to test the functionality of the ensemble run
       !and can serve as an example of how to use the ensemble in 
       !real life data assimilations.

       !allocate ens array(Temp(;,;,;,Nens))
       !define domain for ens members(Send_domain(Nens))
       !define domain for filter (filter_domain)
       !broadcast_domain(Send_domain)
       !do i=1,Nens
       !   redist(Send_domain, TEMP, filter_domain, Temp(:,:,:,i)
       !enddo
       !do i=1,Nens
       !   redist(filter_domain, temp(:,:,:,i), Send_domain, Temp2)
       !enddo
       !is(Temp .eq. Temp2)?
       
       !
       !Get ensemble size and id
       !npes_pm (# of pe's per member) are assumed to be the same for all members
       
       ensemble_id = get_ensemble_id() 
       ensemble_siz = get_ensemble_size()   
       ensemble_size = ensemble_siz(1)      
       total_npes_pm = ensemble_siz(2)               
       ocean_npes_pm = ensemble_siz(3) 
       atmos_npes_pm = ensemble_siz(4) 

       !First save the current pelist
       allocate(original_pelist(ocean_npes_pm)) !can the argument be mpp_npes() ?
       call mpp_get_current_pelist(original_pelist)

       !Get ocean_pelist
       allocate(ocean_pelist(1:ensemble_size,1:ocean_npes_pm))   
       call get_ensemble_pelist(ocean_pelist,'ocean') 

       !Set pelist to ocean_pelist for this ensemble member
       call mpp_set_current_pelist(ocean_pelist(ensemble_id,:))

       !
       !Define "send" domains for individual ensemble members
       !

       call mpp_define_domains((/1,ni,1,nj/), Domain%layout, Send_domain(ensemble_id),&
            xflags=Domain%xflags, yflags=Domain%yflags, xhalo = Domain%xhalo, &
            yhalo = Domain%yhalo, name='send')

       allocate(put_array(isd:ied,jsd:jed,nk))
       allocate(got_array(isd:ied,jsd:jed,nk))

       id_put_array = register_diag_field('oda','oda_put' ,&
         Grid%tracer_axes(1:3),Time,'oda_put_array','none',&
         missing_value=-1.0e+10 )

       id_got_array = register_diag_field('oda','oda_got' ,&
         Grid%tracer_axes(1:3),Time,'oda_got_array','none',&
         missing_value=-1.0e+10 )
       
       !
       !Define the "filter" domain which is the domain for doing 
       !ensemble calculations (avereages, assimilations, ...)
       !
       !The "filter" domain should be defined on a pelist which
       !is the union of Ocean PE's for all ensemble members.
       !
       allocate(ocean_pelist_filter(ensemble_size * ocean_npes_pm))
       call get_ensemble_filter_pelist(ocean_pelist_filter,'ocean') 

       call mpp_set_current_pelist(ocean_pelist_filter) 
       
       call mpp_define_domains((/1,ni,1,nj/), assim_layout, Filter_domain, &
            xflags=Domain%xflags, yflags=Domain%yflags,&
            xhalo = filter_halo_x, yhalo = filter_halo_y, name='filter')

       call mpp_get_data_domain(Filter_domain, isd_filt,ied_filt,jsd_filt,jed_filt)
       call mpp_get_compute_domain(Filter_domain, isc_filt,iec_filt,jsc_filt,jec_filt)       

       do m = 1, ensemble_size
          allocate(T_ens(m)%data(isd_filt:ied_filt,jsd_filt:jed_filt,nk))
       enddo

       !return to original pelist
       call mpp_set_current_pelist(original_pelist) 

       return

    case (EAKF)

        call mpp_error(FATAL,'Not currently supported')
        
    case (Var2d)

        call mpp_error(FATAL,'Not currently supported')
        
!        xhalo=1;yhalo=1

!        Asm_domain => Domain%domain2d

!        allocate(First_guess%temp%ex(isd:ied,jsd:jed,nk))

!        First_guess%temp%ex = 0.0
!        First_guess%temp%grid => T_grid
       
!        ! assign tracer values on computational domain
!        ! assume for now that the ocean decomposition is identical
!        ! to assim decompositionn with the exception of halo sizes

!        allocate(Analysis%temp%ex(isd:ied,jsd:jed,nk))

!        Analysis%temp%grid => First_guess%temp%grid

!        call oda_core_init(Domain%domain2d, T_grid,localize=.true.)
!        call dr_oi_init(Domain, Grid, Time, nk_asm)

   case (NO_ASSIM)

       xhalo=1;yhalo=1

       Asm_domain => Domain%domain2d

       allocate(First_guess%temp%ex(isd:ied,jsd:jed,nk))

       First_guess%temp%ex = 0.0
       First_guess%temp%grid => T_grid
       
       ! assign tracer values on computational domain
       ! assume for now that the ocean decomposition is identical
       ! to assim decompositionn with the exception of halo sizes

       call oda_core_init(Domain%domain2d, T_grid,localize=.true.)
      
    end select

    id_omf = register_diag_field('oda','omf_temp',&
         Grid%tracer_axes(1:3),Time,'data-first guess misfit (temp)','deg_C',&
         missing_value=-20.,range=(/-20.,20./))

    
    id_oma = register_diag_field('oda','oma_temp',&
         Grid%tracer_axes(1:3),Time,'data-analysis misfit (temp)','deg_C',&
         missing_value=-20.,range=(/-20.,20./))
    
    id_amf = register_diag_field('oda','amf_temp',&
         Grid%tracer_axes(1:3),Time,'correction for temperature','W/m^2',&
         missing_value=-2.e3,range=(/-2.e3,2.e3/))

    allocate(Profiles(max_profiles))
    allocate(model_obs(max_profiles), model_obsd(max_profiles))    
    allocate(Sfc_obs(max_sfc_obs))

    call write_ocean_data_init()
    
    return

  end subroutine init_oda
! </SUBROUTINE> NAME="init_oda"


!#######################################################################
! <SUBROUTINE NAME="oda">
!
! <DESCRIPTION>
! Request ocean state increments
! </DESCRIPTION>
! <IN NAME="Ocn_Time" TYPE="ocean_type_type"></IN>
! <INOUT NAME="T_prog" TYPE="ocean_prog_tracer_type(:)"></INOUT>
! <INOUT NAME="Thickness" TYPE="ocean_thickness_type"></INOUT>
! <INOUT NAME="Dens" TYPE="ocean_density_type"></INOUT>
  
  subroutine oda(Ocn_Time, T_prog, Thickness, Dens)

!  use eakf_mod, only : ensemble_filter
    use diag_manager_mod, only : send_data

    type(ocean_time_type), target :: Ocn_Time
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    type(ocean_thickness_type), intent(inout), optional :: Thickness  
    type(ocean_density_type), intent(inout), optional :: Dens
    
    type(time_type) :: time

! snz adds a time-step index as the analysis time-step using time_step
    integer :: time_step_index = 0

    integer :: i,j, k, m, n, nprof, nsfc, itt, tau, taup1, &
         nprof_g
    real :: t, u, v, tfg, z1, z2
    logical :: assim_time, used

  integer :: stdoutunit 
  stdoutunit=stdout() 

    if (.not.use_this_module) return
    
    time = Ocn_time%model_time

    itt = Ocn_time%itt

    tau = Ocn_time%tau
           
    if (mod(itt-1,assim_tsteps) .eq. 0 .or. itt .eq. 0) then
        assim_time = .true.
        call get_obs(time, Profiles, Sfc_obs, nprof, nsfc)
        amf = 0.0;corr=0.0
    else
        assim_time = .false.
    endif
    
    

    select case (asm_method)

    case(TEST_ENSEMBLE_ENGINE)

       !First save the current pelist
       call mpp_get_current_pelist(original_pelist)

       !Set pelist to all filter PE's
       call get_ensemble_filter_pelist(ocean_pelist_filter,'ocean') 

       call mpp_set_current_pelist(ocean_pelist_filter) 
!       call mpp_set_current_pelist()

       put_array = T_prog(1)%field(:,:,:,1)
       got_array = -1.0

       do m = 1, ensemble_size
          call mpp_broadcast_domain(Send_domain(m))
       enddo
       
       !Distribute each ensemble member array to filter domain
       do m = 1, ensemble_size
          T_ens(m)%data(:,:,:) = 0.0
          call mpp_redistribute(Send_domain(m), put_array(:,:,:), &
               Filter_domain, T_ens(m)%data(:,:,:))          
       enddo

       !Here are filter calculations
       do m = 1, ensemble_size
          call mpp_update_domains(T_ens(m)%data(:,:,:), Filter_domain)
       enddo

       !Distribute the array back to ensemble members domains
       do m = 1, ensemble_size
          call mpp_redistribute(Filter_domain, T_ens(m)%data(:,:,:), &
               Send_domain(m), got_array(:,:,:))
       enddo

       !Since there was no real filter calculation the two arrays must match on each ensemble member
       do k=1,nk
          do j=jsc,jec; do i=isc,iec
             if (got_array(i,j,k)-put_array(i,j,k) .ne. 0.0) &
                  call mpp_error(FATAL,'Ensemble loopback test failed!')
          enddo;enddo
       enddo
       
       if (id_put_array > 0) used = send_data(id_put_array,put_array(isc:iec,jsc:jec,:),Time)
       if (id_got_array > 0) used = send_data(id_got_array,got_array(isc:iec,jsc:jec,:)&
                                                          -put_array(isc:iec,jsc:jec,:),Time)


       !return to original pelist
       call mpp_set_current_pelist(original_pelist) 
       

       return

    case (EAKF)

        call mpp_error(FATAL,'Not currently supported')

    case (Var2d)

       ! set first guess field to tau tracer

       First_guess%temp%ex(:,:,:) = &
            T_prog(index_temp)%field(:,:,:,tau)
       call mpp_update_domains( First_guess%temp%ex, Asm_domain)

       tau = Ocn_time%tau
       taup1 = Ocn_time%taup1       

       if (assim_time) then
           write(stdoutunit,*) 'itt = ',itt
           write(stdoutunit,*) 'Calculating analysis using observed data'
!           call dr_oi(First_guess,Profiles(1:nprof), Sfc_obs(1:nsfc), Analysis, Time)
           amf(:,:,:) = (Analysis%temp%ex(:,:,:) - &
                 First_guess%temp%ex(:,:,:))
           if (do_iau) then
! sub-divide increments equally between analysis intervals               
               corr = amf/assim_tsteps
           else
               corr = amf
           endif
       else
           if (do_iau) then
               write(stdoutunit,*) 'itt = ',itt
               write(stdoutunit,*) 'Applying analysis increment to ocean state'
           else
               corr = 0.0
           endif
       endif

       if (do_iau.or.assim_time) then
           wrk1(:,:,:) = T_prog(index_temp)%field(:,:,:,taup1)
           T_prog(index_temp)%field(:,:,:,taup1) = &
                wrk1(:,:,:) + corr(:,:,:)

           where (T_prog(index_temp)%field(:,:,:,taup1) .lt. -1.9)
               T_prog(index_temp)%field(:,:,:,taup1) = -1.9
           end where
           
           if (present(Dens).and.do_convect_adjust) then
               call convection(Ocn_Time, Thickness, T_prog, Dens)
               call mpp_update_domains(T_prog(index_temp)%field(:,:,:,taup1),Dom%domain2d)
               call mpp_update_domains(T_prog(index_salt)%field(:,:,:,taup1),Dom%domain2d)
           else
               write(stdoutunit,*) "WARNING:  Not doing convection in oda"
               call mpp_update_domains(T_prog(index_temp)%field(:,:,:,taup1),Dom%domain2d)               
           endif

           if (PRESENT(Thickness)) then
               amf(:,:,:) = dtimer*rho_cp*Thickness%dzt(:,:,:)*(T_prog(index_temp)%field(:,:,:,taup1) - wrk1(:,:,:))
           else
! NOTE: Warning!! this diagnostic is inaccurate
               do k=1,nk
                  amf(:,:,k) = dtimer*rho_cp*Grd%dzt(k)*(T_prog(index_temp)%field(:,:,k,taup1) - wrk1(:,:,k))
               enddo
           endif
       else
           amf = 0.0
       endif

   end select

   if (assim_time) then
       if (save_omf_snaps .or. save_oma_snaps) call open_profile_snaps(Time,nprof)

       if (asm_method == NO_ASSIM) then
           tau = Ocn_time%tau
           taup1 = Ocn_time%taup1                  
           First_guess%temp%ex(:,:,:) = &
                T_prog(index_temp)%field(:,:,:,taup1)
           call mpp_update_domains( First_guess%temp%ex, Asm_domain)
       endif
       

       if (nprof .gt. 0 .and. save_omf_snaps) then
   
   
           write(stdoutunit,*) 'itt = ',itt
           write(stdoutunit,*) 'Calculating first guess differences using observed data'
           call copy_obs(Profiles(1:nprof),model_obs(1:nprof)) ! copy profile data to separate array for holding model guess
           call copy_obs(Profiles(1:nprof),model_obsd(1:nprof))
           call forward_obs( model_obs(1:nprof), First_guess%temp%ex) ! put first guess in observation data (forward operator gets defined here)
           call assign_forward_model(model_obs(1:nprof),model_obsd(1:nprof))
           call diff_obs(model_obs(1:nprof),Profiles(1:nprof),model_obsd(1:nprof)) ! difference obs - model
           if (id_omf .gt. 0 ) then
               call backward_obs(model_obsd(1:nprof),omf)  ! D^T To
           endif

       endif

       if (save_omf_snaps) then
           call mpp_set_current_pelist((/mpp_pe()/))
       
           do i=1,nprof
              if (Profiles(i)%accepted) then
                  call write_profile(iunit_fg, model_obs(i))
                  call write_profile(iunit_omf, model_obsd(i))        
              endif
           enddo
           
           call mpp_set_current_pelist()
       endif
       
   endif



   
   if (asm_method /= NO_ASSIM .and. assim_time) then

       if (nprof .gt. 0 .and. save_oma_snaps) then
           write(stdoutunit,*) 'itt = ',itt
           write(stdoutunit,*) 'Calculating analysis differences using observed data'
           call forward_obs( model_obs(1:nprof), Analysis%temp%ex) ! put first guess in observation data (forward operator gets defined here)
           call diff_obs(model_obs(1:nprof),Profiles(1:nprof),model_obsd(1:nprof)) ! difference obs - model

       endif

       if (save_oma_snaps) then
           call mpp_set_current_pelist((/mpp_pe()/))
           
           do i=1,nprof
              if (Profiles(i)%accepted) then
                  
                  call write_profile(iunit_anal, model_obs(i))
                  call write_profile(iunit_oma,model_obsd(i))        
              endif
           enddo
  
           call mpp_set_current_pelist()
       endif

       
       if (id_oma .gt. 0 ) then
           call backward_obs(model_obsd(1:nprof),oma)  ! D^T To
       endif

   endif

       
           
   if (assim_time .and. (save_omf_snaps .or. save_oma_snaps)) call close_profile_snaps(nprof)
   
    if (id_omf > 0) used = send_data(id_omf,omf(isc:iec,jsc:jec,:),&
         &Time,rmask=Grd%tmask(isc-Dom%ioff:iec-Dom%ioff,&
         jsc-Dom%joff:jec-Dom%joff,:))

    if (id_oma > 0) used = send_data(id_oma,oma(isc:iec,jsc:jec,:),&
         &Time,rmask=Grd%tmask(isc-Dom%ioff:iec-Dom%ioff,&
         jsc-Dom%joff:jec-Dom%joff,:))

    if (id_amf > 0) used = send_data(id_amf,amf(isc:iec,jsc:jec,:),&
         &Time,rmask=Grd%tmask(isc-Dom%ioff:iec-Dom%ioff,&
         jsc-Dom%joff:jec-Dom%joff,:))
    

    return

  end subroutine oda
! </SUBROUTINE> NAME="init_oda"


  subroutine open_profile_snaps(Time,nprof)

    type(time_type), intent(in) :: Time
    integer, intent(in) :: nprof
    
    character(len=2) :: cmon,cday,chr
    character(len=4) :: cyr, cpe
    integer :: ipe
    integer :: year,month,day,hour,minute,second
  
  
    ipe = mpp_pe()

!    if (nprof.eq.0) return
    
    if (ipe .lt. 10) then
        write(cpe,'(a3,i1)') '000',ipe
    else if (ipe .lt. 100) then
        write(cpe,'(a2,i2)') '00',ipe
    else if (ipe .lt. 1000) then
        write(cpe,'(a1,i3)') '0',ipe
    else
        write(cpe,'(i4)') ipe
    endif
  
    call get_date(Time,year,month,day,hour,minute,second)
  
    if (year .lt. 10) then
        write(cyr,'(a3,i1)') '000',year
    else if (year .lt. 100) then
        write(cyr,'(a2,i2)') '00',year
    else if (year .lt. 1000) then
        write(cyr,'(a1,i3)') '0',year
    else
        write(cyr,'(i4)') year
    endif
  
    if (month .lt. 10) then
        write(cmon,'(a1,i1)') '0',month
    else
        write(cmon,'(i2)') month
    endif
    
    if (day .lt. 10) then
        write(cday,'(a1,i1)') '0',day
    else
        write(cday,'(i2)') day
    endif
    
    if (hour .lt. 10) then
        write(chr,'(a1,i1)') '0',hour
    else
        write(chr,'(i2)') hour
    endif
    

    if (save_omf_snaps) then
        write(fname_fg,'(a3,a4,a2,a2,a2,a3,a4,a3)') 'fg_',cyr,cmon,cday,chr,'.pe',cpe,'.nc'
        write(fname_omf,'(a4,a4,a2,a2,a2,a3,a4,a3)') 'omf_',cyr,cmon,cday,chr,'.pe',cpe,'.nc'
    endif

    if (save_oma_snaps) then
        write(fname_anal,'(a5,a4,a2,a2,a2,a3,a4,a3)') 'anal_',cyr,cmon,cday,chr,'.pe',cpe,'.nc'
        write(fname_oma,'(a4,a4,a2,a2,a2,a3,a4,a3)') 'oma_',cyr,cmon,cday,chr,'.pe',cpe,'.nc'  
    endif
    

    call mpp_set_current_pelist((/mpp_pe()/))

    if (nprof.ne.0 .and. save_omf_snaps) then
        iunit_fg = open_profile_file(fname_fg,1,&
             grid_lon=Grd%grid_x_t,grid_lat=Grd%grid_y_t,thread=MPP_SINGLE,fset=MPP_SINGLE)
        iunit_omf = open_profile_file(fname_omf,1,&
             grid_lon=Grd%grid_x_t,grid_lat=Grd%grid_y_t,thread=MPP_SINGLE,fset=MPP_SINGLE)
    endif

    if (nprof.ne.0 .and. save_oma_snaps) then    
        iunit_anal = open_profile_file(fname_anal,1,&
             grid_lon=Grd%grid_x_t,grid_lat=Grd%grid_y_t,thread=MPP_SINGLE,fset=MPP_SINGLE)
        iunit_oma = open_profile_file(fname_oma,1,&
             grid_lon=Grd%grid_x_t,grid_lat=Grd%grid_y_t,thread=MPP_SINGLE,fset=MPP_SINGLE)
    endif

    call mpp_set_current_pelist()

  end subroutine open_profile_snaps

  subroutine close_profile_snaps(nprof)

    integer, intent(in) :: nprof

    
    call mpp_set_current_pelist((/mpp_pe()/))

    if (nprof.ne.0 .and. save_omf_snaps) then
        call close_profile_file(iunit_fg)
        call close_profile_file(iunit_omf)
    endif

    if (nprof.ne.0 .and. save_oma_snaps) then    
        call close_profile_file(iunit_anal)
        call close_profile_file(iunit_oma)
    endif
    
    call mpp_set_current_pelist()

  end subroutine close_profile_snaps

#endif 
  
end module oda_driver_mod

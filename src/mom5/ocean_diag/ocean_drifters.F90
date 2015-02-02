module ocean_drifters_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Advect USER supplied drifters using the shared/drifters package.
!</OVERVIEW>
!
!<DESCRIPTION>
! Advect USER supplied drifters using the shared/drifters package.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_drifters_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to run this module. Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="output_interval" TYPE="integer">
!  Interval in timesteps between drifter writes
!  </DATA> 
!
!</NAMELIST>

  use drifters_mod
  use fms_mod
  use mpp_domains_mod
  use mpp_mod,          only : input_nml_file 
  use time_manager_mod, only : time_type, get_time

  use ocean_domains_mod,    only : get_local_indices
  use ocean_parameters_mod, only : rho0r
  use ocean_types_mod,      only : ocean_grid_type, ocean_domain_type
  use ocean_types_mod,      only : ocean_time_type, ocean_velocity_type
  use ocean_types_mod,      only : ocean_adv_vel_type, ocean_prog_tracer_type

  implicit none

  type(drifters_type),save :: drfts

  character(len=128) :: ermesg
  real :: xcmin, xcmax,ycmin,ycmax,xdmin,xdmax,ydmin,ydmax,xmin,xmax
  integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
  integer :: dt ! ocean time step in seconds
  real, dimension(:), allocatable :: x1d, y1d, z1d, dlon, dlat
  real, dimension(:,:,:), allocatable, save :: u, v, w, lon3d, lat3d, depth3d

  ! nml variables 
  logical :: use_this_module=.false.
  integer :: output_interval=1
  
  namelist /ocean_drifters_nml/ use_this_module, output_interval
  
contains

  subroutine ocean_drifters_init(Domain,Grid, Time, T_prog, Velocity, Adv_vel)

    type(ocean_domain_type) :: Domain
    type(ocean_grid_type)   :: Grid
    type(ocean_time_type)   :: Time
    type(ocean_prog_tracer_type), dimension(:) :: T_prog
    type(ocean_velocity_type) :: Velocity
    type(ocean_adv_vel_type) :: Adv_vel
    
    integer :: i,j,k, days, npes, ioun, ierr, io_status
    integer, dimension(:), allocatable :: pelist

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_drifters_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_drifters_nml')
#else
    ioun = open_namelist_file()
    read  (ioun, ocean_drifters_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_drifters_nml')
    call close_file (ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_drifters_nml)
    write (stdlogunit, ocean_drifters_nml)

    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk=Grid%nk

    if(.not. use_this_module) return 


    ! Determine local compute and data domain boundaries
    ! Note: This is only appropriate for spherical grids
    ! will not work for the tripolar grid region, for instance.

    allocate(x1d(isd:ied),y1d(jsd:jed),z1d(nk))
    allocate(u(isd:ied,jsd:jed,nk),v(isd:ied,jsd:jed,nk),w(isd:ied,jsd:jed,nk))
    allocate(lon3d(isd:ied,jsd:jed,nk),lat3d(isd:ied,jsd:jed,nk),depth3d(isd:ied,jsd:jed,nk))    
    allocate(dlon(isd:ied),dlat(jsd:jed))

    u=0.0;v=0.0;w=0.0
    
    do i=isc,iec
       x1d(i)=Grid%grid_x_t(i)
    enddo
    
    if (isc.eq.1) then
        x1d(isd)=Grid%grid_x_t(Grid%ni)-360.0
    else
        x1d(isd)=Grid%grid_x_t(isd)
    endif

    if (iec.eq.Grid%ni) then
        x1d(ied)=Grid%grid_x_t(1)+360.0
    else
        x1d(ied)=Grid%grid_x_t(ied)
    endif

    
    do j=jsc,jec
       y1d(j)=Grid%grid_y_t(j)
    enddo

    if (jsc.eq.1) then
        y1d(jsd) = Grid%grid_y_t(1)-1.e-5
    else
        y1d(jsd) = Grid%grid_y_t(jsd)
    endif

    if (jec.eq.Grid%nj) then
        y1d(jed) = Grid%grid_y_t(Grid%nj)+1.e-5
    else
        y1d(jed) = Grid%grid_y_t(jed)
    endif
    

    do k=1,nk
       z1d(k)=Grid%zt(k)
    enddo
    
    xcmin = x1d(isc)
    xcmax = x1d(iec)
    ycmin = y1d(jsc)
    ycmax = y1d(jec)    

    xdmin = x1d(isd)
    xdmax = x1d(ied)
    ydmin = y1d(jsd)
    ydmax = y1d(jed)    

    xmin = Grid%grid_x_t(Grid%ni)-360.;xmax = Grid%grid_x_t(1)+360.
    
    call get_time(Time%Time_Step,dt,days)
    
    call drifters_new(drfts, &
       & input_file ='INPUT/drifters_inp.nc'  ,  &
       & output_file='DRIFTERS/drifters_out.nc', &
       & ermesg=ermesg)

    npes=mpp_npes()
    allocate(pelist(npes))
    call mpp_get_pelist(Domain%domain2d,pelist)

    drfts%comm%pe_beg=pelist(1)
    drfts%comm%pe_end=pelist(npes)
    
    
  ! set the initial time and dt
    drfts%time = 0.0
    drfts%dt   = float(dt)


  ! set the PE domain boundaries. Xmin_comp/ymin_comp, xmax_comp/ymax_comp
  ! refer to the "compute" domain, which should cover densily the domain:
  ! xcmax[pe] = xcmin[pe_east]
  ! ycmax[pe] = ycmin[pe_north]
  ! Xmin_data/ymin_data, xmax_data/ymax_data refer to the "data" domain, which
  ! should be larger than the compute domain and therefore overlap:
  ! xdmax[pe] > xdmin[pe_east]
  ! ydmax[pe] > ydmin[pe_north]
  ! Particles in the overlap regions are tracked by several PEs. 
  call drifters_set_domain(drfts, &
       & xmin_comp=xcmin, xmax_comp=xcmax, &
       & ymin_comp=ycmin, ymax_comp=ycmax, &
       & xmin_data=xdmin, xmax_data=xdmax, &
       & ymin_data=ydmin, ymax_data=ydmax, &
       & xmin_glob=xmin , xmax_glob=xmax,  & ! this will set periodicity in x 
       & ermesg=ermesg)
  
  ! set neighboring PEs [domain2d is of type(domain2d)]
  call drifters_set_pe_neighbors(drfts, domain=Domain%domain2d, ermesg=ermesg)

  ! set the velocities axes. Each velocity can have different axes.
  call drifters_set_v_axes(drfts, component='u', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)

  call drifters_set_v_axes(drfts, component='v', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)

  call drifters_set_v_axes(drfts, component='w', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)  
  

  ! Distribute the drifters across PEs
  call drifters_distribute(drfts, ermesg)


  ! Push the drifters. u_comp, v_comp etc are provided by the host code

  do j=jsc,jec
     dlat(j) = 0.5*(y1d(j+1)-y1d(j-1))
  enddo
  dlat(jsd)=dlat(jsc)
  dlat(jed)=dlat(jec)

  do i=isc,iec
     dlon(i) = 0.5*(x1d(i+1)-x1d(i-1))
  enddo
  dlon(isd)=dlon(isc)
  dlon(ied)=dlon(iec)
  
  
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           u(i,j,k)=Velocity%u(i,j,k,1,Time%taup1)/Grid%dxt(i,j)*dlon(i)
           v(i,j,k)=Velocity%u(i,j,k,2,Time%taup1)/Grid%dyt(i,j)*dlat(j)
           w(i,j,k)=-1.0*Adv_vel%wrho_bt(i,j,k)*rho0r
           lon3d(i,j,k) = Grid%xt(i,j)
           lat3d(i,j,k) = Grid%yt(i,j)
           depth3d(i,j,k) = Grid%zt(k)                      
        enddo
     enddo
  enddo
  
  call drifters_push(drfts, u=u, v=v, w=w, ermesg=ermesg)


  ! check if RK4 integration is complete
  if(drfts%rk4_completed .and. mod(Time%itt0,output_interval) .eq. 0 ) then

      ! interpolate fields      
      call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, z=z1d, &
           &    data=lon3d, ermesg=ermesg)
      
      call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, z=z1d, &
           &    data=lat3d, ermesg=ermesg)             

      call drifters_set_field(drfts, index_field=3, x=x1d, y=y1d, z=z1d, &
           &    data=depth3d, ermesg=ermesg)
      
      call drifters_set_field(drfts, index_field=4, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(1)%field(:,:,:,Time%taup1), ermesg=ermesg)

      call drifters_set_field(drfts, index_field=5, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(2)%field(:,:,:,Time%taup1), ermesg=ermesg)      

      
      ! save data 
      call drifters_save(drfts, ermesg=ermesg)
      
  endif
  
  end subroutine ocean_drifters_init


  subroutine update_ocean_drifters(Velocity, Adv_vel, T_prog, Grid, Time)

    type(ocean_velocity_type) :: Velocity
    type(ocean_adv_vel_type) :: Adv_vel
    type(ocean_prog_tracer_type), dimension(:) :: T_prog
    type(ocean_grid_type) :: Grid
    type(ocean_time_type) :: Time
    
    integer :: i,j,k

    if(.not. use_this_module) return 

    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             u(i,j,k)=Velocity%u(i,j,k,1,Time%taup1)/Grid%dxt(i,j)*dlon(i)
             v(i,j,k)=Velocity%u(i,j,k,2,Time%taup1)/Grid%dyt(i,j)*dlat(j)
             w(i,j,k)=-1.0*Adv_vel%wrho_bt(i,j,k)*rho0r
          enddo
       enddo
    enddo


    ! push the drifters
    call drifters_push(drfts, u=u, v=v, w=w, ermesg=ermesg)

 
    ! check if RK4 integration is complete
    if(drfts%rk4_completed .and. mod(Time%itt0,output_interval) .eq. 0 ) then

        ! interpolate fields
        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, z=z1d, &
             &    data=lon3d, ermesg=ermesg)

        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, z=z1d, &
             &    data=lat3d, ermesg=ermesg)             

        call drifters_set_field(drfts, index_field=3, x=x1d, y=y1d, z=z1d, &
             &    data=depth3d, ermesg=ermesg)

        call drifters_set_field(drfts, index_field=4, x=x1d, y=y1d, z=z1d, &
	     &    data=T_prog(1)%field(:,:,:,Time%taup1), ermesg=ermesg)

        call drifters_set_field(drfts, index_field=5, x=x1d, y=y1d, z=z1d, &
             &    data=T_prog(2)%field(:,:,:,Time%taup1), ermesg=ermesg)

        call drifters_save(drfts, ermesg=ermesg)

    endif 
    
  end subroutine update_ocean_drifters


  subroutine ocean_drifters_end(Grid)
    type(ocean_grid_type) :: Grid

    if(.not. use_this_module) return 

    ! write restart file, optionally with lon/lat data coordinates (under development)
    call drifters_write_restart(drfts, filename='RESTART/drifters_inp.nc', &
!       & x1=Grid%grid_x_u, y1=Grid%grid_y_u, geolon1=Grid%xt, &
!       & x2=Grid%grid_x_u, y2=Grid%grid_y_u, geolat2=Grid%yt, &
         ermesg=ermesg)  

    ! destroy
    call drifters_del(drfts, ermesg=ermesg)

    deallocate(x1d,y1d,z1d,u,v,w,lon3d,lat3d,depth3d,dlon,dlat)


  end subroutine ocean_drifters_end


end module ocean_drifters_mod

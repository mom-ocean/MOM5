module land_data_mod

use mpp_mod           , only : mpp_get_current_pelist, mpp_pe
use constants_mod     , only : PI
use mpp_domains_mod   , only : domain2d, mpp_get_compute_domain, &
     mpp_define_layout, mpp_define_domains, mpp_define_io_domain, &
     mpp_get_current_ntile, mpp_get_tile_id, CYCLIC_GLOBAL_DOMAIN, &
     mpp_get_io_domain, mpp_get_pelist, mpp_get_layout
use mpp_mod,            only : mpp_chksum
use fms_mod           , only : write_version_number, mpp_npes, &
                               error_mesg, FATAL, stdout
use time_manager_mod  , only : time_type
use tracer_manager_mod, only : register_tracers, get_tracer_index, NO_TRACER
use field_manager_mod , only : MODEL_LAND
use grid_mod          , only : get_grid_ntiles, get_grid_size, get_grid_cell_vertices, &
     get_grid_cell_centers, get_grid_cell_area, get_grid_comp_area, &
     define_cube_mosaic
use land_tile_mod     , only : land_tile_type, land_tile_list_type, &
     land_tile_list_init, land_tile_list_end, nitems

implicit none
private

! ==== public interfaces =====================================================
public :: land_data_init
public :: land_data_end
public :: lnd            ! global data 

public :: atmos_land_boundary_type ! container for information passed from the 
                         ! atmosphere to land
public :: land_data_type ! container for information passed from land to 
                         ! the atmosphere
! both hold information on land grid (that is, after flux exchange translated 
! it from the atmosphere)
public land_data_type_chksum    ! routine to print checksums for land_data_type
public atm_lnd_bnd_type_chksum  ! routine to print checksums for atmos_land_boundary_type

public :: dealloc_land2cplr ! deallocates a land_data_type structure
public :: realloc_land2cplr ! allocates a land_data_type members for current 
                            ! number of tiles
public :: dealloc_cplr2land ! deallocates an atmos_land_boundary_type structure
public :: realloc_cplr2land ! allocates an atmos_land_boundary_type members 
                            ! for current number of tiles
! NOTE: realloc_* procedures can be called regardless of the current state
! of the argument data structures, since they deallocate data first.

public :: land_state_type
! ==== end of public interfaces ==============================================

! ---- module constants ------------------------------------------------------
character(len=*), parameter :: &
     module_name = 'land_data_mod', &
     version     = '$Id: land_data.F90,v 20.0 2013/12/13 23:29:24 fms Exp $', &
     tagname     = '$Name: tikal $'

! init_value is used to fill most of the allocated boundary condition arrays.
! It is supposed to be double-precision signaling NaN, to trigger a trap when
! the program is compiled with trapping non-initialized values.  
! See http://ftp.uniovi.es/~antonio/uned/ieee754/IEEE-754references.html
real, parameter :: init_value = Z'FFF0000000000001'

! ---- types -----------------------------------------------------------------
type :: atmos_land_boundary_type
   ! data passed from the coupler to the surface
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperature of precipitation, degK
   ! components of downward shortwave flux, W/m2  
        sw_flux_down_vis_dir   => NULL(), & ! visible direct 
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
   ! derivatives of the fluxes
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature 
   !
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value, 
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   ! data passed from the surface to the coupler
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   ! NOTE that in contrast to most of the other fields in this structure, the discharges
   ! hold data per-gridcell, rather than per-tile basis. This, and the order of updates,
   ! have implications for the data reallocation procedure.
   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   type(domain2d) :: domain  ! our computation domain
   logical, pointer :: maskmap(:,:) 
   integer, pointer, dimension(:) :: pelist
end type land_data_type


! land_state_type combines the general information about state of the land model:
! domain, coordinates, time steps, etc. There is only one variable of this type,
! and it is public in this module.
type :: land_state_type
   integer        :: is,ie,js,je ! compute domain boundaries
   integer        :: nlon,nlat   ! size of global grid
   integer        :: ntprog      ! number of prognostic tracers
   integer        :: isphum      ! index of specific humidity in tracer table
   integer        :: ico2        ! index of carbon dioxide in tracer table
   type(time_type):: dt_fast     ! fast (physical) time step
   type(time_type):: dt_slow     ! slow time step
   type(time_type):: time        ! current time

   real, pointer  :: lon (:,:), lat (:,:) ! domain grid center coordinates, radian
   real, pointer  :: lonb(:,:), latb(:,:) ! domain grid vertices, radian
   real, pointer  :: area(:,:)  ! land area per grid cell, m2
   real, pointer  :: cellarea(:,:)  ! grid cell area, m2
   real, pointer  :: coord_glon(:), coord_glonb(:) ! longitudes for use in diag axis and such, degrees East
   real, pointer  :: coord_glat(:), coord_glatb(:) ! latitudes for use in diag axis and such, degrees North

   ! map of tiles
   type(land_tile_list_type), pointer :: tile_map(:,:)
   
   type(domain2d) :: domain ! our domain -- should be the last since it simplifies
                            ! debugging in totalview
   integer :: nfaces ! number of mosaic faces
   integer :: face  ! the current mosaic face
   integer, allocatable :: pelist(:) ! list of processors that run land model
   integer, allocatable :: io_pelist(:) ! list of processors in our io_domain
   ! if io_domain was not defined, then there is just one element in this
   ! array, and it's equal to current PE
   integer :: io_id     ! suffix in the distributed files.
   logical :: append_io_id ! if FALSE, io_id is not appended to the file names
                           ! (for the case io_layout = 1,1)
end type land_state_type

! ---- public module variables -----------------------------------------------
type(land_state_type),save :: lnd


! ---- private module variables ----------------------------------------------
logical :: module_is_initialized =.FALSE.


#define __DEALLOC__(x) if (associated(x)) deallocate(x)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



! ============================================================================
subroutine land_data_init(layout, io_layout, time, dt_fast, dt_slow)
  integer, intent(inout) :: layout(2) ! layout of our domains
  integer, intent(inout) :: io_layout(2) ! layout for land model io
  type(time_type), intent(in) :: &
       time,    & ! current model time
       dt_fast, & ! fast (physical) time step
       dt_slow    ! slow time step

  ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid 
  integer :: ntracers, ndiag ! non-optional output from register_tracers
  integer, allocatable :: tile_ids(:) ! mosaic tile IDs for the current PE
  integer :: i,j
  type(domain2d), pointer :: io_domain ! our io_domain
  integer :: n_io_pes(2) ! number of PEs in our io_domain along x and y
  integer :: io_id(1)

  ! write the version and tag name to the logfile
  call write_version_number(version, tagname)

  ! define the processor layout information according to the global grid size 
  call get_grid_ntiles('LND',ntiles)
  lnd%nfaces = ntiles
  call get_grid_size('LND',1,nlon,nlat)
  ! set the size of global grid
  lnd%nlon = nlon; lnd%nlat = nlat
  if( layout(1)==0 .AND. layout(2)==0 ) &
       call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes()/ntiles, layout )
  if( layout(1)/=0 .AND. layout(2)==0 )layout(2) = mpp_npes()/(layout(1)*ntiles)
  if( layout(1)==0 .AND. layout(2)/=0 )layout(1) = mpp_npes()/(layout(2)*ntiles)

  ! define land model domain
  if (ntiles==1) then
     call mpp_define_domains ((/1,nlon, 1, nlat/), layout, lnd%domain, xhalo=1, yhalo=1,&
          xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL')
  else
     call define_cube_mosaic ('LND', lnd%domain, layout, halo=1)
  endif

  ! define io domain
  call mpp_define_io_domain(lnd%domain, io_layout)

  ! set up list of processors for collective io: only the first processor in this
  ! list actually writes data, the rest just send the data to it.
  io_domain=>mpp_get_io_domain(lnd%domain)
  if (associated(io_domain)) then
     call mpp_get_layout(io_domain,n_io_pes)
     allocate(lnd%io_pelist(n_io_pes(1)*n_io_pes(2)))
     call mpp_get_pelist(io_domain,lnd%io_pelist)
     io_id = mpp_get_tile_id(io_domain)
     lnd%io_id = io_id(1)
  else
     allocate(lnd%io_pelist(1))
     lnd%io_pelist(1) = mpp_pe()
     lnd%io_id        = mpp_pe()
  endif
  lnd%append_io_id = (io_layout(1)/=1.or.io_layout(2)/=1)
     

  ! get the domain information
  call mpp_get_compute_domain(lnd%domain, lnd%is,lnd%ie,lnd%js,lnd%je)

  ! get the mosaic tile number for this processor: this assumes that there is only one
  ! mosaic tile per PE.
  allocate(tile_ids(mpp_get_current_ntile(lnd%domain)))
  tile_ids = mpp_get_tile_id(lnd%domain)
  lnd%face = tile_ids(1)
  deallocate(tile_ids)

  allocate(lnd%tile_map (lnd%is:lnd%ie, lnd%js:lnd%je))

  allocate(lnd%lonb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%latb    (lnd%is:lnd%ie+1, lnd%js:lnd%je+1))
  allocate(lnd%lon     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%lat     (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%area    (lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%cellarea(lnd%is:lnd%ie,   lnd%js:lnd%je))
  allocate(lnd%coord_glon(nlon), lnd%coord_glonb(nlon+1))
  allocate(lnd%coord_glat(nlat), lnd%coord_glatb(nlat+1))

  ! initialize coordinates
  call get_grid_cell_vertices('LND',lnd%face,lnd%coord_glonb,lnd%coord_glatb)
  call get_grid_cell_centers ('LND',lnd%face,lnd%coord_glon, lnd%coord_glat)
  call get_grid_cell_area    ('LND',lnd%face,lnd%cellarea, domain=lnd%domain)
  call get_grid_comp_area    ('LND',lnd%face,lnd%area,     domain=lnd%domain)
  
  ! set local coordinates arrays -- temporary, till such time as the global arrays
  ! are not necessary
  call get_grid_cell_vertices('LND',lnd%face,lnd%lonb,lnd%latb, domain=lnd%domain)
  call get_grid_cell_centers ('LND',lnd%face,lnd%lon, lnd%lat, domain=lnd%domain)
  ! convert coordinates to radian; note that 1D versions stay in degrees
  lnd%lonb = lnd%lonb*pi/180.0 ; lnd%lon = lnd%lon*pi/180.0
  lnd%latb = lnd%latb*pi/180.0 ; lnd%lat = lnd%lat*pi/180.0
  
  ! initialize land tile map
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_init(lnd%tile_map(i,j))
  enddo
  enddo

  ! initialize land model tracers, if necessary
  ! register land model tracers and find specific humidity
  call register_tracers ( MODEL_LAND, ntracers, lnd%ntprog, ndiag )
  lnd%isphum = get_tracer_index ( MODEL_LAND, 'sphum' )
  if (lnd%isphum==NO_TRACER) then
     call error_mesg('land_model_init','no required "sphum" tracer',FATAL)
  endif
  lnd%ico2 = get_tracer_index ( MODEL_LAND, 'co2' )
  ! NB: co2 might be absent, in this case ico2 == NO_TRACER

  ! initialize model's time-related parameters
  lnd%time    = time
  lnd%dt_fast = dt_fast
  lnd%dt_slow = dt_slow

  ! initialize the land model processor list
  allocate(lnd%pelist(0:mpp_npes()-1))
  call mpp_get_current_pelist(lnd%pelist)

end subroutine land_data_init

! ============================================================================
subroutine land_data_end()

  integer :: i,j
  
  module_is_initialized = .FALSE.

  ! deallocate land tile map here. 
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     call land_tile_list_end(lnd%tile_map(i,j))
  enddo
  enddo

  ! deallocate grid data
  deallocate(lnd%lonb, lnd%latb, lnd%lon, lnd%lat,&
       lnd%area, lnd%cellarea,&
       lnd%coord_glonb, lnd%coord_glon, &
       lnd%coord_glatb, lnd%coord_glat, &       
       lnd%tile_map,&
       lnd%pelist, lnd%io_pelist)

end subroutine land_data_end


! ============================================================================
! allocates boundary data for land domain and current number of tiles
subroutine realloc_land2cplr ( bnd )
  type(land_data_type), intent(inout) :: bnd     ! data to allocate

  ! ---- local vars
  integer :: n_tiles

  call dealloc_land2cplr(bnd, dealloc_discharges=.FALSE.)

  bnd%domain = lnd%domain
  n_tiles = max_n_tiles()


  ! allocate data according to the domain boundaries
  allocate( bnd%mask(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  allocate( bnd%tile_size(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_surf(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%t_ca(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%tr(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles,lnd%ntprog) )
  allocate( bnd%albedo(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dir(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%albedo_nir_dif(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_mom(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_heat(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )
  allocate( bnd%rough_scale(lnd%is:lnd%ie,lnd%js:lnd%je,n_tiles) )

  bnd%mask              = .FALSE.
  bnd%tile_size         = init_value
  bnd%t_surf            = init_value
  bnd%t_ca              = init_value
  bnd%tr                = init_value
  bnd%albedo            = init_value
  bnd%albedo_vis_dir    = init_value
  bnd%albedo_nir_dir    = init_value
  bnd%albedo_vis_dif    = init_value
  bnd%albedo_nir_dif    = init_value
  bnd%rough_mom         = init_value
  bnd%rough_heat        = init_value
  bnd%rough_scale       = init_value

  ! in contrast to the rest of the land boundary condition fields, discharges 
  ! are specified per grid cell, not per tile; therefore they should not be 
  ! re-allocated when the number of tiles changes. In fact, they must not be
  ! changed at all here because their values are assigned in update_land_model_fast,
  ! not in update_land_bc_*, and therefore would be lost if re-allocated.
  if (.not.associated(bnd%discharge)) then
     allocate( bnd%discharge          (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_heat     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow     (lnd%is:lnd%ie,lnd%js:lnd%je) )
     allocate( bnd%discharge_snow_heat(lnd%is:lnd%ie,lnd%js:lnd%je) )

     ! discharge and discaharge_snow must be, in contrast to the rest of the boundary
     ! values, filled with zeroes. The reason is because not all of the usable elements
     ! are updated by the land model (only coastal points are).
     bnd%discharge           = 0.0
     bnd%discharge_heat      = 0.0
     bnd%discharge_snow      = 0.0
     bnd%discharge_snow_heat = 0.0
  endif
end subroutine realloc_land2cplr


! ============================================================================
! deallocates boundary data memory
! NOTE that the discharges should be deallocated only at the final clean-up
! stage; during the model run they should be preserved unchanged even when
! other fields are reallocated.
subroutine dealloc_land2cplr ( bnd, dealloc_discharges )
  type(land_data_type), intent(inout) :: bnd  ! data to de-allocate
  logical, intent(in) :: dealloc_discharges

  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%tile_size )
  __DEALLOC__( bnd%t_surf )
  __DEALLOC__( bnd%t_ca )
  __DEALLOC__( bnd%tr )
  __DEALLOC__( bnd%albedo )
  __DEALLOC__( bnd%albedo_vis_dir )
  __DEALLOC__( bnd%albedo_nir_dir )
  __DEALLOC__( bnd%albedo_vis_dif )
  __DEALLOC__( bnd%albedo_nir_dif )
  __DEALLOC__( bnd%rough_mom )
  __DEALLOC__( bnd%rough_heat )
  __DEALLOC__( bnd%rough_scale )
  __DEALLOC__( bnd%mask )

  if (dealloc_discharges) then
     __DEALLOC__( bnd%discharge           )
     __DEALLOC__( bnd%discharge_heat      )
     __DEALLOC__( bnd%discharge_snow      )
     __DEALLOC__( bnd%discharge_snow_heat )
  end if

end subroutine dealloc_land2cplr


! ============================================================================
! allocates boundary data for land domain and current number of tiles;
! initializes data for data override.
! NOTE: previously the body of the procedure was in the flux_exchange_init,
! currently it is called from land_model_init
subroutine realloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  ! ---- local vars
  integer :: kd

  call dealloc_cplr2land(bnd)

  ! allocate data according to the domain boundaries
  kd = max_n_tiles()

  allocate( bnd%t_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%lprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%fprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%tprec(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%dhdq(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%drdt(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%p_surf(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%tr_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )
  allocate( bnd%dfdtr(lnd%is:lnd%ie,lnd%js:lnd%je,kd,lnd%ntprog) )

  allocate( bnd%lwdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%swdn_flux(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dir(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%sw_flux_down_total_dif(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_t(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%cd_m(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%bstar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%ustar(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%wind(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )
  allocate( bnd%z_bot(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )

  allocate( bnd%drag_q(lnd%is:lnd%ie,lnd%js:lnd%je,kd) )

  bnd%t_flux                 = init_value
  bnd%lw_flux                = init_value
  bnd%sw_flux                = init_value
  bnd%lprec                  = init_value
  bnd%fprec                  = init_value
  bnd%tprec                  = init_value
  bnd%dhdt                   = init_value
  bnd%dhdq                   = init_value
  bnd%drdt                   = init_value
  bnd%p_surf                 = init_value
  bnd%tr_flux                = init_value
  bnd%dfdtr                  = init_value

  bnd%lwdn_flux              = init_value
  bnd%swdn_flux              = init_value
  bnd%sw_flux_down_vis_dir   = init_value
  bnd%sw_flux_down_total_dir = init_value
  bnd%sw_flux_down_vis_dif   = init_value
  bnd%sw_flux_down_total_dif = init_value
  bnd%cd_t                   = init_value
  bnd%cd_m                   = init_value
  bnd%bstar                  = init_value
  bnd%ustar                  = init_value
  bnd%wind                   = init_value
  bnd%z_bot                  = init_value

  bnd%drag_q                 = init_value

end subroutine realloc_cplr2land


! ============================================================================
subroutine dealloc_cplr2land( bnd )
  type(atmos_land_boundary_type), intent(inout) :: bnd

  __DEALLOC__( bnd%t_flux )
  __DEALLOC__( bnd%lw_flux )
  __DEALLOC__( bnd%sw_flux )
  __DEALLOC__( bnd%lprec )
  __DEALLOC__( bnd%fprec )
  __DEALLOC__( bnd%dhdt )
  __DEALLOC__( bnd%dhdq )
  __DEALLOC__( bnd%drdt )
  __DEALLOC__( bnd%p_surf )
  __DEALLOC__( bnd%lwdn_flux )
  __DEALLOC__( bnd%swdn_flux )
  __DEALLOC__( bnd%sw_flux_down_vis_dir )
  __DEALLOC__( bnd%sw_flux_down_total_dir )
  __DEALLOC__( bnd%sw_flux_down_vis_dif )
  __DEALLOC__( bnd%sw_flux_down_total_dif )
  __DEALLOC__( bnd%cd_t )
  __DEALLOC__( bnd%cd_m )
  __DEALLOC__( bnd%bstar )
  __DEALLOC__( bnd%ustar )
  __DEALLOC__( bnd%wind )
  __DEALLOC__( bnd%z_bot )
  __DEALLOC__( bnd%tr_flux )
  __DEALLOC__( bnd%dfdtr )

end subroutine dealloc_cplr2land

! ============================================================================
! get max number of tiles in the domain
function max_n_tiles() result(n)
  integer :: n
  integer :: i,j

  n=1
  do j=lnd%js,lnd%je
  do i=lnd%is,lnd%ie
     n=max(n, nitems(lnd%tile_map(i,j)))
  enddo
  enddo

end function 


! ===========================================================================
!  Prints checksums of the various fields in the atmos_land_boundary_type.
subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)
    character(len=*), intent(in) :: id  ! Label to differentiate where this 
                      ! routine is being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which 
                      ! timestep this routine is being called for.
    type(atmos_land_boundary_type), intent(in) :: albt
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_land_boundary_type):: ', id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%swdn_flux               ', mpp_chksum( albt%swdn_flux)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%dhdq                  ', mpp_chksum( albt%dhdq)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)
    do n = 1,size(albt%tr_flux,4)
    write(outunit,100) 'albt%tr_flux               ', mpp_chksum( albt%tr_flux(:,:,:,n))
    enddo
    do n = 1,size(albt%dfdtr,4)
    write(outunit,100) 'albt%dfdtr                 ', mpp_chksum( albt%dfdtr(:,:,:,n))
    enddo

100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum


! ===========================================================================
!  Prints checksums of the various fields in the land_data_type.
subroutine land_data_type_chksum(id, timestep, land)
    character(len=*), intent(in) :: id ! Label to differentiate where this 
        ! routine in being called from.
    integer         , intent(in) :: timestep ! An integer to indicate which
        ! timestep this routine is being called for.
    type(land_data_type), intent(in) :: land
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(land_data_type):: ', id, timestep
    write(outunit,100) 'land%tile_size         ',mpp_chksum(land%tile_size)
    write(outunit,100) 'land%t_surf            ',mpp_chksum(land%t_surf)
    write(outunit,100) 'land%t_ca              ',mpp_chksum(land%t_ca)
    write(outunit,100) 'land%albedo            ',mpp_chksum(land%albedo)
    write(outunit,100) 'land%albedo_vis_dir    ',mpp_chksum(land%albedo_vis_dir)
    write(outunit,100) 'land%albedo_nir_dir    ',mpp_chksum(land%albedo_nir_dir)
    write(outunit,100) 'land%albedo_vis_dif    ',mpp_chksum(land%albedo_vis_dif)
    write(outunit,100) 'land%albedo_nir_dif    ',mpp_chksum(land%albedo_nir_dif)
    write(outunit,100) 'land%rough_mom         ',mpp_chksum(land%rough_mom)
    write(outunit,100) 'land%rough_heat        ',mpp_chksum(land%rough_heat)
    write(outunit,100) 'land%rough_scale       ',mpp_chksum(land%rough_scale)

    do n = 1, size(land%tr,4)
    write(outunit,100) 'land%tr                ',mpp_chksum(land%tr(:,:,:,n))
    enddo
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)
    write(outunit,100) 'land%discharge_heat    ',mpp_chksum(land%discharge_heat)


100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine land_data_type_chksum

end module land_data_mod

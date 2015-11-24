#include <fms_platform.h>

! ============================================================================
! an "umbrella" module for all land-related modules
! ============================================================================

module land_model_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Christopher Milly
! </CONTACT> 

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Elena Shevliakova
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Sergey Malyshev
! </REVIEWER>

! <OVERVIEW>
!   This module contains calls to the sub-models to calculate quantities
!   on the fast and slow time scales and update the boundary conditions and
!   tiling structure where necessary.
! </OVERVIEW>

! <DESCRIPTION>
!   The sub-models (soil and vegetation) are first called in top-bottom order
!   to evaluate the derivatives. Then they are called in bottom-top order to
!   finish the implicit calculations. On the upward pass, they have a chance
!   to update boundary values returned to the atmosphere. These calls are done
!   on the fast time scale to calculate the fluxes.
!
!   On the slow time scale, the land is updated and the boundary conditions
!   may be updated. The land boundary data is updated on the fast time scale 
!   without updating the tiling structure. The land boundary data for the
!   atmosphere is updated on the slow time scale and changes the tiling
!   structure, if necessary, as well as the albedo, drag coefficient and such.
!
!   The diagnostic horizontal axes for the land grid is initialized. All 
!   sub-models use them instead of creating their own. Also, NetCDF library 
!   messages are printed out, including file names and line numbers.
! </DESCRIPTION>


  ! things from other modules which are used in the interface part
  ! of the land module
  use time_manager_mod, only: time_type

  use mpp_domains_mod,  only: domain1D, domain2d, mpp_get_layout, mpp_define_layout, &
                              mpp_define_domains, mpp_get_compute_domain,            &
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domains,         &
                              mpp_get_domain_components, mpp_get_pelist,             &
                              mpp_define_mosaic

  use fms_mod,          only: write_version_number, error_mesg, FATAL, NOTE, &
                              mpp_pe, mpp_npes, mpp_root_pe, &
                              mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                              clock_flag_default, CLOCK_COMPONENT, file_exist, &
                              field_size, read_data, field_exist, get_mosaic_tile_grid, &
                              set_domain, nullify_domain

  use vegetation_mod,   only: vegetation_type, vegetation_init, vegetation_end,       &
                              update_vegetation_fast_down, update_vegetation_fast_up, &
                              update_vegetation_slow, update_vegetation_bnd_fast,     &
                              update_vegetation_bnd_slow, vegetation_stock_pe

  use soil_mod,         only: soil_type, soil_init, soil_end, update_soil_fast, &
                              update_soil_slow, update_soil_bnd_fast,           &
                              update_soil_bnd_slow, soil_stock_pe

  use land_types_mod,   only: land_data_type, allocate_boundary_data, &
                              atmos_land_boundary_type, deallocate_boundary_data, &
                              land_data_type_chksum, atm_lnd_bnd_type_chksum 

  use numerics_mod,     only: is_latlon

  use constants_mod,   only:  radius, hlf, hlv, tfreeze, pi 

  use diag_manager_mod, only: diag_axis_init

  use field_manager_mod,  only: MODEL_LAND
  use tracer_manager_mod, only: register_tracers, get_tracer_index, NO_TRACER

  use mosaic_mod,         only: get_mosaic_ntiles, calc_mosaic_grid_area, &
                                get_mosaic_xgrid_size, get_mosaic_xgrid

implicit none
private


! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! dummy routines
public update_land_model_fast   ! fast time-scale update of the land
public update_land_model_slow   ! slow time-scale update of the land

public Lnd_stock_pe             ! calculate and return total amount of requested quantitiy per PE
public land_data_type_chksum, atm_lnd_bnd_type_chksum
! ==== end of public interfaces ==============================================

! ==== public data type =====================================================
public land_type, land_data_type, atmos_land_boundary_type

! ==== generic interfaces ====================================================
! <INTERFACE NAME="land_model_init">
!   <OVERVIEW>
!     Initializes the state of the land model.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Initializes the land state in three ways. 1) The grid description file may
!      be used as an input. The land grid boundaries and area of land are read from
!      a grid description file. 2) A logical land mask may be used to calculate the
!      area of the grid cells where the mask is true and calls init_land_with_area.
!      Naturally, in this case there are no partly covered land cells -- every cell
!      is either land or ocean. 3) The land state may be initialized from the land area
!      for each of the grid points.
!   </DESCRIPTION>
interface land_model_init
   module procedure init_land_with_area
   module procedure init_land_with_xgrid
end interface
! </INTERFACE>


! <TYPE NAME="land_type">
!   <DESCRIPTION>
!     The type describing the state of the land model. It is private to this
!     module and is basically a remnant of the previous design. There is only one
!     variable of this type, called "theLand", which is used in the model to
!     hold the data.
!     Strictly speaking, this type is not necessary, but there is no harm in
!     having it and it is possible to image a situation where a single entity
!     describing the land-surface model state would be useful.
!   </DESCRIPTION>
type land_type

   private  ! it's an opaque type: all data members are private
   
!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     Domain of calculations
!   </DATA>
!   <DATA NAME="is" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="ie" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="js" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="je" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="n_tiles" TYPE="integer">
!     Domain bound
!   </DATA>
   type(domain2d)         :: domain     ! domain of calculations
   integer :: is,ie,js,je,n_tiles       ! domain bounds, for convenience

!   <DATA NAME="soil" TYPE="soil_type">
!     Soil component data
!   </DATA>
!   <DATA NAME="vegetation" TYPE="vegetation_type">
!     Vegetation component data
!   </DATA>
   type(soil_type)        :: soil       ! soil component data
   type(vegetation_type)  :: vegetation ! vegetation component data

   ! the values below is just a quick fix to start running Land with original
   ! exchange grid, since it assumes something like that to be present. It also 
   ! assumes the number of tiles to be fixed.
   !   Of course, in general case there is no such thing as "latitude/longitude"
   ! boundaries of the cell, since cell boundaries do not have to be parallel
   ! to the lat/lon coordinate axes

!   <DATA NAME="blon" TYPE="real, pointer" DIM="2">
!     Longitude domain corners
!   </DATA>
!   <DATA NAME="blat" TYPE="real, pointer" DIM="2">
!     Latitude domain corners
!   </DATA>
!   <DATA NAME="mask" TYPE="logical, pointer" DIM="2">
!     Land mask, true where there is land
!   </DATA>
   real, _ALLOCATABLE          :: blon(:,:) _NULL  ! longitude corners of our domain
   real, _ALLOCATABLE          :: blat(:,:) _NULL  ! latitude corners of our domain 
   logical, _ALLOCATABLE       :: mask(:,:) _NULL  ! land mask (true where ther is some land)

end type land_type
! </TYPE>


! ==== some names, for information only ======================================
logical                       :: module_is_initialized = .FALSE.
character(len=*),   parameter :: module_name = 'land_mod'
character(len=128), parameter :: version     = '$Id: land_model.F90,v 19.0 2012/01/06 20:38:57 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'


! ==== local module variables ================================================
integer                 :: n_tiles = 1  ! number of tiles
type(land_type),save    :: theLand  ! internal state of the model
integer :: landClock
integer :: isphum ! number of specific humidity in the tracer array, or NO_TRACER

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! <SUBROUTINE NAME="init_land_with_xgrid" INTERFACE="land_model_init">
!   <OVERVIEW>
!     Initializes the land model
!   </OVERVIEW>
!   <DESCRIPTION>
!     The state of the land model is initialized using the grid description
!     file as an input. This routine reads land grid corners and area of
!     land from a grid description file.
!   </DESCRIPTION>
!   <PUBLICROUTINE INTERFACE="land_model_init">
subroutine init_land_with_xgrid &
     (atmos2land, bnd, time_init, time, dt_fast, dt_slow, &
#ifdef LAND_GRID_FROM_ATMOS
      glonb, glatb, &
#endif
      atmos_domain)

  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: bnd ! land boundary data
  type(time_type), intent(in)    :: time_init ! initial time of simulation (?)
  type(time_type), intent(in)    :: time      ! current time
  type(time_type), intent(in)    :: dt_fast   ! fast time step
  type(time_type), intent(in)    :: dt_slow   ! slow time step
#ifdef LAND_GRID_FROM_ATMOS
  real,            intent(in)    :: glonb(:,:), glatb(:,:)
  type(domain2d),  intent(in), target :: atmos_domain ! domain of computations
#else
  type(domain2d),  intent(in), optional :: atmos_domain ! domain of computations
#endif
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
#ifndef LAND_GRID_FROM_ATMOS
  real, allocatable :: glonb(:,:)   ! lon corners of global grid
  real, allocatable :: glatb(:,:)   ! lat corners of global grid
#endif
  real, allocatable :: glon(:,:)    ! lon centers of global grid
  real, allocatable :: glat(:,:)    ! lat centers of global grid
  real, allocatable :: gfrac    (:,:) ! fraction of land in cells
  real, allocatable :: gcellarea(:,:) ! area of land cells
  real, allocatable :: tmpx(:,:), tmpy(:,:)
  real, allocatable :: tmp_latb(:), tmp_lat(:)
  real, allocatable :: geo_lonv(:,:)
  real, allocatable :: geo_latv(:,:)
  real, allocatable :: xgrid_area(:)
  integer, allocatable :: i1(:), j1(:)
  integer, allocatable :: i2(:), j2(:)

  integer :: nlonb, nlatb  ! number of the cell bounds in lon and lat directions
  integer :: nlon,  nlat   ! number of cells in lon and lat directions
  integer :: siz(4)        ! used to store field size
  integer,save :: ntiles=1 ! number of tiles in land mosaic, should be 1 or 6.
  integer :: nfile_axl     ! number of atmosXland exchange grid file. Can be more than one.
  integer :: nxgrid        ! number of exchange grid cell.
  integer :: i, j, n, m, tile, ind, digit

  character(len=256) :: grid_file   = "INPUT/grid_spec.nc"
  character(len=256) :: tile_file
  character(len=256) :: land_mosaic   ! land mosaic file
  character(len=256) :: axl_file      ! atmosXland exchange grid file

  ! [1] get land grid corners and area of land from the grid description file
  ! open grid description file
  ! for mosaic grid, grid area will be calculated based on exchange grid and grid cell.
  
#ifndef LAND_GRID_FROM_ATMOS
  if( field_exist(grid_file, 'AREA_LND') ) then
     call field_size( grid_file, 'AREA_LND', siz)
     nlon = siz(1)
     nlat = siz(2)
     nlonb = nlon + 1
     nlatb = nlat + 1
     ! allocate data for longitude and latitude bounds
     allocate( glonb(nlonb,nlatb), glatb(nlonb,nlatb), glon(nlon,nlat), glat(nlon,nlat),&
               gcellarea(nlon,nlat), gfrac(nlon,nlat), tmp_lat(nlat), tmp_latb(nlatb) )
     ! read coordinates of grid cell vertices
     call read_data(grid_file, "xbl", glonb(:,1), no_domain=.true.)
     call read_data(grid_file, "ybl", tmp_latb, no_domain=.true.)
     do j = 2, nlatb
        glonb(:,j) = glonb(:,1)
     enddo
     do i = 1, nlonb
        glatb(i,:) = tmp_latb(:)
     enddo

     ! read coordinates of grid cell centers
     call read_data(grid_file, "xtl", glon(:,1), no_domain=.true.)
     call read_data(grid_file, "ytl", tmp_lat(:), no_domain=.true.)
     do j = 2, nlat
        glon(:,j) = glon(:,1)
     enddo
     do i = 1, nlon
        glat(i,:) = tmp_lat(:)
     enddo

     ! read land area, land cell area and calculate the fractional land coverage 
     call read_data(grid_file, "AREA_LND_CELL", gcellarea, no_domain=.true.)
     call read_data(grid_file, "AREA_LND",      gfrac,     no_domain=.true.)
     ! calculate land fraction
     gfrac     = gfrac/gcellarea
     ! convert relative area to absolute value, m2
     gcellarea = gcellarea*4*pi*radius**2
  else if( field_exist(grid_file, 'lnd_mosaic_file') ) then ! read from mosaic file
     call read_data(grid_file, "lnd_mosaic_file", land_mosaic)     
     land_mosaic = "INPUT/"//trim(land_mosaic)
     ntiles = get_mosaic_ntiles(land_mosaic)
     tile = 1
     if (ntiles .EQ. 1) then ! lat-lon case
        call read_data (land_mosaic, "gridfiles", tile_file)
        tile_file = 'INPUT/'//trim(tile_file)
     else if (ntiles .EQ. 6) then ! cubed-sphere case
        if (present(atmos_domain)) then
           ! assumes the atmos domain has 6 tiles (probably need to check)
           call get_mosaic_tile_grid (tile_file, land_mosaic, atmos_domain)
        else
           if (mod(mpp_npes(),ntiles) .NE. 0) call error_mesg('land_model_init', &
                       'Number of processors must be a multiple of 6', FATAL)
           tile = ntiles*mpp_pe()/mpp_npes() + 1  ! assumption
           call read_data (land_mosaic, "gridfiles", tile_file, level=tile)
           tile_file = 'INPUT/'//trim(tile_file)
        endif
     else
        call error_mesg('land_model_init',  &
           ' ntiles should be 1 or 6 for land mosaic, contact developer', FATAL)
     endif
     call field_size(tile_file, "x", siz)
     if( mod(siz(1)-1,2) .NE. 0) call error_mesg("land_model_init", "size(x,1) - 1 should be divided by 2", FATAL);
     if( mod(siz(2)-1,2) .NE. 0) call error_mesg("land_model_init", "size(x,2) - 1 should be divided by 2", FATAL);
     nlon = (siz(1)-1)/2
     nlat = (siz(2)-1)/2
     nlonb = nlon + 1
     nlatb = nlat + 1
     allocate( glonb(nlonb,nlatb), glatb(nlonb,nlatb), glon(nlon,nlat), glat(nlon,nlat),&
               gcellarea(nlon,nlat), gfrac(nlon,nlat) )
     !--- read the grid information on supergrid.
     allocate( tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
     call read_data(tile_file, "x", tmpx, no_domain=.TRUE.)
     call read_data(tile_file, "y", tmpy, no_domain=.TRUE.)
     !--- make sure the grid is regular lat-lon grid.

     do j = 1, nlatb
     do i = 1, nlonb
        glonb(i,j) = tmpx(2*i-1,2*j-1)
        glatb(i,j) = tmpy(2*i-1,2*j-1)
     end do     
     end do     
     do j = 1, nlat
     do i = 1, nlon
        glon(i,j) = tmpx(2*i,2*j)
        glat(i,j) = tmpy(2*i,2*j)
     end do     
     end do     
     !--- land_cell_area will be calculated using the same way to calculate the area of xgrid.
     allocate(geo_lonv(nlonb,nlatb), geo_latv(nlonb,nlatb))
     do j = 1, nlatb
        do i = 1, nlonb
           geo_lonv(i,j) = glonb(i,j)*pi/180.0
           geo_latv(i,j) = glatb(i,j)*pi/180.0
        end do
     end do

     call calc_mosaic_grid_area(geo_lonv, geo_latv, gcellarea)
     deallocate(tmpx, tmpy, geo_lonv, geo_latv)

     !--- land area will be calculated based on exchange grid area.
     call field_size(grid_file, "aXl_file", siz)
     nfile_axl = siz(2)
     gfrac = 0
     do n = 1, nfile_axl
        call read_data(grid_file, "aXl_file", axl_file, level=n)
        if(ntiles .EQ. 6) then
          ind = index(axl_file, "land_mosaic_tile")
          if(ind ==0) call error_mesg('land_model_init', &
              'string land_mosaic_tile should contains in the field aXl_file in file '//trim(grid_file), FATAL )
          ind = ind + 16
          digit = ichar(axl_file(ind:ind)) - ichar('0')
          if( digit > 9 .OR. digit < 1) call error_mesg('land_model_init', &
              'The character following land_mosaic_tile in field aXl_file of file ' &
              //trim(grid_file)//' should be a positive integer', FATAL)
          if( tile .NE. digit ) cycle
        endif
        axl_file = 'INPUT/'//trim(axl_file)
        nxgrid = get_mosaic_xgrid_size(axl_file)
        allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
        call get_mosaic_xgrid(aXl_file, i1, j1, i2, j2, xgrid_area)
        do m = 1, nxgrid
           i = i2(m); j = j2(m)
           gfrac(i,j) = gfrac(i,j) + xgrid_area(m)
        end do
        deallocate(i1, j1, i2, j2, xgrid_area)
     end do     
     gfrac = gfrac*4*PI*RADIUS*RADIUS/gcellarea
  else
     call error_mesg('land_model_init','both AREA_LND and lnd_mosaic_file do not exist in file '//trim(grid_file),FATAL)
  end if

  ! [2] convert lon and lat coordinates from degrees to radian
  glonb = glonb*pi/180.0
  glatb = glatb*pi/180.0
  glon  = glon*pi/180.0
  glat  = glat*pi/180.0

#else

! initialize grid info from input arguments
  nlonb = size(glonb,1)
  nlatb = size(glonb,2)
  nlon = nlonb-1
  nlat = nlatb-1
  allocate( glon(nlon,nlat), glat(nlon,nlat),&
            gcellarea(nlon,nlat), gfrac(nlon,nlat) )

! approximate coordinates of grid cell centers
  glon = (glonb(1:nlon,1:nlat)+glonb(2:nlon+1,1:nlat)+glonb(1:nlon,2:nlat+1)+glonb(2:nlon+1,2:nlat+1))/4.
  glat = (glatb(1:nlon,1:nlat)+glatb(2:nlon+1,1:nlat)+glatb(1:nlon,2:nlat+1)+glatb(2:nlon+1,2:nlat+1))/4.

! read land fraction from file (land=1)
! otherwise initialize with no land
  if (file_exist('INPUT/land_mask.nc')) then
     call read_data ('INPUT/land_mask.nc', 'land_mask', gfrac, no_domain=.true.) !, atmos_domain)
     where (gfrac > 0.5)
        gfrac = 1.0
     elsewhere
        gfrac = 0.0
     endwhere
  else
     gfrac = 0.0 ! aqua-planet
  endif

! uniform area
  gcellarea = 1./real(6*(nlonb-1)*(nlatb-1))
  ! convert relative area to absolute value, m2
  gcellarea = gcellarea*4*pi*radius**2

#endif

  ! [3] initialize the land model using the global grid we just obtained
  call init_land_with_area &
       ( atmos2land, bnd, glonb, glatb, gcellarea, gfrac, time, dt_fast, dt_slow, &
       atmos_domain, glon=glon, glat=glat, numtiles=ntiles )

  ! [4] deallocate memory that we are not going to use any more
#ifndef LAND_GRID_FROM_ATMOS
  deallocate ( glonb, glatb, glon, glat, gfrac, gcellarea )
#else
  deallocate ( glon, glat, gfrac, gcellarea )
#endif

!   <NOTE>
!     Theoretically, the grid description file can specify any regular
!     rectangular grid for land, not just lon/lat grid. Therefore the variables
!     "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
!     boundaries of the grid.
!     However, at this time the module land_properties assumes that grid _is_
!     lon/lat and therefore the entire module also have to assume that the
!     land grid is lon/lat. lon/lat grid is also assumed for the diagnostics, 
!     but this is probably not so critical. 
!   </NOTE>
!Balaji: looks like this is the only public procedure under land_model_init
!        called from coupler_init
  landClock = mpp_clock_id( 'Land', flags=clock_flag_default, grain=CLOCK_COMPONENT )
!should also initialize subcomponent clocks for veg, hydro, etc. here
end subroutine init_land_with_xgrid
! </SUBROUTINE>


! <SUBROUTINE NAME="init_land_with_area" INTERFACE="land_model_init">
!   <OVERVIEW>
!     Initializes the land model
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initializes the land state using land area for each of the grid points.
!   </DESCRIPTION>
!   <PUBLICROUTINE INTERFACE="land_model_init">
subroutine init_land_with_area &
     ( atmos2land, bnd, gblon, gblat, gcellarea, gfrac, time, dt_fast, dt_slow, domain, &
     glon, glat, numtiles )

  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: bnd        ! state to update
  real,              intent(in) :: gblon(:,:)! lon corners of the grid cells
  real,              intent(in) :: gblat(:,:)! lat corners of the grid cells
  real,              intent(in) :: gcellarea(:,:) ! full area of the grid cells
  real,              intent(in) :: gfrac(:,:)     ! fraction of land in the grid cell
  type(time_type),   intent(in) :: time    ! current time
  type(time_type),   intent(in) :: dt_fast ! fast time step
  type(time_type),   intent(in) :: dt_slow ! slow time step
  type(domain2d),    intent(in), optional :: domain  ! domain of computations
  real,              intent(in), optional :: glon(:,:), glat(:,:) ! centers
  integer,           intent(in), optional :: numtiles
                          ! of the grid cells
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: layout(2) ! layout of our domains
  integer :: nlon, nlat
  integer :: is,ie,js,je,k
  integer :: id_lon, id_lat  ! IDs of land diagnostic axes
  integer :: ntrace, ntprog, ntdiag, ntfamily ! numbers of tracers
  integer :: ntiles ! number of tiles for land mosaic

  module_is_initialized = .TRUE.

  ! write the version and tagname to the logfile
  call write_version_number(version, tagname)

! initialize tracers
#ifdef LAND_BND_TRACERS
  ! register land model tracers and find specific humidity
  call register_tracers ( MODEL_LAND, ntrace, ntprog, ntdiag, ntfamily )
  isphum = get_tracer_index ( MODEL_LAND, 'sphum' )
  if (isphum==NO_TRACER) then
     call error_mesg('land_model_init','no required "sphum" tracer',FATAL)
  endif
#else
  isphum = NO_TRACER
#endif

  ! get the size of the global grid
  nlon = size(gfrac, 1)
  nlat = size(gfrac, 2)
  ! number of land mosaic tiles
  ntiles = 1
  if (present(numtiles)) ntiles = numtiles
  if (mod(mpp_npes(),ntiles) .NE. 0) call error_mesg('land_model_init', &
      'Number of processors must be a multiple of number of mosaic tiles', FATAL)
  ! get the processor layout information, either from upper-layer module
  ! decomposition, or define according to the grid size
  if(present(domain)) then
     call mpp_get_layout (domain, layout)
  else
     call mpp_define_layout &
          ((/1,nlon, 1,nlat/),mpp_npes()/ntiles,layout)
  endif

  ! create our computation domain according to obtained processor layout 
  if (ntiles .EQ. 1) then
     call mpp_define_domains ( &
          (/1,nlon, 1, nlat/),           &  ! global grid size
          layout,                        &  ! layout for our domains
          theLand%domain,                &  ! domain to define
          xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL' )
  else if (ntiles .EQ. 6) then
     call define_cube_mosaic ( theLand%domain, nlon, nlat, layout )
  else
     call error_mesg('land_model_init',  &
        'number of mosaic tiles should be 1 or 6', FATAL)
  endif

  ! get the size of our computation domain
  call mpp_get_compute_domain ( theLand%domain, is,ie,js,je )

  theLand%is=is; theLand%ie=ie; theLand%js=js; theLand%je=je
  theLand%n_tiles=n_tiles

  ! NOTE: the piece of code from here to the end of the subroutine should
  ! probably be revised to accommodate tiling structure. The problem is that
  ! it is theoretically possible that the soil and vegetation have different
  ! tiling, or there may be an additional sub-model with yet different tiling. 
  ! It is case it is not quite clear what does the "land tile" mean -- it may 
  ! be the tile of the uppermost land sub-model (since this is what atmosphere 
  ! sees) or it may be decart product of all the sub-model tiling structures, 
  ! just to name two possibilities.

  ! allocate fractional area for soil

  ! allocate storage for the boundary data 
  call allocate_boundary_data ( atmos2land, bnd, theLand%domain, n_tiles, ntprog )

  ! set up boundary values that we know already
  bnd % domain = theLand % domain
  do k = 1, size(bnd%tile_size,3)
     where (gfrac(is:ie,js:je)>0)
        bnd % tile_size (:,:,k) = 1.0/n_tiles
        bnd % mask      (:,:,k) = .true.
     elsewhere
        bnd % tile_size (:,:,k) = 0.0
        bnd % mask      (:,:,k) = .false.
     endwhere
  enddo

  ! init grid variables: this is a quick fix since it is not yet
  ! clear how to initialize land properties with arbitrary grid
  allocate(theLand%blon(is:ie+1,js:je+1), theLand%blat(is:ie+1,js:je+1), theLand%mask(is:ie,js:je))
  theLand%blon(:,:)   = gblon(is:ie+1,js:je+1)
  theLand%blat(:,:)   = gblat(is:ie+1,js:je+1)
  theLand%mask(:,:)   = any( bnd % mask, DIM=3 )

  ! initialize land diagnostics -- create diagnostic axes to be used by all
  ! diagnostic routine in submodules
  call init_land_diag ( gblon, gblat, theLand%domain, id_lon, id_lat, &
       glon=glon, glat=glat )
  bnd%axes = (/id_lon, id_lat/)

  ! initialize all submodules
  call soil_init ( theLand%soil, gblon, gblat, gcellarea, gfrac, &
       time, dt_fast, dt_slow, theLand%domain, bnd % tile_size, bnd % mask, &
       id_lon, id_lat )

  call vegetation_init ( &
       theLand%vegetation, gblon, gblat, gcellarea, gfrac, &
       time, dt_fast, dt_slow, theLand%domain, bnd % tile_size, bnd % mask, &
       id_lon, id_lat )

  ! update boundary conditions 
  call update_land_bnd_slow ( bnd )
  call update_land_bnd_fast ( bnd )

!   <NOTE>
!     NOTE: we assume every sub-model has the same domain layout and is on
!     the same grid, so we do not need to do general-case flux exchange
!     between them.
!   </NOTE>

end subroutine init_land_with_area
! </SUBROUTINE>



! <SUBROUTINE NAME="land_model_end">

!   <OVERVIEW>
!     Destructs the land model.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Destructs the land model and its sub-models. Deallocates the arrays
!     and the boundary exchange data 
!   </DESCRIPTION>

!   <TEMPLATE>
!     call land_model_end ( bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine land_model_end ( atmos2land, bnd )
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd        ! state to update
!   </PUBLICROUTINE>

  module_is_initialized = .FALSE.
  call soil_end ( theLand%soil )
  call vegetation_end ( theLand%vegetation )

  ! deallocate variables
  deallocate ( theLand%blon, theLand%blat, theLand%mask )

  ! deallocate boundary exchange data
  call deallocate_boundary_data ( atmos2land, bnd )
  
end subroutine land_model_end
! </SUBROUTINE>

subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name

    call error_mesg ('land_model_restart in land_model_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)

end subroutine land_model_restart

! <SUBROUTINE NAME="update_land_model_fast">

!   <OVERVIEW>
!     Updates the state of the land model on the fast time scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates the state of the land model on the fast time scale. To implement
!     implicit calculation scheme, first the land sub-models are called in
!     top-bottom general order to evaluate the derivatives; then they are
!     called in bottom-top order to finish implicit calculations. On the
!     upward pass they have a chance to update boundary values returned to
!     the atmosphere. After calculations, update the boundary conditions
!     visible to the atmosphere.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_land_model_fast( atmos2land, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="">
subroutine update_land_model_fast ( atmos2land, bnd )

  type(atmos_land_boundary_type), intent(inout)    :: atmos2land ! quantities exchanged between
                                                              ! the atmosphere and the land
  type(land_data_type),  intent(inout) :: bnd ! state to update
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  real, dimension(theLand%is:theLand%ie, theLand%js:theLand%je, theLand%n_tiles) :: &
       q_flux, dedt     ! slm Mar 19 2002 

  call mpp_clock_begin(landClock)
  call set_domain (theLand%domain)

  call update_vegetation_fast_down (      &
       theLand%vegetation, theLand%soil,  &
#ifdef LAND_BND_TRACERS
       atmos2land%tr_flux(:,:,:,isphum),  &
       atmos2land%dfdtr(:,:,:,isphum),    &
#else
       atmos2land%q_flux,                 &
       atmos2land%dedq,                   &
#endif
       atmos2land%drag_q,                 &
       atmos2land%p_surf,                 &
       q_flux,  &
       dedt     )
  
  call update_soil_fast ( theLand%soil, &
       atmos2land%sw_flux, &
       atmos2land%lw_flux, &
       atmos2land%t_flux,  &
       q_flux,             &
       atmos2land%dhdt,    &
       dedt,               &
       atmos2land%drdt,    &
       atmos2land%lprec,   &
       atmos2land%fprec    )
  
  call update_vegetation_fast_up (        &
       theLand%vegetation, theLand%soil,  &
       atmos2land%drag_q,                 &
#ifdef LAND_BND_TRACERS
       atmos2land%tr_flux(:,:,:,isphum),  &
       atmos2land%dfdtr(:,:,:,isphum)     )
#else
       atmos2land%q_flux,                 &
       atmos2land%dedq                    )
#endif
  ! after calculations, update boundary conditions visible to the atmosphere
  call update_land_bnd_fast ( bnd )

  call nullify_domain()
  call mpp_clock_end(landClock)
end subroutine update_land_model_fast
! </SUBROUTINE>



! <SUBROUTINE NAME="update_land_model_slow">
!   <OVERVIEW>
!     Updates the state of the land model on the slow time scale.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Updates land, and boundary conditions if
!     necessary, on the slow time scale. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call update_land_model_slow( bnd )
!   </TEMPLATE>
!   <PUBLICROUTINE>
subroutine update_land_model_slow ( atmos2land, bnd )
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd        ! state to update
!   </PUBLICROUTINE>

  call mpp_clock_begin(landClock)
  call set_domain (theLand%domain)
  call update_vegetation_slow(theLand%vegetation)
  call update_soil_slow(theLand%soil, theLand%blon, theLand%blat)

  ! update boundary conditions on slow time-scale, if necessary
  call update_land_bnd_slow ( bnd )
  call nullify_domain()
  call mpp_clock_end(landClock)

end subroutine update_land_model_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="update_land_bnd_fast">
!   <OVERVIEW>
!     Updates land boundary data.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Updates land boundary data (the ones that atmosphere sees) on the
!     fast time scale. This routine does not update the tiling structure
!     because it is assumed that the tiling does not change on fast time
!     scale.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call update_land_bnd_fast ( bnd )
!   </TEMPLATE>
!   <PUBLICROUTINE>
subroutine update_land_bnd_fast ( bnd )

  type(land_data_type), intent(inout) :: bnd        ! state to update
!   </PUBLICROUTINE>

  call update_soil_bnd_fast ( theLand%soil, bnd )
  call update_vegetation_bnd_fast ( theLand%vegetation, bnd )

end subroutine update_land_bnd_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="update_land_bnd_slow">
!   <OVERVIEW>
!     Updates the land boundary data for the atmosphere on the slow time scale.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Updates the land boundary data for the atmosphere on the slow time scale.
!     This subroutine is responsible for changing of the tiling structure, if
!     necessary, as well as for changing albedo, drag coefficients and such.     
!   </DESCRIPTION>
!   <TEMPLATE>
!     call update_land_bnd_slow ( bnd )
!   </TEMPLATE>
!   <PUBLICROUTINE>
subroutine update_land_bnd_slow ( bnd )

  type(land_data_type), intent(inout) :: bnd        ! state to update
!   </PUBLICROUTINE>

  call update_soil_bnd_slow ( theLand%soil, bnd )
  call update_vegetation_bnd_slow ( theLand%vegetation, bnd )

end subroutine update_land_bnd_slow
!   <NOTE>
!     If the tiling structure has been modified, then probably the
!     distribution of other boundary values, such as temperature and
!     surface humidity, should be modified too, not yet clear how.
!   </NOTE>
! </SUBROUTINE>

! ============================================================================
! calculate and return total amount of requested quantitiy per PE
subroutine Lnd_stock_pe(bnd,index,value)
  type(land_data_type), intent(in)  :: bnd
  integer             , intent(in)  :: index ! ID of the stock to calculate
  real                , intent(out) :: value ! calculated value of the stock

  real :: soil_value, vegn_value

  if(bnd%pe) then
    call       soil_stock_pe(theLand%soil,       index, soil_value)
    call vegetation_stock_pe(theLand%vegetation, index, vegn_value)
    value = soil_value+vegn_value
  else
    value = 0.0
  endif

end subroutine Lnd_stock_pe


! <SUBROUTINE NAME="init_land_diag">
!   <OVERVIEW>
!     Initializes the horizontal axes for the land grid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialize the horizontal axes for the land grid so that all
!     submodules use them, instead of creating their own.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_land_diag(glonb, glatb, domain, id_lon, id_lat)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine init_land_diag(glonb, glatb, domain, id_lon, id_lat, glon, glat)

  real, intent(in)           :: glonb(:,:), glatb(:,:) ! longitude and latitude corners of grid
                                                       ! cells, specified for the global grid
                                                       ! (not only domain)
  type(domain2d), intent(in) :: domain ! the domain of operations
  integer, intent(out)       :: id_lon, id_lat         ! IDs of respective diag. manager axes
  real, intent(in),optional  :: glon(:,:), glat(:,:)   ! coordinates of grid cell centers
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: id_lonb, id_latb ! IDs for cell boundaries
  integer :: nlon, nlat       ! sizes of respective axes, just for convenience
  real    :: rad2deg          ! conversion factor radian -> degrees
  real    :: lon(size(glonb,1)-1), lat(size(glatb,2)-1)
  integer :: i, j

  rad2deg = 180./pi
  nlon = size(glonb,1)-1
  nlat = size(glatb,2)-1

  if (is_latlon(glonb,glatb)) then

 !----- lat/lon grid -----

     if(present(glon)) then
        lon = glon(:,1)
     else
        lon = (glonb(1:nlon,1)+glonb(2:nlon+1,1))/2
     endif
     if(present(glat)) then
        lat = glat(1,:)
     else
        lat = (glatb(1,1:nlat)+glatb(1,2:nlat+1))/2
     endif

     ! define longitude axes and its edges
     id_lonb = diag_axis_init ( &
          'lonb', glonb(:,1)*rad2deg, 'degrees_E', 'X', 'longitude edges', &
          set_name='land', domain2=domain )
     id_lon  = diag_axis_init (                                                &
          'lon', lon(:)*rad2deg, 'degrees_E', 'X',  &
          'longitude', set_name='land',  edges=id_lonb, domain2=domain )

     ! define latitude axes and its edges
     id_latb = diag_axis_init ( &
          'latb', glatb(1,:)*rad2deg, 'degrees_N', 'Y', 'latitude edges',  &
          set_name='land',  domain2=domain   )
     id_lat = diag_axis_init (                                                &
          'lat', lat(:)*rad2deg, 'degrees_N', 'Y', &
          'latitude', set_name='land', edges=id_latb, domain2=domain   )

!----- cubed sphere grid -----

  else
     do i = 1, nlon
       lon(i) = real(i)
     enddo
     do j = 1, nlat
       lat(j) = real(j)
     enddo

     id_lon = diag_axis_init ( 'grid_xt', lon, 'degrees_E', 'X', &
            'T-cell longitude', set_name='land',  domain2=domain )
     id_lat = diag_axis_init ( 'grid_yt', lat, 'degrees_N', 'Y', &
             'T-cell latitude', set_name='land',  domain2=domain )

  endif


end subroutine init_land_diag
! </SUBROUTINE>

! <SUBROUTINE NAME="define_cube_mosaic">
subroutine define_cube_mosaic ( Domain, nx, ny, layout )
type(domain2d), intent(inout) :: Domain
integer, intent(in) :: nx, ny, layout(2)
integer, parameter :: ntiles = 6, num_contact = 12, ng = 0
integer :: global_indices(4,ntiles), layout_2d(2,ntiles)
integer :: pe_start(ntiles), pe_end(ntiles)
integer :: east(4), west(4), north(4), south(4)
integer :: tile1(num_contact), tile2(num_contact),  &
           index1(num_contact,4), index2(num_contact,4)
integer :: n

    do n = 1, ntiles
       global_indices(:,n) = (/ 1, nx, 1, ny /)
       layout_2d     (:,n) = layout
       pe_start        (n) = (n-1)*layout(1)*layout(2)
       pe_end          (n) =     n*layout(1)*layout(2) - 1
    enddo

  ! istart, iend, jstart, jend for each face
    west =  (/  1,  1,  1, ny /)
    east =  (/ nx, nx,  1, ny /)
    south = (/  1, nx,  1,  1 /)
    north = (/  1, nx, ny, ny /)

  ! define the 12 contact lines bewteen the faces
    tile1( 1) = 1; index1( 1,:) = east;   tile2( 1) = 2; index2( 1,:) = west
    tile1( 2) = 1; index1( 2,:) = north;  tile2( 2) = 3; index2( 2,:) = west
    tile1( 3) = 1; index1( 3,:) = west;   tile2( 3) = 5; index2( 3,:) = north
    tile1( 4) = 1; index1( 4,:) = south;  tile2( 4) = 6; index2( 4,:) = north
    tile1( 5) = 2; index1( 5,:) = north;  tile2( 5) = 3; index2( 5,:) = south
    tile1( 6) = 2; index1( 6,:) = east;   tile2( 6) = 4; index2( 6,:) = south
    tile1( 7) = 2; index1( 7,:) = south;  tile2( 7) = 6; index2( 7,:) = east
    tile1( 8) = 3; index1( 8,:) = east;   tile2( 8) = 4; index2( 8,:) = west
    tile1( 9) = 3; index1( 9,:) = north;  tile2( 9) = 5; index2( 9,:) = west
    tile1(10) = 4; index1(10,:) = north;  tile2(10) = 5; index2(10,:) = south
    tile1(11) = 4; index1(11,:) = east;   tile2(11) = 6; index2(11,:) = south
    tile1(12) = 5; index1(12,:) = east;   tile2(12) = 6; index2(12,:) = west

  ! create the domain2d variable
    call mpp_define_mosaic ( global_indices, layout_2d, Domain,                  &
                             ntiles, num_contact, tile1, tile2,                  &
                             index1(:,1), index1(:,2), index1(:,3), index1(:,4), &
                             index2(:,1), index2(:,2), index2(:,3), index2(:,4), &
                             pe_start=pe_start, pe_end=pe_end, symmetry=.true.,  &
                             shalo = ng, nhalo = ng, whalo = ng, ehalo = ng,     &
                             name = "Land Model Cubic-Grid" )

end subroutine define_cube_mosaic
! </SUBROUTINE>

! <DIAGFIELDS>
!   <NETCDF NAME="lonb" UNITS="degrees_E">
!     Longitude cell boundaries
!   </NETCDF>
!   <NETCDF NAME="lon" UNITS="degrees_E">
!     Diagnostic manager longitude axis
!   </NETCDF>
!   <NETCDF NAME="latb" UNITS="degrees_N">
!     Latitude cell boundaries
!   </NETCDF>
!   <NETCDF NAME="lat" UNITS="degrees_N">
!     Diagnostic manager latitude axis
!   </NETCDF>
! </DIAGFIELDS>

end module land_model_mod
   

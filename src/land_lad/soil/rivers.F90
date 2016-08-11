#include <fms_platform.h>

! ============================================================================
! River-related processes
! ============================================================================

! ---- useful macro definitions ----------------------------------------------
#define __ERROR__(message) \
call print_error_mesg(mod_name,message,FATAL,file_name,__LINE__)
#define __NOTE__(message) \
call print_error_mesg(mod_name,message,NOTE,file_name,__LINE__)
#define __ASSERT__(x, message) \
if(.NOT.(x))call print_error_mesg(mod_name,message,FATAL,file_name,__LINE__)

module rivers_mod

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
!   Module containing processes relating to the rivers.
! </OVERVIEW>

! <DESCRIPTION>
!   Allocates river data and sets river destination indices. Includes tests
!   for checking whether points are coastal. Updates the state of the rivers
!   and runoff on the fast and slow time-scales. River boundary data is
!   updated on the slow time-scale. Given the local runoff field, calculates
!   the local portion of discharge.
!
!   Initializes diagnostics and sends the diagnostic field output for the
!   slow time-scale river fields and static fields to the diagnostic manager.
! </DESCRIPTION>

  use time_manager_mod,    only: time_type, get_time, increment_time

  use diag_manager_mod,    only: diag_axis_init, register_diag_field,      &
                                 register_static_field, send_data

  use mpp_mod,             only: mpp_sum, mpp_send, mpp_recv, mpp_npes,    &
                                 mpp_sync_self

  use mpp_domains_mod,     only: domain2d, mpp_get_compute_domain,         &
                                 mpp_get_global_domain, mpp_global_field
  use mpp_io_mod,          only: mpp_open, MPP_RDONLY

  use fms_mod,             only: open_namelist_file, stdlog, error_mesg, FATAL, NOTE, &
                                 file_exist, close_file, check_nml_error,mpp_pe,      &
                                 mpp_root_pe, write_version_number, stderr

  use constants_mod,       only: pi

  use land_types_mod,      only: land_data_type
 
  use land_properties_mod, only: regrid_discrete_field

  use numerics_mod,        only: bisect, is_latlon


implicit none
private


! ==== public interfaces =====================================================
public :: rivers_type        

public :: rivers_init        ! initialize rivers module
public :: rivers_end         ! finishing actions 

public :: update_rivers_fast ! fast time-scale update of state
public :: update_rivers_slow ! slow time-scale update of state

public :: update_rivers_bnd_slow
! ==== end of public interfaces ==============================================

! <TYPE NAME="rivers_type">

!   <DESCRIPTION>
!     Describes the domain and physical properties of the rivers.
!   </DESCRIPTION>

type rivers_type
   private

!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     Computational domain
!   </DATA>

!   <DATA NAME="is" TYPE="integer">
!     Computational domain bounds
!   </DATA>

!   <DATA NAME="ie" TYPE="integer">
!     Computational domain bounds
!   </DATA>

!   <DATA NAME="js" TYPE="integer">
!     Computational domain bounds
!   </DATA>

!   <DATA NAME="je" TYPE="integer">
!     Computational domain bounds
!   </DATA>

!   <DATA NAME="gnlon" TYPE="integer">
!     Size of global grid
!   </DATA>

!   <DATA NAME="gnlat" TYPE="integer">
!     Size of global grid
!   </DATA>

   type(domain2d)   :: domain      ! our domain
   integer :: is,ie,js,je          ! computational domain bounds
   integer :: gnlon, gnlat         ! size of global grid

!   <DATA NAME="i_dest" TYPE="integer, pointer" DIM="2">
!     Longitude destination indices
!   </DATA>

!   <DATA NAME="j_dest" TYPE="integer, pointer" DIM="2">
!     Latitude destination indices
!   </DATA>

!   <DATA NAME="warea" UNITS="m2" TYPE="real, pointer" DIM="2">
!     Area covered by water per grid cell
!   </DATA>

!   <DATA NAME="larea" UNITS="m2" TYPE="real, pointer" DIM="2">
!     Area of land per grid cell
!   </DATA>

   integer, _ALLOCATABLE :: i_dest(:,:) _NULL ! lon dst indices
   integer, _ALLOCATABLE :: j_dest(:,:) _NULL  ! lat dst indices
   real,    _ALLOCATABLE :: warea(:,:) _NULL  ! area covered by water per grid cell
   real,    _ALLOCATABLE :: larea(:,:) _NULL  ! area of land per grid cell

!   <DATA NAME="discharge" UNITS="kg/m2/s" TYPE="real, pointer" DIM="2">
!     Water discharge
!   </DATA>

!   <DATA NAME="discharge_snow" UNITS="kg/m2/s" TYPE="real, pointer" DIM="2">
!     Snow discharge
!   </DATA>

   real,    _ALLOCATABLE :: discharge     (:,:) _NULL ! water discharge
   real,    _ALLOCATABLE :: discharge_snow(:,:) _NULL ! snow discharge
   
!   <DATA NAME="time" UNITS="s" TYPE="time_type">
!     Current time
!   </DATA>

!   <DATA NAME="dt" UNITS="s" TYPE="real">
!     Fast time step
!   </DATA>

!   <DATA NAME="dt_slow" UNITS="s" TYPE="real">
!     Slow time step
!   </DATA>

   type(time_type)  :: time         ! current time
   real             :: dt           ! fast time step, s
   real             :: dt_slow      ! slow time step, s
end type rivers_type
! </TYPE>

! ---- module constants -----------------------------------------------------
! some names, for information only
character(len=*), parameter :: mod_name = 'rivers'
character(len=*), parameter :: file_name = __FILE__
character(len=128), parameter :: version  = '$Id: rivers.F90,v 19.0 2012/01/06 20:39:34 fms Exp $'
character(len=128), parameter :: tagname      = '$Name: tikal $'


! <NAMELIST NAME="rivers_nml">

!   <DATA NAME="use_single_basin" TYPE="logical" DEFAULT=".false.">
!     If true, all river discharge goes to a single point
!   </DATA>

!   <DATA NAME="i_dest0" TYPE="integer" DEFAULT="1">
!     Destination point longitude index
!   </DATA>

!   <DATA NAME="j_dest0" TYPE="integer" DEFAULT="1">
!     Destination point latitude index
!   </DATA>

!   <DATA NAME="min_water_frac" TYPE="real" DEFAULT="1e-4">
!     Minimum water area fraction for the point to be a
!     valid discharge destination
!   </DATA>

!   <DATA NAME="min_land_frac" TYPE="real" DEFAULT="1e-4">
!     "Land" for discharge destination calculations
!   </DATA>

logical :: use_single_basin = .false. ! if true, all river discharge goes to 
                                      ! a single point
integer :: i_dest0 = 1 ! dst point lon index
integer :: j_dest0 = 1 ! dst point lat index
real    :: min_water_frac = 1e-4 ! min water area fraction for the point to be valid 
                          ! discharge destination
real    :: min_land_frac  = 1e-4 ! min land fraction for the cell to be considered 
                          ! "land" for discharge destination calculations
! </NAMELIST>


namelist/rivers_nml/ use_single_basin, i_dest0, j_dest0, &
     min_land_frac, min_water_frac

! ---- diagnostic field ids --------------------------------------------------
integer :: id_discharge, id_discharge_snow
integer :: id_adis_water,id_adis_snow ! area-weighted discharges
integer :: id_dest,  id_basin
logical :: module_is_initialized =.FALSE.

contains  ! ******************************************************************


! <SUBROUTINE NAME="rivers_init">

!   <OVERVIEW>
!     Initializes rivers data instance.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Initializes rivers data instance.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call rivers_init ( Rivers, gblon, gblat, garea, gfrac, time, dt_fast,
!        dt_slow, domain, id_lon, id_lat )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine rivers_init &
     ( Rivers, gblon, gblat, garea, gfrac, time, dt_fast, dt_slow, domain, &
     id_lon, id_lat )

  type(rivers_type), intent(inout) :: Rivers     ! data to initialize
  real,            intent(in)      :: gblon(:,:) ! lon corners of the grid cells
  real,            intent(in)      :: gblat(:,:) ! lat corners of the grid cells
  real,            intent(in)      :: garea(:,:) ! total area of each grid cell
  real,            intent(in)      :: gfrac(:,:) ! fraction of the cell covered by land
  type(time_type), intent(in)      :: time       ! current time
  type(time_type), intent(in)      :: dt_fast    ! fast time step
  type(time_type), intent(in)      :: dt_slow    ! slow time step
  type(domain2d),  intent(in)      :: domain     ! our domain
  integer,         intent(in)      :: id_lon     ! ID of land longitude (X) diag axis
  integer,         intent(in)      :: id_lat     ! ID of land latitude (Y) diag axis
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: unit                 ! unit number for i/o
  integer :: io                   ! status code for i/o
  integer :: ierr                 ! error code
  integer :: sec, day             ! components of time
  integer :: logunit

  ! copy domain
  Rivers%domain = domain

  ! get boundaries of our domain
  call mpp_get_compute_domain & 
       ( Rivers%domain, Rivers%is,Rivers%ie,Rivers%js,Rivers%je )

  ! and the size of global domain
  Rivers%gnlon = size(gblon,1)-1
  Rivers%gnlat = size(gblat,2)-1

  ! set up time-related values
  rivers % time = time
  call get_time(dt_fast, sec, day); rivers%dt      = day*86400.0+sec
  call get_time(dt_slow, sec, day); rivers%dt_slow = day*86400.0+sec

  ! read namelist
  if (file_exist('input.nml')) then
     unit = open_namelist_file ()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=rivers_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'rivers_nml')
     enddo
10   continue
     call close_file (unit)
  endif

  ! write the namelist to the log file
  if ( mpp_pe() == mpp_root_pe() )  then
     call write_version_number(version, tagname)
     logunit = stdlog()
     write (logunit, nml=rivers_nml)
  endif

  ! allocate data
  allocate( &
       Rivers%i_dest(size(garea,1),size(garea,2)), &
       Rivers%j_dest(size(garea,1),size(garea,2)), &
       Rivers%larea (size(garea,1),size(garea,2)), &
       Rivers%warea (Rivers%is:Rivers%ie,Rivers%js:Rivers%je), &
       Rivers%discharge     (Rivers%is:Rivers%ie,Rivers%js:Rivers%je), &
       Rivers%discharge_snow(Rivers%is:Rivers%ie,Rivers%js:Rivers%je), &
       stat = ierr)
  __ASSERT__(ierr==0,"Cannot allocate river data")

  Rivers%larea = garea*gfrac
  Rivers%warea = garea(Rivers%is:Rivers%ie,Rivers%js:Rivers%je)  & 
       * (1-gfrac(Rivers%is:Rivers%ie,Rivers%js:Rivers%je))

  if(use_single_basin) then
     Rivers%i_dest = i_dest0
     Rivers%j_dest = j_dest0
  else
     call init_routing ( gblon, gblat, gfrac, &
           1,size(gfrac,1), 1,size(gfrac,2), &
           Rivers%i_dest, Rivers%j_dest )
  endif

  call init_rivers_diag ( id_lon, id_lat, Time )

  call diag_static ( Rivers )

  Rivers%discharge      = 0;
  Rivers%discharge_snow = 0;
       
  module_is_initialized =.TRUE. 

!   <ERROR MSG="Cannot allocate river data" STATUS="FATAL">
!   </ERROR>

end subroutine rivers_init
! </SUBROUTINE>


! <SUBROUTINE NAME="rivers_end">

!   <OVERVIEW>
!     Deallocates rivers data instance.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Deallocates rivers data instance.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call rivers_end(Rivers)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine rivers_end(Rivers)

  type(rivers_type), intent(inout) :: Rivers ! data to deallocate/finish
!   </PUBLICROUTINE>

  ! insert saving of restart file here, if necessary
  module_is_initialized =.FALSE.
  ! release memory
  deallocate (               &
       Rivers%i_dest,        &
       Rivers%j_dest,        &
       Rivers%warea,         &
       Rivers%larea,         &
       Rivers%discharge,     &
       rivers%discharge_snow  )

end subroutine rivers_end
! </SUBROUTINE>


! <SUBROUTINE NAME="init_routing">

!   <OVERVIEW>
!     Sets river destination indices.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Sets river destination indices.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call init_routing ( glonb, glatb, gfrac, is,ie, js,je, i_dest,j_dest ) 
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine init_routing ( glonb, glatb, gfrac, is,ie, js,je, i_dest,j_dest )

  real,    intent(in)  :: glonb(:,:)          ! lon corners of the global grid
  real,    intent(in)  :: glatb(:,:)          ! lat corners of the global grid
  real,    intent(in)  :: gfrac(:,:)          ! global array of land fractional area
  integer, intent(in)  :: is,ie,js,je         ! boundaries of our domain

  integer, intent(out) :: i_dest(is:ie,js:je) ! lon index of dest points
  integer, intent(out) :: j_dest(is:ie,js:je) ! lat index of dest points
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: unit               ! unit number for i/o
  integer :: ierr               ! error code
  character(80) :: input_format ! input format of the data
  integer :: n_lon_in, n_lat_in ! number of lat and lon in input grid
  real    :: wb_degrees         ! western boundary of input greed, degrees

  real    :: wb, sb             ! western and southern bounds of input data
  real    :: dx, dy             ! size of input grid cells

  integer :: i, j, n, ii        ! various iterators
  integer :: n_dest             ! number of dest points in the list
  logical :: match              ! true if dest point has been found in the list
  integer :: &
       i_dest_bad, j_dest_bad, &! dest point which falls on land
       i_dest_new, j_dest_new   ! redirected dest point

  integer :: next_i, next_j     ! coordinates of next redir. point to try
  integer :: gnlon, gnlat       ! size of the global grid
  real    :: lon, lat           ! work values for destination calculations

  integer, allocatable :: input_1(:,:), input_2(:,:) ! input buffers
  integer, allocatable :: route_input(:,:)
  integer, allocatable :: route_compact(:,:), route_list(:)

  integer, allocatable :: route(:,:)

  ! path to find redirection points
  integer, parameter :: n_didj = 80
  integer, parameter :: di(n_didj) = &
      (/ -4, 4,-4, 4,-4, 4,-3, 3, 4,-4, 3,-3,-2, 2,-4, 4,-4, 4,-2, 2, &
          3,-3, 3,-3,-1, 1, 4,-4,-4, 4,-1, 1, 0, 0,-4, 4,             &
          2,-3, 3,-3, 3,-3, 2,-2,                                     &
          1,-1, 3,-3,-3, 3, 1,-1, 3,-3, 0, 0,                         &
         -2, 2,-2, 2,-1, 1, 2,-2,-2, 2,-1, 1, 0, 0,-2, 2,             &
          1,-1, 1,-1, 1,-1, 0, 0/)

  integer, parameter :: dj(n_didj) = &
      (/  4,-4,-4, 4, 3,-3,-4, 4, 3,-3,-4, 4, 4,-4,-2, 2,-2,-2,-4, 4, &
          3,-3,-3, 3,-4, 4,-1, 1,-1, 1, 4,-4,-4, 4, 0, 0,             &
          3,-2,-2, 2, 2,-2,-3, 3,                                     &
         -3, 3, 1,-1, 1,-1, 3,-3, 0, 0, 3,-3,                         &
          2,-2,-2, 2,-2, 2,-1, 1,-1, 1, 2,-2,-2, 2, 0, 0,             &
          1,-1,-1, 1, 0, 0, 1,-1/)

  ! diagnostics-only variables
  integer, dimension(:,:), allocatable :: moved_to ! moved river destination points
  integer, dimension(:,:), allocatable :: bad_dest ! array of bad river destination points -- the ones that
                                                   ! was impossible to reroute to coastal points
                                                   ! was impossible to reroute to coastal points
  integer, dimension(size(gfrac,1),size(gfrac, 2)) :: &
       g_moved_to       ! same as above, but gathered from all 
                        ! processors (only used on root processor)
  integer :: pe         ! processor index, for diagnostic gathering
  integer :: ip, npes, root_pe
  logical :: periodic
  integer :: logunit, errunit

  ! size of the global grid
  gnlon = size(gfrac, 1)
  gnlat = size(gfrac, 2)

  allocate(moved_to(gnlon,gnlat),bad_dest(gnlon,gnlat))  ! Want SGI symmetric memory for send buffers

  ! first read and re-grid pointers as with other fields
  if (.not.file_exist("INPUT/river_destination_field")) &
     call error_mesg('init_routing', &
     "Cannot find file INPUT/river_destination_field: provide this file or set use_single_basin = .true. in rivers_nml", FATAL)
! <ERROR MSG="Cannot find file INPUT/river_destination_field: provide this file or set use_single_basin = .true. in rivers_nml" STATUS="FATAL">
!   The namelist settings request the river destination field to be read, but the file 
!   cannot be found. Either provide this file, or change module namelist parameters so 
!   it is not necessary. To do the latter, set use_single_basin = .true. in the namelist rivers_nml.
! </ERROR>
  call mpp_open (unit, 'INPUT/river_destination_field', action = MPP_RDONLY)
  read (unit,*) n_lon_in, n_lat_in, wb_degrees
  read (unit,*) input_format
  allocate ( input_1     (n_lon_in, n_lat_in), &
             input_2     (n_lon_in, n_lat_in), &
             route_input (n_lon_in, n_lat_in), &
             stat = ierr )
  __ASSERT__(ierr==0, "Cannot allocate buffers")
!   <ERROR MSG="Cannot allocate buffers" STATUS="FATAL">
!   </ERROR>

  do j=1,n_lat_in
     read(unit,input_format) (input_1(i,j),i=1,n_lon_in)
  end do
  do j=1,n_lat_in
     read(unit,input_format) (input_2(i,j),i=1,n_lon_in)
  end do
  call close_file(unit)

  ! __NOTE__("Data input DONE")

  route_input = input_1 + n_lon_in*(input_2-1)
  deallocate ( input_1, input_2 )

  ! set the parameters of input grid
  wb = pi*wb_degrees/180.    ! western boundary
  sb = -pi/2                 ! southern boundary
  
  dx = 2.*pi/float(n_lon_in) ! lon size of the cell
  dy = pi/float(n_lat_in)    ! lat size of the cell

  ! compact the range of indices appearing in the routing array, for
  ! quicker processing by regrid_discrete_field. (outer 'if' skips ocean cells,
  ! which just point to themselves, and which are many in number.)
  
  allocate ( route_compact ( n_lon_in,n_lat_in ), &
             route_list    ( n_lon_in*n_lat_in ), &
             stat = ierr)
  __ASSERT__(ierr==0,"Cannot allocate buffers for route compacting")
!   <ERROR MSG="Cannot allocate buffers for route compacting" STATUS="FATAL">
!   </ERROR>

  n_dest        = 0
  route_compact = 0
  route_list    = 0
  do i = 1, n_lon_in
     do j = 1, n_lat_in
        if ( route_input(i,j) /= i+n_lon_in*(j-1) ) then
           match = .FALSE.; n = 0
           do while ( (.not.match) .and. (n < n_dest) )
              n = n+1
              match = (route_input(i,j)==route_list(n))
           enddo
           if (.not.match) then
              n_dest = n_dest + 1
              route_list(n_dest) = route_input(i,j)
              where (route_input==route_input(i,j))  &
                   route_compact = n_dest
           endif
        endif
     enddo
  enddo

  ! write(*,*) "PE=",mpp_pe(),"N_DEST=",n_dest  
   
  allocate ( route (is:ie,js:je), &
             stat = ierr )

  __ASSERT__(ierr==0,"Cannot allocate routing table")
!   <ERROR MSG="Cannot allocate routing table" STATUS="FATAL">    
!   </ERROR>

! NOTE regrid_discrete_field also needs to be changed to use gfrac instead of gland

  route = 0
  call regrid_discrete_field ( &
       route_compact,  wb, sb, dx, dy, &
       glonb(is:ie+1,js:je+1), glatb(is:ie+1,js:je+1), &
       n_dest, (gfrac(is:ie,js:je)>0), route )

  do i = 1, n_lon_in
     do j = 1, n_lat_in
        __ASSERT__(route_compact(i,j)>=0, "Route_compact index below zero")
!   <ERROR MSG="Route_compact index below zero" STATUS="FATAL">
!   </ERROR>
        __ASSERT__(route_compact(i,j)<=n_dest, "Route_compact index above n_dest")
!   <ERROR MSG="Route_compact index above n_dest" STATUS="FATAL">
!   </ERROR>

     enddo
  enddo
  deallocate (route_input,          &
              route_compact         )

  ! __NOTE__("Regridding DONE")

  !-------- then adjust values of pointers to be consistent with new grid
  ! (next code inefficiently fails to take full advantage of new compact array,
  ! but it should do the job for now. we could operate on the short list
  ! of destinations instead of the full array of destinations.)

  ! extract destination info for land cells
  do i = is,ie
     do j = js,je
        if (route(i,j)/=0) then
           __ASSERT__(route(i,j)>=0,      "Route index below zero")
!   <ERROR MSG="Route index below zero" STATUS="FATAL">
!   </ERROR>
           __ASSERT__(route(i,j)<=n_dest, "Route index above n_dest")
!   <ERROR MSG="Route index above n_dest" STATUS="FATAL">
!   </ERROR>
           i_dest(i,j) = mod(route_list(route(i,j))-1,n_lon_in)+1
           j_dest(i,j) = 1+(route_list(route(i,j))-i_dest(i,j))/n_lon_in
        else
           i_dest(i,j) = i
           j_dest(i,j) = j
        endif
     enddo
  enddo
  deallocate (route_list)

  ! __NOTE__("destinations CALCULATED")

  ! convert input-grid indices to model-grid indices (valid values only
  ! for land points)
  periodic = is_latlon(glonb,glatb)
  do i = is,ie
     do j = js,je
        ! For each point in the local domain (nlon, nlat), find the 
        ! corresponding index in the global domain by searching over the 
        ! entire global domain. This is necessary since the destination point 
        ! for a point in a local domain may fall outside of the local domain 
        ! (but not outside the global domain).
        lon = wb + (float(i_dest(i,j))-0.5)*dx ! grid is 0-centered
        lat = sb + (float(j_dest(i,j))-0.5)*dy

        i_dest(i,j) = bisect ( glonb(:,1), lon, periodic=periodic)
        j_dest(i,j) = bisect ( glatb(1,:), lat )
     enddo
  enddo

  ! __NOTE__("Conversion to dest grid DONE")

  ! clean up diagnostic arrays
  bad_dest = 0
  moved_to = 0

  ! clean up where last step made some destinations be land points (and
  ! assign destinations for ocean cells, though this may never be needed).
  ! Use gland rather than land where i_dest, j_dest are subscripts:
  ! Use gnlon, gnlat when checking whether destination points are in valid range.
  ! For ocean points, use is+i-1 rather than i, etc.

  do i = is,ie
     do j = js,je
        if (gfrac(i,j).lt.(1.0-min_water_frac)) then
           ! (i,j) is a (partly) ocean cell, direct own runoff to self
           i_dest(i,j) = i
           j_dest(i,j) = j
        else if ( .not.is_coastal(i_dest(i,j),j_dest(i,j),gfrac) ) then
           ! (i,j) is an all-land cell whose runoff is directed to an all-land 
           !   cell or to a non-coastal all-ocean cell. this is no good, so we 
           !   need to change it
           i_dest_bad = i_dest(i,j)
           j_dest_bad = j_dest(i,j)
           i_dest_new = i_dest_bad
           j_dest_new = j_dest_bad
           ! find closest coastal cell to redirect runoff to. to avoid branching 
           ! out of loop, currently search all possibilities, over-writing from 
           ! most distant to least distant
           do ii = 1, size(di(:))
              next_i = modulo ( i_dest(i,j)+di(ii)+gnlon-1, gnlon ) + 1
              next_j = max ( min ( j_dest(i,j)+dj(ii),gnlat ), 1 )
              if ( is_coastal(next_i,next_j,gfrac) ) then
                 ! found a new destination that is OK
                 i_dest_new = next_i
                 j_dest_new = next_j
              endif
           enddo
           ! change runoff destination to new location
           where (i_dest == i_dest(i,j) .and. j_dest == j_dest(i,j))
              i_dest = i_dest_new
              j_dest = j_dest_new
           endwhere
           ! check if we indeed found a coastal-cell destination for river in the end
           if ( .not. is_coastal(i_dest(i,j), j_dest(i,j), gfrac) ) then
              bad_dest(i_dest(i,j),j_dest(i,j)) = bad_dest(i_dest(i,j),j_dest(i,j))+1
           else
              moved_to(i_dest_bad, j_dest_bad) = j_dest_new*size(gfrac,1)+i_dest_new-1
           endif
        endif
     enddo
  enddo

  deallocate (route)


  ! --- re-routing problem diagnostic section --------------------------------
 
  ! gather bad destination points information from all processors
  call mpp_sum(bad_dest, size(bad_dest(:,:)))

  npes    = mpp_npes()
  root_pe = mpp_root_pe()
  if ( mpp_pe() == root_pe ) then
     ! gather re-routing diagnostics from all processors
     g_moved_to = moved_to
     do ip = 0, npes-1
        pe = ip + root_pe
        if (pe/=mpp_pe()) then 
           ! Force use of "scalar", integer pointer mpp interface
           call mpp_recv(moved_to(1,1),glen=size(moved_to,1)*size(moved_to,2),from_pe=pe)
           where (g_moved_to == 0) g_moved_to = moved_to
        endif
     enddo

     ! print out bad destination points
     logunit = stdlog()
     errunit = stderr()
     do j = 1, size(bad_dest,2)
        do i = 1, size(bad_dest,1)
           if(bad_dest(i,j)/=0) then
              write(logunit,"(a,i4,a,i3,x,i3,a)")             &
                   'ERROR :: RUNOFF SINK FOR ',bad_dest(i,j), &
                   ' CELLS (AT',i,j,') IS NOT NEAR COAST'
              write(errunit,"(a,i4,a,i3,x,i3,a)") &
                   'ERROR :: RUNOFF SINK FOR ',bad_dest(i,j),  &
                   ' CELLS (AT',i,j,') IS NOT NEAR COAST'
!   <ERROR MSG="Runoff sink is not near coast" STATUS="ERROR">
!   </ERROR>
           endif
        enddo
     enddo

     ! print out re-routing information
     do j = 1, size(moved_to,2)
        do i = 1, size(moved_to,1)
           if(g_moved_to(i,j)/=0) &
              write(logunit,"(a,i3,x,i3,a,i3,x,i3,a)")                         &
              'NOTE :: RUNOFF SINK RE-ROUTED FROM LAND OR OCEAN CELL (',        &
              i,j, ') TO COASTAL CELL (',                                       &
              mod(g_moved_to(i,j),size(gfrac,1))+1, g_moved_to(i,j)/size(gfrac,1), ')'
        enddo
     enddo

!   <ERROR MSG="Runoff sink re-routed from land or ocean cell to coastal cell" STATUS="NOTE">
!   </ERROR>

  else
     ! Force use of "scalar", integer pointer mpp interface
     call mpp_send(moved_to(1,1),plen=size(moved_to,1)*size(moved_to,2),to_pe=root_pe)
     call mpp_sync_self()
  endif

  deallocate(moved_to,bad_dest)

  ! __NOTE__("init_routing DONE")

end subroutine init_routing
! </SUBROUTINE>


! <FUNCTION NAME="is_coastal">

!   <OVERVIEW>
!     Tests whether a point is a coastal point.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Returns true if and only if point (i,j) is a coastal point on the map gfrac.
!   </DESCRIPTION>

!   <TEMPLATE>
!     value=is_coastal(i, j, gfrac)
!   </TEMPLATE>

!   <PUBLICROUTINE>
logical function is_coastal(i, j, gfrac)

  integer, intent(in) :: i, j       ! coordinates of the point in question
  real,    intent(in) :: gfrac(:,:) ! fractional area of the land
!   </PUBLICROUTINE>

  ! ---- local vars ---------------------------------------------------------
  integer ip,in, jp,jn ! coordinates of surrounding points

  ! find coordinates of surrounding points
  ip = i-1; if (ip<1) ip = size(gfrac, 1)
  in = i+1; if (in>size(gfrac,1)) in = 1
  jp = max(1, j-1); 
  jn = min(size(gfrac,2), j+1)

  is_coastal = gfrac(i,j) < (1.0-min_water_frac) .and.( &
       gfrac(i ,j ) > min_land_frac .or. &
       gfrac(in,j ) > min_land_frac .or. &
       gfrac(ip,j ) > min_land_frac .or. &
       gfrac(i ,jn) > min_land_frac .or. &
       gfrac(i ,jp) > min_land_frac      )
           
end function is_coastal
! </FUNCTION>


! <SUBROUTINE NAME="init_rivers_diag">

!   <OVERVIEW>
!     Initializes diagnostics for rivers module.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Initializes diagnostics for rivers module.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call init_rivers_diag ( id_lon, id_lat, Time)
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="">
subroutine init_rivers_diag ( id_lon, id_lat, Time)

  integer,         intent(in) :: id_lon  ! ID of land longitude (X) diag axis
  integer,         intent(in) :: id_lat  ! ID of land latitude (Y) diag axis
  type(time_type), intent(in) :: Time    ! current time
!   </PUBLICROUTINE>

  ! regular diagnostic fields
  id_discharge     = register_diag_field (                        &
       mod_name, 'discharge', (/id_lon, id_lat/),                 &
       Time, 'discharge', 'kg/s', missing_value=-999.0            )
  id_discharge_snow = register_diag_field (                       &
       mod_name, 'discharge_snow',  (/id_lon, id_lat/),           &
       Time,  'discharge of snow',   'kg/s', missing_value=-999.0 )
  id_adis_water = register_diag_field (                              &
       mod_name, 'discharge_aw', (/id_lon, id_lat/),  Time,         &
       'water discharge per unit ocean area',                     &
       'kg/(m2 s)', missing_value=-999.0       )
  id_adis_snow = register_diag_field (                              &
       mod_name, 'discharge_snow_aw',  (/id_lon, id_lat/), Time,     &
       'snow discharge per unit ocean area',                      &
       'kg/(m2 s)', missing_value=-999.0 )

  ! static fields
  id_dest = register_static_field ( &
       mod_name, 'dest', (/id_lon, id_lat/), &
       'destination points', 'none', missing_value=0.0 )
  id_basin = register_static_field ( &
       mod_name, 'basin', (/id_lon, id_lat/), &
       'river basins', 'none', missing_value=0.0 )

end subroutine init_rivers_diag
! </SUBROUTINE>


! <DIAGFIELDS>
!   <NETCDF NAME="discharge" UNITS="kg/s">
!     Discharge
!   </NETCDF>
!   <NETCDF NAME="discharge_snow" UNITS="kg/s">
!     Snow Discharge
!   </NETCDF>
!   <NETCDF NAME="discharge_aw" UNITS="kg/(m2 s)">
!     Water discharge per unit ocean area
!   </NETCDF>
!   <NETCDF NAME="discharge_snow_aw" UNITS="kg/(m2 s)">
!     Snow discharge per unit ocean area
!   </NETCDF>
!   <NETCDF NAME="dest" UNITS="dimensionless">
!     Destination points
!   </NETCDF>
!   <NETCDF NAME="basin" UNITS="dimensionless">
!     River basins
!   </NETCDF>
! </DIAGFIELDS>


! <SUBROUTINE NAME="update_rivers_fast">

!   <OVERVIEW>
!     Fast time-scale update of state.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates state of the rivers and runoff on the fast time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_rivers_fast ( Rivers )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_rivers_fast ( Rivers )

  type(rivers_type), intent(inout) :: Rivers ! data to update
!   </PUBLICROUTINE>

  !  increment time
  rivers%Time = increment_time(rivers%Time, int(rivers%dt), 0)

end subroutine update_rivers_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="calc_discharge">

!   <OVERVIEW>
!     Given local runoff field, calculates local portion of discharge.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Given local runoff field, calculates local portion of discharge. Both
!     runoff and discharge fields are of local domain size.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call calc_discharge ( rivers, runoff, discharge )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine calc_discharge ( rivers, runoff, discharge )

  type(rivers_type), intent(in)  :: rivers                           ! data to update
  real,              intent(in)  :: runoff(rivers%is:,rivers%js:)    ! runoff, kg/m2/s
  real,              intent(out) :: discharge(rivers%is:,rivers%js:) ! resulting discharge, kg/s
!   </PUBLICROUTINE>

  ! ---- local vars --------------------------------------------------------
  integer :: i,j  ! iterators
  ! global runoff fields of water and snow, respectively
! real    :: g_runoff(rivers%gnlon,rivers%gnlat) 
  real,allocatable,save :: g_runoff(:,:) 
  logical,save :: first_time=.true.
  
  if(first_time)then
    first_time = .false.
    ALLOCATE(g_runoff(rivers%gnlon,rivers%gnlat))
  endif
  ! clean up the fields and set up local valued
  discharge = 0
  g_runoff  = 0
  g_runoff ( rivers%is:rivers%ie,rivers%js:rivers%je ) = runoff
  ! obtain the global field
  call mpp_sum ( g_runoff, rivers%gnlon*rivers%gnlat )
  ! 
  do i = lbound(g_runoff,1),ubound(g_runoff,1)
     do j = lbound(g_runoff,2),ubound(g_runoff,2)
        if ( rivers%is<=rivers%i_dest(i,j).and.rivers%i_dest(i,j)<=rivers%ie.and.&
             rivers%js<=rivers%j_dest(i,j).and.rivers%j_dest(i,j)<=rivers%je ) &

             discharge(rivers%i_dest(i,j),rivers%j_dest(i,j)) = &
             discharge(rivers%i_dest(i,j),rivers%j_dest(i,j)) + g_runoff(i,j)*rivers%larea(i,j)
     enddo
  enddo
  
end subroutine calc_discharge
! </SUBROUTINE>


! <SUBROUTINE NAME="update_rivers_slow">

!   <OVERVIEW>
!     Slow time-scale update of state.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates state of the rivers and runoff on the slow time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_rivers_slow ( rivers, runoff_w, runoff_s )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_rivers_slow ( rivers, runoff_w, runoff_s )

  type(rivers_type), intent(inout) :: rivers         ! river data to update
  real,              intent(in)    :: runoff_w(:,:)  ! runoff of liquid water, kg/(m2 s)
  real,              intent(in)    :: runoff_s(:,:)  ! snow runoff,  kg/(m2 s)
!   </PUBLICROUTINE>  

  call calc_discharge ( rivers, runoff_w, rivers%discharge )
  call calc_discharge ( rivers, runoff_s, rivers%discharge_snow )
  ! may be it is more effective to exchange both fields at once: less memory-
  ! efficient, though

  call diag_slow ( rivers )
   
end subroutine update_rivers_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="update_rivers_bnd_slow">

!   <OVERVIEW>
!     Updates the river boundary data on the slow time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Updates the river boundary data on the slow time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_rivers_bnd_slow ( rivers, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_rivers_bnd_slow ( rivers, bnd )

  type(rivers_type),    intent(in)    :: rivers  ! river data to update
  type(land_data_type), intent(inout) :: bnd     ! land boundary data
!   </PUBLICROUTINE>

  bnd % discharge      = rivers % discharge
  bnd % discharge_snow = rivers % discharge_snow
  
  where ( rivers%warea > 0.0 )
     ! adjust partial-ocean-cell discharge values for area fraction
     ! and recalculate to kg/(m2 s)
     bnd % discharge      = bnd % discharge      / rivers%warea
     bnd % discharge_snow = bnd % discharge_snow / rivers%warea
  elsewhere
     bnd % discharge      = 0
     bnd % discharge_snow = 0
  endwhere

end subroutine update_rivers_bnd_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="diag_slow">

!   <OVERVIEW>
!     Diagnostic field output of the slow time-scale variables.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Diagnostic field output of the slow time-scale variables.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_slow ( rivers )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_slow ( rivers )

  type(rivers_type), intent(in) :: rivers  ! river data to update
!   </PUBLICROUTINE>

  ! ---- local vars ---------------------------------------------------------
  logical :: used
  real    :: x(rivers%is:rivers%ie, rivers%js:rivers%je)

  if (id_discharge > 0) &
     used = send_data (id_discharge, rivers%discharge, rivers%Time)
  if (id_discharge_snow > 0) &
     used = send_data (id_discharge_snow, rivers%Discharge_snow, rivers%Time)
  if (id_adis_water > 0) then
     x=0
     where ( rivers%warea > 0.0 ) &
          x = rivers%discharge      / rivers%warea
     used = send_data (id_adis_water, x, rivers%Time)
  endif
  if (id_adis_snow  > 0) then
     x=0
     where ( rivers%warea > 0.0 ) &
          x = rivers%discharge_snow / rivers%warea
     used = send_data (id_adis_snow, x, rivers%Time)
  endif

end subroutine diag_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="diag_static">

!   <OVERVIEW>
!     Static diagnostic fields output.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Static diagnostic fields output.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call diag_static ( rivers )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine diag_static ( rivers )

  type(rivers_type), intent(in) :: rivers  ! river data
!   </PUBLICROUTINE>

  ! ---- local vars ---------------------------------------------------------
  logical :: used
  integer :: i,j

  integer, allocatable :: &
       dpoint(:,:), &  ! destination point map
       basin (:,:)     ! basin number

  integer :: nbasins   ! current number of basins
  integer :: id, jd    ! indices of current dest point

  ! + slm Apr 09 2002
  if( id_dest > 0.or.id_basin>0 ) then
     ! allocate global destination points arrays
     allocate( &
          dpoint( lbound(rivers%i_dest,1):ubound(rivers%i_dest,1),    &
                  lbound(rivers%i_dest,2):ubound(rivers%i_dest,2)  ), &
          basin ( lbound(rivers%i_dest,1):ubound(rivers%i_dest,1),    &
                  lbound(rivers%i_dest,2):ubound(rivers%i_dest,2)  )  )

     ! clean up the arrays of destination points, basins, and the number of
     ! river basins
     dpoint  = 0
     basin   = 0
     nbasins = 0

     do j = lbound(rivers%i_dest,2),ubound(rivers%i_dest,2)
        do i = lbound(rivers%i_dest,1),ubound(rivers%i_dest,1)
           id = rivers%i_dest(i,j); jd = rivers%j_dest(i,j)
           if (id /= i.or.jd/=j) then  ! point does not discharge to itself => it is not ocean
              if(dpoint(id,jd).eq.0) then
                 ! we have not encountered this destination point before,
                 ! add a basin
                 nbasins       = nbasins+1
                 dpoint(id,jd) = nbasins
              endif
              basin(i,j) = dpoint(id,jd)
           endif
        enddo
     enddo
     if (id_dest > 0) &
          used = send_data (id_dest, &
          dpoint(rivers%is:rivers%ie, rivers%js:rivers%je)*1.0, rivers%time)
     if (id_basin > 0) &
          used = send_data (id_basin, &
          basin(rivers%is:rivers%ie, rivers%js:rivers%je)*1.0, rivers%time)
     ! clean up after ourselves
     deallocate(dpoint, basin )
  endif
  ! - slm Apr 09 2002

end subroutine diag_static
! </SUBROUTINE>


! <SUBROUTINE NAME="print_error_mesg">

!   <OVERVIEW>
!     Reports error, including file name and line.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Reports error, including file name and line.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call print_error_mesg(mod_name, message, mode, file, line)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine print_error_mesg(mod_name, message, mode, file, line)

  character(len=*), intent(in) :: mod_name   ! module name
  character(len=*), intent(in) :: message    ! error message
  integer,          intent(in) :: mode       ! error mode
  character(len=*), intent(in) :: file       ! file containing error
  integer,          intent(in) :: line       ! line number or error
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  character(len=512) :: mesg
  write(mesg,'("File ",a," Line ",i4.4," :: ",a)')&
       file, line, trim(message)
  call error_mesg(mod_name, mesg, mode)

end subroutine
! </SUBROUTINE>


end module rivers_mod

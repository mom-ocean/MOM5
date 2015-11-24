!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module hgrid_mod
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  !
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>


  !<OVERVIEW>
  ! <TT>hgrid_mod</TT> Generate horizontal grid. The horizontal grid can be conventional lon-lat 
  ! spherical grid or a reprojected rotated tripolar grid (R. Murray, "Explicit generation of 
  ! orthogonal grids for ocean models", 1996, J.Comp.Phys., v. 126, p. 251-273.). 
  !</OVERVIEW>

  !<DESCRIPTION>
  ! There are four subgrids, labeled T (for tracer), C (corner of T), N (north of T) and E (east of T). 
  ! The following schematic describes the grid cell notation.
  !
  !
  !<PRE>
  !                          Ni,j
  !               +----------+-----------+Ci,j
  !               |                      |     
  !               |                      |
  !               |                      |
  !               +          +Ti,j       +Ei,j
  !               |                      |
  !               |                      |
  !               +----------+-----------+
  !
  !                          
  !</PRE>
  !
  ! The grid_spec file would contains all of the following information on each subgrid. 
  ! The following example is for T subgrid. Repeated for E, C, and N subgrids.
  !
  !
  !<PRE>
  !   x_T, y_T           = Geographic location of T-cell center
  !   x_vert_T, y_vert_T = Geographic location of T-cell vertices(each cell has 4 vertices)
  !   area_T             = area of T-cell
  !   angle_T            = Angle clockwise between logical and geographic east of T-cell
  !   ds_00_02_T         = Length of western face of T-cell
  !   ds_20_22_T         = Length of eastern face of T-cell
  !   ds_02_22_T         = Length of northern face of T-cell
  !   ds_00_20_T         = Length of southern face of T-cell
  !   ds_00_01_T         = Distance from southwest corner to western face center of T-cell
  !   ds_01_02_T         = Distance from northwest corner to western face center of T-cell
  !   ds_02_12_T         = Distance from northwest corner to northern face center of T-cell
  !   ds_12_22_T         = Distance from northeast corner to northern face center of T-cell
  !   ds_21_22_T         = Distance from northeast corner to eastern face center of T-cell
  !   ds_20_21_T         = Distance from southeast corner to eastern face center of T-cell
  !   ds_10_20_T         = Distance from southeast corner to southern face center of T-cell
  !   ds_00_10_T         = Distance from southwest corner to southern face center of T-cell
  !   ds_01_11_T         = Distance from center to western face of T-cell
  !   ds_11_12_T         = Distance from center to northern face of T-cell
  !   ds_11_21_T         = Distance from center to eastern face of T-cell
  !   ds_10_11_T         = Distance from center to southern face of T-cell
  !   ds_01_21_T         = width of T-cell
  !   ds_10_12_T         = height of T-cell
  !
  !  Distances between points are described in the following schematics (for T-cell).
  !   
  !
  !
  !               +<----ds_02_12_T---->+<----ds_12_22_T---->+
  !               ^                    ^                    ^
  !               |                    |                    |
  !               |                    |                    |
  !          ds_01_02_T           ds_11_12_T           ds_21_22_T
  !               |                    |                    |
  !               |                    |                    |
  !               v                    v                    v
  !               +<----ds_01_11_T---->+<----ds_11_21_T---->+
  !               ^                    ^                    ^
  !               |                    |                    |
  !               |                    |                    |
  !          ds_00_01_T           ds_10_11_T           ds_20_21_T
  !               |                    |                    |
  !               |                    |                    |
  !               v                    v                    v
  !               +<----ds_00_10_T---->+<----ds_10_20_T---->+
  !
  !
  !
  !               <-------------- ds_02_22_T---------------->
  !             ^ +--------------------+--------------------+ ^
  !             | |                    ^                    | |
  !             | |                    |                    | |
  !             | |                    |                    | |
  !             | |                    |                    | |
  !             | |                    |                    | |
  !     ds_00_02_T|<-------------------+--ds_01_21_T------->| ds_20_22_T              
  !             | |                    |                    | |
  !             | |               ds_10_12_T                | |
  !             | |                    |                    | |
  !             | |                    |                    | |
  !             | |                    |                    | | 
  !             | |                    v                    | |
  !             v +--------------------+--------------------+ v
  !               <-------------- ds_00_20_T---------------->
  !
  !
  !   The other three subgrids (E, N, C subgrids) have similiar name but replacing T.
  !
  ! Axis specifications involve specifying the number of regions for varying
  ! resolution, the bondaries of said regions and the nominal resolution in 
  ! the respective regions.  
  !
  ! For instance, for longitude axis specification:
  !
  !        dx_lon(1) = 4                         dx_lon(2) = 6
  ! |<----|----|----|----|----|----|------|------|------|------|------|------>|
  ! |                              |                                          |
  ! x_lon(1)                    x_lon(2)                                   x_lon(3)
  !
  ! Grid cells are constructed such that
  ! 
  ! dxt(i) = 0.5*(dxu(i-1)+dxu(i))
  !
  !</PRE>
  !
  !</DESCRIPTION>

  use mpp_mod,         only : mpp_pe, mpp_root_pe, mpp_npes, mpp_error, mpp_chksum, FATAL, NOTE
  use mpp_mod,         only : uppercase, lowercase
  use mpp_io_mod,      only : MPP_RDONLY, MPP_ASCII, MPP_MULTI, MPP_SINGLE, MPP_NETCDF   
  use mpp_io_mod,      only : mpp_write_meta, mpp_write, mpp_read, mpp_open, mpp_close 
  use mpp_io_mod,      only : mpp_get_id, mpp_get_info, axistype, fieldtype
  use mpp_io_mod,      only : mpp_get_atts, mpp_get_axes, mpp_get_fields, mpp_get_axis_data
  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_define_domains, mpp_global_field 
  use fms_mod,         only : write_version_number, open_namelist_file, file_exist
  use fms_mod,         only : close_file, check_nml_error, stdout, string
  use constants_mod,   only : radius, pi, RADIAN, epsln
  use axis_utils_mod,  only : nearest_index, lon_in_range, get_axis_cart
  use grids_type_mod,  only : hgrid_data_type, cell_type 
  use grids_util_mod,  only : make_axis, gcell, get_file_unit, write_field_meta 
  use grids_util_mod,  only : write_field_data, set_grid
  implicit none
  private

  integer, parameter :: maxlen=10000,maxbounds=1000

  !------ namelist interface ---------------------------------------------
  !------ specify a spherical grid resolution in latitude, longitude, and depth
  !-----------------------------------------------------------------------
  integer                    :: nxlons = 0     
  integer                    :: nylats = 0     
  real, dimension(maxbounds) :: x_lon, dx_lon, y_lat, dy_lat
  logical                    :: square_grid=.false. 
  logical                    :: extend_square_grid=.false.
  logical                    :: cyclic_x=.true.
  logical                    :: cyclic_y=.false.
  logical                    :: tripolar_grid = .false. 
  real                       :: lat_join = 65.0         ! requested latitude for rotated grid
  logical                    :: read_my_grid = .false.
  character(len=128)         :: my_grid_file = 'my_hgrid'
  logical                    :: debug = .false.
  logical                    :: f_plane             = .false.
  logical                    :: beta_plane          = .false.
  logical                    :: simple_cartesian    = .false.
  real                       :: simple_cartesian_dx = 0
  real                       :: simple_cartesian_dy = 0
  real                       :: f_plane_latitude  =  30.0
  character(len=24) :: lon_axis_t   = 'GRID_X_T'
  character(len=24) :: lat_axis_t   = 'GRID_Y_T'
  character(len=24) :: lon_axis_u   = 'none'
  character(len=24) :: lat_axis_u   = 'none'

  !<NAMELIST NAME="hgrid_nml">
  ! <DATA NAME="nxlons" TYPE="integer">
  ! number of zonal regions for varying resolution
  ! </DATA>
  ! <DATA NAME="nylats" TYPE="integer">
  ! number of latitude regions for varying resolution
  ! </DATA>
  ! <DATA NAME="x_lon" TYPE="real" DIM="(nxlons)" UNITS="degrees">
  ! boundaries for defining zonal regions of varying resolution. When tripolar_grid 
  ! is .true., x_lon also defines the longitude of the two new poles. 
  ! lon_start = x_lon(1) and lon_end = x_lon(nxlons) are longitude of the two new 
  ! poles. In this case, the program will ignore the value x_lon(2:nxlons-1)
  ! and set grid resolution to dx_lon(1). When tripolar_grid is true, you
  ! need to be careful about your choice of x_lon, because there might be ocean
  ! at the grid singularity. The recommended choice of x_lon is x_lon = -280,80,
  ! this will put the singularity over land.
  ! </DATA>
  ! <DATA NAME="dx_lon" TYPE="real" DIM="(nxlons)" UNITS="degrees">
  ! nominal resolution of zonal regions
  ! </DATA>
  ! <DATA NAME="cyclic_x" TYPE="logical">
  ! True if grid is connected in i-direction 
  ! </DATA>
  ! <DATA NAME="cyclic_y" TYPE="logical">
  ! True if grid is connected in j-direction 
  ! </DATA>
  ! <DATA NAME="y_lat" TYPE="real" DIM="(nylats)" UNITS="degrees">
  ! boundaries for defining meridional regions of varying resolution
  ! </DATA>
  ! <DATA NAME="dy_lat" TYPE="real" DIM="(nxlons)" UNITS="degrees">
  ! nominal resolution of meridional regions
  ! </DATA>
  !<DATA NAME="tripolar_grid" TYPE="logical">
  ! convert portion of spherical grid north of lat_join to a bipolar rotated grid
  !</DATA>
  !<DATA NAME="square_grid" TYPE="logical">
  ! latitudinal grid spacing matches convergence of meridians
  !</DATA>
  !<DATA NAME="extend_square_grid" TYPE="logical">
  ! extend square grid to poles
  !</DATA>
  !<DATA NAME="lat_join" TYPE="real">
  ! requested latitude for joining spherical and rotated bipolar grid
  !</DATA>
  !<DATA NAME="read_my_grid" TYPE="logical">
  ! read ASCII grid information for supplying user-defined grids. 
  !</DATA>
  !<DATA NAME="my_grid_file" TYPE="character(len=128)">
  ! Name of ASCII user grid file
  !</DATA>
  !<DATA NAME="f_plane" TYPE="logical">
  !  For setting geometric fractors according to f-plane. 
  !</DATA> 
  !<DATA NAME="beta_plane" TYPE="logical">
  !  For setting geometric fractors according to beta plane. 
  !</DATA> 
  !<DATA NAME="f_plane_latitude" TYPE="real">
  !  Central latitude to define f_plane and beta_plane.
  !</DATA> 
  !<DATA NAME="simple_cartesian" TYPE="logical">
  !  For setting simple cartesian grid. When set true, simple_cartesian_dx 
  !  and simple_cartesian_dy need to be set. The grid box length in x-direction 
  !  will be uniform to be simple_cartesian_dx. The grid box length in y-direction 
  !  will be uniform to be simple_cartesian_dy.
  !</DATA> 
  !<DATA NAME="simple_cartesian_dx" TYPE="real">
  !  uniform grid length in x-direction ( units is meter).
  !</DATA> 
  !<DATA NAME="debug" TYPE="logical">
  !  uniform grid length in y-direction ( units is meter).
  !</DATA>
  !</NAMELIST>

  namelist /hgrid_nml/ nxlons, nylats, x_lon, dx_lon, y_lat, dy_lat, square_grid,      &
       extend_square_grid, cyclic_x, cyclic_y, tripolar_grid, lat_join, read_my_grid,  &
       my_grid_file, debug, f_plane, beta_plane, f_plane_latitude,                     &
       lon_axis_t, lat_axis_t, lon_axis_u, lat_axis_u, simple_cartesian,               &
       simple_cartesian_dx, simple_cartesian_dy

  !-----------------------------------------------------------------------
  !--------private data---------------------------------------------------
  integer                         :: ni                         ! number of grid cells in zonal direction
  integer                         :: nj                         ! number of grid cells in meridional direction
  integer                         :: isc, iec, jsc, jec         ! compute domain
  integer                         :: isd, ied, jsd, jed         ! data domain
  real, dimension(:), allocatable :: xt0, xu0, yt0, yu0
  real                            :: x_period = 0.0             ! period in i-direction
  real                            :: y_period = 0.0             ! period in j-direction
  real                            :: lon_start=0.0, lon_end=0.0 ! beginning and ending longitudes for tripolar grid
  real                            :: lon_bpeq, lon_bpnp, lon_bpsp, lam0, rp
  real                            :: D2R, hpi, tpi
  logical                         :: module_is_initialized = .false.
  type(domain2d),save             :: Domain

  !---------version information-------------------------------------------
  character(len=128) :: version = '$Id: hgrid.f90,v 14.0 2007/03/15 22:46:29 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 

  !---------public interface----------------------------------------------
  public :: generate_hgrid, hgrid_init, hgrid_end, write_hgrid_global_meta
  public :: write_hgrid_field_meta, write_hgrid_data

contains

  !#######################################################################
  ! <SUBROUTINE NAME="hgrid_init" >
  !   <OVERVIEW>
  !    Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    Read namelist, write out version and namelist informaiton, generate longitude
  !     and latitude resolution.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call hgrid_init ( )
  !   </TEMPLATE>
  ! </SUBROUTINE>

  subroutine hgrid_init 

    !--- local variables -------------------------------------------------
    integer                     :: i, j, pe, unit, ierr, io, npes
    integer                     :: ndim, nvar, natt, ntime, len
    character(len=128)          :: txt, name
    real, allocatable           :: data2d(:,:), data3d(:,:,:)
    type(axistype), allocatable :: axes(:)
    type(fieldtype),allocatable :: flds(:)
    logical                     :: found_x_t, found_y_t, found_x_vert_t, found_y_vert_t
    logical                     :: found_x_u, found_y_u

    D2R = PI/180.
    hpi = PI*0.5
    tpi = PI*2.0

    x_lon(:)    = 0.0
    dx_lon(:)   = 0.0
    y_lat(:)    = 0.0
    dy_lat(:)   = 0.0

    !---- read namelist --------------------------------------------------
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
       read  (unit, nml=hgrid_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'hgrid_nml')  ! also initializes nml error codes
    enddo
10  call close_file (unit)

    !--- write version info and namelist to logfile ----------------------

    call write_version_number(version, tagname)
    write (stdout(), nml=hgrid_nml)

    module_is_initialized = .true.

    ! --- namelist check -----------------------------------------------
    if(.not. read_my_grid) then
       if(nxlons .lt. 2 .or. nylats .lt. 2) call mpp_error(FATAL, &
            'hgrid_mod: nml "nxlons"= '//trim(string(nxlons))//' and "nylats"= '// &
             trim(string(nylats))//' should be no less than 2')
       if (nxlons .gt. maxbounds )  &
            call mpp_error(FATAL,'hgrid_mod: nml "nxlons"= '//trim(string(nxlons))// &
              ' is greater than "maxbounds"= '//trim(string(maxbounds)) )
       if (nylats .gt. maxbounds )  &
            call mpp_error(FATAL,'hgrid_mod: nml "nylats"= '//trim(string(nylats))// &
              ' is greater than "maxbounds"= '//trim(string(maxbounds)) )
       if ( y_lat(nylats) - y_lat(1) > 180.0) write(stdout(),*)   &
            "=>Warning: Latitudinal domain width exceeds 180 deg.", &
            "Restricting last bounds to ",dy_lat(nylats) +180.0

       !--- if cyclic condition, set the period
       if(cyclic_x) then
          x_period = x_lon(nxlons) - x_lon(1)
          if(abs(x_period - 360) < epsln ) x_period = 360.0
          write(stdout(),*) "x-boundary is cyclic with period = ", x_period
          if(dx_lon(nxlons) .NE. dx_lon(1)) call mpp_error(FATAL, &
               'hgrid_mod: dx_lon(1) must equal dx_lon(last), modify nml "dx_lon" ' )
       end if
       if(cyclic_y) then
          y_period = y_lat(nylats) - y_lat(1)
          write(stdout(),*) "y-boundary is cyclic with period = ", y_period
          if(dy_lat(nylats) .NE. dy_lat(1)) call mpp_error(FATAL, &
               'hgrid_mod: dy_lat(1) must equal dy_lat(last), modify nml "dy_lat" ' )
       end if

       !--- for tripolar grid, nxlons should be 2.
       if(tripolar_grid) then
          !--- tripolar grid can not be cyclic_y
          if(cyclic_y) call mpp_error(FATAL, "hgrid_mod:nml cyclic_y and tripolar can not both be true")
          if(nxlons .NE. 2) call mpp_error(FATAL, "hgrid_mod: nml nxlons should be set to 2 for tripolar grid")
          !-- for tripolar grid, the period should be 360.
          if( x_period .NE. 360.0 ) call mpp_error(FATAL, &
               "hgrid_mod: for tripolar grid, period in i-direction should be 360")
          if(x_lon(1) .ne. -280 .and. x_lon(2) .ne. 80 ) then
             write(stdout(),*) ' '
             write(stdout(),*)' WARNING from hgrid_mod: the grid is tripolar grid, but ', &
                  'the longitude of the two poles are not -280 and 80, you might ', &
                  'put the singularity over ocean with your choice of longitude of ', &
                  'the two poles: ', x_lon(1), ', ', x_lon(2)
             write(stdout(),*) ' '
          endif
          if ( y_lat(nylats) < 90.) call mpp_error(FATAL, &
             "hgrid_mod: the last latitude bound is less than 90 degree and is tripolar grid")

          if(dx_lon(1) .NE. dx_lon(2)) call mpp_error(FATAL,  &
               'hgrid_mod: dx_lon(1) should equal dx_lon(2) for tripolar grid ')
          if(mod(x_period,dx_lon(1)) .ne. 0) call mpp_error(FATAL,  &
               'hgrid_mod: non-integral number of longitude cells. modify nml "dx_lon" ')
       end if

       !--- at most one of f_plane and beta_plane can be true
       if(f_plane .and. beta_plane) call mpp_error(FATAL, "hgrid_mod: f_plane and beta_plane can not both be true")
       if(f_plane .or. beta_plane) then
          if(f_plane_latitude > 90 .or. f_plane_latitude < -90.) call mpp_error(FATAL, &
                "hgrid_mod: nml f_plane_latitude should be between -90 and 90.")
          if(f_plane_latitude>y_lat(nylats) .or. f_plane_latitude<y_lat(1)) then
              write(stdout(),*)" Warning from hgrid_mod: nml f_plane_latitude is not inside the latitude range of the grid"
          end if
          write(stdout(),*)" setting geometric factor according to f-plane with f_plane_latitude = ", f_plane_latitude
          !--- for f_plane or beta_plane, can not be tripolar grid.
          if(tripolar_grid) call mpp_error(FATAL, "hgrid_mod: for f_plane or beta_plane, can not be tripolar grid.")
       end if
    end if

    !--- when simple_cartesian is set to true, f_plane, beta_plane and tripolar_grid all should be false
    if(simple_cartesian) then
       if(f_plane .or. beta_plane .or. tripolar_grid) call mpp_error(FATAL,  &
           "hgrid_mod: when simple_cartesian is set to true, _plane, beta_plane and tripolar_grid all should be false")
       if(simple_cartesian_dx == 0 .or. simple_cartesian_dy == 0 ) call mpp_error(FATAL,  &
           "hgrid_mod: both imple_cartesian_dx and simple_cartesian_dy need to be set")
    end if

    allocate (xt0(0:maxlen), xu0(0:maxlen), yt0(0:maxlen), yu0(0:maxlen))

    if (read_my_grid) then

       if (tripolar_grid) call mpp_error(FATAL,&
               'hgrid_mod: nml "read_my_grid"= true contradicts with nml "tripolar_grid"=true')
       if (f_plane .or. beta_plane ) call mpp_error(FATAL,&
               'hgrid_mod: nml "read_my_grid"= true contradicts with nml "f_plane"/"beta_plane" =true')
       !
       ! read grid information from an ASCII file.
       ! this is an alternative to namelist grid specification.
       ! For instance, to generate an R30 ocean grid.
       ! we are expecting t-cell centers and uv-cell centers (upper right).
       ! 
       ni=0;nj=0
       if(.not. file_exist(trim(my_grid_file))) &
            call mpp_error(FATAL,'hgrid_mod: file '//trim(my_grid_file)//' does not exist')
       !--- my_grid_file can be in netcdf format or ascii format. The grid should be 
       !--- spherical grid.
       len = len_trim(my_grid_file)
       if(my_grid_file(len-2:len) == '.nc') then
          call mpp_open(unit,trim(my_grid_file),action=MPP_RDONLY,form=MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
          call mpp_get_info(unit, ndim, nvar, natt, ntime)
          allocate(axes(ndim))
          call mpp_get_axes(unit, axes)
          found_x_t = .false.; found_x_vert_t = .false.
          found_y_t = .false.; found_y_vert_t = .false.
          found_x_u = .false.; found_y_u = .false.
          do i=1, ndim
             call mpp_get_atts(axes(i),name=name,len=len)  
             if (trim(lowercase(name)) .eq. trim(lowercase(lat_axis_t))) then
               found_y_t = .true.
               nj = len
               call mpp_get_axis_data(axes(i), yt0(1:nj) )
             endif  
             if (trim(lowercase(name)) .eq. trim(lowercase(lon_axis_t))) then
               found_x_t = .true.
               ni = len
               call mpp_get_axis_data(axes(i), xt0(1:ni) )
             endif  
             if (trim(lowercase(name)) .eq. trim(lowercase(lat_axis_u))) then
               found_y_u = .true.
               nj = len
               call mpp_get_axis_data(axes(i), yu0(1:nj) )
               yu0(0) = 2.0*yu0(1)-yu0(2)
             endif  
             if (trim(lowercase(name)) .eq. trim(lowercase(lon_axis_u))) then
               found_x_u = .true.
               ni = len
               call mpp_get_axis_data(axes(i), xu0(1:ni) )
               xu0(0) = 2.0*xu0(1)-xu0(2)
             endif  
          enddo
          if (ni == 0 .or. nj == 0 ) &
             call mpp_error(FATAL,'Error reading grid size from grid file '//trim(my_grid_file) )
          allocate(flds(nvar))
          call mpp_get_fields(unit, flds)
          do i=1, nvar
             call mpp_get_atts(flds(i),name=name)
             select case(trim(uppercase(name)))
             case('X_T')
                found_x_t = .true.
                allocate(data2d(ni,nj))
                call mpp_read(unit,flds(i), data2d(:,:))
                xt0(1:ni) = data2d(1:ni,1)
                deallocate(data2d)
             case('Y_T')
                found_y_t = .true.
                allocate(data2d(ni,nj))
                call mpp_read(unit,flds(i), data2d(:,:))
                yt0(1:nj) = data2d(1,1:nj)
                deallocate(data2d)
             case('X_U')
                found_x_u = .true.
                allocate(data2d(ni,nj))
                call mpp_read(unit,flds(i), data2d(:,:))
                xu0(1:ni) = data2d(1:ni,1)
                xu0(0) = 2.0*xu0(1)-xu0(2)
                deallocate(data2d)
             case('Y_U')
                found_y_u = .true.
                allocate(data2d(ni,nj))
                call mpp_read(unit,flds(i), data2d(:,:))
                yu0(1:nj) = data2d(1,1:nj)
                yu0(0) = 2.0*yu0(1)-yu0(2)
                deallocate(data2d)
             case('X_VERT_T')
                found_x_vert_t = .true.
                allocate(data3d(ni,nj,4))
                call mpp_read(unit,flds(i), data3d)
                xu0(0)    = data3d(1,1,1)
                xu0(1:ni) = data3d(1:ni,1,3)
                deallocate(data3d)
             case('Y_VERT_T')
                found_y_vert_t = .true.
                allocate(data3d(ni,nj,4))
                call mpp_read(unit,flds(i), data3d)
                yu0(0)    = data3d(1,1,1)
                yu0(1:nj) = data3d(1,1:nj,3)
                deallocate(data3d)
             end select
          enddo
          if(.not. found_x_t) call mpp_error(FATAL,'field x_T does not exist in file '//trim(my_grid_file) )
          if(.not. found_y_t) call mpp_error(FATAL,'field y_T does not exist in file '//trim(my_grid_file) )
          if(.not. found_x_vert_t .and. .not. found_x_u) call mpp_error(FATAL, &
             'neither field x_vert_U nor field x_U do not exist in file '//trim(my_grid_file) )
          if(.not. found_y_vert_t .and. .not. found_y_u) call mpp_error(FATAL, &
             'neither field y_vert_U nor field y_U do not exist in file '//trim(my_grid_file) )
          if(found_x_vert_t .and. found_x_u) call mpp_error(FATAL, &
             'field x_vert_U and field x_U exist in file '//trim(my_grid_file) )
          if(found_y_vert_t .and. found_y_u) call mpp_error(FATAL, &
             'field y_vert_U adn field y_U exist in file '//trim(my_grid_file) )
          deallocate(axes, flds)
       else
          call mpp_open(unit,trim(my_grid_file),action=MPP_RDONLY,form=MPP_ASCII,threading=MPP_MULTI,fileset=MPP_SINGLE)
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) ni
          read(unit,*) xt0(1:ni)
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) nj
          read(unit,*) yt0(1:nj)     
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) 
          read(unit,*) xu0(0:ni)
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) 
          read(unit,*) yu0(0:nj)     
          call mpp_close(unit)

          if (ni == 0 .or. nj == 0 ) call mpp_error(FATAL,'error reading my_grid file')
       endif
       if (cyclic_x) then  
          x_period = xu0(ni) - xu0(0)
          xt0(0) = xt0(ni)-x_period
          xt0(ni+1) = xt0(1)+x_period
          xu0(ni+1) = xu0(0)+x_period
       else
          xt0(0) = 2.0*xt0(1)-xt0(2)
          xt0(ni+1) = 2.0*xt0(ni)-xt0(ni-1)
          xu0(ni+1) = 2.0*xu0(ni)-xu0(ni-1)
       endif
       if (cyclic_y) then  
          y_period    = yu0(nj) - yu0(0)
          yt0(0)    = yt0(nj)-y_period
          yt0(nj+1) = yt0(1)+y_period
          yu0(nj+1) = yu0(0)+y_period
       else
          yt0(0) = 2.0*yt0(1)-yt0(2)
          yt0(nj+1) = 2.0*yt0(nj)-yt0(nj-1)
          yu0(nj+1) = 2.0*yu0(nj)-yu0(nj-1)
       end if
    else
       write (stdout(),'(//,36x,a,/)') 'H G R I D   G E N E R A T I O N'

       ! Tile full longitudinal domain with i=1,ni cells. 

       write (stdout(),'(/a)')'Generating longitudinal resolution for the computational domain:'
       call make_axis('X',maxlen, nxlons,x_lon,dx_lon,xu0,xt0,ni,square_grid, extend_square_grid)
       write (stdout(),'(/a)')'Generating latitudinal resolution for the computational domain:'
       call make_axis('Y',maxlen, nylats,y_lat,dy_lat,yu0,yt0,nj,square_grid, extend_square_grid )
    endif

    ! Print all grid coordinates

    write (stdout(),'(//,40x,a,//,a,g14.7,a/,a,g14.7,a/,a,g14.7,a/)') &
         ' Grid Point Coordinate details: ',' The western edge of the T-cell at i=1 is at:',&
         xu0(0),' (deg) longitude',' The southern edge of the T-cell at j=1 is at:',yu0(0),'(deg) latitude'

    write (stdout(),'(a,i4,a,g14.7,a/)') ' The eastern edge of the T-cell at i=',ni,' is at:',&
         xu0(ni),' (deg) longitude',' The northern edge of the T-cell at j=',nj,' is at:',&
         yu0(nj),' (deg) latitude'

    write (stdout(),'(/,a,g14.7,a/,a,g14.7,a/,a,g14.7,a/)') ' The western edge of the C-cell at i=1 is at:', xt0(1),&
         '(deg) longitude',' The southern edge of the C cell at j=1 is at:', yt0(1),&
         '(deg) latitude'

    write (stdout(),'(a,i4,a,g14.7,a/)') ' The eastern edge of the C-cell at i=',ni,' is at:',&
         xt0(ni+1),'(deg) longitude',' The northern edge of the C-cell at j=',nj,' is at:',&
         yt0(nj+1),'(deg) latitude'

    write (stdout(),9105) nj
    write (stdout(),9001) (yt0(j),j=1,nj)
    write (stdout(),9106) nj
    write (stdout(),9001) (yu0(j),j=1,nj)
    write (stdout(),9107) ni
    write (stdout(),9001) (xt0(i),i=1,ni)
    write (stdout(),9108) ni
    write (stdout(),9001) (xu0(i),i=1,ni)
    write (stdout(),'(/)')
    ! Compute a grid checksum
    if(debug) then     
       pe = mpp_pe()
       write (stdout(),*) 'Grid checksum = ',mpp_chksum(xt0(1:ni),(/pe/)) + mpp_chksum(yt0(1:nj),(/pe/)) + &
            mpp_chksum(xu0(1:ni),(/pe/)) + mpp_chksum(yu0(1:nj),(/pe/)) 
       write (stdout(),'(/)')
    endif

    !--- domain decomposition --------------------------------------------
    npes = mpp_npes()
    call mpp_define_domains((/1,ni,1,nj/),(/1,npes/), Domain, xhalo=1, yhalo=1)
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)

    return

9001 format (1x,10f10.4)
9002 format (1x,10f10.2)
9101 format (/,  a,i4,' in units of ',a,' as follows:')
9105 format (/,' Latitude of T points (deg): yt(j) j=1,',i4)
9106 format (/,' Latitude of C points (deg): yu(j) j=1,',i4)
9107 format (/,' Longitude of T points (deg): xt(i) i=1,',i4)
9108 format (/,' Longitude of C points (deg): xu(i) i=1,',i4)

  end subroutine hgrid_init

  !#######################################################################
  ! <SUBROUTINE NAME="generate_hgrid" >
  !   <OVERVIEW>
  !    Generate horizontal grid.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    Define geographical locations of center and vertices of 
  !    T, C, E, N-cell and also calculate the area, orientation, cell size, face lengths, 
  !    half face lengths and center to face size of each T, E, C, N cell.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call generate_hgrid (Hgrid)
  !   </TEMPLATE>
  !   <INOUT NAME="Hgrid" TYPE="hgrid_data_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </INOUT>
  ! </SUBROUTINE>
  subroutine generate_hgrid (Hgrid)
    type(hgrid_data_type), intent(inout) :: Hgrid
    integer :: i, j

    Hgrid%ni            = ni
    Hgrid%nj            = nj
    Hgrid%Domain        = Domain
    Hgrid%tripolar_grid = tripolar_grid
    Hgrid%cyclic_x      = cyclic_x
    Hgrid%cyclic_y      = cyclic_y

    call allocate_hgrid(Hgrid)

    do j=0,nj+1
       Hgrid%T%x(0:ni+1,j) = xt0(0:ni+1)
    enddo
    do i=0,ni+1
       Hgrid%T%y(i,0:nj+1) = yt0(0:nj+1)
       if (Hgrid%T%y(i,nj+1) .gt. 90.) then
          Hgrid%T%y(i,nj+1) = 180.-Hgrid%T%y(i,nj+1)
       endif
    enddo
    do j=0,nj+1
       Hgrid%C%x(0:ni+1,j) = xu0(0:ni+1)
    enddo
    do i=0,ni+1
       Hgrid%C%y(i,0:nj+1) = yu0(0:nj+1)
       if (Hgrid%C%y(i,nj+1) .gt. 90.) then
          Hgrid%C%y(i,nj+1) = 180.-Hgrid%C%y(i,nj+1)     
       endif
    enddo

    Hgrid%E%x(:,:) = Hgrid%C%x(:,:)
    Hgrid%E%y(:,:) = Hgrid%T%y(:,:)
    Hgrid%N%x(:,:) = Hgrid%T%x(:,:)
    Hgrid%N%y(:,:) = Hgrid%C%y(:,:)  

    ! write spherical grid information to Hgrid (center of T,C,E and N-cell) and 
    ! set up the variables to calculate the half length of each cell. 
    if(f_plane .or. beta_plane) then
       call f_plane_grid_init(Hgrid)
    else if(simple_cartesian) then
       call simple_cartesian_grid_init(Hgrid)
    else
       call spherical_grid_init(Hgrid)
    end if

    !--- call set_grid for setting up io
    call set_grid(xt0(1:ni), yt0(1:nj), xu0(1:ni),yu0(1:nj))

    ! if tripolar_grid, then overload all the grid information above the lat_join,
    ! and the variables to calculate half length of cell.
    if (tripolar_grid) call tripolar_grid_init(Hgrid)

    ! write grid information of vertices of each cell to the Hgrid.
    call define_vertices(Hgrid)

    if(debug) then
       call cell_chksum(Hgrid%T, 'T')
       call cell_chksum(Hgrid%E, 'E')
       call cell_chksum(Hgrid%N, 'N')
       call cell_chksum(Hgrid%C, 'C')
    endif

    return

  end subroutine generate_hgrid

  !#######################################################################
  ! <SUBROUTINE NAME="write_hgrid_data">
  !   <OVERVIEW>
  !     write the Hgrid data to netcdf file
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_hgrid_data (unit,Hgrid)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Hgrid" TYPE="hgrid_data_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_hgrid_data (file,Hgrid)
    character(len=*),      intent(in) :: file
    type(hgrid_data_type), intent(in) :: Hgrid

    !--- write out Hgrid T-cell data -------------------------------
    call write_cell_data(file, Hgrid%T, 'T') 

    !--- write out Hgrid E-cell data -------------------------------
    call write_cell_data(file, Hgrid%E, 'E') 

    !--- write out Hgrid N-cell data -------------------------------
    call write_cell_data(file, Hgrid%N, 'N') 

    !--- write out Hgrid C-cell data -------------------------------
    call write_cell_data(file, Hgrid%C, 'C') 

    return

  end subroutine write_hgrid_data

  !#####################################################################

  subroutine write_hgrid_global_meta(file)
    character(len=*), intent(in) :: file
    integer                      :: unit

    if(mpp_pe() .ne. mpp_root_pe() ) return

    unit = get_file_unit(file)

    call mpp_write_meta(unit,'xname', cval='longitude')
    call mpp_write_meta(unit,'yname', cval='latitude')
    call mpp_write_meta(unit,'vertex_convention', cval='SWCCW')

    if (tripolar_grid) then
       call mpp_write_meta(unit,'join_lat',rval=lat_join)
       call mpp_write_meta(unit, 'y_boundary_type',cval='fold_north_edge')
    else if(cyclic_y) then
       call mpp_write_meta(unit, 'y_boundary_type',cval='cyclic')
    else
       call mpp_write_meta(unit,'y_boundary_type',cval='solid_walls')
    endif

    if (cyclic_x) then
       call mpp_write_meta(unit,'x_boundary_type',cval='cyclic')
    else
       call mpp_write_meta(unit,'x_boundary_type',cval='solid_walls')
    endif

    if(f_plane) call mpp_write_meta(unit,'f_plane', cval='y')
    if(beta_plane) call mpp_write_meta(unit,'beta_plane', cval='y')
    if(f_plane .or. beta_plane) call mpp_write_meta(unit,'f_plane_latitude',rval=f_plane_latitude)

  end subroutine write_hgrid_global_meta

  !#######################################################################
  ! <SUBROUTINE NAME="write_hgrid_meta">

  !   <OVERVIEW>
  !     Write out horizontal grid meta data.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_hgrid_meta(unit, Hgrid, axis_x, axis_y)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Hgrid" TYPE="hgrid_data_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </IN>
  !   <OUT NAME="axis_x, axis_y" TYPE="type(axistype), optional">
  !     axis of T-cell center
  !   </OUT>
  ! </SUBROUTINE>

  subroutine write_hgrid_field_meta(file)
    character(len=*), intent(in) :: file

    if(mpp_pe() .ne. mpp_root_pe() ) return

    !--- write out T-cell field meta --------------------------------
    call write_cell_meta(file, 'T' )

    !--- write out E-cell field meta --------------------------------
    call write_cell_meta(file, 'E' )

    !--- write out N-cell field meta --------------------------------
    call write_cell_meta(file, 'N' )

    !--- write out C-cell field meta --------------------------------
    call write_cell_meta(file, 'C' )

    return

  end subroutine write_hgrid_field_meta

  !#######################################################################
  ! <SUBROUTINE NAME="hgrid_end">

  !   <OVERVIEW>
  !     Destruction routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "hgrid_data_type" variables.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call hgrid_end ( Hgrid )
  !   </TEMPLATE>
  !   <INOUT NAME="Hgrid" TYPE="hgrid_data_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </INOUT>
  ! </SUBROUTINE>

  subroutine hgrid_end(Hgrid)

    type(hgrid_data_type), intent(inout) :: Hgrid
    !--- release memory of module variables ------------------------------

    deallocate(xt0, xu0, yt0, yu0)
    !--- release memory of Hgrid T-cell data -------------------------------
    call deallocate_hgrid_cell(Hgrid%T)

    !--- release memory of Hgrid E-cell data -------------------------------
    call deallocate_hgrid_cell(Hgrid%E)

    !--- release memory of Hgrid N-cell data -------------------------------
    call deallocate_hgrid_cell(Hgrid%N)

    !--- release memory of Hgrid C-cell data -------------------------------
    call deallocate_hgrid_cell(Hgrid%C)

    module_is_initialized = .false.

    return

  end subroutine hgrid_end

  !#######################################################################
  !allocate memory to the hgrid data type
  subroutine allocate_hgrid(Hgrid)
    type(hgrid_data_type), intent(inout) :: Hgrid

    !--- allocate memory to Hgrid T-cell data variables ------------ 
    call allocate_hgrid_cell(Hgrid%T)
    !--- allocate memory to Hgrid E-cell data variables ------------ 
    call allocate_hgrid_cell(Hgrid%E)
    !--- allocate memory to Hgrid N-cell data variables ------------ 
    call allocate_hgrid_cell(Hgrid%N)
    !--- allocate memory to Hgrid C-cell data variables ------------ 
    call allocate_hgrid_cell(Hgrid%C)

    return
  end subroutine allocate_hgrid

  !#######################################################################
  !--- generate spherical grid
  subroutine spherical_grid_init(Hgrid)

    type(hgrid_data_type), intent(inout) :: Hgrid

    !---------------------------------------------------------------------
    !--- calculate face length, half length area and rotation angle of T-cell
    call spherical_cell_length_area(Hgrid%T, Hgrid%T%x(isc:iec,jsc:jec), Hgrid%T%y(isc:iec,jsc:jec),             &
         Hgrid%E%x(isd:iec,jsc:jec), Hgrid%E%y(isd:iec,jsc:jec), Hgrid%N%x(isc:iec,jsd:jec),     &
         Hgrid%N%y(isc:iec,jsd:jec), Hgrid%C%x(isd:iec,jsd:jec), Hgrid%C%y(isd:iec,jsd:jec)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of E-cell

    call spherical_cell_length_area(Hgrid%E, Hgrid%E%x(isc:iec,jsc:jec), Hgrid%E%y(isc:iec,jsc:jec),             &
         Hgrid%T%x(isc:ied,jsc:jec), Hgrid%T%y(isc:ied,jsc:jec), Hgrid%C%x(isc:iec,jsd:jec),     &
         Hgrid%C%y(isc:iec,jsd:jec), Hgrid%N%x(isc:ied,jsd:jec), Hgrid%C%y(isc:ied,jsd:jec)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of N-cell
    call spherical_cell_length_area(Hgrid%N, Hgrid%N%x(isc:iec,jsc:jec), Hgrid%N%y(isc:iec,jsc:jec),             &
         Hgrid%C%x(isd:iec,jsc:jec), Hgrid%C%y(isd:iec,jsc:jec), Hgrid%T%x(isc:iec,jsc:jed),     &
         Hgrid%T%y(isc:iec,jsc:jed), Hgrid%E%x(isd:iec,jsc:jed), Hgrid%E%y(isd:iec,jsc:jed)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of C-cell
    call spherical_cell_length_area(Hgrid%C, Hgrid%C%x(isc:iec,jsc:jec), Hgrid%C%y(isc:iec,jsc:jec),             &
         Hgrid%N%x(isc:ied,jsc:jec), Hgrid%N%y(isc:ied,jsc:jec), Hgrid%E%x(isc:iec,jsc:jed),     &
         Hgrid%E%y(isc:iec,jsc:jed), Hgrid%T%x(isc:ied,jsc:jed), Hgrid%T%y(isc:ied,jsc:jed)  )

    !--- for spherical grid, rotation angle always is 0 ------------------
    Hgrid%T%angle = 0.0
    Hgrid%E%angle = 0.0
    Hgrid%N%angle = 0.0
    Hgrid%C%angle = 0.0  

    return
  end subroutine spherical_grid_init


  !#######################################################################
  !--- generate f_plane/beta_plane grid
  subroutine f_plane_grid_init(Hgrid)
    type(hgrid_data_type), intent(inout) :: Hgrid

    !---------------------------------------------------------------------
    !--- calculate face length, half length area and rotation angle of T-cell
    call f_plane_cell_length_area(Hgrid%T, Hgrid%T%x(isc:iec,jsc:jec), Hgrid%T%y(isc:iec,jsc:jec),   &
         Hgrid%E%x(isd:iec,jsc:jec), Hgrid%E%y(isd:iec,jsc:jec), Hgrid%N%x(isc:iec,jsd:jec),         &
         Hgrid%N%y(isc:iec,jsd:jec), Hgrid%C%x(isd:iec,jsd:jec), Hgrid%C%y(isd:iec,jsd:jec)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of E-cell

    call f_plane_cell_length_area(Hgrid%E, Hgrid%E%x(isc:iec,jsc:jec), Hgrid%E%y(isc:iec,jsc:jec),   &
         Hgrid%T%x(isc:ied,jsc:jec), Hgrid%T%y(isc:ied,jsc:jec), Hgrid%C%x(isc:iec,jsd:jec),         &
         Hgrid%C%y(isc:iec,jsd:jec), Hgrid%N%x(isc:ied,jsd:jec), Hgrid%C%y(isc:ied,jsd:jec)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of N-cell
    call f_plane_cell_length_area(Hgrid%N, Hgrid%N%x(isc:iec,jsc:jec), Hgrid%N%y(isc:iec,jsc:jec),   &
         Hgrid%C%x(isd:iec,jsc:jec), Hgrid%C%y(isd:iec,jsc:jec), Hgrid%T%x(isc:iec,jsc:jed),         &
         Hgrid%T%y(isc:iec,jsc:jed), Hgrid%E%x(isd:iec,jsc:jed), Hgrid%E%y(isd:iec,jsc:jed)  )

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of C-cell
    call f_plane_cell_length_area(Hgrid%C, Hgrid%C%x(isc:iec,jsc:jec), Hgrid%C%y(isc:iec,jsc:jec),   &
         Hgrid%N%x(isc:ied,jsc:jec), Hgrid%N%y(isc:ied,jsc:jec), Hgrid%E%x(isc:iec,jsc:jed),         &
         Hgrid%E%y(isc:iec,jsc:jed), Hgrid%T%x(isc:ied,jsc:jed), Hgrid%T%y(isc:ied,jsc:jed)  )

    !--- for spherical grid, rotation angle always is 0 ------------------
    Hgrid%T%angle = 0.0
    Hgrid%E%angle = 0.0
    Hgrid%N%angle = 0.0
    Hgrid%C%angle = 0.0  

    return
  end subroutine f_plane_grid_init

 !#######################################################################
  !--- generate simple cartesian grid
  subroutine simple_cartesian_grid_init(Hgrid)
    type(hgrid_data_type), intent(inout) :: Hgrid

    !---------------------------------------------------------------------
    !--- calculate face length, half length area and rotation angle of T-cell
    call cartesian_cell_length_area(Hgrid%T)

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of E-cell

    call cartesian_cell_length_area(Hgrid%E)

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of N-cell
    call cartesian_cell_length_area(Hgrid%N)

    !---------------------------------------------------------------------
    !--- calculate face length, half length, area and rotation angle of C-cell
    call cartesian_cell_length_area(Hgrid%C)

    !--- for spherical grid, rotation angle always is 0 ------------------
    Hgrid%T%angle = 0.0
    Hgrid%E%angle = 0.0
    Hgrid%N%angle = 0.0
    Hgrid%C%angle = 0.0  

    return
  end subroutine simple_cartesian_grid_init


  !#######################################################################
  !--- generate tripolar grid
  subroutine tripolar_grid_init(Hgrid)

    type(hgrid_data_type), intent(inout) :: Hgrid

    !--- local variables--------------------------------------------------
    integer :: i, j
    real :: lon_last, lat_join_old
    !---------------------------------------------------------------------
    lon_start = x_lon(1)
    lon_end   = x_lon(2)
    ! recompute join latitude based on nearest spherical latitude
    lat_join_old = lat_join
    lat_join = Hgrid%C%y(1,nearest_index(lat_join, yu0(1:nj)))
    lon_bpeq = lon_start + 90.
    lon_bpnp = lon_start 
    lon_bpsp = lon_start+180.

    if(lat_join_old .ne. lat_join) then
       write(stdout(),*) 'NOTE: Change join latitude from ', lat_join_old, 'to ', lat_join
    endif

    !---------------------------------------------------------------------
    ! Transform from bipolar grid coordinates (bp_lon, bp_lat) to
    ! geographic coordinates (geolon_t, geolat_t) following R. Murray,
    ! "Explicit generation of orthogonal grids for ocean models",
    ! 1996, J.Comp.Phys., v. 126, p. 251-273.  All equation
    ! numbers refer to the Murray paper.
    !---------------------------------------------------------------------
    lam0 = mod(lon_bpeq*D2R + tpi,tpi)
    rp   = tan((hpi-lat_join*D2R)/2.)             ! eqn. 2

    !--- update T-cell length at bi-polar region ---------------------------
    call update_cell_length(Hgrid%T, Hgrid%T%y(1,jsc:jec), Hgrid%T%x(isc:iec,jsc:jec), Hgrid%T%y(isc:iec,jsc:jec), &
         Hgrid%E%x(isd:iec,jsc:jec), Hgrid%E%y(isd:iec,jsc:jec), Hgrid%N%x(isc:iec,jsd:jec),                       &
         Hgrid%N%y(isc:iec,jsd:jec), Hgrid%C%x(isd:iec,jsd:jec), Hgrid%C%y(isd:iec,jsd:jec)  )

    !--- update E-cell length at bi-polar region ---------------------------
    call update_cell_length(Hgrid%E, Hgrid%T%y(1,jsc:jec), Hgrid%E%x(isc:iec,jsc:jec), Hgrid%E%y(isc:iec,jsc:jec), &
         Hgrid%T%x(isc:ied,jsc:jec), Hgrid%T%y(isc:ied,jsc:jec), Hgrid%C%x(isc:iec,jsd:jec),                       &
         Hgrid%C%y(isc:iec,jsd:jec), Hgrid%N%x(isc:ied,jsd:jec), Hgrid%C%y(isc:ied,jsd:jec)  )

    !--- update N-cell length at bi-polar region ---------------------------
    call update_cell_length(Hgrid%N, Hgrid%T%y(1,jsc:jec), Hgrid%N%x(isc:iec,jsc:jec), Hgrid%N%y(isc:iec,jsc:jec), &
         Hgrid%C%x(isd:iec,jsc:jec), Hgrid%C%y(isd:iec,jsc:jec), Hgrid%T%x(isc:iec,jsc:jed),                       &
         Hgrid%T%y(isc:iec,jsc:jed), Hgrid%E%x(isd:iec,jsc:jed), Hgrid%E%y(isd:iec,jsc:jed)  )

    !--- update C-cell length at bi-polar region ---------------------------
    call update_cell_length(Hgrid%C, Hgrid%T%y(1,jsc:jec), Hgrid%C%x(isc:iec,jsc:jec), Hgrid%C%y(isc:iec,jsc:jec), &
         Hgrid%N%x(isc:ied,jsc:jec), Hgrid%N%y(isc:ied,jsc:jec), Hgrid%E%x(isc:iec,jsc:jed),                       &
         Hgrid%E%y(isc:iec,jsc:jed), Hgrid%T%x(isc:ied,jsc:jed), Hgrid%T%y(isc:ied,jsc:jed)  )

    !--- update the length of north to make sure the last grid box is folded.
    call update_north_length(Hgrid)

    ! define corner locations
    do j=0,nj+1
       do i=0,ni+1
          lon_last = Hgrid%C%x(max(i-1,0),j)
          if (Hgrid%T%y(i,min(nj+1,j+1)) >= lat_join) then
             call tp_trans(Hgrid%C%x(i,j),Hgrid%C%y(i,j),lon_last) 
             ! geographical position of bipolar points
          end if
       end do
    end do

    ! define tracer, east and north locations
    do j=0,nj+1
       do i=0,ni+1
          lon_last = Hgrid%T%x(max(i-1,0),j)
          if (Hgrid%T%y(1,j) >= lat_join) then
             call tp_trans(Hgrid%T%x(i,j),Hgrid%T%y(i,j), lon_last) 
          end if
       end do
    end do

    ! define tracer locations
    do j=0,nj+1
       do i=0,ni+1
          lon_last = Hgrid%E%x(max(i-1,0),j)
          if (Hgrid%T%y(1,j) >= lat_join) then
             call tp_trans(Hgrid%E%x(i,j),Hgrid%E%y(i,j), lon_last) 
          end if
       end do
    end do

    ! define north locations
    do j=0,nj+1
       do i=0,ni+1
          lon_last = Hgrid%N%x(max(i-1,0),j)
          if (Hgrid%T%y(1,min(nj+1,j+1)) >= lat_join) then
             call tp_trans(Hgrid%N%x(i,j),Hgrid%N%y(i,j), lon_last) 
          end if
       end do
    end do

    !--- update T-cell area at tripolar region ---------------------------
    call update_cell_area(Hgrid%T,  Hgrid%T%y(1,jsc:jec),                                           &
         Hgrid%T%x(isc:iec,jsc:jec), Hgrid%T%y(isc:iec,jsc:jec), Hgrid%E%x(isd:iec,jsc:jec), Hgrid%E%y(isd:iec,jsc:jec),   &
         Hgrid%N%x(isc:iec,jsd:jec), Hgrid%N%y(isc:iec,jsd:jec), Hgrid%C%x(isd:iec,jsd:jec), Hgrid%C%y(isd:iec,jsd:jec)  )

    !--- update E-cell area at bi-polar region ---------------------------
    call update_cell_area(Hgrid%E,  Hgrid%T%y(1,jsc:jec),                                               &
         Hgrid%E%x(isc:iec,jsc:jec), Hgrid%E%y(isc:iec,jsc:jec), Hgrid%T%x(isc:ied,jsc:jec), Hgrid%T%y(isc:ied,jsc:jec),   &
         Hgrid%C%x(isc:iec,jsd:jec), Hgrid%C%y(isc:iec,jsd:jec), Hgrid%N%x(isc:ied,jsd:jec), Hgrid%N%y(isc:ied,jsd:jec)  )

    !--- update N-cell area at bi-polar region ---------------------------
    call update_cell_area(Hgrid%N,  Hgrid%T%y(1,jsc:jec),                                                  &
         Hgrid%N%x(isc:iec,jsc:jec), Hgrid%N%y(isc:iec,jsc:jec), Hgrid%C%x(isd:iec,jsc:jec), Hgrid%C%y(isd:iec,jsc:jec),    &
         Hgrid%T%x(isc:iec,jsc:jed), Hgrid%T%y(isc:iec,jsc:jed), Hgrid%E%x(isd:iec,jsc:jed), Hgrid%E%y(isd:iec,jsc:jed) )

    !--- update C-cell area at bi-polar region ---------------------------
    call update_cell_area(Hgrid%C,  Hgrid%T%y(1,jsc:jec),                                                    &
         Hgrid%C%x(isc:iec,jsc:jec), Hgrid%C%y(isc:iec,jsc:jec), Hgrid%N%x(isc:ied,jsc:jec), Hgrid%N%y(isc:ied,jsc:jec),    &
         Hgrid%E%x(isc:iec,jsc:jed), Hgrid%E%y(isc:iec,jsc:jed), Hgrid%T%x(isc:ied,jsc:jed), Hgrid%T%y(isc:ied,jsc:jed) )  

    !--- update rotation angle of T-cell -----------------------------
    call update_cell_angle(Hgrid%T%angle, Hgrid%T%y(isc:iec,jsc:jec), &
         Hgrid%E%x(isd:iec,jsc:jec), Hgrid%E%y(isd:iec,jsc:jec) )

    !--- update rotation angle of E-cell -----------------------------
    call update_cell_angle(Hgrid%E%angle, Hgrid%E%y(isc:iec,jsc:jec), &
         Hgrid%T%x(isc:ied,jsc:jec), Hgrid%T%y(isc:ied,jsc:jec) )

    !--- update rotation angle of N-cell -----------------------------
    call update_cell_angle(Hgrid%N%angle, Hgrid%N%y(isc:iec,jsc:jec), &
         Hgrid%C%x(isd:iec,jsc:jec), Hgrid%C%y(isd:iec,jsc:jec) )

    !--- update rotation angle of C-cell -----------------------------
    call update_cell_angle(Hgrid%C%angle, Hgrid%C%y(isc:iec,jsc:jec), &
         Hgrid%N%x(isc:ied,jsc:jec), Hgrid%N%y(isc:ied,jsc:jec) ) 
    return

  end subroutine tripolar_grid_init

  !#######################################################################
  !--- define vertices for each cell
  subroutine define_vertices(Hgrid)

    type(hgrid_data_type), intent(inout) :: Hgrid

    !--- define vertices of T-cell
    call define_cell_vertices(Hgrid%T, Hgrid%C, 0, 0)
    !--- define vertices of E-cell
    call define_cell_vertices(Hgrid%E, Hgrid%N, 1, 0)
    !--- define vertices of N-cell
    call define_cell_vertices(Hgrid%N, Hgrid%E, 0, 1)
    !--- define vertices of C-cell
    call define_cell_vertices(Hgrid%C, Hgrid%T, 1, 1)

    return
  end subroutine define_vertices

  !#######################################################################
  !--- memory allocation
  subroutine allocate_hgrid_cell(cell)

    type(cell_type), intent(inout) :: cell

    allocate(cell%x(0:ni+1,0:nj+1), cell%y(0:ni+1,0:nj+1),  &
         cell%x_vert(isc:iec,jsc:jec,4),  cell%y_vert(isc:iec,jsc:jec,4),   &
         cell%area(isc:iec,jsc:jec),      cell%angle(isc:iec,jsc:jec),      &
         cell%ds_00_02(isc:iec,jsc:jec),  cell%ds_20_22(isc:iec,jsc:jec),   &
         cell%ds_02_22(isc:iec,jsc:jec),  cell%ds_00_20(isc:iec,jsc:jec),   &
         cell%ds_00_01(isc:iec,jsc:jec),  cell%ds_01_02(isc:iec,jsc:jec),   & 
         cell%ds_02_12(isc:iec,jsc:jec),  cell%ds_12_22(isc:iec,jsc:jec),   &
         cell%ds_21_22(isc:iec,jsc:jec),  cell%ds_20_21(isc:iec,jsc:jec),   &
         cell%ds_10_20(isc:iec,jsc:jec),  cell%ds_00_10(isc:iec,jsc:jec),   &
         cell%ds_01_11(isc:iec,jsc:jec),  cell%ds_11_12(isc:iec,jsc:jec),   &
         cell%ds_11_21(isc:iec,jsc:jec),  cell%ds_10_11(isc:iec,jsc:jec),   &
         cell%ds_01_21(isc:iec,jsc:jec),  cell%ds_10_12(isc:iec,jsc:jec) )
    return

  end subroutine allocate_hgrid_cell

  !#######################################################################
  !--- release memory
  subroutine deallocate_hgrid_cell(cell)
    type(cell_type), intent(inout) :: cell

    deallocate(cell%x, cell%y, cell%x_vert, cell%y_vert, cell%area, cell%angle, &
         cell%ds_00_02, cell%ds_20_22, cell%ds_02_22, cell%ds_00_20,      &
         cell%ds_00_01, cell%ds_01_02, cell%ds_02_12, cell%ds_12_22,      &
         cell%ds_21_22, cell%ds_20_21, cell%ds_10_20, cell%ds_00_10,      &
         cell%ds_01_11, cell%ds_11_12, cell%ds_11_21, cell%ds_10_11,      &
         cell%ds_01_21,  cell%ds_10_12 ) 
    return

  end subroutine deallocate_hgrid_cell

  !#######################################################################
  !--- write meta data 
  subroutine write_cell_meta(file, cell_id )

    character(len=*), intent(in) :: file
    character(len=1), intent(in) :: cell_id
    character(len=1)             :: x_pos, y_pos

    select case(cell_id)
    case('T')
       x_pos='T'; y_pos='T'
    case('E')
       x_pos='C'; y_pos='T'
    case('N')
       x_pos='T'; y_pos='C'
    case('C')
       x_pos='C'; y_pos='C'
    end select

    call write_field_meta(file, 'x_'//cell_id, 'degree_east',  &
         'Geographic longitude of '//cell_id//'_cell centers', 2, x_pos, y_pos)
    call write_field_meta(file, 'y_'//cell_id, 'degree_north',  &
         'Geographic latitude of '//cell_id//'_cell centers', 2, x_pos, y_pos)
    call write_field_meta(file,  'x_vert_'//cell_id,'degree_east', 'Geographic longitude of ' &
         //cell_id//'_cell vertices begin southwest counterclockwise', 3, x_pos, y_pos)
    call write_field_meta(file,  'y_vert_'//cell_id,'degree_north', 'Geographic latitude of ' &
         //cell_id //'_cell vertices begin southwest counterclockwise', 3, x_pos, y_pos)
    call write_field_meta(file,  'area_'//cell_id, 'm2','Area of '//cell_id//'_cell ',2,x_pos, y_pos)
    call write_field_meta(file,  'angle_'//cell_id, 'degree','Angle clockwise between logical ' &
         //'and geographic east of '// cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file,  'ds_00_02_'//cell_id, 'm',  &
         'Length of western face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_20_22_'//cell_id, 'm',  &
         'Length of eastern face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_02_22_'//cell_id, 'm',  &
         'Length of northern face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_00_20_'//cell_id, 'm',  &
         'Length of southern face of '//cell_id//'_cell', 2,x_pos, y_pos)  
    call write_field_meta(file, 'ds_00_01_'//cell_id, 'm',  &
         'Distance from southwest corner to western face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_01_02_'//cell_id, 'm',  &
         'Distance from northwest corner to western face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_02_12_'//cell_id, 'm',  &
         'Distance from northwest corner to northern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_12_22_'//cell_id, 'm',  &
         'Distance from northeast corner to northern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_21_22_'//cell_id, 'm',  &
         'Distance from northeast corner to eastern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_20_21_'//cell_id, 'm',  &
         'Distance from southeast corner to eastern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_10_20_'//cell_id, 'm',  &
         'Distance from southeast corner to southern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_00_10_'//cell_id, 'm',  &
         'Distance from southwest corner to southern face center of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_01_11_'//cell_id, 'm',  &
         'Distance from center to western face of '//cell_id//'_cell', 2, x_pos, y_pos)
    call write_field_meta(file, 'ds_11_12_'//cell_id, 'm',  &
         'Distance from center to northern face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_11_21_'//cell_id, 'm',  &
         'Distance from center to eastern face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_10_11_'//cell_id, 'm',  &
         'Distance from center to southern face of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_01_21_'//cell_id, 'm', 'width of '//cell_id//'_cell', 2,x_pos, y_pos)
    call write_field_meta(file, 'ds_10_12_'//cell_id, 'm','height of '//cell_id//'_cell', 2,x_pos, y_pos)
    return

  end subroutine write_cell_meta

  !#######################################################################
  !--- write out data
  subroutine write_cell_data(file, cell, cell_id)

    character(len=*), intent(in) :: file
    type(cell_type),  intent(in) :: cell
    character(len=1), intent(in) :: cell_id
    real, allocatable :: tmp_2d(:,:), tmp_3d(:,:,:)

    allocate(tmp_2d(ni,nj), tmp_3d(ni,nj,4))

    !--- write out center of cell
    tmp_2d = cell%x(1:ni,1:nj)
    call write_field_data(file, 'x_'//cell_id, tmp_2d)
    tmp_2d = cell%y(1:ni,1:nj)
    call write_field_data(file, 'y_'//cell_id, tmp_2d)  

    !--- write out vertices cell 
    call mpp_global_field(Domain, cell%x_vert, tmp_3d)
    call write_field_data(file, 'x_vert_'//cell_id, tmp_3d)
    call mpp_global_field(Domain, cell%y_vert, tmp_3d)
    call write_field_data(file, 'y_vert_'//cell_id, tmp_3d)

    !--- write out area of cell 
    call mpp_global_field(Domain, cell%area, tmp_2d)
    call write_field_data(file, 'area_'//cell_id, tmp_2d)

    !--- write out angle between i-unit and x-unit vector of cell 
    call mpp_global_field(Domain, cell%angle, tmp_2d)
    call write_field_data(file, 'angle_'//cell_id, tmp_2d)

    !--- write out half and face length of cell
    call mpp_global_field(Domain, cell%ds_00_02, tmp_2d) 
    call write_field_data(file, 'ds_00_02_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_20_22, tmp_2d) 
    call write_field_data(file, 'ds_20_22_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_02_22, tmp_2d) 
    call write_field_data(file, 'ds_02_22_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_00_20, tmp_2d) 
    call write_field_data(file, 'ds_00_20_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_00_01, tmp_2d) 
    call write_field_data(file, 'ds_00_01_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_01_02, tmp_2d) 
    call write_field_data(file, 'ds_01_02_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_02_12, tmp_2d) 
    call write_field_data(file, 'ds_02_12_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_12_22, tmp_2d) 
    call write_field_data(file, 'ds_12_22_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_21_22, tmp_2d) 
    call write_field_data(file, 'ds_21_22_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_20_21, tmp_2d) 
    call write_field_data(file, 'ds_20_21_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_10_20, tmp_2d) 
    call write_field_data(file, 'ds_10_20_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_00_10, tmp_2d) 
    call write_field_data(file, 'ds_00_10_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_01_11, tmp_2d) 
    call write_field_data(file, 'ds_01_11_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_11_12, tmp_2d) 
    call write_field_data(file, 'ds_11_12_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_11_21, tmp_2d) 
    call write_field_data(file, 'ds_11_21_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_10_11, tmp_2d) 
    call write_field_data(file, 'ds_10_11_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_01_21, tmp_2d) 
    call write_field_data(file, 'ds_01_21_'//cell_id, tmp_2d)
    call mpp_global_field(Domain, cell%ds_10_12, tmp_2d) 
    call write_field_data(file, 'ds_10_12_'//cell_id, tmp_2d)

    deallocate(tmp_2d, tmp_3d)

    return

  end subroutine write_cell_data

  !#######################################################################
  !--- define cell vertices
  subroutine define_cell_vertices(cell_out, cell_in, ioff, joff )
    type(cell_type), intent(inout) :: cell_out
    type(cell_type),   intent(in)  :: cell_in
    integer,           intent(in)  :: ioff, joff

    cell_out%x_vert(isc:iec,jsc:jec,1) = cell_in%x(isc-1+ioff:iec-1+ioff, jsc-1+joff:jec-1+joff )
    cell_out%x_vert(isc:iec,jsc:jec,2) = cell_in%x(isc+ioff  :iec+ioff,   jsc-1+joff:jec-1+joff )
    cell_out%x_vert(isc:iec,jsc:jec,3) = cell_in%x(isc+ioff  :iec+ioff,   jsc+joff  :jec+joff   )
    cell_out%x_vert(isc:iec,jsc:jec,4) = cell_in%x(isc-1+ioff:iec-1+ioff, jsc+joff  :jec+joff   )
    cell_out%y_vert(isc:iec,jsc:jec,1) = cell_in%y(isc-1+ioff:iec-1+ioff, jsc-1+joff:jec-1+joff )
    cell_out%y_vert(isc:iec,jsc:jec,2) = cell_in%y(isc+ioff  :iec+ioff,   jsc-1+joff:jec-1+joff )
    cell_out%y_vert(isc:iec,jsc:jec,3) = cell_in%y(isc+ioff  :iec+ioff,   jsc+joff  :jec+joff   )
    cell_out%y_vert(isc:iec,jsc:jec,4) = cell_in%y(isc-1+ioff:iec-1+ioff, jsc+joff  :jec+joff   )

    return

  end subroutine define_cell_vertices

  !#######################################################################
  !--- distance (in degrees) between points on lat. circle
  function lat_dist(x1,x2) 
    real :: x1, x2, lat_dist

    lat_dist=min(mod(x1-x2+720,360.),mod(x2-x1+720,360.))

  end function lat_dist

  !#######################################################################
  !--- find bipolar grid longitude given geo. coordinates
  function bp_lam(x,y)
    real :: x, y, bp_lam
    !  bp_lam = ((90-y)/(90-lat_join))*90
    ! invert eqn. 5 with phic=0 to place point at specified geo. lat
    bp_lam = 2.*atan(tan((hpi-y*D2R)/2)/rp)/D2R
    if (lat_dist(x,lon_bpeq)<90.) bp_lam = -bp_lam
  end function bp_lam

  !#######################################################################
  !--- find bipolar grid latitude given geo. coordinates
  function bp_phi(x,y) ! find bipolar grid latitude given geo. coordinates
    real :: x, y, bp_phi

    if (lat_dist(x,lon_bpsp)<90.) then
       bp_phi = (-90+lat_dist(x,lon_bpsp))
    else
       bp_phi = ( 90-lat_dist(x,lon_bpnp))
    end if
  end function bp_phi

  !#######################################################################
  !--- calculate tripolar grid
  subroutine tp_trans(lon,lat,lon_ref)

    real, intent(inout) :: lon, lat
    real,    intent(in) :: lon_ref

    real :: lamc, phic, lams, chic, phis

    lamc = bp_lam(lon,lat)*D2R
    phic = bp_phi(lon,lat)*D2R


    if (abs(lat-90.) < 1.e-4) then
       if (phic > 0) then
          lon=lon_in_range(lon_start,lon_ref)
       else
          lon=lon_start+180.
       endif
       chic = acos(cos(lamc)*cos(phic))                          ! eqn. 6
       phis = pi/2-2*atan(rp*tan(chic/2))                        ! eqn. 5
       lat = phis/D2R
       return
    endif

    if (abs(lamc) < 1.e-4 .and. abs(phic) < 1.e-4) then
       lat=90.;lon=lon_ref
    else
       lams = mod(lam0+pi+pi/2-atan2(sin(lamc),tan(phic)),2*pi)  ! eqn. 5
       chic = acos(cos(lamc)*cos(phic))                          ! eqn. 6
       phis = pi/2-2*atan(rp*tan(chic/2))                        ! eqn. 5
       lon = lams/D2R
       lon = lon_in_range(lon,lon_ref) 
       lat = phis/D2R
    endif

    return

  end subroutine tp_trans

  !#######################################################################
  !--- distant on the earth
  function dist(a, b, met1, met2)

    real, intent(in) :: a, b, met1, met2
    real             :: dist

    dist = abs(a-b)/radian*(met1+met2)/2.

  end function dist

  !#######################################################################
  !--- distance between spherical grid on the earth
  function spherical_dist(x1,y1,x2,y2)
    real, intent(in) :: x1, y1, x2, y2
    real             :: spherical_dist, h1, h2

    spherical_dist = 0.0
    if(x1 == x2) then
       h1 = radius
       h2 = radius
       spherical_dist = dist(y1,y2,h1,h2)
    else if(y1 == y2) then
       h1 = radius * cos(y1/RADIAN)
       h2 = radius * cos(y2/RADIAN) 
       spherical_dist = dist(x1,x2,h1,h2)
    else 
       call mpp_error(FATAL,'hgrid_mod: This is not rectangular grid')
    endif
    return

  end function spherical_dist

 !#######################################################################
  !--- distance between cartesian grid on the earth
  function cartesian_dist(x1,y1,x2,y2)
    real, intent(in) :: x1, y1, x2, y2
    real             :: cartesian_dist

    cartesian_dist = 0.0
    if(x1 == x2) then
       cartesian_dist = abs(y2-y1)/radian*radius
    else if(y1 == y2) then
       cartesian_dist = abs(x2-x1)/radian*radius*cos(f_plane_latitude*D2R)
    else 
       call mpp_error(FATAL,'hgrid_mod( from cartesian_dist): This is not rectangular grid')
    endif
    return

  end function cartesian_dist

  !#######################################################################
  !--- distance of bipolar grids
  function bipolar_dist(x1,y1,x2,y2,num_div)
    real,     intent(in) :: x1, y1, x2, y2
    integer, intent(out) :: num_div
    real                 :: bipolar_dist

    real, allocatable, dimension(:) :: x, y, h1, h2 , bp_lon, bp_lat, metric 
    real :: chic
    integer :: n
    !--- right now suppose each line is divied into one segment ----------
    !--- higher order approximation will be implemented later ------------
    num_div = 1

    allocate(x(num_div+1),y(num_div+1), h1(num_div+1), h2(num_div+1), &
         bp_lon(num_div+1), bp_lat(num_div+1), metric(num_div+1) )
    do n = 1, num_div+1
       x(n) = x1 + real(n-1)*(x2-x1)/real(num_div)
       y(n) = y1 + real(n-1)*(y2-y1)/real(num_div)
    enddo

    !--- get the bipolar grid and metric term ----------------------------
    do n =1, num_div+1
       bp_lon(n) = bp_lam(x(n),y(n)) ! longitude (degrees) in bipolar grid system
       bp_lat(n) = bp_phi(x(n),y(n)) ! latitude (degrees) in bipolar grid system
       h1(n) = radius*cos(bp_lat(n)*D2R)
       h2(n) = radius
       metric(n) = 1.0
       if (abs(y(n)-90.0) < 1.e-4 .or. abs(bp_lon(n)*D2R) .ge. 1.e-4 .or. abs(bp_lat(n)*D2R) .ge. 1.e-4) then
          chic = acos(cos(bp_lon(n)*D2R)*cos(bp_lat(n)*D2R))  ! eqn. 6
          metric(n) = rp*(1/cos(chic/2)**2)/(1+(rp**2)*(tan(chic/2)**2)) ! eq 3
       endif
    enddo

    !--- then calculate the distance -------------------------------------
    bipolar_dist = 0.0
    do n = 1, num_div
       if(x1 == x2) then
          bipolar_dist = bipolar_dist+dist(bp_lon(n),bp_lon(n+1),metric(n)*h1(n),metric(n+1)*h1(n+1))
       else if(y1 == y2) then
          bipolar_dist = bipolar_dist+dist(bp_lat(n),bp_lat(n+1),metric(n)*h2(n),metric(n+1)*h2(n+1))
       endif
    enddo

    deallocate(x, y, h1, h2, bp_lon, bp_lat, metric)

    return

  end function bipolar_dist

  !#######################################################################
  !--- rectangular grid box area
  function spherical_box_area(x1,y1,x2,y2,x3,y3,x4,y4)

    real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
    real :: spherical_box_area

    real, dimension(4) :: x_vert,y_vert
    integer :: i, ip, n
    real    :: lat1, lat2, dx

    x_vert(1) = x1; y_vert(1) = y1
    x_vert(2) = x2; y_vert(2) = y2
    x_vert(3) = x3; y_vert(3) = y3
    x_vert(4) = x4; y_vert(4) = y4

    spherical_box_area = 0.0
    n = size(x_vert)
    do i = 1, n
       ip = i+1
       if(ip .gt. n) ip = ip - n
       dx = x_vert(ip) - x_vert(i)
       if(cyclic_x) then
          if(dx > x_period/2) dx = dx - x_period
          if(dx < -x_period/2) dx = dx + x_period
       end if
       dx = dx*D2R
       if(dx==0.0) cycle
       lat1 = y_vert(ip)*D2R
       lat2 = y_vert(i)*D2R

       if (lat1 == lat2) then ! cheap area calculation along latitude 
          spherical_box_area = spherical_box_area - dx*sin(lat1)
       else
          spherical_box_area = spherical_box_area - dx*(sin(lat1)+sin(lat2))/2   !  TRAPEZOID_RULE
       endif
    enddo

    spherical_box_area = spherical_box_area * radius * radius

    return

  end function spherical_box_area

  !#######################################################################
  !--- rectangular grid box area for cartesian grid
  function cartesian_box_area(x1,y1,x2,y2)

    real, intent(in) :: x1, y1, x2, y2
    real :: cartesian_box_area

    real    :: dx, dy

    dx = x2-x1
    if(cyclic_x) then
       if(dx > x_period/2)  dx = dx - x_period
       if(dx < -x_period/2) dx = dx + x_period
    end if
    dx = dx * D2R

    dy = y2-y1
    if(cyclic_y) then
       if(dy > y_period/2)  dy = dy - y_period
       if(dy < -y_period/2) dy = dy + y_period
    end if
    dy = dy * D2R
    cartesian_box_area = dx*dy*cos(f_plane_latitude*D2R)*radius*radius

    return

  end function cartesian_box_area

  !#######################################################################
  !--- bipolar grid area
  function bipolar_area(x1,y1,x2,y2,x3,y3,x4,y4 )
    real, intent(in)    :: x1,y1,x2,y2,x3,y3,x4,y4

    real :: bipolar_area

    real, dimension(8) :: x,y
    integer :: i, ip, n
    real    :: lat1, lat2, dx

    x(1) = x1; y(1) = y1
    x(2) = x2; y(2) = y2
    x(3) = x3; y(3) = y3
    x(4) = x4; y(4) = y4

    !--- first fix the longitude at the pole -----------------------------
    call lon_fix(x, y, 4, n, 180.)

    !--- calculate the area ----------------------------------------------  
    bipolar_area = 0.0

    do i = 1, n
       ip = i+1
       if(ip .gt. n) ip = 1
       dx = (x(ip) - x(i))*D2R
       lat1 = y(ip)*D2R
       lat2 = y(i)*D2R
       if(dx==0.0) cycle
       if(dx > pi)  dx = dx - tpi
       if(dx < -pi) dx = dx + tpi

       if (lat1 == lat2) then ! cheap area calculation along latitude 
          bipolar_area = bipolar_area - dx*sin(lat1)
       else
          bipolar_area = bipolar_area - dx*(sin(lat1)+sin(lat2))/2   !  TRAPEZOID_RULE
       endif
    enddo

    bipolar_area = bipolar_area * radius * radius

  end function bipolar_area

  !#######################################################################
  !calculate face length and half length of each cell for spherical grid.
  subroutine spherical_cell_length_area(cell, x, y, x_e, y_e, x_n, y_n, x_c, y_c)
    type(cell_type),                intent(inout) :: cell
    real, dimension(isc:iec,jsc:jec), intent(in)  :: x, y     ! center of the cell.
    real, dimension(isd:iec,jsc:jec), intent(in)  :: x_e, y_e ! east of the cell.
    real, dimension(isc:iec,jsd:jec), intent(in)  :: x_n, y_n ! north of the cell.
    real, dimension(isd:iec,jsd:jec), intent(in)  :: x_c, y_c ! four vertices of cell

    integer :: i, j
    real :: area1, area2, area3, area4

    do j = jsc,jec
       do i = isc,iec
          !--- define half length ------------------------------------------- 
          cell%ds_00_01(i,j) = spherical_dist(x_e(i-1,j), y_e(i-1,j), x_c(i-1,j-1), y_c(i-1,j-1) )
          cell%ds_01_02(i,j) = spherical_dist(x_c(i-1,j), y_c(i-1,j), x_e(i-1,j),   y_e(i-1,j)   )
          cell%ds_02_12(i,j) = spherical_dist(x_n(i,j),   y_n(i,j),   x_c(i-1,j),   y_c(i-1,j)   )
          cell%ds_12_22(i,j) = spherical_dist(x_c(i,j),   y_c(i,j),   x_n(i,j),     y_n(i,j)     )
          cell%ds_21_22(i,j) = spherical_dist(x_c(i,j),   y_c(i,j),   x_e(i,j),     y_e(i,j)     )
          cell%ds_20_21(i,j) = spherical_dist(x_e(i,j),   y_e(i,j),   x_c(i,j-1),   y_c(i,j-1)   )
          cell%ds_10_20(i,j) = spherical_dist(x_c(i,j-1), y_c(i,j-1), x_n(i,j-1),   y_n(i,j-1)   )
          cell%ds_00_10(i,j) = spherical_dist(x_n(i,j-1), y_n(i,j-1), x_c(i-1,j-1), y_c(i-1,j-1) )
          cell%ds_01_11(i,j) = spherical_dist(x(i,j),     y(i,j),     x_e(i-1,j),   y_e(i-1,j)   )
          cell%ds_11_12(i,j) = spherical_dist(x_n(i,j),   y_n(i,j),   x(i,j),       y(i,j)       )
          cell%ds_11_21(i,j) = spherical_dist(x_e(i,j),   y_e(i,j),   x(i,j),       y(i,j)       )
          cell%ds_10_11(i,j) = spherical_dist(x(i,j),     y(i,j),     x_n(i,j-1),   y_n(i,j-1)   )  
          !--- define full length -------------------------------------------
          cell%ds_00_02(i,j) = cell%ds_00_01(i,j) + cell%ds_01_02(i,j)
          cell%ds_20_22(i,j) = cell%ds_20_21(i,j) + cell%ds_21_22(i,j)
          cell%ds_00_20(i,j) = cell%ds_00_10(i,j) + cell%ds_10_20(i,j)
          cell%ds_02_22(i,j) = cell%ds_02_12(i,j) + cell%ds_12_22(i,j)
          cell%ds_01_21(i,j) = cell%ds_01_11(i,j) + cell%ds_11_21(i,j)
          cell%ds_10_12(i,j) = cell%ds_10_11(i,j) + cell%ds_11_12(i,j)
       enddo
    enddo

    !--- we divided the cell into four parts and the area of the cell ----
    !--- is the sum of areas of the four parts ---------------------------
    do j=jsc,jec
       do i=isc,iec
          area1 = spherical_box_area(x_c(i-1,j-1), y_c(i-1,j-1), x_n(i,j-1), y_n(i,j-1), &
               x(i,j), y(i,j), x_e(i-1,j), y_e(i-1,j)  )
          area2 = spherical_box_area(x_n(i,j-1), y_n(i,j-1), x_c(i,j-1), y_c(i,j-1), &
               x_e(i,j), y_e(i,j), x(i,j), y(i,j)  )
          area3 = spherical_box_area(x(i,j), y(i,j), x_e(i,j), y_e(i,j), &
               x_c(i,j), y_c(i,j), x_n(i,j), y_n(i,j)  )
          area4 = spherical_box_area(x_e(i-1,j), y_e(i-1,j), x(i,j), y(i,j), &
               x_n(i,j), y_n(i,j), x_c(i-1,j), y_c(i-1,j) )
          cell%area(i,j) = area1+area2+area3+area4
       enddo
    enddo

    return

  end subroutine spherical_cell_length_area

  !#######################################################################
  !calculate face length and half length of each cell for f_plane/beta_plane grid.
  subroutine f_plane_cell_length_area(cell, x, y, x_e, y_e, x_n, y_n, x_c, y_c)
    type(cell_type),                intent(inout) :: cell
    real, dimension(isc:iec,jsc:jec), intent(in)  :: x, y     ! center of the cell.
    real, dimension(isd:iec,jsc:jec), intent(in)  :: x_e, y_e ! east of the cell.
    real, dimension(isc:iec,jsd:jec), intent(in)  :: x_n, y_n ! north of the cell.
    real, dimension(isd:iec,jsd:jec), intent(in)  :: x_c, y_c ! four vertices of cell

    integer :: i, j

    do j = jsc,jec
       do i = isc,iec
          !--- define half length ------------------------------------------- 
          cell%ds_00_01(i,j) = cartesian_dist(x_e(i-1,j), y_e(i-1,j), x_c(i-1,j-1), y_c(i-1,j-1) )
          cell%ds_01_02(i,j) = cartesian_dist(x_c(i-1,j), y_c(i-1,j), x_e(i-1,j),   y_e(i-1,j)   )
          cell%ds_02_12(i,j) = cartesian_dist(x_n(i,j),   y_n(i,j),   x_c(i-1,j),   y_c(i-1,j)   )
          cell%ds_12_22(i,j) = cartesian_dist(x_c(i,j),   y_c(i,j),   x_n(i,j),     y_n(i,j)     )
          cell%ds_21_22(i,j) = cartesian_dist(x_c(i,j),   y_c(i,j),   x_e(i,j),     y_e(i,j)     )
          cell%ds_20_21(i,j) = cartesian_dist(x_e(i,j),   y_e(i,j),   x_c(i,j-1),   y_c(i,j-1)   )
          cell%ds_10_20(i,j) = cartesian_dist(x_c(i,j-1), y_c(i,j-1), x_n(i,j-1),   y_n(i,j-1)   )
          cell%ds_00_10(i,j) = cartesian_dist(x_n(i,j-1), y_n(i,j-1), x_c(i-1,j-1), y_c(i-1,j-1) )
          cell%ds_01_11(i,j) = cartesian_dist(x(i,j),     y(i,j),     x_e(i-1,j),   y_e(i-1,j)   )
          cell%ds_11_12(i,j) = cartesian_dist(x_n(i,j),   y_n(i,j),   x(i,j),       y(i,j)       )
          cell%ds_11_21(i,j) = cartesian_dist(x_e(i,j),   y_e(i,j),   x(i,j),       y(i,j)       )
          cell%ds_10_11(i,j) = cartesian_dist(x(i,j),     y(i,j),     x_n(i,j-1),   y_n(i,j-1)   )  
          !--- define full length -------------------------------------------
          cell%ds_00_02(i,j) = cell%ds_00_01(i,j) + cell%ds_01_02(i,j)
          cell%ds_20_22(i,j) = cell%ds_20_21(i,j) + cell%ds_21_22(i,j)
          cell%ds_00_20(i,j) = cell%ds_00_10(i,j) + cell%ds_10_20(i,j)
          cell%ds_02_22(i,j) = cell%ds_02_12(i,j) + cell%ds_12_22(i,j)
          cell%ds_01_21(i,j) = cell%ds_01_11(i,j) + cell%ds_11_21(i,j)
          cell%ds_10_12(i,j) = cell%ds_10_11(i,j) + cell%ds_11_12(i,j)
       enddo
    enddo

    !--- we divided the cell into four parts and the area of the cell ----
    !--- is the sum of areas of the four parts ---------------------------
    do j=jsc,jec
       do i=isc,iec
          cell%area(i,j) = cartesian_box_area(x_c(i-1,j-1), y_c(i-1,j-1), x_c(i,j), y_c(i,j))
       enddo
    enddo

    return

  end subroutine f_plane_cell_length_area

  !#######################################################################
  !calculate face length and half length of each cell for cartesian grid.
  subroutine cartesian_cell_length_area(cell)
    type(cell_type),                intent(inout) :: cell

    integer :: i, j

    do j = jsc,jec
       do i = isc,iec
          !--- define half length ------------------------------------------- 
          cell%ds_00_01(i,j) = 0.5*simple_cartesian_dy
          cell%ds_01_02(i,j) = 0.5*simple_cartesian_dy
          cell%ds_02_12(i,j) = 0.5*simple_cartesian_dx
          cell%ds_12_22(i,j) = 0.5*simple_cartesian_dx
          cell%ds_21_22(i,j) = 0.5*simple_cartesian_dy
          cell%ds_20_21(i,j) = 0.5*simple_cartesian_dy
          cell%ds_10_20(i,j) = 0.5*simple_cartesian_dx
          cell%ds_00_10(i,j) = 0.5*simple_cartesian_dx
          cell%ds_01_11(i,j) = 0.5*simple_cartesian_dx
          cell%ds_11_12(i,j) = 0.5*simple_cartesian_dy
          cell%ds_11_21(i,j) = 0.5*simple_cartesian_dx
          cell%ds_10_11(i,j) = 0.5*simple_cartesian_dy
          !--- define full length -------------------------------------------
          cell%ds_00_02(i,j) = simple_cartesian_dy
          cell%ds_20_22(i,j) = simple_cartesian_dy
          cell%ds_00_20(i,j) = simple_cartesian_dx
          cell%ds_02_22(i,j) = simple_cartesian_dx
          cell%ds_01_21(i,j) = simple_cartesian_dx
          cell%ds_10_12(i,j) = simple_cartesian_dy
       enddo
    enddo

    !--- we divided the cell into four parts and the area of the cell ----
    !--- is the sum of areas of the four parts ---------------------------
    do j=jsc,jec
       do i=isc,iec
          cell%area(i,j) = simple_cartesian_dx * simple_cartesian_dy
       enddo
    enddo

    return

  end subroutine cartesian_cell_length_area

  !#######################################################################
  !--- calculate tripolar grid distance
  subroutine update_cell_length(cell, lat, x, y, x_e, y_e, x_n, y_n, x_c, y_c)
    type(cell_type),        intent(inout) :: cell
    real, dimension(jsc:jec),         intent(in)  :: lat
    real, dimension(isc:iec,jsc:jec), intent(in)  :: x, y     ! center of the cell.
    real, dimension(isd:iec,jsc:jec), intent(in)  :: x_e, y_e ! east of the cell.
    real, dimension(isc:iec,jsd:jec), intent(in)  :: x_n, y_n ! north of the cell.
    real, dimension(isd:iec,jsd:jec), intent(in)  :: x_c, y_c ! four vertices of cell

    integer :: i, j
    integer, allocatable :: num_div(:,:,:)

    allocate(num_div(ni,nj,12))

    !---update length-----------------------------------------------------
    do j = jsc,jec
       do i = isc,iec
          if(lat(j) .ge. lat_join) then
             !--- define half length ----------------------------------------
             cell%ds_00_01(i,j) = bipolar_dist(x_e(i-1,j),y_e(i-1,j),x_c(i-1,j-1),y_c(i-1,j-1), num_div(i,j,8) )
             cell%ds_01_02(i,j) = bipolar_dist(x_c(i-1,j), y_c(i-1,j), x_e(i-1,j), y_e(i-1,j), num_div(i,j,7))
             cell%ds_02_12(i,j) = bipolar_dist(x_n(i,j), y_n(i,j), x_c(i-1,j), y_c(i-1,j), num_div(i,j,6))
             cell%ds_12_22(i,j) = bipolar_dist(x_c(i,j), y_c(i,j), x_n(i,j), y_n(i,j), num_div(i,j,5))
             cell%ds_21_22(i,j) = bipolar_dist(x_c(i,j), y_c(i,j), x_e(i,j), y_e(i,j), num_div(i,j,4))
             cell%ds_20_21(i,j) = bipolar_dist(x_e(i,j), y_e(i,j), x_c(i,j-1), y_c(i,j-1), num_div(i,j,3))
             cell%ds_10_20(i,j) = bipolar_dist(x_c(i,j-1), y_c(i,j-1), x_n(i,j-1), y_n(i,j-1), num_div(i,j,2))
             cell%ds_00_10(i,j) = bipolar_dist(x_n(i,j-1), y_n(i,j-1), x_c(i-1,j-1), y_c(i-1,j-1), num_div(i,j,1))
             cell%ds_01_11(i,j) = bipolar_dist(x(i,j), y(i,j), x_e(i-1,j), y_e(i-1,j), num_div(i,j,12))
             cell%ds_11_12(i,j) = bipolar_dist(x_n(i,j), y_n(i,j), x(i,j), y(i,j), num_div(i,j,11) )
             cell%ds_11_21(i,j) = bipolar_dist(x(i,j), y(i,j),x_e(i,j), y_e(i,j), num_div(i,j,10) )
             cell%ds_10_11(i,j) = bipolar_dist(x(i,j), y(i,j), x_n(i,j-1), y_n(i,j-1) , num_div(i,j,9)) 
             !--- define full length ----------------------------------------
             cell%ds_00_02(i,j) = cell%ds_00_01(i,j) + cell%ds_01_02(i,j)
             cell%ds_20_22(i,j) = cell%ds_20_21(i,j) + cell%ds_21_22(i,j)
             cell%ds_02_22(i,j) = cell%ds_02_12(i,j) + cell%ds_12_22(i,j)
             cell%ds_00_20(i,j) = cell%ds_00_10(i,j) + cell%ds_10_20(i,j)
             cell%ds_01_21(i,j) = cell%ds_01_11(i,j) + cell%ds_11_21(i,j)
             cell%ds_10_12(i,j) = cell%ds_10_11(i,j) + cell%ds_11_12(i,j)
          endif
       enddo
    enddo

    deallocate(num_div)

    return

  end subroutine update_cell_length

  !#######################################################################
  !--- calculate tripolar grid area
  subroutine update_cell_area(cell, lat, x, y, x_e, y_e, x_n, y_n, x_c, y_c)
    type(cell_type),                intent(inout) :: cell
    real, dimension(jsc:jec),         intent(in)  :: lat
    real, dimension(isc:iec,jsc:jec), intent(in)  :: x, y     ! center of the cell.
    real, dimension(isd:iec,jsc:jec), intent(in)  :: x_e, y_e ! east of the cell.
    real, dimension(isc:iec,jsd:jec), intent(in)  :: x_n, y_n ! north of the cell.
    real, dimension(isd:iec,jsd:jec), intent(in)  :: x_c, y_c ! four vertices of cell

    integer :: i, j
    real :: area1, area2, area3, area4

    !--- update area -----------------------------------------------------
    do j = jsc,jec
       do i = isc,iec
          if(lat(j) .ge. lat_join) then
             area1 = bipolar_area(x_c(i-1,j-1), y_c(i-1,j-1), x_n(i,j-1), y_n(i,j-1), &
                  x(i,j), y(i,j), x_e(i-1,j), y_e(i-1,j)   )
             area2 = bipolar_area(x_n(i,j-1), y_n(i,j-1), x_c(i,j-1), y_c(i,j-1), &
                  x_e(i,j), y_e(i,j), x(i,j), y(i,j) )
             area3 = bipolar_area(x(i,j), y(i,j), x_e(i,j), y_e(i,j),          &
                  x_c(i,j), y_c(i,j), x_n(i,j), y_n(i,j)  )
             area4 = bipolar_area(x_e(i-1,j), y_e(i-1,j), x(i,j), y(i,j),         &
                  x_n(i,j), y_n(i,j), x_c(i-1,j), y_c(i-1,j) )
             if(area1 .lt. 0 .or. area2 .lt. 0 .or. area3 .lt. 0 .or. area4 .lt. 0) then
                area1 = abs(area1); area2 = abs(area2) 
                area3 = abs(area3); area4 = abs(area4)
             endif
             cell%area(i,j) = area1+area2+area3+area4
          endif
       enddo
    enddo

  end subroutine update_cell_area

  !#######################################################################
  !--- update length at the folded region, for this reason, 
  !--- the domain layout always is (1,npes). in this case
  !--- isc = 1, iec = ni
  subroutine update_north_length(Hgrid)
    type(hgrid_data_type), intent(inout) :: Hgrid

    integer :: i

    if(jec == nj) then

       do i = 1, ni/2
          !--- update most northern T-cell length ---------------------------
          Hgrid%T%ds_02_22(i,nj) = Hgrid%T%ds_02_22(ni-i+1,nj)
          Hgrid%T%ds_02_12(i,nj) = Hgrid%T%ds_12_22(ni-i+1,nj)
          Hgrid%T%ds_12_22(i,nj) = Hgrid%T%ds_02_12(ni-i+1,nj)
          !--- update most southern E-cell length ---------------------------
          Hgrid%E%ds_02_22(i,nj) = Hgrid%E%ds_02_22(ni-i,nj)
          Hgrid%E%ds_02_12(i,nj) = Hgrid%E%ds_12_22(ni-i,nj)
          Hgrid%E%ds_12_22(i,nj) = Hgrid%E%ds_02_12(ni-i,nj)
          !--- update most southern N-cell length ---------------------------
          Hgrid%N%ds_02_22(i,nj) = Hgrid%N%ds_00_20(ni-i+1,nj)
          Hgrid%N%ds_00_20(i,nj) = Hgrid%N%ds_02_22(ni-i+1,nj)
          Hgrid%N%ds_00_02(i,nj) = Hgrid%N%ds_20_22(ni-i+1,nj)
          Hgrid%N%ds_20_22(i,nj) = Hgrid%N%ds_00_02(ni-i+1,nj)
          Hgrid%N%ds_02_12(i,nj) = Hgrid%N%ds_10_20(ni-i+1,nj)
          Hgrid%N%ds_12_22(i,nj) = Hgrid%N%ds_00_10(ni-i+1,nj)
          Hgrid%N%ds_00_10(i,nj) = Hgrid%N%ds_12_22(ni-i+1,nj)
          Hgrid%N%ds_10_20(i,nj) = Hgrid%N%ds_02_12(ni-i+1,nj)
          Hgrid%N%ds_00_01(i,nj) = Hgrid%N%ds_21_22(ni-i+1,nj)
          Hgrid%N%ds_01_02(i,nj) = Hgrid%N%ds_20_21(ni-i+1,nj)
          Hgrid%N%ds_20_21(i,nj) = Hgrid%N%ds_01_02(ni-i+1,nj)
          Hgrid%N%ds_21_22(i,nj) = Hgrid%N%ds_00_01(ni-i+1,nj)
          Hgrid%N%ds_01_11(i,nj) = Hgrid%N%ds_11_21(ni-i+1,nj)
          Hgrid%N%ds_11_12(i,nj) = Hgrid%N%ds_10_11(ni-i+1,nj)
          Hgrid%N%ds_11_21(i,nj) = Hgrid%N%ds_01_11(ni-i+1,nj)
          Hgrid%N%ds_10_11(i,nj) = Hgrid%N%ds_11_12(ni-i+1,nj)
          Hgrid%N%ds_01_21(i,nj) = Hgrid%N%ds_01_21(ni-i+1,nj)
          Hgrid%N%ds_10_12(i,nj) = Hgrid%N%ds_10_12(ni-i+1,nj)
          !--- update most southern C-cell length ---------------------------
          Hgrid%C%ds_02_22(i,nj) = Hgrid%C%ds_00_20(ni-i,nj)
          Hgrid%C%ds_00_20(i,nj) = Hgrid%C%ds_02_22(ni-i,nj)
          Hgrid%C%ds_00_02(i,nj) = Hgrid%C%ds_20_22(ni-i,nj)
          Hgrid%C%ds_20_22(i,nj) = Hgrid%C%ds_00_02(ni-i,nj)
          Hgrid%C%ds_02_12(i,nj) = Hgrid%C%ds_10_20(ni-i,nj)
          Hgrid%C%ds_12_22(i,nj) = Hgrid%C%ds_00_10(ni-i,nj)
          Hgrid%C%ds_00_10(i,nj) = Hgrid%C%ds_12_22(ni-i,nj)
          Hgrid%C%ds_10_20(i,nj) = Hgrid%C%ds_02_12(ni-i,nj)
          Hgrid%C%ds_00_01(i,nj) = Hgrid%C%ds_21_22(ni-i,nj)
          Hgrid%C%ds_01_02(i,nj) = Hgrid%C%ds_20_21(ni-i,nj)
          Hgrid%C%ds_20_21(i,nj) = Hgrid%C%ds_01_02(ni-i,nj)
          Hgrid%C%ds_21_22(i,nj) = Hgrid%C%ds_00_01(ni-i,nj)
          Hgrid%C%ds_01_11(i,nj) = Hgrid%C%ds_11_21(ni-i,nj)
          Hgrid%C%ds_11_12(i,nj) = Hgrid%C%ds_10_11(ni-i,nj)
          Hgrid%C%ds_11_21(i,nj) = Hgrid%C%ds_01_11(ni-i,nj)
          Hgrid%C%ds_10_11(i,nj) = Hgrid%C%ds_11_12(ni-i,nj)
          Hgrid%C%ds_01_21(i,nj) = Hgrid%C%ds_01_21(ni-i,nj)
          Hgrid%C%ds_10_12(i,nj) = Hgrid%C%ds_10_12(ni-i,nj)
       enddo

    endif

    return

  end subroutine update_north_length

  !#######################################################################
  !calculate the rotation angle of the cell
  subroutine update_cell_angle(angle, lat, lone, late)
    real, dimension(isc:iec,jsc:jec), intent(in)    :: lat
    real, dimension(isc:ied,jsc:jec), intent(in)    :: lone, late
    real, dimension(isc:iec,jsc:jec), intent(inout) :: angle
    integer :: i, j
    real    :: lon_scale

    do j = jsc,jec
       do i = isc,iec
          lon_scale = cos(lat(i,j)*D2R)
          angle(i,j) = (atan2(late(i+1,j) - late(i,j), (lone(i+1,j) - lone(i,j))*lon_scale))/D2R
       enddo
    enddo

    return
  end subroutine update_cell_angle

  !#######################################################################
  !--- write out chksum
  subroutine cell_chksum(cell, id)
    type(cell_type), intent(in) :: cell
    character(len=1), intent(in) :: id

    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_x_'//id//' = ', mpp_chksum(cell%x(isc:iec,jsc:jec))
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_y_'//id//' = ', mpp_chksum(cell%y(isc:iec,jsc:jec))
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_x_vert_'//id//'= ', mpp_chksum(cell%x_vert)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_y_vert_'//id//'= ', mpp_chksum(cell%y_vert)  
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_area_'//id//' = ',  mpp_chksum(cell%area)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_angle_'//id//' = ', mpp_chksum(cell%angle)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_00_02_'//id//' = ', mpp_chksum(cell%ds_00_02)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_20_22_'//id//' = ', mpp_chksum(cell%ds_20_22)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_02_22_'//id//' = ', mpp_chksum(cell%ds_02_22)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_00_20_'//id//' = ', mpp_chksum(cell%ds_00_20)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_01_21_'//id//' = ', mpp_chksum(cell%ds_01_21)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_10_12_'//id//' = ', mpp_chksum(cell%ds_10_12)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_00_01_'//id//' = ', mpp_chksum(cell%ds_00_01)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_01_02_'//id//' = ', mpp_chksum(cell%ds_01_02)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_02_12_'//id//' = ', mpp_chksum(cell%ds_02_12)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_12_22_'//id//' = ', mpp_chksum(cell%ds_12_22)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_21_22_'//id//' = ', mpp_chksum(cell%ds_21_22)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_20_21_'//id//' = ', mpp_chksum(cell%ds_20_21)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_10_20_'//id//' = ', mpp_chksum(cell%ds_10_20)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_00_10_'//id//' = ', mpp_chksum(cell%ds_00_10)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_01_11_'//id//' = ', mpp_chksum(cell%ds_01_11)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_11_12_'//id//' = ', mpp_chksum(cell%ds_11_12)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_11_21_'//id//' = ', mpp_chksum(cell%ds_11_21)
    write (stdout(),'(/(a,I20/))')'hgrid chksum: chksum_ds_10_11_'//id//' = ', mpp_chksum(cell%ds_10_11)

  end subroutine cell_chksum

  !#######################################################################
  !--- fix longitude at pole.
  subroutine lon_fix(x,y,n_in,n_out,tlon)
    real, dimension(:), intent(inout) :: x, y
    integer,               intent(in) :: n_in
    integer,              intent(out) :: n_out
    real,                  intent(in) :: tlon

    integer :: i, ip, im
    real :: dx, x_sum

    n_out = n_in
    i = 1
    do
       if(i .gt. n_out) exit

       if(abs(y(i)) .eq. 90.) then
          im= i - 1
          if(im .le. 0) im = im + n_out
          ip = i + 1
          if(ip .gt. n_out) ip = ip - n_out
          !--- all pole points must be paired ----------------------------
          if(y(im) .eq. y(i) .and. y(ip) .eq. y(i) ) then
             call vtx_delete(x,y, n_out, i)
             i = i - 1
          else if(y(im) .ne. y(i) .and. y(ip) .ne. y(i) ) then
             call vtx_insert(x,y,n_out,i)
             i = i + 1
          endif
       endif
       i = i + 1
    enddo

    !--- first of pole pair has longitude of previous vertex -------------
    !--- second of pole pair has longitude of subsequent vertex ----------
    do i = 1, n_out
       if(abs(y(i)) .eq. 90.) then
          im= i - 1
          if(im .le. 0) im = im + n_out
          ip = i + 1
          if(ip .gt. n_out) ip = ip - n_out

          if(y(im) .ne. y(i)) x(i) = x(im)
          if(y(ip) .ne. y(i)) x(i) = x(ip)
       endif
    enddo

    if(n_out == 0) return

    x_sum = x(1)
    do i = 2, n_out
       dx = x(i) - x(i-1)
       if(dx < -180) then
          dx = dx + 360
       else if (dx >  180) then
          dx = dx - 360
       endif
       x(i) = x(i-1) + dx
       x_sum = x_sum + x(i)
    enddo

    dx = x_sum/real(n_out) - tlon
    if (dx < -180.) then
       do i = 1, n_out
          x(i) = x(i) + 360.
       enddo
    else if (dx > 180.) then
       do i = 1, n_out
          x(i) = x(i) - 360.
       enddo
    endif

    return

  end subroutine lon_fix

  !#######################################################################
  ! delete vertex
  subroutine vtx_delete(x,y, n, n_del)
    real, dimension(:), intent(inout) :: x,y
    integer,            intent(inout) :: n
    integer,            intent(in)    :: n_del

    integer :: i

    do i = n_del, n-1
       x(i) = x(i+1)
       y(i) = y(i+1)
    enddo
    n = n -1

    return 

  end subroutine vtx_delete

  !#######################################################################
  ! insert vertex
  subroutine vtx_insert(x,y, n, n_ins)
    real, dimension(:), intent(inout) :: x,y
    integer,            intent(inout) :: n
    integer,            intent(in)    :: n_ins

    integer :: i

    do i = n, n_ins, -1
       x(i+1) = x(i)
       y(i+1) = y(i)
    enddo
    n = n + 1

    return 

  end subroutine vtx_insert

  !#######################################################################

end module hgrid_mod

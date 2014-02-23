module atmos_grid_mod
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
  ! <TT>atmos_grid_mod</TT> Generate horizontal grid ( either bgrid or spectral grid )
  !     for atmosphere and land model. 
  !</OVERVIEW>

  !<DESCRIPTION>
  ! There are four subgrids, labeled T (for tracer), C (corner of T), N (north of T) 
  ! and E (east of T). The following schematic describes the grid cell notation.
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
  ! The grid_spec file would contains all of the following information.
  !
  !<PRE>
  !   x_T, y_T           = Geographic location of T-cell center
  !   x_vert_T, y_vert_T = Geographic location of T-cell vertices(each cell has 4 vertices)
  !   x_C, y_C           = Geographic location of C-cell center
  !
  !</PRE>
  !
  !</DESCRIPTION>

  use mpp_mod,        only : mpp_pe, mpp_root_pe, mpp_chksum, mpp_error, FATAL, NOTE
  use mpp_io_mod,     only : mpp_write_meta, mpp_write, axistype, fieldtype
  use fms_mod,        only : write_version_number, file_exist, stdout, field_size, string
  use fms_mod,        only : open_namelist_file, close_file, check_nml_error, read_data
  use transforms_mod, only : transforms_init, get_grid_boundaries, get_deg_lon, get_deg_lat 
  use constants_mod,  only : PI, RADIUS

  implicit none
  private

  !------ namelist interface ---------------------------------------------
  !<NAMELIST NAME="atmos_grid_nml">
  !<DATA NAME="grid_type" TYPE="character(len=24)" >
  !  type of grid. Its value can be "bgrid" or "spectral". Its default value is "bgrid". 
  !</DATA>
  !<DATA NAME="debug" TYPE="logical">
  ! control standard output.
  !</DATA>
  !<DATA NAME="num_lon" TYPE="integer">
  ! number of longitude points.
  !</DATA>
  !<DATA NAME="num_lat" TYPE="integer">
  ! number of latitude points.
  ! </DATA>
  !<DATA NAME="num_fourier" TYPE="integer">
  ! the fourier wavenumbers in the truncation are set equal to fourier_inc*m, where m = 0, 
  ! num_fourier therefore, the total number of fourier modes is num_fourier +1. 
  ! This namelist option is only for grid_type = 'spectral'
  ! </DATA>
  !<DATA NAME="num_spherical" TYPE="integer">
  ! the wavenumber increment (see num_fourier above). This namelist option is only
  ! for grid_type = 'spectral'
  ! </DATA>
  !<DATA NAME="fourier_inc" TYPE="integer">
  ! the maximum meridional wavenumber Retained meridional wavewnumbers are n = 0,
  ! num_spherical The total spherical wavenumber is L = n+m. This namelist option is only
  ! for grid_type = 'spectral'
  ! </DATA>
  !<DATA NAME="lon_begin, lon_end" TYPE="real" >
  !  range of the longitude. Default value is : lon_begin=0, lon_end=360. If you want to
  !  generate regional grid, you may need to set these namelists. Otherwise use default value. 
  !  This namelist option is only for grid_type = 'bgrid'
  ! </DATA>
  !<DATA NAME="lat_begin, lat_end" TYPE="real" >
  !  range of the latitude. Default value is : lat_begin=-90, lat_end=90. If you want to
  !  generate regional grid, you may need to set these namelists. Otherwise use default value. 
  !  This namelist option is only for grid_type = 'bgrid'
  ! </DATA>
  !</NAMELIST>
  character(len=24) :: grid_type = 'bgrid'
  logical           :: debug = .false.
  integer           :: num_lon = 0
  integer           :: num_lat = 0
  integer           :: num_fourier = 0
  integer           :: num_spherical = 0
  integer           :: fourier_inc = 1
  real              :: lon_begin = 0.0
  real              :: lon_end   = 360.0
  real              :: lat_begin = -90.0
  real              :: lat_end   = 90.0

  namelist /atmos_grid_nml/ debug, grid_type, num_lon, num_lat, num_fourier, num_spherical, &
       fourier_inc, lon_begin, lon_end, lat_begin, lat_end

  !<PUBLICTYPE >
  type atmos_grid_type
     real, dimension(:,:), pointer   :: x_T => NULL()      ! geographical longitude of T-cell center
     real, dimension(:,:), pointer   :: y_T => NULL()      ! geographical latitude of T-cell center 
     real, dimension(:,:), pointer   :: x_C => NULL()      ! geographical longitude of C-cell center
     real, dimension(:,:), pointer   :: y_C => NULL()      ! geographical latitude of C-cell center 
     real, dimension(:,:,:), pointer :: x_vert_T => NULL() ! geographical longitude of T-cell vertices
     real, dimension(:,:,:), pointer :: y_vert_T => NULL() ! geographical latitude of T-cell vertices
  end type atmos_grid_type
  !</PUBLICTYPE>

  !-----------------------------------------------------------------------
  !--------private data---------------------------------------------------
  logical :: module_is_initialized = .false.

  !---  axis and field data type for output file -----------------------
  type(axistype),save  :: axis_xt, axis_yt, axis_xc, axis_yc, axis_v
  type(fieldtype),save :: fld_x_T, fld_y_T, fld_x_vert_T, fld_y_vert_T, fld_x_C, fld_y_C

  !---------version information-------------------------------------------
  character(len=128) :: version = '$Id: atmos_grid.f90,v 10.0 2003/10/24 22:01:43 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 
  !---------public interface----------------------------------------------
  public :: generate_atmos_grid, atmos_grid_init, atmos_grid_end, atmos_grid_type
  public :: write_atmos_grid_meta, write_atmos_grid_data

contains

  !#######################################################################
  ! <SUBROUTINE NAME="atmos_grid_init" >
  !   <OVERVIEW>
  !    Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    Read namelist, write out version and namelist informaiton.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call atmos_grid_init ( )
  !   </TEMPLATE>
  ! </SUBROUTINE>

  subroutine atmos_grid_init 

    !--- local variables -------------------------------------------------
    integer ::  unit, ierr, io

    !---- read namelist --------------------------------------------------
    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, nml=atmos_grid_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'atmos_grid_nml')  ! also initializes nml error codes
       enddo
10     call close_file (unit)
    else
       call mpp_error(FATAL,'grid_generator: file input.nml does not exist')
    endif
    !--- write version info and namelist to logfile ----------------------

    call write_version_number(version,tagname)
    write (stdout(), nml=atmos_grid_nml)

    module_is_initialized = .true.

    !--- check namelist value
    if(trim(grid_type) .ne. 'bgrid' .and. trim(grid_type) .ne. 'spectral') &
         call mpp_error(FATAL,'atmos_grid_mod: '//trim(grid_type)//' is not a valid option of nml "grid_type" ')

    if(num_lon .le. 0 .or. num_lat .le. 0) then
       call mpp_error(FATAL, 'atmos_grid_mod: nml num_lon= '//trim(string(num_lon))//' and num_lat= ' &
            //trim(string(num_lat))//' should be positive number')
    endif

    if(trim(grid_type) == 'spectral') then
       if(num_fourier == 0 .or. num_spherical == 0) &
            call mpp_error(FATAL, 'atmos_grid_mod: nml "num_fourier"= '//trim(string(num_fourier))// &
            ' and "num_spherical"= '// trim(string(num_spherical))//' shoulbe be positive number')
    else   ! grid_type == bgrid
       if(lon_begin .ge. lon_end ) &
            call mpp_error(FATAL, 'atmos_grid_mod: nml lon_begin= '//trim(string(lon_begin))// &
            ' should be less than nml "lon_end"= '//trim(string(lon_end)) )
       if(lat_begin .ge. lat_end ) &
            call mpp_error(FATAL, 'atmos_grid_mod: nml lat_begin= '//trim(string(lat_begin))// &
            ' should be less than nml "lat_end"= '//trim(string(lat_end)) )
    endif

  end subroutine atmos_grid_init


  !#######################################################################
  ! <SUBROUTINE NAME="generate_atmos_grid" >
  !   <OVERVIEW>
  !    Generate horizontal grid.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    Calculate geographic locations of T and C-cell center. 
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call generat_atmos_grid (Hgrid)
  !   </TEMPLATE>
  !   <INOUT NAME="Hgrid" TYPE="atmos_grid_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </INOUT>
  ! </SUBROUTINE>

  subroutine generate_atmos_grid (Hgrid)

    type(atmos_grid_type), intent(inout) :: Hgrid

    select case(trim(grid_type))
    case('bgrid')
       call generate_bgrid (Hgrid)
    case('spectral')
       call generate_spectral_grid (Hgrid)
    end select

  end subroutine generate_atmos_grid

  !#######################################################################
  !--- generate spectral horizontal grid
  subroutine generate_spectral_grid (Hgrid)

    type(atmos_grid_type), intent(inout) :: Hgrid

    real    :: dlon, dlat, R2D 
    integer :: i, j, pe, root_pe
    integer :: chksum_x_t, chksum_y_t, chksum_x_c, chksum_y_c, chksum_x_vert_t, chksum_y_vert_t
    real, allocatable, dimension(:) :: lonb, latb, deg_lon, deg_lat

    R2D = 180./PI
    pe       = mpp_pe()
    root_pe  = mpp_root_pe()

    call transforms_init(RADIUS, num_lat, num_lon, num_fourier, fourier_inc, num_spherical)

    allocate(lonb(num_lon+1), latb(num_lat+1))
    allocate(deg_lon(num_lon), deg_lat(num_lat))
    call get_deg_lon(deg_lon)
    call get_deg_lat(deg_lat)
    call get_grid_boundaries(lonb, latb, global=.true.)
    allocate(hgrid%x_C(num_lon,num_lat), hgrid%y_C(num_lon,num_lat), &
         hgrid%x_T(num_lon,num_lat), hgrid%y_T(num_lon,num_lat), &
         hgrid%x_vert_T(num_lon,num_lat,4), hgrid%y_vert_T(num_lon,num_lat,4))

    do i=1,num_lon
       hgrid%x_C(i,:) = lonb(i+1)*R2D
       hgrid%x_T(i,:) = deg_lon(i)
       hgrid%x_vert_T(i,:,1) = lonb(i)*R2D
       hgrid%x_vert_T(i,:,2) = lonb(i+1)*R2D
       hgrid%x_vert_T(i,:,3) = lonb(i+1)*R2D
       hgrid%x_vert_T(i,:,4) = lonb(i)*R2D
    enddo

    do j=1,num_lat
       hgrid%y_C(:,j) = latb(j+1)*R2D
       hgrid%y_T(:,j) = deg_lat(j)
       hgrid%y_vert_T(:,j,1) = latb(j)*R2D
       hgrid%y_vert_T(:,j,2) = latb(j)*R2D
       hgrid%y_vert_T(:,j,3) = latb(j+1)*R2D
       hgrid%y_vert_T(:,j,4) = latb(j+1)*R2D
    enddo

    ! Print all grid coordinates

    write (stdout(),'(//,36x,a,/)') 'H G R I D   G E N E R A T I O N'
    write (stdout(),'(/)')
    write (stdout(),9105) num_lat
    write (stdout(),9001) (hgrid%y_T(1,j),j=1,num_lat)
    write (stdout(),9106) num_lat
    write (stdout(),9001) (hgrid%y_C(1,j),j=1,num_lat)
    write (stdout(),9107) num_lon
    write (stdout(),9001) (hgrid%x_T(i,1),i=1,num_lon)
    write (stdout(),9108) num_lon
    write (stdout(),9001) (hgrid%x_T(i,1),i=1,num_lon)
    write (stdout(),'(/)')

    if(debug) then
       chksum_x_t = mpp_chksum(Hgrid%x_T,(/pe/))
       chksum_y_t = mpp_chksum(Hgrid%y_T,(/pe/))
       chksum_x_c = mpp_chksum(Hgrid%x_C,(/pe/))
       chksum_y_c = mpp_chksum(Hgrid%y_C,(/pe/))
       chksum_x_vert_t = mpp_chksum(Hgrid%x_vert_T,(/pe/))
       chksum_y_vert_t = mpp_chksum(Hgrid%y_vert_T,(/pe/))
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_x_T = ', chksum_x_t
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_y_T = ', chksum_y_t
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_x_C = ', chksum_x_c
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_y_c = ', chksum_y_c
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_x_vert_T = ', chksum_x_vert_t
       write (stdout(),'(/(a,I20/))')'spectral_hgrid chksum: chksum_y_vert_T = ', chksum_y_vert_t
    endif

    deallocate(lonb, latb, deg_lon, deg_lat)
    return

9001 format (1x,10f10.4)
9105 format (/,' Latitude of T points (deg): yt(j) j=1,',i4)
9106 format (/,' Latitude of C points (deg): yu(j) j=1,',i4)
9107 format (/,' Longitude of T points (deg): xt(i) i=1,',i4)
9108 format (/,' Longitude of C points (deg): xu(i) i=1,',i4)
  end subroutine generate_spectral_grid

  !#####################################################################
  !--- generate horizontal bgrid
  subroutine generate_bgrid (hgrid)

    type(atmos_grid_type), intent(inout) :: hgrid

    real    :: dlon, dlat
    integer :: i, j, pe
    integer :: chksum_x_t, chksum_y_t, chksum_x_c, chksum_y_c, chksum_x_vert_t, chksum_y_vert_t

    pe       = mpp_pe()

    allocate(hgrid%x_C(num_lon,num_lat), hgrid%y_C(num_lon,num_lat), &
         hgrid%x_T(num_lon,num_lat), hgrid%y_T(num_lon,num_lat), &
         hgrid%x_vert_T(num_lon,num_lat,4), hgrid%y_vert_T(num_lon,num_lat,4))

    dlon = (lon_end-lon_begin)/num_lon
    dlat = (lat_end-lat_begin)/num_lat

    do i=1,num_lon
       hgrid%x_C(i,:)        = lon_begin + dlon*real(i)
       hgrid%x_T(i,:)        = lon_begin + dlon*(i-.5)
       hgrid%x_vert_T(i,:,1) = lon_begin + dlon*(i-1.0)
       hgrid%x_vert_T(i,:,2) = lon_begin + dlon*real(i)
       hgrid%x_vert_T(i,:,3) = lon_begin + dlon*real(i)
       hgrid%x_vert_T(i,:,4) = lon_begin + dlon*(i-1.0)
    enddo

    do j=1,num_lat
       hgrid%y_C(:,j)        = lat_begin + dlat*real(j)
       hgrid%y_T(:,j)        = lat_begin + dlat*(j-.5) 
       hgrid%y_vert_T(:,j,1) = lat_begin + dlat*real(j-1.0)
       hgrid%y_vert_T(:,j,2) = lat_begin + dlat*real(j-1.0)
       hgrid%y_vert_T(:,j,3) = lat_begin + dlat*real(j)
       hgrid%y_vert_T(:,j,4) = lat_begin + dlat*real(j)

    enddo

    ! Print all grid coordinates

    write (stdout(),'(//,36x,a,/)') 'H G R I D   G E N E R A T I O N'
    write (stdout(),'(/)')
    write (stdout(),9105) num_lat
    write (stdout(),9001) (hgrid%y_T(1,j),j=1,num_lat)
    write (stdout(),9106) num_lat
    write (stdout(),9001) (hgrid%y_C(1,j),j=1,num_lat)
    write (stdout(),9107) num_lon
    write (stdout(),9001) (hgrid%x_T(i,1),i=1,num_lon)
    write (stdout(),9108) num_lon
    write (stdout(),9001) (hgrid%x_C(i,1),i=1,num_lon)
    write (stdout(),'(/)')

    if(debug ) then
       chksum_x_t = mpp_chksum(Hgrid%x_T,(/pe/))
       chksum_y_t = mpp_chksum(Hgrid%y_T,(/pe/))
       chksum_x_c = mpp_chksum(Hgrid%x_C,(/pe/))
       chksum_y_c = mpp_chksum(Hgrid%y_C,(/pe/))
       chksum_x_vert_t = mpp_chksum(Hgrid%x_vert_T,(/pe/))
       chksum_y_vert_t = mpp_chksum(Hgrid%y_vert_T,(/pe/))
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_x_T = ', chksum_x_t
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_y_T = ', chksum_y_t
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_x_C = ', chksum_x_c
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_y_c = ', chksum_y_c
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_x_vert_T = ', chksum_x_vert_t
       write (stdout(),'(/(a,I20/))')'bgrid_hgrid chksum: chksum_y_vert_T = ', chksum_y_vert_t
    endif

    return
9001 format (1x,10f10.4)
9105 format (/,' Latitude of T points (deg): yt(j) j=1,',i4)
9106 format (/,' Latitude of C points (deg): yu(j) j=1,',i4)
9107 format (/,' Longitude of T points (deg): xt(i) i=1,',i4)
9108 format (/,' Longitude of C points (deg): xu(i) i=1,',i4)

  end subroutine generate_bgrid

  !#######################################################################
  ! <SUBROUTINE NAME="write_atmos_grid_data">
  !   <OVERVIEW>
  !     write the Hgrid data to netcdf file
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_atmos_grid_data (unit,Hgrid)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Hgrid" TYPE="atmos_grid_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_atmos_grid_data (unit,Hgrid)
    integer,               intent(in) :: unit
    type(atmos_grid_type), intent(in) :: Hgrid

    !--- write out axis data
    call mpp_write(unit, axis_xt)
    call mpp_write(unit, axis_yt)
    call mpp_write(unit, axis_xc)
    call mpp_write(unit, axis_yc)
    call mpp_write(unit, axis_v)

    call mpp_write(unit, fld_x_T, hgrid%x_T)
    call mpp_write(unit, fld_y_T, hgrid%y_T)
    call mpp_write(unit, fld_x_C, hgrid%x_C)
    call mpp_write(unit, fld_y_C, hgrid%y_C)
    call mpp_write(unit, fld_x_vert_T, hgrid%x_vert_T)
    call mpp_write(unit, fld_y_vert_T, hgrid%y_vert_T)  

    return

  end subroutine write_atmos_grid_data

  !#######################################################################
  ! <SUBROUTINE NAME="write_atmos_grid_meta">

  !   <OVERVIEW>
  !     Write out horizontal grid meta data.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_atmos_grid_meta(unit, Hgrid)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Hgrid" TYPE="atmos_grid_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_atmos_grid_meta(unit, Hgrid )
    integer,               intent(in) :: unit
    type(atmos_grid_type), intent(in) :: Hgrid

    !--- write out axis meta ----------------------------------------
    call mpp_write_meta(unit, axis_xt,'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = hgrid%x_T(:,1) )
    call mpp_write_meta(unit, axis_yt,'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = hgrid%y_T(1,:) )
    call mpp_write_meta(unit, axis_xc,'grid_x_C','degree_east','Nominal Longitude of C-cell center', &
         cartesian ='X', data = hgrid%x_C(:,1) )
    call mpp_write_meta(unit, axis_yc,'grid_y_C','degree_north','Nominal Latitude of C-cell center', &
         cartesian ='Y', data = hgrid%y_C(1,:) )
    call mpp_write_meta(unit, axis_v, 'vertex', 'none ', 'Vertex position from southwest couterclockwise', &
         data = (/1.,2.,3.,4./) )
    call mpp_write_meta(unit, fld_x_T, (/axis_xt, axis_yt/), 'x_T', 'degree_east',  &
         'Geographic longitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_T, (/axis_xt, axis_yt/), 'y_T', 'degree_north',  &
         'Geographic latitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_x_C, (/axis_xc, axis_yc/), 'x_C', 'degree_east',  &
         'Geographic longitude of C_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_C, (/axis_xc, axis_yc/), 'y_C', 'degree_north',  &
         'Geographic latitude of C_cell centers', pack=1)
    call mpp_write_meta(unit, fld_x_vert_T, (/axis_xt, axis_yt, axis_v/), 'x_vert_T', &
         'degree_east', 'Geographic longitude of T_cell vertices begin southwest counterclockwise', pack=1)
    call mpp_write_meta(unit, fld_y_vert_T, (/axis_xt, axis_yt, axis_v/), 'y_vert_T', &
         'degree_north', 'Geographic latitude of T_cell vertices begin southwest counterclockwise', pack=1)

    return
  end subroutine write_atmos_grid_meta

  !#######################################################################
  ! <SUBROUTINE NAME="atmos_grid_end">
  !   <OVERVIEW>
  !     Destruction routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "atmos_grid_type" variables.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call atmos_grid_end ( Hgrid )
  !   </TEMPLATE>
  !   <INOUT NAME="Hgrid" TYPE="atmos_grid_type">
  !     A derived-type variable that contains horizontal grid information.
  !   </INOUT>
  ! </SUBROUTINE>

  subroutine atmos_grid_end(Hgrid)

    type(atmos_grid_type), intent(inout) :: Hgrid
    !--- release memory of module variables ------------------------------

    deallocate(Hgrid%x_C, Hgrid%y_C, Hgrid%x_T, Hgrid%y_T, Hgrid%x_vert_T, Hgrid%y_vert_T)

    module_is_initialized = .false.

    return

  end subroutine atmos_grid_end

  !#######################################################################

end module atmos_grid_mod

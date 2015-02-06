program grid_transfer
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
  !  Converts the grid specification netcdf file using old name convention to a netcdf 
  !  file with new grid name convention. This program can be applied on any ocean grid 
  !  or exchange grid. 
  !</OVERVIEW>
  use fms_mod,       only : fms_init, fms_end, file_exist, close_file, open_namelist_file
  use fms_mod,       only : check_nml_error, write_version_number, stdout
  use mpp_mod,       only : mpp_error, FATAL, mpp_npes
  use constants_mod, only : constants_init, PI

  implicit none
#include <netcdf.inc>

  !--- namelist interface
  !
  !<NAMELIST NAME="grid_transfer_nml">
  !<DATA NAME="old_grid" TYPE="character(len=128)">
  !  name of input grid file to be converted.
  !</DATA>
  !<DATA NAME="new_grid" TYPE="character(len=128)">
  !  name of output grid file converted from old_grid.
  !</DATA>
  !</NAMELIST>
  character(len=128) :: old_grid = 'old_grid.nc'
  character(len=128) :: new_grid = 'new_grid.nc'
  namelist /grid_transfer_nml/ old_grid, new_grid

  !---------version information-------------------------------------------
  character(len=128) :: version = '$Id: grid_transfer.F90,v 13.0 2006/03/28 21:44:17 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 

  !---------------------------------------------------------------------
  integer                             :: nlon, nlat, nk
  real, dimension(:), allocatable     :: grid_x_T, grid_y_T, grid_x_c, grid_y_c, zb
  real, dimension(:,:), allocatable   :: x_T, y_T, x_c, y_c, angle_c, depth_t, num_levels
  real, dimension(:,:), allocatable   :: ds_11_21_t, ds_01_11_t, ds_11_12_t, ds_10_11_t 
  real, dimension(:,:), allocatable   :: ds_01_21_t, ds_10_12_t, ds_20_22_t, ds_02_22_t
  real, dimension(:,:), allocatable   :: ds_00_20_c, ds_00_02_c, ds_11_21_c, ds_01_11_c
  real, dimension(:,:), allocatable   :: ds_11_12_c, ds_10_11_c, ds_01_21_c, ds_10_12_c
  real, dimension(:,:), allocatable   :: ds_02_22_c, ds_20_22_c, area_t
  real, dimension(:,:), allocatable   :: ds_01_21_e, ds_10_12_n
  real, dimension(:,:,:), allocatable :: x_vert_T, y_vert_T


  !--- begin of the program 
  call fms_init
  call constants_init

  !--- read the namelist
  call grid_transfer_init

  !--- read the old_grid file.
  call read_old_grid

  !--- write the new_grid file
  call write_new_grid

  call grid_transfer_end

  call fms_end

contains

  !#####################################################################
  !--- initializtion routine
  !--- Read namelist and write out namelist and version information
  subroutine grid_transfer_init

    integer :: unit, ierr, io_status

    if(mpp_npes() .ne. 1)&
         call mpp_error(FATAL,' grid_transfer: do not support parallel, number of pes should be 1')

    if(file_exist('input.nml')) then
       unit = open_namelist_file()
       read (unit,grid_transfer_nml,IOSTAT=io_status)
       write (stdout(),'(/)')
       write (stdout(),grid_transfer_nml)  
       ierr = check_nml_error(io_status, 'grid_transfer_nml')
       call close_file(unit)
    else
       call mpp_error(FATAL, 'grid_transfer: file input.nml does not exist' )
    endif

    !--- write version information
    call write_version_number(version, tagname)

  end subroutine grid_transfer_init

  !#####################################################################
  !--- read_old_grid
  subroutine read_old_grid
    integer :: ncid, rcode, start(4), nread(4), dims(4)
    integer :: id_geolon_t, id_geolat_t, id_geolon_c, id_geolat_c
    integer :: id_geolon_vert_t, id_geolat_vert_t, id_sin_rot, id_ht, id_kmt
    integer :: id_zw, id_gridlon_t, id_gridlat_t, id_gridlon_c, id_gridlat_c
    integer :: id_dte, id_dtw, id_dtn, id_dts, id_dxt, id_dyt, id_dyte, id_dxtn
    integer :: id_dxte, id_dytn, id_due, id_duw, id_dun, id_dus, id_dxu, id_dyu
    integer :: id_dxun, id_dyue
    real, dimension(:,:), allocatable :: geolon_vert_t, geolat_vert_t, sin_rot  

    rcode = nf_open(trim(old_grid), NF_NOWRITE, ncid)
    if(rcode /= 0) call mpp_error(FATAL, 'Cannot open file '//trim(old_grid))
    !--- first get the dimension of the grid
    rcode = nf_inq_varid(ncid, 'geolon_t', id_geolon_t )
    if(rcode /= 0) call mpp_error(FATAL, ' field geolon_t is not in the file' //trim(old_grid) )
    rcode = nf_inq_vardimid(ncid, id_geolon_t, dims)
    rcode = nf_inq_dimlen(ncid, dims(1), nlon)
    rcode = nf_inq_dimlen(ncid, dims(2), nlat)
    rcode = nf_inq_varid(ncid, 'zw', id_zw )
    if(rcode /= 0) call mpp_error(FATAL, 'zw is not in the file' //trim(old_grid) )
    rcode = nf_inq_vardimid(ncid, id_zw, dims)
    rcode = nf_inq_dimlen(ncid, dims(1), nk)    
    !--- get the axis data
    allocate(grid_x_t(nlon), grid_y_t(nlat), grid_x_c(nlon), grid_y_c(nlat), zb(nk) )
    rcode = nf_inq_varid(ncid, 'gridlon_t', id_gridlon_t )  
    if(rcode /= 0) call mpp_error(FATAL, 'gridlon_t is not in the file' //trim(old_grid) )  
    rcode = nf_inq_varid(ncid, 'gridlat_t', id_gridlat_t )
    if(rcode /= 0) call mpp_error(FATAL, 'gridlat_t is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'gridlon_c', id_gridlon_c )
    if(rcode /= 0) call mpp_error(FATAL, 'gridlon_c is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'gridlat_c', id_gridlat_c )
    if(rcode /= 0) call mpp_error(FATAL, 'gridlat_c is not in the file' //trim(old_grid) )
    start = 1; nread = 1
    nread(1) = nlon
    rcode = nf_get_vara_double(ncid, id_gridlon_t, start, nread,grid_x_t )
    rcode = nf_get_vara_double(ncid, id_gridlon_c, start, nread,grid_x_c )
    nread(1) = nlat
    rcode = nf_get_vara_double(ncid, id_gridlat_t, start, nread,grid_y_t )
    rcode = nf_get_vara_double(ncid, id_gridlat_c, start, nread,grid_y_c )
    nread(1) = nk
    rcode = nf_get_vara_double(ncid, id_zw, start, nread,zb )

    !--- get the field data

    allocate(x_t(nlon,nlat), y_t(nlon,nlat), x_c(nlon,nlat), y_c(nlon,nlat))
    allocate(x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4))
    allocate(geolon_vert_t(nlon+1,nlat+1), geolat_vert_t(nlon+1,nlat+1) )
    allocate(ds_11_21_t(nlon,nlat), ds_01_11_t(nlon,nlat), ds_11_12_t(nlon,nlat), ds_10_11_t(nlon,nlat) )
    allocate(ds_01_21_t(nlon,nlat), ds_10_12_t(nlon,nlat), ds_20_22_t(nlon,nlat), ds_02_22_t(nlon,nlat) )
    allocate(ds_00_20_c(nlon,nlat), ds_00_02_c(nlon,nlat), ds_11_21_c(nlon,nlat), ds_01_11_c(nlon,nlat) )
    allocate(ds_11_12_c(nlon,nlat), ds_10_11_c(nlon,nlat), ds_01_21_c(nlon,nlat), ds_10_12_c(nlon,nlat) )
    allocate(ds_02_22_c(nlon,nlat), ds_20_22_c(nlon,nlat), area_t(nlon,nlat) )
    allocate(ds_10_12_n(nlon,nlat), ds_01_21_e(nlon,nlat) )
    allocate(angle_c(nlon,nlat), sin_rot(nlon,nlat), depth_t(nlon,nlat), num_levels(nlon,nlat)  )

    rcode = nf_inq_varid(ncid, 'geolat_t', id_geolat_t )
    if(rcode /= 0) call mpp_error(FATAL, 'geolat_t is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'geolon_c', id_geolon_c )
    if(rcode /= 0) call mpp_error(FATAL, 'geolon_c is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'geolat_c', id_geolat_c )
    if(rcode /= 0) call mpp_error(FATAL, 'geolat_c is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'geolon_vert_t', id_geolon_vert_t )
    if(rcode /= 0) call mpp_error(FATAL, 'geolon_vert_t is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'geolat_vert_t', id_geolat_vert_t )
    if(rcode /= 0) call mpp_error(FATAL, 'geolat_vert_t is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'sin_rot', id_sin_rot )
    if(rcode /= 0) call mpp_error(FATAL, 'sin_rot is not in the file' //trim(old_grid) )

    nread(1) = nlon; nread(2) = nlat
    rcode = nf_get_vara_double(ncid, id_geolon_t, start, nread, x_t)
    rcode = nf_get_vara_double(ncid, id_geolat_t, start, nread, y_t)
    rcode = nf_get_vara_double(ncid, id_geolon_c, start, nread, x_c)
    rcode = nf_get_vara_double(ncid, id_geolat_c, start, nread, y_c)
    rcode = nf_get_vara_double(ncid, id_sin_rot, start, nread, sin_rot)
    nread(1) = nlon+1; nread(2) = nlat+1
    rcode = nf_get_vara_double(ncid, id_geolon_vert_t, start, nread, geolon_vert_t )
    rcode = nf_get_vara_double(ncid, id_geolat_vert_t, start, nread, geolat_vert_t )

    !--- define x_vert_t and y_vert_t
    x_vert_t(:,:,1) = geolon_vert_t(1:nlon,   1:nlat  )
    x_vert_t(:,:,2) = geolon_vert_t(2:nlon+1, 1:nlat   )
    x_vert_t(:,:,3) = geolon_vert_t(2:nlon+1, 2:nlat+1)
    x_vert_t(:,:,4) = geolon_vert_t(1:nlon,   2:nlat+1)
    y_vert_t(:,:,1) = geolat_vert_t(1:nlon,   1:nlat  )
    y_vert_t(:,:,2) = geolat_vert_t(2:nlon+1, 1:nlat  )
    y_vert_t(:,:,3) = geolat_vert_t(2:nlon+1, 2:nlat+1)
    y_vert_t(:,:,4) = geolat_vert_t(1:nlon,   2:nlat+1)

    !--- define rotation angle
    angle_c = asin(sin_rot)*180./PI

    !--- get the grid resolution variables
    rcode = nf_inq_varid(ncid, 'dte', id_dte )
    if(rcode /= 0) call mpp_error(FATAL, 'dte is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dtw', id_dtw )
    if(rcode /= 0) call mpp_error(FATAL, 'dtw is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dtn', id_dtn )
    if(rcode /= 0) call mpp_error(FATAL, 'dtn is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dts', id_dts )
    if(rcode /= 0) call mpp_error(FATAL, 'dts is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxt', id_dxt )
    if(rcode /= 0) call mpp_error(FATAL, 'dxt is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dyt', id_dyt )
    if(rcode /= 0) call mpp_error(FATAL, 'dyt is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dyte', id_dyte )
    if(rcode /= 0) call mpp_error(FATAL, 'dyte is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxtn', id_dxtn )
    if(rcode /= 0) call mpp_error(FATAL, 'dxtn is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxte', id_dxte )
    if(rcode /= 0) call mpp_error(FATAL, 'dxte is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dytn', id_dytn )
    if(rcode /= 0) call mpp_error(FATAL, 'dytn is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'due', id_due )
    if(rcode /= 0) call mpp_error(FATAL, 'due is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'duw', id_duw )
    if(rcode /= 0) call mpp_error(FATAL, 'duw is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dun', id_dun )
    if(rcode /= 0) call mpp_error(FATAL, 'dun is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dus', id_dus )
    if(rcode /= 0) call mpp_error(FATAL, 'dus is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxu', id_dxu )
    if(rcode /= 0) call mpp_error(FATAL, 'dxu is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dyu', id_dyu )
    if(rcode /= 0) call mpp_error(FATAL, 'dyu is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxun', id_dxun )
    if(rcode /= 0) call mpp_error(FATAL, 'dxun is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dyue', id_dyue )
    if(rcode /= 0) call mpp_error(FATAL, 'dyue is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dxte', id_dxte )
    if(rcode /= 0) call mpp_error(FATAL, 'dxte is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'dytn', id_dytn )
    if(rcode /= 0) call mpp_error(FATAL, 'dytn is not in the file' //trim(old_grid) )

    !--- get the grid resolution data
    start=1; nread=1
    nread(1)=nlon; nread(2) = nlat
    rcode = nf_get_vara_double(ncid, id_dte, start, nread, ds_11_21_t)
    rcode = nf_get_vara_double(ncid, id_dtw, start, nread, ds_01_11_t)
    rcode = nf_get_vara_double(ncid, id_dtn, start, nread, ds_11_12_t)
    rcode = nf_get_vara_double(ncid, id_dts, start, nread, ds_10_11_t)
    rcode = nf_get_vara_double(ncid, id_dxt, start, nread, ds_01_21_t)
    rcode = nf_get_vara_double(ncid, id_dyt, start, nread, ds_10_12_t)
    rcode = nf_get_vara_double(ncid, id_dyte, start, nread, ds_20_22_t)
    rcode = nf_get_vara_double(ncid, id_dxtn, start, nread, ds_02_22_t)
    rcode = nf_get_vara_double(ncid, id_dxte, start, nread, ds_00_20_c)
    rcode = nf_get_vara_double(ncid, id_dytn, start, nread, ds_00_02_c)
    rcode = nf_get_vara_double(ncid, id_due, start, nread, ds_11_21_c)
    rcode = nf_get_vara_double(ncid, id_duw, start, nread, ds_01_11_c)
    rcode = nf_get_vara_double(ncid, id_dun, start, nread, ds_11_12_c)
    rcode = nf_get_vara_double(ncid, id_dus, start, nread, ds_10_11_c)
    rcode = nf_get_vara_double(ncid, id_dxu, start, nread, ds_01_21_c)
    rcode = nf_get_vara_double(ncid, id_dyu, start, nread, ds_10_12_c)
    rcode = nf_get_vara_double(ncid, id_dxun, start, nread, ds_02_22_c)
    rcode = nf_get_vara_double(ncid, id_dyue, start, nread, ds_20_22_c)
    rcode = nf_get_vara_double(ncid, id_dxte, start, nread, ds_01_21_e)
    rcode = nf_get_vara_double(ncid, id_dytn, start, nread, ds_10_12_n)

    !--- calculate area of T-cell, which is dxt*dyt
    area_t =  ds_01_21_t * ds_10_12_t

    !--- get the topography data
    rcode = nf_inq_varid(ncid, 'ht', id_ht )
    if(rcode /= 0) call mpp_error(FATAL, 'ht is not in the file' //trim(old_grid) )
    rcode = nf_inq_varid(ncid, 'kmt', id_kmt )
    if(rcode /= 0) call mpp_error(FATAL, 'kmt is not in the file' //trim(old_grid) )

    rcode = nf_get_vara_double(ncid, id_ht, start, nread, depth_t)
    rcode = nf_get_vara_double(ncid, id_kmt, start, nread, num_levels)


  end subroutine read_old_grid

  !#####################################################################
  !--- write out the new_grid
  subroutine write_new_grid

    integer :: ncid, rcode, dims(10), start(4), nwrite(4), dim(4)
    integer :: id_grid_x_t, id_grid_y_t, id_grid_x_c, id_grid_y_c, id_zb, id_vertex
    integer :: id_x_t, id_y_t, id_x_c, id_y_c, id_x_vert_t, id_y_vert_t
    integer :: id_ds_11_21_t, id_ds_01_11_t, id_ds_11_12_t, id_ds_10_11_t
    integer :: id_ds_01_21_t, id_ds_10_12_t, id_ds_20_22_t, id_ds_02_22_t
    integer :: id_ds_00_20_c, id_ds_00_02_c, id_ds_11_21_c, id_ds_01_11_c
    integer :: id_ds_11_12_c, id_ds_10_11_c, id_ds_01_21_c, id_ds_10_12_c
    integer :: id_ds_02_22_c, id_ds_20_22_c, id_area_t, id_ds_01_21_e, id_ds_10_12_n
    integer :: id_angle_c, id_depth_t, id_num_levels


    call system("cp "//trim(old_grid)//" "//trim(new_grid) )
    rcode = nf_open(trim(new_grid), NF_WRITE, ncid)
    rcode = nf_redef(ncid)
    rcode = nf_def_dim(ncid, 'grid_x_T', nlon, dims(1))
    rcode = nf_def_dim(ncid, 'grid_y_T', nlat, dims(2))
    rcode = nf_def_dim(ncid, 'vertex',      4, dims(3))
    rcode = nf_def_dim(ncid, 'zb',      nk, dims(4))
    rcode = nf_def_dim(ncid, 'grid_x_C', nlon, dims(5))
    rcode = nf_def_dim(ncid, 'grid_y_C', nlat, dims(6))
    rcode = nf_def_var(ncid, "grid_x_T", NF_FLOAT, 1, dims(1), id_grid_x_t)
    rcode = nf_put_att_text(ncid,id_grid_x_t, 'cartesian_axis', 1,'X')
    rcode = nf_put_att_text(ncid,id_grid_x_t, 'units', 11,'degree_east')
    rcode = nf_def_var(ncid, "grid_y_T", NF_FLOAT, 1, dims(2), id_grid_y_t)
    rcode = nf_put_att_text(ncid,id_grid_y_t, 'cartesian_axis', 1,'Y')
    rcode = nf_put_att_text(ncid,id_grid_y_t, 'units', 12,'degree_north')
    rcode = nf_def_var(ncid, "grid_x_C", NF_FLOAT, 1, dims(5), id_grid_x_c)
    rcode = nf_put_att_text(ncid,id_grid_x_c, 'cartesian_axis', 1,'X')
    rcode = nf_put_att_text(ncid,id_grid_x_c, 'units', 11,'degree_east')
    rcode = nf_def_var(ncid, "grid_y_C", NF_FLOAT, 1, dims(6), id_grid_y_c)
    rcode = nf_put_att_text(ncid,id_grid_y_c, 'cartesian_axis', 1,'Y')
    rcode = nf_put_att_text(ncid,id_grid_y_c, 'units', 12,'degree_north')
    rcode = nf_def_var(ncid, "vertex", NF_FLOAT, 1, dims(3), id_vertex)
    rcode = nf_def_var(ncid, "zb", NF_FLOAT, 1, dims(4), id_zb)
    rcode = nf_put_att_text(ncid,id_zb, 'cartesian_axis', 1,'Z')
    rcode = nf_put_att_text(ncid,id_zb, 'units', 6,'meters')

    rcode = nf_def_var(ncid, 'x_T', NF_DOUBLE, 2, dims(1:2), id_x_t)
    rcode = nf_def_var(ncid, 'y_T', NF_DOUBLE, 2, dims(1:2), id_y_t)
    rcode = nf_def_var(ncid, 'x_C', NF_DOUBLE, 2, dims(5:6), id_x_c)
    rcode = nf_def_var(ncid, 'y_C', NF_DOUBLE, 2, dims(5:6), id_y_c)
    rcode = nf_def_var(ncid, 'x_vert_T', NF_DOUBLE, 3, dims(1:3), id_x_vert_t)
    rcode = nf_def_var(ncid, 'y_vert_T', NF_DOUBLE, 3, dims(1:3), id_y_vert_t)
    rcode = nf_def_var(ncid, 'area_T', NF_DOUBLE, 2, dims(1:2), id_area_t)
    rcode = nf_def_var(ncid, 'ds_11_21_T', NF_DOUBLE, 2, dims(1:2), id_ds_11_21_t)
    rcode = nf_def_var(ncid, 'ds_01_11_T', NF_DOUBLE, 2, dims(1:2), id_ds_01_11_t)
    rcode = nf_def_var(ncid, 'ds_11_12_T', NF_DOUBLE, 2, dims(1:2), id_ds_11_12_t)
    rcode = nf_def_var(ncid, 'ds_10_11_T', NF_DOUBLE, 2, dims(1:2), id_ds_10_11_t)
    rcode = nf_def_var(ncid, 'ds_01_21_T', NF_DOUBLE, 2, dims(1:2), id_ds_01_21_t)
    rcode = nf_def_var(ncid, 'ds_10_12_T', NF_DOUBLE, 2, dims(1:2), id_ds_10_12_t)
    rcode = nf_def_var(ncid, 'ds_20_22_T', NF_DOUBLE, 2, dims(1:2), id_ds_20_22_t)
    rcode = nf_def_var(ncid, 'ds_02_22_T', NF_DOUBLE, 2, dims(1:2), id_ds_02_22_t)
    rcode = nf_def_var(ncid, 'ds_00_20_C', NF_DOUBLE, 2, dims(5:6), id_ds_00_20_c)
    rcode = nf_def_var(ncid, 'ds_00_02_C', NF_DOUBLE, 2, dims(5:6), id_ds_00_02_c)
    rcode = nf_def_var(ncid, 'ds_11_21_C', NF_DOUBLE, 2, dims(5:6), id_ds_11_21_c)
    rcode = nf_def_var(ncid, 'ds_01_11_C', NF_DOUBLE, 2, dims(5:6), id_ds_01_11_c)
    rcode = nf_def_var(ncid, 'ds_11_12_C', NF_DOUBLE, 2, dims(5:6), id_ds_11_12_c)
    rcode = nf_def_var(ncid, 'ds_10_11_C', NF_DOUBLE, 2, dims(5:6), id_ds_10_11_c)
    rcode = nf_def_var(ncid, 'ds_01_21_C', NF_DOUBLE, 2, dims(5:6), id_ds_01_21_c)
    rcode = nf_def_var(ncid, 'ds_10_12_C', NF_DOUBLE, 2, dims(5:6), id_ds_10_12_c)
    rcode = nf_def_var(ncid, 'ds_02_22_C', NF_DOUBLE, 2, dims(5:6), id_ds_02_22_c)
    rcode = nf_def_var(ncid, 'ds_20_22_C', NF_DOUBLE, 2, dims(5:6), id_ds_20_22_c)
    dim(1) = dims(5); dim(2) = dims(2)
    rcode = nf_def_var(ncid, 'ds_01_21_E', NF_DOUBLE, 2, dim,       id_ds_01_21_e)
    dim(1) = dims(1); dim(2) = dims(6)
    rcode = nf_def_var(ncid, 'ds_10_12_N', NF_DOUBLE, 2, dim,       id_ds_10_12_n)
    rcode = nf_def_var(ncid, 'angle_C', NF_DOUBLE, 2, dims(5:6), id_angle_c)
    rcode = nf_def_var(ncid, 'depth_t', NF_DOUBLE, 2, dims(1:2), id_depth_t)
    rcode = nf_def_var(ncid, 'num_levels', NF_DOUBLE, 2, dims(1:2), id_num_levels)

    rcode = nf_enddef(ncid)

    start = 1; nwrite = 1;
    nwrite(1) = nlon
    rcode = nf_put_vara_double(ncid, id_grid_x_t, start, nwrite, grid_x_t)
    rcode = nf_put_vara_double(ncid, id_grid_x_c, start, nwrite, grid_x_c)
    nwrite(1) = nlat
    rcode = nf_put_vara_double(ncid, id_grid_y_t, start, nwrite, grid_y_t)
    rcode = nf_put_vara_double(ncid, id_grid_y_c, start, nwrite, grid_y_c)
    nwrite(1) = nk
    rcode = nf_put_vara_double(ncid, id_zb, start, nwrite, zb)

    nwrite(1) = nlon; nwrite(2) = nlat
    rcode = nf_put_vara_double(ncid, id_x_t, start, nwrite, x_t)
    rcode = nf_put_vara_double(ncid, id_y_t, start, nwrite, y_t) 
    rcode = nf_put_vara_double(ncid, id_x_c, start, nwrite, x_c)
    rcode = nf_put_vara_double(ncid, id_y_c, start, nwrite, y_c) 
    nwrite(3) = 4   
    rcode = nf_put_vara_double(ncid, id_x_vert_t, start, nwrite, x_vert_t)
    rcode = nf_put_vara_double(ncid, id_y_vert_t, start, nwrite, y_vert_t)

    nwrite = 1; nwrite(1) = nlon; nwrite(2) = nlat
    rcode = nf_put_vara_double(ncid, id_area_t, start, nwrite, area_t)
    rcode = nf_put_vara_double(ncid, id_ds_11_21_t, start, nwrite, ds_11_21_t)
    rcode = nf_put_vara_double(ncid, id_ds_01_11_t, start, nwrite, ds_01_11_t)
    rcode = nf_put_vara_double(ncid, id_ds_11_12_t, start, nwrite, ds_11_12_t)
    rcode = nf_put_vara_double(ncid, id_ds_10_11_t, start, nwrite, ds_10_11_t)
    rcode = nf_put_vara_double(ncid, id_ds_01_21_t, start, nwrite, ds_01_21_t)
    rcode = nf_put_vara_double(ncid, id_ds_10_12_t, start, nwrite, ds_10_12_t)
    rcode = nf_put_vara_double(ncid, id_ds_20_22_t, start, nwrite, ds_20_22_t)
    rcode = nf_put_vara_double(ncid, id_ds_02_22_t, start, nwrite, ds_02_22_t)
    rcode = nf_put_vara_double(ncid, id_ds_00_20_c, start, nwrite, ds_00_20_c)
    rcode = nf_put_vara_double(ncid, id_ds_00_02_c, start, nwrite, ds_00_02_c)
    rcode = nf_put_vara_double(ncid, id_ds_11_21_c, start, nwrite, ds_11_21_c)
    rcode = nf_put_vara_double(ncid, id_ds_01_11_c, start, nwrite, ds_01_11_c)
    rcode = nf_put_vara_double(ncid, id_ds_11_12_c, start, nwrite, ds_11_12_c)
    rcode = nf_put_vara_double(ncid, id_ds_10_11_c, start, nwrite, ds_10_11_c)
    rcode = nf_put_vara_double(ncid, id_ds_01_21_c, start, nwrite, ds_01_21_c)
    rcode = nf_put_vara_double(ncid, id_ds_10_12_c, start, nwrite, ds_10_12_c)
    rcode = nf_put_vara_double(ncid, id_ds_02_22_c, start, nwrite, ds_02_22_c)
    rcode = nf_put_vara_double(ncid, id_ds_20_22_c, start, nwrite, ds_20_22_c)
    rcode = nf_put_vara_double(ncid, id_ds_01_21_e, start, nwrite, ds_01_21_e)
    rcode = nf_put_vara_double(ncid, id_ds_10_12_n, start, nwrite, ds_10_12_n)
    rcode = nf_put_vara_double(ncid, id_angle_c, start, nwrite, angle_c )
    rcode = nf_put_vara_double(ncid, id_depth_t, start, nwrite, depth_t)
    rcode = nf_put_vara_double(ncid, id_num_levels, start, nwrite, num_levels)

    rcode = nf_close(ncid)

  end   subroutine write_new_grid

  !#####################################################################
  !--- release the memory
  subroutine grid_transfer_end

    deallocate( grid_x_T, grid_y_T, grid_x_c, grid_y_c, zb, x_T, y_T, x_c, y_c )
    deallocate( x_vert_T, y_vert_T, angle_c, depth_t, num_levels )
    deallocate( ds_11_21_t, ds_01_11_t, ds_11_12_t, ds_10_11_t )
    deallocate( ds_01_21_t, ds_10_12_t, ds_20_22_t, ds_02_22_t )
    deallocate( ds_00_20_c, ds_00_02_c, ds_11_21_c, ds_01_11_c )
    deallocate( ds_11_12_c, ds_10_11_c, ds_01_21_c, ds_10_12_c )
    deallocate( ds_01_21_e, ds_10_12_n )

  end subroutine grid_transfer_end

  !#####################################################################

end program grid_transfer

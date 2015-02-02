program create_grid
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
  !
  ! <OVERVIEW>
  !   This program works together with program make_xgrids, runoff_regrid to 
  !   remap runoff data from a spherical grid onto any grid. This program can 
  !   generate a spherical grid netcdf file that can be used to generate exchange 
  !   grid grid_spec.nc with the destination grid file ( ocean grid file). 
  !   
  ! </OVERVIEW>

  !<DESCRIPTION>
  ! This program expects to get grid data from a netcdf file, which is specfied
  ! by the namelist variable "src_data". The name of the runoff field is specified 
  ! by the namelist variable "src_fld_name". The output file is a netcdf file with
  ! name "atmos_grid.nc". 
  !</DESCRIPTION>
  !

  implicit none

#include "netcdf.inc"

  !--- namelist interface ----------------------------------------------
  !<NAMELIST NAME="runoff_regrid_nml">
  ! <DATA NAME="src_data" TYPE="character(len=128)" >
  !  Name of input file containing runoff data.
  ! </DATA>
  ! <DATA NAME="src_fld_name" TYPE="character(len=128),dimension(2)" >
  !  Name of runoff field.
  ! </DATA>
  !</NAMELIST>


  character(len=128) :: src_fld_name = 'runoff'      ! name of the runoff field
  character(len=128) :: src_data     = 'src_data.nc' ! name of the runoff source data file

  namelist /create_grid_nml/ src_fld_name, src_data

  !--- version information ---------------------------------------------
  character(len=128)   :: version= '$ID$'
  character(len=128)   :: tagname= '$Name: tikal $'

  !--- other variables -------------------------------------------------
  real,    allocatable :: lon(:), lat(:)
  integer              :: ni, nj
  integer              :: stdout = 6 

  !--- initialize routine 
  call create_grid_init()

  !--- read source data file to get longitude and latitude -------------
  call read_src_data()

  !--- generate the simple latlon atmosphere grid
  !--- the source data should be a simple latlon grid
  call generate_grid()

contains

  !#######################################################################
  !--- initialization routine : read namelist
  subroutine create_grid_init
    integer :: unit=7, iostat
    logical :: opened

    !--- read namelist 
    do 
       inquire( unit, opened = opened )
       if( .not. opened ) exit
       unit = unit + 1
       if( unit.EQ.100 )call error_handler( 'Error: unable to locate unit number.' )
    enddo

    open (unit, file='input.nml', iostat = iostat )
    if(iostat .ne. 0) call error_handler('error in open input.nml file')

    read (unit, create_grid_nml, iostat=iostat)
    if(iostat .ne. 0) call error_handler('error in while reading create_grid_nml ')

    close(unit)

    !--- write out namelist and version information
    write( stdout,'(/a)' )'create_grid '//trim(version)//trim(tagname)
    write( stdout, create_grid_nml)

  end subroutine create_grid_init

  !#######################################################################
  !--- read the longitude and latitude from src_data
  subroutine read_src_data

    integer           :: dims(4)
    integer           :: rcode, ncid, nvar, id_runoff, i
    character(len=64) :: lon_name, lat_name, name

    rcode = nf_open(trim(src_data), NF_NOWRITE, ncid)
    call error_handler('error in opening file '//trim(src_data), rcode)

    rcode = nf_inq_nvars(ncid, nvar)
    call error_handler('error in inquring nvars', rcode)  

    rcode = nf_inq_varid(ncid, trim(src_fld_name), id_runoff )
    call error_handler('error in inquring variable id', rcode)

    rcode = nf_inq_vardimid(ncid, id_runoff, dims)
    call error_handler('error in inquring dimension id', rcode)

    rcode = nf_inq_dimlen(ncid, dims(1), ni)
    call error_handler('error in inquring dimension length ni', rcode)

    rcode = nf_inq_dimlen(ncid, dims(2), nj)
    call error_handler('error in inquring dimension length nj', rcode)

    allocate(lon(ni), lat(nj) )

    rcode = nf_inq_dimname(ncid, dims(1), lon_name)
    call error_handler('error in inquring lon_name', rcode)

    rcode = nf_inq_dimname(ncid, dims(2), lat_name)
    call error_handler('error in inquring lat_name', rcode)

    do i = 1, nvar
       rcode = nf_inq_varname(ncid,i,name)
       if(trim(name) == trim(lon_name)) then
          rcode = nf_get_var_double(ncid,i,lon)
       else if(trim(name) == trim(lat_name)) then
          rcode = nf_get_var_double(ncid,i,lat)
       endif
    enddo

  end subroutine read_src_data

  !#####################################################################
  !--- generate the spherical grid
  subroutine generate_grid

    integer           :: i, j, rcode, ncid, start(4), nwrite(4), dims(3), id_vertex 
    integer           :: id_x_t, id_y_t, id_x_vert_t, id_y_vert_t, id_grid_x_t, id_grid_y_t
    real, allocatable :: x_t(:,:), y_t(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
    real, allocatable :: lonb(:), latb(:)
    integer            ::fsize = 65536, inital = 0


    allocate(x_t(ni, nj), x_vert_t(ni, nj,4) )
    allocate(y_t(ni, nj), y_vert_t(ni, nj,4) )
    allocate(lonb(ni+1), latb(nj+1))

    do i = 1, ni
       x_t(i,:) = lon(i)
    enddo

    do j = 1, nj
       y_t(:,j) = lat(j)
    enddo

    lonb(1) = (3.0*lon(1) - lon(2))*0.5
    lonb(ni+1) = (3.0*lon(ni) - lon(ni-1))*0.5
    do i = 2, ni
       lonb(i) = 0.5*(lon(i-1) + lon(i))
    enddo

    latb(1) = (3.0*lat(1) - lat(2)) * 0.5
    latb(nj+1) = (3.0*lat(nj) - lat(nj-1)) *0.5
    do j = 2, nj
       latb(j) = 0.5*(lat(j-1) + lat(j))
    enddo

    do i = 1, ni
       x_vert_t(i,:,1) = lonb(i)
       x_vert_t(i,:,2) = lonb(i+1)
       x_vert_t(i,:,3) = lonb(i+1)
       x_vert_t(i,:,4) = lonb(i)
    enddo

    do j = 1, nj
       y_vert_t(:,j,1) = latb(j)
       y_vert_t(:,j,2) = latb(j)
       y_vert_t(:,j,3) = latb(j+1)
       y_vert_t(:,j,4) = latb(j+1)
    enddo

#ifdef use_netCDF3
    rcode = NF__CREATE('atmos_grid.nc', NF_CLOBBER, inital, fsize, ncid)
#else  
    rcode = NF__CREATE('atmos_grid.nc', IOR(NF_NETCDF4,NF_CLASSIC_MODEL), inital, fsize, ncid )
#endif
    call error_handler('error in create file atmos_grid.nc', rcode)
    rcode = nf_def_dim(ncid, 'grid_x_T', ni, dims(1))
    call error_handler('error in defining x dimension', rcode)
    rcode = nf_def_dim(ncid, 'grid_y_T', nj, dims(2))
    call error_handler('error in defining y dimension', rcode)
    rcode = nf_def_dim(ncid, 'vertex',      4, dims(3))
    call error_handler('error in defining vertex dimension', rcode)
    rcode = nf_def_var(ncid, 'grid_x_T', NF_FLOAT, 1, dims(1), id_grid_x_t)
    call error_handler('error in defining axis grid_x_T', rcode)
    rcode = nf_def_var(ncid, 'grid_y_T', NF_FLOAT, 1, dims(2), id_grid_y_t)
    call error_handler('error in defining axis grid_y_T', rcode)
    rcode = nf_def_var(ncid, 'vertex',   NF_FLOAT, 1, dims(3), id_vertex)
    call error_handler('error in defining axis vertex', rcode)
    rcode = nf_def_var(ncid, 'x_T', NF_DOUBLE, 2, dims(1:2), id_x_t)
    call error_handler('error in defining field x_T', rcode)
    rcode = nf_def_var(ncid, 'y_T', NF_DOUBLE, 2, dims(1:2), id_y_t)
    call error_handler('error in defining field y_T', rcode)
    rcode = nf_def_var(ncid, 'x_vert_T', NF_DOUBLE, 3, dims(1:3), id_x_vert_t)
    call error_handler('error in defining field x_vert_T', rcode)
    rcode = nf_def_var(ncid, 'y_vert_T', NF_DOUBLE, 3, dims(1:3), id_y_vert_t)
    call error_handler('error in defining field y_vert_T', rcode)
    rcode = nf_enddef(ncid)

    start = 1; nwrite = 1
    nwrite(1) = ni
    rcode = nf_put_vara_double(ncid, id_grid_x_t, start, nwrite, lon)
    call error_handler('error in putting the data for axis grid_x_T', rcode)
    nwrite(1) = nj
    rcode = nf_put_vara_double(ncid, id_grid_y_t, start, nwrite, lat)
    call error_handler('error in putting the data for axis grid_y_T', rcode)
    nwrite(1) = ni; nwrite(2) = nj
    rcode = nf_put_vara_double(ncid, id_x_t, start, nwrite, x_t)
    call error_handler('error in putting the data for field x_T', rcode)
    rcode = nf_put_vara_double(ncid, id_y_t, start, nwrite, y_t) 
    call error_handler('error in putting the data for field y_T', rcode)
    nwrite(3) = 4   
    rcode = nf_put_vara_double(ncid, id_x_vert_t, start, nwrite, x_vert_t)
    call error_handler('error in putting the data for field x_vert_T', rcode)
    rcode = nf_put_vara_double(ncid, id_y_vert_t, start, nwrite, y_vert_t)
    call error_handler('error in put the data for field y_vert_T', rcode)
    rcode = nf_close(ncid)

    deallocate(x_t, y_t, x_vert_t, y_vert_t, lonb, latb, lon, lat)

  end subroutine generate_grid


  !#####################################################################
  ! error handling routine.
  subroutine error_handler(mesg, status)
    character(len=*),  intent(in) :: mesg
    integer, optional, intent(in) :: status
    character(len=128) :: msg

    if(present(status)) then
       if(status == 0) return
       write(msg,*) mesg, ', status code is ', status
    else
       write(msg,*) mesg
    endif

    write(stdout,*) msg

    call ABORT()

  end subroutine error_handler

  !#######################################################################

end program create_grid

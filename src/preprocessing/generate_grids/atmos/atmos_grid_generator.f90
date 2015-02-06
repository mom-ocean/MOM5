program grid_generator
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
  !  Generate horizontal grid for atmosphere or land model. 
  !  The grid is either bgrid or spectral grid.
  !</OVERVIEW>

  use fms_mod,         only : check_nml_error, open_namelist_file, close_file, stdout
  use fms_mod,         only : write_version_number, fms_init, fms_end, file_exist
  use mpp_mod,         only : mpp_pe, mpp_root_pe, mpp_error, FATAL
  use mpp_io_mod,      only : mpp_open, MPP_OVERWR,MPP_NETCDF, MPP_SINGLE, axistype
  use atmos_grid_mod,  only : generate_atmos_grid, atmos_grid_init, atmos_grid_end
  use atmos_grid_mod,  only : write_atmos_grid_meta, write_atmos_grid_data, atmos_grid_type
  use constants_mod,   only : constants_init

  implicit none

  character(len=128), parameter :: version= '$Id: atmos_grid_generator.f90,v 10.0 2003/10/24 22:01:43 fms Exp $'
  character(len=128), parameter :: tagname='$Name: tikal $'

  type(atmos_grid_type) :: Hgrid  
  integer               :: unit   ! output_file io unit
  !--- namelist interface
  character(len=128) :: output_file = "atmos_grid.nc"
  !
  !<NAMELIST NAME="atmos_grid_generator_nml">
  !<DATA NAME="output_file" TYPE="character(len=128)">
  !  name of grid file to be created. Default value is "atmos_grid.nc".
  !</DATA>
  !</NAMELIST>
  namelist / atmos_grid_generator_nml / output_file
  integer :: io, ierr

  !---------------------------------------------------------------------------

  call fms_init()
  call constants_init()

  !--- read namelist information ---------------------------------------
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr=1
     do while (ierr /= 0)
        read  (unit, nml=atmos_grid_generator_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'atmos_grid_generator_nml')  ! also initializes nml error codes
     enddo
10   call close_file (unit)
  else
     call mpp_error(FATAL,'atmos_grid_generator: file input.nml does not exist')
  endif
  !--- write version info and namelist to logfile ----------------------
  write (stdout(), nml=atmos_grid_generator_nml)
  call write_version_number(version,tagname)

  !--- generate data ---------------------------------------------------
  call atmos_grid_init
  call generate_atmos_grid(Hgrid)

  !--- write out data-----------------------------------------------------

  if(mpp_pe() == mpp_root_pe()) then
     call mpp_open(unit,trim(output_file),MPP_OVERWR,MPP_NETCDF, &
          threading=MPP_SINGLE,  fileset=MPP_SINGLE)
     call write_atmos_grid_meta(unit, Hgrid)
     call write_atmos_grid_data(unit, Hgrid)
  endif

  !--- release the memory
  call atmos_grid_end(Hgrid)
  call fms_end()

end program grid_generator

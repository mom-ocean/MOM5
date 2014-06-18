program ocean_grid_generator
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>

  !<OVERVIEW>
  !  Generate a grid specification data file for ocean.
  !</OVERVIEW>

  !<DESCRIPTION>
  !  This program can generate horizontal grid, vertical grid or grids with topography and land/sea mask. .
  !  The namelist option grid_type control the type of grid created. when grid_type equal
  !  <PRE>
  ! 1. "hgrid": only horizontal grid will be created.
  ! 2. "vgrid": only vertical grid will be created
  ! 4. "hgrid_topog" : horizontal grid, topography and land/sea mask will be created. 
  !    Topography is obtained by remapping onto current grid from some topography source data. 
  !    land/sea mask is determined by the topography. In this case, topography does not depend 
  !    on vertical grid and no vertical grid will be created. You need to set topog_depend_on_vgrid of 
  !    topog_nml to .false. .
  ! 5. "hgrid_vgrid_topog": horizontal grid, vertical grid, topography and land/sea mask 
  !    will be created. The topography is mom4-specific topography (could be idealized or from 
  !    some source file), which depends on vertical grid. The land/sea mask is determined by 
  !    topography. In this case, you need to set topog_depend_on_vgrid of topog_nml to .true. .
  !  </PRE>
  !</DESCRIPTION>
  use fms_mod,        only : check_nml_error, open_namelist_file, close_file, stdout
  use fms_mod,        only : write_version_number, fms_init, fms_end, file_exist
  use mpp_mod,        only : mpp_pe, mpp_root_pe, mpp_error, FATAL
  use mpp_io_mod,     only : mpp_open, MPP_OVERWR, MPP_NETCDF, MPP_SINGLE, axistype
  use grids_type_mod, only : hgrid_data_type, vgrid_data_type, topog_data_type
  use hgrid_mod,      only : generate_hgrid, hgrid_init, hgrid_end, write_hgrid_global_meta 
  use hgrid_mod,      only : write_hgrid_field_meta, write_hgrid_data
  use vgrid_mod,      only : generate_vgrid, vgrid_init, vgrid_end, write_vgrid_meta, write_vgrid_data
  use topog_mod,      only : generate_topog, topog_init, topog_end, write_topog_global_meta
  use topog_mod,      only : write_topog_field_meta, write_topog_data
  use constants_mod,  only : constants_init

  implicit none

  character(len=128) :: version= '$Id: ocean_grid_generator.f90,v 13.0 2006/03/28 21:44:58 fms Exp $'
  character(len=128) :: tagname='$Name: tikal $'

  type(hgrid_data_type) :: Hgrid
  type(vgrid_data_type) :: Vgrid
  type(topog_data_type) :: Topog
  integer               :: unit, ierr, io

  !--- namelist interface
  character(len=128) :: output_file = "ocean_grid.nc"
  character(len=64)  :: grid_type   = "hgrid_vgrid_topog"
  !
  !<NAMELIST NAME="ocean_grid_generator_nml">
  !<DATA NAME="grid_type" TYPE="character(len=64)">
  !  Control the type of grid will be created. Its value can be hgrid, vgrid, hgrid_mask, 
  !  hgrid_topog_mask, hgrid_vgrid_topog_mask. Default value is hgrid_vgrid_topog. 
  !  See module description for details.
  !</DATA>
  !<DATA NAME="output_file" TYPE="character(len=128)">
  !  name of grid file to be created. Default value is "ocean_grid.nc".
  !</DATA>
  !</NAMELIST>
  namelist / ocean_grid_generator_nml / output_file, grid_type
  !---------------------------------------------------------------------------

  call fms_init()
  call constants_init()

  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr=1
     do while (ierr /= 0)
        read(unit, nml=ocean_grid_generator_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'ocean_grid_generator_nml')
     enddo
10   call close_file(unit)
  else
     call mpp_error(FATAL,'ocean_grid_generator: file input.nml does not exist')
  endif

  call write_version_number(version,tagname)
  if (mpp_pe() == mpp_root_pe()) write(stdout(), nml=ocean_grid_generator_nml)

  !--- generate data ---------------------------------------------------
  select case(trim(grid_type))
  case ('hgrid')
     call hgrid_init
     call generate_hgrid(Hgrid)
  case ('vgrid')
     call vgrid_init
     call generate_vgrid(Vgrid)
  case ('hgrid_topog')
     call hgrid_init
     call generate_hgrid(Hgrid)
     call topog_init(Topog, Hgrid)
     call generate_topog(Topog, Hgrid)
  case ('hgrid_vgrid_topog')
     call hgrid_init
     call generate_hgrid(Hgrid)
     call vgrid_init
     call generate_vgrid(Vgrid)
     call topog_init(Topog, Hgrid)
     call generate_topog(Topog, Hgrid, Vgrid)
  case default
     call mpp_error(FATAL, 'ocean_grid_generator: '//trim(grid_type)//' is not a valid option of nml "grid_type"')
  end select

  !--- write out data-----------------------------------------------------

  !--- write meta
  select case(trim(grid_type))
  case ('hgrid')
     call write_hgrid_global_meta(output_file)
     call write_hgrid_field_meta(output_file)
     call write_hgrid_data(output_file, Hgrid)
  case ('vgrid') 
     call write_vgrid_meta(output_file, Vgrid)
     call write_vgrid_data(output_file)
  case ('hgrid_topog')
     call write_hgrid_global_meta(output_file)
     call write_topog_global_meta(output_file)
     call write_hgrid_field_meta(output_file)
     call write_topog_field_meta(output_file)
     call write_hgrid_data(output_file, Hgrid)
     call write_topog_data(output_file, Topog)
  case ('hgrid_vgrid_topog')
     call write_hgrid_global_meta(output_file)
     call write_topog_global_meta(output_file)
     call write_vgrid_meta(output_file, Vgrid)
     call write_hgrid_field_meta(output_file)
     call write_topog_field_meta(output_file)
     call write_vgrid_data(output_file)
     call write_hgrid_data(output_file, Hgrid)
     call write_topog_data(output_file, Topog)
  end select

  !--- release the memory
  if(trim(grid_type) .ne. 'vgrid') call hgrid_end(Hgrid)
  if(trim(grid_type) .eq. 'vgrid' .or. (trim(grid_type) .eq. 'hgrid_vgrid_topog' )) &
       call vgrid_end(Vgrid)
  if(trim(grid_type) .eq. 'hgrid_topog' .or. trim(grid_type) .eq.'hgrid_vgrid_topog')  &
       call topog_end(Topog)

  call fms_end()

end program ocean_grid_generator

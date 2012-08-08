module grids_type_mod
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov" >Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>


  ! <OVERVIEW>
  !   Define grid types, including cell type, horizontal grid, vertical grid and topography type
  ! </OVERVIEW>

  !<DESCRIPTION>
  !    From discussions with Robert Hallberg, as a prototype for ESMF
  !</DESCRIPTION>

  use mpp_domains_mod, only: domain2d

  implicit none
  private

  public :: hgrid_data_type, vgrid_data_type, topog_data_type, cell_type

  !<PUBLICTYPE >
  type cell_type
     real, dimension(:,:), pointer   :: x => NULL()        ! geographical longitude of cell center
     real, dimension(:,:), pointer   :: y => NULL()        ! geographical latitude of cell center 
     real, dimension(:,:,:), pointer :: x_vert => NULL()   ! geographical longitude of cell vertices
     real, dimension(:,:,:), pointer :: y_vert => NULL()   ! geographical latitude of cell vertices
     real, dimension(:,:), pointer   :: area => NULL()     ! cell area
     real, dimension(:,:), pointer   :: angle => NULL()    ! Angle between i-unit and x-unit vector of cell
     real, dimension(:,:), pointer   :: ds_00_02 => NULL() ! Length of western face of cell
     real, dimension(:,:), pointer   :: ds_20_22 => NULL() ! Length of eastern face of cell
     real, dimension(:,:), pointer   :: ds_02_22 => NULL() ! Length of northern face of cell
     real, dimension(:,:), pointer   :: ds_00_20 => NULL() ! Length of southern face of cell
     real, dimension(:,:), pointer   :: ds_01_21 => NULL() ! width of cell
     real, dimension(:,:), pointer   :: ds_10_12 => NULL() ! height of cell
     real, dimension(:,:), pointer   :: ds_00_01 => NULL() ! Distance from southwest corner to western face center of cell
     real, dimension(:,:), pointer   :: ds_01_02 => NULL() ! Distance from northwest corner to western face center of cell
     real, dimension(:,:), pointer   :: ds_02_12 => NULL() ! Distance from northwest corner to northern face center of cell
     real, dimension(:,:), pointer   :: ds_12_22 => NULL() ! Distance from northeast corner to northern face center of cell
     real, dimension(:,:), pointer   :: ds_21_22 => NULL() ! Distance from northeast corner to eastern face center of cell
     real, dimension(:,:), pointer   :: ds_20_21 => NULL() ! Distance from southeast corner to eastern face center of cell
     real, dimension(:,:), pointer   :: ds_10_20 => NULL() ! Distance from southeast corner to southern face center of cell
     real, dimension(:,:), pointer   :: ds_00_10 => NULL() ! Distance from southwest corner to southern face center of cell
     real, dimension(:,:), pointer   :: ds_01_11 => NULL() ! Distance from center to western face of cell
     real, dimension(:,:), pointer   :: ds_11_12 => NULL() ! Distance from center to northern face of cell
     real, dimension(:,:), pointer   :: ds_11_21 => NULL() ! Distance from center to eastern face of cell 
     real, dimension(:,:), pointer   :: ds_10_11 => NULL() ! Distance from center to southern face of cell
  end type cell_type
  !</PUBLICTYPE>

  !<PUBLICTYPE >
  type hgrid_data_type
     type(cell_type) :: T             ! T-cell
     type(cell_type) :: E             ! E-cell
     type(cell_type) :: N             ! N-cell
     type(cell_type) :: C             ! C-cell
     logical         :: tripolar_grid ! true means tripolar grid
     logical         :: cyclic_x      ! true means cyclic in i-direction
     logical         :: cyclic_y      ! true means cyclic in j-direction
     integer         :: ni, nj        ! grid size
     type(domain2d)  :: Domain
  end type hgrid_data_type
  !</PUBLICTYPE>

  !<PUBLICTYPE >
  type vgrid_data_type
     real, dimension(:), pointer :: zt => NULL()   ! vertical level at T-cell center 
     real, dimension(:), pointer :: zb => NULL()   ! vertical level at T-cell boundary
  end type vgrid_data_type
  !</PUBLICTYPE>

  !<PUBLICTYPE >
  type topog_data_type
     real, dimension(:,:), pointer :: depth_t => NULL()    ! topographic depth of T-cell
     real, dimension(:,:), pointer :: num_levels => NULL() ! number of vertical T-cells
     real, dimension(:,:), pointer :: wet => NULL()        ! wet mask of T-cell
     real, dimension(:,:), pointer :: depth_c => NULL()    ! topographic depth of C-cell
     real, dimension(:,:), pointer :: num_levels_c => NULL() ! number of vertical C-cells
     real, dimension(:,:), pointer :: wet_c => NULL()        ! wet mask of C-cell
  end type topog_data_type
  !</PUBLICTYPE>

end module grids_type_mod

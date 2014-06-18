module river_type_mod
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
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Kirsten Findell </CONTACT> 
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT> 

  use time_manager_mod, only : time_type

  implicit none
  private

!--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river_type.F90,v 20.0 2013/12/13 23:29:45 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

!--- public interface ------------------------------------------------
  public :: river_type, Leo_Mad_trios
  integer, public, parameter :: NO_RIVER_FLAG = -9999

!--- public data type ------------------------------------------------

  type river_type
     real, dimension(:),        pointer :: lon_1d        => NULL()  ! in degree
     real, dimension(:),        pointer :: lat_1d        => NULL()  ! in degree
     real, dimension(:,:),      pointer :: lon           => NULL()  ! in radians
     real, dimension(:,:),      pointer :: lat           => NULL()  ! in radians
     real, dimension(:,:),      pointer :: reach_length  => NULL()
     real, dimension(:,:),      pointer :: landfrac      => NULL()
     logical, dimension(:,:),   pointer :: mask          => NULL()
     real, dimension(:,:),      pointer :: land_area     => NULL()
     integer, dimension(:,:),   pointer :: basinid       => NULL() 
     integer, dimension(:,:),   pointer :: tocell        => NULL()
     integer, dimension(:,:),   pointer :: travel        => NULL()
     integer, dimension(:,:),   pointer :: i_tocell      => NULL()
     integer, dimension(:,:),   pointer :: j_tocell      => NULL()
     real, dimension(:,:),      pointer :: storage       => NULL()     
     real, dimension(:,:),      pointer :: stordis       => NULL()     
     real, dimension(:,:),      pointer :: run_stor      => NULL()  ! runoff storage
     real, dimension(:,:,:),    pointer :: storage_c     => NULL()     
     real, dimension(:,:,:),    pointer :: stordis_c     => NULL()     
     real, dimension(:,:,:),    pointer :: run_stor_c    => NULL()  ! tracer runoff storage
     real, dimension(:,:),      pointer :: inflow        => NULL()
     real, dimension(:,:,:),    pointer :: inflow_c      => NULL()
     real, dimension(:,:),      pointer :: infloc        => NULL()
     real, dimension(:,:,:),    pointer :: infloc_c      => NULL()
     real, dimension(:,:),      pointer :: outflow       => NULL()
     real, dimension(:,:,:),    pointer :: outflow_c     => NULL()
     real, dimension(:,:),      pointer :: lake_outflow   => NULL()
     real, dimension(:,:,:),    pointer :: lake_outflow_c => NULL()
     real, dimension(:,:),      pointer :: disw2o        => NULL()
     real, dimension(:,:),      pointer :: disw2l        => NULL()
     real, dimension(:,:,:),    pointer :: disc2o        => NULL()
     real, dimension(:,:,:),    pointer :: disc2l        => NULL()
     real, dimension(:,:),      pointer :: melt          => NULL()
     real, dimension(:,:,:),    pointer :: removal_c     => NULL()
     real, dimension(:,:),      pointer :: outflowmean   => NULL()
     real, dimension(:,:),      pointer :: o_coef        => NULL()
     real, dimension(:,:),      pointer :: d_coef        => NULL()
     real, dimension(:,:),      pointer :: w_coef        => NULL()
     real, dimension(:,:,:),    pointer :: source_conc   => NULL()
     real, dimension(:,:,:),    pointer :: source_flux   => NULL()
     real, dimension(:,:),      pointer :: So            => NULL()
     real, dimension(:,:),      pointer :: depth         => NULL()
     real, dimension(:,:),      pointer :: width         => NULL()
     real, dimension(:,:),      pointer :: vel           => NULL()
     real, dimension(:),        pointer :: t_ref         => NULL()
     real, dimension(:),        pointer :: vf_ref        => NULL()
     real, dimension(:),        pointer :: q10           => NULL()
     real, dimension(:),        pointer :: kinv          => NULL()
     real                               :: o_exp
     real                               :: d_exp
     real                               :: w_exp
     real                               :: channel_tau
     type (time_type)                   :: Time
     integer                            :: dt_fast, dt_slow
     integer                            :: nlon, nlat, num_species, num_c
     integer                            :: num_phys
     logical                            :: do_age
  end type river_type

  type Leo_Mad_trios
     real :: on_V, on_d, on_w
  end type Leo_Mad_trios


end module river_type_mod

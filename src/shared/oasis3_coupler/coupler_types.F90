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
module coupler_types_mod  
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov">
! Stephen Griffies 
! </CONTACT>
!
!<OVERVIEW>
! This module contains simplified type declarations for the coupler
! of an ocean-only model. With mom4, this is used with ocean_solo_mod
! as the driver, and with HIM it is used with HIM_driver.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains simplified type declarations for the coupler
! of an ocean-only model. With mom4, this is used with ocean_solo_mod
! as the driver, and with HIM it is used with HIM_driver.
!</DESCRIPTION>
!

use field_manager_mod, only: fm_field_name_len, fm_string_len

implicit none ; private

type, public    :: coupler_2d_values_type
  character(len=fm_field_name_len)                      :: name = ' '
  real, pointer, dimension(:,:)                         :: values => NULL()
  logical                                               :: mean = .true.
  logical                                               :: override = .false.
  integer                                               :: id_diag = 0
  character(len=fm_string_len)                          :: long_name = ' '
  character(len=fm_string_len)                          :: units = ' '
end type coupler_2d_values_type

type, public    :: coupler_2d_field_type  !{
  character(len=fm_field_name_len)                      :: name = ' '
  integer                                               :: num_fields = 0
  type(coupler_2d_values_type), pointer, dimension(:)   :: field => NULL()
  character(len=fm_string_len)                          :: flux_type = ' '
  character(len=fm_string_len)                          :: implementation = ' '
  real, pointer, dimension(:)                           :: param => NULL()
  logical, pointer, dimension(:)                        :: flag => NULL()
  integer                                               :: atm_tr_index = 0
  character(len=fm_string_len)                          :: ice_restart_file = ' '
  character(len=fm_string_len)                          :: ocean_restart_file = ' '
  logical                                               :: use_atm_pressure
  logical                                               :: use_10m_wind_speed
  logical                                               :: pass_through_ice
  real                                                  :: mol_wt = 0.0
end type coupler_2d_field_type

type, public :: coupler_2d_bc_type 
  integer    :: num_bcs = 0
  type(coupler_2d_field_type), pointer, dimension(:)    :: bc => NULL()
end type coupler_2d_bc_type

end module coupler_types_mod 

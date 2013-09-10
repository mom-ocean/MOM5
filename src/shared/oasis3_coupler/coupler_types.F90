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

implicit none ; private

type, public :: coupler_2d_bc_type 
  integer    :: num_bcs = 0
end type coupler_2d_bc_type

end module coupler_types_mod 

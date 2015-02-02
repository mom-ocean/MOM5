!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
!-----------------------------------------------------------------------

module plev_constants_mod
implicit none
private

real, public, parameter :: GRAV   = 9.80
real, public, parameter :: RDGAS  = 287.04
real, public, parameter :: RVGAS  = 461.50
real, public, parameter :: TFREEZE  = 273.16

end module plev_constants_mod

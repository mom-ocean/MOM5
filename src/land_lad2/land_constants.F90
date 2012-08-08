module land_constants_mod

use constants_mod, only : rdgas, rvgas, wtmair

implicit none
private

! ==== public interfaces =====================================================
integer, public, parameter :: &
     NBANDS   = 2, & ! number of spectral bands for short-wave radiation calculations
     BAND_VIS = 1, & ! visible radiation (wavelenght range?)
     BAND_NIR = 2    ! near infra-red radiation (wavelenght range?)

real, public, parameter :: d622 = rdgas/rvgas
real, public, parameter :: d378 = 1.0-d622
real, public, parameter :: d608 = d378/d622

real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1

real, public, parameter :: seconds_per_year = 86400.0*365.0
real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
real, public, parameter :: mol_air = wtmair/1000.0 ! molar mass of air, kg
real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: land_constants.F90,v 17.0 2009/07/21 03:02:18 fms Exp $', &
     tagname = '$Name: siena_201207 $'

end module

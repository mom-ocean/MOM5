MODULE aerosol_params_mod

use   fms_mod,        only: write_version_number

implicit none
private

!-------------------------------------------------------------------------
!--interfaces-------------------------------------------------------------
public  aerosol_params_init

!------------------------------------------------------------------------
!----version number------------------------------------------------------
Character(len=128) :: Version = '$Id: aerosol_params.F90,v 19.0 2012/01/06 20:31:39 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!------------------------------------------------------------------------
!--namelist--------------------------------------------------------------

!------------------------------------------------------------------------
!--public data-----------------------------------------------------------

!!REAL, PARAMETER :: rho_sulf = 1770. !kg/m3 ammonium sulfate
REAL, PARAMETER, public :: rho_sulf = 1840.   ! kg/m3 sulfuric acid
REAL, PARAMETER, public :: sigma_sulf = 2.

REAL, PARAMETER, public :: rho_bc = 1000.  ! Hess et al., bams 98
REAL, PARAMETER, public :: sigma_bc = 2.

!------------------------------------------------------------------------
! note:
! dust number calculation is  hard-coded for 5 dust bins 
! with constant dV/dlnr 
! and the following bin bounds:
! [0.1,1.], [1., 2.], [2., 3.], [3.,6.], [6.,10.] micrometer
! and with densities
! rho_dust = 2500., 2650.,2650., 2650.,2650. kg/m^3
!
! The number in each bin (N_i) is (see also bin.pro)
! N = 2* M/ (rho_dust * pi) *  (1./d1^3 - 1./d2^3) / (ln(d2) - ln(d1) )
!   = M * Nfact_du
! M total mass in bin i
! d1 lower bound 
! d2 upper bound 

! the bin bounds are from atmos_tracer_driver, 
! arguments to atmos_dust_sourcesink
! ideally, the bounds would be specified 
! for example in a namelist and the Nfacts be calculated
! in an initialization routine... 
!------------------------------------------------------------------------
 REAL, PARAMETER, public ::  Nfact_du1 = 1.105e17
 REAL, PARAMETER, public ::  Nfact_du2 = 3.033e14
 REAL, PARAMETER, public ::  Nfact_du3 = 5.212e13
 REAL, PARAMETER, public ::  Nfact_du4 = 1.123e13
 REAL, PARAMETER, public ::  Nfact_du5 = 1.707e12

!------------------------------------------------------------------------
!mean dust radius for bin [D1,D2]: 1.5*(D1^-2 - D2^-2)/(D1^-3 - D2^-3) !
!------------------------------------------------------------------------
REAL, PARAMETER, public :: rbar_du1 = 7.432e-8
REAL, PARAMETER, public :: rbar_du2 = 6.429e-7
REAL, PARAMETER, public :: rbar_du3 = 1.184e-6
REAL, PARAMETER, public :: rbar_du4 = 1.929e-6
REAL, PARAMETER, public :: rbar_du5 = 3.674e-6

!------------------------------------------------------------------------
! similarly for sea salt 
! [0.1, 0.5],[0.5,1.],[1.,2.5],[2.5,5.],[5.,10.]
! rho_sea_salt = 2200 kg/m^3 for all bins
!------------------------------------------------------------------------
REAL, PARAMETER, public ::  Nfact_ss1 = 1.105e17
REAL, PARAMETER, public ::  Nfact_ss2 = 2.922e15
REAL, PARAMETER, public ::  Nfact_ss3 = 2.956e14
REAL, PARAMETER, public ::  Nfact_ss4 = 2.338e13
REAL, PARAMETER, public ::  Nfact_ss5 = 2.922e12


logical :: module_is_initialized




CONTAINS



!#########################################################################

subroutine aerosol_params_init


      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    write version number to output file.
!-------------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------------
!    declare module to be initialized.
!-------------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------------

end subroutine aerosol_params_init


!#########################################################################


END MODULE aerosol_params_mod

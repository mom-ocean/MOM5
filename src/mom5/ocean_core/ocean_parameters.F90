module ocean_parameters_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module contains some parameters used in MOM.
!</OVERVIEW>
!
!<DESCRIPTION>
! The parameter settings for numerical and/or physical schemes. 
! Also some physical constants, whose values can be modified via namelist.  
!</DESCRIPTION>
!
!
!<NAMELIST NAME="ocean_parameters_nml">
!
!  <DATA NAME="cp_ocean" UNITS="J/(kg degC)" TYPE="real">
!  Specific heat capacity J/(kg degC) for liquid seawater.
!  Values are taken from from Jackett etal (2006) for preTEOS10 and 
!  from TEOS-10 manual for TEOS10 value. The default values differ
!  from that in shared/constants since the MOM defaults are more updated.  
!  Note that there is a check inside of ocean_tempsalt.F90 to ensure that
!  cp_ocean=cp_ocean_teos10 if using the teos10 recommendations, 
!  and cp_ocean=cp_ocean_preteos10 for cases not using teos10.  
!  </DATA> 
!
!  <DATA NAME="cp_solid_runoff" UNITS="J/(kg degC)" TYPE="real">
!  Specific heat capacity J/(kg degC) for solid water runoff via calving land ice.
!  Default cp_solid_runoff = 2106.0 is consistent with that used in the 
!  GFDL land model. 
!  </DATA> 
!
!  <DATA NAME="cp_liquid_runoff" UNITS="J/(kg degC)" TYPE="real">
!  Specific heat capacity J/(kg degC) for liquid water runoff from land. 
!  Default cp_liquid_runoff = 4218.0 is consistent with that used in the 
!  GFDL land model. 
!  </DATA> 
!
!  <DATA NAME="rho0" UNITS="kg/m^3" TYPE="real">
!  Boussinesq reference density.  Default rho0=1035.0
!  corresponds to the value in Gill (page 47), where he notes 
!  that the ocean density typically deviates less than 2 per cent
!  from this value. But if using the Boussinesq approximation for 
!  other water bodies, such as the Baltic, then may wish to change
!  rho0 to a more appropriate value. 
!  </DATA> 
!
!  <DATA NAME="tfreeze" UNITS="Kelvin" TYPE="real">
!  freezing point of fresh water at standard atmos pressure.
!  Default tfreeze=273.15
!  </DATA> 
!
!  <DATA NAME="omega_earth" UNITS="radians per second" TYPE="real">
!  rotation of earth in radians per second
!  Default omega_earth= 7.2921e-5, as per equation (4.1) in Griffies (2004).
!  </DATA> 
!
!  <DATA NAME="grav" UNITS="m/s^2" TYPE="real">
!  Gravitational acceleration at earth surface. Assumed to be constant
!  throughout the ocean domain.  Default grav=9.8 corresponds to the
!  "grav" parameter from shared/constants.F90. 
!  </DATA> 
!
!  <DATA NAME="von_karman" UNITS="dimensionless" TYPE="real">
!  Von Karman constant for law of wall turbulence. 
!  Default von_karman=0.4. 
!  Note: due to answer changes on some compilers, von_karman 
!  has been removed from nml and is now set as a hard-value
!  to be consistent with earlier simulations (06mar2012).  
!  </DATA> 
!
!</NAMELIST>

  use mpp_mod,  only: stdlog, stdout, input_nml_file
  use fms_mod,  only: open_namelist_file, check_nml_error, close_file

  implicit none
  private

  public ocean_parameters_init
  public ocean_parameters_end 


  ! some numerical constants
  real, parameter, public  :: missing_value=-1.e20
  real, parameter, public  :: onehalf      = 1.0/2.0 
  real, parameter, public  :: onethird     = 1.0/3.0
  real, parameter, public  :: twothirds    = 2.0/3.0
  real, parameter, public  :: onefourth    = 1.0/4.0 
  real, parameter, public  :: onesixth     = 1.0/6.0 
  real, parameter, public  :: oneeigth     = 1.0/8.0 
  real, parameter, public  :: onesixteenth = 1.0/16.0 
  real, parameter, public  :: sec_in_day   = 86400.0 
  real, parameter, public  :: sec_in_day_r = 1.0/86400.0 
  real, parameter, public  :: sec_in_yr    = 86400.0*365.25 
  real, parameter, public  :: sec_in_yr_r  = 1.0/(86400.0*365.25)
  real, parameter, public  :: m3toSv       = 1.e-6   ! volume based Sv (10^6 m^3/s)
  real, parameter, public  :: kgtoSv       = 1.e-9   ! mass   based Sv (10^9 kg/s)

  ! parameters for choosing the time tendency calculation
  integer, parameter, public :: TWO_LEVEL   = 2
  integer, parameter, public :: THREE_LEVEL = 3

  ! parameters for vertical coordinates 
  integer, parameter, public :: GEOPOTENTIAL = 1
  integer, parameter, public :: ZSTAR        = 2
  integer, parameter, public :: ZSIGMA       = 3
  integer, parameter, public :: PRESSURE     = 4
  integer, parameter, public :: PSTAR        = 5 
  integer, parameter, public :: PSIGMA       = 6 

  integer, parameter, public :: DEPTH_BASED    = 1
  integer, parameter, public :: PRESSURE_BASED = 2

  integer, parameter, public :: QUASI_HORIZONTAL  = 1
  integer, parameter, public :: TERRAIN_FOLLOWING = 2

  ! parameters for choosing method of computing vertical cell thickness 
  integer, parameter, public :: ENERGETIC   =1
  integer, parameter, public :: FINITEVOLUME=2

  ! parameters for choosing the temperature variable
  integer, parameter, public :: CONSERVATIVE_TEMP = 1
  integer, parameter, public :: POTENTIAL_TEMP    = 2

  ! parameters for choosing the salinity variable
  integer, parameter, public :: PREFORMED_SALT = 1
  integer, parameter, public :: PRACTICAL_SALT = 2

  ! parameters for choosing the vertical mixing scheme 
  integer, parameter, public :: VERTMIX_CONST      = 1
  integer, parameter, public :: VERTMIX_PP         = 2
  integer, parameter, public :: VERTMIX_CHEN       = 3
  integer, parameter, public :: VERTMIX_GOTM       = 4
  integer, parameter, public :: VERTMIX_KPP_MOM4P0 = 5
  integer, parameter, public :: VERTMIX_KPP_MOM4P1 = 6
  integer, parameter, public :: VERTMIX_KPP_TEST   = 7

  ! parameters for choosing the horizontal layout of variables 
  integer, parameter, public :: MOM_BGRID = 1
  integer, parameter, public :: MOM_CGRID = 2

  integer, parameter, public :: TEMP_ID = 1
  integer, parameter, public :: SALT_ID = 2

  ! parameters for tracer advection 
  integer, parameter, public :: ADVECT_UPWIND            = 1
  integer, parameter, public :: ADVECT_2ND_ORDER         = 2
  integer, parameter, public :: ADVECT_4TH_ORDER         = 3
  integer, parameter, public :: ADVECT_6TH_ORDER         = 4 
  integer, parameter, public :: ADVECT_QUICKER           = 5
  integer, parameter, public :: ADVECT_QUICKMOM3         = 6
  integer, parameter, public :: ADVECT_MDFL_SUP_B        = 7
  integer, parameter, public :: ADVECT_MDPPM             = 8
  integer, parameter, public :: ADVECT_MDFL_SWEBY        = 9
  integer, parameter, public :: ADVECT_DST_LINEAR        = 10
  integer, parameter, public :: ADVECT_PSOM              = 11
  integer, parameter, public :: ADVECT_MDFL_SWEBY_TEST   = 12
  integer, parameter, public :: ADVECT_MDPPM_TEST        = 13
  integer, parameter, public :: ADVECT_DST_LINEAR_TEST   = 14
  integer, parameter, public :: ADVECT_MDMDT_TEST        = 15
 
  ! parameters for linear mometum  advection 
  integer, parameter, public :: VEL_ADVECT_UPWIND     = 1
  integer, parameter, public :: VEL_ADVECT_2ND_ORDER  = 2


  !---------------------------------------------------------
  ! the following are some physical parameters; some 
  ! can be modified through namelist settings. 

  ! acceleration due to gravity (m/s^2), assumed constant 
  ! throughout the ocean domain.  
  real, public :: grav = 9.80

  ! specific heat capacity J/(kg degC) for seawater 
  real, public            :: cp_ocean           = 3992.10322329649d0  
  real, parameter, public :: CP_OCEAN_PRETEOS10 = 3992.10322329649d0  
  real, parameter, public :: CP_OCEAN_TEOS10    = 3991.86795711963d0 

  ! specific heat capacity J/(kg degC) for calving land ice. 
  ! this value is consistent with that used in the GFDL land model. 
  real, public :: cp_solid_runoff = 2106.

  ! specific heat capacity J/(kg degC) for liquid runoff from land.
  ! this value is consistent with that used in the GFDL land model. 
  real, public :: cp_liquid_runoff = 4218.

  ! Boussinesq reference density (kg/m3) and its inverse.  Ideally, this 
  ! density should be set equal to the domain averaged density 
  ! in the model.   
  real, public :: rho0  = 1035.0
  real, public :: rho0r = 1.0/1035.0 

  ! product of rho0*cp_ocean
  ! (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C)
  real, public :: rho_cp = 1035.0 * 3992.10322329649d0 

  ! freezing point of fresh water at standard atmos pressure 
  real, public :: tfreeze  = 273.15

  ! rotation of earth in radians per second 
  ! from equation (4.1) in Griffies (2004)
  real, public :: omega_earth = 7.2921e-5 

  ! von Karman constant (dimensionless)
  ! 06mar2012: reset as real, parameter to keep
  ! consistent with earlier compilers 
  real, parameter, public :: von_karman=0.4


  character(len=128) :: version = &
     '$Id: ocean_parameters.F90,v 20.0 2013/12/14 00:10:55 fms Exp $'
  character (len=128) :: tagname = &
     '$Name: tikal $'

  namelist /ocean_parameters_nml/ cp_ocean, cp_liquid_runoff, cp_solid_runoff, &
                                  rho0, tfreeze, omega_earth, grav
  
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_parameters_init">
!
! <DESCRIPTION>
! Initialize the parameter module, passing cp_ocean back to ocean_model.F90
! and setting all other parameters that will be used throughout the model
! simulation. 
!
! Note: we do not enable check_nml_error, since the default settings 
! are generally those recommended for simulations.
!
! </DESCRIPTION>
!
subroutine ocean_parameters_init()

    integer :: stdoutunit, stdlogunit
    integer :: ioun, io_status

    stdoutunit = stdout()
    stdlogunit = stdlog()

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_parameters_nml, iostat=io_status)
!    ierr = check_nml_error(io_status,'ocean_parameters_nml')
#else
    ioun = open_namelist_file()
    read  (ioun, ocean_parameters_nml,iostat=io_status)
!    ierr = check_nml_error(io_status, 'ocean_parameters_nml')
    call close_file(ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_parameters_nml)
    write (stdlogunit, ocean_parameters_nml)
    write( stdlogunit,'(/a/)') trim(version)

    rho0r  = 1.0/rho0
    rho_cp = rho0*cp_ocean 

    return

  end subroutine ocean_parameters_init
! </SUBROUTINE> NAME="ocean_parameters_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_parameters_end">
!
! <DESCRIPTION>
! Summarize the basic physical parameters used in the simulation. 
! </DESCRIPTION>
!
subroutine ocean_parameters_end()

    integer :: stdoutunit
    stdoutunit=stdout()
 
    write (stdoutunit,'(a)')  ' '
    write (stdoutunit,'(2x,a)') 'Some of the basic physical parameters used in the simulation'
    write (stdoutunit,'(2x,a,e20.10)')  'Reference density for the Boussinesq approximation (kg/m^3) = ',rho0
    write (stdoutunit,'(2x,a,e20.10)')  'Gravitational acceleration (m/s^2)                          = ',grav
    write (stdoutunit,'(2x,a,e20.10)')  'Rotational rate of the earth (radians/sec)                  = ',omega_earth
    write (stdoutunit,'(2x,a,e20.10)')  'Specific heat capacity of seawater      J/(kg*degC)         = ',cp_ocean
    write (stdoutunit,'(2x,a,e20.10)')  'Specific heat capacity of liquid runoff J/(kg*degC)         = ',cp_liquid_runoff
    write (stdoutunit,'(2x,a,e20.10)')  'Specific heat capacity of solid  runoff J/(kg*degC)         = ',cp_liquid_runoff
    write (stdoutunit,'(2x,a,e20.10)')  'Von Karman dimensionless constant for turbulent mixing      = ',von_karman
    write (stdoutunit,'(a)')  ' '


  end subroutine ocean_parameters_end
! </SUBROUTINE> NAME="ocean_parameters_end"


end module ocean_parameters_mod



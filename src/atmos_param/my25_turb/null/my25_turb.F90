  MODULE MY25_TURB_MOD

!=======================================================================
!   MELLOR-YAMADA LEVEL 2.5 TURBULENCE CLOSURE SCHEME - GFDL VERSION   !
!=======================================================================

 use Fms_Mod, ONLY: ERROR_MESG, FATAL, mpp_pe, mpp_root_pe, write_version_number

!---------------------------------------------------------------------
 implicit none
 private
!---------------------------------------------------------------------

 public :: MY25_TURB, MY25_TURB_INIT, MY25_TURB_END, TKE_SURF, get_tke
 public :: my25_turb_restart

!---------------------------------------------------------------------

 character(len=128) :: version = '$Id: my25_turb.F90,v 17.0 2009/07/21 02:55:45 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'

 contains

!#######################################################################

 SUBROUTINE MY25_TURB( is, js, delt, fracland, phalf, pfull, theta, &   
                       um,   vm,       zhalf, zfull, z0,    &
                       el0,      el,    akm,   akh,   &
                       mask, kbot,     ustar, bstar, h    )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       delt     -  Time step in seconds
!       fracland -  Fractional amount of land beneath a grid box
!       phalf    -  Pressure at half levels
!       pfull    -  Pressure at full levels
!       theta    -  Potential temperature
!       um, vm   -  Wind components
!       zhalf    -  Height at half levels
!       zfull    -  Height at full levels
!       z0       -  Roughness length
!       mask     -  OPTIONAL; floating point mask (0. or 1.) designating
!                   where data is present
!       kbot     -  OPTIONAL;lowest model level index (integer);
!                    at levels > kbot, mask = 0.
!       ustar    -  OPTIONAL:friction velocity (m/sec)
!       bstar    -  OPTIONAL:buoyancy scale (m/sec**2)
!---------------------------------------------------------------------
  integer, intent(in)                   :: is, js
  real,    intent(in)                   :: delt 
  real,    intent(in), dimension(:,:)   :: fracland, z0
  real,    intent(in), dimension(:,:,:) :: phalf, pfull, zhalf, zfull
  real,    intent(in), dimension(:,:,:) :: um, vm, theta

  integer, intent(in), OPTIONAL, dimension(:,:)   :: kbot
  real,    intent(in), OPTIONAL, dimension(:,:,:) :: mask
  real,    intent(in), OPTIONAL, dimension(:,:)   :: ustar, bstar

!---------------------------------------------------------------------
! Arguments (Intent out)
!       el0  -  characteristic length scale
!       el   -  master length scale
!       akm  -  mixing coefficient for momentum
!       akh  -  mixing coefficient for heat and moisture
!         h  -  OPTIONAL, diagnosed depth of planetary boundary 
!                         layer (m)
!---------------------------------------------------------------------
  real, intent(out), dimension(:,:)   :: el0
  real, intent(out), dimension(:,:,:) :: akm, akh, el
  real, intent(out), OPTIONAL, dimension(:,:) :: h

call error_mesg('MY25_TURB', &
 'This module is not supported as part of the public release', FATAL)

!====================================================================
  end SUBROUTINE MY25_TURB 

!#######################################################################
subroutine get_tke(is, ie, js, je, tke_out)
integer, intent(in) :: is, ie, js, je
real, intent(out), dimension(:,:,:) :: tke_out

call error_mesg('get_tke', &
 'This module is not supported as part of the public release', FATAL)

end subroutine get_tke
!#######################################################################

  SUBROUTINE MY25_TURB_INIT( ix, jx, kx )

!=======================================================================
! ***** INITIALIZE MELLOR-YAMADA
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     ix, jx  - Horizontal dimensions for global storage arrays
!     kx      - Number of vertical levels in model
!---------------------------------------------------------------------
 integer, intent(in) :: ix, jx, kx

 call error_mesg('MY25_TURB_INIT', &
 'This module is not supported as part of the public release', FATAL)
 
  end SUBROUTINE MY25_TURB_INIT

!#######################################################################

  SUBROUTINE MY25_TURB_END
!=======================================================================

 call error_mesg('MY25_TURB_END', &
 'This module is not supported as part of the public release', FATAL)
 
!=====================================================================

  end SUBROUTINE MY25_TURB_END

!#######################################################################

 SUBROUTINE TKE_SURF ( is, js, u_star, kbot )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       u_star -  surface friction velocity (m/s)
!       kbot   -  OPTIONAL;lowest model level index (integer);
!                 at levels > Kbot, Mask = 0.
!---------------------------------------------------------------------
  integer, intent(in) :: is, js
  real, intent(in), dimension(:,:)   :: u_star

  integer, intent(in), OPTIONAL, dimension(:,:) :: kbot

 call error_mesg('TKE_SURF', &
 'This module is not supported as part of the public release', FATAL)

!=======================================================================
 end SUBROUTINE TKE_SURF 

!#######################################################################
! <SUBROUTINE NAME="my25_turb_restart">
!
! <DESCRIPTION>
!   dummy interface
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine my25_turb_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

end subroutine my25_turb_restart
! </SUBROUTINE> NAME="my25_turb_restart"

!#######################################################################
  end MODULE MY25_TURB_MOD

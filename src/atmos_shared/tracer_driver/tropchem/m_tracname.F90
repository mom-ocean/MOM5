
      module m_tracname_mod
!-----------------------------------------------------------
!       ... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use mo_grid_mod,   only : pcnst
      use chem_mods_mod, only : grpcnt

      implicit none

character(len=128), parameter :: version     = '$Id: m_tracname.F90,v 19.0 2012/01/06 20:32:14 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      save

      character(len=8) :: tracnam(pcnst)          ! species names
      character(len=8) :: natsnam(max(1,grpcnt))  ! names of non-advected trace species

      end module m_tracname_mod

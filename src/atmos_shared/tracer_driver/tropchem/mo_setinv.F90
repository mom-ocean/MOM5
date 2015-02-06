      module MO_SETINV_MOD

implicit none
character(len=128), parameter :: version     = '$Id: mo_setinv.F90,v 19.0 2012/01/06 20:34:06 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine setinv( invariants, tfld, h2ovmr, pmid, inv_data, &
                         do_interactive_h2o, plonl )
!-----------------------------------------------------------------
!        ... Set the invariant densities (molecules/cm**3)
!-----------------------------------------------------------------
      
      implicit none

!-----------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)  :: plonl
      real,    intent(in)  :: tfld(:,:)           ! temperature
      real,    intent(in)  :: h2ovmr(:,:)         ! water vapor vmr
      real,    intent(in)  :: pmid(:,:)           ! pressure
      real,    intent(in)  :: inv_data(:,:,:)     ! invariants
      logical, intent(in)  :: do_interactive_h2o  ! include h2o sources/sinks?
      real,    intent(out) :: invariants(:,:,:)   ! invariant array
!-----------------------------------------------------------------
!        .. Local variables
!-----------------------------------------------------------------
      real, parameter ::  boltz = 1.38044e-16      ! erg/K

      integer :: k,j,n
      integer :: plev
      
      plev = SIZE(tfld,2)

!-----------------------------------------------------------------
!        NOTE: Invariants are in cgs density units.
!              The pmid array is in pascals and must be
!              mutiplied by 10. to yield dynes/cm**2.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... Set M, N2, O2, and H2O densities
!-----------------------------------------------------------------
      do k = 1,plev
         invariants(:,k,1) = 10. * pmid(:,k) / (boltz*tfld(:,k))
         invariants(:,k,2) = .79 * invariants(:,k,1)
         invariants(:,k,3) = .21 * invariants(:,k,1)
         if (do_interactive_h2o) then
            n=3
         else
            n=4
            invariants(:,k,n) = h2ovmr(:,k) * invariants(:,k,1)
         end if
         do j = 1, size(invariants,3)-n
            invariants(:,k,j+n) = inv_data(:,k,j) * invariants(:,k,1)
         enddo
      end do

      end subroutine SETINV

      end module MO_SETINV_MOD

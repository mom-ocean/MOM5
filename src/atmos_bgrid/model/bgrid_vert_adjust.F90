
module bgrid_vert_adjust_mod

!-----------------------------------------------------------------------
! This module computes the following terms:
!    * surface pressure tendency
!    * vertical part of thermodynamic term (aka, omega-alpha term)
!    * vertical mass flux (sigma-dot)
!-----------------------------------------------------------------------

use bgrid_vert_mod, only:  vert_grid_type
use        fms_mod, only:  write_version_number,  &
                           mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                           MPP_CLOCK_SYNC, CLOCK_MODULE

implicit none
private

!-----------------------------------------------------------------------
!---- public interfaces -----

public :: vert_adjust, vert_adjust_init

!-----------------------------------------------------------------------

character(len=128) :: version='$Id: bgrid_vert_adjust.F90,v 11.0 2004/09/28 19:07:04 fms Exp $'
character(len=128) :: tagname='$Name: tikal $'

integer :: id_clock
logical :: initialized = .false.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_adjust ( Vgrid, res, div, wta, wtb, mask,   &
                         omgalf, etadot, psdt)

!-----------------------------------------------------------------------
! INPUT
!   Vgrid   = vertical grid constants
!   res     = reciprocal of eta at the surface (=1 for sigma)
!   div     = mass divergence (Pa/s)
!   wta,wtb = weights for obtaining full level values
!   mask    = grid box mask for eta coordinate mountains
! INPUT/OUTPUT
!   omgalf  = thermodynamic (omega-alpha) term (Pa/s)
!              (input: horiz part; output: horiz+vert part)
!   etadot  = vertical mass flux (a summation) (Pa/s)
! OUTPUT
!   psdt    = surface pressure tendency (Pa/s)
!-----------------------------------------------------------------------

type(vert_grid_type), intent(in)      :: Vgrid
real, intent(in),    dimension(:,:)   :: res
real, intent(in),    dimension(:,:,:) :: div, wta, wtb, mask
real, intent(inout), dimension(:,:,:) :: omgalf
real, intent(inout), dimension(:,:,:) :: etadot
real, intent(out),   dimension(:,:)   :: psdt

! sdiv = vertical integral of div
real, dimension(size(div,1),size(div,2),size(div,3)) ::  sdiv
!-----------------------------------------------------------------------
   if (.not.initialized) call vert_adjust_init
   call mpp_clock_begin (id_clock)

!----- vertical adjustments ----

   call ps_tendency (div, psdt, sdiv)
   call vert_omgalf (Vgrid, res, wta, wtb, sdiv, mask, omgalf)

!----- compute "eta-dot", the vertical velocity -----

   call vert_velocity (res, sdiv, mask, Vgrid%eta, etadot)

   call mpp_clock_end (id_clock)
!-----------------------------------------------------------------------

end subroutine vert_adjust

!#######################################################################

subroutine vert_adjust_init

   call write_version_number(version, tagname)

   id_clock = mpp_clock_id ('BGRID: vert_adjust', flags=MPP_CLOCK_SYNC, &
                                                  grain=CLOCK_MODULE)
   initialized=.true.

end subroutine vert_adjust_init

!#######################################################################

subroutine ps_tendency (div, psdt, sdiv)

!-----------------------------------------------------------------------
!
!          Routine for calculating the vertical integral
!          of divergence and surface pressure tendency
!
!-----------------------------------------------------------------------

   real, intent(in)  :: div(:,:,:)
   real, intent(out) :: psdt(:,:), sdiv(:,:,:)

   integer  k, kdim
!-----------------------------------------------------------------------

      kdim = size(div,3)

!---------integrate divergence to get surface pressure tendency---------

      sdiv(:,:,1) = div(:,:,1)
   do k = 2, kdim
      sdiv(:,:,k) = sdiv(:,:,k-1) + div(:,:,k)
   enddo
      psdt(:,:) = -sdiv(:,:,kdim)

!-----------------------------------------------------------------------

end subroutine ps_tendency

!#######################################################################

subroutine vert_velocity (res, sdiv, mask, eta, etadot)

!-----------------------------------------------------------------------
!
!    Routine for calculating the vertical velocity ("eta-dot")
!
!-----------------------------------------------------------------------

   real, intent(in)  :: res(:,:), sdiv(:,:,:),  &
                        mask(:,:,:), eta(:)
   real, intent(inout) :: etadot(:,:,:)

   real, dimension(size(res,1),size(res,2)) :: pret
   integer  k, kdim

!---- computation of etadot (add onto previous value) ----

   kdim = size(sdiv,3)

   pret(:,:) = -sdiv(:,:,kdim)*res(:,:)

   do k = 2, kdim
      etadot(:,:,k) = etadot(:,:,k)  &
                      - (pret(:,:)*eta(k)+sdiv(:,:,k-1))*mask(:,:,k)
   enddo

!-----------------------------------------------------------------------

end subroutine vert_velocity

!#######################################################################

subroutine vert_omgalf (Vgrid, res, wta, wtb, sdiv, mask, omgalf)

!-----------------------------------------------------------------------
!
!   Routine for calculating the vertical part of the omega-alpha term
!
!-----------------------------------------------------------------------

   type(vert_grid_type), intent(in)       :: Vgrid
   real, intent(in),    dimension(:,:)    :: res
   real, intent(in),    dimension(:,:,:)  :: wta, wtb, sdiv, mask
   real, intent(inout), dimension(:,:,:)  :: omgalf

   real, dimension(size(res,1),size(res,2)) :: pret
   integer  k, kdim

!-----------------------------------------------------------------------
!-----kinetic energy generation terms in the thermodynamic equation-----

   kdim = size(sdiv,3)

      omgalf(:,:,1) = omgalf(:,:,1) - wtb(:,:,1)*sdiv(:,:,1)*mask(:,:,1)

!!!do k=2,kdim-1
   do k=2,kdim
      omgalf(:,:,k) = omgalf(:,:,k) - mask(:,:,k)* &
                    (wta(:,:,k)*sdiv(:,:,k-1)+wtb(:,:,k)*sdiv(:,:,k))
   enddo

!  --- may need to add/modify for eta coordinate? ---
!  if (kdim > 1) then
!     pret(:,:) = -sdiv(:,:,kdim)*res(:,:)
!     omgalf(:,:,kdim) = omgalf(:,:,kdim) + mask(:,:,kdim)* &
!                   (wtb(:,:,k)*pret(:,:)-wta(:,:,k)*sdiv(:,:,kdim-1))
!  endif

!-----------------------------------------------------------------------

end subroutine vert_omgalf

!#######################################################################

end module bgrid_vert_adjust_mod



module bgrid_masks_mod

!-----------------------------------------------------------------------
!
!    allocates memory and initializes masks and vertical indexing
!    associated with the step mountain vertical coordinate.
!
!-----------------------------------------------------------------------

use bgrid_horiz_mod, only:  horiz_grid_type
use  bgrid_halo_mod, only:  update_halo, NOPOLE, UWND
use  bgrid_vert_mod, only:  vert_grid_type
use         fms_mod, only:  write_version_number, stdlog, stdout
use         mpp_mod, only:  mpp_pe, mpp_root_pe, mpp_max

implicit none
private

!-----------------------------------------------------------------------
!    ---- public interfaces ----

   public :: grid_masks_init, grid_masks_end

!-----------------------------------------------------------------------
!    ---- public data types ----

   public ::      mask_type
   public :: grid_mask_type

   type mask_type
      real,    pointer, dimension(:,:,:) :: mask =>NULL()
      integer, pointer, dimension(:,:)   :: kbot =>NULL()
      integer                            :: kbotmin
   end type mask_type

!  mask  = step-mountain topography mask (0. or 1.)
!          mask = 0. for grid boxes that form the step-mountain
!          mask = 1. for grid boxes above ground
!  kbot  = lowest model level above ground
!  kbotmin = smallest value of kbot across all processors
!
!  NOTE:  For the sigma coordinate, mask = 1.0 everywhere, and
!         kbot = number of vertical model levels

   type grid_mask_type
      type(mask_type) :: Tmp, Vel
      logical :: sigma
   end type grid_mask_type

!  Tmp = grid masking values for the temperature/mass grid
!  Vel = grid masking values for the velocity grid
!  sigma = logical flag that specifies whether the vertical coordinate
!          is the eta/step-mountain coordinate or sigma coordinate
!
!-----------------------------------------------------------------------
!  private interfaces and data

   private :: compute_mass_mask, compute_vel_mask, compute_lowest_level

   logical :: sigma  ! local variable

!  version number info

   character(len=128) :: version='$Id: bgrid_masks.F90,v 19.0 2012/01/06 19:55:15 fms Exp $'
   character(len=128) :: tagname='$Name: tikal $'
   logical :: do_log = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine grid_masks_init (Hgrid, Vgrid, res, Masks)

!-----------------------------------------------------------------------
! input arguments
!     Hgrid = horizontal grid constants
!     Vgrid = vertical grid constants
!     res   = 1/eta(surface); i.e., reciprical of eta at the surface
! output
!     Masks = grid box masks for the eta coordinate
!-----------------------------------------------------------------------

   type (horiz_grid_type), intent(inout) :: Hgrid
   type  (vert_grid_type), intent(in)    :: Vgrid
   real,                   intent(in)    :: res(Hgrid%ilb:,Hgrid%jlb:)
   type  (grid_mask_type), intent(inout) :: Masks

!-----------------------------------------------------------------------

   logical :: sigma   ! flag for sigma(T) or eta(F)
   integer :: k, outunit, logunit
   real    :: maxres
   real    :: aeta(Vgrid%nlev)

   if (do_log) then
      call write_version_number(version, tagname)
      do_log = .false.
   endif
   logunit = stdlog()
   outunit = stdout()

   ! allocate global storage

 allocate ( Masks%Tmp%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, Vgrid%nlev), &
            Masks%Vel%mask (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub, Vgrid%nlev), &
            Masks%Tmp%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub), &
            Masks%Vel%kbot (Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub)  )

   ! determine if this is a sigma or eta
   ! vertical coordinate from the value of res

      maxres = maxval(res)
      call mpp_max (maxres)

      if (maxres > 1.0001) then
         Masks % sigma=.false.
         sigma=.false.
         if (mpp_pe() == mpp_root_pe()) then
            if (Vgrid%hybrid) then
               write (logunit,100) 'eta/hybrid'
            else
               write (logunit,100) 'eta'
            endif
             ! print gaudy banner to standard output for untested eta option
               write (outunit,'(4(a/))') &
                      '*******************************************************', &
                      'WARNING: The eta coordinate option has not been tested.', &
                      '         Proceed with caution.',                          &
                      '*******************************************************'
         endif
      else
         Masks % sigma=.true.
         sigma=.true.
         if (mpp_pe() == mpp_root_pe()) then
            if (Vgrid%hybrid) then
               write (logunit,100) 'sigma/hybrid'
            else
               write (logunit,100) 'sigma'
            endif
         endif
      endif

 100  format (/,'B-grid dynamical core has been initialized with the ',a,' vertical coordinate.')

!--------------topography masks ----------------------------------------

    ! average eta at model levels
      do k = 1, Vgrid%nlev
        aeta(k) = (Vgrid%eta(k)+Vgrid%eta(k+1))*0.5
      enddo

      Masks % Tmp % mask = compute_mass_mask (res, aeta)
      Masks % Vel % mask = compute_vel_mask  (res, aeta)
      call update_halo (Hgrid, UWND, Masks%Vel%mask, flags=NOPOLE)
!!!!! call update_halo (Hgrid, UWND, Masks%Vel%mask)   ! sets mask=0 at poles

!------------- compute the lowest model level --------------------------

      Masks % Tmp % kbot = compute_lowest_level (Masks % Tmp % mask)
      Masks % Vel % kbot = compute_lowest_level (Masks % Vel % mask)

!     ----- global values -----

      Masks % Tmp % kbotmin = minval(Masks % Tmp % kbot)
      Masks % Vel % kbotmin = minval(Masks % Vel % kbot)

!-----------------------------------------------------------------------

end subroutine grid_masks_init

!#######################################################################

subroutine grid_masks_end ( Masks )

type  (grid_mask_type), intent(inout) :: Masks

! release memory
  deallocate ( Masks%Tmp%mask, Masks%Vel%mask, &
               Masks%Tmp%kbot, Masks%Vel%kbot  )
  Masks%Tmp%kbotmin = 0  ! set unrealistic value
  Masks%Vel%kbotmin = 0

end subroutine grid_masks_end

!#######################################################################

   function compute_mass_mask (res, aeta) result (mask)

   real, intent(in) :: res(:,:), aeta(:)
   real, dimension(size(res,1),size(res,2),size(aeta,1)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=1,size(res,2); do i=1,size(res,1)
         do k=1,size(aeta(:))
            if (aeta(k) > (1.0/res(i,j))) mask(i,j,k) = 0.0
         enddo; enddo; enddo
      endif

   end function compute_mass_mask

!#######################################################################

   function compute_vel_mask (res, aeta) result (mask)

   real, intent(in) :: res(:,:), aeta(:)
   real, dimension(size(res,1),size(res,2),size(aeta,1)) :: mask
   integer  i, j, k

      mask = 1.0

      if (.not.sigma) then
         do j=2,size(res,2); do i=2,size(res,1)
         do k=1,size(aeta(:))
            if (aeta(k) > (1.0/res(i,j))) then
                mask(i-1,j-1,k) = 0.0
                mask(i  ,j-1,k) = 0.0
                mask(i-1,j  ,k) = 0.0
                mask(i  ,j  ,k) = 0.0
            endif
         enddo; enddo; enddo
      endif

   end function compute_vel_mask

!#######################################################################

   function compute_lowest_level (mask) result (kbot)

   real, intent(in) :: mask(:,:,:)
   integer, dimension(size(mask,1),size(mask,2)) :: kbot
   integer   i, j, k, kdim

      kdim = size(mask,3)
      kbot = kdim

      if (.not.sigma) then
         do j=1,size(mask,2); do i=1,size(mask,1)
         do k=1,kdim
            if (mask(i,j,k) < 0.50) then
               kbot(i,j)=k-1; exit
            endif
         enddo; enddo; enddo
       ! must not be zero
         kbot = max(kbot,1)
      endif

   end function compute_lowest_level

!#######################################################################

end module bgrid_masks_mod


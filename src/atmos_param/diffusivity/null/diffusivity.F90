
module diffusivity_mod

!=======================================================================
!
!                          DIFFUSIVITY MODULE
!
!     Routines for computing atmospheric diffusivities in the 
!       planetary boundary layer and in the free atmosphere
!
!=======================================================================


use constants_mod, only : grav, vonkarm, cp_air, rdgas, rvgas

use       fms_mod, only:  error_mesg, FATAL, file_exist,   &
                          check_nml_error, open_namelist_file,      &
                          mpp_pe, mpp_root_pe, close_file, &
                          write_version_number, stdlog

use monin_obukhov_mod, only : mo_diff

implicit none
private

! public interfaces
!=======================================================================

 public diffusivity, pbl_depth, molecular_diff, &
        diffusivity_init, diffusivity_end

!
!=======================================================================


!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: diffusivity.F90,v 10.0 2003/10/24 22:00:28 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
logical            :: module_is_initialized   = .false.

!=======================================================================

contains

!=======================================================================

subroutine diffusivity_init


!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized  = .true.

!---------------------------------------------------------------------

      call error_mesg('diffusivity_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diffusivity_init

!=======================================================================

subroutine diffusivity_end

      module_is_initialized  = .false.

!---------------------------------------------------------------------

      call error_mesg('diffusivity_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diffusivity_end

!=======================================================================

subroutine diffusivity(t, q, u, v, p_full, p_half, z_full, z_half,  &
                       u_star, b_star, h, k_m, k_t, kbot)

real,    intent(in),           dimension(:,:,:) :: t, q, u, v
real,    intent(in),           dimension(:,:,:) :: p_full, p_half
real,    intent(in),           dimension(:,:,:) :: z_full, z_half
real,    intent(in),           dimension(:,:)   :: u_star, b_star
real,    intent(inout),        dimension(:,:,:) :: k_m, k_t
real,    intent(out),          dimension(:,:)   :: h
integer, intent(in), optional, dimension(:,:)   :: kbot

! input:
  
!        t     : real, dimension(:,:,:) -- (:,:,pressure), third index running
!                          from top of atmosphere to bottom
!                 temperature (K)
!
!        q     : real, dimension(:,:,:)
!                 water vapor specific humidity (nondimensional)
!
!        u     : real, dimension(:,:)
!                 zonal wind (m/s)
!
!        v     : real, dimension(:,:,:) 
!                 meridional wind (m/s) 
!
!        z_full  : real, dimension(:,:,: 
!                 height of full levels (m)
!                 1 = top of atmosphere; size(p_half,3) = surface
!                 size(z_full,3) = size(t,3)
!
!        z_half  : real, dimension(:,:,:)
!                 height of  half levels (m)
!                 size(z_half,3) = size(t,3) +1
!              z_half(:,:,size(z_half,3)) must be height of surface!
!                                  (if you are not using eta-model)
!
!        u_star: real, dimension(:,:)
!                friction velocity (m/s)
!
!        b_star: real, dimension(:,:)
!                buoyancy scale (m/s**2)

!   (u_star and b_star can be obtained by calling 
!     mo_drag in monin_obukhov_mod)

! output:

!        h     : real, dimension(:,:,) 
!                 depth of planetary boundary layer (m)
!
!        k_m   : real, dimension(:,:,:)
!                diffusivity for momentum (m**2/s)
!
!                defined at half-levels
!                size(k_m,3) should be at least as large as size(t,3)
!                only the returned values at 
!                      levels 2 to size(t,3) are meaningful
!                other values will be returned as zero
!
!        k_t   : real, dimension(:,:,:)
!                diffusivity for temperature and scalars (m**2/s)
!
!---------------------------------------------------------------------

      call error_mesg('diffusivity', &
      'This module is not supported as part of the public release', FATAL)

end subroutine diffusivity

!=======================================================================

subroutine pbl_depth(t, u, v, z, u_star, b_star, h, kbot)


real,   intent(in) ,           dimension(:,:,:) :: t, u, v, z
real,   intent(in) ,           dimension(:,:)   :: u_star,b_star
real,   intent(out),           dimension(:,:)   :: h
integer,intent(in) , optional, dimension(:,:)   :: kbot

!---------------------------------------------------------------------

      call error_mesg('pbl_depth', &
      'This module is not supported as part of the public release', FATAL)

end subroutine pbl_depth

!=======================================================================

subroutine molecular_diff ( temp, press, k_m, k_t)

real, intent(in),    dimension (:,:,:)  ::  temp, press
real, intent(inout), dimension (:,:,:)  ::  k_m, k_t    

!---------------------------------------------------------------------

      call error_mesg('molecular_diff', &
      'This module is not supported as part of the public release', FATAL)

end subroutine molecular_diff 

!=======================================================================

end module diffusivity_mod


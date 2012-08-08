! The include file fms_platform.h will handle the conversion of POINTER to _ALLOCATABLE arrays
! for derived type members. The conversion affects performance only and should not change 
! any numeric result. It is limited to member arrays that are used within MOM only
! and to arrays that are never associated (=>) with another array.

! Fortran 90 requires array members of derived type to have the POINTER attribute. 
! However, most Fortran 95 compilers now also support _ALLOCATABLE array components 
! (a Fortran 2003 feature). This avoids the aliasing problem afflicting pointers. 
! Some compilers may require an additional switch (e.g. -fno-alias) to fully exploit
! the performance benefit of the conversion.
!
! Macros used from fms_platform.h:
! __ALLOCATABLE maps to either POINTER  or _ALLOCATABLE
! _NULL        maps to either =>NULL() or "nothing"
#include <fms_platform.h>

module wave_types_mod

#include <ocean_memory.h>

  type, public :: ocean_wave_type

#ifndef MOM_STATIC_ARRAYS
    real, dimension(:,:,:),_ALLOCATABLE:: xmom, ymom  _NULL !x & y-component of local wave momentum (scaled by 1/rho/g -> m s)
    real, dimension(:,:),  _ALLOCATABLE:: wave_k      _NULL !wave number (1/m)
    real, dimension(:,:),  _ALLOCATABLE:: height      _NULL !significant wave height (m)
    real, dimension(:,:),  _ALLOCATABLE:: wave_p      _NULL !peak frequency at time level taup1
#else
    real, dimension(isd:ied,jsd:jed,2) :: xmom, ymom	    !x & y-component of local wave momentum
    real, dimension(isd:ied,jsd:jed)   :: wave_k	    !wave number (1/m)
    real, dimension(isd:ied,jsd:jed)   :: height	    !significant wave height (m)
    real, dimension(isd:ied,jsd:jed)   :: wave_p	    !peak frequency at time level taup1
#endif

  end type ocean_wave_type

end module wave_types_mod

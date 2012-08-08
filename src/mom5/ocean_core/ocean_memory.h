! This file contains some commonly used indices 
! needed for setting up the various domains in MOM.  
!
! The MOM_STATIC_ARRAYS preprocessor option sets these 
! indices according to values specified at compile time, 
! whereas dynamic memory determines the indices at runtime.  
!
! Dynimical allocation is more conveneient as it allows one to 
! make changes to the model grid configuration without recompiling.
! However, on some platforms, use of static allocated arrays,
! whose size is set at runtime, provides for a more efficient
! executable.  Tests should be run to determine if significant
! differences exist on your specific computer and model 
! configuration.  
!

#ifdef MOM_STATIC_ARRAYS

  ! computational domain indices (excludes halos)
  integer, parameter :: isc = 1
  integer, parameter :: iec = NI_LOCAL_
  integer, parameter :: jsc = 1
  integer, parameter :: jec = NJ_LOCAL_

  ! data domain indices (includes halo size of 1) 
  integer, parameter :: isd = 0
  integer, parameter :: ied = NI_LOCAL_ + 1 
  integer, parameter :: jsd = 0
  integer, parameter :: jed = NJ_LOCAL_ + 1

  ! global domain indices
  integer, parameter :: isg = 1
  integer, parameter :: ieg = NI_
  integer, parameter :: jsg = 1
  integer, parameter :: jeg = NJ_
  integer, parameter :: ni = NI_
  integer, parameter :: nj = NJ_
  integer, parameter :: nk  = NK_

#else

  integer :: isd
  integer :: ied
  integer :: jsd
  integer :: jed
  integer :: isc
  integer :: iec
  integer :: jsc
  integer :: jec
  integer :: isg
  integer :: ieg
  integer :: jsg
  integer :: jeg
  integer :: ni
  integer :: nj
  integer :: nk

#endif

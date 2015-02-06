module leapfrog_mod

use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

!===================================================================================================
implicit none
private
!===================================================================================================

interface leapfrog
  module procedure leapfrog_2d_complex, leapfrog_3d_complex
end interface

interface leapfrog_2level_A
  module procedure leapfrog_2level_A_2d_complex, &
                   leapfrog_2level_A_3d_complex, &
                   leapfrog_2level_A_3d_real
end interface

interface leapfrog_2level_B
  module procedure leapfrog_2level_B_2d_complex, &
                   leapfrog_2level_B_3d_complex, &
                   leapfrog_2level_B_3d_real
end interface

character(len=128), parameter :: version = '$Id leapfrog.f90 $'
character(len=128), parameter :: tagname = '$Name: tikal $'

public :: leapfrog, leapfrog_2level_A, leapfrog_2level_B

logical :: entry_to_logfile_done = .false.

contains

!================================================================================

subroutine leapfrog_2level_A_3d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff)

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

return
end subroutine leapfrog_2level_A_3d_complex

!================================================================================

subroutine leapfrog_2level_B_3d_complex (a, current, future, robert_coeff)

complex, intent(inout), dimension(:,:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future)

return
end subroutine leapfrog_2level_B_3d_complex

!================================================================================

subroutine leapfrog_2level_A_3d_real (a, dt_a, previous, current, future, delta_t, robert_coeff)

real, intent(inout), dimension(:,:,:,:) :: a
real, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

return
end subroutine leapfrog_2level_A_3d_real

!================================================================================

subroutine leapfrog_2level_B_3d_real (a, current, future, robert_coeff)

real, intent(inout), dimension(:,:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future)

return
end subroutine leapfrog_2level_B_3d_real


!================================================================================

subroutine leapfrog_2level_A_2d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff)

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(in),    dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff

complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_2level_A_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff)
a = a_3d(:,:,1,:)

end subroutine leapfrog_2level_A_2d_complex
!================================================================================

subroutine leapfrog_2level_B_2d_complex (a, current, future, robert_coeff)

complex, intent(inout), dimension(:,:,:) :: a
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff
complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
call leapfrog_2level_B_3d_complex (a_3d, current, future, robert_coeff)
a = a_3d(:,:,1,:)

end subroutine leapfrog_2level_B_2d_complex
!================================================================================

subroutine leapfrog_3d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current) + a(:,:,:,future ))
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*a(:,:,:,future)
endif

return
end subroutine leapfrog_3d_complex

!================================================================================

subroutine leapfrog_2d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff)

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(in),    dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff

complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff)
a = a_3d(:,:,1,:)

return
end subroutine leapfrog_2d_complex

!================================================================================

end module leapfrog_mod

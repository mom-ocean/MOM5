      subroutine upper(string,length)

!***********************************************************************
!
!     upper.f - change lower case letter to upper case letter          *
!                                                                      *
!     George Lai Tue Jun 28 16:37:00 1994                              *
!                                                                      *
!***********************************************************************

      implicit         none

!      character string(length)
!      character*20 string
!      character, dimension(length) :: string
!      character (len=*), intent(inout) ::  string
!      character (len=*) ::  string
!      character (len=1), intent(inout) ::  string(20)
!ok      character (len=20), intent(inout) ::  string
      character (len=*), intent(inout) ::  string
      character char1
      integer,   intent(in)    ::  length
      integer i
      integer a, z, dist
      a = ichar('a')
      z = ichar('z')
      dist = ichar('A') - a

      do i = 1,length
        char1=string(i:i)
        if (ichar(char1) .ge. a .and.       &
            ichar(char1) .le. z) then
          string(i:i) = char(ichar(char1)+dist)
        endif
      end do

      return
      end


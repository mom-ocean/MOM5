! ============================================================================
! module numerics: a collection of useful general-purpose routines
! ============================================================================
#define __ERROR__(message) \
call my_error(mod_name,message,FATAL,__FILE__,__LINE__)
#define __NOTE__(message) \
call my_error(mod_name,message,NOTE,__FILE__,__LINE__)
#define __ASSERT__(x, message) \
if(.NOT.(x))call my_error(mod_name,message,FATAL,__FILE__,__LINE__)




module numerics_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Christopher Milly
! </CONTACT> 

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Elena Shevliakova
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Sergey Malyshev
! </REVIEWER>

! <OVERVIEW>
!   A collection of useful general-purpose numerical routines, including
!   a bisect function and linear interpolation routines.
! </OVERVIEW>

  ! thing we use repeatedly inside the module
  use fms_mod, only: error_mesg, FATAL, NOTE, write_version_number
  use constants_mod, only : PI
implicit none
private

! ==== public interfaces =====================================================
public :: bisect    ! finds a position of point in array of bounds
public :: lin_int   ! linear interpolation
public :: numerics_init
public :: is_latlon, expand_cell
! ==== end of public interfaces ==============================================

! <INTERFACE NAME="lin_int">

!   <OVERVIEW>
!     Linear interpolation.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Linear interpolation. 
!   </DESCRIPTION>

interface lin_int
   module procedure lin_int0
   module procedure lin_int1
   module procedure lin_int2
   module procedure lin_int1m
   module procedure lin_int2m
end interface
! </INTERFACE>

logical :: module_is_initialized =.FALSE.
! module constants
character(len=*), parameter :: mod_name = "Numerics_mod"
character(len=128) :: version = '$Id: numerics.F90,v 15.0 2007/08/14 04:00:08 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

contains


! ============================================================================

! <SUBROUTINE NAME="numerics_init">

!   <OVERVIEW>
!     Initializes the numerics module.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Initializes the numerics module.
!   </DESCRIPTION>

subroutine numerics_init()

  module_is_initialized =.TRUE. 
  call write_version_number(version, tagname)

end subroutine numerics_init
! </SUBROUTINE>


! <FUNCTION NAME="bisect">

!   <OVERVIEW>
!     Finds a position of point in array of bounds.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Finds a position of point in array of bounds. Returns i, so that x is
!     between xx(i) and xx(i+1).
!   </DESCRIPTION>

!   <TEMPLATE>
!     value=bisect( xx, x1, periodic )
!   </TEMPLATE>

!   <PUBLICROUTINE>

function bisect(xx, x1, periodic)

  real, intent(in)              :: xx(:)     ! array of boundaries
  real, intent(in)              :: x1        ! point to locate
  logical, intent(in), optional :: periodic  ! if present and true, the data
                                             ! domain is assumed to be periodic

!   </PUBLICROUTINE>

  ! ---- result type ---------------------------------------------------
  integer bisect

  ! ---- local vars ----------------------------------------------------
  real    :: x              ! duplicate of input value
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: ascending      ! if true, the coordinates are in ascending order

  n = size(xx(:))
  x = x1

  ! bring the point inside bounds of the periodic domain
  if (present(periodic)) then
     if(periodic) then
       __ASSERT__(xx(n)-xx(1)/=0,"periodic bisect: period equal to zero")
       x = modulo(x-min(xx(1),xx(n)),abs(xx(n)-xx(1)))+min(xx(1),xx(n))
    endif
  endif

!   <ERROR MSG="Periodic bisect: period equal to zero" STATUS="FATAL">
!     Period is equal to zero, [ (xx(n)-xx(1)) / = 0 ]
!   </ERROR>

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     ascending = xx(n) > xx(1)
     do while (high-low > 1)
        mid = (low+high)/2
        if (ascending.eqv.xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else
     bisect = -1
  endif

end function bisect
! </FUNCTION>


! <SUBROUTINE NAME="lin_int0" INTERFACE="lin_int">

!   <OVERVIEW>
!     Linearly interpolates 1-D data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Linearly interpolates 1-D data.
!   </DESCRIPTION>

!   <PUBLICROUTINE INTERFACE="lin_int">
subroutine lin_int0(data, xx, x, res)

  real, intent(in) :: data(:)    ! data to interpolate
  real, intent(in) :: xx(:)      ! coord. corresponding to data
  real, intent(in) :: x          ! coord to interpolate to
  real, intent(inout) :: res     ! result of interpolation
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx(:)),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(i1)*f1+data(i2)*f2

end subroutine lin_int0
! </SUBROUTINE>


! <SUBROUTINE NAME="lin_int1" INTERFACE="lin_int">

!   <OVERVIEW>
!     Linearly interpolates 1-D data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Linearly interpolates 1-D data.
!     
!   </DESCRIPTION>

!   <PUBLICROUTINE INTERFACE="lin_int">
subroutine lin_int1(data, xx, x, res)

  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coord. corresponding to data
  real, intent(in) :: x            ! coord to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1

!   <ERROR MSG="Coordinate is out of range" STATUS="FATAL">
!     Coordinate corresponding to data is out of range.
!   </ERROR>

  __ASSERT__(i1>0.AND.i1<size(xx(:)),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(:,i1)*f1+data(:,i2)*f2

end subroutine lin_int1
! </SUBROUTINE>



! <SUBROUTINE NAME="lin_int2" INTERFACE="lin_int">

!   <OVERVIEW>
!     Interpolates prescribed over time.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Interpolates prescribed over time.
!   </DESCRIPTION>

!   <PUBLICROUTINE INTERFACE="lin_int">
subroutine lin_int2(data, tt, t, res)

  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt(:)),"Coordinate is out of range")

!   <ERROR MSG="Coordinate is out of range" STATUS="FATAL">
!     Time moments corresponding to data points is out of range.
!   </ERROR>

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  res = data(:,:,i1)*f1+data(:,:,i2)*f2

end subroutine lin_int2
! </SUBROUTINE>


! <SUBROUTINE NAME="lin_int1m" INTERFACE="lin_int">

!   <OVERVIEW>
!     Linearly interpolates 1-D data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Linearly interpolates 1-D data.
!   </DESCRIPTION>

!   <PUBLICROUTINE INTERFACE="lin_int">
subroutine lin_int1m(data, xx, x, res, mask)

  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coord. corresponding to data
  real, intent(in) :: x            ! coord to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation
  logical, intent(in) :: mask(:)   ! valid data mask
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx(:)),"Coordinate is out of range")

!   <ERROR MSG="Coordinate is out of range" STATUS="FATAL">
!     Coordinate corresponding to data is out of range.
!   </ERROR>

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! finally, update the result
  where (mask) 
     res = data(:,i1)*f1+data(:,i2)*f2
  endwhere

end subroutine lin_int1m
! </SUBROUTINE>

! <SUBROUTINE NAME="lin_int2m" INTERFACE="lin_int">

!   <OVERVIEW>
!     Interpolates prescribed over time.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Interpolates prescribed over time.
!   </DESCRIPTION>

!   <PUBLICROUTINE INTERFACE="lin_int">
subroutine lin_int2m(data, tt, t, res, mask)

  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result
  logical, intent(in) :: mask(:,:) ! interpolation mask
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt(:)),"Coordinate is out of range")

!   <ERROR MSG="Coordinate is out of range" STATUS="FATAL">
!     Time moments corresponding to data points is out of range.
!   </ERROR>

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  where (mask) 
     res = data(:,:,i1)*f1+data(:,:,i2)*f2
  endwhere
end subroutine lin_int2m
! </SUBROUTINE>

!######################################################################
!   Routines added to deal with cubed sphere geometry
!######################################################################

function is_latlon ( lon, lat )
real, intent(in) :: lon(:,:), lat(:,:)

! Determines if the latitude/longitude values form a
! regular lat/lon grid (or something else).
!

logical :: is_latlon
integer :: i, j 

  is_latlon = .true.

  do j = 2, size(lon,2)
  do i = 1, size(lon,1)
     if (lon(i,j) .ne. lon(i,1)) then
        is_latlon = .false. 
        return  
     endif   
  enddo
  enddo

  do j = 1, size(lat,2)
  do i = 2, size(lat,1)
     if (lat(i,j) .ne. lat(1,j)) then
        is_latlon = .false. 
        return  
     endif   
  enddo
  enddo

end function is_latlon

!=================================================

 subroutine expand_cell(lon, lat, fac)
! Util for land model (for BW)
!
!        4----3
!        |  . |
!        1----2
!
      real, intent(inout):: lon(2,2), lat(2,2)
      real, intent(in):: fac    ! expansion factor: outside: > 1
                                ! fac = 1: qq1 returns q1
                                ! fac = 0: qq1_4 retruns the center position
! Local
      real qq1(3), qq2(3), qq3(3), qq4(3)
      real p1(3), p2(3), p3(3), p4(3)
      real ec(3)
      real dd, d1, d2, d3, d4
      integer k

! Transform to (x,y,z)
      call latlon2xyz((/lon(1,1),lat(1,1)/), p1)
      call latlon2xyz((/lon(2,1),lat(2,1)/), p2)
      call latlon2xyz((/lon(1,2),lat(1,2)/), p3)
      call latlon2xyz((/lon(2,2),lat(2,2)/), p4)

! Get center position:
      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd   ! cell center position
      enddo

! Perform the "extrapolation" to outside if fac > 1
      do k=1,3
         qq1(k) = ec(k) + fac*(p1(k)-ec(k))
         qq2(k) = ec(k) + fac*(p2(k)-ec(k))
         qq3(k) = ec(k) + fac*(p3(k)-ec(k))
         qq4(k) = ec(k) + fac*(p4(k)-ec(k))
      enddo

!--------------------------------------------------------
! Force the points to be on the sphere via normalization
!--------------------------------------------------------
      d1 = sqrt( qq1(1)**2 + qq1(2)**2 + qq1(3)**2 )
      d2 = sqrt( qq2(1)**2 + qq2(2)**2 + qq2(3)**2 )
      d3 = sqrt( qq3(1)**2 + qq3(2)**2 + qq3(3)**2 )
      d4 = sqrt( qq4(1)**2 + qq4(2)**2 + qq4(3)**2 )
      do k=1,3
         qq1(k) = qq1(k) / d1
         qq2(k) = qq2(k) / d2
         qq3(k) = qq3(k) / d3
         qq4(k) = qq4(k) / d4
      enddo

!----------------------------------------
! Transform back to lat-lon coordinates:
!----------------------------------------

      call cart_to_latlon1(qq1, lon(1,1), lat(1,1))
      call cart_to_latlon1(qq2, lon(2,1), lat(2,1))
      call cart_to_latlon1(qq3, lon(1,2), lat(1,2))
      call cart_to_latlon1(qq4, lon(2,2), lat(2,2))

 end subroutine expand_cell

!=================================================

 subroutine latlon2xyz(p, e)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real, intent(in) :: p(2)
 real, intent(out):: e(3)

 integer n
 real :: q(2)
 real :: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz

!=================================================

 subroutine cart_to_latlon1(q, xs, ys)
  real, intent(inout):: q(3)
  real, intent(inout):: xs, ys
! local
  real, parameter:: esl=1.e-10
  real :: p(3)
  real :: pi, dist, lat, lon
  integer :: k

  pi = 4.*atan(1.)

     do k=1,3
        p(k) = q(k)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0.
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = 2.*pi + lon
     lat = asin(p(3))
     
     xs = lon
     ys = lat
! q Normalized:
     do k=1,3
        q(k) = p(k)
     enddo

 end  subroutine cart_to_latlon1

!######################################################################

! <SUBROUTINE NAME="my_error">

!   <OVERVIEW>
!     Reports error, including file name and line.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Reports error, including file name and line.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call my_error(mod_name, message, mode, file, line)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine my_error(mod_name, message, mode, file, line)

  character(len=*), intent(in) :: mod_name
  character(len=*), intent(in) :: message
  integer,          intent(in) :: mode
  character(len=*), intent(in) :: file
  integer,          intent(in) :: line
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  character(len=512) :: mesg
  write(mesg,'("File ",a," Line ",i4.4," :: ",a)')&
       file, line, trim(message)
  call error_mesg(mod_name, mesg, mode)

end subroutine
! </SUBROUTINE>

end module numerics_mod

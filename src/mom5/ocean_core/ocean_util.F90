module ocean_util_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Tim.Leslie@gmail.com"> Tim Leslie 
!</CONTACT>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! This module contains many routines of use for MOM. 
!</OVERVIEW>
!
!<DESCRIPTION>
! A utility module for MOM. 
!</DESCRIPTION>
!

#define COMP isc:iec,jsc:jec

use platform_mod,        only: i8_kind
use constants_mod,       only: epsln
use diag_manager_mod,    only: send_data, register_diag_field 
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_sum, NON_BITWISE_EXACT_SUM
use mpp_mod,             only: stdout, stdlog, FATAL
use mpp_mod,             only: mpp_error, mpp_pe, mpp_root_pe, mpp_min, mpp_max, mpp_chksum
use time_manager_mod,    only: time_type, get_date

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_density_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_time_type

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()

public iplot
public matrix
public ocean_util_init
public invtri
public invtri_bmf
public write_timestamp
public write_note
public write_warning
public write_line
public write_chksum_header
public write_chksum_2d
public write_chksum_2d_int
public write_chksum_3d
public write_chksum_3d_int
public check_restart
public write_summary

public diagnose_2d
public diagnose_2d_u
public diagnose_2d_en
public diagnose_2d_int
public diagnose_2d_comp
public diagnose_3d
public diagnose_3d_u
public diagnose_3d_en
public diagnose_3d_int
public diagnose_3d_comp
public diagnose_sum
public register_2d_t_field
public register_3d_t_field
public register_3d_xte_field
public register_3d_ytn_field
public register_3d_ztb_field



contains


!#######################################################################
! <SUBROUTINE NAME="ocean_util_init">
!
! <DESCRIPTION>
! Initialize MOM utilities.
! </DESCRIPTION>
!
subroutine ocean_util_init (Domain, Grid)
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_grid_type),   intent(in)         :: Grid
  integer :: stdlogunit

  stdlogunit=stdlog()

  write( stdlogunit,'(/a/)') trim(version)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grid%nk
#endif

  Dom => Domain

end subroutine ocean_util_init
! </SUBROUTINE> NAME="ocean_util_init">


!#######################################################################
! <SUBROUTINE NAME="invtri">
!
! <DESCRIPTION>
! Solve the vertical diffusion equation implicitly using the
! method of inverting a tridiagonal matrix as described in
! Numerical Recipes in Fortran, The art of Scientific Computing,
! Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
! pages 42,43.
!
! This routine assumes that the variables are defined at grid points,
! and the top and bottom b.c. are flux conditions.
!
! inputs:
!
! z         = right hand side terms
!
! nk        = number of vertical levels
!
! topbc     = top boundary condition
!
! botbc     = bottom boundary condition
!
! dcb       = vertical mixing coeff at base of cell
!
! tdt       = timestep over which do implicit update
!
! kmz       = level indicator
!
! mask      = land/sea mask
!
! outputs:
!
! z         = returned solution
!
! </DESCRIPTION>
!
subroutine invtri (z, topbc, botbc, dcb, tdt, kmz, mask, dh, dhw, aidif, nk)

  integer, intent(in)                             :: nk
  real, intent(inout), dimension(isd:,jsd:,:)     :: z
  real, intent(in), dimension(isd:,jsd:,:)        :: dcb
  real, intent(in), dimension(isd:,jsd:,:)        :: dh
  real, intent(in), dimension(isd:,jsd:,:)        :: mask
  real, intent(in), dimension(isd:,jsd:,0:)       :: dhw
  real, intent(in), dimension(isd:,jsd:)          :: topbc
  real, intent(in), dimension(isd:,jsd:)          :: botbc
  real, intent(in)                                :: tdt
  integer, intent(in), dimension(isd:ied,jsd:jed) :: kmz
  real, intent(in)                                :: aidif

  real, dimension(isd:ied,0:nk) :: a, b, c, e, f
  real, dimension(isd:ied) :: bet

  integer :: i, j, k, km1, kp1
  real :: eps, factu, factl, tdt_aidif

  eps = 1.e-30
  tdt_aidif = tdt*aidif
  do j=jsc,jec
    do k=1,nk
      km1   = max(1,k-1)
      kp1   = min(k+1,nk)
      do i=isc,iec
        factu  = tdt_aidif/(dhw(i,j,k-1)*dh(i,j,k))
        factl  = tdt_aidif/(dhw(i,j,k)*dh(i,j,k))
        a(i,k) = -dcb(i,j,km1)*factu*mask(i,j,k)
        c(i,k) = -dcb(i,j,k)*factl*mask(i,j,kp1)
        f(i,k) = z(i,j,k)*mask(i,j,k) 
        b(i,k) = 1.0 - a(i,k) - c(i,k)
      enddo
    enddo

    do i=isc,iec
      a(i,1)  = 0.0
      c(i,nk) = 0.0
      b(i,1)  = 1.0 - a(i,1) - c(i,1)
      b(i,nk) = 1.0 - a(i,nk) - c(i,nk)

      ! top and bottom b.c.
      f(i,1)  = z(i,j,1) + topbc(i,j)*tdt_aidif*mask(i,j,1)/dh(i,j,1)
      k = max(2,kmz(i,j))
      f(i,k)  = z(i,j,k) - botbc(i,j)*tdt_aidif*mask(i,j,k)/dh(i,j,k)
    enddo

    ! decomposition and forward substitution
    do i=isc,iec
      bet(i)   = mask(i,j,1)/(b(i,1) + eps)
      z(i,j,1) = f(i,1)*bet(i)
    enddo
    do k=2,nk
      do i=isc,iec
        e(i,k)   = c(i,k-1)*bet(i)
        bet(i)   = mask(i,j,k)/(b(i,k) - a(i,k)*e(i,k) + eps)
        z(i,j,k) = (f(i,k) - a(i,k)*z(i,j,k-1))*bet(i)
      enddo
    enddo

    ! back substitution
    do k=nk-1,1,-1
      do i=isc,iec
        z(i,j,k) = z(i,j,k) - e(i,k+1)*z(i,j,k+1)
      enddo
    enddo
  enddo

end subroutine invtri
! </SUBROUTINE> NAME="invtri">


!#######################################################################
! <SUBROUTINE NAME="invtri_bmf">
!
! <DESCRIPTION>
! Solve the vertical friction equation implicitly using the
! method of inverting a tridiagonal matrix as described in
! Numerical Recipes in Fortran, The art of Scientific Computing,
! Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
! pages 42,43.
!
! This routine assumes that the variables are defined at grid points,
! and the top b.c. is a flux condition.  The bottom b.c. is assumed
! to be a bottom drag which is implemented implicitly, thus allowing
! for large values of the bottom drag coefficient. 
!
! NOTE: This routine is generally only called when doing the bmf 
! implicitly in time.  The original invtri is used for explicit
! bmf.  
!
! inputs:
!
! z         = right hand side terms
!
! nk        = number of vertical levels
!
! topbc     = top boundary condition
!
! botbc     = time explicit bottom boundary condition (zero in this routine, since bmf is implicit) 
!
! gamma     = botttom drag factor scaling the u(taup1) contribution to bottom drag 
! 
! dcb       = vertical mixing coeff at base of cell
!
! tdt       = timestep over which do implicit update
!
! kmz       = level indicator
!
! mask      = land/sea mask
!
! outputs:
!
! z         = returned solution
!
! </DESCRIPTION>
!
subroutine invtri_bmf (z, topbc, gamma, dcb, tdt, kmz, mask, dh, dhw, aidif, nk)

  integer, intent(in)                             :: nk
  real, intent(inout), dimension(isd:,jsd:,:)     :: z
  real, intent(in), dimension(isd:,jsd:)          :: gamma
  real, intent(in), dimension(isd:,jsd:,:)        :: dcb
  real, intent(in), dimension(isd:,jsd:,:)        :: dh
  real, intent(in), dimension(isd:,jsd:,:)        :: mask
  real, intent(in), dimension(isd:,jsd:,0:)       :: dhw
  real, intent(in), dimension(isd:,jsd:)          :: topbc
  real, intent(in)                                :: tdt
  integer, intent(in), dimension(isd:ied,jsd:jed) :: kmz
  real, intent(in)                                :: aidif

  real, dimension(isd:ied,0:nk) :: a, b, c, e, f
  real, dimension(isd:ied) :: bet

  integer :: i, j, k, km1, kp1
  real :: eps, factu, factl, tdt_aidif

  eps = 1.e-30
  tdt_aidif = tdt*aidif
  do j=jsc,jec
    do k=1,nk
      km1   = max(1,k-1)
      kp1   = min(k+1,nk)
      do i=isc,iec
        factu  = tdt_aidif/(dhw(i,j,k-1)*dh(i,j,k))
        factl  = tdt_aidif/(dhw(i,j,k)*dh(i,j,k))
        a(i,k) = -dcb(i,j,km1)*factu*mask(i,j,k)
        c(i,k) = -dcb(i,j,k)*factl*mask(i,j,kp1)
        f(i,k) = z(i,j,k)*mask(i,j,k) 
        b(i,k) = 1.0 - a(i,k) - c(i,k)
      enddo
    enddo

    do i=isc,iec

      a(i,1)  = 0.0
      c(i,nk) = 0.0
      b(i,1)  = 1.0 - a(i,1)  - c(i,1)
      b(i,nk) = 1.0 - a(i,nk) - c(i,nk)

      ! top and bottom b.c.
      f(i,1)  = z(i,j,1)  + topbc(i,j)*tdt_aidif*mask(i,j,1)/dh(i,j,1)
      k = max(2,kmz(i,j))
      f(i,k)  = z(i,j,k)
      c(i,k)  = 0.0
      b(i,k)  = 1.0 - a(i,k) - c(i,k) + gamma(i,j)*tdt_aidif*mask(i,j,k)/dh(i,j,k)
    enddo

    ! decomposition and forward substitution
    do i=isc,iec
      bet(i)   = mask(i,j,1)/(b(i,1) + eps)
      z(i,j,1) = f(i,1)*bet(i)
    enddo
    do k=2,nk
      do i=isc,iec
        e(i,k)   = c(i,k-1)*bet(i)
        bet(i)   = mask(i,j,k)/(b(i,k) - a(i,k)*e(i,k) + eps)
        z(i,j,k) = (f(i,k) - a(i,k)*z(i,j,k-1))*bet(i)
      enddo
    enddo

    ! back substitution
    do k=nk-1,1,-1
      do i=isc,iec
        z(i,j,k) = z(i,j,k) - e(i,k+1)*z(i,j,k+1)
      enddo
    enddo
  enddo

end subroutine invtri_bmf
! </SUBROUTINE> NAME="invtri_bmf">


!#######################################################################
! <SUBROUTINE NAME="iplot">
!
! <DESCRIPTION>
!
!  map integer array "iarray" into characters for printing with
!  format (a1) to provide a contour map of the integer field.
!  note: max number of unique characters = 80
!
! inputs:
!
! iarray = integer array to be plotted
!
! is     = starting index along inner dimension of "iarray"
!
! ie     = ending index along inner dimension of "iarray"
!
! js     = starting index along outer dimension of "iarray"
!
! je     = ending index along outer dimension of "iarray"
!
! output: prints contour map of "iarray"
!
! </DESCRIPTION>
!
subroutine iplot (iarray, is, ie, js, je, ni, nj)
  
  integer, intent(in) :: is, ie, js, je, ni, nj
  integer, intent(in), dimension(ni,nj) ::  iarray
  character*80 levels
  character*80 lev1
  integer :: i, j, il, jl, l, incr, ii, jj, inc, last, jinc, maxint, minint
  save levels
  integer :: stdoutunit 
  stdoutunit=stdout() 

  write (stdoutunit,*) ' '

  ! set character markers
  lev1(1:51) = '.+*ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuv'
  levels = lev1(1:51)//'wxyz0123456789-=!@#$%<>[]{}()'

  ! find range of integers
  maxint = iarray(is,js)
  minint = iarray(is,js)
  do j=js,je
    do i=is,ie
      maxint = max(maxint,iarray(i,j))
      minint = min(minint,iarray(i,j))
    enddo
  enddo

  ! show mapping of integers into characters
  write (stdoutunit,*) ' '
  write (stdoutunit,*) ' "iplot" mapping of integers to characters is as follows:'
  inc  = 3
  last = min(minint+80-1,maxint) 
  do i=minint,last,inc
    ii = i-minint+1
    if (i+inc <= last) then
      jinc = inc
    else
      jinc = last-i+1
    endif
    write (stdoutunit,'(6(1x,i6,a,a,3x))')  (j+minint-1,' is printed as ',levels(j:j),j=ii,ii+jinc-1)
  enddo
  write (stdoutunit,*) ' '

  if (maxint - minint + 1 > 80) then
    write (stdoutunit,*) ' Note: there are ',maxint-minint+1,' integers in the field'
    write (stdoutunit,*) '       "iplot" cannot uniquely assign more than 80 characters for plotting symbols.'
    write (stdoutunit,*) '       therefore integers are represented by cyclicly reusing the list of plotting symbols'
    write (stdoutunit,*) ' '
  endif

  ! print character representation of integers
  inc=124
  il = ie-is+1
  jl = je-js+1
  do l=0,il,inc
    incr = min(inc,il-l)
    write (stdoutunit,8800) (l+i,i=1,incr,4)
    do jj=1,jl
      j = jl+1-jj
      write (stdoutunit,8900) j, (levels(mod(iarray(l+i+is-1,j+js-1)-minint+1-1,80)+1:&
             mod(iarray(l+i+is-1,j+js-1)-minint+1-1,80)+1),i=1,incr) 
    enddo
  enddo
  8800  format (/, 2x, 31i4)
  8900  format (1x,i3,1x, 124a1)

end subroutine iplot
! </SUBROUTINE> NAME="iplot">


!#######################################################################
! <SUBROUTINE NAME="matrix">
!
! <DESCRIPTION>
! matrix is a general two-dimensional array printing routine,
! input:
!
! array = the array to be printed
!
! istrt = the 1st element of the 1st dimension to be printed
!
! im    = the last element of the 1st dimension to be printed
!
! jstrt = the 1st element of the 2nd dimension to be printed
!
! jm    = the last element of the 2nd dimension to be printed
!         the 2nd dimension is printed in reverse order if both
!         jstrt and jm are negative
!
! scale = a scaling factor by which array is divided before
!         printing.  (if this is zero, no scaling is done.)
!         if scale=0, 10 columns are printed across in e format
!         if scale>0, 20 columns are printed across in f format
!
! output: print "array" as a matrix
!
! </DESCRIPTION>
!
subroutine matrix (array, istrt, im, jstrt, jm, scale)

  integer :: i, l, istrt, im, jstrt, jm, is, ie, js, je, jinc, unit
  real, dimension(istrt:im,abs(jstrt):abs(jm)) ::  array
  real :: scale, scaler

  unit = 6
  if (jstrt*jm < 0) then
    write (unit,999)  jstrt, jm
    call mpp_error(FATAL,'==>Error in ocean_util_mod (matrix): jstrt*jm < 0 found in matrix')
  endif

  ! allow for inversion of 2nd dimension
  if (jm < 0) then
    js   = -jm
    je   = -jstrt
    jinc = -1
  else
    js   = jstrt
    je   = jm
    jinc = 1
  endif

  if (scale == 0.0) then

    do is=istrt,im,10
      ie = min(is + 9,im)
      write (unit,9001) (i, i=is,ie)
      do l=js,je,jinc
        write (unit,9002) l, (array(i,l),i=is,ie)
      enddo
      write (unit,'(/)')
    enddo
  else
    scaler = 1.0/scale
    do is=istrt,im,20
      ie = min(is + 19,im)
      write (unit,9003) (i, i=is,ie)
      do l=js,je,jinc
        write (unit,9004) l, (array(i,l)*scaler,i=is,ie)
      enddo
      write (unit,'(/)')
    enddo
  endif
  999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
  9001  format(10i13)
  9002  format(i3,10(1pe13.5))
  9003  format(3x,20i6)
  9004  format(1x,i3,1x,20f6.2)

end subroutine matrix
! </SUBROUTINE> NAME="matrix">


!#######################################################################
! <SUBROUTINE NAME="write_timestamp">
!
! <DESCRIPTION>
! Write the time stamp.
! </DESCRIPTION>
!
subroutine write_timestamp(Time)

  type(time_type), intent(in) :: Time
  integer :: yr, mon, day, hr, min, sec
  integer :: stdoutunit 

  stdoutunit=stdout() 

  call get_date(Time, yr, mon, day, hr, min, sec)

  write(stdoutunit,'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)' ) &
       'yyyy/mm/dd hh:mm:ss = ', yr, '/',mon,'/',day, hr, ':',min,':', sec

end subroutine write_timestamp
! </SUBROUTINE> NAME="write_timestamp">

!#######################################################################
! <SUBROUTINE NAME="write_chksum_header">
!
! <DESCRIPTION>
! Write the checksum header.
! </DESCRIPTION>
!
subroutine write_chksum_header(filename, desc, model_time)

  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: desc
  type(time_type), intent(in) :: model_time
  integer :: stdoutunit 

  stdoutunit=stdout() 
  write(stdoutunit,'(a)') ' ' 
  write(stdoutunit,'(a)') '[CHKSUM BLOCK] '//filename//': '//desc
  call write_timestamp(model_time)

end subroutine
! </SUBROUTINE> NAME="write_chksum_header">


!#######################################################################
! <SUBROUTINE NAME="write_note">
!
! <DESCRIPTION>
! Write a note.
! </DESCRIPTION>
!
subroutine write_note(filename, msg)    

  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: msg
  integer :: stdoutunit

  stdoutunit=stdout()
  write(stdoutunit, '(/a)') "[Note] "//filename//": "// msg

end subroutine write_note
! </SUBROUTINE> NAME="write_note">


!#######################################################################
! <SUBROUTINE NAME="write_warning">
!
! <DESCRIPTION>
! Write a warning.
! </DESCRIPTION>
!
subroutine write_warning(filename, msg)    

  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: msg
  integer :: stdoutunit

  stdoutunit=stdout()
  write(stdoutunit, '(/a)') "[Warning] "//filename//": "// msg

end subroutine write_warning
! </SUBROUTINE> NAME="write_warning">


!#######################################################################
! <SUBROUTINE NAME="write_line">
!
! <DESCRIPTION>
! Write a message.
! </DESCRIPTION>
!
subroutine write_line(msg)

  character(len=*), intent(in) :: msg
  integer :: stdoutunit

  stdoutunit=stdout()
  write(stdoutunit, '(a)')  "          "// msg

end subroutine write_line
! </SUBROUTINE> NAME="write_line">


!#######################################################################
! <SUBROUTINE NAME="write_chksum_3d">
!
! <DESCRIPTION>
! Write a 3d checksum.
! </DESCRIPTION>
!
subroutine write_chksum_3d(name, data, chksum)
  character(len=*), intent(in) :: name
  real, dimension(isc:iec,jsc:jec,nk), intent(in) :: data
  integer(i8_kind), optional, intent(inout) :: chksum
  integer(i8_kind)                          :: chk_sum

  integer :: stdoutunit
  character(len=40) :: c
  c = '[chksum] '//name
  stdoutunit=stdout()
  chk_sum = mpp_chksum(data(isc:iec,jsc:jec,:nk))
  write(stdoutunit, '(a40,i30)') c, chk_sum
  if (present(chksum)) chksum = chk_sum

end subroutine write_chksum_3d
! </SUBROUTINE> NAME="write_chksum_3d">

!#######################################################################
! <SUBROUTINE NAME="write_chksum_3d_int">
!
! <DESCRIPTION>
! Write a 3d integer checksum.
! </DESCRIPTION>
!
subroutine write_chksum_3d_int(name, data)
  character(len=*), intent(in) :: name
  integer, dimension(isc:iec,jsc:jec,nk) :: data
  integer :: stdoutunit
  character(len=40) :: c
  c = '[chksum] '//name

  stdoutunit=stdout()
  write(stdoutunit, '(a40,i30)') c, mpp_chksum(data(isc:iec,jsc:jec,:nk))

end subroutine write_chksum_3d_int
! </SUBROUTINE> NAME="write_chksum_3d_int">


!#######################################################################
! <SUBROUTINE NAME="write_chksum_2d">
!
! <DESCRIPTION>
! Write a 2d checksum.
! </DESCRIPTION>
!
subroutine write_chksum_2d(name, data)

  character(len=*), intent(in) :: name
  real, dimension(isc:iec,jsc:jec) :: data
  integer :: stdoutunit
  character(len=40) :: c
  c = '[chksum] '//name

  stdoutunit=stdout()
  write(stdoutunit, '(a40,i30)') c, mpp_chksum(data(isc:iec,jsc:jec))


end subroutine write_chksum_2d
! </SUBROUTINE> NAME="write_chksum_2d">


!#######################################################################
! <SUBROUTINE NAME="write_chksum_2d_int">
!
! <DESCRIPTION>
! Write a 2d integer checksum.
! </DESCRIPTION>
!
subroutine write_chksum_2d_int(name, data)

  character(len=*), intent(in) :: name
  integer, dimension(isc:iec,jsc:jec) :: data
  integer :: stdoutunit
  character(len=40) :: c
  c = '[chksum] '//name

  stdoutunit=stdout()
  write(stdoutunit, '(a40,i30)') c, mpp_chksum(data(isc:iec,jsc:jec))

end subroutine write_chksum_2d_int
! </SUBROUTINE> NAME="write_chksum_2d_int">


!#######################################################################
! <SUBROUTINE NAME="check_restart">
!
! <DESCRIPTION>
! Write a note regarding the restart setup.
! </DESCRIPTION>
!
subroutine check_restart(filename, write_a_restart)
  character(len=*), intent(in) :: filename
  logical, intent(in) :: write_a_restart

  if(.not. write_a_restart) then 
     call write_note(filename, 'running with write_a_restart=.false.')
     call write_line('Will NOT write restart file, and so cannot restart the run.')
  endif

end subroutine check_restart
! </SUBROUTINE> NAME="check_restart">


!#######################################################################
! <SUBROUTINE NAME="write_summary">
!
! <DESCRIPTION>
! Write a summary note.
! </DESCRIPTION>
!
subroutine write_summary(filename, msg, model_time)

  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: msg
  type(time_type),  intent(in) :: model_time  
  integer :: stdoutunit

  stdoutunit=stdout()
  write(stdoutunit, '(/a)') '[SUMMARY] '//filename//': '//msg
  call write_timestamp(model_time)
  
end subroutine write_summary
! </SUBROUTINE> NAME="write_summary">


!#######################################################################
! <SUBROUTINE NAME="diagnose_3d">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the grid tmask.
! </DESCRIPTION>
!
subroutine diagnose_3d(Time, Grid, id_name, data, nk_lim, use_mask, abs_max, abs_min)
    type(ocean_time_type),        intent(in) :: Time
    type(ocean_grid_type),        intent(in) :: Grid
    integer,                      intent(in) :: id_name
    real, dimension(isd:,jsd:,:), intent(in) :: data
    integer,                      intent(in), optional :: nk_lim
    logical,                      intent(in), optional :: use_mask
    real,                         intent(in), optional :: abs_max
    real,                         intent(in), optional :: abs_min

    call diagnose_3d_mask(Time, Grid%tmask(:,:,:), id_name, data, nk_lim, use_mask, abs_max, abs_min)

end subroutine diagnose_3d
! </SUBROUTINE> NAME="diagnose_3d"

!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_u">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the grid umask.
! </DESCRIPTION>
!
subroutine diagnose_3d_u(Time, Grid, id_name, data, nk_lim, use_mask, abs_max, abs_min)
    type(ocean_time_type),        intent(in) :: Time
    type(ocean_grid_type),        intent(in) :: Grid
    integer,                      intent(in) :: id_name
    real, dimension(isd:,jsd:,:), intent(in) :: data
    integer,                      intent(in), optional :: nk_lim
    logical,                      intent(in), optional :: use_mask
    real,                         intent(in), optional :: abs_max
    real,                         intent(in), optional :: abs_min

    call diagnose_3d_mask(Time, Grid%umask(:,:,:), id_name, data, nk_lim, use_mask, abs_max, abs_min)

end subroutine diagnose_3d_u
! </SUBROUTINE> NAME="diagnose_3d_u"

!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_en">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the grid tmasken.
! </DESCRIPTION>
!
subroutine diagnose_3d_en(Time, Grid, id_name1, id_name2, data, nk_lim, use_mask, abs_max, abs_min)
    type(ocean_time_type),        intent(in) :: Time
    type(ocean_grid_type),        intent(in) :: Grid
    integer,                      intent(in) :: id_name1
    integer,                      intent(in) :: id_name2
    real, dimension(isd:,jsd:,:,:), intent(in) :: data
    integer,                      intent(in), optional :: nk_lim
    logical,                      intent(in), optional :: use_mask
    real,                         intent(in), optional :: abs_max
    real,                         intent(in), optional :: abs_min

    call diagnose_3d_mask(Time, Grid%tmasken(:,:,:,1), id_name1, data(:,:,:,1), nk_lim, use_mask, abs_max, abs_min)
    call diagnose_3d_mask(Time, Grid%tmasken(:,:,:,2), id_name2, data(:,:,:,2), nk_lim, use_mask, abs_max, abs_min)

end subroutine diagnose_3d_en
! </SUBROUTINE> NAME="diagnose_3d_en"

!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_mask">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the given mask.
! </DESCRIPTION>
!
subroutine diagnose_3d_mask(Time, mask, id_name, data, nk_lim, use_mask, abs_max, abs_min)
    type(ocean_time_type),        intent(in) :: Time
    real, dimension(isd:,jsd:,:), intent(in) :: mask
    integer,                      intent(in) :: id_name
    real, dimension(isd:,jsd:,:), intent(in) :: data
    integer,                      intent(in), optional :: nk_lim
    logical,                      intent(in), optional :: use_mask
    real,                         intent(in), optional :: abs_max
    real,                         intent(in), optional :: abs_min

    logical :: use_mask_, used
    integer :: nk_lim_
    real, dimension(isd:ied,jsd:jed,nk) :: threshold_mask

    nk_lim_ = nk
    if (present(nk_lim)) then
       nk_lim_ = nk_lim
    endif

    use_mask_ = .true.
    if (present(use_mask)) then
       use_mask_ = use_mask
    endif

    threshold_mask(:,:,:) = 1.0
    if (present(abs_max)) then
       where (abs(data(COMP,:)) > abs_max) threshold_mask(COMP,:) = 0.0
    endif
    if (present(abs_min)) then
       where (abs(data(COMP,:)) < abs_min) threshold_mask(COMP,:) = 0.0
    endif

    if (id_name > 0) then
       if (use_mask_) then
          if (nk_lim_ /= nk) then
             call mpp_error(FATAL, &
                  '==> Error from ocean_nphysics_util_new (diagnose_3d): nk_lim must equal nk.')
          endif
          threshold_mask(:,:,:) = 1.0
          if (present(abs_max)) then
             where (abs(data(COMP,:)) > abs_max) threshold_mask(COMP,:) = 0.0
          endif
          if (present(abs_min)) then
             where (abs(data(COMP,:)) < abs_min) threshold_mask(COMP,:) = 0.0
          endif
          used = send_data(id_name, data(:,:,:),                              &
               Time%model_time, rmask=mask(:,:,:)*threshold_mask(:,:,:),&
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk_lim_)
       else
          used = send_data(id_name, data(:,:,:), &
               Time%model_time,                  &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk_lim_)
       endif
    endif

end subroutine diagnose_3d_mask
! </SUBROUTINE> NAME="diagnose_3d_mask"


!#######################################################################
! <SUBROUTINE NAME="diagnose_2d">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using the grid tmask.
! </DESCRIPTION>
!
subroutine diagnose_2d(Time, Grid, id_name, data, abs_max, abs_min)
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_grid_type),      intent(in) :: Grid
    integer,                    intent(in) :: id_name
    real, dimension(isd:,jsd:), intent(in) :: data
    real,                       intent(in), optional :: abs_max
    real,                       intent(in), optional :: abs_min

    call diagnose_2d_mask(Time, Grid%tmask(:,:,1), id_name, data, abs_max, abs_min)

end subroutine diagnose_2d
! </SUBROUTINE> NAME="diagnose_2d"

!#######################################################################
! <SUBROUTINE NAME="diagnose_2d_u">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using the grid umask.
! </DESCRIPTION>
!
subroutine diagnose_2d_u(Time, Grid, id_name, data, abs_max, abs_min)
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_grid_type),      intent(in) :: Grid
    integer,                    intent(in) :: id_name
    real, dimension(isd:,jsd:), intent(in) :: data
    real,                       intent(in), optional :: abs_max
    real,                       intent(in), optional :: abs_min

    call diagnose_2d_mask(Time, Grid%umask(:,:,1), id_name, data, abs_max, abs_min)

  end subroutine diagnose_2d_u
! </SUBROUTINE> NAME="diagnose_2d_u"

!#######################################################################
! <SUBROUTINE NAME="diagnose_2d_en">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using the grid tmasken.
! </DESCRIPTION>
!
subroutine diagnose_2d_en(Time, Grid, id_name1, id_name2, data, abs_max, abs_min)
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_grid_type),      intent(in) :: Grid
    integer,                    intent(in) :: id_name1
    integer,                    intent(in) :: id_name2
    real, dimension(isd:,jsd:,:), intent(in) :: data
    real,                       intent(in), optional :: abs_max
    real,                       intent(in), optional :: abs_min

    call diagnose_2d_mask(Time, Grid%tmasken(:,:,1,1), id_name1, data(:,:,1), abs_max, abs_min)
    call diagnose_2d_mask(Time, Grid%tmasken(:,:,1,2), id_name2, data(:,:,2), abs_max, abs_min)

  end subroutine diagnose_2d_en
! </SUBROUTINE> NAME="diagnose_2d_en"

!#######################################################################
! <SUBROUTINE NAME="diagnose_2d_mask">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using a given mask.
! </DESCRIPTION>
!
subroutine diagnose_2d_mask(Time, mask, id_name, data, abs_max, abs_min)
    type(ocean_time_type),      intent(in) :: Time
    real, dimension(isd:,jsd:), intent(in) :: mask
    integer,                    intent(in) :: id_name
    real, dimension(isd:,jsd:), intent(in) :: data
    real,                       intent(in), optional :: abs_max
    real,                       intent(in), optional :: abs_min

    logical :: used
    real, dimension(isd:ied,jsd:jed) :: threshold_mask

    threshold_mask(:,:) = 1.0
    if (present(abs_max)) then
       where (abs(data(COMP)) > abs_max) threshold_mask(COMP) = 0.0
    endif
    if (present(abs_min)) then
       where (abs(data(COMP)) < abs_min) threshold_mask(COMP) = 0.0
    endif

    if (id_name > 0) used = send_data(id_name, data(:,:),             &
         Time%model_time, rmask=mask(:,:)*threshold_mask(:,:),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

end subroutine diagnose_2d_mask
! </SUBROUTINE> NAME="diagnose_2d_mask"

!#######################################################################
! <SUBROUTINE NAME="diagnose_2d_int">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using the grid tmask.
! </DESCRIPTION>
!
subroutine diagnose_2d_int(Time, Grid, id_name, data)
    type(ocean_time_type),         intent(in) :: Time
    type(ocean_grid_type),         intent(in) :: Grid
    integer,                       intent(in) :: id_name
    integer, dimension(isd:,jsd:), intent(in) :: data

    real, dimension(isd:ied,jsd:jed) :: data_

    data_ = 1.0*data
    call diagnose_2d(Time, Grid, id_name, data_)

end subroutine diagnose_2d_int
! </SUBROUTINE> NAME="diagnose_2d_int"

!#######################################################################
! <SUBROUTINE NAME="diagnose_2d_comp">
!
! <DESCRIPTION>
! Helper function for diagnosting 2D data using the grid tmask on the 
! computational domain.
! </DESCRIPTION>
!
subroutine diagnose_2d_comp(Time, Grid, id_name, data)
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_grid_type),      intent(in) :: Grid
    integer,                    intent(in) :: id_name
    real, dimension(isc:,jsc:), intent(in) :: data

    logical :: used

    if (id_name > 0) used = send_data(id_name, data, &
         Time%model_time, rmask=Grid%tmask(isc:iec,jsc:jec,1))

end subroutine diagnose_2d_comp
! </SUBROUTINE> NAME="diagnose_2d_comp"

!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_int">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the grid tmask.
! </DESCRIPTION>
!
subroutine diagnose_3d_int(Time, Grid, id_name, data)
    type(ocean_time_type),           intent(in) :: Time
    type(ocean_grid_type),           intent(in) :: Grid
    integer,                         intent(in) :: id_name
    integer, dimension(isd:,jsd:,:), intent(in) :: data

    real, dimension(isd:ied,jsd:jed,nk) :: data_

    data_ = 1.0*data
    call diagnose_3d(Time, Grid, id_name, data_)

end subroutine diagnose_3d_int
! </SUBROUTINE> NAME="diagnose_3d_int"

!#######################################################################
! <SUBROUTINE NAME="diagnose_3d_comp">
!
! <DESCRIPTION>
! Helper function for diagnosting 3D data using the grid tmask on the 
! computational domain.
! </DESCRIPTION>
!
subroutine diagnose_3d_comp(Time, Grid, id_name, data)
    type(ocean_time_type),      intent(in) :: Time
    type(ocean_grid_type),      intent(in) :: Grid
    integer,                    intent(in) :: id_name
    real, dimension(isc:,jsc:,:), intent(in) :: data

    logical :: used

    if (id_name > 0) used = send_data(id_name, data, &
         Time%model_time, rmask=Grid%tmask(isc:iec,jsc:jec,:))

end subroutine diagnose_3d_comp
! </SUBROUTINE> NAME="diagnose_3d_comp"


subroutine diagnose_sum(Time, Grid, Dom, id_name, data, factor)

  type(ocean_time_type), intent(in) :: Time
  type(ocean_grid_type), intent(in) :: Grid
  type(ocean_domain_type), intent(in) :: Dom
  integer, intent(in) :: id_name
  real, dimension(isd:,jsd:), intent(in) :: data
  real, intent(in) :: factor

  real, dimension(isd:ied,jsd:jed) :: work
  real :: total
  logical :: used

  if(id_name > 0) then 
     work(:,:) = Grid%tmask(:,:,1)*Grid%dat(:,:)*data(:,:)
     total = mpp_global_sum(Dom%domain2d, work(:,:), NON_BITWISE_EXACT_SUM)*factor
     used = send_data (id_name, total, Time%model_time)
  endif

end subroutine diagnose_sum

function register_2d_t_field(Grid, Time, name, desc, units, range)

    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    character(len=*),      intent(in) :: name
    character(len=*),      intent(in) :: desc
    character(len=*),      intent(in) :: units
    real, dimension(2),    intent(in) :: range

    integer :: register_2d_t_field

    register_2d_t_field = register_diag_field ('ocean_model', name, &
         Grid%tracer_axes(1:2), Time%model_time, desc, units,       &
         missing_value=missing_value, range=range)

end function register_2d_t_field

function register_3d_t_field(Grid, Time, name, desc, units, range)

    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    character(len=*),      intent(in) :: name
    character(len=*),      intent(in) :: desc
    character(len=*),      intent(in) :: units
    real, dimension(2),    intent(in) :: range

    integer :: register_3d_t_field

    register_3d_t_field = register_diag_field ('ocean_model', name, &
         Grid%tracer_axes(1:3), Time%model_time, desc, units,       &
         missing_value=missing_value, range=range)

end function register_3d_t_field

function register_3d_xte_field(Grid, Time, name, desc, units, range)

    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    character(len=*),      intent(in) :: name
    character(len=*),      intent(in) :: desc
    character(len=*),      intent(in) :: units
    real, dimension(2),    intent(in) :: range

    integer :: register_3d_xte_field

    register_3d_xte_field = register_diag_field ('ocean_model', name, &
         Grid%tracer_axes_flux_x(1:3), Time%model_time, desc, units,  &
         missing_value=missing_value, range=range)

end function register_3d_xte_field

function register_3d_ytn_field(Grid, Time, name, desc, units, range)

    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    character(len=*),      intent(in) :: name
    character(len=*),      intent(in) :: desc
    character(len=*),      intent(in) :: units
    real, dimension(2),    intent(in) :: range

    integer :: register_3d_ytn_field

    register_3d_ytn_field = register_diag_field ('ocean_model', name, &
         Grid%tracer_axes_flux_y(1:3), Time%model_time, desc, units,  &
         missing_value=missing_value, range=range)

end function register_3d_ytn_field

function register_3d_ztb_field(Grid, Time, name, desc, units, range)

    type(ocean_grid_type), intent(in) :: Grid
    type(ocean_time_type), intent(in) :: Time
    character(len=*),      intent(in) :: name
    character(len=*),      intent(in) :: desc
    character(len=*),      intent(in) :: units
    real, dimension(2),    intent(in) :: range

    integer :: register_3d_ztb_field

    register_3d_ztb_field = register_diag_field ('ocean_model', name, &
         Grid%tracer_axes_wt(1:3), Time%model_time, desc, units,      &
         missing_value=missing_value, range=range)

end function register_3d_ztb_field



end module ocean_util_mod

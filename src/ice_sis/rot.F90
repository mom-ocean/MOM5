!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! rot_mod - rotates a spherical coordinate system, from M. Eby via MOM         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module rot_mod

  use mpp_mod,        only: stdout
  use axis_utils_mod, only: nearest_index

  implicit none
  private
  public rot_intrp_sclr, rot_intrp_vctr, rotate, psir, thetar, phir
  ! rotation angles set in ice_grid

  real               :: psir = 0.0, thetar = 0.0, phir = 0.0

  integer            :: warns  = 0 ! > 0 allows #warns warnings from interpolator


contains

  subroutine rot_intrp_vctr (g, xg, yg, ig, jg, r, xr, yr, ir, jr, &
       psir, thetar, phir)
    !
    !=======================================================================
    !     interpolate vector data from an geographic data grid to a
    !     rotated model grid
    !
    !     input
    !     psir, thetar, phir = Euler angles defining rotation
    !     g  = vector on geographic data grid
    !     xg = longitude of data points on geographic data grid
    !     yg = latitude of data points on geographic data grid
    !     ig = number of longitudes in geographic data grid 
    !     jg = number of latitudes in geographic data grid 
    !     xr = longitude of points on rotated model grid
    !     yr = latitude of points on rotated model grid
    !     ir = number of longitudes in rotated model grid 
    !     jr = number of latitudes in rotated model grid 
    !
    !     output
    !     r  = vector on rotated model grid
    !
    !     internal
    !     (rln,rlt) = (longitude,latitude) in rotated coordinates
    !     (gln,glt) = (longitude,latitude) in geographic coordinates
    !     xg(iw) = point on the geographic grid to the west of (gln,glt)
    !     xg(ie) = point on the geographic grid to the east of (gln,glt)
    !     yg(js) = point on the geographic grid to the south of (gln,glt)
    !     yg(jn) = point on the geographic grid to the north of (gln,glt)
    !
    !=======================================================================
    !
  integer                  :: ig, jg, ir, jr
  real, dimension(ig,jg,2) ::  g
  real, dimension(ig     ) :: xg
  real, dimension(   jg  ) :: yg
    real, dimension(ir,jr,2) ::  r
  real, dimension(ir     ) :: xr
  real, dimension(   jr  ) :: yr
  real                     :: psir, thetar, phir

    integer                  :: i, j
    real                     :: rad, a, angle, vmag, glt, gln, rlt, rln
    !
    rad = acos(-1.)/180.
    !
    ! interpolate vector components as scalers on rotated model grid 
    !
    call rot_intrp_sclr (g(1,1,1), xg, yg, ig, jg, r(1,1,1), xr, yr, ir, jr, &
         psir, thetar, phir)
    call rot_intrp_sclr (g(1,1,2), xg, yg, ig, jg, r(1,1,2), xr, yr, ir, jr, &
         psir, thetar, phir)
    !
    ! correct vector direction
    !
    do j=1,jr
       do i=1,ir
          vmag = sqrt(r(i,j,1)**2 + r(i,j,2)**2)
          if (vmag .gt. 0.) then 
             a = r(i,j,1)/vmag
             a = min(a, 1.)
             a = max(a, -1.)
             a = acos(a)
             if (r(i,j,2) .lt. 0.) a = -a
             call rotate (yr(j), xr(i), -psir, -thetar, -phir, glt, gln)
             call rotvec(glt, gln, phir, thetar, psir, angle)
             a = a + angle*rad
             r(i,j,1) = vmag*cos(a)
             r(i,j,2) = vmag*sin(a)
          else
             r(i,j,1) = 0.
             r(i,j,2) = 0.
          endif
       enddo
    enddo
    !
    return
  end subroutine rot_intrp_vctr


  subroutine rot_intrp_sclr (g, xg, yg, ig, jg, r, xr, yr, &
       ir, jr, psir, thetar, phir    )
    !
    !=======================================================================
    !     interpolate scaler data from an geographic data grid to a 
    !     rotated model grid
    !
    !     input
    !     psir, thetar, phir = Euler angles defining rotation
    !     g  = scaler on geographic data grid
    !     xg = longitude of data points on geographic data grid
    !     yg = latitude of data points on geographic data grid
    !     ig = number of longitudes in on geographic data grid 
    !     jg = number of latitudes in on geographic data grid 
    !     xr = longitude of points on rotated model grid
    !     yr = latitude of points on rotated model grid
    !     ir = number of longitudes in rotated model grid 
    !     jr = number of latitudes in rotated model grid 
    !
    !     output
    !     r  = scaler on rotated model grid
    !
    !     internal
    !     (rln,rlt) = (longitude,latitude) in rotated coordinates
    !     (gln,glt) = (longitude,latitude) in geographic coordinates
    !     xg(iw) = point on the geographic grid to the west of (gln,glt)
    !     xg(ie) = point on the geographic grid to the east of (gln,glt)
    !     yg(js) = point on the geographic grid to the south of (gln,glt)
    !     yg(jn) = point on the geographic grid to the north of (gln,glt)
    !
    !=======================================================================
    !
    integer                :: ig, jg, ir, jr, outunit
    real                   :: psir, thetar, phir
    real, dimension(ig,jg) :: g
    real, dimension(ig   ) :: xg
    real, dimension(   jg) :: yg
    real, dimension(ir,jr) :: r
    real, dimension(ir   ) :: xr
    real, dimension(   jr) :: yr

    real    :: rln, rlt, gln, glt
    integer :: i, j, iw, ie, js, jn, istrt, jstrt, iend, jend, ln_err, lt_err
    real    :: epsln, glt_min, glt_max, gln_min, gln_max, del, wtw, wte, wts, wtn
    !
    epsln = 1.e-10
    glt_min = 90.
    glt_max = -90.
    gln_min = 360.
    gln_max = -360.
    ln_err = 0
    lt_err = 0
    !
    !     find longitude points of data within interval [0., 360.]
    istrt = 1
    do i=2,ig
       if (xg(i-1) .lt. 0. .and. xg(i) .ge. 0.) istrt = i
    enddo
    iend = ig
    do i=2,ig
       if (xg(i-1) .lt. 360. .and. xg(i) .ge. 360.) iend = i
    enddo
    !
    !     find latitude points of data within interval [-90., 90.]
    jstrt = 1
    do j=2,jg
       if (yg(j-1) .lt. -90. .and. yg(j) .ge. -90.) jstrt = j
    enddo
    jend = jg
    do j=2,jg
       if (yg(j-1) .lt. 90. .and. yg(j) .ge. 90.) jend = j
    enddo
    !
    !     interpolate data to model grid 
    !
    do j=1,jr
       do i=1,ir
          call rotate (yr(j), xr(i), -psir, -thetar, -phir, glt, gln)
          if (gln .lt. 0.) gln = gln + 360.
          if (gln .ge. 360.) gln = gln - 360.
          glt_min = min(glt,glt_min)
          glt_max = max(glt,glt_max)
          gln_min = min(gln,gln_min)
          gln_max = max(gln,gln_max)
          !
          iw = nearest_index (gln, xg(istrt:iend) ) + istrt - 1
          if (xg(iw) .gt. gln) iw = iw - 1
          ie = iw + 1
          if (iw .ge. istrt .and. ie .le. iend) then
             del = xg(ie) - xg(iw)
             wtw = (xg(ie) - gln)/del
          else
             !     east or west of the last data value. this could be because a
             !     cyclic condition is needed or the dataset is too small. in either 
             !     case apply a cyclic condition
             ln_err = 1
             iw = iend
             ie = istrt
             del = xg(ie) + 360. + epsln - xg(iw) 
             if (xg(ie) .ge. gln) then
                wtw = (xg(ie) - gln)/del
             else
                wtw = (xg(ie) + 360. + epsln - gln)/del
             endif
          endif
          wte = 1. - wtw
          !
          js = nearest_index (glt, yg(jstrt:jend) ) + jstrt - 1
          if (yg(js) .gt. glt) js = max(js - 1,jstrt)
          jn = min(js + 1,jend)
          if (yg(jn) .ne. yg(js) .and. yg(js) .le. glt) then
             wts = (yg(jn) - glt)/(yg(jn) - yg(js))
          else
             !     north or south of the last data value. this could be because a
             !     pole is not included in the data set or the dataset is too small.
             !     in either case extrapolate north or south
             lt_err = 1
             wts = 1.
          endif
          wtn = 1. - wts
          !
          r(i,j) = g(ie,jn)*wte*wtn + g(ie,js)*wte*wts &
               + g(iw,jn)*wtw*wtn + g(iw,js)*wtw*wts
          !
       enddo
    enddo
    !
    outunit = stdout()
    if (ln_err .eq. 1 .and. warns > 0) then
       write (outunit,'(/,(1x,a))')                                      &
            '==> Warning: the geographic data set does not extend far   ', &
            '             enough east or west - a cyclic boundary       ', &
            '             condition was applied. check if appropriate   '
       write (outunit,'(/,(1x,a,2f8.2))')                                &
            '    data required between longitudes:', gln_min, gln_max,     &
            '      data set is between longitudes:', xg(istrt), xg(iend)
       warns = warns - 1
    endif
    !
    if (lt_err .eq. 1 .and. warns > 0) then
       write (outunit,'(/,(1x,a))')                                     &
            '==> Warning: the geographic data set does not extend far   ',&
            '             enough north or south - extrapolation from    ',&
            '             the nearest data was applied. this may create ',&
            '             artificial gradients near a geographic pole   ' 
       write (outunit,'(/,(1x,a,2f8.2))')                             &
            '    data required between latitudes:', glt_min, glt_max,   &
            '      data set is between latitudes:', yg(jstrt), yg(jend)
       warns = warns - 1
    endif
    !
    return
  end subroutine rot_intrp_sclr


  subroutine rotate (glt, gln, phir, thetar, psir, rlt, rln)
    real :: glt, gln, phir, thetar, psir, rlt, rln
    !
    !=======================================================================
    !     subroutine rotate takes a geographic latitude and longitude and 
    !     finds the the equivalent latitude and longitude on a rotated grid.
    !     when going from a geographic grid to a rotated grid, all of the 
    !     defined rotation angles given to rotate by the calling program 
    !     are positive, but when going from a rotated grid back to the 
    !     geographic grid, the calling program must reverse the angle order 
    !     (phir and psir are switched) and all of the angles made negative.
    !
    !     the first rotation angle phir is defined as a rotation about the
    !     original z axis. the second rotation angle thetar is defined as a
    !     rotation about the new x axis. the final rotation angle psir is
    !     defined as a rotation about the new z axis. these rotation angles
    !     are just the Euler angles as defined in "classical mechanics"
    !     Goldstein (1951).
    !
    !     author:   M. Eby            e-mail eby@uvic.ca
    !=======================================================================
    !
    !     g...  = geographic value
    !     r...  = rotated value
    !     ...lt = latitude (or equivalent spherical coordinate)
    !     ...ln = longitude (or equivalent spherical coordinate)
    !     ...x  = x coordinate
    !     ...y  = y coordinate
    !     ...z  = z coordinate
    !     psir, thetar, phir = Euler angles defining rotation
    !
    real rad, gx, gy, gz, phis, thetas, rx, ry, rz
    !     define rad for conversion to radians.
    rad = acos(-1.)/180.
    !
    !     convert latitude and longitude to spherical coordinates
    thetas = gln
    if (thetas .gt. 180.) thetas = thetas - 360.
    if (thetas .lt. -180.) thetas = thetas + 360.
    phis = (90. - glt)*rad
    thetas = thetas*rad
    !
    !     translate point into Cartesian coordinates for rotation.
    gx = sin(phis)*cos(thetas)
    gy = sin(phis)*sin(thetas)
    gz = cos(phis)
    !
    !     rotate the point (gx, gy, gz) about the z axis by phir then the x
    !     axis by thetar and finally about the z axis by psir.
    ! 
    rx = gx*(cos(psir)*cos(phir) - cos(thetar)*sin(phir)*sin(psir)) + &
         gy*(cos(psir)*sin(phir) + cos(thetar)*cos(phir)*sin(psir)) + &
         gz*sin(psir)*sin(thetar)
    !
    ry = gx*(-sin(psir)*cos(phir) - cos(thetar)*sin(phir)*cos(psir)) + &
         gy*(-sin(psir)*sin(phir) + cos(thetar)*cos(phir)*cos(psir)) + &
         gz*(cos(psir)*sin(thetar))
    !
    rz = gx*(sin(thetar)*sin(phir)) + gy*(-sin(thetar)*cos(phir)) +  &
         gz*(cos(thetar))
    !
    !     convert rotated point back to spherical coordinates
    !
    !     check for rounding error (arccos(x): abs(x) must be .le. 1)
    rz = min(rz, 1.)
    rz = max(rz, -1.)
    rlt = acos(rz)
    !     if point is at a pole set rotated longitude equal to initial.
    if (rlt .le. 0. .or. rlt .ge. 180.*rad) then
       rln = thetas
    else
       !     if rln lies between -135 and -45 or between 45 and 135 degrees
       !     it is more accurate to use an arccos calculation.
       if (abs(rx/sin(rlt)) .lt. cos(45.*rad)) then
          rln = rx/sin(rlt)
          !     check for rounding error (arccos(x): abs(x) must be .le. 1)
          rln = min(rln, 1.)
          rln = max(rln, -1.)
          rln = acos(rln)
          !     arccos will give rln between 0 and 180 degrees.  if the point
          !     is negative in y, rln must be equal to negative rln.
          if (ry .lt. 0.) rln = -rln
       else
          !     if rln lies between -45 and 45 or between 135 and -135 degrees
          !     it is more accurate to use an arcsin calculation.
          rln = ry/sin(rlt)
          !     check for rounding error (arcsin(x): abs(x) must be .le. 1)
          rln = min(rln, 1.)
          rln = max(rln, -1.)
          rln = asin(rln)
          !     arcsin will give rln between -90 and 90 degrees. if the point
          !     is negative in x, rln must be equal to 180 degrees minus rln.
          if (rx .lt. 0.) rln = 180.*rad - rln
       endif
    endif
    !
    !     convert back to degrees of latitude and longitude.
    rlt = 90. - rlt/rad
    rln = rln/rad
    if (rln .gt. 180.) rln = rln - 360.
    if (rln .le. -180.) rln = rln + 360.
    !
    return
  end subroutine rotate

  subroutine rotvec (glt, gln, phir, thetar, psir, angle)
    real :: glt, gln, phir, thetar, psir, angle

    real rad, delta, rlt, rln, rlth, rlnh, glth, glnh, dst, t
    !
    !=======================================================================
    !     subroutine rotvec takes a geographic latitude and longitude and 
    !     finds the the vector rotation angle angle (in degrees) for a 
    !     vector on the rotated grid. when going from the geographic to a 
    !     rotated grid, all of the defined rotation angles given to rotvec 
    !     by the calling program are positive, but when going from a 
    !     rotated grid back to the geographic grid, the calling program 
    !     must reverse the angle order (phir and psir are switched) and all 
    !     of the angles made negative. if a pole is detected then an angle 
    !     of zero is returned.
    !
    !     rotvec rotates the point defining the head of a very short north
    !     or south pointing geographic vector. the angle between this 
    !     vector and a similar direction vector defined in the new grid is 
    !     calculated using the law of cosines. the accuracy of this 
    !     calculation depends on the size of the direction vector (delta) 
    !     and the precision of the computation. the smaller the vector the 
    !     more accurate the calculation but then the more precision 
    !     required. double precision is strongly recommended.
    !
    !     author:   M. Eby            e-mail eby@uvic.ca
    !=======================================================================
    !
    ! include "stdunits.h"
    !     if gx and gy are vector components at glt and gln, the corrected 
    !     vector components rx and ry at rlt and rln can be calculated as 
    !     follows:
    !
    !  rad = acos(-1.)/180.
    !  call rotate (rlt, rln, -psir, -thetar, -phir, glt, gln)
    !  call get_vector (glt, gln, gx, gy)   ! some routine to get vector
    !  r = sqrt(gx**2 + gy**2)
    !  if (r .gt. 0.) then 
    !    a = gx/r
    !    a = min(a, 1.)
    !    a = max(a, -1.)
    !    a = acos(a)
    !    if (gy .lt. 0.) a = -a
    !    call rotvec(glt, gln, phir, thetar, psir, angle)
    !    a = a + angle*rad
    !    rx = r*cos(a)
    !    ry = r*sin(a)
    !  endif
    !
    ! g...   = geographic value
    ! r...   = rotated value
    ! ...lt  = latitude
    ! ...ln  = longitude
    ! ...lth = latitude of head of vector
    ! ...lnh = longitude of head of vector
    ! dst    = distance between heads of vectors
    ! delta  = length of vector
    ! angle  = angle to rotate vectors back to original orientation
    ! psir, thetar, phir = Euler angles defining rotation
    !
    ! define multiplier rad for conversion to radians.
    rad = acos(-1.)/180.
    !
    ! define length of direction vector (single precision may require a 
    ! longer and thus less accurate vector length).
    delta = 0.001
    !
    ! if the base is in the north of the geographic grid use a south 
    ! pointing vector to avoid any possibility of going over the pole.
    if (glt .ge. 0.) delta = -delta
    !
    ! find the base of the direction vectors in the rotated grid.
    call rotate (glt, gln, phir, thetar, psir, rlt, rln)
    !
    ! if base in the rotated grid is near a pole return an angle of zero
    if (abs(rlt) .ge. 90.-abs(delta)) then
       angle = 0.
       return
    endif
    !
    ! find the head of the geographic grid direction vector in the 
    ! rotated grid.
    call rotate (glt+delta, gln, phir, thetar, psir, glth, glnh)
    !
    ! if the base is in opposite hemispheres switch the vector
    ! direction for better accuracy.
    if (glt*rlt .lt. 0) delta = -delta
    !
    ! find the head of the rotated grid direction vector.
    rlth = rlt + delta
    rlnh = rln
    !
    ! find the distance between the heads of the direction vectors.
    call dist (glth, glnh, rlth, rlnh, dst)
    !
    ! find the angle between direction vectors with the law of cosines.
    delta = abs(delta)
    angle = (cos(dst)-cos(delta)**2)/(sin(delta)**2)
    angle = min(angle, 1.)
    angle = max(angle, -1.)
    angle = acos(angle)/rad
    t = abs(delta)
    !
    ! adjust the angle if the direction vectors are opposite.
    if (glt*rlt .lt. 0) angle = 180. - angle
    ! determine the sign of the angle by checking the offset longitudes.
    if (glnh+360. .gt. rlnh+360.) angle = -angle
    ! change sign if the original direction vector was pointing south.
    if (glt .ge. 0.) angle = -angle
    !
    return
  end subroutine rotvec


  subroutine dist (lat1, lng1, lat2, lng2, dst)
    !
    !=======================================================================
    !     subroutine dist calculates the arc distance between two 
    !     points given their latitudes and longitudes
    !=======================================================================
    !
    real lat1, lng1, lat2, lng2, dst, rad
    !
    !     define multiplier rad for conversion to radians.
    rad = acos(-1.)/180.
    !
    !     check input.
    lat1 = min(lat1, 90.)
    lat1 = max(lat1, -90.)
    if (lng1 .lt. 0.) lng1 = lng1 + 360.
    if (lng1 .gt. 360.) lng1 = lng1 - 360.
    lat2 = min(lat2, 90.)
    lat2 = max(lat2, -90.)
    if (lng2 .lt. 0.) lng2 = lng2 + 360.
    if (lng2 .gt. 360.) lng2 = lng2 - 360.
    !
    dst = sin(lat1*rad)*sin(lat2*rad)+cos(lat1*rad) &
         *cos(lat2*rad)*cos((lng1-lng2)*rad)
    dst = min(dst, 1.)
    dst = max(dst, -1.)
    dst = (acos(dst)/rad)
    return
  end subroutine dist

end module rot_mod

      real function gmean(im, jm, jfirst, jlast, q)
      use fv_pack, only: grid_weight
      implicit none

! !INPUT PARAMETERS:
      integer  im, jm                        ! Horizontal dimensions
      integer  jfirst, jlast                 ! Latitude strip
      real q(im,jfirst:jlast)            ! 2D field 

      integer i, j
      real xsum(jfirst:jlast)

          xsum = 0.
      do j=jfirst,jlast
        do i=1,im
          xsum(j) = xsum(j) + q(i,j)
        enddo
          xsum(j) = xsum(j)*grid_weight(j)
      enddo

      call par_vecsum( jm, jfirst, jlast, xsum, gmean)
      end

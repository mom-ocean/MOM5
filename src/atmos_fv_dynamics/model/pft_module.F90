module pft_module
#ifndef USE_LIMA
  use ecmfft
  implicit none
   private
   public pft2d, pft2d_phys, pft_init, pft_end
!
   integer ifax(13)             ! FFT internal information
   real, allocatable :: trigs(:)
   real, allocatable :: sc(:), se(:)
   real, allocatable :: damp(:,:,:)
   integer js2g0

! if "SGI_FFT" is defined SGI's FFT library will be used (instead of the
! default FFT f77 code from ECMWF)
! (Need to link to SGI library with  "-lcomplib.sgimath" )

! if "ALT_PFT" is defined the prognostic winds on the D grid will be directly
! filtered (instead of the tendencies), and there will be 2 less calls
! to pft2d. However, the faster algebraic 3-pt filter will not be used.

! Programmer: S.-J. Lin, NOAA/GFDL, Princeton, NJ
 
CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft_init --- Calculate algebraic and FFT polar filters
!
! !INTERFACE: 
 subroutine pft_init(im, jm, beglat, endlat, ighost, cosp, cose, print_msg)

! !USES:

! !INPUT PARAMETERS:
      integer, intent(in):: im        ! Global X dimension
      integer, intent(in):: jm        ! Global Y dimension
      integer, intent(in):: beglat, endlat
      integer, intent(in):: ighost    ! Normally 0
      real, intent(in)::   cosp(jm)            ! cosine array
      real, intent(in)::   cose(jm)            ! cosine array

      logical, intent(in)::  print_msg

! !DESCRIPTION:
!
!   Compute coefficients for the 3-point algebraic and the FFT
!   polar filters.
!
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipft1, ipft2
      integer jn1g1                   ! j north limit ghosted 1 (starts jm)
      integer i, j
      real    ycrit               ! critical value
      real    lat_s
      real  dl, coszc, cutoff, phi, dmp, rat
      real  pi

      cutoff = 1.e-20
      pi = 4.d0 * datan(1.d0)

      ipft1 = max(2, beglat-ighost)
      js2g0 = max(2, beglat)

      ipft2 = min(jm-1, endlat+ighost)
      jn1g1 = min(jm,   endlat+1)

      allocate( sc(ipft1:ipft2) )
      allocate( se(js2g0:jn1g1) )
      allocate( damp(im,ipft1:jn1g1,2) )
      allocate( trigs(3*im/2+1) )

      call fftfax(im, trigs)

!---------------------------------------------
! Determine ycrit such that effective DX >= DY
!---------------------------------------------
      rat = float(im)/float(2*(jm-1))
      ycrit = acos( min(0.81, rat) ) * (180./pi)

      if(print_msg) then
         write(6,*) 'Initializing polar filter .........'
#ifdef ALT_PFT
         write(6,*) 'Starting latitude for FFT pft=', ycrit
#else
         write(6,*) 'Starting latitude for algebraic pft=', ycrit
#endif
      endif

      coszc = cos(ycrit*pi/180.)

! INIT fft polar coefficients:
      dl = pi/dble(im)

      do j=ipft1, ipft2
         do i=1,im
            damp(i,j,1) = 1.
         enddo
      enddo

      do j=js2g0,jn1g1
         do i=1,im
            damp(i,j,2) = 1.
         enddo
      enddo

!************
! Cell center
!************
      do j=ipft1, ipft2
         sc(j) = (coszc/cosp(j))**2

#ifdef ALT_PFT
         if(sc(j) > 1.) then
#else
         if(sc(j) > 1. .and. sc(j) <= 2.0) then
            sc(j) =  1. +  (sc(j)-1.)/(sc(j)+1.)
         elseif(sc(j) > 2. .and. sc(j) <= 4.) then
            sc(j) =  1. +  sc(j)/(8.-sc(j))
            sc(j) = min(2., sc(j))
            lat_s = acos(cosp(j+1))*180/pi
         elseif(sc(j) > 4. ) then
#endif

! FFT filter
         do i=1,im/2
            phi = dl * i
            dmp = min((cosp(j)/coszc)/sin(phi), 1.)**2
            if(dmp < cutoff) dmp = 0.
            damp(2*i-1,j,1) = dmp
            damp(2*i  ,j,1) = dmp
         enddo

         endif
      enddo

#if ( !defined ALT_PFT )
!     if(print_msg) write(6,*) 'Switching to FFT at latitude =', lat_s
#endif

!************
! Cell edges
!************
      do j=js2g0,jn1g1
         se(j) = (coszc/cose(j))**2

#ifdef ALT_PFT
         if( se(j) > 1. ) then
#else
         if(se(j) > 1. .and. se(j) <= 2.0 ) then
            se(j) =  1. +  (se(j)-1.)/(se(j)+1.)
         elseif(se(j) > 2. .and. se(j) <= 4.) then
            se(j) =  1. +  se(j)/(8.-se(j))
            se(j) = min(2., se(j))
         elseif(se(j) > 4. ) then
#endif
! FFT
            do i=1,im/2
               phi = dl * i
               dmp = min((cose(j)/coszc)/sin(phi), 1.)**2
               if(dmp < cutoff) dmp = 0.
               damp(2*i-1,j,2) = dmp
               damp(2*i  ,j,2) = dmp
            enddo
         endif
      enddo
!EOC
 end subroutine pft_init
!-----------------------------------------------------------------------


 subroutine pft_end

   deallocate( trigs )
   deallocate( sc )
   deallocate( se )
   deallocate( damp )

 end subroutine pft_end


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft2d --- Two-dimensional fast Fourier transform
!
! !INTERFACE: 
 subroutine pft2d(p, im, jp, q1, q2, igrid)

! !USES:

! !INPUT PARAMETERS:
      integer, intent(in):: im     ! Total X dimension
      integer, intent(in):: jp     ! Total Y dimension
      integer, intent(in):: igrid  ! 1: cell center; 2: edge

! !INPUT/OUTPUT
      real, intent(inout)::  p(im,jp)           ! Array to be polar filtered
      real, intent(inout):: q1( im+2, *)        ! Work array
      real, intent(inout):: q2(*)               ! Work array

! !DESCRIPTION:
!
!   Perform a two-dimensional fast Fourier transformation.
!
! !REVISION HISTORY:
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real s(jp)             ! 3-point algebraic filter
      real ptmp(0:im+1)
      integer jf(jp)
      real     rsc, bt 
      integer i, j, n, nj


      if ( igrid == 1 ) then
         do j=1,jp
            s(j) = sc(j+js2g0-1)
         enddo
      else
         do j=1,jp
            s(j) = se(j+js2g0-1)
         enddo
      endif

      nj = 0       ! Number of latitudes to be filtered by zonal FFT

      do j=1,jp

#ifdef ALT_PFT
      if(s(j) > 1.01) then
#else
      if(s(j) > 1.01 .and. s(j) <= 2.) then

         rsc = 1./s(j)
         bt  = 0.5*(s(j)-1.)

         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = ptmp(im)
           ptmp(im+1) = ptmp(1)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo

      elseif(s(j) > 2.) then
#endif

! Packing for FFT 
           nj  = nj + 1
         jf(nj) = j

         do i=1,im
            q1(i,nj) = p(i,j)
         enddo
            q1(im+1,nj) = 0.
            q1(im+2,nj) = 0.

      endif
      enddo

      if( nj == 0) return

      call rfftmlt(q1, q2, 1, im+2, im, nj, -1)

      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n)+js2g0-1, igrid)
         enddo
      enddo

      call rfftmlt(q1, q2, 1, im+2, im, nj, 1)

      do n=1,nj
         do i=1,im
            p(i,jf(n)) = q1(i,n)
         enddo
      enddo
!EOC
 end subroutine pft2d

 subroutine pft2d_phys(p, im, jp)
! This is a very weak zonal filter for filtering the 
! tendencies from the physics

! !INPUT PARAMETERS:
    integer, intent(in):: im     ! Total X dimension
    integer, intent(in):: jp     ! Total Y dimension

! !INPUT/OUTPUT PARAMETERS:
    real, intent(inout):: p(im,jp)  ! Array to be polar filtered

! !LOCAL VARIABLES:
    real ptmp(0:im+1)
    real rsc, bt 
    integer  i, j

    do j=1,jp
      if(sc(js2g0+j-1) > 1.01 ) then
         rsc =  1./ min(2., sc(js2g0+j-1))
         bt  = 0.5*(min(2., sc(js2g0+j-1))-1.)
         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = ptmp(im)
           ptmp(im+1) = ptmp(1)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo
      endif
    enddo

 end subroutine pft2d_phys
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: rfftmlt --- Apply fast Fourier transform
!
! !INTERFACE: 
 subroutine rfftmlt(a,work,inc,jump,n,lot,isign)


! !DESCRIPTION:
!
!   Apply the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer jump, inc, n, lot, isign
#ifdef SGI_FFT
! Need to link to SGI library with  "-lcomplib.sgimath"
!
      real   a(jump,lot)
      real   work(1)                       ! Not used; here for plug reason
! Local
      integer*4 iisign,in,iinc,ijump,ilot
      integer i, j
      real scale

!-----convert to i4
      iisign = isign
      iinc = inc
      ijump = jump
      in = n
      ilot = lot

      if( iisign < 0 ) then
!-----forward
          call dzfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)

      elseif( iisign > 0 ) then
!-----backward
          call zdfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)

          scale = 1.0/float(n)
          do j=1,lot
             do i=1,jump
                a(i,j) = scale*a(i,j)
             enddo
          enddo
      endif
#else
!
! Default f77 version
!
      real a(jump*lot)
      real work((n+1)*lot)
      integer nfax, nx, ink, ibase, i, j, k, l, ia, ib, nb, m
      integer la, nh, jbase,igo
 
!     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
 
!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2
 
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
 
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
 
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
 
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
 
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
 
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign==+1) go to 30
 
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      igo=50
      if (mod(nfax,2)==1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
 
      igo=60
      go to 40
 
!     PREPROCESSING (ISIGN=+1)
 
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
 
!     COMPLEX TRANSFORM
 
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo==60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,     &
                  ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,     &
                  2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
 
      if (isign==-1) go to 130
 
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      if (mod(nfax,2)==1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
 
!     FILL IN ZEROS AT END
  110 continue
      ib=n*inc+1
!DIR$ IVDEP
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
 
!     POSTPROCESSING (ISIGN=-1):
 
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
 
  140 continue
#endif
!EOC
 end subroutine rfftmlt



 subroutine fftfax (n, trigs)

 integer n
#ifdef SGI_FFT
 real    trigs(1)
! local
 integer*4 nn

   nn=n
   call dzfftm1dui (nn,trigs)
#else
 real trigs(3*n/2+1)
 integer:: mode=3
 integer i

    call fax (ifax, n, mode)
    i = ifax(1)
    if (ifax(i+1) > 5 .or. n <= 4) ifax(1) = -99
    if (ifax(1) <= 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop 'set99'
    endif
    call fftrig (trigs, n, mode)
#endif
 end subroutine fftfax

!
! Add Lima code here.
!
#else
  use ecmfft
!BOP
!
! !PUBLIC MEMBER FUNCTIONS:
   public pft2d, pft2d_phys, pft_init, rfftmlt, fftfax 
!
! !DESCRIPTION:
!
!      This module provides fast-Fourier transforms
!
!      \begin{tabular}{|l|l|} \hline \hline
!         pft2d     &  \\ \hline
!         pft\_cf   &  \\ \hline
!         rfftmlt   &  \\ \hline
!         fftfax    &  \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.30   Lin        Integrated into this module
!   01.03.26   Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft2d --- Two-dimensional fast Fourier transform
!
! !INTERFACE: 
 subroutine pft2d(p, s, damp, im,  jp, ifax, trigs, q1, q2)
  implicit none

! !USES:

! !INPUT PARAMETERS:
      integer im                   ! Total X dimension
      integer jp                   ! Total Y dimension
      integer ifax(13)             ! FFT internal information
      real   s(jp)             ! 3-point algebraic filter
      real  damp(im,jp)        ! FFT damping coefficients
      real trigs(3*im/2+1)

! !INPUT/OUTPUT PARAMETERS:
      real q1( im+2, *)        ! Work array
      real q2(*)               ! Work array
      real  p(im,jp)           ! Array to be polar filtered

! !DESCRIPTION:
!
!   Perform a two-dimensional fast Fourier transformation.
!
! !REVISION HISTORY:
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!   02.04.05   Sawyer       Integrated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real rsc, bt 
      integer  i, j, n, nj

!Local Auto arrays:
      real ptmp(0:im+1)
!!!      real q1(  im+2, jp)
!!!      real q2( (im+1)*jp )
      integer  jf(jp)

      nj = 0

      do 200 j=1,jp

#ifdef ALT_PFT
      if(s(j) > 1.01) then
#else
      if(s(j) > 1.01 .and. s(j) <= 2.) then

         rsc = 1./s(j)
         bt  = 0.5*(s(j)-1.)

         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = p(im,j)
           ptmp(im+1) = p(1 ,j)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo

      elseif(s(j) > 4.) then
#endif

! Packing for FFT 
           nj  = nj + 1
         jf(nj) = j

         do i=1,im
            q1(i,nj) = p(i,j)
         enddo
            q1(im+1,nj) = 0.
            q1(im+2,nj) = 0.

      endif
200   continue

      if( nj == 0) return
      call rfftmlt(q1,  q2, trigs, ifax, 1, im+2, im, nj, -1)

      do n=1,nj
         do i=5,im+2
            q1(i,n) = q1(i,n) * damp(i-2,jf(n))
         enddo
      enddo

      call rfftmlt(q1, q2, trigs, ifax, 1, im+2, im, nj, 1)

      do n=1,nj
         do i=1,im
            p(i,jf(n)) = q1(i,n)
         enddo
      enddo
!EOC
 end subroutine pft2d

 subroutine pft2d_phys(p, s, im,  jp)
  implicit none
! This is a very weak zonal filter for filtering the 
! tendencies from the physics

! !INPUT PARAMETERS:
    integer, intent(in):: im     ! Total X dimension
    integer, intent(in):: jp     ! Total Y dimension
    real, intent(in):: s(jp)     ! 3-point algebraic filter

! !INPUT/OUTPUT PARAMETERS:
    real, intent(inout):: p(im,jp)  ! Array to be polar filtered

! !LOCAL VARIABLES:
    real ptmp(0:im+1)
    real rsc, bt 
    integer  i, j

    do j=1,jp
      if(s(j) > 1.01 ) then
         rsc = 1./ min(2., s(j))
         bt  = 0.5*(min(2.,s(j))-1.)
         do i=1,im
            ptmp(i) = p(i,j)
         enddo
           ptmp(   0) = p(im,j)
           ptmp(im+1) = p(1 ,j)

         do i=1,im
            p(i,j) = rsc * ( ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)) )
         enddo
      endif
    enddo

 end subroutine pft2d_phys
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: pft_init --- Calculate algebraic and FFT polar filters
!
! !INTERFACE: 
 subroutine pft_init(im, jm, js2g0, jn1g1, sc, se, dc, de,          &
                     cosp, cose, ycrit, ipft1, ipft2)
  implicit none

! !USES:

! !INPUT PARAMETERS:
      integer im                      ! Total X dimension
      integer jm                      ! Total Y dimension
      integer js2g0                   ! j south limit ghosted 0 (SP: from 2)
      integer jn1g1                   ! j north limit ghosted 1 (starts jm)
      integer ipft1, ipft2
      real    cosp(jm)            ! cosine array
      real    cose(jm)            ! cosine array
      real    ycrit               ! critical value

! !OUTPUT PARAMETERS:
      real    sc(ipft1:ipft2)     ! Algebric filter at center
      real    se(js2g0:jn1g1)     ! Algebric filter at edge

      real    dc(im,ipft1:ipft2)  ! FFT filter at center
      real    de(im,js2g0:jn1g1)  ! FFT filter at edge

! !DESCRIPTION:
!
!   Compute coefficients for the 3-point algebraic and the FFT
!   polar filters.
!
! !REVISION HISTORY:
!
!   99.01.01   Lin          Creation
!   99.08.20   Sawyer/Lin   Changes for SPMD mode
!   01.01.30   Lin          Put into this module
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j
      real  dl, coszc, cutoff, phi, damp
      real  pi

      pi = 4.d0 * datan(1.d0)

      coszc = cos(ycrit*pi/180.)

! INIT fft polar coefficients:
      dl = pi/dble(im)
      cutoff = 1.e-20

      do j=ipft1, ipft2
         do i=1,im
            dc(i,j) = 1.
         enddo
      enddo

      do j=js2g0,jn1g1
         do i=1,im
            de(i,j) = 1.
         enddo
      enddo

!     write(6,*) 'Initializing 3-point polar filter coefficients:'

!************
! Cell center
!************
      do j=ipft1, ipft2
            sc(j) = (coszc/cosp(j))**2

#ifdef ALT_PFT
         if(sc(j) > 1.) then
#else
         if(sc(j) > 1. .and. sc(j) <= 2.0) then
            sc(j) =  1. +  (sc(j)-1.)/(sc(j)+1.)
         elseif(sc(j) > 2. .and. sc(j) <= 4.) then
            sc(j) =  1. +  sc(j)/(8.-sc(j))
            sc(j) = min(2., sc(j))
         elseif(sc(j) > 4. ) then
#endif

! FFT filter
         do i=1,im/2
            phi = dl * i
            damp = min((cosp(j)/coszc)/sin(phi), 1.)**2
            if(damp < cutoff) damp = 0.
            dc(2*i-1,j) = damp
            dc(2*i  ,j) = damp
         enddo

         endif
      enddo

!************
! Cell edges
!************
      do j=js2g0,jn1g1
            se(j) = (coszc/cose(j))**2

#ifdef ALT_PFT
         if( se(j) > 1. ) then
#else
         if(se(j) > 1. .and. se(j) <= 2.0 ) then
            se(j) =  1. +  (se(j)-1.)/(se(j)+1.)
         elseif(se(j) > 2. .and. se(j) <= 4.) then
            se(j) =  1. +  se(j)/(8.-se(j))
            se(j) = min(2., se(j))
         elseif(se(j) > 4. ) then
#endif
! FFT
            do i=1,im/2
               phi = dl * i
               damp = min((cose(j)/coszc)/sin(phi), 1.)**2
               if(damp < cutoff) damp = 0.
               de(2*i-1,j) = damp
               de(2*i  ,j) = damp
            enddo
         endif
      enddo
!EOC
 end subroutine pft_init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: rfftmlt --- Apply fast Fourier transform
!
! !INTERFACE: 
 subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isign)


! !DESCRIPTION:
!
!   Apply the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC


#ifdef SGI_FFT
!
! WS 99.11.24 : Added SGI wrappers
!
      integer inc, jump, n, lot, isign
      real   a(jump,lot)
      real   trigs(1)
      real   work(1)                       ! Not used; here for plug reason
      integer ifax(*)
! Local
      integer*4 iisign,in,iinc,ijump,ilot
      integer i, j
      real scale

!-----convert to i4
      iisign = isign
      iinc = inc
      ijump = jump
      in = n
      ilot = lot

      if( iisign < 0 ) then
!-----forward
          call dzfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)
       endif

      if( iisign > 0 ) then
!-----backward
          call zdfftm1du (iisign,in,ilot,a,iinc,ijump,trigs)

          scale = 1.0/float(n)
          do j=1,lot
             do i=1,jump
                a(i,j) = scale*a(i,j)
             enddo
          enddo
       endif
#else
!
! Default f77 version
!
      real a(jump*lot)
      real work((n+1)*lot)
      real trigs(3*n/2+1)
      integer ifax(13)
 
!     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
 
!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2
 
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
 
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
 
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
 
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
 
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
 
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign==+1) go to 30
 
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      igo=50
      if (mod(nfax,2)==1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
 
      igo=60
      go to 40
 
!     PREPROCESSING (ISIGN=+1)
 
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
 
!     COMPLEX TRANSFORM
 
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo==60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,     &
                  ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,     &
                  2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
 
      if (isign==-1) go to 130
 
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      if (mod(nfax,2)==1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
!DIR$ IVDEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
 
!     FILL IN ZEROS AT END
  110 continue
      ib=n*inc+1
!DIR$ IVDEP
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
 
!     POSTPROCESSING (ISIGN=-1):
 
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
 
  140 continue
#endif
!EOC
 end subroutine rfftmlt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fftfax --- Initialize FFT
!
! !INTERFACE: 
 subroutine fftfax (n, ifax, trigs)

! !USES:


! !DESCRIPTION:
!
!   Initialize the fast Fourier transform.  If CPP token SGI_FFT is
!   set, SGI libraries will be used.  Otherwise the Fortran code
!   is inlined.
!
! !REVISION HISTORY:
!
!   99.11.24   Sawyer       Added wrappers for SGI
!   01.03.26   Sawyer       Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef SGI_FFT
      real    trigs(1)
      integer ifax(*)
      integer n
! local
      integer*4 nn

      nn=n
      call dzfftm1dui (nn,trigs)
#else
       integer ifax(13)
       real trigs(3*n/2+1)
 
! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.
 
      data mode /3/
      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) > 5 .or. n <= 4) ifax(1) = -99
      if (ifax(1) <= 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop 'set99'
      endif
      call fftrig (trigs, n, mode)
#endif
!EOC
 end subroutine fftfax
!-----------------------------------------------------------------------
#endif

end module pft_module

! End of Lima pft_module

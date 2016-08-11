! driver for old supersource MCA and large-scale condensation routines
! Version $Id: mcm_mca_lsc.F90,v 11.0 2004/09/28 19:28:56 fms Exp $  

! TK mod 5/22/02: use mixing ratio, not spec humidity
! MS mod: use FMS values for constants rgas, grav, hl, d622

! THE FOLLOWING NOTES ARE OBSOLETE 5/22/02
! ----------------------------------------

!  TK modified to run with specific humidity as input. 
!    1) Make local copy of specific humidity input and convert to 
!         mixing ratio. 
!    2) In the precip accumulation step, make sure to
!         vertically integrate dq rather than dr.
!         (i.e., find the change in q at each level and integrate).
!    These were done for the tk1_1A test (8/15/01).

!  TK modified:
!    The mcm_mca_lsc routine outputs relative humidity (based on
!    mixing ratio) to use in cloud prediction scheme.  
!    (This replaces rh_calc functionality in the default FMS, 
!    which is based on specific humidity.)  
! ----------------------------------------


      module mcm_mca_lsc_mod

      use     fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL,   &
                             write_version_number, stdlog, close_file, &
                             open_namelist_file, check_nml_error

      use constants_mod, only: rdgas, rvgas, grav, cp_air, hlf, hlv, tfreeze

      use spectral_dynamics_mod, only: get_pk_bk, get_num_levels

      use spec_mpp_mod, only: get_grid_domain

      implicit none
      private

      public :: mcm_mca_lsc

      character(len=128) :: version = &
      '$Id: mcm_mca_lsc.F90,v 11.0 2004/09/28 19:28:56 fms Exp $'

      character(len=128) :: tagname = &
      '$Name: tikal $'

      real,allocatable,dimension(:) :: qmh,q,dqph,dq,tauran,tausno

      real,allocatable,dimension(:,:) :: psp, convq

      real,allocatable,dimension(:,:,:) :: ta,ra

!     TK mod: (analogous to ta,ra)
      real,allocatable,dimension(:,:,:) :: rh_local

      real,dimension(5102) :: etabl

      real :: rswt1,rswt2,rix

      real, parameter :: joules_per_calorie = 4.186
      real, parameter :: gm_per_kg = 1.e3
      real, parameter :: cm_per_meter = 100.

!     cpp = .24000000000000002  calories/(gram*degree K)
!     ara = .068571428571428575 calories/(gram*degree K)
      real, parameter :: cpp =    cp_air/(joules_per_calorie*gm_per_kg)
      real, parameter :: ara = rdgas/(joules_per_calorie*gm_per_kg)

!     Units of og are seconds**2/cm
      real, parameter :: og  = 1.0/(grav*cm_per_meter)

!     hlvv = 597.2288580984233  calories/gm
!     hlss = 677.01863354037266 calories/gm
      real, parameter :: hlvv =  hlv/(gm_per_kg*joules_per_calorie)
      real, parameter :: hlss = (hlv+hlf)/(gm_per_kg*joules_per_calorie)

!     d622 = .62197183098591557
      real, parameter :: d622 = rdgas/rvgas

      integer :: kx,ix,krs1,krs2,kp,km

!     nhem and ih are legacy variables from supersource.
!     I use them to provide a point of reference for those familiar with supersource.
      integer :: nhem = 1,ih=1
      logical :: do_init = .true.
      logical :: use_mixing_ratio = .false.

      namelist / mcm_mca_lsc_nml / use_mixing_ratio

      contains

!     ---------------------------------------------------------------------------------

      subroutine mcm_mca_lsc(tin,qin,p_half,tdel,qdel,rain,snow, rh_out)
      real,intent(in),dimension(:,:,:) :: tin,qin,p_half
      real,intent(out),dimension(:,:,:) :: tdel,qdel
      real,intent(out),dimension(:,:) :: rain,snow

!     TK mod:
      real,intent(out),dimension(:,:,:) :: rh_out

      integer j
      
      if (do_init) call init_mcm_mca_lsc
      do j = 1,size(p_half,2)
!        convert temperature from deg K to deg C and select slab
         ta(:,:,ih) = tin(:,j,:) - tfreeze
!        convert pressure from pascal to dynes/cm^2 and select jrow
         psp(:,ih) = p_half(:,j,size(p_half,3))*10
         
         if(use_mixing_ratio) then
           ra(:,:,ih) = qin(:,j,:)
         else
           ra(:,:,ih) = qin(:,j,:)/(1.-qin(:,j,:))
         endif
            
         call convad
         
!        calculate temperature difference (convert deg C to deg K)
         tdel(:,j,:) = ta(:,:,ih)+tfreeze-tin(:,j,:)
         if(use_mixing_ratio) then
           qdel(:,j,:) = ra(:,:,ih)-qin(:,j,:)
         else
           qdel(:,j,:) = ra(:,:,ih)/(1.+ra(:,:,ih))-qin(:,j,:)
         endif
         
!        TK mod:  Return rh_local to calling routine           
         rh_out(:,j,:) = rh_local(:,:,ih)        
         
!        use rain and snow amounts (convert cm to kg/m^2): factor of 10
!        factor of 2 comes in because tauran and tausno were multiplied by 0.5 
!        (taupre is calculated over 2dt time intervale in supersource)
         rain(:,j) = tauran(:)*10*2
         snow(:,j) = tausno(:)*10*2
      end do
      end subroutine mcm_mca_lsc
! ---------------------------------------------------------------------------------

      subroutine init_mcm_mca_lsc
      real, allocatable, dimension(:) :: pk
      integer :: is, ie, js, je, unit, io, ierr

      unit = open_namelist_file()
      ierr=1
      do while (ierr /= 0)
        read(unit, nml=mcm_mca_lsc_nml, iostat=io, end=20)
        ierr = check_nml_error (io, 'mcm_mca_lsc_nml')
      enddo
20    call close_file (unit)

      call write_version_number(version, tagname)
      if(mpp_pe() == mpp_root_pe()) write (stdlog(),nml=mcm_mca_lsc_nml)

!  ix and kx stay fixed througout run
      call get_grid_domain(is, ie, js, je)
      ix = ie - is + 1
      call get_num_levels(kx)

      if(.not.use_mixing_ratio) then
        allocate(convq(ix,kx))
      endif
      
      kp = kx + 1
      km = kx - 1
      allocate(pk(0:kx), qmh(0:kx), q(kx), dq(kx), dqph(km))
!     assume sigma levels are being used
      call get_pk_bk(pk,qmh)
      if(sum(pk) /= 0.) then
        call error_mesg('init_mcm_mca_lsc','Must be pure sigma', FATAL)
      endif
      allocate(psp(ix,1),ta(ix,kx,1),ra(ix,kx,1))
      
!     TK Mod:
      allocate(rh_local(ix,kx,1))
      
      allocate(tauran(ix),tausno(ix))
!     parmtr sets the following variables:
!     rswt1,rswt2,krs1,krs2
!     q,qmh,dq,dqph,etabl
      call parmtr
      deallocate(pk)
      do_init = .false.

      return
      end subroutine init_mcm_mca_lsc
! ---------------------------------------------------------------------------------

      subroutine convad
 
! this subroutine controls the computation of convective adjustment.
! the input is pressure, temperature and mixing ratio.
! the output are the new adjusted temps and mixing ratios.
 
      real :: taupre(ix)
      real :: tcc   (ix)
      real ::  hl   (ix)
      real :: r0    (ix,kx)
      real :: t0    (ix,kx)
      real :: temp  (ix,kx)

! to match fms diagnostics, convt and convr no longer have units of "amount"/s
! instead, they have units of "amount"
      real,dimension(ix,kx) :: convt,convr

!     TK mod:  8/15/01
!     real,dimension(ix,kx) :: convq

      integer :: ixxx,i,k
      real :: fact
 
!   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
 
! test temperatures for range between -173.16 and +102 degrees c.
! if any temperature is outside of this range, program stops.
 
      fact = 1.0
!      if (ktau .eq. 0) fact = 0.5
      do 70 k=1,kx
       do 70 i=1,ix
        temp(i,k) = 102.0 - ta(i,k,ih)
        temp(i,k) = temp(i,k) * (ta(i,k,ih)+173.16)
   70 continue
      ixxx = 0
      do 71 k=1,kx
       do 71 i=1,ix
        if (temp(i,k) .lt. 0.0) go to 771
   71 continue
      ixxx = 0
      go to 772
  771 continue
      ixxx = 1
  772 continue
      if (ixxx .ne. 0)  then
!        write(6,9030) jrow, ktau
        do 7770 k=1,kx
!          write(6,9031) k, (ta(i,k,ih),i=1,ix)
          
          write(6,*) "bad temperature: k, ta(:,k,ih)",k, &
                (ta(i,k,ih),i=1,ix)
7770    continue
!        call aquit(20)
        call error_mesg('convad in mcm_mca_lsc','bad temperature',FATAL)
      endif
 
! save initial mixing ratio and temperature.
 
      do 110 k=1,kx
        do 110 i=1,ix
          r0 (i,k) = ra (i,k,ih)
          t0 (i,k) = ta (i,k,ih)
  110 continue
 
! test to determine whether precipitation occurs as rain or snow.
! heat of condensation is adjusted accordingly.
 
! tcc(i)=+1.0 when the boundary layer temperature .gt.0 deg.c
! tcc(i)= 0.0 when the boundary layer temperature .le.0 deg.c
 
! height of freezing level is taken to be 350 meters based on
! empirical analysis by nmc.  krs1, krs2, rswt1, and rswt2 are
! determined in subroutine parmtr.
 
      do 74 i=1,ix
       tcc(i) = -rswt1*ta(i,krs1,ih) - rswt2*ta(i,krs2,ih)
   74 continue
      do 75 i=1,ix
       if (tcc(i) .lt. 0.0) then
         hl (i) = hlvv
         tcc(i) =   1.0
       else
         hl (i) = hlss
         tcc(i) =   0.0
       endif
   75 continue
 
!  call mstcnv which computes the moist adjustment.
 
!  TK Mod:
!      call mstcnv(hl)
      call mstcnv(hl, r0)
 
! compute change of temperature and mixing ratio due to condensation,
! convective and dry adjustment.
 
      do 531 k=1,kx
        do 531 i=1,ix
          convt(i,k) = ta(i,k,ih) - t0(i,k)
          convr(i,k) = ra(i,k,ih) - r0(i,k)
          
!         TK mod: 8/15/01
!           Compute the convective adjustment of water
!           in terms of specific humidity in preparation
!           for vertical integration to compute the precipitation.

          if(.not.use_mixing_ratio) then
            convq(i,k) =  ra(i,k,ih) / (ra(i,k,ih) + 1.) &
                       - (r0(i,k) / (r0(i,k) + 1.))             
          endif
          
  531 continue
 
! compute total precipitation and heating due to all adjustments
! in column.
 
! convert convt and convr to rates of change.
 
      do 532 i=1,ix
!        cmcah (i) =  dq(1) * convt(i,1) * psp(i,ih)
         if(use_mixing_ratio) then
           taupre(i) = -dq(1) * convr(i,1) ! TK mod 5/22/02
         else
           taupre(i) = -dq(1) * convq(i,1) ! TK mod: 8/15/01
         endif
  532 continue
      do 534 k=2,kx
        do 534 i=1,ix
!          cmcah (i) = cmcah (i) + dq(k) * convt(i,k) * psp(i,ih)
! TK mod: 8/15/01.  Vertically integratge the change in q:
           if(use_mixing_ratio) then
             taupre(i) = taupre(i) - dq(k) * convr(i,k) ! TK mod 5/22/02
           else
             taupre(i) = taupre(i) - dq(k) * convq(i,k)
           endif
  534 continue
      do 550 i=1,ix
        taupre(i) = taupre(i) * psp(i,ih) * og
  550 continue
!      rdt = 1.0 / float ( idt )
!      rdt = 1.0 / dt
!      do 535 k=1,kx
!        do 535 i=1,ix
!          convt(i,k) = rdt * convt(i,k)
!          convr(i,k) = rdt * convr(i,k)
!  535 continue

! set total precipitation as either rain or snow depending upon
! criterion tcc determined near beginning of subroutine.

! taupre is computed over a 2dt time interval and, therefore,
! must be divided in half.

      do 600 i=1,ix
        taupre(i) = taupre(i) * 0.5 * fact
        tauran(i) = taupre(i) *        tcc(i)
        tausno(i) = taupre(i) * (1.0 - tcc(i))
  600 continue

      return
      end subroutine convad
! ---------------------------------------------------------------------------------

      subroutine parmtr
!     ========== ======
!  this subroutine is called at the beginning of each run from atmos.
!  it calculates  constants, including those in the es table.
!  entry points:       cpoly, lgndre1, lgndre, lgndxy
!  subroutines called : o3init, gaussg

      real :: table(5102),esnew(2551),esupc(200)

      integer :: k,ii,i
      real :: qrs,esbasw,tbasw,esbasi,tbasi,tem,aa,b,c,d, &
           e,esh20,wice,wh2o

! compute sigma values at model full levels
      do 41 k=1,kx
        q(k) = 0.5*(qmh(k)+qmh(k-1))
   41 continue

!  set fundamental atmospheric constants

!  load common block /param/. for additional information,
!  see june 26,1981 notes.
      rix = 1.0/float(ix)

! Determine the levels and weights to be used for discriminating
! between rain and snow.  The 350m level is used as the critical
! level, and q=0.9575 is used as an approximation to this level
! based on the U. S. standard atmosphere.
      qrs = 0.9575
      do 71 k=1,kx-1
        if (q(k).le.qrs .and. q(k+1).gt.qrs) then
          krs1 = k
          krs2 = k+1
          rswt1 = (q(k+1)-qrs)/(q(k+1)-q(k))
          rswt2 = 1.0 - rswt1
        endif
   71 continue

!  compute level dependent constants
!  load common block /sigwgt/.  for additional information,
!  see october 19,1981 notes.
      do 10 k=1,kx
      dq    (k) = qmh(k)-qmh(k-1)
 10   continue
      do 20 k=1,km
      dqph  (k  ) = q(k+1)-q(k)
 20   continue

!  ================================================================
!  +                                                              +
!  +             construction of the es table                     +
!  +                                                              +
!  + this table is constructed from es equations from the smithsonian
!  + tables. the es input  is computed from values  (in           +
!  + one-tenth of a degree increments) of es over ice from  -153c +
!  + to 0c  and values of es over water from 0c to 102c.   output +
!  + table contains these data interleaved with their derivatives +
!  + with   respect  to  temperature  except  between -20c and 0c +
!  + where  blended  (over water  and over ice)   es  values  and +
!  + derivatives are calculated.                                  +
!  note: all es computation is done in microbars and not millibars.
!  ================================================================

!#ifdef moist
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16

!  compute es over ice between -153 c and 0 c.
!  see smithsonian meteorological tables page 350.
      do 120 i=1,1530
      tem = 120.16+0.1*float(i-1)
      aa  = -9.09718 *(tbasi/tem-1.0)
      b   = -3.56654 *alog10(tbasi/tem)
      c   =  0.876793*(1.0-tem/tbasi)
      e   = alog10(esbasi)
      table(i)=10**(aa+b+c+e)
 120  continue

!  compute es over water between -20c and freezing.
!  see smithsonian meteorological tables page 350.
      do 128 i=1,1221
      tem = 253.16+0.1*float(i-1)
      aa  = -7.90298*(tbasw/tem-1)
      b   =  5.02808*alog10(tbasw/tem)
      c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
      d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
      e   = alog10(esbasw)
      esh20  = 10**(aa+b+c+d+e)
      if (i.le.200) then
        esupc(i) = esh20
      else
        table(i+1330) = esh20
      endif
  128 continue

!  derive blended es over ice and supercooled water between -20c and 0c
      do 130 i=1,200
      tem  = 253.16+0.1*float(i-1)
      wice = 0.05*(273.16-tem)
      wh2o = 0.05*(tem-253.16)
      table(i+1330) = wice*table(i+1330)+wh2o*esupc(i)
  130 continue

!  estimate es derivative in microbars per degree celsius by
!    differencing es values in the table
      do 134 i=2,2550
      esnew(i) = 5.0*(table(i+1)-table(i-1))
  134 continue
      esnew(   1) = esnew(   2)
      esnew(2551) = esnew(2550)

!  interleave es and its derivative into one table of alternating
!    variables (es entries are odd numbered, derivatives are even)
      do 135 i=1,2551
                  ii = (i-1)*2
      table(5101-ii) = table(2552-i)
  135 continue
      do 136 i=1,2551
            ii    = (i-1)*2
      table(ii+2) = esnew(i)
  136 continue
      do 140 i=1,5102
      etabl(i)=table(i)
  140 continue
!#endif

      return

      end subroutine parmtr
! ---------------------------------------------------------------------------------

      subroutine mstcnv(hl, r0)

!     input : ta and ra
!     output: ta and ra

!      TK Mod:  additional input:  r0
!      TK Mod:  additional output: rh_local

! this subroutine computes precipitation due to both large scale
! condensation and small scale moist convection.

      integer,parameter :: iters = 1
      integer :: iter

      real,intent(in),dimension(:) ::     hl
      real,dimension(size(hl,1)) :: tauran,tausno
      real,dimension(ix,kx) ::  wk, alrm
      real,dimension(ix,kx) :: c
      real,dimension(ix,kx) :: rsmin

      real,dimension(ix,kx) :: es,des,psq
      real,dimension(ix,kx) :: rs
      integer,dimension(ix,kx) :: mvf
      real,dimension(ix,kx) :: psqm
      integer,dimension(ix,kx) :: it
      real, dimension(ix,kx) :: temnew,ratnew,temavg
      integer,dimension(ix,kx) :: itest,iwk
      real,dimension(ix,kx) :: term2
      real,dimension(ix,kx) ::  wt1, wt3,delta
      real,dimension(ix,kx) :: sum0,term1,test
      real,dimension(ix,kp) :: sum2, sum1
      real :: delmin,esmin
      integer :: i,k,iflag1,iflag2
      
!  TK mod:
      real,dimension(ix,kx) :: r0
      
      real :: hc=1.00

! initialize dummy arrays to zero for moist adiabatic adjustment
! computation.

      do 72 i=1,ix
       sum0(i,kx) = 0.0
       sum1(i, 1) = 0.0
       sum2(i, 1) = 0.0
   72 continue

! compute pressure at middle of layers for use in moist adiabatic
! lapse rate evaluation.

      do 130 k=1,km
        do 130 i=1,ix
          psqm(i,k) = psp(i,ih) * qmh(k)
  130 continue

! compute pressure at individual levels.

      do 140 k=1,kx
        do 140 i=1,ix
          psq(i,k) = psp(i,ih) * q(k)
  140 continue

! beginning of code for determining precipitation by both large
! and small scale processes.  compute es, c, and rs at the
! individual levels which are
!   es = saturation vapor pressure.
!    c = derivative of saturation mixing ratio with respect to
!        temperature.
!   rs = saturation mixing ratio.
!   rsmin = minimum saturation mixing ratio on es table
!   esmin = minimum saturation vapor pressure on es table
!   delmin = minimum p-es allowed, this is done to keep p-es
!            positive.

! note; repeated iteration through the column is no longer performed
! (i.e. iters=1).

      delmin = 1.0e-10
      esmin  = etabl(1)
      do 519 k=1,kx
       do 519 i=1,ix
!       rsmin(i,k) = .622*esmin / (psq(i,k)-esmin) ! MS mod 5/22/02
        rsmin(i,k) = d622*esmin / (psq(i,k)-esmin) ! MS mod 5/22/02
  519 continue

      do 520 iter=1,iters

! es values are taken from a table(120k-375k) at .1 degree intervals

        do 149 k=1,kx
          do 149 i=1,ix
            it(i,k) = 10.0 * (ta(i,k,ih) + 153.05)
  149   continue
        do 80 k=1,kx
         do 80 i=1,ix
          it(i,k) = 2*it(i,k) + 1
          it(i,k) = max(it(i,k), 1)
          it(i,k) = min(it(i,k), 5102-3)
   80   continue
        do 81 k=1,kx
         do 81 i=1,ix
           es(i,k) = etabl(it(i,k)  )
          des(i,k) = etabl(it(i,k)+1)
   81   continue
        do 150 k=1,kx
          do 150 i=1,ix
!cc         wk(i,k)    = .622 * hc *  es(i,k) ! MS mod 5/22/02
            wk(i,k)    = d622 * hc *  es(i,k) ! MS mod 5/22/02

!cc         des(i,k)   = .622 * hc * des(i,k) ! MS mod 5/22/02
            des(i,k)   = d622 * hc * des(i,k) ! MS mod 5/22/02

            delta(i,k) = psq(i,k) - hc *  es(i,k)
  150   continue

!  make sure delta, (press - es press) does not become negative
        do 82 k=1,kx
         do 82 i=1,ix
          delta(i,k) = max(delta(i,k), delmin)
   82   continue
        do 151 k=1,kx
          do 151 i=1,ix
             c(i,k) = psq(i,k) * des(i,k) / (delta(i,k) * delta(i,k))
            rs(i,k) = wk(i,k) / delta(i,k)
  151   continue
        do 83 k=1,kx
         do 83 i=1,ix
          rs(i,k) = max(rs(i,k), rsmin(i,k) )
          
!         TK Mod:  Determine relative humidity for cloud determination 
!                 (later in rh_clouds).  Note that this is within an interation 
!                 loop where temperature could change, but the number of
!                 iterations is currently set to 1.

          rh_local(i,k,ih) = r0(i,k) / rs(i,k)
          
   83   continue

! compute mean temperature of layer for use in moist adiabatic
! lapse rate evaluation.

        do 160 k=1,km
          do 160 i=1,ix
            temavg(i,k) = (ta(i,k,ih) + ta(i,k+1,ih)) * 0.5
  160   continue

! compute es and des (derivative of saturation vapor pressure)
! at the mid-point of the layer for use in moist adiabatic lapse
! rate evaluation.

        do 169 k=1,km
          do 169 i=1,ix
            it(i,k) = 10.0 * (temavg(i,k) + 153.05)
  169   continue
        do 85 k=1,kx
         do 85 i=1,ix
          it(i,k) = 2*it(i,k) + 1
          it(i,k) = max(it(i,k), 1)
          it(i,k) = min(it(i,k), 5102-3)
   85   continue
        do 86 k=1,kx
         do 86 i=1,ix
           es(i,k) = etabl(it(i,k)  )
          des(i,k) = etabl(it(i,k)+1)
   86   continue
        do 170 k=1,km
          do 170 i=1,ix
!cc         es (i,k) = .622 * hc * es (i,k) ! MS mod 5/22/02
            es (i,k) = d622 * hc * es (i,k) ! MS mod 5/22/02

!cc         des(i,k) = .622 * hc * des(i,k) ! MS mod 5/22/02
            des(i,k) = d622 * hc * des(i,k) ! MS mod 5/22/02

  170   continue

! compute alrm (moist adiabatic lapse rate) for each layer.

        do 171 k=1,km
          do 171 i=1,ix
            alrm(i,k) = (dqph(k) * (ara * psqm(i,k) * &
                        (temavg(i,k) + tfreeze) + hl(i) * es(i,k))) &
                      / (qmh(k) * (psqm(i,k) * cpp &
                                     + hl(i) * des(i,k)))
  171   continue

! tests are now performed to determine occurrence and type of
! precipitation.  for small scale, moist convection, layer must be
! both super-moist adiabatic and supersaturated.  for large scale
! condensation, just the individual layer must be supersaturated.

! here begins the computation of the integer arrays iwk and itest
! which will determine the levels where the moist convective and
! large scale adjustments occur, respectively.  mvf is the same as
! iwk except it applies to the layers rather than the levels
! involved.

        if (iter .le. 1) then

!         test(i,k) = 1.0 - ra(i,k,ih) / rs(i,k)

! set test to negative number for supersaturation.

          do 88 k=1,kx
           do 88 i=1,ix
            test(i,k) = wk(i,k) - ra(i,k,ih)*delta(i,k)
   88     continue
          do 89 k=1,kx
           do 89 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   89     continue

          do 181 i=1,ix
            mvf(i,kx) = 0
  181     continue

          do 90 k=1,km
           do 90 i=1,ix
            mvf(i,k) = itest(i,k) * itest(i,k+1)
   90     continue

          do 91 k=1,km
           do 91 i=1,ix
            test(i,k) = alrm(i,k) + 0.01 + ta(i,k,ih) - ta(i,k+1,ih)
   91     continue

          do 92 k=1,km
           do 92 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   92     continue

          do 93 k=1,km
           do 93 i=1,ix
            mvf(i,k) = mvf(i,k) * itest(i,k)
   93     continue

!     mvf(i,k)=1 where moist adiabatic adjustment is required
!     mvf(i,k)=0 where either condensation or nothing is required

          do 94 i=1,ix
           iwk(i,1) = mvf(i,1)
   94     continue

          do 95 k=2,kx
           do 95 i=1,ix
            iwk(i,k) = max(mvf(i,k),mvf(i,k-1))
   95     continue

          do 96 k=1,kx
           do 96 i=1,ix
            test(i,k) = rs(i,k) - ra(i,k,ih)
   96     continue
          do 97 k=1,kx
           do 97 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   97     continue

          do 98 k=1,kx
           do 98 i=1,ix
            itest(i,k) = itest(i,k) * (1-iwk(i,k))
   98     continue

          do 99 k=1,kx
           do 99 i=1,ix
            wt1(i,k) = mvf(i,k)
            wt3(i,k) = iwk(i,k)
   99     continue

! iwk now contains the convection flagbits and itest the
! condensation flagbits: iwk=1 if and only if convection.
! itest=1 if and only if condensation, otherwise they equal zero.
! the value of zero implies no adjustment.

! solution of simultaneous equations for the new temperatures and
! mixing ratios as a result of both condensation and convective
! adjustments.  here
!       ta = temperature after dry adjustment.
!       ra = old mixing ratio
!   temnew = new temperature
!   ratnew = new mixing ratio

          iflag1 = 0
          iflag2 = 0
          do 188 k=1,kx
           do 188 i=1,ix
            iflag1 = iflag1 + itest(i,k)
            iflag2 = iflag2 +   mvf(i,k)
  188     continue
          if (iflag1 .eq. 0 .and. iflag2 .eq. 0) go to 525
          if (iflag1 .eq. 0) then
            do 3376 k=1,kx
              do 3376 i=1,ix
                temnew(i,k) = ta(i,k,ih)
                ratnew(i,k) = ra(i,k,ih)
 3376       continue
          else

! large scale condensation adjustment.

            do 3333 k=1,kx
              do 3333 i=1,ix
                test(i,k) = itest(i,k)
 3333       continue

! **    itest now has a "1" where condensation is allowed to occur,
!        i.e., where moist adjustment did not occur and where the
!        humidity (ratnew/rs) is greater than 1.0, and "0" otherwise.

            do 3004 k=1,kx
              do 3004 i=1,ix
                temnew(i,k) = ta(i,k,ih) + &
                            test(i,k) * (hl(i) * (ra(i,k,ih) - rs(i,k)) &
                              / (cpp + hl(i) * c(i,k)))
                ratnew(i,k) = test(i,k) * (rs(i,k) + c(i,k) &
                                           * (temnew(i,k) - ta(i,k,ih)) &
                                           - ra(i,k,ih)) + ra(i,k,ih)
 3004       continue
          endif

        endif

        if (iter.gt.1) then
          do 3006 k=1,kx
            do 3006 i=1,ix
              temnew(i,k) = ta(i,k,ih) &
                          + test(i,k) * (hl(i) * (ra(i,k,ih) - rs(i,k)) &
                            / (cpp + hl(i) * c(i,k)))
              ratnew(i,k) = test(i,k) * (rs(i,k) + c(i,k) &
                                         * (temnew(i,k) - ta(i,k,ih)) &
                                         - ra(i,k,ih)) + ra(i,k,ih)
 3006     continue
        endif

! small scale moist convective adjustment.

        if (iflag2 .ne. 0) then
          do 700 k=1,km
            do 700 i=1,ix
              sum0(i,kx-k) = wt1(i,kx-k) * (sum0(i,kx+1-k) + alrm(i,kx-k))
  700     continue
          do 701 k=1,kx
            do 701 i=1,ix
              term2(i,k) = dq(k) * (cpp + hl(i) * c(i,k))
              term1(i,k) = term2(i,k) * (ta(i,k,ih) + sum0(i,k)) &
                        + dq(k) * hl(i) * (ra(i,k,ih) - rs(i,k))
  701     continue
          do 702 k=1,kx
            do 706 i=1,ix
              sum1(i,k+1) = wt3(i,k) * (sum1(i,k) + term1(i,k))
              sum2(i,k+1) = wt3(i,k) * (sum2(i,k) + term2(i,k))
  706       continue
            do 172 i=1,ix
             if (wt1(i,k) .ne. wt3(i,k)) then
               temnew(i,k  ) = sum1(i,k+1) / sum2(i,k+1)
                 sum1(i,k+1) = 0.0
                 sum2(i,k+1) = 0.0
             endif
  172       continue
  702     continue
          do 703 k=1,km
            do 703 i=1,ix
              temnew(i,kx-k) = (1.0 - wt1(i,kx-k)) * temnew(i,kx  -k) &
                             +        wt1(i,kx-k)  * temnew(i,kx+1-k)
  703     continue
          do 704 k=1,km
            do 704 i=1,ix
              temnew(i,k) = temnew(i,k) - sum0(i,k)
  704     continue
          do 705 k=1,kx
            do 705 i=1,ix
              ratnew(i,k) = (rs(i,k) + c(i,k) &
                             * (temnew(i,k) - ta(i,k,ih))) * wt3(i,k) &
                          + (1.0 - wt3(i,k)) * ratnew(i,k)
  705     continue
        endif

! end of condensation and convective adjustment algorithm.  replace
! previous temperatures and mixing ratios with new ones.

        do 510 k=1,kx
          do 510 i=1,ix
            ta(i,k,ih) = temnew(i,k)
            ra(i,k,ih) = ratnew(i,k)
  510   continue
  520 continue
  525 continue

      return
      end subroutine mstcnv
! ---------------------------------------------------------------------------------

      end module mcm_mca_lsc_mod

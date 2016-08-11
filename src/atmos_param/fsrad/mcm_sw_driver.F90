      MODULE MCM_SW_DRIVER_MOD
!
!
!   TK modified from Kerr code 
!     /t90/tk/climatemodel/r30/atm_leg1/atm5a/atmlib_mods/no_ifdefs
!
!   Added interface routine (mcm_shortwave_driver) which is called by
!     fsrad and which calls mcm_swnew in this
!     module after constructing the appropriate inputs.
!
!   TK mod (6/01/01) added test to make sure CosZ is greater than
!     zero, otherwise skip sw calculation and set output to zero
!     for that point...
!!
!
!

!      USE   Constants_Mod, ONLY: grav, tfreeze

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

      use mcm_swnew_mod, only: mcm_swnew

implicit none
      private

      integer :: kx, kp
      character(len=128) :: version = '$Id: mcm_sw_driver.F90,v 10.0 2003/10/24 22:00:31 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical :: module_is_initialized = .false.

      public :: mcm_shortwave_driver, mcm_sw_driver_init, &
                mcm_sw_driver_end

!     -------------------------------------------------
! TK NOTE: not ready for this yet...      implicit none




contains

      subroutine mcm_shortwave_driver(                                 &
                     Nclds, KtopSW, KbtmSW, Press, Rh2o, Qo3, CldAmt, &
                     CUVRF, CIRRF, CIRAB, Rco2, CosZ, SSolar,         &
                     Albedo, FSW, DFSW, UFSW, TdtSW, Phalf)


!     TK Original calling parameters in Kerr code:
!     &     ccosz, ttauda, rco2, solc, p, pp, qmh,
!     &     etabl, gvbrps, gqmix, grad, slwdp, rclim, dduo3n, r, 
!     &     flx, heat,
!     &     radian, radang,
!     &     zsnit, zsnitc, zsng, zsngc, zplalb, zplalbc, icf,
!
!     &     q, calb, iflg, ncv,
!
!     &     ipower, kflags,
!
!     &     tfra, alb, ts, zin, ca, oas,
!     &     snwdpt,
!     &     pkice,
!
!     &     kp, kx, ng, lx, nc, il, ix, jrow, iy, ls, ih, nhem)


      integer, intent (in), dimension(:,:)     :: Nclds
      integer, intent (in), dimension(:,:,:)      :: KtopSW, KbtmSW
      real, intent (in)   , dimension(:,:,:)      :: Press, Phalf
      real, intent (in)   , dimension(:,:,:)      :: CldAmt, CUVRF,&
                                                  &  CIRRF, CIRAB
      real, intent (in)   , dimension(:,:,:)      :: Rh2o, Qo3

      real, intent (in)                           :: Rco2
      real, intent (in)   , dimension(:,:)        :: CosZ
      real, intent (in)   , dimension(:,:)        :: SSolar
      real, intent (in)   , dimension(:,:)        :: Albedo

      REAL,   INTENT(OUT), DIMENSION(:,:,:)       :: FSW, DFSW, UFSW
      REAL,   INTENT(OUT), DIMENSION(:,:,:)       :: TdtSW


!----------------LOCAL ARRAY STORAGE------------------------------------

      INTEGER, DIMENSION(kx+2)   :: kthsw, kbhsw
      REAL,    DIMENSION(kx+2)   :: cwca, cwcb, coca, cloudy
      REAL,    DIMENSION(kx+1)   :: flx, heat, UF, DF
      REAL                       :: grdflx
      real,   dimension(0:kx)   :: pp
      real,   dimension(1:kx+1) :: pr2
      real,   dimension(size(press,1), size(press,2), kx) :: sigma_level

      integer ncv
      real pchg

      INTEGER :: ix, jx
      integer :: ipt, jrow, k

      if(.not.module_is_initialized) then
        call error_mesg('mcm_shortwave_driver', &
                        'shortwave_driver_mod is not initialized',FATAL)
      endif

      ix = SIZE(Press,1)
      jx = SIZE(Press,2)

!     Compute sigma levels, which are needed for an mcm_swnew routine input:
      do k=1,kx
        sigma_level(:,:,k) = Press(:,:,k) / Phalf(:,:,kp)
      enddo

!     Loop over all points in the horizontal array, computing
!       radiation for one column at a time in mcm_swnew...

      do jrow = 1, jx

         do ipt = 1, ix

!           Create necessary column input arrays for mcm_swnew routine:

            do k = 1, kp

               kthsw (k) = KtopSW (ipt,jrow,k)
               kbhsw (k) = KbtmSW (ipt,jrow,k)

               cwca (k) = CIRRF (ipt,jrow,k)
               cwcb (k) = CIRAB (ipt,jrow,k)
               coca (k) = CUVRF (ipt,jrow,k)
               cloudy(k) = CldAmt (ipt,jrow,k)

            end do

            ncv = Nclds(ipt,jrow) + 2

!           Load surface albedo into 2 input arrays in location ncv:
            cwca (ncv) = Albedo (ipt,jrow)           
            coca (ncv) = Albedo (ipt,jrow) 

!           Set ncv element of cloudy to 1.0 to match supersource:
            cloudy(ncv) = 1.0

!           Also initialize ncv+1 element of the cloud arrays (not
!            sure if this is actually done in SS but trying to
!            match the input arrays for the sample point -- TK):

            cwca (ncv+1) = Albedo (ipt,jrow)           
            coca (ncv+1) = Albedo (ipt,jrow) 
            cloudy(ncv+1) = 1.0
            cwcb(ncv+1) = 0.0

!          Initialize element 1 of the following arrays exactly as
!           in supersource:

            kbhsw (  1) = 1
            cwca  (  1) = 0.0
            cwcb  (  1) = 0.0
            coca  (  1) = 0.0
            cloudy(  1) = 0.0


!           Compute pp, pr2 input arrays for mcm_swnew:

            do k = 0, kx
               pp(k) = Phalf(ipt,jrow,k+1) * 10
            end do

            pchg = sqrt (Phalf(ipt,jrow,kp)/101325.0)
            do k = 1, kx
               pr2(k) = sqrt(sigma_level(ipt,jrow,k)) * pchg
            end do  

!!           TK Print diagnostics for auxilary inputs:
!            if ((ipt .eq. 49) .and. (jrow .eq. 27)) then
!               print *
!               print *, 'TK chkpt 1 in mcm_sw_mod.F.i=49,j= ', jrow
!               print *
!               print *, ' Cosz (FMS computed) = ', Cosz(ipt,jrow)
!
!               print *, ' rco2 = ', rco2
!               print *, ' rh2o = ', Rh2o(ipt,jrow,:)
!               print *, ' ro3 = ', Qo3(ipt,jrow,:)
!               print *, ' pp[0:kx] = ', pp
!               print *, ' Ssolar (calculated) = ', Ssolar(ipt,jrow)
!               print *, ' pr2[1:kx+1] = ', pr2
!               print *, ' kx = ', kx
!               print *
!               print *, ' Cloud data for shortwave: '
!               print *, '      (use only k = 2 to nc1 for clouds)'
!               print *, '      (use k = ncv for surface albedo)'
!               print *, ' ncv = ', ncv  
!               print *, ' k cloudy kthsw kbhsw      cwca ',       &
!     &             '             cwcb               coca    '
!               print *
!               do 325 k = 1, kx
!                  write(6,323) k, cloudy(k), kthsw(k), kbhsw(k),  &
!     &               cwca(k), cwcb(k), coca(k)
!323               format (i3, 2x, f4.2, 2x, i3, 2x, i3, 3(2x, e15.6))
!325            continue
!               print *
!               do k = 0, kx
!                  print *, ' k, pp(k) = ', k, pp(k)
!               end do
!               do k = 1, kx
!                  print *, ' k, pr2(k) = ', k, pr2(k)
!               end do
!               print *
!            end if

!         Only call shortwave routine if CosZ is greater than 0.0

            if (CosZ(ipt,jrow) > 0.0) then       

               call mcm_swnew (CosZ(ipt,jrow),                &
                 rco2, Rh2o(ipt,jrow,:), Qo3(ipt,jrow,:), pp, &
                 cwca, cwcb, coca, cloudy, kthsw, kbhsw,      &
                 Ssolar(ipt,jrow), pr2,                       &
                 flx, heat, grdflx, ncv, kx, UF, DF)

!              Fill the appropriate output arrays with the column output:

               do k = 1, kp
                  FSW(ipt,jrow,k) = flx(k)
                  DFSW(ipt,jrow,k) = DF(k)
                  UFSW(ipt,jrow,k) = UF(k)
               end do

               do k = 1, kx
                  TdtSw(ipt,jrow,k) = heat(k)
               end do

            else
!              CosZ is  0.0, so zero sw output arrays at that point:
                  
               do k = 1, kp
                  FSW(ipt,jrow,k) = 0.
                  DFSW(ipt,jrow,k) = 0.
                  UFSW(ipt,jrow,k) = 0.
               end do

               do k = 1, kx
                  TdtSw(ipt,jrow,k) = 0.
               end do

            end if

!!           TK Print diagnostics:
!            if ((ipt .eq. 49) .and. (jrow .eq. 27)) then
!               print *
!               print *, ' mcm_swnew output arrays, etc: '
!               print *
!               print *, ' k      flx (ergs/cm2/sec)      heat (K/day) '
!               do 328 k = 1, kx+1
!                  write(6,327) k, flx(k), heat(k)
!327               format (i3, 7x, f14.6, 3x, f15.9)
!328            continue
!               print *
!               print *, ' grdflx [ergs/cm2/sec] = ', grdflx
!               print *
!            end if

         end do
      end do

      end subroutine mcm_shortwave_driver
! ---------------------------------------------------------------------------------------
      subroutine mcm_sw_driver_init(kx_in)
      integer, intent(in) :: kx_in

      kx = kx_in
      kp = kx + 1
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

      return
      end subroutine mcm_sw_driver_init
! ---------------------------------------------------------------------------------------

      subroutine mcm_sw_driver_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      end subroutine mcm_sw_driver_end
! ---------------------------------------------------------------------------------------
      end module MCM_SW_DRIVER_MOD

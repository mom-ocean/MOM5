      module mcm_lw_mod

!   TK modified from ajb code 
!     /net/ajb/radiation_code/updates/v4.3/lw_mod.F 

!   Added interface routine (lw_rad_ss) which is called by
!     fsrad and which calls lwcool in this
!     module after constructing the appropriate inputs.

      USE   Constants_Mod, ONLY: grav, tfreeze

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

      implicit none
      private

      integer, parameter :: nb_lw=19, ng=3, ngp=ng+1
      integer :: ix, jx, kx, kp, km
!------------ VERSION NUMBER ----------------

      character(len=128) :: version = '$Id: mcm_lw.F90,v 13.0 2006/03/28 21:09:36 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical :: module_is_initialized = .false.

      public :: MCM_LW_RAD, mcm_lw_init, mcm_lw_end

!     -------------------------------------------------
! TK NOTE: not ready for this yet...      implicit none

!-----------------------------------------------------------------------
!--------------------- G L O B A L   D A T A ---------------------------
!-----------------------------------------------------------------------

! #include "parm.h"

! TK note: (ix,jx are set up below -- depend on domain decomposition)
!      integer, parameter :: ix = 96
!      integer, parameter :: jx = 80

! TK note: (kx,kp,km defined in rdparm)
!      integer, parameter :: kx = 14
!      integer, parameter :: kp = kx + 1
!      integer, parameter :: km = kx - 1

! TK note: (nb,ng,ngp are defined in rdparm)
!      integer, parameter :: nb=19
!      integer, parameter :: ng=3
!      integer, parameter :: ngp=ng +1

!   nb=number of spectral bands
!   ng=number of absorbing gases

!    TK: This is a name change to avoid having to change nb locally:
      INTEGER, PARAMETER :: nb = nb_lw

      INTEGER :: LMAX, i

      real :: grav_accel
      real :: pi_alpha

      real :: co2_mixrat
      real :: t_freeze

      real :: h2o_line_width(nb)
      real :: h2o_line_strength(nb)
      real :: h2o_corr_a(nb,2)
      real :: h2o_corr_b(nb,2)
      real :: absorp_table(30,12)
      real :: d_absorp_dt(7,12)

      data h2o_line_width/ &
         .93000e-01,  .18200e-00,  .94000e-01,  .79700e-01,  .73300e-01, &
         .52000e-01,  .67000e-01,  .45900e-01,  .10000e+01,  .10000e+01, &
         .89000e-01,  .23000e-00,  .32000e-00,  .29600e-00,  .45200e-00, &
         .35900e+00,  .16500e+00,  .10400e+00,  .11600e+00 /

      data h2o_line_strength/ &
         .57975e+03,  .72103e+04,  .60248e+04,  .14039e+04,  .78952e+02, &
         .42700e+01,  .71500e-01,  .27500e-01,  .00000e+00,  .00000e+00, &
         .12650e+02,  .13440e+03,  .63290e+03,  .33120e+03,  .43410e+03, &
         .13600e+03,  .35650e+02,  .90150e+01,  .15290e+01 /

      data (h2o_corr_a(i,1),i=1,19)/ &
        -.87600e-02, -.26800e-02,  .20300e-02,  .96100e-02,  .15820e-01, &
         .01705e+00,  .02620e+00,  .02918e+00,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0         /
      data (h2o_corr_a(i,2),i=1,19)/ &
        -.67500e-02, -.29300e-02,  .14300e-02,  .98400e-02,  .13710e-01, &
         .01579e+00,  .02410e+00,  .02596e+00,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0         /

      data (h2o_corr_b(i,1),i=1,19)/ &
         .14150e-04,  .15700e-05, -.10300e-04, -.43140e-04, -.37440e-04, &
        -.51440e-04, -.74100e-04, -.80760e-04,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0         /
      data (h2o_corr_b(i,2),i=1,19)/ &
         .85500e-05,  .20100e-05, -.13000e-04, -.40810e-04, -.16150e-04, &
        -.44510e-04, -.40300e-04, -.66720e-04,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0        ,  .0        , &
         .0        ,  .0        ,  .0        ,  .0         /

      data (absorp_table(i,12),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0004,.0006,.0011,.0020,.0034,.0058, &
      .0097,.0156,.0248,.0386,.0586,.0861,.1205,.1607,.2040,.2480, &
      .2918,.3362,.3812,.4250,.4665,.5045,.5386,.5707,.6037,.6395/
      data (absorp_table(i,11),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0004,.0006,.0011,.0020,.0034,.0057, &
      .0093,.0149,.0231,.0351,.0520,.0744,.1032,.1388,.1803,.2255, &
      .2715,.3157,.3580,.3993,.4404,.4812,.5200,.5558,.5899,.6248/
      data (absorp_table(i,10),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0004,.0006,.0011,.0019,.0033,.0054, &
      .0086,.0131,.0194,.0280,.0398,.0557,.0771,.1047,.1385,.1777, &
      .2206,.2648,.3080,.3497,.3910,.4328,.4748,.5155,.5535,.5886/
      data (absorp_table(i,9),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0004,.0006,.0011,.0019,.0031,.0049, &
      .0075,.0110,.0157,.0221,.0309,.0430,.0591,.0798,.1056,.1366, &
      .1728,.2133,.2565,.3004,.3436,.3856,.4271,.4680,.5079,.5456/
      data (absorp_table(i,8),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0004,.0006,.0010,.0017,.0027,.0041, &
      .0060,.0086,.0122,.0170,.0236,.0323,.0437,.0585,.0772,.1005, &
      .1287,.1622,.2005,.2426,.2865,.3301,.3721,.4127,.4525,.4915/
      data (absorp_table(i,7),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0003,.0006,.0009,.0015,.0023,.0034, &
      .0049,.0070,.0099,.0138,.0189,.0256,.0342,.0454,.0597,.0778, &
      .1003,.1278,.1603,.1977,.2390,.2824,.3257,.3673,.4073,.4466/
      data (absorp_table(i,6),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0003,.0005,.0008,.0013,.0019,.0027, &
      .0039,.0056,.0079,.0109,.0149,.0200,.0265,.0349,.0457,.0595, &
      .0769,.0986,.1250,.1565,.1928,.2331,.2758,.3187,.3603,.4003/
      data (absorp_table(i,5),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0003,.0005,.0007,.0010,.0014,.0021, &
      .0029,.0042,.0058,.0079,.0107,.0143,.0189,.0247,.0321,.0415, &
      .0535,.0688,.0879,.1113,.1395,.1727,.2104,.2515,.2943,.3366/
      data (absorp_table(i,4),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0003,.0004,.0006,.0008,.0011,.0015, &
      .0021,.0029,.0039,.0052,.0069,.0091,.0118,.0152,.0194,.0248, &
      .0316,.0402,.0510,.0647,.0819,.1031,.1290,.1598,.1954,.2351/
      data (absorp_table(i,3),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0003,.0004,.0005,.0007,.0010,.0013, &
      .0018,.0024,.0032,.0042,.0054,.0069,.0087,.0109,.0137,.0171, &
      .0214,.0267,.0333,.0416,.0520,.0652,.0818,.1025,.1277,.1579/
      data (absorp_table(i,2),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0002,.0004,.0005,.0007,.0009,.0013, &
      .0017,.0022,.0029,.0037,.0046,.0058,.0071,.0086,.0105,.0127, &
      .0154,.0186,.0225,.0273,.0331,.0404,.0496,.0613,.0762,.0949/
      data (absorp_table(i,1),i=1,30)/ &
      .0000,.0001,.0001,.0002,.0002,.0004,.0005,.0007,.0009,.0013, &
      .0017,.0022,.0028,.0036,.0044,.0054,.0066,.0079,.0094,.0111, &
      .0131,.0154,.0180,.0212,.0249,.0295,.0351,.0421,.0511,.0627/

      data d_absorp_dt/ &
      0.,.0006,.0017,.0036,.0061,.0092,.0158,0.,.0006,.0017,.0035,.0065, &
      .0114,.0230,0.,.0005,.0016,.0038,.0081,.0169,.0376,0.,.0005,.0017, &
      .0045,.0111,.0256,.0568,0.,.0004,.0021,.0065,.0178,.0421,.0864,    &
      0.,.0004,.0011,.0088,.0248,.0579,.1046,0.,.0004,.0032,.0112,.0316, &
      .0725,.1160,0.,.0004,.0036,.0141,.0392,.0880,.1256,0.,.0004,.0038, &
      .0187,.0508,.1045,.1345,0.,.0003,.0036,.0220,.0650,.1138,.1380,0., &
      .0002,.0028,.0234,.0859,.1180,.1410,0.,.0001,.0022,.0257,.0936,    &
      .1183,.1444/

      contains

!#######################################################################
      subroutine mcm_lw_init(ix_in, jx_in, kx_in)
      integer, intent(in) :: ix_in, jx_in, kx_in

      ix = ix_in
      jx = jx_in
      kx = kx_in
      kp = kx + 1
      km = kx -1

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

      return
      end subroutine mcm_lw_init
!#######################################################################

subroutine mcm_lw_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine mcm_lw_end

!#######################################################################

      SUBROUTINE MCM_LW_RAD (KTOP,KBTM,NCLDS,EMCLD, &
                      PRES,TEMP,RH2O,QO3,CAMT, &
                      RRVCO2,  HEATRA,GRNFLX,TOPFLX, phalf)

!   This illustrates the intended sizes of arrays:
!      INTEGER, INTENT(IN), DIMENSION(ix,jx,kp) :: KTOP,KBTM
!      INTEGER, INTENT(IN), DIMENSION(ix,jx)    :: NCLDS
!      REAL,    INTENT(IN), DIMENSION(ix,jx,kp) :: EMCLD
!      REAL,    INTENT(IN), DIMENSION(ix,jx,kp) :: PRES,TEMP
!      REAL,    INTENT(IN), DIMENSION(ix,jx,kx) :: RH2O,QO3
!      REAL,    INTENT(IN), DIMENSION(ix,jx,kp) :: CAMT
!      REAL,    INTENT(IN)                      :: RRVCO2

!      REAL,   INTENT(OUT), DIMENSION(ix,jx,kx) :: HEATRA
!      REAL,   INTENT(OUT), DIMENSION(ix,jx)    :: GRNFLX,TOPFLX

!     TK mod:
!      REAL,    INTENT(IN), DIMENSION(ix,jx,kp) :: phalf

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)    :: NCLDS
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRES,TEMP
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: CAMT
      REAL,    INTENT(IN)                      :: RRVCO2
 
      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: HEATRA
      REAL,   INTENT(OUT), DIMENSION(:,:)    :: GRNFLX,TOPFLX

!     TK mod:
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: phalf

!----------------LOCAL ARRAY STORAGE------------------------------------
      real, dimension(SIZE(PRES,1)) :: dummy
      REAL, DIMENSION(SIZE(PRES,1),SIZE(PRES,1),0:kp) :: cloud_cover

      REAL, DIMENSION(ix,jx,1:kx) :: sigma_level
      REAL, DIMENSION(ix,jx,0:kx) :: sigma_half_level
      REAL, DIMENSION(ix,jx,1:kx) :: sigma_thick
      REAL, DIMENSION(kx) :: o3_mixrat
 
      INTEGER :: i,j, n, ipr

!     TK mods:
      integer :: iindex, jindex, klev, kpr

      real :: pi

      if(.not.module_is_initialized) then
        call error_mesg('MCM_LW_RAD','module is not initialized.',FATAL)
      endif

!     ix = SIZE(PRES,1)
!     jx = SIZE(PRES,2)
      lmax = ix * (kx + 1)

!     Gather/calculate other needed parameters for ss longwave code:

!     Convert grav to grav_accel in cgs units:
      grav_accel = grav * 100.

      pi = 4. * atan(1.)
      pi_alpha = 0.28 * pi

      t_freeze = tfreeze

!     Multiply CO2 volume mixing ratio by 1.5194 to convert to mass mixing
!     ratio... TK -- need to check accuracy of this...

      co2_mixrat = RRVCO2 * 1.5194

!     Compute sigma levels, which are needed by ss routine:
!     This can be done on any x,y point in the domain.
      do klev = 1, kx
        do j=1,jx
          do i=1,ix
            sigma_level(i,j,klev) = pres(i,j,klev) / phalf(i,j,kp)
          enddo
        enddo
      enddo

!     Compute half sigma levels, which are needed by ss routine:
      do klev = 1, kp
        do j=1,jx
         do i=1,ix
          sigma_half_level(i,j,klev-1) = phalf(i,j,klev) / phalf(i,j,kp)
         enddo
        enddo
      enddo

!     Compute delta sigma of layers, which are needed by ss routine:
      do klev = 1, kx
        do j=1,jx
          do i=1,ix
            sigma_thick(i,j,klev) = sigma_half_level(i,j,klev) - &
                                    sigma_half_level(i,j,klev-1) 
          enddo
        enddo
      enddo

!     Create cloud_cover input field for ss longwave code:

!     Expand_cloud routine creates a Manabe Climate Model cloud_cover
!     field from the FMS arrays: nclds, ktop, kbtm, camt, and emcld.

!     This is patterned after the expand_cloud subroutine in
!      the clouds.f90 module.  See also the supersource code
!      in impt.f 
!     Accounts for cloud emmissivity...

      cloud_cover = 0.0
      do jindex=1,size(nclds,2)
      do iindex=1,size(nclds,1)
         do n=2,nclds(iindex,jindex)+1
            cloud_cover(iindex,jindex,  &
                ktop(iindex,jindex,n):kbtm(iindex,jindex,n)) =   &
                camt(iindex,jindex,n) * emcld(iindex,jindex,n)
         enddo
      enddo
      enddo
      
!     TK Mod to mimic supersource.  This is a temporary patch
!     which needs to be cleaned up.  A search is made for locations
!     where the LW cloud emissivity has been set to 0.6.
!     If that location is an isolated cloud with the layer above
!     and below being cloud free, the cloud_cover value is left
!     as 0.6, otherwise it is set to 1.0.   This is the
!     cirrus1L test.   7/09/01

      do jindex=1,size(nclds,2)
      do iindex=1,size(nclds,1)
!        For supersource runs, there should not be any cloud
!        at the top model level (k=1), nor should the high cloud
!        issue come into play the bottom level (k=kx).      

         do klev = 2, kx-1
            if (abs(cloud_cover(iindex,jindex,klev)-0.6) .lt. 0.001) then
 
!             The cloud at level klev has been previously reset to 0.6
 
               if ((cloud_cover(iindex,jindex,klev-1) .gt. 0.0) .or. &
                   (cloud_cover(iindex,jindex,klev+1) .gt. 0.0)) then
!             The cloud with 0.6 LW emiss has a neighboring cloud
!              so reset the LW emiss (cloud_cover) to 1.0 as follows:
                  cloud_cover(iindex,jindex,klev) = 1.0
               end if
 
            end if
         enddo
      enddo
      enddo

!     TK  Can Artificially set cloud_cover to 1 at the test point:
!         This method of doing this is now obsolete.  See radiation_driver.
!      print *, '***** TK Artificially set cloud_cover(49,27,14) = 1 ***'
!      cloud_cover(49,27,14) = 1.0

      DO j=1,jx

!        Specify ozone for the latitude (zonal mean used):
         do klev=1,kx
            o3_mixrat(klev) = QO3(1,j,klev)
         end do
 
!        TK Printout diagnostics:
        if ((j .eq. 27)) then
!        if ((j .eq. 27) .or. (j .eq. 69)) then
!        if ((j .eq. 10000) .or. (j .eq. 10000)) then

           ipr = 48

           print *
           print *, 'TK chkpt 1 in mcm_lw_mod.f.  Trap input to '
           print *, '  lwcool for i,j = ', ipr, j

           print *, 'k     cloud_cover     PRESSURE'
           do kpr = 1, 14
              write (6, 120) kpr, cloud_cover (ipr,j,kpr), &
                    PRES(ipr,j,kpr) 
120           format (i2, 2x, e10.3, 3x, e10.3)
           end do
           print *
           print *, 'Surface pressure = ', PRES(ipr,j,15)
        
           print *
           print *, 'k   TEMP        RH2O        QO3         CAMT', &
           '      RRVCO2'
           do kpr = 1, 14
              write (6, 125) kpr, TEMP(ipr,j,kpr),RH2O(ipr,j,kpr), &
                  QO3(ipr,j,kpr), CAMT(ipr,j,kpr), RRVCO2 
125           format (i2, 2x, f10.6, 2x, 2(e10.5,2x), 2(e8.3,2x))
           end do
           print *
           print *, 'Surface temperature = ', TEMP(ipr,j,15)
           print *
        end if
        
!       TK note that lw_down_sfc is just a throwaway variable in fms, 
!           so a dummy variable is used.  Some extra input
!           arguments are sent to lwcool.  Send i=1 values of 
!           ozone since zonally symmetric.
!           Use air temps in deg C as input parameter.
!           Convert surface pressures from Pa to dynes/cm2 for
!           input to the SS-derived routine...
!           The extra parameter j is for debugging only...

         call lwcool ( cloud_cover(:,j,:), temp(:,j,1:kx)-t_freeze, &
            RH2O(:,j,1:kx), phalf(:,j,kp)*10., temp(:,j,kp), &
            dummy, HEATRA(:,j,:), TOPFLX(:,j), GRNFLX(:,j), j, &
            sigma_level(:,j,:), sigma_half_level(:,j,:),  &
            sigma_thick(:,j,:), o3_mixrat)

!        TK Printout diagnostics:
         if ((j .eq. 27) ) then
!         if ((j .eq. 27) .or. (j .eq. 69)) then
!         if ((j .eq. 10000) .or. (j .eq. 10000)) then

           ipr = 48

           print *
           print *, 'TK chkpt 2b in mcm_lw_mod.f.  Trap output from '
           print *, '  lwcool for i,j =  ', ipr, j

!          Convert from ly/min to W/m2:
           print *, 'k    HEATRA(K/d)   GRNFLX      TOPFLX  (W/m2) '
           do kpr = 1, 14
              write (6, 130) kpr, HEATRA(ipr,j,kpr)*24*3600,  &
                 GRNFLX(ipr,j)*697., TOPFLX(ipr,j)*697. 
130           format (i2, 2x, f10.6, 3x, f10.6, 3x, f10.6)
           end do
           print *
         end if
      ENDDO

!      Convert top of atm and sfc fluxes into units FMS is
!      expecting as output: (milliWatts/m2, K/day).
!      TK mod 5/14/01: use identical conversion factor to SS: 697.496.
!      Note that 4186/6 = 697.66667.

       GRNFLX = GRNFLX * 1.0E3 * 697.496
       TOPFLX = TOPFLX * 1.0E3 * 697.496

!       GRNFLX = GRNFLX * 1.0E3 * 4186/6.
!       TOPFLX = TOPFLX * 1.0E3 * 4186/6.

       HEATRA = HEATRA * 24 * 3600


      END SUBROUTINE MCM_LW_RAD

!#######################################################################

      subroutine lwcool(cloud_cover, tj, rj, psj, t_sfc, &
       lw_down_sfc, lw_cooling_rate, lw_up_toa, lw_net_sfc, jrow_flag, &
       sigma_level, sigma_half_level, sigma_thick, o3_mixrat)
 

!  compute net longwave cooling rates

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  routine :   lwcool
!  called by : impt
!  calls :     lwtran
!  purpose :   compute net longwave cooling rates

!  inputs:
!     tj(i,k)             temperature (deg C) at each grid point
!     rj(i,k)             mixing ratio of h2o (g/g) at each grid point
!     psj(i)              surface pressure (dynes/cm**2) at each grid point
!     t_sfc(i)            surface temperature (deg K) at each grid point
!     cloud_cover(i,k)    fractional cloud cover

!  outputs:
!     lw_up_toa(i)           upward flux at top of atmosphere (ly/min)
!     lw_net_sfc(i)          net flux at ground (ly/min)
!     lw_down_sfc(i)         downward flux at ground (ly/min)
!     lw_cooling_rate(i,k)   cooling rate at each grid point (deg/s)

!  procedure:
!     the temperature at each grid point is converted from celsius to
!     absolute (kelvin); from the temperatures planck blackbody
!     radiances are computed for each grid point.  the longwave
!     transmission coefficients (recomputed every ldisk timesteps) are
!     then combined with the blackbody radiances to determine clear
!     fluxes (i.e., the flux which would occur in the absence of
!     clouds) between all pairs of levels.  dqx(i,k,l), for example, is
!     the clear-path downward flux from vertical level l to level k.
!     then, from the known distribution of clouds, the actual upward and
!     downward fluxes at each level are computed.  these are used to
!     compute cooling rates as well as some diagnostic quantities.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      real, parameter :: sigma = 5.673e-5, c1 = 3.740e-5, c2 = 1.4385 
      real, parameter :: cp1 = 0.24, cp2 = 4.1867e7 * cp1 

      real,    intent(in)  :: cloud_cover(ix,0:kp)
      real,    intent(in)  :: tj(ix,kx), rj(ix,kx), psj(ix)
      real,    intent(in)  :: t_sfc(ix)
      real,    intent(in)  :: sigma_level(ix,kx), sigma_half_level(ix,0:kx)
      real,    intent(in)  :: sigma_thick(ix,kx), o3_mixrat(kx)
      integer, intent(in)  :: jrow_flag

      real,    intent(out) :: lw_up_toa(ix)
      real,    intent(out) :: lw_net_sfc(ix)
      real,    intent(out) :: lw_cooling_rate(ix,kx)
      real,    intent(out) :: lw_down_sfc(ix)


!----------------LOCAL ARRAY STORAGE------------------------------------
      real :: dqx (ix,0:kp,0:kp), uqx (ix,0:kp,0:kp)
      real :: dfxc(ix,0:kp),      ufxc(ix,0:kp)
      real :: psigl(ix,2:kx)

      integer kloud (ix,0:kp)
      integer kldtop(ix,2:kx), kldmid(ix,2:kx), kldbot(ix,2:kx)

      integer :: numtop(2:kx), nummid(2:kx), numbot (2:kx)
      real :: ufxkb (ix),   uqxkt (ix),   ufxb(ix,2:kx)
      real :: dfxkt (ix),   dqxkb (ix),   dqxb(ix,2:kx)
      real :: ppfkb (ix),   ppfkt (ix),   ppfb(ix,2:kx)

      real :: cdxa  (kx),   cuxa(0:kx)
      real :: dfx(ix,0:kp), ufx(ix,0:kp)

      real temp_kelvin(ix,0:kp)
      real planck_func(ix,0:kp)
      real lw_trans_coeff(ix,0:kx,0:kx), lw_toa_correct(ix,0:kx)
      real lw_abs_quarter(ix,0:kx,2)
      real dpbb(ix,0:kx)
      integer :: k, l, kb, kt
      real :: cof

!    TK: I took out these equivalences but needed to modify the
!         code below so that dfx and ufx were initialized properly.
!         Printout answers for test points reproduced.  5/03/01 
!       equivalence (dfx, dqx(1,0,0))
!       equivalence (ufx, uqx(1,0,kp))


!      data cdxa / km * 0.5, 1.0 /, cuxa / -1.0, kx * -0.5 /

      do k=1,km
        cdxa(k) = 0.5
      end do
      cdxa(kx) = 1.0

      cuxa(0) = -1.0
      do k=1,kx
        cuxa(k) = -0.5
      end do

! ----------------------------------------------------------------------
!  compute clear (cloud-free) fluxes between all pairs of levels
! ----------------------------------------------------------------------

! ****** compute pbb (planck blackbody radiance), and d(pbb)/dk

      do k=1,kx
        do i=1,ix
          temp_kelvin(i,k) = tj(i,k) + t_freeze
        end do
      end do

      do i=1,ix
        temp_kelvin(i, 0) = temp_kelvin(i,1)
        temp_kelvin(i,kp) = t_sfc(i)
        end do
      end do

      do k=0,kp
        do i=1,ix
          planck_func(i,k) = sigma*temp_kelvin(i,k)**4
        end do
      end do

      do k=0,kx
        do i=1,ix
          dpbb(i,k) = planck_func(i,k+1) - planck_func(i,k)
        end do
      end do

! ----------------------------------------------------------------------
!  Obtain transmission functions
! ----------------------------------------------------------------------
      call lwtran (planck_func, temp_kelvin, rj, psj,  &
        lw_trans_coeff, lw_toa_correct, lw_abs_quarter, jrow_flag,  &
        sigma_level, sigma_half_level, sigma_thick, o3_mixrat)

! ****** compute upward and downward fluxes at "flux" levels

      do i=1,ix
        dqx(i,0,0) = 0.0
      end do

      do k=1,kx
        do i=1,ix
          dqx(i,k,k) = cdxa(k) * dpbb(i,k) * lw_abs_quarter(i,k,2) + planck_func(i,k)
        end do
      end do

      do l=kx-1,0,-1
        do k=l+1,kx
          do i=1,ix
            dqx(i,k,l) = dqx(i,k,l+1) - lw_trans_coeff(i,k,l)*dpbb(i,l)
          end do
        end do
      end do

      do k=1,kx
        do i=1,ix
          dqx(i,k,0) = dqx(i,k,0) - lw_toa_correct(i,k) * planck_func(i,0)
        end do
      end do

      do k=0,kx
        do i=1,ix
          uqx(i,k,k+1) = cuxa(k) * dpbb(i,k) * lw_abs_quarter(i,k,1) + planck_func(i,k+1)
        end do
      end do

      do l=2,kp
        do k=0,l-2
          do i=1,ix
            uqx(i,k,l) = uqx(i,k,l-1) + lw_trans_coeff(i,k,l-1) * dpbb(i,l-1)
          end do
        end do
      end do

!  TK Add code to initialize ufx to take place of equivalence:
      do k = 0, kp
         do i = 1, ix
            ufx(i,k) = uqx(i,k,kp)
        end do
      end do

!  TK Add code to initialize dfx to take place of equivalence:
      do k = 0, kp
         do i = 1, ix
            dfx(i,k) = dqx(i,k,0)
        end do
      end do

! ----------------------------------------------------------------------
!  compute actual fluxes from cloud distribution
! ----------------------------------------------------------------------

      do k=0,kp
       do i=1,ix
        if (cloud_cover(i,k) .ne. 0) then
          kloud(i,k) = 1
        else
          kloud(i,k) = 0
        endif
       end do
      end do
      do k=2,kx
       do i=1,ix
        kldtop(i,k) = kloud(i,k) * (1 - kloud(i,k-1))
        kldmid(i,k) = kloud(i,k) *      kloud(i,k+1)
        kldbot(i,k) = kloud(i,k) * (1 - kloud(i,k+1))
       end do
      end do

      do kb=2,kx
        do k=kb,kx
         do i=1,ix
          if (kldbot(i,kb) .eq. 1) &
            dfx(i,k) = dfx(i,k) * (1.0-cloud_cover(i,kb)) &
                     + dqx(i,k,kb) * cloud_cover(i,kb)
        end do
       end do
      end do

      do kt=kx,2,-1
        do k=0,kt-1
         do i=1,ix
          if (kldtop(i,kt) .eq. 1) &
            ufx(i,k) = ufx(i,k) * (1.0-cloud_cover(i,kt)) &
                     + uqx(i,k,kt) * cloud_cover(i,kt)
        end do
       end do
      end do

      do k=2,kx
        numtop(k) = 0
        numbot(k) = 0
        nummid(k) = 0
      end do
      do k=2,kx
       do i=1,ix
        numtop(k) = numtop(k) + kldtop(i,k)
        numbot(k) = numbot(k) + kldbot(i,k)
        nummid(k) = nummid(k) + kldmid(i,k)
       end do
      end do

      do k=kx,2,-1
       do i=1,ix
        if (kldbot(i,k) .eq. 1) then
          ufxkb(i) = ufx(i,k)
          dqxkb(i) = dqx(i,k,k)
          ppfkb(i) = sigma_half_level(i,k)
        endif
       end do
       if (nummid(k) .eq. 0) goto 40
       do i=1,ix
        if (kldmid(i,k) .eq. 1) then
          ufxb(i,k) = ufxkb(i)
          dqxb(i,k) = dqxkb(i)
          ppfb(i,k) = ppfkb(i)
        endif
       end do
  40  continue

      do k=2,kx
       if (nummid(k) .eq. 0) goto 50
         do i=1,ix
          if (kldtop(i,k) .eq. 1) then
            dfxkt(i) = dfx(i,k-1)
            uqxkt(i) = uqx(i,k-1,k)
            ppfkt(i) = sigma_half_level(i,k-1)
          endif
        end do
         do i=1,ix
          if (kldmid(i,k) .eq. 1) then
            psigl(i,k) = (sigma_half_level(i,k) - ppfkt(i)) / &
                            (ppfb(i,k) - ppfkt(i))
            dfxc (i,k) = dfxkt(i) + psigl(i,k)*(dqxb(i,k) - dfxkt(i))
            ufxc (i,k) = uqxkt(i) + psigl(i,k)*(ufxb(i,k) - uqxkt(i))
          endif
        end do
      end do

      do k=2,kx
       do i=1,ix
        if (kldmid(i,k) .eq. 1) then
          dfx(i,k) = dfx(i,k)*(1.0-cloud_cover(i,k)) + dfxc(i,k)*cloud_cover(i,k)
          ufx(i,k) = ufx(i,k)*(1.0-cloud_cover(i,k)) + ufxc(i,k)*cloud_cover(i,k)
        endif
      end do
    end do

! ----------------------------------------------------------------------
!  compute net cooling rates from fluxes
! ----------------------------------------------------------------------

! 1.43306e-6 converts ergs / (cm**2*sec) to ly/min

   80 continue

        do i=1,ix
         lw_up_toa(i)     = ufx(i,0) * 1.43306e-06
         lw_net_sfc(i)     = (ufx(i,kx) - dfx(i,kx)) * 1.43306e-06
         lw_down_sfc(i) = -dfx(i,kx) * 1.43306e-6
        end do

! 1.1574e-5 = 1.0 / 86400.0

        cof = grav_accel * 86400.0 / cp2
        cof = cof * 1.1574e-5

        do k=1,kx
         do i=1,ix
          lw_cooling_rate(i,k) = cof/(sigma_thick(i,k)*psj(i))
         end do
        end do

        do k=1,kx
         do i=1,ix
          lw_cooling_rate(i,k) = lw_cooling_rate(i,k) &
             * (ufx(i,k) - ufx(i,k-1) + dfx(i,k-1) - dfx(i,k))
         end do
        end do

      return
      end subroutine lwcool

function  expx(x)
real, intent(in) :: x
real :: expx

! ----------------------------------------------------------------------
!  statement function
! ----------------------------------------------------------------------

!  the approximation to exp(x) which follows is appropriate only given
!  the following assumptions:
!    1) since the correction coefficients (h2o_corr_a & h2o_corr_b) 
!       are known only to about four digits of accuracy, 6-7 digits 
!       of accuracy for the approximation is ample.
!    2) the input range for x is [-6.735,2.289].  this depends on:
!       a) the valid temperature range, now [-173.16,102] (celsius)
!       b) the values of the coefficients h2o_corr_a & h2o_corr_b
!

      real, parameter :: c0e = 9.999999810120017e-1 
      real, parameter :: c1e = 9.999999599140620e-1/16.**1 
      real, parameter :: c2e = 5.000056658124412e-1/16.**2 
      real, parameter :: c3e = 1.666837874572320e-1/16.**3 
      real, parameter :: c4e = 4.146341773620255e-2/16.**4 
      real, parameter :: c5e = 7.279860860195404e-3/16.**5 

      expx = ((( (((((c5e*x+c4e)*x+c3e)*x+c2e)*x+c1e)*x+c0e) &
                   **2)**2)**2)**2
end function expx


      subroutine lwtran (planck_func, temp_kelvin, rj, psj, &
        lw_trans_coeff, lw_toa_correct, lw_abs_quarter, jrow_flag, &
        sigma_level, sigma_half_level, sigma_thick, o3_mixrat)

!  compute longwave atmospheric transmission coefficients

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!  routine :   lwtran
!  called by : lwcool
!  calls :     lwco2
!  purpose :   compute longwave atmospheric transmission coefficients

!  inputs:
!     rj(i,k)             mixing ratio of h2o at each grid point (g/g)
!     psj(i,k)            surface presure at each ground point (dynes/cm**2)
!     temp_kelvin(i,k)    temperature (deg K)
!     planck_func(i,k)    blackbody function (total radiance; ergs/(cm**2)*s)

!  outputs:
!     lw_trans_coeff         longwave transmission coefficient (dimensionless)
!     lw_toa_correct         correction to above for top level (dimensionless)
!     lw_abs_quarter         quarter-level absorption (dimensionless)

!  procedure:
!     the main purpose of this routine is to compute the band-integrated
!     transmission coefficients transm(i,k,l) for radiation passing from
!     flux level l to flux level k.  in addition, it computes two
!     corrections. the first, transz(i,k), is for radiation from the top
!     level (split off from transm to improve code efficiency, i think).
!     the second, absrbx(i,k,kq), is a half-layer (flux level to data-
!     input level, or vice versa) absorption correction; kq=1 is the
!     correction for upward flux and kq=2 is for downward flux.

!     this routine splits the longwave radiation spectrum into nb bands,
!     ranging from about 20 to 2000 cm-1.  transmission coefficients are
!     computed for each absorber contributing to absorption in that
!     band.  the single-absorber transmission coefficients are
!     multiplied together to produce the total single-band transmission,
!     and then subtracted from 1 to produce an absorption.  the single-
!     band absorptions are weighted by the relative fraction of
!     blackbody radiation in that band (for a particular temperature),
!     and summed to produce the total band-integrated absorption
!     coefficient, which is converted back to a transmission.

!     transmissions for h2o are computed using a goody random model.
!     the function form and coefficients come from a paper by rogers and
!     walshaw (rogers, c.d. and c.d. walshaw, the computation of infra-
!     red cooling in planetary atmospheres, quart. j. royal meteorol.
!     soc., v.92, pp.67-92, 1966).

!        the water vapor continuum is now computed by using the
!     formulation given in roberts(first fit). reference:
!     roberts,r.e.,et al,applied optics,v. 15,pp 2085-2090,1976.

!     the co2 transmission table was generated by a program obtained
!     privately from s.r. drayson, university of michigan.  this program
!     assumed a homogeneous path between flux levels; a version which
!     did not assume homogeneous paths was later published by drayson
!     (drayson, s.r., atmospheric transmission in the co2 bands between
!     12 microns and 18 microns, appl. opt., v.5, pp.385-391, 1973).

!     the o3 transmission computation, changed from a table lookup in
!     the original longwave radiation code, also uses a goody random
!     model.  the functional form and coefficients were published by
!     rogers (rogers, c.d., some extensions and applications of the new
!     random model for molecular band transmission, q.j.r.m.s, v.94,
!     pp.99-102, 1968).

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, intent(in) :: planck_func(ix,0:kp)
      real, intent(in) :: temp_kelvin(ix,0:kp)
      real, intent(in) :: rj(ix,kx)
      real, intent(in) :: psj(ix)

      real, intent(in) :: sigma_level(ix,kx), sigma_half_level(ix,0:kx)
      real, intent(in) :: sigma_thick(ix,kx), o3_mixrat(kx)
      integer, intent(in)  :: jrow_flag

      real, intent(out) :: lw_trans_coeff(ix,0:kx,0:kx)
      real, intent(out) :: lw_toa_correct(ix,0:kx)
      real, intent(out) :: lw_abs_quarter(ix,0:kx,2)

      real, parameter :: sigma = 5.673e-5, c1 = 3.740e-5, c2 = 1.4385 
      real, parameter :: po = 1013.25e3, opo = 1.0 / po, o2po = 0.5 * opo 
!     parameter ( opo = 9.8692271e-7, o2po = 4.9346136e-7 )
      real, parameter :: ak1 = 223.0, ak2=10.0 
!     parameter ( alpha = 0.28 )
!     parameter ( pi = 3.1415926535898 )
!     parameter ( pialph = pi * alpha )
      real, parameter :: ak1x4 = 4.0 * ak1 
      real, parameter :: ak2x4 = 4.0 * ak2 
!     parameter ( pialph = 0.8796459, ak1x4 = 892.0, ak2x4 = 40.0 )
      real, parameter ::            azint = 1.66                 
!     parameter ( g = 980.6, azint = 1.66, g1 = azint / g )
!     parameter ( g1h = 0.5 * g1 )
!     parameter ( g1 = 1.6938774e-3, g1h = 0.5 * g1 )
      real, parameter :: c3 = 1.0 / 120.0, a1 = 75.48, a2 = 23.44 
!     parameter ( c3 = 8.333333e-3, a1 = 75.48, a2 = 23.44 )
      real, parameter :: c3a1 = c3 * a1, c3a2 = c3 * a2 
      real, parameter :: c0e = 9.999999810120017e-1 
      real, parameter :: c1e = 9.999999599140620e-1/16.**1 
      real, parameter :: c2e = 5.000056658124412e-1/16.**2 
      real, parameter :: c3e = 1.666837874572320e-1/16.**3 
      real, parameter :: c4e = 4.146341773620255e-2/16.**4 
      real, parameter :: c5e = 7.279860860195404e-3/16.**5 

      real :: bw(ix,0:kx,nb), bwz(ix,nb), expzs(ix,nb)
      real :: trans(ix,0:kx,0:kx), dbbdt(ix,0:kx)

      real :: vcube(nb), c(nb), v(nb), dv(nb)
      real :: pdi   (ix,0:kp), pfl   (ix,0:kx), dpdi (ix,  kx)
      real :: rh2o  (ix,0:kp), ro3   (ix,0:kp), qdi  (ix,  kx)
      real :: duh2o (ix,  kx), duco2 (ix,  kx), duo3 (ix,  kx)
      real :: duqh2o(ix,  kx), duqfac(ix,  kx), dufac(ix,  kx)
      real :: fbignl(ix,  kx), bignel(ix,0:kx), bign (ix,  kx)
      real :: tdpfac(ix,  kx), tfac  (ix,  kx), teff (ix,  kx)
      real :: sfac  (ix,0:kx), pfac  (ix,0:kx), peff (ix,0:kx)
      real :: tcof  (ix,0:kx), fac1  (ix,0:kx), fac2 (ix,0:kx)
      real :: t260  (ix,0:kx), tadj  (ix,0:kx,2)
      real :: ulog  (ix,0:kx), plog  (ix,0:kx)


!  quarter-layer variables

      real :: txdegk(ix,0:kx), qx    (ix,0:kx), dpfac(ix,0:kx)
      real :: duxh2o(ix,0:kx), duxco2(ix,0:kx), duxo3(ix,0:kx)
      real :: g1, g1h, bc,  bandc1, bandc2, cco2, asodsq, diag, &
              factr1, factr2, diago3
      integer :: lb, k, l, len, kq

!c  TK: Removed unnecessary equivalences:
!      equivalence (ulog,  fac1  ), (plog,  fac2 ), (teff, peff )
!      equivalence (txdegk,dpdi  ), (qx,    qdi  )
!      equivalence (duxh2o,duh2o ), (duxco2,duco2), (duxo3,duo3 )
!      equivalence (tadj,  duqfac), (tadj(1,0,2),dufac,tdpfac,dpfac)

!  constants for h2o bignell continuum absorption

!        these are the roberts continuum coefficients
      data  c / 4*0.0, 149.566, 38.103, 11.312, 7.9193, 6.0741, 4.9179, 9*0.0 /

!  longwave frequency bands and band widths (inverse cm)

      data v /   20.0,  100.0,  220.0,  340.0,  470.0, &
                670.0,  850.0,  930.0, 1020.0, 1140.0, &
               1275.0, 1400.0, 1500.0, 1600.0, 1700.0, &
               1800.0, 1900.0, 2000.0, 2125.0 /
      data dv /  40.0, 120.0, 120.0, 120.0, 140.0, &
                260.0, 100.0,  60.0, 120.0, 120.0, &
                150.0, 100.0, 100.0, 100.0, 100.0, &
                100.0, 100.0, 100.0, 150.0 /

! ----------------------------------------------------------------------
!  description of variables and constants used in lwtran
! ----------------------------------------------------------------------

!  variables centered on data input levels:
!     sigma_level             1-d "sigma" (pressure/po) 
!     sigma_thick             d(sigma)/dk
!                             (not to be confused with planck bb sigma)
!     pdi, dpdi, qdi          pressure, d(pdi)/dk, and pdi/po
!     rh2o, co2_mixrat, ro3         mixing ratios for the various absorbers
!     duh2o, duco2, duo3      absorber amounts from top of atmosphere
!        (gm/cm**2)           (also called unadjusted optical depth)
!     temp_kelvin                   temperature (kelvin)

!  variables centered on flux levels (see last paragraph, below):
!     sigma_half_level, pfl                1-d "sigma", 2-d pressure

!  assorted additional variables and constants:
!     po, opo, o2po           average pressure at sea level, 1/po, 2/po
!     planck_func             planck blackbody function
!     dbbdt                   d(planck_func)/dt
!     sigma, c1, c2           constants in planck blackbody equation
!     v, dv                   frequencies and band widths (1/cm)
!                             (wavelength range is 5-500 microns)
!     bw, bwz                 weights for freq. bands derived from planck_func
!     alpha, c3, a1, a2,      empirical constants used in computing
!       ak1, ak2                ozone absorption
!     h2o_line_width, 
!     h2o_line_strength       empirical band constants for water lines
!       h2o_line_width          coefficient for h2o line width
!       h2o_line_strength       coefficient for h2o line strength
!     sfac, pfac, tfac        intensity, pressure, temperature factors
!     peff, teff              effective pressure, temperature
!     ulog, plog              logs of co2 amount, effective pressure
!     tadj                    temperature adjustments for h2o lines
!     fbignl, bignel, bign    factors associated with bignell continuum
!     tcof                    transmission coefficient, single absorber
!     trans                   transmission matrix, single frequency band

!  output variables:
!     transm, transz          band-integrated transmission coefficients
!     absrbx                  quarter-layer correction to above (absrp.)

!  quarter layer variables (x usually indicates quarter-layer):
!     txdegk, qx              quarter-layer temperature, "sigma"
!     duxh2o, duxco2, duxo3   quarter-layer absorber amounts

!  miscellaneous constants:
!     azint                   azimuthal integral factor for diffuse rad.
!     g1, g1h                 convert dp*mix-ratio to absorber amount

!  miscellaneous temporaries:
!     expzs, vcube, bc, bandc1, bandc2
!     t260, fac1, fac2, dufac, dpfac, duqfac, duqh2o, cco2

!     vertical levels in the model are numbered from the top of the
!  atmosphere to the bottom; thus vertical level 1 is in the
!  stratosphere and level kx is just above the ground.  in the
!  radiation code, additional "pseudo-levels" are tacked onto the top
!  (level 0) and the bottom (level kp).  these levels make it easier
!  to handle incoming solar radiation and longwave radiation emanating
!  from the ground.

!     an additional complication to the level numbering system for the
!  radiation code is the use of "flux" levels and "data-input" levels.
!  the simplest numbering scheme for radiation calculations focuses
!  on even and odd levels: the grid variables (pressure, temperature,
!  etc.) reside at even levels, and flux is computed across them, i.e.,
!  from odd level to odd level.  this a pedagogically sound approach,
!  but it poses problems for computers, particularly vector types.
!  thus this code numbers what would be odd and even levels separately.
!  "data-input" levels are normal model vertical levels; "flux" levels
!  are halfway between them.  flux level k is one-half level closer to
!  the ground than data-input level k.  in the diagram below, data-input
!  levels are shown as "++++++" and flux levels are shown as "------".

!           description          level  index
!    incoming solar radiation      0      0   +++++++
!                                 1/2     0   -------
!             top model level      1      1   +++++++
!                                  .
!                                  .
!                                  .
!          bottom model level     kx     kx   +++++++
!                               kx+1/2   kx   -------
!                ground level     kp     kp   +++++++

! ----------------------------------------------------------------------
!  Compute various constants involving gravity.
! ----------------------------------------------------------------------
      g1  = azint / grav_accel
      g1h = 0.5 * g1

! ----------------------------------------------------------------------
!  compute band weights
! ----------------------------------------------------------------------

      do lb=1,nb
        vcube(lb) = v(lb) * v(lb) ** 2
      end do

      do lb=1,nb
        do i=1,ix
          expzs(i,lb) = exp(c2 / temp_kelvin(i,0) * v(lb)) - 1.0
        end do
      end do

      do lb=1,nb
        bc = c1 * vcube(lb) * dv(lb)
        do i=1,ix
          bwz(i,lb) = bc / (expzs(i,lb) * planck_func(i,0))
        end do
      end do

      do k=0,kx
        do i=1,ix
          dbbdt(i,k) = 4.0 * sigma * temp_kelvin(i,k)**3
        end do
      end do

      do lb=1,nb
        bandc1 = c1 * vcube(lb)
        bandc2 = c2 * v(lb)
        bc     = bandc1 * bandc2 * dv(lb)
        do k=0,kx
          do i=1,ix
            bw(i,k,lb) = exp(bandc2 / temp_kelvin(i,k))
            bw(i,k,lb) = bc * bw(i,k,lb) &
                         / ((bw(i,k,lb) - 1.0) * temp_kelvin(i,k))**2
            bw(i,k,lb) = bw(i,k,lb) / dbbdt(i,k)
          end do
        end do
      end do

! ----------------------------------------------------------------------
!  compute field arrays needed for trans calculation
!    (note that temp_kelvin has already been computed and stored in flux)
! ----------------------------------------------------------------------

      do i=1,ix
        pdi(i, 0) = 0.0
        pfl(i, 0) = 0.0
        pdi(i,kp) = psj(i)
      end do

      do k=1,kx
        do i=1,ix
          pdi (i,k) = psj(i) * sigma_level(i,k)
          pfl (i,k) = psj(i) * sigma_half_level(i,k)
          dpdi(i,k) = psj(i) * sigma_thick(i,k)
          ro3 (i,k) = o3_mixrat(k)
        end do
      end do

      do k=1,kx
        do i=1,ix
          rh2o  (i,k) = rj(i,k)
          qdi   (i,k) = (pfl(i,k) + pfl(i,k-1)) * o2po
          t260  (i,k) = temp_kelvin(i,k) - 260.0
          fbignl(i,k) = exp(1800.0/temp_kelvin(i,k)-6.081081081)
!        the kernel is      1800.0*(1.0/temp-1.0/296.0)
        end do
      end do

      do i=1,ix
        rh2o(i, 0) = 0.0
        ro3 (i, 0) = 0.0
        rh2o(i,kp) = rh2o(i,kx)
        ro3 (i,kp) = ro3 (i,kx)
      end do

      do k=1,kp
       do i=1,ix
         rh2o(i,k) = max(rh2o(i,k),3.0e-6)
       end do
      end do

      cco2 = co2_mixrat * g1
      do k=1,kx
        do i=1,ix
          duh2o (i,k) = dpdi  (i,k) * rh2o  (i,k) * g1
          duco2 (i,k) = dpdi  (i,k) * cco2
          duo3  (i,k) = dpdi  (i,k) * ro3   (i,k) * g1
          duqh2o(i,k) = duh2o (i,k) * qdi   (i,k)
          bignel(i,k) = (fbignl(i,k)*duqh2o(i,k)*rh2o(i,k)) / (0.622+rh2o(i,k))
        end do
      end do

! ----------------------------------------------------------------------
!  compute transm and transz

!  transm(i,k,l) is the transmission coefficient for longwave radiation
!    passing from vertical level l to level k
!  transz is a correction for radiation from the top of the atmosphere
! ----------------------------------------------------------------------

      do l=0,kx
       do k=0,kx
        do i=1,ix
         trans (i,k,l) = 0.0
         lw_trans_coeff(i,k,l) = 0.0
        end do
       end do
      end do
      do k=0,kx
       do i=1,ix
        lw_toa_correct(i,k) = 0.0
       end do
      end do

!  ***** begin loop on longwave frequency bands *****

      do lb= 1, nb

! ****** h2o contribution

! ******  compute temperature adjustment to line strength and line width

        do l=1,2
          do k=1,kx
            do i=1,ix
              tadj(i,k,l) = (h2o_corr_a(lb,l) + h2o_corr_b(lb,l) * t260(i,k)) * t260(i,k)
              tadj(i,k,l) = expx(tadj(i,k,l))
            end do
          end do
        end do

        do k=0,kx
         do i=1,ix
          sfac(i,k) = 0.0
          pfac(i,k) = 0.0
         end do
        end do
        do k=1,kx
         do i=1,ix
          bign(i,k) = 0.0
         end do
        end do

!     compute h2o contribution to lower triangular half of trans matrix
!     note that a non-physical factor of sod has been folded into duqfac
!     to cancel the factor that arises when sfac is squared in 1550 loop

        asodsq = h2o_line_strength(lb) * h2o_line_width(lb)
        if (h2o_line_strength(lb) .eq. 0.0) asodsq = h2o_line_width(lb)
        do k=1,kx
          do i=1,ix
            dufac (i,k) = h2o_line_strength(lb) * duh2o (i,k) * tadj(i,k,2)
            duqfac(i,k) = asodsq  * duqh2o(i,k) * tadj(i,k,1)
          end do
        end do

        do l=kx-1,0,-1
          do k=l+1,kx
            do i=1,ix
              sfac(i,k) = sfac(i,k) + dufac (i,l+1)
              pfac(i,k) = pfac(i,k) + duqfac(i,l+1)
              bign(i,k) = bign(i,k) + bignel(i,l+1)
            end do
          end do
        end do

          do k=l+1,kx
            do i=1,ix
              trans(i,k,l) = exp(-c(lb) * bign(i,k) - sfac(i,k) / sqrt(1.0 + sfac(i,k)**2 / pfac(i,k)))
            end do
          end do

! ****** compute h2o contribution to diagonal elements of trans matrix

        diag = exp(-c(lb) - h2o_line_strength(lb)/sqrt(1.0 + h2o_line_strength(lb)/h2o_line_width(lb)))

! ****** co2 contribution

        if (lb .ne. 6) goto 2500

        do k=0,kx
         do i=1,ix
          sfac(i,k) = 0.0
          pfac(i,k) = 0.0
         end do
        end do
        do k=1,kx
         do i=1,ix
          tfac(i,k) = 0.0
         end do
        end do

! **** compute co2 contribution to lower triangular half of trans matrix

        do k=1,kx
          do i=1,ix
            duqfac(i,k) = duco2(i,k) * qdi (i,k)
            tdpfac(i,k) = temp_kelvin(i,k) * dpdi(i,k)
          end do
        end do

        do l=kx-1,0,-1
          do k=l+1,kx
            do i=1,ix
              sfac(i,k) = sfac(i,k) + duco2 (i,l+1)
              pfac(i,k) = pfac(i,k) + duqfac(i,l+1)
              tfac(i,k) = tfac(i,k) + tdpfac(i,l+1)
              teff(i,k) = tfac(i,k) / (pfl(i,k) - pfl(i,l))
            end do
          end do
          do k=l+1,kx
            do i=1,ix
              ulog(i,k) = alog10(sfac(i,k))
              plog(i,k) = alog10(pfac(i,k)) - ulog(i,k)
              ulog(i,k) = ulog(i,k) + 2.70668
            end do
          end do

!  ---------BEGIN TK MODIFICATIONS

          if (l .eq. kx-1) then
            i = ix
            ulog(i,kx-1) = 2.70668
            plog(i,kx-1) = 0.0
            teff(i,kx-1) = 0.0
            len = ix + 1

! ******    call lwco2 to compute tcof values

            call lwco2(ulog(i,kx-1),  plog(i,kx-1), teff(i,kx-1), &
                       len, tcof(i,kx-1))

! ******    compute co2 contribution to diagonal elements of trans matrix

            diag = tcof(i,kx-1) * diag

           else
             i = 1
             len = ix * (kx - l)

! ******    call lwco2 to compute tcof values

            call lwco2(ulog(i,l+1),  plog(i,l+1), teff(i,l+1), &
                       len, tcof(i,l+1))
          endif

!  -----------END OF TK MODIFICATIONS -------------------

! ****** accumulate co2 tcof into trans (lower triangular half)

          do k=l+1,kx
            do i=1,ix
              trans(i,k,l) = tcof(i,k) * trans(i,k,l)
            end do
          end do
        end do

! ****** ozone contribution

 2500   if (lb .ne. 9) go to 2855

        do k=0,kx
         do i=1,ix
          sfac(i,k) = 0.0
          pfac(i,k) = 0.0
         end do
        end do

! **** compute o3 contribution to lower triangular half of trans matrix

        do k=1,kx
          do i=1,ix
            duqfac(i,k) = pi_alpha * duo3(i,k) * qdi(i,k)
          end do
        end do

        do l=kx-1,0,-1
          do k=l+1,kx
            do i=1,ix
              sfac(i,k) = sfac(i,k) + duo3  (i,l+1)
              pfac(i,k) = pfac(i,k) + duqfac(i,l+1)
            end do
          end do
          do k=l+1,kx
            do i=1,ix
              peff (i,k  ) = pfac(i,k) / sfac(i,k)
              fac1 (i,k  ) = sfac(i,k) / peff(i,k)
              fac2 (i,k  ) = 1.0 - exp(-5.0 * peff(i,k) * (sqrt(1.0 + ak2x4 * fac1(i,k)) - 1.0))
              fac1 (i,k  ) = 1.0 - exp(-5.0 * peff(i,k) * (sqrt(1.0 + ak1x4 * fac1(i,k)) - 1.0))
              tcof (i,k  ) = 1.0 - (c3a1 * fac1(i,k) + c3a2 * fac2(i,k))
              trans(i,k,l) = tcof(i,k) * trans(i,k,l)
            end do
          end do
        end do

! ****** compute o3 contribution to diagonal elements of trans matrix

        factr1 = 1.0 - exp(-5.0 * pi_alpha * (sqrt(1.0 + ak1x4 / pi_alpha) - 1.0))
        factr2 = 1.0 - exp(-5.0 * pi_alpha * (sqrt(1.0 + ak2x4 / pi_alpha) - 1.0))
        diago3 = 1.0 - (c3a1 * factr1 + c3a2 * factr2)
        diag   = diago3 * diag

! ****** apply band weights to trans; accumulate
! ****** as absorp. in transm & transz

! ****** first convert trans and diag to absorption

 2855   continue
        do l=0,kx
         do k=0,kx
          do i=1,ix
           trans(i,k,l) = 1.0 - trans(i,k,l)
          end do
         end do
        end do

        diag  = 1.0 - diag

        do l=0,kx-1
          do k=l+1,kx
            do i=1,ix
              lw_trans_coeff(i,k,l) = lw_trans_coeff(i,k,l) + trans(i,k,l) * bw(i,l,lb)
            end do
          end do
        end do

        do l=1,kx
          do k=0,l-1
            do i=1,ix
              lw_trans_coeff(i,k,l) = lw_trans_coeff(i,k,l) + trans(i,l,k) * bw(i,l,lb)
            end do
          end do
        end do

        do l=0,kx
          do i=1,ix
            lw_trans_coeff(i,l,l) = lw_trans_coeff(i,l,l) + diag * bw(i,l,lb)
          end do
        end do

        do l=1,kx
          do i=1,ix
            lw_toa_correct(i,l) = lw_toa_correct(i,l) + trans(i,l,0) * bwz(i,lb)
          end do
        end do

        do i=1,ix
          lw_toa_correct(i,0) = lw_toa_correct(i,0) + diag * bwz(i,lb)
        end do

      end do

!  ***** end loop on longwave frequency bands *****

!  convert transm & transz from absorption to transmission coefficients

      do l=0,kx
       do k=0,kx
        do i=1,ix
         lw_trans_coeff(i,k,l) = 1.0 - lw_trans_coeff(i,k,l)
        end do
       end do
      end do
      do k=0,kx
       do i=1,ix
        lw_toa_correct(i,k) = 1.0 - lw_toa_correct(i,k)
       end do
      end do

! ----------------------------------------------------------------------
!  compute absrbx

!  absrbx is a correction for the finite thickness of the vertical
!    levels, and represents within-layer absorption of radiation
! ----------------------------------------------------------------------

!     absrbx = 0.0
      do l=1,2
       do k=0,kx
        do i=1,ix
         lw_abs_quarter(i,k,l) = 0.0
        end do
       end do
      end do

!  kq=1 is the upward-flux correction; kq=2 is for downward flux

      do kq=1,2

!  compute coefficients needed for quarter-layer correction

        if (kq.eq.2) goto 3050

        do k=1,kx-1
          do i=1,ix
            dpfac(i,k) = (pdi(i,k+1) - pfl(i,k)) * g1h
          end do
        end do

        do k=1,kx-1
          do i=1,ix
            txdegk(i,k) = 0.375 * temp_kelvin(i,k) + 0.625 * temp_kelvin(i,k+1)
            qx    (i,k) = (0.75 * pfl(i,k) + 0.25 * pdi(i,k+1)) * opo
            duxh2o(i,k) = (0.375 * rh2o(i,k) + 0.625 * rh2o(i,k+1)) * dpfac(i,k)
            duxco2(i,k) = co2_mixrat * dpfac(i,k)
            duxo3 (i,k) = (0.375 * ro3(i,k) + 0.625 * ro3(i,k+1)) * dpfac(i,k)
            bignel(i,k) = (0.375 * rh2o(i,k) + 0.625 * rh2o(i,k+1))**2 &
                          * dpfac(i,k) * fbignl(i,k+1) / &
                          (0.622 + 0.375*rh2o(i,k) + 0.625*rh2o(i,k+1))
          end do
        end do

        do i=1,ix
          txdegk(i, 0) = 0.75 * temp_kelvin(i,0) + 0.25 * temp_kelvin(i,1)
          txdegk(i,kx) = temp_kelvin(i,kp)
          qx    (i, 0) = (0.75 * pdi(i,0) + 0.25 * pdi(i,1)) * opo
          qx    (i,kx) = psj(i) * opo
          duxh2o(i, 0) = (0.25 * rh2o(i,0) + 0.75 * rh2o(i,1)) * (pdi(i,1) - pdi(i,0)) * g1h
          duxh2o(i,kx) = 0.0
          duxco2(i, 0) = co2_mixrat * (pdi(i,1) - pdi(i,0)) * g1h
          duxco2(i,kx) = 0.0
          duxo3 (i, 0) = (0.25 * ro3(i,0) + 0.75 * ro3(i,1)) * (pdi(i,1) - pdi(i,0)) * g1h
          duxo3 (i,kx) = 0.0
          bignel(i, 0) = (0.75 * rh2o(i,0) + 0.25 * rh2o(i,1))**2 &
                       * (pdi(i,1)-pdi(i,0)) * g1h * fbignl(i,1) / &
                         (0.622 + 0.75*rh2o(i,0) + 0.25*rh2o(i,1))
          bignel(i,kx) = 0.0
        end do

        do k=0,kx
          do i=1,ix
            t260(i,k) = txdegk(i,k) - 260.0
          end do
        end do

        go to 3100

 3050   do k=1,kx-1
          do i=1,ix
            dpfac(i,k) = (pfl(i,k) - pdi(i,k)) * g1h
          end do
        end do

        do k=1,kx-1
          do i=1,ix
            txdegk(i,k) = 0.625 * temp_kelvin(i,k) + 0.375 * temp_kelvin(i,k+1)
            qx    (i,k) = (0.75 * pfl(i,k) + 0.25 * pdi(i,k)) * opo
            duxh2o(i,k) = (0.625 * rh2o(i,k) + 0.375 * rh2o(i,k+1)) * dpfac(i,k)
            duxco2(i,k) = co2_mixrat * dpfac(i,k)
            duxo3 (i,k) = (0.625 * ro3(i,k) + 0.375 * ro3(i,k+1)) * dpfac(i,k)
            bignel(i,k) = (0.625*rh2o(i,k) + 0.375*rh2o(i,k+1))**2 &
                          * dpfac(i,k) * fbignl(i,k) / &
                          (0.622 + 0.625*rh2o(i,k) + 0.375*rh2o(i,k+1))
          end do
        end do

        do i=1,ix
          txdegk(i, 0) = temp_kelvin(i,0)
          txdegk(i,kx) = 0.75 * temp_kelvin(i,kp) + 0.25 * temp_kelvin(i,kx)
          qx    (i, 0) = 1.0
          qx    (i,kx) = (0.75 * pdi(i,kp) + 0.25 * pdi(i,kx)) * opo
          duxh2o(i, 0) = 0.0
          duxh2o(i,kx) = (0.75 * rh2o(i,kp) + 0.25 * rh2o(i,kx)) * (pdi(i,kp) - pdi(i,kx)) * g1h
          duxco2(i, 0) = 0.0
          duxco2(i,kx) = co2_mixrat * (pdi(i,kp) - pdi(i,kx)) * g1h
          duxo3 (i, 0) = 0.0
          duxo3 (i,kx) = (0.75 * ro3(i,kp) + 0.25 * ro3(i,kx)) * (pdi(i,kp) - pdi(i,kx)) * g1h
          bignel(i, 0) = 0.
          bignel(i,kx) = (0.75 * rh2o(i,kp) + 0.25 * rh2o(i,kx))**2 &
                       * (pdi(i,kp)-pdi(i,kx)) * g1h * fbignl(i,kx) / &
                         (0.622 + 0.75*rh2o(i,kp) + 0.25*rh2o(i,kx))
        end do

        do k=0,kx
          do i=1,ix
            t260(i,k) = txdegk(i,k) - 260.0
          end do
        end do

!  ***** begin loop on longwave frequency bands *****

 3100   do lb=1,nb

! ****** quarter-layer h2o absorption

! ****** compute temperature adjustment to line strength and line width

          do k=0,kx
            do i=1,ix
              tadj(i,k,1) = (h2o_corr_a(lb,1) + h2o_corr_b(lb,1) * t260(i,k)) * t260(i,k)
              tadj(i,k,2) = (h2o_corr_a(lb,2) + h2o_corr_b(lb,2) * t260(i,k)) * t260(i,k)
              tadj(i,k,1) = tadj(i,k,1) - tadj(i,k,2)
            end do
          end do

          do l=1,2
            do k=0,kx
              do i=1,ix
                tadj(i,k,l) = expx(tadj(i,k,l))
              end do
            end do
          end do

! ****** compute h2o contribution to quarter-layer absorption matrix

          do k=0,kx
            do i=1,ix
              sfac (i,k   ) = h2o_line_strength(lb) * duxh2o(i,k) * tadj(i,k,2)
              pfac (i,k   ) = h2o_line_width(lb) * qx(i,k) * tadj(i,k,1)
              trans(i,k,kq) = exp(-c(lb) * bignel(i,k) -sfac(i,k)/sqrt(1.0+sfac(i,k)/pfac(i,k)))
            end do
          end do

! ****** quarter-layer co2 absorption

          if (lb .ne. 6) go to 3700

! ****** compute co2 contribution to quarter-layer absorption matrix

          len = ix * kp
          do k=0,kx
           do i=1,ix
            sfac(i,k) = duxco2(i,k)
           end do
          end do
          do k=0,kx
           do i=1,ix
            duxco2(i,k) = max(duxco2(i,k),1.0e-6)
           end do
          end do
          do k=0,kx
            do i=1,ix
              ulog(i,k) = alog10(duxco2(i,k)) + 2.70668
              plog(i,k) = alog10(qx(i,k))
            end do
          end do

! ****** call lwco2 to compute tcof values

          call lwco2(ulog(1,0), plog(1,0), txdegk(1,0), len, tcof(1,0))

! ****** accumulate co2 ccp into trans

          do k=0,kx
            do i=1,ix
              tcof(i,k) = tcof(i,k) * trans(i,k,kq)
            end do
          end do
          do k=0,kx
           do i=1,ix
            if(sfac(i,k) .gt. 0.0) then
              trans(i,k,kq) = tcof(i,k)
            endif
           end do
          end do

! ****** quarter-layer o3 absorption

 3700     if (lb .ne. 9) go to 3900

          do k=0,kx
            do i=1,ix
              peff (i,k   ) = pi_alpha * qx(i,k)
              fac2 (i,k   ) = duxo3(i,k) / peff(i,k)
              fac1 (i,k   ) = 1.0 - exp(-5.0 * peff(i,k) * (sqrt(1.0 + ak1x4 * fac2(i,k)) - 1.0))
              fac2 (i,k   ) = 1.0 - exp(-5.0 * peff(i,k) * (sqrt(1.0 + ak2x4 * fac2(i,k)) - 1.0))
              tcof (i,k   ) = 1.0 - (c3a1 * fac1(i,k) + c3a2 *fac2(i,k))
              trans(i,k,kq) = tcof(i,k) * trans(i,k,kq)
            end do
          end do
 
! ****** total quarter-layer absorption
 
 3900     do k=0,kx
            do i=1,ix
              lw_abs_quarter(i,k,kq) = lw_abs_quarter(i,k,kq) + (1.0 - trans(i,k,kq)) * bw(i,k,lb)
            end do
          end do
      end do
 
!  ***** end of principal quarter-layer loop *****
 
      return
      end subroutine lwtran
 
      subroutine lwco2 (ulog,plog,teff,ln,tco2)

!  table lookup routine for co2 transmission
 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
!  routine:   lwco2
!  called by: lwtran
!  calls:     none
!  purpose:   compute longwave co2 transmission
 
!  inputs:
!     ulog           log of co2 amount
!     plog           log of effective pressure
!     teff           effective temperature
!     ln             length of call-argument arrays
 
!  output:
!     tco2           co2 transmission coefficient
 
!  side effects:
!     none
!     (i.e., no input variables or common blocks modified, etc.)
 
!  procedure:
!     a 2d table lookup of co2 absorption is performed.  the axes are
!     log of the co2 amount and log of the effective pressure.  an
!     effective temperature of 300k is assumed for this first step.
!     a correction for temperatures in the range 200k-300k is made by
!     multiplying a temperature correction factor by the difference
!     between the effective temperature and 300k, and adding the
!     result to the previously calculated absorption.  the temperature
!     correction factor is found by 2d table lookup, using the same
!     axes as were used to obtain the 300k absorption.  no temperature
!     correction is applied if the temperature exceeds 300k; if the
!     temperature falls below 200k, the 200k correction is applied.
!     if the co2 concentration is too low, no temperature correction
!     is made.
 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
!  TK Mod:  set up in parent routine:
!      parameter ( lmax = ix * (kx + 1) )
 
      integer, intent(in) :: ln
      real, intent(in) :: ulog(ln) 
      real, intent(in) :: plog(ln) 
      real, intent(in) :: teff(ln)
 
      real, intent(out) :: tco2(ln) 
 
      real :: ahi   (lmax), alo   (lmax), cct   (lmax), cdt   (lmax), fulog (lmax)
      integer :: nplog (lmax), iulog (lmax), index (lmax)
      real :: plnum (lmax), plden (lmax), plfac (lmax), tcor  (lmax)
      real :: alogoo(lmax), alogpo(lmax), alogop(lmax), alogpp(lmax)
      real :: pltab (  12)
      integer :: iu, l, n
 
!c TK: Removed unnecessary equivalences:
!c TK:      equivalence (alogoo, plnum), (alogpo, plden, alo, tcor)
!c TK:      equivalence (alogop, iulog), (alogpp, index, ahi, cdt )
 
      data pltab / &
           -4.50000, -4.00000, -3.40000, -2.88081, -2.27875, -1.88081, &
           -1.57978, -1.27875, -0.88081, -0.48287,  0.     ,  0.29528 /

! ****** check for out-of-range ulog input values

      iu = 0
      do l=1,ln
        if ( ulog(l).lt.-4.25 ) iu=iu+1
      end do
      if ( iu.ne.0 ) go to 1200

! ****** compute pressure range for tables,
! ****** and pressure interpolation factor

      do l=1,ln
        nplog(l) = 1
      end do

      do n=2,11
        do l=1,ln
          if ( plog(l).gt.pltab(n) ) then
            nplog(l) = n
            plnum(l) = plog(l) - pltab(n)
            plden(l) = pltab(n+1) - pltab(n)
          endif
        end do
      end do

      do l=1,ln
        plfac(l) = plnum(l) / plden(l)
      end do

! ****** compute co2 absorption for an effective temperature of 300k

      do l=1,ln
        fulog(l) = 4 * ulog(l) + 18
        iulog(l) = int(fulog(l))

!        TK diagnostic printout:
!        print *, 'TK chkpt 1 in lwco2:  l, ulog(l), fulog(l), ', &
!                   'int(fulog(l)), iulog(l) = ', &
!                  l, ulog(l), fulog(l), int(fulog(l)), iulog(l)

        if ( iulog(l).gt.29 ) iulog(l) = 29
        fulog(l) = fulog(l) - real(iulog(l))
      end do

      do l = 1, ln
        index(l) = iulog(l) + 30 * (nplog(l) - 1)
      end do

      do l=1,ln
        alogoo(l) = absorp_table(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + 1
        alogpo(l) = absorp_table(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + (30 - 1)
        alogop(l) = absorp_table(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + 1
        alogpp(l) = absorp_table(index(l),1)
      end do

      do l = 1, ln
        ahi(l) = alogop(l) + fulog(l) * (alogpp(l) - alogop(l))
        alo(l) = alogoo(l) + fulog(l) * (alogpo(l) - alogoo(l))
        cct(l) = alo   (l) + plfac(l) * (ahi   (l) - alo   (l))
      end do

! ****** compute temperature correction coefficient for co2 absorption

      do l=1,ln
        fulog(l) = ulog(l) + 4
        iulog(l) = int(fulog(l))
        if ( iulog(l).gt.6 ) iulog(l) = 6
        fulog(l) = fulog(l) - real(iulog(l))
      end do

      do l = 1, ln
        index(l) = iulog(l) + 7 * (nplog(l) - 1)
      end do

      do l=1,ln
        alogoo(l) = d_absorp_dt(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + 1
        alogpo(l) = d_absorp_dt(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + (7 - 1)
        alogop(l) = d_absorp_dt(index(l),1)
      end do

      do l=1,ln
        index (l) = index(l) + 1
        alogpp(l) = d_absorp_dt(index(l),1)
      end do

      do i = 1, ln
        ahi(i) = alogop(i) + fulog(i) * (alogpp(i) - alogop(i))
        alo(i) = alogoo(i) + fulog(i) * (alogpo(i) - alogoo(i))
        cdt(i) = alo   (i) + plfac(i) * (ahi   (i) - alo   (i))
      end do

! ****** restrict effective temperature input values to proper range

      do l=1,ln
        tcor(l) = teff(l)
        if ( tcor(l).gt.320.0 ) tcor(l) = 320
        if ( tcor(l).lt.180.0 ) tcor(l) = 180
      end do

! ****** apply temperature correction to co2
! ****** absorption if ulog is in range

      do l=1,ln
        cdt (l) = cct(l) + 0.01 * (cdt(l) * (tcor(l) - 300))
        tco2(l) = 1 - 1.3077 * cdt(l)
      end do

      do l=1,ln
        if ( ulog(l).lt.-3 ) tco2(l) = 1 - 1.3077 * cct(l)
      end do

      return

! ****** fatal error: ulog out-of-range

 1200 write(6,1201) ulog
      call Error_Mesg ( 'lwco2 in MCM_LW_MOD','ulog out of range.', FATAL)

 1201 format('1ulog > 3. or < -4.25; execution terminated' / ' the ulog vector is:' / (8e16.8))

        end subroutine lwco2

      end module mcm_lw_mod

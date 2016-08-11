
module fs_profile_mod

use fms_mod, only : mpp_pe, mpp_root_pe, write_version_number, &
                    error_mesg, FATAL

implicit none
private

!-----------------------------------------------------------------------
! **      THIS PROGRAM CALCULATES TEMPERATURES ,H2O MIXING RATIOS     **
! **      AND O3 MIXING RATIOS BY USING AN ANALYTICAL                 **
! **      FUNCTION WHICH APPROXIMATES                                 **
! **      THE US STANDARD (1976).  THIS IS                            **
! **      CALCULATED IN FUNCTION 'ANTEMP', WHICH IS CALLED BY THE     **
! **      MAIN PROGRAM.  THE FORM OF THE ANALYTICAL FUNCTION WAS      **
! **      SUGGESTED TO ME IN 1971 BY RICHARD S. LINDZEN.              **
!
!*****THIS VERSION IS ONLY USABLE FOR 1976 US STD ATM AND OBTAINS
!     QUANTITIES FOR CO2 INTERPOLATION AND INSERTION INTO OPERA-
!     TIONAL RADIATION CODES
!
!    definitions:       
!    -----------
!      pd,pd8: pressures (mb) for data levels. pd is for the case where
!              p(sfc)=1013.25 mb; pd8 applies when p(sfc)=810.6 mb.
!              in either case, index (nlev+1) is at the sfc.
!      press:  same as pd, but with indices reversed,index 1 at the
!              surface, and index (nlev+1) at the top (nonzero) data
!              level.
!-----------------------------------------------------------------------

 Public   fs_profile, fs_profile_init, fs_profile_end

 Real, Public, Allocatable, Dimension(:) :: pd1013,plm1013,pd810,plm810

!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: fs_profile.F90,v 10.0 2003/10/24 22:00:30 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
 logical            :: module_is_initialized = .false.

CONTAINS

!#######################################################################

  subroutine fs_profile (Pref,stemp,gtemp)

!-----------------------------------------------------------------------
  Real, Intent(IN) , Dimension(:,:) :: Pref
  Real, Intent(OUT), Dimension(:)   :: stemp,gtemp
!-----------------------------------------------------------------------

  Real,Dimension(Size(Pref,1)) :: press

  Real :: PSMAX = 1013.250

  Real     DELZAP,R,G0,ZMASS,HT,DZ,ZNINT,HTA,RK1,RK2,RK3,RK4,  &
           DLOGP,PSTAR
  Integer  n,m,nlev,NINT
!-----------------------------------------------------------------------

      nlev=Size(Pref,1)-1

      DELZAP=0.5
      R=8.31432
      G0=9.80665
      ZMASS=28.9644

      stemp(nlev+1)=ANTEMP(6,0.0)

!*******DETERMINE THE PRESSURES (press)
      PSTAR=PSMAX

      CALL SIGP (Pref,gtemp)    

!----- convert to mb -----
      Do n=1,nlev+1
         press(n)=Pref(nlev+2-n,1)*0.01
      EndDo
         press(1)=PSTAR

!    *** CALCULATE TEMPS ***

         HTA=0.0

         Do n=1,nlev

! **      ESTABLISH COMPUTATATIONAL LEVELS BETWEEN USER LEVELS AT     **
! **      INTERVALS OF APPROXIMATELY 'DELZAP' KM.                     **

            DLOGP=7.0*LOG(press(n)/press(n+1))
            NINT=DLOGP/DELZAP
            NINT=NINT+1
            ZNINT=NINT
            DZ=R*DLOGP/(7.0*ZMASS*G0*ZNINT)
            HT=HTA

! **      CALCULATE HEIGHT AT NEXT USER LEVEL BY MEANS OF             **
! **                RUNGE-KUTTA INTEGRATION.                          **

            Do m=1,NINT
               RK1=ANTEMP(6,HT)*DZ
               RK2=ANTEMP(6,HT+0.5*RK1)*DZ
               RK3=ANTEMP(6,HT+0.5*RK2)*DZ
               RK4=ANTEMP(6,HT+RK3)*DZ
               HT=HT+0.16666667*(RK1+RK2+RK2+RK3+RK3+RK4)
            EndDo
            HTA=HT
            stemp(nlev+1-n)=ANTEMP(6,HT)

         EndDo

!-----------------------------------------------------------------------

  end subroutine fs_profile

!#######################################################################

  Subroutine SigP (Pref,gtemp)

!-----------------------------------------------------------------------
  Real, Intent(IN) , Dimension(:,:) :: Pref
  Real, Intent(OUT), Dimension(:)   :: gtemp
!-----------------------------------------------------------------------

  Real, Dimension(Size(Pref,1)) :: pd, pd8, plm, plm8
  REAL     PSS
  Integer  k,nlev

!     PSS = surface pressure (specified)
!-----------------------------------------------------------------------

      nlev = Size(Pref,1)-1

      If (Allocated(pd1013)) DeAllocate (pd1013, plm1013, pd810, plm810)

      Allocate (pd1013(nlev+1), plm1013(nlev+1),  &
                pd810 (nlev+1), plm810 (nlev+1))

!-----------------------------------------------------------------------
!-------- first pass: PSS=1013.25 MB --------

      PSS = 1013250.

      pd(:)=Pref(:,1)*10.
      pd(nlev+1)=PSS

      plm(1)=0.
      Do k=1,nlev-1
         plm(k+1)=0.5*(pd(k)+pd(k+1))
      EndDo
      plm(nlev+1)=PSS

      Do k=1,nlev
         gtemp(k)=pd(k)**0.2*(1.+pd(k)/30000.)**0.8/1013250.
      EndDo
      gtemp(nlev+1)=0.

!--- pd1013, plm1013 are used by the co2 interpolation prgm (PS=1013mb)
!    THE FOLLOWING PUTS P-DATA INTO MB

       pd1013(:)= pd(:)*1.e-3
      plm1013(:)=plm(:)*1.e-3

!-----------------------------------------------------------------------
!-------- second pass: PSS=810MB, gtemp NOT COMPUTED --------

      PSS=0.8*1013250.

      pd8(:)=Pref(:,2)*10.
      pd8(nlev+1)=PSS

      plm8(1)=0.
      Do k=1,nlev-1
         plm8(k+1)=0.5*(pd8(k)+pd8(k+1))
      EndDo
      plm8(nlev+1)=PSS

!--- pd810, plm810 are used by the co2 interpolation prgm (PS=810mb)
!    THE FOLLOWING PUTS P-DATA INTO MB

       pd810(:)= pd8(:)*1.e-3
      plm810(:)=plm8(:)*1.e-3

!-----------------------------------------------------------------------

  End Subroutine SigP

!#######################################################################

  FUNCTION ANTEMP (L,Z)

!-----------------------------------------------------------------------
  Integer, Intent(IN)  :: L
  Real,    Intent(IN)  :: Z
  Real                 :: ANTEMP
!-----------------------------------------------------------------------
  Real    ZB(10,7),C(11,7),DELTA(10,7),TSTAR(7)
  Real    temp,expo,x,y,zlog,expp,faclog
  Integer nlast,n

!--------------- TROPICAL SOUNDING -------------------------------------
      Data (ZB(n,1),n=1,10)/   2.0,   3.0,   16.5,  21.5,  45., &
                              51.0,  70.0,  100.,  200.,  300.  /
      Data (C(n,1),n=1,11)/   -6.0,  -4.0,  -6.7,   4.0,   2.2,  &
                               1.0,  -2.8,  -.27,   0.0,   0.0,  0.0 /
      Data (DELTA(n,1),n=1,10)/ .5,    .5,    .3,    .5,   1.0,  &
                               1.0,   1.0,   1.0,   1.0,   1.0   /
!--------------- SUB-TROPICAL SUMMER -----------------------------------
      Data (ZB(n,2),n=1,10)/  1.5,   6.5,  13.0,  18.0,  26.0,  &
                             36.0,  48.0,  50.0, 70.0,  100./
      Data (C(n,2),n=1,11)/  -4.0,  -6.0,  -6.5,   0.0,   1.2,  &
                              2.2,   2.5,   0.0,  -3.0,  -0.25,  0.0/
      Data (DELTA(n,2),n=1,10)/ .5,  1.0,    .5,    .5,   1.0,  &
                               1.0,  2.5,    .5,   1.0,   1.0/
!--------------- SUB-TROPICAL WINTER -----------------------------------
      Data (ZB(n,3),n=1,10)/ 3.0,  10.0,  19.0,  25.0,  32.0,  &
                              44.5, 50.0,  71.0,  98.0,  200.0/
      Data (C(n,3),n=1,11)/  -3.5,  -6.0,  -0.5,  0.0,   0.4,  &
                              3.2,   1.6,  -1.8, -0.7,   0.0,   0.0/
      Data (DELTA(n,3),n=1,10)/ .5,   .5,  1.0,   1.0,   1.0,  &
                               1.0,  1.0,  1.0,   1.0,   1.0/
!--------------- SUB-ARCTIC SUMMER -------------------------------------
      Data (ZB(n,4),n=1,10)/ 4.7, 10.0,  23.0,  31.8,  44.0,  &
                              50.2, 69.2, 100.0, 102.0, 103.0/
      Data (C(n,4),n=1,11)/  -5.3, -7.0,   0.0,  1.4,   3.0,  &
                              0.7, -3.3,  -0.2,  0.0,   0.0,  0.0/
      Data (DELTA(n,4),n=1,10)/ .5,   .3,  1.0,   1.0,   2.0,  &
                               1.0,  1.5,  1.0,   1.0,   1.0/
!------------- SUB-ARCTIC WINTER ---------------------------------------
      Data (ZB(n,5),n=1,10)/ 1.0,   3.2,   8.5,   15.5,   25.0,  &
                              30.0,  35.0,  50.0,  70.0,  100.0/
      Data (C(n,5),n=1,11)/  3.0,  -3.2,  -6.8,  0.0,  -0.6,  &
                             1.0,   1.2,   2.5, -0.7,  -1.2,  0.0/
      Data (DELTA(n,5),n=1,10)/ .4,   1.5,    .3 ,   .5,   1.0,  &
                               1.0,   1.0,   1.0,   1.0,   1.0/
!------------- US STANDARD 1976 ----------------------------------------
      Data (ZB(n,6),n=1,10)/ 11.0,  20.0,  32.0,  47.0,  51.0,  &
                             71.0,  84.8520,  90.0,  91.0,  92.0/
      Data (C(n,6),n=1,11)/ -6.5,   0.0,   1.0,   2.80,  0.0,  &
                            -2.80, -2.00,  0.0,   0.0,   0.0,  0.0/
      Data (DELTA(n,6),n=1,10)/ 0.3,   1.0,   1.0,   1.0,   1.0,  &
                                1.0,   1.0,   1.0,   1.0,   1.0/

!------------- ENLARGED US STANDARD 1976 -------------------------------
      Data (ZB(n,7),n=1,10)/ 11.0,  20.0,  32.0,  47.0,  51.0,  &
                             71.0,  84.8520,  90.0,  91.0,  92.0/
      Data (C(n,7),n=1,11)/ -6.5,   0.0,   1.0,   2.80,  0.0,  &
                            -2.80, -2.00,  0.0,   0.0,   0.0,  0.0/
      Data (DELTA(n,7),n=1,10)/ 0.3,   1.0,   1.0,   1.0,   1.0,  &
                                1.0,   1.0,   1.0,   1.0,   1.0/


      Data TSTAR / 300.0,  294.0,  272.2,  287.0,  257.1, 2*288.15/

!-----------------------------------------------------------------------

      nlast=10
      temp=TSTAR(L)+C(1,L)*Z

      Do n=1,nlast
         expo=(Z-ZB(n,L))/DELTA(n,L)
         If (abs(expo) <= 60.) Then
            x=exp(expo)
            y=x+1.0/x
            zlog=log(y)
         Else
            zlog=abs(expo)
         EndIf
         expp=ZB(n,L)/DELTA(n,L)
         If (abs(expp) <= 60.) Then
            x=exp(expp)
            y=x+1.0/x
            faclog=log(y)
         Else
            faclog=abs(expp)
         EndIf
         temp=temp+(C(n+1,L)-C(n,L))*0.5*(Z+DELTA(n,L)*(zlog-faclog))
      EndDo

      ANTEMP=temp

!-----------------------------------------------------------------------

  END FUNCTION ANTEMP

!#######################################################################

      Subroutine fs_profile_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      End Subroutine fs_profile_init

!#######################################################################

      Subroutine fs_profile_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      End Subroutine fs_profile_end

!#######################################################################

end module fs_profile_mod


      module mcm_swnew_mod

      use mcm_swtbls_mod, only: aaa, aab
      Use       Fms_Mod, ONLY: write_version_number, mpp_pe, mpp_root_pe, &
                               error_mesg, FATAL

!implicit none 
private 

      character(len=128) :: version = '$Id: mcm_swnew.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
      character(len=128) :: tagname = '$Name: tikal $'
      logical            :: module_is_initialized = .false.

public mcm_swnew, mcm_swnew_init, mcm_swnew_end

      contains

      subroutine mcm_swnew( cosz, rco2, rh2o, ro3, pp, &
              cwca, cwcb, coca, cloudy, kthsw, kbhsw, ssolar, pr2, &
              flx, heat, grdflx, ncv, kx, UF, DF)

!          TK Original array parameters:
!              cosz, tauda, rco2, rh2o, ro3, pp, &
!              cwca, cwcb, coca, cloudy, kthsw, kbhsw, solrsw, pr2, &
!              flx, heat, grdflx, ncv, kx)

!     TK Modified from Kerr's re-write of SS shortwave routine
!     swnew.F.  FMS developmental version 5/29/01.

!     This routine performs SS shortwave radiation on a single
!     column (of dimension kx).  For kx use RDPARM_MOD.

!     rh2o and ro3 have been modified to have dimension kx, rather
!     than kx+1.

!     New output arrays UF, DF have been added as needed by FMS.

!     cosz, tauda, solrsw  trio of input variables replaced by
!     cosz and ssolar as in FMS.  Code modified below.


      parameter (nbsw=9)

! TK      real   , intent (in)  :: cosz, tauda, rco2, solrsw
      real   , intent (in)  :: cosz, rco2, Ssolar

! TK      real   , intent (in), dimension(kx+1)      :: rh2o, ro3
      real   , intent (in), dimension(kx)        :: rh2o, ro3
      real   , intent (in), dimension(0:kx)      :: pp
      real   , intent (in), dimension(1:kx+2)    :: cwca, cwcb, coca, cloudy
      real   , intent (in), dimension(1:kx+1)    :: pr2

      integer, intent (in), dimension(1:kx+2)    :: kthsw, kbhsw
      integer, intent (in)                       :: ncv, kx

      real   , intent (out)                      :: grdflx
! TK      real   , intent (out), dimension(1:kx+1)   :: flx, heat
      real   , intent (out), dimension(1:kx+1)   :: flx, heat, UF, DF

      dimension abcff (nbsw)
      dimension absdo3(kx+1)
      dimension absuo3(kx+1)
      dimension aco2  ((kx+1)*2)
      dimension adco2 (kx+1)
      dimension alfa  (nbsw)
      dimension alfat (nbsw,kx+1)
      dimension alfau (nbsw,kx+1)
      dimension ao3   ((kx+1)*2)
      dimension auco2 (kx+1)
      dimension axx   ((kx+1)*2)
      dimension ayy   ((kx+1)*2)
      dimension azz   ((kx+1)*2)
      dimension cca   (nbsw,kx+2)
      dimension ccb   (nbsw,kx+2)
      dimension cr    (nbsw,kx+2)
      dimension cro3  (kx+2)
      dimension ct    (nbsw,kx+2)
      dimension cto3  (kx+2)
      dimension ddrv  (nbsw,kx)
! TK      dimension df    (kx+1)
      dimension dfn   (nbsw,kx+1)
      dimension dpcld (kx)
      dimension dpsw  (kx+1)
      dimension duco2 (kx+1)
      dimension duo3  (kx+1)
      dimension dusw  (kx+1)
      dimension ff    (kx+1)
      dimension ffco2 (kx+1)
      dimension ffo3  (kx+1)
      dimension in    ((kx+1)*2)
      dimension ixx   ((kx+1)*2)
      dimension ppress(kx+1)
      dimension pwts  (nbsw)
      dimension t1    (nbsw,kx+1)
      dimension t2    (nbsw,kx+1)
      dimension tcld  (nbsw,kx+1)
      dimension tclu  (nbsw,kx+1)
      dimension tdh2o (nbsw,kx+1)
      dimension ttd   (nbsw,kx+1)
      dimension ttu   (nbsw,kx+1)
      dimension tuh2o (nbsw,kx+1)
      dimension uco2  ((kx+1)*2)
      dimension ud    (kx+1)
      dimension udco2 (kx+1)
      dimension udo3  (kx+1)
      dimension udrv  (nbsw,kx)
! TK      dimension uf    (kx+1)
      dimension ufn   (nbsw,kx+1)
      dimension uo3   ((kx+1)*2)
      dimension ur    (kx+1)
      dimension urco2 (kx+1)
      dimension uro3  (kx+1)
      dimension vv    (nbsw)
      dimension vvd   (nbsw)
      dimension vvu   (nbsw)

      data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190./
      data pwts/.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
      data cfco2,cfo3/508.96,466.64/
      data reflo3/1.9/
      data rrayav/0.144/
      data g1/1.020408e-3/

      kp    = kx + 1
      kpx2  = kp*2
      kpx2m = kpx2-1
      nc0   = ncv - 2
      nc1   = nc0 + 1
      nc2   = nc0 + 2
      nc3   = nc0 + 3

      ppress(1:kp) = pp(0:kx)

      do 351 i=1,kp
      ff(i)=1.66
      ffco2(i)=1.66
      ffo3(i)=1.90
351   continue
      do 339 k=1,nc2
      cro3(k)=coca(k)*cloudy(k)
339   continue
      do 334 kk=1,nc2
      cto3(kk)=1.-cro3(kk)
334   continue
      do 335 n=2,nbsw
      do 335 k=1,nc2
      cca(n,k)=cwca(k)
      ccb(n,k)=cwcb(k)
335   continue
      do 336 n=2,nbsw
      do 333 kk=1,nc2
      cr(n,kk)=cca(n,kk)*cloudy(kk)
      ct(n,kk)=cloudy(kk)*(1.-(cca(n,kk)+ccb(n,kk)))+1.-cloudy(kk)
333   continue
336   continue
      do 337 kk=1,nc2
      cr(1,kk)=cro3(kk)
      ct(1,kk)=cto3(kk)
337   continue
!     kbhsw,kthsw are cloud level indices for bottom and top of cloud
!     cr,ct are cloud reflectivity,transmissivity with index 2
!     representing the topmost cloud
!     ***********************************************************
!     calculate initial ozone reflectivity
      rray=0.219/(1.+0.816*cosz)
      rg=cro3(nc2)
      refl=rray+(1.-rray)*(1.-rrayav)*rg/(1.-rg*rrayav)
      cr(1,nc2)=refl
      ltop=kthsw(2)
      secz=1./cosz
      do 2 i=1,kx
2     dpsw(i)=ppress(i+1)-ppress(i)
      do 3 i=1,kx
      duo3(i)=ro3(i)*dpsw(i)*g1
      duco2(i)=rco2*dpsw(i)*g1
3     dusw(i)=rh2o(i)*dpsw(i)*g1
      do 131 i=1,ltop
      ffo3(i)=secz
      ffco2(i)=secz
      ff(i)=secz
131   continue
!     calculate pressure-weighted optical paths in units of g/cm2.
!     pressure weighting is by p**0.5
!     ud is the downward path,ur the upward path,
!     and the calculation is made by taking a path with an angle
!     of (secz) from the top of the atmosphere to the topmost cloud,
!     then using the diffusivity factor (1.66) to the surface and for
!     reflected radiation. the code below reflects this.
      ud(1)=0.
      udco2(1)=0.
      udo3(1)=0.
      do 4 i=2,kp
      ud(i)=ud(i-1)+dusw(i-1)*pr2(i-1)*ff(i)
      udco2(i)=udco2(i-1)+duco2(i-1)*pr2(i-1)*ffco2(i)*cfco2
      udo3(i)=udo3(i-1)+duo3(i-1)*ffo3(i)*cfo3
4     continue
!     udo3,uro3 are in units of cm. cfco2,cfo3 is the conversion
!     factor from gm/cm2 to cm.
      ur(kp)=ud(kp)
      urco2(kp)=udco2(kp)
      uro3(kp)=udo3(kp)
      do 5 i=1,kx
      ur(i)=ud(kp)+1.66*(ud(kp)/ff(kp)-ud(i)/ff(i+1)) &
                +ud(ltop)*(ff(kp)-ff(i+1))/ff(i+1)
      urco2(i)=urco2(kp)+1.66*(udco2(kp)/ffco2(kp)-udco2(i)/ffco2(i+1)) &
         +udco2(ltop)*(ffco2(kp)-ffco2(i+1))/ffco2(i+1)
      uro3(i)=uro3(kp)+reflo3*(udo3(kp)/ffo3(kp)-udo3(i)/ffo3(i+1)) &
           +udo3(ltop)*(ffo3(kp)-ffo3(i+1))/ffo3(i+1)
5     continue
 
!     maximize the size of o3 and co2 path lengths to avoid going
!     off tables

      uo3(1:kp)    = udo3(1:kp)
      uo3(kx+2:kp*2) = uro3(1:kp)

      uco2(1:kp)    = udco2(1:kp)
      uco2(kx+2:kp*2) = urco2(1:kp)

      do 621 k=1,kpx2
       uco2(k) =amin1(uco2(k),998.9)
       uo3(k) =amin1(uo3(k),3.998)
621   continue
 
!     calculate entering flux at the top
      do 6 n=1,nbsw
! TK      dfn(n,1)=solrsw*6.97667e5*cosz*tauda*pwts(n)
! TK      Replace 
      dfn(n,1)=Ssolar*6.97667e5*pwts(n)
6     continue
!     calculate water vapor transmission functions for bands 2-9;
!     t.f. for band 1= t.f for band 2
      do 7 k=1,kp
      do 7 n=2,nbsw
      t1(n,k)=amin1(abcff(n)*ud(k),50.)
      t2(n,k)=amin1(abcff(n)*ur(k),50.)
7     continue
      do 8 n=2,nbsw
      tdh2o(n,1)=1.
8     continue
      do 9 n=2,nbsw
      do 9 k=2,kp
      tdh2o(n,k)=exp(-t1(n,k))
9     continue
      do 10 n=2,nbsw
      tuh2o(n,kp)=tdh2o(n,kp)
10    continue
      do 11 n=2,nbsw
      do 11 k=1,kx
      tuh2o(n,k)=exp(-t2(n,k))
11    continue
      do 12 k=1,kp
      tdh2o(1,k)=tdh2o(2,k)
      tuh2o(1,k)=tuh2o(2,k)
12    continue
!     calculate co2 absorptions . they will be used in bands 2-9.
!     since these occupy 50 percent of the solar spectrum the
!     absorptions will be multiplied by 2.
!     the absorptions are obtained by table lookup in array aab,
!     common block swtabl,and then interpolation.
      do 614 k=1,kpx2
      ixx(k)=uco2(k)+1.
614   continue
      do 615 k=1,kpx2
      axx(k)=uco2(k)-ifix(uco2(k))
615   continue
      call mcm_sif1d ( ixx, kpx2m, ayy, azz )
      do 617 k=1,kpx2
      aco2(k)=ayy(k)+axx(k)*(azz(k)-ayy(k))
617   continue

      adco2(1:kp) = aco2(1:kp)
      auco2(1:kp) = aco2(kx+2:kp*2)

      do 26 k=1,kp
      adco2(k)=2.*adco2(k)
      auco2(k)=2.*auco2(k)
26    continue
!     now calculate ozone absorptions. these will be used in
!     band 1. as this occupies 50 percent of the solar spectrum
!     the ozone absorptions will be multiplied by 2.
!     the ozone absorptions are obtained by table lookup from
!     array aaa, common block swtabl. no interpolation is done.
      do 603 k=1,kpx2
      in(k)=500.*uo3(k)+1
603   continue
      call mcm_sif1 ( in, kpx2m, ao3 )

      absdo3(1:kp) = ao3(1:kp)
      absuo3(1:kp) = ao3(kx+2:kp*2)

      do 33 k=1,kp
      absdo3(k)=2.*absdo3(k)
      absuo3(k)=2.*absuo3(k)
33    continue
!     combine absorptions and transmissions to obtain a
!     transmission function for each of the 9 bands.
      do 41 n=1,nbsw
      ttd(n,1)=1.
41    continue
      do 42 k=2,kp
      ttd(1,k)=tdh2o(1,k)*(1.-absdo3(k))
42    continue
      do 43 k=1,kx
      ttu(1,k)=tuh2o(1,k)*(1.-absuo3(k))
43    continue
      do 44 n=2,nbsw
      do 44 k=2,kp
      ttd(n,k)=tdh2o(n,k)*(1.-adco2(k))
44    continue
      do 45 n=1,nbsw
      ttu(n,kp)=ttd(n,kp)
45    continue
      do 46 n=2,nbsw
      do 46 k=1,kx
      ttu(n,k)=tuh2o(n,k)*(1.-auco2(k))
46    continue
!     the following calculation is for alfat: the ratio between
!     the downward flux at the top of the highest cloud (if
!     present) or the ground to the upward flux at the same
!     level, taking into account multiple reflections from
!     clouds, if present
      do 53 n=1,nbsw
      do 53 nn=1,nc1
      alfat(n,nn)=cr(n,nc3-nn)
53    continue
      do 55 nn=1,nc1
      do 55 n=1,nbsw
      alfau(n,nn)=0.
55    continue
      if ( nc0 .eq. 0 )  go to 58
      do 51 ll=1,nc1
      kt1=kthsw(nc3-ll)
      kt2=kthsw(nc2-ll)
      kt3=kbhsw(nc2-ll)
!     tclu(ll) is transmission function from top of next lower
!     cloud to top of upper cloud. tcld(ll) is t.f. from top
!     of next lower cloud to bottom of upper cloud. ll=1 is
!     the lowest bottom cloud (the ground) ; ll=nc1 is the
!     highest upper cloud. this calculation is used for the
!     alfat calculation only if cloud is present.
      do 52 n=1,nbsw
      tclu(n,ll)=ttd(n,kt1)/ttd(n,kt2)
      tcld(n,ll)=ttd(n,kt1)/ttd(n,kt3)
52    continue
51    continue
      do 56 nn=1,nc0
      do 57 n=1,nbsw
      alfau(n,nn+1)=(cr(n,nc3-nn)+alfau(n,nn))*ct(n,nc2-nn)* &
           ct(n,nc2-nn)*tclu(n,nn)*tclu(n,nn)/ &
           (1.-(cr(n,nc3-nn)+alfau(n,nn))*cr(n,nc2-nn)* &
           tcld(n,nn)*tcld(n,nn))
57    continue
56    continue
58    continue
      do 59 n=1,nbsw
      alfa(n)=alfat(n,nc1)+alfau(n,nc1)
59    continue
!     downward flux above topmost cloud
      do 61 n=1,nbsw
      vv(n)=dfn(n,1)
61    continue
      do 62 n=1,nbsw
      do 62 k=2,ltop
      dfn(n,k)=vv(n)*ttd(n,k)
62    continue
!     upward flux above topmost cloud
      ltopm=ltop-1
      do 63 n=1,nbsw
      ufn(n,ltop)=alfa(n)*dfn(n,ltop)
63    continue
      do 64 n=1,nbsw
      vv(n)=ufn(n,ltop)/ttu(n,ltop)
64    continue
      do 65 n=1,nbsw
      do 65 k=1,ltopm
      ufn(n,k)=vv(n)*ttu(n,k)
65    continue
      if ( nc0 .eq. 0 )  go to 91
!     calculate ufn at cloud tops and dfn at cloud bottoms
      do 71 ll=2,nc1
      do 72 n=1,nbsw
      ufn(n,kthsw(ll+1))=alfau(n,nc3-ll)*dfn(n,kbhsw(ll-1))* &
           tcld(n,nc3-ll)/(tclu(n,nc2-ll)*ct(n,ll))
72    continue
      do 73 n=1,nbsw
      dfn(n,kbhsw(ll))=dfn(n,kbhsw(ll-1))*tcld(n,nc3-ll)*tclu(n,nc2-ll)* &
           ct(n,ll)/tcld(n,nc2-ll)+ufn(n,kthsw(ll+1))*tcld(n,nc2-ll)* &
           cr(n,ll)
73    continue
71    continue
!     now obtain dfn and ufn for levels between the clouds
      do 74 ll=1,nc0
      ltop=kbhsw(ll+1)
      ltopp=ltop+1
      lbot=kthsw(ll+2)
      lbotm=lbot-1
      if (ltop.eq.lbot) go to 74
      do 75 n=1,nbsw
      vvu(n)=ufn(n,lbot)/ttu(n,lbot)
      vvd(n)=dfn(n,ltop)/ttd(n,ltop)
75    continue
      do 76 k=ltop,lbotm
      do 76 n=1,nbsw
      ufn(n,k)=vvu(n)*ttu(n,k)
76    continue
      do 77 k=ltopp,lbot
      do 77 n=1,nbsw
      dfn(n,k)=vvd(n)*ttd(n,k)
77    continue
74    continue
!     now obtain downward and upward fluxes for levels,if any,
!     between the tops and bottoms of clouds. the assumption of
!     constant heating rate in these regions is used.
      do 78 ll=1,nc0
      ltop=kthsw(ll+1)
      lbot=kbhsw(ll+1)
      if((lbot-ltop).le.1) go to 78
      ltopp=ltop+1
      lbotm=lbot-1
      dpcld(ll)=ppress(ltop)-ppress(lbot)
      do 79 n=1,nbsw
      udrv(n,ll)=(ufn(n,ltop)-ufn(n,lbot))/dpcld(ll)
      ddrv(n,ll)=(dfn(n,ltop)-dfn(n,lbot))/dpcld(ll)
79    continue
      do 80 n=1,nbsw
      vvu(n)=ufn(n,ltop)
      vvd(n)=dfn(n,ltop)
80    continue
      do 81 k=ltopp,lbotm
      do 81 n=1,nbsw
      ufn(n,k)=vvu(n)+udrv(n,ll)*(ppress(k)-ppress(ltop))
      dfn(n,k)=vvd(n)+ddrv(n,ll)*(ppress(k)-ppress(ltop))
81    continue
78    continue
91    continue
!     sum over bands
 
      do 14 k=1,kp
      df(k)=0.
      uf(k)=0.
      do 15 n=1,nbsw
      df(k)=df(k)+dfn(n,k)
      uf(k)=uf(k)+ufn(n,k)
15    continue
14    continue
      do 16 k=1,kp
      flx(k)=uf(k)-df(k)
16    continue
! TK REMOVE: cdir$ novector
      do 17 k=1,kx
      heat(k)=8.42668*(flx(k+1)-flx(k))/dpsw(k)
17    continue
! TK REMOVE: cdir$ vector
!     8.42668=g(cgs units)*(no.sec/da)/cp(cgs units)
      grdflx=(1.-refl)*dfn(1,kp)
      do 19 n=2,nbsw
      grdflx=grdflx+ct(1,nc2)*dfn(n,kp)
19    continue
      return
      end subroutine mcm_swnew
! --------------------------------------------
      subroutine mcm_sif1(ndx,len,tgt)
      dimension ndx(1),tgt(1)
      leng=len+1
      do 10 j=1,leng
        jj=ndx(j)
        tgt(j)=aaa(jj)
  10  continue
      return
      end subroutine mcm_sif1
! --------------------------------------------
      subroutine mcm_sif1d(ndx,len,tgt1,tgt2)
      dimension ndx(1),tgt1(1),tgt2(1)
      leng=len+1
      do 10 j=1,leng
        jj=ndx(j)
        tgt1(j)=aab(jj)
        tgt2(j)=aab(jj+1)
  10  continue
      return
      end subroutine mcm_sif1d

! ---------------------------------------------------------------------------------------
      subroutine mcm_swnew_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      end subroutine mcm_swnew_init
! ---------------------------------------------------------------------------------------
      subroutine mcm_swnew_end

      module_is_initialized = .false.

!---------------------------------------------------------------------
      end subroutine mcm_swnew_end


      end module mcm_swnew_mod

module MO_FASTJX_MOD
!----------------------------------------------------------------------
!       FAST-JX Photolysis rate 
!       Revised by Junfeng Liu
!       
!       Feb 19, 2010
!       Junfeng.Liu
!----------------------------------------------------------------------
!         Below is the original information from FASTJX code
!
!  >>>>>>>>>>>>>>>>current code revised to JX ver 6.4 (8/08)<<<<<<<<<<<<
!                                                                       
! version 6.4                                                           
! >>>> allows for shortened, speeded up troposphere versions<<<<<<      
!     STD:           W_=18                                              
!        identical results to v-6.2 if cloud OD is consistent           
!     TROP-ONLY:     W_=12                                              
!        collapses the wavelength bins from 18 to 12 (5-6-7-8 & 11-18)  
!        drops many 'stratospheric' cross-sections (denoted by 'x' in 2n
!        allows use of single standard spectral data set:  FJX_spec.dat 
!        results close to W_=18, largest difference is J-O2 (<1% in 13-1
!        This is recommended as accurate for troposphere only calculatio
!     TROP-QUICK:    W_=8                                               
!        reverts to original fast-J 7-bins (12-18) plus 1 scaled UV (5) 
!        errors in 12-18 km range for J-O2, high sun are 10%, worse abov
!     ***Photolysis of O2 in the upper tropical troposphere is an import
!        source of O3.  It needs to be included in tropospheric runs.   
!        TROP-ONLY is recommended, W_=8 is a quick fix if speed essentia
!                                                                       
!     Major rewrite of code to minimize calls and allow better vector-ty
!     loop over wavelengths internal to Mie soln.                       
!     Driven by profiling of CTM code, may still need optimization.     
!     Wavelengths can be contracted to W_=12 (trop only) and strat-only 
!        X-sections are dropped.  With parm W_=18, the std fast-JX is re
!     Many call eliminated and summations in BLKSLV and GEN_ID are expli
!     GEN_ID replaces GEN and calculates all matrix coeff's (1:L_) at on
!     RD_XXX changed to collapse wavelengths & x-sections to Trop-only: 
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengt
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).     
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...
!                                                                       
! version 6.3                                                           
!     revise cloud/aerosol OD & wavelength properties for CTM link:     
!         OPTICL is new sub for cloud optical properties, but it        
!              now starts with cloud OD @ 600 nm and cloud NDX          
!              fast-JX now uses cloud NDX to scale OD to other wavelengt
!         OPTICA & OPTICM are new subs to convert aerosol path (g/m2) to
!              A is std UCI scat data                                   
!              M is U Michigan data tables for aerosols, includes Rel Hu
!     drop sub GAUSSP and put into Parameter statement (parm_MIE.f)     
!                                                                       
! version 6.2                                                           
!     corrects a long-standing problem at SZA > 89 degrees.             
!     In prior versions the ray-tracing of the path (and air-mass functi
!     back to the sun was done at the edges of the CTM layers (it was de
!     for the grid-point J-value code at Harvard/GISS/UCI).  This left t
!     interpolation to the mid-layer (needed for J's) open.  The prior m
!     gave irregular fluctuations in the direct solar beam at mid-layer 
!     large SZA > 88.  This is now corrected with exact ray-tracing from
!     the mid-pt of each CTM layer.  For small SZA, there is no effectiv
!     difference, for large SZA, results could be erratic.              
!                                                                       
!   v-6.2 fix should be easy if you have migrated to v6.1, else some min
!      caution may be needed:                                           
!      replace sub SPHERE with SPHERE2, AMF2 report factors for mid and 
!      replace sub OPMIE with new OPMIE, this uses the new AMF2 correctl
!      replace sub PHOTOJ with new PHOTOJ, this just hands off AMF2 from
!            SPHERE2 to OPMIE.                                          
!                                                                       
! version 6.1 adds                                                      
!      6.1b simplifies calling sequences feeds solar factor, albedo, to 
!         and read LAT, LNG directly.  No substantive changes.          
!      new read-in of scat data for clouds/aerosols to allow for UMich d
!      This has required substantial rewrite of some of the core subrout
!         OPMIE is now called for each wavelength and without aersol/clo
!              all subs below OPMIE are unchanged                       
!         OPTICD & OPTICM are new subs to convert path (g/m2) to OD and 
!              D is std UCI scat data (re-ordered for clouds 1st)       
!              M is U Michigan data tables for aerosols, includes Rel Hu
!         PHOTOJ now assembles the aerosol data (better for CTM implemen
!      This version can reproduce earlier versions exactly, but the test
!         is changed from OD and NDX to PATH (g/m2) and NDX.            
!                                                                       
! version 6.0 adds                                                      
!      new 200-nm scattering data so that stratospheric aerosols can be 
! version 5.7                                                           
!     adds the new flux diagnostics (including heating rates)           
!        accurate fluxes for spherical atmos and SZA > 90 !             
!     recommend geometric delta-tau factor from 1.18 to 1.12 for more ac
!        heating rates (but more layers!)                               
!     tuned and corrected to be almost flux conserving (1.e-5), except  
!        deep clouds, where diffusive flux is created (1.e-4)           
!     still needs to return to the original 1970-code for the block-tri 
!        after extensive profiling with F95 and 'modern' versions       
!        it was found that they are much more expensive!!!              
!     corrects typo in JAC(2000) fast-J paper on I+ (reflected from l.b.
!        I+(lb) = refl/(1+refl) * (4*Integ[j(lb)*mu*dmu] + mu0*Fdirect(l
! version 5.6 adds                                                      
!      clean up problems with thick clouds does correct solar attenuatio
!        into cloud sub-layers and into the mid-point of the CTM level  
!      New calculated upward and downward FLUXES at each wavelength at T
!      Correct deposition of solar flux in each CTM layer (spherical)   
!        awaits new diagnostics of the h's for heating rates.           
!      back to old matrix solver (UCI blocksolver and matinv-4)         
! version 5.5 adds                                                      
!      new code for generating and solving the block tri-diagonal scatte
!           problem.  Uses single call to GEM and general 4x4 block-tri 
! version 5.3c adds                                                     
!      calculates reflected UV-vis solar energy (relative to 1.0)       
!      new solar spectrum (J-O2 increases in strat by 10%, J-NO by 15+%)
!                                                                       
! version 5.3b changes include:                                         
!      new data files for specral Xsection and mie-scattering.          
!      add sub-layers (JXTRA) to thick cloud/aerosol layers,            
!           sets up log-spaced sub-layers of increasing thickness ATAU  
!      correction 'b' does massive clean up of the linking code,        
!           now the only subroutine that has access to CTM arrays is PHO
!           Also, the access to the cmn_JVdat.f is 'read-only' after ini
!           This should enable safe openMP/MPI coding.                  
!                                                                       
! common files and what they mean:                                      
!   parm_CTM.f  dimensions & params for code (CTM and fast-JX)          
!   parm_MIE.f  dimensions for mie code variables.                      
!   cmn_metdat.f  CTM 3-D arrays, time of day, grid,  etc.              
!   cmn_JVdat.f   Xsects, Mie, etc., (initialized and then read-only)   
!                                                                       
!<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
! subroutines:                                                          
!                                                                       
!   SET_ATM(GMTAU)                                                      
!           set ups atmosphere (p,T,O3,airmass, etc) for time GMTAU     
!              COMMON BLOCKS: cmn_metdat.f                              
!                                                                       
!   SET_AER(GMTAU)                                                      
!              set ups aerosols for time GMTAU = DUMMY                  
!                                                                       
!   SET_CLD(GMTAU)                                                      
!           set ups clouds for time GMTAU = DUMMY                       
!                                                                       
!   INPHOT                                                              
!           Init. photolysis rate, called once by CHMSET                
!              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f                 
!              Input files: ECT42_grid.dat                              
!                                                                       
!   RD_JS(NJ1,NAMFIL)                                                   
!           Read labels of photo. rates, called once by INPHOT.         
!              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f                 
!              Input files: chem_Js.dat                                 
!                                                                       
!   RD_PROF(NJ2,NAMFIL)                                                 
!           Read T & O3 climatology, called once by INPHOT.             
!              COMMON BLOCKS: cmn_metdat.f                              
!              Input files: atmos_std.dat                               
!                                                                       
!   SET_CLD0(TINIT,CLDP,NCLD,ALBEDO)                                    
!              Initialize cloud and surface properties, called by MAIN. 
!              COMMON BLOCKS: cmn_metdat.f                              
!                                                                       
!   SET_AER0 (AER1,AER2, NAA1,NAA2)                                     
!              Iniitalize (climatology) aerosol OD and types (3 arrays) 
!              COMMON BLOCKS: cmn_metdat.f                              
!                                                                       
!   SET_ATM0                                                            
!             Initialize climatologies for T & O3, set up atmospheric pr
!              COMMON BLOCKS: cmn_metdat.f                              
!                                                                       
!<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
!                                                                       
!   PHOTOJ(UTIME,IDAY,ILNG,JLAT,SOLF,SZA,U0,FREFL,ZPJ)                  
!              Gateway to fast-JX, Update the photolysis rates          
!              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f                 
!                                                                       
!<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
!  N.B. all these need access to cmn_JVdat.f, but do NOT write into it. 
!           also have no need to access cmn_metdat.f                    
!                                                                       
!   OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)                               
!        UCI CLOUD data sets, calculate scattering properties,          
!           scales OD at 600 nm to JX wavelength                        
!              COMMON BLOCKS: cmn_JVdat.f                               
!                                                                       
!   OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L)                               
!   OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)                               
!         UC Irvine & U Michigan aerosol data sets, generates fast-JX da
!              COMMON BLOCKS: cmn_JVdat.f                               
!                                                                       
!   JRATET(PPJ,TTJ,FFF, VALJL)                                          
!           Calculate J-value, called by PTOTOJ.                        
!              COMMON BLOCKS: cmn_JVdat.f                               
!                                                                       
!   JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,DTAUX,POMEGAX,JXTRA)                
!              print out atmosphere used in J-value calc.               
!              >>>will be superseded by CTM routines                    
!              COMMON BLOCKS: cmn_JVdat.f                               
!                                                                       
!   RD_XXX(NJ1,NAMFIL)                                                  
!           Read wavelength bins, solar fluxes, Rayleigh, Temp-dep X-sec
!           >>>v-6.4 now collapse to Trop-only Wl's and Xs's if W_=12   
!              called once by INPHOT.                                   
!              COMMON BLOCKS: cmn_JVdat.f                               
!              Input files: FJX_spec.dat                                
!                                                                       
!   RD_MIE(NJ1,NAMFIL)                                                  
!           Set aerosols/cloud scattering, called once by INPHOT        
!              COMMON BLOCKS: cmn_JVdat.f                               
!              Input files: FJX_scat.dat                                
!                                                                       
!   RD_UM(NJ1,NAMFIL)                                                   
!           UMich aerosol optical data, called once by INPHOT           
!              COMMON BLOCKS: cmn_JVdat.f                               
!              Input files: FJX_UMaer.dat                               
!                                                                       
!   SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)                   
!     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)                 
!              calc SZA and Solar Flux factor for given lat/lon/UT      
!                                                                       
!   SPHERE2(GMU,RAD,ZHL,ZZHT,AMF2,L1_)                                  
!              calculate spherical geometry, air-mass factors (v 6.2)   
!                                                                       
!   EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)                   
!              add sub-layers (JXTRA) to thick cloud/aerosol layers     
!                                                                       
!   OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,FJACT,FJTOP,FJBOT,FSBOT,FJFLX
!              calculate mean intensity (actinic) at each CTM levels    
!              calculate fluxes and deposition (heating rates)          
!              COMMON BLOCKS: cmn_JVdat.f                               
!                                                                       
!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!                                                                       
!   MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)                  
!            include 'parm_MIE.f' = dimension parameters                
!                                                                       
!   LEGND0 (X,PL,N)                                                     
!                                                                       
!   BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)          
!                                                                       
!   GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)            
!                                                                       
!   MATINW (B,A)                                                        
!                                                                       
!-----------------------------------------------------------------------
      use              fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       mpp_clock_begin, mpp_clock_end, &
                                       mpp_clock_id, CLOCK_MODULE, &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL
      use           mpp_io_mod, only : mpp_open, mpp_close, MPP_RDONLY, &
                                       MPP_ASCII, MPP_SEQUENTIAL,   &
                                       MPP_MULTI, MPP_SINGLE 

      implicit none
      
      private
      public :: fastjx_init, fastjx_end, fastjx_photo, JVN_

!-----------------------------------------------------------------------
!     version number and tagname.
!-----------------------------------------------------------------------
      character(len=128)            :: version     = '$Id: mo_fastjx.F90,v 20.0 2013/12/13 23:24:59 fms Exp $'
      character(len=128)            :: tagname     = '$Name: tikal $'

!    include 'parm_CTM.f'  for fast-JX code v5.3+ (prather 6/05)
!
!     I_ = longitude dim of CTM grid
!     J_ = latitude  dim of CTM grid
!     L_ = altitude(levels) dim of CTM grid
!     LWE_ = altitude(level) dim for trop processes (clouds, rain)
!     JVL_ = vertical(levels) dim for J-values
!     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_
!     JVN_ =  no. of J-values
!     W_   = dim = no. of Wavelength bins   (now in parm_MIE.f)
!     WX_  = dim = no. of wavelengths in input file
!     X_   = dim = no. of X-section data sets (input data)
!     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
!     MX   = no. of aerosol/cloud types supplied from CTM
!     NTR_ = no. of CTM tracers
!     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
!
!-----------------------------------------------------------------------
!      integer, parameter               ::      I_      =128            !longitude dim of CTM grid
!      integer, parameter               ::      J_      =64             !latitude  dim of CTM grid
      integer, parameter                ::      L_      =48             !altitude(levels) dim of CTM grid
!      integer, parameter               ::      LWE_    =37             !altitude(level) dim for trop processes (clouds, rain)
      integer, parameter                ::      JVL_    =48             !vertical(levels) dim for J-values
      integer, parameter                ::      L1_     =L_+1           !assuming another layer above the top layer
      integer, parameter                ::      L2_     =2*L_+2         !2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_ 
!      integer, parameter               ::      MX      =4              !no. of aerosol/cloud types supplied from CTM, no use
!      integer, parameter               ::      NTR_    =1              !no. of CTM tracers, is O3 here
      real*8,  parameter                ::      SZAMAX  =98.0d0         !Solar zenith angle cut-off, above which to skip calculation
      integer, parameter                ::      WX_     =18             !dim = no. of wavelengths in input file
      integer, parameter                ::      X_      =64             !dim = no. of X-section data sets (input data)
      integer, parameter                ::      A_      =40             !dim = no. of Aerosol/cloud Mie sets (input data)
      integer, parameter                ::      A_AM3   =1000           !dim = no. of Aerosol/cloud Mie sets FOR AM3(input data)     
      integer, parameter                ::      JVN_    =62             !no. of J-values
!----------------------------------------------------------------------- 
!    include 'parm_mie.f'  for fast-JX code v5.3+ (prather 6/05)
!
!     N_  = no. of levels in Mie scattering arrays
!         = 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
!     M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
!     M2_ = 2*M_ = 8, replaces MFIT
!     W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
!
!-----------------------------------------------------------------------
      integer, parameter                ::      N_      =501            !no. of levels in Mie scattering arrays
                                                                        != 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
      integer, parameter                ::      M_      =4              !no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter                ::      M2_     =2*M_           !2*M_ = 8, replaces MFIT
      integer, parameter                ::      W_      =18             !dim = no. of Wavelength bins:  =18 std, =12 trop only
!-----------------------------------------------------------------------
!    4 Gauss pts = 8-stream
      real*8, dimension(M_), parameter  ::      EMU = [.06943184420297d0, &
                                                        .33000947820757d0, .66999052179243d0, .93056815579703d0]
      real*8, dimension(M_), parameter  ::      WT  = [.17392742256873d0, &
                                                        .32607257743127d0, .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------
!    include 'cmn_JVdat.f'  for fast-JX code v 6.0+ (prather 4/07)
! 
! NB - ALL of these common variables are set paramters, 
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
!-----INPHOT------------------------------------------------------------ 
      real*8                            ::      RAD                     !Radius of Earth (cm)  
      real*8                            ::      ZZHT                    !Effective scale height above top of atmosphere (cm)
!-----RD_XXX  <- FJX_spec.dat-------------------------------------------
      real*8                            ::      WBIN(WX_+1)             !Boundaries of wavelength bins !!!NO USE
      integer                           ::      NJVAL                   !Number of J values
      integer                           ::      NW1                     ! is 1
      integer                           ::      NW2                     ! is 18 wavelines
      real*8                            ::      WL(WX_)                 !Centres of wavelength bins - 'effective wavelength' (nm)
      real*8                            ::      FL(WX_)                 !Solar flux incident on top of atmosphere (cm-2.s-1)
      real*8                            ::      QRAYL(WX_+1)            !Rayleigh parameters (effective cross-section) (cm2) 
      real*8                            ::      QO2(WX_,3)              !O2 cross-sections 
      real*8                            ::      QO3(WX_,3)              !O3 cross-sections
      real*8                            ::      Q1D(WX_,3)              !O3 => O(1D) quantum yield
      character*7                       ::      TITLEJ(X_)              !Species name in FJX_spec.dat
      character*7                       ::      TITLEJ2                 !Temp var 
      character*7                       ::      TITLEJ3                 !Temp var      
      real*8                            ::      TQQ(3,X_)               !Temperature for supplied cross sections(K)
      real*8                            ::      QQQ(WX_,2,X_)           !Supplied cross sections in each wavelength bin (cm2)
!-----RD_MIE  <- FJX_scat.dat-------------------------------------------
      integer                           ::      NAA                     !Number of categories for scattering phase functions, is 30
      integer                           ::      NAA_AM3                 !Number of categories for scattering phase functions for am3, is 880
      real*8                            ::      MEE_AM3(5,A_AM3)        !Aerosol mass extinction efficiency, MEE*colume mass =od
      real*8                            ::      SAA_AM3(5,A_AM3)        !Single scattering albedo
      real*8                            ::      PAA_AM3(8,5,A_AM3)      !Phase function: first 8 terms of expansion      
      real*8                            ::      WAA_AM3(5,A_AM3)        !Wavelengths for the NK supplied phase functions(nm)
      character*78                      ::      TITLE0                  !Title of the first line in FJX_spec.dat
      character*78                      ::      TITLE0_AM3              !Title of the first line in am3_scat.dat
      integer                           ::      JTAUMX                  !Max # is 250
      real*8                            ::      ATAU                    ! is 1.120
      real*8                            ::      ATAU0                   ! is 0.010
      character*20                      ::      TITLAA(A_)              ! A_ =40
      real*8                            ::      RAA(A_)                 !Effective radius associated with aerosol type
      real*8                            ::      DAA(A_)                 ! rho
      real*8                            ::      WAA(5,A_)               !Wavelengths for the NK supplied phase functions(nm)
      real*8                            ::      QAA(5,A_)               !Aerosol scattering phase functions
      real*8                            ::      SAA(5,A_)               !Single scattering albedo
      real*8                            ::      PAA(8,5,A_)             !Phase function: first 8 terms of expansion      
!-----RD_UM  <-  FJX_UMaer.dat------------------------------------------  
      character*20                      ::      TITLUM(33)              !UM aerosol type names
      real*8                            ::      UMAER(3,5,21,33)        !UM aerosol data
!-----RD_JS  <-  chem_Js.dat--------------------------------------------
      real*8                            ::      JFACTA(JVN_)            !Quantum yield (or multiplication factor) for photolysis
      character*7                       ::      JLABEL(JVN_)            !Reference label identifying appropriate J-value to use
      integer                           ::      JIND(JVN_)              !Set index arrays that map Jvalue(j) onto rates
      integer                           ::      NRATJ                   !number of J species
      logical                           ::      module_is_initialized = .false.      
      real*8,    parameter              ::      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)         !used to calc ZH
      real*8,    parameter              ::      pfactor1 = 1.e3/9.8/28.9644*6.02e23/1.e4                !pa => molec/cm2
      real*8, dimension(L1_)            ::      ZH   !
      

      CONTAINS
! <SUBROUTINE NAME="fastjx_init">
!   <OVERVIEW>
!     Initialize data for fastjx calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the calculation of photolysis rates
!     for FAST-JX calculation
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fastjx_init
!   </TEMPLATE>                                                                        
!-----------------------------------------------------------------------
      subroutine fastjx_init 
!-----------------------------------------------------------------------
!  Routine to initialise photolysis rate data, called directly from the 
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction in
!-----------------------------------------------------------------------
!                                                                       
!     IPH       Channel number for reading all data files               
!     RAD       Radius of Earth (cm)                                    
!     ZZHT      Effective scale height above top of atmosphere (cm)     
!     DATUMX    Maximum opt.depth above which sub layers should be inser
!     SZAMAX    Solar zenith angle cut-off, above which to skip calculat
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_metdat.f' 
!      include 'cmn_JVdat.f' 
                                                                        
       integer  I, J, K 
!      real*8, parameter::   PI = 3.141592653589793d0 
!      real*8, parameter::  PI180=3.141592653589793d0/180.d0 
                                                                        
      if (module_is_initialized) return
                           
!----------------------------------------------------------------------
!     output version number and tagname to logfile.
!----------------------------------------------------------------------
      call write_version_number (version, tagname)

!                                                                       
! Defaults & constants                                                  
      RAD  = 6375.d5   !cm earth radii
      ZZHT = 5.d5 
!      STT(:,:,:,:) = 1.d6 
!      TCNAME(:)    = 'CO2' 
!                                                                       
! Read in annual mean pressure field                                    
!RSH if this code is ever reactivated, change to use mpp_open
!      open (IPH,file='ECT42_grid.dat',form='formatted',status='old') 
!      read (IPH,*) 
!      read (IPH,'(16F8.2)') (XDGRD(I),I=1,I_) 
!      read (IPH,*) 
!      read (IPH,'(16F8.2)') (YDGRD(J),J=1,J_) 
!      read (IPH,*) 
!      read (IPH,'(16F8.2)') (AREAXY(1,J),J=1,J_) 
!      read (IPH,*) 
!      read (IPH,'(16F8.1)') ((PMEAN(I,J),I=1,I_),J=1,J_) 
!      close(IPH) 
!RSH   call mpp_close (IPH)
                                                                        
!      do J=1,J_ 
!        AREAXY(1,J) = AREAXY(1,J)*1.d9 
!        do I=2,I_ 
!        AREAXY(I,J) = AREAXY(1,J) 
!        enddo 
!      enddo 
                                                                        
!      do I = 1,I_ 
!        XGRD(I) = XDGRD(I)*PI180 
!      enddo 
!      do J = 1,J_ 
!        YGRD(J) = YDGRD(J)*PI180 
!      enddo 
!                                                                       
! Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX  
      call RD_XXX('INPUT/FJX_spec.dat') 
                                                                        
! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX  
      call RD_MIE('INPUT/FJX_scat.dat') 

! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX  
      call RD_MIE_AM3('INPUT/am3_scat.dat') 
      
! Read in UMich aerosol optical data    <<<<<<<<<<<<<<<<<< new fast-JX 6
!      call RD_UM (IPH,'INPUT/FJX_UMaer.dat') 
                                                                        
! Read in labels of photolysis rates required   >>>>> keyed to users che
!   this is a tranfer map from the J's automatically calculated in fast-
!   onto the names and order in the users chemistry code                
      call RD_JS('INPUT/chem_Js.dat') 
!                                                                       
! Read in T & O3 climatology                    >>>> general backup clim
!      call RD_PROF(IPH,'atmos_std.dat')      

!-----------------------------------------------------------------------
!   mark module as initialized.
!-----------------------------------------------------------------------
         module_is_initialized = .true.                                                                     
    end subroutine fastjx_init                                         
 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    include 'cmn_metdat.f'  for fast-JX code v5.3+ (prather 6/05)
!
!        needs 'parm_ctm.f' for dimensions
!        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJX
!        >>>>this is for standalone fast-JX ver 5.3   (6.05)
!-----------------------------------------------------------------------
!
!      real*8  P(I_,J_)         !  Surface pressure
!      real*8  T(I_,J_,L_)      !  Temperature profile
!
!      real*8  XGRD(I_)         !  Longitude (midpoint, radians)
!      real*8  XDGRD(I_)
!      real*8  YGRD(J_)         !  Latitude  (midpoint, radians)
!!      real*8  YDGRD(J_) 
!      real*8  ETAA(L1_)        !  Eta(a) value for level boundaries
!      real*8  ETAB(L1_)        !  Eta(b) value for level boundaries
!      real*8  AREAXY(I_,J_)    !  area (m^2)
!      integer  MONTH
!      integer  NSLAT         ! Latitude(J) index of current column
!      integer  NSLON         ! Longitude(I) index of current column

!   Note that the climatology for O3, T, clouds, etc include a layer above
!        the CTM (L1_=L_+1) for calculating the J's
!      real*8, dimension(L1_)       :: RELH   !  Rel Hum (0.00 to 1.00)
!      real*8, dimension(I_,J_,L1_) :: TJ, DM, DO3, ZH
!      real*8, dimension(I_,J_,L1_) :: CLDLWP,AER1P, AER2P
!      integer,dimension(I_,J_,L1_) :: CLDNDX,AER1N, AER2N
!      real*8, dimension(I_,J_)     :: PMEAN, SA



!     PJ       Pressure at boundaries of model levels (hPa)             
!     MASFAC   Conversion factor for pressure to column density         
!                                                                       
!     TJ       Temperature profile on model grid                        
!     DM       Air column for each model level (molecules.cm-2)         
!     DO3      Ozone column for each model level (molecules.cm-2)       
!     ZH       Altitude of boundaries of model levels (cm)              
!     PSTD     Approximate pressures of levels for supplied climatology 

!
!      real*8   STT(I_,J_,L_,NTR_)
!      real*8   TREF(51,18,12),OREF(51,18,12)
!      character*10  TCNAME(NTR_)
!
!      common/metdat/P,T, STT,  &
!     &  XGRD,YGRD,XDGRD,YDGRD,ETAA,ETAB,AREAXY,  &
!     &  TREF,OREF,PMEAN,SA, MONTH,NSLAT,NSLON, TCNAME
!
!      common /jvdatIJ/TJ,DM,DO3,ZH,RELH, CLDLWP,AER1P,AER2P, &
!     &                  CLDNDX,AER1N,AER2N
!-----------------------------------------------------------------------


      
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)       
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<< 
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     NJ1      Channel number for reading data file                     
!     NAA      Number of categories for scattering phase functions      
!     QAA      Aerosol scattering phase functions                       
!     NK       Number of wavelengths at which functions supplied (set as
!     WAA      Wavelengths for the NK supplied phase functions          
!     PAA      Phase function: first 8 terms of expansion               
!     RAA      Effective radius associated with aerosol type            
!     SAA      Single scattering albedo                                 
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J
!     NJ1      Channel number for reading data file                     
!                                                                       
!     NJVAL    Number of species to calculate J-values for              
!     NWWW     Number of wavelength bins, from 1:NWWW                   
!     WBIN     Boundaries of wavelength bins                            
!     WL       Centres of wavelength bins - 'effective wavelength'      
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)      
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)      
!     QO2(18,3)      O2 cross-sections                                        
!     QO3(18,3)      O3 cross-sections                                        
!     Q1D(18,3)      O3 => O(1D) quantum yield                                
!     TQQ(3,64)      Temperature for supplied cross sections                  
!     QQQ      Supplied cross sections in each wavelength bin (cm2)   
!  
!-----------------------------------------------------------------------
!      common /jvchem/JFACTA,JIND,NRATJ,JLABEL
!
!      common /jvdat/WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,  &
!     &        ZZHT,ATAU,ATAU0, WAA,QAA,PAA,SAA,RAA,DAA,RAD, UMAER,  &
!     &        JTAUMX, NJVAL,NW1,NW2,NAA ,TITLE0,TITLEJ,TITLAA,TITLUM
!-----------------------------------------------------------------------



!     fastjx_photo should be called by fphoto, its function is to provide 
!                  the calculated J-values to fphoto, based on clould and 
!                  aerosol information
!
                  
      subroutine fastjx_photo( U0, SOLF, phalf1, zhalf1, pfull1, tfull, &
                              XO3, pwt, lwc, cloudsf, clouds_ndx ,qfld, &
                              srf_alb, aerop, aeron, o3_column_top, ZPJ)                   

!-----------------------------------------------------------------------
!                                                                       
!          PHOTOJ is the gateway to fast-JX calculations:                       
!        only access to CTM 3-D GLOBAL arrays                           
!        sets up the 1-D column arrays for calculating J's              
! v6.1   calculates the optical properties for each wavelength          
!          and passes single column, single wavelength data to OPMIE    
!        OPMIE no longer needs aerosol data, just the overall properties
!           DTAUX() = opt depth of each layer,                          
!           POMEGAX(8,) = scat phase fn, includes s-s albedo factor.    
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_metdat.f' 
!      include 'cmn_JVdat.f' 
      real*8,           intent(in)      ::  U0,      &                  !cos(SZA)
                                            SOLF                        !solar-earth distance factor 
                                                                        !In AM3, the first layer is top of ATM, 
                                                                        !but in FASTJ, the first layer is surface           
      real*8,           intent(in)      ::  phalf1(L1_),&               !Pressure at boundaries (Pa) 49 layer
                                            zhalf1(L1_),&               !height at layer boundaries      
                                            pfull1(L_),&                !Pressure at mid-layer (Pa)  48
                                            tfull(L_)                   !Temperature at mid-Layer (K)           
      real*8,           intent(in)      ::  XO3(L_)                     !O3 mixing ratio (VMR)              
      real*8,           intent(in)      ::  pwt(L_)                     !Column air density (Kg/m2)              
      real*8,           intent(in)      ::  lwc(:,:)                    !(L_,ncloud=2)Liquid water content (Kg/Kg)
      real*8,           intent(in)      ::  cloudsf(:,:)                !Clouds fraction
      integer,          intent(in)      ::  clouds_ndx(:,:)             ! clouds type index =2: 1 water, 2 ice
      real*8,           intent(in)      ::  qfld(L_)                    !Specific humidity (kg/kg)
      real*8,           intent(in)      ::  srf_alb                     !Surface albedo               
      real*8,           intent(in)      ::  aerop(:,:)                  !Aerosol mass path (g/m2) (L_,AERONUM)
      integer,          intent(in)      ::  aeron(:,:)                  !(AERONUM) Aerosol type index, see FJX_scat.dat 
                                                                        !for aerosol type info 
      real*8,           intent(in)      ::  o3_column_top               !O3 column density, a small number
      real*8,           intent(out)     ::  ZPJ(JVL_,JVN_)              !ZPJ(JVL_,JVN_) 2-D array of J-values                                                  
!-----------------------------------------------------------------------
      real*8                            ::  SZA                         !solar zenith angle
      real*8                            ::  FREFL                       !fraction of energy refle                                            
      real*8                            ::  ZH(L1_)                     !height of each layer (cm)
      real*8                            ::  SCALEH                      !scale height (cm)
      real*8                            ::  pz, qs, wrk, es             !parameters to calc RH
      real*8                            ::  phalf(L1_), pfull(L_)       !pressures at boundaries and mid-layer (hPa)
      real*8, dimension(L1_+1)          ::  TTJ,&                       !Temperature used by FAST-J(K)
                                            DDJ,&                       !Air column density (molecs/cm2) for each layer
                                            ZZJ,&                       !O3 column density (molecs/cm2)
                                            ZHL                         !Distance of each layer to surface (cm), same to ZH,
                                                                        !but have another top layer, used to calculate cloud density
      real*8, dimension(L1_+1)          ::  PPJ                         !pressure at boundaries (hPa)
!--------key amtospheric data needed to solve plane-parallel J--------- 
      integer,dimension(L2_+1)          ::  JXTRA
      real*8, dimension(W_)             ::  FJTOP,&
                                            FJBOT,&
                                            FSBOT,&
                                            FLXD0,&
                                            RFL 
      real*8, dimension(L_, W_)         ::  AVGF,&
                                            FJFLX 
      real*8, dimension(L1_,W_)         ::  DTAUX,&
                                            FLXD 
      real*8, dimension(8,L1_,W_)       ::  POMEGAX 
      real*8, dimension(W_,L1_)         ::  FFX 
      real*8, dimension(W_,8)           ::  FFXNET 
!---flux/heating arrays (along with FJFLX,FLXD,FLXD0)                   
      real*8                            ::  FLXJ(L1_),&
                                            FFX0,&
                                            FXBOT,&
                                            FABOT                                                                         
      real*8                            ::  ODABS,&
                                            ODRAY 
      real*8                            ::  RFLECT,&
                                            FREFS,&
                                            FREFI !SOLF
      real*8                            ::  AMF2(2*L1_+1,2*L1_+1) 
!------------key SCATTERING arrays for clouds+aerosols------------------
      real*8                            ::  OPTX(5),&
                                            SSAX(5),&
                                            SLEGX(8,5) 
      real*8                            ::  OD(5,L1_),&
                                            SSA(5,L1_),&
                                            SLEG(8,5,L1_) 
      real*8                            ::  OD600(L1_) 
      real*8                            ::  PATH,PATH1,PATH2,&
                                            DENSWP,&
                                            DENS,&
                                            RH,&
                                            RHO,&
                                            REFF,&
                                            XTINCT,&
                                            ODCLD,ODCLD1,ODCLD2                         !optical depth of cloud
!------------key arrays AFTER solving for J's---------------------------
      real*8                            ::  FFF(W_,JVL_),&
                                            VALJ(X_) 
      real*8                            ::  FLXUP(W_),&
                                            FLXDN(W_),&
                                            DIRUP(W_),&
                                            DIRDN(W_) 
!------------2-D array of J_s returned by JRATET------------------------    
      real*8                            ::  VALJL(JVL_,NJVAL)                                                                         
      integer                           ::  I,J,K,&
                                            L,M,&
                                            KMIE,KW,&
                                            NAER,NDCLD,NDCLD1,NDCLD2, &
                                            RATIO(W_) 
      real*8                            ::  XQO3,XQO2,&
                                            DTAUC,&
                                            WAVE, &
!                                           FLINT, &
                                            TTT 
      logical                           ::  test_J = .false.

!-------------check model initialization -------------------------------
      if (.not. module_is_initialized) &
            call error_mesg ('ATMOS: fastjx_init:','fastjx_init must be called first.', FATAL)
!-----------------------------------------------------------------------
      ZPJ(:,:)  = 0.d0 
      FFF(:,:) = 0.d0 
      FREFI = 0.d0 
      FREFL = 0.d0 
      FREFS = 0.d0       

!-----------------------------------------------------------------------
!     call SOLARZ(UTIME,IDAY,YGRD,XGRD1, SZA,U0,SOLF) 
!-----------------------------------------------------------------------
!-----SOLF = 1.d0   ! this needs to be dropped to include 6.7% annual cy
 
      SZA    = acos(U0)*180./ 3.141592653589793d0
                                                                                                                                                
!-----check for dark conditions SZA > 98.0 deg => tan ht = 63 km          
!                        or         99.                  80 km          
      if (SZA .gt. SZAMAX) goto 99 
                                                                        
!-----load the amtospheric column data     
      phalf=phalf1/100.
      pfull=pfull1/100.  
      
!-----calculate effective altitude of each CTM level edge (cm)
!      ZH(1) = 16d5*log10(1000.d0/phalf(L1_))
!      do L = 1,L_
!         SCALEH  = 1.3806d-19*MASFAC*tfull(L1_-L)
!         ZH(L+1) = ZH(L) -(log(phalf(L1_-L)/phalf(L1_-L+1))*SCALEH)
!      enddo
      do L=1,L1_
         ZH(L)= zhalf1(L1_+1-L)*100.            ! height at boundaries (cm)
      end do
!            if (mpp_pe() == mpp_root_pe() ) then 
!              write(*,*) 'fastjx: zh(cm)=', zh    
!           end if                                                         
!---calculate spherical weighting functions (AMF: Air Mass Factor)      
      do L = 1,L1_ 
         ZHL(L) =   ZH(L) 
      enddo 
      ZHL(L1_+1) = ZHL(L1_) + ZZHT 
                                  
      do L = 1,L1_                !L=1 surface  AM3 L=1 TOA
!        PPJ(L) = ETAA(L) + ETAB(L)*P(ILNG,JLAT) 
!        TTJ(L) = TJ(ILNG,JLAT,L) 
!        DDJ(L) = DM(ILNG,JLAT,L) 
!        ZZJ(L) = DO3(ILNG,JLAT,L) 
         PPJ(L) = phalf(L1_+1-L)   ! pressure at boundaries (hPa)
         if(L <L1_)then
            TTJ(L) = tfull(L1_ -L) !TJ(ILNG,JLAT,L) 
            DDJ(L) =   pwt(L1_ -L)*6.023e22/28.92   !DM(ILNG,JLAT,L)  air den moleculars/cm2 
            ZZJ(L) =   XO3(L1_ -L)*DDJ(L)           !O3 column density (molecs/cm2) per layer
         else
            TTJ(L) = TTJ(L1_-1) + (TTJ(L1_-1)-TTJ(L1_-2))*(ZHL(L1_+1)-ZHL(L1_-1))/(ZHL(L1_)-ZHL(L1_-2))                         !K
            DDJ(L) = phalf(1)/pfactor1*100.     !hPa=>molec/cm2
            ZZJ(L) = o3_column_top                !molec/cm2 
         end if 
      enddo 
      PPJ(L1_+1) = 0.d0        
                                                            
!-----------------------------------------------------------------------
      call SPHERE2 (U0,RAD,ZHL,ZZHT,AMF2,L1_)                           !calculate how much sunlight is blocked by earth
!                                                                       !at a few lowest layers when zenith angle is high
!-----------------------------------------------------------------------   
    
                                                                     
!---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8)
!---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      do L = 1,L1_ 
         OD600(L) = 0.D0 
         do K=1,5 
            OD(K,L)  = 0.d0 
            SSA(K,L) = 0.d0 
            do I=1,8 
               SLEG(I,K,L) = 0.d0 
            enddo 
         enddo 
      enddo 
      
!-----sst and phase-fn are weighted by optical depth                                                                                   
      do L = 1,L1_   
   !      do M = 1,size(lwc,2)      
                                                                            
!-----cloud in layer:  if valid cloud index (now 4:13)        
            if(L<L1_)then            
               NDCLD1 = clouds_ndx(L1_-L,1)                                     !CLDNDX(ILNG,JLAT,L) 
               NDCLD2 = clouds_ndx(L1_-L,2)                                     !CLDNDX(ILNG,JLAT,L) 
               PATH1  = lwc(L1_-L,1)*pwt(L1_-L)*1.e3                            !(g/m2)  CLDLWP(ILNG,JLAT,L) liquid cloud 
               PATH2  = lwc(L1_-L,2)*pwt(L1_-L)*1.e3                            !(g/m2)  CLDLWP(ILNG,JLAT,L)  ice cloud
            else
               NDCLD1 = 0.                                              !CLDNDX(ILNG,JLAT,L) 
               NDCLD2 = 0.                                              !CLDNDX(ILNG,JLAT,L) 
               PATH1  = 0.                                              !(g/m2)  CLDLWP(ILNG,JLAT,L)    
               PATH2  = 0.                                              !(g/m2)  CLDLWP(ILNG,JLAT,L)    
            endif
!               if (mpp_pe() == mpp_root_pe() ) then 
!                 write(*,*)'fastjx: cloud PATH(g/m2): L(m)=', ZHL(L)/100.,'    PATH1=', PATH1, '       PATH2=',PATH2
!               endif
   
  
!            DENS  = PATH / ZH(L)
!            if (mpp_pe() == mpp_root_pe() ) then 
!              write(*,*) 'fastjx: dens(org)=', dens    
!           end if
!            if (mpp_pe() == mpp_root_pe() ) then 
!              write(*,*) 'fastjx: dens(new)=', dens    
!           end if       
            if (PATH1+PATH2 .gt. 0.d0) then 
               if (NDCLD1 .lt.4 .or. NDCLD1.gt.13) then 
                   NDCLD1 = 9 
               endif 
!---now calculate cloud OD at 600 nm from mixed formulae:               
               REFF = RAA(NDCLD1) 
               RHO  = DAA(NDCLD1) 
      
      
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(
               if(L<L1_)then
                  ODCLD1 = PATH1*0.75d0*QAA(4,NDCLD1)/(REFF*RHO) * cloudsf(L1_-L,1)**0.5                ! assuming approximate random overlap
               else
                  ODCLD1 = 0.
               end if  
       
       
               REFF = RAA(NDCLD2) 
               RHO  = DAA(NDCLD2)
!               if (NDCLD2 .eq. 12 .or. NDCLD2 .eq. 13) then 
               DENS  = PATH2 /(ZHL(L+1)-ZHL(L)) * 100.                  ! !!!JUL CHANGED DENSITY HERE(g/m2)/cm       
!---Ice Water Content (density at g/m3) varies from 0.0001 to 0.1       
!---correct R-eff for ice clouds, increase Reff (microns) with density  
               REFF = 50.d0 * (1.d0 + 8.333d0*DENS)   ! shoudl check the unit of DENS jul !!!!!!!!!!!!!!
!               endif 
               if(L<L1_)then
                  ODCLD2 = PATH2*0.75d0*QAA(4,NDCLD2)/(REFF*RHO) * cloudsf(L1_-L,2)**0.5                ! assuming approximate random overlap
               else
                  ODCLD2 = 0.
               end if  
               ODCLD = ODCLD1 + ODCLD2   
               if(ODCLD1 .ge. ODCLD2)then
                  NDCLD = NDCLD1
               else
                  NDCLD = NDCLD2
               endif
!               if (mpp_pe() == mpp_root_pe() ) then
!                 write(*,*)'fastjx: L(m)=', ZHL(L)/100., '     NDDLD=',NDCLD,  'ODCLD1=',ODCLD1, '     ODCLD2=',ODCLD2, '      DENS2=',DENS
!              endif
               call OPTICL (OPTX,SSAX,SLEGX,  ODCLD,NDCLD) 
               do K=1,5 
                  OD(K,L)  = OD(K,L)  + OPTX(K) 
                  SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K) 
                  do I=1,8 
                     SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K) 
                  enddo 
               enddo 
            endif !if path >0
  !       end do !M
                                                                       
!---use OD of clouds (not aerosols) at 600 nm to determine added layers 
         OD600(L) = OD(4,L)                                                                         
!---aerosols in layer: check aerosol index                              
!---this uses data from climatology OR from current CTM (STT of aerosols                                                                       
!---FIND useful way to sum over different aerosol types!                
         do M = 1,size(aerop,2)
            if(L<L1_)then
               NAER = aeron(L1_-L,M) 
               PATH = aerop(L1_-L,M) !g/m2
            else
               NAER = 0 
               PATH = 0.D0
            endif                                                                      
!---subroutines OPTICA & OPTICM return the same information:            
!---  optical depth (OPTX), single-scat albedo (SSAX) and phase fn (SLEG
!---subs have slightly different inputs:                                
!---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index    
!---  UMich aerosols use relative humidity (RH)                         
            if (PATH .gt. 0.d0 .and. NAER .gt. 0) then
                call OPTICAM3 (OPTX,SSAX,SLEGX,  PATH, NAER)
!               if (NAER .gt.0) then 
!                   call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER) 
!               else 
!                   NAER = -NAER
!                   call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,NAER) 
!               endif
                do K=1,5 
                   OD(K,L)  = OD(K,L)  + OPTX(K) 
                   SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K) 
                   do I=1,8 
                      SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K) 
                   enddo !I
                enddo !K
            endif
         enddo !M
                                                                        
         do K=1,5 
            if (OD(K,L) .gt. 0.d0) then 
               SSA(K,L) = SSA(K,L)/OD(K,L) 
               do I=1,8 
                  SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L) 
               enddo 
            endif 
         enddo                                                                        
      enddo !L
                                                                        
!---can add aerosol OD at 600 nm to determine added layers, but not done yet
!---    OD600(L) = OD(4,L) + aerosol OD's                               
                                                                                                                                                
!---when combining with Rayleigh and O2-O3 abs, remember the SSA and    
!---  phase fn SLEG are weighted by OD and OD*SSA, respectively.        
                                                                        
!---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
!       additonal levels at top of clouds (now uses log spacing)        
!-----------------------------------------------------------------------
      call EXTRAL(OD600,L1_,L2_,N_,JTAUMX,ATAU,ATAU0, JXTRA) 
!-----------------------------------------------------------------------
                                                                        
!---set surface reflectance                                             
      RFLECT = srf_alb !SA(ILNG,JLAT) 
      RFL(:) = max(0.0,min(1.0,RFLECT))                                                                                                                                                 
!-----------------------------------------------------------------------
!---Loop over all wavelength bins to calc mean actinic flux AVGF(L)     
!-----------------------------------------------------------------------                                                                       
      do K = 1,W_ 
                                                                        
         WAVE = WL(K) 
!---Pick nearest Mie wavelength to get scattering properites------------                                           
                                KMIE=1 ! use 200 nm prop for <255 nm                                        
         if( WAVE .gt. 255.d0 ) KMIE=2 ! use 300 nm prop for 255-355 nm                                        
         if( WAVE .gt. 355.d0 ) KMIE=3 ! use 400 nm prop for 355-500 nm
         if( WAVE .gt. 500.d0 ) KMIE=4 
         if( WAVE .gt. 800.d0 ) KMIE=5                                                                                                                                                 
!---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
!---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
         do L = 1,L1_ 
            TTT     = TTJ(L) 
            XQO3 = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2)                    &   ! O3 cross-section intepolation
                           ,QO3(K,1),QO3(K,2),QO3(K,3))                
                                                                        
            XQO2 = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1)                    &   ! O2 cross-section intepolation
                           ,QO2(K,1),QO2(K,2),QO2(K,3))                
                                                                        
            ODABS = XQO3*ZZJ(L) + XQO2*DDJ(L)*0.20948d0         !o3+o2
            ODRAY = DDJ(L)*QRAYL(K)                             !air
                                                                        
            DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY             !OD: cloud+aerosol
                                                                        
            do I=1,8 
               POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L) 
            enddo 
            POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY 
            POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY 
            do I=1,8 
               POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K) 
            enddo 
         enddo !L                                                                       
      enddo !K                                                                       
!-----------------------------------------------------------------------
                                                                        
      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,                      &
             AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)                  
                                                                        
!-----------------------------------------------------------------------                                                                       
      do K = 1,W_ 
                                                                        
!----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convection)
!----     also at bottom (DN), does not include diffuse reflected flux. 
         FLXUP(K) =  FJTOP(K) 
         DIRUP(K) = -FLXD0(K) 
         FLXDN(K) = -FJBOT(K) 
         DIRDN(K) = -FSBOT(K) 
                                                                        
         do L = 1,JVL_ 
            FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K) 
         enddo 
         FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K) 
         FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K) 
         FREFS = FREFS + SOLF*FL(K)/WL(K) 
                                                                        
!---for each wavelength calculate the flux budget/heating rates:        
!  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(
!            but for spherical atmosphere!                              
!  FJFLX(L) = diffuse flux across top of layer L                        
                                                                        
!---calculate divergence of diffuse flux in each CTM layer (& t-o-a)    
!---     need special fix at top and bottom:                            
!---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.     
         FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K)) 
         FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K)) 
         FLXJ(1) = FJFLX(1,K) - FXBOT 
         do L=2,L_ 
            FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K) 
         enddo 
         FLXJ(L_+1) = FJTOP(K) - FJFLX(L_,K) 
!---calculate net flux deposited in each CTM layer (direct & diffuse):  
         FFX0 = 0.d0 
         do L=1,L1_ 
            FFX(K,L) = FLXD(L,K) - FLXJ(L) 
            FFX0 = FFX0 + FFX(K,L) 
         enddo 
                                                                        
!  NB: the radiation level ABOVE the top CTM level is included in these 
!      these are the flux budget/heating terms for the column:          
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spheric
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)  
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface       
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos         
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos            
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)  
!       these are surface fluxes to compare direct vs. diffuse:         
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for sr
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)        
                                                                        
         FFXNET(K,1) = FLXD0(K) 
         FFXNET(K,2) = FSBOT(K) 
         FFXNET(K,3) = FLXD0(K) + FSBOT(K) 
         FFXNET(K,4) = FJTOP(K) 
         FFXNET(K,5) = FFX0 
         FFXNET(K,6) = FABOT 
         FFXNET(K,7) = FSBOT(K) 
         FFXNET(K,8) = FJBOT(K)                                                                         
!-----------------------------------------------------------------------                                         
      enddo ! end loop over wavelength K 
!-----------------------------------------------------------------------                                                                                                           
      FREFL = FREFL/FREFS !calculate reflected flux (energy wei
      FREFI = FREFI/FREFS 
                                                                        
!---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin
                                                                         
!-----------------------------------------------------------------------
      call JRATET(pfull,tfull,FFF, VALJL) 
!-----------------------------------------------------------------------                                                                         
!---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)    

      test_j =.false.
      do L = 1,JVL_ 
         do J = 1,NRATJ 
            if (JIND(J).gt.0) then 
               ZPJ(JVL_+1-L,J) = VALJL(L,JIND(J))*JFACTA(J) 
               if(ZPJ(JVL_+1-L,J) .lt. 0.0D0) then
                 test_j = .true.
               endif
            else 
               ZPJ(JVL_+1-L,J) = 0.d0 
            endif 
         enddo 
      enddo 
      
!      if(test_j)then
!            if (mpp_pe() == mpp_root_pe() ) then
!              write(*,*) 'fastjx: test_negative: (1) input variables'
!              write(*,*) 'fastjx: solar zenith angly',SZA, '  COSZ', U0   
!              write(*,*) 'fastjx:(1): ZHL', ZHL
!              write(*,*) 'fastjx:(1): PPJ', PPJ  
!              write(*,*) 'fastjx:(1): TTJ', TTJ  
!              write(*,*) 'fastjx:(1): DDJ', DDJ  
!              write(*,*) 'fastjx:(1): ZZJ', ZZJ  
!              write(*,*) 'fastjx:(1): AMF2', AMF2  
!              write(*,*) 'fastjx:(1): FFF', FFF  
!              write(*,*) 'fastjx:(1): AVGF', AVGF  
          
!            end if       
!      endif
      
                                                                        
!---diagnostics that are NOT returned to the CTM code                   
                                                                        
!-----------------------------------------------------------------------
!      write(6,*)'fast-JX-(6.4)----PHOTOJ internal print: Atmosphere----' 
!---used last called values of DTAUX and POMEGAX, should be 600 nm      
                                                                        
!      call JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,                             &
!     &            DTAUX(1,W_),POMEGAX(1,1,W_),JXTRA)                    
                                                                        
!---PRINT SUMMARY of mean intensity, flux, heating rates:               
!      write(6,*) 
!      write(6,*)'fast-JX(6.4)----PHOTOJ internal print: Mean Intens----' 
!      write(6,'(a,5f10.4)')                                             &
!     & ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',                &
!     &  RFLECT,SZA,U0,FREFI,FREFL                                       
                                                                        
!      write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1) 
!      write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1) 
!      write(6,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin' 
!      do L = JVL_,1,-1 
!       do K=NW1,NW2 
!        RATIO(K) = (1.d5*FFF(K,L)/FL(K)) 
!       enddo 
!        write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1) 
!      enddo 
                                                                        
!      write(6,*) 
!      write(6,*)'fast-JX(6.4)----PHOTOJ internal print: Net Fluxes----' 
!      write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
!      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
!      write(6,*) ' ---NET FLUXES--- ' 
!      write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1) 
!      write(6,*) ' ---SRF FLUXES--- ' 
!      write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1) 
!      write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1) 
!      write(6,'(2a)') '  ---NET ABS per layer:       10000=Fsolar',     &
!     & '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
!      do L = JVL_,1,-1 
!       do K=NW1,NW2 
!        RATIO(K) = 1.d5*FFX(K,L) 
!       enddo 
!        write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1) 
!      enddo                                                                         
!-----------------------------------------------------------------------
                                                                        
   99 continue 
      return 
      end subroutine fastjx_photo                                          
                                                                        
!<<<<<<<<<<<<<<<<<<<<<<<end CTM-fastJX linking subroutines<<<<<<<<<<<<<<
 

      subroutine fastjx_end
         implicit none
         module_is_initialized = .false.
      end subroutine fastjx_end

                                                                        
!-----------------------------------------------------------------------
!      subroutine SET_ATM(GMTAU) 
!-----------------------------------------------------------------------
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'cmn_metdat.f' 
                                                                        
!      real*8 MASFAC,SCALEH,GMTAU, PJ(L1_+1) 
!      integer I,J,L,N 
!                                                                       
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0) 
!                                                                        
!      do J = 1,J_ 
!      do I = 1,I_ 
!        P(I,J) = PMEAN(I,J) 
!        do L = 1,L1_ 
!          PJ(L) = ETAA(L) + ETAB(L)*P(I,J) 
!        enddo 
!          PJ(L1_+1) = 0.d0 
!        do L = 1,L1_ 
!          DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC 
!        enddo 
!        do L = 1,L_ 
!          TJ(I,J,L)  = T(I,J,L) 
!        enddo 
!-------calculate effective altitude of each CTM level edge             
!          ZH(I,J,1) = 16d5*log10(1000.d0/P(I,J)) 
!        do L = 1,L_ 
!          SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L) 
!          ZH(I,J,L+1) = ZH(I,J,L) -(log(PJ(L+1)/PJ(L))*SCALEH) 
!        enddo 
!      enddo 
!      enddo 
                                                                        
!---load O3 from CTM is being calculated:                               
!      do N = 1,NTR_ 
!        if (TCNAME(N) .eq. 'O3') then 
!          do J = 1,J_ 
!          do I = 1,I_ 
!                              ! unit of DO3:  # molecules/cm^2          
!          do L = 1,L_ 
!             DO3(I,J,L) = 6.023d26*STT(I,J,L,N)/48.d0                   &
!     &                    *1.d-4 /AREAXY(I,J)                           
!            enddo 
!          enddo 
!          enddo 
!        endif 
!      enddo 
                                                                        !
!      return 
!      END                                           
                                                                        
!      subroutine SET_AER(GMTAU) 
!      return 
!      END                                           
                                                                        
!      subroutine SET_CLD(GMTAU) 
!      return 
!      END                                           
                                                                        
      subroutine RD_JS(NAMFIL) 
!-----------------------------------------------------------------------
!  Reread the chem_Js.dat file to map photolysis rate to reaction       
!  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'             
!-----------------------------------------------------------------------
!                                                                       
!     jfacta    Quantum yield (or multiplication factor) for photolysis 
!     jlabel    Reference label identifying appropriate J-value to use  
!     ipr       Photolysis reaction counter - should total 'JVN_'       
!                                                                       
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_metdat.f' 
!      include 'cmn_JVdat.f' 
!                                                                       
      character(*), intent(in) ::  NAMFIL 
                                                                        
      integer  NJ1 
      integer  IPR, I, J, K 
      character*120 CLINE 
      character*7  T_JX, T_CHEM 
                                                                        
! Reread the chem_Js.dat file to map photolysis rate to reaction        
!                     Read in quantum yield jfacta and fastj2 label jlab
      IPR = 0 
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
   10 read (NJ1,'(A)',err=20)  CLINE 
      if (IPR .eq. JVN_) goto 20 
                                                                        
      if (CLINE(2:5).eq.'9999') then 
        go to 20 
      elseif (CLINE(1:1).eq.'#') then 
        go to 10 
      elseif (CLINE(5:5).eq.'$') then 
        go to 10 
      else 
        IPR = IPR+1 
        read (CLINE(79:83),'(F5.1)') JFACTA(IPR) 
        read (CLINE(86:92),'(A7)')   JLABEL(IPR) 
        JFACTA(IPR) = JFACTA(IPR)/100.d0 
        go to 10 
      endif 
   20 call mpp_close(NJ1) 
                                                                        
      NRATJ = IPR 
                                                                        
!-----------------------------------------------------------------------
!  compare Xsections titles with J-values listed in chem code (jratd.dat
!  map J-values needed for chemistry (chem_Js.dat) onto the fast-JX rate
!  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<
!          >>>this must now follow the read in of Xsects, etc<<<        
!-----------------------------------------------------------------------
                                                                        
!---Zero / Set index arrays that map Jvalue(j) onto rates               
      do J = 1,JVN_ 
        JIND(J) = 0 
      enddo 
      do J = 1,NJVAL 
        T_JX = TITLEJ(J) 
      do K = 1,NRATJ 
        T_CHEM = JLABEL(K) 
        if (T_CHEM(1:6) .eq. T_JX(1:6)) JIND(K)=J 
      enddo 
      enddo 
                                                                        
!      write(6,'(a,i4,a)') ' Photochemistry Scheme with ',IPR,' J-values' 
      do K=1,NRATJ 
        J = JIND(K) 
        if (J.eq.0) then 
!         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),         &
!     &         ' has no mapping onto onto fast-JX'                      
        else 
!         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),         &
!     &         ' mapped onto fast-JX:',J,TITLEJ(J)                      
        endif 
      enddo 
                                                                        
      return 
      END  subroutine RD_JS                                         
                                                                                                                                                
!-----------------------------------------------------------------------
!      subroutine RD_PROF(NJ2,NAMFIL) 
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles                         
!-----------------------------------------------------------------------
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'cmn_metdat.f' 
                                                                        
!      integer, intent(in) ::  NJ2 
!      character(*), intent(in) ::  NAMFIL 
!                                                                       
!      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216 
!      real*8  OFAC, OFAK 
!      character*72 TITLE0 
!                                                                       
!RSH if this code is ever reactivated, change to use mpp_open:
!      open (NJ2,file=NAMFIL,status='old',form='formatted') 
!!      read (NJ2,'(A)') TITLE0 
!      read (NJ2,'(2I5)') NTLATS,NTMONS 
!      write(6,'(1X,A)') TITLE0                                         
!      write(6,1000) NTLATS,NTMONS 
!      N216  = min(216, NTLATS*NTMONS) 
!      do IA = 1,N216 
!        read (NJ2,'(1X,I3,3X,I2)') LAT, MON 
!        M = min(12, max(1, MON)) 
!        L = min(18, max(1, (LAT+95)/10)) 
!        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41) 
!        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31) 
!      enddo 
!      close (NJ2) 
!RSH   call mpp_close (NJ2)
                                                                        
!  Extend climatology to 100 km                                         
!      OFAC = exp(-2.d5/5.d5) 
!      do I = 32,51 
!        OFAK = OFAC**(I-31) 
!        do M = 1,NTMONS 
!        do L = 1,NTLATS 
!          OREF(I,L,M) = OREF(31,L,M)*OFAK 
!        enddo 
!        enddo 
!      enddo 
!      do L = 1,NTLATS 
!      do M = 1,NTMONS 
!      do I = 42,51 
!        TREF(I,L,M) = TREF(41,L,M) 
!      enddo 
!      enddo 
!      enddo 
                                                                        
!      return 
! 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon') 
!      END                                           
                                                                        
!-----------------------------------------------------------------------
!      subroutine SET_CLD0(TINIT,CLDP,NCLD,ALBEDO) 
!-----------------------------------------------------------------------
!  Routine to set cloud and surface properties: now loads the input the 
!     input profiles (FJX_test.dat) INCLUDING Temperatures and Surf Albe
!                                                                       
!  >>>>>> this subroutine will need to be customized                    
!  >>>>>> it is separate from aerosols since it comes from met fields   
!-----------------------------------------------------------------------
!      implicit none 
!!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_metdat.f' 
!      include 'cmn_JVdat.f' 
                                                                        
!      integer, intent(in) :: NCLD(L_) 
!      real*8,  intent(in) :: TINIT(L_),CLDP(L_),ALBEDO 
!      integer  I, J, L 
!                                                                       
!      do J = 1,J_ 
!      do I = 1,I_ 
!        do L = 1,L_ 
!          CLDLWP(I,J,L) = CLDP(L) 
!          CLDNDX(I,J,L) = NCLD(L) 
!          T(I,J,L)      = TINIT(L) 
!        enddo 
!          SA(I,J) = ALBEDO 
!      enddo 
!      enddo 
!                                                                       
!      return 
!      END                                           
                                                                        
                                                                        
!-----------------------------------------------------------------------
!      subroutine SET_AER0 (AER1,AER2, NAA1,NAA2) 
!-----------------------------------------------------------------------
!  Set up aerosols >>>>customize for climatology or CTM as source.      
!     this needs to load in the aerosol mass/column and type in each lay
!  In a CTM, this should be done as a climatology or in real time using 
!     computed STT(kg) of the aerosol and calculating PATH=1000.*STT/ARE
!-----------------------------------------------------------------------
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_metdat.f' 
!      include 'cmn_JVdat.f' 
                                                                        
!      integer, intent(in) :: NAA1(L_),NAA2(L_) 
!      real*8,  intent(in) :: AER1(L_),AER2(L_) 
                                                                        
!      integer I, J, L, K 
                                                                        
!          AER1P(:,:,:) = 0.d0 
!!          AER2P(:,:,:) = 0.d0 
!          AER1N(:,:,:) = 0 
!          AER2N(:,:,:) = 0 
!      do J = 1,J_ 
!      do I = 1,I_ 
!        do L = 1,L_ 
!          AER1P(I,J,L) = AER1(L) 
!          AER2P(I,J,L) = AER2(L) 
!          AER1N(I,J,L) = NAA1(L) 
!          AER2N(I,J,L) = NAA2(L) 
!        enddo 
!      enddo 
!      enddo 
                                                                        
!      END                                           
                                                                        
                                                                        
!-----------------------------------------------------------------------
!      subroutine SET_ATM0 
!-----------------------------------------------------------------------
!  Routine to set up atmospheric profiles required by Fast-J2 using a   
!  doubled version of the level scheme used in the CTM. First pressure  
!  and z* altitude are defined, then O3 and T are taken from the supplie
!  climatology and integrated to the CTM levels (may be overwritten with
!  values directly from the CTM, if desired).                           
!                                       Oliver (04/07/99) & MJP (7/05)  
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!     PJ       Pressure at boundaries of model levels (hPa)             
!     MASFAC   Conversion factor for pressure to column density         
!                                                                       
!     TJ       Temperature profile on model grid                        
!     DM       Air column for each model level (molecules.cm-2)         
!     DO3      Ozone column for each model level (molecules.cm-2)       
!     ZH       Altitude of boundaries of model levels (cm)              
!     PSTD     Approximate pressures of levels for supplied climatology 
!-----------------------------------------------------------------------
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'cmn_metdat.f' 
                                                                        
!      integer  I, J, K, L, M, N 
!      real*8   PSTD(52),OREF2(51),TREF2(51),PJ(L1_+1) 
!      real*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH 
!                                                                       
!  Select appropriate month                                             
!      M = max(1,min(12,MONTH)) 
                                                                        
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0) 
                                                                        
!-----------------------------------------------------------------------
!      do J = 1,J_ 
                                                                        
!  Select appropriate latitudinal profiles                              
!          N = max(1, min(18, (int(YDGRD(J))+99)/10 )) 
!  Temporary zonal arrays for climatology data                          
!          do K = 1,51 
!            OREF2(K) = OREF(K,N,M) 
!            TREF2(K) = TREF(K,N,M) 
!          enddo 
!                                                                        
!        do I = 1,I_ 
                                                                        
!  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
!  with mass (pressure) weighting, assuming constant mixing ratio and   
!  temperature half a layer on either side of the point supplied.       
!  L1_ = L_+1:                                                          
!       PJ(L=1:L1_) = pressure at CTM layer edge, PJ(L1_+1)=0 (top-of-at
                                                                        
!           PJ(L1_+1) = 0.d0 
!         do K = 1,L1_ 
!           PJ(K) = ETAA(K) + ETAB(K)*PMEAN(I,J) 
!         enddo 
                                                                        
!  Set up pressure levels for O3/T climatology - assume that value      
!  given for each 2 km z* level applies from 1 km below to 1 km above,  
!  so select pressures at these boundaries. Surface level values at     
!  1000 mb are assumed to extend down to the actual P(ILNG,JLAT).       
!                                                                       
!           PSTD(1) = max(PJ(1),1000.d0) 
!           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0) 
!           DLOGP   = 10.d0**(-2.d0/16.d0) 
!         do K = 3,51 
!           PSTD(K) = PSTD(K-1)*DLOGP 
!         enddo 
!           PSTD(52)  = 0.d0 
                                                                        
!         do L = 1,L1_ 
!           F0 = 0.d0 
!           T0 = 0.d0 
!           do K = 1,51 
!             PC   = min(PJ(L),PSTD(K)) 
!             PB   = max(PJ(L+1),PSTD(K+1)) 
!             if (PC .gt. PB) then 
!               XC = (PC-PB)/(PJ(L)-PJ(L+1)) 
!               F0 = F0 + OREF2(K)*XC 
!               T0 = T0 + TREF2(K)*XC 
!             endif !
!           enddo 
!           TJ(I,J,L)  = T0 
!           DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC 
!           DO3(I,J,L) = F0*1.d-6*DM(I,J,L) 
!         enddo 
                                                                        
!  Calculate effective altitudes using scale height at each level       
!           ZH(I,J,1) = 0.d0 
!         do L = 1,L_ 
!           SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L) 
!           ZH(I,J,L+1) = ZH(I,J,L) -( LOG(PJ(L+1)/PJ(L)) * SCALEH ) 
!         enddo 
                                                                        
!        enddo 
!      enddo 
                                                                        
!      return 
!      END                                           
                                                                        
!<<<<<<<<<<<<<<<<<<<<<<<<end CTM-specific subroutines<<<<<<<<<<<<<<<<<<<
                                                                        
                                                                        
!<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
                                                                        
!-----------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD) 
!-----------------------------------------------------------------------
!---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-60
!---   this now requires the cloud layer OD rather than water path      
!---   assumed CLOUD Optical Depth given with type for optical propertie
!                                                                       
!  04 W_H01   (H1/Deir)GAMMA:r-m=0.1/alf=2 n=1.335   reff=0.250___G=.094
!  05 W_H04   (H1/Deir)GAMMA:r-m=0.4/alf=2 n=1.335   reff=1.000___G=1.50
!  06 W_H40   (H1/Deir)GAMMA:r-m=4.0/alf=2 n=1.335   reff=10.00___G=146.
!  07 W_C02   (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.5
!  08 W_C04   (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.1
!  09 W_C08   (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.
!  10 W_C13   (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.
!  11 W_L06   (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.
!  12 Ice-Hexagonal (Mishchencko)                    reff=50.00___G=999.
!  13 Ice-Irregular (Mishchencko)                    reff=50.00___G=999.
                                                                        
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
                                          ! optical depth of layer      
      real*8, intent(out)::    OPTD(5) 
                                          ! single-scattering albedo    
      real*8, intent(out)::    SSALB(5) 
                                          ! scatt phase fn (Leg coeffs) 
      real*8, intent(out)::    SLEG(8,5) 
                                          ! optical depth of cloud layer
      real*8, intent(in)::     ODCLD 
                                          ! index of cloud layer:  4:13 
      integer,intent(inout)::  NDCLD 
                                                                        
      integer I,J 
      real*8  XTINCT, REFF,RHO 
                                                                        
!---default cloud type C1, Reff = 12 microns                            
      if (NDCLD .gt. 13 .or. NDCLD .lt.4) then 
         NDCLD = 9 
      endif 
                                                                        
!--rescale OD by Qext at 600 nm (J=4)                                   
      do J=1,5 
         OPTD(J) = ODCLD * QAA(J,NDCLD)/QAA(4,NDCLD) 
         SSALB(J) = SAA(J,NDCLD) 
        do I=1,8 
         SLEG(I,J) =  PAA(I,J,NDCLD) 
        enddo 
      enddo 
                                                                        
      return 
      END  subroutine OPTICL                                         
                                                                        
 
                                                                        
!-----------------------------------------------------------------------
      subroutine OPTICAM3 (OPTD,SSALB,SLEG, PATH,L) 
!-----------------------------------------------------------------------
!---for the AM3 aerosol data sets, calculates optical properties at fast-JX
!              std 5 wavelengths:200-300-400-600-999nm                  :                                    
!                                                                         
!                    
!                                                                        
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
                                          ! optical depth of layer      
      real*8, intent(out)::    OPTD(5) 
                                          ! single-scattering albedo    
      real*8, intent(out)::    SSALB(5) 
                                          ! scatt phase fn (Leg coeffs) 
      real*8, intent(out)::    SLEG(8,5) 
                                          ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     PATH 
                                          ! index of cloud/aerosols     
      integer,intent(inout)::     L 
                                                                        
      integer I,J 
!      real*8  XTINCT, REFF,RHO 
                                                                        
      if (L .gt. NAA_AM3 ) then 
         write(6,*) 'ATMOS:fastjx_photo:OPTICAM3:FATAL: aerosol index out-of-range: L/NAA',L,NAA_AM3
         call error_mesg ('ATMOS: fastjx_init:OPTICAM3','aerosol index out-of-range.', FATAL)
      endif 
                                                                        
!          MEE1 = MEE(L)                                                          
!         REFF = RAA(L) 
!         RHO = DAA(L) 
      do J=1,5 
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(
!         XTINCT = 0.75d0*QAA(J,L)/(REFF*RHO) 
         OPTD(J) = PATH*MEE_AM3(J,L) 
         SSALB(J) = SAA_AM3(J,L) 
       do I=1,8 
         SLEG(I,J) =  PAA_AM3(I,J,L) 
       enddo 
      enddo 
                                                                        
      return 
      END  subroutine OPTICAM3                                         
 
                                                                        
!-----------------------------------------------------------------------
!      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L) 
!-----------------------------------------------------------------------
!---for the UCI aerosol data sets, calculates optical properties at fast
!              std 5 wavelengths:200-300-400-600-999nm                  
!---UCI aersols optical data  v-6.1:                                    
!                                                                       
!   14 S-Bkg LOGN:r=.090 s=.600 reff=.221  n=1.514/1.473/1.459/1.448/1.4
!   15 S-Vol LOGN:r=.080 s=.800 reff=.386  n=1.514/1.473/1.459/1.448/1.4
!   16 UT-sulfate LOGN:r=0.05 s=.693 n=1.44           reff=0.166___rho=1
!   17 UT-sulfate LOGN:r=0.05 s=.693 n=1.46           reff=0.166___rho=1
!   18 UM-SULFate LOGN:r=.050 s=.642 n=1.53           reff=0.140___rho=1
!   19 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i     reff=0.140___rho=1
!   20 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i     reff=0.150___rho=1
!   21 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i  reff=0.149___rho=1
!   22 UM-FF04(%BC) LOGN:r=.050 s=.642 n=1.541+0.02i  reff=0.140___rho=1
!   23 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i  reff=0.140___rho=1
!   24 Mdust .15 (R.V. Martin generated phase fns)                      
!   25 Mdust .25 (R.V. Martin generated phase fns)                      
!   26 Mdust 0.4 (R.V. Martin generated phase fns)                      
!   27 Mdust 0.8 (R.V. Martin generated phase fns)                      
!   28 Mdust 1.5 (R.V. Martin generated phase fns)                      
!   29 Mdust 2.5 (R.V. Martin generated phase fns)                      
!   30 Mdust 4.0 (R.V. Martin generated phase fns)                      
                                                                        
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
                                          ! optical depth of layer      
!      real*8, intent(out)::    OPTD(5) 
                                          ! single-scattering albedo    
!      real*8, intent(out)::    SSALB(5) 
                                          ! scatt phase fn (Leg coeffs) 
!      real*8, intent(out)::    SLEG(8,5) 
                                          ! path (g/m2) of aerosol/cloud
!      real*8, intent(in)::     PATH 
                                          ! relative humidity (0.00->1.0
!      real*8, intent(in)::     RELH 
                                          ! index of cloud/aerosols     
!      integer,intent(inout)::     L 
                                                                        
!      integer I,J 
!      real*8  XTINCT, REFF,RHO 
                                                                        
!      if (L .gt. NAA .or. L .lt. 4) then 
!         write(6,*) 'ATMOS:fastjx_photo: aerosol index out-of-range: L/NAA',L,NAA 
!         L = 18 
!      endif 
                                                                        
!      if (L .lt. 14) then 
!         write(6,*) 'ATMOS:fastjx_photo: aerosol as cloud: L/NAA',L,NAA 
!      endif 
                                                                        
!         REFF = RAA(L) 
!         RHO = DAA(L) 
!      do J=1,5 
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(
!         XTINCT = 0.75d0*QAA(J,L)/(REFF*RHO) 
!         OPTD(J) = PATH*XTINCT 
!         SSALB(J) = SAA(J,L) 
!       do I=1,8 
!         SLEG(I,J) =  PAA(I,J,L) 
!       enddo 
!      enddo 
                                                                        
!      return 
!      END  subroutine OPTICA                                         
                                                                        
                                                                        
!-----------------------------------------------------------------------
!      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L) 
!-----------------------------------------------------------------------
!---for the U Michigan aerosol data sets, this generate fast-JX data for
!---NB Approximates the Legendre expansion(L) of the scattering phase fn
!---           as (2*L+1)*g**L                                          
!---UMAER(I,J,K,L):                                                     
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]                                  
!   J=1:6 = [200, 300, 400, 550, 600 , 1000 nm]                         
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]                     
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
!      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
                                          ! optical depth of layer      
!      real*8, intent(out)::    OPTD(5) 
                                          ! single-scattering albedo    
!      real*8, intent(out)::    SSALB(5) 
                                          ! scatt phase fn (Leg coeffs) 
!      real*8, intent(out)::    SLEG(8,5) 
                                          ! path (g/m2) of aerosol/cloud
!      real*8, intent(in)::     PATH 
                                          ! relative humidity (0.00->1.0
!      real*8, intent(in)::     RELH 
                                          ! index of cloud/aerosols     
!      integer,intent(inout)::     L 
                                                                        
!      integer KR,J 
!      real*8  R,FRH, GCOS, XTINCT 
                                                                        
!---calculate fast-JX properties at the std 5 wavelengths:200-300-400-60
!---interpolate in Relative Humidity                                    
!---extrapolate pahse fn from first term (g)                            
!      if (L .gt. 33) then 
!         write(6,*) 'ATMOS:fastjx_init: UM aer index too large: L',L 
!         L = 1 
!      endif 
                                                                        
!      R = 100.d0*min(1.d0, max(.01d0, RELH)) 
!         KR = (R/5.d0) 
!         KR = max(0, min(19, KR)) + 1 
!        if (KR.lt.20) then 
!         FRH = 0.20d0*(R - 5.d0*float(KR-1)) 
!        else 
!         FRH = 0.25d0*(R - 5.d0*float(KR-1)) 
!        endif 
                                                                        
!      do J=1,5 
!       SSALB(J) = UMAER(1,J,KR,L)*(1.d0-FRH) + UMAER(1,J,KR+1,L)*FRH 
                                                                        
!       XTINCT = UMAER(3,J,KR,L)*(1.d0-FRH) + UMAER(3,J,KR+1,L)*FRH 
!       OPTD(J) = PATH*XTINCT 
                                                                        
!       GCOS   = UMAER(2,J,KR,L)*(1.d0-FRH) + UMAER(2,J,KR+1,L)*FRH 
!       SLEG(1,J) =  1.d0 
!       SLEG(2,J) =  3.d0*GCOS 
!       SLEG(3,J) =  5.d0*GCOS**2 
!       SLEG(4,J) =  7.d0*GCOS**3 
!       SLEG(5,J) =  9.d0*GCOS**4 
!       SLEG(6,J) = 11.d0*GCOS**5 
!       SLEG(7,J) = 13.d0*GCOS**6 
!       SLEG(8,J) = 15.d0*GCOS**7 
!      enddo 
!      return 
!      END subroutine OPTICM                                          
                                                                        
!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL) 
!-----------------------------------------------------------------------
! in:                                                                   
!        PPJ(L1_+1) = pressure profile at edges    (hPa)                     
!        TTJ(L1_) = = temperatures at mid-level                         
!        FFF(K=1:NW, L=1:JVL_) = mean actinic flux                      
! out:                                                                  
!        VALJ(JVL_,NJVAL)  JVL_ = no of levels                          
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
!      real*8, intent(in)  ::  PPJ(L1_+1),TTJ(L1_) 
      real*8, intent(in)  ::  PPJ(L_),TTJ(L_)  !set to pfull (hPa) and tfull from AM3
      real*8, intent(in)  ::  FFF(W_,JVL_) 
      real*8, intent(out) ::  VALJL(JVL_,NJVAL) 
                                                                        
                                ! temp for call J's at one L            
      real*8  VALJ(X_) 
      real*8  QO2TOT(W_), QO3TOT(W_), QO31DY(W_), QO31D, QQQT, TFACT 
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B 
      integer J,K,L, IV 
                                                                        
                       ! master loop over layer = L                     
      do L = 1,JVL_ 
                                                                        
!---need temperature, pressure, and density at mid-layer (for some quant
          TT   = TTJ(L1_-L) 
!         if (L .eq. 1) then 
          PP = PPJ(L1_-L)  !PPJ(1) 
!         else 
!          PP  = (PPJ(L)+PPJ(L+1))*0.5d0 
!         endif 
          DD = 7.24e18*PP/TT   !molec/cm3
                                                                        
        do J = 1,NJVAL 
          VALJ(J) = 0.d0 
        enddo !j
                                                                        
        if (TT .le. TQQ(2,1))  then 
          if (TT .le. TQQ(1,1))  then 
            TFACT = 0.d0 
          else 
            TFACT = (TT -TQQ(1,1))/(TQQ(2,1) -TQQ(1,1)) 
          endif 
          do K = 1,W_ 
           QO2TOT(K) = QO2(K,1) + (QO2(K,2) - QO2(K,1))*TFACT 
          enddo 
        else 
          if (TT .ge. TQQ(3,1))  then 
            TFACT = 1.d0 
          else 
            TFACT = (TT -TQQ(2,1))/(TQQ(3,1) -TQQ(2,1)) 
          endif 
          do K = 1,W_ 
           QO2TOT(K) = QO2(K,2) + (QO2(K,3) - QO2(K,2))*TFACT 
          enddo 
        endif 

                                                                        
        if (TT .le. TQQ(2,2))  then 
          if (TT .le. TQQ(1,2))  then 
            TFACT = 0.d0 
          else 
            TFACT = (TT -TQQ(1,2))/(TQQ(2,2) -TQQ(1,2)) 
          endif 
          do K = 1,W_ 
           QO3TOT(K) = QO3(K,1) + (QO3(K,2) - QO3(K,1))*TFACT 
          enddo 
        else 
          if (TT .ge. TQQ(3,2))  then 
            TFACT = 1.d0 
          else 
            TFACT = (TT -TQQ(2,2))/(TQQ(3,2) -TQQ(2,2)) 
          endif 
          do K = 1,W_ 
           QO3TOT(K) = QO3(K,2) + (QO3(K,3) - QO3(K,2))*TFACT 
          enddo 
        endif 

                                                                        
        if (TT .le. TQQ(2,3))  then 
          if (TT .le. TQQ(1,3))  then 
            TFACT = 0.d0 
          else 
            TFACT = (TT -TQQ(1,3))/(TQQ(2,3) -TQQ(1,3)) 
          endif 
          do K = 1,W_ 
           QO31DY(K) = Q1D(K,1) + (Q1D(K,2) - Q1D(K,1))*TFACT 
          enddo 
        else 
          if (TT .ge. TQQ(3,3))  then 
            TFACT = 1.d0 
          else 
            TFACT = (TT -TQQ(2,3))/(TQQ(3,3) -TQQ(2,3)) 
          endif 
          do K = 1,W_ 
           QO31DY(K) = Q1D(K,2) + (Q1D(K,3) - Q1D(K,2))*TFACT 
          enddo 
        endif 

                                                                        
        do K = 1,W_ 
           QO31D  = QO31DY(K)*QO3TOT(K) 
          VALJ(1) = VALJ(1) + QO2TOT(K)*FFF(K,L) 
          VALJ(2) = VALJ(2) + (QO3TOT(K) -QO31D) *FFF(K,L)  !changed by jul
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L) 
        enddo 

                                                                        
        do J = 4,NJVAL 
                                                                        
          if (TQQ(2,J) .gt. TQQ(1,J)) then 
           TFACT = max(0.0,min(1.0,(TT-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)))) 
          else 
           TFACT = 0.d0 
          endif 
                                                                        
          do K = 1,W_ 
            QQQT    = QQQ(K,1,J) + (QQQ(K,2,J) - QQQ(K,1,J))*TFACT 
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L) 
          enddo 
                                                                        
!>>>>Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])      
          if (TITLEJ(J).eq.'MeVK  ') then 
            VALJ(J) = VALJ(J)/(1.0 + 1.67e-19*DD) 
          endif 
!>>>>Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])      
          if (TITLEJ(J).eq.'MEKeto') then 
            VALJ(J) = VALJ(J)/(1.0 + 0.80E-19*DD) 
          endif 
!>>>>Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])     
          if (TITLEJ(J).eq.'MGlyxl') then 
            VALJ(J) = VALJ(J)/(1.0 + 1.66e-19*DD) 
          endif 
                                                                        
        enddo 

                                                                        
      if (TITLEJ(NJVAL-1).eq.'Acet-a') then 
!--------------J-ref v8.3 includes Blitz ACETONE q-yields-------------- 
!---Acetone is a special case:   (as per Blitz et al GRL, 2004)         
!---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3                   
!---     62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3                
          VALJ(NJVAL-1) = 0.d0 
          VALJ(NJVAL)   = 0.d0 
!---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp fa
        IV    = NJVAL-1 
        TFACA = (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) 
        TFACA = max(0.0, min(1.0, TFACA)) 
!---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic 
!---      but force to Q2=0 by 210K                                     
        IV    = NJVAL 
        TFAC0 = ( (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) )**2 
        if (TT .lt. TQQ(1,IV)) then 
          TFAC0 = (TT - 210.d0)/(TQQ(1,IV)-210.d0) 
        endif 
        TFAC0 = max(0.0, min(1.0, TFAC0)) 
!---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-30
        IV    = NJVAL+1 
        TT200 = min(300.0, max(200.0, TT)) 
        TFAC1 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) 
!---IV=NJVAL+2 = Q1B for Acetone => (1)                                 
        IV    = NJVAL+2 
        TFAC2 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) 
                                                                        
!---now integrate over wavelengths                                      
        do K = 1,W_ 
!---NJVAL-1 = Xsect (total abs) for Acetone                             
          IV   = NJVAL-1 
          QQQA = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFACA 
!---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic i
          IV   = NJVAL 
          QQ2  = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC0 
          if (TT .lt. TQQ(1,IV)) then 
            QQ2 = QQQ(K,1,IV)*TFAC0 
          endif 
!---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K 
          IV   = NJVAL+1 
          QQ1A = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC1 
!---NJVAL+2 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19           
          IV   = NJVAL+2 
          QQ1B = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC2 
          QQ1B = QQ1B*4.d-20 
!---J(61)                                                               
          VALJ(NJVAL-1) = VALJ(NJVAL-1)                                 &
     &         + FFF(K,L)*QQQA*(1.d0-QQ2)/(QQ1A + QQ1B*DD)              
!---J(62)                                                               
          VALJ(NJVAL) = VALJ(NJVAL) + FFF(K,L)*QQQA*QQ2 
                                                                        
                 !K                                                     
        enddo 
!-----------end v-8.3 includes Blitz ACETONE q-yields--------------     
      endif 
                                                                        
!----Load array of J-values in native order, need to be indexed/scaled  
!    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ
        do J=1,NJVAL 
          VALJL(L,J) = VALJ(J) 
!         if(VALJ(J) .lt. 0.D0)then
!            write(*,*)'fastjx: !!!! negative J: J=',J, ' TITLE=',TITLEJ(J),'  VALJ(J)=',VALJ(J),' L=',L
!         endif   
        enddo 
                                                                        
                                          
      enddo ! master loop over L=1,JVL_   
      return 
      END subroutine JRATET                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
!      subroutine JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,DTAUX,POMEGAX,JXTRA) 
!-----------------------------------------------------------------------
!      implicit none 
!      include 'parm_CTM.f' 
!-----------------------------------------------------------------------
!--------key amtospheric data needed to solve plane-parallel J----------
!      real*8, dimension(L1_+1) :: TTJ,DDJ,ZZJ,ZHL 
!      real*8, dimension(L1_+1) :: PPJ 
!      integer,dimension(L2_+1) :: JXTRA 
!      real*8                 :: ZZHT 
!      real*8                 :: DTAUX(L1_),POMEGAX(8,L1_) 
!-----------------------------------------------------------------------
!      integer  I,J,K,L 
!      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP 
                                                                        
!      write(6,'(4a)') '   L z(km)     p      T   ',                     &
!     & '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',    &
!     & '  g(cos) CTM lyr=>'                                             
                                                                        
!          COLO2 = 0.d0 
!          COLO3 = 0.d0 
!          ZTOP = ZHL(L1_) + ZZHT 
                                                                        
!        do L = L1_,1,-1 
!          COLO2 = COLO2 + DDJ(L)*0.20948d0 
!          COLO3 = COLO3 + ZZJ(L) 
!          DELZ = ZTOP-ZHL(L) 
!          ZTOP = ZHL(L) 
!          ZKM = ZHL(L)*1.d-5 
                                                                        
!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!     &      L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ,                &
!     &      COLO2,COLO3,DTAUX(L),POMEGAX(1,L),POMEGAX(2,L)/3.d0,        &
!     &      JXTRA(L+L),JXTRA(L+L-1)                                     
                                                                        
                                                                        
!        enddo 
                                                                        
!      return 
!      END  subroutine JP_ATM                                         
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine RD_XXX(NAMFIL) 
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.   
!                                                                       
!>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-onl
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengt
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).     
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) s
!                                                                       
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J
!     NJ1      Channel number for reading data file                     
!                                                                       
!     NJVAL    Number of species to calculate J-values for              
!     NWWW     Number of wavelength bins, from 1:NWWW                   
!     WBIN     Boundaries of wavelength bins                            
!     WL       Centres of wavelength bins - 'effective wavelength'      
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)      
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)      
!     QO2(18,3)      O2 cross-sections                                        
!     QO3(18,3)      O3 cross-sections                                        
!     Q1D(18,3)      O3 => O(1D) quantum yield                                
!     TQQ(3,64)      Temperature for supplied cross sections                  
!     QQQ      Supplied cross sections in each wavelength bin (cm2)   
!  
!-----------------------------------------------------------------------
!jul++
! RD_XXX  
! will access WL(18), FL(18), QRAYL(18)
!             QO2(18,3)      O2 cross-sections                                        
!             QO3(18,3)      O3 cross-sections                                        
!             Q1D(18,3)      O3 => O(1D) quantum yield   
!             TITLEJ(64),   SPECIES NAME  O2 O3
!             TQQ(3, 64),   TEMPERATURE 180k, 260k, 300K
!             QQQ(18,2,64)
!jul--
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
      character(*), intent(in) ::  NAMFIL 
                                                                        
      integer  NJ1 
      integer  I, J, JJ, K, IW, NQQQ, NWWW, NQRD 
                                                                        
      TQQ(:,:) = 0.d0 
                                                                        
!----------spectral data----set for new format data J-ver8.3------------
!         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects rea
!         for 2005a data, NJVAL = 62 (including a spare XXXX) and       
!              NQQQ = 64 so that 4 wavelength datasets read in for aceto
!         note NQQQ is not used outside this subroutine!                
                                                                        
! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-s
                                                                        
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
      read (NJ1,100) TITLE0 
      read (NJ1,101) NJVAL,NQRD, NWWW 
         NW1 = 1 
         NW2 = NWWW 
      if (NJVAL.gt.X_ .or. NQRD.gt.X_) then 
        write(6,201) 'ATMOS:fastjx_init: RD_XXX stop',NJVAL,X_ 
        stop 
      endif 
!      write(6,'(1X,A)') TITLE0 
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects                       
      read (NJ1,102) (WL(IW),IW=1,NWWW) 
      read (NJ1,102) (FL(IW),IW=1,NWWW) 
      read (NJ1,102) (QRAYL(IW),IW=1,NWWW) 
                                                                        
!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps
      read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW) 
      read (NJ1,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW) 
      read (NJ1,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW) 
                                                                        
      read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW) 
      read (NJ1,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW) 
      read (NJ1,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW) 
                                                                        
      read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW) 
      read (NJ1,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW) 
      read (NJ1,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW) 
                                                                        
      do J = 1,3 
              if (mpp_pe() == mpp_root_pe()) then 
                write(*,*)'ATMOS:fastjx_init: RD_XXX:'
                write(*,200) J,TITLEJ(J),(TQQ(I,J),I=1,3) 
              end if
      enddo 
                                                                        
!---Read remaining species:  X-sections at 2 T_s                        
        JJ = 4 
      do J = 4,NQRD 
        read (NJ1,103) TITLEJ(JJ),TQQ(1,JJ),(QQQ(IW,1,JJ),IW=1,NWWW) 
        read (NJ1,103) TITLEJ2,  TQQ(2,JJ),(QQQ(IW,2,JJ),IW=1,NWWW) 
                                                                        
        if (W_.eq.18 .or. TITLEJ2(7:7).ne.'x') then 
!---include stratospheric J's (this also includes Cl and Br compounds!) 
          if (mpp_pe() == mpp_root_pe()) then 
            write(*,200) JJ,TITLEJ(JJ),(TQQ(I,JJ),I=1,2) 
          end if  
          JJ = JJ+1 
        endif 
                                                                        
      enddo 
       NQQQ = JJ-1 
       NJVAL = NJVAL + (NQQQ - NQRD) 
                                                                        
!---truncate number of wavelengths to do troposphere-only               
      if (W_ .ne. WX_) then 
!---TROP-ONLY                                                           
       if (W_ .eq. 12) then 
!        write(6,'(a)')                                                  &
!        ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'   
        NW2 = 12 
        do IW = 1,4 
          WL(IW) = WL(IW+4) 
          FL(IW) = FL(IW+4) 
          QRAYL(IW) = QRAYL(IW+4) 
         do K = 1,3 
          QO2(IW,K) = QO2(IW+4,K) 
          QO3(IW,K) = QO3(IW+4,K) 
          Q1D(IW,K) = Q1D(IW+4,K) 
         enddo 
         do J = 4,NQQQ 
          QQQ(IW,1,J) = QQQ(IW+4,1,J) 
          QQQ(IW,2,J) = QQQ(IW+4,2,J) 
         enddo 
        enddo 
        do IW = 5,12 
          WL(IW) = WL(IW+6) 
          FL(IW) = FL(IW+6) 
          QRAYL(IW) = QRAYL(IW+6) 
         do K = 1,3 
          QO2(IW,K) = QO2(IW+6,K) 
          QO3(IW,K) = QO3(IW+6,K) 
          Q1D(IW,K) = Q1D(IW+6,K) 
         enddo 
         do J = 4,NQQQ 
          QQQ(IW,1,J) = QQQ(IW+6,1,J) 
          QQQ(IW,2,J) = QQQ(IW+6,2,J) 
         enddo 
        enddo 
!---TROP-QUICK  (must scale solar flux for W=5)                         
       elseif (W_ .eq. 8) then 
!        write(6,'(a)')                                                  &
!        ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'   
        NW2 = 8 
        do IW = 1,1 
          WL(IW) = WL(IW+4) 
          FL(IW) = FL(IW+4)  * 2.d0 
          QRAYL(IW) = QRAYL(IW+4) 
         do K = 1,3 
          QO2(IW,K) = QO2(IW+4,K) 
          QO3(IW,K) = QO3(IW+4,K) 
          Q1D(IW,K) = Q1D(IW+4,K) 
         enddo 
         do J = 4,NQQQ 
          QQQ(IW,1,J) = QQQ(IW+4,1,J) 
          QQQ(IW,2,J) = QQQ(IW+4,2,J) 
         enddo 
        enddo 
        do IW = 2,8 
          WL(IW) = WL(IW+10) 
          FL(IW) = FL(IW+10) 
          QRAYL(IW) = QRAYL(IW+10) 
         do K = 1,3 
          QO2(IW,K) = QO2(IW+10,K) 
          QO3(IW,K) = QO3(IW+10,K) 
          Q1D(IW,K) = Q1D(IW+10,K) 
         enddo 
         do J = 4,NQQQ 
          QQQ(IW,1,J) = QQQ(IW+10,1,J) 
          QQQ(IW,2,J) = QQQ(IW+10,2,J) 
         enddo 
        enddo 
                                                                        
       else 
         write(6,*) 'ATMOS:fastjx_init: number of used wavelengths wrong:',W_ 
         stop 
       endif 
      endif 
                                                                        
!  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s       
!   61: C3H6O  = Acet-a     (CH3CO + CH3)                               
!   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)                            
                                                                        
      TITLEJ(NJVAL-1) = 'Acet-a' 
      TITLEJ(NJVAL)   = 'Acet-b' 
                                                                        
      call mpp_close (NJ1)
                                                                        
  100 format(a) 
  101 format(10x,5i5) 
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3)) 
  103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3)) 
  200 format(1x,' x-sect:',i3,a10,3(3x,f6.2)) 
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/,        &
            ' Maximum number allowed (X_) only set to: ',i3,           &
            ' - increase in cmn_jv.f')                                 
                                                                        
      return 
      END subroutine RD_XXX                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine RD_MIE(NAMFIL) 
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)       
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<< 
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     NJ1      Channel number for reading data file                     
!     NAA      Number of categories for scattering phase functions      
!     QAA      Aerosol scattering phase functions                       
!     NK       Number of wavelengths at which functions supplied (set as
!     WAA      Wavelengths for the NK supplied phase functions          
!     PAA      Phase function: first 8 terms of expansion               
!     RAA      Effective radius associated with aerosol type            
!     SAA      Single scattering albedo                                 
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
      character(*), intent(in) ::  NAMFIL 
      integer  I, J, K, NJ1 
!jul++
!   NAA,TITLE0
!   JTAUMX,ATAU,ATAU0 
!   TITLAA(J),RAA(J),DAA(J) 
!   WAA(K,J),QAA(K,J),SAA(K,J),PAA(I,K,J)
!
!jul--                                                                        
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
                                                                        
      read (NJ1,'(i2,a78)') NAA,TITLE0 
        if (NAA .gt. A_) then 
          write(*,*) 'ATMOS:fastjx_init: too many scat-data sets:', NAA, A_ 
          stop 
        endif 
      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0 
      if (mpp_pe() == mpp_root_pe()) then 
               write(*,*)'ATMOS:fastjx_init: RD_MIE:'
               write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX 
      end if
        
      read (NJ1,*) 
      do J = 1,NAA 
          read (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)')                       &
     &           TITLAA(J),RAA(J),DAA(J)                                
                       ! ver 6.0 extend to 5 ref wavelengths for mie-sca
        do K = 1,5 
          read (NJ1,'(f4.0,f7.4,f7.4,7f6.3,1x,f7.3,f8.4)')              &
     &  WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)                   
          PAA(1,K,J) = 1.d0 
        enddo 
      enddo 
                                                                        
      call mpp_close (NJ1)
      if (mpp_pe() == mpp_root_pe()) then 
        write(*,'(a,9f8.1)') ' ATMOS:fastjx_init: RD_MIE: Aerosol optical: r-eff/rho/Q(@wavel):'     &
                  ,(WAA(K,1),K=1,5)                                    
        write(*,*) TITLE0 
      end if                                                                        

      do J=1,NAA 
        if (mpp_pe() == mpp_root_pe()) then 
           write(*,*) ' ATMOS:fastjx_init: RD_MIE:'
           write(*,'(i3,1x,a8,7f8.3)')                                 &
                   J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)   
        end if          
      enddo 
                                                                        
      return 
      END  subroutine RD_MIE                                         
 
 
 
 
                                                                       
!-----------------------------------------------------------------------
      subroutine RD_MIE_AM3(NAMFIL) 
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)       
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<< 
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     NJ1      Channel number for reading data file                     
!     NAA      Number of categories for scattering phase functions      
!     QAA      Aerosol scattering phase functions                       
!     NK       Number of wavelengths at which functions supplied (set as
!     WAA      Wavelengths for the NK supplied phase functions          
!     PAA      Phase function: first 8 terms of expansion               
!     RAA      Effective radius associated with aerosol type            
!     SAA      Single scattering albedo                                 
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
      character(*), intent(in) ::  NAMFIL 
      integer  I, J, K, NJ1 
      character*90 :: str1
!jul++
!   NAA,TITLE0
!   JTAUMX,ATAU,ATAU0 
!   TITLAA(J),RAA(J),DAA(J) 
!   WAA(K,J),QAA(K,J),SAA(K,J),PAA(I,K,J)
!
!jul--                                                                        
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
                                                                        
      read (NJ1,'(i4,a78)') NAA_AM3, TITLE0_AM3 
      if (NAA_AM3 .gt. A_AM3) then 
          write(*,*) 'ATMOS:fastjx_init: too many scat-data sets for AM3:', NAA_AM3, A_AM3 
          stop 
      endif 
!      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0 
!      if (mpp_pe() == mpp_root_pe()) then 
!               write(*,*)'ATMOS:fastjx_init: RD_MIE_AM3:'
!               write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX 
!      end if
        
!      read (NJ1,*) 
      do J = 1,NAA_AM3
          read (NJ1,'(a90)') str1
          if (mpp_pe() == mpp_root_pe()) then 
             write (*,'(a12,a90)')              &
                   'fastjx_init:',str1               
          endif       
!          read (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)')                       &
!     &           TITLAA(J),RAA(J),DAA(J)                                
                       ! ver 6.0 extend to 5 ref wavelengths for mie-sca
        do K = 1,5 
          read (NJ1,'(f4.0,2e14.6,7f6.3,1x,f7.3,f8.4)')              &
                   WAA_AM3(K,J),MEE_AM3(K,J),SAA_AM3(K,J),(PAA_AM3(I,K,J),I=2,8)                   
          PAA_AM3(1,K,J) = 1.d0 
  
          if (mpp_pe() == mpp_root_pe()) then 
!!!              write (*,'(a12,i4,1x,f4.0,2e14.6,8f6.3,1x,f7.3,f8.4)')              &
!!!                  'fastjx_init:',J,WAA_AM3(K,J),MEE_AM3(K,J),SAA_AM3(K,J),(PAA_AM3(I,K,J),I=1,8)                   
          endif
        enddo 
      enddo 
                                                                        
       call mpp_close (NJ1)
!      if (mpp_pe() == mpp_root_pe()) then 
!        write(*,'(a,9f8.1)') ' ATMOS:fastjx_init: RD_MIE: Aerosol optical: r-eff/rho/Q(@wavel):'     &
!                  ,(WAA_AM3(K,1),K=1,5)                                    
!!        write(*,*) TITLE0_AM3 
!      end if                                                                        
!        if (mpp_pe() == mpp_root_pe()) then 

!      do J=1,NAA_AM3 
!          write(*,*) ' ATMOS:fastjx_init: RD_MIE_AM3:'
!           write(*,'(i3,1x,a8,7f8.3)')                                 &
!                   J,TITLAA_AM3(J),(QAA(K,J),K=1,5)   
!        end if          
!      enddo 
                                                                        
      return 
      END  subroutine RD_MIE_AM3                                         
 
 
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine RD_UM(NAMFIL) 
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)               
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)        
!     NJ1      Channel number for reading data file                     
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 

!jul++
!  TITLUM(L)
!  UMAER
!jul--                                                                        
      character(*), intent(in) ::  NAMFIL 
                                                                        
      integer  I, J, K, L, NJ1 
                                                                        
      call mpp_open (NJ1, trim(NAMFIL), MPP_RDONLY, MPP_ASCII,  &
                     MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)
                                                                        
      read (NJ1,'(a78)') TITLE0 
!        write(6,*) 'UMichigan Aerosol optical data' 
!        write(6,*) TITLE0 
                                                                        
!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, 
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC) 
      do L=1,33 
          read(NJ1,'(a4)') TITLUM(L) 
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%                  
        do K=1,21 
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 4=600nm, 5=10
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext                       
          read(NJ1,'(9f9.5,27x,6f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,5) 
        enddo 
      enddo 
                                                                        
      call mpp_close (NJ1)
                                                                        
!        write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33) 
                                                                        
      return 
      END subroutine RD_UM                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      real*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3) 
!-----------------------------------------------------------------------
!  Three-point linear interpolation function                            
!-----------------------------------------------------------------------
      real*8  TINT,T1,T2,T3,F1,F2,F3 
      if (TINT .le. T2)  then 
        if (TINT .le. T1)  then 
          FLINT = F1 
        else 
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1) 
        endif 
      else 
        if (TINT .ge. T3)  then 
          FLINT = F3 
        else 
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2) 
        endif 
      endif 
      return 
      END  FUNCTION FLINT                                         
                                                                        
                                                                        
!-----------------------------------------------------------------------
!      subroutine SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX) 
!-----------------------------------------------------------------------
!     GMTIME = UT for when J-values are wanted                          
!           (for implicit solver this is at the end of the time step)   
!     NDAY   = integer day of the year (used for solar lat and declin)  
!     YGRDJ  = laitude (radians) for grid (I,J)                         
!     XGDRI  = longitude (radians) for grid (I,J)                       
!                                                                       
!     SZA = solar zenith angle in degrees                               
!     COSSZA = U0 = cos(SZA)                                            
!-----------------------------------------------------------------------
!      implicit none 
!      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI 
!      integer, intent(in) ::  NDAY 
!      real*8, intent(out) ::  SZA,COSSZA,SOLFX 
                                                                        
!      real*8  PI, PI180, LOCT 
!      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ 
!                                                                       
!      PI     = 3.141592653589793d0 
!      PI180  = PI/180.d0 
!      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180) 
!      SOLDEK = asin(SINDEC) 
!      COSDEC = cos(SOLDEK) 
!      SINLAT = sin(YGRDJ) 
!      SOLLAT = asin(SINLAT) 
!      COSLAT = cos(SOLLAT) 
!                                                                       
!      LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI 
!      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT 
!      SZA    = acos(COSSZA)/PI180 
!       if (mpp_pe() == mpp_root_pe()) then 
!         write(*,*) 'FASTJ:SOLARZ: XGRDI,YGRDJ',XGRDI,YGRDJ 
!         write(*,*) 'FASTJ:SOLARZ: LOCT (rad)',LOCT 
!         write(*,*) 'FASTJ:SOLARZ: SINDEC,COSDEC', SINDEC,COSDEC 
!         write(*,*) 'FASTJ:SOLARZ: SINLAT,COSLAT', SINLAT,COSLAT 
!         write(*,*) 'FASTJ:SOLARZ: COS, SZA',COSSZA,SZA 
!       end if                                                                        
                                                                        
!      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0)) 
                                                                        
!      return 
!      END   subroutine SOLARZ                                         
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine SPHERE2(GMU,RAD,ZHL,ZZHT,AMF2,L1_) 
!-----------------------------------------------------------------------
!----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90   
!  This new AMF2 does each of the half-layers of the CTM separately,    
!     whereas the original, based on the pratmo code did the whole layer
!     and thus calculated the ray-path to the CTM layre edges, NOT the m
!  Since fast-JX is meant to calculate the intensity at the mid-layer, t
!     solar beam at low sun (interpolated between layer edges) was incor
!  This new model does make some approximations of the geometry of the l
!     the CTM layer is split evenly in mass (good) and in height (approx
!                                                                       
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when          
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent  
!  beam (where tangent height is below altitude J-value desired at).    
!-----------------------------------------------------------------------
! in:                                                                   
!     GMU     = MU0 = cos(solar zenith angle)                           
!     RAD     radius of Earth mean sea level (cm)                       
!     ZHL(L)  height (cm) of the bottome edge of CTM level L            
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1)        
!     L1_     dimension of CTM = levels +1 (L+1 = above-CTM level)      
! out:                                                                  
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching 
!-----------------------------------------------------------------------
      implicit none 
      integer, intent(in) ::   L1_ 
      real*8, intent(in)  ::   GMU,RAD,ZHL(L1_+1),ZZHT 
      real*8, intent(out) ::   AMF2(2*L1_+1,2*L1_+1) 
                                                                        
!     RZ      Distance from centre of Earth to each point (cm)          
!     RQ      Square of radius ratios                                   
!     SHADHT  Shadow height for the current SZA                         
!     XL      Slant path between points                                 
                                                                        
      integer  I, J, K, II, L2 
      real*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(L1_+1) 
      real*8   RZ2(2*L1_+1),RQ2(2*L1_+1) 
!                                                                       
!--- must have top-of-atmos (NOT top-of-CTM) defined                    
!      ZHL(L1_+1) = ZHL(L1_) + ZZHT                                     
                                                                        
        RZ(1) = RAD + ZHL(1) 
      do II = 2,L1_+1 
        RZ(II)   = RAD + ZHL(II) 
      enddo 
                                                                        
!---calculate heights for edges of split CTM-layers                     
      L2 = 2*L1_ 
      do II = 2,L2,2 
        I = II/2 
        RZ2(II-1) = RZ(I) 
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1)) 
      enddo 
        RZ2(L2+1) = RZ(L1_+1) 
      do II = 1,L2 
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2 
      enddo 
                                                                        
                                                                        
!---shadow height for SZA > 90                                          
      if (GMU .lt. 0.0d0)  then 
        SHADHT = RZ2(1)/dsqrt(1.0d0-GMU**2) 
      else 
        SHADHT = 0.d0 
      endif 
                                                                        
!---up from the surface calculating the slant paths between each level  
!---  and the level above, and deriving the appropriate Air Mass Factor 
         AMF2(:,:) = 0.d0 
                                                                        
      do 16 J = 1,2*L1_+1 
                                                                        
!  Air Mass Factors all zero if below the tangent height                
        if (RZ2(J) .lt. SHADHT) goto 16 
!  Ascend from layer J calculating AMF2s                                
        XMU1 = abs(GMU) 
        do I = J,2*L1_ 
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2)) 
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1 
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I)) 
          XMU1     = XMU2 
        enddo 
!--fix above top-of-atmos (L=L1_+1), must set DTAU(L1_+1)=0             
          AMF2(2*L1_+1,J) = 1.d0 
!                                                                       
!  Twilight case - Emergent Beam, calc air mass factors below layer     
        if (GMU .ge. 0.0d0) goto 16 
                                                                        
!  Descend from layer J                                                 
          XMU1       = abs(GMU) 
         do II = J-1,1,-1 
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II) 
                                                ! filter                
          if (II.eq.1)  DIFF = max(DIFF,0.0) 
!  Tangent height below current level - beam passes through twice       
          if (DIFF .lt. 0.0d0)  then 
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II)) 
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2) 
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II)) 
            XMU1      = XMU2 
!  Lowest level intersected by emergent beam                            
          else 
            XL        = RZ2(II+1)*XMU1*2.0d0 
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II)) 
            goto 16 
          endif 
         enddo 
                                                                        
   16 continue 
      return 
      END  subroutine SPHERE2                                         
                                                                        
                                                                                                                                                
!-----------------------------------------------------------------------
      subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA) 
!-----------------------------------------------------------------------
!                                                                       
!    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol laye
!    this version sets up log-spaced sub-layers of increasing thickness ATAU
!                                                                       
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)   
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L                 
!        Set for log-spacing of tau levels, increasing top-down.        
!                                                                       
!     N.B. the TTAU, etc calculated here are NOT used elsewhere         
                                                                        
!---The log-spacing parameters have been tested for convergence and chosen
!---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0 
!---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!-----------------------------------------------------------------------
!                                                                       
      implicit none 
                                                
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol                                               
      integer, intent(in) ::  NX              !Mie scattering array size                                              
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0                                               
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!                                                                       
      integer JTOTL,I,L,L2 
      real*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1 
!                                                                       
!---Reinitialize arrays                                                 
      TTAU(:)  = 0.d0 
      JXTRA(:) = 0 
!                                                                       
!---combine these edge- and mid-layer points into grid of size:         
!---              L2X+1 = 2*L1X+1 = 2*L_+3                              
!---calculate column optical depths above each level, TTAU(1:L2X+1)     
!---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD                  
!                                                                       
!---Divide thick layers to achieve better accuracy in the scattering code
!---In the original fast-J, equal sub-layers were chosen, this is wasteful
!---and this new code (ver 5.3) uses log-scale:                         
!---        Each succesive layer (down) increase thickness by ATAU > 1  
!---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!---        4 sub-layers with ODs = 1 - 2 - 4 - 8                       
!---The key parameters are:                                             
!---        ATAU = factor increase from one layer to the next           
!---        ATAUMN = the smallest OD layer desired                      
!---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN
!---These are hardwired below, can be changed, but have been tested/optimized
                                                                        
      ATAU1  = ATAU - 1.d0 
      ATAULN = log(ATAU) 
        TTAU(L2X+1)  = 0.0d0 
      do L2 = L2X,1,-1 
        L         = (L2+1)/2 
        DTAUJ     = 0.5d0 * DTAUX(L) 
        TTAU(L2)  = TTAU(L2+1) + DTAUJ 
!---Now compute the number of log-spaced sub-layers to be added in      
!---   the interval TTAU(L2) > TTAU(L2+1)                               
!---The objective is to have successive TAU-layers increasing by factor ATAU>1
!---the number of sub-layers + 1                                        
        if (TTAU(L2) .lt. ATAU0) then 
          JXTRA(L2) = 0 
        else 
          ATAUM    = max(ATAU0, TTAU(L2+1)) 
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN 
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0))) 
        endif 
      enddo 
                                                                        
!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2 
      do L2 = L2X,1,-1 
        JTOTL  = JTOTL + JXTRA(L2) 
        if (JTOTL .gt. NX/2)  then 
          write(6,'(A,2I5,F9.2)') 'ATMOS:fastjx_photo: N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2 
          do L = L2,1,-1 
            JXTRA(L) = 0 
          enddo 
          go to 10 
        endif 
      enddo 
   10 continue 
                                                                        
      return 
      END subroutine EXTRAL                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,                &
     &                  FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)       
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_CTM.f' 
!      include 'parm_MIE.f' 
!      include 'cmn_JVdat.f' 
                                                                        
      real*8, intent(in)        ::      DTAUX(L1_,W_),&         ! optical depth
                                        POMEGAX(8,L1_,W_)       ! scattering phase fn
      real*8, intent(in)        ::      AMF2(2*L1_+1,2*L1_+1) 
      real*8, intent(in)        ::      U0,RFL(W_) 
      integer, intent(in)       ::      JXTRA(L2_+1) 
      real*8, intent(out)       ::      FJACT(L_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_) 
      real*8, intent(out)       ::      FJFLX(L_,W_),FLXD(L1_,W_),FLXD0(W_) 
!                                                                       
      integer JNDLEV(L_),JNELEV(L1_) 
      integer JADDLV(L2_+1),JADDTO(L2_+1),L2LEV(L2_+1) 
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND 
      integer LZ0,LZ1,LZMID 
      real*8   SUMT,SUMJ 
                                                                        
      real*8  DTAU(L1_+1,W_),POMEGAJ(M2_,L2_+1,W_),TTAU(L2_+1,W_) 
      real*8  FTAU2(L2_+1,W_),POMEGAB(M2_,W_) 
      real*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0 
      real*8, dimension(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX 
!--- variables used in mie code-----------------------------------------
      real*8, dimension(W_)         :: FJT,FJB 
      real*8, dimension(N_,W_)      :: FJ,FZ,ZTAU 
      real*8, dimension(M2_,N_,W_) :: POMEGA 
      real*8, dimension(2*L1_,W_)   :: FLXD2 
                                                                        
!  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts     
!                                                                       
! in:                                                                   
!     DTAUX(1:L1_,1:W_) = optical depth of each layer                   
!     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s abledo)
!     U0  = cos (SZA)                                                   
!     RFL(1:W_) = Lambertian albedo of surface                          
!     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!        AMF2 now does both edges and middle of CTM layers              
!     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
! out:                                                                  
!     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
!  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
!     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)        
!     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)  
!     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L       
!        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
!        this should take into account sphericity, and is not just = mu0
!     FLXD0(1:W_) = sum of solar flux deposited in atmos                
!        does NOT include flux on lower surface, does NOT mean absorbed!
!-----------------------------------------------------------------------
!                                                                       
!     DTAU     Local optical depth of each CTM level                    
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU2     Attenuation of solar beam                               
!     POMEGAJ  Scattering phase function                                
!                                                                       
!---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the  
!   factor increase from sub-layer to sub-layer                         
!                                                                       
!---------------------SET UP FOR MIE CODE-------------------------------
                                                                        
!-----------------wavelength independent--------------------------------
!                                                                       
!  Transpose the ascending TTAU grid to a descending ZTAU grid.         
!  Double the resolution - TTAU points become the odd points on the     
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.  
!  Odd point added at top of grid for unattenuated beam   (Z='inf')     
!                                                                       
!  The following mapping holds for JADDLV=0                             
!        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)                        
!        Top:       TTAU(L2_)  ==> ZTAU(3)                              
!        Infinity:     0.0     ==> ZTAU(1)                              
!        index: 2*(L2_+1-L2)+1 ==> LZ                                   
!                                                                       
!  Mie scattering code only used from surface to level L2_              
!-----------------------------------------------------------------------
!                                                                       
!-----------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere  
!  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!  to be incremented linearly, and the flux fz to be attenuated top-down
!    (avoiding problems where lower level fluxes are zero).             
!-----------------------------------------------------------------------
!                                                                       
!  Ascend through atmosphere transposing grid and adding extra points   
!  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!    because we need to insert the intermediate layers (even LZ) for the
!    asymmetric scattering code.                                        
                                                                        
                                                                        
!  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse 
!    order, expanded, doubled-level scatter grid.                       
!    Note that we need to deal with the expansion by JADD levels (L2L). 
!      These JADDLV levels are skipped and need to be interpolated later
!    Note that only odd LZ levels are filled,                           
                                                                        
!----------------------re-grid data-------------------------------------
!  Calculate cumulative total and define levels we want J-values at.    
!  Sum upwards for levels, and then downwards for Mie code readjustments
!                                                                       
!     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)    
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDLV(L2)  Number of new levels actually added at each wavelength
!            where JADDLV = 0 when there is effectively no FTAU2        
!     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!     JNDLEV(L) = L2 index that maps on CTM mid-layer L                 
!                                                                       
!---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!---    JADDLV is taken from JXTRA, which is based on visible OD.       
!---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!---these should be fixed for all wavelengths to lock-in the array sizes
      do L2 = 1,L2_,1 
        JADDLV(L2) = JXTRA(L2) 
      enddo 
        JADDTO(L2_+1) = 0 
      do L2 = L2_,1,-1 
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2) 
      enddo 
                                                                        
!---expanded grid now included CTM edge and mid layers plus expanded    
!---    grid to allow for finer delta-tau at tops of clouds.            
!---    DIM of new grid = L2_ + JADDTO(1) + 1                           
                                                                        
!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV) 
!     in absence of JADDLV, L2LEV(L2) = L2                              
        L2LEV(1)  = 1 
      do L2 = 2,L2_+1 
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1) 
      enddo 
                                                                        
!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L      
!---JNELEV(L=1:L_) = L2-index for top of layer L                        
      do L = 1,L_ 
        JNDLEV(L) = L2LEV(2*L) 
        JNELEV(L) = L2LEV(2*L+1) 
      enddo 
                          !need to set this to top-of-atmosphere        
        JNELEV(L_+1) = 0 
                                                                        
      ND = 2*L2_ + 2*JADDTO(1) + 1 
                                                                        
      if(ND .gt. N_) then 
        write(6,'(a,2i9)') 'ATMOS:fastjx_photo: overflow of scatter arrays:',ND,N_ 
        stop 
      endif 
                                                                        
!----------------begin wavelength dependent set up----------------------
                                                                        
!---Reinitialize arrays                                                 
      ZTAU(:,:)     = 0.d0 
      FZ(:,:)       = 0.d0 
      POMEGA(:,:,:) = 0.d0 
                                                                        
      do K=1,W_ 
                                                                        
!---Set up optical depth DTAU(L)                                        
       do L = 1,L1_ 
        DTAU(L,K) = DTAUX(L,K) 
       enddo 
        DTAU(L1_+1,K) = 0.d0 
                                                                        
!---Define the total scattering phase fn for each CTM layer L=1:L_+1    
!---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh               
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8     
       do L = 1,L1_ 
        do I = 1,M2_ 
          POMEGAJ(I,L,K) = POMEGAX(I,L,K) 
        enddo 
       enddo 
                                                                        
!---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
!---      at the middle & edges of the CTM layers L=1:2*L1_+1           
!---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0       
!---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0.d0 
        FTAU2(L2_+1,:) = 1.0d0 
       do LL = 1,2*L1_+1 
         L = (LL+1)/2 
        if (AMF2(LL,LL) .gt. 0.0d0) then 
           XLTAU = 0.0d0 
         do II = 1,2*L1_+1 
           I = (II+1)/2 
           XLTAU = XLTAU + 0.5d0*DTAU(I,K)*AMF2(II,LL) 
         enddo 
                                      ! zero out flux at 1e-33          
         if (XLTAU .lt. 76.d0) then 
          FTAU2(LL,K) = exp(-XLTAU)    !!!!!!!!good place to
         endif 
        endif 
       enddo 
                                                                        
!---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)   
          FLXD2(:,:) = 0.d0 
       do LL = 1,2*L1_ 
        if (AMF2(LL,LL) .gt. 0.d0) then 
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL) 
        endif 
       enddo 
        if (AMF2(1,1) .gt. 0.d0) then 
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1) 
        else 
          FSBOT(K) = 0.d0 
        endif 
                                                                        
       do LL = 2,2*L1_,2 
         L=LL/2 
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K) 
       enddo 
                                                                        
!---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0 
       if (AMF2(2*L1_,2*L1_) .gt. 0.d0) then 
        do L=1,L1_ 
         FLXD0(K) = FLXD0(K) + FLXD(L,K) 
        enddo 
       endif 
                                                                        
!-----------------------------------------------------------------------
!  Take optical properties on CTM layers and convert to a photolysis    
!  level grid corresponding to layer centres and boundaries. This is    
!  required so that J-values can be calculated for the centre of CTM    
!  layers; the index of these layers is kept in the JNDLEV array.       
!-----------------------------------------------------------------------
!---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer     
!---    points (1:L_) plus 1 for the mid point of added top layer.      
!---combine these edge- and mid-layer points into grid of size:         
!---              L2_+1 = 2*L1_+1 = 2*L_+3                              
!---calculate column optical depths above each level, TTAU(1:L2_+1)     
!---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD                  
                                                                        
        TTAU(L2_+1,K) = 0.0d0 
       do L2 = L2_,1,-1 
        L          = (L2+1)/2 
        DTAUJ      = 0.5d0 * DTAU(L,K) 
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ 
       enddo 
                                                                        
!----solar flux incident on lower boundary & Lambertian reflect factor: 
       if (FSBOT(K) .gt. 0.d0) then 
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K)) 
       else 
        ZFLUX(K) = 0.d0 
       endif 
                                                                        
!  Calculate scattering properties, level centres then level boundaries 
!>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upw
       do L2 = L2_,2,-2 
        L   = L2/2 
        do I = 1,M2_ 
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K) 
        enddo 
       enddo 
!---lower boundary value is set (POMEGAJ(I,1), but set upper:           
       do I = 1,M2_ 
         POMEGAJ(I,L2_+1,K) = POMEGAJ(I,L2_,K) 
       enddo 
!---now have POMEGAJ filled at even points from L2=3:L2_-1              
!---use inverse interpolation for correct tau-weighted values at edges  
       do L2 = 3,L2_-1,2 
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K) 
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K) 
        do I = 1,M2_ 
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN +                  &
     &           POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)               
        enddo 
       enddo 
                                                                        
!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)              
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface     
! Ftau2 indicate how much of light is still alive at that layer                                                                        
                                
       do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers                               
        L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)                                     
        LZ  = ND + 2 - 2*L2L    ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K) 
          FZ(LZ,K)   = FTAU2(L2,K) 
        do I=1,M2_ 
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K) 
        enddo 
       enddo 
                                                                        
!   Now go thru the pairs of L2 levels to see if we need JADD levels    
                                 ! L2 = index of CTM edge- and mid-layer
       do L2 = 1,L2_ 
                                 ! L2L = index for L2 in expanded scale(
          L2L = L2LEV(L2) 
                                ! LZ = index for L2 in scatt arrays     
          LZ  = ND + 2 - 2*L2L 
                                             ! L22 = 0 if no added level
          L22 = L2LEV(L2+1) - L2LEV(L2) - 1 
                                                                        
          if (L22 .gt. 0) then 
             TAUBTM(K) = TTAU(L2,K) 
             TAUTOP(K) = TTAU(L2+1,K) 
             FBTM(K)   = FTAU2(L2,K) 
             FTOP(K)   = FTAU2(L2+1,K) 
             do I = 1,M2_ 
                POMEGAB(I,K) = POMEGAJ(I,L2,K) 
             enddo 
                                                                        
!---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU 
!---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM       
             ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1)) 
                                ! add odd levels between L2LEV(L2) & L2L
             do L = 1,L22 
                               ! LZZ = index(odd) of added level in scat
                LZZ = LZ - 2*L 
                ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ 
                                                                        
!---fraction from TAUBTM=>TAUTOP                                        
                ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K)) 
!---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg),   
!---else scale by TAU                                                   
                if (U0 .gt. 0.02d0) then 
                    FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0) 
                else 
                    if (FBTM(K) .lt. 1.d-32) then 
                        FZ(LZZ,K) = 0.d0 
                    else 
                        FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA 
                    endif 
                endif 
                do I = 1,M2_ 
                   POMEGA(I,LZZ,K) = POMEGAB(I,K) +                            &
     &               ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))             
                enddo 
                TAUBTM(K)    = ZTAU(LZZ,K) 
                FBTM(K)      = FZ(LZZ,K) 
                do I = 1,M2_ 
                   POMEGAB(I,K) = POMEGA(I,LZZ,K) 
                enddo 
             enddo !Loop L
          endif 
       enddo ! Loop L2
                                                                        
!   Now fill in the even points with simple interpolation in scatter arrays
       do LZ = 2,ND-1,2 
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K)) 
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K)) 
        do I=1,M2_ 
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K)) 
        enddo 
       enddo 
                                                                        
                                                     
      enddo ! wavelength loop! 
                                                                        
!-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND ) 
!-----------------------------------------------------------------------
                                                                        
!---Move mean intensity from scatter array FJ(LZ=1:ND)                  
!---              to CTM mid-level array FJACT(L=1:L_)                  
                                                                        
      do K=1,W_ 
                                                                        
!---mean intensity:  4*<I> + solar at mid-layer                         
       do L = 1,L_ 
        L2L = JNDLEV(L) 
        LZ  = ND+2 - 2*L2L 
        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K) 
       enddo 
       
 
                                                                        
!---mean diffuse flux:  4<I*mu> (not solar) at top of layer L           
!---      average (tau-wtd) the h's just above and below the L-edge     
       do L = 1,L_ 
        L2L = JNELEV(L) 
        LZ  = ND+2 - 2*L2L 
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K)) 
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0)) 
       enddo 
                                                                        
!---NB if one needs the mean intensity throughout layer L (instead of mi
!---   then average (tau-weighted) the odd-points from: NELEV(L-1) to NE
!---NB This is NOT now used.                                            
!      real*8   FJACT2(L_,W_)                                           
!-----LZ are indices to the 1:ND array     LZ1 > LZMID > LZ0            
!       LZ1 = ND                                                        
!      do L = 1,L_                                                      
!       LZMID = ND+2-2*JNDLEV(L)                                        
!       LZ0 = ND+2-2*JNELEV(L)                                          
!         SUMT = 0.d0                                                   
!         SUMJ = 0.d0                                                   
!       do L2 = LZ0,LZ1-2,2                                             
!         SUMT = SUMT + ZTAU(L2+2,K)-ZTAU(L2,K)                         
!         SUMJ = SUMJ + (ZTAU(L2+2,K)-ZTAU(L2,K))*(FJ(L2,K)+FJ(L2+2,K)) 
!       enddo                                                           
!         FJACT2(L,K) = 2.d0*SUMJ/SUMT + FZ(LZMID,K)                    
!         LZ1 = LZ0                                                     
!      enddo                                                            
                                                                        
!---diffuse fluxes reflected at top, incident at bottom                 
         FJTOP(K) = FJT(K) 
         FJBOT(K) = FJB(K) 
                                                                        
             ! wavelength loop!                                         
      enddo 
 
!if( any(FJACT(:,:) < 0.0D0 )) then
!  write(*,*) 'Mie: find negative FJACT'
!  write(*,*) 'Mie: U0=', U0, ' ND=',ND
!  if (mpp_pe() == mpp_root_pe() ) then   
!   write(*,*) 'Mie: FJ diffusive', FJ
!   write(*,*) 'Mie: FZ solar at mid-layer', FZ
!   write(*,*) 'Mie: ZTAU', ZTAU
!   write(*,*) 'Mie: RFL', RFL
!   write(*,*) 'Mie: FJT', FJT
!   write(*,*) 'Mie: FJB', FJB
!   write(*,*) 'Mie: POMEGA', POMEGA
!   write(*,*) 'Mie: ZFLUX', ZFLUX
!   write(*,*) 'Mie: AMF2', AMF2
!   endif
!end if   
 
                                                                        
      return 
      END  subroutine OPMIE                                         
                                                                        
                                                                        
!<<<<<<<<<<<<<<<<<<<<<<<end core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
                                                                        
                                                                        
!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND) 
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_MIE.f' 
                                                                        
      integer, intent(in) ::  ND 
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_)   &
     &                       ,RFL(W_),U0,ZFLUX(W_)                      
      real*8, intent(out) ::  FJ(N_,W_),FJT(W_),FJB(W_) 
                                                                        
      real*8  PM(M_,M2_),PM0(M2_)          
      integer I, IM  ,K 
                                                                        
!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.                        
!         Sol_n of inhomogeneous Rayleigh scattering atmosphere.        
!         (original Rayleigh w/ polarization)                           
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.                   
!         Raman scattering in the atmospheres of the major planets.     
!         (first use of anisotropic code)                               
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002
!         Chemistry of a polluted cloudy boundary layer,                
!         (documentation of extension to anisotropic scattering)        
!                                                                       
!    takes atmospheric structure and source terms from std J-code       
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)  
!-----------------------------------------------------------------------
      do I = 1,M_                         ! M_ =4 fixed
       call LEGND0 (EMU(I),PM0,M2_)       ! EMU =gauss weighting, M2_=8, return PM0
       do IM = 1,M2_ 
         PM(I,IM) = PM0(IM) 
       enddo 
      enddo 
                                                                        
       call LEGND0 (-U0,PM0,M2_)         ! -cosz , 
       do IM=1,M2_ 
         PM0(IM) = 0.25d0*PM0(IM) 
       enddo 
                                                                        
!---BLKSLV now called with all the wavelength arrays (K=1:W_)           
                                                                        
      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND) 
                                                                        
      return 
      END subroutine MIESCT                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N) 
!-----------------------------------------------------------------------
!---Calculates ORDINARY Legendre fns of X (real)                        
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)            
      implicit none 
      integer, intent(in) :: N 
      real*8, intent(in)  :: X 
      real*8, intent(out) :: PL(N) 
      integer I 
      real*8  DEN 
!---Always does PL(2) = P[1]                                            
        PL(1) = 1.d0 
        PL(2) = X 
        do I = 3,N 
         DEN = (I-1) 
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN) 
        enddo 
      return 
      END subroutine LEGND0                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine BLKSLV                                                 &
     &     (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)          
!-----------------------------------------------------------------------
!  Sets up and solves the block tri-diagonal system:                    
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)            
!  This goes back to the old, dumb, fast version 5.3                    
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_MIE.f' 
!--- expect parameters M_, N_ in parm_MIE.f-----------------------------
                                                                        
      integer, intent(in) ::  ND                   !total Tau layers < N_
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),&  !phase fn  
                              FZ(N_,W_),&          !total attenuation of sun light from top layer to layer 
                              ZTAU(N_,W_)   &      !optical depth
     &                       ,PM(M_,M2_),&         !legender coef. at gauss points
                              PM0(M2_)&            !legendar coef at cosz
     &                       ,RFL(W_),&            !surface albedo
                              ZFLUX(W_)            !refected solar flux at surface layer
                        
      real*8, intent(out) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_) 
                                                                        
      real*8, dimension(M_,N_,W_)    ::  A,C,H,   RR 
                                                                        
      real*8, dimension(M_,M_,N_,W_) ::  B,AA,CC,  DD 
      real*8, dimension(W_,M_,M_) ::  F,E 
      real*8  SUMB,SUMBX,SUMT 
      integer I, J, K, L 
                                                                        
                                                                        
      do K=1,W_ 
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K),    &
     &     PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),                  &
     &             A(1,1,K),H(1,1,K),C(1,1,K), ND)                      
      enddo 
                                                                        
!-----------UPPER BOUNDARY L=1                                          
      do K = 1,W_ 
       do J = 1,M_ 
        do I = 1,M_ 
         F(K,I,J) = B(I,J,1,K) 
        enddo 
       enddo 
      enddo 
                                                                        
      call MATINW (F,E) 
                                                                        
      do K = 1,W_ 
       do J = 1,M_ 
        do I = 1,M_ 
         DD(I,J,1,K) = -E(K,I,1)*CC(1,J,1,K)-E(K,I,2)*CC(2,J,1,K)       &
     &                 -E(K,I,3)*CC(3,J,1,K)-E(K,I,4)*CC(4,J,1,K)       
        enddo 
         RR(J,1,K) = E(K,J,1)*H(1,1,K)+E(K,J,2)*H(2,1,K)                &
     &              +E(K,J,3)*H(3,1,K)+E(K,J,4)*H(4,1,K)                
       enddo 
      enddo 
                                                                        
!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1            
      do L = 2,ND-1 
                                                                        
       do K = 1,W_ 
        do J = 1,M_ 
         do I = 1,M_ 
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K) 
         enddo 
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K) 
        enddo 
                                                                        
        do J = 1,M_ 
         do I = 1,M_ 
          F(K,I,J) = B(I,J,L,K) 
         enddo 
        enddo 
       enddo 
                                                                        
       call MATINW (F,E) 
                                                                        
       do K = 1,W_ 
        do J = 1,M_ 
         do I = 1,M_ 
          DD(I,J,L,K) = - E(K,I,J)*C(J,L,K) 
         enddo 
          RR(J,L,K) = E(K,J,1)*H(1,L,K)+E(K,J,2)*H(2,L,K)               &
     &              + E(K,J,3)*H(3,L,K)+E(K,J,4)*H(4,L,K)               
        enddo 
       enddo 
                                                                        
      enddo 
                                                                        
!---------FINAL DEPTH POINT: L=ND                                       
      L = ND 
                                                                        
      do K = 1,W_ 
        do J = 1,M_ 
         do I = 1,M_ 
          B(I,J,L,K) = B(I,J,L,K)                                       &
     &     + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K)      &
     &     + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)      
         enddo 
          H(J,L,K) = H(J,L,K)                                           &
     &     - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K)          &
     &     - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)          
        enddo 
        do J = 1,M_ 
         do I = 1,M_ 
          F(K,I,J) = B(I,J,L,K) 
         enddo 
        enddo 
      enddo 
                                                                        
      call MATINW (F,E) 
                                                                        
       do K = 1,W_ 
        do J = 1,M_ 
         RR(J,L,K) = E(K,J,1)*H(1,L,K)+E(K,J,2)*H(2,L,K)                &
     &              +E(K,J,3)*H(3,L,K)+E(K,J,4)*H(4,L,K)                
        enddo 
       enddo 
!-----------BACK SOLUTION                                               
      do L = ND-1,1,-1 
       do K = 1,W_ 
        do J = 1,M_ 
         RR(J,L,K) = RR(J,L,K)                                          &
     &    + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K)           &
     &    + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)           
        enddo 
       enddo 
      enddo 
                                                                        
!----------mean J & H                                                   
        FJ(:,:) = 0.d0 
      do L = 1,ND,2 
       do K = 1,W_ 
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2)                     &
     &          + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)                     
       enddo 
      enddo 
      do L = 2,ND,2 
       do K = 1,W_ 
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2)       &
     &          + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)       
       enddo 
      enddo 
                                                                        
!---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)     
!---FJBOT = scaled diffuse flux onto surface:                           
!---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)        
!---SUMBX = flux from Lambert reflected I+                              
      do K = 1,W_ 
                                                                        
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2)         &
     &      + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)         
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2)         &
     &      + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)         
       SUMBX = 4.d0*SUMB*RFL(K)/(1.0d0 + RFL(K)) + ZFLUX(K) 
                                                                        
       FJTOP(K) = 4.d0*SUMT 
       FJBOT(K) = 4.d0*SUMB - SUMBX 
                                                                        
      enddo 
                                                                        
      return 
      END subroutine BLKSLV                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0                 &
     &              ,B,CC,AA,A,H,C,  ND)                                
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:    
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)            
!-----------------------------------------------------------------------
      implicit none 
!      include 'parm_MIE.f' 
!--- expect parameters M_, N_ in parm_MIE.f-----------------------------
      integer, intent(in) ::  ND 
      real*8, intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_) 
      real*8, intent(in)  ::  ZFLUX,RFL 
      real*8, intent(in),dimension(N_) :: FZ,ZTAU 
                                                                        
      real*8, intent(out),dimension(M_,M_,N_) ::  B,AA,CC 
      real*8, intent(out),dimension(M_,N_) ::  A,C,H 
                                                                        
      integer I, J, K, L1,L2,LL 
      real*8  SUM0, SUM1, SUM2, SUM3 
      real*8  DELTAU, D1, D2, SURFAC 
!                                                                       
      real*8, dimension(M_,M_) :: S,T,U,V,W 
!---------------------------------------------                          
                                                                        
!---------upper boundary:  2nd-order terms                              
       L1 = 1 
       L2 = 2 
       do I = 1,M_ 
        SUM0 =                                                          &
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)      &
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)      
        SUM2 =                                                          &
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)      &
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)      
        SUM1 =                                                          &
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)      &
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)      
        SUM3 =                                                          &
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)      &
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)      
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2)) 
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2)) 
       enddo 
                                                                        
       do I = 1,M_ 
        do J = 1,I 
         SUM0 =                                                         &
     &   POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)    &
     & + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)    
         SUM2 =                                                         &
     &   POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)    &
     & + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)    
         SUM1 =                                                         &
     &   POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)    &
     & + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)    
         SUM3 =                                                         &
     &   POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)    &
     & + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)    
         S(I,J) = - SUM2*WT(J) 
         S(J,I) = - SUM2*WT(I) 
         T(I,J) = - SUM1*WT(J) 
         T(J,I) = - SUM1*WT(I) 
         V(I,J) = - SUM3*WT(J) 
         V(J,I) = - SUM3*WT(I) 
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J) 
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I) 
        enddo 
       enddo 
                                                                        
       do I = 1,M_ 
         S(I,I)   = S(I,I)   + 1.0d0 
         T(I,I)   = T(I,I)   + 1.0d0 
         V(I,I)   = V(I,I)   + 1.0d0 
         B(I,I,L1)= B(I,I,L1) + 1.0d0 
                                                                        
         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2)         &
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)         
       enddo 
                                                                        
       do I = 1,M_ 
        do J = 1,M_ 
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2)           &
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)           
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)           &
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)           
        enddo 
       enddo 
!-------------upper boundary, 2nd-order, C-matrix is full (CC)          
         DELTAU = ZTAU(L2) - ZTAU(L1) 
         D2 = 0.25d0*DELTAU 
       do I = 1,M_ 
        do J = 1,M_ 
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) 
         CC(I,J,L1) = D2*U(I,J) 
        enddo 
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1) 
         A(I,L1) = 0.0d0 
       enddo 
       do I = 1,M_ 
        D1 = EMU(I)/DELTAU 
        B(I,I,L1)  = B(I,I,L1) + D1 
        CC(I,I,L1) = CC(I,I,L1) - D1 
       enddo 
                                                                        
!------------intermediate points:  can be even or odd, A & C diagonal   
!---mid-layer h-points, Legendre terms 2,4,6,8                          
       do LL=2,ND-1,2 
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1) 
        do I = 1,M_ 
          A(I,LL) = EMU(I)/DELTAU 
          C(I,LL) = -A(I,LL) 
          H(I,LL) = FZ(LL)*(                                            &
     &     POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4)    &
     &   + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))   
        enddo 
        do I = 1,M_ 
         do J=1,I 
          SUM0 =                                                        &
     &     POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4)  &
     &    +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)  
          B(I,J,LL) =  - SUM0*WT(J) 
          B(J,I,LL) =  - SUM0*WT(I) 
         enddo 
        enddo 
        do I = 1,M_ 
          B(I,I,LL) = B(I,I,LL) + 1.0d0 
        enddo 
       enddo 
                                                                        
!---odd-layer j-points, Legendre terms 1,3,5,7                          
       do LL=3,ND-2,2 
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1) 
        do I = 1,M_ 
          A(I,LL) = EMU(I)/DELTAU 
          C(I,LL) = -A(I,LL) 
          H(I,LL) = FZ(LL)*(                                            &
     &     POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3)    &
     &   + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))   
        enddo 
        do I = 1,M_ 
         do J=1,I 
          SUM0 =                                                        &
     &     POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3)  &
     &    +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)  
          B(I,J,LL) =  - SUM0*WT(J) 
          B(J,I,LL) =  - SUM0*WT(I) 
         enddo 
        enddo 
        do I = 1,M_ 
          B(I,I,LL) = B(I,I,LL) + 1.0d0 
        enddo 
       enddo 
                                                                        
!---------lower boundary:  2nd-order terms                              
       L1 = ND 
       L2 = ND-1 
       do I = 1,M_ 
        SUM0 =                                                          &
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)      &
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)      
        SUM2 =                                                          &
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)      &
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)      
        SUM1 =                                                          &
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)      &
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)      
        SUM3 =                                                          &
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)      &
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)      
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2)) 
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2)) 
       enddo 
                                                                        
       do I = 1,M_ 
        do J = 1,I 
         SUM0 =                                                         &
     &    POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)   &
     &  + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)   
         SUM2 =                                                         &
     &    POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)   &
     &  + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)   
         SUM1 =                                                         &
     &    POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)   &
     &  + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)   
         SUM3 =                                                         &
     &    POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)   &
     &  + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)   
         S(I,J) = - SUM2*WT(J) 
         S(J,I) = - SUM2*WT(I) 
         T(I,J) = - SUM1*WT(J) 
         T(J,I) = - SUM1*WT(I) 
         V(I,J) = - SUM3*WT(J) 
         V(J,I) = - SUM3*WT(I) 
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J) 
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I) 
        enddo 
       enddo 
                                                                        
       do I = 1,M_ 
         S(I,I)   = S(I,I)   + 1.0d0 
         T(I,I)   = T(I,I)   + 1.0d0 
         V(I,I)   = V(I,I)   + 1.0d0 
         B(I,I,L1)= B(I,I,L1) + 1.0d0 
                                                                        
         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2)         &
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)         
       enddo 
                                                                        
       do I = 1,M_ 
        do J = 1,M_ 
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2)           &
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)           
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)           &
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)           
        enddo 
       enddo 
                                                                        
!------------lower boundary, 2nd-order, A-matrix is full (AA)           
         DELTAU = ZTAU(L1) - ZTAU(L2) 
         D2 = 0.25d0*DELTAU 
         SURFAC = 4.0d0*RFL/(1.0d0 + RFL) 
       do I = 1,M_ 
          D1 = EMU(I)/DELTAU 
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4)) 
          SUM1 = SURFAC*SUM0 
        do J = 1,M_ 
         AA(I,J,L1) = - D2*U(I,J) 
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J) 
        enddo 
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX 
       enddo 
                                                                        
       do I = 1,M_ 
          D1 = EMU(I)/DELTAU 
        AA(I,I,L1) = AA(I,I,L1) + D1 
        B(I,I,L1)  = B(I,I,L1) + D1 
        C(I,L1) = 0.0d0 
       enddo 
                                                                        
      return 
      END subroutine GEN_ID                                          
                                                                        
                                                                        
!-----------------------------------------------------------------------
      subroutine MATINW (B,A) 
!-----------------------------------------------------------------------
!  invert 4x4 matrix B(4,4) with L-U decomposition, return as A(4,4) (mjp, old)
      implicit none 
!      include 'parm_MIE.f' 
      real*8, intent(in)   ::  B(W_,4,4) 
      real*8, intent(out)  ::  A(W_,4,4) 
                                                                        
      A(:,:,:) = B(:,:,:) 
!---SETUP L AND U                                                       
      A(:,2,1) = A(:,2,1)/A(:,1,1) 
      A(:,2,2) = A(:,2,2)-A(:,2,1)*A(:,1,2) 
      A(:,2,3) = A(:,2,3)-A(:,2,1)*A(:,1,3) 
      A(:,2,4) = A(:,2,4)-A(:,2,1)*A(:,1,4) 
      A(:,3,1) = A(:,3,1)/A(:,1,1) 
      A(:,3,2) = (A(:,3,2)-A(:,3,1)*A(:,1,2))/A(:,2,2) 
      A(:,3,3) = A(:,3,3)-A(:,3,1)*A(:,1,3)-A(:,3,2)*A(:,2,3) 
      A(:,3,4) = A(:,3,4)-A(:,3,1)*A(:,1,4)-A(:,3,2)*A(:,2,4) 
      A(:,4,1) = A(:,4,1)/A(:,1,1) 
      A(:,4,2) = (A(:,4,2)-A(:,4,1)*A(:,1,2))/A(:,2,2) 
      A(:,4,3) = (A(:,4,3)-A(:,4,1)*A(:,1,3)-A(:,4,2)*A(:,2,3))/A(:,3,3) 
      A(:,4,4) =                                                        &
     &    A(:,4,4)-A(:,4,1)*A(:,1,4)-A(:,4,2)*A(:,2,4)-A(:,4,3)*A(:,3,4)
!---INVERT L                                                            
      A(:,4,3) = -A(:,4,3) 
      A(:,4,2) = -A(:,4,2)-A(:,4,3)*A(:,3,2) 
      A(:,4,1) = -A(:,4,1)-A(:,4,2)*A(:,2,1)-A(:,4,3)*A(:,3,1) 
      A(:,3,2) = -A(:,3,2) 
      A(:,3,1) = -A(:,3,1)-A(:,3,2)*A(:,2,1) 
      A(:,2,1) = -A(:,2,1) 
!---INVERT U                                                            
      A(:,4,4) = 1.d0/A(:,4,4) 
      A(:,3,4) = -A(:,3,4)*A(:,4,4)/A(:,3,3) 
      A(:,3,3) = 1.d0/A(:,3,3) 
      A(:,2,4) = -(A(:,2,3)*A(:,3,4)+A(:,2,4)*A(:,4,4))/A(:,2,2) 
      A(:,2,3) = -A(:,2,3)*A(:,3,3)/A(:,2,2) 
      A(:,2,2) = 1.d0/A(:,2,2) 
      A(:,1,4) =                                                        &
     & -(A(:,1,2)*A(:,2,4)+A(:,1,3)*A(:,3,4)+A(:,1,4)*A(:,4,4))/A(:,1,1)
      A(:,1,3) = -(A(:,1,2)*A(:,2,3)+A(:,1,3)*A(:,3,3))/A(:,1,1) 
      A(:,1,2) = -A(:,1,2)*A(:,2,2)/A(:,1,1) 
      A(:,1,1) = 1.d0/A(:,1,1) 
!---MULTIPLY (:,U-INVERSE)*(:,L-INVERSE)                                
      A(:,1,1) =                                                        &
     &   A(:,1,1)+A(:,1,2)*A(:,2,1)+A(:,1,3)*A(:,3,1)+A(:,1,4)*A(:,4,1) 
      A(:,1,2) = A(:,1,2)+A(:,1,3)*A(:,3,2)+A(:,1,4)*A(:,4,2) 
      A(:,1,3) = A(:,1,3)+A(:,1,4)*A(:,4,3) 
      A(:,2,1) = A(:,2,2)*A(:,2,1)+A(:,2,3)*A(:,3,1)+A(:,2,4)*A(:,4,1) 
      A(:,2,2) = A(:,2,2)+A(:,2,3)*A(:,3,2)+A(:,2,4)*A(:,4,2) 
      A(:,2,3) = A(:,2,3)+A(:,2,4)*A(:,4,3) 
      A(:,3,1) = A(:,3,3)*A(:,3,1)+A(:,3,4)*A(:,4,1) 
      A(:,3,2) = A(:,3,3)*A(:,3,2)+A(:,3,4)*A(:,4,2) 
      A(:,3,3) = A(:,3,3)+A(:,3,4)*A(:,4,3) 
      A(:,4,1) = A(:,4,4)*A(:,4,1) 
      A(:,4,2) = A(:,4,4)*A(:,4,2) 
      A(:,4,3) = A(:,4,4)*A(:,4,3) 
                                                                        
      return 
      END  subroutine MATINW                                          
end module MO_FASTJX_MOD

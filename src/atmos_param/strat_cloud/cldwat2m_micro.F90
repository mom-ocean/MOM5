#ifdef CAM_COMPATIBLE_MICROP 
#undef GFDL_COMPATIBLE_MICROP
#else
#define GFDL_COMPATIBLE_MICROP
#endif
#ifdef GFDL_COMPATIBLE_MICROP 
module cldwat2m_micro_mod
#else
module cldwat2m_micro
#endif

!-------------------------------------------------------------------------
! Purpose:
!   CAM Interface for microphysics
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres 
!                                                                (G2010)  
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!-------------------------------------------------------------------------
! modification for sub-columns, HM, (orig 8/11/10)
! This is done using the logical 'sub_column' set to .true. = subcolumns
!-------------------------------------------------------------------------
#ifdef GFDL_COMPATIBLE_MICROP
! Interfaces, diagnostics, used modules, constants modified for use in GFDL
!based models
! Huan Guo and Rick Hemler, 2010 - 2012
!--------------------------------------------------------------------------
#endif

#ifdef GFDL_COMPATIBLE_MICROP
  use constants_mod,             only: grav, rdgas, rvgas, cp_air, hlv, &
                                       hlf, hls, tfreeze
  use strat_cloud_utilities_mod, only: diag_id_type, diag_pt_type,   &
                                       strat_nml_type
  use mpp_mod,                   only: input_nml_file
  use fms_mod,                   only: mpp_pe, file_exist, error_mesg,  &
                                       open_namelist_file, FATAL, &
                                       stdlog, write_version_number, &
                                       check_nml_error, close_file, &
                                       mpp_root_pe
  use simple_pdf_mod,            only: simple_pdf
#else
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o
  use physconst,      only: latvap, latice
  use wv_saturation,  only: cp, polysvp, epsqs, vqsatd_water
  use cam_history,    only: addfld, add_default, phys_decomp, outfld, &
                            fillvalue
  use cam_logfile,    only: iulog
  use phys_control,   only: phys_getopts
  use cldwat2m_macro, only: rhmini
#endif

  implicit none
  private
  public :: ini_micro, mmicro_pcond, gamma, mmicro_end
  save
 
#ifdef GFDL_COMPATIBLE_MICROP
!------------------------------------------------------------------------
!--version number--------------------------------------------------------
 
character(len=128) :: Version = '$Id: cldwat2m_micro.F90,v 20.0 2013/12/13 23:21:51 fms Exp $'
character(len=128) :: Tagname = '$Name: tikal $'

 
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
INTEGER, PARAMETER :: r8 = dp
#endif
!#else

! logical, public :: liu_in = .true.   ! True = Liu et al 2007 Ice nucleation 
!                                      ! False = cooper fixed ice nucleation (MG2008)


!#ifdef GFDL_COMPATIBLE_MICROP
!INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
!INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
!INTEGER, PARAMETER :: r8 = dp
!#endif

!constants remapped
real(r8), private::  g              !gravity
real(r8), private::  r              !Dry air Gas constant
real(r8), private::  rv             !water vapor gas contstant
real(r8), private::  cpp            !specific heat of dry air
real(r8), private::  rhow           !density of liquid water
real(r8), private::  xxlv           ! latent heat of vaporization
real(r8), private::  xlf            !latent heat of freezing
real(r8), private::  xxls           !latent heat of sublimation

#ifndef GFDL_COMPATIBLE_MICROP
real(r8), private:: rhosn  ! bulk density snow
#endif
real(r8), private:: rhoi   ! bulk density ice

real(r8), private:: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
real(r8), private:: ci,di    !ice mass-diameter relation parameters
real(r8), private:: cs,ds    !snow mass-diameter relation parameters
real(r8), private:: cr,dr    !drop mass-diameter relation parameters
real(r8), private:: f1s,f2s  !ventilation param for snow
real(r8), private:: Eii      !collection efficiency aggregation of ice
real(r8), private:: Ecc      !collection efficiency
real(r8), private:: Ecr      !collection efficiency cloud droplets/rain
real(r8), private:: f1r,f2r  !ventilation param for rain
#ifndef GFDL_COMPATIBLE_MICROP
real(r8), private:: DCS      !autoconversion size threshold
#endif
real(r8), private:: qsmall   !min mixing ratio 
real(r8), private:: bimm,aimm !immersion freezing
real(r8), private:: rhosu     !typical 850mn air density
real(r8), private:: mi0       ! new crystal mass
real(r8), private:: rin       ! radius of contact nuclei
real(r8), private:: qcvar     ! 1/relative variance of sub-grid qc
real(r8), private:: pi       ! pi

! Additional constants to help speed up code

real(r8), private:: cons1
real(r8), private:: cons2
real(r8), private:: cons3
real(r8), private:: cons4
real(r8), private:: cons5
real(r8), private:: cons6
real(r8), private:: cons7
real(r8), private:: cons8
real(r8), private:: cons9
real(r8), private:: cons10
real(r8), private:: cons11
real(r8), private:: cons12
real(r8), private:: cons13
real(r8), private:: cons14
real(r8), private:: cons15
real(r8), private:: cons16
real(r8), private:: cons17
real(r8), private:: cons18
real(r8), private:: cons19
real(r8), private:: cons20
real(r8), private:: cons21
real(r8), private:: cons22
real(r8), private:: cons23
real(r8), private:: cons24
real(r8), private:: cons25
real(r8), private:: cons27
real(r8), private:: cons28

real(r8), private:: lammini
real(r8), private:: lammaxi
real(r8), private:: lamminr
real(r8), private:: lammaxr
real(r8), private:: lammins
real(r8), private:: lammaxs

real(r8), private:: tmax_fsnow ! max temperature for transition to 
                               ! convective snow
real(r8), private:: tmin_fsnow ! min temperature for transition to 
                               ! convective snow

real(r8), private:: tt0       ! Freezing temperature

real(r8), private:: csmin,csmax,minrefl,mindbz

!real(r8), private:: Berg_factor = 1._r8


#ifdef GFDL_COMPATIBLE_MICROP

real       :: dcs = 400.e-6_r8    !autoconversion size threshold
real       :: min_diam_ice = 10.e-6_r8    
logical    :: allow_all_cldtop_collection = .false.
logical    :: rho_factor_in_max_vt = .true.
real       :: max_rho_factor_in_vt = 1.0
real       :: lowest_temp_for_sublimation = 180._r8
real       :: rhosn = 100._r8
!--> cjg: modifications incorporated from Huan's code
logical    :: allow_rain_num_evap = .false.
real       :: accretion_scale = 1.0_r8
!<--cjg

!-->cjg: imposed cloud or ice number
logical,  private::             nccons = .false. ! nccons = true to specify constant cloud droplet number
logical,  private::             nicons = .false. ! nicons = true to specify constant cloud ice number
real(r8), private::             ncnst  = 100.e6  ! specified value (m-3) droplet num concentration (in-cloud not grid-mean) 
real(r8), private::             ninst  = 0.1e6   ! specified value (m-3) ice num concentration (in-cloud not grid-mean)
!<--cjg

logical           :: use_qcvar_in_accretion = .false.
real(r8), private :: qcvar_min4accr         = 0.1
real(r8), private :: qcvar_max4accr         = 0.5
real(r8), private :: accretion_scale_max    = 2.0
real(r8), private :: accr_scale

! <---h1g, 2012-06-12
logical           :: liu_in = .false. ! True = Liu et al 2007 Ice nucleation 
                                   ! False = cooper fixed ice nucleation 
                                   !         (MG2008)
namelist / cldwat2m_micro_nml /   &
                 dcs, min_diam_ice,  &
                 allow_all_cldtop_collection, &
                 max_rho_factor_in_vt, &
                rho_factor_in_max_vt, lowest_temp_for_sublimation, &
!                                     lowest_temp_for_sublimation, &
!--> cjg: modifications incorporated from Huan's code
                 allow_rain_num_evap, accretion_scale, &
!<--cjg
                rhosn,               &  ! h1g
       nccons, ncnst,                &  ! cjg
       nicons, ninst,                &  ! cjg
       liu_in,                       &  ! h1g
       use_qcvar_in_accretion,       &  ! h1g
       qcvar_min4accr,               &  ! h1g
       qcvar_max4accr,               &  ! h1g 
       accretion_scale_max              ! h1g


real(r8), private, parameter :: tmelt  = tfreeze
real(r8), private::             rhmini = 0.80_r8  ! Minimum rh for ice 
                                                  ! cloud fraction
real(r8), private, parameter :: epsqs = rdgas/rvgas
!real(r8),parameter           :: d378 = 1. - epsqs      
real(r8), private, parameter :: d378 = 1. - epsqs      
#else
logical, public :: liu_in = .true. ! True = Liu et al 2007 Ice nucleation 
                                   ! False = cooper fixed ice nucleation 
                                   !         (MG2008)
#endif

contains

!=========================================================================

#ifdef GFDL_COMPATIBLE_MICROP
subroutine ini_micro (qcvar_in)

real,  intent(in) :: qcvar_in

#else
subroutine ini_micro
#endif

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize constants for the morrison microphysics
! called from stratiform.F90
! 
! Author: Andrew Gettelman Dec 2005
! 
!-----------------------------------------------------------------------


#ifdef GFDL_COMPATIBLE_MICROP
   INTEGER   :: unit, io, ierr, logunit
#endif

#ifndef GFDL_COMPATIBLE_MICROP
   integer k

   integer l,m, iaer
   real(r8) surften       ! surface tension of water w/respect to air (N/m)
   real(r8) arg

   character(len=16) :: eddy_scheme = ' '
   logical           :: history_microphysics   ! output variables for 
                                               ! microphysics diagnostics 
                                               ! package

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out = eddy_scheme,           &
                     history_microphysics_out = history_microphysics   )

   ! diagnostic precip
   call addfld ('QRAIN   ','kg/kg   ',pver, 'A',  &
                'Diagnostic grid-mean rain mixing ratio', phys_decomp)
   call addfld ('QSNOW   ','kg/kg   ',pver, 'A',  &
                'Diagnostic grid-mean snow mixing ratio', phys_decomp)
   call addfld ('NRAIN   ','m-3     ',pver, 'A',  &
                'Diagnostic grid-mean rain number conc' ,phys_decomp)
   call addfld ('NSNOW   ','m-3     ',pver, 'A',  &
                'Diagnostic grid-mean snow number conc' ,phys_decomp)

   ! size of precip
   call addfld ('RERCLD   ','m      ',pver, 'A',  &
                'Diagnostic effective radius of Liquid Cloud and Rain' , &
                phys_decomp)

   call addfld ('DSNOW   ','m       ',pver, 'A', &
                'Diagnostic grid-mean snow diameter',phys_decomp)

   ! diagnostic radar reflectivity, cloud-averaged 
   call addfld ('REFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity', &
                phys_decomp)
   call addfld ('AREFL  ','DBz  ',pver, 'A',  &
                'Average 94 GHz radar reflectivity',phys_decomp)
   call addfld ('FREFL  ','fraction  ',pver, 'A',  &
                'Fractional occurance of radar reflectivity' ,phys_decomp)
   
   call addfld ('CSRFL  ','DBz  ',pver, 'A',   &
                '94 GHz radar reflectivity (CloudSat thresholds)' ,  &
                phys_decomp)
   call addfld ('ACSRFL  ','DBz  ',pver, 'A',   &
                'Average 94 GHz radar reflectivity (CloudSat thresholds)',&
                phys_decomp)
   call addfld ('FCSRFL  ','fraction  ',pver, 'A',      &
                'Fractional occurance of radar reflectivity &
                &(CloudSat thresholds)' ,phys_decomp)
 
   call addfld ('AREFLZ ','mm^6/m^3 ',pver, 'A',  &
                'Average 94 GHz radar reflectivity',phys_decomp)

   ! Aerosol information
    call addfld ('NCAL    ','#/m3   ',pver, 'A',  &
                 'Number Concentation Activated for Liquid',phys_decomp)
    call addfld ('NCAI    ','#/m3   ',pver, 'A',    &
                 'Number Concentation Activated for Ice',phys_decomp)


   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call addfld ('AQRAIN   ','kg/kg   ',pver, 'A',  &
                'Average rain mixing ratio' ,phys_decomp)
   call addfld ('AQSNOW   ','kg/kg   ',pver, 'A',    &
                'Average snow mixing ratio' ,phys_decomp)
   call addfld ('ANRAIN   ','m-3     ',pver, 'A',   &
                'Average rain number conc' ,phys_decomp)
   call addfld ('ANSNOW   ','m-3     ',pver, 'A',    &
                'Average snow number conc' ,phys_decomp)
   call addfld ('ADRAIN   ','Micron  ',pver, 'A',   &
                'Average rain effective Diameter',phys_decomp)
   call addfld ('ADSNOW   ','Micron  ',pver, 'A',    &
                'Average snow effective Diameter' ,phys_decomp)
   call addfld ('FREQR  ','fraction  ',pver, 'A',  &
                'Fractional occurance of rain' ,phys_decomp)
   call addfld ('FREQS  ','fraction  ',pver, 'A',   &
                'Fractional occurance of snow' ,phys_decomp)

   if ( history_microphysics) then
      call add_default ('AQSNOW   ', 1, ' ')
      call add_default ('FREQR    ', 1, ' ')
      call add_default ('FREQS    ', 1, ' ')
      call add_default ('AQRAIN   ', 1, ' ')
      call add_default ('AQSNOW   ', 1, ' ')
      call add_default ('ANRAIN   ', 1, ' ')
      call add_default ('ANSNOW   ', 1, ' ')
   end if

#endif

!declarations for morrison codes (transforms variable names):

#ifdef GFDL_COMPATIBLE_MICROP

   qcvar = qcvar_in

!obtain constants from constants_mod when in FMS
   g = grav          !gravity
   r = rdgas         !Dry air Gas constant: 
                     !        note units(phys_constants are in J/K/kmol)
   rv= rvgas         !water vapor gas contstant
   cpp=cp_air        !specific heat of dry air
   rhow = 1000.0_r8  !density of liquid water

! latent heats
   xxlv = hlv        !latent heat vaporization
   xlf  = hlf        !latent heat freezing
   xxls = hls        !latent heat of sublimation

!---------------------------------------------------------------
!     process namelist
!---------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cldwat2m_micro_nml, iostat=io)
      ierr = check_nml_error(io,'cldwat2m_micro_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cldwat2m_micro_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cldwat2m__micro_nml')
        enddo
10      call close_file (unit)
      endif
#endif
 
!-----------------------------------------------------------------------
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                      write (logunit, nml=cldwat2m_micro_nml)
#else
   g= gravit         !gravity
   r= rair           !Dry air Gas constant: 
                     !         note units(phys_constants are in J/K/kmol)
   rv= rh2o          !water vapor gas contstant
   cpp = cpair       !specific heat of dry air
   rhow = rhoh2o     !density of liquid water
   rhosn = 250._r8   !bulk density snow  (++ ceh)

! latent heats

   xxlv = latvap     ! latent heat vaporization
   xlf = latice      ! latent heat freezing
   xxls = xxlv + xlf ! latent heat of sublimation
#endif


! parameters for snow/rain fraction for convective clouds

   tmax_fsnow = tmelt
   tmin_fsnow = tmelt - 5._r8

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

   rhoi = 500._r8     ! bulk density ice
   rhow = 1000._r8    ! bulk density liquid


! fall speed parameters, V = aD^b
! V is in m/s

! droplets
   ac = 3.e7_r8
   bc = 2._r8

! snow
   as = 11.72_r8
   bs = 0.41_r8

! cloud ice
   ai = 700._r8
   bi = 1._r8

! rain
   ar = 841.99667_r8
   br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

   pi= 3.1415927_r8

! cloud ice mass-diameter relationship

   ci = rhoi*pi/6._r8
   di = 3._r8

! snow mass-diameter relationship

   cs = rhosn*pi/6._r8
   ds = 3._r8

! drop mass-diameter relationship

   cr = rhow*pi/6._r8
   dr = 3._r8

! ventilation parameters for snow
! hall and prupacher

   f1s = 0.86_r8
   f2s = 0.28_r8

! collection efficiency, aggregation of cloud ice and snow

   Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain

   Ecr = 1.0_r8

! ventilation constants for rain

   f1r = 0.78_r8
   f2r = 0.32_r8

! autoconversion size threshold for cloud ice to snow (m)

#ifndef GFDL_COMPATIBLE_MICROP
   Dcs = 400.e-6_r8
#endif

! smallest mixing ratio considered in microphysics

   qsmall = 1.e-14_r8  

! immersion freezing parameters, bigg 1953

   bimm = 100._r8
   aimm = 0.66_r8

! typical air density at 850 mb
   rhosu = 85000._r8/(r * tmelt)
!#ifdef GFDL_COMPATIBLE_MICROP
!        rhosu = 85000._r8/(rdgas * tmelt)
!#else
!        rhosu = 85000._r8/(rair * tmelt)
!#endif

! mass of new crystal due to aerosol freezing and growth (kg)

#ifndef GFDL_COMPATIBLE_MICROP
   mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)
#else
   mi0 = 4._r8/3._r8*pi*rhoi*(min_diam_ice)*(min_diam_ice)*(min_diam_ice)
#endif

! radius of contact nuclei aerosol (m)

   rin = 0.1e-6_r8

! 1 / relative variance of sub-grid cloud water distribution
! see morrison and gettelman, 2007, J. Climate for details

#ifndef GFDL_COMPATIBLE_MICROP
   qcvar = 2._r8
#endif

! freezing temperature
   tt0=tmelt     

#ifndef GFDL_COMPATIBLE_MICROP
   pi=4._r8*atan(1.0_r8)
#endif

!Range of cloudsat reflectivities (dBz) for analytic simulator
   csmin= -30._r8
   csmax= 26._r8
   mindbz = -99._r8
!      minrefl = 10._r8**(mindbz/10._r8)
   minrefl = 1.26e-10_r8

! Define constants to help speed up code (limit calls to gamma function)

   cons1=gamma(1._r8+di)
   cons2=gamma(qcvar+2.47_r8)
   cons3=gamma(qcvar)
   cons4=gamma(1._r8+br)
   cons5=gamma(4._r8+br)
   cons6=gamma(1._r8+ds)
   cons7=gamma(1._r8+bs)     
   cons8=gamma(4._r8+bs)     
   cons9=gamma(qcvar+2._r8)     
   cons10=gamma(qcvar+1._r8)
   cons11=gamma(3._r8+bs)
   cons12=gamma(qcvar+1.15_r8)
   cons13=gamma(5._r8/2._r8+br/2._r8)
   cons14=gamma(5._r8/2._r8+bs/2._r8)
   cons15=gamma(qcvar+bc/3._r8)
   cons16=gamma(1._r8+bi)
   cons17=gamma(4._r8+bi)
   cons18=qcvar**2.47_r8
   cons19=qcvar**2
   cons20=qcvar**1.15_r8
   cons21=qcvar**(bc/3._r8)
   cons22=(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
   cons23=dcs**3
   cons24=dcs**2
   cons25=dcs**bs
   cons27=xxlv**2
   cons28=xxls**2

#ifdef GFDL_COMPATIBLE_MICROP
   lammaxi = 1._r8/min_diam_ice 
   lammaxs = 1._r8/min_diam_ice 
#else
   lammaxi = 1._r8/10.e-6_r8
   lammaxs = 1._r8/10.e-6_r8
#endif
   lammini = 1._r8/(2._r8*dcs)
   lammaxr = 1._r8/20.e-6_r8
   lamminr = 1._r8/500.e-6_r8
   lammins = 1._r8/2000.e-6_r8

   return

end subroutine ini_micro

!=========================================================================
!microphysics routine for each timestep goes here...


#ifdef GFDL_COMPATIBLE_MICROP
subroutine mmicro_pcond (dqa_activation, total_activation,    &
                         tiedtke_macrophysics, sub_column,&
                         j, jdim, pver, pcols, ncol, deltatin, tn, &
                         qn, qc_in, qi_in, nc_in, ni_in, p, pdel, cldn,   &
                         liqcldf, icecldf, delta_cf,                   &
                         D_eros_l4, nerosc4, D_eros_i4, nerosi4, dqcdt, &
                         dqidt, naai, npccnin, rndst, nacon,       &
                         tlat, qvlat, qctend, qitend, nctend, nitend,&
                         prect, preci, rflx, sflx,        &
                         qrout,qsout, lsc_rain_size, lsc_snow_size,   &
                         f_snow_berg, Nml, qa0, gamma_mg, SA_0, SA, &
                         ssat_disposal, n_diag_4d, diag_4d, diag_id, &
                         diag_pt, do_clubb, qcvar_clubb)

logical,  intent(in) :: dqa_activation
logical,  intent(in) :: total_activation
logical,  intent(in) :: tiedtke_macrophysics
logical,  intent(in) :: sub_column ! True = configure for sub-columns 
                                   ! False = use w/o sub-columns (standard)
integer,  intent(in) :: j, jdim
integer,  intent(in) :: pver, pcols
integer,  intent(in) :: ncol
real(r8), intent(in) :: deltatin        ! time step (s)
real(r8), intent(in) :: tn(pcols,pver)  ! input temperature (K)
real(r8), intent(in) :: qn(pcols,pver)  ! input h20 vapor mixing ratio 
                                           ! (kg/kg)
! note: all input cloud variables are grid-averaged
real(r8), intent(in) :: qc_in(pcols,pver) ! cloud water mixing ratio 
                                          ! (kg/kg)  
real(r8), intent(in) :: qi_in(pcols,pver) ! cloud ice mixing ratio (kg/kg)
real(r8), intent(in) :: nc_in(pcols,pver) ! grid-average cloud water number
                                          ! conc (#/kg)
real(r8), intent(in) :: ni_in(pcols,pver) ! grid-average cloud ice number 
                                          ! conc (#/kg)
real(r8), intent(in) :: p(pcols,pver)     ! air pressure (pa)
real(r8), intent(in) :: pdel(pcols,pver)  ! pressure difference across 
                                          ! level (pa)
real(r8), intent(inout) :: cldn(pcols,pver)  
                                          ! cloud fraction
real(r8), intent(in) :: liqcldf(pcols,pver)  ! liquid cloud fraction
real(r8), intent(in) :: icecldf(pcols,pver)  ! ice cloud fraction   
real(r8), intent(in) :: delta_cf(pcols,pver) ! 
real(r8), intent(inout) :: D_eros_l4(pcols,pver) ! 
real(r8), intent(in) :: nerosc4 (pcols,pver) ! 
real(r8), intent(inout) :: D_eros_i4(pcols,pver) ! 
real(r8), intent(in) :: nerosi4 (pcols,pver) ! 
real(r8), intent(in) :: dqcdt   (pcols,pver) ! 
real(r8), intent(in) :: dqidt   (pcols,pver) ! 
real(r8), intent(in) :: naai(pcols,pver)     ! ice nulceation number 
real(r8), intent(in) :: npccnin(pcols,pver)  ! ccn activated number 
real(r8), intent(in) :: rndst(pcols,pver,4)  ! radius of 4 dust bins for 
                                             ! contact freezing 
real(r8), intent(in) :: nacon(pcols,pver,4)  ! number in 4 dust bins for 
                                             ! contact freezing 
real(r8), intent(out) :: tlat(pcols,pver)    ! latent heating rate (W/kg)
real(r8), intent(out) :: qvlat(pcols,pver)   ! microphysical tendency qv 
                                             ! (1/s)
real(r8), intent(out) :: qctend(pcols,pver)  ! microphysical tendency qc  
                                             ! (1/s)
real(r8), intent(out) :: qitend(pcols,pver)  ! microphysical tendency qi 
                                             ! (1/s)
real(r8), intent(out) :: nctend(pcols,pver)  ! microphysical tendency nc 
                                             ! (1/s)
real(r8), intent(out) :: nitend(pcols,pver)  ! microphysical tendency ni 
                                             ! (1/s)
real(r8), intent(out) :: prect(pcols)        ! surface precip rate (m/s)
real(r8), intent(out) :: preci(pcols)        ! cloud ice/snow precip rate 
                                             !                        (m/s)
real(r8), intent(out) :: rflx(pcols,pver+1)  ! grid-box average rain flux 
                                             ! (kg m^-2 s^-1)
real(r8), intent(out) :: sflx(pcols,pver+1)  ! grid-box average snow flux 
                                             ! (kg m^-2 s^-1)
real(r8), intent(out) :: qrout(pcols,pver)   ! grid-box average rain 
                                             ! mixing ratio (kg/kg)
real(r8), intent(out) :: qsout(pcols,pver)   ! snow mixing ratio (kg/kg)
real(r8), intent(out) :: lsc_rain_size(pcols,pver) ! raindrop effective  
                                                   ! size (diameter) 
                                                   ! (micron)
real(r8), intent(out) :: lsc_snow_size(pcols,pver) ! snow flake effective 
                                                   ! size (diameter)  
                                                   ! (micron)
real(r8), intent(out) :: f_snow_berg  (pcols,pver) ! ratio of bergeron 
                                                   ! production of qi to 
                                                   ! sum of bergeron, 
                                                   ! riming and freezing
type(strat_nml_type), intent(in) :: Nml
real(r8), intent(in) :: qa0(pcols,pver)       ! 
real(r8), intent(in) :: gamma_mg(pcols,pver)       ! 
real(r8), intent(in) :: SA_0 (pcols,pver)       ! 
real(r8), intent(inout) :: SA   (pcols,pver)       ! 
real(r8), intent(out) :: ssat_disposal(pcols,pver) 
                                 ! disposition of supersaturation at end 
                                 ! of step; 0.= no ssat, 1.= liq, 2.=ice)
INTEGER,INTENT(IN) :: n_diag_4d
REAL, dimension( ncol, jdim, pver, 0:n_diag_4d ), INTENT(INOUT) ::  diag_4d
TYPE(diag_id_type),INTENT(IN) :: diag_id
TYPE(diag_pt_type),INTENT(INout) :: diag_pt

! --> h1g, 2012-10-05
   integer,  intent(in), optional :: do_clubb
   real(r8), intent(in), optional :: qcvar_clubb(pcols,pver)
! <-- h1g, 2012-10-05

#else

subroutine mmicro_pcond ( sub_column,       &
   lchnk, ncol, deltatin, tn,               &
   qn, qc, qi,                              &
   nc, ni, p, pdel, cldn,                   &
   liqcldf, icecldf,                        &
   cldo,                                    &
   rate1ord_cw2pr_st,                       &   
   naai, npccnin, rndst,nacon,              &
   tlat, qvlat,        &
   qctend, qitend, nctend, nitend, effc,    &
   effc_fn, effi, prect, preci,             &  
   nevapr, evapsnow,      &
   prain, prodsnow, cmeout, deffi, pgamrad, &
   lamcrad,qsout,dsout, &
   rflx,sflx, qrout,reff_rain,reff_snow,  &
   qcsevap,qisevap,qvres,cmeiout, &
   vtrmc,vtrmi,qcsedten,qisedten, &
   prao,prco,mnuccco,mnuccto,msacwio,psacwso,&
   bergso,bergo,melto,homoo,qcreso,prcio,praio,qireso,&
   mnuccro,pracso,meltsdt,frzrdt,mnuccdo)

   logical,  intent(in) :: sub_column  ! True = configure for sub-columns 
                                       ! False = use w/o sub-columns 
                                       !                        (standard)
   integer,  intent(in) :: lchnk
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: deltatin        ! time step (s)
   real(r8), intent(in) :: tn(pcols,pver)  ! input temperature (K)
   real(r8), intent(in) :: qn(pcols,pver)  ! input h20 vapor mixing ratio 
                                           ! (kg/kg)
   real(r8), intent(inout) :: qc(pcols,pver) ! cloud water mixing ratio
                                             ! (kg/kg)
   real(r8), intent(inout) :: qi(pcols,pver) ! cloud ice mixing ratio 
                                             ! (kg/kg)
   real(r8), intent(inout) :: nc(pcols,pver) ! cloud water number conc 
                                             ! (1/kg)
   real(r8), intent(inout) :: ni(pcols,pver) ! cloud ice number conc 
                                             ! (1/kg)
   real(r8), intent(in) :: p(pcols,pver)     ! air pressure (pa)
   real(r8), intent(in) :: pdel(pcols,pver)  ! pressure difference across 
                                             ! level (pa)
   real(r8), intent(in) :: cldn(pcols,pver)  ! cloud fraction
   real(r8), intent(in) :: icecldf(pcols,pver) ! ice cloud fraction   
   real(r8), intent(in) :: liqcldf(pcols,pver) ! liquid cloud fraction
   real(r8), intent(inout) :: cldo(pcols,pver) ! old cloud fraction
   real(r8), intent(out) :: rate1ord_cw2pr_st(pcols,pver) 
                                ! 1st order rate for direct cw to 
                                ! precip conversion 
   real(r8), intent(in) :: naai(pcols,pver)    ! ice nulceation number 
                                               ! (from microp_aero_ts) 
   real(r8), intent(in) :: npccnin(pcols,pver) ! ccn activated number 
                                               ! (from microp_aero_ts)
   real(r8), intent(in) :: rndst(pcols,pver,4) ! radius of 4 dust bins for
                                               ! contact freezing (from 
                                               ! microp_aero_ts)
   real(r8), intent(in) :: nacon(pcols,pver,4) ! number in 4 dust bins for 
                                               ! contact freezing  (from 
                                               ! microp_aero_ts)
   real(r8), intent(out) :: tlat(pcols,pver)   ! latent heating rate  
                                               !     (W/kg)
   real(r8), intent(out) :: qvlat(pcols,pver)  ! microphysical tendency qv
                                               ! (1/s)
   real(r8), intent(out) :: qctend(pcols,pver) ! microphysical tendency qc 
                                               ! (1/s) 
   real(r8), intent(out) :: qitend(pcols,pver) ! microphysical tendency qi 
                                               ! (1/s)
   real(r8), intent(out) :: nctend(pcols,pver) ! microphysical tendency nc
                                               ! (1/(kg*s))
   real(r8), intent(out) :: nitend(pcols,pver) ! microphysical tendency ni 
                                               ! (1/(kg*s))
   real(r8), intent(out) :: effc(pcols,pver)   ! droplet effective radius 
                                               ! (micron)
   real(r8), intent(out) :: effc_fn(pcols,pver)! droplet effective radius, 
                                               ! assuming nc = 1.e8 kg-1
   real(r8), intent(out) :: effi(pcols,pver)   ! cloud ice effective radius
                                               ! (micron)
   real(r8), intent(out) :: prect(pcols)       ! surface precip rate (m/s)
   real(r8), intent(out) :: preci(pcols)       ! cloud ice/snow precip rate
                                               ! (m/s)
   real(r8), intent(out) :: nevapr(pcols,pver)  ! evaporation rate of rain 
                                                ! + snow
   real(r8), intent(out) :: evapsnow(pcols,pver)! sublimation rate of snow
   real(r8), intent(out) :: prain(pcols,pver)   ! production of rain + snow
   real(r8), intent(out) :: prodsnow(pcols,pver)! production of snow
   real(r8), intent(out) :: cmeout(pcols,pver)  ! evap/sub of cloud
   real(r8), intent(out) :: deffi(pcols,pver)   ! ice effective diameter 
                                                ! for optics (radiation)
   real(r8), intent(out) :: pgamrad(pcols,pver) ! ice gamma parameter for 
                                                ! optics (radiation)
   real(r8), intent(out) :: lamcrad(pcols,pver) ! slope of droplet 
                                                ! distribution for optics 
                                                ! (radiation)
   real(r8), intent(out) :: qsout(pcols,pver)   ! snow mixing ratio (kg/kg)
   real(r8), intent(out) :: dsout(pcols,pver)   ! snow diameter (m)
   real(r8), intent(out) :: rflx(pcols,pver+1)  ! grid-box average rain 
                                                ! flux (kg m^-2 s^-1)
   real(r8), intent(out) :: sflx(pcols,pver+1)  ! grid-box average snow 
                                                ! flux (kg m^-2 s^-1)
   real(r8), intent(out) :: qrout(pcols,pver)   ! grid-box average rain 
                                                ! mixing ratio (kg/kg)
   real(r8), intent(out) :: reff_rain(pcols,pver) ! rain effective radius 
                                                  ! (micron)
   real(r8), intent(out) :: reff_snow(pcols,pver) ! snow effective radius 
                                                  ! (micron)
   real(r8), intent(out) :: qcsevap(pcols,pver) ! cloud water evaporation 
                                                ! due to sedimentation
   real(r8), intent(out) :: qisevap(pcols,pver) ! cloud ice sublimation due
                                                ! to sublimation
   real(r8), intent(out) :: qvres(pcols,pver)   ! residual condensation 
                                                ! term to ensure RH < 100%
   real(r8), intent(out) :: cmeiout(pcols,pver) ! grid-mean cloud ice 
                                                ! sub/dep
   real(r8), intent(out) :: vtrmc(pcols,pver)   ! mass-weighted cloud water
                                                ! fallspeed
   real(r8), intent(out) :: vtrmi(pcols,pver)   ! mass-weighted cloud ice 
                                                ! fallspeed
   real(r8), intent(out) :: qcsedten(pcols,pver)! qc sedimentation tendency
   real(r8), intent(out) :: qisedten(pcols,pver)! qi sedimentation tendency

! microphysical process rates for output (mixing ratio tendencies)
   real(r8), intent(out) :: prao(pcols,pver) ! accretion of cloud by rain 
   real(r8), intent(out) :: prco(pcols,pver) ! autoconversion of cloud to 
                                             ! rain
   real(r8), intent(out) :: mnuccco(pcols,pver) ! mixing rat tend due to 
                                                ! immersion freezing
   real(r8), intent(out) :: mnuccto(pcols,pver) ! mixing ratio tend due to
                                                ! contact freezing
   real(r8), intent(out) :: msacwio(pcols,pver) ! mixing ratio tend due to 
                                                ! H-M splintering
   real(r8), intent(out) :: psacwso(pcols,pver) ! collection of cloud water
                                                ! by snow
   real(r8), intent(out) :: bergso(pcols,pver)  ! bergeron process on snow
   real(r8), intent(out) :: bergo(pcols,pver)   ! bergeron process on 
                                                ! cloud ice
   real(r8), intent(out) :: melto(pcols,pver)   ! melting of cloud ice
   real(r8), intent(out) :: homoo(pcols,pver)   ! homogeneos freezing 
                                                ! cloud water
   real(r8), intent(out) :: qcreso(pcols,pver)  ! residual cloud conden-
                                                ! sation due to removal of 
                                                ! excess supersat
   real(r8), intent(out) :: prcio(pcols,pver)   ! autoconversion of cloud 
                                                ! ice to snow
   real(r8), intent(out) :: praio(pcols,pver)   ! accretion of cloud ice by
                                                ! snow
   real(r8), intent(out) :: qireso(pcols,pver)  ! residual ice deposition 
                                                ! due to removal of excess 
                                                ! supersat
   real(r8), intent(out) :: mnuccro(pcols,pver) ! mixing ratio tendency due
                                                ! to heterogeneous freezing
                                                ! of rain to snow (1/s)
   real(r8), intent(out) :: pracso (pcols,pver) ! mixing ratio tendency due
                                                ! to accretion of rain by 
                                                ! snow (1/s)
   real(r8), intent(out) :: meltsdt(pcols,pver) ! latent heating rate due 
                                                ! to melting of snow (W/kg)
   real(r8), intent(out) :: frzrdt (pcols,pver) ! latent heating rate due 
                                                ! to homogeneous freezing 
                                                ! of rain (W/kg)
   real(r8), intent(out) :: mnuccdo(pcols,pver) ! mass tendency from ice 
                                                ! nucleation
#endif

!Author: Hugh Morrison, Andrew Gettelman, NCAR
! e-mail: morrison@ucar.edu, andrew@ucar.edu

#ifdef GFDL_COMPATIBLE_MICROP
!   these variables are output by NCAR, but not GFDL, so need to be 
!   declared as local for GFDL implementation

   real(r8) :: qc(pcols,pver)      ! cloud water mixing ratio (kg/kg)
   real(r8) :: qi(pcols,pver)      ! cloud ice mixing ratio (kg/kg)
   real(r8) :: nc(pcols,pver)      ! cloud water number conc (1/kg)
   real(r8) :: ni(pcols,pver)      ! cloud ice number conc (1/kg)
   real(r8) :: rate1ord_cw2pr_st(pcols,pver) 
                                   ! 1st order rate for direct cw to 
                                   ! precip conversion used for scavenging
   real(r8) :: effc(pcols,pver)    ! droplet effective radius (micron)
   real(r8) :: effc_fn(pcols,pver) ! droplet effective radius, 
                                   ! assuming nc = 1.e8 kg-1
   real(r8) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
   real(r8) :: nevapr(pcols,pver)  ! evaporation rate of rain + snow
   real(r8) :: evapsnow(pcols,pver)! sublimation rate of snow
   real(r8) :: prain(pcols,pver)   ! production of rain + snow
   real(r8) :: prodsnow(pcols,pver)! production of snow
   real(r8) :: cmeout(pcols,pver)  ! evap/sub of cloud
   real(r8) :: deffi(pcols,pver)   ! ice effective diameter for optics 
                                   ! (radiation)
   real(r8) :: pgamrad(pcols,pver) ! ice gamma parameter for optics 
                                   ! (radiation)
   real(r8) :: lamcrad(pcols,pver) ! slope of droplet distribution for 
                                   ! optics (radiation)
   real(r8) :: dsout(pcols,pver)   ! snow diameter (m)
   real(r8) :: qcsevap(pcols,pver) ! cloud water evaporation due to 
                                   ! sedimentation
   real(r8) :: qisevap(pcols,pver) ! cloud ice sublimation due to 
                                   ! sublimation
   real(r8) :: qvres(pcols,pver)   ! residual condensation term to ensure 
                                   ! RH < 100%
   real(r8) :: cmeiout(pcols,pver) ! grid-mean cloud ice sub/dep
   real(r8) :: vtrmc(pcols,pver)   ! mass-weighted cloud water fallspeed
   real(r8) :: vtrmi(pcols,pver)   ! mass-weighted cloud ice fallspeed
   real(r8) :: qcsedten(pcols,pver)! qc sedimentation tendency
   real(r8) :: qisedten(pcols,pver)! qi sedimentation tendency
   real(r8) :: prao(pcols,pver)    ! accretion of cloud by rain 
   real(r8) :: prco(pcols,pver)    ! autoconversion of cloud to rain
   real(r8) :: mnuccco(pcols,pver) ! mixing rat tend due to immersion 
                                   ! freezing
   real(r8) :: mnuccto(pcols,pver) ! mixing ratio tend due to contact 
                                   ! freezing
   real(r8) :: msacwio(pcols,pver) ! mixing ratio tend due to H-M 
                                   ! splintering
   real(r8) :: psacwso(pcols,pver) ! collection of cloud water by snow
   real(r8) :: bergso(pcols,pver)  ! bergeron process on snow
   real(r8) :: bergo(pcols,pver)   ! bergeron process on cloud ice
   real(r8) :: melto(pcols,pver)   ! melting of cloud ice
   real(r8) :: homoo(pcols,pver)   ! homogeneos freezign cloud water
   real(r8) :: qcreso(pcols,pver)  ! residual cloud condensation due to 
                                   ! removal of excess supersat
   real(r8) :: prcio(pcols,pver)   ! autoconversion of cloud ice to snow
   real(r8) :: praio(pcols,pver)   ! accretion of cloud ice by snow
   real(r8) :: qireso(pcols,pver)  ! residual ice deposition due to removal
                                   ! of excess supersat
   real(r8) :: mnuccro(pcols,pver) ! mixing ratio tendency due to 
                                   ! heterogeneous freezing of rain to 
                                   ! snow (1/s)
   real(r8) :: pracso (pcols,pver) ! mixing ratio tendency due to accretion
                                   ! of rain by snow (1/s)
   real(r8) :: meltsdt(pcols,pver) ! latent heating rate due to melting of 
                                   ! snow  (W/kg)
   real(r8) :: frzrdt (pcols,pver) ! latent heating rate due to homogeneous
                                   ! freezing of rain (W/kg)
   real(r8) :: mnuccdo(pcols,pver) ! mass tendency from ice nucleation

!  these variables are only used in the GFDL implementation

   real(r8) :: cmelo(pcols,pver)   ! liquid condensation           
   real(r8) :: eroslo(pcols,pver)  ! liquid erosion                
   real(r8) :: erosio(pcols,pver)  ! ice erosion                
   real(r8) :: preo(pcols,pver)    ! rain evaporation 
   real(r8) :: prdso(pcols,pver)   ! snow sublimation 
   logical  :: do_berg1
   logical  :: limit_berg = .false.
   real(r8) :: berg_lim = 1.0e-6_r8
   real(r8) :: dum3                ! temporary dummy variable

! droplet number
   real(r8) :: nucclim(pver)
   real(r8) :: nucclimo(pcols,pver)
   real(r8) :: npccno(pcols,pver)
   real(r8) :: nnuccco(pcols,pver)
   real(r8) :: nnuccto(pcols,pver)
   real(r8) :: npsacwso(pcols,pver)
   real(r8) :: nsubco(pcols,pver)
   real(r8) :: nerosco(pcols,pver)
   real(r8) :: nprao(pcols,pver)
   real(r8) :: nprc1o(pcols,pver)
   real(r8) :: nerosc(pcols,pver)
   real(r8) :: D_eros_l(pcols,pver)

! cloud ice number
   real(r8) :: nucclim1i(pver)
   real(r8) :: nucclim1io(pcols,pver)
   real(r8) :: nnuccdo(pcols,pver)
   real(r8) :: nsacwio(pcols,pver)
   real(r8) :: nsubio(pcols,pver)
   real(r8) :: nerosio(pcols,pver)
   real(r8) :: nprcio(pcols,pver)
   real(r8) :: npraio(pcols,pver)
   real(r8) :: nerosi(pcols,pver)
   real(r8) :: D_eros_i(pcols,pver)


   real(r8) :: cmel     (pcols,pver)
   real(r8) :: cmel_orig(pcols,pver)
   real(r8) :: cmei_orig(pcols,pver)
   real(r8) :: berg_orig(pcols,pver)

   real(r8) :: sum_freeze(pcols,pver)
   real(r8) :: sum_freeze2(pcols,pver)
   real(r8) :: sum_rime  (pcols,pver)
   real(r8) :: sum_berg  (pcols,pver)
   real(r8) :: sum_ice_adj(pcols,pver)
   real(r8) :: sum_bergs (pcols,pver)
   real(r8) :: sum_cond  (pcols,pver)
   real(r8) :: sum_splinter(pcols,pver)
   real(r8) :: qldt_sum
   real(r8) :: eslt, esit, rhi, qs_d, tc
   real(r8) :: qs2d(pcols,pver)
   real(r8) :: qtot(pcols,pver)

   logical  :: lflag = .false.

#endif

#ifndef GFDL_COMPATIBLE_MICROP
!  these variables are not used in the GFDL implementation

   real(r8) :: qcld              ! total cloud water
   real(r8) :: lcldn(pcols,pver) ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver) ! fractional coverage of old liquid cloud
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: arg               ! argument of erfc
   real(r8) :: drout(pcols,pver) ! rain diameter (m)
#endif

! local workspace
! all units mks unless otherwise stated

! temporary variables for sub-stepping 
        real(r8) :: t1(pcols,pver)
        real(r8) :: q1(pcols,pver)
        real(r8) :: qc1(pcols,pver)
        real(r8) :: qi1(pcols,pver)
        real(r8) :: nc1(pcols,pver)
        real(r8) :: ni1(pcols,pver)
        real(r8) :: tlat1(pcols,pver)
        real(r8) :: qvlat1(pcols,pver)
        real(r8) :: qctend1(pcols,pver)
        real(r8) :: qitend1(pcols,pver)
        real(r8) :: nctend1(pcols,pver)
        real(r8) :: nitend1(pcols,pver)
        real(r8) :: prect1(pcols)
        real(r8) :: preci1(pcols)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(r8) :: deltat  ! sub-time step (s)
        real(r8) :: omsm    ! number near unity for round-off issues
        real(r8) :: dto2    ! dt/2 (s)
        real(r8) :: mincld  ! minimum allowed cloud fraction
        real(r8) :: q(pcols,pver) ! water vapor mixing ratio (kg/kg)
        real(r8) :: t(pcols,pver) ! temperature (K)
        real(r8) :: rho(pcols,pver) ! air density (kg m-3)
        real(r8) :: dv(pcols,pver)  ! diffusivity of water vapor in air
        real(r8) :: mu(pcols,pver)  ! viscocity of air
        real(r8) :: sc(pcols,pver)  ! schmidt number
        real(r8) :: kap(pcols,pver) ! thermal conductivity of air
        real(r8) :: rhof(pcols,pver) ! air density correction factor for 
                                     ! fallspeed
        real(r8) :: cldmax(pcols,pver) ! precip fraction assuming maximum 
                                       ! overlap
        real(r8) :: cldm(pcols,pver)   ! cloud fraction
        real(r8) :: icldm(pcols,pver)  ! ice cloud fraction
        real(r8) :: lcldm(pcols,pver)  ! liq cloud fraction
        real(r8) :: icwc(pcols)    ! in cloud water content (liquid+ice)
        real(r8) :: calpha(pcols)  ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cbeta(pcols)   ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cbetah(pcols)  ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cgamma(pcols)  ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cgamah(pcols)  ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: rcgama(pcols)  ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cmec1(pcols)   ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cmec2(pcols)   ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cmec3(pcols)   ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: cmec4(pcols)   ! parameter for cond/evap 
                                   ! (Zhang et al. 2003)
        real(r8) :: qtmp           ! dummy qv 
        real(r8) :: dum            ! temporary dummy variable

!       real(r8) :: cme(pcols,pver)  ! total (liquid+ice) cond/evap rate 
                                     ! of cloud
        real(r8) :: cmei(pcols,pver) ! dep/sublimation rate of cloud ice
        real(r8) :: cwml(pcols,pver) ! cloud water mixing ratio
        real(r8) :: cwmi(pcols,pver) ! cloud ice mixing ratio
        real(r8) :: nnuccd(pver)     ! ice nucleation rate from 
                                     ! deposition/cond.-freezing
        real(r8) :: mnuccd(pver)     ! mass tendency from ice nucleation

! for calculation of rate1ord_cw2pr_st:
        real(r8) :: qcsinksum_rate1ord(pver) ! sum over iterations of cw 
                                             ! to precip sink
        real(r8) :: qcsum_rate1ord(pver)     ! sum over iterations of 
                                             ! cloud water       
        real(r8) :: alpha

        real(r8) :: dum1,dum2   !general dummy variables

        real(r8) :: npccn(pver)      ! droplet activation rate
        real(r8) :: qcic(pcols,pver) ! in-cloud cloud liquid mixing ratio
        real(r8) :: qiic(pcols,pver) ! in-cloud cloud ice mixing ratio
        real(r8) :: qniic(pcols,pver)! in-precip snow mixing ratio
        real(r8) :: qric(pcols,pver) ! in-precip rain mixing ratio
        real(r8) :: ncic(pcols,pver) ! in-cloud droplet number conc
        real(r8) :: niic(pcols,pver) ! in-cloud cloud ice number conc
        real(r8) :: nsic(pcols,pver) ! in-precip snow number conc
        real(r8) :: nric(pcols,pver) ! in-precip rain number conc
        real(r8) :: lami(pver)       ! slope of cloud ice size distr
        real(r8) :: n0i(pver)        ! intercept of cloud ice size distr
        real(r8) :: lamc(pver)       ! slope of cloud liquid size distr
        real(r8) :: n0c(pver)        ! intercept of cloud liquid size distr
        real(r8) :: lams(pver)       ! slope of snow size distr
        real(r8) :: n0s(pver)        ! intercept of snow size distr
        real(r8) :: lamr(pver)       ! slope of rain size distr
        real(r8) :: n0r(pver)        ! intercept of rain size distr
        real(r8) :: cdist1(pver)     ! size distr parameter to calculate 
                                     ! droplet freezing
        real(r8) :: rercld(pcols,pver) ! effective radius calculation for 
                                       ! rain + cloud (combined size of 
                                       ! precip & cloud drops)
        real(r8) :: arcld(pcols,pver)  ! averaging control flag
        real(r8) :: Actmp              ! area cross section of drops
        real(r8) :: Artmp              ! area cross section of rain

        real(r8) :: pgam(pver)         ! spectral width parameter of 
                                       ! droplet size distr
        real(r8) :: lammax             ! maximum allowed slope of size 
                                       ! distr
        real(r8) :: lammin             ! minimum allowed slope of size 
                                       ! distr
        real(r8) :: nacnt              ! number conc of contact ice nuclei
        real(r8) :: mnuccc(pver)       ! mixing ratio tendency due to 
                                       ! freezing of cloud water
        real(r8) :: nnuccc(pver)       ! number conc tendency due to 
                                       ! freezing of cloud water

        real(r8) :: mnucct(pver)       ! mixing ratio tendency due to 
                                       ! contact freezing of cloud water
        real(r8) :: nnucct(pver)       ! number conc tendency due to 
                                       ! contact freezing of cloud water
        real(r8) :: msacwi(pver)       ! mixing ratio tendency due to HM 
                                       ! ice multiplication
        real(r8) :: nsacwi(pver)       ! number conc tendency due to HM 
                                       ! ice multiplication

        real(r8) :: prc(pver)          ! qc tendency due to autoconversion 
                                       ! of cloud droplets
        real(r8) :: nprc(pver)         ! number conc tendency due to 
                                       ! autoconversion of cloud droplets
        real(r8) :: nprc1(pver)        ! qr tendency due to autoconversion 
                                       ! of cloud droplets
        real(r8) :: nsagg(pver)        ! ns tendency due to 
                                       ! self-aggregation of snow
        real(r8) :: dc0                ! mean size droplet size distr
        real(r8) :: ds0                ! mean size snow size distr 
                                       ! (area weighted)
        real(r8) :: eci                ! collection efficiency for riming 
                                       ! of snow by droplets
        real(r8) :: psacws(pver)       ! mixing rat tendency due to 
                                       ! collection of droplets by snow
        real(r8) :: npsacws(pver)      ! number conc tendency due to 
                                       ! collection of droplets by snow
        real(r8) :: uni                ! number-weighted cloud ice 
                                       ! fallspeed
        real(r8) :: umi                ! mass-weighted cloud ice fallspeed
        real(r8) :: uns(pver)          ! number-weighted snow fallspeed
        real(r8) :: ums(pver)          ! mass-weighted snow fallspeed
        real(r8) :: unr(pver)          ! number-weighted rain fallspeed
        real(r8) :: umr(pver)          ! mass-weighted rain fallspeed
        real(r8) :: unc                ! number-weighted cloud droplet 
                                       ! fallspeed
        real(r8) :: umc                ! mass-weighted cloud droplet 
                                       ! fallspeed
        real(r8) :: pracs(pver)        ! mixing rat tendency due to 
                                       ! collection of rain by snow
        real(r8) :: npracs(pver)       ! number conc tendency due to 
                                       ! collection of rain by snow
        real(r8) :: mnuccr(pver)       ! mixing rat tendency due to 
                                       ! freezing of rain
        real(r8) :: nnuccr(pver)       ! number conc tendency due to 
                                       ! freezing of rain
        real(r8) :: pra(pver)          ! mixing rat tendnency due to 
                                       ! accretion of droplets by rain
        real(r8) :: npra(pver)         ! nc tendnency due to accretion of 
                                       ! droplets by rain
        real(r8) :: nragg(pver)        ! nr tendency due to 
                                       ! self-collection of rain
        real(r8) :: prci(pver)         ! mixing rat tendency due to auto-
                                       ! conversion of cloud ice to snow
        real(r8) :: nprci(pver)        ! number conc tendency due to auto-
                                       ! conversion of cloud ice to snow
        real(r8) :: prai(pver)         ! mixing rat tendency due to 
                                       ! accretion of cloud ice by snow
        real(r8) :: nprai(pver)        ! number conc tendency due to 
                                       ! accretion of cloud ice by snow
        real(r8) :: qvs                ! liquid saturation vapor mixing 
                                       ! ratio
        real(r8) :: qvi                ! ice saturation vapor mixing ratio
        real(r8) :: dqsdt              ! change of sat vapor mixing ratio 
                                       ! with temperature
        real(r8) :: dqsidt             ! change of ice sat vapor mixing 
                                       ! ratio with temperature
        real(r8) :: ab                 ! correction factor for rain evap 
                                       ! to account for latent heat
        real(r8) :: qclr               ! water vapor mixing ratio in clear
                                       ! air
        real(r8) :: abi                ! correction factor for snow 
                                       ! sublimation to account for 
                                       ! latent heat
        real(r8) :: epss               ! 1/ sat relaxation timescale for 
                                       ! snow
        real(r8) :: epsr               ! 1/ sat relaxation timescale for 
                                       ! rain
        real(r8) :: pre(pver)          ! rain mixing rat tendency due to 
                                       ! evaporation
        real(r8) :: prds(pver)         ! snow mixing rat tendency due to 
                                       ! sublimation
        real(r8) :: qce                ! dummy qc for conservation check
        real(r8) :: qie                ! dummy qi for conservation check
        real(r8) :: nce                ! dummy nc for conservation check
        real(r8) :: nie                ! dummy ni for conservation check
        real(r8) :: ratio              ! parameter for conservation check
        real(r8) :: dumc(pcols,pver)   ! dummy in-cloud qc
        real(r8) :: dumnc(pcols,pver)  ! dummy in-cloud nc
        real(r8) :: dumi(pcols,pver)   ! dummy in-cloud qi
        real(r8) :: dumni(pcols,pver)  ! dummy in-cloud ni
        real(r8) :: dums(pcols,pver)   ! dummy in-cloud snow mixing rat
        real(r8) :: dumns(pcols,pver)  ! dummy in-cloud snow number conc
        real(r8) :: dumr(pcols,pver)   ! dummy in-cloud rain mixing rat
        real(r8) :: dumnr(pcols,pver)  ! dummy in-cloud rain number conc

! these are parameters for cloud water and cloud ice sedimentation 
! calculations:
        real(r8) :: fr(pver)
        real(r8) :: fnr(pver)
        real(r8) :: fc(pver)
        real(r8) :: fnc(pver)
        real(r8) :: fi(pver)
        real(r8) :: fni(pver)
        real(r8) :: fs(pver)
        real(r8) :: fns(pver)
        real(r8) :: faloutr(pver)
        real(r8) :: faloutnr(pver)
        real(r8) :: faloutc(pver)
        real(r8) :: faloutnc(pver)
        real(r8) :: falouti(pver)
        real(r8) :: faloutni(pver)
        real(r8) :: falouts(pver)
        real(r8) :: faloutns(pver)
        real(r8) :: faltndr
        real(r8) :: faltndnr
        real(r8) :: faltndc
        real(r8) :: faltndnc
        real(r8) :: faltndi
        real(r8) :: faltndni
        real(r8) :: faltnds
        real(r8) :: faltndns
        real(r8) :: faltndqie
        real(r8) :: faltndqce



        real(r8) :: relhum(pcols,pver) ! relative humidity
        real(r8) :: csigma(pcols)      ! parameter for cond/evap of 
                                       ! cloud water/ice
        real(r8) :: rgvm               ! max fallspeed for all species
        real(r8) :: arn(pcols,pver)    ! air density corrected rain 
                                       ! fallspeed parameter
        real(r8) :: asn(pcols,pver)    ! air density corrected snow 
                                       ! fallspeed parameter
        real(r8) :: acn(pcols,pver)    ! air density corrected cloud 
                                       ! droplet fallspeed parameter
        real(r8) :: ain(pcols,pver)    ! air density corrected cloud ice 
                                       ! fallspeed parameter
        real(r8) :: nsubi(pver)        ! evaporation of cloud ice number
        real(r8) :: nsubc(pver)        ! evaporation of droplet number
        real(r8) :: nsubs(pver)        ! evaporation of snow number
        real(r8) :: nsubr(pver)        ! evaporation of rain number
        real(r8) :: mtime              ! factor to account for droplet 
                                       ! activation timescale
        real(r8) :: dz(pcols,pver)     ! height difference across model 
                                       ! vertical level

        real(r8) :: nfice(pcols,pver)  ! fice variable

! precip flux variables for sub-stepping:
        real(r8) :: rflx1(pcols,pver+1)
        real(r8) :: sflx1(pcols,pver+1)

! returns from function/subroutine calls:
        real(r8) :: tsp(pcols,pver)    ! saturation temp (K)
        real(r8) :: qsp(pcols,pver)    ! saturation mixing ratio (kg/kg)
        real(r8) :: qsphy(pcols,pver)  ! saturation mixing ratio (kg/kg): 
                                       ! hybrid rh
        real(r8) :: qs(pcols)          ! liquid-ice weighted sat mixing 
                                       ! rat (kg/kg)
        real(r8) :: es(pcols)          ! liquid-ice weighted sat vapor 
                                       ! press (pa)
        real(r8) :: esl(pcols,pver)    ! liquid sat vapor pressure (pa)
        real(r8) :: esi(pcols,pver)    ! ice sat vapor pressure (pa)
        real(r8) :: gammas(pcols)      ! parameter for cond/evap of cloud 
                                       ! water

! sum of source/sink terms for diagnostic precip:
        real(r8) :: qnitend(pcols,pver)! snow mixing ratio source/sink term
        real(r8) :: nstend(pcols,pver) ! snow number concentration 
                                       ! source/sink term
        real(r8) :: qrtend(pcols,pver) ! rain mixing ratio source/sink term
        real(r8) :: nrtend(pcols,pver) ! rain number concentration 
                                       ! source/sink term
        real(r8) :: qrtot              ! vertically-integrated rain mixing 
                                       ! rat source/sink term
        real(r8) :: nrtot              ! vertically-integrated rain number 
                                       ! conc source/sink term
        real(r8) :: qstot              ! vertically-integrated snow mixing 
                                       ! rat source/sink term
        real(r8) :: nstot              ! vertically-integrated snow number 
                                       ! conc source/sink term

! new terms for Bergeron process
        real(r8) :: dumnnuc            ! provisional ice nucleation rate 
                                       ! (for calculating bergeron)
        real(r8) :: ninew              ! provisional cloud ice number conc 
                                       ! (for calculating bergeron)
        real(r8) :: qinew              ! provisional cloud ice mixing ratio
                                       ! (for calculating bergeron)
        real(r8) :: qvl                ! liquid sat mixing ratio   
        real(r8) :: epsi               ! 1/ sat relaxation timecale for 
                                       ! cloud ice
        real(r8) :: prd                ! provisional deposition rate of 
                                       ! cloud ice at water sat 
        real(r8) :: berg(pcols,pver)   ! mixing rat tendency due to 
                                       ! bergeron process for cloud ice
        real(r8) :: bergs(pver)        ! mixing rat tendency due to 
                                       ! bergeron process for snow

!bergeron terms
        real(r8) :: bergtsf            ! bergeron timescale to remove all 
                                       ! liquid
        real(r8) :: rhin               !modified RH for vapor deposition

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!
        real(r8) :: nrout(pcols,pver)  ! rain number concentration (1/m3)
        real(r8) :: nsout(pcols,pver)  ! snow number concentration (1/m3)

!averaged rain/snow for history
        real(r8) :: qrout2(pcols,pver)
        real(r8) :: qsout2(pcols,pver)
        real(r8) :: nrout2(pcols,pver)
        real(r8) :: nsout2(pcols,pver)
        real(r8) :: freqs(pcols,pver)
        real(r8) :: freqr(pcols,pver)
        real(r8) :: dumfice
        real(r8) :: drout2(pcols,pver) ! mean rain particle diameter (m)
        real(r8) :: dsout2(pcols,pver) ! mean snow particle diameter (m)

!ice nucleation, droplet activation
        real(r8) :: dum2i(pcols,pver)  ! number conc of ice nuclei 
                                       ! available (1/kg)
        real(r8) :: dum2l(pcols,pver)  ! number conc of CCN (1/kg)
        real(r8) :: ncmax
        real(r8) :: nimax

!output fields for number conc
        real(r8) :: ncai(pcols,pver)   ! output number conc of ice nuclei 
                                       ! available (1/m3)
        real(r8) :: ncal(pcols,pver)   ! output number conc of CCN (1/m3)

! loop array variables
        integer i, k, nstep, n, l
        integer ii, kk, m

! loop variables for sub-step solution
        integer iter, it, ltrue(pcols)

! used in contact freezing via dust particles
        real(r8)  tcnt, viscosity, mfp
        real(r8)  slip1, slip2, slip3, slip4
        real(r8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
        real(r8)  nslip1, nslip2, nslip3, nslip4

! used in ice effective radius
        real(r8)  bbi, cci, ak, iciwc, rvi

! used in Bergeron processe and water vapor deposition
        real(r8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
        real(r8)  cldmw(pcols,pver)

! used in secondary ice production
        real(r8) ni_secp

! variables to check for RH after rain evap
        real(r8) :: esn
        real(r8) :: qsn
        real(r8) :: ttmp


        real(r8) :: refl(pcols,pver)   ! analytic radar reflectivity     
        real(r8) :: rainrt(pcols,pver) ! rain rate for reflectivity 
                                       ! calculation
        real(r8) :: rainrt1(pcols,pver)
        real(r8) :: csrfl(pcols,pver)  ! cloudsat reflectivity 
        real(r8) :: arefl(pcols,pver)  ! average reflectivity will zero 
                                       ! points outside valid range
        real(r8) :: acsrfl(pcols,pver) ! cloudsat average
        real(r8) :: frefl(pcols,pver)
        real(r8) :: fcsrfl(pcols,pver)
        real(r8) :: areflz(pcols,pver) !average reflectivity in z.
        real(r8) :: tmp

        real(r8) dmc,ssmc,dstrn        ! variables for modal scheme.
  
        real(r8), parameter :: cdnl    = 0.e6_r8    ! cloud droplet number limiter


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize  output fields for number conc qand ice nucleation
    ncai(1:ncol,1:pver)=0._r8 
    ncal(1:ncol,1:pver)=0._r8  

!Initialize rain size
    rercld(1:ncol,1:pver)=0._r8
    arcld(1:ncol,1:pver)=0._r8

!initialize radiation output variables
    pgamrad(1:ncol,1:pver)=0._r8 ! liquid gamma parameter for optics 
    lamcrad(1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics
    deffi  (1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics

!initialize water vapor tendency term output
    qcsevap(1:ncol,1:pver)=0._r8 
    qisevap(1:ncol,1:pver)=0._r8 
    qvres  (1:ncol,1:pver)=0._r8 
    cmeiout (1:ncol,1:pver)=0._r8
    vtrmc (1:ncol,1:pver)=0._r8
    vtrmi (1:ncol,1:pver)=0._r8
    qcsedten (1:ncol,1:pver)=0._r8
    qisedten (1:ncol,1:pver)=0._r8    

!initialize arrays which accumulate tendencies across sub-steps
    prao(1:ncol,1:pver)=0._r8 
    prco(1:ncol,1:pver)=0._r8 
    mnuccco(1:ncol,1:pver)=0._r8 
    mnuccto(1:ncol,1:pver)=0._r8 
    msacwio(1:ncol,1:pver)=0._r8 
    psacwso(1:ncol,1:pver)=0._r8 
    bergso(1:ncol,1:pver)=0._r8 
    bergo(1:ncol,1:pver)=0._r8 
    melto(1:ncol,1:pver)=0._r8 
    homoo(1:ncol,1:pver)=0._r8 
    qcreso(1:ncol,1:pver)=0._r8 
    prcio(1:ncol,1:pver)=0._r8 
    praio(1:ncol,1:pver)=0._r8 
    qireso(1:ncol,1:pver)=0._r8 
    mnuccro(1:ncol,1:pver)=0._r8 
    pracso (1:ncol,1:pver)=0._r8 
    meltsdt(1:ncol,1:pver)=0._r8
    frzrdt (1:ncol,1:pver)=0._r8
    mnuccdo(1:ncol,1:pver)=0._r8
#ifdef GFDL_COMPATIBLE_MICROP
    preo(1:ncol,1:pver) =0._r8
    prdso(1:ncol,1:pver)=0._r8
    cmelo(1:ncol,1:pver) =0._r8
    eroslo(1:ncol,1:pver) =0._r8
    erosio(1:ncol,1:pver) =0._r8
!droplet number
    nucclimo(1:ncol,1:pver)   = 0._r8
    npccno(1:ncol,1:pver)     = 0._r8
    nnuccco(1:ncol,1:pver)    = 0._r8
    nnuccto(1:ncol,1:pver)    = 0._r8
    npsacwso(1:ncol,1:pver)   = 0._r8
    nsubco(1:ncol,1:pver)     = 0._r8
    nerosco(1:ncol,1:pver)    = 0._r8
    nprao(1:ncol,1:pver)      = 0._r8
    nprc1o(1:ncol,1:pver)     = 0._r8
!ice number
    nucclim1io(1:ncol,1:pver) = 0._r8
    nnuccdo(1:ncol,1:pver)    = 0._r8
    nsacwio(1:ncol,1:pver)    = 0._r8
    nsubio(1:ncol,1:pver)     = 0._r8
    nerosio(1:ncol,1:pver)    = 0._r8
    nprcio(1:ncol,1:pver)     = 0._r8
    npraio(1:ncol,1:pver)     = 0._r8
#endif



! assign variable deltat for sub-stepping...
        deltat=deltatin

! parameters for scheme

        omsm=0.99999_r8
        dto2=0.5_r8*deltat
        mincld=0.0001_r8

! initialize multi-level fields
        q(1:ncol,1:pver)=qn(1:ncol,1:pver)
        t(1:ncol,1:pver)=tn(1:ncol,1:pver)

#ifdef GFDL_COMPATIBLE_MICROP
        qc(1:ncol,1:pver) = qc_in(1:ncol,1:pver)
        qi(1:ncol,1:pver) = qi_in(1:ncol,1:pver)
        nc(1:ncol,1:pver) = nc_in(1:ncol,1:pver)
        ni(1:ncol,1:pver) = ni_in(1:ncol,1:pver)
        if (PRESENT(do_clubb)) lflag=(do_clubb>0)
#endif

! initialize time-varying parameters

        do k=1,pver
          do i=1,ncol
            rho(i,k) = p(i,k)/(r*t(i,k))
            dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
            mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/(t(i,k) + 120._r8)
            sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
            kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)**1.5_r8/   &
                                                      (t(i,k) + 120._r8) 

! air density adjustment for fallspeed parameters
! includes air density correction factor to the
! power of 0.54 following Heymsfield and Bansemer 2007

            rhof(i,k) = (rhosu/rho(i,k))**0.54_r8
#ifdef GFDL_COMPATIBLE_MICROP
            rhof(i,k) = MIN (rhof(i,k), max_rho_factor_in_vt)
#endif
            arn(i,k) = ar*rhof(i,k)
            asn(i,k) = as*rhof(i,k)
            acn(i,k) = ac*rhof(i,k)
            ain(i,k) = ai*rhof(i,k)
   
!#ifdef GFDL_COMPATIBLE_MICROP
!           if (.not. rho_factor_in_max_vt) rhof(i,k) = 1.0
!           rhof(i,k) = MIn(rhof(i,k), 1.6)
!#endif

! get dz from dp and hydrostatic approx
! keep dz positive (define as layer k-1 - layer k)

            dz(i,k) = pdel(i,k)/(rho(i,k)*g)
          end do
        end do

! initialization -- these variables retain the input fields during 
!                   sub-stepping
        t1(1:ncol,1:pver) = t(1:ncol,1:pver)
        q1(1:ncol,1:pver) = q(1:ncol,1:pver)
        qc1(1:ncol,1:pver) = qc(1:ncol,1:pver)
        qi1(1:ncol,1:pver) = qi(1:ncol,1:pver)
        nc1(1:ncol,1:pver) = nc(1:ncol,1:pver)
        ni1(1:ncol,1:pver) = ni(1:ncol,1:pver)

! initialize tendencies to zero
        tlat1(1:ncol,1:pver)=0._r8
        qvlat1(1:ncol,1:pver)=0._r8
        qctend1(1:ncol,1:pver)=0._r8
        qitend1(1:ncol,1:pver)=0._r8
        nctend1(1:ncol,1:pver)=0._r8
        nitend1(1:ncol,1:pver)=0._r8

! initialize precip output
        qrout(1:ncol,1:pver)=0._r8
        qsout(1:ncol,1:pver)=0._r8
        nrout(1:ncol,1:pver)=0._r8
        nsout(1:ncol,1:pver)=0._r8
        dsout(1:ncol,1:pver)=0._r8

#ifdef GFDL_COMPATIBLE_MICROP
!  initialize bergeron fraction arrays
        sum_freeze(1:ncol,1:pver) = 0._r8
        sum_freeze2(1:ncol,1:pver) = 0._r8
        sum_rime(1:ncol,1:pver)   = 0._r8
        sum_berg(1:ncol,1:pver)  = 0._r8
        sum_ice_adj(1:ncol,1:pver)  = 0._r8
        sum_bergs(1:ncol,1:pver)  = 0._r8
        sum_cond (1:ncol,1:pver)  = 0._r8
        sum_splinter(1:ncol,1:pver)  = 0._r8
#endif

#ifndef GFDL_COMPATIBLE_MICROP
        drout(1:ncol,1:pver)=0._r8
!! initialize as fillvalue to avoid Floating Exceptions
        reff_rain(1:ncol,1:pver)=fillvalue
        reff_snow(1:ncol,1:pver)=fillvalue
#endif

! initialize variables for trop_mozart
        nevapr(1:ncol,1:pver)   = 0._r8
        evapsnow(1:ncol,1:pver) = 0._r8
        prain(1:ncol,1:pver)    = 0._r8
        prodsnow(1:ncol,1:pver) = 0._r8
        cmeout(1:ncol,1:pver)   = 0._r8

! for refl calc
        rainrt1(1:ncol,1:pver) = 0._r8

! initialize precip fraction and output tendencies
        cldmax(1:ncol,1:pver) = mincld

!initialize aerosol number
        dum2l(1:ncol,1:pver)=0._r8
        dum2i(1:ncol,1:pver)=0._r8

! initialize avg precip rate
        prect1(1:ncol)=0._r8
        preci1(1:ncol)=0._r8

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Get humidity and saturation vapor pressures

      do k=1,pver

! find wet bulk temperature and saturation value for provisional t and q 
! without condensation

        call vqsatd_water (t(1,k),p(1,k),es,qs,gammas,ncol) ! use rhw

        do i=1,ncol

          esl(i,k) = polysvp(t(i,k),0)
          esi(i,k) = polysvp(t(i,k),1)

! hm fix, make sure when above freezing that esi=esl, not active yet
          if (t(i,k).gt.tmelt) esi(i,k) = esl(i,k)

          relhum(i,k) = q(i,k)/qs(i)

! get cloud fraction, check for minimum
          cldm(i,k)=max(cldn(i,k),mincld)
          cldmw(i,k)=max(cldn(i,k),mincld)

          icldm(i,k)=max(icecldf(i,k),mincld)
          lcldm(i,k)=max(liqcldf(i,k),mincld)

! subcolumns, set cloud fraction variables to one
! if cloud water or ice is present, if not present
! set to mincld (mincld used instead of zero, to prevent
! possible division by zero errors

          if (sub_column) then
            cldm(i,k)=mincld
            cldmw(i,k)=mincld
            icldm(i,k)=mincld
            lcldm(i,k)=mincld
            if (qc(i,k).ge.qsmall) then
              lcldm(i,k)=1.           
              cldm(i,k)=1.
              cldmw(i,k)=1.
            end if
            if (qi(i,k).ge.qsmall) then             
              cldm(i,k)=1.
              icldm(i,k)=1.
            end if
          end if    ! sub-column

! calculate nfice based on liquid and ice mmr (no rain and snow mmr 
! available yet)

          nfice(i,k)=0._r8
          dumfice=qc(i,k)+qi(i,k)
          if (dumfice.gt.qsmall .and. qi(i,k).gt.qsmall) then
            nfice(i,k)=qi(i,k)/dumfice
          endif

! determine number of activated ice nuclei on this step

#ifdef GFDL_COMPATIBLE_MICROP

          if (t(i,k).lt.tmelt - 5._r8) then

! if aerosols interact with ice set number of activated ice nuclei
            if ( liu_in ) then 
              dum2=naai(i,k)
              dumnnuc = (dum2 - ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
            elseif ( Nml%do_ice_nucl_wpdf ) THEN             
              if (total_activation) then
                dum2 = naai(i,k)
                if (Nml%activate_all_ice_always) then
                   dumnnuc = (dum2 - ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
                else
                  if (delta_cf(i,k) .gt. 0._r8) then
                    dumnnuc = (dum2 - ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
                  else
                    dumnnuc = 0._r8
                  endif
                endif
              else if (dqa_activation) then
                dum2 = naai(i,k)
                dumnnuc = max(delta_cf(i,k), 0._r8)*dum2/deltat
              endif
            else
! default when aerosol field not used for nuclei source
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
              dum2=0.005_r8*exp(0.304_r8*(tmelt -t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
              dum2=min(dum2,208.9e3_r8)/rho(i,k) ! convert from m-3 to kg-1
              dumnnuc = (dum2 - ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
            endif

            dumnnuc = max(dumnnuc, 0._r8)

! get provisional ni and qi after nucleation in order to calculate
! Bergeron process below
            ninew = ni(i,k) + dumnnuc*deltat
!What is proper if test here --  dqa or tiedtke or  ????
! or is this incorrect and qi SHOULD increase here (and so also below where
! mnuccd is defined)
            if ( tiedtke_macrophysics .or. dqa_activation) then
              qinew = qi(i,k)
            else
              qinew = qi(i,k) + dumnnuc*deltat*mi0
            endif
          else   ! T>268
            ninew=ni(i,k)
            qinew=qi(i,k)
          end if ! T>268
#else
          if (t(i,k).lt.tmelt - 5._r8) then

! if aerosols interact with ice set number of activated ice nuclei
            if (liu_in) then 
              dum2=naai(i,k)

            else
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
              dum2=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
              dum2=min(dum2,208.9e3_r8)/rho(i,k) ! convert from m-3 to kg-1
            endif

            dumnnuc = (dum2 - ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
            dumnnuc = max(dumnnuc, 0._r8)

! get provisional ni and qi after nucleation in order to calculate
! Bergeron process below
            ninew = ni(i,k) + dumnnuc*deltat
            qinew = qi(i,k) + dumnnuc*deltat*mi0
          else   !  T>268
            ninew = ni(i,k)
            qinew = qi(i,k)
          end if    !  T>268
#endif

! get in-cloud qi and ni after nucleation
          if (icldm(i,k) .gt. 0._r8) then 
            qiic(i,k) = qinew/icldm(i,k)
            niic(i,k) = ninew/icldm(i,k)
          else
            qiic(i,k) = 0._r8 
            niic(i,k) = 0._r8
          endif

!-->cjg
! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
        niic(i,k)=ninst/rho(i,k)
        end if
!<--cjg

!-------------------------------------------------------------------
!Bergeron process

! initialize bergeron process terms to zero
          cmei(i,k)= 0._r8
          berg(i,k) = 0._r8
          prd       = 0._r8

#ifdef GFDL_COMPATIBLE_MICROP
!  define the large-scale cloud water and ice condensation/evaporation from
!  values passed into routine, and their sum.
          cmel(i,k) = dqcdt(i,k)
          cmei(i,k) = dqidt(i,k)
          dum2 = cmel(i,k) + cmei(i,k)

!  if bergeron process to be active for any non-zero condensate, set flag
!  so indicating.
!  current setting is limit_berg = .false. (controlled by nml)
          IF ( .NOT.  limit_berg ) THEN
            if (dum2 .ge. 0._r8) then
              do_berg1 = .true.
            else
              do_berg1 = .false.
            end if
          ELSE
! GFDL has option to not allow bergeron process when cloud ice is less than
! berg_lim, even if have positive condensate. set flag appropriately. 
            if (dum2 .ge. 0._r8 .and. qinew .gt. berg_lim ) then
              do_berg1 = .true.
            else
              do_berg1 = .false.
            end if
          END If
 
          if (do_berg1 .or. lflag ) THEN
#endif

!  calculate bergeron term.
!  temp must be cold enough
            if (t(i,k).lt.tmelt) then
!  ice must exist
              if (qi(i,k).gt.qsmall) then
                bergtsf = 0._r8 ! bergeron time scale (fraction of 
                                !                               timestep)
                qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8 - epsqs)*esi(i,k))
                qvl = epsqs*esl(i,k)/(p(i,k) - (1._r8 - epsqs)*esl(i,k))

!LIMITS  RSH 8/14/12: probably not needed here since liquid not likely to 
!                     be present at these pressures, but not guaranteed. 
              if( .not. lflag ) then
                qvi = MAX(0._r8, MIN (qvi,1.0_r8))
                qvl = MAX(0._r8, MIN (qvl,1.0_r8))
              endif
                dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                abi = 1._r8 + dqsidt*xxls/cpp

! get ice size distribution parameters
                if (qiic(i,k).ge.qsmall) then
                  lami(k) = (cons1*ci* &
                             niic(i,k)/qiic(i,k))**(1._r8/di)
                  n0i(k) = niic(i,k)*lami(k)

! check for slope
! adjust vars
                  if (lami(k).lt.lammini) then
                    lami(k) = lammini
                    n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
                  else if (lami(k).gt.lammaxi) then
                    lami(k) = lammaxi
                    n0i(k) = lami(k)**(di + 1._r8)*qiic(i,k)/(ci*cons1)
                  end if
                  epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

! liquid must exist  
#ifdef GFDL_COMPATIBLE_MICROP
                  if (qc(i,k) + dqcdt(i,k)*deltat .gt. qsmall) then
#else
                  if (qc(i,k) .gt. qsmall) then
#endif

! calculate Bergeron process
                    prd = epsi*(qvl - qvi)/abi
                  else
                    prd = 0._r8
                  end if

! multiply by cloud fraction
                  prd = prd*min(icldm(i,k), lcldm(i,k))

! transfer of existing cloud liquid to ice
                  berg(i,k) = max(0._r8, prd)
                end if  !end qiic   exists bergeron

                if (berg(i,k).gt.0._r8) then
#ifdef GFDL_COMPATIBLE_MICROP
                 if( lflag ) then
                   bergtsf = max(0._r8, (qc(i,k)/berg(i,k))/deltat) 
                   if (bergtsf.lt.1._r8) berg(i,k) = max(0._r8,   &
                                                         qc(i,k)/deltat)
                 else
                  bergtsf = max(0._r8,   &
                                 ((dqcdt(i,k) + qc(i,k)/deltat)/berg(i,k)))

                  if (bergtsf.lt.1._r8) berg(i,k) =    &
                                   max(0._r8, dqcdt(i,k) + qc(i,k)/deltat)
                 endif
#else
                  bergtsf = max(0._r8, (qc(i,k)/berg(i,k))/deltat) 
                  if (bergtsf.lt.1._r8) berg(i,k) = max(0._r8,   &
                                                         qc(i,k)/deltat)
#endif
                endif
#ifdef GFDL_COMPATIBLE_MICROP
! Marc includes a restriction on berg at T < -40C
                if (t(i,k) < tmelt    - 40._r8 .and.  &
        ( .not. lflag ) ) then
                  berg(i,k) = 0._r8
                endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef GFDL_COMPATIBLE_MICROP
!  As per Marc, the  following is inconsistent with the Tiedtke asumption 
!  of in-cloud RH = 1., so it is excluded in the Tiedtke case.
                if ( .not. tiedtke_macrophysics) then
#endif
                  if (bergtsf.lt.1._r8.or.icldm(i,k).gt.lcldm(i,k)) then
                    if (qiic(i,k).ge.qsmall) then

! first case is for case when liquid water is present, but is completely depleted in time step, i.e., bergrsf > 0 but < 1
                      if (qc(i,k).ge.qsmall) then
                        rhin  = (1.0_r8 + relhum(i,k)) / 2._r8
                        if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
                          prd = epsi*(rhin*qvl-qvi)/abi

! multiply by cloud fraction assuming liquid/ice maximum overlap
                          prd = prd*min(icldm(i,k),lcldm(i,k))

! add to cmei
                          cmei(i,k) = cmei(i,k) + (prd * (1._r8- bergtsf))

                        end if ! rhin 
                      end if ! qc > qsmall

! second case is for pure ice cloud, either no liquid, or icldm > lcldm

                      if (qc(i,k).lt.qsmall.or.icldm(i,k).gt.lcldm(i,k))  &
                                                                    then

! note: for case of no liquid, need to set liquid cloud fraction to zero
! store liquid cloud fraction in 'dum'
                        if (qc(i,k).lt.qsmall) then 
                          dum=0._r8 
                        else
                          dum=lcldm(i,k)
                        end if

! set RH to grid-mean value for pure ice cloud
                        rhin = relhum(i,k)
                        if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
                          prd = epsi*(rhin*qvl-qvi)/abi

! multiply by relevant cloud fraction for pure ice cloud
! assuming maximum overlap of liquid/ice
                          prd = prd*max((icldm(i,k)-dum),0._r8)
                          cmei(i,k) = cmei(i,k) + prd
                        end if ! rhin
                      end if ! qc or icldm > lcldm
                    end if   ! qiic .ge.qsmall
                  end if    ! bergtsf or icldm > lcldm

!  if deposition, it should not reduce grid mean rhi below 1.0
                  if (cmei(i,k) > 0.0_r8 .and.     &
                           (relhum(i,k)*esl(i,k)/esi(i,k)) > 1._r8 ) &
                        cmei(i,k) = min(cmei(i,k),   &
                              (q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat)

#ifdef GFDL_COMPATIBLE_MICROP
                endif  ! (tiedtke_macrophysics)
#endif
              end if            !end ice exists loop qi(i,k).gt.qsmall
            end if  ! t(i,k).lt.tmelt  

#ifdef GFDL_COMPATIBLE_MICROP
          endif ! (do_berg1 .or. do_clubb>0)
#endif

!!!!!  END OF BERGERON CALCULATION

#ifdef GFDL_COMPATIBLE_MICROP
! evaporation should not exceed available water
          if( lflag ) then
             if ((-berg(i,k)).lt.-qc(i,k)/deltat) &
                                 berg(i,k) = max(qc(i,k)/deltat, 0._r8)
          endif 
#else 
          if ((-berg(i,k)).lt.-qc(i,k)/deltat) &
                                 berg(i,k) = max(qc(i,k)/deltat, 0._r8)
#endif

#ifdef GFDL_COMPATIBLE_MICROP
! if Tiedtke scheme, this already supplied in input args -- calculated 
! in nc_cond.F90 for Tiedtke scheme

          if (.not. tiedtke_macrophysics) then
#endif

!  sublimation process...

            if ((relhum(i,k)*esl(i,k)/esi(i,k)).lt.1._r8 .and.   &
                                              qiic(i,k).ge.qsmall ) then
              qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8-epsqs)*esi(i,k))
              qvl = epsqs*esl(i,k)/(p(i,k) - (1._r8-epsqs)*esl(i,k))
              dqsidt =  xxls*qvi/(rv*t(i,k)**2)
              abi = 1._r8 + dqsidt*xxls/cpp

! get ice size distribution parameters
              lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
              n0i(k) = niic(i,k)*lami(k)

! check for slope
! adjust vars
              if (lami(k).lt.lammini) then         
                lami(k) = lammini
                n0i(k) = lami(k)**(di + 1._r8)*qiic(i,k)/(ci*cons1)
              else if (lami(k).gt.lammaxi) then
                lami(k) = lammaxi
                n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
              end if
      
              epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

! modify for ice fraction below
              prd = epsi*(relhum(i,k)*qvl-qvi)/abi*icldm(i,k)
              cmei(i,k) = min(prd, 0._r8)
            endif 

! sublimation should not exceed available ice
            if (cmei(i,k).lt.-qi(i,k)/deltat)  cmei(i,k) = -qi(i,k)/deltat

! sublimation should not increase grid mean rhi above 1.0 
            if (cmei(i,k) < 0.0_r8 .and.   &
                 (relhum(i,k)*esl(i,k)/esi(i,k)) < 1._r8 ) &
                          cmei(i,k) = min(0._r8, max(cmei(i,k),  &
                           (q(i,k) - qs(i)*esi(i,k)/esl(i,k))/abi/deltat))

#ifdef GFDL_COMPATIBLE_MICROP
          endif   ! tiedtke_macrophysics
#endif

          cmei(i,k) = cmei(i,k)*omsm

#ifdef GFDL_COMPATIBLE_MICROP
       if( .not. lflag ) &  
          cmel(i,k) = cmel(i,k)*omsm
#endif 

! calculate ice nucleation dum2i

#ifdef GFDL_COMPATIBLE_MICROP
          if (t(i,k).lt.(tmelt - 5._r8)) then 
            if ( liu_in) then
              dum2i(i,k) = naai(i,k)
            elseif ( Nml%do_ice_nucl_wpdf ) THEN
! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
! ice nucleation rate (dum2) has already been calculated and read in (naai)

! if aerosols interact with ice set number of activated ice nuclei
              if (total_activation) then
                dum2i(i,k) = naai(i,k)
              else if (dqa_activation) then
                dum2i(i,k) = naai(i,k)
              endif
            else
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
              dum2i(i,k)=0.005_r8*exp(0.304_r8*(tmelt    -t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
              dum2i(i,k)=min(dum2i(i,k),208.9e3_r8)/rho(i,k) ! convert from
                                                             ! m-3 to kg-1
            endif
          else
            dum2i(i,k)=0._r8
          end if  ! t(i,k).lt.(tmelt - 5._r8)
#else
          if (t(i,k) .lt. (tmelt - 5._r8)) then 
            if (liu_in) then
! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
! ice nucleation rate (dum2) has already been calculated and read in (naai)

! if aerosols interact with ice set number of activated ice nuclei
              dum2i(i,k) = naai(i,k)
            else
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
              dum2i(i,k)=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
              dum2i(i,k)=min(dum2i(i,k),208.9e3_r8)/rho(i,k) ! convert from
                                                             !  m-3 to kg-1
            endif
          else
            dum2i(i,k)=0._r8
          end if  ! t(i,k).lt.(tmelt - 5._r8)
#endif

        end do ! i loop
      end do ! k loop

#ifndef GFDL_COMPATIBLE_MICROP
       cldo(:ncol,:)=cldn(:ncol,:)
#endif


#ifdef GFDL_COMPATIBLE_MICROP
! code for pdf cloud option -- THIS HAS NOT BEEN TESTED AT ALL !!

!  re-calculate cloud fraction
       IF (Nml%do_pdf_clouds .AND.   &
            (Nml%super_ice_opt .EQ. 1 .OR. Nml%super_ice_opt .EQ. 2)) THEN
         IF (Nml%super_ice_opt .EQ. 1 ) THEN

           DO k=1,pver
             DO i= 1,ncol
               ttmp = t(i,k) 
               IF (ttmp .LT. tmelt - 40._r8 .OR.  (ttmp .LE. tmelt .AND. &
                     qc(i,k) + (cmel(i,k) - berg(i,k))/deltat .LT.  &
                                                     3._r8*Nml%qmin)) THEN 
                 eslt = polysvp(ttmp,1)
               ELSE
                 eslt = polysvp(ttmp,0)
               END IF
               qs_d = p(i,k) - d378*eslt
               qs_d = max(qs_d,eslt)
               qs2d(i,k) = epsqs*eslt/qs_d 
             END DO
           END DO
         END IF

         IF  ( Nml%super_ice_opt .EQ. 2 ) THEN
           DO k=1,pver
             DO i= 1,ncol
               ttmp = t(i,k) 
               IF (ttmp .LT. tmelt - 40._r8 .OR. (ttmp .LE. tmelt .AND. &
                      qc(i,k) + ( cmel(i,k) - berg(i,k))/deltat .LT.    &
                                                    3._r8*Nml%qmin) ) THEN 
                 eslt = polysvp(ttmp,1)
                 tc=ttmp-tmelt   
!!!              rhi=MIN( max_super_ice , 0.000195*tc**2+0.00266*tc+1.005)
                 rhi = 0.000195_r8*tc**2+0.00266_r8*tc+1.005_r8
               ELSE
                 eslt = polysvp(ttmp,0)
                 rhi = 1._r8
               END IF
               qs_d = p(i,k) - d378*eslt
               qs_d = max(qs_d,eslt)
               qs2d(i,k)= rhi * epsqs*eslt/qs_d 
             END DO
           END DO
         END IF
         qtot = qn+qc_in+qi_in 

         IF ( Nml%pdf_org )   call error_mesg ( 'cldwat2m_micro', &
                                         'ERROR 1 simple_pdf ', FATAL)
         CALL  simple_pdf(j, ncol, jdim, pver, Nml%qmin, qa0, qtot,    &
                          qs2d,  gamma_mg, Nml%qthalfwidth, Nml%betaP,  &
                          1._r8/deltat, SA_0, n_diag_4d, diag_4d,  &
                          diag_id, diag_pt, SA, cldn)

         do k=1,pver
           do i=1,ncol
             cldm(i,k)=max(cldn  (i,k),mincld)
             lcldm(i,k)=cldm(i,k)
             icldm(i,k)=cldm(i,k)
           end do
         end do 
       END IF  
#endif

!! initialize sub-step precip flux variables
!! flux is zero at top interface.
!      do i=1,ncol
!        rflx1(i,1)=0._r8
!        sflx1(i,1)=0._r8
!      end do 
       do k=1,pver+1
         do i=1,ncol
           rflx1(i,k)=0._r8
           sflx1(i,k)=0._r8
         end do 
       end do 

!! initialize final precip flux variables.
!! flux is zero at top interface, so these should stay as 0.
!      do i=1,ncol
!        rflx(i,1)=0._r8
!        sflx(i,1)=0._r8
       do k=1,pver+1
         do i=1,ncol
           rflx(i,k)=0._r8
           sflx(i,k)=0._r8
         end do 
       end do 

! skip microphysical calculations if no cloud water
       do i=1,ncol
         ltrue(i)=0
         do k=1,pver

#ifndef GFDL_COMPATIBLE_MICROP
           if (qc(i,k).ge.qsmall .or. qi(i,k).ge.qsmall.or.   &
                                       cmei(i,k).ge.qsmall ) ltrue(i) = 1
#endif
#ifdef GFDL_COMPATIBLE_MICROP
           if( lflag ) then
             if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
           else
             if (qc(i,k).ge.qsmall .or. qi(i,k).ge.qsmall .or.  &
                cmei(i,k).ge.qsmall .or. cmel(i,k).ge.qsmall) ltrue(i) = 1
!cms also skip if total water amount is negative anywhere within the column
           if (  qc(i,k) + qi(i,k) + qn(i,k)  .lt. -1.e-9_r8 .OR.   &
                                         qn(i,k)  .lt. -1.e-9_r8 ) then
             ltrue(i)=0
             end if ! (qc(i,k) + qi(i,k) + qn(i,k)  .lt. -1.e-9_r8 .OR. ... )
           end if ! ( lflag )
#endif
         end do
       end do

! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in MG2008

      iter = 2

! get sub-step time step
      deltat = deltat/real(iter)

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid 
! vertical velocity
! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 
! 20 min bulk
! note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8

      mtime=1._r8
      rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01

!!!! skip calculations if no cloud water
!  define output fields if doing so.
      do i=1,ncol
        if (ltrue(i).eq.0) then
          tlat(i,1:pver)=0._r8
          qvlat(i,1:pver)=0._r8
          qctend(i,1:pver)=0._r8
          qitend(i,1:pver)=0._r8
          qnitend(i,1:pver)=0._r8
          qrtend(i,1:pver)=0._r8
          nctend(i,1:pver)=0._r8
          nitend(i,1:pver)=0._r8
          nrtend(i,1:pver)=0._r8
          nstend(i,1:pver)=0._r8
          prect(i)=0._r8
          preci(i)=0._r8
          qniic(i,1:pver)=0._r8
          qric(i,1:pver)=0._r8
          nsic(i,1:pver)=0._r8
          nric(i,1:pver)=0._r8
          rainrt(i,1:pver)=0._r8
!6/6/12
#ifdef GFDL_COMPATIBLE_MICROP
          ssat_disposal(i,1:pver) =0._r8
#endif
          goto 300
        end if

        qcsinksum_rate1ord(1:pver) = 0._r8 
        qcsum_rate1ord(1:pver) = 0._r8 

#ifdef GFDL_COMPATIBLE_MICROP
        cmel_orig(i,1:pver) = cmel(i,1:pver)
        cmei_orig(i,1:pver) = cmei(i,1:pver)
         

        berg_orig(i,1:pver) = berg(i,1:pver)
#endif


!!!!!!!!! BEGIN SUB-STEP LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.......................................................................
        do it=1,iter

! initialize sub-step microphysical tendencies
          tlat(i,1:pver)=0._r8
          qvlat(i,1:pver)=0._r8
          qctend(i,1:pver)=0._r8
          qitend(i,1:pver)=0._r8
          qnitend(i,1:pver)=0._r8
          qrtend(i,1:pver)=0._r8
          nctend(i,1:pver)=0._r8
          nitend(i,1:pver)=0._r8
          nrtend(i,1:pver)=0._r8
          nstend(i,1:pver)=0._r8

! initialize diagnostic precipitation to zero

          qniic(i,1:pver)=0._r8
          qric(i,1:pver)=0._r8
          nsic(i,1:pver)=0._r8
          nric(i,1:pver)=0._r8
   
          rainrt(i,1:pver)=0._r8


! initialize vertically-integrated rain and snow tendencies

          qrtot = 0._r8
          nrtot = 0._r8
          qstot = 0._r8
          nstot = 0._r8

! initialize precip at surface

          prect(i)=0._r8
          preci(i)=0._r8

!  begin new i,k loop. 
          do k=1,pver

! set cwml and cwmi to current qc and qi
            cwml(i,k) = qc(i,k)
            cwmi(i,k) = qi(i,k)

! initialize precip fallspeeds to zero
            ums(k)=0._r8 
            uns(k)=0._r8 
            umr(k)=0._r8 
            unr(k)=0._r8

#ifdef GFDL_COMPATIBLE_MICROP
!  set erosion, bergeron and condensation fields to input values
        if( .not. lflag ) then 
            nerosi(i,k) = nerosi4(i,k)
            nerosc(i,k) = nerosc4(i,k)
            D_eros_l(i,k) = D_eros_l4(i,k)
            D_eros_i(i,k) = D_eros_i4(i,k)
            cmel(i,k) = cmel_orig(i,k)
            cmei(i,k) = cmei_orig(i,k)        
            berg(i,k) = berg_orig(i,k)       
        endif
#endif

!  calculate new cldmax after adjustment to cldm above
! calculate precip fraction based on maximum overlap assumption

! for sub-columns cldm has already been set to 1 if cloud
! water or ice is present, so cldmax will be correctly set below
! and nothing extra needs to be done here

            if (k.eq.1) then
              cldmax(i,k)=cldm(i,k)
            else
! if rain or snow mix ratio is smaller than
! threshold, then set cldmax to cloud fraction at current level
              if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
                cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
              else
                cldmax(i,k)=cldm(i,k)
              end if
            end if

#ifdef GFDL_COMPATIBLE_MICROP
!           nsubc(k) = 0._r8
!           nsubi(k) = 0._r8
            if ( .not. tiedtke_macrophysics) then
! decrease in number concentration due to sublimation/evap
! divide by cloud fraction to get in-cloud decrease
! don't reduce Nc due to bergeron process

              if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and.   &
                                                 cldm(i,k) > mincld) then
                nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
              else
                nsubi(k) = 0._r8
              end if

!!SHOULD NSUBC be nonzero ?? NCAR says no, MG includes code but it is 
!! not activated with Tiedtke. Should it be activated in general for 
!! non-Tiedtke case ??
              if( lflag ) then
                 nsubc(k) = 0._r8  ! (do_clubb >0) 
              else 
              if (cmel(i,k) < 0._r8  .AND. qc(i,k) .ge. qsmall     &
                                     .and. cldm(i,k) > mincld )      then
                nsubc(k) = cmel(i,k)/qc(i,k)*nc(i,k)/cldm(i,k)
              else
                nsubc(k) = 0._r8
              end if ! (cmel(i,k) < 0._r8  .AND. qc(i,k) .ge. qsmall
              endif  ! (do_clubb <=0) 
            else  ! tiedtke_macrophysics
              nsubc(k) = 0._r8
              nsubi(k) = 0._r8
            endif ! tiedtke_macrophysics
#else
! decrease in ice number concentration due to sublimation/evap
! divide by cloud fraction to get in-cloud decrease
! don't reduce Nc due to bergeron process

            if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and.   &
                                                cldm(i,k) > mincld) then
              nsubi(k) = cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
            else
              nsubi(k) = 0._r8
            end if
            nsubc(k) = 0._r8
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  calculate ice nucleation
!  ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
!  note: this is gridbox averaged

#ifndef GFDL_COMPATIBLE_MICROP
 
            if (dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
                 relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini + 0.05_r8) then
 
!if NCAI > 0. then set numice = ncai (as before)
              nnuccd(k) = (dum2i(i,k) - ni(i,k)/icldm(i,k))/deltat*  &
                                                                icldm(i,k)
              nnuccd(k) = max(nnuccd(k), 0._r8)
              nimax = dum2i(i,k)*icldm(i,k)

!Calc mass of new particles using new crystal mass...
!also this will be multiplied by mtime as nnuccd is...
              mnuccd(k) = nnuccd(k)*mi0
 
!  add mnuccd to cmei....
              cmei(i,k) = cmei(i,k) + mnuccd(k)*mtime

!  limit cmei
              qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8-epsqs)*esi(i,k))
              dqsidt =  xxls*qvi/(rv*t(i,k)**2)
              abi = 1._r8 + dqsidt*xxls/cpp
              cmei(i,k)  =min(cmei(i,k), (q(i,k) - qvi)/abi/deltat)

! limit for roundoff error
              cmei(i,k) = cmei(i,k)*omsm

            else
              nnuccd(k) =0._r8
              nimax = 0._r8
              mnuccd(k) = 0._r8
            end if  
#endif

#ifdef GFDL_COMPATIBLE_MICROP
!  Note that ice nucleation calculated later on in code for 
!  Tiedtke macrophysics case.
            if (.not. tiedtke_macrophysics) then
              if (dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
                    relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini+0.05_r8) then
                if( liu_in ) then
                  nnuccd(k) = (dum2i(i,k) - ni(i,k)/icldm(i,k))/deltat*  &
                                                                icldm(i,k)
                  nnuccd(k) = max(nnuccd(k), 0._r8)
                elseif (total_activation) then
                  nnuccd(k) = (dum2i(i,k) - ni(i,k)/icldm(i,k))/deltat*   &
                                                                icldm(i,k)
                  nnuccd(k) = max(nnuccd(k), 0._r8)
                else if (dqa_activation) then
                  nnuccd(k) = max(delta_cf(i,k),0._r8) *dum2i(i,k)/deltatin
                endif
                nimax = dum2i(i,k)*icldm(i,k)

                if (.not. dqa_activation) then

!Calc mass of new particles using new crystal mass...
!also this will be multiplied by mtime as nnuccd is...
                  mnuccd(k) = nnuccd(k) * mi0

!  add mnuccd to cmei....
                  cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime

!  limit cmei
                  qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8-epsqs)*esi(i,k))
                  dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                  abi = 1._r8 + dqsidt*xxls/cpp
                  cmei(i,k) = min(cmei(i,k),(q(i,k) - qvi)/abi/deltat)

! limit for roundoff error
                  cmei(i,k)=cmei(i,k)*omsm
                else
                  mnuccd(k) = 0.
                endif ! (dqa_activation)
              else
                nnuccd(k)=0._r8
                nimax = 0._r8
                mnuccd(k) = 0._r8
              end if  
            endif  ! (tiedtke_macrophyics)
#endif

!c........................................................................
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! obtain in-cloud values of cloud water/ice mixing ratios and number 
! concentrations for microphysical process calculations
! units are kg/kg for mixing ratio, 1/kg for number conc

! limit in-cloud values to 0.005 kg/kg

            qcic(i,k) = min(cwml(i,k)/lcldm(i,k), 5.e-3_r8)
            qiic(i,k) = min(cwmi(i,k)/icldm(i,k), 5.e-3_r8)
            ncic(i,k) = max(nc(i,k)/lcldm(i,k),0._r8)
            niic(i,k) = max(ni(i,k)/icldm(i,k),0._r8)

!-->cjg
! hm add 6/2/11 specify droplet concentration
           if (nccons) then
           ncic(i,k)=ncnst/rho(i,k)
           end if

! hm add 6/2/11 switch for specification of cloud ice number
           if (nicons) then
           niic(i,k)=ninst/rho(i,k)
           end if
!<--cjg

!  adjust previously calculated tendencies to avoid creating negative
!  water species

#ifdef GFDL_COMPATIBLE_MICROP
          if( lflag ) then
            if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
              qcic(i,k)=0._r8
              ncic(i,k)=0._r8
              if (qc(i,k)-berg(i,k)*deltat.lt.0._r8) then
                 berg(i,k)=qc(i,k)/deltat*omsm
              end if
            end if

            if (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat.lt.qsmall) then
              qiic(i,k)=0._r8
              niic(i,k)=0._r8
              if (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat.lt.0._r8) then
                cmei(i,k) = (-qi(i,k)/deltat - berg(i,k))*omsm
              end if
            end if

          else
            if (qc(i,k) + (cmel(i,k) + D_eros_l(i,k) -    &
                                         berg(i,k))*deltat.lt.qsmall) then
              qcic(i,k)=0._r8
              ncic(i,k)=0._r8
              if (qc(i,k) + (cmel(i,k) + D_eros_l(i,k) -   &
                                          berg(i,k))*deltat.lt.0._r8) then
                if (cmel(i,k).lt.0._r8) then
 !++ first only scale cmel, d_eros
                  dum = -cmel(i,k) - D_eros_l(i,k)
                  if (dum .gt. 1.e-30_r8) then
                    dum3 = qc(i,k)/deltat/dum*omsm
                  else
                    dum3 = 0._r8
                  end if
                  cmel(i,k) = dum3*cmel(i,k)
                  D_eros_l(i,k) = dum3*D_eros_l(i,k)
                  dum = -cmel(i,k) - D_eros_l(i,k) + berg(i,k)
                  if (dum .gt. 1.e-30_r8) then
                    dum3 = qc(i,k)/deltat/dum*omsm
                  else
                    dum3 = 0._r8
                  end if
                  cmel(i,k) = dum3*cmel(i,k)
                  D_eros_l(i,k) = dum3*D_eros_l(i,k)
                  berg(i,k) = dum3*berg(i,k)
                else
                  dum = -D_eros_l(i,k) + berg(i,k)
                  if (dum .gt. 1.e-30_r8) then
                    dum3 = ( qc(i,k)/deltat +  cmel(i,k) ) / dum * omsm
                  else
                    dum3 = 0._r8
                  end if
                  D_eros_l(i,k) = D_eros_l(i,k)*dum3
                  berg(i,k) = berg(i,k)*dum3
                endif
              endif
            end if

            if (qi(i,k) + (cmei(i,k) + D_eros_i(i,k) + berg(i,k))*   &
                                                    deltat.lt.qsmall) then
              qiic(i,k)=0._r8
              niic(i,k)=0._r8
              if (qi(i,k) + (cmei(i,k) + D_eros_i(i,k) + berg(i,k))*  &
                                                     deltat.lt.0._r8) then
                if (cmei(i,k).lt.0._r8) then
                  dum = - cmei(i,k) - D_eros_i(i,k)
                  if (dum .gt. 1.e-30_r8) then
                    dum3 = (qi(i,k)/deltat + berg(i,k))/dum*omsm
                  else
                    dum3 = 0._r8
                  end if
                  cmei(i,k) = dum3 * cmei(i,k)
                  D_eros_i(i,k) = dum3 *  D_eros_i(i,k)
                else
                  dum = - D_eros_i(i,k)
                  if (dum .gt. 1.e-30_r8) then
                    dum3 = (qi(i,k)/deltat + cmei(i,k) + berg(i,k))/  &
                                                                 dum*omsm
                  else
                    dum3 = 0._r8
                  end if
                  D_eros_i(i,k) = dum3*D_eros_i(i,k)
                end if
              end if
            end if

          endif
#else
            if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
              qcic(i,k) = 0._r8
              ncic(i,k) = 0._r8
              if (qc(i,k) - berg(i,k)*deltat.lt.0._r8) then
                berg(i,k) = qc(i,k)/deltat*omsm
              end if
            end if

            if (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat.lt.qsmall) then
              qiic(i,k)=0._r8
              niic(i,k)=0._r8
              if (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat.lt.0._r8) then
                cmei(i,k) = (-qi(i,k)/deltat - berg(i,k))*omsm
              end if
            end if
#endif 

! add to cme output
            cmeout(i,k) = cmeout(i,k) + cmei(i,k)

#ifdef GFDL_COMPATIBLE_MICROP
!  calculate ice nuclei activation for tiedtke macrophysics
!  note: this is gridbox averaged
            if ( tiedtke_macrophysics) then
              if (qiic(i,k).ge.qsmall .and. t(i,k).lt.tmelt - 5._r8) then
                if (total_activation) then
                  nnuccd(k) = (dum2i(i,k) - ni(i,k)/icldm(i,k))/deltat*  &
                                                                icldm(i,k)
                  nnuccd(k) = max(nnuccd(k), 0._r8)
                else if (dqa_activation) then
                  nnuccd(k) = max(delta_cf(i,k),0._r8)*dum2i(i,k)/deltatin
                endif
                nimax = dum2i(i,k)*icldm(i,k)
!               if (.not. dqa_activation) then

!Calc mass of new particles using new crystal mass...
!also this will be multiplied by mtime as nnuccd is...
!                 mnuccd(k) = nnuccd(k) * mi0

!  add mnuccd to cmei....
!                 cmei(i,k) = cmei(i,k) + mnuccd(k)*mtime

!  limit cmei
!                 qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8-epsqs)*esi(i,k))
!                 dqsidt =  xxls*qvi/(rv*t(i,k)**2)
!                 abi = 1._r8 + dqsidt*xxls/cpp
!                 cmei(i,k) = min(cmei(i,k), (q(i,k)-qvi)/abi/deltat)

! limit for roundoff error
!                 cmei(i,k) = cmei(i,k)*omsm
!               else
!                 mnuccd(k) = 0.
!               endif   !(dqa_activation)
                mnuccd(k) = 0._r8
              else
                nnuccd(k)=0._r8
                nimax = 0._r8
                mnuccd(k) = 0._r8
              end if 
            endif  ! (tiedtke_macrophysics)
#endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! droplet activation
! calculate potential for droplet activation if cloud water is present
! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
! number (npccnin) is read in from companion routine

! assume aerosols already activated are equal to number of existing droplets for simplicity
! multiply by cloud fraction to obtain grid-average tendency

#ifdef GFDL_COMPATIBLE_MICROP
            if (qcic(i,k).ge.qsmall) then   
              if( lflag ) then
                npccn(k) = max(0._r8, npccnin(i,k))
                dum2l(i,k) = (nc(i,k) + npccn(k)*deltat)/cldm(i,k)
                dum2l(i,k) = max(dum2l(i,k), cdnl/rho(i,k)) ! sghan minimum
                                                       ! in #/cm3  
                ncmax = dum2l(i,k)*cldm(i,k)

              else
              IF ( total_activation) THEN
                dum2l(i,k) = max(0._r8, npccnin(i,k))  
                npccn(k) = ((dum2l(i,k) - nc(i,k)/cldm(i,k))/deltat)* &
                                                                 cldm(i,k)
                npccn(k) = max(0._r8,npccn(k))
                dum2l(i,k) = (nc(i,k) + npccn(k)*deltat)/cldm(i,k)
                dum2l(i,k) = max(dum2l(i,k), cdnl/rho(i,k)) ! sghan minimum
                                                            ! in #/cm3  
              ELSE IF   ( dqa_activation    ) THEN
!delta_cf:  A_dt * (1.-qabar)   where A_dt = A*dt , A source rate
! Eq. 7 of Yi's 2007 paper
!dum2l has already been multiplied by 1.e6/airdens(i,k)
                npccn(k) = max (delta_cf(i,k), 0._r8)*npccnin(i,k)/deltatin
                dum2l(i,k) = (nc(i,k) + npccn(k)*deltat)/cldm(i,k)
              END IF
              ncmax = npccnin(i,k)*cldm(i,k)
              endif
            else
              npccn(k)=0._r8
              ncmax = 0._r8
              dum2l(i,k) = 0._r8
            end if
#else
            if (qcic(i,k).ge.qsmall) then   
              npccn(k) = max(0._r8, npccnin(i,k))
              dum2l(i,k) = (nc(i,k) + npccn(k)*deltat)/cldm(i,k)
              dum2l(i,k) = max(dum2l(i,k),cdnl/rho(i,k)) ! sghan minimum  
                                                       ! in #/cm3  
              ncmax = dum2l(i,k)*cldm(i,k)
            else
              npccn(k)=0._r8
              dum2l(i,k)=0._r8
              ncmax = 0._r8
            end if
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get size distribution parameters based on in-cloud cloud water/ice 
! the calculations also ensure consistency between number and mixing ratio
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!......................................................................
! cloud ice

            if (qiic(i,k).ge.qsmall) then

! impose upper limit on in-cloud number concentration to prevent numerical 
! error
              niic(i,k) = min(niic(i,k), qiic(i,k)*1.e20_r8)
              lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
              n0i(k) = niic(i,k)*lami(k)
! check for slope
! adjust vars
              if (lami(k).lt.lammini) then
                lami(k) = lammini
                n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
                niic(i,k) = n0i(k)/lami(k)
              else if (lami(k).gt.lammaxi) then
                lami(k) = lammaxi
                n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
                niic(i,k) = n0i(k)/lami(k)
              end if
            else
              lami(k) = 0._r8
              n0i(k) = 0._r8
            end if  !qiic(i,k).ge.qsmall 

            if (qcic(i,k).ge.qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical 
! error
              ncic(i,k) = min(ncic(i,k), qcic(i,k)*1.e20_r8)
              ncic(i,k) = max(ncic(i,k), cdnl/rho(i,k)) ! sghan minimum 
                                                        ! in #/cm  

! get pgam from fit to observations of martin et al. 1994
              pgam(k) = 0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k)) +   &
                                                                 0.2714_r8
              pgam(k) = 1._r8/(pgam(k)**2) - 1._r8
              pgam(k) = max(pgam(k), 2._r8)
              pgam(k) = min(pgam(k), 15._r8)

! calculate lamc
              lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k) + 4._r8)/ &
                        (qcic(i,k)*gamma(pgam(k) + 1._r8)))**(1._r8/3._r8)

! lammin, 50 micron diameter max mean size
              lammin = (pgam(k) + 1._r8)/50.e-6_r8
              lammax = (pgam(k) + 1._r8)/2.e-6_r8

              if (lamc(k).lt.lammin) then
                lamc(k) = lammin
                ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)*  &
                                 gamma(pgam(k)+1._r8)/ &
                                    (pi*rhow*gamma(pgam(k) + 4._r8))
              else if (lamc(k).gt.lammax) then
                lamc(k) = lammax
                ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                                    gamma(pgam(k)+1._r8)/ &
                                        (pi*rhow*gamma(pgam(k) + 4._r8))
              end if

! parameter to calculate droplet freezing
              cdist1(k) = ncic(i,k)/gamma(pgam(k) + 1._r8) 
            else
              lamc(k) = 0._r8
              cdist1(k) = 0._r8
            end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin micropysical process calculations 
!.................................................................
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

            if (qcic(i,k).ge.1.e-8_r8) then

! nprc is increase in rain number conc due to autoconversion
! nprc1 is decrease in cloud droplet conc due to autoconversion
! assume exponential sub-grid distribution of qc, resulting in additional
! factor related to qcvar below

 ! hm switch for sub-columns, don't include sub-grid qc
              if (sub_column) then
                prc(k) = 1350._r8*qcic(i,k)**2.47_r8* &
                      (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
                nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
                nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
              else
#ifdef GFDL_COMPATIBLE_MICROP
                if( present(qcvar_clubb) ) then
                  prc(k) = gamma(qcvar_clubb(i,k)+2.47_r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k)**2.47_r8)*1350._r8*qcic(i,k)**2.47_r8* &
                               (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
                else
                  prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8* &
                               (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
                endif
#else
                prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8*&
                               (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
#endif
                nprc(k) = prc(k)/cons22
                nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
              end if               ! sub-column switch
            else
              prc(k)=0._r8
              nprc(k)=0._r8
              nprc1(k)=0._r8
            end if
 
!  add autoconversion to precip from above to get provisional rain mixing 
!  ratio and number concentration (qric and nric)
!  0.45 m/s is fallspeed of new rain drop (80 micron diameter)
            dum = 0.45_r8
            dum1 = 0.45_r8

            if (k.eq.1) then
              qric(i,k) = prc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
              nric(i,k) = nprc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            else
              if (qric(i,k-1).ge.qsmall) then
                dum = umr(k-1)
                dum1 = unr(k-1)
              end if

#ifdef GFDL_COMPATIBLE_MICROP
              if (allow_all_cldtop_collection) then
! NCAR allows no autoconversion of rain number if rain/snow falling from 
! above. this assumes that new drizzle drops formed by autoconversion are 
! rapidly collected by the existing rain/snow particles falling from above.
! Marc's code allowed autoconversion to change rain number, so  variable 
! allow_all_cldtop_collection  introduced, which when .true. would turn off
! this effect. By default, it is .false. for GFDL (as in MG), in contrast 
! to what NCAR does (ifndef GFDL_COMPATIBLE_MICROP). 
                if (qric(i,k-1).ge.1.e-9_r8 .or.    &
                                        qniic(i,k-1).ge.1.e-9_r8) then
                  nprc(k)=0._r8
                end if
              endif  !  allow_all_cldtop_collection
#else
! no autoconversion of rain number if rain/snow falling from above
! this assumes that new drizzle drops formed by autoconversion are rapidly
! collected by the existing rain/snow particles from above
              if (qric(i,k-1).ge.1.e-9_r8 .or.     &
                                         qniic(i,k-1).ge.1.e-9_r8) then
                nprc(k)=0._r8
              end if
#endif

              qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                          (rho(i,k)*dz(i,k)*((pra(k-1) + prc(k))*   &
                           lcldm(i,k) + (pre(k-1) - pracs(k-1) -   &
                                  mnuccr(k-1))*cldmax(i,k))))/    &
                                              (dum*rho(i,k)*cldmax(i,k))
              nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                          (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k) +  &
                             (nsubr(k-1) - npracs(k-1) - nnuccr(k-1)   &
                              + nragg(k-1))*cldmax(i,k))))   &
                                         /(dum1*rho(i,k)*cldmax(i,k))
            end if ! k > 1

!.......................................................................
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)
! note: assumes autoconversion timescale of 180 sec

            if (t(i,k).le.tmelt    .and.qiic(i,k).ge.qsmall) then
               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
                            (cons23/lami(k) + 3._r8*cons24/lami(k)**2 + &
                                     6._r8*dcs/lami(k)**3 +    &
                                        6._r8/lami(k)**4)*exp(-lami(k)*dcs)
            else
              prci(k)=0._r8
              nprci(k)=0._r8
            end if

! add autoconversion to flux from level above to get provisional 
! snow mixing ratio and number concentration (qniic and nsic)
            dum = (asn(i,k)*cons25)
            dum1 = (asn(i,k)*cons25)
            if (k.eq.1) then
              qniic(i,k) = prci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
              nsic(i,k) = nprci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            else
              if (qniic(i,k-1).ge.qsmall) then
                dum = ums(k-1)
                dum1 = uns(k-1)
              end if
              qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*  &
                                                      cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*((prci(k) + prai(k-1) +   &
                               psacws(k-1) + bergs(k-1))*icldm(i,k) +  &
                                 (prds(k-1) + pracs(k-1) + mnuccr(k-1))   &
                                        *cldmax(i,k))))&
                                             /(dum*rho(i,k)*cldmax(i,k))
              nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                              (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k) +   &
                               (nsubs(k-1) + nsagg(k-1) + nnuccr(k-1))*  &
                                 cldmax(i,k))))/(dum1*rho(i,k)*cldmax(i,k))
            end if  ! k > 1

! if precip mix ratio is zero so should number concentration
            if (qniic(i,k).lt.qsmall) then
              qniic(i,k)=0._r8
              nsic(i,k)=0._r8
            end if

            if (qric(i,k).lt.qsmall) then
              qric(i,k)=0._r8
              nric(i,k)=0._r8
            end if

! make sure number concentration is a positive number to avoid
! taking root of negative later

            nric(i,k) = max(nric(i,k), 0._r8)
            nsic(i,k) = max(nsic(i,k), 0._r8)

!.......................................................................
! get size distribution parameters for precip
!......................................................................
! rain
            if (qric(i,k).ge.qsmall) then
              lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
              n0r(k) = nric(i,k)*lamr(k)

! check for slope
! adjust vars
              if (lamr(k).lt.lamminr) then
                lamr(k) = lamminr
                n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                nric(i,k) = n0r(k)/lamr(k)
              else if (lamr(k).gt.lammaxr) then
                lamr(k) = lammaxr
                n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                nric(i,k) = n0r(k)/lamr(k)
              end if

! provisional rain number and mass weighted mean fallspeed (m/s)
              unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
              umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),   &
                                                        9.1_r8*rhof(i,k))
            else
              lamr(k) = 0._r8
              n0r(k) = 0._r8
              umr(k) = 0._r8
              unr(k) = 0._r8
            end if  ! qric(i,k).ge.qsmall

!......................................................................
! snow

            if (qniic(i,k).ge.qsmall) then
              lams(k) = (cons6*cs*nsic(i,k)/ &
              qniic(i,k))**(1._r8/ds)
              n0s(k) = nsic(i,k)*lams(k)

! check for slope
! adjust vars
              if (lams(k).lt.lammins) then
                lams(k) = lammins
                n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
                nsic(i,k) = n0s(k)/lams(k)
              else if (lams(k).gt.lammaxs) then
                lams(k) = lammaxs
                n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
                nsic(i,k) = n0s(k)/lams(k)
              end if

! provisional snow number and mass weighted mean fallspeed (m/s)

              ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),   &
                                                         1.2_r8*rhof(i,k))
              uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))
            else
              lams(k) = 0._r8
              n0s(k) = 0._r8
              ums(k) = 0._r8
              uns(k) = 0._r8
            end if  ! qniic(i,k).ge.qsmall

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! heterogeneous freezing of cloud water
            if (qcic(i,k).ge.qsmall .and. t(i,k).lt.tmelt - 4._r8 ) then  

! immersion freezing (Bigg, 1953)
! subcolumns
              if (sub_column) then
                mnuccc(k) = pi*pi/36._r8*rhow* &
                   cdist1(k)*gamma(7._r8+pgam(k))* &
                       bimm*(exp(aimm*(273.15_r8 - t(i,k))) - 1._r8)/ &
                                                     lamc(k)**3/lamc(k)**3
                nnuccc(k) = &
                    pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8)*bimm* &
                          (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
              else
! ---> h1g, this is the MG 2011-02 version 
              if( present(qcvar_clubb) ) then
                mnuccc(k) = gamma(qcvar_clubb(i,k)+2._r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k)**2)* &
                            pi*pi/36._r8*rhow* &
                            cdist1(k)*gamma(7._r8+pgam(k))* &
                            bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                            lamc(k)**3/lamc(k)**3
                nnuccc(k) = gamma(qcvar_clubb(i,k)+1._r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k))* &
                            pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                            *bimm* &
                           (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
              else
                mnuccc(k) = cons9/(cons3*cons19)* &
                       pi*pi/36._r8*rhow*cdist1(k)*gamma(7._r8+pgam(k))* &
                        bimm*(exp(aimm*(tmelt - t(i,k)))-1._r8)/ &
                                                    lamc(k)**3/lamc(k)**3
                nnuccc(k) = cons10/(cons3*qcvar)* &
                         pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8)*bimm* &
                          (exp(aimm*(tmelt - t(i,k))) - 1._r8)/lamc(k)**3
                endif
! <--- h1g, the MG 2011-02 version
! ---> h1g, this is the MG 2010-09 version
               ! mnuccc(k) = cons9/(cons3*cons19)* &
               !            pi*pi/36._r8*rhow* &
               !        cdist1(k)*gamma(7._r8+pgam(k))* &
               !         bimm*exp(aimm*(273.15_r8-t(i,k)))/ &
               !         lamc(k)**3/lamc(k)**3
               !
               ! nnuccc(k) = cons10/(cons3*qcvar)* &
               !  pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
               !      *bimm* &
               !       exp(aimm*(273.15_r8-t(i,k)))/lamc(k)**3
! <--- h1g, the MG 2010-09 version
              end if           ! sub-columns

#ifdef GFDL_COMPATIBLE_MICROP
!   contact freezing not currently availablein GFDL. Need to get proper
!   dust input fields.
         if( .not. lflag ) then
              mnucct(k) = 0._r8
              nnucct(k) = 0._r8
         else
! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
! dust size and number in 4 bins are read in from companion routine
           tcnt=(270.16_r8-t(i,k))**1.3_r8
           viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
           mfp=2.0_r8*viscosity/(p(i,k)  &                   ! Mean free path (m)
               *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))           

           nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
           nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
           nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
           nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,4)/mfp))))

           ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*rndst(i,k,1))  ! aerosol diffusivity (m2/s)
           ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*rndst(i,k,2))
           ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*rndst(i,k,3))
           ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*rndst(i,k,4))

           if (sub_column) then
 
               mnucct(k) = &
                        (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                        cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
 
               nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                        cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
 
           else

! ---> h1g, this is the MG 2011-02 version
             if( present(qcvar_clubb) ) then
               mnucct(k) = gamma(qcvar_clubb(i,k)+4._r8/3._r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k)**(4._r8/3._r8))*  &
                       (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                       cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
               nnucct(k) =  gamma(qcvar_clubb(i,k)+1._r8/3._r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k)**(1._r8/3._r8))*  &
                         (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                       cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
             else
               mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
                       (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                       cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
               nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
                         (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                       cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
             endif
! <--- h1g, the MG 2011-02 version


! ---> h1g, this is the MG 2010-09 version
           ! mnucct(k) = cons10/(cons3*qcvar)*  &
           !            (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow*&
           !            cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
           ! nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*&
           !            cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
! <--- h1g, the MG 2010-09 version

           end if      ! sub-column switch
         endif
#endif
#ifndef GFDL_COMPATIBLE_MICROP

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated 
! dust.  dust size and number in 4 bins are read in from companion routine
              tcnt = (tmelt - 3._r8 - t(i,k))**1.3_r8
              viscosity = 1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8 ! Viscosity
                                                                !  (kg/m/s)
              mfp = 2.0_r8*viscosity/(p(i,k)  &      ! Mean free path (m)
                      *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k)))) 
              nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*  &
                  Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
              nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*  &
                                        Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
              nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*  &
                                        Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
              nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*  &
                                         Exp(-(1.1_r8*rndst(i,k,4)/mfp))))

              ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*  &
                                rndst(i,k,1))  ! aerosol diffusivity (m2/s)
              ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*  &
                                                       rndst(i,k,2))
              ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*  &
                                                             rndst(i,k,3))
              ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*  &
                                                              rndst(i,k,4))

              if (sub_column) then
                mnucct(k) = &
                       (ndfaer1*(nacon(i,k,1)*tcnt)+  &
                         ndfaer2*(nacon(i,k,2)*tcnt)+  &
                         ndfaer3*(nacon(i,k,3)*tcnt)+  &
                         ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                                cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
 
                nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+   &
                              ndfaer2*(nacon(i,k,2)*tcnt)+   &
                              ndfaer3*(nacon(i,k,3)*tcnt)+   &
                              ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                                cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
              else

! ---> h1g, this is the MG 2011-02 version
                mnucct(k) = gamma(qcvar+4._r8/3._r8)/    &
                                      (cons3*qcvar**(4._r8/3._r8))*  &
                       (ndfaer1*(nacon(i,k,1)*tcnt)+    &
                        ndfaer2*(nacon(i,k,2)*tcnt)+   &
                        ndfaer3*(nacon(i,k,3)*tcnt)+   &
                        ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                              cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
                nnucct(k) =  gamma(qcvar+1._r8/3._r8)/   &
                                       (cons3*qcvar**(1._r8/3._r8))*  &
                       (ndfaer1*(nacon(i,k,1)*tcnt)+   &
                        ndfaer2*(nacon(i,k,2)*tcnt)+   &
                        ndfaer3*(nacon(i,k,3)*tcnt)+   &
                        ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                                cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
! <--- h1g, the MG 2011-02 version


! ---> h1g, this is the MG 2010-09 version
              ! mnucct(k) = cons10/(cons3*qcvar)*  &
              !            (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow*&
              !            cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4
              ! nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*&
              !            cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)
! <--- h1g, the MG 2010-09 version

              end if      ! sub-column switch
#endif

! make sure number of droplets frozen does not exceed available ice nuclei 
! concentration 
! this prevents 'runaway' droplet freezing

              if (nnuccc(k)*lcldm(i,k).gt.nnuccd(k)) then
                dum = (nnuccd(k)/(nnuccc(k)*lcldm(i,k)))

! scale mixing ratio of droplet freezing with limit
                mnuccc(k)=mnuccc(k)*dum
                nnuccc(k)=nnuccd(k)/lcldm(i,k)
              end if
            else
              mnuccc(k) = 0._r8
              nnuccc(k) = 0._r8
              mnucct(k) = 0._r8
              nnucct(k) = 0._r8
            end if  ! qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8

!.......................................................................
! snow self-aggregation from passarelli, 1978, used by reisner, 1998
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

            if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt) then
              nsagg(k) =    &
                -1108._r8*asn(i,k)*Eii*pi**((1._r8-bs)/3._r8)*  &
                rhosn**((-2._r8-bs)/3._r8)*rho(i,k)**((2._r8+bs)/3._r8)*&
                qniic(i,k)**((2._r8+bs)/3._r8)* &
                     (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
                                                 (4._r8*720._r8*rho(i,k))
            else
              nsagg(k) = 0._r8
            end if

!.......................................................................
! accretion of cloud droplets onto snow/graupel
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

! ignore collision of snow with droplets above freezing
 
            if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt .and. &
                                            qcic(i,k).ge.qsmall) then

! put in size dependent collection efficiency
! mean diameter of snow is area-weighted, since
! accretion is function of crystal geometric area
! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

              dc0 = (pgam(k)+1._r8)/lamc(k)
              ds0 = 1._r8/lams(k)
              dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
              eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

              eci = max(eci,0._r8)
              eci = min(eci,1._r8)


! no impact of sub-grid distribution of qc since psacws
! is linear in qc

              psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
                                   n0s(k)*Eci*cons11/ lams(k)**(bs+3._r8)
              npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
                                   n0s(k)*Eci*cons11/lams(k)**(bs+3._r8)
            else
              psacws(k) = 0._r8
              npsacws(k) = 0._r8
            end if

#ifdef GFDL_COMPATIBLE_MICROP
            if (Nml%do_hallet_mossop .or. lflag ) then
#endif
! add secondary ice production due to accretion of droplets by snow 
! (Hallet-Mossop process) (from Cotton et al., 1986)
              if ((t(i,k).lt.tmelt - 3._r8) .and.     &
                                      (t(i,k).ge.tmelt - 5._r8)) then
                ni_secp = 3.5e8_r8*(tmelt - 3._r8-t(i,k))/2.0_r8*psacws(k)
                nsacwi(k) = ni_secp
                msacwi(k) = min(ni_secp*mi0, psacws(k))
              else if((t(i,k).lt.tmelt -5._r8) .and.   &
                                      (t(i,k).ge.tmelt - 8._r8)) then
                ni_secp = 3.5e8_r8*(t(i,k) - (tmelt -8._r8))/   &
                                                       3.0_r8*psacws(k)
                nsacwi(k) = ni_secp
                msacwi(k) = min(ni_secp*mi0, psacws(k))
              else
                ni_secp   = 0.0_r8
                nsacwi(k) = 0.0_r8
                msacwi(k) = 0.0_r8
              endif
              psacws(k) = max(0.0_r8, psacws(k) - ni_secp*mi0)
#ifdef GFDL_COMPATIBLE_MICROP
            else
              msacwi(k) = 0.0_r8
              nsacwi(k) = 0.0_r8
            endif ! (do_hallet_mossop or do_clubb>0)
#endif

!.......................................................................
! accretion of rain water by snow
! formula from ikawa and saito, 1991, used by reisner et al., 1998

            if (qric(i,k).ge.1.e-8_r8 .and.   &
                       qniic(i,k).ge.1.e-8_r8 .and. & 
                                        t(i,k).le.tmelt) then

              pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k) -    &
                                                 0.95_r8*ums(k))**2 + &
                          0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
                          n0r(k)*n0s(k)* &
                               (5._r8/(lamr(k)**6*lams(k)) + &
                                     2._r8/(lamr(k)**5*lams(k)**2) + &
                                        0.5_r8/(lamr(k)**4*lams(k)**3)))

              npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*   &
                                   (unr(k) - uns(k))**2 + &
                             0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                                  (1._r8/(lamr(k)**3*lams(k)) + &
                                     1._r8/(lamr(k)**2*lams(k)**2) + &
                                              1._r8/(lamr(k)*lams(k)**3))

            else
              pracs(k) = 0._r8
              npracs(k) = 0._r8
            end if

!.......................................................................
! heterogeneous freezing of rain drops
! follows from Bigg (1953)

            if (t(i,k).lt.tmelt -4._r8 .and. qric(i,k).ge.qsmall) then

! ---> h1g, this is the MG 2011-02 version
              mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
                 (exp(aimm*(tmelt - t(i,k))) - 1._r8)/lamr(k)**3 &
                                                          /lamr(k)**3
              nnuccr(k) = pi*nric(i,k)*bimm* &
                             (exp(aimm*(tmelt - t(i,k)))-1._r8)/lamr(k)**3
! ---> h1g, the MG 2011-02 version

! ---> h1g, this is the MG 2010-09 version
          !    mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
          !           exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3 &
          !                                                 /lamr(k)**3
          !    nnuccr(k) = pi*nric(i,k)*bimm* &
          !            exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3
! ---> h1g, the MG 2010-09 version

            else
              mnuccr(k) = 0._r8
              nnuccr(k) = 0._r8
            end if

!.......................................................................
! accretion of cloud liquid water by rain
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

            if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then

! include sub-grid distribution of cloud water

! add sub-column switch
 
              if (sub_column) then
                pra(k) = &
                            67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
                npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
              else
!--> cjg: modifications incorporated from Huan's code
#ifdef GFDL_COMPATIBLE_MICROP
                if( present( qcvar_clubb)) then
                  if( use_qcvar_in_accretion ) then
                    if ( qcvar_clubb(i,k) > qcvar_max4accr ) then
                        accr_scale = 1.0
                    elseif( qcvar_clubb(i,k) < qcvar_min4accr ) then
                        accr_scale = accretion_scale_max
                    else
                        accr_scale =   (accretion_scale_max-1.0)/(1.0/qcvar_min4accr - 1.0/qcvar_max4accr) &
                               * (1.0/qcvar_clubb(i,k) - 1.0/qcvar_max4accr) + 1.0
                    endif 
                  else
                    accr_scale = accretion_scale                  
                  endif
                  pra(k) = accr_scale* gamma(qcvar_clubb(i,k)+1.15_r8)/(gamma(qcvar_clubb(i,k))*qcvar_clubb(i,k)**1.15_r8)* &
                      67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
                else
                pra(k) = accretion_scale*cons12/(cons3*cons20)* &
                      67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
                endif
#else
                pra(k) = cons12/(cons3*cons20)* &
                      67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
#endif
!<--cjg
                npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
              end if               ! sub-column switch
            else
              pra(k) = 0._r8
              npra(k) = 0._r8
            end if

!.......................................................................
! Self-collection of rain drops
! from Beheng(1994)

            if (qric(i,k).ge.qsmall) then
              nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
            else
              nragg(k) = 0._r8
            end if

!.......................................................................
! Accretion of cloud ice by snow
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

            if (qniic(i,k).ge.qsmall .and. qiic(i,k).ge.qsmall .and.   &
                                                    t(i,k).le.tmelt) then
              prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
                                 n0s(k)*Eii*cons11/lams(k)**(bs+3._r8)
              nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
                             rho(i,k)*n0s(k)*Eii*cons11/lams(k)**(bs+3._r8)
            else
              prai(k) = 0._r8
              nprai(k) = 0._r8
            end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate evaporation/sublimation of rain and snow
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

! initialize evap/sub tendncies
            pre(k)=0._r8
            prds(k)=0._r8

! evaporation of rain
! only calculate if there is some precip fraction > cloud fraction
!!!NOTE    NOTE   NOTE 
!RSH 7/31/12: The use of 1.e-6 as the limit here results in snow sublim-
!              ation  at high levels (p=1.3 hPa) and a resultant heating
!              term which becomes excessive and causes model blowup 
!              (likely associated with assumed cloud fraction).
!          I have therefore changed the if test to be < qsmall, as is done
!       in all other tests against qcic, qiic, etc in this routine.
!  Marc's code used this qsmall test until 11/27/08 when he changed it --
!   I will also check there to see if it solves problem. 
!  8/7: In MG code, switching to use qsmall also eliminates model 
!       blowups previously encountered. 
!END NOTE

            if (qcic(i,k) + qiic(i,k) .lt. qsmall .or.    &
                                      cldmax(i,k) .gt. lcldm(i,k)) then

! set temporary cloud fraction to zero if cloud water + ice is very small
! this will ensure that evaporation/sublimation of precip occurs over
! entire grid cell, since min cloud fraction is specified otherwise
              if (qcic(i,k) + qiic(i,k) .lt. qsmall) then
                dum = 0._r8
              else
                dum = lcldm(i,k)
              end if

! saturation vapor pressure
              esn = polysvp(t(i,k),0)
              qsn = min(epsqs*esn/(p(i,k) - (1._r8 - epsqs)*esn), 1._r8)
!RSH 8/1/12: Need to prevent negative values which may occur at low p
#ifdef GFDL_COMPATIBLE_MICROP
              if( .not. lflag ) &
              qsn = max(qsn, 0._r8)
#endif
! recalculate saturation vapor pressure for liquid and ice
              esl(i,k) = esn
              esi(i,k) = polysvp(t(i,k),1)
! hm fix, make sure when above freezing that esi=esl, not active yet
              if (t(i,k) .gt. tmelt) esi(i,k) = esl(i,k)

! calculate q for out-of-cloud region
              qclr = (q(i,k) - dum*qsn)/(1._r8 - dum)

              if (qric(i,k).ge.qsmall) then
                qvs = epsqs*esl(i,k)/(p(i,k) - (1._r8 - epsqs)*esl(i,k))
                dqsdt = xxlv*qvs/(rv*t(i,k)**2)
                ab = 1._r8 + dqsdt*xxlv/cpp
                epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
                                                 (f1r/(lamr(k)*lamr(k)) + &
                          f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                                   sc(i,k)**(1._r8/3._r8)*cons13/ &
                                      (lamr(k)**(5._r8/2._r8 + br/2._r8)))
                pre(k) = epsr*(qclr - qvs)/ab

! only evaporate in out-of-cloud region
! and distribute across cldmax
                pre(k) = min(pre(k)*(cldmax(i,k)-dum), 0._r8)
                pre(k) = pre(k)/cldmax(i,k)
              end if

! sublimation of snow
              if (qniic(i,k).ge.qsmall) then
                qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8 - epsqs)*esi(i,k))
                dqsidt =  xxls*qvi/(rv*t(i,k)**2)
                abi = 1._r8 + dqsidt*xxls/cpp
                epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                                           (f1s/(lams(k)*lams(k)) + &
                          f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                                 sc(i,k)**(1._r8/3._r8)*cons14/ &
                                      (lams(k)**(5._r8/2._r8 + bs/2._r8)))
                prds(k) = epss*(qclr - qvi)/abi

! only sublimate in out-of-cloud region and distribute over cldmax
                prds(k) = min(prds(k)*(cldmax(i,k)-dum), 0._r8)
                prds(k) = prds(k)/cldmax(i,k)
              end if

! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
! get updated RH at end of time step based on cloud water/ice condensation/evap

#ifdef GFDL_COMPATIBLE_MICROP
            if( .not. lflag ) then
              qtmp = q(i,k) -     &
                       (D_eros_l(i,k) + D_eros_i(i,k) + cmel(i,k) +   &
                        cmei(i,k) + (pre(k) + prds(k))*cldmax(i,k))*deltat
              ttmp = t(i,k) +(  &
                       (D_eros_l(i,k) + cmel(i,k) + pre(k)*cldmax(i,k))*  &
                                                              xxlv + &
                       (D_eros_i(i,k) + cmei(i,k) + prds(k)*cldmax(i,k))* &
                                                         xxls)*deltat/cpp
            else
              qtmp = q(i,k) -    &
                        (cmei(i,k) + (pre(k) + prds(k))*cldmax(i,k))*deltat
              ttmp = t(i,k) +    &
                      ((pre(k)*cldmax(i,k))*xxlv + &
                        (cmei(i,k) + prds(k)*cldmax(i,k))*xxls)*deltat/cpp
            endif 
#else
              qtmp = q(i,k) -    &
                        (cmei(i,k) + (pre(k) + prds(k))*cldmax(i,k))*deltat
              ttmp = t(i,k) +    &
                      ((pre(k)*cldmax(i,k))*xxlv + &
                        (cmei(i,k) + prds(k)*cldmax(i,k))*xxls)*deltat/cpp
#endif

!limit range of temperatures!
#ifdef GFDL_COMPATIBLE_MICROP
              ttmp = max(lowest_temp_for_sublimation, min(ttmp, 323._r8))
#else
              ttmp = max(180._r8, min(ttmp, 323._r8))
#endif
              esn = polysvp(ttmp,0) ! use rhw to allow ice supersaturation
              qsn = min(epsqs*esn/(p(i,k) - (1._r8 - epsqs)*esn), 1._r8)
#ifdef GFDL_COMPATIBLE_MICROP
!RSH 8/1/12: Need to prevent negative values which may occur at low p
            if( .not. lflag ) &
              qsn = max(qsn, 0._r8)
#endif

! modify precip evaporation rate if q > qsat
              if (qtmp .gt. qsn) then
                if (pre(k) + prds(k) .lt. -1.e-20_r8) then
                  dum1 = pre(k)/(pre(k) + prds(k))
! recalculate q and t after cloud water cond but without precip evap
#ifdef GFDL_COMPATIBLE_MICROP
            if( .not. lflag ) then
                  qtmp = q(i,k) -    &
                               (D_eros_l(i,k) + D_eros_i(i,k) +   &
                                          cmel(i,k) + cmei(i,k))*deltat
                  ttmp = t(i,k)+    &
                            ((D_eros_l(i,k) + cmel(i,k))*xxlv + &
                               (D_eros_i(i,k) + cmei(i,k))*xxls)*deltat/cpp
            else
                  qtmp = q(i,k) - (cmei(i,k))*deltat
                  ttmp = t(i,k) + (cmei(i,k)*xxls)*deltat/cpp
            endif 
#else
                  qtmp = q(i,k) - (cmei(i,k))*deltat
                  ttmp = t(i,k) + (cmei(i,k)*xxls)*deltat/cpp
#endif
                  esn = polysvp(ttmp,0) ! use rhw to allow ice 
                                        ! supersaturation
                  qsn = min(epsqs*esn/(p(i,k) - (1._r8 - epsqs)*esn),  &
                                                                     1._r8)
!RSH 8/1/12:  Need to prevent negative values which may occur at low p
#ifdef GFDL_COMPATIBLE_MICROP
            if( .not. lflag ) &
                  qsn = max(qsn, 0._r8) 
#endif

                  dum = (qtmp - qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
                  dum = min(dum, 0._r8)

! modify rates if needed, divide by cldmax to get local (in-precip) value
                  pre(k) = dum*dum1/deltat/cldmax(i,k)

! do separately using RHI for prds....
                  esn = polysvp(ttmp,1) ! use rhi to allow ice 
                                        ! supersaturation
                  qsn = min(epsqs*esn/(p(i,k) - (1._r8 - epsqs)*esn),  &
                                                                    1._r8)
                  dum = (qtmp - qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
                  dum = min(dum, 0._r8)

! modify rates if needed, divide by cldmax to get local (in-precip) value
                  prds(k) = dum*(1._r8 - dum1)/deltat/cldmax(i,k)
                end if  ! pre(k)+prds(k).lt.-1.e-20_r8 
              end if  !qtmp.gt.qsn 
            end if   ! qcic(i,k)+qiic(i,k).lt.qsmall .or.cldmax(i,k)...

! bergeron process - evaporation of droplets and deposition onto snow

            if (qniic(i,k) .ge. qsmall .and.   &
                 qcic(i,k) .ge. qsmall .and.   &
                                           t(i,k) .lt. tmelt) then
              qvi = epsqs*esi(i,k)/(p(i,k) - (1._r8 - epsqs)*esi(i,k))
              qvs = epsqs*esl(i,k)/(p(i,k) - (1._r8 - epsqs)*esl(i,k))

#ifdef GFDL_COMPATIBLE_MICROP
!8/1/12 RSH:   Need to prevent negative values which may occur at low p
            if( .not. lflag ) then
              qvi = MAX(0._r8, MIN (qvi,1.0_r8))
              qvs = MAX(0._r8, MIN (qvs,1.0_r8))
            endif
#endif
              dqsidt =  xxls*qvi/(rv*t(i,k)**2)
              abi = 1._r8 + dqsidt*xxls/cpp
              epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                                          (f1s/(lams(k)*lams(k)) + &
                             f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                                   sc(i,k)**(1._r8/3._r8)*cons14/ &
                                      (lams(k)**(5._r8/2._r8 + bs/2._r8)))
              bergs(k) = epss*(qvs - qvi)/abi
            else
              bergs(k) = 0._r8
            end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! conservation to ensure no negative values of cloud water/precipitation
! in case microphysical process rates are large

! make sure and use end-of-time step values for cloud water, ice, due
! condensation/deposition

! note: for check on conservation, processes are multiplied by omsm
! to prevent problems due to round off error

! include mixing timescale  (mtime)

#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
            qce = (qc(i,k) +    &
                            (D_eros_l(i,k) + cmel(i,k) - berg(i,k))*deltat)
          else
            qce = (qc(i,k) - berg(i,k)*deltat)
          endif 
#else
            qce = (qc(i,k) - berg(i,k)*deltat)
#endif
            nce = (nc(i,k) + npccn(k)*deltat*mtime)
#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
            qie = (qi(i,k) +    &
                           (D_eros_i(i,k) + cmei(i,k) + berg(i,k))*deltat)
          else
            qie = (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat)
          endif
#else
            qie = (qi(i,k) + (cmei(i,k) + berg(i,k))*deltat)
#endif
            nie = (ni(i,k) + nnuccd(k)*deltat*mtime)

! conservation of qc

            dum = (prc(k) + pra(k) + mnuccc(k) + mnucct(k) + msacwi(k) + &
                                 psacws(k) + bergs(k))*lcldm(i,k)*deltat
            if (dum .gt. qce) then
              ratio = qce/deltat/lcldm(i,k)/    &
                          (prc(k) + pra(k) + mnuccc(k) + mnucct(k) +  &
                                    msacwi(k) + psacws(k) + bergs(k))*omsm 
              prc(k) = prc(k)*ratio
              pra(k) = pra(k)*ratio
              mnuccc(k) = mnuccc(k)*ratio
              mnucct(k) = mnucct(k)*ratio
              msacwi(k) = msacwi(k)*ratio
              psacws(k) = psacws(k)*ratio
              bergs(k) = bergs(k)*ratio
            end if

! conservation of nc

#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
            dum = (nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) + &
                    npsacws(k) - nsubc(k) - nerosc(i,k))*lcldm(i,k)*deltat
          else
            dum = (nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) + &
                                npsacws(k) - nsubc(k))*lcldm(i,k)*deltat
          endif
#else
            dum = (nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) + &
                                npsacws(k) - nsubc(k))*lcldm(i,k)*deltat
#endif 

            if (dum .gt. nce) then
#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
              ratio = nce/deltat/   &
                     ((nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) +  &
                     npsacws(k) - nsubc(k) - nerosc(i,k))*lcldm(i,k))*omsm
          else
              ratio = nce/deltat/   &
                       ((nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) +   &
                                   npsacws(k) - nsubc(k))*lcldm(i,k))*omsm
          endif
#else
              ratio = nce/deltat/   &
                       ((nprc1(k) + npra(k) + nnuccc(k) + nnucct(k) +   &
                                   npsacws(k) - nsubc(k))*lcldm(i,k))*omsm
#endif
              nprc1(k) = nprc1(k)*ratio
              npra(k) = npra(k)*ratio
              nnuccc(k) = nnuccc(k)*ratio
              nnucct(k) = nnucct(k)*ratio  
              npsacws(k) = npsacws(k)*ratio
              nsubc(k)=nsubc(k)*ratio
#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag )  &    
              nerosc(i,k)=nerosc(i,k)*ratio
#endif
            end if

! conservation of qi

            dum = ((-mnuccc(k) - mnucct(k) - msacwi(k))*lcldm(i,k) +   &
                                  (prci(k) + prai(k))*icldm(i,k))*deltat
            if (dum .gt. qie) then
              ratio = (qie/deltat +    &
                        (mnuccc(k) + mnucct(k) + msacwi(k))*lcldm(i,k))/  &
                                     ((prci(k) + prai(k))*icldm(i,k))*omsm 
              prci(k) = prci(k)*ratio
              prai(k) = prai(k)*ratio
            end if

! conservation of ni

#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
            dum = ((-nnucct(k) - nsacwi(k))*lcldm(i,k) + (nprci(k) + &
                     nprai(k) - nsubi(k) - nerosi(i,k))*icldm(i,k))*deltat
          else
            dum = ((-nnucct(k) - nsacwi(k))*lcldm(i,k) + (nprci(k) + &
                                 nprai(k) - nsubi(k))*icldm(i,k))*deltat
          endif 
#else
            dum = ((-nnucct(k) - nsacwi(k))*lcldm(i,k) + (nprci(k) + &
                                 nprai(k) - nsubi(k))*icldm(i,k))*deltat
#endif
            if (dum .gt. nie) then

#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
              ratio = (nie/deltat + (nnucct(k) + nsacwi(k))*lcldm(i,k))/ &
                         ((nprci(k) + nprai(k) - nsubi(k) - nerosi(i,k))* &
                                                          icldm(i,k))*omsm
          else
              ratio = (nie/deltat + (nnucct(k) + nsacwi(k))*lcldm(i,k))/ &
                        ((nprci(k) + nprai(k) - nsubi(k))*icldm(i,k))*omsm
          endif       
#else
              ratio = (nie/deltat + (nnucct(k) + nsacwi(k))*lcldm(i,k))/ &
                        ((nprci(k) + nprai(k) - nsubi(k))*icldm(i,k))*omsm
#endif
              nprci(k) = nprci(k)*ratio
              nprai(k) = nprai(k)*ratio
              nsubi(k) = nsubi(k)*ratio
#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag )  &  
              nerosi(i,k) = nerosi(i,k)*ratio
#endif
            end if

! for preciptiation conservation, use logic that vertical integral 
! of tendency from current level to top of model (i.e., qrtot) cannot be negative

! conservation of rain mixing rat

            if (((prc(k) + pra(k))*lcldm(i,k) +   &
                    (-mnuccr(k) + pre(k) - pracs(k))*cldmax(i,k))*  &
                                 dz(i,k)*rho(i,k) + qrtot .lt. 0._r8) then
              if (-pre(k) + pracs(k) + mnuccr(k) .ge. qsmall) then
                ratio = (qrtot/(dz(i,k)*rho(i,k)) +   &
                                       (prc(k) + pra(k))*lcldm(i,k))/&
                       ((-pre(k) + pracs(k) + mnuccr(k))*cldmax(i,k))*omsm 
                pre(k) = pre(k)*ratio
                pracs(k) = pracs(k)*ratio
                mnuccr(k) = mnuccr(k)*ratio
              end if  ! -pre(k)+pracs(k)+mnuccr(k).ge.qsmall
            end if

! conservation of nr
! for now neglect evaporation of nr
            nsubr(k)=0._r8
!--> cjg: modifications incorporated from Huan's code
            if (allow_rain_num_evap)  then
              if (qric(i,k) .ge. qsmall) nsubr(k)= max(pre(k)/qric(i,k)*nric(i,k), -nric(i,k)/deltat)
            endif
!<--cjg

            if ((nprc(k)*lcldm(i,k) + (-nnuccr(k) + nsubr(k) -   &
                   npracs(k) + nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k) +  &
                                                    nrtot .lt. 0._r8) then
              if (-nsubr(k) - nragg(k) + npracs(k) +   &
                                               nnuccr(k) .ge.qsmall) then
                ratio = (nrtot/(dz(i,k)*rho(i,k)) + nprc(k)*lcldm(i,k))/&
                            ((-nsubr(k) - nragg(k) + npracs(k) +    &
                                               nnuccr(k))*cldmax(i,k))*omsm
                nsubr(k) = nsubr(k)*ratio
                npracs(k) = npracs(k)*ratio
                nnuccr(k) = nnuccr(k)*ratio
                nragg(k) = nragg(k)*ratio
              end if  
            end if

! conservation of snow mix ratio

            if (((bergs(k) + psacws(k))*lcldm(i,k) +    &
                          (prai(k) + prci(k))*icldm(i,k) + (pracs(k) + &
                            mnuccr(k) + prds(k))*cldmax(i,k))*dz(i,k)*  &
                                         rho(i,k) + qstot .lt. 0._r8) then
              if (-prds(k) .ge. qsmall) then
                ratio = (qstot/(dz(i,k)*rho(i,k)) +    &
                            (bergs(k) + psacws(k))*lcldm(i,k) +    &
                                        (prai(k) + prci(k))*icldm(i,k) + &
                                 (pracs(k) + mnuccr(k))*cldmax(i,k))/   &
                                               (-prds(k)*cldmax(i,k))*omsm
                prds(k) = prds(k)*ratio
              end if  ! -prds(k).ge.qsmall
            end if

! conservation of ns

! calculate loss of number due to sublimation
! for now neglect sublimation of ns
            nsubs(k)=0._r8
            if ((nprci(k)*icldm(i,k) + (nnuccr(k) + nsubs(k) +   &
                                           nsagg(k))*cldmax(i,k))*&
                                dz(i,k)*rho(i,k) + nstot .lt. 0._r8) then
              if (-nsubs(k) - nsagg(k) .ge. qsmall) then
                ratio = (nstot/(dz(i,k)*rho(i,k)) + nprci(k)*   &
                             icldm(i,k) + nnuccr(k)*cldmax(i,k))/   &
                               ((-nsubs(k) - nsagg(k))*cldmax(i,k))*omsm
                nsubs(k) = nsubs(k)*ratio
                nsagg(k) = nsagg(k)*ratio
              end if
            end if

! get tendencies due to microphysical conversion processes
! note: tendencies are multiplied by appropaiate cloud/precip 
! fraction to get grid-scale values
! note: cmei is already grid-average values

#ifdef GFDL_COMPATIBLE_MICROP
        if( .not. lflag ) then
            qvlat(i,k) = qvlat(i,k) - &
                          (pre(k) + prds(k))*cldmax(i,k) - cmel(i,k) -  &
                           cmei(i,k) - D_eros_l(i,k) - D_eros_i(i,k) 
            tlat(i,k) = tlat(i,k) + ((pre(k)*cldmax(i,k) + cmel(i,k) +  &
                          D_eros_l(i,k))*xxlv + (prds(k)*cldmax(i,k) +  &
                          cmei(i,k) + D_eros_i(i,k))*xxls + &
                          ((bergs(k) + psacws(k) + mnuccc(k) +   &
                            mnucct(k) + msacwi(k))*lcldm(i,k) +    &
                             (mnuccr(k) + pracs(k))*cldmax(i,k) +   &
                                                           berg(i,k))*xlf)
            qctend(i,k) = qctend(i,k) + &
                            (-pra(k) - prc(k) - mnuccc(k) - mnucct(k) -   &
                              msacwi(k) - psacws(k) - bergs(k))*   &
                                                           lcldm(i,k) +   &
                                    cmel(i,k) - berg(i,k) + D_eros_l(i,k)
            qitend(i,k) = qitend(i,k) + &
                            (mnuccc(k) + mnucct(k) + msacwi(k))*  &
                                                           lcldm(i,k) +   &
                             (-prci(k) - prai(k))*icldm(i,k) +   &
                                     cmei(i,k) + berg(i,k) + D_eros_i(i,k)
        else
            qvlat(i,k) = qvlat(i,k) -   &
                            (pre(k) + prds(k))*cldmax(i,k) - cmei(i,k) 
            tlat(i,k) = tlat(i,k) + ((pre(k)*cldmax(i,k))*xxlv +  &
                           (prds(k)*cldmax(i,k) + cmei(i,k))*xxls + &
                          ((bergs(k) + psacws(k) + mnuccc(k) +    &
                            mnucct(k) + msacwi(k))*lcldm(i,k) +   &
                           (mnuccr(k) + pracs(k))*cldmax(i,k) +   &
                                                          berg(i,k))*xlf)
            qctend(i,k) = qctend(i,k) + &
                           (-pra(k) - prc(k) - mnuccc(k) - mnucct(k) -  &
                             msacwi(k) - psacws(k) - bergs(k))*  &
                                                   lcldm(i,k) - berg(i,k) 
            qitend(i,k) = qitend(i,k) + &
                           (mnuccc(k) + mnucct(k) + msacwi(k))*    &
                                                           lcldm(i,k) +    &
                           (-prci(k) - prai(k))*icldm(i,k)    &
                                                   + cmei(i,k) + berg(i,k)

        endif
#else
            qvlat(i,k) = qvlat(i,k) -   &
                            (pre(k) + prds(k))*cldmax(i,k) - cmei(i,k) 
            tlat(i,k) = tlat(i,k) + ((pre(k)*cldmax(i,k))*xxlv +  &
                           (prds(k)*cldmax(i,k) + cmei(i,k))*xxls + &
                          ((bergs(k) + psacws(k) + mnuccc(k) +    &
                            mnucct(k) + msacwi(k))*lcldm(i,k) +   &
                           (mnuccr(k) + pracs(k))*cldmax(i,k) +   &
                                                          berg(i,k))*xlf)
            qctend(i,k) = qctend(i,k) + &
                           (-pra(k) - prc(k) - mnuccc(k) - mnucct(k) -  &
                             msacwi(k) - psacws(k) - bergs(k))*  &
                                                   lcldm(i,k) - berg(i,k) 
            qitend(i,k) = qitend(i,k) + &
                           (mnuccc(k) + mnucct(k) + msacwi(k))*    &
                                                           lcldm(i,k) +   &
                           (-prci(k) - prai(k))*icldm(i,k)    &
                                                   + cmei(i,k) + berg(i,k)
#endif
            qrtend(i,k) = qrtend(i,k) + &
                           (pra(k) + prc(k))*lcldm(i,k) + (pre(k) -   &
                           pracs(k) - mnuccr(k))*cldmax(i,k)
            qnitend(i,k) = qnitend(i,k) + &
                            (prai(k) + prci(k))*icldm(i,k) +   &
                              (psacws(k) + bergs(k))*lcldm(i,k) +   &
                               (prds(k) + pracs(k) + mnuccr(k))*cldmax(i,k)
! add output for cmei (accumulate)
            cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)

! assign variables for trop_mozart, these are grid-average
! evaporation/sublimation is stored here as positive term

            evapsnow(i,k) = evapsnow(i,k) - prds(k)*cldmax(i,k)
            nevapr(i,k) = nevapr(i,k) - pre(k)*cldmax(i,k)

! change to make sure prain is positive: do not remove snow from
! prain used for wet deposition
            prain(i,k) = prain(i,k) + (pra(k) + prc(k))*lcldm(i,k) +  &
                                      (-pracs(k) - mnuccr(k))*cldmax(i,k)
            prodsnow(i,k) = prodsnow(i,k) + (prai(k) +    &
                             prci(k))*icldm(i,k) +    &
                                   (psacws(k) + bergs(k))*lcldm(i,k) + &
                                        (pracs(k) + mnuccr(k))*cldmax(i,k)

! following are used to calculate 1st order conversion rate of cloud water
! to rain and snow (1/s), for later use in aerosol wet removal routine
! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may 
! be smaller than the qc used to calculate pra, prc, ... in this routine
! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of 
!                      cloud water to rain & snow }
!                      (no cloud ice or bergeron terms)
! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the 
!                      transfer terms }

            qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) +    &
                                        (pra(k) + prc(k) + psacws(k))*  &
                                                                 lcldm(i,k)
            qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k) 

! microphysics output, note this is grid-averaged
#ifdef GFDL_COMPATIBLE_MICROP
            preo(i,k) = preo(i,k) + pre(k)*cldmax(i,k)
            prdso(i,k) = prdso(i,k) + prds(k)*cldmax(i,k)
            if( .not. lflag ) then
            cmelo(i,k) = cmelo(i,k) + cmel(i,k)
            eroslo(i,k) = eroslo(i,k) + D_eros_l(i,k)
            erosio(i,k) = erosio(i,k) + D_eros_i(i,k)
            endif
#endif
            prao(i,k) = prao(i,k) + pra(k)*lcldm(i,k)
            prco(i,k) = prco(i,k) + prc(k)*lcldm(i,k)
            mnuccco(i,k) = mnuccco(i,k) + mnuccc(k)*lcldm(i,k)
            mnuccto(i,k) = mnuccto(i,k) + mnucct(k)*lcldm(i,k)
            mnuccdo(i,k) = mnuccdo(i,k) + mnuccd(k)*lcldm(i,k)
            msacwio(i,k) = msacwio(i,k) + msacwi(k)*lcldm(i,k)
            psacwso(i,k) = psacwso(i,k) + psacws(k)*lcldm(i,k)
            bergso(i,k) = bergso(i,k) + bergs(k)*lcldm(i,k)
            bergo(i,k) = bergo(i,k) + berg(i,k)
            prcio(i,k) = prcio(i,k) + prci(k)*icldm(i,k)
            praio(i,k) = praio(i,k) + prai(k)*icldm(i,k)
            mnuccro(i,k) = mnuccro(i,k) + mnuccr(k)*cldmax(i,k)
            pracso (i,k) = pracso (i,k) + pracs (k)*cldmax(i,k)

! multiply activation/nucleation by mtime to account for fast timescale

#ifdef GFDL_COMPATIBLE_MICROP
         if( .not. lflag ) then
            nctend(i,k) = nctend(i,k) + npccn(k)*mtime + &
                          (-nnuccc(k) - nnucct(k) - npsacws(k) +   &
                            nsubc(k) + nerosc(i,k) - npra(k) - nprc1(k))* &
                                                                 lcldm(i,k)
            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime + &
                           (nnucct(k) + nsacwi(k))*lcldm(i,k) +   &
                            (nsubi(k) + nerosi(i,k) - nprci(k) - &
                                                       nprai(k))*icldm(i,k)
         else
            nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
                  (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) &
                  -npra(k)-nprc1(k))*lcldm(i,k)

            nitend(i,k) = nitend(i,k)+ nnuccd(k)*mtime+ &
                  (nnucct(k)+nsacwi(k))*lcldm(i,k)+(nsubi(k)-nprci(k)- &
                  nprai(k))*icldm(i,k)

         endif 
#else
!#endif
!#ifndef GFDL_COMPATIBLE_MICROP
            nctend(i,k) = nctend(i,k) + npccn(k)*mtime +  &  
                           (-nnuccc(k) - nnucct(k) - npsacws(k) +  &
                            nsubc(k) - npra(k) - nprc1(k))*lcldm(i,k)
            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime +  &
                           (nnucct(k) + nsacwi(k))*lcldm(i,k) +    &
                            (nsubi(k) - nprci(k) - nprai(k))*icldm(i,k)
#endif

            nstend(i,k) = nstend(i,k) + (nsubs(k) +    &
                           nsagg(k) + nnuccr(k))*cldmax(i,k) +   &
                                                       nprci(k)*icldm(i,k)
            nrtend(i,k) = nrtend(i,k) + &
                           nprc(k)*lcldm(i,k) + (nsubr(k) - npracs(k) -  &
                                          nnuccr(k) + nragg(k))*cldmax(i,k)

#ifdef GFDL_COMPATIBLE_MICROP
! save current tendencies so that upcoming adjustment may be captured
            IF (diag_id%qndt_nucclim +     &
                                      diag_id%qn_nucclim_col  > 0) THEN
              nucclim(k) = nctend(i,k)
            END IF

            IF (diag_id%qnidt_nucclim1 +    &
                                       diag_id%qni_nucclim1_col > 0) THEN
              nucclim1i(k) = nitend(i,k)
            END IF
#endif
! make sure that nc and ni at advanced time step do not exceed
! maximum (existing N + source terms*dt), which is possible due to
! fast nucleation timescale

            if (nctend(i,k) .gt. 0._r8 .and.    &
                    nc(i,k) + nctend(i,k)*deltat .gt. ncmax) then
              nctend(i,k) = max(0._r8, (ncmax - nc(i,k))/deltat)
            end if
            if (nitend(i,k) .gt. 0._r8 .and.    &
                    ni(i,k) + nitend(i,k)*deltat .gt. nimax) then
              nitend(i,k) = max(0._r8, (nimax - ni(i,k))/deltat)
            end if

#ifdef GFDL_COMPATIBLE_MICROP
! complete diagnostic calculation 
            IF (diag_id%qndt_nucclim +     &
                                    diag_id%qn_nucclim_col  > 0) THEN
              nucclim(k) = nctend(i,k) - nucclim(k)
            END IF
            IF (diag_id%qnidt_nucclim1 +     &
                                    diag_id%qni_nucclim1_col > 0) THEN
              nucclim1i(k) = nitend(i,k) - nucclim1i(k)
            END IF
#endif

! get final values for precipitation q and N, based on
! flux of precip from above, source/sink term, and terminal fallspeed
! see eq. 15-16 in MG2008

! rain

            if (qric(i,k) .ge. qsmall) then
              if (k .eq. 1) then
                qric(i,k) = qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
                nric(i,k) = nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
              else
                qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*   &
                                                       cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*qrtend(i,k)))/   &
                                            (umr(k)*rho(i,k)*cldmax(i,k))
                nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*   &
                                                      cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*nrtend(i,k)))/    &
                                            (unr(k)*rho(i,k)*cldmax(i,k))
              end if
            else
              qric(i,k)=0._r8
              nric(i,k)=0._r8
            end if

! snow
            if (qniic(i,k) .ge. qsmall) then
              if (k .eq. 1) then
                qniic(i,k) = qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
                nsic(i,k) = nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
              else
                qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*   &
                                                      cldmax(i,k-1) + &
                              (rho(i,k)*dz(i,k)*qnitend(i,k)))/   &
                                             (ums(k)*rho(i,k)*cldmax(i,k))
                nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*    &
                                                     cldmax(i,k-1) + &
                             (rho(i,k)*dz(i,k)*nstend(i,k)))/     &
                                             (uns(k)*rho(i,k)*cldmax(i,k))
              end if
            else
              qniic(i,k)=0._r8
              nsic(i,k)=0._r8
            end if

! calculate precipitation flux at surface
! divide by density of water to get units of m/s

            prect(i) = prect(i) + (qrtend(i,k)*dz(i,k)*rho(i,k) +   &
                       qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
            preci(i) = preci(i) + qnitend(i,k)*dz(i,k)*rho(i,k)/rhow
!RSH 4/3/12
!       prect(i) = max(prect(i), 0._r8)
!       preci(i) = max(preci(i), 0._r8)

! convert rain rate from m/s to mm/hr
            rainrt(i,k) = qric(i,k)*rho(i,k)*umr(k)/rhow*3600._r8*1000._r8

! vertically-integrated precip source/sink terms (note: grid-averaged)
#ifdef GFDL_COMPATIBLE_MICROP
      if( .not. lflag ) then
!RSH 4/3/12
!       qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!       qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
!       nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
!       nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
            qrtot = (qrtot + qrtend(i,k)*dz(i,k)*rho(i,k))
            qstot = (qstot + qnitend(i,k)*dz(i,k)*rho(i,k))
            nrtot = (nrtot + nrtend(i,k)*dz(i,k)*rho(i,k))
            nstot = (nstot + nstend(i,k)*dz(i,k)*rho(i,k))
      else
        qrtot = max(qrtot + qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
        qstot = max(qstot + qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
        nrtot = max(nrtot + nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
        nstot = max(nstot + nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
      endif
#else
        qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
        qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
        nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
        nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)
#endif
! calculate melting and freezing of precip 
! melt snow at +2 C

            if (t(i,k) + tlat(i,k)/cpp*deltat > tmelt + 2._r8) then
              if (qstot > 0._r8) then

! make sure melting snow doesn't reduce temperature below threshold
                dum = -xlf/cpp*qstot/(dz(i,k)*rho(i,k))
                if (t(i,k) + tlat(i,k)/cpp*deltat     &
                                           + dum .lt. tmelt + 2._r8) then
                  dum = (t(i,k) + tlat(i,k)/cpp*deltat -   &
                                                  (tmelt + 2._r8))*cpp/xlf
                  dum = dum/(xlf/cpp*qstot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8, dum)
                  dum = min(1._r8, dum)
                else
                  dum = 1._r8
                end if
                qric(i,k) = qric(i,k) + dum*qniic(i,k)
                nric(i,k) = nric(i,k) + dum*nsic(i,k)
                qniic(i,k) = (1._r8 - dum)*qniic(i,k)
                nsic(i,k) = (1._r8 - dum)*nsic(i,k)
! heating tendency 
                tmp = -xlf*dum*qstot/(dz(i,k)*rho(i,k))
                meltsdt(i,k) = meltsdt(i,k) + tmp
                tlat(i,k) = tlat(i,k) + tmp

#ifdef GFDL_COMPATIBLE_MICROP
                if (diag_id%snow_melt + diag_id%snow_melt_col > 0) & 
                  diag_4d(i,j,k, diag_pt%snow_melt) =   &
                         diag_4d(i,j,k, diag_pt%snow_melt) +   &
                                 dum*preci(i)*rhow/(rho(i,k)*dz(i,k))
#endif
                qrtot = qrtot + dum*qstot
                nrtot = nrtot + dum*nstot
                qstot = (1._r8 - dum)*qstot
                nstot = (1._r8 - dum)*nstot
                preci(i) = (1._r8 - dum)*preci(i)
              end if
            end if

! freeze all rain at -5C for Arctic

            if (t(i,k) + tlat(i,k)/cpp*deltat < (tmelt - 5._r8)) then
              if (qrtot > 0._r8) then

! make sure freezing rain doesn't increase temperature above threshold
                dum = xlf/cpp*qrtot/(dz(i,k)*rho(i,k))
                if (t(i,k) + tlat(i,k)/cpp*deltat +    &
                                            dum .gt. (tmelt - 5._r8)) then
                  dum = -(t(i,k) + tlat(i,k)/cpp*deltat -    &
                                                 (tmelt - 5._r8))*cpp/xlf
                  dum = dum/(xlf/cpp*qrtot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8, dum)
                  dum = min(1._r8, dum)
                else
                  dum = 1._r8
                end if
            
                qniic(i,k) = qniic(i,k) + dum*qric(i,k)
                nsic(i,k) = nsic(i,k) + dum*nric(i,k)
                qric(i,k) = (1._r8 - dum)*qric(i,k)
                nric(i,k) = (1._r8 - dum)*nric(i,k)
! heating tendency 
                tmp = xlf*dum*qrtot/(dz(i,k)*rho(i,k))
                frzrdt(i,k) = frzrdt(i,k) + tmp

                tlat(i,k) = tlat(i,k) + tmp
#ifdef GFDL_COMPATIBLE_MICROP
                diag_4d(i,j,k, diag_pt%rain_freeze) =  &
                         diag_4d(i,j,k, diag_pt%rain_freeze) +  &
                                            dum*(prect(i) - preci(i))*  &
                                                   rhow /(rho(i,k)*dz(i,k))
#endif
                qstot = qstot + dum*qrtot
                qrtot = (1._r8 - dum)*qrtot
                nstot = nstot + dum*nrtot
                nrtot = (1._r8 - dum)*nrtot
                preci(i) = preci(i) + dum*(prect(i) - preci(i))
              end if  !qrtot > 0._r8 
            end if ! t(i,k)+tlat(i,k)/cp*deltat < (tmelt - 5._r8)

! if rain/snow mix ratio is zero so should number concentration

            if (qniic(i,k) .lt. qsmall) then
              qniic(i,k) = 0._r8
              nsic(i,k) = 0._r8
            end if

            if (qric(i,k) .lt. qsmall) then
              qric(i,k) = 0._r8
              nric(i,k) = 0._r8
            end if

! make sure number concentration is a positive number to avoid 
! taking root of negative

            nric(i,k) = max(nric(i,k), 0._r8)
            nsic(i,k) = max(nsic(i,k), 0._r8)

!.......................................................................
! get size distribution parameters for fallspeed calculations
!......................................................................
! rain

            if (qric(i,k) .ge. qsmall) then
              lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
              n0r(k) = nric(i,k)*lamr(k)

! check for slope
! change lammax and lammin for rain and snow
! adjust vars

              if (lamr(k) .lt. lamminr) then
                lamr(k) = lamminr
                n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                nric(i,k) = n0r(k)/lamr(k)
              else if (lamr(k) .gt. lammaxr) then
                lamr(k) = lammaxr
                n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
                nric(i,k) = n0r(k)/lamr(k)
              end if

! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

              unr(k) = min(arn(i,k)*cons4/lamr(k)**br, 9.1_r8*rhof(i,k))
              umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),  &
                                                          9.1_r8*rhof(i,k))
            else
              lamr(k) = 0._r8
              n0r(k) = 0._r8
              umr(k)=0._r8
              unr(k)=0._r8
            end if

!calculate mean size of combined rain and snow

            if (lamr(k) .gt. 0._r8) then
              Artmp = n0r(k)*pi/(2*lamr(k)**3._r8)
            else 
              Artmp = 0._r8
            endif

            if (lamc(k) .gt. 0._r8) then
              Actmp = cdist1(k)*pi*gamma(pgam(k) + 3._r8)/     &
                                                   (4._r8*lamc(k)**2._r8)
            else 
              Actmp = 0._r8
            endif

!8/15        if (Actmp.gt.0_r8.or.Artmp.gt.0) then
            if (Actmp .gt. 0._r8 .or. Artmp .gt. 0._r8) then
              rercld(i,k) = rercld(i,k) + 3._r8*(qric(i,k) + qcic(i,k))/ &
                                          (4._r8*rhow*(Actmp + Artmp))
              arcld(i,k) = arcld(i,k) + 1._r8
            endif

!......................................................................
! snow

            if (qniic(i,k) .ge. qsmall) then
              lams(k) = (cons6*cs*nsic(i,k)/qniic(i,k))**(1._r8/ds)
              n0s(k) = nsic(i,k)*lams(k)

! check for slope
! adjust vars

              if (lams(k) .lt. lammins) then
                lams(k) = lammins
                n0s(k) = lams(k)**(ds + 1._r8)*qniic(i,k)/(cs*cons6)
                nsic(i,k) = n0s(k)/lams(k)
              else if (lams(k) .gt. lammaxs) then
                lams(k) = lammaxs
                n0s(k) = lams(k)**(ds + 1._r8)*qniic(i,k)/(cs*cons6)
                nsic(i,k) = n0s(k)/lams(k)
              end if

! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

              ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),   &
                                                         1.2_r8*rhof(i,k))
              uns(k) = min(asn(i,k)*cons7/lams(k)**bs, 1.2_r8*rhof(i,k))
            else
              lams(k) = 0._r8
              n0s(k) = 0._r8
              ums(k) = 0._r8
              uns(k) = 0._r8
            end if

!c........................................................................
! sum over sub-step for average process rates

! convert rain/snow q and N for output to history, note, 
! output is for gridbox average

            qrout(i,k) = qrout(i,k) + qric(i,k)*cldmax(i,k)
            qsout(i,k) = qsout(i,k) + qniic(i,k)*cldmax(i,k)
            nrout(i,k) = nrout(i,k) + nric(i,k)*rho(i,k)*cldmax(i,k)
            nsout(i,k) = nsout(i,k) + nsic(i,k)*rho(i,k)*cldmax(i,k)

            tlat1(i,k) = tlat1(i,k) + tlat(i,k)
            qvlat1(i,k) = qvlat1(i,k) + qvlat(i,k)
            qctend1(i,k) = qctend1(i,k) + qctend(i,k)
            qitend1(i,k) = qitend1(i,k) + qitend(i,k)
            nctend1(i,k) = nctend1(i,k) + nctend(i,k)
            nitend1(i,k) = nitend1(i,k) + nitend(i,k)

            t(i,k) = t(i,k) + tlat(i,k)*deltat/cpp
            q(i,k) = q(i,k) + qvlat(i,k)*deltat
            qc(i,k) = qc(i,k) + qctend(i,k)*deltat
            qi(i,k) = qi(i,k) + qitend(i,k)*deltat
            nc(i,k) = nc(i,k) + nctend(i,k)*deltat
            ni(i,k) = ni(i,k) + nitend(i,k)*deltat

            rainrt1(i,k) = rainrt1(i,k) + rainrt(i,k)

! divide rain radius over substeps for average
            if (arcld(i,k) .gt. 0._r8) then
              rercld(i,k) = rercld(i,k)/arcld(i,k)
            end if

! calculate precip fluxes and adding them to summing sub-stepping variables
! flux is zero at top interface
            rflx(i,1) = 0.0_r8
            sflx(i,1) = 0.0_r8

! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*
! airdensity(kg/m3)*massweightedfallspeed(m/s)
!!RSH 11/28/11 FIX OF BUG ??
            rflx(i,k+1)=qrout(i,k)*rho(i,k)*umr(k)
            sflx(i,k+1)=qsout(i,k)*rho(i,k)*ums(k)

#ifdef GFDL_COMPATIBLE_MICROP
      if( .not. lflag ) then
            rflx(i,k+1) = qric(i,k)*cldmax(i,k)*rho(i,k)*umr(k)
            sflx(i,k+1) = qniic(i,k)*cldmax(i,k)*rho(i,k)*ums(k)
      endif
#endif

!  add to summing sub-stepping variable
            rflx1(i,k+1) = rflx1(i,k+1) + rflx(i,k+1)
            sflx1(i,k+1) = sflx1(i,k+1) + sflx(i,k+1)

!c........................................................................

#ifdef GFDL_COMPATIBLE_MICROP
!droplet number
            npccno(i,k)     = npccno(i,k)   + npccn(k)*mtime
            nnuccco(i,k)    = nnuccco(i,k)  - nnuccc(k)*lcldm(i,k)
            nnuccto(i,k)    = nnuccto(i,k)  - nnucct(k)*lcldm(i,k)
            npsacwso(i,k)   = npsacwso(i,k) - npsacws(k)*lcldm(i,k)
            nsubco(i,k)     = nsubco(i,k)   + nsubc(k)*lcldm(i,k)
            if( .not. lflag ) &
            nerosco(i,k)    = nerosco(i,k)  + nerosc(i,k)*lcldm(i,k)
            nprao(i,k)      = nprao(i,k)    - npra(k)*lcldm(i,k)
            nprc1o(i,k)     = nprc1o(i,k)   - nprc1(k)*lcldm(i,k)
            nucclimo(i,k)   = nucclimo(i,k) + nucclim(k)

!ice number
            nnuccdo(i,k)    = nnuccdo(i,k) + nnuccd(k)*mtime
            nsacwio(i,k)    = nsacwio(i,k) + nsacwi(k)*lcldm(i,k)
            nsubio(i,k)     = nsubio(i,k)  + nsubi(k)*icldm(i,k)
            if( .not. lflag ) &
            nerosio(i,k)    = nerosio(i,k) + nerosi(i,k)*icldm(i,k)
            nprcio(i,k)     = nprcio(i,k)  - nprci(k)*icldm(i,k)
            npraio(i,k)     = npraio(i,k)  - nprai(k)*icldm(i,k)
            nucclim1io(i,k) = nucclim1io(i,k) + nucclim1i(k)
#endif

          end do ! k loop
!c........................................................................
          prect1(i) = prect1(i) + prect(i)
          preci1(i) = preci1(i) + preci(i)
!!!!!!!!!!!
!    END OF ITER LOOP
!!!!!!!!!!!
        end do ! it loop, sub-step

!----------------------------------------------------------------------

        do k = 1, pver               
          rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)/  &
                                     max(qcsum_rate1ord(k), 1.0e-30_r8) 
        end do                       

 300    continue  ! continue if no cloud water   ! GO TO 300

!!!!!!!!!!!
!    END OF I LOOP
!!!!!!!!!!!
      end do ! i loop

! convert dt from sub-step back to full time step
      deltat = deltat*real(iter)

!c........................................................................

      do i=1,ncol

! skip all calculations if no cloud water
        if (ltrue(i) .eq. 0) then

          do k=1,pver 
! assign default values for effective radius
            effc(i,k) = 10._r8
            effi(i,k) = 25._r8
            effc_fn(i,k) = 10._r8
            lamcrad(i,k) = 0._r8
            pgamrad(i,k) = 0._r8
            deffi(i,k) = 0._r8
          end do
#ifdef GFDL_COMPATIBLE_MICROP
          if (diag_id%vfall > 0) diag_4d(i,j,:,diag_pt%vfall) = 0.0_r8
#endif
          goto 500   ! EXIT FROM LOOP
        end if

! initialize nstep for sedimentation sub-steps
        nstep = 1

! divide precip rate by number of sub-steps to get average over time step

        prect(i) = prect1(i)/real(iter)
        preci(i) = preci1(i)/real(iter)

#ifdef GFDL_COMPATIBLE_MICROP
! convert unit  from m/s to kg/m2/s by multiply water density 1000 kg/m3
        diag_4d(i,j,:, diag_pt%snow_melt) =  &
                          diag_4d(i,j,:, diag_pt%snow_melt)/real(iter)
        diag_4d(i,j,:, diag_pt%rain_freeze) =  &
                          diag_4d(i,j,:, diag_pt%rain_freeze)/real(iter)
! <--- h1g, 2010-11-15
#endif

        do k=1,pver

! assign variables back to start-of-timestep values before updating 
! after sub-steps 

          t(i,k) = t1(i,k)
          q(i,k) = q1(i,k)
          qc(i,k)= qc1(i,k)
          qi(i,k) = qi1(i,k)
          nc(i,k) = nc1(i,k)
          ni(i,k) = ni1(i,k)

! divide microphysical tendencies by number of sub-steps to get average 
! over time step

          tlat(i,k) = tlat1(i,k)/real(iter)
          qvlat(i,k) = qvlat1(i,k)/real(iter)
          qctend(i,k) = qctend1(i,k)/real(iter)
          qitend(i,k) = qitend1(i,k)/real(iter)
          nctend(i,k) = nctend1(i,k)/real(iter)
          nitend(i,k) = nitend1(i,k)/real(iter)
 
          rainrt(i,k) = rainrt1(i,k)/real(iter)

! divide by number of sub-steps to find final values
          rflx(i,k+1) = rflx1(i,k+1)/real(iter)
          sflx(i,k+1) = sflx1(i,k+1)/real(iter)

! divide output precip q and N by number of sub-steps to get average over 
! time step

          qrout(i,k) = qrout(i,k)/real(iter)
          qsout(i,k) = qsout(i,k)/real(iter)
          nrout(i,k) = nrout(i,k)/real(iter)
          nsout(i,k) = nsout(i,k)/real(iter)

! divide trop_mozart variables by number of sub-steps to get average over 
! time step 

          nevapr(i,k) = nevapr(i,k)/real(iter)
          evapsnow(i,k) = evapsnow(i,k)/real(iter)
          prain(i,k) = prain(i,k)/real(iter)
          prodsnow(i,k) = prodsnow(i,k)/real(iter)
          cmeout(i,k) = cmeout(i,k)/real(iter)

          cmeiout(i,k) = cmeiout(i,k)/real(iter)
          meltsdt(i,k) = meltsdt(i,k)/real(iter)
          frzrdt (i,k) = frzrdt (i,k)/real(iter)

! microphysics output
#ifdef GFDL_COMPATIBLE_MICROP
! vapor
          preo(i,k) = preo(i,k)/real(iter)
          prdso(i,k) = prdso(i,k)/real(iter)
          if( .not. lflag ) then
          cmelo(i,k) =  cmelo(i,k)/real(iter)
          eroslo(i,k) = eroslo(i,k)/real(iter)
          erosio(i,k) = erosio(i,k)/real(iter)
          endif
#endif
! liquid
          prao(i,k) = prao(i,k)/real(iter)
          prco(i,k) = prco(i,k)/real(iter)
          mnuccco(i,k) = mnuccco(i,k)/real(iter)
          mnuccto(i,k) = mnuccto(i,k)/real(iter)
          msacwio(i,k) = msacwio(i,k)/real(iter)
          psacwso(i,k) = psacwso(i,k)/real(iter)
          bergso(i,k) = bergso(i,k)/real(iter)
          bergo(i,k) = bergo(i,k)/real(iter)
          prcio(i,k) = prcio(i,k)/real(iter)
          praio(i,k) = praio(i,k)/real(iter)

          mnuccro(i,k) = mnuccro(i,k)/real(iter)
          pracso (i,k) = pracso (i,k)/real(iter)

          mnuccdo(i,k) = mnuccdo(i,k)/real(iter)

! modify to include snow. in prain & evap (diagnostic here: for wet dep)
          nevapr(i,k) = nevapr(i,k) + evapsnow(i,k)
          prain(i,k) = prain(i,k) + prodsnow(i,k)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate sedimentation for cloud water and ice
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! update in-cloud cloud mixing ratio and number concentration 
! with microphysical tendencies to calculate sedimentation, assign to dummy
! vars
! note: these are in-cloud values***, hence we divide by cloud fraction

          dumc(i,k) = (qc(i,k) + qctend(i,k)*deltat)/lcldm(i,k)
          dumi(i,k) = (qi(i,k) + qitend(i,k)*deltat)/icldm(i,k)
          dumnc(i,k) = max((nc(i,k) + nctend(i,k)*deltat)/lcldm(i,k),0._r8)
          dumni(i,k) = max((ni(i,k) + nitend(i,k)*deltat)/icldm(i,k),0._r8)

!-->cjg
! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
        dumnc(i,k)=ncnst/rho(i,k)
        end if

! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
        dumni(i,k)=ninst/rho(i,k)
        end if
!<--cjg

! obtain new slope parameter to avoid possible singularity

          if (dumi(i,k) .ge. qsmall) then
! add upper limit to in-cloud number concentration to prevent 
! numerical error
            dumni(i,k) = min(dumni(i,k), dumi(i,k)*1.e20_r8)

            lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)
            lami(k) = max(lami(k),lammini)
            lami(k) = min(lami(k),lammaxi)
          else
            lami(k)=0._r8
          end if

          if (dumc(i,k) .ge. qsmall) then
! add upper limit to in-cloud number concentration to prevent 
! numerical error
            dumnc(i,k) = min(dumnc(i,k), dumc(i,k)*1.e20_r8)
! add lower limit to in-cloud number concentration
            dumnc(i,k) = max(dumnc(i,k), cdnl/rho(i,k)) ! sghan minimum 
                                                       ! in #/cm3 
#ifdef GFDL_COMPATIBLE_MICROP
         if( .not. lflag ) then
!RSH76     pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
            pgam(k) = 0.0005714_r8*(dumnc(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         else
           pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         endif
#else
           pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
#endif
            pgam(k) = 1._r8/(pgam(k)**2) - 1._r8
            pgam(k) = max(pgam(k), 2._r8)
            pgam(k) = min(pgam(k), 15._r8)

            lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k) + 4._r8)/ &
                       (dumc(i,k)*gamma(pgam(k) + 1._r8)))**(1._r8/3._r8)
            lammin = (pgam(k) + 1._r8)/50.e-6_r8
            lammax = (pgam(k) + 1._r8)/2.e-6_r8
            lamc(k) = max(lamc(k), lammin)
            lamc(k) = min(lamc(k), lammax)
          else
            lamc(k) = 0._r8
          end if

! calculate number and mass weighted fall velocity for droplets
! include effects of sub-grid distribution of cloud water


          if (dumc(i,k).ge.qsmall) then
            unc = &
                   acn(i,k)*gamma(1._r8 + bc+pgam(k))/ &
                          (lamc(k)**bc*gamma(pgam(k) + 1._r8))
            umc = &
                    acn(i,k)*gamma(4._r8 + bc+pgam(k))/ &
                           (lamc(k)**bc*gamma(pgam(k) + 4._r8))
! fallspeed for output
            vtrmc(i,k)=umc
          else
            umc = 0._r8
            unc = 0._r8
          end if

! calculate number and mass weighted fall velocity for cloud ice

          if (dumi(i,k).ge.qsmall) then
            uni =  ain(i,k)*cons16/lami(k)**bi
            umi = ain(i,k)*cons17/(6._r8*lami(k)**bi)
            uni = min(uni, 1.2_r8*rhof(i,k))
            umi = min(umi, 1.2_r8*rhof(i,k))

! fallspeed
            vtrmi(i,k) = umi
          else
            umi = 0._r8
            uni = 0._r8
          end if

#ifdef GFDL_COMPATIBLE_MICROP
          if (diag_id%vfall > 0) diag_4d(i,j,k,diag_pt%vfall) = umi
#endif

          fi(k) = g*rho(i,k)*umi
          fni(k) = g*rho(i,k)*uni
          fc(k) = g*rho(i,k)*umc
          fnc(k) = g*rho(i,k)*unc

! calculate number of split time steps to ensure courant stability criteria
! for sedimentation calculations

          rgvm = max(fi(k), fc(k), fni(k), fnc(k))
          nstep = max(int(rgvm*deltat/pdel(i,k) + 1._r8),nstep)

! redefine dummy variables - sedimentation is calculated over grid-scale
! quantities to ensure conservation

          dumc(i,k) = (qc(i,k) + qctend(i,k)*deltat)
          dumi(i,k) = (qi(i,k) + qitend(i,k)*deltat)
          dumnc(i,k) = max((nc(i,k) + nctend(i,k)*deltat), 0._r8)
          dumni(i,k) = max((ni(i,k) + nitend(i,k)*deltat), 0._r8)

          if (dumc(i,k) .lt. qsmall) dumnc(i,k) = 0._r8
          if (dumi(i,k) .lt. qsmall) dumni(i,k) = 0._r8

        end do       !!! vertical loop

        do n = 1,nstep  !! loop over sub-time step to ensure stability
          do k = 1,pver
             falouti(k) = fi(k)*dumi(i,k)
             faloutni(k) = fni(k)*dumni(i,k)
             faloutc(k) = fc(k)*dumc(i,k)
             faloutnc(k) = fnc(k)*dumnc(i,k)
          end do

! top of model

          k = 1
          faltndi = falouti(k)/pdel(i,k)
          faltndni = faloutni(k)/pdel(i,k)
          faltndc = faloutc(k)/pdel(i,k)
          faltndnc = faloutnc(k)/pdel(i,k)

! add fallout terms to microphysical tendencies

          qitend(i,k) = qitend(i,k) - faltndi/nstep
          nitend(i,k) = nitend(i,k) - faltndni/nstep
          qctend(i,k) = qctend(i,k) - faltndc/nstep
          nctend(i,k) = nctend(i,k) - faltndnc/nstep

! sedimentation tendencies for output
          qcsedten(i,k) = qcsedten(i,k) - faltndc/nstep
          qisedten(i,k) = qisedten(i,k) - faltndi/nstep

#ifdef GFDL_COMPATIBLE_MICROP
          IF ( diag_id%qndt_sedi + diag_id%qn_sedi_col > 0 ) &
                  diag_4d(i,j,k,diag_pt%qndt_sedi) =     &
                       diag_4d(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep

          IF ( diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0 ) &
                 diag_4d(i,j,k,diag_pt%qnidt_sedi) =    &
                      diag_4d(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep
#endif

          dumi(i,k) = dumi(i,k) - faltndi*deltat/nstep
          dumni(i,k) = dumni(i,k) - faltndni*deltat/nstep
          dumc(i,k) = dumc(i,k) - faltndc*deltat/nstep
          dumnc(i,k) = dumnc(i,k) - faltndnc*deltat/nstep

          do k = 2,pver

! for cloud liquid and ice, if cloud fraction increases with height
! then add flux from above to both vapor and cloud water of current level
! this means that flux entering clear portion of cell from above evaporates
! instantly

            dum = lcldm(i,k)/lcldm(i,k-1)
            dum = min(dum, 1._r8)
            dum1 = icldm(i,k)/icldm(i,k-1)
            dum1 = min(dum1, 1._r8)

            faltndqie = (falouti(k) - falouti(k-1))/pdel(i,k)
            faltndi = (falouti(k) - dum1*falouti(k-1))/pdel(i,k)
            faltndni = (faloutni(k) - dum1*faloutni(k-1))/pdel(i,k)
            faltndqce = (faloutc(k) - faloutc(k-1))/pdel(i,k)
            faltndc = (faloutc(k) - dum*faloutc(k-1))/pdel(i,k)
            faltndnc = (faloutnc(k) - dum*faloutnc(k-1))/pdel(i,k)

! add fallout terms to eulerian tendencies

            qitend(i,k) = qitend(i,k) - faltndi/nstep
            nitend(i,k) = nitend(i,k) - faltndni/nstep
            qctend(i,k) = qctend(i,k) - faltndc/nstep
            nctend(i,k) = nctend(i,k) - faltndnc/nstep

! sedimentation tendencies for output
            qcsedten(i,k) = qcsedten(i,k) - faltndc/nstep
            qisedten(i,k) = qisedten(i,k) - faltndi/nstep
 
! add terms to to evap/sub of cloud water

#ifdef GFDL_COMPATIBLE_MICROP
            IF ( diag_id%qndt_sedi  + diag_id%qn_sedi_col > 0 ) &
               diag_4d(i,j,k,diag_pt%qndt_sedi) =    &
                   diag_4d(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep

            IF ( diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0 ) &
                  diag_4d(i,j,k,diag_pt%qnidt_sedi) =    &
                        diag_4d(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep
#endif

            qvlat(i,k) = qvlat(i,k) - (faltndqie-faltndi)/nstep
! for output
            qisevap(i,k) = qisevap(i,k) - (faltndqie - faltndi)/nstep
            qvlat(i,k) = qvlat(i,k) - (faltndqce - faltndc)/nstep
! for output
            qcsevap(i,k) = qcsevap(i,k) - (faltndqce - faltndc)/nstep

            tlat(i,k) = tlat(i,k) + (faltndqie - faltndi)*xxls/nstep
            tlat(i,k) = tlat(i,k) + (faltndqce - faltndc)*xxlv/nstep

            dumi(i,k) = dumi(i,k) - faltndi*deltat/nstep
            dumni(i,k) = dumni(i,k) - faltndni*deltat/nstep
            dumc(i,k) = dumc(i,k) - faltndc*deltat/nstep
            dumnc(i,k) = dumnc(i,k) - faltndnc*deltat/nstep

            Fni(K) = MAX(Fni(K)/pdel(i,K), Fni(K-1)/pdel(i,K-1))*pdel(i,K)
            FI(K) = MAX(FI(K)/pdel(i,K), FI(K-1)/pdel(i,K-1))*pdel(i,K)
            fnc(k) = max(fnc(k)/pdel(i,k), fnc(k-1)/pdel(i,k-1))*pdel(i,k)
            Fc(K) = MAX(Fc(K)/pdel(i,K), Fc(K-1)/pdel(i,K-1))*pdel(i,K)
          end do   !! k loop

! units below are m/s
! cloud water/ice sedimentation flux at surface 
! is added to precip flux at surface to get total precip 
! (cloud + precip water) rate

          prect(i) = prect(i) + (faloutc(pver) + falouti(pver))/  &
                                                       g/nstep/1000._r8
          preci(i) = preci(i) + (falouti(pver)) &
                                                      /g/nstep/1000._r8

#ifdef GFDL_COMPATIBLE_MICROP
!  for the diagnostic of surface precip
          k = 1
          IF ( diag_id%sedi_sfc > 0 ) &
                  diag_4d(i,j,1,diag_pt%sedi_sfc) =    &
                          diag_4d(i,j,1,diag_pt%sedi_sfc) +  &
                                            faloutc(pver)/g/nstep/rhow 
          IF ( diag_id%sedi_ice > 0 ) &
                  diag_4d(i,j,1,diag_pt%sedi_ice) =   &
                           diag_4d(i,j,1,diag_pt%sedi_ice) +   &
                                             falouti(pver)/g/nstep/rhoi
#endif
        end do   !! nstep loop

! end sedimentation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! get new update for variables that includes sedimentation tendency
! note : here dum variables are grid-average, NOT in-cloud

        do k=1,pver

          dumc(i,k) = max(qc(i,k) + qctend(i,k)*deltat, 0._r8)
          dumi(i,k) = max(qi(i,k) + qitend(i,k)*deltat, 0._r8)
          dumnc(i,k) = max(nc(i,k) + nctend(i,k)*deltat, 0._r8)
          dumni(i,k) = max(ni(i,k) + nitend(i,k)*deltat, 0._r8)

!-->cjg
! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
        dumnc(i,k)=ncnst/rho(i,k)*lcldm(i,k)
        end if

! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
        dumni(i,k)=ninst/rho(i,k)*icldm(i,k)
        end if
!<--cjg

          if (dumc(i,k) .lt. qsmall) dumnc(i,k) = 0._r8
          if (dumi(i,k) .lt. qsmall) dumni(i,k) = 0._r8

! calculate instantaneous processes (melting, homogeneous freezing)

          if (t(i,k) + tlat(i,k)/cpp*deltat > tmelt) then
            if (dumi(i,k) > 0._r8) then

! limit so that melting does not push temperature below freezing
              dum = -dumi(i,k)*xlf/cpp
              if (t(i,k) + tlat(i,k)/cpp*deltat + dum .lt. tmelt) then
                dum = (t(i,k) + tlat(i,k)/cpp*deltat - tmelt)*cpp/xlf
                dum = dum/dumi(i,k)*xlf/cpp 
                dum = max(0._r8, dum)
                dum = min(1._r8, dum)
              else
                dum = 1._r8
              end if

              qctend(i,k) = qctend(i,k) + dum*dumi(i,k)/deltat

! for output
              melto(i,k) = dum*dumi(i,k)/deltat

! assume melting ice produces droplet
! mean volume radius of 8 micron

#ifdef GFDL_COMPATIBLE_MICROP
              IF (diag_id%qndt_melt + diag_id%qn_melt_col > 0) &
                   diag_4d(i,j,k,diag_pt%qndt_melt) =  nctend(i,k)
              IF (diag_id%qidt_melt2  + diag_id%qi_melt2_col  > 0) &
                   diag_4d(i,j,k,diag_pt%qidt_melt2) =  qitend(i,k)
              IF (diag_id%qnidt_melt +  diag_id%qni_melt_col > 0) &
                   diag_4d(i,j,k,diag_pt%qnidt_melt) =  nitend(i,k)
#endif

              nctend(i,k) = nctend(i,k) + 3._r8*dum*dumi(i,k)/deltat/ &
                                               (4._r8*pi*5.12e-16_r8*rhow)

              qitend(i,k) = ((1._r8 - dum)*dumi(i,k) - qi(i,k))/deltat
              nitend(i,k) = ((1._r8 - dum)*dumni(i,k) - ni(i,k))/deltat
              tlat(i,k) = tlat(i,k) - xlf*dum*dumi(i,k)/deltat

#ifdef GFDL_COMPATIBLE_MICROP
              IF (diag_id%qndt_melt + diag_id%qn_melt_col > 0) &
                   diag_4d(i,j,k,diag_pt%qndt_melt) =    &
                          nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_melt)
              IF (diag_id%qidt_melt2  + diag_id%qi_melt2_col  > 0) &
                   diag_4d(i,j,k,diag_pt%qidt_melt2) =    &
                          qitend(i,k) - diag_4d(i,j,k,diag_pt%qidt_melt2)
              IF (diag_id%qnidt_melt +  diag_id%qni_melt_col > 0) &
                   diag_4d(i,j,k,diag_pt%qnidt_melt) =    &
                          nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_melt)
#endif
            end if
          end if

! homogeneously freeze droplets at -40 C

          if (t(i,k) + tlat(i,k)/cpp*deltat < tmelt -40._r8) then
            if (dumc(i,k) > 0._r8) then

! limit so that freezing does not push temperature above threshold
              dum = dumc(i,k)*xlf/cpp
              if (t(i,k) + tlat(i,k)/cpp*deltat +     &
                                             dum .gt. tmelt - 40._r8) then
                dum = -(t(i,k) + tlat(i,k)/cpp*deltat - (tmelt -40._r8))* &
                                                                   cpp/xlf
                dum = dum/dumc(i,k)*xlf/cpp
                dum = max(0._r8, dum)
                dum = min(1._r8, dum)
              else
                dum = 1._r8
              end if

              qitend(i,k) = qitend(i,k) + dum*dumc(i,k)/deltat
! for output
              homoo(i,k) = dum*dumc(i,k)/deltat

! assume 25 micron mean volume radius of homogeneously frozen droplets
! consistent with size of detrained ice in stratiform.F90

#ifdef GFDL_COMPATIBLE_MICROP
              IF (diag_id%qldt_freez + diag_id%ql_freez_col > 0) &
                   diag_4d(i,j,k,diag_pt%qldt_freez) = qctend(i,k)
              sum_freeze(i,k) = qctend(i,k)
              IF ( diag_id%qndt_ihom  + diag_id%qn_ihom_col > 0) &
                    diag_4d(i,j,k,diag_pt%qndt_ihom) =  nctend(i,k)
#endif


! 4/24/12: replace the 1.563 with x**3, here and in nCAR routine.
              nitend(i,k) = nitend(i,k) + dum*3._r8*dumc(i,k)/   &
                            (4._r8*3.14_r8*1.563e-14_r8*500._r8)/deltat
              qctend(i,k) = ((1._r8 - dum)*dumc(i,k) - qc(i,k))/deltat
              nctend(i,k) = ((1._r8 - dum)*dumnc(i,k) - nc(i,k))/deltat
              tlat(i,k) = tlat(i,k) + xlf*dum*dumc(i,k)/deltat

#ifdef GFDL_COMPATIBLE_MICROP
              IF (diag_id%qldt_freez + diag_id%ql_freez_col > 0) &
                   diag_4d(i,j,k,diag_pt%qldt_freez) =    &
                         qctend(i,k) - diag_4d(i,j,k,diag_pt%qldt_freez) 
              sum_freeze(i,k) = -(qctend(i,k) - sum_freeze(i,k))
              IF (diag_id%qndt_ihom  + diag_id%qn_ihom_col > 0) &
                   diag_4d(i,j,k,diag_pt%qndt_ihom) =    &
                            nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_ihom)
! 4/24/12: replace the 1.563 with x**3, here and in nCAR routine.
              IF ( diag_id%qnidt_ihom +  diag_id%qni_ihom_col > 0 ) &
                    diag_4d(i,j,k,diag_pt%qnidt_ihom) =    &
                           dum*3._r8*dumc(i,k)/   &
                                (4._r8*3.14_r8*1.563e-14_r8*500._r8)/deltat
#endif
            end if
          end if

! remove any excess over-saturation, which is possible due to non-linearity
! when adding together all microphysical processes
! follow code similar to old CAM scheme

          qtmp = q(i,k) + qvlat(i,k)*deltat
          ttmp = t(i,k) + tlat(i,k)/cpp*deltat

          esn = polysvp(ttmp,0)  ! use rhw to allow ice supersaturation
          qsn = min(epsqs*esn/(p(i,k) - (1._r8 - epsqs)*esn), 1._r8)

#ifdef GFDL_COMPATIBLE_MICROP
          if (qtmp > qsn .and. qsn > 0._r8) then

! expression below is approximate since there may be ice deposition
            dum = (qtmp - qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))/deltat

! add to output cme
            cmeout(i,k) = cmeout(i,k) + dum

! now add to tendencies, partition between liquid and ice based on 
! temperature

         if( .not. lflag ) then
            if (tiedtke_macrophysics) then
              if (ttmp > tmelt - 5._r8) then
                dum1 = 0.0_r8
                ssat_disposal(i,k) = 1._r8
              else if (ttmp < tmelt - 40._r8) then
                dum1 = 1.0_r8
                ssat_disposal(i,k) = 2._r8
              else
                dum1 = 0.0_r8                   
                ssat_disposal(i,k) = 1._r8
              end if  
            else ! (tiedtke_macrophysics)

! for non-tiedtke, need to define how supersaturation removal affects
! particle number, if at all  ??
! for now, assume supersaturation removal does not lead to change in
! particle numbers
! use ice / liq partitioning as in original NCAR
              ssat_disposal(i,k) = 0._r8
              if (ttmp > tmelt - 5._r8) then
                dum1 = 0.0_r8
              else if (ttmp < tmelt - 35._r8) then
                dum1 = 1.0_r8
              else
                dum1 = (tmelt - 5._r8 - ttmp)/30._r8
              end if  
            endif !(tiedtke_macrophysics)

!????????
!RSH 10/05 0553
!           dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
            dum = (qtmp - qsn)/(1._r8 +   &
                               (xxlv*dum1 + xxlv*(1._r8 - dum1))**2* &
                                             qsn/(cpp*rv*ttmp**2))/deltat

         else
           if (ttmp > 268.15_r8) then
              dum1=0.0_r8
! now add to tendencies, partition between liquid and ice based on te
           else if (ttmp < 238.15_r8) then
              dum1=1.0_r8
           else
              dum1=(268.15_r8-ttmp)/30._r8
           end if  
           dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                     *qsn/(cpp*rv*ttmp**2))/deltat
         endif
            qctend(i,k) = qctend(i,k) + dum*(1._r8 - dum1)
! for output
            qcreso(i,k) = dum*(1._r8 - dum1)
            qitend(i,k) = qitend(i,k) + dum*dum1
            qireso(i,k) = dum*dum1
            qvlat(i,k) = qvlat(i,k) - dum
! for output
            qvres(i,k) = -dum
            tlat(i,k) = tlat(i,k) + dum*(1._r8 - dum1)*xxlv + dum*dum1*xxls
          else 
            if( .not. lflag ) &
            ssat_disposal(i,k) = 0._r8
          end if
#else
          if (qtmp > qsn .and. qsn > 0._r8) then

! expression below is approximate since there may be ice deposition
            dum = (qtmp - qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))/deltat

! add to output cme
            cmeout(i,k) = cmeout(i,k) + dum

! now add to tendencies, partition between liquid and ice based on 
! temperature
            if (ttmp > tmelt - 5._r8) then
              dum1 = 0.0_r8
            else if (ttmp < tmelt - 35._r8) then
              dum1 = 1.0_r8
            else
              dum1 = (tmelt - 5._r8 - ttmp)/30._r8
            end if  

            dum = (qtmp - qsn)/(1._r8 +   &
                                  (xxls*dum1 + xxlv*(1._r8 - dum1))**2 *  &
                                             qsn/(cpp*rv*ttmp**2))/deltat
            qctend(i,k) = qctend(i,k) + dum*(1._r8 - dum1)
! for output
            qcreso(i,k) = dum*(1._r8 - dum1)
            qitend(i,k) = qitend(i,k) + dum*dum1
            qireso(i,k) = dum*dum1
            qvlat(i,k) = qvlat(i,k) - dum
! for output
            qvres(i,k) = -dum
            tlat(i,k) = tlat(i,k) + dum*(1._r8 - dum1)*xxlv+dum*dum1*xxls
          endif
#endif
!.........................................................................
! calculate effective radius for pass to radiation code
! if no cloud water, default value is 10 micron for droplets,
! 25 micron for cloud ice

! update cloud variables after instantaneous processes to get effective 
! radius
! variables are in-cloud to calculate size dist parameters

          dumc(i,k) = max(qc(i,k) + qctend(i,k)*deltat, 0._r8)/lcldm(i,k)
          dumi(i,k) = max(qi(i,k) + qitend(i,k)*deltat, 0._r8)/icldm(i,k)
          dumnc(i,k) = max(nc(i,k) + nctend(i,k)*deltat, 0._r8)/lcldm(i,k)
          dumni(i,k) = max(ni(i,k) + nitend(i,k)*deltat, 0._r8)/icldm(i,k)

!-->cjg
! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
        dumnc(i,k)=ncnst/rho(i,k)
        end if

! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
        dumni(i,k)=ninst/rho(i,k)
        end if
!<--cjg

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

          dumc(i,k) = min(dumc(i,k), 5.e-3_r8)
          dumi(i,k) = min(dumi(i,k), 5.e-3_r8)

!...................
! cloud ice effective radius

          if (dumi(i,k) .ge. qsmall) then
#ifdef GFDL_COMPATIBLE_MICROP
            IF (diag_id%qnidt_size_adj +  diag_id%qni_size_adj_col > 0) &
                 diag_4d(i,j,k,diag_pt%qnidt_size_adj ) =  nitend(i,k)
#endif

! add upper limit to in-cloud number concentration to prevent numerical 
! error
            dumni(i,k) = min(dumni(i,k), dumi(i,k)*1.e20_r8)
            lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)

            if (lami(k) .lt. lammini) then
              lami(k) = lammini
              n0i(k) = lami(k)**(di + 1._r8)*dumi(i,k)/(ci*cons1)
              niic(i,k) = n0i(k)/lami(k)

! adjust number conc if needed to keep mean size in reasonable range
              nitend(i,k) = (niic(i,k)*icldm(i,k) - ni(i,k))/deltat
            else if (lami(k) .gt. lammaxi) then
              lami(k) = lammaxi
              n0i(k) = lami(k)**(di + 1._r8)*dumi(i,k)/(ci*cons1)
              niic(i,k) = n0i(k)/lami(k)

! adjust number conc if needed to keep mean size in reasonable range
              nitend(i,k) = (niic(i,k)*icldm(i,k) - ni(i,k))/deltat
            end if

#ifdef GFDL_COMPATIBLE_MICROP
            IF (diag_id%qnidt_size_adj +  diag_id%qni_size_adj_col > 0) &
                 diag_4d(i,j,k,diag_pt%qnidt_size_adj ) =     &
                     nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_size_adj ) 
#endif

            effi(i,k) = 1.5_r8/lami(k)*1.e6_r8
          else
            effi(i,k) = 25._r8
          end if


! cloud droplet effective radius

          if (dumc(i,k) .ge. qsmall) then
#ifdef GFDL_COMPATIBLE_MICROP
            IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col  > 0) &
                 diag_4d(i,j,k,diag_pt%qndt_size_adj ) =  nctend(i,k)
#endif

! add upper limit to in-cloud number concentration to prevent numerical 
! error
            dumnc(i,k) = min(dumnc(i,k), dumc(i,k)*1.e20_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set tendency to ensure minimum droplet concentration after update by 
! microphysics, except when lambda exceeds bounds on 
! mean drop size or if there is no cloud water
            if (dumnc(i,k) .lt. cdnl/rho(i,k)) then   
              nctend(i,k) = (cdnl/rho(i,k)*cldm(i,k) - nc(i,k))/deltat   
            end if
            dumnc(i,k) = max(dumnc(i,k), cdnl/rho(i,k)) ! sghan minimum 
                                                        ! in #/cm3 
!-->cjg
! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
! make sure nc is consistence with the constant N by adjusting tendency, need
! to multiply by cloud fraction
! note that nctend may be further adjusted below if mean droplet size is
! out of bounds

        nctend(i,k)=(ncnst/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat
        end if
!<--cjg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef GFDL_COMPATIBLE_MICROP
      if( .not. lflag ) then
!RSH76      pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
            pgam(k) = 0.0005714_r8*(dumnc(i,k)/1.e6_r8*rho(i,k)) +   &
                                                                 0.2714_r8
      else
           pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
      endif
#else
           pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
#endif
            pgam(k) = 1._r8/(pgam(k)**2) - 1._r8
            pgam(k) = max(pgam(k), 2._r8)
            pgam(k) = min(pgam(k), 15._r8)
           
            lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k) + 4._r8)/ &
                        (dumc(i,k)*gamma(pgam(k) + 1._r8)))**(1._r8/3._r8)
            lammin = (pgam(k) + 1._r8)/50.e-6_r8
            lammax = (pgam(k) + 1._r8)/2.e-6_r8
            if (lamc(k) .lt. lammin) then
              lamc(k) = lammin
              ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)*  &
                                        gamma(pgam(k) + 1._r8)/&
                                           (pi*rhow*gamma(pgam(k) + 4._r8))

! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k) = (ncic(i,k)*lcldm(i,k) - nc(i,k))/deltat

            else if (lamc(k) .gt. lammax) then
              lamc(k) = lammax
              ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                                      gamma(pgam(k) + 1._r8)/ &
                                          (pi*rhow*gamma(pgam(k) + 4._r8))

! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k) = (ncic(i,k)*lcldm(i,k) - nc(i,k))/deltat
            end if

#ifdef GFDL_COMPATIBLE_MICROP
            IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col  > 0) &
                 diag_4d(i,j,k,diag_pt%qndt_size_adj ) =     &
                        nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_size_adj )
#endif

            effc(i,k) = gamma(pgam(k) + 4._r8)/ &
                              gamma(pgam(k) + 3._r8)/lamc(k)/2._r8*1.e6_r8

!assign output fields for shape here
            lamcrad(i,k) = lamc(k)
            pgamrad(i,k) = pgam(k)
          else
            effc(i,k) = 10._r8
            lamcrad(i,k)=0._r8
            pgamrad(i,k)=0._r8
          end if

! ice effective diameter for david mitchell's optics
          deffi(i,k) = effi(i,k)*rhoi/917._r8*2._r8


! recalculate effective radius for constant number, in order to separate
! first and second indirect effects
! assume constant number of 10^8 kg-1

          dumnc(i,k) = 1.e8_r8

          if (dumc(i,k) .ge. qsmall) then
#ifdef GFDL_COMPATIBLE_MICROP
          if( .not. lflag ) then
!RSH76      pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
            pgam(k) = 0.0005714_r8*(dumnc(i,k)/1.e6_r8*rho(i,k)) +   &
                                                                 0.2714_r8
          else
            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
          endif
#else
            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
#endif
            pgam(k) = 1._r8/(pgam(k)**2) - 1._r8
            pgam(k) = max(pgam(k), 2._r8)
            pgam(k) = min(pgam(k), 15._r8)

            lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k) + 4._r8)/ &
                         (dumc(i,k)*gamma(pgam(k) + 1._r8)))**(1._r8/3._r8)
            lammin = (pgam(k) + 1._r8)/50.e-6_r8
            lammax = (pgam(k) + 1._r8)/2.e-6_r8
            if (lamc(k) .lt. lammin) then
              lamc(k) = lammin
            else if (lamc(k) .gt. lammax) then
              lamc(k) = lammax
            end if
            effc_fn(i,k) = gamma(pgam(k) + 4._r8)/ &
                             gamma(pgam(k) + 3._r8)/lamc(k)/2._r8*1.e6_r8
          else
            effc_fn(i,k) = 10._r8
          end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

        end do ! vertical k loop

 500    continue

        do k=1,pver
! if updated q (after microphysics) is zero, then ensure updated n is 
! also zero

#ifdef GFDL_COMPATIBLE_MICROP
          IF (diag_id%qndt_fill2  + diag_id%qn_fill2_col > 0) &
               diag_4d(i,j,k,diag_pt%qndt_fill2 ) =  nctend(i,k)
          IF (diag_id%qnidt_fill2 +  diag_id%qni_fill2_col > 0) &
               diag_4d(i,j,k,diag_pt%qnidt_fill2 ) =  nitend(i,k)
#endif
          if (qc(i,k) + qctend(i,k)*deltat .lt. qsmall) nctend(i,k) =  &
                                                            -nc(i,k)/deltat
          if (qi(i,k) + qitend(i,k)*deltat .lt. qsmall) nitend(i,k)=   &
                                                            -ni(i,k)/deltat
#ifdef GFDL_COMPATIBLE_MICROP
          IF (diag_id%qndt_fill2  + diag_id%qn_fill2_col > 0) &
               diag_4d(i,j,k,diag_pt%qndt_fill2 ) =     &
                          nctend(i,k) - diag_4d(i,j,k,diag_pt%qndt_fill2 ) 
          IF (diag_id%qnidt_fill2 +  diag_id%qni_fill2_col > 0) &
               diag_4d(i,j,k,diag_pt%qnidt_fill2 ) =     &
                         nitend(i,k) - diag_4d(i,j,k,diag_pt%qnidt_fill2 ) 
#endif
        end do ! k loop

      end do ! i loop

! hm add rain/snow mixing ratio and number concentration as diagnostic

#ifndef GFDL_COMPATIBLE_MICROP
      call outfld('QRAIN',qrout,   pcols, lchnk)
      call outfld('QSNOW',qsout,   pcols, lchnk)
      call outfld('NRAIN',nrout,   pcols, lchnk)
      call outfld('NSNOW',nsout,   pcols, lchnk)
#endif

! add snow output
      do i = 1,ncol
        do k=1,pver
          if (qsout(i,k) .gt. 1.e-7_r8 .and. nsout(i,k) .gt. 0._r8) then
            dsout(i,k) = 3._r8*rhosn/917._r8*    &
                           (pi*rhosn*nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
          endif
        end do
      end do

#ifndef GFDL_COMPATIBLE_MICROP
      call outfld('DSNOW',dsout,   pcols, lchnk)
 
! calculate effective radius of rain and snow in microns for COSP using 
! Eq. 9 of COSP v1.3 manual
      do i = 1,ncol
        do k=1,pver
!! RAIN
          if (qrout(i,k) .gt. 1.e-7_r8 .and. nrout(i,k) .gt. 0._r8) then
            reff_rain(i,k) = 1.5_r8*   &
                   (pi*rhow*nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8 
          endif
!! SNOW
          if (qsout(i,k) .gt. 1.e-7_r8 .and. nsout(i,k) .gt. 0._r8) then
            reff_snow(i,k) = 1.5_r8*    &
                   (pi*rhosn*nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8 
          endif
        end do
      end do

#else
! calculate effective radius of rain and snow in microns for COSP using 
! Eq. 9 of COSP v1.3 manual
! convert to diameter to pass out for use in radiation package
! snow_size is not currently used  -- as per mns
      do i = 1,ncol
        do k=1,pver
!! RAIN
          if (qrout(i,k) .gt. 1.e-7_r8 .and. nrout(i,k) .gt. 0._r8) then
            lsc_rain_size(i,k) = 3.0_r8*  &
                    (pi*rhow*nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8

! ---> h1g, hard-write rain effective radius range 30--750 um, 2012-04-25
            if( lflag ) then  
             lsc_rain_size(i,k) = max(  60.0_r8, lsc_rain_size(i,k) )
             lsc_rain_size(i,k) = min(1500.0_r8, lsc_rain_size(i,k) )
            endif
          else
            lsc_rain_size(i,k) = 100._r8
          endif
!! SNOW
         if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
           if( lflag ) &
           lsc_snow_size(i,k) = 3.0_r8*    &
                  (pi*rhosn*nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8
         endif
        end do
      end do
#endif

! analytic radar reflectivity
! formulas from Matthew Shupe, NOAA/CERES
! *****note: radar reflectivity is local (in-precip average)
! units of mm^6/m^3

      do i = 1,ncol
        do k=1,pver
          if (qc(i,k) + qctend(i,k)*deltat .ge. qsmall) then
            dum = ((qc(i,k) + qctend(i,k)*deltat)/lcldm(i,k)  &
                                              *rho(i,k)*1000._r8)**2 /  &
                  (0.109_r8*(nc(i,k) + nctend(i,k)*deltat)/lcldm(i,k)*  &
                                   rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
          else
            dum=0._r8
          end if
          if (qi(i,k) + qitend(i,k)*deltat .ge. qsmall) then
            dum1 = ((qi(i,k) + qitend(i,k)*deltat)*rho(i,k)/   &
                         icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)* &
                                                   icldm(i,k)/cldmax(i,k)
          else 
            dum1 = 0._r8
          end if
         
          if (qsout(i,k) .ge. qsmall) then
            dum1 = dum1 +   &
                     (qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
          end if
            
          refl(i,k) = dum + dum1
 
! add rain rate, but for 37 GHz formulation instead of 94 GHz
! formula approximated from data of Matrasov (2007)
! rainrt is the rain rate in mm/hr
! reflectivity (dum) is in DBz
! don't include rain rate in R calculation for values less than 0.001 mm/hr

          if (rainrt(i,k) .ge. 0.001_r8) then
            dum = log10(rainrt(i,k)**6._r8) + 16._r8
 
! convert from DBz to mm^6/m^3
            dum = 10._r8**(dum/10._r8)
          else
            dum=0._r8
          end if
 
! add to refl
 
          refl(i,k) = refl(i,k) + dum
 
! output reflectivity in Z.
          areflz(i,k) = refl(i,k)
 
! convert back to DBz 
 
          if (refl(i,k) .gt. minrefl) then 
            refl(i,k) = 10._r8*log10(refl(i,k))
          else
            refl(i,k) = -9999._r8
          end if
  
! set averaging flag
          if (refl(i,k).gt.mindbz) then 
            arefl(i,k) = refl(i,k)
            frefl(i,k) = 1.0_r8  
          else
            arefl(i,k) = 0._r8
            areflz(i,k) = 0._r8
            frefl(i,k) = 0._r8
          end if
 
! bound cloudsat reflectivity
 
          csrfl(i,k) = min(csmax, refl(i,k))
 
! set averaging flag
          if (csrfl(i,k) .gt. csmin) then 
            acsrfl(i,k) = refl(i,k)
            fcsrfl(i,k) = 1.0_r8  
          else
            acsrfl(i,k) = 0._r8
            fcsrfl(i,k) = 0._r8
          end if
 
        end do
      end do

#ifndef GFDL_COMPATIBLE_MICROP
      call outfld('REFL',refl,   pcols, lchnk)
      call outfld('AREFL',arefl,   pcols, lchnk)
      call outfld('AREFLZ',areflz,   pcols, lchnk)
      call outfld('FREFL',frefl,   pcols, lchnk)
      call outfld('CSRFL',csrfl,   pcols, lchnk)
      call outfld('ACSRFL',acsrfl,   pcols, lchnk)
      call outfld('FCSRFL',fcsrfl,   pcols, lchnk)

      call outfld('RERCLD',rercld,   pcols, lchnk)
#endif

! averaging for snow and rain number and diameter

      qrout2(:,:) = 0._r8
      qsout2(:,:) = 0._r8
      nrout2(:,:) = 0._r8
      nsout2(:,:) = 0._r8
      drout2(:,:) = 0._r8
      dsout2(:,:) = 0._r8
      freqs(:,:) = 0._r8
      freqr(:,:) = 0._r8
      do i = 1,ncol
        do k=1,pver
          if (qrout(i,k) .gt. 1.e-7_r8 .and. nrout(i,k) .gt. 0._r8) then
            qrout2(i,k) = qrout(i,k)
            nrout2(i,k) = nrout(i,k)
            drout2(i,k) = (pi*rhow*nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
            freqr(i,k) = 1._r8
          endif
          if (qsout(i,k) .gt. 1.e-7_r8 .and. nsout(i,k) .gt. 0._r8) then
            qsout2(i,k) = qsout(i,k)
            nsout2(i,k) = nsout(i,k)
            dsout2(i,k) = (pi*rhosn*nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
            freqs(i,k) = 1._r8
          endif
        end do
      end do

! output activated liquid and ice (convert from #/kg -> #/m3)
      do i = 1,ncol
        do k=1,pver
          ncai(i,k) = dum2i(i,k)*rho(i,k)
          ncal(i,k) = dum2l(i,k)*rho(i,k)
        end do
      end do

#ifndef GFDL_COMPATIBLE_MICROP
      call outfld('NCAL',ncal,    pcols,lchnk)
      call outfld('NCAI',ncai,    pcols,lchnk)

!add averaged output fields.
      call outfld('AQRAIN',qrout2,    pcols,lchnk)
      call outfld('AQSNOW',qsout2,    pcols,lchnk)
      call outfld('ANRAIN',nrout2,    pcols,lchnk)
      call outfld('ANSNOW',nsout2,    pcols,lchnk)
      call outfld('ADRAIN',drout2,    pcols,lchnk)
      call outfld('ADSNOW',dsout2,    pcols,lchnk)
      call outfld('FREQR',freqr,    pcols,lchnk)
      call outfld('FREQS',freqs,    pcols,lchnk)
#endif

!redefine fice here....
      nfice(:,:) = 0._r8
      do k=1,pver
        do i=1,ncol
          dumc(i,k) = (qc(i,k) + qctend(i,k)*deltat)
          dumi(i,k) = (qi(i,k) + qitend(i,k)*deltat)
          dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)  
          if (dumfice .gt. qsmall .and.    &
                             (qsout(i,k) + dumi(i,k) .gt. qsmall)) then
            nfice(i,k) = (qsout(i,k) + dumi(i,k))/dumfice
          endif

          if (nfice(i,k) .gt. 1._r8) then
            nfice(i,k) = 1._r8
          endif

        enddo  ! i loop
      enddo  ! k loop

#ifndef GFDL_COMPATIBLE_MICROP
      call outfld('FICE',nfice,   pcols, lchnk)
#endif

#ifdef GFDL_COMPATIBLE_MICROP
! diagnostics for water tendencies
! water  vapor specific humicity
      if  (diag_id%rain_evap + diag_id%rain_evap_col > 0)  &
              diag_4d(:,j,:,diag_pt%rain_evap)  = -preo( : , : )
      if  (diag_id%qdt_rain_evap > 0)  &
              diag_4d(:,j,:,diag_pt%qdt_rain_evap)  = -preo( : , : )
      if (diag_id%qdt_cond   > 0) &
              diag_4d(:,j,:,diag_pt%qdt_cond)  = -cmelo(:,:)
      if  (diag_id%qdt_snow_sublim > 0 .or.   &
                                     diag_id%q_snow_sublim_col > 0 )  &
              diag_4d(:,j,:,diag_pt%qdt_snow_sublim )  =  - prdso( : , : )
      if  (diag_id%qdt_deposition > 0)  &
              diag_4d(:,j,:,diag_pt%qdt_deposition )  = -cmeiout( : , : )
      if  (diag_id%qdt_sedi_ice2vapor> 0)  &
              diag_4d(:,j,:,diag_pt%qdt_sedi_ice2vapor) = qisevap( : , : )
      if  (diag_id%qdt_sedi_liquid2vapor> 0)  &
            diag_4d(:,j,:,diag_pt%qdt_sedi_liquid2vapor) = qcsevap( : , : )
      if  (diag_id%qdt_super_sat_rm > 0)  &
              diag_4d(:,j,:,diag_pt%qdt_super_sat_rm) = qvres( : , : )

! cloud liquid water
      if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_accr)  = - prao(:,:)
      if (diag_id%qldt_auto  + diag_id%ql_auto_col > 0)&
              diag_4d(:,j,:,diag_pt%qldt_auto)  = -prco(:,:)
      if (diag_id%qldt_freez2 + diag_id%ql_freez2_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_freez2) =   &
                                          -(mnuccco(:,:) + mnuccto(:,:) )
      sum_freeze2(:,:) =  mnuccco(:,:) + mnuccto(:,:)
      if (diag_id%qldt_accrs  + diag_id%ql_accrs_col > 0) & 
              diag_4d(:,j,:,diag_pt%qldt_accrs)  = -psacwso(:,:) 
      sum_rime(:,:) =  psacwso(:,:)
      if (diag_id%qldt_HM_splinter + diag_id%ql_HM_splinter_col > 0)&
              diag_4d(:,j,:,diag_pt%qldt_HM_splinter)  = -msacwio(:,:)
      sum_splinter(:,:) =  msacwio(:,:)
      if (diag_id%qldt_bergs + diag_id%ql_bergs_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_bergs)  = -bergso(:,:)
      sum_bergs(:,:) =  bergso (:,:)
      if (diag_id%qidt_dep + diag_id%qi_dep_col > 0)    &
              diag_4d(:,j,:,diag_pt%qidt_dep)  = max(cmeiout(:,:),0._r8) 
      sum_cond(:,:) = max(cmeiout(:,:),0._r8) 
      if (diag_id%qidt_subl + diag_id%qi_subl_col > 0)  &
              diag_4d(:,j,:,diag_pt%qidt_subl)  =     &
                                         -max(-1._r8*cmeiout(:,:),0._r8) 
      if (diag_id%qldt_cond + diag_id%ql_cond_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_cond)  =  max(cmelo(:,:), 0._r8)
      if (diag_id%qldt_evap  + diag_id%ql_evap_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_evap)  =   &
                                           - max(-1._r8*cmelo(:,:),0._r8)
      if (diag_id%qldt_eros + diag_id%ql_eros_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_eros)  =   eroslo   (:,:)
      if (diag_id%qdt_eros_l                       > 0) &
              diag_4d(:,j,:,diag_pt%qdt_eros_l)  =  -eroslo   (:,:)
      if (diag_id%qidt_eros + diag_id%qi_eros_col > 0) &
              diag_4d(:,j,:,diag_pt%qidt_eros)  =  erosio   (:,:)
      if (diag_id%qdt_eros_i                       > 0) &
              diag_4d(:,j,:,diag_pt%qdt_eros_i)  = -erosio   (:,:)
      if (diag_id%qldt_berg + diag_id%ql_berg_col > 0) &
              diag_4d(:,j,:,diag_pt%qldt_berg)  =  -bergo(:,:)
      sum_berg(:,:) =  bergo(:,:)
      IF ( diag_id%qldt_sedi  + diag_id%ql_sedi_col > 0 ) &
              diag_4d(:,j,1:pver,diag_pt%qldt_sedi) = qcsedten(:,1:pver)
      IF ( diag_id%liq_adj  + diag_id%liq_adj_col > 0 ) &
              diag_4d(:,j,:,diag_pt%liq_adj) = qcreso(:,:)

! cloud ice water
      if (diag_id%qidt_auto + diag_id%qi_auto_col > 0) &
             diag_4d(:,j,:,diag_pt%qidt_auto) = -prcio(:,:)
      if (diag_id%qidt_accr  + diag_id%qi_accr_col > 0) &
             diag_4d(:,j,:,diag_pt%qidt_accr) = -praio(:,:)
      IF ( diag_id%qidt_fall  + diag_id%qi_fall_col > 0 ) &
              diag_4d(:,j,1:pver,diag_pt%qidt_fall) = qisedten(:,1:pver)
      IF ( diag_id%ice_adj  +  diag_id%ice_adj_col > 0 ) &
              diag_4d(:,j,:,diag_pt%ice_adj) = qireso(:,:)
      sum_ice_adj(:,:) = qireso(:,:)

! ---> rain water mixing ratio
      if (diag_id%srfrain_accrs + diag_id%srfrain_accrs_col > 0)    &
             diag_4d(:,j,:,diag_pt%srfrain_accrs)  = -pracso(:,:)
      if (diag_id%srfrain_freez + diag_id%srfrain_freez_col > 0)    &
             diag_4d(:,j,:,diag_pt%srfrain_freez)  = -mnuccro(:,:)

! ---> snow mixing ratio

! --->liquid droplet number
      if (diag_id%qndt_cond + diag_id%qn_cond_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_cond)  = npccno(:,:)  /real(iter)
      if (diag_id%qndt_freez + diag_id%qn_freez_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_freez)  = nnuccco(:,:)  /real(iter)
      if (diag_id%qndt_contact_frz > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_contact_frz)  =    &
                                                 nnuccto(:,:) /real(iter)
      if (diag_id%qndt_sacws + diag_id%qn_sacws_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_sacws)  = npsacwso(:,:) /real(iter)
      if (diag_id%qndt_evap + diag_id%qn_evap_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_evap)  = nsubco(:,:)  /real(iter)
      if (diag_id%qndt_eros + diag_id%qn_eros_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_eros)  = nerosco(:,:)  /real(iter)
      if (diag_id%qndt_pra + diag_id%qn_pra_col > 0)    &
             diag_4d(:,j,:,diag_pt%qndt_pra)  = nprao(:,:)  /real(iter)
      if (diag_id%qndt_auto + diag_id%qn_auto_col > 0)    &
              diag_4d(:,j,:,diag_pt%qndt_auto)  = nprc1o(:,:)  /real(iter)
      if ( diag_id%qndt_nucclim  + diag_id%qn_nucclim_col  > 0 ) &
            diag_4d(:,j,:,diag_pt%qndt_nucclim)  =     &
                                                nucclimo(:,:)  /real(iter)

! ---> ice number 
      if (diag_id%qnidt_nnuccd +  diag_id%qni_nnuccd_col > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nnuccd)  =    &
                                                 nnuccdo(:,:)  /real(iter)
      if (diag_id%qnidt_nsacwi> 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nsacwi)  =     &
                                                 nsacwio(:,:)  /real(iter)
      if (diag_id%qnidt_nsubi  + diag_id%qni_nsubi_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nsubi)  = nsubio(:,:)  /real(iter)
      if (diag_id%qnidt_nerosi  + diag_id%qni_nerosi_col  > 0)    &
           diag_4d(:,j,:,diag_pt%qnidt_nerosi)  = nerosio(:,:)  /real(iter)
      if (diag_id%qnidt_nprci  + diag_id%qni_nprci_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nprci)  = nprcio(:,:)  /real(iter)
      if (diag_id%qnidt_nprai  + diag_id%qni_nprai_col  > 0)    &
             diag_4d(:,j,:,diag_pt%qnidt_nprai)  = npraio(:,:)  /real(iter)
      if (diag_id%qnidt_nucclim1 +  diag_id%qni_nucclim1_col > 0 ) &
             diag_4d(:,j,:,diag_pt%qnidt_nucclim1) =      &
                                              nucclim1io(:,:)  /real(iter)
!RSH:
!   calculate fraction of total ice / snow creation that requires 
!   ice-forming nuclei
      do k=1,pver
        do i=1,ncol
          qldt_sum = sum_cond(i,k) + sum_rime(i,k) + sum_berg(i,k) + &
                     sum_ice_adj(i,k) + MAX(sum_bergs(i,k), 0.0) + &
                     sum_freeze(i,k) + sum_freeze2(i,k) + sum_splinter(i,k)
          if (ABS(qldt_sum) > 0.0            ) then
            f_snow_berg(i,k) = (sum_berg(i,k) + sum_cond(i,k) +   &
                                sum_ice_adj(i,k) +    &
                                MAX( sum_bergs(i,k), 0.0) +     &
                                sum_freeze (i,k))/qldt_sum        
          else
            f_snow_berg(i,k) = 0._r8
          endif
        end do
      end do

#endif

      return

end subroutine mmicro_pcond

!########################################################################

subroutine mmicro_end

return


end subroutine mmicro_end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION GAMMA(X)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!D    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY

      real(r8) gamma
      REAL(r8) &
!D    DOUBLE PRECISION
         C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
          SQRTPI/0.9189385332046727417803297E0_r8/, &
          PI/3.1415926535897932384626434E0_r8/
!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/, &
          XINF/3.4E38_r8/
!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,&
            -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,&
            8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,&
            -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/
      DATA Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
           -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
             2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
           -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
          -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
          -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
           5.7083835261E-03_r8/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I,r8)
!D    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
!D900 DGAMMA = RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END function gamma



#ifdef GFDL_COMPATIBLE_MICROP
!########################################################################
      function polysvp (T,type)
!  Compute saturation vapor pressure by using
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      real(r8) dum

      real(r8) T,polysvp

      integer type

! ice

      if (type.eq.1) then

! Goff Gatch equation (good down to -100 C)

         polysvp = 10._r8**(-9.09718_r8*(273.16_r8/t-1._r8)-3.56654_r8* &
          log10(273.16_r8/t)+0.876793_r8*(1._r8-t/273.16_r8)+ &
          log10(6.1071_r8))*100._r8

      end if

! Goff Gatch equation, uncertain below -70 C

      if (type.eq.0) then
         polysvp = 10._r8**(-7.90298_r8*(373.16_r8/t-1._r8)+ &
             5.02808_r8*log10(373.16_r8/t)- &
             1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/373.16_r8))-1._r8)+ &
             8.1328e-3_r8*(10._r8**(-3.49149_r8*(373.16_r8/t-1._r8))-1._r8)+ &
             log10(1013.246_r8))*100._r8
         end if


      end function polysvp
!#########################################################################


!#########################################################################
subroutine vqsatd_water(t       ,p       ,es      ,qs      ,gam      , &
                        len     )

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: len       ! vector length
   real(r8), intent(in) :: t(len)       ! temperature
   real(r8), intent(in) :: p(len)       ! pressure

!
! Output arguments
!
   real(r8), intent(out) :: es(len)   ! saturation vapor pressure
   real(r8), intent(out) :: qs(len)   ! saturation specific humidity
   real(r8), intent(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
!
!--------------------------Local Variables------------------------------
!
   integer i      ! index for vector calculations
!
   real(r8) omeps     ! 1. - 0.622
   real(r8) hltalt    ! appropriately modified hlat for T derivatives
!
   real(r8) hlatsb    ! hlat weighted in transition region
   real(r8) hlatvp    ! hlat modified for t changes above freezing
   real(r8) desdt     ! d(es)/dT
!
!-----------------------------------------------------------------------
!
   omeps = 1.0_r8 - epsqs
   do i=1,len
      es(i) = polysvp(t(i),0)
!
! Saturation specific humidity
!
      qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
      qs(i) = min(1.0_r8,qs(i))
!
      if (qs(i) < 0.0_r8) then
         qs(i) = 1.0_r8
         es(i) = p(i)
      end if
   end do
!
! No icephs or water to ice transition
!
   do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
      hlatvp = hlv - 2369.0_r8*(t(i)-tmelt)
      hlatsb = hlv
      if (t(i) < tmelt) then
         hltalt = hlatsb
      else
         hltalt = hlatvp
      end if
      desdt  = hltalt*es(i)/(rvgas*t(i)*t(i))
      gam(i) = hltalt*qs(i)*p(i)*desdt/(cpp*es(i)*(p(i) - omeps*es(i)))
      if (qs(i) == 1.0_r8) gam(i) = 0.0_r8
   end do
!
   return
end subroutine vqsatd_water

!#########################################################################
end module cldwat2m_micro_mod

#else

end module cldwat2m_micro

#endif

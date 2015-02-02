                 module microphys_rad_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  smf
! </REVIEWER>
! <OVERVIEW>
!  Code to provide micro physics subroutines for radiation calculation
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!
!  shared modules:

use mpp_mod,               only:  input_nml_file
use fms_mod,               only:  fms_init, open_namelist_file, &
                                  write_version_number, mpp_pe, &
                                  mpp_root_pe, stdlog, file_exist,  &
                                  check_nml_error, error_mesg,   &
                                  FATAL, close_file
use constants_mod,         only:  constants_init, diffac, radian
use time_manager_mod,      only:  time_type
use diag_manager_mod,      only:  register_diag_field, send_data, &
                                  diag_manager_init
use random_numbers_mod,    only:  randomNumberStream,   &
                                  initializeRandomNumberStream, &
                                  getRandomNumbers,             &
                                  constructSeed
use cloud_generator_mod,   only:  cloud_generator_init, &
                                  cloud_generator_end

!  shared radiation package modules:

use rad_utilities_mod,     only:  rad_utilities_init, Lw_control, &
                                  Cldrad_control, thickavg,  &
                                  cld_specification_type, &
                                  microrad_properties_type, &
                                  cldrad_properties_type, &
                                  microphysics_type, Lw_parameters
use longwave_params_mod,   only:  NBLW, longwave_params_init
use esfsw_parameters_mod,  only:  esfsw_parameters_init, Solar_spect

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    microphys_rad_mod produces cloud radiative properties 
!    based upon input microphysical properties.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: microphys_rad.F90,v 20.0 2013/12/13 23:20:12 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
           microphys_rad_init, microphys_sw_driver,     &
           microphys_lw_driver, lwemiss_calc,    &
           comb_cldprops_calc, microphys_rad_end, &
           isccp_microphys_lw_driver, isccp_microphys_sw_driver

private         &
!    called from microphys_sw_driver:
           cloudpar, &
!    called from cloudpar:
           slingo, savijarvi, fu, icesolar, snowsw,  &
!    called from microphys_lw_driver:
           cloud_lwpar, cloud_lwem_oneband,    &
!    caled from cloud_lwpar:
           el, el_dge, cliqlw, furainlw, fusnowlw, &
!    currently not used:
           snowlw


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)   :: lwem_form=' '     ! longwave emissivity param-
                                         ! eterization; either 'fuliou'
                                         ! or 'ebertcurry'
logical       ::  do_orig_donner_stoch = .true.
logical       ::  do_delta_adj = .false.
logical       ::  do_const_asy = .false.
real          ::  val_const_asy = 0.75
real          ::  alpha = 0.1
                  ! frequency-independent parameter for absorption due 
                  ! to cloud drops in the infrared. this value is given 
                  ! in held et al, JAS, 1993. [ m**2 / g ]
logical :: ignore_donner_cells = .false.
                  ! when set to .true., the effects of donner cell clouds 
                  ! in the radiation code are ignored

namelist /microphys_rad_nml /     &
                               lwem_form, &
                               do_orig_donner_stoch, &
                               do_delta_adj, &
                               do_const_asy, val_const_asy, &
                               alpha, ignore_donner_cells

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!----------------------------------------------------------------------
!    parameters and data needed when frequency-dependent longwave 
!    cloud emissivities are activated.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    N_EMISS_BDS : number of  infrared frequency bands over which 
!                  frequency-dependent emissivities are computed 
!                  using appropriate parameterizations. 
!----------------------------------------------------------------------
integer, parameter       :: N_EMISS_BDS = 7

!----------------------------------------------------------------------
!    wavenumber limits for the cloud emissivity frequency bands.
!    THESE MAY BE CHANGED ONLY BY THE KEEPER OF THE RADIATION CODE.
!
!    cldbandlo : low wave number  boundary for the emissivity band
!    cldbandhi : high wave number  boundary for the emissivity band
!----------------------------------------------------------------------
real, dimension (N_EMISS_BDS)     :: cldbandlo, cldbandhi

data cldbandlo /                    &
             0.0, 560.0, 800.0, 900.0, 990.0, 1070.0, 1200.0 /
data cldbandhi /                    &
           560.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /

!----------------------------------------------------------------------
!    note: the cloud properties for wavelengths beyond 1400 wavenumbers
!    are included in the results for the first band, ie, that band
!    actually is 0-560, 1400-2200 cm-1. thus the following indices 
!    include an extra band.
!  
!     istartcldband : starting wave number index for emissivity band
!     iendcldband   : ending wave number index for emissivity band
!----------------------------------------------------------------------
integer, dimension (N_EMISS_BDS + 1) :: istartcldband, iendcldband

data istartcldband /                &
               1,   57,   81,    91,   100,    108,   121,   141 /
data iendcldband /                  &
              56,   80,   90,    99,   107,    120,   140,   220 /

!----------------------------------------------------------------------
!    parameters for Fu and Liou lw snow water parameterization.
!    NBFL : number of frequency bands in the lw snow water para-
!           meterization. corresponds to bands 7-18 in Fu and Liou 
!           (Table 2).
!    NBA  : number of terms in parameterization for series
!           expansion (not counting the 0th power term) for ai
!    NBB  : number of terms in parameterization for series
!           expansion (not counting the 0th power term) for bi
!    NBC  : number of terms in parameterization for series
!           expansion (not counting the 0th power term) for ci
!   NBA2
!   NBB2
!   NBC2
!----------------------------------------------------------------------
integer, parameter  :: NBFL= 12
integer, parameter  :: NBA = 3
integer, parameter  :: NBB = 3
integer, parameter  :: NBC = 3
integer, parameter  :: NBD = 3
 
integer, parameter  :: NBA2 = 5
integer, parameter  :: NBB2 = 5
integer, parameter  :: NBC2 = 5

!----------------------------------------------------------------------
!    wavenumber ranges  for Fu-Liou ice crystal parameterizations
!    these apply to the ice crystal (El), cloud rain (Furainlw)
!    and cloud snow (Fusnowlw) parameterizations. note: the cloud 
!    liquid drop parameterization (Cliqlw) is frequency-independent.
!
!    endfubands : high wavenumber boundary of wavenumber bands used
!                 in Fu-Liou parameterization. since model wavenumber 
!                 bands are in order of increasing wavenumber, the Fu 
!                 coefficients have been reversed; thus bands 1 -> 12
!                 correspond to Fu bands 18 -> 7.
!    iendfubands : index of model 10 cm-1 bands corresponding to
!                  the value of endfubands. computed in the code.
!----------------------------------------------------------------------
real,    dimension (NBFL)        :: endfubands
integer, dimension (NBFL)        :: iendfubands

data endfubands /                  &
                 280,    400,    540,   670,  800,  980,  1100,  1250, &
                1400,   1700,   1900,   2200 /

!----------------------------------------------------------------------
!    weighting factors relating fu lw bands to model lw frequency
!    bands.
!
!    fulwwts         : fraction of total planck function in emissivity
!                      band n that is in fu band ni. for a given emis-
!                      sivity band, the sum of fulwwts over all fu bands
!                      equals 1.0
!    planckivlicecld : value of the planck function in the portion of
!                      the spectrum common to emissivity band n and 
!                      fu band ni. 
!----------------------------------------------------------------------
real,    dimension (N_EMISS_BDS, NBFL)     :: fulwwts
real,    dimension (N_EMISS_BDS + 1, NBFL) :: planckivlicecld

!---------------------------------------------------------------------
!    nivl1lwicecld  :  fu band index corresponding to the lowest wave
!                      number of the emissivity band
!    nivl2lwicecld  :  fu band index corresponding to the highest wave
!                      number of the emissivity band
!    planckcldband  :  sum of the planck function over the given emis-
!                      sivity band
!---------------------------------------------------------------------
integer, dimension (N_EMISS_BDS + 1) :: nivl1lwicecld, nivl2lwicecld
real,    dimension (N_EMISS_BDS)     :: planckcldband

!----------------------------------------------------------------------
!    parameters for shortwave cloud-radiation parameterizations.
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    NICECLDIVLS  : the number of scattering spectral intervals for the
!                   ice crystals used in the fu parameterization (1996).
!    NICESOLARCLDIVLS 
!                 : the number of scattering spectral intervals for the
!                   ice crystals used in the icesolar parameterization
!                   (fu and liou, 1993).
!    NLIQCLDIVLS  : the number of scattering spectral intervals for the 
!                   cloud drops                                    
!    NRAINCLDIVLS : the number of scattering spectral intervals for the 
!                   rain drops                                         
!    NSNOWCLDIVLS : the number of scattering spectral intervals for snow
!----------------------------------------------------------------------
integer, parameter         ::   NICECLDIVLS      = 25 
integer, parameter         ::   NICESOLARCLDIVLS = 6
integer, parameter         ::   NLIQCLDIVLS      = 24 
integer, parameter         ::   NRAINCLDIVLS     = 4 
integer, parameter         ::   NSNOWCLDIVLS     = 6 
 
!---------------------------------------------------------------------
!    define the spectral limits for drop, rain, ice and snow single    
!    scattering properties in shortwave frequency ranges. 
!---------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    wavenumber limits for slingo cloud drop intervals.             
!---------------------------------------------------------------------
integer, dimension (NLIQCLDIVLS)       :: endliqcldwvn
 
data endliqcldwvn /  2924,  3437,  4202,  4695,  6098,  6536,  7813,   &
                     8404,  9091, 10000, 11494, 12821, 13333, 14493,   &
                    15625, 17544, 19231, 20833, 22727, 25000, 27778,   &
                    30303, 33333, 57600 /
 
!-------------------------------------------------------------------
!    wavenumber limits for Savijarvi rain drop intervals.
!-------------------------------------------------------------------
integer, dimension (NRAINCLDIVLS)      :: endraincldwvn
 
data endraincldwvn / 4202, 8403, 14493, 57600 /
 
!----------------------------------------------------------------------
!    wavenumber limits for icesolar ice crystal intervals.
!---------------------------------------------------------------------- 
integer, dimension (NICESOLARCLDIVLS)  :: endicesolcldwvn
 
data endicesolcldwvn / 2857, 4000, 5263, 7692, 14493, 57600 /
 
!---------------------------------------------------------------------
!    wavenumber limits for fu ice crystal intervals.
!--------------------------------------------------------------------
integer, dimension (NICECLDIVLS)       :: endicecldwvn
 
data endicecldwvn /  2000,  2924,  3437,  4202,  4695,  6098,  6536,   &
                     7092,  8404,  9091, 10000, 11494, 12821, 13333,   &
                    14493, 15625, 17544, 19231, 20833, 22727, 25000,   &
                    27778, 30303, 33333, 57600 /
 
!---------------------------------------------------------------------
!    wavenumber limits for the Fu snow intervals                       
!--------------------------------------------------------------------
integer, dimension (NSNOWCLDIVLS)      :: endsnowcldwvn
 
data endsnowcldwvn / 2857, 4000, 5263, 7692, 14493, 57600 /
 
!----------------------------------------------------------------------
!    these arrays define the intersection points of the solar spectral
!    bands and the wavenumber bands for each of the microphysical
!    species. these must be allocated since the number of spectral
!    bands is determined from namelist input.
!
!       nivl1liqcld    :  cloud droplet band index corresponding to the
!                         lowest wave number of spectral band n
!       nivl1icecld    :  ice crystal band index corresponding to the
!                         lowest wave number of spectral band n, 
!                         (fu, 1996)
!       nivl1icesolcld :  ice crystal band index corresponding to the
!                         lowest wave number of spectral band n
!                         (fu and liou, 1993)
!       nivl1raincld   :  rain drop band index corresponding to the
!                         lowest wave number of spectral band n
!       nivl1snowcld   :  snow flake band index corresponding to the
!                         lowest wave number of spectral band n
!       nivl2liqcld    :  cloud droplet band index corresponding to the
!                         highest wave number of spectral band n
!       nivl2icecld    :  ice crystal band index corresponding to the
!                         highest wave number of spectral band n, 
!                         (fu, 1996)
!       nivl2icesolcld :  ice crystal band index corresponding to the
!                         highest wave number of spectral band n
!                         (fu and liou, 1993)
!       nivl2raincld   :  rain drop band index corresponding to the
!                         highest wave number of spectral band n
!       nivl2snowcld   :  snow flake band index corresponding to the
!                         highest wave number of spectral band n
!----------------------------------------------------------------------
integer, dimension(:), allocatable  :: nivl1liqcld,   &
                                       nivl1icecld,   &
                                       nivl1icesolcld,   &
                                       nivl1raincld,  &
                                       nivl1snowcld,  &
                                       nivl2liqcld,   &
                                       nivl2icecld,   &
                                       nivl2icesolcld,   &
                                       nivl2raincld,  &
                                       nivl2snowcld 

!---------------------------------------------------------------------
!    these arrays define the sum of the toa solar flux in the wave
!    number spectrum common to solar spectral band n and particle
!    spectral band ni.
!
!    solivlicecld    : solar flux in solar spectral band n and ice
!                      crystal band ni (fu, 1996)
!    solivlicesolcld : solar flux in solar spectral band n and ice
!                      crystal band ni (fu and liou, 1996)
!    solivlliqcld    : solar flux in solar spectral band n and cloud
!                      droplet band ni (fu, 1996)
!    solivlraincld   : solar flux in solar spectral band n and rain
!                      drop band ni (fu, 1996)
!    solivlsnowcld    : solar flux in solar spectral band n and snow
!                      flake band ni (fu, 1996)
!---------------------------------------------------------------------
real, dimension(:,:), allocatable  :: solivlicecld,   &
                                      solivlicesolcld, &
                                      solivlliqcld, &
                                      solivlraincld , &
                                      solivlsnowcld

real                               :: min_cld_drop_rad, max_cld_drop_rad
real                               :: min_cld_ice_size, max_cld_ice_size

!----------------------------------------------------------------------
!   variables needed for random number seed:
!----------------------------------------------------------------------
real, dimension(:,:), allocatable  :: lats, lons ! lat and lon of columns
                                               ! in this processor's
                                               ! domain [ degrees ]

!----------------------------------------------------------------------
!    diagnostics variables
!----------------------------------------------------------------------
character(len=16)   :: mod_name = 'microphys_rad'
real                :: missing_value = -999.

!-------------------------------------------------------------------
!    logical variables:
!-------------------------------------------------------------------
logical    :: module_is_initialized = .false. ! module is initialized ?

!--------------------------------------------------------------------
!--------------------------------------------------------------------



                    contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="microphys_rad_init">
!  <OVERVIEW>
!   The microphys_rad module constructor
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine initializes micro physics module data,
!   determines micro physics parameterization scheme based
!   on initialization input data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_rad_init(cldhm_abs_in, cldml_abs_in)
!  </TEMPLATE>
!  <IN NAME="cldhm_abs_in" TYPE="cldhm_abs_in">
!   boundaries in sigma pressure level between high and middle
!   clouds
!  </IN>
!  <IN NAME="cldml_abs_in" TYPE="cldml_abs_in">
!   boundaries in sigma pressure level between middle and low
!   clouds
!  </IN>
! </SUBROUTINE>
!
subroutine microphys_rad_init (min_cld_drop_rad_in, max_cld_drop_rad_in, &
                               min_cld_ice_size_in, max_cld_ice_size_in, &
                               axes, Time, lonb, latb) 

!------------------------------------------------------------------
!    subroutine microphys_rad_init is the constructor for 
!    microphys_rad_mod.
!--------------------------------------------------------------------

real,                    intent(in)    :: min_cld_drop_rad_in, &
                                          max_cld_drop_rad_in, &
                                          min_cld_ice_size_in, &
                                          max_cld_ice_size_in     
integer, dimension(4),   intent(in)    ::  axes
type(time_type),         intent(in)    ::  Time
real, dimension(:,:),    intent(in)    ::  lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]

!----------------------------------------------------------------------
! local variables:                                                  

      real, dimension (NBLW)    ::  src1nb
      real                      :: c1, centnb, sc, x, x1
      real                      :: sumsol1, sumsol2, sumsol3, sumsol4, &
                                   sumsol5
      real                      :: sumplanck 
      real                      :: xtemv = 233.15
      integer                   :: unit, ierr, io, logunit
      integer                   :: nivl, nband
      integer                   :: nivl1, nivl2, nivl3, nivl4, nivl5
      integer                   :: n, ib, nw, ni

!---------------------------------------------------------------------
! local variables:                                                  
!
!     src1nb            radiation emitted in a wavelength band  
!                       [ sec (-3) or ergs/(sec*cm**2) ]
!     c1                expression in Planck's radiation law, 
!                       2*pi*speed of light squared* planck's constant/
!                       wavelength**3 [ gm cm / (sec**3) ]
!     centb             wave number at center of spectral interval
!                       [ cm (-1) ]
!     sc                radiation emitted per unit wavelength
!                       [ 1/(sec**3*cm) or ergs/(sec*cm**2) ]
!     x                 expression in planck's law: 
!                       (h*c)/(k*lambda*temp) [ nondimensional ]
!     x1                expression in planck's law:
!                       exp ( (h*c)/(k*lambda*temp) ) [ nondimensional ]
!     sumsol1           scalar used to accumulate sum of toa solar flux
!                       over a spectral interval common to a cloud drop-
!                       let band and a solar spectral band 
!     sumsol2           scalar used to accumulate sum of toa solar flux
!                       over a spectral interval common to a fu (1996)
!                       ice crystal band and a solar spectral band 
!     sumsol3           scalar used to accumulate sum of toa solar flux
!                       over a spectral interval common to a rain drop 
!                       band and a solar spectral band 
!     sumsol4           scalar used to accumulate sum of toa solar flux
!                       over a spectral interval common to a snow flake
!                       band and a solar spectral band 
!     sumsol5           scalar used to accumulate sum of toa solar flux
!                       over a spectral interval common to a fu and liou
!                       (1993) ice crystal band and a solar spectral 
!                       band 
!     sumplanck         sum of the planck function over a spectral
!                       interval common to a cloud emissivity band and
!                       a water substance band
!     xtemv             temperature at which planck function is 
!                       evaluated [  deg k ]
!     unit              io unit for reading nml file and writing logfile
!     io                error status returned from io operation  
!     ierr              error code
!     nivl              fu band index for lw case
!     nband             cloud band index for lw portion of routine,
!                       sw parameterization band index when doing sw
!     nivl1             cloud droplet band index, sw case
!     nivl2             fu (1996) ice band index, sw case
!     nivl3             rain band index, sw case
!     nivl4             snow band index, sw case
!     nivl5             fu and liou (1993) ice crystal band index, sw
!                       case
!     i,j,k,n,ib,nw,ni  do-loop indices
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
        
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call constants_init
      call rad_utilities_init
      call longwave_params_init
      call diag_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=microphys_rad_nml, iostat=io)
      ierr = check_nml_error(io,'microphys_rad_nml')
#else
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=microphys_rad_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'microphys_rad_nml')
        enddo
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write namelist and version number to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                      write (logunit, nml=microphys_rad_nml)

      min_cld_drop_rad  = min_cld_drop_rad_in
      max_cld_drop_rad  = max_cld_drop_rad_in
      min_cld_ice_size = min_cld_ice_size_in
      max_cld_ice_size = max_cld_ice_size_in

!--------------------------------------------------------------------
!    verify that Lw_control%do_lwcldemiss has been initialized.
!--------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss_iz) then
      else
        call error_mesg ('microphys_rad_mod', &
          ' Lw_control%do_lwcldemiss has not been initialized', FATAL)
      endif

!--------------------------------------------------------------------
!    perform consistency checks between lwem_form and desired lwcld 
!    emissivity. do_lwcldemiss being true implies a multi-band emiss-
!    ivity parameterization; do_lw_micro implies a microphysically-
!    based parameterization. ebert-curry is a single band, microphysic-
!    ally based scheme, fuliou is a multi-band microphysically based
!    scheme.
!--------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss)    then
        if (trim(lwem_form) == 'fuliou') then
        else if (trim(lwem_form) == 'ebertcurry') then
          call error_mesg('microphys_rad_mod',  &
              'ebert-curry not implemented for multi-bands', FATAL)
        else
          call error_mesg('microphys_rad_mod',  &
           'incorrect specification of lwem_form for multi-band', FATAL)
        endif
      else
        if (trim(lwem_form) == 'fuliou') then
          call error_mesg('microphys_rad_mod',  &
          'fu parameterization implemented only for multi-band', FATAL)
        else if (trim(lwem_form) == 'ebertcurry') then
        else
          call error_mesg('microphys_rad_mod',  &
         'incorrect specification of lwem_form for single band', FATAL)
        endif
      endif
      if (Cldrad_control%do_strat_clouds_iz) then     
        if (Cldrad_control%do_strat_clouds) then     
          if (Cldrad_control%do_stochastic_clouds_iz) then
            if (Cldrad_control%do_stochastic_clouds) then
              if (trim(lwem_form)  == 'ebertcurry') then
                call error_mesg ('microphys_rad_mod',  &
              'ebert-curry not allowed with stochastic clouds', FATAL)
              endif
            endif
          else
            call error_mesg ('microphys_rad_mod', &
        'Cldrad_control%do_stochastic_clouds has not been initialized',&
                                                                 FATAL)
          endif
        endif
      else
        call error_mesg ('microphys_rad_mod', &
          ' Cldrad_control%do_strat_clouds has not been initialized',  &
                                                                FATAL)
      endif

!--------------------------------------------------------------------
!    compute band-averaged coefficients for h2o forms in infrared 
!    frequency ranges for the fuliou multi-band parameterization.
!    the actual extinction coefficients (to become emissivities)
!    (and other coefficients) will be calculated as time-dependent
!    quantities.
!    at present, the species included are:
!    1)   cloud drops
!    2)   ice crystals
!    3)   cloud rain
!    4)   cloud snow
!--------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss) then

!--------------------------------------------------------------------
!    verify that  Lw_parameters%lw_band_resolution has been initialized.
!--------------------------------------------------------------------
      if (Lw_parameters%lw_band_resolution_iz) then
      else
        call error_mesg ('microphys_rad_mod', &
         'Lw_parameters%lw_band_resolution has not been initialized',&
                                                                FATAL)
      endif

!--------------------------------------------------------------------
!    determine the band indices in the full lw spectrum which corres-
!    pond to the high wavenumber boundary of the bands used in the 
!    Fu-Liou parameterization. the resolution of the full lw spectrum
!    is given by Lw_parameters%lw_band_resolution.
!--------------------------------------------------------------------
        iendfubands(:) = INT((endfubands(:) + 0.01)/  &
                                     Lw_parameters%lw_band_resolution)

!--------------------------------------------------------------------
!    compute weighting function for each wavenumber band. according to 
!    Fu and Liou, this should be the Planck function evaluated at -40C.
!    (src1nb)
!--------------------------------------------------------------------
        do n=1,NBLW 
          centnb    = 5.0 + (n - 1)*Lw_parameters%lw_band_resolution
          c1        = (3.7412E-05)*centnb   **3
          x         = 1.4387E+00*centnb   /xtemv
          x1        = EXP(x   )
          sc        = c1   /(x1    - 1.0E+00)
          src1nb(n) = Lw_parameters%lw_band_resolution*sc
        end do
 
!--------------------------------------------------------------------
!    add the planck function from each full-spectrum band to the
!    proper cloud band sum. 
!--------------------------------------------------------------------
        planckcldband(:) = 0.0E+00
        do n=1,N_EMISS_BDS
          do ib=istartcldband(n),iendcldband(n)
            planckcldband(n) = planckcldband(n) + src1nb(ib)
          end do
        end do

!--------------------------------------------------------------------
!    contributions from the last cloud band (1400-2200 cm-1) are added
!    to the first cloud band.
!--------------------------------------------------------------------
        do ib=istartcldband(N_EMISS_BDS+1),iendcldband(N_EMISS_BDS+1)
          planckcldband(1) = planckcldband(1) + src1nb(ib)
        end do
 
!---------------------------------------------------------------------
!    compute the sum of the planck function over each spectral segment
!    created when the fu bands and the cloud bands are overlapped. nivl
!    is the fu band index while nband is the cloud band index.
!    planckivlicecld(nband, nivl) is the sum of the planck function in 
!    the portion of the wave number spectrum common to cloud band nband 
!    and fu band nivl. nivl1icecld(nband) is the fu band index 
!    corresponding to the lowest wave number of the nbandth cloud band,
!    and nivl2icecld(nband) is the fu band index corresponding to the
!    highest wave number of the nbandth cloud band.
!---------------------------------------------------------------------
        nivl = 1
        sumplanck = 0.0
        nband = 1
        planckivlicecld(:,:) = 0.0
        nivl1lwicecld(1) = 1
 
        do nw = 1,NBLW
          sumplanck = sumplanck + src1nb(nw)
          if (nw == iendfubands(nivl)) then
            planckivlicecld(nband,nivl) = sumplanck
            sumplanck = 0.0
          end if
          if (nw == iendcldband(nband)) then
            if (nw /= iendfubands(nivl)) then
              planckivlicecld(nband,nivl) = sumplanck 
              sumplanck = 0.0
            end if
            nivl2lwicecld(nband) = nivl
            nband = nband + 1
            if (nband <= N_EMISS_BDS+1) then
              if (nw == iendfubands(nivl)) then
                nivl1lwicecld(nband) = nivl + 1
              else
                nivl1lwicecld(nband) = nivl
              end if
            end if
          end if
          if (nw == iendfubands(nivl)) nivl = nivl + 1
          if (nw >= iendcldband(N_EMISS_BDS+1)) exit
        end do

!--------------------------------------------------------------------
!    compute the fraction of the total planck function in cloud band
!    n that is present in fu band ni. the sum of fulwwts over ni should
!    be 1.00.
!--------------------------------------------------------------------
        fulwwts(:,:) = 0.0E+00
        do n=1,N_EMISS_BDS
          do ni=nivl1lwicecld(n),nivl2lwicecld(n)
            fulwwts(n,ni) = planckivlicecld(n,ni)/planckcldband(n)
          end do
        end do

!--------------------------------------------------------------------
!    the contributions from cloud band (N_EMISS_BDS+1) are included
!    in band 1.
!--------------------------------------------------------------------
        do ni=nivl1lwicecld(N_EMISS_BDS+1),nivl2lwicecld(N_EMISS_BDS+1)
          fulwwts(1,ni) = planckivlicecld(N_EMISS_BDS+1,ni)/   &
                          planckcldband(1)
        end do
      endif  ! (do_lwcldemiss)

!--------------------------------------------------------------------
!    verify that Cldrad_control%do_sw_micro has been initialized.
!--------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro_iz) then
      else
        call error_mesg ('microphys_rad_mod', &
         'Cldrad_control%do_sw_micro has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    the following section is executed when the shortwave parameter-
!    ization is based on microphysical information.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then

!--------------------------------------------------------------------
!    make certain esfsw_parameters_mod has been initialized.
!--------------------------------------------------------------------
        call esfsw_parameters_init

!--------------------------------------------------------------------
!    verify consistency between the highest wavenumber in the solar 
!    spectrum and the highest wavenumber in the various particle
!    spectral intervals to assure that all solar spectral bands are 
!    assigned to a particle band.
!--------------------------------------------------------------------
        if (Solar_spect%tot_wvnums /= endliqcldwvn(NLIQCLDIVLS) .or.  &
            Solar_spect%tot_wvnums /= endicecldwvn(NICECLDIVLS) .or.  &
            Solar_spect%tot_wvnums /=      &
                               endicesolcldwvn(NICESOLARCLDIVLS) .or.  &
            Solar_spect%tot_wvnums /= endraincldwvn(NRAINCLDIVLS) .or.  &
            Solar_spect%tot_wvnums /= endsnowcldwvn(NSNOWCLDIVLS) ) then
          call error_mesg ( 'microphys_rad_mod',  &
              ' highest wavenumber in particle spectrum differs '//&
                'from highest wavenumber in solar spectrum ', FATAL)
        endif

!--------------------------------------------------------------------
!    allocate the module variables that are needed for the shortwave
!    parameterization.
!--------------------------------------------------------------------
        allocate ( nivl1liqcld    (Solar_spect%nbands),    &
                   nivl1icecld    (Solar_spect%nbands),   &
                   nivl1icesolcld (Solar_spect%nbands),   &
                   nivl1raincld   (Solar_spect%nbands),  &
                   nivl1snowcld   (Solar_spect%nbands),  &
                   nivl2liqcld    (Solar_spect%nbands),   &
                   nivl2icecld    (Solar_spect%nbands),   &
                   nivl2icesolcld (Solar_spect%nbands),   &
                   nivl2raincld   (Solar_spect%nbands),   & 
                   nivl2snowcld   (Solar_spect%nbands) ) 

        allocate ( solivlicecld   (Solar_spect%nbands, nicecldivls),  &
                   solivlicesolcld                                  &
                               (Solar_spect%nbands, nicesolarcldivls), &
                   solivlliqcld   (Solar_spect%nbands, nliqcldivls), &
                   solivlraincld  (Solar_spect%nbands, nraincldivls), &
                   solivlsnowcld  (Solar_spect%nbands, nsnowcldivls) )

!---------------------------------------------------------------------
!    compute the sum of the toa solar flux over each spectral segment
!    created when the individual particle bands (fu ice, liquid, rain, 
!    snow, icesolar ice) and the parameterization bands are overlapped. 
!    nivlx is the particle band index while nband is the parameteriz-
!    ation band index. thus solivlxxxcld(nband, nivlx) is the sum of the
!    toa solar flux in the portion of the wave number spectrum common 
!    to parameterization band nband and particle band nivl. 
!    nivl1xxxcld(nband) is the particle band index corresponding to the
!    lowest wave number of the nbandth parameterization band, and 
!    nivl2xxxcld(nband) is the particle band index corresponding to the
!    highest wave number of the nbandth parameterization band.
!    the naming convention is as follows:
!        xxx1 refers to liquid cloud particles
!        xxx2 refers to fu ice cloud particles
!        xxx3 refers to rain particles
!        xxx4 refers to snow particles      
!        xxx5 refers to icesolar ice particles
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    initialize the indices and sums.
!--------------------------------------------------------------------
        nivl1 = 1
        nivl2 = 1
        nivl3 = 1
        nivl4 = 1
        nivl5 = 1
        sumsol1 = 0.0
        sumsol2 = 0.0
        sumsol3 = 0.0
        sumsol4 = 0.0
        sumsol5 = 0.0
        nband = 1
        solivlliqcld(:,:) = 0.0
        solivlicecld(:,:) = 0.0
        solivlicesolcld(:,:) = 0.0
        solivlraincld(:,:) = 0.0
        solivlsnowcld(:,:) = 0.0

!---------------------------------------------------------------------
!    the single scattering properties for wavenumbers < the lower limit
!    to which the various parameterizations apply are assigned the
!    values in the lowest interval of the parameterization; thus all
!    solar spectral parameterization bands are assigned to a 
!    water-species parameterization band.
!---------------------------------------------------------------------
        nivl1liqcld(1) = 1
        nivl1icecld(1) = 1
        nivl1icesolcld(1) = 1
        nivl1raincld(1) = 1
        nivl1snowcld(1) = 1
 
!---------------------------------------------------------------------
!    integrate over wavenumber, summing the solar spectral bands con-
!    sistent with the parameterization band structure and the particle 
!    band structure. when the loop is ended, the solar flux in each
!    spectral region generated by overlapping the parameterization
!    spectrum and a particular particle spectrum will be resident in
!    the solivlxxxcld arrays.
!---------------------------------------------------------------------
        do nw = 1,Solar_spect%tot_wvnums           
          sumsol1 = sumsol1 + Solar_spect%solarfluxtoa(nw) 
          sumsol2 = sumsol2 + Solar_spect%solarfluxtoa(nw) 
          sumsol3 = sumsol3 + Solar_spect%solarfluxtoa(nw) 
          sumsol4 = sumsol4 + Solar_spect%solarfluxtoa(nw) 
          sumsol5 = sumsol5 + Solar_spect%solarfluxtoa(nw) 
          if ( nw == endliqcldwvn(nivl1) ) then
            solivlliqcld(nband,nivl1) = sumsol1
            sumsol1 = 0.0
          end if
          if ( nw == endicecldwvn(nivl2) ) then
            solivlicecld(nband,nivl2) = sumsol2
            sumsol2 = 0.0
          end if
          if ( nw == endraincldwvn(nivl3) ) then
            solivlraincld(nband,nivl3) = sumsol3
            sumsol3 = 0.0
          end if
          if ( nw == endsnowcldwvn(nivl4) ) then
            solivlsnowcld(nband,nivl4) = sumsol4
            sumsol4 = 0.0
          end if
          if ( nw == endicesolcldwvn(nivl5) ) then
            solivlicesolcld(nband,nivl5) = sumsol5
            sumsol5 = 0.0
          end if
          if ( nw == Solar_spect%endwvnbands(nband) ) then
            if ( nw /= endliqcldwvn(nivl1) ) then
              solivlliqcld(nband,nivl1) = sumsol1 
              sumsol1 = 0.0
            end if
            if ( nw /= endicecldwvn(nivl2) ) then
              solivlicecld(nband,nivl2) = sumsol2 
              sumsol2 = 0.0
            end if
            if ( nw /= endraincldwvn(nivl3) ) then
              solivlraincld(nband,nivl3) = sumsol3 
              sumsol3 = 0.0
            end if
            if ( nw /= endsnowcldwvn(nivl4) ) then
              solivlsnowcld(nband,nivl4) = sumsol4 
              sumsol4 = 0.0
            end if
            if ( nw /= endicesolcldwvn(nivl5) ) then
              solivlicesolcld(nband,nivl5) = sumsol5 
              sumsol5 = 0.0
            end if
            nivl2liqcld(nband) = nivl1
            nivl2icecld(nband) = nivl2
            nivl2raincld(nband) = nivl3
            nivl2snowcld(nband) = nivl4
            nivl2icesolcld(nband) = nivl5
 
            nband = nband + 1
 
            if ( nband <= Solar_spect%nbands ) then
 
              if ( nw == endliqcldwvn(nivl1) ) then
                nivl1liqcld(nband) = nivl1 + 1
              else
                nivl1liqcld(nband) = nivl1
              end if
              if ( nw == endicecldwvn(nivl2) ) then
                nivl1icecld(nband) = nivl2 + 1
              else
                nivl1icecld(nband) = nivl2
              end if
              if ( nw == endraincldwvn(nivl3) ) then
                nivl1raincld(nband) = nivl3 + 1
              else
                nivl1raincld(nband) = nivl3
              end if
              if ( nw == endsnowcldwvn(nivl4) ) then
                nivl1snowcld(nband) = nivl4 + 1
              else
                nivl1snowcld(nband) = nivl4
              end if
              if ( nw == endicesolcldwvn(nivl5) ) then
                nivl1icesolcld(nband) = nivl5 + 1
              else
                nivl1icesolcld(nband) = nivl5
              end if
            end if
          end if
          if ( nw == endliqcldwvn(nivl1) ) nivl1 = nivl1 + 1
          if ( nw == endicecldwvn(nivl2) ) nivl2 = nivl2 + 1
          if ( nw == endraincldwvn(nivl3) ) nivl3 = nivl3 + 1
          if ( nw == endsnowcldwvn(nivl4) ) nivl4 = nivl4 + 1
          if ( nw == endicesolcldwvn(nivl5) ) nivl5 = nivl5 + 1
        end do
      endif   ! (do_sw_micro)

!-----------------------------------------------------------------
!    mark module as initialized.
!----------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------


end subroutine microphys_rad_init


!###################################################################
! <SUBROUTINE NAME="microphys_sw_driver">
!  <OVERVIEW>
!   Subroutine to deploy micro physics radiation calculation,
!   particularly for cloud parameterizations in the shortwave
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes cloud micro physics parameters and
!   calculate broad band cloud radiation parameters. For example
!   the input parameters are cloud droplet concentration, size,
!   and composition; the output parameters are cloud extinction
!   coefficient, scattering coefficient, and assymmetry parameters.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_sw_driver (is, ie, js, je, Cloud_microphysics,  &
!                                Cloud_rad_props, donner_flag )
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of x dimension in current physics window
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending indice of x dimension in current physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of y dimension in current physics window
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending indice of y dimension in current physics window
!  </IN>
!  <IN NAME="Cloud_microphysics" TYPE="microphysics_type">
!   derived type variable containing cloud 
!                          microphysical specification information 
!  </IN>
!  <INOUT NAME="Cloud_rad_props" TYPE="microrad_properties_type">
!   derived type variable containing the micro-
!                          physically-based sw cloud radiative proper-
!                          ties [ microrad_properties_type ]
!                              the components defined in this routine:
!                            %cldext   parameterization band values of 
!                                      the cloud extinction coefficient 
!                                      [ km**(-1) ]   
!                            %cldsct   parameterization band values of 
!                                      the cloud scattering coefficient 
!                                      [ km**(-1) ] 
!                            %cldasymm parameterization band values of 
!                                      the asymmetry factor 
!                                      [ dimensionless ]
!  </INOUT>
!  <IN NAME="donner_flag" TYPE="logical">
!   OPTIONAL: logical flag which if present indicates
!                           that clouds from donner_deep_mod are being
!                           processed, and that an ice parameterization
!                           associated with that scheme (which differs
!                           from that used by strat_cloud_mod) is to
!                           be used.
!  </IN>
! </SUBROUTINE>
subroutine microphys_sw_driver (is, ie, js, je, Cloud_microphysics,  &
                                Cloud_rad_props, Micro_rad_props, &
                                donner_flag )

!---------------------------------------------------------------------
!    microphys_sw_driver obtains microphysically-based cloud shortwave
!    radiative properties for the cloud field described by Cloud_micro-
!    physics and returnms them in Cloud_rad_props.
!---------------------------------------------------------------------

integer,                        intent(in)      :: is, ie, js, je
type(microphysics_type),        intent(in)      :: Cloud_microphysics
type(microrad_properties_type), intent(inout), optional   :: Micro_rad_props
type(cldrad_properties_type), intent(inout), optional   :: Cloud_rad_props
logical,                        intent(in),                         &
                                       optional :: donner_flag 

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je         starting/ending subdomain i,j indices of data
!                          in the physics_window being integrated
!      Cloud_microphysics  derived type variable containing cloud 
!                          microphysical specification information 
!                          [ microphysics_type ]
!
!   intent(inout) variable:
!
!      Cloud_rad_props     derived type variable containing the micro-
!                          physically-based sw cloud radiative proper-
!                          ties [ microrad_properties_type ]
!                              the components defined in this routine:
!                            %cldext   parameterization band values of 
!                                      the cloud extinction coefficient 
!                                      [ km**(-1) ]   
!                            %cldsct   parameterization band values of 
!                                      the cloud scattering coefficient 
!                                      [ km**(-1) ] 
!                            %cldasymm parameterization band values of 
!                                      the asymmetry factor 
!                                      [ dimensionless ]
!
!    intent(in), optional variable:
!
!      donner_flag          logical flag which if present indicates
!                           that clouds from donner_deep_mod are being
!                           processed, and that an ice parameterization
!                           associated with that scheme (which differs
!                           from that used by strat_cloud_mod) is to
!                           be used. 
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:                                                  

      real, dimension   &
                   (size(Cloud_microphysics%size_drop,1), &
                    size(Cloud_microphysics%size_drop,2), &
                    size(Cloud_microphysics%size_drop,3), &
                                      Solar_spect%nbands) :: &
                       size_drop, size_ice, conc_drop, conc_ice
      logical, dimension &
                   (size(Cloud_microphysics%size_drop,1), &
                    size(Cloud_microphysics%size_drop,2), &
                    size(Cloud_microphysics%size_drop,3), &
                                      Solar_spect%nbands) :: &
                                                        dge_column

      logical       :: do_dge_sw, isccp_call
      integer       :: nbmax
      integer       :: n
      integer       :: nnn, nbprofiles, nonly

!----------------------------------------------------------------------
!  local variables:                                                  
!
!      do_dge_sw     logical flag; when .true., indicates that 
!                    donner_deep_mod clouds are being processed
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!    define variable indicating whether doner_deep_mods clouds or 
!    large-scale clouds are currently being processed.
!---------------------------------------------------------------------
      if (present (donner_flag )) then
        if (donner_flag) then
          do_dge_sw = .true.
        else
          do_dge_sw = .false.
        endif
      else
        do_dge_sw = .false.
     endif
     if (Cldrad_control%using_fu2007) then
       do_dge_sw = .true.
     endif

     dge_column = do_dge_sw
     isccp_call = .false.

!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic sw
!    clouds has been activated, define the microphysical inputs
!    for each sw parameterization band and the number of such fields.
!---------------------------------------------------------------------
      if (.not. present (donner_flag)) then
        if (Cldrad_control%do_stochastic_clouds) then
          size_drop = Cloud_microphysics%sw_stoch_size_drop
          size_ice  = Cloud_microphysics%sw_stoch_size_ice
          conc_drop = Cloud_microphysics%sw_stoch_conc_drop
          conc_ice  = Cloud_microphysics%sw_stoch_conc_ice
          if (Cldrad_control%do_ica_calcs) then
            nbprofiles = Solar_spect%nbands
            nbmax = 1
          else
            nbprofiles = 1
            nbmax = Solar_spect%nbands
          endif

!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic sw
!    clouds are not desired, define the microphysical inputs
!    as obtained previously from the full predicted cloud field.
!---------------------------------------------------------------------
        else
          size_drop(:,:,:,1) = Cloud_microphysics%size_drop
          size_ice(:,:,:,1)  = Cloud_microphysics%size_ice
          conc_drop(:,:,:,1) = Cloud_microphysics%conc_drop
          conc_ice(:,:,:,1)  = Cloud_microphysics%conc_ice
          nbmax = 1
        endif

!---------------------------------------------------------------------
!    call cloudpar to define microphysically-based sw cloud properties.
!---------------------------------------------------------------------
      if (present(Cloud_rad_props)) then
      do nnn=1,nbprofiles ! loop over profiles
      nonly = 0
      call cloudpar                                                 &
                     (nonly, nbmax, nnn, size_drop, size_ice,   &
                      Cloud_microphysics%size_rain,              &
                      conc_drop, conc_ice, &
                      Cloud_microphysics%conc_rain, &
                      Cloud_microphysics%conc_snow, do_dge_sw,   &
                      dge_column, isccp_call, &
                      Cloud_rad_props%cldext(:,:,:,:,nnn),   &
                      Cloud_rad_props%cldsct(:,:,:,:,nnn), &
                      Cloud_rad_props%cldasymm(:,:,:,:,nnn))
       end do
      else
      nnn = 1
      nonly = 0
      call cloudpar                                                 &
                     (nonly, nbmax, nnn, size_drop, size_ice,   &
                      Cloud_microphysics%size_rain,              &
                      conc_drop, conc_ice, &
                      Cloud_microphysics%conc_rain, &
                      Cloud_microphysics%conc_snow, do_dge_sw,   &
                      dge_column, isccp_call, &
                      Micro_rad_props%cldext, Micro_rad_props%cldsct, &
                      Micro_rad_props%cldasymm)
       endif

!-------------------------------------------------------------------
!    if donner cloud fields are being processed, there is currently
!    no stochastic component. use the same properties for each sw
!    parameterization band.
!-------------------------------------------------------------------
      else  ! (donner_flag)
        nnn = 1
        nbprofiles = 1
        nbmax = Solar_spect%nbands
        nonly = 0
        do n=1, Solar_spect%nbands
        size_drop(:,:,:,n) = Cloud_microphysics%size_drop
        size_ice(:,:,:,n)  = Cloud_microphysics%size_ice
        conc_drop(:,:,:,n) = Cloud_microphysics%conc_drop
        conc_ice(:,:,:,n)  = Cloud_microphysics%conc_ice
        end do

        call cloudpar                                 &
                     (nonly, nbmax, nnn, size_drop, size_ice,   &
                      Cloud_microphysics%size_rain,              &
                      conc_drop, conc_ice, &
                      Cloud_microphysics%conc_rain, &
                      Cloud_microphysics%conc_snow, do_dge_sw,   &
                      dge_column, isccp_call, &
                      Micro_rad_props%cldext, Micro_rad_props%cldsct, &
                      Micro_rad_props%cldasymm)
      endif ! (present(donner_flag))
 
!--------------------------------------------------------------------



end subroutine microphys_sw_driver


!#####################################################################
! <SUBROUTINE NAME="microphys_lw_driver">
!  <OVERVIEW>
!   Subroutine to deploy micro physics radiation calculation,
!   particularly for cloud parameterizations in the longwave
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes cloud micro physics parameters and
!   calculate broad band cloud radiation parameters. For example
!   the input parameters are cloud droplet concentration, size,
!   and composition; the output parameters are cloud extinction
!   coefficient, scattering coefficient, and assymmetry parameters.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_lw_driver (is, ie, js, je, Cloud_microphysics,  &
!                                Cloud_rad_props, donner_flag )
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of x dimension in current physics window
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending indice of x dimension in current physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of y dimension in current physics window
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending indice of y dimension in current physics window
!  </IN>
!  <IN NAME="Cloud_microphysics" TYPE="microphysics_type">
!   derived type variable containing cloud 
!                          microphysical specification information 
!  </IN>
!  <INOUT NAME="Cloud_rad_props" TYPE="microrad_properties_type">
!   derived type variable containing the micro-
!                          physically-based sw cloud radiative proper-
!                          ties [ microrad_properties_type ]
!                              the components defined in this routine:
!                            %cldext   parameterization band values of 
!                                      the cloud extinction coefficient 
!                                      [ km**(-1) ]   
!                            %cldsct   parameterization band values of 
!                                      the cloud scattering coefficient 
!                                      [ km**(-1) ] 
!                            %cldasymm parameterization band values of 
!                                      the asymmetry factor 
!                                      [ dimensionless ]
!  </INOUT>
!  <IN NAME="donner_flag" TYPE="logical">
!   OPTIONAL: logical flag which if present indicates
!                           that clouds from donner_deep_mod are being
!                           processed, and that an ice parameterization
!                           associated with that scheme (which differs
!                           from that used by strat_cloud_mod) is to
!                           be used.
!  </IN>
! </SUBROUTINE>
subroutine microphys_lw_driver (is, ie, js, je, Cloud_microphysics,  &
                                Cloud_rad_props, Micro_rad_props, &
                                donner_flag )

!---------------------------------------------------------------------
!    microphys_lw_driver obtains microphysically-based cloud longwave
!    radiative properties for the cloud field described by Cloud_micro-
!    physics and returnms them in Cloud_rad_props.
!---------------------------------------------------------------------

integer,                        intent(in)      :: is, ie, js, je
type(microphysics_type),        intent(in)      :: Cloud_microphysics
type(cldrad_properties_type), intent(inout),optional   :: Cloud_rad_props
type(microrad_properties_type), intent(inout),optional   :: Micro_rad_props
logical,                        intent(in),                         &
                                       optional :: donner_flag 

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je         starting/ending subdomain i,j indices of data
!                          in the physics_window being integrated
!      Cloud_microphysics  derived type variable containing cloud 
!                          microphysical specification information 
!                          [ microphysics_type ]
!
!   intent(inout) variable:
!
!      Cloud_rad_props     derived type variable containing the micro-
!                          physically-based sw cloud radiative proper-
!                          ties [ microrad_properties_type ]
!                              the component defined in this routine:
!                             %abscoeff absorption coefficient for 
!                                       clouds in each of the longwave 
!                                       frequency bands [ km**(-1) ]
!
!    intent(in), optional variable:
!
!      donner_flag          logical flag which if present indicates
!                           that clouds from donner_deep_mod are being
!                           processed, and that an ice parameterization
!                           associated with that scheme (which differs
!                           from that used by strat_cloud_mod) is to
!                           be used. 
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:                                                  

      real, dimension   &
                  (size(Cloud_microphysics%size_drop,1), &
                   size(Cloud_microphysics%size_drop,2), &
                   size(Cloud_microphysics%size_drop,3), &
                   Cldrad_control%nlwcldb) :: &
                        size_drop, size_ice, conc_drop, conc_ice
      logical       :: do_dge_lw, isccp_call
      integer :: nbmax, nbprofiles, nnn, nonly
      logical, dimension   &
                  (size(Cloud_microphysics%size_drop,1), &
                   size(Cloud_microphysics%size_drop,2), &
                   size(Cloud_microphysics%size_drop,3), &
                   Cldrad_control%nlwcldb) ::   dge_column

!----------------------------------------------------------------------
!  local variables:                                                  
!
!      do_dge_lw     logical flag; when .true., indicates that 
!                    donner_deep_mod clouds are being processed
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!    define variable indicating whether donner_deep_mods clouds or 
!    large-scale clouds are currently being processed.
!---------------------------------------------------------------------
      if (present (donner_flag )) then
        if (donner_flag) then
          do_dge_lw = .true.
        else
          do_dge_lw = .false.
        endif
      else
        do_dge_lw = .false.
      endif

      if (Cldrad_control%using_fu2007) then
        do_dge_lw = .true.
      endif

       dge_column = do_dge_lw
       isccp_call = .false.

!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic lw
!    clouds has been activated, define the microphysical inputs
!    for each lw parameterization band and the number of such fields.
!---------------------------------------------------------------------
      if ( .not. present (donner_flag )) then
        if (Cldrad_control%do_stochastic_clouds) then
          size_drop = Cloud_microphysics%lw_stoch_size_drop
          size_ice  = Cloud_microphysics%lw_stoch_size_ice
          conc_drop = Cloud_microphysics%lw_stoch_conc_drop
          conc_ice  = Cloud_microphysics%lw_stoch_conc_ice
          if (Cldrad_control%do_ica_calcs) then
            nbprofiles = Cldrad_control%nlwcldb
            nbmax = 1
          else
            nbprofiles = 1
            nbmax = Cldrad_control%nlwcldb
          endif


!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic lw
!    clouds are not desired, define the microphysical inputs
!    as obtained previously from the full predicted cloud field.
!---------------------------------------------------------------------
        else
          size_drop(:,:,:,1) = Cloud_microphysics%size_drop
          size_ice(:,:,:,1)  = Cloud_microphysics%size_ice
          conc_drop(:,:,:,1) = Cloud_microphysics%conc_drop
          conc_ice(:,:,:,1)  = Cloud_microphysics%conc_ice
          nbmax = 1
        endif

!-------------------------------------------------------------------
!    if the fuliou lw emissivity was selected, call cloud_lwpar to 
!    compute multi-band emissivities based on fu parameterizations for
!    cloud drop, cloud ice, snow and rain. 
!---------------------------------------------------------------------
      if (trim(lwem_form) == 'fuliou') then
      if (present(Cloud_rad_props)) then
      do nnn=1,nbprofiles ! loop over profiles
        nonly = 0
        call cloud_lwpar (nonly, nbmax, nnn, size_drop, size_ice,   &
                          Cloud_microphysics%size_rain,              &
                          conc_drop, conc_ice, &
                          Cloud_microphysics%conc_rain, &
                          Cloud_microphysics%conc_snow, do_dge_lw,   &
                          dge_column, isccp_call, &
                          Cloud_rad_props%abscoeff(:,:,:,:,nnn))
     end do
    else  ! ((present(Cloud_rad_props))
        nnn = 1
        nonly = 0
        call cloud_lwpar (nonly, nbmax, nnn, size_drop, size_ice,   &
                          Cloud_microphysics%size_rain,              &
                          conc_drop, conc_ice, &
                          Cloud_microphysics%conc_rain, &
                          Cloud_microphysics%conc_snow, do_dge_lw,   &
                          dge_column, isccp_call, &
                          Micro_rad_props%abscoeff)
     endif   ! ((present(Cloud_rad_props))

!---------------------------------------------------------------------
!    if the ebert-curry emissivity was selected, call cloud_lwem_oneband
!    to compute a single value for the lw emissivity (including effects
!    of drops and ice) based on the ebert and curry parameterization.
!---------------------------------------------------------------------
      else if (trim(lwem_form) == 'ebertcurry') then
       if (present(Cloud_rad_props)) then
        call cloud_lwem_oneband (Cloud_microphysics%conc_drop,   &
                                 Cloud_microphysics%conc_ice,    &
                                 Cloud_microphysics%size_drop,    &
                                 Cloud_microphysics%size_ice,      &
                                 Cloud_rad_props%abscoeff(:,:,:,1,1))
       else
        call cloud_lwem_oneband (Cloud_microphysics%conc_drop,   &
                                 Cloud_microphysics%conc_ice,    &
                                 Cloud_microphysics%size_drop,    &
                                 Cloud_microphysics%size_ice,      &
                                 Micro_rad_props%abscoeff(:,:,:,1))
       endif
      endif

!-------------------------------------------------------------------
!    if donner cloud fields are being processed, there is currently
!    no stochastic component. use the same properties for each lw
!    parameterization band.
!-------------------------------------------------------------------
    else ! (donner_flag)
      nbmax = 1
      size_drop(:,:,:,1) = Cloud_microphysics%size_drop
      size_ice(:,:,:,1)  = Cloud_microphysics%size_ice
      conc_drop(:,:,:,1) = Cloud_microphysics%conc_drop
      conc_ice(:,:,:,1)  = Cloud_microphysics%conc_ice
 
!---------------------------------------------------------------------
!    if the fuliou lw emissivity was selected, call cloud_lwpar to
!    compute multi-band emissivities based on fu parameterizations for
!    cloud drop, cloud ice, snow and rain.
!---------------------------------------------------------------------
     if (trim(lwem_form) == 'fuliou') then
       nnn = 1
       nonly = 0
       call cloud_lwpar (nonly, nbmax, nnn,  size_drop, size_ice,   &
                         Cloud_microphysics%size_rain,              &
                         conc_drop, conc_ice, &
                         Cloud_microphysics%conc_rain, &
                         Cloud_microphysics%conc_snow, do_dge_lw,   &
                          dge_column, isccp_call, &
                         Micro_rad_props%abscoeff)

!---------------------------------------------------------------------
!    if the ebert-curry emissivity was selected, call cloud_lwem_oneband
!    to compute a single value for the lw emissivity (including effects
!    of drops and ice) based on the ebert and curry parameterization.
!---------------------------------------------------------------------
     else if (trim(lwem_form) == 'ebertcurry') then
       call cloud_lwem_oneband (Cloud_microphysics%conc_drop,   &
                                Cloud_microphysics%conc_ice,    &
                                Cloud_microphysics%size_drop,    &
                                Cloud_microphysics%size_ice,      &
                                Cloud_rad_props%abscoeff(:,:,:,1,1))
     endif
   endif

!--------------------------------------------------------------------



end subroutine microphys_lw_driver  




!###################################################################
! <SUBROUTINE NAME="isccp_microphys_sw_driver">
!  <OVERVIEW>
!   Subroutine to deploy micro physics radiation calculation,
!   particularly for cloud parameterizations in the shortwave
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes cloud micro physics parameters and
!   calculate broad band cloud radiation parameters. For example
!   the input parameters are cloud droplet concentration, size,
!   and composition; the output parameters are cloud extinction
!   coefficient, scattering coefficient, and assymmetry parameters.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_microphys_sw_driver (is, js, iswband,
!                                   Cloud_microphysics,cldext )
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of x dimension in current physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of y dimension in current physics window
!  </IN>
!  <IN NAME="iswband" TYPE="integer">
!   swband whose extinction we desire
!  </IN>
!  <IN NAME="Cloud_microphysics" TYPE="microphysics_type">
!   derived type variable containing cloud 
!                          microphysical specification information 
!  </IN>
!  <OUT NAME="cldext " TYPE="real">
!   derived type variable containing the micro-
!                          physically-based sw cloud radiative proper-
!                          ties [ microrad_properties_type ]
!                              the components defined in this routine:
!                            %cldext   parameterization band values of 
!                                      the cloud extinction coefficient 
!                                      [ km**(-1) ]   
!  </OUT>
! </SUBROUTINE>
subroutine isccp_microphys_sw_driver (is, js, iswband, &
                                Cloud_microphysics, cldext )

!---------------------------------------------------------------------
!    isccp_microphys_sw_driver obtains microphysically-based cloud shortwave
!    radiative properties for the cloud field described by Cloud_micro-
!    physics and returnms them in Cloud_rad_props.
!---------------------------------------------------------------------

integer,                        intent(in)      :: is, js
integer,                        intent(in)      :: iswband
type(microphysics_type),        intent(in)      :: Cloud_microphysics
real, intent(out), dimension(:,:,:,:)           :: cldext

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is, js              starting  subdomain i,j indices of data
!                          in the physics_window being integrated
!      Cloud_microphysics  derived type variable containing cloud 
!                          microphysical specification information 
!                          [ microphysics_type ]
!
!   intent(out) variable:
!
!      cldext   parameterization band values of 
!                                      the cloud extinction coefficient 

!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:                                                  

      real, dimension   &
                   (size(Cloud_microphysics%size_drop,1), &
                    size(Cloud_microphysics%size_drop,2), &
                    size(Cloud_microphysics%size_drop,3), &
                                      Solar_spect%nbands) :: &
                       size_drop, size_ice, conc_drop, conc_ice, &
                       tmpcldext,tmpcldsct,tmpcldasymm
      logical, dimension &
                   (size(Cloud_microphysics%size_drop,1), &
                    size(Cloud_microphysics%size_drop,2), &
                    size(Cloud_microphysics%size_drop,3), &
                                      Solar_spect%nbands) :: dge_column

      
      logical       :: do_dge_sw, isccp_call
      integer       :: nbmax
      integer       :: i,j,k
      integer       :: nnn, nbprofiles, nonly

!----------------------------------------------------------------------
!  local variables:                                                  
!
!      do_dge_sw     logical flag; when .true., indicates that 
!                    donner_deep_mod clouds are being processed
!
!----------------------------------------------------------------------
        
!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!    define variable indicating whether doner_deep_mods clouds or 
!    large-scale clouds are currently being processed.
!   
!    donner clouds will not be accessed here
!---------------------------------------------------------------------
      if (Cldrad_control%using_fu2007) then
        do_dge_sw = .true.
      else
        do_dge_sw = .false.
      endif

      isccp_call = .true.


!---------------------------------------------------------------------
!    loop over profiles
!---------------------------------------------------------------------
     nonly = iswband
     nbmax = Solar_spect%nbands
     nbprofiles = size(Cloud_microphysics%stoch_size_drop,4)
     
     do nnn = 1, nbprofiles

!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic sw
!    clouds has been activated, define the microphysical inputs
!    for each sw parameterization band and the number of such fields.
!---------------------------------------------------------------------
        
      size_drop(:,:,:,nonly) = Cloud_microphysics%stoch_size_drop(:,:,:,nnn)
      size_ice(:,:,:,nonly)  = Cloud_microphysics%stoch_size_ice(:,:,:,nnn)
      conc_drop(:,:,:,nonly) = Cloud_microphysics%stoch_conc_drop(:,:,:,nnn)
      conc_ice(:,:,:,nonly)  = Cloud_microphysics%stoch_conc_ice(:,:,:,nnn)
      if (Cldrad_control%using_fu2007) then
        dge_column = .true.
      else
      do k=1, size(cldext,3)
        do j=1, size(cldext,2)
          do i=1, size(cldext,1)
            if (Cloud_microphysics%stoch_cloud_type(i,j,k,nnn) == 2   &
                        .or. &
                Cloud_microphysics%stoch_cloud_type(i,j,k,nnn) == 3 ) &
                         then 
              dge_column(i,j,k,nonly) = .true.
            else
              dge_column(i,j,k,nonly) = .false.
            endif
          end do
        end do
      end do
      endif
           
!---------------------------------------------------------------------
!    call cloudpar to define microphysically-based sw cloud properties.
!---------------------------------------------------------------------
         call cloudpar                                        &
                     (nonly, nbmax, 1, size_drop, size_ice,   &
                      Cloud_microphysics%size_rain,              &
                      conc_drop, conc_ice, &
                      Cloud_microphysics%conc_rain, &
                      Cloud_microphysics%conc_snow, do_dge_sw,   &
                      dge_column, isccp_call, &
                      tmpcldext, tmpcldsct,tmpcldasymm)
       
!---------------------------------------------------------------------
!   save desired profile
!---------------------------------------------------------------------

         cldext(:,:,:,nnn) = tmpcldext(:,:,:,iswband)

     enddo        !loop over profiles


end subroutine isccp_microphys_sw_driver





!#####################################################################
! <SUBROUTINE NAME="isccp_microphys_lw_driver">
!  <OVERVIEW>
!   Subroutine to deploy micro physics radiation calculation,
!   particularly for cloud parameterizations in the longwave
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes cloud micro physics parameters and
!   calculate broad band cloud radiation parameters. For example
!   the input parameters are cloud droplet concentration, size,
!   and composition; the output parameters are cloud extinction
!   coefficient, scattering coefficient, and assymmetry parameters.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_microphys_lw_driver (is, js, ilwband, &
!                                Cloud_microphysics,abscoeff )
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of x dimension in current physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of y dimension in current physics window
!  </IN>
!  <IN NAME="ilwband" TYPE="integer">
!   lwband whose absorption we desire
!  </IN>
!  <IN NAME="Cloud_microphysics" TYPE="microphysics_type">
!   derived type variable containing cloud 
!                          microphysical specification information 
!  </IN>
!  <OUT NAME="abscoeff" TYPE="real">
!        abscoeff absorption coefficient for 
!        clouds in each of the longwave frequency bands [ km**(-1) ]
!  </OUT>
! </SUBROUTINE>
subroutine isccp_microphys_lw_driver (is, js, ilwband, &
                                Cloud_microphysics, abscoeff )

!---------------------------------------------------------------------
!    isccp_microphys_lw_driver obtains microphysically-based cloud longwave
!    radiative properties for the cloud field described by Cloud_micro-
!    physics and returnms them in Cloud_rad_props.
!---------------------------------------------------------------------

integer,                        intent(in)      :: is, js
integer,                        intent(in)      :: ilwband
type(microphysics_type),        intent(in)      :: Cloud_microphysics
real, intent(out), dimension(:,:,:,:)           :: abscoeff

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je         starting/ending subdomain i,j indices of data
!                          in the physics_window being integrated
!      Cloud_microphysics  derived type variable containing cloud 
!                          microphysical specification information 
!                          [ microphysics_type ]
!
!   intent(out) variable:
!
!      abscoeff absorption coefficient for 
!                                       clouds in each of the longwave 
!                                       frequency bands [ km**(-1) ]
!
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:                                                  

      real, dimension   &
                  (size(Cloud_microphysics%size_drop,1), &
                   size(Cloud_microphysics%size_drop,2), &
                   size(Cloud_microphysics%size_drop,3), &
                   Cldrad_control%nlwcldb) :: &
                        size_drop, size_ice, conc_drop, conc_ice, &
                        tmpabscoeff
      logical, dimension &
                  (size(Cloud_microphysics%size_drop,1), &
                   size(Cloud_microphysics%size_drop,2), &
                   size(Cloud_microphysics%size_drop,3), &
                   Cldrad_control%nlwcldb) :: dge_column
      logical       :: do_dge_lw, isccp_call
      integer :: nbmax, nbprofiles, nnn, nonly
      integer   :: i, j, k

!----------------------------------------------------------------------
!  local variables:                                                  
!
!      do_dge_lw     logical flag; when .true., indicates that 
!                    donner_deep_mod clouds are being processed
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!    define variable indicating whether doner_deep_mods clouds or 
!    large-scale clouds are currently being processed.
!
!     Donner is never on for this loop
!---------------------------------------------------------------------
      if (Cldrad_control%using_fu2007) then
        do_dge_lw = .true.
      else
        do_dge_lw = .false.
      endif
      isccp_call = .true.
      

!---------------------------------------------------------------------
!    loop over profiles
!---------------------------------------------------------------------
     nonly = ilwband
     nbmax = Cldrad_control%nlwcldb
     nbprofiles = size(Cloud_microphysics%stoch_size_drop,4)
     
     do nnn = 1, nbprofiles

!---------------------------------------------------------------------
!    if large-scale clouds are being processed and stochastic sw
!    clouds has been activated, define the microphysical inputs
!    for each sw parameterization band and the number of such fields.
!---------------------------------------------------------------------
     size_drop(:,:,:,nonly) = Cloud_microphysics%stoch_size_drop(:,:,:,nnn)
     size_ice(:,:,:,nonly)  = Cloud_microphysics%stoch_size_ice(:,:,:,nnn)
     conc_drop(:,:,:,nonly) = Cloud_microphysics%stoch_conc_drop(:,:,:,nnn)
     conc_ice(:,:,:,nonly)  = Cloud_microphysics%stoch_conc_ice(:,:,:,nnn)
        
      if (Cldrad_control%using_fu2007) then
        dge_column = .true.
      else
      do k=1, size(abscoeff,3)
        do j=1, size(abscoeff,2)
          do i=1, size(abscoeff,1)
            if (Cloud_microphysics%stoch_cloud_type(i,j,k,nnn) == 2   &
                        .or. &
                Cloud_microphysics%stoch_cloud_type(i,j,k,nnn) == 3 ) &
                         then 
              dge_column(i,j,k,nonly) = .true.
            else
              dge_column(i,j,k,nonly) = .false.
            endif
          end do
        end do
      end do
      endif
!-------------------------------------------------------------------
!    if the fuliou lw emissivity was selected, call cloud_lwpar to 
!    compute multi-band emissivities based on fu parameterizations for
!    cloud drop, cloud ice, snow and rain. 
!---------------------------------------------------------------------
      if (trim(lwem_form) == 'fuliou') then
      
         
        
        call cloud_lwpar (nonly, nbmax, 1, size_drop, size_ice,   &
                          Cloud_microphysics%size_rain,              &
                          conc_drop, conc_ice, &
                          Cloud_microphysics%conc_rain, &
                          Cloud_microphysics%conc_snow, do_dge_lw,   &
                           dge_column, isccp_call, &
                          tmpabscoeff)
        abscoeff(:,:,:,nnn)=tmpabscoeff(:,:,:,ilwband)                  
!---------------------------------------------------------------------
!    if the ebert-curry emissivity was selected, call cloud_lwem_oneband
!    to compute a single value for the lw emissivity (including effects
!    of drops and ice) based on the ebert and curry parameterization.
!---------------------------------------------------------------------
      else if (trim(lwem_form) == 'ebertcurry') then
        call cloud_lwem_oneband (Cloud_microphysics%conc_drop,   &
                                 Cloud_microphysics%conc_ice,    &
                                 Cloud_microphysics%size_drop,    &
                                 Cloud_microphysics%size_ice,      &
                                 tmpabscoeff(:,:,:,1))
                                 
        abscoeff(:,:,:,nnn)=tmpabscoeff(:,:,:,1)                         
      endif


  enddo        !loop over profiles


!--------------------------------------------------------------------



end subroutine isccp_microphys_lw_driver  



!####################################################################
! <SUBROUTINE NAME="lwemiss_calc">
!  <OVERVIEW>
!   Subroutine to compute infrared emissivity from the absorption
!   coefficient
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute infrared emissivity from the absorption
!   coefficient
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lwemiss_calc(     deltaz, abscoeff, cldemiss)
!  </TEMPLATE>
!  <IN NAME="deltaz" TYPE="real">
!   Pressure layer thickness
!  </IN>
!  <IN NAME="abscoeff" TYPE="real">
!   Absorption coefficient
!  </IN>
!  <OUT NAME="cldemiss" TYPE="real">
!   Emissivity calculated from absorption coefficient
!  </OUT>
! </SUBROUTINE>
!
subroutine lwemiss_calc (deltaz, abscoeff, cldemiss)

!---------------------------------------------------------------------
!    lwemiss_calc computes the infrared emissivity from the absorption 
!    coefficient.
!---------------------------------------------------------------------
 
real, dimension(:,:,:),   intent(in)   :: deltaz  
real, dimension(:,:,:,:,:), intent(in)   :: abscoeff
real, dimension(:,:,:,:,:), intent(out)  :: cldemiss

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      deltaz    depth of the model layer [ meters ]
!      abscoeff  lw absorption coefficient for each of the nlwcldb
!                bands [ km**(-1) ]
!                                                                   
!   intent(out) variables:                                     
!                                                                   
!      cldemiss  the infrared cloud emissivity for each of the nlwcldb
!                bands [ dimensionless ]
!                                                                   
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      integer           :: n, np     ! do-loop index

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!----------------------------------------------------------------------
!    define the cloud emissivity. over a single frequency band 
!    (see goody and yung, eq. 6.72), the emissivity in a layer is 
!    defined as (1 - T(f))    where T(f) is the flux transmissivity, 
!    which may be computed as exp(-(1.66)*(abs. coeff)*
!    (layer thickness)) where the factor 1.66 is the diffusivity factor
!    (diffac). 1.0E-3 is conversion factor from (m) to (km), needed 
!    because abscoeff is in [ km**(-1) ] and deltaz is in [ m ].
!----------------------------------------------------------------------
      do np =1, size(cldemiss,5)
      do n=1,Cldrad_control%nlwcldb
        cldemiss(:,:,:,n,np)  = 1.0E+00 -                           &
                 exp(-diffac*abscoeff(:,:,:,n,np)*deltaz(:,:,:)*1.0E-03)
      end do
      end do
 
!---------------------------------------------------------------------


end subroutine lwemiss_calc



!#################################################################
! <SUBROUTINE NAME="comb_cldprops_calc">
!  <OVERVIEW>
!   Subroutine to define the total-cloud radiative
!    properties
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to define the total-cloud radiative
!    properties to be seen by the radiation package, obtained by the 
!    appropriate combination of the large-scale, donner mesoscale and 
!    cell-scale and uw shallow clouds present in a grid box.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comb_cldprops_calc ( is, js, Rad_time, deltaz,     &
!                               cldext, cldsct, cldasymm, abscoeff, &
!                               Lsc_microphys, Cell_microphys,    &
!                               Meso_microphys, Shallow_microphys, &
!                               Lscrad_props,    &
!                               Cellrad_props, Mesorad_props,  &
!                               Shallowrad_props)
!  </TEMPLATE>
!  <INOUT NAME="cldext" TYPE="real">
!   parameterization band values of the cloud      
!                      extinction coefficient [ km**(-1) ]  
!  </INOUT>
!  <INOUT NAME="cldsct" TYPE="real">
!   parameterization band values of the cloud      
!                      scattering coefficient [ km**(-1) ] 
!  </INOUT>
!  <INOUT NAME="cldasymm" TYPE="real">
!   parameterization band values of the asymmetry  
!                      factor [ dimensionless ]
!  </INOUT>
!  <INOUT NAME="abscoeff" TYPE="real">
!   combined absorption coefficient for clouds in 
!                      each of the longwave frequency bands [ km**(-1) ]
!  </INOUT>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                      clouds
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell 
!                      clouds associated with donner convection
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                      clouds assciated with uw shallow convection
!  </IN>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale 
!                      clouds associated with donner convection
!  </IN>
!  <IN NAME="Shallowrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the  
!                      clouds associated with uw shallow convection
!  </IN>
! </SUBROUTINE>
! 
subroutine comb_cldprops_calc ( is, js, Rad_time, Time_next, deltaz,  &
                               stoch_cloud_type, &
                               cldext, cldsct, cldasymm, abscoeff, &
                               Lsc_microphys, Cell_microphys,    &
                               Meso_microphys, &
                               Shallow_microphys, &
                               Lscrad_props, Cellrad_props,  &
                               Mesorad_props, Shallowrad_props)

!---------------------------------------------------------------------
!    subroutine comb_cldprops_calc defines the total-cloud radiative
!    properties to be seen by the radiation package, obtained by the 
!    appropriate combination of the large-scale, donner mesoscale and 
!    cell-scale, and uw shallow clouds present in a grid box.
!---------------------------------------------------------------------

integer,                intent(in)  :: is, js
type(time_type),        intent(in)  :: Rad_time, Time_next
real, dimension(:,:,:), intent(in)  :: deltaz
integer, dimension(:,:,:,:), intent(in)  :: stoch_cloud_type
real, dimension(:,:,:,:,:),       intent(inout)       ::  cldext,     &
                                                        cldsct, &
                                                        cldasymm
real, dimension(:,:,:,:,:),       intent(inout)       ::  abscoeff
type(microphysics_type),        intent(in), optional :: Lsc_microphys, &
                                                        Cell_microphys,&
                                                        Meso_microphys,&
                                                    Shallow_microphys
type(microrad_properties_type), intent(in), optional :: Lscrad_props, &
                                                        Cellrad_props, &
                                                        Mesorad_props, &
                                                    Shallowrad_props

!--------------------------------------------------------------------
!   intent(inout) variables:
!
!       cldext         parameterization band values of the cloud      
!                      extinction coefficient [ km**(-1) ]   
!       cldsct         parameterization band values of the cloud      
!                      scattering coefficient [ km**(-1) ] 
!       cldasymm       parameterization band values of the asymmetry  
!                      factor [ dimensionless ]
!       abscoeff       combined absorption coefficient for clouds in 
!                      each of the longwave frequency bands [ km**(-1) ]
!
!   intent(in), optional variables:
!
!       Lsc_microphys  microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!       Meso_microphys microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!       Cell_microphys microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!       Shallow_microphys 
!                      microphysical specification for 
!                      clouds associated with uw shallow convection
!                      [ microphysics_type ]
!       Lscrad_props   cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!       Mesorad_props  cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!       Cellrad_props  cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!       Shallowrad_props  cloud radiative properties for 
!                      clouds associated with uw shallow convection  
!                      [ microrad_properties_type ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
 
      real, dimension (size(cldext,1), size(cldext,2),              &
                                    size(cldext,3))    :: cldsum
      real :: cltau,cldextdu
      integer :: i, j, k, n
     
      integer :: nn
  
!------------------------------------------------------------------
!    diagnostics
!------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!        cldsum         total cloud amount in grid box [ dimensionless ]
!        cltau
!        cldextu
!        i,j,k,n        do-loop indices
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!    make sure that appropriate combinations of optional arguments are
!    present.
!---------------------------------------------------------------------
      if (present(Lscrad_props) .or.   &
          present(Lsc_microphys))  then     
        if (.not. present(Lscrad_props) .or.   &
            .not. present(Lsc_microphys))  then     
          call error_mesg ('microphys_rad_mod', &
              ' both Lscrad_props and Lsc_microphys must be present '//&
               'when one is', FATAL)
        endif
      endif
      if (present (Cellrad_props) .or. present (Cell_microphys) .or. &
          present (Mesorad_props) .or. present (Meso_microphys)) then
        if ( .not. present (Cellrad_props) .or.     &
             .not. present (Cell_microphys) .or. &
             .not. present (Mesorad_props) .or.   &
             .not. present (Meso_microphys)) then
          call error_mesg ('microphys_rad_mod', &
            ' either all or none of the cell-scale and meso-scale '//&
              'cloud arguments must be present.', FATAL)
        endif
      endif
      if (present(Shallowrad_props) .or.   &
          present(Shallow_microphys))  then
       if (.not. present(Shallowrad_props) .or.   &
           .not. present(Shallow_microphys))  then
         call error_mesg ('microphys_rad_mod', &
               ' both Shallowrad_props and Shallow_microphys   &
               &must be present when one is', FATAL)
       endif
     endif

!---------------------------------------------------------------------
!    define appropriately-weighted total-cloud radiative properties
!    when large-scale, donner meso-scale and cell-scale, and uw shallow
!    clouds may be present.
!----------------------------------------------------------------------
      if (present(Lscrad_props) .and. present(Cellrad_props) .and. &
          present(Shallowrad_props)) then

!---------------------------------------------------------------------
!    define total cloud fraction.
!---------------------------------------------------------------------
        if (.not. Cldrad_control%do_stochastic_clouds) then
          cldsum = Lsc_microphys%cldamt + Cell_microphys%cldamt +   &
                     Meso_microphys%cldamt + Shallow_microphys%cldamt

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present, 
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*  &
                                     Lscrad_props%cldext(i,j,k,n) + &
                                     Cell_microphys%cldamt(i,j,k)*   &
                                     Cellrad_props%cldext(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*   &
                                     Mesorad_props%cldext(i,j,k,n) + &
                                   Shallow_microphys%cldamt(i,j,k)* &
                                   Shallowrad_props%cldext(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldext(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                exp(-Shallowrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
                  if (cldext(i,j,k,n,1) .gt. cldextdu)   &
                               cldext(i,j,k,n,1)=cldextdu

                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*   &
                                     Lscrad_props%cldsct(i,j,k,n) +   &
                                     Cell_microphys%cldamt(i,j,k)*  &
                                     Cellrad_props%cldsct(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*  &
                                     Mesorad_props%cldsct(i,j,k,n) + &
                                   Shallow_microphys%cldamt(i,j,k)*  &
                              Shallowrad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldsct(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                exp(-Shallowrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                 cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                  cldasymm(i,j,k,n,1) =        &
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n)* &
                           Lscrad_props%cldasymm(i,j,k,n) +&
                           Cell_microphys%cldamt(i,j,k)*  &
                           Cellrad_props%cldsct(i,j,k,n)* &
                           Cellrad_props%cldasymm(i,j,k,n) + &
                           Meso_microphys%cldamt(i,j,k)*   &
                           Mesorad_props%cldsct(i,j,k,n)*  &
                           Mesorad_props%cldasymm(i,j,k,n) + &
                        Shallow_microphys%cldamt(i,j,k)*   &
                        Shallowrad_props%cldsct(i,j,k,n)*  &
                       Shallowrad_props%cldasymm(i,j,k,n) ) /&
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n) +        &
                           Cell_microphys%cldamt(i,j,k)*    &
                           Cellrad_props%cldsct(i,j,k,n) +          &
                           Meso_microphys%cldamt(i,j,k)*    &
                           Mesorad_props%cldsct(i,j,k,n) + &
                        Shallow_microphys%cldamt(i,j,k)*    &
                        Shallowrad_props%cldsct(i,j,k,n) )
                endif
              end do
            end do
          end do
        end do

!------------------------------------------------------------
!    define the total-cloud lw emissivity when large-scale, meso-scale
!    and cell-scale clouds may be present.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%nlwcldb
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu          =                            &
                             (Lsc_microphys%cldamt(i,j,k)*  &
                              Lscrad_props%abscoeff(i,j,k,n) +&
                              Cell_microphys%cldamt(i,j,k)* &
                              Cellrad_props%abscoeff(i,j,k,n) +    &
                              Meso_microphys%cldamt(i,j,k)*   &
                              Mesorad_props%abscoeff(i,j,k,n) +  &
                              Shallow_microphys%cldamt(i,j,k)*   &
                              Shallowrad_props%abscoeff(i,j,k,n)) /   &
                              cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%abscoeff(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Shallowrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
                endif
              end do
            end do
          end do
        end do
      else

!------------------------------------------------------------------------
!    stochastic clouds are being used. we compare the cell and meso-scale
!    cloud amounts to a random number, and replace the large-scale clouds
!    and clear sky in each subcolum with the properties of the cell or 
!    meso-scale clouds when the number is less than the cloud fraction. 
!    we use the maximum overlap assumption. we treat the random number 
!    as the location with the PDF of total water. cells are at the top 
!    of the PDF; then meso-scale anvils, then large-scale clouds and 
!    clear sky. 
!------------------------------------------------------------
      
!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
!       nSubCols = size(Lsc_microphys%stoch_cldamt, 4)
      
!----------------------------------------------------------------------
!    shortwave cloud properties, band by band
!----------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3) ! Levels
            do j=1,size(cldext,2) ! Lons
              do i=1,size(cldext,1) ! Lats
                if (stoch_cloud_type(i,j,k,n) == 3) then
!----------------------------------------------------------------------
!    it's a cell.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Cellrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Cellrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n, 1) = Cellrad_props%cldasymm(i,j,k,n)
                else if (stoch_cloud_type(i,j,k,n) == 2) then
                 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Mesorad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Mesorad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Mesorad_props%cldasymm(i,j,k,n)
                else if (stoch_cloud_type(i,j,k,n) == 4) then
                 
!----------------------------------------------------------------------
!    it's a uw shallow cloud.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Shallowrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Shallowrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Shallowrad_props%cldasymm(i,j,k,n)
                else if (stoch_cloud_type(i,j,k,n) == 1) then
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Lscrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Lscrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Lscrad_props%cldasymm(i,j,k,n)
                else
                  cldext(i,j,k,n,1) = 0. 
                  cldsct(i,j,k,n,1) = 0. 
                  cldasymm(i,j,k,n,1) = 1. 
                endif
              end do 
            end do 
          end do 
        end do 

!----------------------------------------------------------------------
!    longwave cloud properties, band by band
!----------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        nn = Solar_spect%nbands + n
        do k=1,size(cldext,3) ! Levels
          do j=1,size(cldext,2) ! Lons
            do i=1,size(cldext,1) ! Lats
              if (stoch_cloud_type(i,j,k,nn) == 3) then
                 
!----------------------------------------------------------------------
!    it's a cell.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Cellrad_props%abscoeff(i,j,k,n)
              else if (stoch_cloud_type(i,j,k,nn) == 2) then
                 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Mesorad_props%abscoeff(i,j,k,n)
              else if (stoch_cloud_type(i,j,k,nn) == 4) then
                 
!----------------------------------------------------------------------
!    it's a uw shallow.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Shallowrad_props%abscoeff(i,j,k,n)
              else if (stoch_cloud_type(i,j,k,nn) == 1) then
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Lscrad_props%abscoeff(i,j,k,n)
              else
                abscoeff(i,j,k,n,1) = 0. 
              endif
            end do 
          end do 
        end do 
      end do 
      
      

     endif  ! (do_stochastic)

!---------------------------------------------------------------------
!    define appropriately-weighted total-cloud radiative properties
!    when large-scale, meso-scale and cell-scale clouds may be present.
!----------------------------------------------------------------------
      else if (present(Lscrad_props) .and. present(Cellrad_props)) then

!---------------------------------------------------------------------
!    define total cloud fraction.
!---------------------------------------------------------------------
        if (.not. Cldrad_control%do_stochastic_clouds) then
          cldsum = Lsc_microphys%cldamt + Cell_microphys%cldamt +   &
                     Meso_microphys%cldamt

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present, 
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*  &
                                     Lscrad_props%cldext(i,j,k,n) + &
                                     Cell_microphys%cldamt(i,j,k)*   &
                                     Cellrad_props%cldext(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*   &
                                     Mesorad_props%cldext(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldext(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
                  if (cldext(i,j,k,n,1) .gt. cldextdu)   &
                               cldext(i,j,k,n,1)=cldextdu

                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*   &
                                     Lscrad_props%cldsct(i,j,k,n) +   &
                                     Cell_microphys%cldamt(i,j,k)*  &
                                     Cellrad_props%cldsct(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*  &
                                     Mesorad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldsct(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                  cldasymm(i,j,k,n,1) =        &
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n)* &
                           Lscrad_props%cldasymm(i,j,k,n) +&
                           Cell_microphys%cldamt(i,j,k)*  &
                           Cellrad_props%cldsct(i,j,k,n)* &
                           Cellrad_props%cldasymm(i,j,k,n) + &
                           Meso_microphys%cldamt(i,j,k)*   &
                           Mesorad_props%cldsct(i,j,k,n)*  &
                           Mesorad_props%cldasymm(i,j,k,n)) /&
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n) +        &
                           Cell_microphys%cldamt(i,j,k)*    &
                           Cellrad_props%cldsct(i,j,k,n) +          &
                           Meso_microphys%cldamt(i,j,k)*    &
                           Mesorad_props%cldsct(i,j,k,n) )
                endif
              end do
            end do
          end do
        end do

!------------------------------------------------------------
!    define the total-cloud lw emissivity when large-scale, meso-scale
!    and cell-scale clouds may be present.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%nlwcldb
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu          =                            &
                             (Lsc_microphys%cldamt(i,j,k)*  &
                              Lscrad_props%abscoeff(i,j,k,n) +&
                              Cell_microphys%cldamt(i,j,k)* &
                              Cellrad_props%abscoeff(i,j,k,n) +    &
                              Meso_microphys%cldamt(i,j,k)*   &
                              Mesorad_props%abscoeff(i,j,k,n)) /   &
                              cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%abscoeff(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
                endif
              end do
            end do
          end do
        end do
      else

        if (do_orig_donner_stoch) then

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present,
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
       do n=1,Solar_spect%nbands
         cldsum = Lsc_microphys%sw_stoch_cldamt(:,:,:,n) +   &
                  Cell_microphys%cldamt + Meso_microphys%cldamt
         do k=1,size(cldext,3)
           do j=1,size(cldext,2)
             do i=1,size(cldext,1)
               if (cldsum(i,j,k) > 0.0) then
                 cldextdu        = (  &
                              Lsc_microphys%sw_stoch_cldamt(i,j,k,n)*  &
                                  Lscrad_props%cldext(i,j,k,n) +   &
                                     Cell_microphys%cldamt(i,j,k)*   &
                                      Cellrad_props%cldext(i,j,k,n) +  &
                                      Meso_microphys%cldamt(i,j,k)*   &
                                      Mesorad_props%cldext(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                            exp(-Cellrad_props%cldext(i,j,k,n)*     &
                           deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                          cldsum(i,j,k))* &
                    exp(-Mesorad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                               1000.)
                 cltau=cltau+(Lsc_microphys%sw_stoch_cldamt(i,j,k,n)/  &
                            cldsum(i,j,k))* &
                     exp(-Lscrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
                 if (cldext(i,j,k,n,1) .gt. cldextdu)   &
                               cldext(i,j,k,n,1)=cldextdu

          cldextdu        = (Lsc_microphys%sw_stoch_cldamt(i,j,k,n)*   &
                                    Lscrad_props%cldsct(i,j,k,n) +   &
                                     Cell_microphys%cldamt(i,j,k)*  &
                                      Cellrad_props%cldsct(i,j,k,n) +  &
                                      Meso_microphys%cldamt(i,j,k)*  &
                                    Mesorad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                        exp(-Cellrad_props%cldsct(i,j,k,n)*     &
                           deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                        cldsum(i,j,k))* &
                  exp(-Mesorad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                             1000.)
               cltau=cltau+(Lsc_microphys%sw_stoch_cldamt(i,j,k,n)/  &
                          cldsum(i,j,k))* &
                   exp(-Lscrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                          1000.)
                cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                 cldasymm(i,j,k,n,1) =        &
                          (Lsc_microphys%sw_stoch_cldamt(i,j,k,n)*  &
                           Lscrad_props%cldsct(i,j,k,n)* &
                          Lscrad_props%cldasymm(i,j,k,n) +&
                           Cell_microphys%cldamt(i,j,k)*  &
                          Cellrad_props%cldsct(i,j,k,n)* &
                        Cellrad_props%cldasymm(i,j,k,n) + &
                           Meso_microphys%cldamt(i,j,k)*   &
                          Mesorad_props%cldsct(i,j,k,n)*  &
                                Mesorad_props%cldasymm(i,j,k,n)) /&
                       (Lsc_microphys%sw_stoch_cldamt(i,j,k,n)*  &
                            Lscrad_props%cldsct(i,j,k,n) +        &
                           Cell_microphys%cldamt(i,j,k)*    &
                           Cellrad_props%cldsct(i,j,k,n) +          &
                           Meso_microphys%cldamt(i,j,k)*    &
                            Mesorad_props%cldsct(i,j,k,n) )
                endif
             end do
            end do
          end do
        end do

!---------------------------------------------------------------------
!    define the total-cloud lw emissivity when large-scale, meso-scale
!    and cell-scale clouds may be present.
!---------------------------------------------------------------------
       do n=1,Cldrad_control%nlwcldb
         cldsum = Lsc_microphys%lw_stoch_cldamt(:,:,:,n) +   &
                  Cell_microphys%cldamt + Meso_microphys%cldamt
         do k=1,size(cldext,3)
           do j=1,size(cldext,2)
             do i=1,size(cldext,1)
               if (cldsum(i,j,k) > 0.0) then
                 cldextdu          =                            &
                      (Lsc_microphys%lw_stoch_cldamt(i,j,k,n)*  &
                               Lscrad_props%abscoeff(i,j,k,n) +&
                               Cell_microphys%cldamt(i,j,k)* &
                              Cellrad_props%abscoeff(i,j,k,n) +    &
                              Meso_microphys%cldamt(i,j,k)*   &
                             Mesorad_props%abscoeff(i,j,k,n)) /   &
                               cldsum(i,j,k)
                 cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                   exp(-Cellrad_props%abscoeff(i,j,k,n)*             &
                       deltaz(i,j,k)/1000.)
                cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                          cldsum(i,j,k))* &
         exp(-Mesorad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/          &
                           1000.)
      cltau=cltau+(Lsc_microphys%lw_stoch_cldamt(i,j,k,n)/  &
                           cldsum(i,j,k))* &
           exp(-Lscrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/          &
                          1000.)
              abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
       if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
               endif
          end do
           end do
         end do
       end do 
       
        else  ! (using new donner-stochastic connection)
!------------------------------------------------------------------------
!    stochastic clouds are being used. we compare the cell and meso-scale
!    cloud amounts to a random number, and replace the large-scale clouds
!    and clear sky in each subcolum with the properties of the cell or 
!    meso-scale clouds when the number is less than the cloud fraction. 
!    we use the maximum overlap assumption. we treat the random number 
!    as the location with the PDF of total water. cells are at the top 
!    of the PDF; then meso-scale anvils, then large-scale clouds and 
!    clear sky. 
!------------------------------------------------------------
      
!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
!       nSubCols = size(Lsc_microphys%stoch_cldamt, 4)

!----------------------------------------------------------------------
!    shortwave cloud properties, band by band
!----------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3) ! Levels
            do j=1,size(cldext,2) ! Lons
              do i=1,size(cldext,1) ! Lats
                if ( stoch_cloud_type(i,j,k,n) == 3) then 
!----------------------------------------------------------------------
!    it's a cell.
!----------------------------------------------------------------------
              IF (ignore_donner_cells) then
                  cldext(i,j,k,n,1) = 0. 
                  cldsct(i,j,k,n,1) = 0. 
                  cldasymm(i,j,k,n,1) = 1. 
              ELSE
                  cldext(i,j,k,n,1) = Cellrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Cellrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n, 1) = Cellrad_props%cldasymm(i,j,k,n)
              ENDIF
                else if ( stoch_cloud_type(i,j,k,n) == 2) then 
                 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Mesorad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Mesorad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Mesorad_props%cldasymm(i,j,k,n)
                else if ( stoch_cloud_type(i,j,k,n) == 1) then 
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Lscrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Lscrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Lscrad_props%cldasymm(i,j,k,n)
                else
                  cldext(i,j,k,n,1) = 0. 
                  cldsct(i,j,k,n,1) = 0. 
                  cldasymm(i,j,k,n,1) = 1. 
                endif
              end do 
            end do 
          end do 
        end do 

!----------------------------------------------------------------------
!    longwave cloud properties, band by band
!----------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        nn = Solar_spect%nbands + n
        do k=1,size(cldext,3) ! Levels
          do j=1,size(cldext,2) ! Lons
            do i=1,size(cldext,1) ! Lats
                if ( stoch_cloud_type(i,j,k,nn) == 3) then 
              IF (ignore_donner_cells) then
                abscoeff(i,j,k,n,1) = 0.                              
              ELSE
                 
!----------------------------------------------------------------------
!    it's a cell.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Cellrad_props%abscoeff(i,j,k,n)
               ENDIF
                else if ( stoch_cloud_type(i,j,k,nn) == 2) then 
                 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Mesorad_props%abscoeff(i,j,k,n)
                else if ( stoch_cloud_type(i,j,k,nn) == 1) then 
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Lscrad_props%abscoeff(i,j,k,n)
              else
                abscoeff(i,j,k,n,1) = 0. 
              endif
            end do 
          end do 
        end do 
      end do 
      

      endif ! (do_orig_donner_stoch)
     endif  ! (do_stochastic)

!---------------------------------------------------------------------
!    define appropriately-weighted total-cloud radiative properties
!    when large-scale, and uw shallow clouds may be present.
!----------------------------------------------------------------------
     else if (present(Lscrad_props) .and. present(Shallowrad_props)) then

!---------------------------------------------------------------------
!    define total cloud fraction.
!---------------------------------------------------------------------
        if (.not. Cldrad_control%do_stochastic_clouds) then
          cldsum = Lsc_microphys%cldamt + Shallow_microphys%cldamt 

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present, 
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*  &
                                     Lscrad_props%cldext(i,j,k,n) + &
                                   Shallow_microphys%cldamt(i,j,k)*   &
                                   Shallowrad_props%cldext(i,j,k,n) )/ &
                                     cldsum(i,j,k)
                cltau=(Shallow_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                           exp(-Shallowrad_props%cldext(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
                  if (cldext(i,j,k,n,1) .gt. cldextdu)   &
                               cldext(i,j,k,n,1)=cldextdu

                  cldextdu        = (Lsc_microphys%cldamt(i,j,k)*   &
                                     Lscrad_props%cldsct(i,j,k,n) +   &
                                  Shallow_microphys%cldamt(i,j,k)*  &
                                  Shallowrad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                cltau=(Shallow_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                            exp(-Shallowrad_props%cldsct(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                  cldasymm(i,j,k,n,1) =        &
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n)* &
                           Lscrad_props%cldasymm(i,j,k,n) +&
                        Shallow_microphys%cldamt(i,j,k)*   &
                        Shallowrad_props%cldsct(i,j,k,n)*  &
                        Shallowrad_props%cldasymm(i,j,k,n)) /&
                          (Lsc_microphys%cldamt(i,j,k)*  &
                           Lscrad_props%cldsct(i,j,k,n) +        &
                        Shallow_microphys%cldamt(i,j,k)*    &
                        Shallowrad_props%cldsct(i,j,k,n) )
                endif
              end do
            end do
          end do
        end do

!------------------------------------------------------------
!    define the total-cloud lw emissivity when large-scale, meso-scale
!    and cell-scale clouds may be present.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%nlwcldb
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu          =                            &
                             (Lsc_microphys%cldamt(i,j,k)*  &
                              Lscrad_props%abscoeff(i,j,k,n) +&
                        Shallow_microphys%cldamt(i,j,k)*   &
                           Shallowrad_props%abscoeff(i,j,k,n)) /   &
                              cldsum(i,j,k)
                cltau=(Shallow_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                         exp(-Shallowrad_props%abscoeff(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Lsc_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Lscrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
                endif
              end do
            end do
          end do
        end do
      else

!------------------------------------------------------------------------
!    stochastic clouds are being used. we compare the cell and meso-scale
!    cloud amounts to a random number, and replace the large-scale clouds
!    and clear sky in each subcolum with the properties of the cell or 
!    meso-scale clouds when the number is less than the cloud fraction. 
!    we use the maximum overlap assumption. we treat the random number 
!    as the location with the PDF of total water. cells are at the top 
!    of the PDF; then meso-scale anvils, then large-scale clouds and 
!    clear sky. 
!------------------------------------------------------------
      
!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
!       nSubCols = size(Lsc_microphys%stoch_cldamt, 4)

!----------------------------------------------------------------------
!    shortwave cloud properties, band by band
!----------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3) ! Levels
            do j=1,size(cldext,2) ! Lons
              do i=1,size(cldext,1) ! Lats
                if ( stoch_cloud_type(i,j,k,n) == 4) then 
!----------------------------------------------------------------------
!    it's a uw shallow
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Shallowrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Shallowrad_props%cldsct(i,j,k,n)
              cldasymm(i,j,k,n, 1) = Shallowrad_props%cldasymm(i,j,k,n)
                else if ( stoch_cloud_type(i,j,k,n) == 1) then 
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                  cldext(i,j,k,n,1) = Lscrad_props%cldext(i,j,k,n)
                  cldsct(i,j,k,n,1) = Lscrad_props%cldsct(i,j,k,n)
                  cldasymm(i,j,k,n,1) = Lscrad_props%cldasymm(i,j,k,n)
                else
                  cldext(i,j,k,n,1) = 0. 
                  cldsct(i,j,k,n,1) = 0. 
                  cldasymm(i,j,k,n,1) = 1. 
                endif
              end do 
            end do 
          end do 
        end do 

!----------------------------------------------------------------------
!    longwave cloud properties, band by band
!----------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        nn = Solar_spect%nbands + n
        do k=1,size(cldext,3) ! Levels
          do j=1,size(cldext,2) ! Lons
            do i=1,size(cldext,1) ! Lats
                if ( stoch_cloud_type(i,j,k,nn) == 4) then 
                 
!----------------------------------------------------------------------
!    it's a uw shallow.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Shallowrad_props%abscoeff(i,j,k,n)
                else if ( stoch_cloud_type(i,j,k,nn) == 1) then 
                 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                abscoeff(i,j,k,n,1) = Lscrad_props%abscoeff(i,j,k,n)
              else
                abscoeff(i,j,k,n,1) = 0. 
              endif
            end do 
          end do 
        end do 
      end do 
      

     endif  ! (do_stochastic)

!---------------------------------------------------------------------
!    define appropriately-weighted total-cloud radiative properties
!    when donner deep convective clouds and uw shallow clouds
!    may be present.
!----------------------------------------------------------------------
     else if (present(Cellrad_props) .and. present(Shallowrad_props)) then

!---------------------------------------------------------------------
!    define total cloud fraction.
!---------------------------------------------------------------------
        cldsum = Shallow_microphys%cldamt + Cell_microphys%cldamt +   &
                     Meso_microphys%cldamt

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present, 
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu        = (Shallow_microphys%cldamt(i,j,k)*  &
                                   Shallowrad_props%cldext(i,j,k,n) + &
                                     Cell_microphys%cldamt(i,j,k)*   &
                                     Cellrad_props%cldext(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*   &
                                     Mesorad_props%cldext(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldext(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                 exp(-Shallowrad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
                  if (cldext(i,j,k,n,1) .gt. cldextdu)   &
                               cldext(i,j,k,n,1)=cldextdu

               cldextdu        = (Shallow_microphys%cldamt(i,j,k)*   &
                                  Shallowrad_props%cldsct(i,j,k,n) +   &
                                     Cell_microphys%cldamt(i,j,k)*  &
                                     Cellrad_props%cldsct(i,j,k,n) +  &
                                     Meso_microphys%cldamt(i,j,k)*  &
                                     Mesorad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldsct(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Shallowrad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                  cldasymm(i,j,k,n,1) =        &
                          (Shallow_microphys%cldamt(i,j,k)*  &
                           Shallowrad_props%cldsct(i,j,k,n)* &
                           Shallowrad_props%cldasymm(i,j,k,n) +&
                           Cell_microphys%cldamt(i,j,k)*  &
                           Cellrad_props%cldsct(i,j,k,n)* &
                           Cellrad_props%cldasymm(i,j,k,n) + &
                           Meso_microphys%cldamt(i,j,k)*   &
                           Mesorad_props%cldsct(i,j,k,n)*  &
                           Mesorad_props%cldasymm(i,j,k,n)) /&
                          (Shallow_microphys%cldamt(i,j,k)*  &
                           Shallowrad_props%cldsct(i,j,k,n) +        &
                           Cell_microphys%cldamt(i,j,k)*    &
                           Cellrad_props%cldsct(i,j,k,n) +          &
                           Meso_microphys%cldamt(i,j,k)*    &
                           Mesorad_props%cldsct(i,j,k,n) )
                endif
              end do
            end do
          end do
        end do

!------------------------------------------------------------
!    define the total-cloud lw emissivity when large-scale, meso-scale
!    and cell-scale clouds may be present.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%nlwcldb
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu          =                            &
                             (Shallow_microphys%cldamt(i,j,k)*  &
                              Shallowrad_props%abscoeff(i,j,k,n) +&
                              Cell_microphys%cldamt(i,j,k)* &
                              Cellrad_props%abscoeff(i,j,k,n) +    &
                              Meso_microphys%cldamt(i,j,k)*   &
                              Mesorad_props%abscoeff(i,j,k,n)) /   &
                              cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%abscoeff(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                   cltau=cltau+(Shallow_microphys%cldamt(i,j,k)/  &
                           cldsum(i,j,k))* &
                    exp(-Shallowrad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                           1000.)
                 abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
                endif
              end do
            end do
          end do
        end do
!--------------------------------------------------------------------
!    define the total-cloud radiative properties when only meso-scale 
!    and cell-scale clouds may be present.
!--------------------------------------------------------------------
      else if (present(Cellrad_props)) then

!---------------------------------------------------------------------
!    define total cloud fraction.
!---------------------------------------------------------------------
        cldsum = Cell_microphys%cldamt + Meso_microphys%cldamt

!---------------------------------------------------------------------
!     define the cloud scattering, cloud extinction and cloud asymmetry
!     factor in each of the spectral bands. if cloud is not present, 
!     values remain at the non-cloudy initialized values.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu        = (Cell_microphys%cldamt(i,j,k)* &
                                     Cellrad_props%cldext(i,j,k,n) +   &
                                     Meso_microphys%cldamt(i,j,k)*  &
                                     Mesorad_props%cldext(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldext(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldext(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                 cldext(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldext(i,j,k,n,1) .gt. cldextdu) cldext(i,j,k,n,1)=cldextdu

                  cldextdu        = (Cell_microphys%cldamt(i,j,k)*  &
                                     Cellrad_props%cldsct(i,j,k,n) +   &
                                     Meso_microphys%cldamt(i,j,k)*   &
                                     Mesorad_props%cldsct(i,j,k,n)) / &
                                     cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%cldsct(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%cldsct(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                 cldsct(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (cldsct(i,j,k,n,1) .gt. cldextdu) cldsct(i,j,k,n,1)=cldextdu

                  cldasymm(i,j,k,n,1) =                           &
                          (Cell_microphys%cldamt(i,j,k)* &
                           Cellrad_props%cldsct(i,j,k,n)* &
                           Cellrad_props%cldasymm(i,j,k,n) + &
                           Meso_microphys%cldamt(i,j,k)*   &
                           Mesorad_props%cldsct(i,j,k,n)*  &
                           Mesorad_props%cldasymm(i,j,k,n)) /&
                          (Cell_microphys%cldamt(i,j,k)*  &
                           Cellrad_props%cldsct(i,j,k,n) +     &
                           Meso_microphys%cldamt(i,j,k)*   &
                           Mesorad_props%cldsct(i,j,k,n) )
                endif     
              end do
            end do
          end do
        end do

!---------------------------------------------------------------------
!    define the total-cloud lw emissivity when only meso-scale and
!    cell-scale clouds may be present.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%nlwcldb
          do k=1,size(cldext,3)
            do j=1,size(cldext,2)
              do i=1,size(cldext,1)
                if (cldsum(i,j,k) > 0.0) then
                  cldextdu          =                           &
                       (Cell_microphys%cldamt(i,j,k)*  &
                        Cellrad_props%abscoeff(i,j,k,n) +   &
                        Meso_microphys%cldamt(i,j,k)*   &
                        Mesorad_props%abscoeff(i,j,k,n)) /   &
                        cldsum(i,j,k)
                   cltau=(Cell_microphys%cldamt(i,j,k)/cldsum(i,j,k))* &
                               exp(-Cellrad_props%abscoeff(i,j,k,n)*     &
                          deltaz(i,j,k)/1000.)
                   cltau=cltau+(Meso_microphys%cldamt(i,j,k)/  &
                         cldsum(i,j,k))* &
                   exp(-Mesorad_props%abscoeff(i,j,k,n)*deltaz(i,j,k)/  &
                              1000.)
                 abscoeff(i,j,k,n,1)=-1000.*alog(cltau)/deltaz(i,j,k)
        if (abscoeff(i,j,k,n,1) .gt. cldextdu) abscoeff(i,j,k,n,1)=cldextdu
                endif
              end do
            end do
          end do
        end do
      endif


!---------------------------------------------------------------------



end subroutine comb_cldprops_calc




!#################################################################
! <SUBROUTINE NAME="microphys_rad_end">
!  <OVERVIEW>
!   microphys_rad_end is the destructor for microphys_rad_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   microphys_rad_end is the destructor for microphys_rad_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_rad_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine microphys_rad_end

!-------------------------------------------------------------------
!    microphys_rad_end is the destructor for microphys_rad_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg('microphys_rad_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    deallocate module variables.
!--------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        deallocate ( nivl1liqcld    ,   &
                     nivl1icecld    ,   &
                     nivl1icesolcld ,   &
                     nivl1raincld   ,   &
                     nivl1snowcld   ,   &
                     nivl2liqcld    ,   &
                     nivl2icecld    ,   &
                     nivl2icesolcld ,   &
                     nivl2raincld   ,   & 
                     nivl2snowcld   ,   &
                     solivlicecld   ,   &
                     solivlicesolcld,   &
                     solivlliqcld   ,   &
                     solivlraincld  ,   &
                     solivlsnowcld  )
      endif

!--------------------------------------------------------------------
!    mark the module as no longer being initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine microphys_rad_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!#################################################################
! <SUBROUTINE NAME="cloudpar">
!  <OVERVIEW>
!   Subroutine to determine cloud single scattering parameters
!  </OVERVIEW>
!  <DESCRIPTION>
!   determine the parameterization band values of the single scattering   
! parameters (extinction coefficient, scattering coefficient and   
! asymmetry factor) for clouds from the size and/or concentration of    
! each constituent (cloud drops, rain drops, ice crystals and snow)     
! present.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudpar (   nonly, nbmax, nnn, &
!                     size_drop, size_ice, size_rain, conc_drop,   &
!                     conc_ice, conc_rain, conc_snow, do_dge_sw,   &
!                     cldext, cldsct, cldasymm)
!  </TEMPLATE>
!  <IN NAME="nonly" TYPE="integer">
!   The single band for calculations.  Note that this is used
!   only in the case of a call from cloudrad_diagnostics to 
!   do isccp simulator work.  For all other calls, nonly should
!   be 0 and will have no effect on the calculations below
!  </IN>
!  <IN NAME="nbmax" TYPE="integer">
!   The number of individual bands to do calculations over. Note
!   that for normal GCM calls this will be 1.  For calls using
!   stochastic clouds with or without the isccp simulator this will
!   be equal to the number of shortwave bands
!  </IN>
!  <IN NAME="nnn" TYPE="integer">
!   This integer controls which cloud state to access for radiation
!   calculations.  For normal GCM applications this will be 1. For
!   Full Independent Column Approximation calculations with stochast-
!   ic clouds this will be the profile number being accessed. 
!  </IN>
!  <IN NAME="conc_drop" TYPE="real">
!   total cloud droplet concentration
!  </IN>
!  <IN NAME="conc_ice" TYPE="real">
!   ice crystal concentration
!  </IN>
!  <IN NAME="conc_rain" TYPE="real">
!   rain droplet concetration
!  </IN>
!  <IN NAME="conc_snow" TYPE="real">
!   snow concentration
!  </IN>
!  <IN NAME="size_drop" TYPE="real">
!   cloud droplet size distribution
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   ice crystal size distribution
!  </IN>
!  <IN NAME="size_rain" TYPE="real">
!   rain drop size distribution
!  </IN>
!  <IN NAME="do_dge_sw" TYPE="logical">
!   use sw parameterizations using generalized effective 
!                  size developed by Fu et al (1998) (if true). 
!                  otherwise use parameterizations by Fu et al using 
!                  effective size.
!  </IN>
!  <OUT NAME="cldext" TYPE="real">
!   the parameterization band values of the cloud      
!               extinction coefficient in kilometer**(-1)
!  </OUT>
!  <OUT NAME="cldsct" TYPE="real">
!   the parameterization band values of the cloud      
!               scattering coefficient in kilometer**(-1)
!  </OUT>
!  <OUT NAME="cldasymm" TYPE="real">
!   the parameterization band values of the asymmetry  
!               factor
!  </OUT>
! </SUBROUTINE>
subroutine cloudpar (nonly, nbmax, nnn, size_drop, size_ice, size_rain, & 
                     conc_drop, conc_ice, conc_rain, conc_snow, do_dge_sw, &
                     dge_column, isccp_call, cldext, cldsct, cldasymm)
 
!----------------------------------------------------------------------
!    subroutine cloudpar determines the parameterization band values of
!    the single scattering parameters (extinction coefficient, scatter-
!    ing coefficient and asymmetry factor) for clouds from the size and/
!    or concentration of each constituent (cloud drops, rain drops, ice
!    crystals and snow) present.                          
!----------------------------------------------------------------------

integer,                   intent(in)     ::  nonly, nbmax, nnn
real, dimension (:,:,:),   intent(in)     ::  size_rain, conc_rain,   &
                                              conc_snow
real, dimension (:,:,:,:), intent(in)     ::  size_drop, size_ice,    &
                                              conc_drop, conc_ice
logical, dimension (:,:,:,:), intent(in)     ::  dge_column
logical,                   intent(in)     ::  do_dge_sw
logical,                   intent(in)     ::  isccp_call
real, dimension (:,:,:,:), intent(inout)  ::  cldext, cldsct, cldasymm
 
!-------------------------------------------------------------------
! intent(in) variables:                                            
!                                                                  
!       size_drop  the cloud drop effective diameter [ microns ]    
!       size_ice   the ice crystal effective size  [ microns ]     
!       size_rain  the rain drop effective diameter [ microns ]    
!       conc_drop  the cloud drop liquid water concentration 
!                  [ grams / meter**3 ]                            
!       conc_ice   the ice water concentation 
!                  [ grams / meter**3 ]                            
!       conc_rain  the rain drop water concentration 
!                  [ grams / meter**3 ]                            
!       conc_snow  the snow concentration 
!                  [ grams / meter**3 ]                            
!       do_dge_sw  use sw parameterizations using generalized effective 
!                  size developed by Fu et al (1998) (if true). 
!                  otherwise use parameterizations by Fu et al using 
!                  effective size.
! 
! intent(inout) variables:                                           
!                                                                  
!       cldext     the parameterization band values of the cloud      
!                  extinction coefficient [ kilometer**(-1) ]          
!       cldsct     the parameterization band values of the cloud      
!                  scattering coefficient [ kilometer**(-1) ]         
!       cldasymm   the parameterization band values of the asymmetry  
!                  factor  [ dimensionless ]                    
!----------------------------------------------------------------------
 
!---------------------------------------------------------------------
! local variables:                                                   
!--------------------------------------------------------------------
      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), NLIQCLDIVLS)  ::   &
                    cldextivlliq, cldssalbivlliq, cldasymmivlliq

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), NRAINCLDIVLS)  ::   &
                    cldextivlrain, cldssalbivlrain, cldasymmivlrain

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), NSNOWCLDIVLS)  ::   &
                    cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), NICECLDIVLS)  ::   &
                    cldextivlice, cldssalbivlice, cldasymmivlice

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), NICESOLARCLDIVLS)  ::   &
                    cldextivlice2, cldssalbivlice2, cldasymmivlice2

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), Solar_spect%nbands)  ::    &
                    cldextbandliq, cldssalbbandliq, cldasymmbandliq,&
                    cldextbandice, cldssalbbandice, cldasymmbandice,&
                    cldextbandrain, cldssalbbandrain, cldasymmbandrain,&
                    cldextbandsnow, cldssalbbandsnow, cldasymmbandsnow

      logical, dimension (size(conc_drop,1), size(conc_drop,2), &
                          size(conc_drop,3))  ::   maskl, anymask, &
                                                   maskr, maski, masks,&
                                                   maskif, maskis
      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                          size(conc_drop,3))  ::   tempext, tempext2, &
                                     tempssa, tempssa2, &
                                             tempasy, tempasy2

      integer  :: nb
      integer  :: i,j,k
      real :: sum, sum2, sum3

!----------------------------------------------------------------------
!   local variables:
!
!      cldextivlliq     cloud extinction coefficient over the spectral
!                       intervals relevant to cloud droplets 
!                       [ km**(-1)]
!      cldssalbivlliq   cloud single scattering albedo over the spectral
!                       intervals relevant to cloud droplets 
!                       [ non-dimensional ]
!      cldasymmivlliq   cloud asymmetry factor over the spectral
!                       intervals relevant to cloud droplets 
!                       [ non-dimensional ]
!      cldextivlrain    cloud extinction coefficient over the spectral
!                       intervals relevant to rain drops [ km**(-1)]  
!      cldssalbivlrain  cloud single scattering albedo over the spectral
!                       intervals relevant to rain drops 
!                       [ non-dimensional ]
!      cldasymmivlrain  cloud asymmetry factor over the spectral
!                       intervals relevant to rain drops         
!                       [ non-dimensional ]
!      cldextivlsnow    cloud extinction coefficient over the spectral
!                       intervals relevant to snow flakes  [ km**(-1)] 
!      cldssalbivlsnow  cloud single scattering albedo over the spectral
!                       intervals relevant to snow flakes 
!                       [ non-dimensional ]
!      cldasymmivlsnow  cloud asymmetry factor over the spectral
!                       intervals relevant to snow flakes       
!                       [ non-dimensional ]
!      cldextivlice     cloud extinction coefficient over the spectral
!                       intervals relevant to fu (1996) ice crystals 
!                       [ km**(-1)]
!      cldssalbivlice   cloud single scattering albedo over the spectral
!                       intervals relevant to fu(1996) ice crystals 
!                       [ non-dimensional ]
!      cldasymmivlice   cloud asymmetry factor over the spectral
!                       intervals relevant to fu(1996) ice crystals
!                       [ non-dimensional ]
!      cldextivlice2    cloud extinction coefficient over the spectral
!                       intervals relevant to fu and liou(1993) ice
!                       crystals  [ km**(-1)]
!      cldssalbivlice2  cloud single scattering albedo over the spectral
!                       intervals relevant to fu and liou(1993) ice
!                       crystals  [ non-dimensional ]
!      cldasymmivlice2  cloud asymmetry factor over the spectral
!                       intervals relevant to fu and liou(1993) ice
!                       crystals
!                       [ non-dimensional ]
!      cldextbandliq    cloud extinction coefficient for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud droplets  [ km**(-1)]  
!      cldssalbbandliq  cloud single scattering albedo for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud droplets
!                       [ non-dimensional ]
!      cldasymmbandliq  cloud asymmetry factor for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud droplets [ non-dimensional ]
!      cldextbandice    cloud extinction coefficient for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud ice  [ km**(-1)]  
!      cldssalbbandice  cloud single scattering albedo for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud ice 
!                       [ non-dimensional ]
!      cldasymmbandice  cloud asymmetry factor for each spectral
!                       parameterization bands resulting from the 
!                       presence of cloud ice [ non-dimensional ]
!      cldextbandrain   cloud extinction coefficient for each spectral
!                       parameterization bands resulting from the 
!                       presence of rain drops  [ km**(-1)]  
!      cldssalbbandrain cloud single scattering albedo for each spectral
!                       parameterization bands resulting from the 
!                       presence of rain drops  
!                       [ non-dimensional ]
!      cldasymmbandrain cloud asymmetry factor for each spectral
!                       parameterization bands resulting from the 
!                       presence of rain drops [ non-dimensional ]
!      cldextbandsnow   cloud extinction coefficient for each spectral
!                       parameterization bands resulting from the 
!                       presence of snow flakes  [ km**(-1)]  
!                       [ non-dimensional ]
!      cldssalbbandsnow cloud single scattering albedo for each spectral
!                       parameterization bands resulting from the 
!                       presence of snow flakes 
!                       [ non-dimensional ]
!      cldasymmbandsnow cloud asymmetry factor for each spectral
!                       parameterization bands resulting from the 
!                       presence of snow flakes [ non-dimensional ]
!      nb               do-loop index
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!
!  NOTE THE FOLLOWING LOGICAL TO THE LOOPS BELOW
!
!        do nb=1,nbmax
!             if (nbmax==1) then
!                  call slingo(conc_drop(:,:,:,nnn).....)
!             else
!                  call slingo(conc_drop(:,:,:,nb),....)
!             end if
!
!        enddo           
!             
!        Note that nbmax = 1 in the following cases:
!
!                 (a) standard GCM applications which do not use
!                     McICA  
!                 (b) Full Independent Column Approximation 
!                     calculations
!
!        Note that nbmax = Solar_spect%nbands in the following cases
!               
!                 (c) McICA calculations where nb = the cloud 
!                     profile being used
!                 (d) ISCCP simulator calls from cloudrad_diagnostics
!                     where "nonly" will be used
!
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!    call slingo to define the single scattering parameters for cloud 
!    droplets for each of the slingo cloud droplet spectral intervals.
!----------------------------------------------------------------------
      do nb=1,nbmax
        if (nbmax == 1) then
          call slingo (conc_drop(:,:,:,nnn), size_drop(:,:,:,nnn),  &
                       cldextivlliq, cldssalbivlliq, cldasymmivlliq, &
                       maskl)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for cloud
!    droplets that were calculated for each cloud droplet spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
          call thickavg (nivl1liqcld, nivl2liqcld, NLIQCLDIVLS,     &
                         Solar_spect%nbands, cldextivlliq,   &
                         cldssalbivlliq, cldasymmivlliq, solivlliqcld, &
                         Solar_spect%solflxbandref, maskl, &
                         cldextbandliq, &
                         cldssalbbandliq, cldasymmbandliq)

!----------------------------------------------------------------------
!    call savijarvi to define the single scattering parameters for 
!    rain drops for each of the savijarvi rain drop spectral intervals. 
!----------------------------------------------------------------------
          call savijarvi (conc_rain, size_rain, cldextivlrain,   &
                      cldssalbivlrain, cldasymmivlrain,maskr) 

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for rain 
!    drops that were calculated for each rain drop spectral interval 
!    to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
          call thickavg (nivl1raincld, nivl2raincld, NRAINCLDIVLS,  &
                     Solar_spect%nbands, cldextivlrain,    &
                     cldssalbivlrain , cldasymmivlrain,  &
                     solivlraincld, Solar_spect%solflxbandref, maskr,   &
                     cldextbandrain, cldssalbbandrain,  &
                     cldasymmbandrain)
 
!---------------------------------------------------------------------
!    on calls from the isccp simulator with stochastic clouds, call all 
!    parameterizations, since some subcolumns may have cloud types 
!    using that parameterization. calculation control is provided 
!    within the parameterization subroutine.
!---------------------------------------------------------------------
         if (isccp_call) then

!----------------------------------------------------------------------
!    define the single scattering parameters for ice crystals.        
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    call the ice crystal parameterization scheme of fu et al(1998) 
!    using generalized effective size. call subroutine fu to calculate 
!    the single scattering parameters.
!----------------------------------------------------------------------
            call fu (conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),   &
                     dge_column(:,:,:,nnn), &
                     cldextivlice, cldssalbivlice,  &
                     cldasymmivlice, maski)
 
!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
            call thickavg (nivl1icecld, nivl2icecld, NICECLDIVLS,   &
                           Solar_spect%nbands, cldextivlice,     &
                           cldssalbivlice, cldasymmivlice,     &
                           solivlicecld, Solar_spect%solflxbandref,   &
                           maski, &
                           cldextbandice, cldssalbbandice,   &
                           cldasymmbandice)

!----------------------------------------------------------------------
!    call the ice crystal parameterization scheme of fu et al(1993) 
!    using effective size. call subroutine icesolar to calculate
!    the single scattering parameters.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1993) using
!    effective size is to be used, call subroutine icesolar to calculate
!    the single scattering parameters.
!----------------------------------------------------------------------
            call icesolar (conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn), &
                           dge_column(:,:,:,nnn), &
                           cldextivlice2, cldssalbivlice2,   &
                           cldasymmivlice2, maski)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
            call thickavg (nivl1icesolcld, nivl2icesolcld,    &
                           NICESOLARCLDIVLS, Solar_spect%nbands, &
                           cldextivlice2, cldssalbivlice2,  &
                           cldasymmivlice2, solivlicesolcld,  &
                           Solar_spect%solflxbandref, maski,  &
                           cldextbandice,  cldssalbbandice,  &
                           cldasymmbandice)

!--------------------------------------------------------------------
!    if this is not an isccp call with activated stochastic clouds, all
!    grid columns will use the same parameterization.
!--------------------------------------------------------------------
         else  !(isccp_call)
           if (do_dge_sw) then

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1998) using
!    generalized effective size is to be used, call subroutine fu 
!    to calculate the single scattering parameters.
!----------------------------------------------------------------------
            call fu (conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),   &
                     dge_column(:,:,:,nnn), &
                     cldextivlice, cldssalbivlice,  &
                     cldasymmivlice, maski)
 
!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
            call thickavg (nivl1icecld, nivl2icecld, NICECLDIVLS,   &
                           Solar_spect%nbands, cldextivlice,     &
                           cldssalbivlice, cldasymmivlice,     &
                           solivlicecld, Solar_spect%solflxbandref,   &
                           maski, &
                           cldextbandice, cldssalbbandice,   &
                           cldasymmbandice)

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1993) using
!    effective size is to be used, call subroutine icesolar to calculate
!    the single scattering parameters.
!----------------------------------------------------------------------
          else

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1993) using
!    effective size is to be used, call subroutine icesolar to calculate
!    the single scattering parameters.
!----------------------------------------------------------------------
            call icesolar (conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn), &
                           dge_column(:,:,:,nnn), &
                           cldextivlice2, cldssalbivlice2,   &
                           cldasymmivlice2, maski)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
            call thickavg (nivl1icesolcld, nivl2icesolcld,    &
                           NICESOLARCLDIVLS, Solar_spect%nbands, &
                           cldextivlice2, cldssalbivlice2,  &
                           cldasymmivlice2, solivlicesolcld,  &
                           Solar_spect%solflxbandref, maski,  &
                           cldextbandice,  cldssalbbandice,  &
                           cldasymmbandice)
          endif
        endif !(isccp_call)

!----------------------------------------------------------------------
!    define the single scattering parameters for snow.                
!----------------------------------------------------------------------
          call snowsw (conc_snow, cldextivlsnow, cldssalbivlsnow,    &
                       cldasymmivlsnow, masks)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for each snow flake spectral interval 
!    to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
          call thickavg (nivl1snowcld, nivl2snowcld, NSNOWCLDIVLS,  &
                         Solar_spect%nbands, cldextivlsnow,     &
                         cldssalbivlsnow , cldasymmivlsnow,  &
                         solivlsnowcld, Solar_spect%solflxbandref,  &
                         masks, &
                         cldextbandsnow, cldssalbbandsnow,  &
                         cldasymmbandsnow)
 
        else
          if (nonly.eq.0   .or. nonly.eq.nb ) then
            call slingo (conc_drop(:,:,:,nb), size_drop(:,:,:,nb),  &
                       cldextivlliq, cldssalbivlliq, cldasymmivlliq,  &
                       maskl,  &
                       starting_band = nivl1liqcld(nb), &
                       ending_band = nivl2liqcld(nb))

!----------------------------------------------------------------------
!    call savijarvi to define the single scattering parameters for 
!    rain drops for each of the savijarvi rain drop spectral intervals. 
!----------------------------------------------------------------------
            call savijarvi (conc_rain, size_rain, cldextivlrain,   &
                            cldssalbivlrain, cldasymmivlrain, &
                            maskr,  &
                            starting_band = nivl1raincld(nb), &
                            ending_band = nivl2raincld(nb)) 

!---------------------------------------------------------------------
!    on calls from the isccp simulator with stochastic clouds, call all 
!    parameterizations, since some subcolumns may have cloud types 
!    using that parameterization. calculation control is provided 
!    within the parameterization subroutine.
!---------------------------------------------------------------------
         if (isccp_call) then

!----------------------------------------------------------------------
!    define the single scattering parameters for ice crystals.        
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    call the ice crystal parameterization scheme of fu et al(1998) 
!    using generalized effective size is to be used, call subroutine fu 
!    to calculate the single scattering parameters.
!----------------------------------------------------------------------
              call fu (conc_ice(:,:,:,nb), size_ice(:,:,:,nb),   &
                        dge_column(:,:,:,nb), &
                       cldextivlice, cldssalbivlice,  &
                     cldasymmivlice,   maskif, &
                     starting_band = nivl1icecld(nb), &
                     ending_band = nivl2icecld(nb))

!----------------------------------------------------------------------
!    call the ice crystal parameterization scheme of fu et al(1993) 
!    using effective size is to be used, call subroutine icesolar to 
!    calculate the single scattering parameters.
!----------------------------------------------------------------------
              call icesolar (conc_ice(:,:,:,nb), size_ice(:,:,:,nb), &
                           dge_column(:,:,:,nb), &
                           cldextivlice2,     &
                           cldssalbivlice2, cldasymmivlice2, &
                           maskis, &
                           starting_band = nivl1icesolcld(nb), &
                           ending_band = nivl2icesolcld(nb))

!--------------------------------------------------------------------
!    if this is not an isccp call with activated stochastic clouds, all
!    grid columns will use the same parameterization.
!--------------------------------------------------------------------
        else   ! (isccp_call)
            if (do_dge_sw) then

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1998) using
!    generalized effective size is to be used, call subroutine fu 
!    to calculate the single scattering parameters.
!----------------------------------------------------------------------
              maskis = .false.
              call fu (conc_ice(:,:,:,nb), size_ice(:,:,:,nb),   &
                        dge_column(:,:,:,nb), &
                       cldextivlice, cldssalbivlice,  &
                     cldasymmivlice,   maskif, &
                     starting_band = nivl1icecld(nb), &
                     ending_band = nivl2icecld(nb))
            else

!----------------------------------------------------------------------
!    if the ice crystal parameterization scheme of fu et al(1993) using
!    effective size is to be used, call subroutine icesolar to calculate
!    the single scattering parameters.
!----------------------------------------------------------------------
              maskif = .false.
              call icesolar (conc_ice(:,:,:,nb), size_ice(:,:,:,nb), &
                           dge_column(:,:,:,nb), &
                           cldextivlice2,     &
                           cldssalbivlice2, cldasymmivlice2, &
                           maskis, &
                           starting_band = nivl1icesolcld(nb), &
                           ending_band = nivl2icesolcld(nb))
            endif
         endif  ! (isccp_call)

!----------------------------------------------------------------------
!    define the single scattering parameters for snow.                
!----------------------------------------------------------------------
            call snowsw (conc_snow, cldextivlsnow, cldssalbivlsnow,    &
                         cldasymmivlsnow, &
                           masks, &
                         starting_band = nivl1snowcld(nb), &
                         ending_band = nivl2snowcld(nb)) 
          
!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for cloud
!    droplets that were calculated for each cloud droplet spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
            if (nonly == 0 ) then
              call thickavg (nb, nivl1liqcld(nb), nivl2liqcld(nb),  &
                            NLIQCLDIVLS,   &
                        Solar_spect%nbands, cldextivlliq,    &
                        cldssalbivlliq, cldasymmivlliq, solivlliqcld, &
                        Solar_spect%solflxbandref(nb),  maskl, &
                        cldextbandliq(:,:,:,nb),   &
                        cldssalbbandliq(:,:,:,nb),  &
                        cldasymmbandliq(:,:,:,nb))

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for rain 
!    drops that were calculated for each rain drop spectral interval 
!    to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              call thickavg ( nb, nivl1raincld(nb), nivl2raincld(nb), &
                              NRAINCLDIVLS,  &
                     Solar_spect%nbands, cldextivlrain,    &
                     cldssalbivlrain , cldasymmivlrain,  &
                     solivlraincld, Solar_spect%solflxbandref(nb),  &
                     maskr,&
                     cldextbandrain(:,:,:,nb),   &
                     cldssalbbandrain(:,:,:,nb), &
                     cldasymmbandrain(:,:,:,nb))

!---------------------------------------------------------------------
!    on calls from the isccp simulator with stochastic clouds, call all 
!    parameterizations, since some subcolumns may have cloud types 
!    using that parameterization. calculation control is provided 
!    within the parameterization subroutine.
!---------------------------------------------------------------------
           if (isccp_call) then

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
                call thickavg (nb, nivl1icecld(nb), nivl2icecld(nb), &
                               NICECLDIVLS, &
                           Solar_spect%nbands, cldextivlice,     &
                           cldssalbivlice, cldasymmivlice,     &
                           solivlicecld,   &
                           Solar_spect%solflxbandref(nb),  &
                           maskif, &
                           tempext, tempssa, tempasy)   
                call thickavg (nb, nivl1icesolcld(nb),  &
                              nivl2icesolcld(nb),    &
                           NICESOLARCLDIVLS, Solar_spect%nbands, &
                           cldextivlice2, cldssalbivlice2,  &
                           cldasymmivlice2,&
                           solivlicesolcld,  &
                           Solar_spect%solflxbandref(nb),&
                           maskis, &
                           tempext2, tempssa2, tempasy2)   

              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    if (maskif(i,j,k)) then
                      cldextbandice(i,j,k,nb) = tempext(i,j,k)
                       cldssalbbandice(i,j,k,nb) = tempssa(i,j,k)
                       cldasymmbandice(i,j,k,nb) = tempasy(i,j,k)
                    else if (maskis(i,j,k)) then
                      cldextbandice(i,j,k,nb) = tempext2(i,j,k)
                       cldssalbbandice(i,j,k,nb) = tempssa2(i,j,k)
                       cldasymmbandice(i,j,k,nb) = tempasy2(i,j,k)
                    endif
                  end do
                  end do
                  end do
!--------------------------------------------------------------------
!    if this is not an isccp call with activated stochastic clouds, all
!    grid columns will use the same parameterization.
!--------------------------------------------------------------------
  else ! (isccp_call)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for ice
!    crystals that were calculated for each ice crystal spectral
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              if (do_dge_sw) then
                maskis = .false.
                call thickavg (nb, nivl1icecld(nb), nivl2icecld(nb), &
                               NICECLDIVLS, &
                           Solar_spect%nbands, cldextivlice,     &
                           cldssalbivlice, cldasymmivlice,     &
                           solivlicecld,   &
                           Solar_spect%solflxbandref(nb),  &
                           maskif, &
                           cldextbandice(:,:,:,nb),   &
                            cldssalbbandice(:,:,:,nb),   &
                           cldasymmbandice(:,:,:,nb))
              else
                maskif = .false.
                call thickavg (nb, nivl1icesolcld(nb),  &
                              nivl2icesolcld(nb),    &
                           NICESOLARCLDIVLS, Solar_spect%nbands, &
                           cldextivlice2, cldssalbivlice2,  &
                           cldasymmivlice2,&
                           solivlicesolcld,  &
                           Solar_spect%solflxbandref(nb),&
                           maskis, &
                           cldextbandice(:,:,:,nb),  &
                           cldssalbbandice(:,:,:,nb),  &
                           cldasymmbandice(:,:,:,nb))
              endif

  endif ! (isccp_call)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for each snow flake spectral interval 
!    to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              call thickavg (nb, nivl1snowcld(nb), nivl2snowcld(nb), &
                            NSNOWCLDIVLS,  &
                     Solar_spect%nbands, cldextivlsnow,     &
                     cldssalbivlsnow , cldasymmivlsnow,  &
                     solivlsnowcld,  &
                     Solar_spect%solflxbandref(nb), masks, &
                     cldextbandsnow(:,:,:,nb),  &
                     cldssalbbandsnow(:,:,:,nb),  &
                     cldasymmbandsnow(:,:,:,nb))

              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    anymask(i,j,k) = maskl(i,j,k) .or. maskr(i,j,k)  &
                                    .or. maskif(i,j,k) .or.  &
                                        maskis(i,j,k) .or. masks(i,j,k)
                  end do
                end do
              end do

              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    if (anymask(i,j,k)) then
                      sum =0.
                      sum2 =0.
                      sum3 =0.
                      if (maskl(i,j,k)) then
                        sum = sum + cldextbandliq(i,j,k,nb)
                        sum2 = sum2 + cldextbandliq(i,j,k,nb)* &
                                      cldssalbbandliq(i,j,k,nb)
                        sum3 = sum3 + (cldextbandliq(i,j,k,nb)* &
                                       cldssalbbandliq(i,j,k,nb))* &
                                       cldasymmbandliq(i,j,k,nb)
                      endif
                      if (maskr(i,j,k)) then
                        sum = sum + cldextbandrain(i,j,k,nb)
                        sum2 = sum2 + cldextbandrain(i,j,k,nb)* &
                                      cldssalbbandrain(i,j,k,nb)
                        sum3 = sum3 + (cldextbandrain(i,j,k,nb)*&
                                       cldssalbbandrain(i,j,k,nb))* &
                                       cldasymmbandrain(i,j,k,nb)
                      endif
                      if (maskis(i,j,k) .or. maskif(i,j,k)) then
                        sum = sum + cldextbandice(i,j,k,nb)
                        sum2 = sum2 + cldextbandice(i,j,k,nb)* &
                                      cldssalbbandice(i,j,k,nb)
                        sum3 = sum3 + (cldextbandice(i,j,k,nb)*&
                                       cldssalbbandice(i,j,k,nb))*&
                                       cldasymmbandice(i,j,k,nb)
                      endif
                      if (masks(i,j,k)) then
                        sum = sum +  cldextbandsnow(i,j,k,nb)
                        sum2 = sum2 + cldextbandsnow(i,j,k,nb)* &
                                      cldssalbbandsnow(i,j,k,nb)
                        sum3 = sum3 + (cldextbandsnow(i,j,k,nb)*&
                                       cldssalbbandsnow(i,j,k,nb))* &
                                       cldasymmbandsnow(i,j,k,nb)
                      endif
                      cldext(i,j,k,nb) = sum
                      cldsct(i,j,k,nb) = sum2
                      cldasymm(i,j,k,nb) = sum3/ (cldsct(i,j,k,nb) +  &
                                                  1.0e-100)
                    else
                      cldext(i,j,k,nb) = 0.0
                      cldsct(i,j,k,nb) = 0.0
                      cldasymm(i,j,k,nb) = 0.0
                    endif
                  end do
                end do
              end do

            else

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for the desired cloud drop spectral 
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              call thickavg (nb, nivl1liqcld(nb), nivl2liqcld(nb),&
                            cldextivlliq,    &
                            solivlliqcld, &
                           Solar_spect%solflxbandref(nb), maskl,  &
                               cldextbandliq(:,:,:,nb))

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for the desired rain drop spectral 
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              call thickavg (nonly, nivl1raincld(nonly),  &
                             nivl2raincld(nonly),  &
                                         cldextivlrain,    &
                     solivlraincld, Solar_spect%solflxbandref(nonly), &
                     maskr, cldextbandrain(:,:,:,nonly))

!---------------------------------------------------------------------
!    on calls from the isccp simulator with stochastic clouds, call all 
!    parameterizations, since some subcolumns may have cloud types 
!    using that parameterization. calculation control is provided 
!    within the parameterization subroutine.
!---------------------------------------------------------------------
           if (isccp_call) then

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for the desired ice crystal spectral 
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
                call thickavg (nb, nivl1icecld(nb), nivl2icecld(nb),  &
                                               cldextivlice,     &
                          solivlicecld,   &
                          Solar_spect%solflxbandref(nb),    &
                           maskif, tempext                )
                call thickavg (nb, nivl1icesolcld(nb),  &
                               nivl2icesolcld(nb),    &
                           cldextivlice2,   &
                       solivlicesolcld,  &
                            Solar_spect%solflxbandref(nb),    &
                           maskis, tempext2               )
              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    if (maskif(i,j,k)) then
                      cldextbandice(i,j,k,nb) = tempext(i,j,k)
                    else if (maskis(i,j,k)) then
                      cldextbandice(i,j,k,nb) = tempext2(i,j,k)
                    endif
                  end do
                  end do
                  end do
                    
!--------------------------------------------------------------------
!    if this is not an isccp call with activated stochastic clouds, all
!    grid columns will use the same parameterization.
!--------------------------------------------------------------------
        else  ! (isccp_call)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for the desired ice crystal spectral 
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              if (do_dge_sw) then
                maskis = .false.
                call thickavg (nb, nivl1icecld(nb), nivl2icecld(nb),  &
                                               cldextivlice,     &
                          solivlicecld,   &
                          Solar_spect%solflxbandref(nb),    &
                           maskif, cldextbandice(:,:,:,nb))
              else
                maskif = .false.
                call thickavg (nb, nivl1icesolcld(nb),  &
                               nivl2icesolcld(nb),    &
                           cldextivlice2,   &
                       solivlicesolcld,  &
                            Solar_spect%solflxbandref(nb),    &
                           maskis, cldextbandice(:,:,:,nb))
              endif
  endif ! (isccp_call)

!----------------------------------------------------------------------
!    call thickavg to map the single-scattering properties for snow 
!    flakes that were calculated for the desired snow flake spectral 
!    interval to the sw parameterization band spectral intervals.
!----------------------------------------------------------------------
              call thickavg (nonly, nivl1snowcld(nonly),   &
                             nivl2snowcld(nonly), &
                                         cldextivlsnow,     &
                  solivlsnowcld, Solar_spect%solflxbandref(nonly),   &
                     masks, cldextbandsnow(:,:,:,nonly))

              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    anymask(i,j,k) = maskl(i,j,k) .or.  &
                                     maskr(i,j,k) .or.  &
                                      maskif(i,j,k) .or. &
                                      maskis(i,j,k) .or. masks(i,j,k)
                  end do
                end do
              end do

              do k=1, size(cldextbandliq,3)
                do j=1, size(cldextbandliq,2)
                  do i=1, size(cldextbandliq,1)
                    if (anymask(i,j,k)) then
                      sum =0.
                      if (maskl(i,j,k)) then
                        sum = sum + cldextbandliq(i,j,k,nonly)
                      endif
                      if (maskr(i,j,k)) then
                        sum = sum + cldextbandrain(i,j,k,nonly)
                      endif
                      if (maskif(i,j,k) .or. maskis(i,j,k)) then
                        sum = sum + cldextbandice(i,j,k,nonly)
                      endif
                      if (masks(i,j,k)) then
                        sum = sum +  cldextbandsnow(i,j,k,nonly)
                      endif
                      cldext(i,j,k,nonly) = sum
                    else
                      cldext(i,j,k,nonly) = 0.0
                    endif
                  end do
                end do
              end do
            endif
          endif !for nonly              
        endif !for (nbmax == 1)
      end do

!----------------------------------------------------------------------
!    combine the contribution to the single-scattering properties from
!    each of the individual constituents to define the overall cloud 
!    values in each sw parameterization band.                
!----------------------------------------------------------------------
      if (nbmax == 1) then
        cldext   =  cldextbandliq + cldextbandrain +  &
                    cldextbandice + cldextbandsnow
        cldsct   =  cldssalbbandliq*cldextbandliq  + &
                    cldssalbbandrain*cldextbandrain  + &
                    cldssalbbandice*cldextbandice +  &
                    cldssalbbandsnow*cldextbandsnow
        cldasymm = (cldasymmbandliq*(cldssalbbandliq*cldextbandliq) + &
                    cldasymmbandrain*                                 &
                              (cldssalbbandrain*cldextbandrain) +     &
                    cldasymmbandice*(cldssalbbandice*cldextbandice) + &
                    cldasymmbandsnow*                                 &
                              (cldssalbbandsnow*cldextbandsnow))/     &
                                                     (cldsct + 1.0e-100)
      endif 



!---------------------------------------------------------------------


end subroutine cloudpar



!#####################################################################
! <SUBROUTINE NAME="slingo">
!  <OVERVIEW>
!   Subroutine to determine single scattering parameters for clouds
!  </OVERVIEW>
!  <DESCRIPTION>
!   define the single scattering parameters for cloud drops using the     
! Slingo parameterization for his spectral intervals.
!   slingo, a., a gcm parameterization of the shortwave properties of     
!      water clouds., j. atmos. sci.,46, 1419-1427, 1989.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call slingo                                               &
!                    (conc_drop   , size_drop     ,                    &
!                     cldextivlliq, cldssalbivlliq, cldasymmivlliq, &
!                     starting_band, ending_band )
!  </TEMPLATE>
!  <IN NAME="conc_drop" TYPE="real">
!   the cloud drop liquid water concentration in grams meter**3 
!  </IN>
!  <IN NAME="size_drop" TYPE="real">
!   the cloud drop effective diameter in microns
!  </IN>
!  <OUT NAME="cldextivlliq" TYPE="real">
!   The specified spectral values of the extinction      
!   coefficient in kilometer**(-1) for drops
!  </OUT>
!  <OUT NAME="cldssalbivlliq" TYPE="real">
!   the specified spectral values of the single-scattering albedo 
!   for drops
!  </OUT>
!  <OUT NAME="cldasymmivlliq" TYPE="real">
!   the specified spectral values of the asymmetry factor for drops
!  </OUT>
!  <IN NAME="starting_band">
!
!  </IN>
!  <IN NAME="ending_band">
!
!  </IN>
! </SUBROUTINE>
!
subroutine slingo (conc_drop, size_drop, cldextivlliq, cldssalbivlliq, &
                   cldasymmivlliq, mask, starting_band, ending_band)
 
!----------------------------------------------------------------------
!    subroutine slingo defines the single scattering parameters for 
!    cloud droplets using the Slingo parameterization for his spectral 
!    intervals. references:                                      
!    slingo, a., a gcm parameterization of the shortwave properties of 
!                water clouds., j. atmos. sci.,46, 1419-1427, 1989.   
!----------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)     ::   conc_drop, size_drop
real, dimension (:,:,:,:), intent(inout)  ::   cldextivlliq,  &
                                               cldssalbivlliq,   &
                                               cldasymmivlliq
logical, dimension(:,:,:), intent(out)    ::   mask
integer, intent(in), optional             ::   starting_band,  &
                                               ending_band

!-------------------------------------------------------------------
!   intent(in) variables:                                              
!                                                                       
!        conc_drop        cloud drop liquid water concentration 
!                         [ grams / meter**3 ]                       
!        size_drop        cloud drop effective diameter [ microns ]    
!
!    intent(out) variables:                                     
!                                                                       
!        cldextivlliq     extinction coefficient in each spectral
!                         interval of the slingo cloud droplet param-
!                         eterization resulting from the presence of 
!                         cloud droplets  [ km **(-1) ]
!        cldssalbivlliq   single scattering albedo in each spectral
!                         interval of the slingo cloud droplet param-
!                         eterization resulting from the presence of 
!                         cloud droplets  [ dimensionless ]
!        cldasymmivlliq   asymmetry factor in each spectral
!                         interval of the slingo cloud droplet param-
!                         eterization resulting from the presence of 
!                         cloud droplets  [ dimensionless ]
!
!   intent(in), optional variables:
!
!        starting_band    the index of the first droplet spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!        ending_band      the index of the last droplet spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!
!----------------------------------------------------------------------
 
!---------------------------------------------------------------------
! local variables:                                                   

      real, dimension (size(conc_drop,1), size(conc_drop,2),  &
                       size(conc_drop,3))   ::  size_d

      real, dimension (NLIQCLDIVLS)         ::  a, b, c, d, e, f
 
      data a /-1.023E+00, 1.950E+00, 1.579E+00, 1.850E+00, 1.970E+00, &
               2.237E+00, 2.463E+00, 2.551E+00, 2.589E+00, 2.632E+00, &
               2.497E+00, 2.622E+00, 2.650E+00, 3.115E+00, 2.895E+00, &
               2.831E+00, 2.838E+00, 2.672E+00, 2.698E+00, 2.668E+00, &
               2.801E+00, 3.308E+00, 2.944E+00, 3.094E+00 /
      data b / 1.933E+00, 1.540E+00, 1.611E+00, 1.556E+00, 1.501E+00, &
               1.452E+00, 1.420E+00, 1.401E+00, 1.385E+00, 1.365E+00, &
               1.376E+00, 1.362E+00, 1.349E+00, 1.244E+00, 1.315E+00, &
               1.317E+00, 1.300E+00, 1.320E+00, 1.315E+00, 1.307E+00, &
               1.293E+00, 1.246E+00, 1.270E+00, 1.252E+00 /
      data c / 2.500E-02, 4.490E-01, 1.230E-01, 1.900E-04, 1.200E-03, &
               1.200E-04, 2.400E-04, 6.200E-05,-2.800E-05,-4.600E-05, &
               9.800E-06, 3.300E-06, 2.300E-06,-2.700E-07,-1.200E-07, &
              -1.200E-06, 0.000E+00, 0.000E+00, 1.000E-06, 0.000E+00, &
               1.000E-06,-3.000E-07,-6.500E-07, 7.900E-07 /
      data d / 1.220E-02, 1.540E-03, 9.350E-03, 2.540E-03, 2.160E-03, &
               6.670E-04, 8.560E-04, 2.600E-04, 8.000E-05, 5.000E-05, &
               2.100E-05, 2.800E-06, 1.700E-06, 1.400E-06, 4.400E-07, &
               4.000E-07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
               0.000E+00, 2.360E-07, 4.330E-07, 3.690E-07 /
      data e / 7.260E-01, 8.310E-01, 8.510E-01, 7.690E-01, 7.400E-01, &
               7.490E-01, 7.540E-01, 7.730E-01, 7.800E-01, 7.840E-01, &
               7.830E-01, 8.060E-01, 8.090E-01, 8.040E-01, 8.180E-01, &
               8.280E-01, 8.250E-01, 8.280E-01, 8.200E-01, 8.400E-01, &
               8.360E-01, 8.390E-01, 8.410E-01, 8.440E-01 /
      data f / 6.652E+00, 6.102E+00, 2.814E+00, 5.171E+00, 7.469E+00, &
               6.931E+00, 6.555E+00, 5.405E+00, 4.989E+00, 4.745E+00, &
               5.035E+00, 3.355E+00, 3.387E+00, 3.520E+00, 2.989E+00, &
               2.492E+00, 2.776E+00, 2.467E+00, 3.004E+00, 1.881E+00, &
               2.153E+00, 1.946E+00, 1.680E+00, 1.558E+00 /

      integer   :: nistart, niend
      integer   :: i, j, k, ni
 
!---------------------------------------------------------------------
! local variables:                                                   
!
!      size_d       droplet effective radius [ microns ]
!      a            slingo parameterization coefficient for cloud
!                   extinction coefficient [ m**2 / g ]
!      b            slingo parameterization coefficient for cloud
!                   extinction coefficient [ micron*m**2 / g ]
!      c            slingo parameterization coefficient for cloud
!                   single scattering albedo [ nondimensional ]
!      d            slingo parameterization coefficient for cloud
!                   single scattering albedo [ micron **(-1) ]
!                   asymmetry factor [ nondimensional ]
!      f            slingo parameterization coefficient for cloud
!                   asymmetry factor [ micron **(-1) ]
!      nistart      first droplet parameterization band to be processed
!      niend        last droplet parameterization band to be processed
!      i,j,k,ni     do-loop indices
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the starting and ending droplet parameterization bands to
!    be processed during this call.
!--------------------------------------------------------------------
      if (present(starting_band)) then
        nistart = starting_band
      else
        nistart = 1
      endif
      if (present(ending_band)) then
        niend = ending_band
      else
        niend = NLIQCLDIVLS
      endif
 
!---------------------------------------------------------------------
      do k=1,size(conc_drop,3)
        do j=1,size(conc_drop,2)
          do i=1,size(conc_drop,1)

!------------------------------------------------------------------
!    bypass calculations if no drops are present. values are set to
!    zero in all spectral bands.
!-----------------------------------------------------------------
            if (conc_drop(i,j,k) == 0.0) then
              mask(i,j,k) = .false.

!--------------------------------------------------------------------
!    convert input variable size from diameter to radius for use in the
!    slingo formulae.
!--------------------------------------------------------------------
            else
              mask(i,j,k) = .true.
              size_d(i,j,k) = 0.5*size_drop(i,j,k)

!----------------------------------------------------------------------
!    the cloud drop effective radius must be between 4.2 and 16.6 
!    microns.                               
!----------------------------------------------------------------------
              if (size_d(i,j,k) <  min_cld_drop_rad) then   
                size_d(i,j,k) =  min_cld_drop_rad
              else if (size_d(i,j,k) > max_cld_drop_rad) then 
                size_d(i,j,k) = max_cld_drop_rad             
              endif

!---------------------------------------------------------------------
!    define values of extinction coefficient, single-scattering albedo
!    and asymmetry factor for each of the slingo parameterization 
!    spectral bands. these values are a function of droplet concen-
!    tration and droplet effective radius. the extinction coefficient 
!    is converted to kilometer**(-1).     
!---------------------------------------------------------------------
                do ni=nistart, niend
                  cldextivlliq(i,j,k,ni) = 1.0E+03*conc_drop(i,j,k)* &
                                           (1.0E-02*a(ni) + (b(ni)/  &
                                           size_d(i,j,k)            ) )
                  cldssalbivlliq(i,j,k,ni) = 1.0 - ( c(ni) + d(ni)* &
                                             size_d(i,j,k) )
                  cldasymmivlliq(i,j,k,ni) = e(ni) + 1.0E-03*f(ni)*  &
                                             size_d(i,j,k)
                end do
            endif     
          end do
        end do
      end do

!-------------------------------------------------------------------


end subroutine slingo




!#####################################################################
! <SUBROUTINE NAME="savijarvi">
!  <OVERVIEW>
!   Subroutine to define the single scattering parameters for rain drop
!  </OVERVIEW>
!  <DESCRIPTION>
!   define the single scattering parameters for rain drops using the      
! Savijarvi parameterization for his spectral intervals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call savijarvi                                          &
!                    (conc_rain    , size_rain      ,                &
!                     cldextivlrain, cldssalbivlrain, cldasymmivlrain)
!  </TEMPLATE>
!  <IN NAME="conc_rain" TYPE="real">
!   the rain drop water concentration in grams / meter**3
!  </IN>
!  <IN NAME="size_rain" TYPE="real">
!   the rain drop effective diameter in microns
!  </IN>
!  <OUT NAME="cldextivlrain" TYPE="real">
!   the specified spectral values of the extinction     
!                   coefficient for rain in kilometers**(-1)
!  </OUT>
!  <OUT NAME="cldssalbivlrain" TYPE="real">
!   the specified spectral values of the single-        
!                   scattering albedo for rain
!  </OUT>
!  <OUT NAME="cldasymmivlrain" TYPE="real">
!   the specified spectral values of the asymmetry      
!                   factor for rain
!  </OUT>
! </SUBROUTINE>
!
subroutine savijarvi (conc_rain, size_rain, cldextivlrain,    &
                      cldssalbivlrain, cldasymmivlrain, &
                      mask, starting_band, ending_band)      
 
!----------------------------------------------------------------------
!    subroutine savijarvi defines the single scattering parameters for 
!    rain drops using the Savijarvi parameterization for his spectral 
!    intervals. references:                                     
!    savijarvi, h., shortwave optical properties of rain., tellus, 49a, 
!    177-181, 1997.                                                   
!---------------------------------------------------------------------- 

real, dimension (:,:,:),   intent(in)   ::  conc_rain, size_rain
real, dimension (:,:,:,:), intent(out)  ::  cldextivlrain,            &
                                            cldssalbivlrain,   &
                                            cldasymmivlrain
logical, dimension(:,:,:), intent(out)    ::   mask
integer, intent(in), optional             ::   starting_band,  &
                                               ending_band

!---------------------------------------------------------------------
!  intent(in) variables:
!
!        conc_rain        rain drop water concentration [ grams / m**3 ]
!        size_rain        rain drop effective diameter [ microns ]     
!
!  intent(out) variables:
!
!        cldextivlrain    extinction coefficient in each spectral
!                         interval of the savijarvi rain drop param-
!                         eterization resulting from the presence of 
!                         rain drops  [ km **(-1) ]
!        cldssalbivlrain  single scattering albedo in each spectral
!                         interval of the savijarvi rain drop param-
!                         eterization resulting from the presence of 
!                         rain drops  [ dimensionless ]
!        cldasymmivlrain  asymmetry factor in each spectral
!                         interval of the savijarvi rain drop param-
!                         eterization resulting from the presence of 
!                         rain drops  [ dimensionless ]
!
!---------------------------------------------------------------------
 
!---------------------------------------------------------------------- 
! local variables:                                                      
 
      real, dimension (size(conc_rain,1), size(conc_rain,2),       &
                       size(conc_rain,3) )       ::                &
                                                      rcap, size_d
 
      real, dimension (NRAINCLDIVLS)          ::  a, b, asymm

      data a     / 4.65E-01, 2.64E-01, 1.05E-02, 8.00E-05 /
      data b     / 1.00E-03, 9.00E-02, 2.20E-01, 2.30E-01 /
      data asymm / 9.70E-01, 9.40E-01, 8.90E-01, 8.80E-01 /

      integer   ::  i, j, k, ni
      integer   ::  nistart, niend

!---------------------------------------------------------------------
!   local variables:
! 
!         rcap       drop size function used in savijarvi parameter-
!                    ization : (drop radius/500.)**4.348
!                    [ dimensionless ]
!         size_d     rain drop effective radius [ microns ]
!         a          interval-dependent parameter in savijarvi single 
!                    scattering albedo formula
!                    [ dimensionless ]
!         b          interval-dependent parameter in savijarvi single 
!                    scattering albedo formula
!                    [ dimensionless ]
!         asymm      asymmetry factor for each savijarvi spectral band
!         i,j,k,ni   do-loop indices
!--------------------------------------------------------------------
!---------------------------------------------------------------------
!    define the starting and ending droplet parameterization bands to
!    be processed during this call.
!--------------------------------------------------------------------
      if (present(starting_band)) then
        nistart = starting_band
      else
        nistart = 1
      endif
      if (present(ending_band)) then
        niend = ending_band
      else
        niend = NRAINCLDIVLS
      endif

!--------------------------------------------------------------------
      do k=1,size(conc_rain,3)
        do j=1,size(conc_rain,2)
          do i=1,size(conc_rain,1)
 
!-----------------------------------------------------------------
!    if no rain is present in a grid box, set the scattering parameters
!    to so indicate.
!-----------------------------------------------------------------
            if (conc_rain(i,j,k) == 0.0) then
              mask(i,j,k) = .false.

!----------------------------------------------------------------------
!    convert input size from drop diameter to drop radius. 
!----------------------------------------------------------------------
            else
              mask(i,j,k) = .true.
              size_d(i,j,k) = 0.5*size_rain(i,j,k) 

!---------------------------------------------------------------------
!    the rain drop effective radius must be between 16.6 and 5000    
!    microns. compute the rcap function, used in the savijarvi formula.
!---------------------------------------------------------------------
              if (size_d(i,j,k) > 16.6 .and.              &
                  size_d(i,j,k) <= 5000. ) then                       
                rcap(i,j,k) = (size_d(i,j,k)/500.) ** 4.348E+00

!--------------------------------------------------------------------
!    compute values for each of the savijarvi rain drop spectral
!    intervals. the extinction coefficient is converted to km**(-1).    
!--------------------------------------------------------------------
!               do ni = 1,NRAINCLDIVLS
                do ni = nistart, niend 
                  cldextivlrain(i,j,k,ni) = 1.00E+03*1.505E+00*     &
                                            conc_rain(i,j,k)/     &
                                            size_d(i,j,k)  
                  cldssalbivlrain(i,j,k,ni) = 1.0E+00 - (a(ni)*    &
                                              (rcap(i,j,k)**b(ni)))
                  cldasymmivlrain(i,j,k,ni) = asymm(ni)
                end do
              else
                call error_mesg ('microphys_rad_mod', &
                          'rain drop size out of range', FATAL)
              endif
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------
 


end subroutine savijarvi


!####################################################################
! <SUBROUTINE NAME="fu">
!  <OVERVIEW>
!   Subroutine to define the single scattering parameters for ice crystals
!  </OVERVIEW>
!  <DESCRIPTION>
!   define the single scattering parameters for ice crystals using the    
! Fu parameterization for his spectral intervals.
!                                                                       
! references:                                                           
!                                                                       
! fu, q., an accurate parameterization of the solar radiative           
!      properties of cirrus clouds for climate models., j. climate,     
!      9, 2058-2082, 1996.                                              
!                                                                       
! notes: the ice crystal effective size (D^sub^ge in his paper) can     
!        only be 18.6 <= D^sub^ge <= 130.2 microns.                     
!                                                                       
!        the single scattering properties for wavenumbers < 2000 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division   
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call fu                                                 &
!                    (conc_ice    , size_ice      ,                  &
!                     cldextivlice, cldssalbivlice, cldasymmivlice)
!  </TEMPLATE>
!  <IN NAME="conc_ice" TYPE="real">
!   the ice water concentation in grams / meter**3
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   the ice crystal effective size in microns
!  </IN>
!  <OUT NAME="cldextivlice" TYPE="real">
!   the specified spectral values of the extinction      
!                  coefficient for ice particles in kilometers**(-1)
!  </OUT>
!  <OUT NAME="cldssalbivlice" TYPE="real">
!   the specified spectral values of the single-         
!                  scattering albedo for ice particles
!  </OUT>
!  <OUT NAME="cldasymmivlice" TYPE="real">
!   the specified spectral values of the asymmetry       
!                  factor for ice particles
!  </OUT>
! </SUBROUTINE>
!
subroutine fu (conc_ice, size_ice, dge_column, cldextivlice,  &
               cldssalbivlice, cldasymmivlice, mask, starting_band, &
               ending_band)
 
!----------------------------------------------------------------------
!    subroutine fu defines the single scattering parameters for ice
!    crystals using the Fu parameterization for his spectral intervals.
!    references:                                                
!    fu, q., an accurate parameterization of the solar radiative    
!    properties of cirrus clouds for climate models., j. climate,     
!    9, 2058-2082, 1996.                                              
!---------------------------------------------------------------------- 
                                                                        
real, dimension (:,:,:),   intent(in)    ::   conc_ice, size_ice
logical, dimension (:,:,:),   intent(in)    ::   dge_column          
real, dimension (:,:,:,:), intent(inout)   ::  cldextivlice,      &
                                               cldssalbivlice,    &
                                               cldasymmivlice
logical, dimension(:,:,:), intent(out)    ::   mask
integer,     intent(in), optional          ::  starting_band, &
                                               ending_band

!----------------------------------------------------------------------
!  intent(in) variables:                                              
!                                                                       
!        conc_ice         ice water concentation [ grams / meter**3 ]
!        size_ice         ice crystal effective size [ microns ]  
!
!  intent(out) variables:                                               
!                                                                       
!        cldextivlice     extinction coefficient in each spectral
!                         interval of the fu ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ km **(-1) ]
!        cldssalbivlice   single scattering albedo in each spectral
!                         interval of the fu ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ dimensionless ]
!        cldasymmivlice   asymmetry factor in each spectral
!                         interval of the fu ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ dimensionless ]
!
!   intent(in), optional variables:
!
!        starting_band    the index of the first ice crystal spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!        ending_band      the index of the last ice crystal spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------c
! local variables:                                                     c

      real, dimension (size(conc_ice,1), size(conc_ice,2),  &
                       size(conc_ice,3))   ::  size_i

      real, dimension (NICECLDIVLS) ::  a0fu, a1fu,             &
                                        b0fu, b1fu, b2fu, b3fu,       &
                                        c0fu, c1fu, c2fu, c3fu
 
      data a0fu / -2.54823E-04, 1.87598E-04, 2.97295E-04, 2.34245E-04, &
                   4.89477E-04,-8.37325E-05, 6.44675E-04,-8.05155E-04, &
                   6.51659E-05, 4.13595E-04,-6.14288E-04, 7.31638E-05, &
                   8.10443E-05, 2.26539E-04,-3.04991E-04, 1.61983E-04, &
                   9.82244E-05,-3.03108E-05,-9.45458E-05, 1.29121E-04, &
                  -1.06451E-04,-2.58858E-04,-2.93599E-04,-2.66955E-04, &
                  -2.36447E-04 /
      data a1fu /  2.52909E+00, 2.51396E+00, 2.48895E+00, 2.48573E+00, &
                   2.48776E+00, 2.52504E+00, 2.47060E+00, 2.57600E+00, &
                   2.51660E+00, 2.48783E+00, 2.56520E+00, 2.51051E+00, &
                   2.51619E+00, 2.49909E+00, 2.54412E+00, 2.50746E+00, &
                   2.50875E+00, 2.51805E+00, 2.52061E+00, 2.50410E+00, &
                   2.52684E+00, 2.53815E+00, 2.54540E+00, 2.54179E+00, &
                   2.53817E+00 /
      data b0fu /  2.60155E-01, 1.96793E-01, 4.64416E-01, 9.05631E-02, &
                   5.83469E-04, 2.53234E-03, 2.01931E-03,-2.85518E-05, &
                  -1.48012E-07, 6.47675E-06,-9.38455E-06,-2.32733E-07, &
                  -1.57963E-07,-2.75031E-07, 3.12168E-07,-7.78001E-08, &
                  -8.93276E-08, 9.89368E-08, 5.08447E-07, 7.10418E-07, &
                   3.25057E-08,-1.98529E-07, 1.82299E-07,-1.00570E-07, &
                  -2.69916E-07 /
      data b1fu/   5.45547E-03, 5.75235E-03, 2.04716E-05, 2.93035E-03, &
                   1.18127E-03, 1.75078E-03, 1.83364E-03, 1.71993E-03, &
                   9.02355E-05, 2.18111E-05, 1.77414E-05, 6.41602E-06, &
                   1.72475E-06, 9.72285E-07, 4.93304E-07, 2.53360E-07, &
                   1.14916E-07, 5.44286E-08, 2.73206E-08, 1.42205E-08, &
                   5.43665E-08, 9.39480E-08, 1.12454E-07, 1.60441E-07, &
                   2.12909E-07 /
      data b2fu / -5.58760E-05,-5.29220E-05,-4.60375E-07,-1.89176E-05, &
                  -3.40011E-06,-8.00994E-06,-7.00232E-06,-7.43697E-06, &
                  -1.98190E-08, 1.83054E-09,-1.13004E-09, 1.97733E-10, &
                   9.02156E-11,-2.23685E-10, 1.79019E-10,-1.15489E-10, &
                  -1.62990E-10,-1.00877E-10, 4.96553E-11, 1.99874E-10, &
                  -9.24925E-11,-2.54540E-10,-1.08031E-10,-2.05663E-10, &
                  -2.65397E-10 /
      data b3fu /  1.97086E-07, 1.76618E-07, 2.03198E-09, 5.93361E-08, &
                   8.78549E-09, 2.31309E-08, 1.84287E-08, 2.09647E-08, &
                   4.01914E-11,-8.28710E-12, 2.37196E-12,-6.96836E-13, &
                  -3.79423E-13, 5.75512E-13,-7.31058E-13, 4.65084E-13, &
                   6.53291E-13, 4.56410E-13,-1.86001E-13,-7.81101E-13, &
                   4.53386E-13, 1.10876E-12, 4.99801E-13, 8.88595E-13, &
                   1.12983E-12 /
      data c0fu /  7.99084E-01, 7.59183E-01, 9.19599E-01, 8.29283E-01, &
                   7.75916E-01, 7.58748E-01, 7.51497E-01, 7.52528E-01, &
                   7.51277E-01, 7.52292E-01, 7.52048E-01, 7.51715E-01, &
                   7.52318E-01, 7.51779E-01, 7.53393E-01, 7.49693E-01, &
                   7.52131E-01, 7.51135E-01, 7.49856E-01, 7.48613E-01, &
                   7.47054E-01, 7.43546E-01, 7.40926E-01, 7.37809E-01, &
                   7.33260E-01 /
      data c1fu /  4.81706E-03, 4.93765E-03, 5.03025E-04, 2.06865E-03, &
                   1.74517E-03, 2.02709E-03, 2.05963E-03, 1.95748E-03, &
                   1.29824E-03, 1.14395E-03, 1.12044E-03, 1.10166E-03, &
                   1.04224E-03, 1.03341E-03, 9.61630E-04, 1.05446E-03, &
                   9.37763E-04, 9.09208E-04, 8.89161E-04, 8.90545E-04, &
                   8.86508E-04, 9.08674E-04, 8.90216E-04, 8.97515E-04, &
                   9.18317E-04 /
      data c2fu / -5.13220E-05,-4.84059E-05,-5.74771E-06,-1.59247E-05, &
                  -9.21314E-06,-1.17029E-05,-1.12135E-05,-1.02495E-05, &
                  -4.99075E-06,-3.27944E-06,-3.11826E-06,-2.91300E-06, &
                  -2.26618E-06,-2.13121E-06,-1.32519E-06,-2.32576E-06, &
                  -9.72292E-07,-6.34939E-07,-3.49578E-07,-3.44038E-07, &
                  -2.59305E-07,-4.65326E-07,-1.87919E-07,-2.17099E-07, &
                  -4.22974E-07 /
      data c3fu /  1.84420E-07, 1.65801E-07, 2.01731E-08, 5.01791E-08, &
                   2.15003E-08, 2.95195E-08, 2.73998E-08, 2.35479E-08, &
                   6.33757E-09,-2.42583E-10,-5.70868E-10,-1.37242E-09, &
                  -3.68283E-09,-4.24308E-09,-7.17071E-09,-3.58307E-09, &
                  -8.62063E-09,-9.84390E-09,-1.09913E-08,-1.10117E-08, &
                  -1.13305E-08,-1.05786E-08,-1.16760E-08,-1.16090E-08, &
                  -1.07976E-08 /
 
      integer     :: nistart, niend
      integer     :: i, j, k, ni
      real        :: fd, f, fw

!---------------------------------------------------------------------
!   local variables:
!
!        a0fu     interval-dependent parameter used to define extinction
!                 coefficient due to ice crystals in the fu 
!                 parameterization
!        a1fu     interval-dependent parameter used to define extinction
!                 coefficient due to ice crystals in the fu 
!                 parameterization
!        b0fu     interval-dependent parameter used to define single-
!                 scattering albedo due to ice crystals in the fu 
!                 parameterization
!        b1fu     interval-dependent parameter used to define single-
!                 scattering albedo due to ice crystals in the fu 
!                 parameterization
!        b2fu     interval-dependent parameter used to define single-
!                 scattering albedo due to ice crystals in the fu 
!                 parameterization
!        b3fu     interval-dependent parameter used to define single-
!                 scattering albedo due to ice crystals in the fu 
!                 parameterization
!        c0fu     interval-dependent parameter used to define asymmetry 
!                 factor due to ice crystals in the fu parameterization
!        c1fu     interval-dependent parameter used to define asymmetry
!                 factor due to ice crystals in the fu parameterization
!        c2fu     interval-dependent parameter used to define asymmetry
!                 factor due to ice crystals in the fu parameterization
!        c3fu     interval-dependent parameter used to define asymmetry 
!                 factor due to ice crystals in the fu parameterization
!        nistart  first ice crystal parameterization band to be
!                 processed
!        niend    last ice crystal parameterization band to be processed
!        i,j,k,ni do-loop indices
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the starting and ending droplet parameterization bands to
!    be processed during this call.
!--------------------------------------------------------------------
      if (present(starting_band)) then
        nistart = starting_band
      else
        nistart = 1
      endif
      if (present(ending_band)) then
        niend = ending_band
      else
        niend = NICECLDIVLS
      endif

!----------------------------------------------------------------------
      do k=1,size(conc_ice,3)
        do j=1,size(conc_ice,2)
          do i=1,size(conc_ice,1)

            if (dge_column(i,j,k)) then
!----------------------------------------------------------------------
!    if no ice crystals are present in a grid box, set the scattering
!    parameters to zero.
!----------------------------------------------------------------------
              if (conc_ice (i,j,k) == 0.0) then
                mask(i,j,k) = .false.

!----------------------------------------------------------------------
!    compute the ice crystal scattering parameters.
!----------------------------------------------------------------------
              else !(conc_ice > 0)
                mask(i,j,k) = .true.

!--------------------------------------------------------------------
!    the ice crystal effective size (D^sub^ge in fu's paper) is limited.
!--------------------------------------------------------------------
                if (size_ice(i,j,k) <= min_cld_ice_size) then
                  size_i(i,j,k) = min_cld_ice_size      
                else if (size_ice(i,j,k) > max_cld_ice_size ) then    
                  size_i(i,j,k) = max_cld_ice_size         
                else
                  size_i(i,j,k) = size_ice(i,j,k)          
                endif

!---------------------------------------------------------------------
!     compute the scattering parameters for each of the fu spectral 
!     intervals. the extinction coefficient is converted to km**(-1).  
!---------------------------------------------------------------------
                do ni = nistart, niend
                  cldextivlice(i,j,k,ni) = 1.0E+03*conc_ice(i,j,k)*  &
                                           (a0fu(ni) + (a1fu(ni)/    &
                                            size_i(i,j,k)     ))
                  cldssalbivlice(i,j,k,ni) =  1.0 -                &
                                     ( b0fu(ni)                    +   &
                                       b1fu(ni)*size_i(i,j,k)    +   &
                                       b2fu(ni)*size_i(i,j,k)**2 +   &
                                       b3fu(ni)*size_i(i,j,k)**3 )
                  cldasymmivlice(i,j,k,ni) =                        &
                          c0fu(ni) +                                   &
                          c1fu(ni)*size_i(i,j,k) +                   &
                          c2fu(ni)*size_i(i,j,k)**2 +                &
                          c3fu(ni)*size_i(i,j,k)**3

                  if (do_delta_adj .and. (.not. do_const_asy)) then
                    fd =                                           &
                          1.1572963e-1 +                       &
                          2.5648064e-4*size_ice(i,j,k) +     &
                          1.9131293e-6*size_ice(i,j,k)**2    &
                         -1.2460341e-8*size_ice(i,j,k)**3
                    f = 0.5/cldssalbivlice(i,j,k,ni) + fd
                    fw = f * cldssalbivlice(i,j,k,ni)
                    cldextivlice(i,j,k,ni) =   &
                            cldextivlice(i,j,k,ni) * (1. - fw)
                    cldssalbivlice(i,j,k,ni) =   &
                          cldssalbivlice(i,j,k,ni) * (1. - f)/(1. - fw)
                    cldasymmivlice(i,j,k,ni) =  &
                          (cldasymmivlice(i,j,k,ni) - f)/(1. - f)
                  endif
  
                  if (do_const_asy .and. (.not. do_delta_adj)) then
                    f = 0.5/cldssalbivlice(i,j,k,ni)
                    fw = f * cldssalbivlice(i,j,k,ni)
                    cldextivlice(i,j,k,ni) = cldextivlice(i,j,k,ni) &
                                                        * (1. - fw)
                    cldssalbivlice(i,j,k,ni) =  &
                           cldssalbivlice(i,j,k,ni) * (1. - f)/(1. - fw)
                    cldasymmivlice(i,j,k,ni) = (val_const_asy - f)/ &
                                                            (1. - f)
                  endif

                  if (do_delta_adj .and. do_const_asy) then
                    fd =                                           &
                        1.1572963e-1 +                         &
                        2.5648064e-4*size_ice(i,j,k) +              &
                        1.9131293e-6*size_ice(i,j,k)**2       &
                       -1.2460341e-8*size_ice(i,j,k)**3
                    f = 0.5/cldssalbivlice(i,j,k,ni) + fd
                    fw = f * cldssalbivlice(i,j,k,ni)
                    cldextivlice(i,j,k,ni) = cldextivlice(i,j,k,ni) &
                                                         * (1. - fw)
                    cldssalbivlice(i,j,k,ni) =  &
                          cldssalbivlice(i,j,k,ni) * (1. - f)/(1. - fw)
                    cldasymmivlice(i,j,k,ni) =  &
                                         (val_const_asy - f)/(1. - f)
                  endif
                end do
              endif ! (conc_ice > 0)
            else ! (dge_column)
              mask(i,j,k) = .false.
            endif ! (dge_column)
          end do
        end do
      end do

!---------------------------------------------------------------------
 


end subroutine fu



!#####################################################################
! <SUBROUTINE NAME="icesolar">
!  <OVERVIEW>
!   Subroutine to define the single scattering parameters for ice crystals
!  </OVERVIEW>
!  <DESCRIPTION>
!   define the single scattering parameters for ice crystals using the    
! Fu parameterization for his spectral intervals.                       
!                                                                       
! references:                                                           
!                                                                       
! fu and Liou (1993, JAS) 
!                                                                       
! notes: the ice crystal effective size (D^sub^e in paper) can          
!        only be 18.6 <= D^sub^e <= 130.2 microns.                     
!                                                                       
!        the single scattering properties for wavenumbers < 2000 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division  
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call icesolar                                           &
!                    (conc_ice    , size_ice      ,                  &
!                     cldextivlice, cldssalbivlice, cldasymmivlice)
!  </TEMPLATE>
!  <IN NAME="conc_ice" TYPE="real">
!   the ice water concentation in grams / meter**3
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   the ice crystal effective size in microns                  
! Corresponds to minimum dimension of hexagonal crystal.
!  </IN>
!  <OUT NAME="cldextivlice" TYPE="real">
!   the specified spectral values of the extinction      
!                  coefficient for ice particles in kilometers**(-1)
!  </OUT>
!  <OUT NAME="cldssalbivlice" TYPE="real">
!   the specified spectral values of the single-         
!                  scattering albedo for ice particles
!  </OUT>
!  <OUT NAME="cldasymmivlice" TYPE="real">
!   the specified spectral values of the asymmetry       
!                  factor for ice particles
!  </OUT>
! </SUBROUTINE>
!
subroutine icesolar (conc_ice, size_ice, dge_column, cldextivlice,    &
                     cldssalbivlice, cldasymmivlice, &
                     mask, &
                     starting_band, ending_band)
 
!---------------------------------------------------------------------- 
!    subroutine icesolar defines the single scattering parameters for 
!    ice crystals using the fu and liou (1993) parameterization for 
!    their spectral intervals. references:                      
!    Fu and Liou (1993, JAS) 
!----------------------------------------------------------------------
                                                                   
real, dimension (:,:,:),   intent(in)   ::   conc_ice, size_ice
logical, dimension (:,:,:),   intent(in)   ::   dge_column        
real, dimension (:,:,:,:), intent(inout)  ::   cldextivlice,           &
                                               cldssalbivlice,         &
                                               cldasymmivlice
logical, dimension(:,:,:), intent(out)    ::   mask
integer,  intent(in), optional            ::   starting_band, &
                                               ending_band

!---------------------------------------------------------------------- 
! intent(in) variables:                                                 
!                                                                       
!        conc_ice         ice water concentation [ grams / meter**3 ]  
!        size_ice         ice crystal effective size. this corresponds 
!                         to the minimum dimension of hexagonal crystal.
!                         [ microns ]
!                                                                       
! intent(out) variables:                                                
!                                                                       
!        cldextivlice     extinction coefficient in each spectral
!                         interval of the fu and liou ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ km **(-1) ]
!        cldssalbivlice   single scattering albedo in each spectral
!                         interval of the fu and liou ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ dimensionless ]
!        cldasymmivlice   asymmetry factor in each spectral interval of
!                         the fu iand liou ice crystal param-
!                         eterization resulting from the presence of 
!                         ice crystals [ dimensionless ]
!
!   intent(in), optional variables:
!
!        starting_band    the index of the first ice crystal spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!        ending_band      the index of the last ice crystal spectral
!                         band contained in the sw parameterization
!                         band(s) being processed
!
!---------------------------------------------------------------------
 
!---------------------------------------------------------------------- 
! local variables:                                                      

      real, dimension (size(conc_ice,1), size(conc_ice,2),  &
                       size(conc_ice,3))   ::  size_i
      real, dimension (1:NICESOLARCLDIVLS, 0:NBB) :: b
      real, dimension (1:NICESOLARCLDIVLS, 0:NBC) :: c
      real, dimension (1:NICESOLARCLDIVLS, 0:NBD) :: d

      data b     /                                           &
 .10998e-5,  .20208e-4, .1359e-3,  -.16598e-2,  .4618,      .42362e-1, &
-.26101e-7,  .96483e-5, .73453e-3,  .20933e-2,  .24471e-3,  .86425e-2, &
 .10896e-8,  .83009e-7, .28281e-5, -.13977e-5, -.27839e-5, -.75519e-4, &
-.47387e-11,-.32217e-9,-.18272e-7, -.18703e-7,  .10379e-7,  .24056e-6/

      data c     /                                            &
 2.211,      2.2151,    2.2376,    2.3012,    2.7975,    1.9655,      &
 -.10398e-2, -.77982e-3, .10293e-2, .33854e-2, .29741e-2, .20094e-1,  &
  .65199e-4,  .6375e-4,  .50842e-4, .23528e-4,-.32344e-4,-.17067e-3,  &
 -.34498e-6, -.34466e-6,-.30135e-6,-.20068e-6, .11636e-6, .50806e-6 /

      data d     /                                           & 
  .12495,    .12363,    .12117,    .11581,   -.15968e-3, .1383,       &
 -.43582e-3,-.44419e-3,-.48474e-3,-.55031e-3, .10115e-4,-.18921e-2,   &
  .14092e-4, .14038e-4, .12495e-4, .98776e-5,-.12472e-6, .1203e-4,    &
 -.69565e-7,-.68851e-7,-.62411e-7,-.50193e-7, .48667e-9,-.31698e-7 /

      real    :: a0 = -6.656e-03
      real    :: a1 =  3.686
      real    :: fgam2, fdel2

      integer :: nistart, niend
      integer :: i, j, k, ni

!----------------------------------------------------------------------
!   local variables:
!
!        b             coefficients in fu and liou expression for
!                      single scattering albedo; when second dimension
!                      has values of 0 -> 3, units are [ dimensionless,
!                      microns**(-1), microns**(-2) and microns**(-3) ]
!                      respectively.
!        c             coefficients in fu and liou expression for
!                      asymmetry factor; when second dimension
!                      has values of 0 -> 3, units are [ dimensionless,
!                      microns**(-1), microns**(-2) and microns**(-3) ]
!                      respectively.
!        d             coefficients in fu and liou expression for
!                      asymmetry factor; when second dimension
!                      has values of 0 -> 3, units are [ dimensionless, 
!                      microns**(-1), microns**(-2) and microns**(-3) ]
!                      respectively.
!        a0             parameter in fu and liou extinction formula 
!                       [ m**2 / g ]
!        a1             parameter in fu and liou extinction formula
!                       [ (m**2)*microns / g ]
!        fgam2          intermediate expression used in evaluating
!                       asymmetry factor
!        fdel2          intermediate expression used in evaluating
!                       asymmetry factor
!        nistart        first ice crystal parameterization band to be
!                       processed
!        niend          last ice crystal parameterization band to be
!                       processed
!        i,j,k,ni       do-loop indices
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the starting and ending droplet parameterization bands to
!    be processed during this call.
!--------------------------------------------------------------------
      if (present(starting_band)) then
        nistart = starting_band
      else
        nistart = 1
      endif
      if (present(ending_band)) then
        niend = ending_band
      else
        niend = NICESOLARCLDIVLS
      endif

!-----------------------------------------------------------------
!    compute scattering parameters for ice crystals. 
!-----------------------------------------------------------------
      do k=1,size(conc_ice,3)
        do j=1,size(conc_ice,2)
          do i=1,size(conc_ice,1)

            if (.not. dge_column(i,j,k)) then
!---------------------------------------------------------------------
!    bypass calculations if no crystals are present. set scattering
!    parameters to values comatible with the absence of cloud.
!-----------------------------------------------------------------
            if (conc_ice (i,j,k) == 0.0) then
              mask(i,j,k) = .false.

!--------------------------------------------------------------------
!    the ice crystal effective size (D^sub^ge in fu's paper) is limited
!    to the range of 18.6 to 130.2 microns.                     
!--------------------------------------------------------------------
            else
              mask(i,j,k) = .true.

              if (size_ice(i,j,k) < min_cld_ice_size) then           
                size_i(i,j,k) = min_cld_ice_size             
              else if(size_ice(i,j,k) > max_cld_ice_size ) then  
                size_i(i,j,k) = max_cld_ice_size       
              else
                size_i(i,j,k) = size_ice(i,j,k)
              endif

!---------------------------------------------------------------------
!     compute the scattering parameters for each of the fu spectral 
!     intervals. the extinction coefficient is converted to km**(-1).  
!---------------------------------------------------------------------
                do ni = nistart,niend
                  cldextivlice(i,j,k,ni) = 1.0E+03*       &
                         conc_ice(i,j,k)*(a0 + (a1/size_i(i,j,k))) 
                  cldssalbivlice(i,j,k,ni) = 1.0 -           &
                             (b(7-ni,0) +                    &
                              b(7-ni,1)*size_i(i,j,k) +    &
                              b(7-ni,2)*size_i(i,j,k)**2 + &
                              b(7-ni,3)*size_i(i,j,k)**3 )
                  fgam2  =                                   &
                              c(7-ni,0) +                    &
                              c(7-ni,1)*size_i(i,j,k) +    &
                              c(7-ni,2)*size_i(i,j,k)**2 + &
                              c(7-ni,3)*size_i(i,j,k)**3
                  fdel2  =                                   &
                              d(7-ni,0) +                    &
                              d(7-ni,1)*size_i(i,j,k) +    &
                              d(7-ni,2)*size_i(i,j,k)**2 + &
                              d(7-ni,3)*size_i(i,j,k)**3
                  cldasymmivlice(i,j,k,ni) =                 &
                              ((1. - fdel2)*fgam2 + 3.*fdel2)/3.
                end do
            endif
           else
             mask(i,j,k) = .false.
           endif
          end do
        end do
      end do

!------------------------------------------------------------------


end subroutine icesolar



!######################################################################
! <SUBROUTINE NAME="snowsw">
!  <OVERVIEW>
!   Subroutine to define the single scattering parameters for snow
!  </OVERVIEW>
!  <DESCRIPTION>
!   define the single scattering parameters for snow using the Fu         
!   parameterization for his spectral intervals.
! author: leo donner, gfdl, 11 Sept 98                                  
!                                                                       
! references:                                                           
!                                                                       
! fu, q., et al., (See notes from Kuo-Nan Liou, 1 Sept 98). (SNOW)      
!                                                                       
! notes: the single scattering properties for wavenumbers < 2500 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is in units of kilometer**(-1)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call snowsw                                            &
!                    (conc_snow,                                    &
!                     cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow)
!  </TEMPLATE>
!  <IN NAME="conc_snow" TYPE="real">
!   the snow concentration in grams / meter**3
!  </IN>
!  <OUT NAME="cldextivlsnow" TYPE="real">
!   the specified spectral values of the extinction     
!                   coefficient for snow in kilometers**(-1)
!  </OUT>
!  <OUT NAME="cldssalbivlsnow" TYPE="real">
!   the specified spectral values of the single-        
!                   scattering albedo for snow
!  </OUT>
!  <OUT NAME="cldasymmivlsnow" TYPE="real">
!   the specified spectral values of the asymmetry      
!                   factor for snow
!  </OUT>
! </SUBROUTINE>
!
subroutine snowsw (conc_snow, cldextivlsnow, cldssalbivlsnow,    &
                   cldasymmivlsnow, mask, starting_band, ending_band)      
 
!----------------------------------------------------------------------
!    subroutine snowsw defines the single scattering parameters for snow
!    flakes  using the Fu parameterization for his spectral intervals. 
!    author: leo donner, gfdl, 11 Sept 98                     
!    references:                                     
!    fu, q., et al., (See notes from Kuo-Nan Liou, 1 Sept 98). (SNOW)  
!----------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)     ::  conc_snow
real, dimension (:,:,:,:), intent(out)    ::  cldextivlsnow,       &
                                              cldssalbivlsnow,     &
                                              cldasymmivlsnow
logical, dimension(:,:,:), intent(out)    ::   mask
integer,  intent(in), optional            ::   starting_band, &
                                               ending_band

!----------------------------------------------------------------------
!   intent(in) variables:                                
!                                                                       
!        conc_snow        snow concentration [ grams / meter**3 ]     
!
!  intent(out) variables:                                               
!                                                                       
!        cldextivlsnow    extinction coefficient in each spectral
!                         interval of the fu snow flake param-
!                         eterization resulting from the presence of 
!                         snow flakes [ km **(-1) ]
!        cldssalbivlsnow  single scattering albedo in each spectral
!                         interval of the fu snow flake param-
!                         eterization resulting from the presence of 
!                         snow flakes [ dimensionless ]
!        cldasymmivlsnow  asymmetry factor in each spectral
!                         interval of the fu snow flake param-
!                         eterization resulting from the presence of 
!                         snow flakes [ dimensionless ]
!
!---------------------------------------------------------------------
 
!---------------------------------------------------------------------
! local variables:                                                     
      real, dimension (NSNOWCLDIVLS)  ::  asymm, ext, ssalb
 
      data asymm / 9.6373E-01,9.8141E-01,9.7816E-01,9.6820E-01,      &
                   8.9940E-01,8.9218E-01 /
      data ext   / 8.3951E-01,8.3946E-01,8.3941E-01,8.3940E-01,      &
                   8.3940E-01,8.3939E-01 /
      data ssalb / 5.3846E-01,5.2579E-01,5.3156E-01,5.6192E-01,      &
                   9.7115E-01,9.99911E-01 /

      real          ::  conc_ref=0.5
      integer       ::  i, j, k, ni
      integer      :: nistart, niend

!---------------------------------------------------------------------
! local variables:                                                     
!
!       asymm        asymmetry factor due to snow flakes in each of the 
!                    snow spectral bands [ dimensionless ]
!       ext          extinction coefficient due to snow flakes in each
!                    of the snow spectral bands, relevant for a snow
!                    concentration of conc_ref [ km**(-1) ]
!       ssalb        single scattering albedo due to snow flakes in each
!                    of the snow flake spectral bands [ dimensionless ]
!       conc_ref     reference snow flake concentration for which the
!                    values of ext apply [ g / m**3 ]
!       i,j,k,ni     do-loop indices
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the starting and ending droplet parameterization bands to
!    be processed during this call.
!--------------------------------------------------------------------
      if (present(starting_band)) then
        nistart = starting_band
      else
        nistart = 1
      endif
      if (present(ending_band)) then
        niend = ending_band
      else
        niend = NSNOWCLDIVLS
      endif

!----------------------------------------------------------------------
      do k=1,size(conc_snow,3)
        do j=1,size(conc_snow,2)
          do i=1,size(conc_snow,1)

!-----------------------------------------------------------------
!    if no snow is present in the box, set the scattering parameters
!    to so indicate.
!-----------------------------------------------------------------
            if (conc_snow(i,j,k) .le. 1.e-5) then
              mask(i,j,k) = .false.

!---------------------------------------------------------------------
!    if snow is present, calculate the scattering parameters over each
!    of the snow spectral intervals. the extinction coefficient is 
!    in units of km**(-1).
!---------------------------------------------------------------------
            else
              mask(i,j,k) = .true.
              do ni = nistart, niend
                cldextivlsnow(i,j,k,ni) = ext(ni)*conc_snow(i,j,k)/   &
                                          conc_ref
                cldssalbivlsnow(i,j,k,ni) = ssalb(ni)
                cldasymmivlsnow(i,j,k,ni) = asymm(ni)
              end do
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
 

end subroutine snowsw


!######################################################################
! <SUBROUTINE NAME="cloud_lwpar">
!  <OVERVIEW>
!   Subroutine to determine cloud infrared emissivity
!  </OVERVIEW>
!  <DESCRIPTION>
! determine the infrared cloud emissivities for specified wavenumber    
! bands from parameterizations for absorption coefficients due to       
! cloud drops, cloud ice crystals, rain and snow. conceptually one      
! could have separate concentrations and sizes for "thin" or randomly   
! overlapped and for maximally overlapped clouds. for now, there is     
! one concentration and size, therefore the two emissivities are set    
! equal.
!  </DESCRIPTION>
!  <TEMPLATE>
!   subroutine cloud_lwpar  (nonly, nbmax, nnn,                  &
!                     size_drop, size_ice, size_rain,            &
!                     conc_drop, conc_ice, conc_rain, conc_snow, &
!                     do_dge_lw, abscoeff)
!  </TEMPLATE>
!  <IN NAME="nonly" TYPE="integer">
!   The single band for calculations.  Note that this is used
!   only in the case of a call from cloudrad_diagnostics to 
!   do isccp simulator work.  For all other calls, nonly should
!   be 0 and will have no effect on the calculations below
!  </IN>
!  <IN NAME="nbmax" TYPE="integer">
!   The number of individual bands to do calculations over. Note
!   that for normal GCM calls this will be 1.  For calls using
!   stochastic clouds with or without the isccp simulator this will
!   be equal to the number of longwave bands
!  </IN>
!  <IN NAME="nnn" TYPE="integer">
!   This integer controls which cloud state to access for radiation
!   calculations.  For normal GCM applications this will be 1. For
!   Full Independent Column Approximation calculations with stochast-
!   ic clouds this will be the profile number being accessed. 
!  </IN>
!  <IN NAME="conc_drop" TYPE="real">
!   total cloud droplet concentration
!  </IN>
!  <IN NAME="conc_ice" TYPE="real">
!   ice cloud droplet concentration
!  </IN>
!  <IN NAME="conc_rain" TYPE="real">
!   rain droplet concetration
!  </IN>
!  <IN NAME="conc_snow" TYPE="real">
!   snow concentration
!  </IN>
!  <IN NAME="size_drop" TYPE="real">
!   cloud droplet size distribution
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   ice droplet size distribution
!  </IN>
!  <IN NAME="size_rain" TYPE="real">
!   rain droplet size distribution
!  </IN>
!  <OUT NAME="abscoeff" TYPE="real">
!   cloud absorption coefficient
!  </OUT>
!  <IN NAME="do_dge_lw" TYPE="logical">
!   flag for using dge longwave parameterization
!  </IN>
! </SUBROUTINE>
subroutine cloud_lwpar (nonly, nbmax, nnn, size_drop, size_ice,  &
                        size_rain, conc_drop, conc_ice, conc_rain,  &
                        conc_snow, do_dge_lw, dge_column, isccp_call, &
                        abscoeff)
 
!----------------------------------------------------------------------
!    cloud_lwpar determines the absorption coefficients due to cloud 
!    drops, cloud ice crystals, rain and snow over appropriate spectral
!    intervals and then maps these into the appropriate lw radiation
!    bands. the contributions from each of the water species are then
!    consolidated into a single absorption coefficient which is output
!    to the calling routine.
!----------------------------------------------------------------------

integer,                   intent(in)   ::   nonly, nbmax, nnn
real, dimension (:,:,:,:), intent(in)   ::   size_drop, size_ice,    &
                                             conc_drop, conc_ice
logical, dimension (:,:,:,:), intent(in)   :: dge_column
real, dimension (:,:,:),   intent(in)   ::   size_rain, conc_rain,   &
                                             conc_snow
logical,                   intent(in)   ::   do_dge_lw
logical,                   intent(in)   ::   isccp_call
real, dimension (:,:,:,:), intent(out)  ::   abscoeff 
 
!
!----------------------------------------------------------------------
! intent(in) variables:                                             
!                                                                   
!       size_drop    cloud drop effective diameter [ microns ]          
!       size_ice     ice crystal effective size [ microns ]    
!       size_rain    rain drop effective diameter [ microns ] 
!       conc_drop    cloud drop liquid water concentration 
!                    [ grams / meter**3 ]                               
!       conc_ice     ice water concentation [ grams / meter**3 ] 
!       conc_rain    rain drop water concentration [ grams / meter**3 ] 
!       conc_snow    snow concentration [ grams / meter**3 ]  
!       do_dge_lw    if true, use parameterization using generalized 
!                    effective size developed by Fu et al (1998); other-
!                    wise use parameterization by Fu et al using 
!                    effective size.
!                                    
! intent(out) variable:                                             
!                                                                   
!       abscoeff     infrared absorption coefficient. [ km**(-1) ]
!                                                                   
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
! local variables:                                                  

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3), N_EMISS_BDS)  ::    &
             cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw, &
             cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw, &
             cldextbndicelw,  cldssalbbndicelw,  cldasymmbndicelw,    &
             cldextbnddroplw
      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3)             )  ::    &
             cldext, cldssa, cldasy, cldext2, cldssa2, cldasy2
      logical, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3)             )  ::    &
                 maskf, maski
      

      integer  :: i,j,k,n

!----------------------------------------------------------------------
! local variables:                                                  
!
!       cldextbndrainlw    values of the extinction coefficient for 
!                          rain water over the wavenumber bands used by 
!                          the radiation code [ km**(-1) ]
!       cldssalbbndrainlw  values of the single-scattering albedo for 
!                          rain water over the wavenumber bands used by
!                          the radiation code  [ dimensionless ] 
!       cldasymmbndrainlw  values of the asymmetry factor for rain water
!                          over the wavenumber bands used by the 
!                          radiation code  
!       cldextbndsnowlw    values of the extinction coefficient for 
!                          snow flakes over the wavenumber bands used by
!                          the radiation code [ km**(-1) ]
!       cldssalbbndsnowlw  values of the single-scattering albedo for 
!                          snow flakes over the wavenumber bands used by
!                          the radiation code  [ dimensionless ] 
!       cldasymmbndsnowlw  values of the asymmetry factor for snow
!                          flakes over the wavenumber bands used by the 
!                          radiation code  
!       cldextbndicelw     values of the extinction coefficient for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code [ km**(-1) ]
!       cldssalbbndicelw   values of the single-scattering albedo for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code  [ dimensionless ] 
!       cldasymmbndicelw   values of the asymmetry factor for ice 
!                          crystals over the wavenumber bands used by 
!                          the radiation code  
!       cldextbnddroplw    values of the extinction coefficient for 
!                          cloud droplets over the wavenumber bands used
!                          by the radiation code [ km**(-1) ]
!       n                  do-loop index
!                                                                   
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    call furainlw to compute the extinction coefficient, single 
!    scattering coefficient and asymmetry parameter for rain.
!-------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        if (nonly == 0 .or. n == nonly) then
          call furainlw (n, conc_rain, cldextbndrainlw(:,:,:,n),   &
                         cldssalbbndrainlw(:,:,:,n), &
                         cldasymmbndrainlw(:,:,:,n))
        endif
      end do

!----------------------------------------------------------------------
!    call fusnowlw to compute the extinction coefficient, single 
!    scattering coefficient and asymmetry parameter for snow.
!----------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        if (nonly == 0 .or. n == nonly) then
          call fusnowlw (n, conc_snow, cldextbndsnowlw(:,:,:,n),   &
                         cldssalbbndsnowlw(:,:,:,n), &
                         cldasymmbndsnowlw(:,:,:,n))
        endif
      end do

!----------------------------------------------------------------------
!
!  NOTE THE FOLLOWING LOGICAL TO THE LOOPS BELOW
!
!        do n=1,Cldrad_control%nlwcldb
!             if (nbmax==1) then
!                  call cliqlw(conc_drop(:,:,:,nnn).....)
!             else
!                  call slingo(conc_drop(:,:,:,n),....)
!             end if
!
!        enddo           
!             
!        Note that nbmax = 1 in the following cases:
!
!                 (a) standard GCM applications which do not use
!                     McICA  
!                 (b) Full Independent Column Approximation 
!                     calculations
!
!        Note that nbmax = Cldrad_control%nlwcldb in the following cases
!               
!                 (c) McICA calculations where n = the cloud 
!                     profile being used
!                 (d) ISCCP simulator calls from cloudrad_diagnostics
!                     where "nonly" will be used
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    call cliqlw to compute the extinction coefficient for cloud drops.
!----------------------------------------------------------------------
      do n=1,Cldrad_control%nlwcldb
        if (nbmax == 1) then
          call cliqlw (conc_drop(:,:,:,nnn), cldextbnddroplw(:,:,:,n))
        else
          if (nonly.eq.0   .or. nonly.eq.n ) then            
            call cliqlw (conc_drop(:,:,:,n), cldextbnddroplw(:,:,:,n))
          endif 
        endif 
      end do
 
!----------------------------------------------------------------------
!    compute the extinction coefficient, single scattering coefficient
!    and asymmetry parameter for cloud ice crystals. if the generalized
!    effectiuve radius parameterization is to be used call subroutine 
!    el_dge; if it is not, call subroutine el.
!----------------------------------------------------------------------
      if (isccp_call) then
        do n=1,Cldrad_control%nlwcldb
          if (nbmax == 1) then
            call el_dge (n, conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),  &
                          
                  dge_column(:,:,:,nnn), maskf, cldext, cldssa, cldasy)
           else
             if (nonly.eq.0   .or. nonly.eq.n ) then            
             call el_dge (n, conc_ice(:,:,:,n), size_ice(:,:,:,n),  &
                   dge_column(:,:,:,n), maskf, cldext, cldssa, cldasy)
             end if ! for nonly
           endif ! for nbmax == 1
          if (nbmax == 1) then
            call el (n, conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),  &
               dge_column(:,:,:,nnn), maski, cldext2, cldssa2, cldasy2)
          else
            if (nonly.eq.0   .or. nonly.eq.n ) then            
            call el (n, conc_ice(:,:,:,n), size_ice(:,:,:,n),  &
                 dge_column(:,:,:,n), maski, cldext2, cldssa2, cldasy2)
            endif ! for nonly 
          endif ! for nbmax == 1
 
              do k=1, size(conc_drop,3)
                do j=1, size(conc_drop    ,2)
                  do i=1, size(conc_drop,1)
                    if (maskf(i,j,k)) then
                     cldextbndicelw(i,j,k,n) = cldext(i,j,k)
                     cldssalbbndicelw(i,j,k,n) = cldssa(i,j,k)
!                    cldasymmbndicelw(i,j,k,n) = cldasy(i,j,k)
                    else if (maski(i,j,k)) then
                     cldextbndicelw(i,j,k,n) = cldext2(i,j,k)
                     cldssalbbndicelw(i,j,k,n) = cldssa2(i,j,k)
!                    cldasymmbndicelw(i,j,k,n) = cldasy2(i,j,k)
                    else
                     cldextbndicelw(i,j,k,n) = 0.0             
                     cldssalbbndicelw(i,j,k,n) = 0.0              
!                    cldasymmbndicelw(i,j,k,n) = cldasy2(i,j,k)
                    endif
                  end do
                  end do
                  end do
        end do

   else    ! (isccp_call)
      if (do_dge_lw) then
          maski = .false.
        do n=1,Cldrad_control%nlwcldb
          maski = .false.
          if (nbmax == 1) then
            call el_dge (n, conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),  &
                          
                         dge_column(:,:,:,nnn), maskf, &
                         cldextbndicelw(:,:,:,n),     &
                         cldssalbbndicelw(:,:,:,n),   &
                         cldasymmbndicelw(:,:,:,n))
           else
             if (nonly.eq.0   .or. nonly.eq.n ) then            
             call el_dge (n, conc_ice(:,:,:,n), size_ice(:,:,:,n),  &
                         dge_column(:,:,:,n), maskf, &
                          cldextbndicelw(:,:,:,n),     &
                          cldssalbbndicelw(:,:,:,n),    &
                          cldasymmbndicelw(:,:,:,n))
             end if ! for nonly
           endif ! for nbmax == 1
         end do
      else
          maskf = .false.
        do n=1,Cldrad_control%nlwcldb
          if (nbmax == 1) then
            call el (n, conc_ice(:,:,:,nnn), size_ice(:,:,:,nnn),  &
                         dge_column(:,:,:,nnn), maski,  &
                     cldextbndicelw(:,:,:,n),     &
                     cldssalbbndicelw(:,:,:,n),    &
                     cldasymmbndicelw(:,:,:,n))
          else
            if (nonly.eq.0   .or. nonly.eq.n ) then            
            call el (n, conc_ice(:,:,:,n), size_ice(:,:,:,n),  &
                         dge_column(:,:,:,n), maski,  &
                     cldextbndicelw(:,:,:,n),     &
                     cldssalbbndicelw(:,:,:,n),  &
                     cldasymmbndicelw(:,:,:,n))
            endif ! for nonly 
          endif ! for nbmax == 1
        end do
      endif
 
    endif   ! (isccp_call)

!----------------------------------------------------------------------
!    compute absorption coefficient for each species as the product of 
!    the extinction coefficient and (1 - single scattering albedo). 
!    the total absorption coefficient is the sum of the species 
!    absorption coefficients. 
!----------------------------------------------------------------------
      if (nonly == 0) then
        do n=1,Cldrad_control%nlwcldb
          abscoeff(:,:,:,n) = cldextbndicelw(:,:,:,n)*       &
                              (1.0E+00 - cldssalbbndicelw(:,:,:,n)) +  &
                              cldextbnddroplw(:,:,:,n)              +  &
                              cldextbndsnowlw(:,:,:,n)*               &
                              (1.0E+00 - cldssalbbndsnowlw(:,:,:,n)) + &
                              cldextbndrainlw(:,:,:,n)*                &
                              (1.0E+00 - cldssalbbndrainlw(:,:,:,n))
        end do
      else 
        abscoeff(:,:,:,nonly) = cldextbndicelw(:,:,:,nonly)*       &
                           (1.0E+00 - cldssalbbndicelw(:,:,:,nonly)) + &
                           cldextbnddroplw(:,:,:,nonly)          +    &
                           cldextbndsnowlw(:,:,:,nonly)*           &
                           (1.0E+00 - cldssalbbndsnowlw(:,:,:,nonly)) +&
                           cldextbndrainlw(:,:,:,nonly)*           &
                           (1.0E+00 - cldssalbbndrainlw(:,:,:,nonly))
      endif
 
!---------------------------------------------------------------------


end subroutine cloud_lwpar



!######################################################################
! <SUBROUTINE NAME="cloud_lwem_oneband">
!  <OVERVIEW>
!   Subroutine to determine a single broadband cloud infrared emissivity
!  </OVERVIEW>
!  <DESCRIPTION>
! determine the infrared cloud emissivities for a single broadband
! from parameterizations for absorption coefficients due to
! cloud drops and cloud ice crystals. the parameterization comes from
! from the 5-band formulation given by Ebert and Curry (1992,
! J. Geophys. Res., vol. 97, pp. 3831-3836). S. Klein  derived the
! coefficients for the 1-band version used here.
!  </DESCRIPTION>
!  <TEMPLATE>
!   subroutine cloud_lwem_oneband                              &
!                    (conc_drop, conc_ice, size_drop, size_ice, &
!                     abscoeff)
!  </TEMPLATE>
!  <IN NAME="conc_drop" TYPE="real">
!   total cloud droplet concentration
!  </IN>
!  <IN NAME="conc_ice" TYPE="real">
!   ice cloud droplet concentration
!  </IN>
!  <IN NAME="size_drop" TYPE="real">
!   cloud droplet size distribution
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   ice droplet size distribution
!  </IN>
!  <OUT NAME="abscoeff" TYPE="real">
!   cloud absorption coefficient
!  </OUT>
! </SUBROUTINE>
 subroutine cloud_lwem_oneband (conc_drop, conc_ice, size_drop,    &
                                size_ice, abscoeff)

!----------------------------------------------------------------------
!    subroutine cloud_lwem_oneband determines the infrared cloud emis-
!    sivities for a single broadband from parameterizations for absor-
!    ption coefficients due to cloud drops and cloud ice crystals. the 
!    parameterization comes from the 5-band formulation given by Ebert 
!    and Curry (1992, J. Geophys. Res., vol. 97, pp. 3831-3836). 
!    S. Klein  derived the coefficients for the 1-band version used 
!    here.
!----------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)     ::   conc_drop, conc_ice, &
                                               size_drop, size_ice
!real, dimension (:,:,:,:), intent(out)    ::   abscoeff         
real, dimension (:,:,:), intent(out)    ::   abscoeff         

!----------------------------------------------------------------------
!
!   intent(in) variables:
!
!      conc_drop      cloud drop liquid water concentration 
!                     [ grams / m**3 ]                           
!      conc_ice       ice water concentation 
!                     [ grams / m**3 ]                           
!      size_drop      cloud drop effective diameter [ microns ]
!      size_ice       ice crystal effective size [ microns ]
!
!   intent(out) variables:
!
!      abscoeff       one-band infrared absorption coefficient.
!                     [ kilometer**(-1) ]
!
!----------------------------------------------------------------------
 

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(conc_drop,1), size(conc_drop,2), &
                       size(conc_drop,3))                    ::   &
                                     reff_ice, k_liq, k_ice, size_i
      integer  :: i, j, k
 
!---------------------------------------------------------------------
!    local variables:
!
!         reff_ice       effective diameter of ice crystals. [ microns ]
!         k_liq          liquid cloud mass absorption coefficient for 
!                        longwave portion of spectrum 
!                        [ m**2 / kg condensate ]
!         k_ice          ice cloud mass absorption coefficient for 
!                        longwave portion of spectrum 
!                        [ m**2 / kg condensate ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    reff_ice is the effective diameter obtained using an expression 
!    provided by S. Klein. 
!---------------------------------------------------------------------
      
      do k=1,size(conc_ice,3)
        do j=1,size(conc_ice,2)
          do i=1,size(conc_ice,1)
            if (size_ice(i,j,k) < min_cld_ice_size) then         
              size_i(i,j,k) = min_cld_ice_size             
            else if(size_ice(i,j,k) > max_cld_ice_size) then    
              size_i(i,j,k) = max_cld_ice_size       
            else
              size_i(i,j,k) = size_ice(i,j,k)
            endif
          end do
         end do
      end do
      reff_ice = ((0.1033741*size_i*size_i +      &
                   0.2115169*(size_i**2.272))**0.5)
  
!---------------------------------------------------------------------
!    define the mass absorption coefficients for liquid and ice
!    clouds for the longwave spectrum.
!---------------------------------------------------------------------
      k_liq(:,:,:) = 140.
      where (size_i(:,:,:) /= 0.0) 
        k_ice(:,:,:) = 4.83591 + 1758.511/reff_ice(:,:,:)
      elsewhere
        k_ice(:,:,:) = 0.0                                  
      end where

!--------------------------------------------------------------------- 
!    compute the absorption coefficient. the division by the diffusivity
!    coefficient (diffac) is due to the assumption that the effect of 
!    angular integration is accounted for on the values of k_liq and 
!    k_ice. since the dimensions of k_liq and k_ice are [m**2/Kg] 
!    and conc_ice and conc_drop is in [g/m**3] the unit conversion
!    factor to obtain [km**-1] is unity.
!--------------------------------------------------------------------- 
 
!      abscoeff(:,:,:,1) = ( k_liq(:,:,:)*conc_drop(:,:,:) +       &
       abscoeff(:,:,:  ) = ( k_liq(:,:,:)*conc_drop(:,:,:) +       &
                             k_ice(:,:,:)*conc_ice (:,:,:))/diffac
      
!---------------------------------------------------------------------



end subroutine cloud_lwem_oneband




!######################################################################
! <SUBROUTINE NAME="el">
!  <OVERVIEW>
!   Subroutine to calculates total optical depth and scattering 
!   optical depth
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine calculates total optical depth and scattering 
!     optical depth
!     for infrared radiation using Fu and Liou (1993,
!     JAS). To be used for crystal effective sizes from 20 to 130 um.
!     limits changed to 18.6 to 130.2 um on 2 April 1999 to
!     match shortwave limits.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call el  (conc_ice    , size_ice      ,                   &
!                cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
!  </TEMPLATE>
!  <IN NAME="conc_ice" TYPE="real">
!   the ice crystal concentation in grams / meter**3
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   the ice crystal effective size in microns
!  </IN>
!  <OUT NAME="cldextbndicelw" TYPE="real">
!   the specified values of the extinction          
!                 coefficient for ice particles in kilometers**(-1)
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldssalbbndicelw" TYPE="real">
!   the specified values of the single-           
!                  scattering albedo for ice particles             
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldasymmbndicelw" TYPE="real">
!   the specified values of the asymmetry         
!                  factor for ice particles                        
!                 over wavenumber bands used by the radiation code
!  </OUT>
! </SUBROUTINE>
!
subroutine el  (nb, conc_ice, size_ice, dge_column, mask,  &
                cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
 
!-----------------------------------------------------------------------
!    subroutine el calculates total optical depth and scattering optical
!    depth for infrared radiation using Fu and Liou (1993, JAS).  the
!    parameterization will be used for crystal effective sizes between
!    18.6 and 130.2 microns.
!    Leo Donner, GFDL, 29 Aug 98
!-----------------------------------------------------------------------

integer,                 intent(in)    ::  nb
real, dimension (:,:,:),   intent(in)    ::  conc_ice, size_ice
logical, dimension (:,:,:),   intent(in )    ::  dge_column        
logical, dimension (:,:,:),   intent(out)    ::  mask              
real, dimension (:,:,:  ), intent(out)   ::  cldextbndicelw,   &
                                             cldssalbbndicelw,    &
                                             cldasymmbndicelw

!----------------------------------------------------------------------
!  intent(in) variables:                                            
!                                                                  
!       conc_ice           ice crystal concentation [ grams / meter**3 ]
!       size_ice           ice crystal effective size [ microns ]      
!                                                                  
!  intent(out) variables:                                           
!                                                                  
!       cldextbndicelw     values of the extinction coefficient for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code [ km**(-1) ]
!       cldssalbbndicelw   values of the single-scattering albedo for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code  [ dimensionless ] 
!       cldasymmbndicelw   values of the asymmetry factor for ice 
!                          crystals over the wavenumber bands used by 
!                          the radiation code  
!
!--------------------------------------------------------------------- 
 
                                                                   
!---------------------------------------------------------------------
! local variables:                                                   
      real, dimension (size(conc_ice,1), size(conc_ice,2), &
                       size(conc_ice,3))                    ::   &
                                                    size_i
      integer     :: n
      integer     :: i, j,k
      real  ::          cldextivlice, cldssalbivlice
!     real  ::          cldasymmivlice

      real     ::               sumext, sumssalb
!     real     ::               sumasymm

      real, dimension (NBFL)  ::   a0, a1, a2
 
      data a0 /                                                      &
          -7.752E-03,  -1.741E-02,  -1.704E-02,  -1.151E-02,         &
          -1.026E-02,  -8.294E-03,  -1.153E-02,  -9.609E-03,         &
          -9.061E-03,  -8.441E-03,  -8.088E-03,  -7.770E-03/
      data a1 /                                                      &
           4.624,   5.541,   4.830,   4.182,   4.105,   3.925,       &
           4.109,   3.768,   3.741,   3.715,   3.717,   3.734/
      data a2 /                                                      &
         -42.010, -58.420,  16.270,  31.130,  16.360,   1.315,       &
          17.320,  34.110,  26.480,  19.480,  17.170,  11.850/

      real, dimension (1:NBFL,0:NBB)     ::   b

      data (b(n,0),n=1,NBFL) /                                     &
          0.8079,   0.3964,   0.1028,   0.3254,   0.5207,   0.5631,  &
          0.2307,   0.2037,   0.3105,   0.3908,   0.3014,   0.1996/
      data (b(n,1),n=1,NBFL) /                                     &
         -0.7004E-02, -0.3155E-02,  0.5019E-02,  0.3434E-02,         &
         -0.9778E-03, -0.1434E-02,  0.3830E-02,  0.4247E-02,         &
          0.2603E-02,  0.1272E-02,  0.2639E-02,  0.3780E-02/
      data (b(n,2),n=1,NBFL) /                                     &
          0.5209E-04,  0.6417E-04, -0.2024E-04, -0.3081E-04,         &
          0.3725E-05,  0.6298E-05, -0.1616E-04, -0.1810E-04,         &
          -0.1139E-04, -0.5564E-05, -0.1116E-04, -0.1491E-04/
      data (b(n,3),n=1,NBFL) /                                     &
         -0.1425E-06, -0.2979E-06,  0.0000E+00,  0.9143E-07,         &
          0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,         &
          0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00/

!     real, dimension (1:NBFL,0:NBC)     ::   cpr

!     data (cpr(n,0),n=1,NBFL) /                                   &
!         0.2292,   0.7024,   0.7290,   0.7678,   0.8454,   0.9092,  &
!         0.9167,   0.8815,   0.8765,   0.8915,   0.8601,   0.7955/
!     data (cpr(n,1),n=1,NBFL) /                                   &
!          1.724E-02,   4.581E-03,   2.132E-03,   2.571E-03,         &
!          1.429E-03,   9.295E-04,   5.499E-04,   9.858E-04,         &
!          1.198E-03,   1.060E-03,   1.599E-03,   2.524E-03/
!     data (cpr(n,2),n=1,NBFL) /                                   &
!         -1.573E-04,  -3.054E-05,  -5.584E-06,  -1.041E-05,         &
!         -5.859E-06,  -3.877E-06,  -1.507E-06,  -3.116E-06,         &
!         -4.485E-06,  -4.171E-06,  -6.465E-06,  -1.022E-05/
!     data (cpr(n,3),n=1,NBFL) /                                   &
!          4.995E-07,   6.684E-08,   0.000E+00,   0.000E+00,         &
!          0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,         &
!          0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/


!---------------------------------------------------------------------
!   local variables:
!
!      cldextivlice     cloud extinction coefficient over the spectral
!                       intervals relevant to fu and liou (1993) ice 
!                       crystals. 
!                       [ km**(-1)]
!      cldssalbivlice   cloud single scattering albedo over the spectral
!                       intervals relevant to fu and liou (1993) ice 
!                       crystals. 
!                       [ non-dimensional ]
!      cldasymmivlice   cloud asymmetry factor over the spectral
!                       intervals relevant to fu and liou (1993) ice 
!                       crystals. not currently used.
!                       [ non-dimensional ]
!      sumext           weighted sum of extinction coefficient over fu 
!                       bands to produce value for lw radiation bands
!                       [ kilometers**(-1) ]
!      sumssalb         sum of single scattering albedo over fu bands 
!                       to produce value for lw radiation bands 
!                       [ dimensionless ]
!      sumasymm         sum of asymmetry factor over fu bands 
!                       to produce value for lw radiation bands 
!                       [ dimensionless ]
!      a0,a1,a2         empirical coefficients for extinction 
!                       coefficient parameterization
!      b                empirical coefficients for single scattering 
!                       albedo parameterization
!      cpr              empirical coefficients for asymmetry parameter
!                       parameterization
!      n,ni             do-loop indices
!
!---------------------------------------------------------------------
      do k=1,size(conc_ice,3)
        do j=1,size(conc_ice,2)
          do i=1,size(conc_ice,1)
            if (size_ice(i,j,k) < min_cld_ice_size) then        
              size_i(i,j,k) = min_cld_ice_size             
            else if(size_ice(i,j,k) > max_cld_ice_size) then    
              size_i(i,j,k) = max_cld_ice_size       
            else
              size_i(i,j,k) = size_ice(i,j,k)
            endif
          end do
         end do
      end do

      do k=1, size(conc_ice,3)
        do j=1, size(conc_ice,2)
          do i=1, size(conc_ice,1)
            if ( .not. dge_column(i,j,k)) then
            sumext = 0.
            sumssalb = 0.
            if (conc_ice(i,j,k) /= 0.0) then

              mask(i,j,k) = .true.
!-----------------------------------------------------------------------
!    calculate extinction coefficient [ km**(-1) ] over the wavenumber
!    bands of the Fu-Liou parameterization (not the radiation code
!    wavenumber bands).
!-----------------------------------------------------------------------
              do n=1,NBFL                                               
                cldextivlice  = 1.0E+03*conc_ice(i,j,k)*        &
                                (a0(n) +                              &
                                 a1(n)/size_i(i,j,k) +              &
                                 a2(n)/size_i(i,j,k)**2)

!-----------------------------------------------------------------------
!    calculate single-scattering albedo and asymmetry parameter.
!    the asymmetry parameter is not currently used in the infrared 
!    code. therefore its calculation is commented out.
!-----------------------------------------------------------------------
                cldssalbivlice  = 1.0E+00 -                           &
                                  (b(n,0) +                           &
                                   b(n,1)*size_i(i,j,k) +           &
                                   b(n,2)*size_i(i,j,k)**2 +        &
                                   b(n,3)*size_i(i,j,k)**3)
!               cldasymmivlice  =                            &
!                                  cpr(n,0) +                         &
!                                  cpr(n,1)*size_i(:,:,:) +         &
!                                  cpr(n,2)*size_i(:,:,:)**2 +      &
!                                  cpr(n,3)*size_i(:,:,:)**3
 
!-----------------------------------------------------------------------
!    use the band weighting factors computed in microphys_rad_init
!    to define the radiation band values for the scattering parameters.
!-----------------------------------------------------------------------
                sumext    = sumext + cldextivlice*fulwwts(nb,n )
                sumssalb  = sumssalb + cldssalbivlice*fulwwts(nb,n )
!               sumasymm  = sumasymm + cldasymmivlice*fulwwts(nb,n )
              end do
            else 
              mask(i,j,k) = .false.
        cldextbndicelw(i,j,k) = 0.0        
        cldssalbbndicelw(i,j,k) = 0.0          
!       cldasymmbndicelw(:,:,:,n) = sumasymm(:,:,:)
            endif
            cldextbndicelw(i,j,k)   = sumext       
            cldssalbbndicelw(i,j,k) = sumssalb       
!           cldasymmbndicelw(i,j,k,n) = sumasymm        
             else
                mask(i,j,k) = .false.
            cldextbndicelw(i,j,k)   = 0.0          
            cldssalbbndicelw(i,j,k) = 0.0            
!           cldasymmbndicelw(i,j,k,n) = sumasymm        
             endif
          end do
        end do
      end do

!--------------------------------------------------------------------
 


end subroutine el




!#####################################################################
! <SUBROUTINE NAME="el_dge">
!  <OVERVIEW>
!   Subroutine to calculate total optical depth and scattering optical
!   depth in infrared.
!  </OVERVIEW>
!  <DESCRIPTION>
!   calculates total optical depth and scattering optical depth
!     for infrared radiation using Fu et al (J. Clim., 11,2223 (1998)).
!     To be used for crystal generalized effective diameters from 
!     18.6 to 130.2 um to match shortwave limits.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call el_dge( conc_ice    , size_ice      ,                &
!                cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
!  </TEMPLATE>
!  <IN NAME="conc_ice" TYPE="real">
!   the ice crystal concentation in grams / meter**3
!  </IN>
!  <IN NAME="size_ice" TYPE="real">
!   the ice crystal effective size in microns
!  </IN>
!  <OUT NAME="cldextbndicelw" TYPE="real">
!   the specified values of the extinction          
!                 coefficient for ice particles in kilometers**(-1)
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldssalbbndicelw" TYPE="real">
!   the specified values of the single-           
!                  scattering albedo for ice particles             
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldasymmbndicelw" TYPE="real">
!   the specified values of the asymmetry         
!                  factor for ice particles                        
!                 over wavenumber bands used by the radiation code
!  </OUT>
! </SUBROUTINE>
!  
 
subroutine el_dge (nb, conc_ice, size_ice, dge_column, mask, &
                   cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
 
!-----------------------------------------------------------------------
!    subroutine el_dge calculates the total optical depth and scattering
!    optical depth for infrared radiation using Fu et al (J. Clim., 11,
!    2223 (1998)). this parameterization will be used for crystal 
!    generalized effective diameters from 18.6 to 130.2 um to match 
!    shortwave limits.
!    Leo Donner, GFDL, 29 Aug 98
!    Dan Schwarzkopf, GFDL 31 July 2001
!-----------------------------------------------------------------------

integer,                   intent(in)    ::  nb
real, dimension (:,:,:),   intent(in)    ::  conc_ice, size_ice
logical, dimension (:,:,:),   intent(in)    ::  dge_column
logical, dimension (:,:,:),   intent(out)    ::  mask       
real, dimension (:,:,:  ), intent(out)   ::  cldextbndicelw,   &
                                             cldssalbbndicelw,    &
                                             cldasymmbndicelw

!----------------------------------------------------------------------
!   intent(in) variables:                                            
!                                                                  
!       conc_ice           ice crystal concentation [ grams / meter**3 ]
!       size_ice           ice crystal effective size [ microns ]      
!                                                                  
!  intent(out) variables:                                           
!                                                                  
!       cldextbndicelw     values of the extinction coefficient for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code [ km**(-1) ]
!       cldssalbbndicelw   values of the single-scattering albedo for 
!                          ice crystals over the wavenumber bands used 
!                          by the radiation code  [ dimensionless ] 
!       cldasymmbndicelw   values of the asymmetry factor for ice 
!                          crystals over the wavenumber bands used by 
!                          the radiation code  
!
!--------------------------------------------------------------------- 
 
!---------------------------------------------------------------------
!  local variables:                                                   
      real, dimension (size(conc_ice,1), size(conc_ice,2), &
                       size(conc_ice,3))                    ::   &
                                                    size_i
      integer     :: n, m
      integer     :: i,j,k
 
      real                        ::    cldextivlice, cldabsivlice, &
                                        cldssalbivlice 
!     real                        ::    cldasymmivlice 

      real    ::        sumext, sumssalb
!     real    ::        sumasymm
      real    ::        fw1,fw2,fw3,fw4

      real, dimension (NBFL)  ::   a0, a1, a2
 
      data a0 /                                                       &
           4.919685E-03,  3.325756E-03, -1.334860E-02, -9.524174E-03, &
          -4.159424E-03, -1.691632E-03, -8.372696E-03, -8.178608E-03, &
          -4.936610E-03, -3.034573E-03, -2.465236E-03, -2.308881E-03/
      data a1 /                                                     &
           2.327741, 2.601360, 4.043808, 3.587742, 3.047325, 2.765756, &
           3.455018, 3.401245, 3.087764, 2.900043, 2.833187, 2.814002  /
      data a2 /                                                    &
          -1.390858E+01, -1.909602E+01, -2.171029E+01, -1.068895E+01, &
          -5.061568E+00, -8.331033E+00, -1.516692E+01, -8.812820E+00, &
          -3.884262E+00, -1.849911E+00, -4.227573E-01,  1.072211E+00/

      real, dimension (1:NBFL,0:NBB)     ::   b
      data (b(n,0),n=1,NBFL) /                                     &
          8.869787E-01,  2.005578E-01,  3.003701E-01,  9.551440E-01, &
          1.466481E+00,  1.195515E+00,  5.409536E-01,  5.874323E-01, &
          7.152274E-01,  8.862434E-01,  7.428957E-01,  4.346482E-01/
      data (b(n,1),n=1,NBFL) /                                     &
          2.118409E-02,  2.132614E-02,  2.051529E-02,  1.309792E-02, &
         -2.129226E-03,  3.350616E-03,  1.949649E-02,  1.876628E-02, &
          1.621734E-02,  1.226538E-02,  1.279601E-02,  1.721457E-02/
      data (b(n,2),n=1,NBFL) /                                     &
         -2.781429E-04, -1.751052E-04, -1.931684E-04, -1.793694E-04, &
         -1.361630E-05, -5.266996E-05, -2.050908E-04, -2.045834E-04, &
         -1.868544E-04, -1.523076E-04, -1.391803E-04, -1.623227E-04/
      data (b(n,3),n=1,NBFL) /                                     &
          1.094562E-06,  5.355885E-07,  6.583031E-07,  7.313392E-07, &
          1.193649E-07,  2.233377E-07,  7.364680E-07,  7.510080E-07, &
          7.078738E-07,  6.000892E-07,  5.180104E-07,  5.561523E-07/

!     real, dimension (1:NBFL,0:NBC)     ::   cpr
!     data (cpr(n,0),n=1,NBFL) /                                   &
!         4.949276E-01,  6.891414E-01,  7.260484E-01,  7.363466E-01, &
!         7.984021E-01,  8.663385E-01,  8.906280E-01,  8.609604E-01, &
!         8.522816E-01,  8.714665E-01,  8.472918E-01,  7.962716E-01/
!     data (cpr(n,1),n=1,NBFL) /                                   &
!         1.186174E-02,  6.192281E-03,  2.664334E-03,  4.798266E-03, &
!         3.977117E-03,  2.797934E-03,  1.903269E-03,  2.200445E-03, &
!         2.523627E-03,  2.455409E-03,  2.559953E-03,  3.003488E-03/
!     data (cpr(n,2),n=1,NBFL) /                                   &
!        -1.267629E-04, -6.459514E-05, -1.251136E-05, -4.513292E-05, &
!        -4.471984E-05, -3.187011E-05, -1.733552E-05, -1.748105E-05, &
!        -2.149196E-05, -2.456935E-05, -2.182660E-05, -2.082376E-05/
!     data (cpr(n,3),n=1,NBFL) /                                   &
!         4.603574E-07,  2.436963E-07,  2.243377E-08,  1.525774E-07, &
!         1.694919E-07,  1.217209E-07,  5.855071E-08,  5.176616E-08, &
!         6.685067E-08,  8.641223E-08,  6.879977E-08,  5.366545E-08/


!+yim small dge
      real, dimension (NBA2,NBFL)     ::   aa
      real, dimension (NBB2,NBFL)     ::   bb
!     real, dimension (NBC2,NBFL)     ::   cc
 
        data aa / &
      -2.005187e+00, 1.887870e+00, -2.394987e-01, 1.226004e-02, &
      -2.176747e-04, &
      -1.221428e+00, 1.190519e+00, -1.081918e-01, 3.207774e-03, &
      -7.790185e-06, &
      -5.522210e-01, 5.556264e-01, 1.350808e-02, -5.182503e-03, &
       1.854760e-04, &
      -2.411192e-01, 2.109769e-01, 7.588264e-02, -9.103300e-03, &
       2.678349e-04, &
      -1.485194e-02, 4.630892e-03, 8.989527e-02, -8.569112e-03, &
       2.290338e-04, &
       4.292661e-02, -7.619363e-04, 5.089112e-02, -4.101744e-03, &
       9.917537e-05, &
      -1.257657e-03, 3.840350e-01, -2.336758e-02, 5.263245e-04, &
       9.536367e-07, &
      -2.482977e-01, 5.149985e-01, -1.086854e-02, -1.909389e-03, &
       8.220600e-05, &
       1.130811e-01, -7.663294e-02, 9.961269e-02, -8.920452e-03, &
       2.325299e-04, &
       1.477471e-01, -1.276555e-01, 5.889066e-02, -3.637540e-03, &
       7.242738e-05, &
       2.778228e-02, 9.410452e-03, 7.771632e-03, -1.847559e-05, &
      -7.178001e-06, &
       2.954018e-03, 1.254725e-01, -3.265442e-03, 2.270727e-04, &
      -6.365789e-06 /
        data bb / &
      -8.768658e-03, 8.493330e-02, -3.632126e-03, 6.987413e-05, &
       2.703965e-07, &
      -7.762272e-03, 1.653825e-01, -1.242696e-02, 4.813596e-04, &
      -6.987702e-06, &
      -1.103846e-02, 1.880946e-01, -1.320780e-02, 4.530029e-04, &
      -5.384886e-06, &
      -1.240034e-02, 1.353184e-01, -6.773254e-03, 1.353446e-04, &
       4.783046e-07, &
      -9.834148e-03, 1.045283e-01, -3.714625e-03, 9.185834e-06, &
       2.434297e-06, &
      -4.989783e-03, 9.761852e-02, -3.464011e-03, 1.681863e-05, &
       1.990612e-06, &
       5.524896e-02, 3.828618e-01, -4.868927e-02, 2.788080e-03, &
      -5.893696e-05, &
      -1.102297e-01, 4.983548e-01, -5.947312e-02, 3.147713e-03, &
      -6.196981e-05, &
      -3.705134e-02, 1.612865e-01, -4.132244e-03, -2.863781e-04, &
       1.374847e-05, &
       5.730367e-03, 3.433887e-02, 3.147511e-03, -3.044807e-04, &
       7.929481e-06, &
       3.126666e-03, 3.533685e-02, 5.299923e-04, -6.735890e-05, &
       1.687872e-06, &
       9.549627e-03, 1.140347e-01, 1.223725e-03, -4.282989e-04, &
       1.343652e-05 /
!       data cc / &
!     -1.592086e-01, 5.165795e-01, -8.889665e-02, 6.133364e-03, &
!     -1.466832e-04, &
!     -2.780309e-01, 5.589181e-01, -9.294043e-02, 6.316572e-03, &
!     -1.501642e-04, &
!     -4.146218e-01, 6.015844e-01, -9.714942e-02, 6.513667e-03, &
!     -1.539503e-04, &
!     -4.644106e-01, 5.861063e-01, -9.087172e-02, 5.917403e-03, &
!     -1.371181e-04, &
!     -4.848736e-01, 5.552414e-01, -8.183047e-02, 5.147920e-03, &
!     -1.164374e-04, &
!     -5.056360e-01, 5.240870e-01, -7.278649e-02, 4.395703e-03, &
!     -9.639759e-05, &
!     -4.991806e-01, 4.601579e-01, -5.805338e-02, 3.236112e-03, &
!     -6.615910e-05, &
!     -4.382576e-01, 3.812485e-01, -4.268756e-02, 2.088357e-03, &
!     -3.689533e-05, &
!     -3.094784e-01, 2.406058e-01, -1.477957e-02, 2.970087e-05, &
!      1.421683e-05, &
!     -9.731071e-02, 4.088258e-02, 2.106015e-02, -2.364895e-03, &
!      6.892137e-05, &
!      7.192725e-02, -8.649291e-02, 3.621089e-02, -2.888238e-03, &
!      7.087982e-05, &
!      6.792641e-02, -5.575384e-02, 1.319878e-02, -4.919461e-04, &
!      8.543384e-07 /



!---------------------------------------------------------------------
!   local variables:
!
!      cldextivlice     cloud extinction coefficient over the spectral
!                       intervals relevant to fu and liou (1998) ice 
!                       crystals. 
!                       [ km**(-1)]
!      cldabsivlice     absorption coefficient over the spectral
!                       intervals relevant to fu and liou (1998) ice 
!                       crystals. 
!                       [ km**(-1)]
!      cldssalbivlice   cloud single scattering albedo over the spectral
!                       intervals relevant to fu and liou (1998) ice 
!                       crystals. 
!                       [ non-dimensional ]
!      cldasymmivlice   cloud asymmetry factor over the spectral
!                       intervals relevant to fu and liou (1998) ice 
!                       crystals. not currently used.
!                       [ non-dimensional ]
!      sumext           weighted sum of extinction coefficient over fu 
!                       bands to produce value for lw radiation bands
!                       [ kilometers**(-1) ]
!      sumssalb         sum of single scattering albedo over fu bands 
!                       to produce value for lw radiation bands 
!                       [ dimensionless ]
!      sumasymm         sum of asymmetry factor over fu bands 
!                       to produce value for lw radiation bands 
!                       [ dimensionless ]
!      a0,a1,a2         empirical coefficients for extinction coef-
!                       ficient parameterization
!      b                empirical coefficients for single scattering 
!                       albedo parameterization
!      cpr              empirical coefficients for asymmetry parameter
!                       parameterization
!      n,ni             do-loop indices
! 
!----------------------------------------------------------------------

      do k=1,size(conc_ice,3)
        do j=1,size(conc_ice,2)
          do i=1,size(conc_ice,1)
            if (size_ice(i,j,k) < min_cld_ice_size) then       
              size_i(i,j,k) = min_cld_ice_size             
            else if(size_ice(i,j,k) > max_cld_ice_size ) then  
              size_i(i,j,k) = max_cld_ice_size       
            else
              size_i(i,j,k) = size_ice(i,j,k)
            endif
          end do
         end do
      end do

      do k=1, size(conc_ice,3)
        do j=1, size(conc_ice,2)
          do i=1, size(conc_ice,1)
            if (dge_column(i,j,k)) then
              sumext = 0.
              sumssalb = 0.
              if (conc_ice(i,j,k) /= 0.0) then
                mask (i,j,k)= .true.

                if ((Cldrad_control%using_fu2007 .and.    &
                    size_i(i,j,k) > 15.0) .or. &
                    .not. Cldrad_control%using_fu2007) then

!-----------------------------------------------------------------------
!    calculate extinction coefficient [km**(-1)] for each wavenumber
!    band of the Fu-Liou parameterization (not the radiation
!    code wavenumber bands).
!-----------------------------------------------------------------------
                  do n=1,NBFL                                   
                    cldextivlice          = 1.0E+03*conc_ice(i,j,k)*  &
                                  (a0(n) +                           &
                                   a1(n)*(1.0/size_i(i,j,k)) +     &
                                   a2(n)*(1.0/size_i(i,j,k)**2))

!-----------------------------------------------------------------------
!    calculate the absorption coefficient. convert to units of 
!    [ km**(-1) ].
!-----------------------------------------------------------------------
                    cldabsivlice          = 1.0E+03*conc_ice(i,j,k)*   &
                                  (1.0/size_i(i,j,k))*        &       
                                   (b(n,0) +                    &
                                    b(n,1)*size_i(i,j,k) +    &
                                    b(n,2)*size_i(i,j,k)**2 + &
                                    b(n,3)*size_i(i,j,k)**3)

!---------------------------------------------------------------------
!    compute the single-scattering albedo. the asymmetry parameter is
!    not currently used in the infrared code, so its calculation is
!    commented out.
!-----------------------------------------------------------------------
                    if (cldextivlice /= 0.0) then
                      cldssalbivlice = 1.0E+00 -    &
                                             cldabsivlice/cldextivlice
                    else
                      cldssalbivlice = 0.0
                    endif 
 
!                   do n=1,NBFL                                        
!                     cldasymmivlice(:,:,:,n) = cpr(n,0) +        &
!                                 cpr(n,1)*size_i(:,:,:) +        &
!                                 cpr(n,2)*size_i(:,:,:)**2 +     &
!                                 cpr(n,3)*size_i(:,:,:)**3
!                   end do
 
!-----------------------------------------------------------------------
!    use the band weighting factors computed in microphys_rad_init
!    to define the values of these parameters for each lw radiation
!    band.
!-----------------------------------------------------------------------
!!                  sumasymm(:,:,:) = 0.
                    sumext        = sumext +         &
                          cldextivlice*fulwwts(nb,n)
                    sumssalb = sumssalb +     &
                            cldssalbivlice*fulwwts(nb,n)
!!                  sumasymm(:,:,:) = sumasymm(:,:,:) +     &
!!                          cldasymmivlice(:,:,:,ni)*fulwwts(n,ni)
                  end do
                else ! ( using_fu2007 .and. size_i > 15.0 .or. not
                     !   using_fu2007)

 !+yim small dge
!-----------------------------------------------------------------------
!    calculate extinction coefficient [km**(-1)] for each wavenumber
!    band of the Fu-Liou parameterization (not the radiation
!    code wavenumber bands).
!-----------------------------------------------------------------------
                  do n=1,NBFL                          
                    m = NBFL - n + 1
                    fw1 = size_i(i,j,k)
                    fw2 = fw1*size_i(i,j,k)
                    fw3 = fw2*size_i(i,j,k)
                    fw4 = fw3*size_i(i,j,k)
       
                    cldextivlice   = 1.0E+03*conc_ice(i,j,k)/fw1*  &
                                 (aa(1,m) +       &
                                  aa(2,m)*fw1 +   &
                                  aa(3,m)*fw2 +   &
                                  aa(4,m)*fw3 +   &
                                  aa(5,m)*fw4 )

!-----------------------------------------------------------------------
!    calculate the absorption coefficient. convert to units of 
!    [ km**(-1) ].
!-----------------------------------------------------------------------
                    cldabsivlice    = 1.0E+03*conc_ice(i,j,k)*      &
                                 (1.0/fw1)*        &
                                  (bb(1,m) +                    &
                                   bb(2,m)*fw1 +  &
                                   bb(3,m)*fw2 +  &
                                   bb(4,m)*fw3 +  &
                                   bb(5,m)*fw4 )

!---------------------------------------------------------------------
!    compute the single-scattering albedo. the asymmetry parameter is
!    not currently used in the infrared code, so its calculation is
!    commented out.
!-----------------------------------------------------------------------
                    if (cldextivlice /= 0.0) then
                      cldssalbivlice = 1.0E+00 -    &
                                            cldabsivlice/cldextivlice
                    else
                      cldssalbivlice = 0.0
                    endif

!                   do n=1,NBFL                           
!                     m = NBFL - n + 1
!                     cldasymmivlice(:,:,:,n) = cc(1,m) +        &
!                                 cc(2,m)*fw1 +     &
!                                 cc(3,m)*fw2 +     &
!                                 cc(4,m)*fw3 +     &
!                                 cc(5,m)*fw4
!                   end do
 
!-----------------------------------------------------------------------
!    use the band weighting factors computed in microphys_rad_init
!    to define the values of these parameters for each lw radiation
!    band.
!-----------------------------------------------------------------------
!!                  sumasymm(:,:,:) = 0.
                    sumext  = sumext + cldextivlice*fulwwts(nb,n)
                    sumssalb = sumssalb + cldssalbivlice*fulwwts(nb,n)
!!                  sumasymm(:,:,:) = sumasymm(:,:,:) +     &
!!                               cldasymmivlice(:,:,:,ni)*fulwwts(n,ni)
                  end do
                endif ! (size > 15)    
                cldextbndicelw(i,j,k) = sumext
                cldssalbbndicelw(i,j,k) = sumssalb
!               cldasymmbndicelw(:,:,:,n) = sumasymm(:,:,:)
              else  ! (conc_ice > 0)
                mask(i,j,k) = .false.
                cldextbndicelw(i,j,k) = 0.0        
                cldssalbbndicelw(i,j,k) = 0.0          
!               cldasymmbndicelw(:,:,:,n) = sumasymm(:,:,:)
              endif ! (convc_ice > 0)
              cldextbndicelw(i,j,k) = sumext
              cldssalbbndicelw(i,j,k) = sumssalb
!             cldasymmbndicelw(:,:,:,n) = sumasymm(:,:,:)
            else  ! (dge_column)
              mask(i,j,k) = .false.
              cldextbndicelw(i,j,k) = 0.0        
              cldssalbbndicelw(i,j,k) = 0.0          
!             cldasymmbndicelw(:,:,:,n) = sumasymm(:,:,:)
            endif !(dge_column)
          end do
        end do
      end do

!---------------------------------------------------------------------


end subroutine el_dge



!#####################################################################
! <SUBROUTINE NAME="cliqlw">
!  <OVERVIEW>
!   Subroutine to calculate longwave absorption optical depth for liquid
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculates longwave absorption optical depth for liquid.
!     Follows Held et al. (J. Atmos. Sci., 1993).
!     
!     Leo Donner, GFDL, 1 Feb 1999
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cliqlw (conc_drop, cldextbnddroplw)
!  </TEMPLATE>
!  <IN NAME="conc_drop" TYPE="real">
!   the cloud drop concentration in grams / meter**3
!  </IN>
!  <OUT NAME="cldextbnddroplw" TYPE="real">
!   the specified values of the extinction         
!                 coefficient for cloud drops in kilometers**(-1)  
!                 over wavenumber bands used by the radiation code
!  </OUT>
! </SUBROUTINE>

subroutine cliqlw (conc_drop, cldextbnddroplw)
 
!-----------------------------------------------------------------------
!    subroutine cliqlw calculates the longwave absorption optical depth
!    for cloud drops. the parameterization follows Held et al. 
!    (J. Atmos. Sci., 1993).
!     Leo Donner, GFDL, 1 Feb 1999
!-----------------------------------------------------------------------

real, dimension (:,:,:),   intent(in)   ::  conc_drop
real, dimension (:,:,:), intent(out)  ::  cldextbnddroplw

!---------------------------------------------------------------------
!   intent(in) variable:                                            
!                                                                  
!     conc_drop        cloud drop concentration [ grams / meter**3 ]
!                                                                  
!   intent(out) variable:                                            
!                                                                  
!     cldextbnddroplw  values of the extinction coefficient for 
!                      cloud droplets over the wavenumber bands 
!                      used by the radiation code [ km**(-1) ]
!                                                                  
!-----------------------------------------------------------------------

 
!---------------------------------------------------------------------
!  local variables:                                                   

!     real        ::   alpha = 0.1   (now is module variable set in nml)

!---------------------------------------------------------------------
!  local variables:                                                   
!
!     n           do-loop index
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the cloud droplet extinction coefficient. convert to
!    units of [ km**(-1) ].
!--------------------------------------------------------------------
        cldextbnddroplw(:,:,:  ) = 1.0E+03*alpha*conc_drop(:,:,:)
 
!---------------------------------------------------------------------

end subroutine cliqlw


!#####################################################################
! <SUBROUTINE NAME="furainlw">
!  <OVERVIEW>
!   Subroutine to calculate total optical depth and scattering optical
!   depth in infrared for cloud rain water
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculates absorption coefficient for cloud rain water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for rain water with radii between 60 um and 1.8 mm.
!      See also notes from Q. Fu (4 Sept 98)
!      note: the size of the rain water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!  </DESCRIPTION>
!  <TEMPLATE>
!   subroutine furainlw                                         &
!                (conc_rain    ,                                   &
!                 cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw)
!  </TEMPLATE>
!  <IN NAME="conc_rain" TYPE="real">
!   the rain drop water concentration in grams / meter**3
!  </IN>
!  <OUT NAME="cldextbndrainlw" TYPE="real">
!   the specified values of the extinction          
!                 coefficient for rain water in kilometers**(-1)
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldssalbbndrainlw" TYPE="real">
!   the specified values of the single-           
!                  scattering albedo for rain water             
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldasymmbndrainlw" TYPE="real">
!   the specified values of the asymmetry         
!                  factor for rain water                       
!                 over wavenumber bands used by the radiation code
!  </OUT>
! </SUBROUTINE>
!  

subroutine furainlw (nb, conc_rain, cldextbndrainlw, cldssalbbndrainlw,&
                     cldasymmbndrainlw)
 
!----------------------------------------------------------------------
!    subroutine furainlw calculates the absorption coefficient for 
!    rain water for longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!    it is designed for use with rain drops with radii between 60 um 
!    and 1.8 mm. see also notes from Q. Fu (4 Sept 98). note that the 
!    size of the rain water from the microphysical model (if existing) 
!    does not enter into the calculations.
!    Leo Donner, GFDL, 20 Mar 99
!----------------------------------------------------------------------

integer,                   intent(in)    ::   nb
real, dimension (:,:,:),   intent(in)    ::   conc_rain
real, dimension (:,:,:  ), intent(out)   ::   cldextbndrainlw,    &
                                              cldssalbbndrainlw,   &
                                              cldasymmbndrainlw

!----------------------------------------------------------------------
!   intent(in) variables:                                             
!
!       conc_rain           rain drop concentration [ grams / meter**3 ]
!
!   intent(out) variables:                                            
!
!       cldextbndrainlw     values of the extinction coefficient for 
!                           rain drops over the wavenumber bands used 
!                           by the radiation code [ km**(-1) ]
!       cldssalbbndrainlw   values of the single-scattering albedo for 
!                           rain drops over the wavenumber bands used 
!                           by the radiation code  [ dimensionless ] 
!       cldasymmbndrainlw   values of the asymmetry factor for rain
!                           drops over the wavenumber bands used by 
!                           the radiation code  
!
!--------------------------------------------------------------------- 


!----------------------------------------------------------------------
!  local variables:                                                  

      real ::          cldextivlrain, cldssalbivlrain, cldasymmivlrain
 
      real ::          sumext, sumssalb, sumasymm

      real, dimension (NBFL)    :: brn, wrnf, grn

      data brn /                                                     &
             1.6765,  1.6149,  1.5993,  1.5862,  1.5741,  1.5647,    &
             1.5642,  1.5600,  1.5559,  1.5512,  1.5478,  1.5454/
      data wrnf /                                                    &
              .55218,  .55334,  .55488,  .55169,  .53859,  .51904,   &
              .52321,  .52716,  .52969,  .53192,  .52884,  .53233/
      data grn /                                                     &
              .90729,  .92990,  .93266,  .94218,  .96374,  .98584,   &
              .98156,  .97745,  .97467,  .97216,  .97663,  .97226/

      real      :: rwc0 = 0.5
      integer   :: n
      integer   :: i,j,k

!----------------------------------------------------------------------
!  local variables:                                                  
!
!       cldextivlrain     the specified spectral values of the 
!                         extinction coefficient for rain water 
!                         [ kilometers**(-1) ] 
!       cldssalbivlrain   the specified spectral values of the single-  
!                         scattering albedo for rain water     
!                         [ dimensionless ]
!       cldasymmivlrain   the specified spectral values of the 
!                         asymmetry factor for rain water    
!                         [ dimensionless ]
!       sumext            weighted sum of extinction coefficient over fu
!                         bands to produce value for lw radiation bands
!                         [ kilometers**(-1) ]
!       sumssalb          sum of single scattering albedo over fu bands 
!                         to produce value for lw radiation bands 
!                         [ dimensionless ]
!       sumasymm          sum of asymmetry factor over fu bands 
!                         to produce value for lw radiation bands 
!                         [ dimensionless ]
!        brn              empirical coefficients for extinction 
!                         coefficient parameterization [ (km**-1) ]
!        wrnf             empirical coefficients for single scattering 
!                         albedo parameterization
!                         [ dimensionless ]
!        grn              empirical coefficients for asymmetry parameter
!                         parameterization
!                         [ dimensionless ]
!        rwc0             rain water content used to obtain above
!                         empirical coefficients.  [ g / m**3 ] 
!        n, ni            do-loop indices         
!
!---------------------------------------------------------------------
      do k=1, size(conc_rain,3)
        do j=1, size(conc_rain,2)
          do i=1, size(conc_rain,1)
            sumext = 0.
            sumssalb = 0.
            sumasymm = 0.
            if (conc_rain(i,j,k) /= 0.0) then
 
!-----------------------------------------------------------------------
!    calculate extinction coefficient (km**(-1)) over the wavenumber
!    bands of the Fu-Liou parameterization (not the radiation code
!    wavenumber bands). define the single-scattering albedo for each
!    band. the asymmetry factor is not currently used, so the code
!    defining it is commented out.
!-----------------------------------------------------------------------
              do n=1,NBFL
                cldextivlrain   = brn(n)*conc_rain(i,j,k)/rwc0
                cldssalbivlrain = wrnf(n)
                cldasymmivlrain = grn(n)
 
!-----------------------------------------------------------------------
!    use the band weighting factors computed in microphys_rad_init
!    to define the values of these parameters for each lw radiation
!    band.
!-----------------------------------------------------------------------
                sumext   = sumext + cldextivlrain*fulwwts(nb,n )
                sumssalb = sumssalb + cldssalbivlrain*fulwwts(nb,n )
                sumasymm = sumasymm + cldasymmivlrain*fulwwts(nb,n)
              end do
            endif
            cldextbndrainlw  (i,j,k  ) = sumext
            cldssalbbndrainlw(i,j,k  ) = sumssalb
            cldasymmbndrainlw(i,j,k  ) = sumasymm
          end do
        end do
      end do

!---------------------------------------------------------------------
 
 
end subroutine furainlw



!#####################################################################
! <SUBROUTINE NAME="fusnowlw">
!  <OVERVIEW>
!   Subroutine to calculate total optical depth and scattering optical
!   depth in infrared for cloud rain water
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculates absorption coefficient for cloud rain water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for rain water with radii between 60 um and 1.8 mm.
!      See also notes from Q. Fu (4 Sept 98)
!      note: the size of the rain water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!  </DESCRIPTION>
!  <TEMPLATE>
!   subroutine fusnowlw                                         &
!                (conc_snow    ,                                   &
!                 cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw)
!  </TEMPLATE>
!  <IN NAME="conc_snow" TYPE="real">
!   the snow drop water concentration in grams / meter**3
!  </IN>
!  <OUT NAME="cldextbndsnowlw" TYPE="real">
!   the specified values of the extinction          
!                 coefficient for snow drop water in kilometers**(-1)
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldssalbbndsnowlw" TYPE="real">
!   the specified values of the single-           
!                  scattering albedo for snow drop water             
!                 over wavenumber bands used by the radiation code
!  </OUT>
!  <OUT NAME="cldasymmbndsnowlw" TYPE="real">
!   the specified values of the asymmetry         
!                  factor for snow drop water                       
!                 over wavenumber bands used by the radiation code
!  </OUT>
! </SUBROUTINE>
!
subroutine fusnowlw (nb, conc_snow, cldextbndsnowlw, cldssalbbndsnowlw,&
                     cldasymmbndsnowlw)
 
!-----------------------------------------------------------------------
!    subroutine fusnowlw calculates the absorption coefficient for cloud
!    snow water for longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!    it is relevant for snow flakes with radii between 60 um and 5.0 mm.
!    see also notes from Q. Fu (4 Sept 98). note that the size of the 
!    snow flakes from the microphysical model (if existing) does not
!    enter into the calculations.
!    Leo Donner, GFDL, 20 Mar 99
!
!-----------------------------------------------------------------------
                                                                    
integer,                   intent(in)     ::   nb
real, dimension (:,:,:),   intent(in)     ::   conc_snow
real, dimension (:,:,:),   intent(out)    ::   cldextbndsnowlw,    &
                                               cldssalbbndsnowlw,  &
                                               cldasymmbndsnowlw

!----------------------------------------------------------------------
!   intent(in) variables:                                             
!
!       conc_snow           snow concentration [ grams / meter**3 ]
!
!   intent(out) variables:                                            
!
!       cldextbndsnowlw     values of the extinction coefficient for 
!                           snow flakes over the wavenumber bands used 
!                           by the radiation code [ km**(-1) ]
!       cldssalbbndsnowlw   values of the single-scattering albedo for 
!                           snow flakes over the wavenumber bands used 
!                           by the radiation code  [ dimensionless ] 
!       cldasymmbndsnowlw   values of the asymmetry factor for snow
!                           flakes over the wavenumber bands used by 
!                           the radiation code  
!
!--------------------------------------------------------------------- 
 
!----------------------------------------------------------------------
!  local variables:                                                  

     real  ::          cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow

      real ::          sumext, sumssalb, sumasymm

      real, dimension (NBFL)     :: brn, wrnf, grn
      data brn /                                                     &
              .87477,  .85421,  .84825,  .84418,  .84286,  .84143,   &
              .84097,  .84058,  .84029,  .83995,  .83979,  .83967/
      data wrnf /                                                    &
              .55474,  .53160,  .54307,  .55258,  .54914,  .52342,   &
              .52446,  .52959,  .53180,  .53182,  .53017,  .53296/
      data grn /                                                     &
              .93183,  .97097,  .95539,  .94213,  .94673,  .98396,   &
              .98274,  .97626,  .97327,  .97330,  .97559,  .97173/
                                                                    
      real      :: swc0 = 0.5
      integer   :: n
      integer   :: i,j,k
 
!---------------------------------------------------------------------
!   local variables:
!
!       cldextivlsnow     the specified spectral values of the 
!                         extinction coefficient for snow flakes
!                         [ kilometers**(-1) ] 
!       cldssalbivlsnow   the specified spectral values of the single-  
!                         scattering albedo for snow flakes    
!                         [ dimensionless ]
!       cldasymmivlsnow   the specified spectral values of the 
!                         asymmetry factor for snow flakes   
!                         [ dimensionless ]
!       sumext            weighted sum of extinction coefficient over fu
!                         bands to produce value for lw radiation bands
!                         [ kilometers**(-1) ]
!       sumssalb          sum of single scattering albedo over fu bands 
!                         to produce value for lw radiation bands 
!                         [ dimensionless ]
!       sumasymm          sum of asymmetry factor over fu bands 
!                         to produce value for lw radiation bands 
!                         [ dimensionless ]
!       brn               empirical coefficients for extinction 
!                         coefficient parameterization (km**-1)
!       wrnf              empirical coefficients for single scattering 
!                         albedo parameterization [ dimensionless ]
!       grn               empirical coefficients for asymmetry factor   
!                         parameterization [ dimensionless ]
!       swc0              snow water content used to obtain above
!                         empirical coefficients  [ g / m**3 ]
!       n,ni              do-loop indices
!
!---------------------------------------------------------------------

      do k=1, size(conc_snow,3)
        do j=1, size(conc_snow,2)
          do i=1, size(conc_snow,1)
            sumext = 0.
            sumssalb = 0.
            sumasymm = 0.
            if (conc_snow(i,j,k) >= 1.e-5) then

!-----------------------------------------------------------------------
!    calculate the extinction coefficient over the wavenumber bands of 
!    the Fu-Liou parameterization (not the radiation code wavenumber 
!    bands). define the single scattering albedo for each band. the
!    asymmetry factor is not currently used, so the code defining it
!    is commented out.
!-----------------------------------------------------------------------
              do n=1,NBFL
                cldextivlsnow   = brn(n)*conc_snow(i,j,k)/swc0
                cldssalbivlsnow = wrnf(n)
!               cldasymmivlsnow = grn(n)
 
!-----------------------------------------------------------------------
!    use the band weighting factors computed in microphys_rad_init     
!    to define the appropriate values for the scattering parameters for
!    each lw radiation band.
!-----------------------------------------------------------------------
                sumext     = sumext + cldextivlsnow*fulwwts(nb,n )
                sumssalb   = sumssalb + cldssalbivlsnow*fulwwts(nb,n )
!               sumasymm   = sumasymm + cldasymmivlsnow*fulwwts(nb,n)
              end do
            endif
            cldextbndsnowlw(i,j,k)   = sumext
            cldssalbbndsnowlw(i,j,k) = sumssalb         
            cldasymmbndsnowlw(i,j,k) = sumasymm
          end do
        end do
      end do

!----------------------------------------------------------------------
 
 
end subroutine fusnowlw



!####################################################################

!!! THIS SUBOUTINE IS CURRENTLY NOT USED.
subroutine snowlw
!subroutine snowlw(riwp, tau)
!
!-----------------------------------------------------------------------
!
!     Calculates emissivity
!     for longwave radiation using Fu et al. (1995,
!     JAS). (See notes from Kuo-Nan Liou, 1 Sept 98).
!
!-----------------------------------------------------------------------
!
!     On Input:
!
!        riwp   snow path  (g/(m**2))
!
!     On Output:
!
!        tau   absorption optical depth
!
!        Leo Donner, GFDL, 11 Sept 98
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!     Calculate extinction coefficient (m**(-1)) and optical depth
!     for absorption. Extinction coefficient is taken as product of
!     total extinction coefficient (.84 km**(-1)) and single-scattering
!     albedo. See Kuo-Nan Liou notes 
!     (1 Sept 98).
!
!-----------------------------------------------------------------------
!
!     riwc0=0.5
!     ext=.4
!     if (riwp .eq. 0.) then
!       tau=0.
!     else
!       tau=ext*.001*riwp/riwc0
!     end if
!     emis=1.-exp(-tau)

!-----------------------------------------------------------------------


end subroutine snowlw



!#################################################################




                   end module microphys_rad_mod


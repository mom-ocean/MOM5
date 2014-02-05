                     module esfsw_driver_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  smf
! </REVIEWER>
!
!
! <OVERVIEW>
!  Code that initializes and calculates shortwave radiative quantities
!  such as flux and heating rate.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes the necessary shortwave radiative parameters
!  in the initialization subroutine. It then uses delta-eddington approximation
!  and doubling and adding technique to calculate solar flux and
!  heating rate.
! </DESCRIPTION>
!

!    shared modules:

use mpp_mod,              only:  input_nml_file
use fms_mod,              only:  open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, close_file
use constants_mod,        only:  PI, GRAV, radcon_mks, o2mixrat, &
                                 rhoair, pstd_mks, WTMAIR, &
                                 constants_init

!  shared radiation package modules:

use esfsw_parameters_mod, only:  Solar_spect, esfsw_parameters_init, &
                                 esfsw_parameters_end
use rad_utilities_mod,    only:  Rad_control, rad_utilities_init, &
                                 cldrad_properties_type, &
                                 cld_specification_type, &
                                 astronomy_type, &
                                 aerosol_diagnostics_type, &
                                 radiative_gases_type, &
                                 aerosol_type, aerosol_properties_type,&
                                 Cldrad_control, &
                                 atmos_input_type, surface_type, &
                                 sw_output_type, Sw_control
use tracer_manager_mod,   only : get_tracer_index,   &
                                 NO_TRACER
use field_manager_mod,    only : MODEL_ATMOS
!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    esfsw_driver_mod is the internal driver for the esf shortwave
!    package.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: esfsw_driver.F90,v 19.0 2012/01/06 20:15:51 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public    &
         esfsw_driver_init, swresf,   &
         esfsw_driver_end

private     &

!   called from swresf:
         adding, deledd


!---------------------------------------------------------------------
!-------- namelist  ---------

 logical     ::  do_ica_calcs=.false.           ! do independent column
                                                ! calculations when sto-
                                                ! chastic clouds are
                                                ! active ?
logical      ::  do_rayleigh_all_bands = .true. ! rayleigh scattering 
                                                ! calculated in all sw 
                                                ! bands ?
logical      ::  do_herzberg = .false.          ! include the herzberg 
                                                ! effect on the o2 
                                                ! optical depth ?
logical      ::  do_quench = .false.            ! include the quenching
                                                ! effect of non-LTE 
                                                ! processes on the co2 
                                                ! optical depth ?
logical      ::  do_ch4_sw_effects = .false.    ! the shortwave effects
                                                ! of ch4 are included ?
logical      ::  do_n2o_sw_effects = .false.    ! the shortwave effects
                                                ! of n2o are included ?
logical      ::  do_coupled_stratozone = .false. ! include the coupled
                                                 ! stratospheric ozone effects?
  

namelist / esfsw_driver_nml /    &
                               do_ica_calcs, do_rayleigh_all_bands, &
                               do_herzberg, do_quench, &
                               do_ch4_sw_effects, do_n2o_sw_effects, &
                               do_coupled_stratozone

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!    variables associated with absorptivity and sw transmission for
!    various gaseous atmospheric components
!
! powph2o      = the scaling factor used in the fit of the h2o          
!                transmission function                                  
!                                                                       
! p0h2o        = the reference pressure (mb) used in the fit of the     
!                h2o transmission function                              
!                                                                       
! c(n)co2(str) = coefficients for the absorptivity expression for co2   
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! c(n)o2(str)  = coefficients for the absorptivity expression for o2    
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! ""(schrun)   = coefficients for the absorptivity expression for the   
!                Schuman-Runge o2 band (non-scaled only)                
!                                                                       
! kh2o         =  the psuedo-absorption coefficients in cm2/gm for h2o  
!                                                                       
! ko3          = the absorption coefficients in cm2/gm for o3           
!                                                                       
! wtfreq       = the weight associated with each exponential term       
! strterm      = logical flag to indicate whether or not a h2o pseudo-  
!                absorption coefficient is assigned a non-scaled        
!                (true) or pressure-scaled (false) gas amount           
!---------------------------------------------------------------------
real, dimension (:), allocatable    :: c1co2, c1co2str, c1o2, c1o2str, &
                                       c2co2, c2co2str, c2o2, c2o2str, &
                                       c3co2, c3co2str, c3o2, c3o2str, &
                                       c4co2, c4co2str, c4o2, c4o2str, &
                                       c1ch4, c1ch4str, c2ch4,         &
                                       c2ch4str, c3ch4, c3ch4str,      &
                                       c4ch4, c4ch4str,                &
                                       c1n2o, c1n2ostr, c2n2o,         &
                                       c2n2ostr, c3n2o, c3n2ostr,      &
                                       c4n2o, c4n2ostr,                &
                                       powph2o, p0h2o
real                                :: c1o2strschrun, c2o2strschrun, &
                                       c3o2strschrun, c4o2strschrun
real, dimension (:), allocatable    :: kh2o, ko3, wtfreq
logical, dimension(:), allocatable  :: strterm

!---------------------------------------------------------------------
!    quantities associated with solar spectral parameterization
!                                                                       
! firstrayband = the first band number where the contribution by        
!                rayleigh scattering is included in the solar           
!                calculations                                           
!                                                                       
! nirbands     = the number of bands in the near-infrared (used in      
!                assigning the value of the surface albedo for the      
!                near-infrared, and the visible and ultraviolet         
!                regions, separately)                                   
! nfreqpts     = the number of pseudo-monochromatic frequencies         
! solflxband   = the solar flux in each parameterization band           
! solflxbandref = the solar flux in each parameterization band, used for
!                 defining band-averaged optical parameters. If the
!                 solar constant is time-invariant, it is also the solar
!                 flux in each parameterization band (solflxband).
! vis_wvnum    = the wavenumber of visible light (corresponds to
!                wavelength of 0.55 microns) [ cm **(-1) ]
!---------------------------------------------------------------------
real                                :: refquanray, solflxtotal
integer                             :: firstrayband, nirbands
integer, dimension (:), allocatable :: nfreqpts
real,    dimension(:), allocatable  :: solflxband
real,    dimension(:), allocatable  :: solflxbandref
real, dimension(:), allocatable     :: wtstr, cosangstr
real, dimension(4)                  :: wtstr_4 =      &
                                       (/0.347854845, 0.652145155,&
                                         0.347854845, 0.652145155/)

integer :: nbands, tot_wvnums, nfrqpts, nh2obands, nstreams
logical :: nstr4 = .false.
real    :: vis_wvnum = 1.0E+04/0.55
real    :: wvnum_870 = 1.0E+04/0.87
real    :: wvnum_340 = 1.0E+04/0.34
real    :: wvnum_380 = 1.0E+04/0.38
real    :: wvnum_440 = 1.0E+04/0.44
real    :: wvnum_670 = 1.0E+04/0.67
real    :: one_micron_wvnum = 1.0E+04/1.00
real    :: onepsix_micron_wvnum = 1.0E+04/1.61
integer :: onepsix_band_indx

!---------------------------------------------------------------------
!    variables associated with rayleigh scattering
!---------------------------------------------------------------------
real, dimension (:), allocatable    :: betaddensitymol

!----------------------------------------------------------------------
!    variables associated with total optical path of species ? - smf
!----------------------------------------------------------------------
real                            :: toto2strmaxschrun
real, dimension(:), allocatable :: totco2max, totco2strmax, &
                                   toto2max, toto2strmax, &
                                   totch4max, totch4strmax, &
                                   totn2omax, totn2ostrmax

!----------------------------------------------------------------------
!    variables associated with the herzberg effect. wtmo2 is the mol-
!    ecular weight of o2. herzberg_fac is a factor used in the last 
!    shortwave band, so that herzberg_fac*wo2 yields a correction for 
!    the o2 optical depth to account for the herzberg o2 heating. this 
!    is done only when do_herzberg is true.
!----------------------------------------------------------------------
real, parameter   :: wtmo2        = 3.19988E+01  
real, parameter   :: herzberg_fac = 9.9488377E-3

!----------------------------------------------------------------------
!    variables associated with the quenching effect. co2_quenchfac is a
!    multiplication factor that reduces the co2 gas optical depth, and 
!    hence, the solar heating in the upper atmosphere, to account for 
!    "quenching" due to non-LTE processes. co2_quenchfac_height is the 
!    reference height for co2_quenchfac [ meters ].
!----------------------------------------------------------------------
real, dimension(30) :: co2_quenchfac
data co2_quenchfac /1.0,.954,.909,.853,.800,.747,.693,.637,.583, .526,&
                    .467,.416,.368,.325,.285,.253,.229,.206,.186,.170,&
                    .163,.156,.151,.144,.138,.132,.127,.124,.068,.037/

real, dimension(30) :: co2_quenchfac_height
data co2_quenchfac_height /67304.,68310.,69303.,70288.,71267.,72245.,&
                           73221.,74195.,75169.,76141.,77112.,78082.,&
                           79051.,80018.,80985.,81950.,82914.,83876.,&
                           84838.,85798.,86757.,87715.,88672.,89627.,&
                           90582.,91535.,92487.,93438.,94387.,106747./



!---------------------------------------------------------------------
!    miscellaneous variables
!---------------------------------------------------------------------
integer, parameter :: NSOLWG = 1
real, dimension(NSOLWG) :: gausswt
logical        :: module_is_initialized = .false.
logical        :: do_esfsw_band_diagnostics = .false.
integer        :: naerosol_optical, naerosoltypes_used


!---------------------------------------------------------------------
!---------------------------------------------------------------------
 


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
! <SUBROUTINE NAME="esfsw_driver_init">
!  <OVERVIEW>
!   Subroutine that defines the time-independent quantities associated
!   with the incoming shortwave radiation in the multiple-band solar
!   radiation parameterization.
!  </OVERVIEW>
!  <DESCRIPTION>
!   It first reads in the input namelist and then allocates gas absorption
!   coefficient variables. It then reads in the shortwave input namelist
!   file and assigns the gas absorption coefficients. Rayleigh scattering
!   coefficient is also calculated based on the temperature and pressure
!   structure of the atmosphere.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call esfsw_driver_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine esfsw_driver_init
 
!---------------------------------------------------------------------- 
!    esfsw_driver_init is the constructor for esfsw_driver_mod. it
!    defines the time-independent quantities associated with the 
!    incoming shortwave radiation in the multiple-band solar radiation 
!    parameterization.                                    
!---------------------------------------------------------------------- 

!---------------------------------------------------------------------
!  local variables:

      real,    dimension(Solar_spect%nbands)    :: freqnu
      real,    dimension(Solar_spect%nstreams)  :: ptstr 
      integer, dimension(0:Solar_spect%nbands) :: endwvnbands

      integer, dimension(:), allocatable  :: nwvnsolar
      real   , dimension(:), allocatable  :: solint   

      character(len=64)    :: file_name
      real, dimension(4)   :: ptstr_4 = (/-0.861136312,&
                                          -0.339981044, &
                                           0.861136312,  &
                                           0.339981044 /)
      real      :: ptstr_1 = 0.2
      real      :: temprefray  = 288.15
      real      :: pressrefray = 101325.       ! MKS units
      real      :: densmolref  = 2.54743E+19
      real      :: convfac     = 1.0E+18
      real      :: corrfac, gamma, f1, f2, f3, pinteg, &
                   twopiesq, densmolrefsqt3, wavelength,  &
                   freqsq, ristdm1, ri
      integer   :: iounit, nband, nf, ni, nw, nw1, nw2, nintsolar
      integer   :: unit, io, ierr, logunit
      integer   :: i
      integer   :: n
      real      :: input_flag = 1.0e-99

!---------------------------------------------------------------------
!  local variables:
!                                                                       
!      freqnu   
!      ptstr          gaussian points and weights for evaluation of the
!                     diffuse beam.
!      nwvnsolar      the number of wavenumbers in each region where the
!                     solar flux is constant                         
!      solint         the solar flux in watts per meter**2 in each      
!                     wavenumber region where it is constant    
!      endwvnbands    the wavenumber value for the band limits   
!      file_name
!      ptstr_4    
!      ptstr_1     
!      temprefray     reference temperature used in defining rayleigh
!                     optical depth
!      pressrefray    reference pressure used in defining rayleigh
!                     optical depth [ Pa ]
!      densmolref     reference density used in defining rayleigh
!                     optical depth
!      convfac     
!      corrfac
!      gamma
!      f1
!      f2
!      f3
!      pinteg
!      twopiesq
!      densmolrefsqt3
!      wavelength
!      freqsq
!      ristdm1
!      ri
!      iounit
!      nband
!      nf
!      ni
!      nw
!      nw1
!      nw2
!      nintsolar      the number of wavenumber regions where the  
!                     solar flux is constant.   
!      unit
!      io
!      ierr
!      i
!                                                                       
!---------------------------------------------------------------------

!---------------------------------------------------------------------
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
      call esfsw_parameters_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=esfsw_driver_nml, iostat=io)
       ierr = check_nml_error(io,'esfsw_driver_nml')
#else   
       if ( file_exist('input.nml')) then
         unit =  open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
         read  (unit, nml=esfsw_driver_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'esfsw_driver_nml')
         end do
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
       call write_version_number (version, tagname)
       logunit = stdlog()
       if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=esfsw_driver_nml)

!---------------------------------------------------------------------
!    define flag indicating if ICA calculations being done.
!---------------------------------------------------------------------
      Cldrad_control%do_ica_calcs = do_ica_calcs
      Cldrad_control%do_ica_calcs_iz = .true.

!---------------------------------------------------------------------
!    allocate module variables
!---------------------------------------------------------------------
      nbands = Solar_spect%nbands       
      tot_wvnums = Solar_spect%tot_wvnums
      nfrqpts = Solar_spect%nfrqpts
      nstreams = Solar_spect%nstreams
      nh2obands = Solar_spect%nh2obands
      allocate ( betaddensitymol (nbands) )
      allocate ( c1co2   (nh2obands),  &
                 c1co2str(nh2obands),  &
                 c1o2    (nh2obands),  &
                 c1o2str (nh2obands),  &
                 c2co2   (nh2obands),  &
                 c2co2str(nh2obands),  &
                 c2o2    (nh2obands),  &
                 c2o2str (nh2obands),  &
                 c3co2   (nh2obands),  &
                 c3co2str(nh2obands),  &
                 c3o2    (nh2obands),  &
                 c3o2str (nh2obands),  &
                 c4co2   (nh2obands),  &
                 c4co2str(nh2obands),  &
                 c4o2    (nh2obands),  &
                 c4o2str (nh2obands),  &
                 c1ch4   (nh2obands),  &
                 c1ch4str(nh2obands),  &
                 c2ch4   (nh2obands),  &
                 c2ch4str(nh2obands),  &
                 c3ch4   (nh2obands),  &
                 c3ch4str(nh2obands),  &
                 c4ch4   (nh2obands),  &
                 c4ch4str(nh2obands),  &
                 c1n2o   (nh2obands),  &
                 c1n2ostr(nh2obands),  &
                 c2n2o   (nh2obands),  &
                 c2n2ostr(nh2obands),  &
                 c3n2o   (nh2obands),  &
                 c3n2ostr(nh2obands),  &
                 c4n2o   (nh2obands),  &
                 c4n2ostr(nh2obands),  &
                 powph2o (nh2obands),  &
                 p0h2o   (nh2obands)    )
      allocate ( nfreqpts        (nbands) )
      allocate ( solflxband      (nbands) )
      allocate ( solflxbandref   (nbands) )
      allocate ( kh2o            (nfrqpts),   & 
                 ko3             (nfrqpts),   &
                 wtfreq          (nfrqpts),   &
                 strterm         (nfrqpts)   )
      allocate ( wtstr           (nstreams),   & 
                 cosangstr       (nstreams)  )
      allocate ( totco2max    (nh2obands),     &
                 totco2strmax (nh2obands),     &
                 toto2max     (nh2obands),     &
                 toto2strmax  (nh2obands),   &
                 totch4max    (nh2obands),   &
                 totch4strmax (nh2obands),   &
                 totn2omax    (nh2obands),   &
                 totn2ostrmax (nh2obands)    )

      betaddensitymol = 0.0 ; c1co2    = 0.0 ; c1co2str = 0.0
      c1o2     = 0.0 ; c1o2str  = 0.0 ; c2co2    = 0.0
      c2co2str = 0.0 ; c2o2     = 0.0 ; c2o2str  = 0.0
      c3co2    = 0.0 ; c3co2str = 0.0 ; c3o2     = 0.0
      c3o2str  = 0.0 ; c4co2    = 0.0 ; c4co2str = 0.0
      c4o2     = 0.0 ; c4o2str  = 0.0 ; powph2o  = 0.0
      p0h2o    = 0.0 ; nfreqpts        = 0.0 ; solflxband      = 0.0
      solflxbandref   = 0.0 ; kh2o            = 0.0 ; ko3             = 0.0
      wtfreq          = 0.0 ; strterm     = .FALSE. ; wtstr           = 0.0
      cosangstr       = 0.0 ; totco2max     = 0.0 ; totco2strmax  = 0.0
      toto2max      = 0.0 ; toto2strmax   = 0.0
      c1ch4    = 0.0 ; c1ch4str = 0.0; c2ch4    = 0.0
      c2ch4str = 0.0 ; c3ch4    = 0.0 ; c3ch4str = 0.0
      c4ch4    = 0.0 ; c4ch4str = 0.0
      totch4max     = 0.0 ; totch4strmax  = 0.0
      c1n2o    = 0.0 ; c1n2ostr = 0.0; c2n2o    = 0.0
      c2n2ostr = 0.0 ; c3n2o    = 0.0 ; c3n2ostr = 0.0
      c4n2o    = 0.0 ; c4n2ostr = 0.0
      totn2omax     = 0.0 ; totn2ostrmax  = 0.0

!---------------------------------------------------------------------
!    allocate local variables.
!---------------------------------------------------------------------
      if (nstreams == 4) then
        ptstr(:) = ptstr_4(:)
        wtstr(:) = wtstr_4(:)
        nstr4 = .true.
      else if (nstreams == 1) then
        ptstr(1) = ptstr_1
        wtstr(1) = 1.0
        nstr4 = .false.
      endif

!---------------------------------------------------------------------
!    read input file for band positions, solar fluxes and band
!    strengths.
!---------------------------------------------------------------------
      if (nbands == 25 .and. nfrqpts == 72) then 
        file_name = 'INPUT/esf_sw_input_data_n72b25'
      else if (nbands == 18 .and. nfrqpts == 38) then
        file_name = 'INPUT/esf_sw_input_data_n38b18'
      else
        call error_mesg ( 'esfsw_driver_mod', &
          'input file for desired bands and frqs is not available', &
                                                               FATAL)
      endif
      iounit = open_namelist_file (file_name)
      read(iounit,101) ( solflxbandref(nband), nband=1,NBANDS )
      read(iounit,102) ( nfreqpts(nband), nband=1,NBANDS )
      read(iounit,103) ( endwvnbands(nband), nband=1,NBANDS )
      read(iounit,103) FIRSTRAYBAND,NIRBANDS
      read(iounit,104) ( powph2o(nband), nband=1,NH2OBANDS )
      read(iounit,104) ( p0h2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1o2str(nband), nband=1,NH2OBANDS ),  &
                         c1o2strschrun
      read(iounit,105) ( c2o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2o2str(nband), nband=1,NH2OBANDS ), &
                         c2o2strschrun
      read(iounit,105) ( c3o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3o2str(nband), nband=1,NH2OBANDS ),  &
                         c3o2strschrun
      read(iounit,105) ( c1ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3ch4(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3ch4str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1n2ostr(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2n2ostr(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3n2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3n2ostr(nband), nband=1,NH2OBANDS )
      do nf = 1,nfrqpts
        read(iounit,106) wtfreq(nf),kh2o(nf),ko3(nf),strterm (nf)
      end do

      read(iounit,107) nintsolar

      allocate ( nwvnsolar (nintsolar) )
      allocate ( solint    (nintsolar) )
 
      do ni = 1,nintsolar
        read(iounit,107) nwvnsolar (ni),solint(ni)
      end do
 
      call close_file (iounit)
 
      if (Solar_spect%tot_wvnums /=       &
                                 endwvnbands(Solar_spect%nbands)) then
        call error_mesg ( 'esfsw_driver_mod', &
         ' inconsistency between highest solar spectrum wavenumber '//&
          'in esfsw_parameters_mod and in esfsw_sriver input file', &
                                                           FATAL)
      endif

!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.34 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_340) then
          Solar_spect%w340_band_indx = ni
          Solar_spect%w340_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to near ultraviolet light
!    (0.38 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_380) then
          Solar_spect%w380_band_indx = ni
          Solar_spect%w380_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to blue light
!    (0.44 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_440) then
          Solar_spect%w440_band_indx = ni
          Solar_spect%w440_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to red light
!    (0.67 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_670) then
          Solar_spect%w670_band_indx = ni
          Solar_spect%w670_band_iz = .true.
          exit
        endif
      end do
!---------------------------------------------------------------------
!    define the band index corresponding to visible light
!    (0.55 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > vis_wvnum) then
          Solar_spect%visible_band_indx = ni
          Solar_spect%visible_band_indx_iz = .true.
          exit
        endif
      end do

!---------------------------------------------------------------------
!    define the band index corresponding to 870nm
!    (0.87 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > wvnum_870) then
          Solar_spect%eight70_band_indx = ni
          Solar_spect%eight70_band_indx_iz = .true.
          exit
        endif
      end do

!---------------------------------------------------------------------
!    define the band index corresponding to near infra red band
!    (1.00 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > one_micron_wvnum) then
          Solar_spect%one_micron_indx = ni
          Solar_spect%one_micron_indx_iz = .true.
          exit
        endif
      end do
 
!---------------------------------------------------------------------
!    define the band index corresponding to               
!    (1.61 microns). note that endwvnbands is in units of (cm**-1).
!---------------------------------------------------------------------
      do ni=1,nbands
        if (endwvnbands(ni) > onepsix_micron_wvnum) then
          onepsix_band_indx = ni
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the wavenumber one solar fluxes.                       
!----------------------------------------------------------------- --
      do ni = 1,nintsolar
        if ( ni.eq.1 ) then
          nw1 = 1
        else
          nw1 = nw1 + nwvnsolar(ni-1)
        end if
        nw2 = nw1 + nwvnsolar(ni) - 1
        do nw = nw1,nw2
          Solar_spect%solarfluxtoa(nw) = solint(ni)
        end do
      end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      deallocate  (solint    )
      deallocate  (nwvnsolar )
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if ( .not. Rad_control%using_solar_timeseries_data) then
        Solar_spect%solflxband = solflxbandref
      endif
      Solar_spect%solflxbandref = solflxbandref
      Solar_spect%endwvnbands = endwvnbands

!--------------------------------------------------------------------
!    override the input file value of firstrayband if the nml control
!    variables indicates rayleigh effects are to be considered in
!    all bands.
!--------------------------------------------------------------------
      if (do_rayleigh_all_bands)  firstrayband = 1

!----------------------------------------------------------------------
!    convert some values to mks to match model units
!--------------------------------------------------------------------
      p0h2o = 1.0E-2/p0h2o  ! invert, and convert mb to mks
      kh2o = kh2o *1.0E-01   ! cgs to mks
      ko3  = ko3 *1.0E-01    ! cgs to mks

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      do n=1,NH2OBANDS
        if (c1co2(n) /= input_flag .and.  &
            c2co2(n) /= input_flag .and.  &
            c3co2(n) /= input_flag ) then
          c4co2(n) = c1co2(n) * c2co2(n) ** c3co2(n)
          c4co2str(n) = c1co2str(n) * c2co2str(n) ** c3co2str(n)
          totco2max(n) = ( (1.0/c1co2(n) ) + c2co2(n) ** c3co2(n) ) ** &
                       (1.0/c3co2(n) ) - c2co2(n)
          if (nbands == 18) then
            if ( n /= 4) then
              totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
            else 
              totco2strmax(n) = HUGE (c4o2strschrun) 
            endif
          else
              totco2strmax(n) = ( (1.0/c1co2str(n) ) + c2co2str(n) ** &
                                 c3co2str(n) ) ** (1.0/c3co2str(n) ) - &
                                c2co2str(n)
          endif
        else
          c4co2(n) = 0.0                              
          c4co2str(n) = 0.0
          totco2max(n) = 0.0                                            
          totco2strmax(n) = 0.0
        endif
        if (c1o2(n) /= input_flag .and.   &
            c2o2(n) /= input_flag .and.   &
            c3o2(n) /= input_flag ) then
          c4o2(n) = c1o2(n) * c2o2(n) ** c3o2(n)
          c4o2str(n) = c1o2str(n) * c2o2str(n) ** c3o2str(n)
          toto2max(n) = ( (1.0/c1o2(n) ) + c2o2(n) ** c3o2(n) ) ** &
                          (1.0/c3o2(n) ) - c2o2(n)
          if (nbands == 18) then
            if ( n /= 4) then
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                c2o2str(n)
            else
              toto2strmax(n) = HUGE (c4o2strschrun) 
            endif
          else
              toto2strmax(n) = ( (1.0/c1o2str(n) ) + c2o2str(n) ** &
                                c3o2str(n) ) ** (1.0/c3o2str(n) ) - &
                                c2o2str(n)
          endif
        else
          c4o2(n) = 0.0                              
          c4o2str(n) = 0.0
          toto2max(n) = 0.0                                            
          toto2strmax(n) = 0.0
        endif
    if (do_ch4_sw_effects) then
        if (c1ch4(n) /= input_flag .and.  &
            c2ch4(n) /= input_flag .and.  &
            c3ch4(n) /= input_flag ) then
          c4ch4(n) = c1ch4(n) * c2ch4(n) ** c3ch4(n)
          c4ch4str(n) = c1ch4str(n) * c2ch4str(n) ** c3ch4str(n)
          totch4max(n) = ( (1.0/c1ch4(n) ) + c2ch4(n) ** c3ch4(n) ) ** &
                           (1.0/c3ch4(n) ) - c2ch4(n)
          totch4strmax(n) = ( (1.0/c1ch4str(n) ) + c2ch4str(n) ** &
                            c3ch4str(n) ) ** (1.0/c3ch4str(n) ) - &
                            c2ch4str(n)
        else
          c4ch4(n) = 0.
          c4ch4str(n) = 0.
          totch4max(n) = 0.
          totch4strmax(n) = 0.     
        endif
     endif
     if (do_n2o_sw_effects) then
        if (c1n2o(n) /= input_flag .and.  &
            c2n2o(n) /= input_flag .and.  &
            c3n2o(n) /= input_flag ) then
          c4n2o(n) = c1n2o(n) * c2n2o(n) ** c3n2o(n)
          c4n2ostr(n) = c1n2ostr(n) * c2n2ostr(n) ** c3n2ostr(n)
          totn2omax(n) = ( (1.0/c1n2o(n) ) + c2n2o(n) ** c3n2o(n) ) ** &
                           (1.0/c3n2o(n) ) - c2n2o(n)
          totn2ostrmax(n) = ( (1.0/c1n2ostr(n) ) + c2n2ostr(n) ** &
                               c3n2ostr(n) ) ** (1.0/c3n2ostr(n) ) - &
                               c2n2ostr(n)
        else
          c4n2o(n) = 0.
          c4n2ostr(n) = 0.
          totn2omax(n) = 0.
          totn2ostrmax(n) = 0.     
        endif
     endif
      end do

      c4o2strschrun = c1o2strschrun * c2o2strschrun ** c3o2strschrun
      toto2strmaxschrun = ( (1.0/c1o2strschrun) + c2o2strschrun ** &
                             c3o2strschrun) ** (1.0/c3o2strschrun) - &
                             c2o2strschrun

!     if (mpp_pe() == 0) then
!       print *, 'c1ch4    ', c1ch4    
!       print *, 'c2ch4    ', c2ch4    
!       print *, 'c3ch4    ', c3ch4    
!       print *, 'c4ch4    ', c4ch4    
!       print *, 'totch4max', totch4max
!       print *, 'c1ch4str ', c1ch4str    
!       print *, 'c2ch4str ', c2ch4str    
!       print *, 'c3ch4str ', c3ch4str    
!       print *, 'c4ch4str ', c4ch4str    
!       print *, 'totch4strmax', totch4strmax
!       print *, 'c1co2    ', c1co2    
!       print *, 'c2ch4    ', c2co2    
!       print *, 'c3ch4    ', c3co2    
!       print *, 'c4ch4    ', c4co2    
!       print *, 'totch4max', totco2max
!       print *, 'c1ch4str ', c1co2str    
!       print *, 'c2ch4str ', c2co2str    
!       print *, 'c3ch4str ', c3co2str    
!       print *, 'c4ch4str ', c4co2str    
!       print *, 'totch4strmax', totco2strmax
!       print *, 'c1n2o    ', c1n2o    
!       print *, 'c2ch4    ', c2n2o    
!       print *, 'c3ch4    ', c3n2o    
!       print *, 'c4ch4    ', c4n2o    
!       print *, 'totch4max', totn2omax
!       print *, 'c1ch4str ', c1n2ostr    
!       print *, 'c2ch4str ', c2n2ostr    
!       print *, 'c3ch4str ', c3n2ostr    
!       print *, 'c4ch4str ', c4n2ostr    
!       print *, 'totch4strmax', totn2ostrmax
!       print *, '1/c1', 1.0/c1co2str
!       print *, '1/c3', 1.0/c3co2str
!       print *, 'c2**c3', c2co2str**c3co2str
!       print *, '(1 / + c2**c3)**1/c3',  &
!                 (1.0/c1co2str + c2co2str**c3co2str)**(1./c3co2str)
!       print *, 'c4o2    ', c4o2    
!       print *, 'c4o2str    ', c4o2str    
!       print *, 'toto2max', toto2max
!       print *, 'toto2strmax', toto2strmax
!     endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if ( .not. Rad_control%using_solar_timeseries_data) then
      solflxtotal = 0.0
      do nband = 1,NBANDS
        solflxtotal = solflxtotal + Solar_spect%solflxbandref(nband)
      end do
      endif
 
!----------------------------------------------------------------------
!    define the wavenumbers to evaluate rayleigh optical depth.      
!-------------------------------------------------------------------
      endwvnbands(0) = 0
      do nband = 1,NBANDS
        freqnu(nband) = 0.5 * ( endwvnbands(nband-1) +   &
                                endwvnbands(nband) )
      end do
 
!---------------------------------------------------------------------
!    define quantities used to determine the rayleigh optical depth. 
!    notes: refquanray is the quantity which multiplies pressure /  
!           temperature to yield the molecular density.                 
!           betaddensitymol is the quantity which multiples the       
!           molecular density to yield the rayleigh scattering      
!           coefficient.                                           
!           1.39E-02 is the depolorization factor.              
!-----------------------------------------------------------------
      refquanray = densmolref * temprefray / pressrefray 
      corrfac = ( 6.0E+00 + 3.0E+00 * 1.39E-02 )/( 6.0E+00 - 7.0E+00 * &
                  1.39E-02 )
      gamma = 1.39E-02 / ( 2.0E+00 - 1.39E-02 )
      f1 = 7.5E-01 / ( gamma * 2.0E+00 + 1.0E+00 )
      f2 = gamma * 3.0E+00 + 1.0E+00 
      f3 = 1.0E+00 - gamma 
      pinteg = 2.0E+00 * PI * ( 2.0E+00 * f1 * f2 * ( 1.0E+00 + f3 / &
               f2 / 3.0E+00 ) )
      twopiesq = 2.0E+00 *  PI ** 2
      densmolrefsqt3 = 3.0E+00 * densmolref ** 2
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do nband = 1,NBANDS
        wavelength = 1.0E+04 / freqnu(nband)
        freqsq = 1.0 / ( wavelength ) ** 2
        ristdm1 = ( 6.4328E+03 + 2.94981E+06 / ( 1.46E+02 - freqsq ) + &
                    2.554E+04 / ( 4.1E+01 - freqsq ) ) * 1.0E-08
        ri = ristdm1 + 1
        betaddensitymol(nband) = twopiesq*( ri ** 2 - 1.0E+00 ) ** 2 * &
                                 corrfac / ( densmolrefsqt3 *  &
                                 wavelength ** 4 ) * pinteg * convfac 
      end do
 
      gausswt(1) = 1.0

!---------------------------------------------------------------------
!    define the gaussian angles for evaluation of the diffuse beam.  
!--------------------------------------------------------------------
      do i = 1,nstreams
        cosangstr(i) = ( ptstr(i) + 1. ) * 5.0E-01
      end do

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------
 101  format( 12f10.4 )
 102  format( 32i4 )
 103  format( 20i6 )
 104  format( 12f10.2 )
 105  format( 1p,16e8.1 )
 106  format( 1p,3e16.6,l16 )
 107  format( i5,1p,e14.5 )
 
!---------------------------------------------------------------------


end subroutine esfsw_driver_init 



!#################################################################
! <SUBROUTINE NAME="swresf">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call swresf(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases,  &
                   Aerosol, Aerosol_props, Astro, Cldrad_props,  &
                   Cld_spec, including_volcanoes, Sw_output,   &
                   Aerosol_diags, r, including_aerosols,   &
                   naerosol_optical) 

!----------------------------------------------------------------------
!    swresf uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

integer,                       intent(in)    :: is, ie, js, je
type(atmos_input_type),        intent(in)    :: Atmos_input
type(surface_type),            intent(in)    :: Surface
type(radiative_gases_type),    intent(in)    :: Rad_gases   
type(aerosol_type),            intent(in)    :: Aerosol     
type(aerosol_properties_type), intent(in)    :: Aerosol_props
type(astronomy_type),          intent(in)    :: Astro
type(cldrad_properties_type),  intent(in)    :: Cldrad_props
type(cld_specification_type),  intent(in)    :: Cld_spec      
logical,                       intent(in)    :: including_volcanoes
type(sw_output_type),          intent(inout) :: Sw_output   
type(aerosol_diagnostics_type),intent(inout) :: Aerosol_diags
real,dimension(:,:,:,:),       intent(inout) :: r
logical,                       intent(in)    :: including_aerosols  
integer,                       intent(in)    :: naerosol_optical


!-------------------------------------------------------------------
!  intent(in) variables:
!
!    is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Surface        surface_type structure, contains variables 
!                     defining the surface characteristics
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Aerosol_props  aerosol radiative property input data for the 
!                     radiation package
!      Astro          astronomy_type structure
!      Cldrad_props   cloud radiative property input data for the 
!                     radiation package
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 

      logical, dimension (size(Atmos_input%temp,1), &
                          size(Atmos_input%temp,2), &
                          size(Atmos_input%temp,3)-1)  :: &
                             cloud

      logical, dimension (size(Atmos_input%temp,1), &
                          size(Atmos_input%temp,2))  ::   &
                              cloud_in_column,   daylight

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1, &
                        Solar_spect%nbands)  :: &
            aeroasymfac,     aerosctopdep,   aeroextopdep, &
            rayopdep

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1, &
                       nfrqpts,NSOLWG)  ::     &
                                                    gasopdep

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1)  :: &
            cloudfrac,      cldfrac_band,    cldfrac,         &
            deltaz,         cloud_deltaz,                     &
            cldext,         cldsct,          cldasymm,         &
            cloudasymfac,   cloudextopdep,                     &
            cloudsctopdep,                   deltap,           &
            densitymol,     extopdepclr,                       &
            extopdepovc,    fclr,            fovc,             &
            gocpdp,         gclr,            gstrclr,          &
            gstrovc,        govc,            omegastrclr,      &
            omegastrovc,    rlayerdif,       rlayerdir,        &
            rlayerdifclr,   rlayerdifovc,    &
            rlayerdirclr,   rlayerdirovc,                     &
            sctopdepclr,    sctopdepovc,      &
            ssalbclr,        ssalbovc,                         &
            taustrclr,      taustrovc,        &
            tlayerde,       tlayerdeclr,      &
            tlayerdeovc,    tlayerdif,                        &
            tlayerdifclr,   tlayerdifovc,     &
            tlayerdir,      tlayerdirclr,     &
            tlayerdirovc
           
      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3))  :: &
              reflectance,  transmittance,  tr_dir            

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3), &
                       Rad_control%nzens)  :: &
            dfswbandclr,     fswbandclr,    ufswbandclr, &
            dfswband,        fswband,       ufswband,  &
            sumtrclr,        sumreclr,      sumtr_dir_clr,  & 
            sumtr,           sumre,         sumtr_dir

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       Rad_control%nzens)  :: sumtr_dir_up

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3))  :: &
            press,           pflux,            pflux_mks, &
            temp,                                       &
            reflectanceclr,  transmittanceclr, tr_dirclr

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)-1 , &
                       Rad_control%nzens)  :: &
                hswbandclr,  hswband

      real, dimension (size(Atmos_input%temp,1), &
                       size(Atmos_input%temp,2))    ::  &
            sfcalb_dir,   sfcalb_dif,  wtfac_p,    &
            fracday_p,    solarflux_p

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),   &
                                         NSOLWG) ::         &
                 cosangsolar_p

      integer :: j, i, k, ng, np, nband, nf, ns, nz
      integer :: nzens

      integer :: nprofile, nprofiles
      real    :: profiles_inverse

      integer :: ix, jx, kx, israd, jsrad, ierad, jerad, ksrad, kerad
      real    :: ssolar  
      real    :: solflxtotal_local




!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_driver_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the solar_constant appropriate at Rad_time when solar
!    input is varying with time.
!---------------------------------------------------------------------
!-- The following fix is for openmp
      if (Rad_control%using_solar_timeseries_data) then
        solflxtotal_local = 0.0
        do nband = 1,NBANDS
          solflxtotal_local = solflxtotal_local + Solar_spect%solflxband(nband)
        end do
      else
        solflxtotal_local = solflxtotal
      endif

!---------------------------------------------------------------------
!
      
      cldfrac = Cld_spec%camtsw

!---------------------------------------------------------------------
!  convert to cgs and then back to mks for consistency with previous 
!---------------------------------------------------------------------
      press(:,:,:) = 0.1*(10.0*Atmos_input%press(:,:,:))
      pflux(:,:,:) =     (10.0*Atmos_input%pflux(:,:,:))
      deltaz(:,:,:) = Atmos_input%deltaz(:,:,:)
      temp  (:,:,:) = Atmos_input%temp  (:,:,:)
      cloud_deltaz = Atmos_input%clouddeltaz

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      ix = size(temp,1)
      jx = size(temp,2)
      kx = size(temp,3) - 1
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = ix
      jerad = jx
      kerad = kx

!----------------------------------------------------------------------c
!    define a flag indicating columns in which there is sunshine during
!    this radiation time step. define a flag indicating points with both
!    sunlight and cloud.      
!----------------------------------------------------------------------c
      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if ( Astro%fracday(i,j) /= 0.0) then
            daylight(i,j) = .true.                 
          else
            daylight(i,j) = .false.                
          endif     
          cloud_in_column(i,j) = .false.
        end do
      end do
        
      call compute_aerosol_optical_props     &
                (Atmos_input, Aerosol, Aerosol_props, &
                 including_volcanoes, Aerosol_diags, r,  &
                 including_aerosols, naerosol_optical, &
                 daylight, aeroextopdep, aerosctopdep, aeroasymfac) 

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ssolar = Sw_control%solar_constant*Astro%rrsun

 
!----------------------------------------------------------------------
!    define a flag indicating points with both sunlight and cloud. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
      if (.not. Cldrad_control%do_stochastic_clouds) then
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (daylight(i,j)) then
              do k=KSRAD,KERAD
                if (cldfrac(i,j,k) > 0.0)  then
                  cloud_in_column(i,j) = .true.
                  cloud(i,j,k) = .true.
                  cloudfrac(i,j,k) = cldfrac(i,j,k)
                else
                  cloud(i,j,k) = .false.
                  cloudfrac(i,j,k) = 0.0
                endif
              end do
            else
              do k=KSRAD,KERAD
                cloud(i,j,k) = .false.
                cloudfrac(i,j,k) = 0.0
              end do
            endif
          end do
        end do
      endif

!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
      pflux_mks = pflux*1.0E-1

      do k = KSRAD+1,KERAD+1
        deltap(:,:,k-1) = pflux_mks(:,:,k) - pflux_mks(:,:,k-1)
        gocpdp(:,:,k-1) = radcon_mks/deltap(:,:,k-1)
      end do
 
      call compute_gas_props (Atmos_input, Rad_gases, Astro,   &
                              daylight, gasopdep)
  
!---------------------------------------------------------------------
!    define the molecular density for use in calculating the           
!    rayleigh optical depth (deltaz is in meters).                     
!--------------------------------------------------------------------
      do k = KSRAD,KERAD
        densitymol(:,:,k) = refquanray * press(:,:,k) / temp(:,:,k)
      end do
 

! assumption is that there is 1 cloud profile for each sw band
      if (do_ica_calcs) then
        nprofiles = nbands
        profiles_inverse = 1.0/nprofiles
      else
        nprofiles = 1
        profiles_inverse = 1.0
      endif

!--------------------------------------------------------------------
!    define the rayleigh optical depths.                                
!---------------------------------------------------------------------
      do nband = 1, NBANDS
        rayopdep(:,:,:,nband) = betaddensitymol(nband)*  &
                                        densitymol(:,:,:)*deltaz(:,:,:)
      end do   ! (nband loop)

      nzens = Rad_control%nzens

!--------------------------------------------------------------------
 
      do nprofile=1, nprofiles
        if (do_ica_calcs) then
          cldfrac_band(:,:,:) = Cld_spec%camtsw_band(:,:,:,nprofile)
        endif

!----------------------------------------------------------------------c
!    np is a counter for the pseudo-monochromatic frequency point 
!    number.
!----------------------------------------------------------------------c
        np = 0
 
        do nz=1,nzens
          if (Rad_control%do_totcld_forcing) then
            dfswbandclr(:,:,:,nz) = 0.0
            fswbandclr(:,:,:,nz) = 0.0
            hswbandclr(:,:,:,nz) = 0.0
            ufswbandclr(:,:,:,nz) = 0.0
          endif
        end do
        reflectanceclr = 0.0
        transmittanceclr = 0.0

!----------------------------------------------------------------------c
!    begin band loop                                                   
!----------------------------------------------------------------------c
        do nband = 1,NBANDS
 
          do nz=1,nzens
            sumtr(:,:,:,nz) = 0.0
            sumtr_dir(:,:,:,nz) = 0.0
            sumtr_dir_up(:,:,nz) = 0.0
            sumre(:,:,:,nz) = 0.0
            if (Rad_control%do_totcld_forcing) then
              sumtrclr(:,:,:,nz) = 0.0
              sumreclr(:,:,:,nz) = 0.0
              sumtr_dir_clr(:,:,:,nz) = 0.0
            endif
          end do
          if (Cldrad_control%do_stochastic_clouds) then
            if (.not. do_ica_calcs) then
              cldfrac_band(:,:,:) = Cld_spec%camtsw_band(:,:,:,nband)
            endif
          endif
 
!----------------------------------------------------------------------
!    if stochastic clouds are activated (cloud fields differ with sw
!    parameterization band), define a flag indicating points with both
!    sunlight and cloud for the current sw parameterization band. set
!    the flag indicating that cloud is present in columns with such
!    points.
!----------------------------------------------------------------------
          if (Cldrad_control%do_stochastic_clouds) then
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                cloud_in_column(i,j) = .false.
                if (daylight(i,j)) then
                  do k = KSRAD,KERAD
                    if (cldfrac_band(i,j,k) > 0.0)  then
                      cloud_in_column(i,j) = .true.
                      cloud(i,j,k) = .true.
                      cloudfrac(i,j,k) = cldfrac_band(i,j,k)
                    else
                      cloud(i,j,k) = .false.
                      cloudfrac(i,j,k) = 0.0
                    endif
                  end do
                else
                  do k = KSRAD,KERAD
                    cloud(i,j,k) = .false.
                    cloudfrac(i,j,k) = 0.0
                  end do
                endif
              end do
            end do
          endif


!---------------------------------------------------------------------
!    obtain cloud properties from the Cldrad_props input variable.
!--------------------------------------------------------------------
          cldext(:,:,:) = Cldrad_props%cldext(:,:,:,nband,nprofile)
          cldsct(:,:,:) = Cldrad_props%cldsct(:,:,:,nband,nprofile)
          cldasymm(:,:,:) = Cldrad_props%cldasymm(:,:,:,nband,nprofile)

          do k = KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD, IERAD
                if (cloud(i,j,k) ) then
                  cloudextopdep(i,j,k) = 1.0E-03*cldext(i,j,k) *    &
                                         cloud_deltaz(i,j,k)
                  cloudsctopdep(i,j,k) = 1.0E-03*cldsct(i,j,k) *    &
                                         cloud_deltaz(i,j,k)
                  cloudasymfac(i,j,k) = cldasymm(i,j,k)
                endif
              end do
            end do
          end do

!-----------------------------------------------------------------
!    define clear sky arrays
!-----------------------------------------------------------------
          if (nband >= firstrayband ) then
            do k=ksrad,kerad
              do j=jsrad,jerad
                do i=israd,ierad
                  if (daylight(i,j) ) then
                    sctopdepclr(i,j,k) = rayopdep(i,j,k,nband) +   &
                                            aerosctopdep(i,j,k,nband)
                    gclr(i,j,k) = aeroasymfac(i,j,k,nband)*  &
                                  aerosctopdep(i,j,k,nband)/&
                                                     sctopdepclr(i,j,k)
                    fclr(i,j,k) = aeroasymfac(i,j,k,nband)*gclr(i,j,k)
                    gstrclr(i,j,k) = ( gclr(i,j,k)  - fclr(i,j,k) )/  &
                                     ( 1.0 - fclr(i,j,k) )
                  endif
                end do
              end do
            end do
          endif

!-----------------------------------------------------------------
!    define cloudy sky arrays
!-----------------------------------------------------------------
          do k=KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (cloud(i,j,k)) then
                  sctopdepovc(i,j,k) = rayopdep(i,j,k,nband) +    &
                                       aerosctopdep(i,j,k,nband) + &
                                       cloudsctopdep(i,j,k) 
                  govc(i,j,k) = ( ( cloudasymfac(i,j,k) *   &
                                    cloudsctopdep(i,j,k) ) +  &
                                  ( aeroasymfac(i,j,k,nband) *   &
                                    aerosctopdep(i,j,k,nband)))/   &
                                    sctopdepovc(i,j,k)
                  fovc(i,j,k) = ( ( cloudasymfac(i,j,k) ** 2 *  &
                                    cloudsctopdep(i,j,k) ) + &
                                  ( aeroasymfac(i,j,k,nband) ** 2 *  &
                                    aerosctopdep(i,j,k,nband) ))/  &
                                    sctopdepovc(i,j,k)
                  gstrovc(i,j,k) = ( govc(i,j,k)  - fovc(i,j,k))/  &
                                   ( 1.0 - fovc(i,j,k) )
                endif
              end do
            end do
          end do

!---------------------------------------------------------------------
!    begin frequency points in the band loop.                          
!--------------------------------------------------------------------
          do nf = 1,nfreqpts(nband)
            np = np + 1
 
!---------------------------------------------------------------------
!    begin gaussian angle loop (ng > 1 only when lswg = true).        
!--------------------------------------------------------------------
            do ng = 1,NSOLWG
 
!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------
              if (nband >= firstrayband )  then
                do k=ksrad,kerad
                  do j=jsrad,jerad
                    do i=israd,ierad
                      if (daylight(i,j) ) then
                        extopdepclr(i,j,k) = gasopdep(i,j,k,np,ng) +  &
                                             rayopdep(i,j,k,nband) +   &
                                             aeroextopdep(i,j,k,nband)
                        ssalbclr(i,j,k) = sctopdepclr(i,j,k)/    &
                                          extopdepclr(i,j,k)
                        taustrclr(i,j,k) = extopdepclr(i,j,k)*( 1.0 -  &
                                           ssalbclr(i,j,k)*fclr(i,j,k) )
                        omegastrclr(i,j,k) =      &
                               ssalbclr(i,j,k)*((1.0 - fclr(i,j,k))/  &
                               (1.0 -  ssalbclr(i,j,k)*fclr(i,j,k)))
                      endif
                    end do
                  end do
                end do
              endif

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                         
!--------------------------------------------------------------------
              do k=KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (cloud(i,j,k) ) then
                      extopdepovc(i,j,k) = gasopdep(i,j,k,np,ng) +    &
                                           rayopdep(i,j,k,nband) +  &
                                           aeroextopdep(i,j,k,nband) + &
                                           cloudextopdep(i,j,k)
                      ssalbovc(i,j,k) = sctopdepovc(i,j,k) /    &
                                        extopdepovc(i,j,k)
                      taustrovc(i,j,k) = extopdepovc(i,j,k)*( 1.0 - &
                                         ssalbovc(i,j,k)*fovc(i,j,k) )
                      omegastrovc(i,j,k) = ssalbovc(i,j,k)*( ( 1.0 - &
                                           fovc(i,j,k) )/( 1.0 -   &
                                           ssalbovc(i,j,k) *   &
                                           fovc(i,j,k) ) )
                    endif
                  end do
                end do
              end do

!---------------------------------------------------------------------
!    do calculations for all desired zenith angles.
!---------------------------------------------------------------------
              do nz = 1, nzens
                if (Rad_control%hires_coszen) then
                  cosangsolar_p(:,:,ng) = Astro%cosz_p(:,:,nz)   
                else
                  cosangsolar_p(:,:,ng) = Astro%cosz(:,:)          
                endif

                where (cosangsolar_p(:,:,:) == 0.0)   &
                                            cosangsolar_p(:,:,:) = 1.0

!---------------------------------------------------------------------
!    clear sky mode                                                    
!    note: in this mode, the delta-eddington method is performed for all
!    spatial columns experiencing sunlight.         
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    calculate the scaled single-scattering quantities for use in the   
!    delta-eddington routine.                                       
!---------------------------------------------------------------------
                if (nband >= firstrayband )  then

!---------------------------------------------------------------------
!    do diffuse calculation only for first zenith angle -- it is 
!    independent of zenith angle.
!---------------------------------------------------------------------
                  if (nz == 1) then
                    call deledd     &
                        (ix, jx, kx, taustrclr, omegastrclr, &
                         gstrclr, cosangsolar_p(:,:,ng), ng, daylight, &
                         rlayerdirclr, tlayerdirclr, tlayerdeclr, &
                         rlayerdif=rlayerdifclr, tlayerdif=tlayerdifclr)
                  else
                    call deledd     &
                        (ix, jx, kx, taustrclr, omegastrclr, &
                         gstrclr, cosangsolar_p(:,:,ng), ng, daylight,&
                         rlayerdirclr, tlayerdirclr, tlayerdeclr)
                  endif

!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1.
!---------------------------------------------------------------------
                  if (ng /= 1) then
                    tlayerdifclr = 0.0       
                    if (NSTREAMS == 1) then
                      tlayerdifclr(:,:,:) =   &
                             exp(-gasopdep(:,:,:,np,ng)/cosangstr(1))
                    else
                      do ns = 1,NSTREAMS
                        tlayerdifclr(:,:,:) =    &
                             tlayerdifclr(:,:,:) +  & 
                                   exp( -gasopdep(:,:,:,np,ng)/&
                                         cosangstr(ns) )*wtstr(ns)* &
                                         cosangstr(ns)
                      end do
                    endif
                    rlayerdifclr = 0.0
                  endif
  
!---------------------------------------------------------------------
!    initialize the layer reflection and transmission arrays for the   
!    non-rayleigh-scattering case.                            
!-------------------------------------------------------------------
                else
!---------------------------------------------------------------------
!    the following needs to be done at daylight points only -- currently
!    this code is not active, since ng == 1, and all bands see rayleigh
!    scattering.
!---------------------------------------------------------------------
                  tlayerdifclr = 0.0       
                  if (NSTREAMS == 1) then
                    tlayerdifclr(:,:,:) =     &
                            exp( -gasopdep(:,:,:,np,ng)/cosangstr(1))
                  else
                    do ns = 1,NSTREAMS
                      tlayerdifclr(:,:,:) =   &
                         tlayerdifclr(:,:,:) +    &
                              exp(-gasopdep(:,:,:,np,ng)/&
                                cosangstr(ns))*wtstr(ns)*cosangstr(ns)
                    end do
                  endif
                  rlayerdirclr(:,:,:) = 0.0
                  do k=KSRAD,KERAD
                    tlayerdirclr(:,:,k) =      &
                                  exp( -gasopdep(:,:,k,np,ng) /   &
                                                 cosangsolar_p(:,:,ng) )
                  end do  
                  tlayerdeclr(:,:,:) = tlayerdirclr(:,:,:)
                  rlayerdifclr = 0.0
                endif

!---------------------------------------------------------------------
!    overcast sky mode                                                  
!    note: in this mode, the delta-eddington method is performed only 
!    for spatial columns containing a cloud and experiencing sunlight. 
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!    calculate the reflection and transmission in the scattering layers 
!    using the delta-eddington method.                                  
!-------------------------------------------------------------------
                if (nz == 1) then
                  call deledd      &
                       (ix, jx, kx, taustrovc, omegastrovc, gstrovc, &
                        cosangsolar_p(:,:,ng), ng, daylight, &
                        rlayerdirovc, tlayerdirovc, tlayerdeovc, &
                        rlayerdif=rlayerdifovc, tlayerdif=tlayerdifovc,&
                        cloud=cloud)
                else
                  call deledd      &
                       (ix, jx, kx, taustrovc, omegastrovc, gstrovc, &
                        cosangsolar_p(:,:,ng), ng, daylight, &
                        rlayerdirovc, tlayerdirovc, tlayerdeovc,   & 
                        cloud=cloud)
                endif
                if (ng /= 1) then
                  tlayerdifovc(:,:,:) = tlayerdifclr(:,:,:)
                  rlayerdifovc(:,:,:) = rlayerdifclr(:,:,:)
                endif
 
!-------------------------------------------------------------------- 
!    weight the reflection and transmission arrays for clear and        
!    overcast sky conditions by the cloud fraction, to calculate the    
!    resultant values.                                                  
!---------------------------------------------------------------------- 
                do k=KSRAD,KERAD
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if ( cloud(i,j,k) ) then
                        rlayerdir(i,j,k) = cloudfrac(i,j,k)*   &
                                           rlayerdirovc(i,j,k) +  &
                                           (1.0 - cloudfrac(i,j,k) )*  &
                                           rlayerdirclr(i,j,k)
                        rlayerdif(i,j,k) = cloudfrac(i,j,k) *  &
                                           rlayerdifovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           rlayerdifclr(i,j,k)
                        tlayerdir(i,j,k) = cloudfrac(i,j,k) *   &
                                           tlayerdirovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdirclr(i,j,k)
                        tlayerdif(i,j,k) = cloudfrac(i,j,k) *   &
                                           tlayerdifovc(i,j,k) +  &
                                           ( 1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdifclr(i,j,k)
                        tlayerde(i,j,k) =  cloudfrac(i,j,k) *   &
                                           tlayerdeovc(i,j,k) +  &
                                           (1.0 - cloudfrac(i,j,k) )* &
                                           tlayerdeclr(i,j,k)
                      else if (daylight(i,j)) then
                        rlayerdir(i,j,k) = rlayerdirclr(i,j,k)
                        tlayerdir(i,j,k) = tlayerdirclr(i,j,k)
                        rlayerdif(i,j,k) = rlayerdifclr(i,j,k)
                        tlayerdif(i,j,k) = tlayerdifclr(i,j,k)
                        tlayerde (i,j,k) = tlayerdeclr (i,j,k)
                      endif
                    end do
                  end do
                end do
 
!---------------------------------------------------------------------
!    define the surface albedo (infrared value for infrared bands,      
!    visible value for the remaining bands).                            
!----------------------------------------------------------------------c
                if (nband <= NIRBANDS ) then
                  sfcalb_dir(:,:) = Surface%asfc_nir_dir(:,:)
                  sfcalb_dif(:,:) = Surface%asfc_nir_dif(:,:)
                else
                  sfcalb_dir(:,:) = Surface%asfc_vis_dir(:,:)
                  sfcalb_dif(:,:) = Surface%asfc_vis_dif(:,:)
                end if
 
!-------------------------------------------------------------------- 
!    calculate the reflection and transmission at flux levels from the  
!    direct and diffuse values of reflection and transmission in the  
!    corresponding layers using the adding method.                      
!---------------------------------------------------------------------
                call adding         &
                    (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif, &
                     tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,    &
                     daylight, reflectance, transmittance, tr_dir)    

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                if (Rad_control%do_totcld_forcing) then
                  call adding       &
                      (ix, jx,  kx, rlayerdirclr, tlayerdirclr,   &
                       rlayerdifclr, tlayerdifclr, tlayerdeclr,   &
                       sfcalb_dir,  sfcalb_dif, cloud_in_column,  &
                       reflectanceclr, transmittanceclr, tr_dirclr)
                endif

!---------------------------------------------------------------------- 
!    weight and sum the reflectance and transmittance to calculate the 
!    band values.                                                     
!-------------------------------------------------------------------
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    wtfac_p(i,j) = wtfreq(np)*gausswt(ng)*   &
                                                  cosangsolar_p(i,j,ng)
                  end do
                end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                do k = KSRAD,KERAD+1
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if (daylight(i,j) ) then
                        sumtr(i,j,k,nz) = sumtr(i,j,k,nz) +    &
                                      transmittance(i,j,k)*wtfac_p(i,j)
                        sumtr_dir(i,j,k,nz) = sumtr_dir(i,j,k,nz) +  &
                                        tr_dir(i,j,k)*wtfac_p(i,j)
                        sumre(i,j,k,nz) = sumre(i,j,k,nz) +     &
                                       reflectance(i,j,k)*wtfac_p(i,j)
                      endif
                    end do
                  end do
                end do
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if (daylight(i,j) ) then
                        sumtr_dir_up(i,j,nz) = sumtr_dir_up(i,j,nz) + &
                          tr_dir(i,j,KERAD+1)*sfcalb_dir(i,j)*wtfac_p(i,j)
                      endif
                    end do
                  end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
                if (Rad_control%do_totcld_forcing) then
                  do k = KSRAD,KERAD+1
                    do j=JSRAD,JERAD
                      do i=ISRAD,IERAD
                        if (cloud_in_column(i,j)) then
                          sumtrclr(i,j,k,nz) =    &
                                sumtrclr(i,j,k,nz) +   &
                                transmittanceclr(i,j,k)* wtfac_p(i,j) 
                          sumtr_dir_clr(i,j,k,nz) = &
                                sumtr_dir_clr(i,j,k,nz) +  &
                                tr_dirclr(i,j,k)*wtfac_p(i,j)
                          sumreclr(i,j,k,nz) = sumreclr(i,j,k,nz) +   &
                                reflectanceclr(i,j,k)*wtfac_p(i,j)
                        else if (daylight(i,j) ) then
                          sumtrclr(i,j,k,nz) = sumtrclr(i,j,k,nz) +   &
                                transmittance(i,j,k)*wtfac_p(i,j)
                          sumtr_dir_clr(i,j,k,nz) =    &
                             sumtr_dir_clr(i,j,k,nz) + tr_dir(i,j,k)*&
                                                            wtfac_p(i,j)
                          sumreclr(i,j,k,nz) = sumreclr(i,j,k,nz) +   &
                                        reflectance(i,j,k)*wtfac_p(i,j)
                        endif
                      end do
                    end do
                  end do
                endif
              end do ! end of nz loop
            end do    ! end of gaussian loop
          end do  ! end of frequency points in the band loop
 
!----------------------------------------------------------------------
!    normalize the solar flux in the band to the appropriate value for  
!    the given total solar insolation.                                 
!---------------------------------------------------------------------
          do nz = 1,nzens
            if (Rad_control%hires_coszen) then
              fracday_p(:,:) = Astro%fracday_p(:,:,nz)
            else
              fracday_p(:,:) = Astro%fracday(:,:)
            endif
            solarflux_p(:,:) = fracday_p(:,:)*  &
                                   Solar_spect%solflxband(nband)*  &
                                                      ssolar/solflxtotal_local
 
            if (nband == Solar_spect%visible_band_indx) then
              Sw_output%bdy_flx(:,:,1,nz) =   &
                  Sw_output%bdy_flx(:,:,1,nz) + sumre(:,:,1,nz)*   &
                                                        solarflux_p(:,:)
              Sw_output%bdy_flx(:,:,3,nz) =    &
                  Sw_output%bdy_flx(:,:,3,nz) + sumtr(:,:,KERAD+1,nz)*&
                                                solarflux_p(:,:) -  &
                                                sumre(:,:,KERAD+1,nz)* &
                                                solarflux_p(:,:) 
            endif
            if (nband == onepsix_band_indx) then
               Sw_output%bdy_flx(:,:,2,nz) =     &
                   Sw_output%bdy_flx(:,:,2,nz) + sumre(:,:,1,nz)*  &
                                                        solarflux_p(:,:)
               Sw_output%bdy_flx(:,:,4,nz) =    &
                   Sw_output%bdy_flx(:,:,4,nz) + sumtr(:,:,KERAD+1,nz)*&
                                                 solarflux_p(:,:) - &
                                                 sumre(:,:,KERAD+1,nz)*&
                                                 solarflux_p(:,:)
            endif
          
!-------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!--------------------------------------------------------------------
            if (do_esfsw_band_diagnostics) then
              do k = KSRAD,KERAD+1
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    dfswband(i,j,k,nz) = sumtr(i,j,k,nz)*   &
                                                       solarflux_p(i,j) 
                    ufswband(i,j,k,nz) = sumre(i,j,k,nz)*      &
                                                       solarflux_p(i,j)
                  end do
                end do
              end do
            endif
 
!----------------------------------------------------------------------
!    sum the band fluxes and heating rates to calculate the total       
!    spectral values.                                                  
!------------------------------------------------------------------
            do k = KSRAD,KERAD+1
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (daylight(i,j) ) then
                    Sw_output%dfsw (i,j,k,nz) =   &
                       Sw_output%dfsw(i,j,k,nz) + sumtr(i,j,k,nz)*&
                                                       solarflux_p(i,j)
                    Sw_output%ufsw (i,j,k,nz) =   &
                       Sw_output%ufsw(i,j,k,nz) + sumre(i,j,k,nz)*  &
                                                       solarflux_p(i,j)
                    fswband(i,j,k,nz) = ((sumre(i,j,k,nz)*  &
                                          solarflux_p(i,j)) - &
                                         (sumtr(i,j,k,nz)*    &
                                                      solarflux_p(i,j)))
                    Sw_output%fsw(i,j,k,nz) = Sw_output%fsw(i,j,k,nz) +&
                                                      fswband(i,j,k,nz)
                  endif
                end do
              end do
            end do
 
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (daylight(i,j) ) then
                  Sw_output%dfsw_dir_sfc(i,j,nz) =   &
                            Sw_output%dfsw_dir_sfc(i,j,nz) +   &
                              sumtr_dir(i,j,KERAD+1,nz)*solarflux_p(i,j)
                  Sw_output%ufsw_dir_sfc(i,j,nz) =   &
                            Sw_output%ufsw_dir_sfc(i,j,nz) +   &
                              sumtr_dir_up(i,j,nz)*solarflux_p(i,j)
                  Sw_output%ufsw_dif_sfc(i,j,nz) =   &
                             Sw_output%ufsw_dif_sfc(i,j,nz) +   &
                                 sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                endif
              end do
            end do

            if (nband > NIRBANDS) then
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (daylight(i,j) ) then
                    Sw_output%dfsw_vis_sfc(i,j,nz) =   &
                          Sw_output%dfsw_vis_sfc(i,j,nz) +   &
                                sumtr(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc(i,j,nz) =   &
                           Sw_output%ufsw_vis_sfc(i,j,nz) +   &
                                 sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%dfsw_vis_sfc_dir(i,j,nz) =   &
                            Sw_output%dfsw_vis_sfc_dir(i,j,nz) +   &
                              sumtr_dir(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc_dir(i,j,nz) =   &
                            Sw_output%ufsw_vis_sfc_dir(i,j,nz) +   &
                              sumtr_dir_up(i,j,nz)*solarflux_p(i,j)
                    Sw_output%ufsw_vis_sfc_dif(i,j,nz) =   &
                             Sw_output%ufsw_vis_sfc_dif(i,j,nz) +   &
                                 sumre(i,j,KERAD+1,nz)*solarflux_p(i,j)
                  endif
                end do
              end do
            endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
            do k = KSRAD,KERAD
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (daylight(i,j) ) then
                    hswband(i,j,k,nz) = (fswband(i,j,k+1,nz) -    &
                                      fswband(i,j,k,nz) )*gocpdp(i,j,k)
                    Sw_output%hsw(i,j,k,nz) =    &
                          Sw_output%hsw(i,j,k,nz) + hswband(i,j,k,nz)
                  endif
                end do
              end do
            end do

!----------------------------------------------------------------------
!    calculate the band fluxes and heating rates.                       
!---------------------------------------------------------------------
            if (nprofile == 1) then  ! clr sky need be done only for 
                                     ! first cloud profile
              if (Rad_control%do_totcld_forcing) then
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (daylight(i,j) ) then
                      Sw_output%dfsw_dir_sfc_clr(i,j,nz) =   &
                         Sw_output%dfsw_dir_sfc_clr(i,j,nz) +   &
                          sumtr_dir_clr(i,j,KERAD+1,nz)*solarflux_p(i,j)
                    endif
                  end do
                end do
                if (nband > NIRBANDS) then
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if (daylight(i,j) ) then
                        Sw_output%dfsw_vis_sfc_clr(i,j,nz) =   &
                           Sw_output%dfsw_vis_sfc_clr(i,j,nz) +   &
                              sumtrclr(i,j,KERAD+1,nz)*solarflux_p(i,j)
                      endif
                    end do
                  end do
                endif
                if (nband == Solar_spect%visible_band_indx) then
                  Sw_output%bdy_flx_clr(:,:,1,nz) =      &
                           sumreclr(:,:,1,nz)*solarflux_p(:,:)
                  Sw_output%bdy_flx_clr(:,:,3,nz) =    &
                          sumtrclr(:,:,KERAD+1,nz)*solarflux_p(:,:) - &
                          sumreclr(:,:,KERAD+1,nz)*solarflux_p(:,:) 
                endif
                if (nband == onepsix_band_indx) then
                  Sw_output%bdy_flx_clr(:,:,2,nz) =    &
                          sumreclr(:,:,1,nz)*solarflux_p(:,:)
                  Sw_output%bdy_flx_clr(:,:,4,nz) =    &
                         sumtrclr(:,:,KERAD+1,nz)*solarflux_p(:,:)  -  &
                         sumreclr(:,:,KERAD+1,nz)*solarflux_p(:,:) 
                endif
          
                if (do_esfsw_band_diagnostics) then
                  do k = KSRAD,KERAD+1
                    do j=JSRAD,JERAD
                      do i=ISRAD,IERAD
                        dfswbandclr(i,j,k,nz) =     &
                                   sumtrclr(i,j,k,nz)*solarflux_p(i,j)
                        ufswbandclr(i,j,k,nz) =    &
                                   sumreclr(i,j,k,nz)*solarflux_p(i,j)
                      end do
                    end do
                  end do
                endif

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total     
!    spectral values.                                                 
!----------------------------------------------------------------------c
                do k = KSRAD,KERAD+1
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      Sw_output%dfswcf(i,j,k,nz) =    &
                              Sw_output%dfswcf(i,j,k,nz) +  &
                                    sumtrclr(i,j,k,nz)*solarflux_p(i,j)
                      Sw_output%ufswcf(i,j,k,nz) =      &
                              Sw_output%ufswcf(i,j,k,nz) +  &
                                    sumreclr(i,j,k,nz)*solarflux_p(i,j)
                      fswbandclr(i,j,k,nz) =    &
                            (sumreclr(i,j,k,nz)*solarflux_p(i,j)) - &
                            (sumtrclr(i,j,k,nz)*solarflux_p(i,j))
                      Sw_output%fswcf(i,j,k,nz) =    &
                                    Sw_output%fswcf(i,j,k,nz) +    &
                                                  fswbandclr(i,j,k,nz)
                    end do
                  end do
                end do

!----------------------------------------------------------------------c
!    sum the band fluxes and heating rates to calculate the total    
!    spectral values.                                               
!----------------------------------------------------------------------c
                do k = KSRAD,KERAD
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      hswbandclr(i,j,k,nz) =    &
                                (fswbandclr(i,j,k+1,nz) -      &
                                    fswbandclr(i,j,k,nz))*gocpdp(i,j,k)
                      Sw_output%hswcf(i,j,k,nz) =   &
                                   Sw_output%hswcf(i,j,k,nz) +    &
                                                   hswbandclr(i,j,k,nz)
                    end do
                  end do
                end do
              endif
            endif ! (nprofile == 1)
          end do ! (nz loop)
        end do   ! (end of band loop)
      end do   ! (end of nprofile loop)

!----------------------------------------------------------------------
!    if the ica calculation was being done, the fluxes and heating rates
!    which have been summed over nprofiles cloud profiles must be 
!    averaged.
!------------------------------------------------------------------
      if (do_ica_calcs) then
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            Sw_output%dfsw_dir_sfc (i,j,:) =   &
                      Sw_output%dfsw_dir_sfc(i,j,:)*profiles_inverse
            Sw_output%ufsw_dif_sfc (i,j,:) =   &
                      Sw_output%ufsw_dif_sfc(i,j,:)*profiles_inverse
            Sw_output%dfsw_vis_sfc (i,j,:) =   &
                      Sw_output%dfsw_vis_sfc(i,j,:)*profiles_inverse
            Sw_output%ufsw_vis_sfc (i,j,:) =   &
                      Sw_output%ufsw_vis_sfc(i,j,:)*profiles_inverse
            Sw_output%dfsw_vis_sfc_dir (i,j,:) =   &
                      Sw_output%dfsw_vis_sfc_dir(i,j,:)*profiles_inverse
            Sw_output%ufsw_vis_sfc_dif (i,j,:) =   &
                      Sw_output%ufsw_vis_sfc_dif(i,j,:)*profiles_inverse
            Sw_output%bdy_flx(i,j,:,:) =  &
                      Sw_output%bdy_flx(i,j,:,:)*profiles_inverse
          end do
        end do
        do k = KSRAD,KERAD+1
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              Sw_output%dfsw (i,j,k,:) = Sw_output%dfsw(i,j,k,:)*  &
                                       profiles_inverse
              Sw_output%ufsw (i,j,k,:) = Sw_output%ufsw(i,j,k,:)*  &
                                       profiles_inverse
              Sw_output%fsw(i,j,k,:) = Sw_output%fsw(i,j,k,:)*  &
                                     profiles_inverse
            end do
          end do
        end do
        do k = KSRAD,KERAD
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              Sw_output%hsw(i,j,k,:) = Sw_output%hsw(i,j,k,:)*  &
                                       profiles_inverse
            end do
          end do
        end do
      endif

      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          if (daylight(i,j) ) then
            Sw_output%dfsw_dif_sfc(i,j,: ) =   &
                              Sw_output%dfsw(i,j,KERAD+1,: ) -   &
                                      Sw_output%dfsw_dir_sfc(i,j,: )
          endif
        end do
      end do

      if (Rad_control%do_totcld_forcing) then
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            if (daylight(i,j) ) then
              Sw_output%dfsw_dif_sfc_clr(i,j,: ) =   &
                                Sw_output%dfswcf(i,j,KERAD+1,:) -   &
                                Sw_output%dfsw_dir_sfc_clr(i,j,:)
            endif
          end do
        end do
      endif

      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          if (daylight(i,j) ) then
            Sw_output%dfsw_vis_sfc_dif(i,j,:) =   &
                                Sw_output%dfsw_vis_sfc(i,j,:) -   &
                                Sw_output%dfsw_vis_sfc_dir(i,j,:) 
          endif
        end do
      end do

!--------------------------------------------------------------------
!    convert sw fluxes to cgs and then back to  mks units.
!---------------------------------------------------------------------
      Sw_output%fsw(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%fsw(:,:,:,:))
      Sw_output%dfsw(:,:,:,:) =    &
                            1.0E-03*(1.0E+03*Sw_output%dfsw(:,:,:,:))
      Sw_output%ufsw(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufsw(:,:,:,:))
      if (Rad_control%do_totcld_forcing) then
        Sw_output%fswcf(:,:,:,:) =   &
                            1.0E-03*(1.0E+03*Sw_output%fswcf(:,:,:,:))
        Sw_output%dfswcf(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%dfswcf(:,:,:,:))
        Sw_output%ufswcf(:,:,:,:) =     &
                            1.0E-03*(1.0E+03*Sw_output%ufswcf(:,:,:,:))
      endif

!---------------------------------------------------------------------


end subroutine swresf



!####################################################################

subroutine esfsw_driver_end

!---------------------------------------------------------------------
!    esfsw_driver_end is the destructor for esfsw_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('esfsw_driver_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules that this module initialized.
!--------------------------------------------------------------------
      call esfsw_parameters_end

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

end subroutine esfsw_driver_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#################################################################
! <SUBROUTINE NAME="compute_aerosol_optical_props">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comput(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine compute_aerosol_optical_props    &
          (Atmos_input, Aerosol, Aerosol_props, including_volcanoes,  &
           Aerosol_diags, r, including_aerosols, naerosol_optical, &
           daylight, aeroextopdep, aerosctopdep, aeroasymfac)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

type(atmos_input_type),        intent(in)    :: Atmos_input
type(aerosol_type),            intent(in)    :: Aerosol     
type(aerosol_properties_type), intent(in)    :: Aerosol_props
logical,                       intent(in)    :: including_volcanoes
type(aerosol_diagnostics_type),intent(inout) :: Aerosol_diags
real,dimension(:,:,:,:),       intent(inout) :: r
logical,                       intent(in)    :: including_aerosols  
integer,                       intent(in)    :: naerosol_optical
logical,dimension(:,:),        intent(in)    :: daylight
real, dimension(:,:,:,:),      intent(out)   :: aeroextopdep, &
                                                aerosctopdep, &
                                                aeroasymfac 


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Aerosol_props  aerosol radiative property input data for the 
!                     radiation package
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 
      real, dimension (size(Atmos_input%temp,3)-1)  :: &
                      arprod, asymm,   arprod2,      deltaz, &     
                      sum_g_omega_tau, sum_ext,      sum_sct

      integer, dimension (size (Atmos_input%press, 3)-1 ) ::   &
                      opt_index_v3, opt_index_v4, &
                      opt_index_v5, opt_index_v6, &
                      opt_index_v7, opt_index_v8, &
                      opt_index_v9, opt_index_v10, &
                      irh

      real, dimension (naerosol_optical)   ::            &
                      aerext,          aerssalb,       aerasymm

      real        :: aerext_i, aerssalb_i, aerasymm_i
      integer     :: j, i, k, nband, nsc
      integer     :: israd, jsrad, ierad, jerad, ksrad, kerad
      integer     :: nextinct  !variable to pass extinction to chemistry

!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size(Atmos_input%temp,1)
      jerad = size(Atmos_input%temp,2)
      kerad = size(Atmos_input%temp,3) - 1

      naerosoltypes_used = size(Aerosol%aerosol,4)

      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if (daylight(i,j) .or. Sw_control%do_cmip_diagnostics) then
            deltaz(:) = Atmos_input%deltaz(i,j,:)

            do nband = 1,Solar_spect%nbands
              if (including_aerosols) then                           
                aerext(:) = Aerosol_props%aerextband(nband,:)
                aerssalb(:) = Aerosol_props%aerssalbband(nband,:)
                aerasymm(:) = Aerosol_props%aerasymmband(nband,:)

!-------------------------------------------------------------------
!    define the local variables for the band values of aerosol and 
!    cloud single scattering parameters.                               
!    note: the unit for the aerosol extinction is kilometer**(-1).     
!--------------------------------------------------------------------
                do k = KSRAD,KERAD
                  irh(k) = MIN(100, MAX( 0,     &
                      NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
                  opt_index_v3(k) = &
                              Aerosol_props%sulfate_index (irh(k), &
                                            Aerosol_props%ivol(i,j,k))
                  opt_index_v4(k) =    &
                              Aerosol_props%omphilic_index( irh(k) )
                  opt_index_v5(k) =    &
                              Aerosol_props%bcphilic_index( irh(k) )
                  opt_index_v6(k) =    &
                              Aerosol_props%seasalt1_index( irh(k) )
                  opt_index_v7(k) =    &
                              Aerosol_props%seasalt2_index( irh(k) )
                  opt_index_v8(k) =    &
                              Aerosol_props%seasalt3_index( irh(k) )
                  opt_index_v9(k) =    &
                              Aerosol_props%seasalt4_index( irh(k) )
                  opt_index_v10(k) =    &
                              Aerosol_props%seasalt5_index( irh(k) )
                end do

!---------------------------------------------------------------------
!    calculate scattering properties for all aerosol constituents 
!    combined.
!---------------------------------------------------------------------
                do k = KSRAD,KERAD
                  sum_g_omega_tau(k) = 0.0
                  sum_ext(k) = 0.
                  sum_sct(k) = 0.
                end do
                do nsc = 1,NAEROSOLTYPES_USED
                  if (Aerosol_props%optical_index(nsc) > 0) then
                    aerext_i =     &
                            aerext(Aerosol_props%optical_index(nsc))
                    aerssalb_i =     &
                            aerssalb(Aerosol_props%optical_index(nsc))
                    aerasymm_i =     &
                            aerasymm(Aerosol_props%optical_index(nsc))
                    do k = KSRAD,KERAD
                      arprod(k) =    &
                            aerext_i*(1.e3*Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb_i*arprod(k)
                      asymm(k)   = aerasymm_i
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + aerssalb_i*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +     &
                                      aerasymm_i*(aerssalb_i*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                                     Aerosol_props%sulfate_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v3(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v3(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v3(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) +    &
                                    aerssalb(opt_index_v3(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                               aerasymm(opt_index_v3(k))*  &
                                   (aerssalb(opt_index_v3(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                                     Aerosol_props%bc_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v3(k)) *    &
                                    (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v3(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v3(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) +       &
                                    aerssalb(opt_index_v3(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) + &
                               aerasymm(opt_index_v3(k)) * &
                                   (aerssalb(opt_index_v3(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                               Aerosol_props%omphilic_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v4(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v4(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v4(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v4(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +   &
                               aerasymm(opt_index_v4(k))*     &
                                  (aerssalb(opt_index_v4(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                                  Aerosol_props%bcphilic_flag) then
                    if (Rad_control%using_im_bcsul) then
                      do k = KSRAD,KERAD
                        arprod(k) = aerext(opt_index_v3(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                        arprod2(k) = aerssalb(opt_index_v3(k))*arprod(k)
                        asymm(k)   = aerasymm(opt_index_v3(k))
                        sum_ext(k) = sum_ext(k) + arprod(k)
                        sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v3(k))*arprod(k)
                        sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                               aerasymm(opt_index_v3(k)) * &
                                  (aerssalb(opt_index_v3(k))*arprod(k))
                      end do
                    else  ! (using_im_bcsul)
                      do k = KSRAD,KERAD
                        arprod(k) = aerext(opt_index_v5(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                        arprod2(k) = aerssalb(opt_index_v5(k))*arprod(k)
                        asymm(k)   = aerasymm(opt_index_v5(k))
                        sum_ext(k) = sum_ext(k) + arprod(k)
                        sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v5(k))*arprod(k)
                        sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                              aerasymm(opt_index_v5(k)) * &
                                 (aerssalb(opt_index_v5(k))*arprod(k))
                      end do
                    endif  !(using_im_bcsul)
                  else if (Aerosol_props%optical_index(nsc) == &
                                Aerosol_props%seasalt1_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v6(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v6(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v6(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v6(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                                aerasymm(opt_index_v6(k)) * &
                                  (aerssalb(opt_index_v6(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                              Aerosol_props%seasalt2_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v7(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v7(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v7(k))
                      sum_ext(k) = sum_ext(k) +  arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v7(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) + &
                               aerasymm(opt_index_v7(k)) * &
                                 (aerssalb(opt_index_v7(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                               Aerosol_props%seasalt3_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v8(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v8(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v8(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v8(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                              aerasymm(opt_index_v8(k)) * &
                                  (aerssalb(opt_index_v8(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                               Aerosol_props%seasalt4_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v9(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v9(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v9(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                     aerssalb(opt_index_v9(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +  &
                                aerasymm(opt_index_v9(k))*  &
                                   (aerssalb(opt_index_v9(k))*arprod(k))
                    end do
                  else if (Aerosol_props%optical_index(nsc) == &
                             Aerosol_props%seasalt5_flag) then
                    do k = KSRAD,KERAD
                      arprod(k) = aerext(opt_index_v10(k)) *    &
                                   (1.e3 * Aerosol%aerosol(i,j,k,nsc))
                      arprod2(k) = aerssalb(opt_index_v10(k))*arprod(k)
                      asymm(k)   = aerasymm(opt_index_v10(k))
                      sum_ext(k) = sum_ext(k) + arprod(k)
                      sum_sct(k) = sum_sct(k) + &
                                    aerssalb(opt_index_v10(k))*arprod(k)
                      sum_g_omega_tau(k) = sum_g_omega_tau(k) +&
                              aerasymm(opt_index_v10(k)) * &
                                 (aerssalb(opt_index_v10(k))*arprod(k))
                    end do
                  endif

                  if (Sw_control%do_cmip_diagnostics) then
                    if (nband == Solar_spect%visible_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,1) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,1) =    &
                                            arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,1) = asymm(:)
                    endif
                    if (nband == Solar_spect%eight70_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,6) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,6) =    &
                                            arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,6) = asymm(:)
                    endif
                    if (nband == Solar_spect%one_micron_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,2) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,2) =    &
                                               arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,2) = asymm(:)
                    endif
                    if (nband == Solar_spect%w340_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,7) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,7) =    &
                                                arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,7) = asymm(:)
                    endif
                    if (nband == Solar_spect%w380_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,8) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,8) =    &
                                                arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,8) = asymm(:)
                    endif
                    if (nband == Solar_spect%w440_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,9) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,9) =    &
                                               arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,9) = asymm(:)
                    endif
                    if (nband == Solar_spect%w670_band_indx) then
                      Aerosol_diags%extopdep(i,j,:,nsc,10) = arprod(:)
                      Aerosol_diags%absopdep(i,j,:,nsc,10) =    &
                                                arprod(:) - arprod2(:)
                      Aerosol_diags%asymdep(i,j,:,nsc,10) = asymm(:)
                    endif
                  endif
                end do

!----------------------------------------------------------------------
!    add the effects of volcanic aerosols, if they are to be included.
!    include generation of diagnostics in the visible (0.55 micron) and
!    nir band (1.0 micron).
!----------------------------------------------------------------------
                if (including_volcanoes) then
                  do k = KSRAD,KERAD
                    sum_ext(k) = sum_ext(k) +    &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                    sum_sct(k) = sum_sct(k) +    &
                                 Aerosol_props%sw_ssa(i,j,k,nband)*  &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                    sum_g_omega_tau(k) =   &
                                 sum_g_omega_tau(k) +&
                                 Aerosol_props%sw_asy(i,j,k,nband)* &
                                 Aerosol_props%sw_ssa(i,j,k,nband)*  &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                    if (Sw_control%do_cmip_diagnostics) then
                      if (nband == Solar_spect%visible_band_indx) then
                           Aerosol_diags%extopdep_vlcno(i,j,k,1) =   &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                           Aerosol_diags%absopdep_vlcno(i,j,k,1) =   &
                            (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
                                Aerosol_props%sw_ext(i,j,k,nband)*  &
                                deltaz(k)
                      endif
                      if (nband == Solar_spect%eight70_band_indx) then
                           Aerosol_diags%extopdep_vlcno(i,j,k,3) =   &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                           Aerosol_diags%absopdep_vlcno(i,j,k,3) =   &
                            (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
                                Aerosol_props%sw_ext(i,j,k,nband)*  &
                                deltaz(k)
                      endif
                      if (nband == Solar_spect%one_micron_indx) then
                           Aerosol_diags%extopdep_vlcno(i,j,k,2) =   &
                                 Aerosol_props%sw_ext(i,j,k,nband)*  &
                                 deltaz(k)
                           Aerosol_diags%absopdep_vlcno(i,j,k,2) =   &
                            (1.0 - Aerosol_props%sw_ssa(i,j,k,nband))*&
                                Aerosol_props%sw_ext(i,j,k,nband)*  &
                                deltaz(k)
                      endif
                    endif
                  end do
                endif   ! (including_volcanoes)
!
!----------------------------------------------------------------------
                do k = KSRAD,KERAD
                  aeroextopdep(i,j,k,nband) = sum_ext(k) 
                  aerosctopdep(i,j,k,nband) = sum_sct(k) 
                  aeroasymfac(i,j,k,nband) = sum_g_omega_tau(k) / &
                                                (sum_sct(k) + 1.0e-30 )
                end do
              else  ! (if not including_aerosols)
                do k = KSRAD,KERAD
                  aeroextopdep(i,j,k,nband) = 0.0                    
                  aerosctopdep(i,j,k,nband) = 0.0                  
                  aeroasymfac(i,j,k,nband) = 0.0                 
                end do
              endif ! (including_aerosols)
            end do ! (nband)
          endif  ! (daylight or cmip_diagnostics)

          if (including_volcanoes) then
            if (do_coupled_stratozone ) then
              nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
              if (nextinct  /= NO_TRACER) &
                    r(i,j,:,nextinct) = Aerosol_props%sw_ext(i,j,:,4)
            endif  
          endif  

        end do ! (i loop)
      end do   ! (j loop)

!---------------------------------------------------------------------
!

!---------------------------------------------------------------------


end subroutine compute_aerosol_optical_props

!#################################################################
! <SUBROUTINE NAME="compute_gas_props">
!  <OVERVIEW>
!   Subroutine that uses the delta-eddington technique in conjunction
!   with a multi-band parameterization for h2o+co2+o2+o3 absorption
!   in the solar spectrum to derive solar fluxes and heating rates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine calculates optical depth, single scattering albedo,
!    asymmetry parameter of a layer based on gaseous absorbers,
!    clouds, aerosols, and rayleigh scattering. It then uses delta-
!    eddington technique to calculate radiative flux at each layer. 
!    Doubling and adding technique is used to combine the layers
!    and calculate flux at TOA and surface and heating rate. This
!    subroutine allocates a substantial amount of memory and deallocates
!    the allocated memory at the end of the subroutine.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call comput(is, ie, js, je, Atmos_input, Surface, Rad_gases, Aerosol, 
!               Astro, &
!               Cldrad_props, Cld_spec, Sw_output)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    Atmos_input_type variable containing the atmospheric
!    input fields on the radiation grid 
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data for shortwave radiation calculation
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!    Astronomy_type variable containing the astronomical
!    input fields on the radiation grid  
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!    Radiative_gases_type variable containing the radiative 
!    gas input fields on the radiation grid 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!    The cloud radiative property input fields on the
!    radiation grid
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!    The shortwave radiation calculation result
!  </INOUT>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface data as boundary condition to radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification data as initial condition to radiation
!  </IN>
! </SUBROUTINE>

subroutine compute_gas_props (Atmos_input, Rad_gases, Astro,   &
                              daylight, gasopdep)

!----------------------------------------------------------------------
!    comput uses the delta-eddington technique in conjunction with a    
!    multiple-band parameterization for h2o+co2+o2+o3 absorption to   
!    derive solar fluxes and heating rates.                             
!    notes: drops are assumed if temp>273.15K, ice crystals otherwise.
!-------------------------------------------------------------------

type(atmos_input_type),        intent(in)    :: Atmos_input
type(radiative_gases_type),    intent(in)    :: Rad_gases   
type(astronomy_type),          intent(in)    :: Astro
logical, dimension(:,:),       intent(in)    :: daylight
real, dimension(:,:,:,:,:),    intent(out)   :: gasopdep              


!-------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Astro          astronomy_type structure
!                                                                 
!   intent(inout) variables:
!
!      Sw_output         shortwave radiation output data
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables:
 

      real, dimension (size(Atmos_input%temp,3)-1)  :: &
                   deltaz,     qo3,         rh2o,                    &
                   efftauo2,   efftauco2,   efftauch4,   efftaun2o, &
                   wh2ostr,    wo3,         wo2,         quenchfac, &
                   opdep,      delpdig,     deltap,      tco2,    &
                   tch4,       tn2o,        to2,         wh2o

           
      real, dimension (size(Atmos_input%temp,3))  :: &
            alphaco2,        alphaco2str,    alphao2,          &
            alphao2str,      alphach4,       alphach4str,      &
            alphan2o,        alphan2ostr,    scale,            &
            scalestr,        totco2,         totco2str,        &
            toto2,           toto2str,       totch4,           &
            totch4str,       totn2o,         totn2ostr,        &
            press,           pflux,          pflux_mks,        &
            temp,            z

      real :: cosangsolar
      real :: denom
      real :: wtquench
      real :: rrvco2 
      real :: rrvch4, rrvn2o

      integer  :: j, i, k, ng, nband, kq
      integer  :: np, nf
      integer  :: israd, jsrad, ierad, jerad, ksrad, kerad


!-----------------------------------------------------------------------
!     local variables:
!
!       aeramt
!       sum_g_omega_tau
!       opt_index_v3
!       irh
!    etc.
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------
      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size(Atmos_input%temp,1)
      jerad = size(Atmos_input%temp,2)
      kerad = size(Atmos_input%temp,3) - 1

!---------------------------------------------------------------------
!    initialize local variables.                                        
!------------------------------------------------------------------
      alphaco2   (1) = 0.0
      alphaco2str(1) = 0.0
      alphao2    (1) = 0.0
      alphao2str (1) = 0.0
      alphach4   (1) = 0.0
      alphach4str(1) = 0.0
      alphan2o   (1) = 0.0
      alphan2ostr(1) = 0.0

      rrvco2 = Rad_gases%rrvco2
      rrvch4 = Rad_gases%rrvch4
      rrvn2o = Rad_gases%rrvn2o

      do j = JSRAD,JERAD
        do i = ISRAD,IERAD
          if ( daylight(i,j) ) then

!---------------------------------------------------------------------
!  convert to cgs and then back to mks for consistency with previous 
!---------------------------------------------------------------------
            do k = KSRAD,KERAD+1
              press(k) = 0.1*(10.0*Atmos_input%press(i,j,k))
              pflux(k) =     (10.0*Atmos_input%pflux(i,j,k))
              temp (k) = Atmos_input%temp  (i,j,k)
            end do
            do k = KSRAD,KERAD
              rh2o  (k) = Atmos_input%rh2o  (i,j,k)
              qo3   (k) = Rad_gases%qo3(i,j,k)
              deltaz(k) = Atmos_input%deltaz(i,j,k)
            end do

!----------------------------------------------------------------------c
!    define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
            pflux_mks(:) = pflux(:)*1.0E-1

            do k = KSRAD+1,KERAD+1
              deltap(k-1) = pflux_mks(k) - pflux_mks(k-1)
              delpdig(k-1) = deltap(k-1)/ GRAV
              scalestr(k) = pflux_mks(k) 
              scale(k) = scalestr(k)*pflux_mks(k)/pstd_mks
            end do
 

            do k = KSRAD,KERAD
              wh2ostr(k) = rh2o(k)*delpdig(k)
              wo3(k)     = qo3(k)*delpdig(k)
              wo2(k) = o2mixrat*(WTMO2/WTMAIR)*delpdig(k)
            end do
  
!---------------------------------------------------------------------
!    if quenching factor effects are desired, calculate the height above
!    the surface of the model flux levels.
!---------------------------------------------------------------------
            if (do_quench) then
              z(KERAD+1) = 0.0
              do k = KERAD,KSRAD,-1
                z(k) = z(k+1) + deltaz(k)
              end do
          
!---------------------------------------------------------------------
!    define the quenching factor for each grid point.
!---------------------------------------------------------------------
              do k = KSRAD,KERAD
                if (z(k) < co2_quenchfac_height(1) ) then
                  quenchfac(k) = 1.0
                else if (z(k) > co2_quenchfac_height(30) ) then 
                  quenchfac(k) = 0.037
                else
                  do kq = 1,29
                    if (z(k) > co2_quenchfac_height(kq) .and. &
                        z(k) <= co2_quenchfac_height(kq+1)) then
                      wtquench = (z(k) - co2_quenchfac_height(kq))/ &
                                 (co2_quenchfac_height(kq+1) - &
                                  co2_quenchfac_height(kq))
                      quenchfac(k) = (1. - wtquench)*   &
                                           co2_quenchfac(kq) +   &
                                      wtquench*co2_quenchfac(kq+1)
                      exit
                    endif
                  end do
                endif
              end do
            else
              quenchfac(:) = 1.0
            endif !(do_quench)


            do ng = 1,NSOLWG
              cosangsolar  = Astro%cosz(i,j)
              if (cosangsolar == 0.0) cosangsolar = 1.0

!----------------------------------------------------------------------c
!    define the scaled and unscaled co2 and o2 pathlengths in 
!    centimeter-atm, and the unscaled h2o and o3 amounts in   
!    kgrams/meter**2. 
!    cm-atm needed as units because of c2co2 having those units.
!----------------------------------------------------------------------c
              denom = 1.0/(GRAV*rhoair*cosangsolar*2.0)
              do k = KSRAD+1,KERAD+1
                totco2(k) = 1.0E+02*rrvco2*scale(k)*denom       
                totco2str(k) = 2.0E+02*rrvco2*scalestr(k)*denom     
                toto2(k) = 1.0E+02*o2mixrat*scale(k)*denom      
                toto2str(k) = 2.0E+02*o2mixrat*scalestr(k)*denom      
                if (do_ch4_sw_effects) then
                  totch4(k) = 1.0E+02*rrvch4*scale(k)*denom     
                  totch4str(k) = 2.0E+02*rrvch4*scalestr(k)*denom     
                endif
                if (do_n2o_sw_effects) then
                  totn2o(k) = 1.0E+02*rrvn2o*scale(k)*denom     
                  totn2ostr(k) = 2.0E+02*rrvn2o*scalestr(k)*denom      
                endif
              end do

              np = 0
              do nband = 1, NBANDS

!-------------------------------------------------------------------
!    define the h2o scaled gas amounts in kgrams/meter**2            
!---------------------------------------------------------------------
                if (nband <= nh2obands) then
                  do k = KSRAD,KERAD
                    wh2o(k) = rh2o(k)*delpdig(k)*   &
                        exp(powph2o(nband)*alog(press(k)*p0h2o(nband)))
                  end do
 
!---------------------------------------------------------------------
!    calculate the "effective" co2, o2, ch4 and n2o gas optical depths 
!    for the appropriate absorbing bands.                               
!    note: for large zenith angles, alpha can exceed 1. In this case,a
!    the optical depths are set to the previous layer values.          
!-------------------------------------------------------------------
                  if ( c1co2(nband).ne.1.0E-99 ) then
                    do k = KSRAD+1,KERAD+1
                      if (totco2(k) < totco2max(nband) .and.  &
                          totco2str(k) < totco2strmax(nband))  then
                        alphaco2(k) =     &
                             c1co2(nband)*exp(c3co2(nband)* &
                                 alog((totco2(k) + c2co2(nband))))  -  &
                                                          c4co2(nband)
                        alphaco2str(k) = &
                          c1co2str(nband)*exp(c3co2str(nband)*  &
                            alog((totco2str(k) + c2co2str(nband)))) - &
                                                        c4co2str(nband)
                        tco2(k-1) =      &
                             (1.0 - alphaco2(k))*   &
                                            (1.0 - alphaco2str(k))/ &
                             ((1.0 - alphaco2(k-1))*    &
                                            (1.0 - alphaco2str(k-1)))
                        efftauco2(k-1) = -cosangsolar*alog( tco2(k-1))
                      else if (k > KSRAD+1) then
                        efftauco2(k-1) = efftauco2(k-2)
                      else
                        efftauco2(k-1) = 0.0
                      end if
                    end do
                  else    !( c1co2(nband).ne.1.0E-99 ) 
                    efftauco2(:) = 0.0
                  end if  !( c1co2(nband).ne.1.0E-99 ) 

                  if (do_ch4_sw_effects) then
                    if (c1ch4(nband).ne.1.0E-99 ) then
                      do k = KSRAD+1,KERAD+1
                        if (totch4(k) < totch4max(nband) .and.  &
                              totch4str(k) < totch4strmax(nband))  then
                           alphach4(k) =    &
                              c1ch4(nband)*exp(c3ch4(nband)*&
                              alog((totch4(k) + c2ch4(nband))))  -   &
                                                           c4ch4(nband)
                           alphach4str(k) = &
                            c1ch4str(nband)*exp(c3ch4str(nband)*  &
                            alog((totch4str(k) + c2ch4str(nband)))) - &
                                                        c4ch4str(nband)
                           tch4(k-1) = &
                                  (1.0 - alphach4(k))*    &
                                           (1.0 - alphach4str(k))/ &
                                   ((1.0 - alphach4(k-1))*   &
                                           (1.0 - alphach4str(k-1)))
                           efftauch4(k-1) = -cosangsolar*alog(tch4(k-1))
                         else if (k > KSRAD+1) then
                           efftauch4(k-1) = efftauch4(k-2)
                         else
                           efftauch4(k-1) = 0.0
                         end if
                       end do
                     else    !( c1ch4(nband).ne.1.0E-99 )
                       efftauch4(:) = 0.0
                     end if  !( c1ch4(nband).ne.1.0E-99 )
                   else    !do_ch4 = .false.
                     efftauch4(:) = 0.0
                   end if

                   if (do_n2o_sw_effects) then
                     if ( c1n2o(nband).ne.1.0E-99 ) then
                       do k = KSRAD+1,KERAD+1
                         if (totn2o(k) < totn2omax(nband) .and.  &
                              totn2ostr(k) < totn2ostrmax(nband)) then
                           alphan2o(k) = &
                                c1n2o(nband)*exp(c3n2o(nband)* &
                                   alog((totn2o(k) +c2n2o(nband)))) -  &
                                                          c4n2o(nband)
                           alphan2ostr(k) = &
                            c1n2ostr(nband)*exp(c3n2ostr(nband)*  &
                            alog((totn2ostr(k) + c2n2ostr(nband)))) -  &
                                                       c4n2ostr(nband)
                           tn2o(k-1) = &
                                    (1.0 - alphan2o(k)) *  &
                                              (1.0 - alphan2ostr(k))/ &
                                    (( 1.0 - alphan2o(k-1)) *  &
                                            (1.0 - alphan2ostr(k-1)))
                           efftaun2o(k-1) = -cosangsolar*alog(tn2o(k-1))
                         else if (k > KSRAD+1) then
                           efftaun2o(k-1) = efftaun2o(k-2)
                         else
                           efftaun2o(k-1) = 0.0
                         end if
                       end do
                     else    !( c1n2o(nband).ne.1.0E-99 )
                       efftaun2o(:) = 0.0
                     end if  !( c1n2o(nband).ne.1.0E-99 )
                   else  !do_n2o = .false.
                     efftaun2o(:) = 0.0
                   end if

                  if ( c1o2(nband).ne.1.0E-99 ) then
                    do k = KSRAD+1,KERAD+1
                      if (toto2(k) .lt. toto2max(nband) .and.   &
                          toto2str(k) .lt. toto2strmax(nband)) then
                        alphao2(k) = c1o2(nband)*exp( c3o2(nband)* &
                                     alog((toto2(k) + c2o2(nband)))) - &
                                                            c4o2(nband)
                        alphao2str( k) = &
                            c1o2str(nband)*exp(c3o2str(nband)*  &
                                alog((toto2str(k) + c2o2str(nband)))) &
                                                        - c4o2str(nband)
                        to2(k-1) = &
                               (1.0 - alphao2(k))*  &
                                    (1.0 - alphao2str(k) )/ &
                                ((1.0 - alphao2(k-1)) *  &
                                            (1.0 - alphao2str(k-1)))
                        efftauo2(k-1) = -cosangsolar*alog(to2(k-1))
                      else if (k.gt.KSRAD+1) then
                        efftauo2(k-1) = efftauo2(k-2)
                      else
                        efftauo2(k-1) = 0.0
                      end if
                    end do
                  else   !  ( c1o2(nband).ne.1.0E-99 ) 
                    efftauo2(:) = 0.0
                  end if  !  ( c1o2(nband).ne.1.0E-99 ) 
                end if  ! (nband <= nh2obands)
 
!---------------------------------------------------------------------
!    calculate the "effective" o2 gas optical depths for the Schuman- 
!    Runge band.                                                        
!-------------------------------------------------------------------
                if ( nband.EQ.NBANDS ) then
                  do k = KSRAD+1,KERAD+1
                    if ( toto2str(k).lt.toto2strmaxschrun) then
                      alphao2str(k) =  &
                           c1o2strschrun*exp( c3o2strschrun*&
                              alog((toto2str(k) + c2o2strschrun))) - &
                                                       c4o2strschrun
                      to2(k-1) = &
                          (1.0 - alphao2str(k))/(1.0 - alphao2str(k-1)) 
                      efftauo2(k-1) =  -cosangsolar*alog(to2(k-1) )
                      if (do_herzberg) then
                        efftauo2(k-1) = efftauo2(k-1) +     &
                                                 wo2(k-1)*herzberg_fac
                      endif
                    else if (k.gt.KSRAD+1) then
                      efftauo2(k-1) = efftauo2(k-2)
                    else
                      efftauo2(k-1) = 0.0
                    end if
                  end do
                end if

                do nf =1,nfreqpts(nband)
                  np = np + 1

!---------------------------------------------------------------------
!    define the h2o + o3 gas optical depths.                           
!--------------------------------------------------------------------
                  if (strterm(np)) then
                    opdep(:) = kh2o(np)*wh2ostr(:) + ko3(np)*wo3(:)
                  else
                    opdep(:) = kh2o(np)*wh2o(:) + ko3(np)*wo3(:)
                  end if

                  gasopdep(i,j,:,np,ng) =    &
                           opdep(:) + quenchfac(:)*efftauco2(:) +   &
                              efftauo2(:) + efftauch4(:) + efftaun2o(:)

                end do  ! (nf loop)
              end do   ! (nband loop)
            end do  ! (ng loop)
          endif  ! (daylight)
        end do ! (i loop)
      end do ! (j loop)

!---------------------------------------------------------------------



end subroutine compute_gas_props


!#####################################################################
!<SUBROUTINE NAME="adding">
! <OVERVIEW>
!  Subroutine that implements doubling and adding technique to combine
!  multiple atmospheric layers to calculate solar fluxes
! </OVERVIEW>
! <DESCRIPTION>
!  This subroutine implements the standard doubling and adding
!  technique to combine reflectance and transmittance of multiple 
!  atmospheric layers to compute solar flux and heating rate.
! </DESCRIPTION>
! <TEMPLATE>
!  call adding ( ix, jx, kx, &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,  &
!                tlayerde, sfcalb, calc_flag, reflectance,   &
!                transmittance)
! </TEMPLATE>
! <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="jx" TYPE="integer">
!  jx is the current latitudinal index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="kx" TYPE="integer">
!  ix is the current vertical index in the physics cell being
!  integrated.
! </IN>
! <IN NAME="rlayerdir" TYPE="real">
!  layer reflectivity to direct incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to direct incident beam
! </IN>
! <IN NAME="rlayerdif" TYPE="real">
!  layer reflectivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerdir" TYPE="real">
!  layer transmissivity to diffuse incident beam
! </IN>
! <IN NAME="tlayerde" TYPE="real">
!  layer diffuse transmissivity to direct incident beam
! </IN>
! <IN NAME="sfcalb" TYPE="real">
!  surface albedo
! </IN>
! <IN NAME="calcflag" TYPE="integer">
!  flag to indicate columns where adding is to be done
! </IN>
! <OUT NAME="reflectance" TYPE="real">
!  diffuse reflectance at a level
! </OUT>
! <OUT NAME="transmittance" TYPE="real">
!  diffuse transmittance at a level
! </OUT>
!</SUBROUTINE>
!

subroutine adding (ix, jx, kx, rlayerdir, tlayerdir, rlayerdif,   &
                   tlayerdif, tlayerde, sfcalb_dir, sfcalb_dif,  &
                   calc_flag, reflectance, transmittance, tr_dir)
 
!-------------------------------------------------------------------
!    adding calculates the reflection and transmission at flux levels 
!    from the direct and diffuse values of reflection and transmission
!    in the corresponding layers using the adding method.           
!    references:                                                        
!    bowen, m.m., and v. ramaswamy, effects of changes in radiatively
!        active species upon the lower stratospheric temperatures.,    
!        j. geophys. res., 18909-18921, 1994.                         
!--------------------------------------------------------------------

integer, intent(in)                    :: ix, jx, kx
real, dimension(:,:,:),   intent(in)   :: rlayerdir, rlayerdif, &
                                          tlayerdir, tlayerdif, & 
                                          tlayerde
real, dimension (:,:),    intent(in)   :: sfcalb_dir, sfcalb_dif
logical, dimension (:,:), intent(in)   :: calc_flag
real, dimension(:,:,:),   intent(out)  :: reflectance, transmittance, &
                                          tr_dir

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx        dimensions of current physics window            
!    rlayerdir       layer reflectivity to a direct incident beam      
!    tlayerdir       layer transmissivity to a direct incident beam   
!    rlayerdif       layer reflectivity to a diffuse incident beam  
!    tlayerdif       layer transmissivity to a diffuse incident beam  
!    tlayerde        layer transmissivity (non-scattered) to the direct 
!                    incident beam                                 
!    sfcalb_dir      surface albedo, direct beam 
!    sfcalb_dif      surface albedo, diffuse beam
!    calc_flag       flag to indicate columns where adding is to be 
!                    done. calculations not done in "dark" columns and 
!                    on clr sky pass in columns without any clouds.
!
!  intent(out) variables:
!
!    reflectance     reflectance of the scattered radiation at a level 
!    transmittance   transmittance of the scattered radiation at a level
!    tr_dir
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:
 
      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)+1) ::  &
                            raddupdif2, raddupdir2

      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)  ) ::  &
                                      radddowndif2,  tadddowndir2

      real :: dm1tl2, dm2tl2, rdm2tl2, dm32, dm3r2, dm3r1p2, alpp2, &
              raddupdif2p, raddupdir2p, tlevel2p, radddowndifm, &
              tadddowndirm
      integer     ::  k, j, i

!-------------------------------------------------------------------
!   local variables:
!
!      raddupdif2
!      raddupdir2
!      tlevel2
!      radddowndif2
!      tadddowndir2
!      tlayerdif2
!      tlayerdir2
!      rlayerdif2
!      rlayerdir2
!      tlayerde2
!      dm1tl2
!      dm2tl2
!      rdm2tl2
!      dm32
!      dm3r2
!      dm3r1p2
!      alpp2
!      raddupdif2
!      raddupdir2p
!      tlevel2p
!      radddowndifm
!      tadddowndirm
!      i,j,k
!
!--------------------------------------------------------------------

!----------------------------------------------------------------------c
!    initialization for the surface layer.                           
!----------------------------------------------------------------------c
      do j=1,jx        
        do i=1,ix
          if (calc_flag(i,j) ) then
 
!------------------------------------------------------------------ 
!    add the inhomogeneous layers upward from the surface to the top of
!    the atmosphere.                                                  
!    radiation incident from above for diffuse beam, reflection of  
!    direct beam and conversion to diffuse.                           
!--------------------------------------------------------------------
            raddupdif2p = sfcalb_dif(i,j)
            raddupdir2p = sfcalb_dir(i,j)
            do k = kx, 1,-1
              dm2tl2    = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*  &
                          raddupdif2p )
              rdm2tl2    = dm2tl2*raddupdif2p     
              raddupdif2(k) = rlayerdif(i,j,k) + tlayerdif(i,j,k)*   &
                              rdm2tl2    
              raddupdir2(k) = rlayerdir(i,j,k) + tlayerde(i,j,k)*   &
                              raddupdir2p* dm2tl2 +   &     
                              (tlayerdir(i,j,k) - tlayerde(i,j,k))*  &
                              rdm2tl2   
              raddupdir2p = raddupdir2(k)
              raddupdif2p = raddupdif2(k)
            end do
 
!---------------------------------------------------------------------
!    define the direct transmittance. add the inhomogeneous layers 
!    downward from the second layer to the surface. radiation incident
!    from below for diffuse beam, transmission of direct beam and 
!    conversion to diffuse.                             
!-------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    initialization for the first scattering layer.                   
!-------------------------------------------------------------------
            tlevel2p         = tlayerde(i,j,1)
            radddowndifm    =  rlayerdif(i,j,1)
            tadddowndirm    =  tlayerdir(i,j,1)
            do k= 2,kx    
              dm1tl2 = tlayerdif(i,j,k)/(1.0 - rlayerdif(i,j,k)*  &
                       radddowndifm)
              radddowndif2(k) = rlayerdif(i,j,k) + radddowndifm* &
                                tlayerdif(i,j,k)*dm1tl2      
              tadddowndir2(k) = tlevel2p*(tlayerdir(i,j,k) + &
                                rlayerdir(i,j,k)*radddowndifm* &
                                dm1tl2) + (tadddowndirm -  &
                                tlevel2p)*dm1tl2           

!---------------------------------------------------------------------
!    add downward to calculate the resultant reflectances and           
!    transmittances at flux levels.                                    
!------------------------------------------------------------------
              dm32  = 1.0/(1.0 - raddupdif2(k)*radddowndifm)
              dm3r2 = dm32*radddowndifm      
              dm3r1p2 = 1.0 + raddupdif2(k)*dm3r2   
              alpp2 = (tadddowndirm - tlevel2p)*dm32   
              transmittance(i,j,k) = (tlevel2p*(1.0 + raddupdir2(k)* &
                                      dm3r2) + alpp2)
              tr_dir(i,j,k) = tlevel2p
              reflectance(i,j,k) = (tlevel2p*raddupdir2(k)*   &
                                    dm3r1p2 + alpp2*   &
                                    raddupdif2(k))
              tlevel2p = tlevel2p*tlayerde (i,j,k) 
              radddowndifm = radddowndif2(k)
              tadddowndirm = tadddowndir2(k)
            end do
!! CORRECT ???
!           dm32  = 1.0/(1.0 - sfcalb(i,j)*radddowndifm)
            dm32          = 1.0/(1.0 - sfcalb_dif(i,j)*   &
                               radddowndifm       )
            dm3r2 = dm32*radddowndifm       
!! CORRECT ???
!           dm3r1p2 = 1.0 + sfcalb(i,j)*dm3r2         
            dm3r1p2          = 1.0 + sfcalb_dif(i,j) * dm3r2
            alpp2 = (tadddowndirm - tlevel2p)*dm32          
            transmittance(i,j,kx+1) = (tlevel2p*(1.0 +   &
!! CORRECT ???
!                                      sfcalb(i,j)* &
!12-08-03:  CHANGE THIS TO _dir as per SMF  sfcalb_dif(i,j)* &
                                       sfcalb_dir(i,j)* &
                                       dm3r2) + alpp2)
            tr_dir(i,j,kx+1) = tlevel2p
            reflectance(i,j,kx+1) = (tlevel2p*  &
!! CORRECT ???
!                                   sfcalb(i,j)*   &
                                    sfcalb_dir(i,j)* &
                                     dm3r1p2 + alpp2* &
!! CORRECT ???
!                                sfcalb(i,j) )
                                    sfcalb_dif(i,j))  
            reflectance(i,j,1) = raddupdir2p         
            transmittance(i,j,1) = 1.0
            tr_dir(i,j,1) = 1.0
          endif
        end do
      end do

!------------------------------------------------------------------


end subroutine adding 



!####################################################################
! <SUBROUTINE NAME="deledd">
!  <OVERVIEW>
!   Subroutine that calculates reflectivity and transmissivity in a
!   scattering layer using delta-eddington method
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine takes layer optical depth, single scattering abledo,
!   and asymmetry parameter, using delta-eddington method, to calculate
!   direct/diffuse reflectivity/transmissivity to direct/diffuse incident
!   radiation. The approximation uses the strong forward scattering of
!   aerosol particles.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deledd (ix, jx, kx,  &
!                taustr, omegastr, gstr, cosang, ng , daylight,  &
!                rlayerdir, tlayerdir, rlayerdif, tlayerdif,   &
!                tlayerde,  cloud)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!  ix is the current longitudinal index in the physics cell being
!  integrated.
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   jx is the current latitudinal index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   ix is the current vertical index in the physics cell being
!   integrated.
!  </IN>
!  <IN NAME="taustr" TYPE="real">
!   the scaled optical depth, true optical depth normalized using
!   delta-eddington approximation
!  </IN>
!  <IN NAME="omegastr" TYPE="real">
!   the scaled single-scattering albedo
!  </IN>
!  <IN NAME="gstr" TYPE="real">
!   the scaled asymmetry factor
!  </IN>
!  <IN NAME="cosang" TYPE="real">
!   cosine of the solar zenith angle
!  </IN>
!  <IN NAME="ng" TYPE="real">
!   the number of gaussian angles to compute the diurnally    
!   averaged solar radiation (=1 unless lswg = true)
!  </IN>
!  <IN NAME="cloud" TYPE="real">
!   flag for existence of a cloud (used only in 'ovc' mode)
!  </IN>
!  <OUT NAME="rlayerdir" TYPE="real">
!   layer reflectivity to direct incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to direct incident beam
!  </OUT>
!  <OUT NAME="rlayerdif" TYPE="real">
!   layer reflectivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerdir" TYPE="real">
!   layer transmissivity to diffuse incident beam
!  </OUT>
!  <OUT NAME="tlayerde" TYPE="real">
!   layer diffuse transmissivity to direct incident beam
!  </OUT>
! </SUBROUTINE>
!
subroutine deledd (ix, jx, kx, taustr, omegastr, gstr, cosang, ng, &
                   daylight, rlayerdir, tlayerdir, tlayerde,   &
                   rlayerdif, tlayerdif, cloud)
 
!---------------------------------------------------------------------- 
!    deledd calculates the reflection and transmission in the 
!    scattering layers using the delta-eddington method.         
!    references:                                                   
!      joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington
!      approximation for radiative flux transfer.,j. atmos. sci.,33,  
!      2452-2459, 1976.                                              
!-------------------------------------------------------------------

integer,                   intent(in)              :: ix, jx, kx
real, dimension(:,:,:),    intent(inout)           :: taustr, omegastr
real, dimension(:,:,:),    intent(in)              :: gstr
real, dimension(:,:),    intent(in)                ::  cosang
integer,                   intent(in)              :: ng
logical, dimension(:,:),   intent(in)              :: daylight
real, dimension(:,:,:),    intent(out)             :: rlayerdir,   &
                                                      tlayerdir,   &
                                                      tlayerde
real, dimension(:,:,:),    intent(inout), optional :: rlayerdif,   &
                                                      tlayerdif
logical, dimension(:,:,:), intent(in), optional    :: cloud          

!----------------------------------------------------------------------
!  intent(in) variables:
!
!    ix,jx,kx
!    gstr        the scaled asymmetry factor                       
!    cosang      the cosine of the solar zenith angle    
!    ng          the number of gaussian angles to compute the diurnally 
!                averaged solar radiation (=1 unless lswg = true)       
!    daylight
!
!  intent(inout) variables:
!
!    taustr      the scaled extinction optical depth                    
!    omegastr    the scaled single-scattering albedo               
!
!  intent(out) variables:
!
!    rlayerdir   the layer reflectivity to a direct incident beam      
!    tlayerdir   the layer transmissivity to a direct incident beam   
!    rlayerdif   the layer reflectivity to a diffuse incident beam   
!    tlayerdif   the layer transmissivity to a diffuse incident beam
!    tlayerde    the layer transmissivity (non-scattered) to the direct 
!                incident beam                                       
!
! intent(in),optional:
!
!    cloud       flag for existence of a cloud (used only in 'ovc' mode)
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real        :: qq(7), rr(5), ss(8), tt(8), ww(4)
      real        :: rsum, tsum
      real        :: onedi3 = 1.0/3.0           
      real        :: twodi3 = 2.0/3.0             
      integer     :: k, i, ns, j, nn, ntot

      integer, dimension(ix, jx, kx) :: cld_index

      real,    dimension(ix)                  ::   &
                                          gstr2, taustr2, omegastr2, &
                                           cosangzk2, rlayerdir2,    &
                                           tlayerde2, tlayerdir2, &
                                           sumr, sumt


!----------------------------------------------------------------------
!  local variables:
!
!      qq
!      rr
!      ss
!      tt
!      ww
!      rsum
!      tsum
!      alpha
!      onedi3
!      twodi3
!      i,j,k
!      ns
!      nn
!      ntot
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=1,kx         
        do j=1,jx         

!---------------------------------------------------------------------
!    overcast sky mode. note: in this mode, the delta-eddington method
!    is performed only for spatial points containing a cloud.   
!-------------------------------------------------------------------
          nn = 0
          if (present(cloud)) then
            do i=1,ix          
              if (cloud(i,j,k) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                gstr2(nn) = gstr(i,j,k)
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j)

!----------------------------------------------------------------------
!    note: the following are done to avoid the conservative scattering 
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                      
!----------------------------------------------------------------------c
                if (omegastr2(nn) >= 1.0) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do

!----------------------------------------------------------------------c
!    clear sky mode. note: in this mode, the delta-eddington method is 
!    performed for all spatial points.                 
!----------------------------------------------------------------------c
          else
            do i=1,ix         
              if (daylight(i,j) ) then
                nn = nn + 1
                cld_index(i,j,k) = nn
                gstr2(nn) = gstr(i,j,k)
                taustr2(nn) = taustr(i,j,k)
                omegastr2(nn) = omegastr(i,j,k)
                cosangzk2(nn) = cosang(i,j   )

!----------------------------------------------------------------------c
!    note: the following are done to avoid the conservative scattering  
!    case, and to eliminate floating point errors in the exponential 
!    calculations, respectively.                    
!----------------------------------------------------------------------c
                if (omegastr2(nn) >= 1.0) omegastr2(nn) = 9.9999999E-01
                if (taustr2(nn) >= 1.0E+02) taustr2(nn) = 1.0E+02
              endif
            end do
          endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          ntot = nn

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          do nn=1,ntot      

!----------------------------------------------------------------------c
!    direct quantities                                            
!----------------------------------------------------------------------c
            ww(1) = omegastr2(nn)
            ww(2) = gstr2(nn)
            ww(3) = taustr2(nn)
            ww(4) = cosangzk2(nn)

            qq(1)     = 3.0 * ( 1.0 - ww(1) )
            qq(2)         = 1.0 - ww(1) * ww(2)
            qq(3)     = qq(1)/qq(2)
            qq(4) = sqrt( qq(1) * qq(2) )
            qq(5) = sqrt (qq(3))
            qq(6) = 1.0 + twodi3 * qq(5)         
            qq(7) = 1.0 - twodi3 * qq(5)       

            rr(1) = 1./qq(6)
            rr(2) = qq(7)*rr(1)
            rr(3) = exp( -ww(3)          * qq(4) )
            rr(4) = 1.0/rr(3)
            rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) * rr(2) )

            ss(1) = 0.75 * ww(1)/(1.0 - (qq(4)*ww(4)      ) ** 2 )
            ss(2) = ss(1)*ww(4)*( 1.0 + ww(2)*qq(1)*onedi3)
            ss(3) = ss(1)*(1.0 + ww(2)*qq(1)*ww(4)** 2 )
            ss(4) = ss(2) - twodi3*ss(3)     
            ss(5) = ss(2) + twodi3 * ss(3)     
            ss(6) = exp( -ww(3)          / ww(4) )
            ss(7) = (ss(4)*ss(6) - ss(5)*rr(3)*rr(2))*rr(5)
            ss(8) = (ss(5) - qq(7)*ss(7))*rr(1)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
            rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - ss(4)
            tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                               qq(7) * rr(4) * ss(7) -  &
                               ss(5) * ss(6) ) + ss(6) )
            tlayerde2(nn) = ss(6)

!----------------------------------------------------------------------c
!    diffuse quantities                                       
!    notes: the number of streams for the diffuse beam is fixed at 4.   
!    this calculation is done only for ng=1.                 
!----------------------------------------------------------------------c
            if (present (tlayerdif) .and. present(rlayerdif)) then
              if ( ng.eq.1 ) then   
                rsum = 0.0
                tsum = 0.0
                do ns = 1,NSTREAMS
                  tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                          cosangstr(ns) ) ** 2 )
                  tt(2) = tt(1) * cosangstr(ns) * ( 1.0 +  &
                          ww(2)        * qq(1) * onedi3 )
                  tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                          cosangstr(ns) ** 2 )
                  tt(4) = tt(2) - twodi3 * tt(3)
                  tt(5) = tt(2) + twodi3 * tt(3)
                  tt(6) = exp( -ww(3)          / cosangstr(ns) )
                  tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                          rr(3) * rr(2)   ) * rr(5)
                  tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)
                  if (nstr4) then
                    rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))* &
                           wtstr(ns)*cosangstr(ns)
                    tsum = tsum + ((rr(3)*qq(6)*tt(8) +   &
                                    qq(7)*rr(4)*tt(7) -   &
                                    tt(5)*tt(6)) + tt(6))*  &
                                    wtstr(ns)*cosangstr(ns)
                  else 
                    rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))
                    tsum = tsum + ( (rr(3)*qq(6)*tt(8) +    &
                                     qq(7)*rr(4)*tt(7) -   &
                                     tt(5)*tt(6)) + tt(6))
                  endif
                end do
                sumr(nn) = rsum
                sumt(nn) = tsum
              endif  !  ng == 1
            endif ! (present (tlayerdiff))
          end do  ! ntot loop

!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------
          if (present(cloud)) then
            do i=1,ix           
              if (cloud(i,j,k) ) then
                rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                if (present(tlayerdif)) then
                  if (ng .eq. 1) then
                    rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
                    tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
                  endif
                endif
              endif
            end do

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
          else
            do i=1,ix            
              if (daylight(i,j) ) then
                rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
                tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
                tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
                if (present(tlayerdif)) then
                  if (ng .eq. 1) then
                    rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
                    tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
                  endif
                endif
              endif
            end do
          endif
        end do
      end do

!---------------------------------------------------------------------
 
end subroutine deledd



!#####################################################################


                   end module esfsw_driver_mod

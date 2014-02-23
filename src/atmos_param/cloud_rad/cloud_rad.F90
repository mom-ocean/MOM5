!FDOC_TAG_GFDL


                 module cloud_rad_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Steve Klein
! </CONTACT>
! <OVERVIEW>
!        The cloud radiation module uses the stored values of the
!     prognostic cloud variables, and computes the cloud albedo and
!     absorption for the two shortwave bands (ultra-violet/visible and
!     near-infrared), the longwave cloud emissivity, and the 
!     fractional areas covered by clouds.
! </OVERVIEW>
! <DESCRIPTION>
!      The cloud radiation module condenses the cloud information 
!     provided by the stratiform cloud scheme and converts it into
!     the areas covered by, the water paths and the effective particle 
!     sizes of liquid and ice. This cloud information is stored into 
!     cloud blocks which are assumed to be randomly overlapped (done 
!     in CLOUD_ORGANIZE subroutine). From these, the single-scattering 
!     albedo, asymmetry parameter, and optical depth for the two short 
!     wave bands and the longwave cloud emissivity for each cloud are 
!     calculated in the subroutine CLOUD_OPTICAL_PROPERTIES. Finally, 
!     the subroutine CLOUD_RAD takes the shortwave cloud properties 
!     and converts them using the Delta-Eddington solution to albedo 
!     and absorption in each of the shortwave bands.
!
!     In CLOUD_OPTICAL_PROPERTIES, the parameterization of Slingo (1989)
!     and Ebert and Curry (1992) are used for the shortwave properties of 
!     liquid and ice clouds, respectively.  For the longwave cloud 
!     emissivity, the empirical observation result of Stephens (1978) is
!     used for liquid clouds whereas the parameterization of Ebert and
!     Curry (1992) is used for ice clouds.
!
!     In CLOUD_ORGANIZE, the effective radius for liquid clouds is 
!     calculated using the parameterization of Martin et al. (1994)
!     whereas the effective radius of ice clouds is parameterized using
!     that of Donner et al. (1997).
!
!  
! </DESCRIPTION>
!

!  <DIAGFIELDS>
!  ************************Note***********************
!   This part of the documentation needs to be updated
!  ***************************************************
!
!  Diagnostic fields may be output to a netcdf file by specifying the
!  module name cloud_rad and the desired field names (given below)
!  in file diag_table. See the documentation for diag_manager.
!  
!  Diagnostic fields for module name: cloud_rad
!  
!     nisccp       frequency of sunlit times at the times of the radiation
!                  calculation at each point {fraction} [real,dimension(:,:)]
!  
!     pc#tau%      where # is a number from 1 to 7
!                  and   % is a number from 1 to 7
!                  {fraction} [real,dimension(:,:)]
!  
!                  Thus there are 49 diagnostic fields of this type.  All
!                  of them are necessary to receive the complete decomposition
!                  of clouds visible from space into the ISCCP categories.
!  
!                  The 7 cloud top pressure ("pc") categories and 7 optical
!                  depth ("tau") categories are defined as:
!  
!                  pc #      pc range (mb)    tau %        tau range
!                  ----    ----------------   -----    ---------------------
!  
!                   1              pc < 180     1     0.0    < tau < taumin 
!                   2        180 < pc < 310     2     taumin < tau < 1.3
!                   3        310 < pc < 440     3     1.3    < tau < 3.6
!                   4        440 < pc < 560     4     3.6    < tau < 9.4
!                   5        560 < pc < 680     5     9.4    < tau < 23
!                   6        680 < pc < 800     6     23     < tau < 60
!                   7        800 < pc                 60     < tau
!  
!                  What is saved in these diagnostics is the time mean of
!                  the area covered by clouds of this type when the sun is
!                  above the horizon. This is done so that the calculation 
!                  will mimic the ISCCP product which is broken down into 
!                  these categories only for sunlit places.
!  
!                  NOTE:  TO DETERMINE THE MEAN AREA COVERED BY A CLOUD TYPE 
!                         WHEN THE SUN IS ABOVE THE HORIZON YOU MUST DIVIDE
!                         BY NISCCP:
!  
!                         area of cloud type pc#tau% =   pc#tau% / nisccp
!  
!     aice         fractional area of sunlit clouds seen from space whose cloud 
!                  top contains ice. {fraction} [real,dimension(:,:)]
!  
!     reffice      time mean ice effective radius of cloud tops visible from
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!  
!                  NOTE:  THUS THE TIME MEAN CLOUD TOP EFFECTIVE RADIUS OF CLOUD 
!                         TOPS WITH ICE VISIBLE FROM SPACE IS:
!  
!                         mean reffice  =    reffice /  aice
!        
!     aliq         fractional area of sunlit clouds seen from space whose cloud 
!                  top contains liquid. {fraction} [real,dimension(:,:)]
!  
!     reffliq      time mean cloud droplet effective radius of cloud tops 
!                  visible from space including areas where there is no such 
!                  cloud {microns} [real,dimension(:,:)]
!     
!                  NOTE:  mean reffliq  =    reffliq / aliq
!  
!     alow         fractional area of sunlit clouds seen from space whose cloud 
!                  tops are low (pc > 680 mb). {fraction} [real,dimension(:,:)]
!  
!     tauicelow    time mean optical depth of ice for cloud tops visible from 
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!    
!     tauliqlow    time mean optical depth of liquid for cloud tops visible from 
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!     
!     tlaylow      time mean of the low level mean temperature (pc > 680 mb) 
!                  when low cloud tops are visible from space including times 
!                  where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!        
!     tcldlow      time mean of the cloud top temperature for cloud tops visible 
!                  from space including times where there is no such cloud 
!                  {microns}  [real,dimension(:,:)]
!        
!                  NOTE:  mean tauicelow  =    tauicelow / alow
!                         mean tauliqlow  =    tauliqlow / alow
!                         mean tlaylow    =    tlaylow   / alow
!                         mean tcldlow    =    tcldlow   / alow
!  
! </DIAGFIELDS>
!  
!  
! <DATASET NAME="">
!
! </DATASET>
!  
!  
! <INFO>
!  
!   <REFERENCE>            
! The shortwave properties of liquid clouds come from:
! 
 !     Slingo, A., 1989: A GCM parameterization for the shortwave 
 !     radiative properties of water clouds. J. Atmos. Sci., vol. 46, 
 !     pp. 1419-1427.
!
!   </REFERENCE>

!   <REFERENCE>            
! The shortwave and longwave properties of ice clouds come from:
! 
 !     Ebert, E. E. and J. A. Curry, 1992: A parameterization of ice cloud
 !     optical properties for climate models. J. Geophys. Res., vol. 97,
 !     D1, pp. 3831-3836.
!
!   </REFERENCE>

!   <REFERENCE>            
! The longwave emissivity parameterization of liquid clouds comes from:
! 
 !     Stephens, G. L., 1978: Radiation profiles in extended water clouds.
 !     II: Parameterization schemes. J. Atmos. Sci., vol. 35, 
 !     pp. 2123-2132.
!
!   </REFERENCE>

!   <REFERENCE>            
! The parameterization of liquid cloud effective radius comes from:
! 
 !     Martin, G. M., D. W. Johnson, and A. Spice, 1994: The measurement 
 !     and parameterization of effective radius of droplets in warm stratocumulus
 !     clouds. J. Atmos. Sci, vol 51, pp. 1823-1842.
!
!   </REFERENCE>

!   <REFERENCE>            
! The parameterization of ice cloud effective radius comes from:
! 
 !     Donner, L. J., C. J. Seman, B. J. Soden, R. S. Hemler, J. C. Warren,
 !     J. Strom, and K.-N. Liou, 1997: Large-scale ice clouds in the GFDL
 !     SKYHI general circulation model. J. Geophys. Res., vol. 102, D18,
 !     pp. 21,745-21,768.
!
!   </REFERENCE>

!   <REFERENCE>            
! The algorithm to reproduce the ISCCP satellite view of clouds comes from:
! 
 !     Klein, S. A., and C. Jakob, 1999: Validation and sensitivities of 
 !     frontal clouds simulated by the ECMWF model. Monthly Weather Review,
 !     127(10),  2514-2531.
!
!   </REFERENCE>

!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!     
!   </NOTE>
!   <FUTURE>The optical depth and particle size for every model level will
!     become a diagnostic output field.               </FUTURE>

! </INFO>

!   shared modules:

use  mpp_mod,             only:  input_nml_file
use  fms_mod,             only:  file_exist, fms_init,       &
                                 stdlog, mpp_pe, mpp_root_pe, &
                                 open_namelist_file, &
                                 write_version_number,  &
                                 error_mesg, FATAL,     &
                                 close_file,  &
                                 check_nml_error
use  constants_mod,       only:  RDGAS, GRAV, TFREEZE, DENS_H2O, &
                                 constants_init, pi
use  gamma_mg_mod,        ONLY : gamma_mg,  gamma_mg_init,  gamma_mg_end
use mg_const_mod,         ONLY : mg_const_init, rhow, di_mg, ci_mg
use  diag_manager_mod,    only:  diag_manager_init,    &
                                 register_diag_field, send_data
use  time_manager_mod,    only:  time_type, time_manager_init

implicit none
private

!---------------------------------------------------------------------
!    cloud_rad_mod does the following:           
!     
!    (a)  subroutine cloud_summary3 returns cloud specification var-
!         iables that are used in calculating the cloud radiative prop-
!         erties. these include cloud locations, water paths and effect-
!         ive particle sizes for use in determining bulk properties and
!         concentrations and drop sizes if microphysically-based prop-
!         erties are desired.
!    (b)  subroutine lw_emissivity returns the long wave cloud emis-
!         sivity and subroutine sw_optical_properties returns the cloud 
!         reflectivities and absorptions in the nir and uv spectral
!         bands when using non-microphysically-based cloud radiative
!         properties.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!------------ version number for this module -------------------------
        
character(len=128) :: version = '$Id: cloud_rad.F90,v 20.0 2013/12/13 23:09:10 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!---------------------------------------------------------------------- 
!--------- interfaces --------

public     &
         cloud_rad_init, &
         cloud_rad_end, &
         cloud_summary, &
         lw_emissivity, &
         sw_optical_properties, &
         cloud_summary3, &
         cloud_rad_k_diag,  &
         cloud_rad, &
         snow_and_rain

!---------------------------------------------------------------------
!    public subroutines:
!
!      cloud_rad_init
!                        Initializes values of qmin, N_land, and 
!                        N_ocean using values from strat_cloud namelist
!                        as well as reads its own namelist variables. 
!                        In addition, it registed diagnostic fields
!                        if needed, and returns the value of the
!                        cloud overlap to strat_cloud.
!      cloud_summary
!                        This is the main driver program of the module
!      cloud_rad   
!                        this solves for the radiative properties of the
!                        clouds given the cloud optical properties
!                        (tau,w0,gg) for each cloud using either a 
!                        Delta-Eddington solution (default) or the
!                        two stream approximation.
!      cloud_optical_properties
!                        for each cloud this calculates the mixed phase
!                        values of the optical depth, the single scat-
!                        tering albedo, and the asymmetry parameter 
!                        (tau, w0,and g) for the visible and near infra-
!                        red bands.  It also computes the longwave 
!                        emissivity of each cloud.
!
!----------------------------------------------------------------------

private     &
         max_rnd_overlap, rnd_overlap, cloud_rad_k
     
!---------------------------------------------------------------------
!    private subroutines:
!
!       max_rnd_overlap
!       rnd_overlap
!       cloud_rad_k
!      cloud_organize
!                        for each cloud this computes the cloud amount,
!                        the liquid and ice water paths, the effective
!                        particle sizes, the tops and bottoms of clouds.
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!-------- namelist  ---------

integer      :: overlap = 1
logical      :: l2strem = .false.
real         :: taucrit = 1.
logical      :: adjust_top = .true.
real         :: scale_factor = 0.85
real         :: qamin = 1.E-2
logical      :: do_brenguier = .true.
real         :: N_min = 1.e6

logical      :: snow_in_cloudrad = .false.
logical      :: rain_in_cloudrad = .false.
real         :: min_diam_ice_mg = 10.e-6
real         :: min_diam_drop_mg = 2.e-6
real         :: max_diam_drop_mg = 50.e-6
real         :: dcs_mg =  200.e-6
real         :: Ni_min = 10.
!--------------------------------------------------------------------
!    namelist variables:
!
!      overlap        integer variable indicating which overlap 
!                     assumption to use:
!                     overlap = 1. means condensate in adjacent levels 
!                                  is treated as part of the same cloud
!                                  i.e. maximum-random overlap
!                     overlap = 2. means condensate in adjacent levels 
!                                  is treated as different clouds
!                                  i.e. random overlap
!      l2strem        logical variable indicating which solution for 
!                     cloud radiative properties is being used.
!                     l2strem = T  2 stream solution
!                     l2strem = F  Delta-Eddington solution
!                     Note that IF l2strem = T then the solution does 
!                     not depend on solar zenith angle
!      taucrit        critical optical depth for switching direct beam 
!                     to diffuse beam for use in Delta-Eddington 
!                     solution [ dimensionless] 
!      adjust_top     logical variable indicating whether or not to use 
!                     the code which places the top and bottom of the 
!                     cloud at the faces which are most in view from
!                     the top and bottom of the cloud block. this is 
!                     done to avoid undue influence of very small cloud
!                     fractions. if true this adjustment of tops is 
!                     performed; if false this is not performed.
!      scale_factor   factor which multiplies actual cloud optical 
!                     depths to account for the plane-parallel homo-
!                     genous cloud bias  (e.g. Cahalan effect).
!                     [ dimensionless] 
!      qamin          minimum permissible cloud fraction 
!                     [ dimensionless] 
!      do_brenguier   should drops at top of stratocumulus clouds be
!                     scaled?
!
!----------------------------------------------------------------------

! <NAMELIST NAME="cloud_rad_nml">
!  <DATA NAME="overlap" UNITS="" TYPE="" DIM="" DEFAULT="">
!integer variable indicating which overlap 
!                     assumption to use:
!                     overlap = 1. means condensate in adjacent levels 
!                                  is treated as part of the same cloud
!                                  i.e. maximum-random overlap
!                     overlap = 2. means condensate in adjacent levels 
!                                  is treated as different clouds
!                                  i.e. random overlap
!  </DATA>
!  <DATA NAME="l2strem" UNITS="" TYPE="" DIM="" DEFAULT="">
!logical variable indicating which solution for 
!                     cloud radiative properties is being used.
!                     l2strem = T  2 stream solution
!                     l2strem = F  Delta-Eddington solution
!                     Note that IF l2strem = T then the solution does 
!                     not depend on solar zenith angle
!  </DATA>
!  <DATA NAME="taucrit" UNITS="" TYPE="" DIM="" DEFAULT="">
! critical optical depth for switching direct beam 
!                     to diffuse beam for use in Delta-Eddington 
!                     solution [ dimensionless] 
!  </DATA>
!  <DATA NAME="adjust_top" UNITS="" TYPE="" DIM="" DEFAULT="">
!logical variable indicating whether or not to use 
!                     the code which places the top and bottom of the 
!                     cloud at the faces which are most in view from
!                     the top and bottom of the cloud block. this is 
!                     done to avoid undue influence of very small cloud
!                     fractions. if true this adjustment of tops is 
!                     performed; if false this is not performed.
!  </DATA>
!  <DATA NAME="scale_factor" UNITS="" TYPE="" DIM="" DEFAULT="">
!factor which multiplies actual cloud optical 
!                     depths to account for the plane-parallel homo-
!                     genous cloud bias  (e.g. Cahalan effect).
!                     [ dimensionless] 
!  </DATA>
!  <DATA NAME="qamin" UNITS="" TYPE="" DIM="" DEFAULT="">
!minimum permissible cloud fraction 
!                     [ dimensionless] 
!  </DATA>
!  <DATA NAME="do_brenguier" UNITS="" TYPE="" DIM="" DEFAULT="">
!should drops at top of stratocumulus clouds be
!                     scaled?
!  </DATA>
! </NAMELIST>
!

namelist /cloud_rad_nml/                                       &
                         overlap, l2strem, taucrit,     &
                         adjust_top, scale_factor, qamin, &
                         do_brenguier, N_min, Ni_min, snow_in_cloudrad, &
                         rain_in_cloudrad, min_diam_ice_mg, &
                         dcs_mg, min_diam_drop_mg, max_diam_drop_mg


!------------------------------------------------------------------
!---- public data ------


!-------------------------------------------------------------------
!---- private data ------

!-------------------------------------------------------------------
!   various physical parameters:
!-------------------------------------------------------------------
real, parameter :: taumin = 1.E-06  ! minimum permissible tau  
                                    ! [ dimensionless ]
real            :: qmin = 1.E-10    ! minimum permissible cloud 
                                    ! condensate [ kg condensate / 
                                    ! kg air ]                 
real            :: N_land = 3.E+08  ! number of cloud droplets in liquid
                                    ! clouds over land  [ m**(-3) ]
real            :: N_ocean = 1.E+08 ! number of cloud droplets in liquid
                                    ! clouds over ocean [ m**(-3) ]
real, parameter :: k_land = 1.143   ! ratio of effective radius to 
                                    ! volume radius for continental
                                    ! air masses  [ dimensionless ]
real, parameter :: k_ocean = 1.077  ! ratio of effective radius to 
                                    ! volume radius for continental
                                    ! air masses  [ dimensionless ]
!       IMPORTANT NOTE qmin, N_land, N_ocean are initialized with a 
!       call from strat_cloud_init to cloud_rad_init.  This guarantees
!       that both strat_cloud and cloud_rad have the exact same values
!       for these parameters.
 
real            :: qcvar

!----------------------------------------------------------------------
!    diagnostics variables.        
!----------------------------------------------------------------------
character(len=8)    :: mod_name = 'cloud_rad'
real                :: missing_value = -999.

integer ::            id_aice, id_reffice, id_aliq, id_reffliq, &
           id_alow, id_tauicelow, id_tauliqlow, id_tlaylow, id_tcldlow


!---------------------------------------------------------------------
!   logical variables:
!--------------------------------------------------------------------
logical   :: module_is_initialized = .false.  ! is module initialized ?
logical   :: do_liq_num          = .false.  ! use prog. droplet number ?
logical   :: do_ice_num          = .false. ! use prog ice crystal number?





!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#######################################################################


! <SUBROUTINE NAME="cloud_rad_init">
!  <OVERVIEW>
!
!   Called once to initialize cloud_rad module.   This routine reads the
!   namelist, registers any requested diagnostic fields, and (when
!   called from strat_cloud_init [standard practice]) returns the
!   overlap assumption to strat_cloud for use in determining cloud and
!   large-scale precipitation overlap. 
!
!  </OVERVIEW>
!  <DESCRIPTION>
!
!   Initializes values of qmin, N_land, and  N_ocean using values from
!   strat_cloud namelist as well as reads its own namelist variables. In
!   addition, it registers diagnostic fields if needed, and returns the 
!   value of the cloud overlap to strat_cloud.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_rad_init (axes, Time, qmin_in, N_land_in, N_ocean_in, &
!                        overlap_out)
!
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer, optional">
!    Axis integers for diagnostics
!  </IN>
!  <IN NAME="Time" TYPE="time_type, optional">
!     Time type variable for diagnostics
!  </IN>
!  <IN NAME="qmin_in" TYPE="real, kg condensate/kg air, optional">
!      Input value of minimum permissible cloud liquid, ice,
!      or fraction                
!  </IN>
!  <IN NAME="N_land_in" TYPE="real, #/(m*m*m), optional">
!    Input value of number of cloud drop per cubic meter
!    over land
!  </IN>
!  <IN NAME="N_ocean_in" TYPE="real, #/(m*m*m), optional">
!    Input value of number of cloud drop per cubic meter
!    over ocean
!  </IN>
!  <OUT NAME="overlap_out" TYPE="integer, optional">
!    Integer indicating the overlap assumption being used 
!                       (1 = maximum-random, 2 = random)
!  </OUT>
!  <ERROR MSG="" STATUS="FATAL">
!   Fatal crashes occur in initialization of the module if:
!
!   1. overlap does not equal 1 or 2
!
!   2. taucrit < 0.
!
!   3. scale_factor < 0.
!
!   4. qamin outside of the range of 0 to 1.
!  </ERROR>
! </SUBROUTINE>
!
subroutine cloud_rad_init (axes, Time, qmin_in, N_land_in, N_ocean_in, &
                           prog_droplet_in, prog_ice_num_in,  qcvar_in, &
                           overlap_out)
                               
!--------------------------------------------------------------------
!    cloud_rad_init is the constructor for cloud_rad_mod.
!--------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine initializes values of qmin, N_land, and 
!       N_ocean using values from the strat_cloud_module, 
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!VARIABLES
!
!
!       --------------
!       OPTIONAL INPUT
!       --------------
!
!
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       axes           axis integers for diagnostics
!
!       Time           time type variable for 
!                      diagnostics
!     
!       qmin_in        input value of minimum per-     kg condensate/ 
!                      missible cloud liquid, ice,     kg air
!                      or fraction                     or fraction
!
!       N_land_in      input value of number of        #/(m*m*m)
!                      of cloud drop per cubic meter
!                      over land
!
!       N_ocean_in     input value of number of        #/(m*m*m)
!                      of cloud drop per cubic meter
!                      over ocean
!
!       ---------------
!       OPTIONAL OUTPUT
!       ---------------
!
!       overlap_out    value of the namelist variable overlap
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       unit,io        namelist integers
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

integer,         intent(in), optional     :: axes(4)
type(time_type), intent(in), optional     :: Time
REAL,     INTENT (IN),  OPTIONAL          :: qmin_in,N_land_in,&
                                             N_ocean_in
LOGICAL,  INTENT (IN), OPTIONAL           :: prog_droplet_in
LOGICAL,  INTENT (IN), OPTIONAL           :: prog_ice_num_in
REAL   ,  INTENT (IN), OPTIONAL           :: qcvar_in
INTEGER,  INTENT (OUT), OPTIONAL          :: overlap_out

!  Internal variables
!  ------------------


INTEGER                                  :: unit,io,ierr, logunit

!-----------------------------------------------------------------------
!       
!       Code
!
    
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call constants_init
      call diag_manager_init

!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloud_rad_nml, iostat=io)
      ierr = check_nml_error(io,'cloud_rad_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
           read  (unit, nml=cloud_rad_nml, iostat=io, end=10)
           ierr = check_nml_error(io,'cloud_rad_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=cloud_rad_nml)

!-----------------------------------------------------------------------
!
!       Prevent unreasonable values

        if (overlap.ne.1 .and. overlap.ne.2) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                                'overlap must be either 1 or 2 ', FATAL)
        if (taucrit .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                  'taucrit must be greater than or equal to 0. ', FATAL)
        if (scale_factor .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                         'scale_factor must be greater than 0. ', FATAL)
        if (qamin .le. 0. .or. qamin .ge. 1.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                               'qamin must be between 0. and 1.', FATAL)
        
!-----------------------------------------------------------------------
!
!       Assign values

        if (present(qmin_in)) then
              qmin = qmin_in
        end if
        if (present(N_land_in)) then
              N_land = N_land_in
        end if
        if (present(N_ocean_in)) then
              N_ocean = N_ocean_in
        end if
        if (present(overlap_out)) then
              overlap_out = overlap
        end if
        if (present(prog_droplet_in)) then
              do_liq_num = prog_droplet_in
        end if
        if (present(prog_ice_num_in)) then
              do_ice_num = prog_ice_num_in
        end if
        if (present(qcvar_in)) then
              qcvar = qcvar_in
        end if
 
        call mg_const_init
        call gamma_mg_init
        
       module_is_initialized = .true.

end subroutine cloud_rad_init

!###################################################################

! <SUBROUTINE NAME="cloud_rad_end">
!  <OVERVIEW>
!
!    A destructor routine for the cloud_rad module.
!
!  </OVERVIEW>
!  <DESCRIPTION>
!
!    A destructor routine for the cloud_rad module.
!
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_rad_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cloud_rad_end

       call gamma_mg_end

       module_is_initialized = .false.

end subroutine cloud_rad_end

!#####################################################################

! <SUBROUTINE NAME="lw_emissivity">
!  <OVERVIEW>
!   
!    Subroutine lw_emissivity computes the longwave cloud emissivity 
!    using the cloud mass absorption coefficient and the water path.
!
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_emissivity (is, js, lwp, iwp, reff_liq, reff_ice,   &
!                       nclds, em_lw)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!     Starting subdomain i index of data 
!     in the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!     Starting subdomain j index of data 
!     in the physics_window being integrated
!  </IN>
!  <IN NAME="lwp" TYPE="real">
!     Liquid water path [ kg / m**2 ]
!  </IN>
!  <IN NAME="iwp" TYPE="real">
!     Ice water path [ kg / m**2 ]
!  </IN>
!  <IN NAME="reff_liq" TYPE="real">
!     Effective cloud drop radius used with
!     bulk cloud physics scheme [ microns ]
!  </IN>
!  <IN NAME="reff_ice" TYPE="real">
!     Effective ice crystal radius used with
!     bulk cloud physics scheme [ microns ]
!  </IN>
!  <IN NAME="nclds" TYPE="integer">
!     Number of random overlapping clouds in column
!  </IN>
!  <OUT NAME="em_lw" TYPE="real">
!     longwave cloud emmissivity [ dimensionless ]
!  </OUT>
! </SUBROUTINE>
!
subroutine lw_emissivity (is, js, lwp, iwp, reff_liq, reff_ice,   &
                          nclds, em_lw)

!---------------------------------------------------------------------
!    subroutine lw_emissivity computes the longwave cloud emissivity 
!    using the cloud mass absorption coefficient and the water path.
!---------------------------------------------------------------------

integer,                 intent(in)   ::  is,js
real, dimension(:,:,:),  intent(in)   ::  lwp, iwp, reff_liq, reff_ice
integer, dimension(:,:), intent(in)   ::  nclds
real, dimension(:,:,:),  intent(out)  ::  em_lw


!--------------------------------------------------------------------
!   intent(in) variables:
!
!        is,js           starting subdomain i,j indices of data 
!                        in the physics_window being integrated     
!        lwp             liquid water path [ kg / m**2 ]
!        iwp             ice water path [ kg / m**2 ]
!        reff_liq        effective cloud drop radius  used with
!                        bulk cloud physics scheme [ microns ]
!        reff_ice        effective ice crystal radius used with
!                        bulk cloud physics scheme [ microns ]
!        nclds           number of random overlapping clouds in column
!
!    intent(out) variables:
!
!        em_lw           longwave cloud emmissivity [ dimensionless ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(em_lw,1), size(em_lw,2),                 &
                                      size(em_lw,3)) ::  k_liq, k_ice

!---------------------------------------------------------------------
!   local variables:
!     
!     k_liq             liquid cloud mass absorption coefficient for 
!                       longwave portion of spectrum 
!                       [ m**2 / kg condensate ]
!     k_ice             ice cloud mass absorption coefficient for 
!                       longwave portion of spectrum 
!                       [ m**2 / kg condensate ]
!     i,j,k             do-loop indices
!
!---------------------------------------------------------------------
              
!----------------------------------------------------------------------
!    compute longwave emissivity, including contributions from both the
!    ice and liquid cloud particles present.
!----------------------------------------------------------------------
      k_liq = 140.
      k_ice = 4.83591 + 1758.511/reff_ice       
 
      em_lw = 1. - exp(-1.*( k_liq*lwp +  k_ice*iwp))

!----------------------------------------------------------------------


    
end subroutine lw_emissivity                   




!######################################################################

! <SUBROUTINE NAME="cloud_summary3">
!  <OVERVIEW>
!
!   cloud_summary3 returns the specification properties of the clouds
!    present in the strat_cloud_mod.
!
!  </OVERVIEW>
!  <DESCRIPTION>
!
!   cloud_summary3 returns the specification properties of the clouds
!    present in the strat_cloud_mod.
!
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_summary3 (is, js, land, ql, qi, qa, qn, pfull, phalf, &
!                        tkel, nclds, cldamt, lwp, iwp, reff_liq,  &
!                        reff_ice, ktop, kbot, conc_drop, conc_ice, &
!                        size_drop, size_ice)
!
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!    Indices for model slab
!  </IN>
!  <IN NAME="land" TYPE="real">
!    Fraction of the grid box covered by land
!                    [ dimensionless ]
!  </IN>
!  <IN NAME="ql" TYPE="real">
!    Cloud liquid condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qi" TYPE="real">
!    Cloud ice condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qa" TYPE="real">
!    Cloud volume fraction [ fraction ]
!  </IN>
!  <IN NAME="qn" TYPE="real">
!    Cloud droplet number [ #/kg air ]
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!    Pressure at full levels [ Pascals ]
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!    Pressure at half levels [ Pascals ]
!    NOTE: it is assumed that phalf(j+1) > phalf(j)
!  </IN>
!  <IN NAME="tkel" TYPE="real">
!    Temperature [ deg. Kelvin ]
!  </IN>
!  <OUT NAME="nclds" TYPE="integer">
!    Number of random-overlap clouds in a column
!  </OUT>
!  <OUT NAME="cldamt" TYPE="real">
!    Cloud amount of condensed cloud
!  </OUT>
!  <OUT NAME="lwp" TYPE="real">
!    Liquid water path
!  </OUT>
!  <OUT NAME="iwp" TYPE="real">
!    Ice water path
!  </OUT>
!  <OUT NAME="reff_liq" TYPE="real">
!    Effective radius of cloud drops
!  </OUT>
!  <OUT NAME="reff_ice" TYPE="real">
!    Effective radius of ice crystals
!  </OUT>
!  <OUT NAME="ktop" TYPE="integer, optional">
!    Integer level for top of cloud, present when 
!    max-random overlap assumption made.
!  </OUT>
!  <OUT NAME="kbot" TYPE="integer, optional">
!    Integer level for bottom of cloud, present when
!    max-random overlap assumption made.
!  </OUT>
!  <OUT NAME="conc_drop" TYPE="real, optional">
!    Liquid cloud droplet mass concentration, present 
!    when microphysically-based cloud radiative
!    properties are desired.
!  </OUT>
!  <OUT NAME="conc_ice" TYPE="real, optional">
!    Ice cloud mass concentration, present when
!    microphysically-based cloud radiative
!    properties are desired
!  </OUT>
!  <OUT NAME="size_drop" TYPE="real, optional">
!    Effective diameter of liquid cloud droplets, 
!    present when microphysically-based cloud radiative
!    properties are desired.
!  </OUT>
!  <OUT NAME="size_ice" TYPE="real, optional">
!     Effective diameter of ice cloud, present when 
!     microphysically-based cloud radiative
!     properties are desired.
!  </OUT>
! </SUBROUTINE>
!
subroutine cloud_summary3 (is, js, land,  use_fu2007, ql, qi, qa, qn, &
                           qni, pfull, phalf, &
                           tkel, nclds, cldamt, lwp, iwp, reff_liq,  &
                           reff_ice, ktop, kbot, conc_drop, conc_ice, &
                           size_drop, size_ice, droplet_number,  &
                           ice_number)
   
!---------------------------------------------------------------------
!    cloud_summary3 returns the specification properties of the clouds
!    present in the strat_cloud_mod.
!---------------------------------------------------------------------
 
integer,                   intent(in)            :: is,js
real, dimension(:,:),      intent(in)            :: land
logical,                   intent(in)             :: use_fu2007
real, dimension(:,:,:),    intent(in)            :: ql, qi, qa, qn, pfull,&
                                                    phalf, tkel, qni
integer, dimension(:,:),   intent(out)           :: nclds          
real, dimension(:,:,:),    intent(out)           :: cldamt, lwp, iwp, &
                                                    reff_liq, reff_ice
integer, dimension(:,:,:), intent(out), optional :: ktop, kbot 
real,    dimension(:,:,:), intent(out), optional :: conc_drop,conc_ice,&
                                                    size_drop,size_ice,&
                                                    droplet_number, &
                                                    ice_number

!---------------------------------------------------------------------
!    intent(in) variables:
!
!       is,js        Indices for model slab
!       land         Fraction of the grid box covered by land
!                    [ dimensionless ]
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       qn           Cloud droplet number [ #/kg air]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels [ Pascals ]
!                    NOTE: it is assumed that phalf(j+1) > phalf(j)
!       tkel         Temperature [ deg. Kelvin ] 
!
!    intent(out) variables:
!
!       nclds        Number of random-overlap clouds in a column
!       cldamt       Cloud amount of condensed cloud
!       lwp          Liquid water path 
!       iwp          Ice water path
!       reff_liq     Effective radius of cloud drops
!       reff_ice     Effective radius of ice crystals
!
!   intent(out), optional variables:
! 
!       ktop         Integer level for top of cloud, present when 
!                    max-random overlap assumption made
!       kbot         Integer level for bottom of cloud, present when
!                    max-random overlap assumption made
!       conc_drop    Liquid cloud droplet mass concentration, present 
!                    when microphysically-based cloud radiative
!                    properties are desired
!       conc_ice     Ice cloud mass concentration, present when
!                    microphysically-based cloud radiative
!                    properties are desired
!       size_drop    Effective diameter of liquid cloud droplets, 
!                    present when microphysically-based cloud radiative
!                    properties are desired
!       size_ice     Effective diameter of ice cloud, present when 
!                    microphysically-based cloud radiative
!                    properties are desired
!       droplet_number
!                    number of cloud droplets [ # / kg(air) ]
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    local variables:

      real,dimension (size(ql,1),size(ql,2),   &
                                 size(ql,3)) :: qa_local, ql_local, &
                                                qi_local, N_drop3D, &
                                                N_ice3D

      real,dimension (size(ql,1),size(ql,2)) :: N_drop2D, k_ratio
      integer  :: i, j, k

!--------------------------------------------------------------------
!    local variables:
!
!       qa_local     local value of qa (fraction)
!       ql_local     local value of ql (kg condensate / kg air)
!       qi_local     local value of qi (kg condensate / kg air)
!       N_drop[23]D  number of cloud droplets per cubic meter
!       k_ratio      ratio of effective radius to mean volume radius
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    create local values of ql and qi. this step is necessary to remove 
!    the values of (qi,ql) which are 0 < (qi,ql) < qmin   or
!    (qi,ql) > qmin and qa <= qamin.
!--------------------------------------------------------------------

      do k=1, size(ql,3)
        do j=1, size(ql,2)
          do i=1, size(ql,1)
            qa_local(i,j,k) = 0.
            if ((qa(i,j,k) > qamin) .and. (ql(i,j,k) > qmin) ) then
              ql_local(i,j,k) = ql(i,j,k)
              qa_local(i,j,k) = qa(i,j,k)
            else
              ql_local(i,j,k) = 0.
            endif      
            if ((qa(i,j,k) > qamin) .and. (qi(i,j,k) > qmin) ) then
              qi_local(i,j,k) = qi(i,j,k)
              qa_local(i,j,k) = qa(i,j,k)
            else
              qi_local(i,j,k) = 0.
            endif       
          end do
        end do
      end do

!--------------------------------------------------------------------
!    define the cloud droplet concentration and the ratio of the 
!    effective drop radius to the mean volume radius.
!--------------------------------------------------------------------
      N_drop2D(:,:)  = N_land*land(:,:) + N_ocean*(1. - land(:,:))
      k_ratio(:,:) = k_land*land(:,:) + k_ocean*(1. - land(:,:))
!yim prognostic droplet number
      if (do_liq_num) then
        N_drop3D=qn
        droplet_number = qn
      else 
        do k=1, size(ql,3)
          do j=1, size(ql,2)
            do i=1, size(ql,1)
              droplet_number(i,j,k) = N_drop2D(i,j)/(pfull(i,j,k)/  &
                                      (RDGAS*tkel(i,j,k)))
            end do
          end do
        end do
      endif    

      if ( do_ice_num ) then
        ice_number = qni
        N_ice3D=qni
      else
        ice_number=-999.
        N_ice3D=-999.
      end if

!--------------------------------------------------------------------
!    execute the following when  the max-random overlap assumption 
!    is being made. 
!--------------------------------------------------------------------
      if (present(ktop) .and. present(kbot)) then    ! max-rnd

!--------------------------------------------------------------------
!    if microphysics output is required, only the random overlap assump-
!    tion is allowed; if max-random overlap is requested, an error
!    message will be issued. if random overlap is requested, call
!    subroutine rnd_overlap to obtain the cloud specification proper-
!    ties, including the microphysical parameters.
!--------------------------------------------------------------------
        if (present (conc_drop) .and.  present (conc_ice ) .and. &
            present (size_ice ) .and.  present (size_drop)) then      
          call error_mesg ( 'cloud_rad_mod', &
       ' max-random overlap not currently available for radiation '//&
              'scheme requiring microphysically-based outputs', FATAL)
     
!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
        else if (present (conc_drop) .or.  present (conc_ice ) .or. &
                 present (size_ice ) .or.  present (size_drop)) then
          call error_mesg ('cloud_rad_mod', &
                ' if any microphysical args present, all must be '//&
                                                    'present', FATAL)

        else
          call  max_rnd_overlap (ql_local, qi_local, qa_local, pfull,  &
                                 phalf, tkel, N_drop3D, N_drop2D, k_ratio, nclds,  &
                                 ktop, kbot, cldamt, lwp, iwp,   &
                                 reff_liq, reff_ice, N_ice3D)
        endif
     
!---------------------------------------------------------------------
!    if only ktop or kbot is present, stop execution; both are needed
!    for max-random overlap and neither are prrmitted when the 
!    random overlap assumption is made.
!---------------------------------------------------------------------
      else if (present(ktop) .or. present(kbot)) then ! error
        call error_mesg ('cloud_rad_mod',  &
                  'kbot and ktop must either both be absent or both '//&
                    'be present', FATAL)

!---------------------------------------------------------------------
!    if neither are present, then random overlap is assumed.
!---------------------------------------------------------------------
      else                 

!---------------------------------------------------------------------
!    if microphysical properties are desired, call subroutine 
!    rnd_overlap to obtain the cloud specification properties, including
!    the microphysical parameters.
!--------------------------------------------------------------------
        if (present (conc_drop) .and.  present (conc_ice ) .and. &
            present (size_ice ) .and.  present (size_drop)) then      
          call rnd_overlap (ql_local, qi_local, qa_local,  &
                            use_fu2007, pfull, phalf, tkel,  &
                            N_drop3D, N_drop2D, N_ice3D, k_ratio, nclds, &
                            cldamt, lwp, iwp, reff_liq, reff_ice,   &
                            conc_drop_org=conc_drop,&
                            conc_ice_org =conc_ice,&
                            size_drop_org=size_drop,&
                            size_ice_org =size_ice)

!--------------------------------------------------------------------
!    account for the plane-parallel homogeneous cloud bias.
!--------------------------------------------------------------------
          conc_drop = scale_factor*conc_drop
          conc_ice  = scale_factor*conc_ice 

!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
        else if (present (conc_drop) .or.  present (conc_ice ) .or. &
                 present (size_ice ) .or.  present (size_drop)) then   
          call error_mesg ('cloud_rad_mod', &
                ' if any microphysical args present, all must '//&
                                                'be present', FATAL)

!----------------------------------------------------------------------
!    if microphysics terms are not required, call rnd_overlap to obtain
!    the cloud specification variables.
!----------------------------------------------------------------------
        else
           call  rnd_overlap (ql_local, qi_local, qa_local,  &
                              use_fu2007, pfull, phalf, tkel,  &
                              N_drop3D, N_drop2D, N_ice3D, &
                              k_ratio, nclds,  &
                              cldamt, lwp, iwp, reff_liq, reff_ice)
        endif
      endif ! (present(ktop and kbot))

!---------------------------------------------------------------------
    


end subroutine cloud_summary3



!######################################################################

! <SUBROUTINE NAME="max_rnd_overlap">
!  <OVERVIEW>
!
!    max_rnd_overlap returns various cloud specification properties
!    obtained with the maximum-random overlap assumption.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!
!    max_rnd_overlap returns various cloud specification properties
!    obtained with the maximum-random overlap assumption.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call max_rnd_overlap (ql, qi, qa, pfull, phalf, tkel, N_drop3D, N_drop2D,  &
!                         k_ratio, nclds, ktop, kbot, cldamt, lwp,  &
!                         iwp, reff_liq, reff_ice)
!
!  </TEMPLATE>
!  <IN NAME="ql" TYPE="real">
!    Cloud liquid condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qi" TYPE="real">
!    Cloud ice condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qa" TYPE="real">
!    Cloud volume fraction [ fraction ]
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!    Pressure at full levels [ Pascals ]
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!    Pressure at half levels, index 1 at model top 
!    [ Pascals ]
!  </IN>
!  <IN NAME="tkel" TYPE="real">
!    Temperature [ deg Kelvin ]
!  </IN>
!  <IN NAME="N_drop[23]D" TYPE="real">
!    Number of cloud droplets per cubic meter (2 and 3 dimensional array)
!  </IN>
!  <IN NAME="k_ratio" TYPE="real">
!    Ratio of effective radius to mean volume radius
!  </IN>
!  <OUT NAME="nclds" TYPE="integer">
!    Number of (random overlapping) clouds in column 
!  </OUT>
!  <OUT NAME="ktop" TYPE="integer">
!    Level of the top of the cloud.
!  </OUT>
!  <OUT NAME="kbot" TYPE="integer">
!    Level of the bottom of the cloud.
!  </OUT>
!  <OUT NAME="cldamt" TYPE="real">
!    Cloud amount of condensed cloud [ dimensionless ]
!  </OUT>
!  <OUT NAME="lwp" TYPE="real">
!    Cloud liquid water path [ kg condensate / m **2 ]
!  </OUT>
!  <OUT NAME="iwp" TYPE="real">
!    Cloud ice path [ kg condensate / m **2 ]
!  </OUT>
!  <OUT NAME="reff_liq" TYPE="real">
!    Effective radius for liquid clouds [ microns ]
!  </OUT>
!  <OUT NAME="reff_ice" TYPE="real">
!    Effective particle size for ice clouds [ microns ]
!  </OUT>
! </SUBROUTINE>
!
subroutine max_rnd_overlap (ql, qi, qa, pfull, phalf, tkel, N_drop3D, N_drop2D,  &
                           k_ratio, nclds, ktop, kbot, cldamt, lwp,  &
                           iwp, reff_liq, reff_ice, N_ice3D)

!----------------------------------------------------------------------
!    max_rnd_overlap returns various cloud specification properties
!    obtained with the maximum-random overlap assumption.
!----------------------------------------------------------------------
 
real,    dimension(:,:,:), intent(in)             :: ql, qi, qa,  &
                                                     pfull, phalf, tkel, &
                                                     N_drop3D, N_ice3D
real,    dimension(:,:),   intent(in)             :: N_drop2D, k_ratio
integer, dimension(:,:),   intent(out)            :: nclds
integer, dimension(:,:,:), intent(out)            :: ktop, kbot
real,    dimension(:,:,:), intent(out)            :: cldamt, lwp, iwp, &
                                                     reff_liq, reff_ice

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels, index 1 at model top 
!                    [ Pascals ]
!       tkel         Temperature [ deg Kelvin ]
!       N_drop       Number of cloud droplets per cubic meter
!       k_ratio      Ratio of effective radius to mean volume radius
!
!   intent(out) variables:
!
!       nclds        Number of (random overlapping) clouds in column 
!       ktop         Level of the top of the cloud
!       kbot         Level of the bottom of the cloud
!       cldamt       Cloud amount of condensed cloud [ dimensionless ]
!       lwp          Cloud liquid water path [ kg condensate / m **2 ]
!       iwp          Cloud ice path [ kg condensate / m **2 ]
!       reff_liq     Effective radius for liquid clouds [ microns ]
!       reff_ice     Effective particle size for ice clouds [ microns ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (size(ql,1), size(ql,2), size(ql,3))  :: &
                    cldamt_cs, lwp_cs, iwp_cs, reff_liq_cs, reff_ice_cs 

      integer    :: kdim
      integer    :: top_t, bot_t
      integer    :: tmp_top, tmp_bot, nlev
      logical    :: already_in_cloud, cloud_bottom_reached
      real       :: sum_liq, sum_ice, maxcldfrac
      real       :: totcld_bot, max_bot
      real       :: totcld_top, max_top, tmp_val
      real       :: reff_liq_local, sum_reff_liq
      real       :: reff_ice_local, sum_reff_ice
      integer    :: i, j, k, kc, t

      real       :: dumc, dumnc, rho, pgam, lamc, lammax, lammin, &
                    dumi, dumni, lami
!--------------------------------------------------------------------
!   local variables:
!
!       kdim              number of model layers
!       top_t             used temporarily as tag for cloud top index
!       bot_t             used temporarily as tag for cloud bottom index
!       tmp_top           used temporarily as tag for cloud top index
!       tmp_bot           used temporarily as tag for cloud bottom index
!       nlev              number of levels in the cloud
!       already_in_cloud  if true, previous layer contained cloud
!       cloud_bottom_reached
!                         if true, the cloud-free layer beneath a cloud
!                         has been reached
!       sum_liq           sum of liquid in cloud 
!                         [ kg condensate / m**2 ]
!       sum_ice           sum of ice in cloud 
!                         [ kg condensate / m**2 ]
!       maxcldfrac        maximum cloud fraction in any layer of cloud
!                         [ fraction ]
!       totcld_bot        total cloud fraction from bottom view
!       max_bot           largest cloud fraction face from bottom view
!       totcld_top        total cloud fraction from top view
!       max_top           largest cloud fraction face from top view
!       tmp_val           temporary number used in the assigning of top 
!                         and bottom
!       reff_liq_local    gridpoint value of reff of liquid clouds 
!                         [ microns ]
!       sum_reff_liq      condensate-weighted sum over cloud of 
!                         reff_liq_local  
!                         [ (kg condensate / m**2) * microns ]
!       reff_ice_local    gridpoint value ofreff of ice clouds  
!                         [ microns ]
!       sum_reff_ice      condensate-weighted sum over cloud of 
!                         reff_ice_local 
!                         [ (kg condensate / m**2) * microns ]
!       i,j,k,kc,t        do-loop indices
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the number of vertical layers in the model. initialize the
!    output fields to correspond to the absence of clouds.
!---------------------------------------------------------------------
      kdim     = size(ql,3)
      nclds    = 0
      ktop     = 1
      kbot     = 0
      cldamt   = 0.
      lwp      = 0.
      iwp      = 0.
      reff_liq = 10.
      reff_ice = 30.

!--------------------------------------------------------------------
!    find the levels with cloud in each column. determine the vertical
!    extent of each individual cloud, treating cloud in adjacent layers
!    as components of a multi-layer cloud, and then calculate appropr-
!    iate values of water paths and effective particle size.
!--------------------------------------------------------------------

      do j=1,size(ql,2)
        do i=1,size(ql,1)

!--------------------------------------------------------------------
!    set a flag indicating that we are searching for the next cloud top.
!--------------------------------------------------------------------
          already_in_cloud  = .false.
          cloud_bottom_reached = .false.

!--------------------------------------------------------------------
!    march down the column.
!--------------------------------------------------------------------
          do k=1,kdim      

!--------------------------------------------------------------------
!    find a layer containing cloud in the column. 
!--------------------------------------------------------------------
            if ( (ql(i,j,k) .gt. qmin) .or. &
                 (qi(i,j,k) .gt. qmin) ) then      

!--------------------------------------------------------------------
!    if the previous layer was not cloudy, then a new cloud has been
!    found. increment the cloud counter, set the flag to indicate the 
!    layer is in a cloud, save its cloud top level, initialize the 
!    values of its ice and liquid contents and fractional area and 
!    effective crystal and drop sizes. 
!--------------------------------------------------------------------
              if (.not. already_in_cloud)  then
                nclds(i,j) = nclds(i,j) + 1
                already_in_cloud = .true.
                cloud_bottom_reached = .false.
                ktop(i,j,nclds(i,j)) = k
                sum_liq          = 0.
                sum_ice          = 0.
                maxcldfrac       = 0.
                sum_reff_liq     = 0.
                sum_reff_ice     = 0.        
              endif

!--------------------------------------------------------------------
!    if liquid is present in the layer, compute the effective drop
!    radius. the following formula, recommended by (Martin et al., 
!    J. Atmos. Sci, vol 51, pp. 1823-1842) is used for liquid droplets:
!    reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = density of cloud droplets (number per cubic meter)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------
do_liq_num_if: if(.not. do_liq_num) then
               if (ql(i,j,k) > qmin) then
                  reff_liq_local = k_ratio(i,j)*620350.49*    &
                                   (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/  &
                                   RDGAS/tkel(i,j,k)/DENS_H2O/  &
                                   N_drop2D(i,j))**(1./3.)
               else
                 reff_liq_local = 0.
               endif
             else ! do_liq_num_if
!cms++
! BETTER SPLIT LOOP AND MOVE MOVE IFs OUTSIDE  and BETTER USE ONE SUBROUTINE
!    INTEAD OF DUPLICATING CODE!!

 mg_if: IF ( .NOT. do_ice_num ) THEN
!cms--

!--------------------------------------------------------------------
! yim: a variant for prognostic droplet number
!    reff (in microns) =  k * 1.E+06 *
!                    (3*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = mixing ratio of cloud droplets (number/kg air)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------
               if (ql(i,j,k) > qmin) then
                 reff_liq_local = k_ratio(i,j)*620350.49*    &
                                  (ql(i,j,k)/DENS_H2O/  &
                                  max(N_drop3D(i,j,k),   &
                                      N_min*max(qa(i,j,k),qmin)/  &
                             (pfull(i,j,k)/RDGAS/tkel(i,j,k))))**(1./3.)
               else
                 reff_liq_local = 0.
               endif

!cms++
ELSE !mg_if

!-----------------------------------------------------------
!
! cloud droplet effective radius
!

! in-cloud mixing ratio and number conc
         dumc = ql(i,j,k)/qa(i,j,k)
         dumnc = N_drop3D(i,j,k)/qa(i,j,k)

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

         dumc =min(dumc,5.e-3)

!cms++2008-11-26

         dumnc =  max(dumnc, N_min)

!cms--
         if ( dumc > qmin ) then


! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc = min(dumnc,dumc*1.e20)

         rho=pfull(i,j,k)/(RDGAS*tkel(i,j,k))

         pgam=0.0005714*(dumnc/1.e6/rho)+0.2714
         pgam=1./(pgam**2)-1.
         pgam=max(pgam,2.)
         pgam=min(pgam,15.)

         lamc = (pi/6.*rhow*dumnc*gamma_mg(pgam+4.)/ &
                 (dumc*gamma_mg(pgam+1.)))**(1./3.)
         lammin = (pgam+1.)/max_diam_drop_mg
         lammax = (pgam+1.)/min_diam_drop_mg
        if (lamc.lt.lammin) then
         lamc = lammin
        else if (lamc.gt.lammax) then
         lamc = lammax
        end if
        reff_liq_local = gamma_mg(qcvar+1./3.)/(gamma_mg(qcvar)*  &
                             qcvar**(1./3.))*gamma_mg(pgam+4.)/ & 
                                          gamma_mg(pgam+3.)/lamc/2.*1.e6
        else
        reff_liq_local = 10.
        end if
        
!ELSE!! qamin_if1
!  reff_ice_local = 0.
!  reff_liq_local = 0.
!      
!END IF   qamin_if1      
END IF mg_if
 endif    do_liq_num_if

!cms--

!----------------------------------------------------------------------
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly 
!    stratified cloud in liquid specific humidity, the cloud top 
!    effective radius is greater than the effective radius of the 
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!    single layer liquid or mixed phase clouds.
!
!---------------------------------------------------------------------- 
              if (do_brenguier) then
              if ( k == 1 ) then
                 if (qa(i,j,2) < qamin) then
                   reff_liq_local = 1.134*reff_liq_local
                 endif
              else if (k == kdim ) then
                 if ( qa(i,j,kdim-1) < qamin) then
                   reff_liq_local = 1.134*reff_liq_local
                 endif
              else if (qa(i,j,k-1) .lt. qamin .and. & 
                       qa(i,j,k+1) .lt. qamin)  then
                reff_liq_local = 1.134*reff_liq_local
              end if
              end if

!cms++

 IF ( do_ice_num ) THEN

! calc. based on Morrison Gettelman code
!  

!qamin_if1: IF ( qa(i,j,k) .GT. qamin ) THEN 

! in-cloud mixing ratio and number conc
         dumi = qi(i,j,k)/qa(i,j,k)
         dumni = N_ice3D(i,j,k)/qa(i,j,k)

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
         dumi =min(dumi ,5.e-3)

!cms++2009-02-19

 dumni =  max(dumni, Ni_min)

!cms--
!...................
! cloud ice effective radius
       if ( dumi > qmin ) then

! add upper limit to in-cloud number concentration to prevent numerical error
       dumni=min(dumni,dumi*1.e20)
       lami = (gamma_mg(1.+di_mg)*ci_mg* &
              dumni/dumi)**(1./di_mg)
!2008-11-18   lammax = 1./10.e-6
        lammax = 1./min_diam_ice_mg
        lammin = 1./(2.*dcs_mg)

        if (lami.lt.lammin) then
        lami = lammin
        else if (lami.gt.lammax) then
        lami = lammax
        end if
        reff_ice_local = 1.5/lami*1.e6
        else
          reff_ice_local = 25.
        end if

!cms++2008-10-31 

       reff_ice_local = 2. *  reff_ice_local !3rd moment/2nd moment as in Fu and Liou paper
!cms--

ELSE
!--------------------------------------------------------------------
!    if ice crystals are present, define their effective size, which
!    is a function of temperature. for ice clouds the effective radius
!    is taken from the formulation in Donner (1997, J. Geophys. Res., 
!    102, pp. 21745-21768) which is based on Heymsfield and Platt (1984)
!    with enhancement for particles smaller than 20 microns.  
!
!              T Range (K)               Reff (microns) 
!     -------------------------------    --------------
!
!     tfreeze-25. < T                       92.46298       
!     tfreeze-30. < T <= Tfreeze-25.        72.35392     
!     tfreeze-35. < T <= Tfreeze-30.        85.19071         
!     tfreeze-40. < T <= Tfreeze-35.        55.65818        
!     tfreeze-45. < T <= Tfreeze-40.        35.29989       
!     tfreeze-50. < T <= Tfreeze-45.        32.89967     
!     Tfreeze-55  < T <= Tfreeze-50         16.60895      
!                   T <= Tfreeze-55.        15.41627    
!
!--------------------------------------------------------------------
              if (qi(i,j,k) > qmin) then
                if (tkel(i,j,k) > TFREEZE - 25. ) then
                  reff_ice_local = 92.46298
                else if (tkel(i,j,k) >  TFREEZE - 30. .and. &
                         tkel(i,j,k) <= TFREEZE - 25.) then
                  reff_ice_local = 72.35392
                else if (tkel(i,j,k) >  TFREEZE - 35. .and. &
                         tkel(i,j,k) <= TFREEZE - 30.) then
                  reff_ice_local = 85.19071 
                else if (tkel(i,j,k) >  TFREEZE - 40. .and. &
                         tkel(i,j,k) <= TFREEZE - 35.) then
                  reff_ice_local = 55.65818
                else if (tkel(i,j,k) >  TFREEZE - 45. .and. &
                         tkel(i,j,k) <= TFREEZE - 40.) then
                  reff_ice_local = 35.29989
                else if (tkel(i,j,k) >  TFREEZE - 50. .and. &
                         tkel(i,j,k) <= TFREEZE - 45.) then
                  reff_ice_local = 32.89967
                else if (tkel(i,j,k) >  TFREEZE - 55. .and. &
                         tkel(i,j,k) <= TFREEZE - 50.) then
                  reff_ice_local = 16.60895
                else
                  reff_ice_local = 15.41627
                end if

              else
                reff_ice_local = 0.
              end if  

END IF
!---------------------------------------------------------------------
!    add this layer's contributions to the current cloud. total liquid
!    content, ice content, largest cloud fraction and condensate-
!    weighted effective droplet and crystal radii are accumulated over 
!    the cloud.
!---------------------------------------------------------------------
              sum_liq = sum_liq + ql(i,j,k)*  &
                        (phalf(i,j,k+1) - phalf(i,j,k))/GRAV
              sum_ice = sum_ice + qi(i,j,k)* &
                        (phalf(i,j,k+1) - phalf(i,j,k))/GRAV
              maxcldfrac = MAX(maxcldfrac,qa(i,j,k))
              sum_reff_liq  = sum_reff_liq + (reff_liq_local*ql(i,j,k)*&
                              (phalf(i,j,k+1) - phalf(i,j,k))/GRAV)
              sum_reff_ice  = sum_reff_ice + &
                              (reff_ice_local * qi(i,j,k) * &
                              (phalf(i,j,k+1) - phalf(i,j,k))/GRAV)
            endif ! (ql > qmin or qi > qmin)

!--------------------------------------------------------------------
!    when the cloud-free layer below a cloud is reached, or if the
!    bottom model level is reached, define the cloud bottom level and
!    set a flag indicating that mean values for the cloud may now be
!    calculated.
!--------------------------------------------------------------------
            if (ql(i,j,k) <= qmin .and. qi(i,j,k) <= qmin .and. &
                already_in_cloud) then                 
              cloud_bottom_reached = .true.
              kbot(i,j,nclds(i,j)) = k - 1
            else if (already_in_cloud .and. k == kdim) then
              cloud_bottom_reached = .true.
              kbot(i,j,nclds(i,j)) = kdim
            endif

!--------------------------------------------------------------------
!    define the cloud fraction as the largest value of any layer in the
!    cloud. define the water paths as the total liquid normalized by the
!    fractional area of the cloud. define the condensate-weighted 
!    effective water and ice radii. 
!--------------------------------------------------------------------
            if (cloud_bottom_reached) then
              cldamt_cs(i,j,nclds(i,j)) = maxcldfrac
              lwp_cs(i,j,nclds(i,j)) = sum_liq/cldamt_cs(i,j,nclds(i,j))
              iwp_cs(i,j,nclds(i,j)) = sum_ice/cldamt_cs(i,j,nclds(i,j))
              if (sum_liq > 0.) then
                reff_liq_cs(i,j,nclds(i,j)) = sum_reff_liq/sum_liq
              else
                reff_liq_cs(i,j,nclds(i,j)) = 10.0
              end if
              if (sum_ice > 0.) then
                reff_ice_cs(i,j,nclds(i,j)) = sum_reff_ice/sum_ice
              else
                reff_ice_cs(i,j,nclds(i,j)) = 30.0
              end if

!----------------------------------------------------------------------
!    if adjust_top is true, the top and bottom indices of multi-layer
!    clouds are adjusted to be those that are the most exposed to top 
!    and bottom view.
!----------------------------------------------------------------------
              if (adjust_top) then
    
!---------------------------------------------------------------------
!    define the cloud thickness.
!---------------------------------------------------------------------
                nlev = kbot(i,j,nclds(i,j)) - ktop(i,j,nclds(i,j)) + 1
                if (nlev > 1) then

!---------------------------------------------------------------------
!    use the current top and bottom as the first guess for the new 
!    values.
!---------------------------------------------------------------------
                  tmp_top = ktop(i,j,nclds(i,j))
                  tmp_bot = kbot(i,j,nclds(i,j))

!--------------------------------------------------------------------
!    initialize local search variables.
!--------------------------------------------------------------------
                  totcld_bot = 0.
                  totcld_top = 0.
                  max_bot    = 0.
                  max_top    = 0.
          
!--------------------------------------------------------------------
!    to find the adjusted cloud top, begin at current top and work 
!    downward. find the layer which is most exposed when viewed from
!    the top; i.e., the cloud fraction increase is largest for that
!    layer. the adjusted cloud base is found equivalently, starting
!    from the actual cloud base and working upwards.
!--------------------------------------------------------------------
                  do t=1,nlev

!--------------------------------------------------------------------
!    find adjusted cloud top.
!--------------------------------------------------------------------
                    top_t   = ktop(i,j,nclds(i,j)) + t - 1
                    tmp_val = MAX(0., qa(i,j,top_t) - totcld_top)
                    if (tmp_val > max_top) then
                      max_top = tmp_val
                      tmp_top = top_t
                    end if
                    totcld_top = totcld_top + tmp_val         
                              
!--------------------------------------------------------------------
!    find adjusted cloud base.
!--------------------------------------------------------------------
                    bot_t   = kbot(i,j,nclds(i,j)) - t + 1
                    tmp_val = MAX(0., qa(i,j,bot_t) - totcld_bot)
                    if (tmp_val > max_bot) then
                      max_bot = tmp_val
                      tmp_bot = bot_t
                    end if
                    totcld_bot = totcld_bot + tmp_val         
                  end do
                       
!--------------------------------------------------------------------
!    assign tmp_top and tmp_bot as the new ktop and kbot.
!--------------------------------------------------------------------
                  ktop(i,j,nclds(i,j)) = tmp_top
                  kbot(i,j,nclds(i,j)) = tmp_bot
                endif  !(nlev > 1)  
              endif  ! (adjust_top)

!---------------------------------------------------------------------
!    reset already_in_cloud and cloud_bottom_reached to indicate that
!    the current cloud has been exited.
!---------------------------------------------------------------------
              already_in_cloud     = .false.
              cloud_bottom_reached = .false.
            endif   ! (cloud_bottom_reached)
          end do
        end do
      end do

!---------------------------------------------------------------------
!    place cloud properties into physical-space arrays for return to
!    calling routine. NOTE THAT ALL LEVELS IN A GIVEN CLOUD ARE
!    ASSIGNED THE SAME PROPERTIES.
!---------------------------------------------------------------------
      do j=1,size(ql,2)
        do i=1,size(ql,1)
          do kc=1, nclds(i,j)
            do k= ktop(i,j,kc), kbot(i,j,kc)
              cldamt(i,j,k)   = cldamt_cs(i,j,kc)
              lwp(i,j,k)      = lwp_cs(i,j,kc)
              iwp(i,j,k)      = iwp_cs(i,j,kc)
              reff_liq(i,j,k) = reff_liq_cs(i,j,kc)
              reff_ice(i,j,k) = reff_ice_cs(i,j,kc)
            end do
          end do
        end do
      end do
     
!---------------------------------------------------------------------

end subroutine max_rnd_overlap




!#####################################################################

! <SUBROUTINE NAME="rnd_overlap">
!  <OVERVIEW>
!    rnd_overlap returns various cloud specification properties, 
!    obtained with the random-overlap assumption. implicit in this
!    assumption is that all clouds are only a single layer thick; i.e.,
!    clouds at adjacent levels in the same column are independent of
!    one another.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    rnd_overlap returns various cloud specification properties, 
!    obtained with the random-overlap assumption. implicit in this
!    assumption is that all clouds are only a single layer thick; i.e.,
!    clouds at adjacent levels in the same column are independent of
!    one another.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rnd_overlap    (ql, qi, qa, pfull, phalf, tkel, N_drop,  &
!                        k_ratio, nclds, cldamt, lwp, iwp, reff_liq, &
!                        reff_ice, conc_drop_org, conc_ice_org,  &
!                        size_drop_org, size_ice_org)
!
!  </TEMPLATE>
!  <IN NAME="ql" TYPE="real">
!    Cloud liquid condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qi" TYPE="real">
!    Cloud ice condensate [ kg condensate/kg air ]
!  </IN>
!  <IN NAME="qa" TYPE="real">
!    Cloud volume fraction [ fraction ]
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!    Pressure at full levels [ Pascals ]
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!    Pressure at half levels, index 1 at model top 
!    [ Pascals ]
!  </IN>
!  <IN NAME="tkel" TYPE="real">
!    Temperature [ deg Kelvin ]
!  </IN>
!  <IN NAME="N_drop" TYPE="real">
!    Number of cloud droplets per cubic meter
!  </IN>
!  <IN NAME="k_ratio" TYPE="real">
!    Ratio of effective radius to mean volume radius
!  </IN>
!  <OUT NAME="nclds" TYPE="integer">
!    Number of (random overlapping) clouds in column 
!  </OUT>
!  <OUT NAME="cldamt" TYPE="real">
!    Cloud amount of condensed cloud [ dimensionless ]
!  </OUT>
!  <OUT NAME="lwp" TYPE="real">
!    Cloud liquid water path [ kg condensate / m **2 ]
!  </OUT>
!  <OUT NAME="iwp" TYPE="real">
!    Cloud ice path [ kg condensate / m **2 ]
!  </OUT>
!  <OUT NAME="reff_liq" TYPE="real">
!    Effective radius for liquid clouds [ microns ]
!  </OUT>
!  <OUT NAME="reff_ice" TYPE="real">
!    Effective particle size for ice clouds [ microns ]
!  </OUT>
!  <OUT NAME="conc_drop_org" TYPE="real, optional">
!    Liquid cloud droplet mass concentration 
!    [ g / m**3 ]
!  </OUT>
!  <OUT NAME="conc_ice_org" TYPE="real, optional">
!    Ice cloud mass concentration [ g / m**3 ]
!  </OUT>
!  <OUT NAME="size_drop_org" TYPE="real, optional">
!    Effective diameter of liquid cloud droplets 
!    [ microns ]
!  </OUT>
!  <OUT NAME="size_ice_org" TYPE="real, optional">
!    Effective diameter of ice clouds [ microns ]
!  </OUT>
! </SUBROUTINE>
!
subroutine rnd_overlap    (ql, qi, qa, use_fu2007, pfull, phalf,   &
                           tkel, N_drop3D, N_drop2D, N_ice3D, &
                           k_ratio, nclds, cldamt, lwp, iwp, reff_liq, &
                           reff_ice, conc_drop_org, conc_ice_org,  &
                           size_drop_org, size_ice_org)

!----------------------------------------------------------------------
!    rnd_overlap returns various cloud specification properties, 
!    obtained with the random-overlap assumption. implicit in this
!    asusmption is that all clouds are only a single layer thick; i.e.,
!    clouds at adjacent levels in the same column are independent of
!    one another.
!----------------------------------------------------------------------
 
real,    dimension(:,:,:), intent(in)             :: ql, qi, qa,  &
                                                     pfull, phalf, tkel, &
                                                     N_drop3D, N_ice3d
logical,                   intent(in)             :: use_fu2007
real,    dimension(:,:),   intent(in)             :: N_drop2D, k_ratio
integer, dimension(:,:),   intent(out)            :: nclds
real,    dimension(:,:,:), intent(out)            :: cldamt, lwp, iwp, &
                                                     reff_liq, reff_ice
real,    dimension(:,:,:), intent(out), optional  :: conc_drop_org,  &
                                                     conc_ice_org,  &
                                                     size_drop_org,  &
                                                     size_ice_org

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels, index 1 at model top 
!                    [ Pascals ]
!       tkel         Temperature [ deg Kelvin ]
!       N_drop       Number of cloud droplets per cubic meter
!       k_ratio      Ratio of effective radius to mean volume radius
!
!   intent(out) variables:
!
!       nclds        Number of (random overlapping) clouds in column 
!       cldamt       Cloud amount of condensed cloud [ dimensionless ]
!       lwp          Cloud liquid water path [ kg condensate / m **2 ]
!       iwp          Cloud ice path [ kg condensate / m **2 ]
!       reff_liq     Effective radius for liquid clouds [ microns ]
!       reff_ice     Effective particle size for ice clouds [ microns ]
!
!    intent(out), optional variables:
!
!       conc_drop_org Liquid cloud droplet mass concentration 
!                     [ g / m**3 ]
!       conc_ice_org  Ice cloud mass concentration [ g / m**3 ]
!       size_drop_org Effective diameter of liquid cloud droplets 
!                     [ microns ]
!       size_ice_org  Effective diameter of ice clouds { microns ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      logical    ::  want_microphysics
      real       ::  reff_liq_local, reff_ice_local
      integer    ::  kdim
      integer    ::  i, j, k
      REAL       ::  dumc, dumnc, rho, pgam, lamc, lammax, lammin, &
                     dumi, dumni, lami

!--------------------------------------------------------------------
!   local variables:
!
!       want_microphysics   logical indicating if microphysical 
!                           parameters are to be calculated
!       reff_liq_local      reff of liquid clouds used locally
!                           [ microns ]
!       reff_ice_local      reff of ice clouds used locally 
!                           [ microns ]
!       kdim                number of vertical layers
!       i,j,k               do-loop indices
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the number of vertical layers in the model. initialize the
!    output fields to correspond to the absence of clouds.
!---------------------------------------------------------------------
      kdim     = size(ql,3)
      nclds    = 0
      cldamt   = 0.
      lwp      = 0.
      iwp      = 0.
      reff_liq = 10.
      reff_ice = 30.

!--------------------------------------------------------------------
!    initialize the optional output arguments, if present. define a
!    logical, indicating their presence.
!--------------------------------------------------------------------
      if (present(conc_drop_org) .and. present(conc_ice_org ) .and. &
          present(size_ice_org ) .and. present(size_drop_org)) then  
        conc_drop_org = 0.
        conc_ice_org  = 0.
        size_drop_org = 20.
        size_ice_org  = 60.
        want_microphysics = .true.

!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
      else if (present(conc_drop_org) .or. present(conc_ice_org ) .or. &
               present(size_ice_org ) .or. present(size_drop_org)) then 
        call error_mesg ('cloud_rad_mod', &
            ' if any microphysical args present, all must be present',&
                                                                FATAL)

!----------------------------------------------------------------------
!    if the optional arguments are not present, set the appropriate
!    flag to indicate that only bulk properties will be calculated.
!---------------------------------------------------------------------
      else
         want_microphysics = .false.
      end if

!--------------------------------------------------------------------
!    find the layers with cloud in each column, starting at model top. 
!--------------------------------------------------------------------
      do k=1,size(ql,3)
        do j=1,size(ql,2)
          do i=1,size(ql,1)
            if (ql(i,j,k) > qmin .or. qi(i,j,k) > qmin) then
               
!---------------------------------------------------------------------
!    when cloud is found, increment the cloud column counter.
!---------------------------------------------------------------------
              nclds(i,j) = nclds(i,j) + 1

!---------------------------------------------------------------------
!    if liquid water is present, compute the liquid water path. 
!---------------------------------------------------------------------
              if (ql(i,j,k) > qmin) then
                lwp(i,j,k) = ql(i,j,k)*    &
                                      (phalf(i,j,k+1) - phalf(i,j,k))/ &
                                      GRAV/qa(i,j,k)

!----------------------------------------------------------------------
!    if microphysical properties are desired, calculate the droplet
!    concentrations. units of concentration are in g / m**3.
!----------------------------------------------------------------------
                if (want_microphysics) then
                  conc_drop_org(i,j,k) =     &
                      1000.*ql(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/&
                      RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/  &
                      MAX(phalf(i,j,k), pfull(i,j,1)))/qa(i,j,k)
                endif  

!---------------------------------------------------------------------
!    compute the effective cloud droplet radius. for liquid clouds the 
!    following formula is used, as recommended by 
!    Martin et al., J. Atmos. Sci, vol 51, pp. 1823-1842:
!
!    reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = density of cloud droplets (number per cubic meter)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!---------------------------------------------------------------------
!       if (.not. do_liq_num) then
     do_liq_num_if2:  if (.not. do_liq_num) then
                reff_liq_local = k_ratio(i,j)* 620350.49 *    &
                                 (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/   & 
                                 RDGAS/tkel(i,j,k)/DENS_H2O/    &
                                 N_drop2D(i,j))**(1./3.)
        else ! do_liq_num_if2

!cms++
! BETTER SPLIT LOOP AND MOVE MOVE IFs OUTSIDE and BETTER USE ONE SUBROUTINE
!    INTEAD OF DUPLICATING CODE!!
 
 mgp_if: IF ( .NOT. do_ice_num ) THEN
!cms--
!--------------------------------------------------------------------
! yim: a variant for prognostic droplet number
!    reff (in microns) =  k * 1.E+06 *
!                    (3*ql/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = mixing ratio of cloud droplets (number/kg air)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------
              if (ql(i,j,k) > qmin) then
                reff_liq_local = k_ratio(i,j)*620350.49*    &
                                 (ql(i,j,k)/DENS_H2O/  &
                                 max(N_drop3D(i,j,k),  &
                                      N_min*max(qa(i,j,k),qmin)/  &
                             (pfull(i,j,k)/RDGAS/tkel(i,j,k))))**(1./3.)
              else
                reff_liq_local = 0.0
              endif
                       

!cms++
ELSE !mg_if

!-----------------------------------------------------------
!
! cloud droplet effective radius
!

! in-cloud mixing ratio and number conc
         dumc = ql(i,j,k)/qa(i,j,k)
         dumnc = N_drop3D(i,j,k)/qa(i,j,k)

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

         dumc =min(dumc,5.e-3)

!cms++2008-11-26

 dumnc =  max(dumnc, N_min)

!cms--

        if ( dumc > qmin ) then


! add upper limit to in-cloud number concentration to prevent numerical error
          dumnc = min(dumnc,dumc*1.e20)

         rho=pfull(i,j,k)/(RDGAS*tkel(i,j,k))

         pgam=0.0005714*(dumnc/1.e6/rho)+0.2714
         pgam=1./(pgam**2)-1.
         pgam=max(pgam,2.)
         pgam=min(pgam,15.)

         lamc = (pi/6.*rhow*dumnc*gamma_mg(pgam+4.)/ &
                 (dumc*gamma_mg(pgam+1.)))**(1./3.)
         lammin = (pgam+1.)/max_diam_drop_mg
         lammax = (pgam+1.)/min_diam_drop_mg
         if (lamc.lt.lammin) then
         lamc = lammin
         else if (lamc.gt.lammax) then
         lamc = lammax
          end if
         reff_liq_local = gamma_mg(qcvar+1./3.)/(gamma_mg(qcvar)*qcvar**(1./3.))* &
             gamma_mg(pgam+4.)/ &
             gamma_mg(pgam+3.)/lamc/2.*1.e6
        else
        reff_liq_local = 10.
        end if
        


!ELSE!! qamin_if
!  reff_ice_local = 0.
!  reff_liq_local = 0.
!      
!END IF   qamin_if      
END IF mgp_if

end if do_liq_num_if2

!cms--

!----------------------------------------------------------------------
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly 
!    stratified cloud in liquid specific humidity, the cloud top 
!    effective radius is greater than the effective radius of the 
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!    single layer liquid or mixed phase clouds.
!----------------------------------------------------------------------
                if (do_brenguier) then
! should this be applied to all clouds ? 
! random overlap ==> all clouds are of 1 layer
                if ( k == 1 ) then
                  if (qa(i,j,2) < qamin) then
                    reff_liq_local = 1.134*reff_liq_local
                  endif
                else if (k == kdim ) then
                  if ( qa(i,j,kdim-1) < qamin) then
                    reff_liq_local = 1.134*reff_liq_local
                  endif
                else if (qa(i,j,k-1) .lt. qamin .and. & 
                         qa(i,j,k+1) .lt. qamin)  then
                  reff_liq_local = 1.134*reff_liq_local
!! ADD for random overlap -- all clouds are 1 layer thick
!!! WAIT FOR SAK REPLY
!               else
!                 reff_liq_local = 1.134*reff_liq_local
                end if
                end if


                reff_liq(i,j,k) =  reff_liq_local
                if (want_microphysics)      &
                   size_drop_org(i,j,k) = 2.*reff_liq_local
              endif  ! (ql > qmin)

!---------------------------------------------------------------------
!    if ice is present, compute the ice water path.
!---------------------------------------------------------------------
              if (qi(i,j,k) .gt. qmin) then
                iwp(i,j,k) = qi(i,j,k)*    &
                                      (phalf(i,j,k+1) - phalf(i,j,k))/ &
                                      GRAV/qa(i,j,k)
                          
!----------------------------------------------------------------------
!    if microphysical properties are desired, calculate the ice con-
!    centration. units of concentration are in g / m**3.
!----------------------------------------------------------------------
                if (want_microphysics) then
                  conc_ice_org (i,j,k) =     &
                      1000.*qi(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                      RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                      MAX(phalf(i,j,k), pfull(i,j,1)))/ qa(i,j,k)
                end if  
                       
!---------------------------------------------------------------------
!    compute the effective ice crystal size. for bulk physics cases, the
!    effective radius is taken from the formulation in Donner 
!    (1997, J. Geophys. Res., 102, pp. 21745-21768) which is based on 
!    Heymsfield and Platt (1984) with enhancement for particles smaller
!    than 20 microns.  
!    if microphysical properties are requested, then the size of the
!    ice crystals comes from the Deff column [ reference ?? ].     
!
!              T Range (K)               Reff (microns)   Deff (microns)
!     -------------------------------    --------------   --------------
!
!     Tfreeze-25. < T                       92.46298         100.6
!     Tfreeze-30. < T <= Tfreeze-25.        72.35392          80.8
!     Tfreeze-35. < T <= Tfreeze-30.        85.19071          93.5
!     Tfreeze-40. < T <= Tfreeze-35.        55.65818          63.9
!     Tfreeze-45. < T <= Tfreeze-40.        35.29989          42.5
!     Tfreeze-50. < T <= Tfreeze-45.        32.89967          39.9
!     Tfreeze-55  < T <= Tfreeze-50         16.60895          21.6
!                   T <= Tfreeze-55.        15.41627          20.2
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    calculate the effective ice crystal size using the data approp-
!    riate for the microphysics case.
!---------------------------------------------------------------------
!cms++
    do_ice_num_3: IF (  do_ice_num ) THEN

! calc. based on Morrison Gettelman code
!     
!qamin_if: IF ( qa(i,j,k) .GT. qamin ) THEN 

! in-cloud mixing ratio and number conc
         dumi = qi(i,j,k)/qa(i,j,k)
         dumni = N_ice3D(i,j,k)/qa(i,j,k)

! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
          dumi =min(dumi ,5.e-3)


!...................
! cloud ice effective radius
         if ( dumi > qmin ) then

! add upper limit to in-cloud number concentration to prevent numerical error
         dumni=min(dumni,dumi*1.e20)
         lami = (gamma_mg(1.+di_mg)*ci_mg* &
              dumni/dumi)**(1./di_mg)
         lammax = 1./ min_diam_ice_mg
         lammin = 1./(2.*dcs_mg)

         if (lami.lt.lammin) then
         lami = lammin
         else if (lami.gt.lammax) then
         lami = lammax
         end if
          reff_ice_local = 1.5/lami*1.e6
        else
         reff_ice_local = 25.
        end if


!cms++2008-10-31 

       reff_ice_local = 2. *  reff_ice_local !3rd moment/2nd moment as in Fu and Liou paper
!cms--



   ! ... assumes spherical ice particles ...
    reff_ice(i,j,k) = reff_ice_local

    size_ice_org(i,j,k) = reff_ice_local

ELSE
!cms--

            if (use_fu2007) then
!+yim Fu's parameterization of dge
              reff_ice_local = 47.05 +   &
                                0.6624*(tkel(i,j,k) - TFREEZE) +&
                               0.001741*(tkel(i,j,k)-TFREEZE)**2
              size_ice_org(i,j,k) = reff_ice_local
            else ! (use_fu2007)
                if (want_microphysics) then
                  if (tkel(i,j,k) > TFREEZE - 25.) then
                    reff_ice_local = 100.6      
                  else if (tkel(i,j,k) >  TFREEZE - 30. .and. &
                           tkel(i,j,k) <= TFREEZE - 25.) then
                    reff_ice_local = 80.8        
                  else if (tkel(i,j,k) >  TFREEZE - 35. .and. &
                           tkel(i,j,k) <= TFREEZE - 30.) then
                    reff_ice_local = 93.5       
                  else if (tkel(i,j,k) >  TFREEZE - 40. .and. &
                           tkel(i,j,k) <= TFREEZE - 35.) then
                    reff_ice_local = 63.9         
                  else if (tkel(i,j,k) >  TFREEZE - 45. .and. &
                           tkel(i,j,k) <= TFREEZE - 40.) then
                    reff_ice_local = 42.5       
                  else if (tkel(i,j,k) >  TFREEZE - 50. .and. &
                           tkel(i,j,k) <= TFREEZE - 45.) then
                    reff_ice_local = 39.9           
                  else if (tkel(i,j,k) >  TFREEZE - 55. .and. &
                           tkel(i,j,k) <= TFREEZE - 50.) then
                    reff_ice_local = 21.6         
                  else
                    reff_ice_local = 20.2             
                  endif

                  size_ice_org(i,j,k) = reff_ice_local
                endif
           endif ! (use_fu2007)

!---------------------------------------------------------------------
!    calculate reff_ice using the bulk physics data.
!---------------------------------------------------------------------
                if (tkel(i,j,k) > TFREEZE - 25.) then
                  reff_ice_local = 92.46298
                else if (tkel(i,j,k) >  TFREEZE - 30. .and. &
                         tkel(i,j,k) <= TFREEZE - 25.) then
                  reff_ice_local = 72.35392
                else if (tkel(i,j,k) >  TFREEZE - 35. .and. &
                         tkel(i,j,k) <= TFREEZE - 30.) then
                  reff_ice_local = 85.19071 
                else if (tkel(i,j,k) >  TFREEZE - 40. .and. &
                         tkel(i,j,k) <= TFREEZE - 35.) then
                  reff_ice_local = 55.65818
                else if (tkel(i,j,k) >  TFREEZE - 45. .and. &
                         tkel(i,j,k) <= TFREEZE - 40.) then
                  reff_ice_local = 35.29989
                else if (tkel(i,j,k) >  TFREEZE - 50. .and. &
                         tkel(i,j,k) <= TFREEZE - 45.) then
                  reff_ice_local = 32.89967
                else if (tkel(i,j,k) >  TFREEZE - 55. .and. &
                         tkel(i,j,k) <= TFREEZE - 50.) then
                  reff_ice_local = 16.60895
                else
                  reff_ice_local = 15.41627
                endif

                reff_ice(i,j,k) = reff_ice_local

  END IF  do_ice_num_3

              end if ! (qi > qmin)                    
                          
!---------------------------------------------------------------------
!    define the cloud fraction.
!---------------------------------------------------------------------
              cldamt(i,j,k) = qa(i,j,k)
            endif   !( ql > 0 or qi > 0)
          end do
        end do
      end do

!-------------------------------------------------------------------
    
end subroutine rnd_overlap   




!#####################################################################

! <SUBROUTINE NAME="CLOUD_RAD_k_diag">
!  <OVERVIEW>
!      This subroutine calculates the following radiative properties
!      for each cloud:
!
!<PRE>               1. r_uv : cloud reflectance in uv band
!               2. r_nir: cloud reflectance in nir band
!               3. ab_uv: cloud absorption in uv band
!               4. ab_nir:cloud absorption in nir band
!</PRE>   
!  </OVERVIEW>
!  <DESCRIPTION>
!      This subroutine calculates the following radiative properties
!      for each cloud:
!
!<PRE>               1. r_uv : cloud reflectance in uv band
!               2. r_nir: cloud reflectance in nir band
!               3. ab_uv: cloud absorption in uv band
!               4. ab_nir:cloud absorption in nir band
!</PRE>
!
!      These quantities are computed by dividing the shortwave
!      spectrum into 4 bands and then computing the reflectance
!      and absorption for each band individually and then setting
!      the uv reflectance and absorption equal to that of band
!      1 and the nir reflectance and absorption equal to the
!      spectrum weighted results of bands 2,3,and 4.  The limits
!      of bands are described in CLOUD_OPTICAL_PROPERTIES.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call CLOUD_RAD_k_diag(tau, direct, w0,gg,coszen,r_uv,r_nir,ab_uv,ab_nir)
!
!  </TEMPLATE>
!  <IN NAME="tau" TYPE="real">
!    Optical depth in 4 bands [ dimensionless ]
!  </IN>
!  <IN NAME="direct" TYPE="logical">
!    Logical variable for each cloud indicating whether
!     or not to use the direct beam solution for the
!     delta-eddington radiation or the diffuse beam
!     radiation solution.
!  </IN>
!  <IN NAME="w0" TYPE="real">
!    Single scattering albedo in 4 bands [ dimensionless ]
!  </IN>
!  <IN NAME="gg" TYPE="real">
!    Asymmetry parameter in 4 bands  [ dimensionless ]
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    Cosine of the zenith angle  [ dimensionless ]
!  </IN>
!  <INOUT NAME="r_uv" TYPE="real">
!    Cloud reflectance in uv band
!  </INOUT>
!  <INOUT NAME="r_nir" TYPE="real">
!    Cloud reflectance in nir band
!  </INOUT>
!  <INOUT NAME="ab_nir" TYPE="real">
!    Cloud absorption in nir band
!  </INOUT>
!  <INOUT NAME="ab_uv" TYPE="real">
!    Cloud absorption in uv band
!  </INOUT>
! </SUBROUTINE>

SUBROUTINE CLOUD_RAD_k_diag(tau, direct, w0,gg,coszen,r_uv,r_nir,ab_uv,ab_nir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following radiative properties
!      for each cloud:
!
!               1. r_uv : cloud reflectance in uv band
!               2. r_nir: cloud reflectance in nir band
!               3. ab_uv: cloud absorption in uv band
!               4. ab_nir:cloud absorption in nir band
!               
!
!      These quantities are computed by dividing the shortwave
!      spectrum into 4 bands and then computing the reflectance
!      and absorption for each band individually and then setting
!      the uv reflectance and absorption equal to that of band
!      1 and the nir reflectance and absorption equal to the
!      spectrum weighted results of bands 2,3,and 4.  The limits
!      of bands are described in CLOUD_OPTICAL_PROPERTIES.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!       tau          Optical depth in 4 bands (dimensionless)
!       direct       Logical variable for each cloud indicating whether
!                      or not to use the direct beam solution for the
!                      delta-eddington radiation or the diffuse beam
!                      radiation solution.
!       w0           Single scattering albedo in 4 bands (dimensionless)
!       gg           Asymmetry parameter in 4 bands (dimensionless)
!       coszen       Cosine of the zenith angle
!
!       ------
!INPUT/OUTPUT:
!       ------
!
!       r_uv         Cloud reflectance in uv band
!       r_nir        Cloud reflectance in nir band
!       ab_uv        Cloud absorption in uv band
!       ab_nir       Cloud absorption in nir band
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!      tau_local    optical depth for the band being solved
!      w0_local     single scattering albedo for the band being solved
!      g_local      asymmetry parameter for the band being solved
!      coszen_3d    3d version of coszen
!      I            looping variable
!      iband        looping variables over band number
!      taucum       cumulative sum of visible optical depth
!      g_prime      scaled g
!      w0_prime     scaled w0
!      tau_prime    scaled tau
!      crit         variable equal to 1./(4 - 3g')
!      AL           variable equal to sqrt(3*(1-w0')*(1-w0'*g'))
!      ALPHV        temporary work variable
!      GAMV         temporary work variable
!      T1V          exp( -1.*AL * tau')
!      trans_dir    direct radiation beam transmittance
!      U            1.5 * (1. - w0'*g')/AL
!      r_diffus     diffuse beam reflection
!      trans_diffus diffuse beam transmission
!
!      r            cloud reflectance for each cloud in each band
!      ab           cloud absorption for each cloud in each band
!      r_dir_uv       direct beam reflection for uv band
!      r_dir_nir      direct beam reflection for nir band
!      trans_dir_uv   direct beam transmission for uv band
!      trans_dir_nir  direct beam transmission for uv band
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:,:,:):: tau,w0,gg
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir
logical,  INTENT (IN   ),DIMENSION(:,:,:):: direct                        

!  Internal variables
!  ------------------

INTEGER                                  :: I,iband
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: coszen_3d
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: w0_local,g_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: g_prime,w0_prime
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_prime,crit,AL
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: ALPHV,GAMV,T1V,U
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: trans_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_dir,trans_dir
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3),SIZE(tau,4)) :: r,ab

!
! Code
! ----

        ! reinitialize variables
        r_uv(:,:,:) = 0.
        r_nir(:,:,:)= 0.
        ab_uv(:,:,:)= 0.
        ab_nir(:,:,:)=0.

        !create 3d zenith angle
        DO I = 1, SIZE(tau,3)
               coszen_3d(:,:,I)=coszen(:,:)
        END DO
        WHERE (coszen_3d(:,:,:) .lt. 1.E-06)
                coszen_3d(:,:,:) = 1.E-06
        END WHERE

        

    !----------------- LOOP OVER BAND -----------------------------!

    DO iband = 1, SIZE(tau,4)

        !-----------------------------------------------------------
        !  assign w0, g, tau to the value appropriate for the band

        w0_local(:,:,:) = w0(:,:,:,iband)
        tau_local(:,:,:)= tau(:,:,:,iband)
        g_local(:,:,:) =  gg(:,:,:,iband)

        !-------------------------------------------------------------------
        ! for delta-Eddington scaled ('prime') g, w0, tau where:
        !
        !               g' = g / (1 + g)
        !              w0' = (1 - g*g) * w0 / (1 - w*g*g)
        !             tau' = (1 - w*g*g) * tau
        !

                 tau_prime(:,:,:) = 1. - &
                      (w0_local(:,:,:)*g_local(:,:,:)*g_local(:,:,:))
                 w0_prime(:,:,:) = w0_local(:,:,:) * &
                   (1. - (g_local(:,:,:)*g_local(:,:,:)))/tau_prime(:,:,:)
                 tau_prime(:,:,:) = tau_prime(:,:,:) * tau_local(:,:,:)
                 g_prime(:,:,:) = g_local(:,:,:) / (1. + g_local(:,:,:))

        !-------------------------------------------------------------------
        ! create other variables
        !
        !        crit = 1./(4 - 3g')
        !
        !      and where w0' < crit set w0' = crit
        !
        !        AL = sqrt( 3. * (1. - w0') * (1. - w0'*g') )
        !

                 crit(:,:,:) = 1./(4.- 3.*g_prime(:,:,:))

                 WHERE (w0_prime(:,:,:) .lt. crit(:,:,:) )
                           w0_prime(:,:,:) = crit(:,:,:)
                 END WHERE

                 AL(:,:,:) =  ( 3. * (1. - w0_prime(:,:,:) ) &
                    * (1. - (w0_prime(:,:,:)*g_prime(:,:,:)))  )**0.5

                 !set up a minimum to AL
                 WHERE (AL(:,:,:) .lt. 1.E-06)
                        AL(:,:,:) = 1.E-06
                 END WHERE


        !-------------------------------------------------------------------
        ! simplifications if not two stream
        !
        !        ALPHV = 0.75*w0'*coszen*(1.+g'(1.-w0'))/
        !                          (1.-(AL*coszen)**2.)
        !        GAMV = 0.5*w0'*(3.*g'*(1.-w0')*coszen*coszen + 1.)/
        !                          (1.-(AL*coszen)**2.)
        !

        IF (.NOT. l2strem) THEN


                ALPHV(:,:,:) = 0.75 * w0_prime(:,:,:)*coszen_3d(:,:,:) * &
                 (1. + (g_prime(:,:,:)*(1. - w0_prime(:,:,:)))) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

                GAMV(:,:,:) =  0.50 * w0_prime(:,:,:) * &
                (  (3.* g_prime(:,:,:) * (1. - w0_prime(:,:,:)) * &
                    coszen_3d(:,:,:) * coszen_3d(:,:,:)) + 1. ) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

        END IF


        !-------------------------------------------------------------------
        ! calculate T1V
        !
        !    T1V = exp (-1* AL * tau' )


                  T1V(:,:,:) = exp( -1.*AL(:,:,:) * tau_prime(:,:,:) )


        !-------------------------------------------------------------------
        !calculate diffuse beam reflection and transmission
        !

        !first calculate U  = 1.5 * (1. - w0'*g')/AL
        U(:,:,:) = 1.5 *(1. - w0_prime(:,:,:)*g_prime(:,:,:))/AL(:,:,:)

        !initialize variables
        r_diffus(:,:,:)= 0.
        trans_diffus(:,:,:) = 1.



        trans_diffus(:,:,:) = 4. * U(:,:,:) * T1V(:,:,:) / &
            ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
              ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                   T1V(:,:,:) *   T1V(:,:,:)   )    )

        r_diffus(:,:,:) =     ((U(:,:,:)*U(:,:,:))-1.) * &
                   ( 1. -   (T1V(:,:,:)*T1V(:,:,:)) ) / &
             ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
               ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                    T1V(:,:,:) *   T1V(:,:,:)   )    )



        !-------------------------------------------------------------------
        ! calculate direct bean transmission
        !
        !
        IF (.NOT. l2strem) THEN


            !initialize variables
            trans_dir(:,:,:) = 1.
            r_dir(:,:,:) = 0.

            r_dir(:,:,:) = ( (ALPHV(:,:,:) - GAMV(:,:,:)) * &
               exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
               trans_diffus(:,:,:) ) +  &
              ( (ALPHV(:,:,:) + GAMV(:,:,:)) * &
              r_diffus(:,:,:) )  -  (ALPHV(:,:,:) - GAMV(:,:,:))

            trans_dir(:,:,:) = &
              ( (ALPHV(:,:,:)+GAMV(:,:,:))*trans_diffus(:,:,:) ) + &
              ( exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
              ( ( (ALPHV(:,:,:)-GAMV(:,:,:))*r_diffus(:,:,:) ) - &
                (ALPHV(:,:,:)+GAMV(:,:,:)) + 1. )   )

        END IF


        !-------------------------------------------------------------------
        ! patch together final solution
        !
        !


        IF (l2strem) THEN

             !two-stream solution
             r(:,:,:,iband) = r_diffus(:,:,:)
             ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) - r_diffus(:,:,:)

        ELSE

             !delta-Eddington solution
             WHERE (.not. direct)

                   r(:,:,:,iband) = r_diffus(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) &
                                     - r_diffus(:,:,:)


             END WHERE

             WHERE (direct)

                   r(:,:,:,iband) = r_dir(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_dir(:,:,:) &
                                     - r_dir(:,:,:)

             END WHERE

        END IF

    !----------------- END LOOP OVER BAND -----------------------------!

    END DO


    !----------------- CREATE SUM OVER BAND ---------------------------!

    r_uv(:,:,:) = r(:,:,:,1)
    ab_uv(:,:,:) = ab(:,:,:,1)

    r_nir(:,:,:) =  (  0.326158 * r(:,:,:,2) + &
                       0.180608 * r(:,:,:,3) + &
                       0.033474 * r(:,:,:,4) ) / 0.540240

    ab_nir(:,:,:) =  (  0.326158 * ab(:,:,:,2) + &
                        0.180608 * ab(:,:,:,3) + &
                        0.033474 * ab(:,:,:,4) ) / 0.540240


        !-------------------------------------------------------------------
        ! guarantee that clouds for tau = 0. have the properties
        ! of no cloud

        
        WHERE(tau(:,:,:,1) .le. 0.)
             r_uv(:,:,:) = 0.
             ab_uv(:,:,:)= 0.                       
        END WHERE
        WHERE((tau(:,:,:,2)+tau(:,:,:,3)+tau(:,:,:,4)) .le. 0.)
             r_nir(:,:,:)= 0.
             ab_nir(:,:,:)=0.
        END WHERE       

        !-------------------------------------------------------------------
        ! guarantee that for coszen lt. or equal to zero that solar
        ! reflectances and absorptances are equal to zero.
        DO I = 1, SIZE(tau,3)
               WHERE (coszen(:,:) .lt. 1.E-06)
                    r_uv(:,:,I) = 0. 
                    ab_uv(:,:,I) = 0.
                    r_nir(:,:,I) = 0.
                    ab_nir(:,:,I) = 0.
               END WHERE
        END DO
        
        !-------------------------------------------------------------------
        ! guarantee that each cloud has some transmission by reducing
        ! the actual cloud reflectance in uv and nir band
        ! this break is necessary to avoid the rest of the
        ! radiation code from breaking up.
        !

        WHERE ( (1. - r_uv(:,:,:) - ab_uv(:,:,:)) .lt. 0.01)
                      r_uv(:,:,:) = r_uv(:,:,:) - 0.01
        END WHERE
        WHERE ( (1. - r_nir(:,:,:) - ab_nir(:,:,:)) .lt. 0.01)
                      r_nir(:,:,:) = r_nir(:,:,:) - 0.01
        END WHERE

        !-------------------------------------------------------------------
        ! guarantee that cloud reflectance and absorption are greater than
        ! or equal to zero

        WHERE (r_uv(:,:,:) .lt. 0.)
               r_uv(:,:,:) = 0.
        END WHERE
        WHERE (r_nir(:,:,:) .lt. 0.)
               r_nir(:,:,:) = 0.
        END WHERE
        WHERE (ab_uv(:,:,:) .lt. 0.)
               ab_uv(:,:,:) = 0.
        END WHERE
        WHERE (ab_nir(:,:,:) .lt. 0.)
               ab_nir(:,:,:) = 0.
        END WHERE


END SUBROUTINE CLOUD_RAD_k_diag



!#####################################################################

! <SUBROUTINE NAME="cloud_rad_k">
!  <OVERVIEW>
!    Subroutine cloud_rad_k calculates the cloud reflectances and
!    absorptions in the uv and nir wavelength bands. These quantities 
!    are computed by dividing the shortwave spectrum into 4 bands and 
!    then computing the reflectance and absorption for each band 
!    individually and then setting the uv reflectance and absorption 
!    equal to that of band 1 and the nir reflectance and absorption 
!    equal to the spectrum-weighted results of bands 2,3,and 4.  The 
!    limits of bands are defined in subroutine sw_optical_properties.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    Subroutine cloud_rad_k calculates the cloud reflectances and
!    absorptions in the uv and nir wavelength bands. These quantities 
!    are computed by dividing the shortwave spectrum into 4 bands and 
!    then computing the reflectance and absorption for each band 
!    individually and then setting the uv reflectance and absorption 
!    equal to that of band 1 and the nir reflectance and absorption 
!    equal to the spectrum-weighted results of bands 2,3,and 4.  The 
!    limits of bands are defined in subroutine sw_optical_properties.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_rad_k (tau, w0, gg, coszen, r_uv, r_nir,    &
!                     ab_nir, ab_uv_out)
!
!  </TEMPLATE>
!  <IN NAME="tau" TYPE="real">
!    Optical depth [ dimensionless ]
!  </IN>
!  <IN NAME="w0" TYPE="real">
!    Single scattering albedo [ dimensionless ]
!  </IN>
!  <IN NAME="gg" TYPE="real">
!    Asymmetry parameter for each band
!    [ dimensionless ]
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    Cosine of zenith angle  [ dimensionless ]
!  </IN>
!  <INOUT NAME="r_uv" TYPE="real">
!    Cloud reflectance in uv band
!  </INOUT>
!  <INOUT NAME="r_nir" TYPE="real">
!    Cloud reflectance in nir band
!  </INOUT>
!  <INOUT NAME="ab_nir" TYPE="real">
!    Cloud absorption in nir band
!  </INOUT>
!  <INOUT NAME="ab_uv_out" TYPE="real, optional">
!    Cloud absorption in uv band
!  </INOUT>
! </SUBROUTINE>
!
subroutine cloud_rad_k (tau, w0, gg, coszen, r_uv, r_nir,    &
                        ab_nir, ab_uv_out)

!----------------------------------------------------------------------
!    subroutine cloud_rad_k calculates the cloud reflectances and
!    absorptions in the uv and nir wavelength bands. these quantities 
!    are computed by dividing the shortwave spectrum into 4 bands and 
!    then computing the reflectance and absorption for each band 
!    individually and then setting the uv reflectance and absorption 
!    equal to that of band 1 and the nir reflectance and absorption 
!    equal to the spectrum-weighted results of bands 2,3,and 4.  The 
!    limits of bands are defined in subroutine sw_optical_properties.
!----------------------------------------------------------------------

real, dimension(:,:,:,:), intent(in)             :: tau, w0, gg
real, dimension(:,:),     intent(in)             :: coszen
real, dimension(:,:,:),   intent(inout)          :: r_uv, r_nir, ab_nir
real, dimension(:,:,:),   intent(inout),optional :: ab_uv_out

!---------------------------------------------------------------------
!   intent(in) variables:
!
!         tau            Optical depth [ dimensionless ]
!         w0             Single scattering albedo [ dimensionless ]
!         gg             Asymmetry parameter for each band
!                        [ dimensionless ]
!         coszen         Cosine of zenith angle  [ dimensionless ]
!
!    intent(inout) variables:
!
!         r_uv            Cloud reflectance in uv band
!         r_nir           Cloud reflectance in nir band
!         ab_nir          Cloud absorption in nir band
!
!    intent(inout), optional variables:
!
!         ab_uv_out       Cloud absorption in uv band
! 
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real,    dimension (size(tau,1), size(tau,2)) :: taucum

      logical, dimension (size(tau,1), size(tau,2),                 &
                                       size(tau,3)) :: direct

      real,    dimension (size(tau,1), size(tau,2),                  &
                                       size(tau,3)) ::              &
                          coszen_3d, tau_local, w0_local, g_local, &
                          g_prime, w0_prime, tau_prime, crit, al,   &
                          alphv, gamv, t1v, u, denom, r_diffus,   &
                          trans_diffus, r_dir, trans_dir, ab_uv

      real, dimension (size(tau,1), size(tau,2),                    &
                       size(tau,3), size(tau,4)) :: r, ab

      integer   :: iband, k

!---------------------------------------------------------------------
!   local variables:
!
!     taucum       cumulative sum of visible optical depth from top
!                  of atmosphere to current layer
!     direct       logical variable for each cloud indicating whether
!                  or not to use the direct beam solution for the
!                  delta-eddington radiation or the diffuse beam
!                  radiation solution.
!      coszen_3d    3d version of coszen
!      tau_local    optical depth for the band being solved
!      w0_local     single scattering albedo for the band being solved
!      g_local      asymmetry parameter for the band being solved
!      g_prime      scaled g
!      w0_prime     scaled w0
!      tau_prime    scaled tau
!      crit         variable equal to 1./(4 - 3g')
!      al           variable equal to sqrt(3*(1-w0')*(1-w0'*g'))
!      alphv        temporary work variable
!      gamv         temporary work variable
!      t1v          exp( -1.*AL * tau')
!      u            1.5 * (1. - w0'*g')/AL
!      denom        [(u+1)*(u+1)] - [(u-1)*(u-1)*t1v*t1v]
!      r_diffus     diffuse beam reflection
!      trans_diffus diffuse beam transmission
!      r_dir        diffuse beam reflection
!      trans_dir    direct radiation beam transmittance
!      ab_uv        cloud absorption in uv band
!      r            cloud reflectance for each cloud in each band
!      ab           cloud absorption for each cloud in each band
!      i,j,k,iband  do-loop indices
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    define some values needed when the delta-eddington approach is
!    being taken.
!-------------------------------------------------------------------
      if (.not. l2strem) then

!--------------------------------------------------------------------
!    create 3d version of zenith angle array. allow it to be no 
!    smaller than 1.0e-06.
!--------------------------------------------------------------------
        do k=1,size(tau,3)
          coszen_3d(:,:,k) = coszen(:,:)
        end do
        where (coszen_3d(:,:,:) < 1.E-06) 
          coszen_3d(:,:,:) = 1.E-06
        end where
        
!--------------------------------------------------------------------
!    initialize taucum and direct. taucum will be the accumulated tau
!    from model top to the current level and will be the basis of
!    deciding whether the incident beam is treated as direct or diffuse.
!--------------------------------------------------------------------
        taucum(:,:) = 0.
        direct(:,:,:) = .true.
        do k=1,size(tau,3)      

!---------------------------------------------------------------------
!    find if taucum from model top to level above has exceeded taucrit.
!----------------------------------------------------------------------
          where (taucum(:,:) > taucrit)
            direct(:,:,k) = .false.
          end where

!---------------------------------------------------------------------
!    increment the cumulative tau.
!---------------------------------------------------------------------
          taucum(:,:) = taucum(:,:) + tau(:,:,k,1)
        end do
      endif 

!---------------------------------------------------------------------
!    loop over the wavelength bands, calculating the reflectivity 
!    and absorption for each band.
!---------------------------------------------------------------------
      do iband=1,size(tau,4)

!---------------------------------------------------------------------
!    assign local values of w0, g and  tau appropriate for the current
!    band under consideration.
!---------------------------------------------------------------------
        w0_local (:,:,:) = w0 (:,:,:,iband)
        tau_local(:,:,:) = tau(:,:,:,iband)
        g_local  (:,:,:) = gg (:,:,:,iband)

!---------------------------------------------------------------------
!    define the scaled (prime) values of g, w0 and tau used in the 
!    delta-Eddington calculation:
!               g' = g / (1 + g)
!              w0' = (1 - g*g) * w0 / (1 - w*g*g)
!             tau' = (1 - w*g*g) * tau
!---------------------------------------------------------------------
        tau_prime(:,:,:) = 1. - (w0_local(:,:,:)*g_local(:,:,:)*   &
                                 g_local(:,:,:))
        w0_prime(:,:,:) = w0_local(:,:,:)*(1. - (g_local(:,:,:)*    &
                          g_local(:,:,:)))/tau_prime(:,:,:)
        tau_prime(:,:,:) = tau_prime(:,:,:)*tau_local(:,:,:)
        g_prime(:,:,:) = g_local(:,:,:)/(1. + g_local(:,:,:))

!-------------------------------------------------------------------
!    define other variables:
!      crit = 1./(4 - 3g')  and where w0' < crit set w0' = crit
!      al = sqrt( 3. * (1. - w0') * (1. - w0'*g') ) and where 
!      al < 1.0e-06, set al = 1.0e-06.
!--------------------------------------------------------------------
        crit(:,:,:) = 1./(4.- 3.*g_prime(:,:,:))
        where (w0_prime(:,:,:) < crit(:,:,:) )
          w0_prime(:,:,:) = crit(:,:,:)
        end where
        al(:,:,:) = (3.*(1. - w0_prime(:,:,:) ) *      &
                    (1. - (w0_prime(:,:,:)*g_prime(:,:,:))) )**0.5
        where (al(:,:,:) < 1.E-06)
          al(:,:,:) = 1.E-06
        end where

!--------------------------------------------------------------------
!    calculate t1v:
!    t1v = exp (-1* al * tau')
!--------------------------------------------------------------------
        t1v(:,:,:) = exp(-1.*al(:,:,:)*tau_prime(:,:,:))

!-------------------------------------------------------------------
!    calculate diffuse beam reflection and transmission. first calculate
!    the factor u  = 1.5 * (1. - w0'*g')/al and the value denom =
!    [ (u+1)*(u+1) - (u-1)*(u-1)*t1v*t1v ].
!---------------------------------------------------------------------
        u(:,:,:) = 1.5*(1. - w0_prime(:,:,:)*g_prime(:,:,:))/al(:,:,:)
        denom(:,:,:) = ( ( (u(:,:,:) + 1.)*(u(:,:,:) + 1.) ) - &
                         ( (u(:,:,:) - 1.)*(u(:,:,:) - 1.)* &
                            t1v(:,:,:)*t1v(:,:,:)   )    )
        trans_diffus(:,:,:) = 4.*u(:,:,:)*t1v(:,:,:)/denom(:,:,:)
        r_diffus(:,:,:) = ((u(:,:,:)*u(:,:,:)) - 1.) * &
                          (1. - (t1v(:,:,:)*t1v(:,:,:)) )/denom(:,:,:)

!-------------------------------------------------------------------
!    calculate direct beam transmission for the delta-eddington case.
!    alphv = 0.75*w0'*coszen*(1.+g'(1.-w0'))/
!            (1.-(al*coszen)**2.)
!    gamv = 0.5*w0'*(3.*g'*(1.-w0')*coszen*coszen + 1.)/
!           (1.-(al*coszen)**2.)
!-------------------------------------------------------------------
        if (.not. l2strem) then
          where (direct)
            alphv(:,:,:) = 0.75 * w0_prime(:,:,:)*coszen_3d(:,:,:)* &
                           (1. + (g_prime(:,:,:)*    &
                           (1. - w0_prime(:,:,:))))/ &
                           (1. - (al(:,:,:)*coszen_3d(:,:,:))**2)
            gamv(:,:,:) =  0.50 * w0_prime(:,:,:)* &
                           (  (3.* g_prime(:,:,:)*(1. -     &
                           w0_prime(:,:,:))*&
                           coszen_3d(:,:,:)*coszen_3d(:,:,:)) + 1. )/ &
                           (1. - (al(:,:,:)*coszen_3d(:,:,:))**2)
            r_dir(:,:,:) = ( (alphv(:,:,:) - gamv(:,:,:)) * &
                           exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:))* &
                           trans_diffus(:,:,:) ) + ((alphv(:,:,:) +   &
                           gamv(:,:,:)) * r_diffus(:,:,:) )  -  &
                           (alphv(:,:,:) - gamv(:,:,:))
            trans_dir(:,:,:) = ( (alphv(:,:,:) + gamv(:,:,:))*  &
                               trans_diffus(:,:,:) ) + &
                               (exp(-1.*tau_prime(:,:,:)/   &
                                coszen_3d(:,:,:)) * &
                                ( ( (alphv(:,:,:) - gamv(:,:,:))*  &
                                r_diffus(:,:,:) ) - (alphv(:,:,:) + &
                                gamv(:,:,:)) + 1. )   )
          end where
        endif 

!-------------------------------------------------------------------
!    patch together final solution.
!-------------------------------------------------------------------
        if (l2strem) then

!---------------------------------------------------------------------
!    two-stream solution.
!---------------------------------------------------------------------
          r(:,:,:,iband)  = r_diffus(:,:,:)
          ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) - r_diffus(:,:,:)
        else

!----------------------------------------------------------------------
!    delta-eddington solution.
!----------------------------------------------------------------------
          where (.not. direct)
            r (:,:,:,iband) = r_diffus(:,:,:)
            ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) - r_diffus(:,:,:)
          end where
          where (direct)
            r (:,:,:,iband) = r_dir(:,:,:)
            ab(:,:,:,iband) = 1. - trans_dir(:,:,:) - r_dir(:,:,:)
          end where
        endif 
      end do  ! (iband loop)


!---------------------------------------------------------------------
!    sum over the apprpriate bands to obtain values for the uv and nir
!    bands.
!----------------------------------------------------------------------

      r_uv (:,:,:) = r (:,:,:,1)
      ab_uv(:,:,:) = ab(:,:,:,1)
      r_nir(:,:,:) = ( 0.326158*r(:,:,:,2) + &
                       0.180608*r(:,:,:,3) + &
                       0.033474*r(:,:,:,4) ) / 0.540240
      ab_nir(:,:,:) =  ( 0.326158*ab(:,:,:,2) + &
                         0.180608*ab(:,:,:,3) + &
                         0.033474*ab(:,:,:,4) ) / 0.540240

!-------------------------------------------------------------------
!    guarantee that when tau = 0., the sw properties are clear-sky 
!    values.
!-------------------------------------------------------------------
      where (tau(:,:,:,1) <= 0.)
        r_uv (:,:,:) = 0.
        ab_uv(:,:,:) = 0.                       
      end where
      where ((tau(:,:,:,2) + tau(:,:,:,3) + tau(:,:,:,4)) <= 0.)
        r_nir (:,:,:) = 0.
        ab_nir(:,:,:) = 0.
      end where       

!-------------------------------------------------------------------
!    guarantee that for coszen .le. 1.0E-6 that solar reflectances and 
!    absorptances are equal to zero.
!-------------------------------------------------------------------
      do k=1,size(tau,3)
        where (coszen(:,:) < 1.E-06)
          r_uv  (:,:,k) = 0. 
          ab_uv (:,:,k) = 0.
          r_nir (:,:,k) = 0.
          ab_nir(:,:,k) = 0.
        end where
      end do
        
!-------------------------------------------------------------------
!    guarantee that each cloud has some transmission by reducing the 
!    actual cloud reflectance. this break is necessary to avoid the 
!    rest of the radiation code from breaking up.
!---------------------------------------------------------------------
      where ( (1. - r_uv(:,:,:) - ab_uv(:,:,:)) < 0.01)
        r_uv(:,:,:) = r_uv(:,:,:) - 0.01
      end where
      where ( (1. - r_nir(:,:,:) - ab_nir(:,:,:)) < 0.01)
        r_nir(:,:,:) = r_nir(:,:,:) - 0.01
      end where

!-------------------------------------------------------------------
!    guarantee that cloud reflectance and absorption are .ge. 0.0.
!-------------------------------------------------------------------
      where (r_uv(:,:,:) < 0.)
        r_uv(:,:,:) = 0.
      end where
      where (r_nir(:,:,:) < 0.)
        r_nir(:,:,:) = 0.
      end where
      where (ab_uv(:,:,:) < 0.)
        ab_uv(:,:,:) = 0.
      end where
      where (ab_nir(:,:,:) < 0.)
        ab_nir(:,:,:) = 0.
      end where

!--------------------------------------------------------------------
!    if ab_uv is desired as an output variable from this subroutine,
!    fill the variable being returned.
!--------------------------------------------------------------------
      if (present (ab_uv_out)) then
        ab_uv_out = ab_uv
      endif

!----------------------------------------------------------------------



end subroutine cloud_rad_k


!#####################################################################

SUBROUTINE CLOUD_RAD(tau,w0,gg,coszen,r_uv,r_nir,ab_uv,ab_nir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following radiative properties
!      for each cloud:
!
!               1. r_uv : cloud reflectance in uv band
!               2. r_nir: cloud reflectance in nir band
!               3. ab_uv: cloud absorption in uv band
!               4. ab_nir:cloud absorption in nir band
!               
!
!      These quantities are computed by dividing the shortwave
!      spectrum into 4 bands and then computing the reflectance
!      and absorption for each band individually and then setting
!      the uv reflectance and absorption equal to that of band
!      1 and the nir reflectance and absorption equal to the
!      spectrum weighted results of bands 2,3,and 4.  The limits
!      of bands are described in CLOUD_OPTICAL_PROPERTIES.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!       tau          optical depth in 4 bands (dimensionless)
!       w0           single scattering albedo in 4 bands (dimensionless)
!       gg           asymmetry parameter in 4 bands (dimensionless)
!       coszen       cosine of the zenith angle
!
!       ------
!INPUT/OUTPUT:
!       ------
!
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!      direct       logical variable for each cloud indicating whether
!                       or not to use the direct beam solution for the
!                       delta-eddington radiation or the diffuse beam
!                       radiation solution.
!      tau_local    optical depth for the band being solved
!      w0_local     single scattering albedo for the band being solved
!      g_local      asymmetry parameter for the band being solved
!      coszen_3d    3d version of coszen
!      I            looping variable
!      iband        looping variables over band number
!      taucum       cumulative sum of visible optical depth
!      g_prime      scaled g
!      w0_prime     scaled w0
!      tau_prime    scaled tau
!      crit         variable equal to 1./(4 - 3g')
!      AL           variable equal to sqrt(3*(1-w0')*(1-w0'*g'))
!      ALPHV        temporary work variable
!      GAMV         temporary work variable
!      T1V          exp( -1.*AL * tau')
!      trans_dir    direct radiation beam transmittance
!      U            1.5 * (1. - w0'*g')/AL
!      r_diffus     diffuse beam reflection
!      trans_diffus diffuse beam transmission
!
!      r            cloud reflectance for each cloud in each band
!      ab           cloud absorption for each cloud in each band
!      r_dir_uv       direct beam reflection for uv band
!      r_dir_nir      direct beam reflection for nir band
!      trans_dir_uv   direct beam transmission for uv band
!      trans_dir_nir  direct beam transmission for uv band
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:,:,:):: tau,w0,gg
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir

!  Internal variables
!  ------------------

INTEGER                                  :: I,iband
REAL,    DIMENSION(SIZE(tau,1),SIZE(tau,2)) :: taucum
LOGICAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: direct
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: coszen_3d
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: w0_local,g_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: g_prime,w0_prime
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_prime,crit,AL
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: ALPHV,GAMV,T1V,U
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: trans_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_dir,trans_dir
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3),SIZE(tau,4)) :: r,ab

!
! Code
! ----

        ! reinitialize variables
        r_uv(:,:,:) = 0.
        r_nir(:,:,:)= 0.
        ab_uv(:,:,:)= 0.
        ab_nir(:,:,:)=0.

        !create 3d zenith angle
        DO I = 1, SIZE(tau,3)
               coszen_3d(:,:,I)=coszen(:,:)
        END DO
        WHERE (coszen_3d(:,:,:) .lt. 1.E-06)
                coszen_3d(:,:,:) = 1.E-06
        END WHERE

        
        !------------------------------------------------------------------
        !do logical variable to determine where total cloud optical depth
        !at uv wavelengths exceeds taucrit
        IF (.NOT. l2strem) THEN

                !---- initialize taucum and direct
                taucum(:,:)=0.
                direct(:,:,:)=.TRUE.

                DO I = 1, SIZE(tau,3)

                      !find if taucum to levels above has exceeded taucrit
                      WHERE (taucum(:,:) .gt. taucrit)
                            direct(:,:,I)=.FALSE.
                      END WHERE

                      !increment cumulative tau
                      taucum(:,:)=taucum(:,:)+tau(:,:,I,1)

                END DO
        END IF

    !----------------- LOOP OVER BAND -----------------------------!

    DO iband = 1, SIZE(tau,4)

        !-----------------------------------------------------------
        !  assign w0, g, tau to the value appropriate for the band

        w0_local(:,:,:) = w0(:,:,:,iband)
        tau_local(:,:,:)= tau(:,:,:,iband)
        g_local(:,:,:) =  gg(:,:,:,iband)

        !-------------------------------------------------------------------
        ! for delta-Eddington scaled ('prime') g, w0, tau where:
        !
        !               g' = g / (1 + g)
        !              w0' = (1 - g*g) * w0 / (1 - w*g*g)
        !             tau' = (1 - w*g*g) * tau
        !

                 tau_prime(:,:,:) = 1. - &
                      (w0_local(:,:,:)*g_local(:,:,:)*g_local(:,:,:))
                 w0_prime(:,:,:) = w0_local(:,:,:) * &
                   (1. - (g_local(:,:,:)*g_local(:,:,:)))/tau_prime(:,:,:)
                 tau_prime(:,:,:) = tau_prime(:,:,:) * tau_local(:,:,:)
                 g_prime(:,:,:) = g_local(:,:,:) / (1. + g_local(:,:,:))

        !-------------------------------------------------------------------
        ! create other variables
        !
        !        crit = 1./(4 - 3g')
        !
        !      and where w0' < crit set w0' = crit
        !
        !        AL = sqrt( 3. * (1. - w0') * (1. - w0'*g') )
        !

                 crit(:,:,:) = 1./(4.- 3.*g_prime(:,:,:))

                 WHERE (w0_prime(:,:,:) .lt. crit(:,:,:) )
                           w0_prime(:,:,:) = crit(:,:,:)
                 END WHERE

                 AL(:,:,:) =  ( 3. * (1. - w0_prime(:,:,:) ) &
                    * (1. - (w0_prime(:,:,:)*g_prime(:,:,:)))  )**0.5

                 !set up a minimum to AL
                 WHERE (AL(:,:,:) .lt. 1.E-06)
                        AL(:,:,:) = 1.E-06
                 END WHERE


        !-------------------------------------------------------------------
        ! simplifications if not two stream
        !
        !        ALPHV = 0.75*w0'*coszen*(1.+g'(1.-w0'))/
        !                          (1.-(AL*coszen)**2.)
        !        GAMV = 0.5*w0'*(3.*g'*(1.-w0')*coszen*coszen + 1.)/
        !                          (1.-(AL*coszen)**2.)
        !

        IF (.NOT. l2strem) THEN


                ALPHV(:,:,:) = 0.75 * w0_prime(:,:,:)*coszen_3d(:,:,:) * &
                 (1. + (g_prime(:,:,:)*(1. - w0_prime(:,:,:)))) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

                GAMV(:,:,:) =  0.50 * w0_prime(:,:,:) * &
                (  (3.* g_prime(:,:,:) * (1. - w0_prime(:,:,:)) * &
                    coszen_3d(:,:,:) * coszen_3d(:,:,:)) + 1. ) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

        END IF


        !-------------------------------------------------------------------
        ! calculate T1V
        !
        !    T1V = exp (-1* AL * tau' )


                  T1V(:,:,:) = exp( -1.*AL(:,:,:) * tau_prime(:,:,:) )


        !-------------------------------------------------------------------
        !calculate diffuse beam reflection and transmission
        !

        !first calculate U  = 1.5 * (1. - w0'*g')/AL
        U(:,:,:) = 1.5 *(1. - w0_prime(:,:,:)*g_prime(:,:,:))/AL(:,:,:)

        !initialize variables
        r_diffus(:,:,:)= 0.
        trans_diffus(:,:,:) = 1.



        trans_diffus(:,:,:) = 4. * U(:,:,:) * T1V(:,:,:) / &
            ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
              ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                   T1V(:,:,:) *   T1V(:,:,:)   )    )

        r_diffus(:,:,:) =     ((U(:,:,:)*U(:,:,:))-1.) * &
                   ( 1. -   (T1V(:,:,:)*T1V(:,:,:)) ) / &
             ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
               ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                    T1V(:,:,:) *   T1V(:,:,:)   )    )



        !-------------------------------------------------------------------
        ! calculate direct bean transmission
        !
        !
        IF (.NOT. l2strem) THEN


            !initialize variables
            trans_dir(:,:,:) = 1.
            r_dir(:,:,:) = 0.

            r_dir(:,:,:) = ( (ALPHV(:,:,:) - GAMV(:,:,:)) * &
               exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
               trans_diffus(:,:,:) ) +  &
              ( (ALPHV(:,:,:) + GAMV(:,:,:)) * &
              r_diffus(:,:,:) )  -  (ALPHV(:,:,:) - GAMV(:,:,:))

            trans_dir(:,:,:) = &
              ( (ALPHV(:,:,:)+GAMV(:,:,:))*trans_diffus(:,:,:) ) + &
              ( exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
              ( ( (ALPHV(:,:,:)-GAMV(:,:,:))*r_diffus(:,:,:) ) - &
                (ALPHV(:,:,:)+GAMV(:,:,:)) + 1. )   )

        END IF


        !-------------------------------------------------------------------
        ! patch together final solution
        !
        !


        IF (l2strem) THEN

             !two-stream solution
             r(:,:,:,iband) = r_diffus(:,:,:)
             ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) - r_diffus(:,:,:)

        ELSE

             !delta-Eddington solution
             WHERE (.not. direct)

                   r(:,:,:,iband) = r_diffus(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) &
                                     - r_diffus(:,:,:)


             END WHERE

             WHERE (direct)

                   r(:,:,:,iband) = r_dir(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_dir(:,:,:) &
                                     - r_dir(:,:,:)

             END WHERE

        END IF

    !----------------- END LOOP OVER BAND -----------------------------!

    END DO


    !----------------- CREATE SUM OVER BAND ---------------------------!

    r_uv(:,:,:) = r(:,:,:,1)
    ab_uv(:,:,:) = ab(:,:,:,1)

    r_nir(:,:,:) =  (  0.326158 * r(:,:,:,2) + &
                       0.180608 * r(:,:,:,3) + &
                       0.033474 * r(:,:,:,4) ) / 0.540240

    ab_nir(:,:,:) =  (  0.326158 * ab(:,:,:,2) + &
                        0.180608 * ab(:,:,:,3) + &
                        0.033474 * ab(:,:,:,4) ) / 0.540240


        !-------------------------------------------------------------------
        ! guarantee that clouds for tau = 0. have the properties
        ! of no cloud

        
        WHERE(tau(:,:,:,1) .le. 0.)
             r_uv(:,:,:) = 0.
             ab_uv(:,:,:)= 0.                       
        END WHERE
        WHERE((tau(:,:,:,2)+tau(:,:,:,3)+tau(:,:,:,4)) .le. 0.)
             r_nir(:,:,:)= 0.
             ab_nir(:,:,:)=0.
        END WHERE       

        !-------------------------------------------------------------------
        ! guarantee that for coszen lt. or equal to zero that solar
        ! reflectances and absorptances are equal to zero.
        DO I = 1, SIZE(tau,3)
               WHERE (coszen(:,:) .lt. 1.E-06)
                    r_uv(:,:,I) = 0. 
                    ab_uv(:,:,I) = 0.
                    r_nir(:,:,I) = 0.
                    ab_nir(:,:,I) = 0.
               END WHERE
        END DO
        
        !-------------------------------------------------------------------
        ! guarantee that each cloud has some transmission by reducing
        ! the actual cloud reflectance in uv and nir band
        ! this break is necessary to avoid the rest of the
        ! radiation code from breaking up.
        !

        WHERE ( (1. - r_uv(:,:,:) - ab_uv(:,:,:)) .lt. 0.01)
                      r_uv(:,:,:) = r_uv(:,:,:) - 0.01
        END WHERE
        WHERE ( (1. - r_nir(:,:,:) - ab_nir(:,:,:)) .lt. 0.01)
                      r_nir(:,:,:) = r_nir(:,:,:) - 0.01
        END WHERE

        !-------------------------------------------------------------------
        ! guarantee that cloud reflectance and absorption are greater than
        ! or equal to zero

        WHERE (r_uv(:,:,:) .lt. 0.)
               r_uv(:,:,:) = 0.
        END WHERE
        WHERE (r_nir(:,:,:) .lt. 0.)
               r_nir(:,:,:) = 0.
        END WHERE
        WHERE (ab_uv(:,:,:) .lt. 0.)
               ab_uv(:,:,:) = 0.
        END WHERE
        WHERE (ab_nir(:,:,:) .lt. 0.)
               ab_nir(:,:,:) = 0.
        END WHERE


END SUBROUTINE CLOUD_RAD

!######################################################################

! <SUBROUTINE NAME="sw_optical_properties">
!  <OVERVIEW>
!    sw_optical_properties computes the needed optical parameters and
!    then calls cloud_rad_k in order to compute the cloud radiative
!    properties.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    sw_optical_properties computes the needed optical parameters and
!    then calls cloud_rad_k in order to compute the cloud radiative
!    properties.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sw_optical_properties (nclds, lwp, iwp, reff_liq, reff_ice, &
!                               coszen, r_uv, r_nir, ab_nir)
!
!  </TEMPLATE>
!  <IN NAME="nclds" TYPE="integer">
!    Number of random overlapping clouds in column





!  </IN>
!  <IN NAME="lwp" TYPE="real">
!    Liquid water path [ kg / m**2 ]
!  </IN>
!  <IN NAME="iwp" TYPE="real">
!    Ice water path [ kg / m**2 ]
!  </IN>
!  <IN NAME="reff_liq" TYPE="real">
!    Effective cloud drop radius  used with
!    bulk cloud physics scheme [ microns ]
!  </IN>
!  <IN NAME="reff_ice" TYPE="real">
!    Effective ice crystal radius used with
!    bulk cloud physics scheme [ microns ]
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    Cosine of zenith angle [ dimensionless ]
!  </IN>
!  <INOUT NAME="r_uv" TYPE="real">
!    Cloud reflectance in uv band
!  </INOUT>
!  <INOUT NAME="r_nir" TYPE="real">
!    Cloud reflectance in nir band
!  </INOUT>
!  <INOUT NAME="ab_nir" TYPE="real">
!    Cloud absorption in nir band
!  </INOUT>
! </SUBROUTINE>
!
subroutine sw_optical_properties (nclds, lwp, iwp, reff_liq, reff_ice, &
                                  coszen, r_uv, r_nir, ab_nir)

!----------------------------------------------------------------------
!    sw_optical_properties computes the needed optical parameters and
!    then calls cloud_rad_k in order to compute the cloud radiative
!    properties.
!---------------------------------------------------------------------
                              
integer, dimension(:,:),   intent(in)      :: nclds
real,    dimension(:,:,:), intent(in)      :: lwp, iwp, reff_liq,   &
                                              reff_ice
real,    dimension(:,:),   intent(in)      :: coszen
real,    dimension(:,:,:), intent(inout)   :: r_uv, r_nir, ab_nir

!----------------------------------------------------------------------
!    intent(in) variables:
!
!        nclds           Number of random overlapping clouds in column
!        lwp             Liquid water path [ kg / m**2 ]
!        iwp             Ice water path [ kg / m**2 ]
!        reff_liq        Effective cloud drop radius  used with
!                        bulk cloud physics scheme [ microns ]
!        reff_ice        Effective ice crystal radius used with
!                        bulk cloud physics scheme [ microns ]
!        coszen          Cosine of zenith angle [ dimensionless ]
!
!    intent(out) variables:
!
!        r_uv            Cloud reflectance in uv band
!        r_nir           Cloud reflectance in nir band
!        ab_nir          Cloud absorption in nir band
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(lwp,1),                                   &
                       size(lwp,2),                                   &
                       size(lwp,3), 4) ::  tau, w0, gg, tau_liq,   &
                                           tau_ice, w0_liq, w0_ice, &
                                           g_liq, g_ice
      integer :: max_cld
      integer :: i,j, k, m

!---------------------------------------------------------------------
!   local variables:
!
!           tau            optical depth [ dimensionless ]
!           w0             single scattering albedo
!           gg             asymmetry parameter for each band
!           tau_liq        optical depth due to cloud liquid 
!                          [ dimensionless ]
!           tau_ice        optical depth due to cloud liquid 
!                          [ dimensionless ]
!           w0_liq         single scattering albedo due to liquid
!           w0_ice         single scattering albedo due to ice
!           g_liq          asymmetry factor for cloud liquid
!           g_ice          asymmetry factor for cloud ice      
!           max_cld        largest number of clouds in any column in
!                          the current physics window
!           i,j,l,m        do-loop indices 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the maximum number of random overlap clouds in any column
!    in the current physics window.
!---------------------------------------------------------------------
      max_cld = maxval (nclds(:,:))

!--------------------------------------------------------------------
!    spatial loop.
!--------------------------------------------------------------------
      do j=1,size(lwp,2)
        do i=1,size(lwp,1)
          do k=1,size(lwp,3)  

!---------------------------------------------------------------------
!    initialize local variables to default values.
!---------------------------------------------------------------------
            tau    (i,j,k ,:) = 0.
            tau_liq(i,j,k ,:) = 0.
            tau_ice(i,j,k ,:) = 0.
            gg     (i,j,k ,:) = 0.85
            g_liq  (i,j,k ,:) = 0.85
            g_ice  (i,j,k ,:) = 0.85
            w0     (i,j,k ,:) = 0.95
            w0_liq (i,j,k ,:) = 0.95
            w0_ice (i,j,k ,:) = 0.95

!--------------------------------------------------------------------
!    compute uv cloud optical depths due to liquid droplets in each 
!    of the wavelength bands. the formulas for optical depth come from 
!    Slingo (1989, J. Atmos. Sci., vol. 46, pp. 1419-1427)
!--------------------------------------------------------------------
            tau_liq(i,j,k,1) = lwp(i,j,k)*1000.* &
                               (0.02817 + (1.305/reff_liq(i,j,k)))
            tau_liq(i,j,k,2) = lwp(i,j,k)*1000.* &
                               (0.02682 + (1.346/reff_liq(i,j,k)))
            tau_liq(i,j,k,3) = lwp(i,j,k)*1000.* &
                               (0.02264 + (1.454/reff_liq(i,j,k)))
            tau_liq(i,j,k,4) = lwp(i,j,k)*1000.* &
                               (0.01281 + (1.641/reff_liq(i,j,k)))
        
!--------------------------------------------------------------------
!    compute uv cloud optical depths due to ice crystals. the ice
!    optical depth is independent of wavelength band. the formulas for 
!    optical depth come from Ebert and Curry (1992, J. Geophys. Res., 
!    vol. 97, pp. 3831-3836. IMPORTANT!!! NOTE WE ARE CHEATING HERE 
!    BECAUSE WE ARE FORCING THE FIVE BAND MODEL OF EBERT AND CURRY INTO
!    THE FOUR BAND MODEL OF SLINGO. THIS IS DONE BY COMBINING BANDS 3 
!    and 4 OF EBERT AND CURRY TOGETHER. EVEN SO THE EXACT BAND LIMITS 
!    DO NOT MATCH.  FOR COMPLETENESS HERE ARE THE BAND LIMITS (MICRONS)
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!--------------------------------------------------------------------
            tau_ice(i,j,k,1) = iwp(i,j,k)*1000.* &
                               (0.003448 + (2.431/reff_ice(i,j,k)))
            tau_ice(i,j,k,2) = tau_ice(i,j,k,1)
            tau_ice(i,j,k,3) = tau_ice(i,j,k,1)
            tau_ice(i,j,k,4) = tau_ice(i,j,k,1)
        
!--------------------------------------------------------------------
!    compute total cloud optical depth as the sum of the liquid and
!    ice components. the mixed phase optical properties are based upon 
!    equation 14 of Rockel et al. 1991, Contributions to Atmospheric 
!    Physics, volume 64, pp.1-12:  
!           tau = tau_liq + tau_ice
!--------------------------------------------------------------------
            tau(i,j,k,:) = tau_liq(i,j,k,:) + tau_ice(i,j,k,:)
        
!---------------------------------------------------------------------
!    compute single-scattering albedo resulting from the presence
!    of liquid droplets.
!---------------------------------------------------------------------
            w0_liq(i,j,k,1) =  5.62E-08 - 1.63E-07*reff_liq(i,j,k)
            w0_liq(i,j,k,2) =  6.94E-06 - 2.35E-05*reff_liq(i,j,k)
            w0_liq(i,j,k,3) = -4.64E-04 - 1.24E-03*reff_liq(i,j,k)
            w0_liq(i,j,k,4) = -2.01E-01 - 7.56E-03*reff_liq(i,j,k)

            w0_liq(i,j,k,:) = w0_liq(i,j,k,:) + 1.

!---------------------------------------------------------------------
!    compute single-scattering albedo resulting from the presence
!    of ice crystals.
!---------------------------------------------------------------------
            w0_ice(i,j,k,1) = -1.00E-05
            w0_ice(i,j,k,2) = -1.10E-04 - 1.41E-05*reff_ice(i,j,k)
            w0_ice(i,j,k,3) = -1.86E-02 - 8.33E-04*reff_ice(i,j,k)
            w0_ice(i,j,k,4) = -4.67E-01 - 2.05E-05*reff_ice(i,j,k)

            w0_ice(i,j,k,:) = w0_ice(i,j,k,:) + 1.
          end do
        end do
      end do

!----------------------------------------------------------------------
!    compute total single scattering albedo. the mixed phase value is
!    obtained from equation 14 of Rockel et al. 1991, Contributions to 
!    Atmospheric Physics, volume 64, pp.1-12:
!           w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                   (          tau_liq  +           tau_ice )
!----------------------------------------------------------------------
      do j=1,size(lwp,2)
        do i=1,size(lwp,1)
          do k=1, size(lwp,3)    
            do m=1,4
              if (tau(i,j,k,m) > 0.0) then
                w0(i,j,k,m) = (w0_liq(i,j,k,m) * tau_liq(i,j,k,m) + &
                               w0_ice(i,j,k,m) * tau_ice(i,j,k,m)) / &
                               tau(i,j,k,m)
              endif
            end do
          end do
        end do
      end do

!---------------------------------------------------------------------
!    compute asymmetry factor.
!---------------------------------------------------------------------
      do j=1,size(lwp,2)
        do i=1,size(lwp,1)
          do k=1, size(lwp,3)  

!---------------------------------------------------------------------
!    compute asymmetry factor resulting from the presence
!    of liquid droplets.
!---------------------------------------------------------------------
            g_liq(i,j,k,1) = 0.829 + 2.482E-03*reff_liq(i,j,k)
            g_liq(i,j,k,2) = 0.794 + 4.226E-03*reff_liq(i,j,k)
            g_liq(i,j,k,3) = 0.754 + 6.560E-03*reff_liq(i,j,k)
            g_liq(i,j,k,4) = 0.826 + 4.353E-03*reff_liq(i,j,k)

!---------------------------------------------------------------------
!    compute asymmetry factor resulting from the presence
!    of ice crystals.
!---------------------------------------------------------------------
            g_ice(i,j,k,1) = 0.7661 + 5.851E-04*reff_ice(i,j,k)
            g_ice(i,j,k,2) = 0.7730 + 5.665E-04*reff_ice(i,j,k)
            g_ice(i,j,k,3) = 0.7940 + 7.267E-04*reff_ice(i,j,k)
            g_ice(i,j,k,4) = 0.9595 + 1.076E-04*reff_ice(i,j,k)
          end do
        end do
      end do

!----------------------------------------------------------------------
!    compute combined asymmetry factor, including effects of both
!    liquid droplets and ice crystals. the mixed phase value is obtained
!    from equation 14 of Rockel et al. 1991, Contributions to 
!    Atmospheric Physics, volume 64, pp.1-12: 
!      g  = ( g_liq * w0_liq * tau_liq + g_ice * w0_ice * tau_ice ) /
!           (         w0_liq * tau_liq +         w0_ice * tau_ice )
!----------------------------------------------------------------------
      do j=1,size(lwp,2)
        do i=1,size(lwp,1)
          do k=1, size(lwp,3)   
            do m=1,4
              if (tau(i,j,k,m) > 0.0) then
                gg(i,j,k,m) = (w0_liq(i,j,k,m)*g_liq(i,j,k,m)*   &
                               tau_liq(i,j,k,m) + w0_ice(i,j,k,m)*  &
                               g_ice(i,j,k,m)*tau_ice(i,j,k,m) ) / &
                               (w0_liq(i,j,k,m)*tau_liq(i,j,k,m) + &
                                w0_ice(i,j,k,m)*tau_ice(i,j,k,m) )

              endif
            end do
          end do
        end do
      end do
        
!--------------------------------------------------------------------
!    apply constraints to the values of these variables.
!--------------------------------------------------------------------
      do j=1,size(lwp,2)
        do i=1,size(lwp,1)
          do k=1, size(lwp,3)
            do m=1,4
              if (tau(i,j,k,m) < taumin .and. tau(i,j,k,m) /= 0.0 ) then
                tau(i,j,k,m) = taumin
              endif
            end do
          end do
        end do
      end do

!---------------------------------------------------------------------
!    account for plane-parallel homogenous cloud bias.
!---------------------------------------------------------------------
      tau(:,:,:,:) = scale_factor*tau(:,:,:,:)

!---------------------------------------------------------------------
!    call cloud_rad_k to calculate the cloud radiative properties.
!---------------------------------------------------------------------
     call cloud_rad_k (tau, w0, gg, coszen, r_uv, r_nir, ab_nir)

!----------------------------------------------------------------------


end subroutine sw_optical_properties


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!         USED WITH ORIGINAL FMS RADIATION:
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!#####################################################################

SUBROUTINE CLOUD_SUMMARY (is,js,                        &
                  LAND,ql,qi,qa,qv,pfull,phalf,TKel,coszen,skt,&
                  nclds,ktop,kbot,cldamt,Time,&
                  r_uv,r_nir,ab_uv,ab_nir,em_lw,&
                  conc_drop,conc_ice,size_drop,size_ice)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following properties of clouds
!
!               1. nclds: # of clouds
!               2. ktop : integer level for top of cloud
!               3. kbot : integer level for bottom of cloud
!               4. cldamt:horizontal cloud amount of every cloud
!
!      Optional arguments
!               5. r_uv : cloud reflectance in uv band
!               6. r_nir: cloud reflectance in nir band
!               7. ab_uv: cloud absorption in uv band
!               8. ab_nir:cloud absorption in nir band
!               9. em_lw :longwave cloud emmissivity
!              10. conc_drop : liquid cloud droplet mass concentration
!              11. conc_ice  : ice cloud mass concentration
!              12. size_drop : effective diameter of liquid cloud droplets
!              13. size_ice  : effective diameter of ice cloud 
!
!      given inputs of ql and qi (liquid and ice condensate),
!      cloud volume fraction, and pressure at the half and full levels
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!       is,js        indices for model slab
!       LAND         fraction of the grid box covered by LAND
!       ql           cloud liquid condensate (kg condensate/kg air)
!       qi           cloud ice condensate (kg condensate/kg air)
!       qa           cloud volume fraction (fraction)
!       qv           water vapor specific humidity (kg vapor/kg air)
!       pfull        pressure at full levels (Pascals)
!       phalf        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       TKel            temperature (Kelvin)
!       coszen       cosine of the zenith angle
!       skt          surface skin temperature (Kelvin)
!
!       -------------
!INPUT/OUTPUT:
!       -------------
!
!       nclds        number of (random overlapping) clouds in column and also
!                        the current # for clouds to be operating on
!       ktop         level of the top of the cloud
!       kbot         level of the bottom of the cloud
!       cldamt       cloud amount of condensed cloud
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!       em_lw        longwave cloud emmissivity
!       conc_drop    liquid cloud droplet mass concentration (g /m3)
!       conc_ice     ice cloud mass concentration (g /m3)
!       size_drop    effective diameter of liquid cloud droplets (microns)
!       size_ice   : effective diameter of ice clouds (microns)
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t      looping variables
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of third dimension points
!       N_drop       number of cloud droplets per cubic meter
!       k_ratio      ratio of effective radius to mean volume radius
!       max_cld      maximum number of clouds in whole array
!       qa_local     local value of qa (fraction)
!       ql_local     local value of ql (kg condensate / kg air)
!       qi_local     local value of qi (kg condensate / kg air)
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!       tau          optical depth in 4 bands (dimensionless)
!       w0           single scattering albedo in 4 bands (dimensionless)
!       gg           asymmetry parameter in 4 bands (dimensionless)
!       rad_prop     logical indicating if you are requesting the
!                    radiative properties of the clouds
!       wat_prop     logical determining if you are requesting the
!                    concentrations and particle sizes of clouds
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

INTEGER,  INTENT (IN)                    :: is,js
REAL,     INTENT (IN), DIMENSION(:,:)    :: LAND,skt
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: ql,qi,qa,qv,pfull,TKel
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
INTEGER,  INTENT (INOUT),DIMENSION(:,:)  :: nclds
INTEGER,  INTENT (INOUT),DIMENSION(:,:,:):: ktop,kbot
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cldamt
type(time_type), intent(in), optional    :: Time
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir,em_lw
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: conc_drop,conc_ice
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: size_drop,size_ice



!  Internal variables
!  ------------------

INTEGER                                           :: i,j,IDIM,JDIM,KDIM,max_cld
LOGICAL                                           :: rad_prop, wat_prop
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2),SIZE(ql,3)) :: qa_local,ql_local,qi_local
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2))            :: N_drop, k_ratio
REAL, DIMENSION(:,:,:), allocatable               :: r_uv_local, r_nir_local
REAL, DIMENSION(:,:,:), allocatable               :: ab_uv_local, ab_nir_local
REAL, DIMENSION(:,:,:), allocatable               :: em_lw_local
REAL, DIMENSION(:,:,:), allocatable               :: conc_drop_local,conc_ice_local
REAL, DIMENSION(:,:,:), allocatable               :: size_drop_local,size_ice_local
REAL, DIMENSION(:,:,:), allocatable               :: LWP,IWP,Reff_liq,Reff_ice
REAL, DIMENSION(:,:,:,:), allocatable             :: tau,w0,gg


!
! Code
! ----

    
    ! reinitialize variables
    IDIM=SIZE(ql,1)
    JDIM=SIZE(ql,2)
    KDIM=SIZE(ql,3)
    nclds(:,:)    = 0
    ktop(:,:,:)   = 0
    kbot(:,:,:)   = 0
    cldamt(:,:,:) = 0.

    rad_prop = .FALSE.
    wat_prop = .FALSE.
    if (PRESENT(r_uv).or.PRESENT(r_nir).or.PRESENT(ab_uv).or.PRESENT(ab_nir).or.&
        PRESENT(em_lw)) then
        rad_prop = .TRUE.
    end if
    if (PRESENT(conc_drop).or.PRESENT(conc_ice).or.PRESENT(size_drop).or.&
        PRESENT(size_ice)) then
        wat_prop = .TRUE.
    end if
    if ((.not.rad_prop).and.(.not.wat_prop)) then
        rad_prop = .TRUE.
    end if

    !create local values of ql and qi
    !this step is necessary to remove the values of (qi,ql) which are
    !           0 < (qi,ql) < qmin   or
    !               (qi,ql) > qmin and qa <= qamin

    ql_local(:,:,:) = 0.
    qi_local(:,:,:) = 0.
    qa_local(:,:,:) = 0.
    WHERE ( (qa(:,:,:) .gt. qamin) .and. (ql(:,:,:) .gt. qmin) )
                  ql_local(:,:,:) = ql(:,:,:)
                  qa_local(:,:,:) = qa(:,:,:)
    END WHERE
    WHERE ( (qa(:,:,:) .gt. qamin) .and. (qi(:,:,:) .gt. qmin) )
                  qi_local(:,:,:) = qi(:,:,:)
                  qa_local(:,:,:) = qa(:,:,:)
    END WHERE

    !compute N_drop and k_ratio
    N_drop(:,:)=N_land*LAND(:,:) + N_ocean*(1.-LAND(:,:))
    k_ratio(:,:)=k_land*LAND(:,:) + k_ocean*(1.-LAND(:,:))

    !do solution for new radiation code
    if (wat_prop) then

         ALLOCATE(conc_drop_local(IDIM,JDIM,KDIM))
         ALLOCATE(conc_ice_local(IDIM,JDIM,KDIM))
         ALLOCATE(size_drop_local(IDIM,JDIM,KDIM))
         ALLOCATE(size_ice_local(IDIM,JDIM,KDIM))
         
         conc_drop_local(:,:,:) = 0.
         conc_ice_local(:,:,:)  = 0.
         size_drop_local(:,:,:) = 20.
         size_ice_local(:,:,:)  = 60.

         call  cloud_organize(ql_local,qi_local,qa_local,&
                  pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  conc_drop_org=conc_drop_local,&
                  conc_ice_org =conc_ice_local,&
                  size_drop_org=size_drop_local,&
                  size_ice_org =size_ice_local)

         !assign to output
         if (PRESENT(conc_drop)) conc_drop = scale_factor*conc_drop_local
         if (PRESENT(conc_ice)) conc_ice = scale_factor*conc_ice_local
         if (PRESENT(size_drop)) size_drop = size_drop_local
         if (PRESENT(size_ice)) size_ice = size_ice_local
    
         DEALLOCATE(conc_drop_local)
         DEALLOCATE(conc_ice_local)
         DEALLOCATE(size_drop_local)
         DEALLOCATE(size_ice_local)
         
    end if
     
         
    !do solution for old radiation code
    if (rad_prop) then

    
         ALLOCATE(r_uv_local(IDIM,JDIM,KDIM))
         ALLOCATE(r_nir_local(IDIM,JDIM,KDIM))
         ALLOCATE(ab_uv_local(IDIM,JDIM,KDIM))
         ALLOCATE(ab_nir_local(IDIM,JDIM,KDIM))
         ALLOCATE(em_lw_local(IDIM,JDIM,KDIM))
         ALLOCATE(LWP(IDIM,JDIM,KDIM))
         ALLOCATE(IWP(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_liq(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_ice(IDIM,JDIM,KDIM))
         ALLOCATE(tau(IDIM,JDIM,KDIM,4))
         ALLOCATE(w0(IDIM,JDIM,KDIM,4))
         ALLOCATE(gg(IDIM,JDIM,KDIM,4))
         
         r_uv_local(:,:,:)   = 0.
         r_nir_local(:,:,:)  = 0.
         ab_uv_local(:,:,:)  = 0.
         ab_nir_local(:,:,:) = 0.
         em_lw_local(:,:,:)  = 0.
         LWP(:,:,:)    = 0.
         IWP(:,:,:)    = 0.
         Reff_liq(:,:,:) = 10.
         Reff_ice(:,:,:) = 30.
         tau(:,:,:,:)    = 0.
         w0(:,:,:,:)     = 0.
         gg(:,:,:,:)     = 0.

    
         call  cloud_organize(ql_local,qi_local,qa_local,&
                  pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  LWP_in=LWP,IWP_in=IWP,Reff_liq_in=Reff_liq,Reff_ice_in=Reff_ice)
         
         !find maximum number of clouds
         max_cld  = MAXVAL(nclds(:,:))
         
         !compute cloud radiative properties
         IF (max_cld .gt. 0) then

              CALL CLOUD_OPTICAL_PROPERTIES(LWP(:,:,1:max_cld),IWP(:,:,1:max_cld),&
                      Reff_liq(:,:,1:max_cld),Reff_ice(:,:,1:max_cld),&
                      tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),gg(:,:,1:max_cld,:),&
                      em_lw_local(:,:,1:max_cld))
              
              !Account for plane-parallel homogenous cloud bias
              tau(:,:,:,:) = scale_factor * tau(:,:,:,:)

              !compute cloud radiative properties
              CALL CLOUD_RAD(tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),&
                       gg(:,:,1:max_cld,:),coszen,&
                       r_uv_local(:,:,1:max_cld),r_nir_local(:,:,1:max_cld),&
                       ab_uv_local(:,:,1:max_cld),ab_nir_local(:,:,1:max_cld))
              
              !assure that zero clouds have properties of zero clouds
              DO i = 1, IDIM
              DO j = 1, JDIM
                       if (nclds(i,j).lt.max_cld) then
                               r_uv_local(i,j,nclds(i,j)+1:max_cld)   = 0.
                               r_nir_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                               ab_uv_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                               ab_nir_local(i,j,nclds(i,j)+1:max_cld) = 0.
                               em_lw_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                       end if
              ENDDO
              ENDDO

         END IF

         if (PRESENT(r_uv)) r_uv = r_uv_local
         if (PRESENT(r_nir)) r_nir = r_nir_local
         if (PRESENT(ab_uv)) ab_uv = ab_uv_local
         if (PRESENT(ab_nir)) ab_nir = ab_nir_local
         if (PRESENT(em_lw)) em_lw = em_lw_local
          
         DEALLOCATE(r_uv_local)
         DEALLOCATE(r_nir_local)
         DEALLOCATE(ab_uv_local)
         DEALLOCATE(ab_nir_local)
         DEALLOCATE(em_lw_local)
         DEALLOCATE(LWP)
         DEALLOCATE(IWP)
         DEALLOCATE(Reff_liq)
         DEALLOCATE(Reff_ice)
         DEALLOCATE(tau)
         DEALLOCATE(w0)
         DEALLOCATE(gg)
         
    end if
    
END SUBROUTINE CLOUD_SUMMARY

SUBROUTINE CLOUD_ORGANIZE(ql,qi,qa,pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  LWP_in,IWP_in,Reff_liq_in,Reff_ice_in, &
                  conc_drop_org,conc_ice_org,size_drop_org,size_ice_org)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following properties of clouds
!
!               1. nclds: # of clouds
!               2. ktop : integer level for top of cloud
!               3. kbot : integer level for bottom of cloud
!               4. cldamt:horizontal cloud amount of every cloud
!               5. LWP :
!
!      given inputs of ql and qi (liquid and ice condensate),
!      cloud volume fraction, and pressure at the half and full levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!       LAND         fraction of the grid box covered by LAND
!       ql           cloud liquid condensate (kg condensate/kg air)
!       qi           cloud ice condensate (kg condensate/kg air)
!       qa           cloud volume fraction (fraction)
!       pfull        pressure at full levels (Pascals)
!       phalf        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       TKel            temperature (Kelvin)
!       coszen       cosine of the zenith angle
!       N_drop       number of cloud droplets per cubic meter
!       k_ratio      ratio of effective radius to mean volume radius
!
!       -------------
!INPUT/OUTPUT:
!       -------------
!
!       nclds        number of (random overlapping) clouds in column and also
!                        the current # for clouds to be operating on
!       ktop         level of the top of the cloud
!       kbot         level of the bottom of the cloud
!       cldamt       cloud amount of condensed cloud
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!       conc_drop_org liquid cloud droplet mass concentration (g /m3)
!       conc_ice_org  ice cloud mass concentration (g /m3)
!       size_drop_org effective diameter of liquid cloud droplets (microns)
!       size_ice_org  effective diameter of ice clouds (microns)
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t      looping variables
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of vertical levels
!       nlev         number of levels in the cloud
!       reff_liq_local   reff of liquid clouds used locally (microns)
!       sum_reff_liq  a sum of reff_liq_local
!       reff_ice_local   reff of ice clouds used locally (microns)
!       sum_reff_ice  a sum of reff_liq_local
!       sum_liq      sum of liquid in cloud (kg condensate per square meter)
!       sum_ice      sum of ice in cloud (kg condensate per square meter)
!       maxcldfrac   maximum cloud fraction in cloud block (fraction)
!       top_t,bot_t  temporary integers used to identify cloud edges
!       totcld_bot   total cloud fraction from bottom view
!       max_bot      largest cloud fraction face from bottom view
!       totcld_top   total cloud fraction from top view
!       max_top      largest cloud fraction face from top view
!       tmp_val      temporary number used in the assigning of top and bottoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      NOTE ON THE FORMULAS FOR EFFECTIVE RADIUS OF LIQUID AND ICE
!      CLOUDS:
!
!
!      FOR LIQUID CLOUDS THE FOLLOWING FORMULA IS USED:
!
!      THIS FORMULA IS THE RECOMMENDATION OF
!      Martin et al., J. Atmos. Sci, vol 51, pp. 1823-1842
!
!
!        reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!       where airdens = density of air in kg air/m3
!                  ql = liquid condensate in kg cond/kg air
!                  qa = cloud fraction
!                  pi = 3.14159
!            Dens_h2o = density of pure liquid water (kg liq/m3) 
!               N_liq = density of cloud droplets (number per cubic meter)
!                   k = factor to account for difference between 
!                       mean volume radius and effective radius
!
!        IN THIS PROGRAM reff_liq is limited to be between 4.2 microns
!        and 16.6 microns, which is the range of validity for the
!        Slingo (1989) radiation.
!
!     For single layer liquid or mixed phase clouds it is assumed that
!     cloud liquid is vertically stratified within the cloud.  Under
!     such situations for observed stratocumulus clouds it is found
!     that the cloud mean effective radius is between 80 and 100% of
!     the cloud top effective radius. (Brenguier et al., Journal of
!     Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  For linearly 
!     stratified cloud in liquid specific humidity, the cloud top 
!     effective radius is greater than the effective radius of the 
!     cloud mean specific humidity by a factor of 2**(1./3.).
!
!     This correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!     single layer liquid or mixed phase clouds.
!
!
!     FOR ICE CLOUDS THE EFFECTIVE RADIUS IS TAKEN FROM THE FORMULATION
!     IN DONNER (1997, J. Geophys. Res., 102, pp. 21745-21768) WHICH IS
!     BASED ON HEYMSFIELD AND PLATT (1984) WITH ENHANCEMENT FOR PARTICLES
!     SMALLER THAN 20 MICRONS.  
!
!              T Range (K)               Reff (microns)   Deff (microns)
!     -------------------------------    --------------   --------------
!
!     Tfreeze-25. < T                       92.46298         100.6
!     Tfreeze-30. < T <= Tfreeze-25.        72.35392          80.8
!     Tfreeze-35. < T <= Tfreeze-30.        85.19071          93.5
!     Tfreeze-40. < T <= Tfreeze-35.        55.65818          63.9
!     Tfreeze-45. < T <= Tfreeze-40.        35.29989          42.5
!     Tfreeze-50. < T <= Tfreeze-45.        32.89967          39.9
!     Tfreeze-55  < T <= Tfreeze-50         16.60895          21.6
!                   T <= Tfreeze-55.        15.41627          20.2
!
!        IN THIS PROGRAM, reff_ice is limited to be between 10 microns
!        and 130 microns, which is the range of validity for the Ebert
!        and Curry (1992) radiation.
!
!        IN THIS PROGRAM, size_ice (i.e. Deff) is limited to be between
!        18.6 microns and 130.2 microns, which is the range of validity 
!        for the Fu Liou JAS 1993 radiation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen, N_drop, k_ratio
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: ql,qi,qa,pfull,TKel
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
INTEGER,  INTENT (INOUT),DIMENSION(:,:)  :: nclds
INTEGER,  INTENT (INOUT),DIMENSION(:,:,:):: ktop,kbot
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cldamt
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: LWP_in,IWP_in,Reff_liq_in,Reff_ice_in
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: conc_drop_org,conc_ice_org
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: size_drop_org,size_ice_org


!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM
INTEGER                                  :: t,top_t,bot_t
INTEGER                                  :: tmp_top,tmp_bot,nlev
LOGICAL                                  :: add_cld,lhsw,sea_esf
REAL                                     :: sum_liq,sum_ice,maxcldfrac
REAL                                     :: totcld_bot,max_bot
REAL                                     :: totcld_top,max_top,tmp_val
REAL                                     :: reff_liq_local,sum_reff_liq
REAL                                     :: reff_ice_local,sum_reff_ice
real, dimension(:,:,:), allocatable :: lwp, iwp, reff_liq, reff_ice

!
! Code
! ----


    ! reinitialize variables
    IDIM=SIZE(ql,1)
    JDIM=SIZE(ql,2)
    KDIM=SIZE(ql,3)
    nclds(:,:)    = 0
    ktop(:,:,:)   = 0
    kbot(:,:,:)   = 0
    cldamt(:,:,:) = 0.

    !decide which type of output is necessary
    lhsw = .FALSE.
    sea_esf = .FALSE.
    if (PRESENT(conc_drop_org).or.PRESENT(conc_ice_org).or.  &
        PRESENT(size_drop_org).or.PRESENT(size_ice_org)) then
        sea_esf = .TRUE.
    end if
    if (PRESENT(LWP_in).or.PRESENT(IWP_in).or.PRESENT(Reff_liq_in).or. &
        PRESENT(Reff_ice_in)) then
        lhsw = .true.
    end if
        allocate ( lwp(idim, jdim,kdim))
        allocate ( iwp(idim, jdim,kdim))
        allocate ( reff_liq(idim, jdim,kdim))
        allocate ( reff_ice(idim, jdim,kdim))
    if ((.not.lhsw).and.(.not.sea_esf)) then
        lhsw = .TRUE.
    end if


    !initialize output fields
         LWP(:,:,:)    = 0.
         IWP(:,:,:)    = 0.
         Reff_liq(:,:,:) = 10.
         Reff_ice(:,:,:) = 30.

    if (sea_esf) then
         conc_drop_org(:,:,:) = 0.
         conc_ice_org (:,:,:) = 0.
         size_drop_org(:,:,:) = 20.
         size_ice_org (:,:,:) = 60.
    end if
        


    !-----------  DETERMINE CLOUD AMOUNT, LWP, IWP FOR EACH CLOUD ------!

    if (overlap .eq. 1) then


         !-----------  prevent user from attempting to do maximum
         !-----------  -random overlap with new radiation code

         if (sea_esf) then

              call error_mesg  ('cloud_rad_organize in cloud_rad module',&
                                'maximum random overlap is not currently '//&
                                'available with sea_esf radiation', FATAL)

        end if

        !---- DO CONDENSING OF CLOUDS ----!

        
        !---loop over vertical levels----!
        DO i = 1, IDIM
        DO j = 1, JDIM
        
        add_cld  = .FALSE.

        DO k = 1, KDIM

                 !identify new cloud tops
                 IF ( ( (ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin) ) .and. &
                           (.NOT. add_cld) ) then
                       nclds(i,j) = nclds(i,j) + 1
                       add_cld = .TRUE.
                       ktop(i,j,nclds(i,j)) = k
                       sum_liq          = 0.
                       sum_ice          = 0.
                       maxcldfrac       = 0.
                       sum_reff_liq     = 0.
                       sum_reff_ice     = 0.        
                 END IF

                 !increment sums where cloud
                 IF (   (ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin)  ) then

                       !compute reff
                       if (ql(i,j,k) .gt. qmin) then
                          reff_liq_local = k_ratio(i,j)* 620350.49 *    &
                             (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/Rdgas/   &
                             TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       else
                          reff_liq_local = 0.
                       end if

                       if ( (k .eq. 1    .and. qa(i,j,2)      .lt. qamin) &
                          .or. &
                          (k .eq. KDIM .and. qa(i,j,KDIM-1) .lt. qamin) &
                          .or. &
                          (k .gt. 1 .and. k .lt. KDIM .and. &
                          qa(i,j,k-1) .lt. qamin .and. &
                          qa(i,j,k+1) .lt. qamin) ) then
                          reff_liq_local = 1.134 * reff_liq_local
                       end if

                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       reff_liq_local = MIN(16.6,reff_liq_local)
                       reff_liq_local = MAX(4.2, reff_liq_local)

                       if (qi(i,j,k) .gt. qmin) then
                          
                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                             reff_ice_local = 92.46298
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                             reff_ice_local = 72.35392
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                             reff_ice_local = 85.19071 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                             reff_ice_local = 55.65818
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                             reff_ice_local = 35.29989
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                             reff_ice_local = 32.89967
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                             reff_ice_local = 16.60895
                          else
                             reff_ice_local = 15.41627
                          end if
                          !limit values to that for which Ebert and
                          !Curry radiation is valid :
                          !  10 microns < reff < 130 microns
                          !
                          reff_ice_local = MIN(130.,reff_ice_local)
                          reff_ice_local = MAX(10.,reff_ice_local)

                       else
                          reff_ice_local = 0.
                       end if   !end if for qi > qmin

                       !increment sums
                       sum_liq = sum_liq + ql(i,j,k)* &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       sum_ice = sum_ice + qi(i,j,k)* &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       maxcldfrac = MAX(maxcldfrac,qa(i,j,k))
                       sum_reff_liq  = sum_reff_liq + &
                          ( reff_liq_local * ql(i,j,k) * &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav )
                       sum_reff_ice  = sum_reff_ice + &
                          ( reff_ice_local * qi(i,j,k) * &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav )

                 END IF


                 !where the first cloud gap exists after a cloud
                 ! or bottom level is reached compute kbot, cldamt,
                 ! LWP, IWP, Reff_liq, and Reff_ice
                 IF (  ( (ql(i,j,k) .le. qmin) .and. &
                         (qi(i,j,k) .le. qmin) .and. &
                         (add_cld) ) .or. &
                         (add_cld .and. k .eq. KDIM)) then

                    !reset add_cld
                    add_cld     = .FALSE.

                    !determine kbot
                    kbot(i,j,nclds(i,j))= k-1
                    if ((ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin)) then
                       kbot(i,j,nclds(i,j)) = k
                    end if

                    cldamt(i,j,nclds(i,j)) = maxcldfrac
                    LWP(i,j,nclds(i,j)) = sum_liq / cldamt(i,j,nclds(i,j))
                    IWP(i,j,nclds(i,j)) = sum_ice / cldamt(i,j,nclds(i,j))
                    if (sum_liq .gt. 0.) then
                       Reff_liq(i,j,nclds(i,j)) = sum_reff_liq / sum_liq
                    end if
                    if (sum_ice .gt. 0.) then
                       Reff_ice(i,j,nclds(i,j)) = sum_reff_ice / sum_ice
                    end if

                    ! If adjust_top is T then
                    !change top and bottom indices to those that
                    !are at the most exposed to top and bottom
                    !view
                    nlev = kbot(i,j,nclds(i,j))-ktop(i,j,nclds(i,j))+1
                    if (adjust_top .and. nlev .gt. 1) then

                       !reset tmp_top,tmp_bot
                       tmp_top = ktop(i,j,nclds(i,j))
                       tmp_bot = kbot(i,j,nclds(i,j))

                       !find top and base of cloud
                       totcld_bot=0.
                       totcld_top=0.
                       max_bot=0.
                       max_top=0.
          
                       DO t = 1,nlev

                          top_t = ktop(i,j,nclds(i,j))+t-1
                          bot_t = kbot(i,j,nclds(i,j))-t+1
                          
                          tmp_val = MAX(0.,qa(i,j,top_t)-totcld_top)
                          if (tmp_val .gt. max_top) then
                             max_top = tmp_val
                             tmp_top = top_t
                          end if
                          totcld_top = totcld_top+tmp_val         
                              
                          tmp_val = MAX(0.,qa(i,j,bot_t)-totcld_bot)
                          if (tmp_val .gt. max_bot) then
                             max_bot = tmp_val
                             tmp_bot = bot_t
                          end if
                          totcld_bot = totcld_bot+tmp_val         
                               
                       END DO
                       
                       !assign tmp_top and tmp_bot to ktop and kbot
                       ktop(i,j,nclds(i,j)) = tmp_top
                       kbot(i,j,nclds(i,j)) = tmp_bot

                    end if  !for adjust_top

                 END IF  !for end of cloud

        END DO
        END DO
        END DO

    else if (overlap .eq. 2) then

           
        !---loop over vertical levels----!
        DO i = 1, IDIM
        DO j = 1, JDIM
        DO k = 1, KDIM
               
                 !where cloud exists compute ktop,kbot, cldamt and LWP and IWP
                 IF ( (ql(i,j,k) .gt. qmin) .or. &
                      (qi(i,j,k) .gt. qmin)  ) then

                    nclds(i,j) = nclds(i,j) + 1
                    ktop(i,j,nclds(i,j)) = k
                    kbot(i,j,nclds(i,j)) = k

                    if (lhsw) then
                       cldamt(i,j,nclds(i,j)) = qa(i,j,k)
                       LWP(i,j,nclds(i,j))    = ql(i,j,k)*    &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav/ &
                          cldamt(i,j,nclds(i,j))
                       IWP(i,j,nclds(i,j))    = qi(i,j,k)*    &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav/ &
                          cldamt(i,j,nclds(i,j))
                    end if  !lhsw if

                    if (sea_esf) then
                       cldamt(i,j,k) = qa(i,j,k)
                       !Note units are in g/m3!
                       conc_drop_org(i,j,k) = 1000.*ql(i,j,k)*                &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Rdgas/TKel(i,j,k)/    &
                          log(phalf(i,j,k+1)/MAX(phalf(i,j,k),pfull(i,j,1)))/ &
                          cldamt(i,j,k)
                       conc_ice_org (i,j,k) = 1000.*qi(i,j,k)*                &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Rdgas/TKel(i,j,k)/    &
                          log(phalf(i,j,k+1)/MAX(phalf(i,j,k),pfull(i,j,1)))/ &
                          cldamt(i,j,k)
                    end if  !sea_esf if

                    !compute reff_liquid
                    if (ql(i,j,k) .gt. qmin) then

                       reff_liq_local = k_ratio(i,j)* 620350.49 *    &
                          (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/Rdgas/   &
                          TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       
                       if ( (k .eq. 1    .and. qa(i,j,2)      .lt. qamin) &
                            .or. &
                            (k .eq. KDIM .and. qa(i,j,KDIM-1) .lt. qamin) &
                            .or. &
                            (k .gt. 1 .and. k .lt. KDIM .and. &
                            qa(i,j,k-1) .lt. qamin .and. &
                            qa(i,j,k+1) .lt. qamin) ) then
                          reff_liq_local = 1.134 * reff_liq_local
                       end if

                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       reff_liq_local = MIN(16.6,reff_liq_local)
                       reff_liq_local = MAX(4.2, reff_liq_local)

                       if (lhsw) Reff_liq(i,j,nclds(i,j)) =  reff_liq_local
                       if (sea_esf) size_drop_org(i,j,k) = 2. * reff_liq_local

                    end if  !ql calculation

                    !compute reff_ice
                    if (qi(i,j,k) .gt. qmin) then
                          
                       if (lhsw) then
                       
                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                              reff_ice_local = 92.46298
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                              reff_ice_local = 72.35392
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                              reff_ice_local = 85.19071 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                              reff_ice_local = 55.65818
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                              reff_ice_local = 35.29989
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                              reff_ice_local = 32.89967
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                              reff_ice_local = 16.60895
                          else
                              reff_ice_local = 15.41627
                          end if

                          !limit values to that for which Ebert and
                          !Curry radiation is valid :
                          !  10 microns < reff < 130 microns
                          !
                          reff_ice_local = MIN(130.,reff_ice_local)
                          Reff_ice(i,j,nclds(i,j)) = MAX(10.,reff_ice_local)                  
                       end if  !end of lhsw if

                       if (sea_esf) then

                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                              reff_ice_local = 100.6
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                              reff_ice_local = 80.8
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                              reff_ice_local = 93.5 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                              reff_ice_local = 63.9
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                              reff_ice_local = 42.5
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                              reff_ice_local = 39.9
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                              reff_ice_local = 21.6
                          else
                              reff_ice_local = 20.2
                          end if

                          !the ice crystal effective size can          
                          !only be 18.6 <= D^sub^e <= 130.2 microns. 
                          !for Fu Liou JAS 1993 code                    
                          reff_ice_local = MIN(130.2,reff_ice_local)
                          size_ice_org(i,j,k) = &
                                           MAX(18.6,reff_ice_local)
                        

                       end if  !end of sea_esf if
                             
                    end if !qi loop                       
                          
                 END IF  !cloud exist if

           END DO
           END DO
           END DO



    end if   !overlap = 2  if

    if (present(lwp_in)) then
      lwp_in = lwp
      iwp_in = iwp
      reff_liq_in = reff_liq
      reff_ice_in = reff_ice
     endif
    
      deallocate (lwp, iwp, reff_liq, reff_ice)

END SUBROUTINE CLOUD_ORGANIZE
SUBROUTINE CLOUD_OPTICAL_PROPERTIES(LWP,IWP,Reff_liq,Reff_ice,&
                        tau,w0,gg,em_lw,tau_ice_diag)

                              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau   :optical depth in each band
!               2. w0    :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw    :longwave cloud emmissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation can be rewritten as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!   Which is what is solved here.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!VARIABLES
!
!       ------
!INPUT:
!       ------
!
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!
!       ------
!INPUT/OUTPUT:
!       ------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emmissivity
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       tau_ice_diag    optical depth in each band
!
!
!       -------------------
!INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq   optical depth            at each band for cloud liquid
!       tau_ice   optical depth            at each band for cloud ice
!       w0_liq    single scattering albedo at each band for cloud liquid
!       w0_ice    single scattering albedo at each band for cloud ice
!       g_liq     asymmetry parameter      at each band for cloud liquid
!       g_ice     asymmetry parameter      at each band for cloud ice
!       k_liq        liquid cloud mass absorption coefficient for longwave
!                         portion of the spectrum (meters**2./kg condensate)
!       k_ice           ice cloud mass absorption coefficient for longwave
!                         portion of the spectrum (meters**2./kg condensate)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------


REAL,     INTENT (IN)   ,DIMENSION(:,:,:)   :: LWP,IWP,Reff_liq,Reff_ice
REAL,     INTENT (INOUT),DIMENSION(:,:,:,:) :: tau,w0,gg
REAL,     INTENT (INOUT),DIMENSION(:,:,:)   :: em_lw
REAL,     INTENT (INOUT), OPTIONAL, DIMENSION(:,:,:,:) :: tau_ice_diag
        
!  Internal variables
!  ------------------
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: tau_liq,tau_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: w0_liq,w0_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: g_liq,g_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3))   :: k_liq,k_ice

!
! Code
! ----
        
        ! reinitialize output variables to default values
        ! (not usually used)
        tau(:,:,:,:)=0.
        gg(:,:,:,:) = 0.85
        w0(:,:,:,:) = 0.95
        em_lw(:,:,:) = 0.
        
        ! reinitialize internal variables (not usually used)
        w0_liq(:,:,:,:) = 0.95
        w0_ice(:,:,:,:) = 0.95
        g_liq(:,:,:,:)  = 0.85
        g_ice(:,:,:,:)  = 0.85
        tau_liq(:,:,:,:)= 0.
        tau_ice(:,:,:,:)= 0.



   !---------------   COMPUTE OPTICAL DEPTH ---------------------------!

        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately

        tau_liq(:,:,:,1) = LWP(:,:,:) * 1000. * &
                           (0.02817 + (1.305/Reff_liq(:,:,:)))
        tau_liq(:,:,:,2) = LWP(:,:,:) * 1000. * &
                           (0.02682 + (1.346/Reff_liq(:,:,:)))
        tau_liq(:,:,:,3) = LWP(:,:,:) * 1000. * &
                           (0.02264 + (1.454/Reff_liq(:,:,:)))
        tau_liq(:,:,:,4) = LWP(:,:,:) * 1000. * &
                           (0.01281 + (1.641/Reff_liq(:,:,:)))
        
        tau_ice(:,:,:,1) = IWP(:,:,:) * 1000. * &
                           (0.003448 + (2.431/Reff_ice(:,:,:)))
        tau_ice(:,:,:,2) = tau_ice(:,:,:,1)
        tau_ice(:,:,:,3) = tau_ice(:,:,:,1)
        tau_ice(:,:,:,4) = tau_ice(:,:,:,1)
        

        ! compute total cloud optical depth
        tau(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)
        

   !---------------   COMPUTE SINGLE SCATTERING ALBEDO ----------------!

        w0_liq(:,:,:,1) =  5.62E-08   - 1.63E-07*Reff_liq(:,:,:)
        w0_liq(:,:,:,2) =  6.94E-06   - 2.35E-05*Reff_liq(:,:,:)
        w0_liq(:,:,:,3) = -4.64E-04   - 1.24E-03*Reff_liq(:,:,:)
        w0_liq(:,:,:,4) = -2.01E-01   - 7.56E-03*Reff_liq(:,:,:)
        w0_liq(:,:,:,:) = w0_liq(:,:,:,:) + 1.

        w0_ice(:,:,:,1) = -1.00E-05
        w0_ice(:,:,:,2) = -1.10E-04   - 1.41E-05*Reff_ice(:,:,:)
        w0_ice(:,:,:,3) = -1.86E-02   - 8.33E-04*Reff_ice(:,:,:)
        w0_ice(:,:,:,4) = -4.67E-01   - 2.05E-05*Reff_ice(:,:,:)
        w0_ice(:,:,:,:) = w0_ice(:,:,:,:) + 1.


        ! compute total single scattering albedo
        WHERE (tau(:,:,:,:) .gt. 0.)
               w0(:,:,:,:) = ( w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                               w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )  /&
                             tau(:,:,:,:)
        END WHERE
        
   !---------------   COMPUTE ASYMMETRY PARAMETER --------------------!


       g_liq(:,:,:,1) = 0.829 + 2.482E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,2) = 0.794 + 4.226E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,3) = 0.754 + 6.560E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,4) = 0.826 + 4.353E-03*Reff_liq(:,:,:)

       g_ice(:,:,:,1) = 0.7661+ 5.851E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,2) = 0.7730+ 5.665E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,3) = 0.7940+ 7.267E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,4) = 0.9595+ 1.076E-04*Reff_ice(:,:,:)

        ! compute  asymmetry parameter
        WHERE (tau(:,:,:,:) .gt. 0. )
              gg(:,:,:,:) = ( &
                 w0_liq(:,:,:,:) * g_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                 w0_ice(:,:,:,:) * g_ice(:,:,:,:) * tau_ice(:,:,:,:) ) &
                       /          (w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                                   w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )
        END WHERE

        
   !---------------   COMPUTE LONGWAVE EMMISSIVITY --------------------!


        k_liq(:,:,:) = 140.
        k_ice(:,:,:) = 4.83591 + 1758.511/Reff_ice(:,:,:)
        
        ! compute combined emmisivity
        em_lw(:,:,:) =  1. - exp( -1. * ( k_liq(:,:,:) * LWP(:,:,:) + &
                                          k_ice(:,:,:) * IWP(:,:,:) ) )

        
   !--------------    RANGE LIMIT QUANTITIES --------------------------!

        WHERE (tau(:,:,:,:) .lt. taumin)
               tau(:,:,:,:) = taumin
        END WHERE

   !----- ---------    EVALUATE TAU_ICE_DIAG (OPTIONAL) ---------------!

        
        if (present(tau_ice_diag)) then
               tau_ice_diag(:,:,:,:) = 0.
               tau_ice_diag(:,:,:,:) = tau_ice(:,:,:,:)
        end if
 

END SUBROUTINE CLOUD_OPTICAL_PROPERTIES


!########################################################################

SUBROUTINE snow_and_rain (qa, pfull, phalf, tkel, cldamt, snow, rain,  &
                          size_snow_in, size_rain_in, conc_rain,  &
                          conc_snow,  size_rain, size_snow )

!----------------------------------------------------------------------
! size_snow currently not used -- 2/10
real, dimension(:,:,:), intent(in)     :: qa, pfull, phalf, tkel
real, dimension(:,:,:), intent(in)     :: snow, rain, size_snow_in,  &
                                          size_rain_in
real, dimension(:,:,:), intent(inout)  :: cldamt
real, dimension(:,:,:), intent(out)    :: conc_rain, conc_snow, &
                                          size_rain, size_snow

      REAL, DIMENSION( size(qa,1), size(qa,2), size(qa,3))  :: cldmax_loc
      REAL :: dum
      integer :: i,j,k


      conc_snow = 0.
      conc_rain = 0.
      size_snow = 1.e-20
      size_rain = 1.e-20

 

      IF (snow_in_cloudrad .OR. rain_in_cloudrad ) THEN

!--------------------------------------------------------------------
!    define the locally seen cloud max above each level.
!--------------------------------------------------------------------
        cldmax_loc = 0.
        do k=1,size(qa,3)
          do j=1,size(qa,2)
            do i=1,size(qa,1)
! this might be o.k. as long as stochastic clouds are used  
              if (k.eq.1) then
                cldmax_loc(i,j,k)  = qa(i,j,k) !max overlap for precip
              else
                cldmax_loc(i,j,k) = max(cldmax_loc(i,j,k-1), qa(i,j,k))
              end if
            end do
          end do
        end do

!-----------------------------------------------------------------------
!    define the snow field to be seen by the radiation code.
!-----------------------------------------------------------------------
        IF (snow_in_cloudrad) THEN
          do k=1,size(qa,3) 
            do j=1,size(qa,2)
              do i=1,size(qa,1)
                IF (cldmax_loc(i,j,k) .GE. qamin) THEN
                  dum =  snow(i,j,k)/cldmax_loc(i,j,k)
                  dum = MIN(dum, 60.e-3)
                  IF (dum .GE. qmin) THEN
                    conc_snow (i,j,k) =     &
                       1000.*dum*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                        RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                        MAX(phalf(i,j,k), pfull(i,j,1)))
!size_snow currently not used
!!                  size_snow (i,j,k) =  MAX( 1.e-20, size_snow_in( i,j,k))
!!                  cldamt(i,j,k) = cldmax_loc(i,j,k)
                  END IF
                END IF 
              end do
            end do
          end do
        END IF  

!-----------------------------------------------------------------------
!    define the rain field to be seen by the radiation code.
!-----------------------------------------------------------------------
        IF (rain_in_cloudrad) THEN
          do k=1,size(qa,3)
            do j=1,size(qa,2)
              do i=1,size(qa,1)
                IF (cldmax_loc(i,j,k) .GE. qamin) THEN
                  dum = rain(i,j,k)/cldmax_loc(i,j,k)
                  dum = MIN(dum, 60.e-3)
                  IF (dum .GE. qmin) THEN
                    conc_rain (i,j,k) =     &
                       1000.*dum*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                       RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                       MAX(phalf(i,j,k), pfull(i,j,1)))
                    size_rain (i,j,k) = MAX( 1.e-20, size_rain_in(i,j,k))
!!                  cldamt(i,j,k) = cldmax_loc(i,j,k)
                  END IF
                END IF 
              end do
            end do
          end do
        END IF 
      END IF !(snow_in_cloudrad .OR. rain_in_cloudrad ) 

!-----------------------------------------------------------------------

END SUBROUTINE snow_and_rain


!########################################################################
      
                  end module cloud_rad_mod

      


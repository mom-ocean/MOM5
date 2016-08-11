                 module aerosolrad_package_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <OVERVIEW>
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!    shared modules:

use mpp_mod,               only: input_nml_file
use fms_mod,               only: open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, close_file
use mpp_io_mod,            only: mpp_open, mpp_close, MPP_RDONLY,   &
                                 MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, &
                                 MPP_SINGLE, mpp_io_init
use time_manager_mod,      only: time_type, time_manager_init,  &
                                 get_date, set_date, operator(+), &
                                 print_date, operator(-), operator(>)
use diag_manager_mod,      only: diag_manager_init, get_base_time
use interpolator_mod,      only: interpolate_type, interpolator_init, &
                                 interpolator, interpolator_end, &
                                 obtain_interpolator_time_slices, &
                                 unset_interpolator_time_flag, &
                                 CONSTANT, INTERP_WEIGHTED_P

! shared radiation package modules:
                                
use rad_utilities_mod,     only: Sw_control, &
                                 Lw_control, &
                                 Rad_control,&
                                 aerosol_type, aerosol_properties_type,&
                                 aerosol_diagnostics_type, &
                                 Lw_parameters, rad_utilities_init, &
                                 thickavg
use esfsw_parameters_mod,  only: Solar_spect, esfsw_parameters_init 
use longwave_params_mod,   only: NBLW, longwave_params_init

!-------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: aerosolrad_package.F90,v 20.0 2013/12/13 23:18:56 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public           &
       aerosolrad_package_init, aerosol_radiative_properties, &
       aerosolrad_package_alloc, aerosolrad_package_time_vary, &
       aerosolrad_package_endts, &
       aerosolrad_package_end,  get_aerosol_optical_info, &
       get_aerosol_optical_index
     
private          &

!  called from aerosolrad_package_init: 
   assign_aerosol_opt_props, read_optical_input_file, &
   sw_aerosol_interaction,  lw_aerosol_interaction
       

!---------------------------------------------------------------------
!-------- namelist  ---------

integer, parameter   ::        &
             MAX_OPTICAL_FIELDS = 1000  ! maximum number of aerosol 
                                       ! optical property types
logical              ::        &
             do_lwaerosol = .false.    ! aerosol efects included in lw
                                       ! radiation ?
logical              ::        &
             do_swaerosol = .false.    ! aerosol effects included in sw
                                       ! radiation ?
logical              ::        &
             force_to_repro_quebec = .false.
                                       ! if true, code sequence is 
                                       ! executed which reproduces
                                       ! quebec+ answers for 
                                       ! AM3p8e in pre_Riga
character(len=48)    ::        &
             aerosol_data_set = ' '    ! source of aerosol data; if 
                                       ! aerosols not desired remains
                                       ! ' ', otherwise is set to either
                                       ! 'shettle_fenn' or 
                                       ! 'Ginoux_Reddy'

!----------------------------------------------------------------------
!    the avaialable aerosol datasets are :
!    1) "shettle_fenn":  
!        Ref: shettle, e.p. and r.w. fenn, models for the aerosols of 
!            the lower atmosphere and the effects of humidity variations
!            on their optical properties,afgl-tr-79-0214,1979,94pp.    
!    2) "Ginoux_Reddy": 3D Aerosol fields are generated online reflecting
!        emissions, transport and deposition:
!        Ref: Reddy et al., 2005, Ginoux et al., 2005
!        
!----------------------------------------------------------------------

character(len=64)    ::        &
             aerosol_optical_names(MAX_OPTICAL_FIELDS) = '  '
                                       ! names associated with the 
                                       ! optical property types that
                                       ! are to be used in this 
                                       ! experiment
character(len=64)    ::        &
             optical_filename = ' '    ! name of file containing the
                                       ! aerosol optical property types
logical              ::        &
             using_volcanic_sw_files = .false.
                                       ! files containing sw aerosol
                                       ! optical properties from vol-
                                       ! canic activity are to be
                                       ! used to supplement those cal-
                                       ! culated by model ?
logical              ::        &
             using_volcanic_lw_files = .false.
                                       ! files containing lw aerosol
                                       ! optical properties from vol-
                                       ! canic activity are to be
                                       ! used to supplement those cal-
                                       ! culated by model ?
character(len=64)    ::        &
              sw_ext_filename = ' '    ! name of file containing the
                                       ! aerosol sw extinction optical
                                       ! depth
character(len=64)    ::        &
              sw_ssa_filename = ' '    ! name of file containing the
                                       ! aerosol sw single scattering 
                                       ! albedo
character(len=64)    ::        &
              sw_asy_filename = ' '    ! name of file containing the
                                       ! aerosol sw asymmetry factor   
character(len=64)    ::        &
              lw_ext_filename = ' '    ! name of file containing the
                                       ! aerosol lw extinction optical
                                       ! depth
character(len=64)    ::        &
              lw_ssa_filename = ' '    ! name of file containing the
                                       ! aerosol lw single scattering 
                                       ! albedo
character(len=64)    ::        &
              lw_asy_filename = ' '    ! name of file containing the
                                       ! aerosol lw asymmetry factor   
                                       ! the supplemental input files
character(len=64)    ::        &
              sw_ext_root                  = '   ' 
                                       ! names given to sw extopdep in
                                       ! input netcdf file
character(len=64)    ::        &
              sw_ssa_root                  = '   ' 
                                       ! name given to sw single scat-
                                       ! tering albedo in input netcdf 
                                       ! file
character(len=64)    ::        &
              sw_asy_root                  = '   ' 
                                       ! name given to sw asymmetry
                                       ! factor in input netcdf file
character(len=64)    ::        &
              lw_ext_root                  = '   ' 
                                       ! name given to lw extopdep in
                                       ! input netcdf file
character(len=64)    ::        &
              lw_ssa_root                  = '   '  
                                       ! name given to lw single scat-
                                       ! tering albedo in input netcdf 
                                       ! file
character(len=64)    ::        &
              lw_asy_root                  = '   '      
                                       ! name given to lw asymmetry
                                       ! factor in input netcdf file
integer, dimension(6) ::       &
              volcanic_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /) 
                                       ! time in volcanic data set
                                       ! corresponding to model
                                       ! initial time 
                                       ! (yr, mo, dy, hr, mn, sc)
logical :: interpolating_volcanic_data = .true.
                                       ! volcanic datasets will be
                                       ! time interpolated rather than
                                       ! held constant for a month ?
logical :: repeat_volcano_year = .false. 
                                      ! the same single year's data from
                                      ! the input data set should be 
                                      ! used for each model year ?
integer :: volcano_year_used = 0      ! year of volcanic data to repeat
                                      ! when repeat_volcano_year is
                                      ! .true.
logical :: using_im_bcsul = .false.   ! bc and sulfate aerosols are 
                                      ! treated as an internal mixture ?
integer, dimension(0:100) ::  omphilic_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  bcphilic_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt1_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt2_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt3_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt4_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  seasalt5_indices = (/        &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
                             0  /)
integer, dimension(0:100) ::  sulfate_indices = (/      &
                           30,30,30,30,30,30,30,30,30,30,30,30,30,30, &
                           30,30,30,30,30,30,30,30,30,30,30,30,30,30, &
                           30,30,30,30,30,35,35,35,35,35,40,40,40,40, &
                           40,45,45,45,45,45,50,50,50,50,50,55,55,55, &
                           55,55,60,60,60,60,60,65,65,65,65,65,70,70, &
                           70,70,70,75,75,75,75,75,80,80,80,80,82,82, &
                           84,84,86,86,88,88,90,91,92,93,94,95,96,97, &
                           98,99,100 /)
!yim
integer, dimension(0:100) ::  sulfate_vol_indices = (/      &
                             100,98,98,96,96,94,94,92,92,90,90,88,88,86,86,84,84,82,82,80, &
                             80,80,80,75,75,75,75,75,70,70,70,70,70,65,65,65, &
                             65,65,60,60,60,60,60,55,55,55,55,55,50,50,50,50,50,45,45,45, &
                             45,45,40,40,40,40,40,35,35,35,35,35,30,30,30,30,30,25,25,25, &
                             25,25,20,20,20,20,20,15,15,15,15,15,10,10,10,10,10,5,5,5, &
                             5,5,0,0,0  /)


namelist / aerosolrad_package_nml /                          &
                                    do_lwaerosol, do_swaerosol, &
                                    force_to_repro_quebec, &
                                    aerosol_data_set, &
                                    aerosol_optical_names, &
                                    sulfate_indices, &
                                    sulfate_vol_indices, &
                                    omphilic_indices, &
                                    bcphilic_indices, &
                                    seasalt1_indices, &
                                    seasalt2_indices, &
                                    seasalt3_indices, &
                                    seasalt4_indices, &
                                    seasalt5_indices, &
                                    optical_filename   , &
                                    using_volcanic_sw_files, &
                                    using_volcanic_lw_files, &
                                    volcanic_dataset_entry, &
                                    interpolating_volcanic_data, &
                                    repeat_volcano_year, &
                                    volcano_year_used, &
                                    using_im_bcsul, &
                                    sw_ext_filename, sw_ssa_filename, &
                                    sw_asy_filename, lw_ext_filename, &
                                    lw_ssa_filename, lw_asy_filename, &
                                    sw_ext_root, sw_ssa_root,   &
                                    sw_asy_root, lw_ext_root,   &
                                    lw_ssa_root, lw_asy_root

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!   the following are interpolate_type variables containing the
!   additional aerosol optical properties that may be included as
!   input to the radiation package.
!---------------------------------------------------------------------
type(interpolate_type), save  :: Sw_aer_extopdep_interp
type(interpolate_type), save  :: Sw_aer_ssalb_interp
type(interpolate_type), save  :: Sw_aer_asymm_interp
type(interpolate_type), save  :: Lw_aer_extopdep_interp
type(interpolate_type), save  :: Lw_aer_ssalb_interp
type(interpolate_type), save  :: Lw_aer_asymm_interp

!---------------------------------------------------------------------
!    the following variables define the number and type of different
!    bands over which the radiation package calculates aerosol 
!    radiative properties.
!---------------------------------------------------------------------
integer, parameter ::    &
         N_AEROSOL_BANDS_FR = 8 ! number of non-continuum ir aerosol
                                ! emissivity bands 
integer, parameter ::     &
         N_AEROSOL_BANDS_CO = 1 ! number of continuum ir aerosol
                                ! emissivity bands  
integer, parameter ::     &
         N_AEROSOL_BANDS_CN = 1 ! number of diagnostic continuum ir 
                                ! aerosol emissivity bands  
integer, parameter ::    &
         N_AEROSOL_BANDS = N_AEROSOL_BANDS_FR + N_AEROSOL_BANDS_CO
                                ! total number of ir aerosol emissivity
                                ! bands 

!--------------------------------------------------------------------
!    num_wavenumbers is the number of wavenumber bands over which the
!    aerosol parameterization provides aerosol radiative property data.
!--------------------------------------------------------------------
integer     ::  num_wavenumbers = 0 ! number of wavenumber bands 
                                    ! present in the aerosol 
                                    ! parameterization

!----------------------------------------------------------------------
!    the following variable defines the number of aerosol property 
!    types that are active.
!----------------------------------------------------------------------
integer     :: naermodels = 0   ! number of aerosol optical properties
                                ! types that are active

!---------------------------------------------------------------------
!    flags indicating an index value characteristic of the optical prop-
!    erties associated with different aerosols
!---------------------------------------------------------------------
integer, PARAMETER ::   SULFATE_FLAG =  0
integer, PARAMETER ::  OMPHILIC_FLAG = -1
integer, PARAMETER ::  BCPHILIC_FLAG = -2
integer, PARAMETER ::  SEASALT1_FLAG = -3
integer, PARAMETER ::  SEASALT2_FLAG = -4
integer, PARAMETER ::  SEASALT3_FLAG = -5
integer, PARAMETER ::  SEASALT4_FLAG = -6
integer, PARAMETER ::  SEASALT5_FLAG = -7
!yim
integer, PARAMETER ::  BC_FLAG = -8
integer, PARAMETER ::  NOT_IN_USE = -2000

!----------------------------------------------------------------------
!    the following index arrays contain the mapping information between 
!    actual model relative humidity and the available enties in the 
!    aerosol optical properties file.
!----------------------------------------------------------------------
integer, dimension(:,:), allocatable :: sulfate_index
integer, dimension(:),   allocatable :: optical_index
integer, dimension(:),   allocatable :: omphilic_index
integer, dimension(:),   allocatable :: bcphilic_index
integer, dimension(:),   allocatable :: seasalt1_index
integer, dimension(:),   allocatable :: seasalt2_index
integer, dimension(:),   allocatable :: seasalt3_index
integer, dimension(:),   allocatable :: seasalt4_index
integer, dimension(:),   allocatable :: seasalt5_index

!---------------------------------------------------------------------
!    the following arrays related to sw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!    here n refers to the bands of the solar parameterization, ni
!    to the bands of the aerosol parameterization, and na to the optical
!    properties type.
!
!      solivlaero(n,ni)  amount of toa incoming solar from solar
!                        spectral band n that is in aerosol parameter-
!                        ization band ni
!      nivl1aero(n)      the aerosol band index corresponding to the 
!                        lowest wave number of spectral band n
!      nivl2aero(n)      the aerosol band index corresponding to the 
!                        highest wave number of spectral band n
!      endaerwvnsf(ni)   ending wave number of aerosol parameterization
!                        band ni
!      aeroextivl(ni,na) extinction coefficient for aerosol parameter-
!                        ization band ni for aerosol optical property 
!                        type na
!      aerossalbivl(ni,na) 
!                        single-scattering albedo for aerosol band 
!                        ni and aerosol optical property type na
!      aeroasymmivl(ni,na)
!                        asymmetry factor for aerosol band ni  and 
!                        aerosol optical property type na
!
!---------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: solivlaero  
integer, dimension(:),   allocatable   :: nivl1aero, nivl2aero
integer, dimension(:),   allocatable   :: endaerwvnsf
real,    dimension(:,:), allocatable   :: aeroextivl, aerossalbivl, &
                                          aeroasymmivl

!---------------------------------------------------------------------
!    sfl following arrays related to lw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!
!    sflwwts(n,ni)     the fraction of the planck function in aerosol 
!                      emissivity band n that is in aerosol param-
!                      eterization band ni
!
!----------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: sflwwts, sflwwts_cn

!--------------------------------------------------------------------
!    logical flags 
!--------------------------------------------------------------------
logical :: module_is_initialized      = .false. ! module has been
                                                ! initialized ?
!logical :: doing_predicted_aerosols   = .false. ! predicted aerosol 
                                                ! scheme being used ?
logical :: band_calculation_completed = .false. ! lw properties have
                                                ! been calculated ?

type(time_type) :: Model_init_time  ! initial calendar time for model  
                                    ! [ time_type ]
type(time_type) :: Volcanic_offset  ! difference between model initial
                                    ! time and volcanic timeseries app-
                                    ! lied at model initial time
                                    ! [ time_type ]
type(time_type) :: Volcanic_entry   ! time in volcanic timeseries which
                                    ! is mapped to model initial time
                                    ! [ time_type ]
logical    :: negative_offset = .false.
                                !  the model initial time is later than
                                !  the volcanic_dataset_entry time  ?
integer :: nfields_sw_ext = 0   ! number of fields contained in 
                                ! supplemental sw_ext file
integer :: nfields_sw_ssa = 0   ! number of fields contained in 
                                ! supplemental sw_ssa file
integer :: nfields_sw_asy = 0   ! number of fields contained in 
                                ! supplemental sw_asy file
integer :: nfields_lw_ext = 0   ! number of fields contained in 
                                ! supplemental lw_ext file
integer :: nfields_lw_ssa = 0   ! number of fields contained in 
                                ! supplemental lw_ssa file
integer :: nfields_lw_asy = 0   ! number of fields contained in 
                                ! supplemental lw_asy file

!-------------------------------------------------------------------
!   arrays holding variable names:
character(len=64), dimension(:), allocatable ::   &
                                sw_ext_name, sw_ssa_name, sw_asy_name, &
                                lw_ext_name, lw_ssa_name, lw_asy_name

!-------------------------------------------------------------------
!    arrays to hold data when not interpolating on every step:
real, dimension(:,:,:,:), allocatable :: sw_ext_save
real, dimension(:,:,:,:), allocatable :: sw_ssa_save
real, dimension(:,:,:,:), allocatable :: sw_asy_save
real, dimension(:,:,:,:), allocatable :: lw_ext_save
real, dimension(:,:,:,:), allocatable :: lw_ssa_save
real, dimension(:,:,:,:), allocatable :: lw_asy_save

!---------------------------------------------------------------------
!   module variables to hold values unchanging in time:
!---------------------------------------------------------------------
real, dimension(:,:), allocatable :: aerextband_MOD
real, dimension(:,:), allocatable :: aerssalbband_MOD
real, dimension(:,:), allocatable :: aerasymmband_MOD
real, dimension(:,:), allocatable :: aerextbandlw_MOD
real, dimension(:,:), allocatable :: aerssalbbandlw_MOD
real, dimension(:,:), allocatable :: aerextbandlw_cn_MOD
real, dimension(:,:), allocatable :: aerssalbbandlw_cn_MOD

!---------------------------------------------------------------------
!    logical variables indicating whether interpolation is currently
!    needed:
logical :: need_sw_ext = .true.
logical :: need_sw_ssa = .true.
logical :: need_sw_asy = .true.
logical :: need_lw_ext = .true.
logical :: need_lw_ssa = .true.
logical :: need_lw_asy = .true.

!---------------------------------------------------------------------
!    logical variables indicating whether the particular radiative 
!    property associated with volcanoes is being supplied:
logical :: using_sw_ext = .false. 
logical :: using_sw_ssa = .false. 
logical :: using_sw_asy = .false. 
logical :: using_lw_ext = .false. 
logical :: using_lw_ssa = .false. 
logical :: using_lw_asy = .false. 

!---------------------------------------------------------------------
!    counters associated with determining when interpolation needs to
!    be done:
logical :: mo_save_set = .false.
integer :: mo_save = 0
integer :: mo_new

type(time_type) :: Volcano_time
integer :: nfields_save
integer :: num_sul, num_bc
integer, dimension(:), allocatable :: sul_ind, bc_ind

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 


                         contains
 

 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_init">
!  <OVERVIEW>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_init (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosolrad_package_init (kmax, aerosol_names, lonb, latb)

!---------------------------------------------------------------------
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!character(len=64), dimension(:), intent(in)  :: aerosol_names
integer,                        intent(in)  :: kmax
character(len=*), dimension(:), intent(in)  :: aerosol_names
real, dimension(:,:),           intent(in)  :: lonb,latb



!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      kmax              number of model levels
!      aerosol_names     the names assigned to each of the activated
!                        aerosol species
!       lonb           2d array of model longitudes at cell corners
!                      [ radians ]
!       latb           2d array of model latitudes at cell corners
!                      [ radians ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer        :: unit, ierr, io, logunit
      integer        :: n
      character(len=16) :: chvers
      character(len=4)  :: chyr  

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
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
      call mpp_io_init
      call fms_init
      call diag_manager_init
      call time_manager_init
      call rad_utilities_init
      call esfsw_parameters_init
      call longwave_params_init

       nfields_save = size(aerosol_names(:))

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=aerosolrad_package_nml, iostat=io)
      ierr = check_nml_error(io,'aerosolrad_package_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosolrad_package_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosolrad_package_nml')
        end do
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=aerosolrad_package_nml)

!---------------------------------------------------------------------
!   exit if aerosols are desired with the lacis-hansen parameterization.
!---------------------------------------------------------------------
     if (Sw_control%do_lhsw_iz) then
       if (Sw_control%do_lhsw .and. do_swaerosol) then
         call error_mesg ('aerosolrad_package_mod', &
         ' cannot activate sw aerosols with lhsw', FATAL)
        endif
     else
       call error_mesg ('aerosolrad_package_mod', &
        'Sw_control%do_lhsw not yet initialized', FATAL)
     endif

!----------------------------------------------------------------------
!    define control variables which indicate whether the impact of 
!    aerosols on radiation is to be included in the sw and lw rad-
!    iation calculations. define a control variable which will be true 
!    if aerosols are included in either the sw or the lw radiation 
!    (Rad_control%do_aerosol).
!----------------------------------------------------------------------
      Sw_control%do_swaerosol = do_swaerosol
      Lw_control%do_lwaerosol = do_lwaerosol
      if (Rad_control%do_lwaerosol_forcing_iz .and. &
          Rad_control%do_swaerosol_forcing_iz) then

      if (do_lwaerosol .or. do_swaerosol  .or.  &
          using_volcanic_lw_files .or. using_volcanic_sw_files .or. &
         Rad_control%do_lwaerosol_forcing .or.  &
         Rad_control%do_swaerosol_forcing)  then
        Rad_control%do_aerosol = .true.
      else
        Rad_control%do_aerosol = .false.
      endif
      else
        call error_mesg ('aerosolrad_package_mod', &
         ' using Rad_control%do_{l,s}waerosol_forcing  before it &
                                             &is defined', FATAL)
      endif

!--------------------------------------------------------------------
!    mark the just defined logicals as initialized.
!--------------------------------------------------------------------
      Sw_control%do_swaerosol_iz = .true.        
      Lw_control%do_lwaerosol_iz = .true.        
      Rad_control%do_aerosol_iz  = .true.
     
!----------------------------------------------------------------------
!    store the control variable indicating whether aerosol internal
!    mixture is being assumed.
!----------------------------------------------------------------------
      Rad_control%using_im_bcsul = using_im_bcsul
      Rad_control%using_im_bcsul_iz = .true.

!---------------------------------------------------------------------
!    exit if an aerosol_data_set is provided when do_aerosol is 
!    .false..
!---------------------------------------------------------------------
      if ( .not. Rad_control%do_aerosol .and.    &
           trim(aerosol_data_set) /= ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are not desired, aerosol_data_set '//&
            'must be set to "   "', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if no aerosol_data_set is provided when do_aerosol  
!    is .true..
!---------------------------------------------------------------------
      if ( Rad_control%do_aerosol .and.    &
           trim(aerosol_data_set) == ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are desired, aerosol_data_set '//&
            'must be non-blank', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if aerosol effects are desired but the aerosol input file
!    provided no aerosol fields.
!---------------------------------------------------------------------
      if (Rad_control%do_aerosol .and. size(aerosol_names(:)) == 0) then
        call error_mesg ('aerosolrad_package_mod', &
          ' aerosols desired  for radiation but no aerosol '//&
            'data_names supplied', FATAL)
      endif


!----------------------------------------------------------------------
!    if aerosol radiative effects are to be included, call 
!    assign_aerosol_opt_props to assign the proper aerosol 
!    properties type to each aerosol type. then call 
!    read_optical_input_file to read the optical input file contain-
!    ing the aerosol parameterization information and data.
!----------------------------------------------------------------------
      if (Rad_control%do_aerosol) then
        call assign_aerosol_opt_props (aerosol_names)
        call read_optical_input_file
      endif
 
!---------------------------------------------------------------------
!    if aerosol effects are to be included in the sw calculation,
!    call sw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using.
!--------------------------------------------------------------------
      if (do_swaerosol .or. Rad_control%do_swaerosol_forcing) then
        call sw_aerosol_interaction                 
      endif

!---------------------------------------------------------------------
!    if aerosol effects are to be included in the lw calculation,
!    call lw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using. if
!    they are not, indicate that this part of the code has been 
!    executed.
!---------------------------------------------------------------------
      if (do_lwaerosol .or. Rad_control%do_lwaerosol_forcing) then
        call lw_aerosol_interaction
      else
        Lw_parameters%n_lwaerosol_bands_iz = .true.
      endif

!---------------------------------------------------------------------
!    make sure consistent nml settings are present. Cannot use volcanic
!    aerosols unless model aerosols are also activated.
!---------------------------------------------------------------------
      if (.not. do_swaerosol  .and.   &
          using_volcanic_sw_files) then
        call error_mesg ('aerosolrad_package_mod', &
         'cant use sw volcanic aerosols without activating standard &
                                               & sw aerosols', FATAL)
      endif
      if (.not. do_lwaerosol  .and.   &
          using_volcanic_lw_files) then
        call error_mesg ('aerosolrad_package_mod', &
         'cant use lw volcanic aerosols without activating standard &
                                               & lw aerosols', FATAL)
      endif

!---------------------------------------------------------------------
!    set the volcanic control variables to .false. when the model 
!    aerosols are not active.
!---------------------------------------------------------------------
      Rad_control%volcanic_lw_aerosols = using_volcanic_lw_files
      Rad_control%volcanic_lw_aerosols_iz = .true.
      Rad_control%volcanic_sw_aerosols = using_volcanic_sw_files
      Rad_control%volcanic_sw_aerosols_iz = .true.

!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the volcanic files are
!    to be used; if not present, the dataset entry point is taken as
!    the model base_time, defined in the diag_table.
!---------------------------------------------------------------------
      Model_init_time = get_base_time()
      if (using_volcanic_sw_files .or.  &
          using_volcanic_lw_files) then
        if (volcanic_dataset_entry(1) == 1 .and. &
            volcanic_dataset_entry(2) == 1 .and. &
            volcanic_dataset_entry(3) == 1 .and. &
            volcanic_dataset_entry(4) == 0 .and. &
            volcanic_dataset_entry(5) == 0 .and. &
            volcanic_dataset_entry(6) == 0 ) then      
          Volcanic_entry = Model_init_time

!----------------------------------------------------------------------
!    define the offset from model base time  (defined in diag_table) 
!    to volcanic_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
      else
        Volcanic_entry  = set_date (volcanic_dataset_entry(1), &
                                    volcanic_dataset_entry(2), &
                                    volcanic_dataset_entry(3), &
                                    volcanic_dataset_entry(4), &
                                    volcanic_dataset_entry(5), &
                                    volcanic_dataset_entry(6))
       endif
     else
       Volcanic_entry  = set_date (volcanic_dataset_entry(1), &
                                   volcanic_dataset_entry(2), &
                                   volcanic_dataset_entry(3), &
                                   volcanic_dataset_entry(4), &
                                   volcanic_dataset_entry(5), &
                                   volcanic_dataset_entry(6))

      endif
      if (using_volcanic_sw_files .or.  &
          using_volcanic_lw_files) then
        if (repeat_volcano_year) then
          if (volcano_year_used == 0) then
            call error_mesg ('aerosolrad_package_init', &
              'valid year must be supplied when &
                       &repeat_volcano_year is .true.', FATAL)
          endif
        endif
        call print_date(Volcanic_entry , str='Data from volcano &
                                           &timeseries at time:')
        call print_date(Model_init_time , str='This data is mapped to &
                                                  &model time:')
        if (repeat_volcano_year) then
          write (chyr, '(i4)') volcano_year_used
          call error_mesg ('aerosolrad_package_init', &
           'volcanic data from dataset year '  // chyr // ' will be &
                                    &used for all model years.', NOTE)
        endif
      endif
      Volcanic_offset = Volcanic_entry - Model_init_time

      if (Model_init_time > Volcanic_entry) then
        negative_offset = .true.
      else
        negative_offset = .false.
      endif

!-----------------------------------------------------------------------
!    if desired, process the sw extinction coefficient file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
      if (using_volcanic_sw_files) then
        if (trim(sw_ext_root) /= ' '  ) then
          nfields_sw_ext = Solar_spect%nbands
          allocate (sw_ext_name (nfields_sw_ext))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_ext_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_ext) )
            sw_ext_save = 0.0
          endif
          do n=1, nfields_sw_ext
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_ext_name(n) = trim(sw_ext_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_ext_name(n) = trim(sw_ext_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_extopdep_interp,  &
                                  sw_ext_filename, lonb, latb,  &
                                  sw_ext_name(:nfields_sw_ext),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_sw_ext = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the sw single scattering albedo file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
        if (trim(sw_ssa_root) /= ' '  ) then
          nfields_sw_ssa = Solar_spect%nbands
          allocate (sw_ssa_name (nfields_sw_ssa))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_ssa_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_ssa) )
            sw_ssa_save = 0.0
          endif
          do n=1, nfields_sw_ssa
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_ssa_name(n) = trim(sw_ssa_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_ssa_name(n) = trim(sw_ssa_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_ssalb_interp,   &
                                  sw_ssa_filename,  lonb, latb,    &
                                  sw_ssa_name(:nfields_sw_ssa),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) ) 
          using_sw_ssa = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the sw asymmetry factor file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
        if (trim(sw_asy_root)    /= ' '  ) then
          nfields_sw_asy = Solar_spect%nbands
          allocate (sw_asy_name (nfields_sw_asy))
          if (.not. interpolating_volcanic_data) then
            allocate (sw_asy_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_sw_asy) )
            sw_asy_save = 0.0
          endif
          do n=1, nfields_sw_asy
            if (n<= 9) then
              write (chvers, '(i1)') n
             sw_asy_name(n) = trim(sw_asy_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              sw_asy_name(n) = trim(sw_asy_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Sw_aer_asymm_interp,   &
                                  sw_asy_filename, lonb, latb,   &
                                  sw_asy_name(:nfields_sw_asy),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_sw_asy = .true.
        endif
      endif

!-----------------------------------------------------------------------
!    if desired, process the lw extinction coefficient file. allocate 
!    space for and define the names of each variable. if not interpol-
!    ating the data, allocate an array to store it between timesteps. 
!    call interpolator_init to initialize the interpolation module for
!    the file.
!-----------------------------------------------------------------------
      if (using_volcanic_lw_files) then
        if (trim(lw_ext_root)    /= ' '  ) then
          nfields_lw_ext = N_AEROSOL_BANDS
          allocate (lw_ext_name (nfields_lw_ext))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_ext_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_ext) )
            lw_ext_save = 0.0
          endif
          do n=1, nfields_lw_ext
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_ext_name(n) = trim(lw_ext_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_ext_name(n) = trim(lw_ext_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_extopdep_interp,   &
                                  lw_ext_filename, lonb,  latb,  &
                                  lw_ext_name(:nfields_lw_ext),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )
          using_lw_ext= .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the lw single scattering albedo file.  it 
!    currently is not needed with the sea lw radiation package. allocate
!    space for and define the names of each variable. call 
!    interpolator_init to initialize the interpolation module for the 
!    file.
!-----------------------------------------------------------------------
        if (trim(lw_ssa_root)    /= ' '  ) then
          nfields_lw_ssa = N_AEROSOL_BANDS
          allocate (lw_ssa_name (nfields_lw_ssa))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_ssa_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_ssa) )
            lw_ssa_save = 0.0
          endif
          do n=1, nfields_lw_ssa
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_ssa_name(n) = trim(lw_ssa_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_ssa_name(n) = trim(lw_ssa_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_ssalb_interp,  &
                                  lw_ssa_filename, lonb, latb,  &
                                  lw_ssa_name(:nfields_lw_ssa),   &
                                  data_out_of_bounds=(/CONSTANT/), &
                                  vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_lw_ssa = .true.
        endif

!--------------------------------------------------------------------
!    if desired, process the lw asymmetry factor file.  it currently is
!    not needed with the sea lw radiation package. allocate space for 
!    and define the names of each variable. call interpolator_init to
!    initialize the interpolation module for the file.
!-----------------------------------------------------------------------
        if (trim(lw_asy_root)    /= ' '  ) then
          nfields_lw_asy = N_AEROSOL_BANDS
          allocate (lw_asy_name (nfields_lw_asy))
          if (.not. interpolating_volcanic_data) then
            allocate (lw_asy_save(size(lonb,1)-1, size(latb,2)-1, &
                                  kmax, nfields_lw_asy) )
            lw_asy_save = 0.0
          endif
          do n=1, nfields_lw_asy
            if (n<= 9) then
              write (chvers, '(i1)') n
             lw_asy_name(n) = trim(lw_asy_root) // '_b0' // trim(chvers)
            else if (n <= 99) then
              write (chvers, '(i2)') n
              lw_asy_name(n) = trim(lw_asy_root) // '_b' // trim(chvers)
            else 
              call error_mesg ('aerosolrad_package_mod', &
                  ' code only handles up to 100 fields', FATAL)
            endif
          end do
          call interpolator_init (Lw_aer_asymm_interp,   &
                                   lw_asy_filename,  lonb, latb,  &
                                   lw_asy_name(:nfields_lw_asy),   &
                                   data_out_of_bounds=(/CONSTANT/), &
                                   vert_interp=(/INTERP_WEIGHTED_P/) )  
          using_lw_asy = .true.
        endif
     endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


end subroutine aerosolrad_package_init


!##########################################################################

subroutine aerosolrad_package_time_vary (Time)

!-------------------------------------------------------------------------
!    aerosolrad_package_time_vary performs time-dependent, space-independent
!    caluclations for this module
!-------------------------------------------------------------------------

type(time_type),         intent(in)   :: Time


      integer  :: yr, mo, dy, hr, mn, sc
      integer  :: na, ni, nw
      type(aerosol_properties_type) :: Aerosol_props_tem
      
    
!---------------------------------------------------------------------
!    define the time for which the volcanic properties will be obtained.
!---------------------------------------------------------------------
        if (using_volcanic_sw_files .or.   &
            using_volcanic_lw_files) then
          if (negative_offset) then
             Volcano_time = Time - Volcanic_offset
          else 
             Volcano_time = Time + Volcanic_offset
          endif
          if (repeat_volcano_year) then
            call get_date (Volcano_time, yr, mo, dy, hr, mn, sc)
            Volcano_time = set_date (volcano_year_used, mo,dy,hr,mn,sc)
          endif

!--------------------------------------------------------------------
!    decide whether the volcanic data must be interpolated on this step.
!    if interpolating_volcanic_data is true, then all variables will
!    always be interpolated. when this is not .true., determine if the
!    month of the data desired has changed from the previous value. if
!    it has set the Volcano_time to 12Z on the 15th of the month, and
!    indicate that new data is needed. On the initial call of the job,
!    one always obtains the data (mo_save_set = .false.).
!--------------------------------------------------------------------
          if (interpolating_volcanic_data) then
            need_sw_ext = .true.
            need_sw_ssa = .true.
            need_sw_asy = .true.
            need_lw_ext = .true.
            need_lw_ssa = .true.
            need_lw_asy = .true.
          else
            call get_date (Volcano_time, yr,mo,dy,hr,mn,sc)
            Volcano_time =  set_date (yr, mo,15,12,0,0)
            if (mo_save_set) then
              if (mo /= mo_save) then
                mo_new = mo       
                need_sw_ext = .true.
                need_sw_ssa = .true.
                need_sw_asy = .true.
                need_lw_ext = .true.
                need_lw_ssa = .true.
                need_lw_asy = .true.
              endif
            else
              need_sw_ext = .true.
              need_sw_ssa = .true.
              need_sw_asy = .true.
              need_lw_ext = .true.
              need_lw_ssa = .true.
              need_lw_asy = .true.
            endif
          endif
        endif ! (using_volcanic_lw or using_volcanic_sw)

!--------------------------------------------------------------------
!    if the volcanic sw aerosol extinction is being supplied, make sure
!    needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_ext) then
          if (need_sw_ext) then
            if (nfields_sw_ext >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_extopdep_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol single scattering albedo is being 
!    supplied, make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_ssa) then
          if (need_sw_ssa) then
            if (nfields_sw_ssa >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_ssalb_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol asymmetry factor is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_sw_asy) then
          if (need_sw_asy) then
            if (nfields_sw_asy >= 1) then
              call obtain_interpolator_time_slices  &
                              (Sw_aer_asymm_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol extinction is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_lw_ext) then
          if (need_lw_ext) then
            if (nfields_lw_ext >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_extopdep_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw single scattering albedo is being supplied, 
!    make sure needed time slices are available.
!--------------------------------------------------------------------
        if (using_lw_ssa) then
          if (need_lw_ssa) then
            if (nfields_lw_ssa >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_ssalb_interp, Volcano_Time)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_asy) then
          if (need_lw_asy) then
            if (nfields_lw_asy >= 1) then
              call obtain_interpolator_time_slices  &
                              (Lw_aer_asymm_interp, Volcano_Time)
            endif
          endif
        endif

          if (force_to_repro_quebec) then
            if (.not. band_calculation_completed) then
              allocate (Aerosol_props_tem%aerextbandlw  &
                                       (N_AEROSOL_BANDS, naermodels))
              allocate (Aerosol_props_tem%aerssalbbandlw  &
                                       (N_AEROSOL_BANDS, naermodels))
              allocate (Aerosol_props_tem%aerextbandlw_cn &
                                       (N_AEROSOL_BANDS_CN, naermodels))
              allocate (Aerosol_props_tem%aerssalbbandlw_cn  &
                                       (N_AEROSOL_BANDS_CN, naermodels))
              aerextbandlw_MOD = 0.0                          
              aerssalbbandlw_MOD = 0.0
              aerextbandlw_cn_MOD = 0.0                              
              aerssalbbandlw_cn_MOD = 0.0
              Aerosol_props_tem%aerextbandlw = 0.0               
              Aerosol_props_tem%aerssalbbandlw = 0.0                 
              Aerosol_props_tem%aerextbandlw_cn = 0.0
              Aerosol_props_tem%aerssalbbandlw_cn = 0.0
              do nw=1,naermodels    
                do na=1,N_AEROSOL_BANDS  
                  do ni=1,num_wavenumbers 
                    Aerosol_props_tem%aerextbandlw(na,nw) =   &
                               Aerosol_props_tem%aerextbandlw(na,nw) + &
                                aeroextivl(ni,nw)*sflwwts(na,ni)*1.0E+03
                    Aerosol_props_tem%aerssalbbandlw(na,nw) =   &
                              Aerosol_props_tem%aerssalbbandlw(na,nw) +   &
                                     aerossalbivl(ni,nw)*sflwwts(na,ni)
                  end do
                end do
              end do
              do nw=1,naermodels    
                do na=1,N_AEROSOL_BANDS_CN
                  do ni=1,num_wavenumbers 
                    Aerosol_props_tem%aerextbandlw_cn(na,nw) = &
                         Aerosol_props_tem%aerextbandlw_cn(na,nw) + &
                            aeroextivl(ni,nw)*sflwwts_cn(na,ni)*1.0E+03
                    Aerosol_props_tem%aerssalbbandlw_cn(na,nw) =    &
                           Aerosol_props_tem%aerssalbbandlw_cn(na,nw) +&
                                  aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
                  end do
                end do
              end do

              aerextbandlw_MOD = Aerosol_props_tem%aerextbandlw
              aerssalbbandlw_MOD = Aerosol_props_tem%aerssalbbandlw
              aerextbandlw_cn_MOD = Aerosol_props_tem%aerextbandlw_cn
              aerssalbbandlw_cn_MOD = Aerosol_props_tem%aerssalbbandlw_cn
              band_calculation_completed = .true.
              deallocate (Aerosol_props_tem%aerextbandlw)
              deallocate (Aerosol_props_tem%aerssalbbandlw)
              deallocate (Aerosol_props_tem%aerextbandlw_cn)
              deallocate (Aerosol_props_tem%aerssalbbandlw_cn)
            endif
          endif

!---------------------------------------------------------------------------


end subroutine aerosolrad_package_time_vary 

!######################################################################

subroutine aerosolrad_package_alloc (ix, jx, kx, Aerosol_props)

integer,                       intent(in) :: ix, jx, kx
type(aerosol_properties_type), intent(inout) :: Aerosol_props


        if (Rad_control%volcanic_sw_aerosols) then
          allocate (Aerosol_props%sw_ext (ix,jx,kx, nfields_sw_ext))
          allocate (Aerosol_props%sw_ssa (ix,jx,kx,nfields_sw_ssa))
          allocate (Aerosol_props%sw_asy (ix,jx,kx,nfields_sw_asy))
        endif
        if (Rad_control%volcanic_lw_aerosols) then
          allocate (Aerosol_props%lw_ext (ix,jx,kx, nfields_lw_ext))
          allocate (Aerosol_props%lw_ssa (ix,jx,kx,nfields_lw_ssa))
          allocate (Aerosol_props%lw_asy (ix,jx,kx,nfields_lw_asy))
        endif
        allocate (Aerosol_props%ivol(ix,jx,kx))
        if (Rad_control%do_swaerosol_forcing .or.  &
                                         Sw_control%do_swaerosol) then
          allocate (Aerosol_props%aerextband   &
                                      (Solar_spect%nbands, naermodels))
          allocate (Aerosol_props%aerssalbband &
                                      (Solar_spect%nbands, naermodels))
          allocate (Aerosol_props%aerasymmband &
                                      (Solar_spect%nbands, naermodels))
        endif
        if (Rad_control%do_lwaerosol_forcing .or.  &
                                         Lw_control%do_lwaerosol) then
          allocate (Aerosol_props%aerextbandlw  &
                                      (N_AEROSOL_BANDS, naermodels))
          allocate (Aerosol_props%aerssalbbandlw  &
                                      (N_AEROSOL_BANDS, naermodels))
          allocate (Aerosol_props%aerextbandlw_cn &
                                      (N_AEROSOL_BANDS_CN, naermodels))
          allocate (Aerosol_props%aerssalbbandlw_cn  &
                                      (N_AEROSOL_BANDS_CN, naermodels))
        endif
        if (Rad_control%using_im_bcsul) then
          allocate (Aerosol_props%sulfate_index (0:100, 0:100))
        else
          allocate (Aerosol_props%sulfate_index (0:100, 0:0))
        endif
        allocate (Aerosol_props%optical_index (nfields_save))
        allocate (Aerosol_props%omphilic_index(0:100))
        allocate (Aerosol_props%bcphilic_index(0:100))
        allocate (Aerosol_props%seasalt1_index(0:100))
        allocate (Aerosol_props%seasalt2_index(0:100))
        allocate (Aerosol_props%seasalt3_index(0:100))
        allocate (Aerosol_props%seasalt4_index(0:100))
        allocate (Aerosol_props%seasalt5_index(0:100))

!--------------------------------------------------------------------


end subroutine aerosolrad_package_alloc


subroutine aerosolrad_package_endts

!---------------------------------------------------------------------
!    when data is not always interpolated, set flags to indicate whether
!    data must be obtained on the next call to this subroutine. if 
!    the current call has obtained data, set the flag indicating that 
!    data is not needed on the next call. also set the flag to indicate 
!    that the initial call has been completed (mo_save_set), and that 
!    the month for which data was obtained has been defined (mo_save). 
!---------------------------------------------------------------------
        if (.not. interpolating_volcanic_data) then
          if (need_sw_ext) then
              need_sw_ext = .false.
              need_sw_ssa = .false.
              need_sw_asy = .false.
              need_lw_ext = .false.
              mo_save_set = .true.
              mo_save = mo_new
          endif
        endif

      if (using_sw_ext) then
        call unset_interpolator_time_flag (Sw_aer_extopdep_interp)
      endif
      if (using_sw_ssa) then
        call unset_interpolator_time_flag (Sw_aer_ssalb_interp)
      endif
      if (using_sw_asy) then
        call unset_interpolator_time_flag (Sw_aer_asymm_interp)
      endif
      if (using_lw_ext) then
        call unset_interpolator_time_flag (Lw_aer_extopdep_interp)
      endif
      if (using_sw_ssa) then
        call unset_interpolator_time_flag (Lw_aer_ssalb_interp)
      endif
      if (using_sw_asy) then
        call unset_interpolator_time_flag (Lw_aer_asymm_interp)
      endif

end subroutine aerosolrad_package_endts


!####################################################################
! <SUBROUTINE NAME="aerosol_radiative_properties">
!  <OVERVIEW>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_radiative_properties (is, ie, js, je, &
!                                         Aerosol, Aerosol_props)
!  </TEMPLATE>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties in radiation package
!  </INOUT>
!  <IN NAME="is, ie" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js, je" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_radiative_properties (is, ie, js, je, &
                                         Time, p_half, Aerosol_diags, &
                                         Aerosol, Aerosol_props)

!---------------------------------------------------------------------
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!---------------------------------------------------------------------

integer,                       intent(in)    :: is, ie, js, je
type(time_type),               intent(in)    :: Time
real, dimension(:,:,:),        intent(in)    :: p_half
type(aerosol_type),            intent(in)    :: Aerosol
type(aerosol_diagnostics_type), intent(inout) :: Aerosol_diags
type(aerosol_properties_type), intent(inout) :: Aerosol_props
 
!----------------------------------------------------------------------
! local variables:                                                     

      integer  :: na, nw, ni       ! do-loop indices
      integer  :: iaer, i, j, k
      real, dimension (size(Aerosol%aerosol,1),  &
                       size(Aerosol%aerosol,2),  &
                       size(Aerosol%aerosol,3))  :: sul, bc
     
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate and initialize arrays to hold aerosol diagnostics.
!---------------------------------------------------------------------
      allocate (Aerosol_diags%extopdep (size(Aerosol%aerosol,1), &
                                        size(Aerosol%aerosol,2), &
                                        size(Aerosol%aerosol,3), &
                                        size(Aerosol%aerosol,4), 10 ))
      Aerosol_diags%extopdep = 0.0
      allocate (Aerosol_diags%absopdep (size(Aerosol%aerosol,1), &
                                        size(Aerosol%aerosol,2), &
                                        size(Aerosol%aerosol,3), &
                                        size(Aerosol%aerosol,4), 10 ))
      Aerosol_diags%absopdep = 0.0
      allocate (Aerosol_diags%asymdep (size(Aerosol%aerosol,1), &
                                       size(Aerosol%aerosol,2), &
                                       size(Aerosol%aerosol,3), &
                                       size(Aerosol%aerosol,4), 10 ))
      Aerosol_diags%asymdep = 0.0
      allocate (Aerosol_diags%extopdep_vlcno    &
                                        (size(Aerosol%aerosol,1), &
                                         size(Aerosol%aerosol,2), &
                                         size(Aerosol%aerosol,3),3))
      Aerosol_diags%extopdep_vlcno = 0.0
      allocate (Aerosol_diags%absopdep_vlcno  &
                                        (size(Aerosol%aerosol,1), &
                                         size(Aerosol%aerosol,2), &
                                         size(Aerosol%aerosol,3),3))
      Aerosol_diags%absopdep_vlcno = 0.0
      allocate (Aerosol_diags%sw_heating_vlcno  &
                                        (size(Aerosol%aerosol,1), &
                                         size(Aerosol%aerosol,2), &
                                         size(Aerosol%aerosol,3), &
                                         Rad_control%nzens))
      Aerosol_diags%sw_heating_vlcno = 0.0
      allocate (Aerosol_diags%lw_extopdep_vlcno    &
                                        (size(Aerosol%aerosol,1), &
                                         size(Aerosol%aerosol,2), &
                                         size(Aerosol%aerosol,3)+1,2))
      Aerosol_diags%lw_extopdep_vlcno = 0.0
      allocate (Aerosol_diags%lw_absopdep_vlcno  &
                                        (size(Aerosol%aerosol,1), &
                                         size(Aerosol%aerosol,2), &
                                         size(Aerosol%aerosol,3)+1,2))
      Aerosol_diags%lw_absopdep_vlcno = 0.0


!--------------------------------------------------------------------
!    if the volcanic sw aerosol extinction is being supplied, obtain
!    the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_ext) then

!---------------------------------------------------------------------
!    if new sw extinction data is needed on this step, call interpolator
!    to obtain it.  if the data is not to be interpolated, save the
!    retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_ext) then
            if (nfields_sw_ext >= 1) then
              call interpolator (Sw_aer_extopdep_interp, Volcano_Time, &
                                 p_half, Aerosol_props%sw_ext,    &
                                 sw_ext_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              sw_ext_save(is:ie,js:je,:,:) = Aerosol_props%sw_ext
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosol_props%sw_ext = sw_ext_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol single scattering albedo is being 
!    supplied, obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_ssa) then

!---------------------------------------------------------------------
!    if new sw single scattering albedo data is needed on this step, 
!    call interpolator to obtain it.  if the data is not to be inter-
!    polated, save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_ssa) then
            if (nfields_sw_ssa >= 1) then
              call interpolator (Sw_aer_ssalb_interp, Volcano_Time, &
                                 p_half, Aerosol_props%sw_ssa,    &
                                 sw_ssa_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              sw_ssa_save(is:ie,js:je,:,:) = Aerosol_props%sw_ssa
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosol_props%sw_ssa = sw_ssa_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic sw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_sw_asy) then

!---------------------------------------------------------------------
!    if new sw asymmetry factor data is needed on this step, call 
!    interpolator to obtain it.  if the data is not to be interpolated,
!    save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_sw_asy) then
            if (nfields_sw_asy >= 1) then
              call interpolator (Sw_aer_asymm_interp, Volcano_Time, &
                                 p_half, Aerosol_props%sw_asy,    &
                                 sw_asy_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              sw_asy_save(is:ie,js:je,:,:) = Aerosol_props%sw_asy
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if (.not. interpolating_volcanic_data) then
              Aerosol_props%sw_asy = sw_asy_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol extinction is being supplied, obtain
!    the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_ext) then

!---------------------------------------------------------------------
!    if new lw extinction data is needed on this step, call interpolator
!    to obtain it.  if the data is not to be interpolated, save the
!    retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_ext) then
            if (nfields_lw_ext >= 1) then
              call interpolator (Lw_aer_extopdep_interp, Volcano_Time, &
                                 p_half, Aerosol_props%lw_ext,    &
                                 lw_ext_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              lw_ext_save(is:ie,js:je,:,:) = Aerosol_props%lw_ext
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosol_props%lw_ext = lw_ext_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw single scattering albedo is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_ssa) then

!---------------------------------------------------------------------
!    if new lw single scattering albedo data is needed on this step, 
!    call interpolator to obtain it.  if the data is not to be inter-
!    polated, save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_ssa) then
            if (nfields_lw_ssa >= 1) then
              call interpolator (Lw_aer_ssalb_interp, Volcano_Time, &
                                 p_half, Aerosol_props%lw_ssa,    &
                                 lw_ssa_name(1), is, js)
            endif
            if ( .not. interpolating_volcanic_data) then
              lw_ssa_save(is:ie,js:je,:,:) = Aerosol_props%lw_ssa
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if ( .not. interpolating_volcanic_data) then
              Aerosol_props%lw_ssa = lw_ssa_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    if the volcanic lw aerosol asymmetry factor is being supplied, 
!    obtain the appropriate data.
!--------------------------------------------------------------------
        if (using_lw_asy) then

!---------------------------------------------------------------------
!    if new lw asymmetry factor data is needed on this step, call 
!    interpolator to obtain it.  if the data is not to be interpolated, 
!    save the retrieved values in a module variable.
!---------------------------------------------------------------------
          if (need_lw_asy) then
            if (nfields_lw_asy >= 1) then
              call interpolator (Lw_aer_asymm_interp, Volcano_Time, &
                                 p_half, Aerosol_props%lw_asy,    &
                                 lw_asy_name(1), is, js)
            endif
            if (.not. interpolating_volcanic_data) then
              lw_asy_save(is:ie,js:je,:,:) = Aerosol_props%lw_asy
            endif

!---------------------------------------------------------------------
!    if new data from the file is not needed on this step, then retrieve
!    the relevant data from the storage variable.
!---------------------------------------------------------------------
          else
            if (.not. interpolating_volcanic_data) then
              Aerosol_props%lw_asy = lw_asy_save(is:ie,js:je,:,:)
            endif
          endif
        endif

!---------------------------------------------------------------------
!    code for treating sulfate and black carbon as an internal aerosol
!    mixture.
!---------------------------------------------------------------------
        if (Rad_control%using_im_bcsul) then
          if (num_sul > 0) then
            sul(:,:,:) = Aerosol%aerosol(:,:,:,sul_ind(1))
            do iaer=2,num_sul
              sul(:,:,:) = sul(:,:,:) +    &
                                Aerosol%aerosol(:,:,:,sul_ind(iaer))
            end do
          else
            sul = 0.
          endif
          if (num_bc > 0) then
            bc(:,:,:) = Aerosol%aerosol(:,:,:,bc_ind(1))
            do iaer=2,num_bc
              bc(:,:,:) = bc(:,:,:) +    &
                                Aerosol%aerosol(:,:,:,bc_ind(iaer))
            end do
          else
            bc = 0.
          endif
          do k = 1,size(Aerosol%aerosol,3)
            do j = 1,size(Aerosol%aerosol,2)
              do i = 1,size(Aerosol%aerosol,1)
                if (bc(i,j,k) > 0 .and. sul(i,j,k) > 0.0) then
                  Aerosol_props%ivol(i,j,k) = 100-MIN(100, MAX( 0,     &
                   NINT(100.*sul(i,j,k)/(sul(i,j,k) +bc(i,j,k)*1.74))))
                else
                  Aerosol_props%ivol(i,j,k) = 0
                end if
              enddo
            end do
          end do
        else
          Aerosol_props%ivol = 0
        endif ! (using_im_bcsul)

!---------------------------------------------------------------------
!    fill the remaining components of the aerosol_properties_type var-
!    iable and return to the calling routine. this variable contains 
!    the aerosol radiative properties for each aerosol properties type 
!    over each solar and aerosol emissivity band.
!---------------------------------------------------------------------
        if (Rad_control%do_swaerosol_forcing .or.  &
                                         Sw_control%do_swaerosol) then
          Aerosol_props%aerextband = aerextband_MOD
          Aerosol_props%aerssalbband = aerssalbband_MOD
          Aerosol_props%aerasymmband = aerasymmband_MOD
        endif

!---------------------------------------------------------------------
!    if longwave aerosol effects are desired, and the following cal-
!    culation has not already been done, calculate the aerosol 
!    properties for each aerosol properties type nw over each aerosol 
!    emissivity band na using the weighted contributions from each
!    aerosol parameterization band ni. mark the calculation as com-
!    pleted.
!
!    the units of extinction coefficient (aeroextivl) are m**2/gm.
!    to make the lw band extinction coefficient (aerextbandlw) have
!    units (m**2/Kg) consistent with the units in FMS models, one
!    must multiply by 1000. this is done below.
!---------------------------------------------------------------------
        if (Rad_control%do_lwaerosol_forcing .or. &
                Lw_control%do_lwaerosol) then
          Aerosol_props%aerextbandlw = aerextbandlw_MOD
          Aerosol_props%aerssalbbandlw = aerssalbbandlw_MOD
          Aerosol_props%aerextbandlw_cn = aerextbandlw_cn_MOD
          Aerosol_props%aerssalbbandlw_cn = aerssalbbandlw_cn_MOD
        endif
        Aerosol_props%sulfate_index = sulfate_index
        Aerosol_props%optical_index = optical_index
        Aerosol_props%omphilic_index =omphilic_index
        Aerosol_props%bcphilic_index =bcphilic_index
        Aerosol_props%seasalt1_index =seasalt1_index
        Aerosol_props%seasalt2_index =seasalt2_index
        Aerosol_props%seasalt3_index =seasalt3_index
        Aerosol_props%seasalt4_index =seasalt4_index
        Aerosol_props%seasalt5_index =seasalt5_index
        Aerosol_props%sulfate_flag =   SULFATE_FLAG
        Aerosol_props%omphilic_flag =  OMPHILIC_FLAG
        Aerosol_props%bcphilic_flag =  BCPHILIC_FLAG
        Aerosol_props%seasalt1_flag =  SEASALT1_FLAG
        Aerosol_props%seasalt2_flag =  SEASALT2_FLAG
        Aerosol_props%seasalt3_flag =  SEASALT3_FLAG
        Aerosol_props%seasalt4_flag =  SEASALT4_FLAG
        Aerosol_props%seasalt5_flag =  SEASALT5_FLAG
        if (Rad_control%using_im_bcsul) then
          Aerosol_props%bc_flag = BC_FLAG
        else
          Aerosol_props%bc_flag = NOT_IN_USE
        endif

end subroutine aerosol_radiative_properties


!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_end">
!  <OVERVIEW>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosolrad_package_end

!--------------------------------------------------------------------
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate module variables.
!---------------------------------------------------------------------
      if (do_swaerosol .or. Rad_control%do_swaerosol_forcing) then      
        deallocate (solivlaero, nivl1aero, nivl2aero, endaerwvnsf, &
                    aeroextivl, aerossalbivl, aeroasymmivl)
      endif
      if (do_lwaerosol .or. Rad_control%do_lwaerosol_forcing) then
        deallocate ( sflwwts)
        deallocate ( sflwwts_cn)
      endif
      
!---------------------------------------------------------------------
!    deallocate elements of the aerosol_properties_type array.
!---------------------------------------------------------------------
      if (allocated(sulfate_index )) deallocate(sulfate_index )
      if (allocated(bcphilic_index)) deallocate(bcphilic_index)
      if (allocated(omphilic_index)) deallocate(omphilic_index)
      if (allocated(seasalt1_index)) deallocate(seasalt1_index)
      if (allocated(seasalt2_index)) deallocate(seasalt2_index)
      if (allocated(seasalt3_index)) deallocate(seasalt3_index)
      if (allocated(seasalt4_index)) deallocate(seasalt4_index)
      if (allocated(seasalt5_index)) deallocate(seasalt5_index)
      if (allocated(optical_index )) deallocate(optical_index )

      
      if (Rad_control%volcanic_lw_aerosols) then
        if (nfields_lw_ext /= 0) then
          call interpolator_end (Lw_aer_extopdep_interp)
        endif
        if (nfields_lw_ssa /= 0) then
        call interpolator_end (Lw_aer_ssalb_interp)
        endif
        if (nfields_lw_asy /= 0) then
        call interpolator_end (Lw_aer_asymm_interp)
        endif
      endif

      if (Rad_control%volcanic_sw_aerosols) then
        if (nfields_sw_ext /= 0) then
        call interpolator_end (Sw_aer_extopdep_interp)
        endif
        if (nfields_sw_ssa /= 0) then
        call interpolator_end (Sw_aer_ssalb_interp)
        endif
        if (nfields_sw_asy /= 0) then
        call interpolator_end (Sw_aer_asymm_interp)
        endif
      endif

      if (.not. interpolating_volcanic_data) then
        deallocate (sw_ext_save, sw_ssa_save, sw_asy_save, &
                    lw_ext_save, lw_ssa_save, lw_asy_save)
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.




end subroutine aerosolrad_package_end


!####################################################################
! <SUBROUTINE NAME="get_aerosol_optical_info">
!  <OVERVIEW>
!    get_aerosol_optical_info accesses data stored by this module.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_aerosol_optical_info accesses data stored by this module.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_aerosol_optical_info( num_categories, nwavenumbers, &
!                                     names, wavenumbers, &
!                                     aer_ext, aer_ss_alb, aer_asymm)
!  </TEMPLATE>
!  <OUT NAME="num_categories" TYPE="integer">
!   number of aerosol properties types
!  </OUT>
!  <OUT NAME="nwavenumbers" TYPE="integer">
!   number of wavenumber bands over which
!                           aerosol properties are defined
!  </OUT>
!  <OUT NAME="names" TYPE="character">
!   names assigned to the optical properties types
!  </OUT>
!  <OUT NAME="wavenumbers" TYPE="real">
!   wavenumber limits for each of the bands for
!                           which aerosol properties are defined
!  </OUT>
!  <OUT NAME="aer_ext, aer_ss_alb, aer_asymm" TYPE="real">
!   Aerosol extinction coefficient, single scattering albedo, and
!   asymmetry parameter
!  </OUT>
! </SUBROUTINE>
!
subroutine get_aerosol_optical_info( num_categories, nwavenumbers, &
                                     names, wavenumbers, &
                                     aer_ext, aer_ss_alb, aer_asymm)

!-----------------------------------------------------------------------
!    get_aerosol_optical_info accesses data stored by this module.
!-----------------------------------------------------------------------

integer,                        intent(out), optional ::       &
                                            num_categories, nwavenumbers
character(len=*), dimension(:), intent(out), optional :: names
integer, dimension(:),          intent(out), optional :: wavenumbers
real, dimension(:,:),           intent(out), optional :: aer_ext, &
                                                         aer_ss_alb, &
                                                         aer_asymm

!----------------------------------------------------------------------
!   intent(out), optional variables:
!
!      num_categories       number of aerosol properties types
!      nwavenumbers         number of wavenumber bands over which
!                           aerosol properties are defined
!      names                names assigned to the optical properties 
!                           types
!      wavenumbers          wavenumber limits for each of the bands for
!                           which aerosol properties are defined
!      aer_ext              extinction coefficient for each aerosol
!                           spectral band and each aerosol optical 
!                           property type
!      aer_ss_ab            single-scattering albedo for each aerosol 
!                           band and each aerosol optical property type 
!      aer_asymm            asymmetry factor for each aerosol band and 
!                           each aerosol optical property type
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the desired output variables.
!---------------------------------------------------------------------
      if( present(num_categories) ) num_categories = naermodels
      if( present(nwavenumbers))    nwavenumbers   = num_wavenumbers
      if( present(names) )          names(:naermodels) =    &
                                     aerosol_optical_names(:naermodels)
      if( present(wavenumbers) )    wavenumbers(:num_wavenumbers) =  &
                                           endaerwvnsf(:num_wavenumbers)
      if( present(aer_ext) )        aer_ext(:,:)    = aeroextivl(:,:)
      if( present(aer_ss_alb) )     aer_ss_alb(:,:) = aerossalbivl(:,:)
      if( present(aer_asymm) )      aer_asymm(:,:)  = aeroasymmivl(:,:)

!---------------------------------------------------------------------

end subroutine get_aerosol_optical_info


!######################################################################
! <SUBROUTINE NAME="get_aerosol_optical_index">
!  <OVERVIEW>
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!  </DESCRIPTION>
!  <TEMPLATE>
!   index = get_aerosol_optical_index( name, naerosol, rh )
!  </TEMPLATE>
!  <IN NAME="name" TYPE="real">
!   aerosol species name for which the optical 
!                      properties index is desired
!  </IN>
!  <IN NAME="naerosol" TYPE="integer">
!   aerosol index of the aerosol for whoch the 
!                      optical properties index is desired
!  </IN>
!  <IN NAME="rh" TYPE="real">
!    relative humidity 
!  </IN>
! </SUBROUTINE>
!
function get_aerosol_optical_index( name, naerosol, rh ) result(index)

!-----------------------------------------------------------------------
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!-----------------------------------------------------------------------

character(len=*),         intent(in) :: name
integer,                   intent(in) :: naerosol
real,                      intent(in) :: rh
integer                               :: index  ! function value

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      name            aerosol species name for which the optical 
!                      properties index is desired
!      naerosol        aerosol index of the aerosol for whoch the 
!                      optical properties index is desired
!      rh              relative humidity 
!
!  function value
!
!      index           returned optical properties index for aerosol
!                      name (aerosol index = naerosol) when the 
!                      relative humidity is rh
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: irh     ! integer value for relative humidity, 
                           ! used as an index
      integer   :: nfields ! total number of active aerosols

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'aerosolrad_package_mod',  &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    be sure the desired aerosol index is valid.
!---------------------------------------------------------------------
      nfields = size(optical_index(:),1)
      if (naerosol > nfields) then
        call error_mesg( 'aerosolrad_package_mod', &
           'aerosol index exceeds number of aerosol fields', FATAL )
      end if

!---------------------------------------------------------------------
!    Determine if the desired aerosol is a hydrophilic or not. 
!    If the aerosol is hydrophilic, then the optical propeties 
!    index will depend on the relative humidity.
!---------------------------------------------------------------------
      if (optical_index(naerosol) == SULFATE_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
!yim no vol info is passed here now. Set 0 for now.
        index = sulfate_index( irh, 0 )
      elseif (optical_index(naerosol) == BC_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
!yim no vol info is passed here now. Set 0 for now.
        index = sulfate_index( irh, 0 )
      elseif (optical_index(naerosol) == &
                                               OMPHILIC_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = omphilic_index( irh )
      elseif (optical_index(naerosol) ==   &
                  BCPHILIC_FLAG .and. Rad_control%using_im_bcsul ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
!yim
        index = sulfate_index( irh, 0 )
      elseif (optical_index(naerosol) == BCPHILIC_FLAG  &
                        .and. .not. Rad_control%using_im_bcsul ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
!yim
        index = bcphilic_index( irh )
      elseif (optical_index(naerosol) ==    &
                                               SEASALT1_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = seasalt1_index( irh )
      elseif (optical_index(naerosol) ==  &            
                                               SEASALT2_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = seasalt2_index( irh )
      elseif (optical_index(naerosol) ==  &            
                                               SEASALT3_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = seasalt3_index( irh )
      elseif (optical_index(naerosol) ==  &           
                                               SEASALT4_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = seasalt4_index( irh )
      elseif (optical_index(naerosol) ==  &           
                                               SEASALT5_FLAG ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = seasalt5_index( irh )
      else
        index = optical_index(naerosol)
      endif

!---------------------------------------------------------------------
!    if no value was obtained for the optical index, stop execution.
!---------------------------------------------------------------------
      if (index == 0 ) then
        call error_mesg ('aerosolrad_package_mod', &
           'Cannot find aerosol optical properties for species = ' // &
                  trim (name), FATAL )
      endif

!----------------------------------------------------------------------


end function get_aerosol_optical_index




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
                                  

!#####################################################################
! <SUBROUTINE NAME="assign_aerosol_opt_props">
!  <OVERVIEW>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </OVERVIEW>
!  <DESCRIPTION>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call assign_aerosol_opt_props (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names associated with each aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine assign_aerosol_opt_props (aerosol_names)

!----------------------------------------------------------------------
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!---------------------------------------------------------------------

!character(len=64), dimension(:), intent(in) :: aerosol_names
character(len=*), dimension(:), intent(in) :: aerosol_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     aerosol_names     names associated with each aerosol species
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      character(len=64) :: name_in, target_name
!yim
      character(len=4)  :: chind, chind2
      integer           :: nfields
!yim
      integer           :: n, noptical, m
      integer           :: ibc, isul

!---------------------------------------------------------------------
!   local variables:
!
!       name_in          variable to hold current aerosol name 
!                        being processed
!       target_name      aerosol_optical_name associated with a given
!                        aerosol species     
!       nfields          number of activated aerosol species
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    count the number of aerosol optical property categories requested
!    via the namelist input.
!----------------------------------------------------------------------
      do n=1,MAX_OPTICAL_FIELDS
        if (aerosol_optical_names(n) /= ' '  ) then
          naermodels = n
        else
          exit
        endif
      end do

!---------------------------------------------------------------------
!    define the number of activated aerosol species.
!---------------------------------------------------------------------
      nfields = size (aerosol_names(:))

      if (Rad_control%using_im_bcsul) then
        allocate (sul_ind(nfields))
        allocate (bc_ind(nfields))
      endif

!---------------------------------------------------------------------
!    allocate components of the aerosol_properties_type module variable
!    which will contain the indices for the different aerosols.
!---------------------------------------------------------------------
      if (Rad_control%using_im_bcsul) then
        allocate (sulfate_index (0:100,0:100))
      else
        allocate (sulfate_index (0:100,0:0  ))
      endif
      allocate (omphilic_index(0:100 ), &
                bcphilic_index(0:100 ), &
                seasalt1_index(0:100 ), &
                seasalt2_index(0:100 ), &
                seasalt3_index(0:100 ), &
                seasalt4_index(0:100 ), &
                seasalt5_index(0:100 ) )
      allocate (optical_index(nfields) )
      optical_index    = 0
      sulfate_index    = 0
      omphilic_index    = 0
      bcphilic_index    = 0
      seasalt1_index    = 0
      seasalt2_index    = 0
      seasalt3_index    = 0
      seasalt4_index    = 0
      seasalt5_index    = 0

!----------------------------------------------------------------------
!    match aerosol optical property indices with aerosol indices.
!    sulfate aerosols are handled separately (below) with RH dependence.
!----------------------------------------------------------------------
      num_sul = 0
      num_bc = 0
      isul = 1
      ibc = 1
      do n=1,nfields
        name_in = trim(aerosol_names(n))
        if (name_in == 'so4' .or. name_in == 'so4_anthro' .or. name_in == 'so4_natural') then
          optical_index(n) = SULFATE_FLAG
          if (Rad_control%using_im_bcsul) then
            num_sul = num_sul +1
            sul_ind(isul) = n
            isul = isul + 1
          endif
        else if (name_in == "omphilic" .or. name_in == "oc_hydrophilic") then
            optical_index(n) = OMPHILIC_FLAG
        else if (name_in == "bcphilic" .or. name_in == "bc_hydrophilic") then
            optical_index(n) = BCPHILIC_FLAG
            if (Rad_control%using_im_bcsul) then
              num_bc = num_bc +1
              bc_ind(ibc) = n
              ibc = ibc + 1
            endif
        else if (name_in == "seasalt1") then
            optical_index(n) = SEASALT1_FLAG
        else if (name_in == "seasalt2") then
            optical_index(n) = SEASALT2_FLAG
        else if (name_in == "seasalt3") then
            optical_index(n) = SEASALT3_FLAG
        else if (name_in == "seasalt4") then
            optical_index(n) = SEASALT4_FLAG
        else if (name_in == "seasalt5") then
            optical_index(n) = SEASALT5_FLAG
!yim
        else if (name_in == "black_carbon" .and.   &
                              Rad_control%using_im_bcsul) then
            optical_index(n) = BC_FLAG
            num_bc = num_bc +1
            bc_ind(ibc) = n
            ibc = ibc + 1
        else 
          select case( name_in )
            case( "anthro_dust_0.1", "natural_dust_0.1" )
              target_name = "dust_0.1"
            case( "anthro_dust_0.2", "natural_dust_0.2" )
              target_name = "dust_0.2"
            case( "anthro_dust_0.4", "natural_dust_0.4" )
              target_name = "dust_0.4"
            case( "anthro_dust_0.8", "natural_dust_0.8" )
              target_name = "dust_0.8"
            case( "anthro_dust_1.0", "natural_dust_1.0" )
              target_name = "dust_1.0"
            case( "anthro_dust_2.0", "natural_dust_2.0" )
              target_name = "dust_2.0"
            case( "anthro_dust_4.0", "natural_dust_4.0" )
              target_name = "dust_4.0"
            case( "anthro_dust_8.0", "natural_dust_8.0" )
              target_name = "dust_8.0"
            case( "black_carbon" )
               target_name = "soot"
            case( "organic_carbon" )
              target_name = "organic_carbon"
            case( "sea_salt" )
              target_name = "sea_salt"
            case( "dust1" )
              target_name = "dust1"
            case( "dust2" )
              target_name = "dust2"
            case( "dust3" )
              target_name = "dust3"
            case( "dust4" )
              target_name = "dust4"
            case( "dust5" )
              target_name = "dust5"
            case( "bcdry" )
              target_name = "bcdry"
            case( "omphobic", "oc_hydrophobic" )
              target_name = "omphobic"
            case( "bcphobic", "bc_hydrophobic" )
              target_name = "bcphobic"
            case DEFAULT
              target_name = name_in
          end select  

!--------------------------------------------------------------------
!    go through the set of aerosol properties types looking for 
!    the target_name defined above. when found, associate the
!    optical properties type index with the current aerosol species.
!--------------------------------------------------------------------
          do noptical=1,naermodels
            if (aerosol_optical_names(noptical) == target_name) then
              optical_index(n) = noptical
              exit
            end if
          end do

!--------------------------------------------------------------------
!    if the target_name is not found, exit with an error message.
!----------------------------------------------------------------------
          if (optical_index(n) == 0 ) then
            call error_mesg( 'aerosolrad_package_mod', &
                'Cannot find aerosol optical model = ' //    &
                                           TRIM( target_name ), FATAL )
          endif
        endif  ! (name_in ==)
      end do  ! (n=1,nfields)

      select case(trim(aerosol_data_set))
        case ('Ginoux_Reddy') 
     if (Rad_control%using_im_bcsul) then

!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
          do m=0,100
            if (sulfate_indices(n) < 10) then
              write (chind, '(i1)') sulfate_indices(n)
            else if (sulfate_indices(n) == 100) then
              write (chind, '(i3)') sulfate_indices(n)
            else
              write (chind, '(i2)') sulfate_indices(n)
            endif
!yim
            if (sulfate_vol_indices(m) < 10) then
              write (chind2, '(i1)') sulfate_vol_indices(m)
            else if (sulfate_vol_indices(m) == 100) then
              write (chind2, '(i3)') sulfate_vol_indices(m)
            else
              write (chind2, '(i2)') sulfate_vol_indices(m)
            endif
!yim format sulfate_10%_10% (RH + volume fraction)
            target_name = 'sulfate_' // trim(chind)  // '%_' // trim(chind2)// '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            sulfate_index(n,m) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (sulfate_index(n,m) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do
      end do
     else
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (sulfate_indices(n) < 10) then
              write (chind, '(i1)') sulfate_indices(n)
            else if (sulfate_indices(n) == 100) then
              write (chind, '(i3)') sulfate_indices(n)
            else
              write (chind, '(i2)') sulfate_indices(n)
            endif
            target_name = 'sulfate_' // trim(chind)  // '%' 

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            sulfate_index(n,0) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (sulfate_index(n,0) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do

     endif

!---------------------------------------------------------------------
!    set up RH-dependent omphilic aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (omphilic_indices(n) < 10) then
              write (chind, '(i1)') omphilic_indices(n)
            else if (omphilic_indices(n) == 100) then
              write (chind, '(i3)') omphilic_indices(n)
            else
              write (chind, '(i2)') omphilic_indices(n)
            endif
            target_name = 'omphilic_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                omphilic_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (omphilic_index(n) == 0 ) then
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do

!---------------------------------------------------------------------
!    set up RH-dependent seasalt1 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (seasalt1_indices(n) < 10) then
              write (chind, '(i1)') seasalt1_indices(n)
            else if (seasalt1_indices(n) == 100) then
              write (chind, '(i3)') seasalt1_indices(n)
            else
              write (chind, '(i2)') seasalt1_indices(n)
            endif
            target_name = 'seasalt1_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                seasalt1_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (seasalt1_index(n) == 0 ) then
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do

!---------------------------------------------------------------------
!    set up RH-dependent seasalt2 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (seasalt2_indices(n) < 10) then
              write (chind, '(i1)') seasalt2_indices(n)
            else if (seasalt2_indices(n) == 100) then
              write (chind, '(i3)') seasalt2_indices(n)
            else
              write (chind, '(i2)') seasalt2_indices(n)
            endif
            target_name = 'seasalt2_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                seasalt2_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (seasalt2_index(n) == 0 ) then
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do

!---------------------------------------------------------------------
!    set up RH-dependent seasalt3 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (seasalt3_indices(n) < 10) then
              write (chind, '(i1)') seasalt3_indices(n)
            else if (seasalt3_indices(n) == 100) then
              write (chind, '(i3)') seasalt3_indices(n)
            else
              write (chind, '(i2)') seasalt3_indices(n)
            endif
            target_name = 'seasalt3_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                seasalt3_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (seasalt3_index(n) == 0 ) then
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do
 
!---------------------------------------------------------------------
!    set up RH-dependent seasalt4 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (seasalt4_indices(n) < 10) then
              write (chind, '(i1)') seasalt4_indices(n)
            else if (seasalt4_indices(n) == 100) then
              write (chind, '(i3)') seasalt4_indices(n)
            else
              write (chind, '(i2)') seasalt4_indices(n)
            endif
            target_name = 'seasalt4_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                seasalt4_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (seasalt4_index(n) == 0 ) then
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do

!-------------------------------------------------------------------
!    set up RH-dependent seasalt5 aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (seasalt5_indices(n) < 10) then
              write (chind, '(i1)') seasalt5_indices(n)
            else if (seasalt5_indices(n) == 100) then
              write (chind, '(i3)') seasalt5_indices(n)
            else
              write (chind, '(i2)') seasalt5_indices(n)
            endif
            target_name = 'seasalt5_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
            do noptical=1,naermodels
              if (aerosol_optical_names(noptical) == target_name ) then
                seasalt5_index(n) = noptical
                exit
              endif
            end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
            if (seasalt5_index(n) == 0 ) then 
              call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
            endif
          end do

          if ( .not. Rad_control%using_im_bcsul) then
!-------------------------------------------------------------------
!    set up RH-dependent bcphilic aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
           do n=0,100
             if (bcphilic_indices(n) < 10) then
               write (chind, '(i1)') bcphilic_indices(n)
             else if (bcphilic_indices(n) == 100) then
               write (chind, '(i3)') bcphilic_indices(n)
             else
               write (chind, '(i2)') bcphilic_indices(n)
             endif
             target_name = 'bcphilic_' // trim(chind)  // '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
             do noptical=1,naermodels
               if (aerosol_optical_names(noptical) == target_name ) then
                 bcphilic_index(n) = noptical
                 exit
               endif
             end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
             if (bcphilic_index(n) == 0 ) then
               call error_mesg( 'aerosolrad_package_mod', &
                  'Cannot find aerosol optical model = ' // &
                                           TRIM( target_name), FATAL )
             endif
           end do
          endif
        case ('shettle_fenn')

   if (Rad_control%using_im_bcsul) then
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
          do m=0,100
            if (sulfate_indices(n) < 10) then
              write (chind, '(i1)') sulfate_indices(n)
            else if (sulfate_indices(n) == 100) then
              write (chind, '(i3)') sulfate_indices(n)
            else
              write (chind, '(i2)') sulfate_indices(n)
            endif
!yim
            if (sulfate_vol_indices(m) < 10) then
              write (chind2, '(i1)') sulfate_vol_indices(m)
            else if (sulfate_vol_indices(m) == 100) then
              write (chind2, '(i3)') sulfate_vol_indices(m)
            else
              write (chind2, '(i2)') sulfate_vol_indices(m)
            endif
!yim format sulfate_10%_10% (RH + volume fraction)
            target_name = 'sulfate_' // trim(chind)  // '%_' // trim(chind2)// '%'

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            sulfate_index(n,m) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (sulfate_index(n,m) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do
      end do
   else 
!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!    define the optical properties type for all possible values of 
!    relative humidity.
!-------------------------------------------------------------------
          do n=0,100
            if (sulfate_indices(n) < 10) then
              write (chind, '(i1)') sulfate_indices(n)
            else if (sulfate_indices(n) == 100) then
              write (chind, '(i3)') sulfate_indices(n)
            else
              write (chind, '(i2)') sulfate_indices(n)
            endif
            target_name = 'sulfate_' // trim(chind)  // '%' 

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            sulfate_index(n,0) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (sulfate_index(n,0) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do

   endif
      end select  ! (aerosol_data_set)  

!---------------------------------------------------------------------



end subroutine assign_aerosol_opt_props



!######################################################################
! <SUBROUTINE NAME="read_optical_input_file">
!  <OVERVIEW>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_optical_input_file
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_optical_input_file

!-----------------------------------------------------------------------
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      real,    dimension(:), allocatable    :: aeroext_in,   &
                                               aerossalb_in,   &
                                               aeroasymm_in
      logical, dimension(:), allocatable    :: found

      integer           :: unit, num_input_categories
      character(len=64) :: name_in
      integer           :: n, noptical

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    open the ASCII input file containing aerosol optical property
!    information.
!----------------------------------------------------------------------
      call mpp_open (unit, 'INPUT/'//optical_filename, MPP_RDONLY,  &
                     MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

!----------------------------------------------------------------------
!    read the dimension information contained in the input file.
!----------------------------------------------------------------------
      read ( unit,* ) num_wavenumbers
      read ( unit,* ) num_input_categories

!----------------------------------------------------------------------
!    read wavenumber limits for aerosol parameterization bands from 
!    the input file.
!----------------------------------------------------------------------
       allocate (endaerwvnsf(num_wavenumbers) )
       read (unit,* )
       read (unit,* ) endaerwvnsf
 
!----------------------------------------------------------------------
!    allocate module arrays to hold the specified sw properties for 
!    each parameterization bnad and each aerosol properties type.
!----------------------------------------------------------------------
      allocate (       &
            aeroextivl   (num_wavenumbers, naermodels),&
            aerossalbivl (num_wavenumbers, naermodels), &
            aeroasymmivl (num_wavenumbers, naermodels) )

!----------------------------------------------------------------------
!    allocate local working arrays.
!----------------------------------------------------------------------
      allocate (aeroext_in   (num_wavenumbers ),             &
                aerossalb_in (num_wavenumbers ),           &
                aeroasymm_in (num_wavenumbers ),           &
                found        (naermodels ) )

!----------------------------------------------------------------------
!    match the names of optical property categories from input file with
!    those specified in the namelist, and store the following data
!    appropriately. indicate that the data has been found.
!----------------------------------------------------------------------
      found(:) = .false.
      do n=1,num_input_categories
        read( unit,* ) name_in
        read( unit,* )
        read( unit,* ) aeroext_in
        read( unit,* )
        read( unit,* ) aerossalb_in
        read( unit,* )
        read( unit,* ) aeroasymm_in
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == name_in) then
            aeroextivl(:,noptical)   = aeroext_in
            aerossalbivl(:,noptical) = aerossalb_in
            aeroasymmivl(:,noptical) = aeroasymm_in
            found( noptical ) = .true.
            exit
          endif
        end do
      end do

!----------------------------------------------------------------------
!    close the ASCII input file.
!----------------------------------------------------------------------
      call mpp_close( unit )

!----------------------------------------------------------------------
!    check to make sure data for all aerosol optical property
!    categories specified in namelist were contained in ASCII
!    input file. if not, exit with a message.
!----------------------------------------------------------------------
      do noptical = 1,naermodels
        if (.not. found( noptical ) ) then
              call error_mesg( 'aerosolrad_package_mod', &
              'Cannot find aerosol optical properties for ' // &
                TRIM(aerosol_optical_names(noptical)),  FATAL )
        endif
      end do

!----------------------------------------------------------------------
!    deallocate local working arrays.
!----------------------------------------------------------------------
      deallocate (aeroext_in, aerossalb_in, aeroasymm_in, found)



end subroutine read_optical_input_file



!#####################################################################
! <SUBROUTINE NAME="sw_aerosol_interaction">
!  <OVERVIEW>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine sw_aerosol_interaction

!-----------------------------------------------------------------------
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      integer           :: nbands, nband, nivl3
      real              :: sumsol3
      integer           :: nw
      integer           :: nmodel

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       nbands           number of bands in solar spectral param-
!                        eterization
!       nband            currently active solar spectrum band 
!       nivl3            currently active aerosol parameterization band
!       sumsol3          sum of solar input in current aerosol param-
!                        eterization band
!       n, nw, noptical  do-loop indices
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!    define the number of bands in the solar spectrum parameterization.
!    allocate space for variables defining the highest and lowest 
!    aerosol parameterization wavenumber in each solar spectral band
!    and the solar flux common to solar spectral band n and aerosol
!    parameterization band ni.
!---------------------------------------------------------------------
      nbands = Solar_spect%nbands 
      allocate ( nivl1aero  (nbands) )
      allocate ( nivl2aero  (nbands) )
      allocate ( solivlaero (nbands, num_wavenumbers))

!---------------------------------------------------------------------
!    define the solar weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the solar
!    spectral intervals and so determine the single-scattering proper-
!    ties on the solar spectral intervals.
!--------------------------------------------------------------------
      nivl3 = 1
      sumsol3 = 0.0
      nband = 1
      solivlaero(:,:) = 0.0
      nivl1aero(1) = 1
      do nw = 1,Solar_spect%endwvnbands(nbands)
        sumsol3 = sumsol3 + Solar_spect%solarfluxtoa(nw)
        if (nw == endaerwvnsf(nivl3) ) then
          solivlaero(nband,nivl3) = sumsol3
          sumsol3 = 0.0
        end if
        if ( nw == Solar_spect%endwvnbands(nband) ) then
          if ( nw /= endaerwvnsf(nivl3) ) then
            solivlaero(nband,nivl3) = sumsol3 
            sumsol3 = 0.0
          end if
          nivl2aero(nband) = nivl3
          nband = nband + 1
          if ( nband <= nbands ) then
            if ( nw == endaerwvnsf(nivl3) ) then
              nivl1aero(nband) = nivl3 + 1
            else
              nivl1aero(nband) = nivl3
            end if
          end if
        end if
        if ( nw == endaerwvnsf(nivl3) ) nivl3 = nivl3 + 1
      end do

!---------------------------------------------------------------------
!    allocate and initialize variables which will hold the aerosol 
!    radiative properties for each solar spectral parameterization band.
!    aerextband     the solar band values of the extinction 
!                   coefficient for aerosols                           
!    aerssalbband   the solar band values of the single-     
!                   scattering albedo for aerosols                      
!    aerasymmband   the solar band values of the asymmetry   
!                   factor for aerosols                                 
!---------------------------------------------------------------------
      allocate     &
        (aerextband_MOD   (Solar_spect%nbands, naermodels), &
         aerssalbband_MOD (Solar_spect%nbands, naermodels), &
         aerasymmband_MOD (Solar_spect%nbands, naermodels) )
      aerextband_MOD   = 0.
      aerssalbband_MOD = 0.
      aerasymmband_MOD = 0.

!--------------------------------------------------------------------
!    if sw aerosol properties are desired and have not yet been calc-
!    ulated, use the thick-averaging technique to define the single-
!    scattering properties for each solar parameterization band n 
!    from the specified properties on the aerosol parameterization 
!    bands ni for each aerosol properties type nmodel. 
! references:                                                          
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!    code I: choosing a configuration for a large-scale model.,     
!    q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: a thin-averaging technique (subroutine thinavg in 
!    rad_utilities_mod) is also available.   
!--------------------------------------------------------------------
      do nmodel=1,naermodels
        call thickavg (nivl1aero, nivl2aero, num_wavenumbers,   &
                       Solar_spect%nbands, aeroextivl(:,nmodel), &
                       aerossalbivl(:,nmodel),    &
                       aeroasymmivl(:,nmodel), solivlaero,   &
                       Solar_spect%solflxbandref,  &
                       aerextband_MOD(:,nmodel),    &
                       aerssalbband_MOD(:,nmodel),   &
                       aerasymmband_MOD(:,nmodel))
      end do

!---------------------------------------------------------------------



end subroutine sw_aerosol_interaction   



!#####################################################################
! <SUBROUTINE NAME="lw_aerosol_interaction">
!  <OVERVIEW>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_aerosol_interaction      

!----------------------------------------------------------------------
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    the following arrays define the wavenumber ranges for the separate
!    aerosol emissivity bands in the model infrared parameterization. 
!    these may be changed only by the keeper of the radiation code.
!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code.
!
!      aerbandlo_fr      low wavenumber limit for the non-continuum 
!                        aerosol emissivity bands
!      aerbandhi_fr      high wavenumber limit for the non-continuum
!                        aerosol emissivity bands
!      istartaerband_fr  starting wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      iendaerband_fr    ending wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      aerbandlo_co      low wavenumber limit for the continuum 
!                        aerosol emissivity bands
!      aerbandhi_co      high wavenumber limit for the continuum
!                        aerosol emissivity bands
!      istartaerband_co  starting wavenumber index for the continuum
!                        aerosol emissivity bands
!      iendaerband_co    ending wavenumber index for the continuum
!                        aerosol emissivity bands
!      aerbandlo         low wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      aerbandhi         high wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      istartaerband     starting wavenumber index for the entire set of
!                        aerosol emissivity bands
!      iendaerband       ending wavenumber index for the entire set of
!                        aerosol emissivity bands
!
!----------------------------------------------------------------------
      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandlo_fr =  &
      (/ 560.0, 630.0, 700.0, 800.0, 900.0,  990.0, 1070.0, 1200.0 /)

      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandhi_fr =  &
      (/ 630.0, 700.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandlo_co =  &
      (/ 560.0 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandhi_co =  &
      (/ 800.0 /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: istartaerband_co =  &
      (/ 57  /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: iendaerband_co =  &
      (/ 80  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: istartaerband_cn =  &
      (/ 81  /)

      integer, dimension (N_AEROSOL_BANDS_CN)  :: iendaerband_cn =  &
      (/ 120 /)

      real,    dimension(N_AEROSOL_BANDS)      :: aerbandlo, aerbandhi
      integer, dimension(N_AEROSOL_BANDS)      :: istartaerband,    &
                                                  iendaerband

!---------------------------------------------------------------------
!    the following arrays define how the ir aerosol band structure 
!    relates to the aerosol parameterization bands.
!
!      nivl1aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl2aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl1aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl2aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl1aer(n)       aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber for the 
!                        ir aerosol emissivity band n
!      nivl2aer(n)       aerosol parameterization band index corres-
!                        ponding to the highest wavenumber for the 
!                        ir aerosol emissivity band n
!      planckaerband(n)  planck function summed over each lw param-
!                        eterization band that is contained in the 
!                        ir aerosol emissivity band n
!
!---------------------------------------------------------------------
      integer, dimension (N_AEROSOL_BANDS_FR)  :: nivl1aer_fr,   &
                                                  nivl2aer_fr
      integer, dimension (N_AEROSOL_BANDS_CO)  :: nivl1aer_co,   &
                                                  nivl2aer_co
      integer, dimension (N_AEROSOL_BANDS_CN)  :: nivl1aer_cn,   &
                                                  nivl2aer_cn
      real,    dimension (N_AEROSOL_BANDS)     :: planckaerband
      real,    dimension (N_AEROSOL_BANDS_CN)  :: planckaerband_cn

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
!        aerextbandlw_fr(n,na)  band averaged extinction coefficient
!                               for non-continuum aerosol emissivity 
!                               band n and aerosol properties type na
!        aerssalbbandlw_fr(n,na)
!                               band averaged single-scattering
!                               coefficient for non-continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        aerextbandlw_co(n,na)  band averaged extinction coefficient
!                               for the continuum aerosol emissivity
!                               band n and aerosol properties type na
!        aerssalbbandlw_co(n,na)
!                               band averaged single-scattering
!                               coefficient for continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        planckivlaer_fr(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity non-
!                               continuum band n and aerosol parameter-
!                               ization band ni
!        planckivlaer_co(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity continuum 
!                               band n and aerosol parameterization 
!                               band ni
!        sflwwts_fr(n,ni)       band weights for the aerosol emissivity
!                               non-continuum band n and the aerosol 
!                               parameterization band ni 
!        sflwwts_co(n,ni)       band weights for the aerosol emissivity
!                               continuum band n and the aerosol 
!                               parameterization band ni 
!        planckivlaer(n,ni)     planck function over the spectral range
!                               common to aerosol emissivity band n and
!                               aerosol parameterization band ni
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real,    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real,    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real,    dimension (N_AEROSOL_BANDS_CN, num_wavenumbers) :: &
                                                  planckivlaer_cn   
      integer, dimension (num_wavenumbers)    ::  iendsfbands

!---------------------------------------------------------------------
!    variables associated with the planck function calculation.
!    the planck function is defined for each of the NBLW longwave 
!    parameterization bands.
!---------------------------------------------------------------------
      real, dimension(NBLW)  :: c1, centnb, sc, src1nb, x, x1
      real                   :: del, xtemv, sumplanck

!---------------------------------------------------------------------
!    miscellaneous variables:

     logical         :: do_band1   !  should we do special calculation 
                                   !  for band 1 ?
     integer         :: ib, nw, nivl, nband, n, ni 
                                   !  do-loop indices and counters
     integer         :: na

!--------------------------------------------------------------------
!    define arrays containing the characteristics of all the ir aerosol
!    emissivity bands, both continuum and non-continuum.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        aerbandlo(n)     = aerbandlo_fr(n)
        aerbandhi(n)     = aerbandhi_fr(n)
        istartaerband(n) = istartaerband_fr(n)
        iendaerband(n)   = iendaerband_fr(n)
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        aerbandlo(n)     = aerbandlo_co     (n - N_AEROSOL_BANDS_FR)
        aerbandhi(n)     = aerbandhi_co     (n - N_AEROSOL_BANDS_FR)
        istartaerband(n) = istartaerband_co (n - N_AEROSOL_BANDS_FR)
        iendaerband(n)   = iendaerband_co   (n - N_AEROSOL_BANDS_FR)
      end do

!---------------------------------------------------------------------
!    define the number of aerosol ir bands to be used in other modules.
!    set the initialization flag to .true.
!---------------------------------------------------------------------
      Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS
      Lw_parameters%n_lwaerosol_bands_iz = .true.

!--------------------------------------------------------------------
!    allocate a module variable which will store the weighting function
!    between the aerosol emissivity bands and the aerosol parameter-
!    ization bands.
!--------------------------------------------------------------------
      allocate (sflwwts (N_AEROSOL_BANDS, num_wavenumbers))
      allocate (sflwwts_cn (N_AEROSOL_BANDS_CN, num_wavenumbers))

!--------------------------------------------------------------------
!    define the ending aerosol band index for each of the aerosol
!    parameterization bands.
!--------------------------------------------------------------------
      iendsfbands(:) = INT((endaerwvnsf(:) + 0.01)/10.0)

!--------------------------------------------------------------------
!    compute the planck function at 10C over each of the longwave
!    parameterization bands to be used as the weighting function. 
!--------------------------------------------------------------------
      do n=1,NBLW 
        del  = 10.0E+00
        xtemv = 283.15
        centnb(n) = 5.0 + (n - 1)*del
        c1(n)     = (3.7412E-05)*centnb(n)**3
        x(n)      = 1.4387E+00*centnb(n)/xtemv
        x1(n)     = EXP(x(n))
        sc(n)     = c1(n)/(x1(n) - 1.0E+00)
        src1nb(n) = del*sc(n)
      end do
 
!--------------------------------------------------------------------
!    sum the weighting function calculated over the longwave param-
!    eterization bands that are contained in each of the aerosol 
!    emissivity bands. 
!--------------------------------------------------------------------
      planckaerband(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS
        do ib = istartaerband(n),iendaerband(n)
          planckaerband(n) = planckaerband(n) + src1nb(ib)
        end do
      end do
      planckaerband_cn(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS_CN
        do ib = istartaerband_cn(n),iendaerband_cn(n)
          planckaerband_cn(n) = planckaerband_cn(n) + src1nb(ib)
        end do
      end do
 
!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the non-
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_fr(:,:) = 0.0
      nivl1aer_fr(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_fr(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_fr(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_fr(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_FR ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_fr(nband) = nivl + 1
            else
              nivl1aer_fr(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband .eq. 1 .and.   &
              iendsfbands(nivl-1) >= istartaerband_fr(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_fr(1)) then
            nivl1aer_fr(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if (nw >= iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_co(:,:) = 0.0
      nivl1aer_co(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_co(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_co(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_co(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CO ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_co(nband) = nivl + 1
            else
              nivl1aer_co(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_co(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_co(1)) then
            nivl1aer_co(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_co(N_AEROSOL_BANDS_CO) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_cn(:,:) = 0.0
      nivl1aer_cn(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_cn(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_cn(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_cn(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_cn(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CN ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_cn(nband) = nivl + 1
            else
              nivl1aer_cn(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_cn(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_cn(1)) then
            nivl1aer_cn(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_cn(N_AEROSOL_BANDS_CN) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the planck-function-weighted band weights for the aerosol
!    parameterization bands onto the non-continuum and continuum ir 
!    aerosol emissivity bands.
!--------------------------------------------------------------------
      sflwwts_fr(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_FR
        do ni=nivl1aer_fr(n),nivl2aer_fr(n)
          sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
        end do
      end do
      sflwwts_co(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CO
        do ni=nivl1aer_co(n),nivl2aer_co(n)
          sflwwts_co(n,ni) = planckivlaer_co(n,ni)/     &
                             planckaerband(N_AEROSOL_BANDS_FR+n)
        end do
      end do
      sflwwts_cn(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CN
        do ni=nivl1aer_cn(n),nivl2aer_cn(n)
          sflwwts_cn(n,ni) = planckivlaer_cn(n,ni)/     &
                             planckaerband_cn(n)
        end do
      end do

!--------------------------------------------------------------------
!    consolidate the continuum and non-continuum weights into an
!    array covering all ir aerosol emissivity bands.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_fr(n,ni)
        end do
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        end do
      end do

!-----------------------------------------------------------------
!    allocate and initialize the arrays in the aerosol_properties_type 
!    variable that will contain the ir aerosol properties for each
!    aerosol optical type over each ir aerosol emissivity band.
!----------------------------------------------------------------
      allocate     &
         (aerextbandlw_MOD   (N_AEROSOL_BANDS, naermodels), &
          aerssalbbandlw_MOD (N_AEROSOL_BANDS, naermodels), &
          aerextbandlw_cn_MOD   (N_AEROSOL_BANDS_CN, naermodels), &
          aerssalbbandlw_cn_MOD (N_AEROSOL_BANDS_CN, naermodels) )
      aerextbandlw_MOD   = 0.0E+00
      aerssalbbandlw_MOD = 0.0E+00
      aerextbandlw_cn_MOD   = 0.0E+00
      aerssalbbandlw_cn_MOD = 0.0E+00

      if (.not. force_to_repro_quebec) then
!---------------------------------------------------------------------
!    if longwave aerosol effects are desired, and the following cal-
!    culation has not already been done, calculate the aerosol 
!    properties for each aerosol properties type nw over each aerosol 
!    emissivity band na using the weighted contributions from each
!    aerosol parameterization band ni. mark the calculation as com-
!    pleted.
!
!    the units of extinction coefficient (aeroextivl) are m**2/gm.
!    to make the lw band extinction coefficient (aerextbandlw) have
!    units (m**2/Kg) consistent with the units in FMS models, one
!    must multiply by 1000. this is done below.
!---------------------------------------------------------------------
        do nw=1,naermodels    
          do na=1,N_AEROSOL_BANDS  
            do ni=1,num_wavenumbers 
              aerextbandlw_MOD(na,nw) = aerextbandlw_MOD(na,nw) + &
                                aeroextivl(ni,nw)*sflwwts(na,ni)*1.0E+03
              aerssalbbandlw_MOD(na,nw) = aerssalbbandlw_MOD(na,nw) + &
                                     aerossalbivl(ni,nw)*sflwwts(na,ni)
            end do
          end do
        end do
        do nw=1,naermodels    
          do na=1,N_AEROSOL_BANDS_CN
            do ni=1,num_wavenumbers 
              aerextbandlw_cn_MOD(na,nw) = aerextbandlw_cn_MOD(na,nw) +&
                            aeroextivl(ni,nw)*sflwwts_cn(na,ni)*1.0E+03
              aerssalbbandlw_cn_MOD(na,nw) =    &
                               aerssalbbandlw_cn_MOD(na,nw) +&
                                  aerossalbivl(ni,nw)*sflwwts_cn(na,ni)
            end do
          end do
        end do
      endif

!----------------------------------------------------------------------



end subroutine lw_aerosol_interaction


!######################################################################



               end module aerosolrad_package_mod






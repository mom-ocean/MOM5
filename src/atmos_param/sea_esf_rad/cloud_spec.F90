                 module cloud_spec_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
! </OVERVIEW>
! <DESCRIPTION>
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
! </DESCRIPTION>

!   shared modules:

use time_manager_mod,         only: time_type, time_manager_init, &
                                    set_time, operator (+)
use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: open_namelist_file, mpp_pe, &
                                    mpp_root_pe, stdlog,  fms_init, &
                                    write_version_number, file_exist, & 
                                    check_nml_error, error_mesg,   &
                                    FATAL, NOTE, close_file
use tracer_manager_mod,       only:         &
!                                   tracer_manager_init,  &
                                    get_tracer_index, NO_TRACER
use field_manager_mod,        only:       &
                                    field_manager_init, &
                                    MODEL_ATMOS
use data_override_mod,        only: data_override
use random_numbers_mod,    only:  randomNumberStream,   &
                                  initializeRandomNumberStream, &
                                  getRandomNumbers,             &
                                  constructSeed
use cloud_generator_mod,   only:  cloud_generator_init, &
                                  cloud_generator_end
use constants_mod,         only : radian, RDGAS

! shared radiation package modules:

use rad_utilities_mod,        only: rad_utilities_init, &
                                    cld_specification_type, &
                                    atmos_input_type, &
                                    surface_type, &
                                    Rad_control, &
                                    microphysics_type,  &         
                                    Cldrad_control
use esfsw_parameters_mod,     only: esfsw_parameters_init, Solar_spect

! interface modules to various cloud parameterizations:

use strat_clouds_W_mod,       only: strat_clouds_W_init,   &
                                    strat_clouds_amt, strat_clouds_W_end
use diag_clouds_W_mod,        only: diag_clouds_W_init,   &
                                    diag_clouds_amt, &
                                    diag_clouds_W_end
use zetac_clouds_W_mod,       only: zetac_clouds_W_init,   &
                                    zetac_clouds_amt, &
                                    zetac_clouds_W_end
use specified_clouds_W_mod,   only: specified_clouds_W_init, &
                                    specified_clouds_amt, &
                                    specified_clouds_W_end
use rh_based_clouds_mod,      only: rh_based_clouds_init,  &
                                    rh_clouds_amt, &
                                    rh_based_clouds_end
use donner_deep_clouds_W_mod, only: donner_deep_clouds_W_init, &
                                    donner_deep_clouds_amt, &
                                    donner_deep_clouds_W_end
use uw_clouds_W_mod,          only: uw_clouds_W_init, &
                                    uw_clouds_amt, &
                                    uw_clouds_W_end
use mgrp_prscr_clds_mod,      only: mgrp_prscr_clds_init, &
                                    prscr_clds_amt,  & 
                                    mgrp_prscr_clds_end 
use standalone_clouds_mod,    only: standalone_clouds_init, &
                                    standalone_clouds_amt, &
                                    standalone_clouds_end
                                 
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present, or from a prescribed formula based on
!    prescribed water paths for high, middle and low clouds.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloud_spec.F90,v 20.0 2013/12/13 23:19:01 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloud_spec_init, cloud_spec,    &
         cloud_spec_dealloc, cloud_spec_end

private    &

!  called from cloud_spec:
         initialize_cldamts, microphys_presc_conc,  &
         combine_cloud_properties


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  ::      &
              cloud_type_form = '     ' ! cloud parameterization being 
                                        ! used; either 'strat', 'rh', 
                                        ! 'deep',  'stratdeep', 'zonal',
                                        ! 'obs', 'prescribed', 'diag', 
                                        ! 'none', 'specified', 'zetac'
                                        ! 'specified_strat', 'stratuw',
                                        ! 'stratdeepuw', 'uw', 'deepuw
                                        ! or 'not_sea_esf'       
real :: wtr_cld_reff=10.                ! assumed cloud drop efective
                                        ! radius [ microns ]  
real :: ice_cld_reff=50.                ! assumed ice cloud effective
                                        ! size [ microns ]
real :: rain_reff=250.                  ! assumed rain drop effective
                                        ! radius [ microns ]
character(len=16) :: overlap_type = 'random'    
                                        ! cloud overlap assumption; 
                                        ! allowable values are 'random'
                                        ! or 'max-random'  
logical :: doing_data_override=.false.
logical :: do_fu2007 = .false.
logical :: do_rain   = .false. !sjl
logical :: do_snow   = .false. !miz
logical :: do_graupel  = .false. !sjl
logical :: force_use_of_temp_for_seed = .false.  
                                        ! if true, when using stochastic 
                                        ! clouds, force the seed to use 
                                        ! top-model-level temps as input to
                                        ! random number generator
                                        ! (needed for some 
                                        ! specialized applications)
logical :: ignore_donner_cells = .false.! when set to .true., the effects 
                                        ! of donner cell clouds in the
                                        ! radiation code are ignored
logical :: do_legacy_seed_generation = .false.
                                        ! setting this variable to .true.
                                        ! (not recommended except to
                                        ! reproduce previous results) will
                                        ! activate the seed generation
                                        ! scheme used previously (through
                                        ! the siena code). this scheme may
                                        ! exhibit flaws at hi-res or when
                                        ! time is held fixed in the 
                                        ! radiation calculation. 

namelist /cloud_spec_nml / cloud_type_form, wtr_cld_reff,   &
                           ice_cld_reff, rain_reff, overlap_type, &
                           doing_data_override, do_fu2007,    &
                           do_rain, do_snow, do_graupel, &
                           do_legacy_seed_generation, &
                           force_use_of_temp_for_seed, ignore_donner_cells

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    assumed water paths.
!--------------------------------------------------------------------
real   ::  lwpath_hi  = 6.313929   ! assumed water path for high clouds
                                   ! [ grams / m**2 ]
real   ::  lwpath_mid = 18.94179   ! assumed water path for middle 
                                   ! clouds [ grams / m**2 ]
real   ::  lwpath_low = 75.76714   ! assumed water path for low clouds
                                   ! [ grams / m**2 ]

!---------------------------------------------------------------------
!    logical  flags.

logical :: module_is_initialized = .false.   ! module initialized ?

!---------------------------------------------------------------------
!    time-step related constants.

integer :: num_pts       !  number of grid columns processed so far that
                         !  have cloud data present (used to identify
                         !  module coldstart condition)
integer :: tot_pts       !  total number of grid columns in the 
                         !  processor's domain

!---------------------------------------------------------------------
!     indices for cloud tracers

integer :: nql           ! tracer index for liquid water
integer :: nqi           ! tracer index for ice water
integer :: nqa           ! tracer index for cloud area
integer :: nqn           ! tracer index for cloud droplet number
integer :: nqni          ! tracer index for ice crystal number
integer :: nqr, nqs, nqg ! tracer index for rainwat, snowwat and graupel           

!----------------------------------------------------------------------
!   variables needed for random number seed:
!----------------------------------------------------------------------
real, dimension(:,:), allocatable  :: lats, lons ! lat and lon of columns
                                               ! in this processor's
                                               ! domain [ degrees ]

!---------------------------------------------------------------------
!     miscellaneous variables:

integer :: num_slingo_bands  ! number of radiative bands over which 
                             ! cloud optical depth is calculated in the
                             ! gordon diag_cloud parameterization
integer :: id, jd, kmax

!----------------------------------------------------------------------
!----------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="cloud_spec_init">
!  <OVERVIEW>
!   Contructor of cloud_spec_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Contructor of cloud_spec_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_init ( pref, lonb, latb, axes, Time)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference pressure levels containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   the longitude array of the model grid box corners
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   the latitude array of the model grid box corners
!  </IN>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
! 
subroutine cloud_spec_init (pref, lonb, latb, axes, Time,  &
                            cloud_type_form_out)

!---------------------------------------------------------------------
!    cloud_spec_init is the constructor for cloud_spec_mod.
!---------------------------------------------------------------------

real, dimension(:,:),     intent(in)   ::  pref        
real, dimension(:,:),     intent(in)   ::  lonb, latb
integer, dimension(4),    intent(in)   ::  axes
type(time_type),          intent(in)   ::  Time
character(len=16),        intent(out)  ::  cloud_type_form_out

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes at cell corners [ radians ]
!       latb      array of model latitudes at cell corners [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
 
      integer   ::   unit, ierr, io, logunit
      integer   ::   ndum, i, j, ii, jj
      

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!      ndum     dummy argument needed for call to field_manager_init
!
!--------------------------------------------------------------------
 
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
      call rad_utilities_init
      call field_manager_init (ndum)
      call esfsw_parameters_init
!  not yet compliant:
!     call tracer_manager_init  ! not public
 
!---------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cloud_spec_nml, iostat=io)
      ierr = check_nml_error(io,"cloud_spec_nml")
#else
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=cloud_spec_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cloud_spec_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=cloud_spec_nml)

      id = size(lonb,1) - 1
      jd = size(latb,2) - 1
      kmax = size(pref,1) - 1

!-----------------------------------------------------------------------
!    define output field.
!-----------------------------------------------------------------------
      cloud_type_form_out = cloud_type_form

!--------------------------------------------------------------------
!    verify a valid type of cloud overlap. set logical variables
!    based on the namelist value.
!--------------------------------------------------------------------
      if (trim(overlap_type) == 'random') then
        Cldrad_control%do_random_overlap = .true.
      else if (trim(overlap_type) == 'max-random') then
        Cldrad_control%do_max_random_overlap = .true.
      else
        call error_mesg ('cloud_spec_mod',  &
         ' invalid specification of overlap_type', FATAL)
      endif

!-------------------------------------------------------------------
!    set the variables indicating that the above control variables have
!    been set.
!--------------------------------------------------------------------
      Cldrad_control%do_random_overlap_iz = .true.
      Cldrad_control%do_max_random_overlap_iz = .true.

!--------------------------------------------------------------------
!    if the sea-esf radiation package is not being used, then 
!    cloud_type_form will have been set to 'not_sea_esf'. in such a
!    case, the clouds will be specified internally within the 
!    radiation_driver_mod, so simply return.
!--------------------------------------------------------------------
      if (trim(cloud_type_form) == 'not_sea_esf')  return
        
!-------------------------------------------------------------------
!    verify that the nml variable cloud_type_form specifies a valid
!    cloud parameterization. set the appropriate logical control
!    variable(s) to .true.. call the constructor modules for the
!    specific cloud scheme(s) requested.
!-------------------------------------------------------------------
      if (trim(cloud_type_form) == 'strat')  then

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the model based on klein 
!    parameterization. strat is an acceptable option both for standalone
!    and gcm applications.
!-------------------------------------------------------------------
        Cldrad_control%do_strat_clouds = .true.
        call strat_clouds_W_init(latb, lonb)

      else if (trim(cloud_type_form) == 'specified_strat')  then
        Cldrad_control%do_specified_strat_clouds = .true.
        Cldrad_control%do_strat_clouds = .true.
        call strat_clouds_W_init(latb, lonb)
        call standalone_clouds_init (pref, lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are diagnosed based on model relative 
!    humidity.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'rh')   then
            Cldrad_control%do_rh_clouds = .true.
            call rh_based_clouds_init 

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the donner deep cloud 
!    (cell cloud, anvil cloud) scheme.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'deep')  then
            Cldrad_control%do_donner_deep_clouds = .true.
            call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)
!RSH 3/6/13: The following call added to allow stochastic clouds to be 
!            run for this case. Ultimately, the do_stochastic_clouds nml 
!            variable should be moved to cloud_spec_nml.
            call strat_clouds_W_init(latb, lonb)

!------------------------------------------------------------------
!    cloud fractions, heights are provided by the uw_conv shallow
!    convection scheme  
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'uw')  then
            Cldrad_control%do_uw_clouds = .true.
            call uw_clouds_W_init (pref, lonb, latb,   &
                                             axes, Time)
!RSH 3/6/13: The following call added to allow stochastic clouds to be 
!            run for this case. Ultimately, the do_stochastic_clouds nml 
!            variable should be moved to cloud_spec_nml.
            call strat_clouds_W_init(latb, lonb)

!-------------------------------------------------------------------
!    cloud fractions, heights are a combination of the donner
!    deep cloud (cell cloud, anvil cloud) and klein large-scale cloud
!    parameterizations.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'stratdeep')  then
            Cldrad_control%do_strat_clouds = .true.
            Cldrad_control%do_donner_deep_clouds = .true.
            call strat_clouds_W_init(latb, lonb)
            call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the donner deep convection
!    (cell cloud, anvil cloud) and uw_conv shallow convection
!    cloud parameterizations.
!-------------------------------------------------------------------
         else if (trim(cloud_type_form) == 'deepuw')  then
           Cldrad_control%do_donner_deep_clouds = .true.
           Cldrad_control%do_uw_clouds = .true.
           call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                             axes, Time)
           call uw_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)
!RSH 3/6/13: The following call added to allow stochastic clouds to be 
!            run for this case. Ultimately, the do_stochastic_clouds nml 
!            variable should be moved to cloud_spec_nml.
            call strat_clouds_W_init(latb, lonb)

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the klein large-scale
!    and uw_conv shallow convection cloud parameterizations.
!-------------------------------------------------------------------
        else if (trim(cloud_type_form) == 'stratuw')  then
          Cldrad_control%do_strat_clouds = .true.
          Cldrad_control%do_uw_clouds = .true.
          call strat_clouds_W_init(latb, lonb)
          call uw_clouds_W_init (pref, lonb, latb,   &
                                             axes, Time)

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the klein large-scale
!    the donner deep convection (cell cloud, anvil cloud) and the
!    uw_conv shallow convection cloud parameterizations.
!-------------------------------------------------------------------
       else if (trim(cloud_type_form) == 'stratdeepuw')  then
         Cldrad_control%do_strat_clouds = .true.
         Cldrad_control%do_donner_deep_clouds = .true.
         Cldrad_control%do_uw_clouds = .true.
         call strat_clouds_W_init(latb, lonb)
         call donner_deep_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)
         call uw_clouds_W_init (pref, lonb, latb,   &
                                            axes, Time)

!---------------------------------------------------------------
!    cloud fractions, heights are prescribed as zonally uniform using
!    the original fms specification.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'zonal')  then
            Cldrad_control%do_zonal_clouds = .true.
            call specified_clouds_W_init (lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are based on observed data set.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'obs')  then
            Cldrad_control%do_obs_clouds = .true.
            call specified_clouds_W_init (lonb, latb)

!-------------------------------------------------------------------
!    cloud fractions, heights are prescribed as zonally invariant, using
!    the formulation from skyhi, with the ability to use a prescribed
!    microphysics. WILL ONLY WORK WITH REGULAR LAT/LON GRID
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'prescribed')  then
            Cldrad_control%do_mgroup_prescribed = .true.
            call mgrp_prscr_clds_init (pref, latb(1,:))

!-------------------------------------------------------------------
!    model is run with gordon diagnostic clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'diag')  then
            Cldrad_control%do_diag_clouds = .true.
            call diag_clouds_W_init (num_slingo_bands)

!-------------------------------------------------------------------
!    model is run with zetac clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'zetac')  then
            Cldrad_control%do_zetac_clouds = .true.
            call zetac_clouds_W_init 
 
!-------------------------------------------------------------------
!    model is run without clouds.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form) == 'none')  then
            Cldrad_control%do_no_clouds = .true.

!-------------------------------------------------------------------
!    model is run with specified clouds and cloud properties.
!-------------------------------------------------------------------
          else if (trim(cloud_type_form)  == 'specified')  then
            Cldrad_control%do_specified_clouds = .true.
            call standalone_clouds_init (pref, lonb, latb)

!-------------------------------------------------------------------
!    failure message if none of the above options was chosen.
!-------------------------------------------------------------------
          else
            call error_mesg ('cloud_spec_mod',  &
              'invalid cloud_type_form specified', FATAL)
      endif  ! (strat)

!--------------------------------------------------------------------
!    define the dimensions of the model subdomain assigned to the 
!    processor.
!--------------------------------------------------------------------
      tot_pts = (size(latb,2)-1)*(size(lonb,1)-1)

!--------------------------------------------------------------------
!    determine if the current run is cold-starting this module. if a 
!    restart file is present, then this is not a coldstart. in that case
!    set num_pts to tot_pts so that if cloud data is not available an 
!    error message can be generated. if this is a coldstart, cloud data
!    will not be available until num_pts equals or exceeds tot_pts, so
!    continue processing without issuing an error message. 
!--------------------------------------------------------------------
      if (file_exist ('INPUT/tracer_cld_amt.res') .or.  &
          file_exist ('INPUT/strat_cloud.res') ) then
        num_pts = tot_pts
      else
        num_pts = 0
      endif

!---------------------------------------------------------------------
!    obtain the tracer indices for the strat_cloud variables when
!    running gcm.
!---------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds .or.  &
            Cldrad_control%do_zetac_clouds) then
          nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
          nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
          nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )

          if (do_rain) then !sjl
             nqr = get_tracer_index ( MODEL_ATMOS, 'rainwat' )
             if (nqr < 0 ) call error_mesg ('cloud_spec_mod', &
                'rainwat tracer not found, but do_rain is true', FATAL)
          end if
          if (do_snow) then !miz
             nqs = get_tracer_index ( MODEL_ATMOS, 'snowwat' )
             if (nqs < 0 ) call error_mesg ('cloud_spec_mod', &
                'snowwat tracer not found, but do_snow is true', FATAL)
          end if
          if (do_graupel) then !sjl
             nqg = get_tracer_index ( MODEL_ATMOS, 'graupel' )
             if (nqg < 0 ) call error_mesg ('cloud_spec_mod', &
                'graupel tracer not found, but do_graupel is true', FATAL)
          end if

          if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)') 'Stratiform cloud tracer ind&
                &ices: nql,nqi,nqa =',nql,nqi,nqa
          if (min(nql,nqi,nqa) <= 0)   &
             call error_mesg ('cloud_spec_mod', &
             'stratiform cloud tracer(s) not found', FATAL)
          if (nql == nqi .or. nqa == nqi .or. nql == nqa)   &
              call error_mesg ('cloud_spec_mod',  &
            'tracers indices cannot be the same (i.e., nql=nqi=nqa).', &
                                                              FATAL)
          nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
          if (nqn /= NO_TRACER)  then
            Cldrad_control%do_liq_num = .true.
          else
            Cldrad_control%do_liq_num = .false.
          endif
          nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
          if (nqni /= NO_TRACER)  then
            Cldrad_control%do_ice_num = .true.
          else
            Cldrad_control%do_ice_num = .false.
          endif
        else
          Cldrad_control%do_liq_num = .false.
          Cldrad_control%do_ice_num = .false.
        endif

        Cldrad_control%do_liq_num_iz = .true.
        Cldrad_control%do_ice_num_iz = .true.

!---------------------------------------------------------------------
!    define the variables indicating that the cloud parameterization
!    control variables have been defined.
!---------------------------------------------------------------------
      Cldrad_control%do_rh_clouds_iz = .true.
      Cldrad_control%do_strat_clouds_iz = .true.
      Cldrad_control%do_zonal_clouds_iz = .true.
      Cldrad_control%do_mgroup_prescribed_iz = .true.
      Cldrad_control%do_obs_clouds_iz = .true.
      Cldrad_control%do_no_clouds_iz = .true.
      Cldrad_control%do_diag_clouds_iz = .true.
      Cldrad_control%do_specified_clouds_iz = .true.
      Cldrad_control%do_specified_strat_clouds_iz = .true.
      Cldrad_control%do_donner_deep_clouds_iz = .true.
      Cldrad_control%do_uw_clouds_iz = .true.
      Cldrad_control%do_zetac_clouds_iz = .true.
      Cldrad_control%do_stochastic_clouds_iz = .true.
 
!--------------------------------------------------------------------
!    if stochastic clouds is active, allocate and define arrays holding
!    the processor's latitudes and longitudes. be sure that the
!    cloud_generator module has been initialized.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds_iz) then
        if (Cldrad_control%do_stochastic_clouds) then
          allocate (lats(size(latb,1),size(latb,2)))
          allocate (lons(size(lonb,1), size(lonb,2)))
          lats(:,:) = latb(:,:)*radian
          lons(:,:) = lonb(:,:)*radian
          call cloud_generator_init

!---------------------------------------------------------------------
!     determine the source of the random number seed generator to be used
!     for the stochastic cloud generation. the legacy scheme may fail to
!     provide spacially unique seeds at hi-res (above c48) or if the time
!     provided the radiation package does not monotonically advance (as in
!     some specialized sensitivity / assessment studies). 
!---------------------------------------------------------------------
          if (do_legacy_seed_generation) then
!---------------------------------------------------------------------
!     if it is desired to force the use of the temperature-based
!     random number seed (as is used when time is not always advancing
!     as seen by the radiation package, or when model resolution is
!     less than 1 degree), set the logical control variable in 
!     Cldrad_control to so indicate. 
!---------------------------------------------------------------------
            if ( force_use_of_temp_for_seed) then
              Cldrad_control%use_temp_for_seed = .true.
              Cldrad_control%use_temp_for_seed_iz = .true.
              call error_mesg ('cloud_spec_init', &
                 'Will use temp as basis for stochastic cloud seed; &
                    &force_use_of_temp_for_seed is set true', NOTE)
            else
              call error_mesg ('cloud_spec_init', &
               ' If model resolution is above c48, it is &
               &HIGHLY RECOMMENDED that you set cloud_spec_nml variable &
               &force_use_of_temp_for_seed to true to assure &
               &reproducibility across pe count and domain layout', NOTE)
              call error_mesg ('cloud_spec_init', &
               'No action is needed at or below c48 resolution.', NOTE)
            endif

!---------------------------------------------------------------------
!     if the latitude and longitude of adjacent points on a pe have the 
!     same integral values (NINT), set the logical control variable in 
!     Cldrad_control to use the model temperature at the top level as the 
!     random number seed to provide spacial uniqueness at the points on the
!     processor. 
!     Note that for model resolutions of ~ 1 degree, some pes may use
!     lat and lon, while others use temperature, and that this may change
!     as the domain decomposition or npes used for the problem are changed.
!     Therefore it is  HIGHLY RECOMMENDED that for  resolutions above c48  
!     that nml variable force_use_of_temp_for_seed be set to .true.; 
!     it may remain the default value of .false. for lower resolution runs
!     or to preserve legacy results, or if reproducibility over npes or 
!     layout is not essential. 
!---------------------------------------------------------------------
            if (.not. Cldrad_control%use_temp_for_seed) then
  jLoop:      do j=1,jd
                do i=1,id
                  do jj=j+1,jd+1
                    do ii=i+1,id+1
                      if (NINT(lats(ii,jj)) == NINT(lats(i,j))) then
                        if (NINT(lons(ii,jj)) == NINT(lons(i,j))) then    
                          Cldrad_control%use_temp_for_seed = .true.
                          Cldrad_control%use_temp_for_seed_iz = .true.
                          call error_mesg ('cloud_spec_init', &
                           'Found grid point within 1 degree of  &
                             &another',NOTE)
                          call error_mesg ('cloud_spec_init', &
                            'if reproducibility across npes and layout is &
                           &desired, you must set cloud_spec_nml variable &
                           &force_use_of_temp_for_seed to true., and &
                           &restart the model.', NOTE)
                          exit jLoop
                        endif
                      endif
                    end do
                  end do
                end do
              end do jLoop
            endif

!-------------------------------------------------------------------------
!    set seed generation source to be the temperature field (current 
!    default).
!-------------------------------------------------------------------------
          else
            Cldrad_control%use_temp_for_seed = .true.
            Cldrad_control%use_temp_for_seed_iz = .true.
          endif
        endif

      else
        call error_mesg ('microphys_rad_mod', &
         ' attempt to use Cldrad_control%do_stochastic_clouds before &
                                                &it is defined', FATAL)
      endif
 
!--------------------------------------------------------------------
!    include do_fu2007 in the cloudrad_control_type variable for use
!    in other modules.
!--------------------------------------------------------------------
      Cldrad_control%using_fu2007 = do_fu2007
      Cldrad_control%using_fu2007_iz = .true.     

!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------



end subroutine cloud_spec_init



!######################################################################
! <SUBROUTINE NAME="cloud_spec">
!  <OVERVIEW>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_radiative_properties defines the cloud radiative properties 
!    appropriate for the radiation options that are active.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time,
!                       Atmos_input, &
!                       Surface, Cld_spec, Lsc_microphys,  &
!                       Meso_microphys, Cell_microphys, lsc_area_in, &
!                       lsc_liquid_in, lsc_ice_in, lsc_droplet_number_in        , r)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation calculation is to apply
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </INOUT>
!  <INOUT NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </INOUT>
!  <IN NAME="lsc_liquid_in" TYPE="real">
!   OPTIONAL: lsc cloud water mixing ratio  present when running 
!    standalone columns or sa_gcm
!  </IN>
!  <IN NAME="lsc_ice_in" TYPE="real">
!   OPTIONAL: cloud ice mixing ratio  present when running 
!    standalone columns or sa_gcm
!  </IN>
!  <IN NAME="lsc_area_in" TYPE="real">
!   OPTIONAL: fractional cloud area, present when running 
!                        standalone columns or sa_gcm
!  </IN>
!  <IN NAME="r" TYPE="real">
!   OPTIONAL: model tracer fields on the current time step
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time, &
                       Atmos_input, Surface, Cld_spec, Lsc_microphys, &
                       Meso_microphys, Cell_microphys,  &
                       Shallow_microphys, lsc_area_in, lsc_liquid_in, &
                       lsc_ice_in, lsc_droplet_number_in,   &
                       lsc_ice_number_in,   &
                       lsc_snow_in, lsc_rain_in,   &
                       lsc_snow_size_in,  lsc_rain_size_in,  r,  &
                       shallow_cloud_area, shallow_liquid, shallow_ice,&
                       shallow_droplet_number, &
                       shallow_ice_number, &
                       cell_cld_frac, cell_liq_amt, cell_liq_size, &
                       cell_ice_amt, cell_ice_size, &
                       cell_droplet_number, &
                       meso_cld_frac, meso_liq_amt, meso_liq_size, &
                       meso_ice_amt, meso_ice_size,  &
                       meso_droplet_number, nsum_out)

!----------------------------------------------------------------------
!    cloud_spec specifies the cloud field seen by the radiation package.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
real, dimension(:,:),         intent(in)             :: lat
real, dimension(:,:,:),       intent(in)             :: z_half, z_full
type(time_type),              intent(in)             :: Rad_time
type(atmos_input_type),       intent(inout)          :: Atmos_input
type(surface_type),           intent(inout)          :: Surface       
type(cld_specification_type), intent(inout)          :: Cld_spec    
type(microphysics_type),      intent(inout)          :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys,&
                                                     Shallow_microphys
real, dimension(:,:,:),       intent(in), optional ::  &
                                      lsc_liquid_in, lsc_ice_in, &
                                      lsc_droplet_number_in, lsc_area_in, &
                                      lsc_ice_number_in
real, dimension(:,:,:),       intent(in), optional ::   &
                                      lsc_snow_in, lsc_rain_in,  &
                                      lsc_snow_size_in,  lsc_rain_size_in
real, dimension(:,:,:,:),     intent(in),   optional :: r
real, dimension(:,:,:),       intent(inout),optional :: &
                      shallow_cloud_area, shallow_liquid, shallow_ice,&
                        shallow_droplet_number, &
                           cell_cld_frac, cell_liq_amt, cell_liq_size, &
                           cell_ice_amt, cell_ice_size, &
                           cell_droplet_number, &
                           meso_cld_frac, meso_liq_amt, meso_liq_size, &
                           meso_ice_amt, meso_ice_size, &
                           meso_droplet_number,  &
                           shallow_ice_number
integer, dimension(:,:),      intent(inout), optional:: nsum_out

!-------------------------------------------------------------------
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      lat               latitude of model points  [ radians ]
!      z_half            height asl at half levels [ m ]
!      z_full            height asl at full levels [ m ]
!      Rad_time          time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ] 
!
!   intent(inout) variables:
!
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ] 
!      Surface           variables defining the surface albedo and land
!                        fraction
!                        [ surface_type ]
!      Cld_spec          variables on the model grid which define all or
!                        some of the following, dependent on the 
!                        specific cloud parameterization: cloud optical 
!                        paths, particle sizes, cloud fractions, cloud 
!                        thickness, number of clouds in a column, 
!                        and /or cloud type (high/mid/low, ice/liq or 
!                        random/max overlap)
!                        [ cld_specification_type ]
!      Lsc_microphys     variables describing the microphysical proper-
!                        ties of the large-scale clouds
!                        [ microphysics_type ]
!      Meso_microphys    variables describing the microphysical proper-
!                        ties of the meso-scale clouds
!                        [ microphysics_type ]
!      Cell_microphys    variables describing the microphysical proper-
!                        ties of the convective cell-scale clouds
!                        [ microphysics_type ]
!
!   intent(in), optional variables:
!
!      lsc_liquid_in     cloud water mixing ratio (or specific humidity 
!                        ????), present when running standalone columns
!                        or sa_gcm
!                        [ non-dimensional ]
!      lsc_ice_in        cloud ice mixing ratio (or specific humidity 
!                         ????), present when running standalone columns
!                        or sa_gcm
!                        [ non-dimensional ]
!      lsc_area_in       fractional cloud area, present when running 
!                        standalone columns or sa_gcm
!                        [ non-dimensional ]
!      r                 model tracer fields on the current time step
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: ix, jx, kx
      integer   :: ierr
      logical   :: override
      type(time_type) :: Data_time
      real, dimension (size (Atmos_input%deltaz,1), &
                       size (Atmos_input%deltaz,2), &
                       size (Atmos_input%deltaz,3)) :: rho

!---------------------------------------------------------------------
!   local variables:
!
!        ix      number of grid points in x direction (on processor)
!        jx      number of grid points in y direction (on processor)
!        kx      number of model layers
!        rho     atmospheric density [ kg / m**3 ]
!        ierr
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    check for the presence of optional input arguments.
!---------------------------------------------------------------------
      ierr = 1
      if (present(lsc_liquid_in) .and.   &
          present(lsc_ice_in) .and. &
          present(lsc_area_in)) then
        ierr = 0
        if (Cldrad_control%do_liq_num) then
          if (.not. present(lsc_droplet_number_in) ) then
            call error_mesg ('cloud_spec_mod', &
               'must input lsc_droplet_number_in when  &
                                   &using do_liq_num', FATAL)
          endif
        endif
      else
        if (present (r)) then
          ierr = 0
        else
          call error_mesg ('cloud_spec_mod', &
              'must input either r or lsc_liquid_in and  &
              &lsc_ice_in when using predicted cloud microphysics', &
                                                               FATAL)
        endif
      endif

      if (Cldrad_control%do_uw_clouds) then
        if ( present (shallow_cloud_area) .and. &
             present (shallow_liquid) .and.  &
             present (shallow_ice)  .and.  &
             present (shallow_droplet_number) .and. &
             present (shallow_ice_number) ) then
        else
          call error_mesg ('cloud_spec_mod',  &
           'optional argument(s) required when shallow clouds &
                                           &are active is missing', FATAL)
        endif
      endif

      if (Cldrad_control%do_donner_deep_clouds) then
        if ( present (cell_cld_frac) .and. &
             present (cell_liq_amt) .and.  &
             present (cell_liq_size)  .and.  &
             present (cell_ice_amt) .and. &
             present (cell_ice_size) .and.  &
             present (cell_droplet_number) .and. &
             present (meso_cld_frac) .and. &
             present (meso_liq_amt) .and.  &
             present (meso_liq_size)  .and.  &
             present (meso_ice_amt) .and. &
             present (meso_ice_size) .and.  &
             present (meso_droplet_number) .and. &
             present (nsum_out) ) then
        else
          call error_mesg ('cloud_spec_mod',  &
           'optional argument(s) required when donner clouds &
                                           &are active is missing', FATAL)
        endif
      endif

!----------------------------------------------------------------------
!    define model dimensions.
!----------------------------------------------------------------------
      ix = size(Atmos_input%deltaz,1)
      jx = size(Atmos_input%deltaz,2)
      kx = size(Atmos_input%deltaz,3)

!----------------------------------------------------------------------
!    call initialize_cldamts to allocate and initialize the arrays
!    contained in the structures used to specify the cloud amounts, 
!    types and locations and the microphysical parameters.
!----------------------------------------------------------------------
      call initialize_cldamts (ix, jx, kx, Lsc_microphys,   &
                               Meso_microphys, Cell_microphys, &
                               Shallow_microphys, Cld_spec)

!---------------------------------------------------------------------
!    define the cloud_water, cloud_ice and cloud_area components of 
!    Cld_spec.
!---------------------------------------------------------------------
      if (present (lsc_ice_in) .and. &
          present (lsc_liquid_in) ) then
        Cld_spec%cloud_ice   = lsc_ice_in
        Cld_spec%cloud_water = lsc_liquid_in
        if (Cldrad_control%do_liq_num) then
          if (present(lsc_droplet_number_in)) then
             Cld_spec%cloud_droplet (:,:,:) = lsc_droplet_number_in
          endif
        endif
        if (Cldrad_control%do_ice_num) then
          if (present(lsc_ice_number_in)) then
            Cld_spec%cloud_ice_num (:,:,:) = lsc_ice_number_in
          end if
        endif
      endif
      if (present (lsc_area_in)  )then
        Cld_spec%cloud_area = lsc_area_in
      endif
      if (present (lsc_snow_in)  )then
        Cld_spec%snow = lsc_snow_in
      endif
      if (present (lsc_rain_in)  )then
        Cld_spec%rain = lsc_rain_in
      endif
      if (present (lsc_snow_size_in)  )then
        Cld_spec%snow_size = lsc_snow_size_in
      endif
      if (present (lsc_rain_size_in)  )then
        Cld_spec%rain_size = lsc_rain_size_in
      endif

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call the appropriate subroutine to define the cloud
!    location, type, amount or whatever other arrays the particular 
!    parameterization uses to specify its clouds. if the model is being
!    run with do_no_clouds = .true., exit from this routine, leaving
!    the cloud specification variables as they were initialized (to a
!    condition of no clouds).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!---------------------------------------------------------------------
!    when running in standalone columns mode, call standalone_clouds_amt
!    to obtain the cloud specification variables. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_specified_clouds .or. &
            Cldrad_control%do_specified_strat_clouds )   then

          call standalone_clouds_amt (is, ie, js, je, lat,     &
                                      Atmos_input%press, Cld_spec)

!---------------------------------------------------------------------
!    if the rh diagnostic cloud scheme is active, call rh_clouds_amt
!    to define the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_rh_clouds) then
          call rh_clouds_amt (is, ie, js, je, Atmos_input%press, lat,  &
                              Cld_spec)

!---------------------------------------------------------------------
!    if either zonal clouds or obs clouds is active, call 
!    specified_clouds_amt to obtain the needed cloud specification
!    variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_zonal_clouds .or. &
                 Cldrad_control%do_obs_clouds) then
          call specified_clouds_amt (is, ie, js, je, Rad_time, lat,    &
                                     Atmos_input%pflux, Cld_spec)

!---------------------------------------------------------------------
!    if mgrp_prscr_clds is active, call prscr_clds_amt to obtain the 
!    needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_mgroup_prescribed) then
          call prscr_clds_amt (is, ie, js, je, Cld_spec)

!----------------------------------------------------------------------
!    if gordon diagnostic clouds are active, call diag_clouds_amt to 
!    obtain the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_diag_clouds) then
          call diag_clouds_amt (is, ie, js, je, lat, Atmos_input%pflux,&
                                Atmos_input%press, Rad_time, Cld_spec, &
                                Lsc_microphys)

!----------------------------------------------------------------------
!    if zetac clouds are active, call zetac_clouds_amt to 
!    obtain the needed cloud specification variables.
!---------------------------------------------------------------------
        else if (Cldrad_control%do_zetac_clouds) then
          if (present (r)) then
            Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
            Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
            Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
!           ierr = 0
!         else
!           call error_mesg ('cloud_spec_mod', &
!             ' must pass tracer array r when using zetac clouds', &
!                                                             FATAL)
          endif
          call zetac_clouds_amt (is, ie, js, je, z_half, z_full, &
                                 Surface%land, Atmos_input%phalf, &
                                 Atmos_input%deltaz, Cld_spec, &
                                 Lsc_microphys)

        endif ! (do_rh_clouds)
!--------------------------------------------------------------------
!    if klein prognostic clouds are active, call strat_clouds_amt to 
!    obtain the needed cloud specification variables.
!--------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then

!---------------------------------------------------------------------
!    if the gcm is being executed, call strat_cloud_avg to obtain the
!    appropriate (either instantaneous or time-averaged) values of
!    cloud water, cloud ice and cloud fraction. if the sa_gcm or the
!    standalone columns mode is being executed with the strat cloud
!    option, then values for the cloud water, cloud ice and when needed
!    cloud area have been input as optional arguments to this sub-
!    routine.
!---------------------------------------------------------------------
          if(present(lsc_liquid_in)) then
            if (Cld_spec%cloud_area(1,1,1) == -99.) then
              if (present (r)) then
                Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
                Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
                Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
                if (Cldrad_control%do_liq_num) then
                  Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
                endif
              else
                call error_mesg ('cloud_spec_mod', &
                   'lsc_area_in, etc not present in restart file (flag &
                    &has been set), and r array is not passed to &
                                 &cloud_spec; cannot proceed', FATAL)
              endif
            endif
          else  ! (present (lsc_liquid_in))
            if (present (r)) then
              Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
              Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
              Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
              if (Cldrad_control%do_liq_num) then
                Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
              endif
            else
              call error_mesg ('cloud_spec_mod', &
                  'neither lsc_area_in, etc nor r array &
                 &has been passed to cloud_spec; cannot proceed', FATAL)
            endif
          endif ! (present(lsc_liquid_in))

          if(present(lsc_ice_number_in)) then
            if (Cld_spec%cloud_ice_num(1,1,1) == -99.) then
              if (present (r)) then
                if (Cldrad_control%do_ice_num) then
                  Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
                endif
              else
                call error_mesg ('cloud_spec_mod', &
                   'lsc_ice_number_in not present in restart file (flag &
                    &has been set), and r array is not passed to &
                                 &cloud_spec; cannot proceed', FATAL)
              endif
            endif
          else  ! (present (lsc_ice_number_in))
            if (Cldrad_control%do_ice_num) then
              if (present (r)) then
                Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
              else
                call error_mesg ('cloud_spec_mod', &
                  'neither lsc_ice_number_in nor r array &
                 &has been passed to cloud_spec; cannot proceed', FATAL)
              endif
            endif
          endif ! (present(lsc_ice_number_in))

          if (present(r)) then
            if (do_rain) then !sjl
              Cld_spec%cloud_water(:,:,:) = Cld_spec%cloud_water(:,:,:)+r(:,:,:,nqr)
            end if
            if (do_snow) then !miz
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqs)
            end if
            if (do_graupel) then !SJL
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqg)
            end if
          endif  

!---------------------------------------------------------------------
!    if the cloud input data is to be overriden, define the time slice
!    of data which is to be used. allocate storage for the cloud data.
!---------------------------------------------------------------------
          if (doing_data_override) then
            Data_time = Rad_time +    &
                          set_time (Rad_control%rad_time_step, 0)
 
!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    water data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qlnew', Cld_spec%cloud_water,   &
                                Data_time, override=override,           &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'ql => cloud_water not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    ice data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current 
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qinew', Cld_spec%cloud_ice,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'qi => cloud_ice   not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    fraction data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qanew', Cld_spec%cloud_area,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'qa => cloud_area not overridden successfully', FATAL)
            endif
            ierr = 0
          endif ! (doing_override)

!---------------------------------------------------------------------
!    if values for the cloud variables have been successfully obtained,
!    call strat_clouds_amt to define the appropriate cloud specification
!    variables.
!---------------------------------------------------------------------
          if (ierr == 0) then
            call strat_clouds_amt (is, ie, js, je, Rad_time, &
                                   Atmos_input%pflux,  &
                                   Atmos_input%press,   &
                                   Atmos_input%cloudtemp, &
                                   Atmos_input%cloudvapor(:,:,:)/  &
                                   (1.0+Atmos_input%cloudvapor(:,:,:)), &
                                   Surface%land,&
                                   Cld_spec, Lsc_microphys)

!----------------------------------------------------------------------
!    if ierr is non-zero, then cloud data was not successfully obtained.
!    if this is not the coldstart step, write an error message and 
!    stop execution.
!----------------------------------------------------------------------
          else 
            if (num_pts >= tot_pts) then
              call error_mesg ('cloud_spec_mod',  &
                     'no strat cloud data available; ierr /= 0', FATAL)

!----------------------------------------------------------------------
!    if this is the coldstart step, retain the input values corres-
!    ponding to no clouds, increment the points counter, and continue. 
!----------------------------------------------------------------------
            else
              num_pts = num_pts + size(Atmos_input%press,1)*   &
                                  size(Atmos_input%press,2)
            endif
          endif
        endif ! (do_strat_clouds)

!--------------------------------------------------------------------
!    since donner_deep_clouds may be active along with strat clouds, 
!    the associated properties are determined outside of the above loop.
!    these properties are placed in Cell_microphys and Meso_microphys.
!----------------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then
          call donner_deep_clouds_amt (is, ie, js, je,  &
                           cell_cld_frac, cell_liq_amt, cell_liq_size, &
                           cell_ice_amt, cell_ice_size, &
                           cell_droplet_number, &
                           meso_cld_frac, meso_liq_amt, meso_liq_size, &
                           meso_ice_amt, meso_ice_size,  &
                           meso_droplet_number,  nsum_out, &
                           Cell_microphys, Meso_microphys)

!---------------------------------------------------------------------
!    convert the cloud and ice amounts from kg(h2o) / kg(air) to 
!    g(h2o) / m**3, as required for use in the microphys_rad routines
!    which compute cloud radiative properties.
!---------------------------------------------------------------------
          rho(:,:,:) = Atmos_input%press(:,:,1:kx)/  &
                       (RDGAS*Atmos_input%temp(:,:,1:kx))
          Cell_microphys%conc_drop = 1.0e03*rho*Cell_microphys%conc_drop
          Cell_microphys%conc_ice  = 1.0e03*rho*Cell_microphys%conc_ice 
          Meso_microphys%conc_drop = 1.0e03*rho*Meso_microphys%conc_drop
          Meso_microphys%conc_ice  = 1.0e03*rho*Meso_microphys%conc_ice 
        endif

!--------------------------------------------------------------------
!    since uw_clouds may be active along with strat clouds and / or 
!    donner deep clouds, the associated properties are determined 
!    outside of the above loop. these properties are placed in  
!    Shallow_microphys.
!----------------------------------------------------------------------
         if (Cldrad_control%do_uw_clouds) then
           call uw_clouds_amt (is, ie, js, je,  &
                            shallow_cloud_area, shallow_liquid, &
                            shallow_ice, shallow_droplet_number, &
                            shallow_ice_number,  &
                            Surface%land, Atmos_input%press,  &
                            Atmos_input%cloudtemp, Shallow_microphys)
        endif

!---------------------------------------------------------------------
!    obtain the microphysical properties (sizes and concentrations) if
!    a prescribed microphysics scheme is active. 
!---------------------------------------------------------------------
        if (Cldrad_control%do_presc_cld_microphys) then
          call microphys_presc_conc (is, ie, js, je,   &
                                     Atmos_input%clouddeltaz,   &
                                     Atmos_input%cloudtemp, &
                                     Cld_spec, Lsc_microphys)
        endif

!---------------------------------------------------------------------
!    call combine_cloud_properties to combine (if necessary) the cloud 
!    properties from multiple cloud types (large-scale, donner deep,
!    uw shallow) into a single set for use by the radiation package. 
!    this is only needed when microphysically-based properties are 
!    present, and when either strat clouds, donner deep and / or uw
!    shallow clouds is activated.
!---------------------------------------------------------------------
        if ( .not. Cldrad_control%do_specified_strat_clouds ) then
          if (Cldrad_control%do_sw_micro  .or.    &
              Cldrad_control%do_lw_micro) then
            if (Cldrad_control%do_strat_clouds .or.    &
                Cldrad_control%do_uw_clouds .or.    &
                Cldrad_control%do_donner_deep_clouds) then
              call combine_cloud_properties ( is, js,  &
                                             Atmos_input%temp(:,:,1), &
                                             Rad_time, &
                                             Lsc_microphys,    &
                                             Meso_microphys,   &
                                             Cell_microphys,   &
                                             Shallow_microphys, &
                                             Cld_spec)
            endif
          endif
        endif
      endif  !  (.not. do_no_clouds)

!--------------------------------------------------------------------
!    if microphysics is active and strat_clouds is not, define the water
!    paths (in units of kg / m**2).  if strat_clouds is active, these 
!    values will have already been defined. when microphysics is active,
!    define the effective sizes for the liquid and ice particles.
!--------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro .or.    &
          Cldrad_control%do_sw_micro)  then
        if (.not. Cldrad_control%do_strat_clouds .and.   &
            .not. Cldrad_control%do_zetac_clouds ) then
          Cld_spec%lwp = 1.0E-03*Lsc_microphys%conc_drop(:,:,:)* &
                         Atmos_input%clouddeltaz(:,:,:)
          Cld_spec%iwp = 1.0E-03*Lsc_microphys%conc_ice(:,:,:)*  &
                         Atmos_input%clouddeltaz(:,:,:)
        endif
        Cld_spec%reff_liq_micro = Lsc_microphys%size_drop
        Cld_spec%reff_ice_micro = Lsc_microphys%size_ice
      endif

!---------------------------------------------------------------------


end subroutine cloud_spec    



!######################################################################
! <SUBROUTINE NAME="cloud_spec_dealloc">
!  <OVERVIEW>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys, Cell_microphys and
!    Shallow_microphys.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys, Cell_microphys and
!    Shallow_microphys.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_dealloc (Cld_spec, Lsc_microphys, Meso_microphys,&
!                               Cell_microphys, Shallow_microphys)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </INOUT>
! </SUBROUTINE>
! 
subroutine cloud_spec_dealloc (Cld_spec, Lsc_microphys, Meso_microphys,&
                               Cell_microphys, Shallow_microphys)

!---------------------------------------------------------------------
!    cloud_spec_dealloc deallocates the component arrays of the 
!    cld_specification_type structure Cld_spec and the microphysics_type
!    structures Lsc_microphys, Meso_microphys, Cell_microphys and
!    Shallow_microphys.
!----------------------------------------------------------------------

type(cld_specification_type), intent(inout) :: Cld_spec
type(microphysics_type),      intent(inout) :: Lsc_microphys,   &
                                               Meso_microphys, &
                                               Cell_microphys, &
                                               Shallow_microphys


!----------------------------------------------------------------------
!    deallocate the array elements of Cld_spec.
!----------------------------------------------------------------------
      deallocate (Cld_spec%camtsw         )  
      deallocate (Cld_spec%cmxolw         )
      deallocate (Cld_spec%crndlw         )
      deallocate (Cld_spec%ncldsw         )
      deallocate (Cld_spec%nmxolw         )
      deallocate (Cld_spec%nrndlw         )
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%camtsw_band    )
        deallocate (Cld_spec%ncldsw_band    )
        deallocate (Cld_spec%cld_thickness_sw_band  )
        deallocate (Cld_spec%lwp_sw_band            )
        deallocate (Cld_spec%iwp_sw_band            )
        deallocate (Cld_spec%reff_liq_sw_band       )
        deallocate (Cld_spec%reff_ice_sw_band       )
        deallocate (Cld_spec%stoch_cloud_type       )
      endif
      if (Cldrad_control%do_stochastic_clouds) then
        deallocate (Cld_spec%crndlw_band    )
        deallocate (Cld_spec%nrndlw_band    )
        deallocate (Cld_spec%cld_thickness_lw_band  )
        deallocate (Cld_spec%lwp_lw_band            )
        deallocate (Cld_spec%iwp_lw_band            )
        deallocate (Cld_spec%reff_liq_lw_band       )
        deallocate (Cld_spec%reff_ice_lw_band       )
      endif
      deallocate (Cld_spec%tau            )
      deallocate (Cld_spec%lwp            )
      deallocate (Cld_spec%iwp            )
      deallocate (Cld_spec%reff_liq       )
      deallocate (Cld_spec%reff_ice       )
      deallocate (Cld_spec%reff_liq_lim   )
      deallocate (Cld_spec%reff_ice_lim   )
      deallocate (Cld_spec%reff_liq_micro )
      deallocate (Cld_spec%reff_ice_micro )
      deallocate (Cld_spec%liq_frac       )
      deallocate (Cld_spec%cld_thickness  )
      deallocate (Cld_spec%hi_cloud       )
      deallocate (Cld_spec%mid_cloud      )
      deallocate (Cld_spec%low_cloud      )
      deallocate (Cld_spec%ice_cloud      )
      deallocate (Cld_spec%cloud_water    )
      deallocate (Cld_spec%cloud_ice      )
      deallocate (Cld_spec%cloud_area     )
      deallocate (Cld_spec%cloud_droplet  )
      deallocate (Cld_spec%cloud_ice_num  )
      deallocate (Cld_spec%snow           )
      deallocate (Cld_spec%rain           )
      deallocate (Cld_spec%snow_size      )
      deallocate (Cld_spec%rain_size      )

!--------------------------------------------------------------------
!    deallocate the elements of Lsc_microphys.
!---------------------------------------------------------------------
      deallocate (Lsc_microphys%conc_drop   )
      deallocate (Lsc_microphys%conc_ice    )
      deallocate (Lsc_microphys%conc_rain   )
      deallocate (Lsc_microphys%conc_snow   )
      deallocate (Lsc_microphys%size_drop   )
      deallocate (Lsc_microphys%size_ice    )
      deallocate (Lsc_microphys%size_rain   )
      deallocate (Lsc_microphys%size_snow   )
      deallocate (Lsc_microphys%cldamt      )
      deallocate (Lsc_microphys%droplet_number )
      deallocate (Lsc_microphys%ice_number  )
      if (Cldrad_control%do_stochastic_clouds) then
        nullify (Lsc_microphys%lw_stoch_conc_drop   )
        nullify (Lsc_microphys%lw_stoch_conc_ice    )
        nullify (Lsc_microphys%lw_stoch_size_drop   )
        nullify (Lsc_microphys%lw_stoch_size_ice    )
        nullify (Lsc_microphys%lw_stoch_cldamt      )
        nullify (Lsc_microphys%lw_stoch_droplet_number)
        nullify (Lsc_microphys%lw_stoch_ice_number  )

        nullify (Lsc_microphys%sw_stoch_conc_drop   )
        nullify (Lsc_microphys%sw_stoch_conc_ice    )
        nullify (Lsc_microphys%sw_stoch_size_drop   )
        nullify (Lsc_microphys%sw_stoch_size_ice    )
        nullify (Lsc_microphys%sw_stoch_cldamt      )
        nullify (Lsc_microphys%sw_stoch_droplet_number)
        nullify (Lsc_microphys%sw_stoch_ice_number  )

        deallocate (Lsc_microphys%stoch_conc_drop   )
        deallocate (Lsc_microphys%stoch_conc_ice    )
        deallocate (Lsc_microphys%stoch_size_drop   )
        deallocate (Lsc_microphys%stoch_size_ice    )
        deallocate (Lsc_microphys%stoch_cldamt      )
        deallocate (Lsc_microphys%stoch_droplet_number )
        deallocate (Lsc_microphys%stoch_ice_number  )
      endif

!--------------------------------------------------------------------
!    deallocate the elements of Cell_microphys.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        deallocate (Cell_microphys%conc_drop   )
        deallocate (Cell_microphys%conc_ice    )
        deallocate (Cell_microphys%conc_rain   )
        deallocate (Cell_microphys%conc_snow   )
        deallocate (Cell_microphys%size_drop   )
        deallocate (Cell_microphys%size_ice    )
        deallocate (Cell_microphys%size_rain   )
        deallocate (Cell_microphys%size_snow   )
        deallocate (Cell_microphys%cldamt      )
        deallocate (Cell_microphys%droplet_number )

!--------------------------------------------------------------------
!    deallocate the elements of Meso_microphys.
!---------------------------------------------------------------------
        deallocate (Meso_microphys%conc_drop   )
        deallocate (Meso_microphys%conc_ice    )
        deallocate (Meso_microphys%conc_rain   )
        deallocate (Meso_microphys%conc_snow   )
        deallocate (Meso_microphys%size_drop   )
        deallocate (Meso_microphys%size_ice    )
        deallocate (Meso_microphys%size_rain   )
        deallocate (Meso_microphys%size_snow   )
        deallocate (Meso_microphys%cldamt      )
        deallocate (Meso_microphys%droplet_number )
      endif

!--------------------------------------------------------------------
!    deallocate the elements of Shallow_microphys.
!---------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds) then
        deallocate (Shallow_microphys%conc_drop   )
        deallocate (Shallow_microphys%conc_ice    )
        deallocate (Shallow_microphys%conc_rain   )
        deallocate (Shallow_microphys%conc_snow   )
        deallocate (Shallow_microphys%size_drop   )
        deallocate (Shallow_microphys%size_ice    )
        deallocate (Shallow_microphys%size_rain   )
        deallocate (Shallow_microphys%size_snow   )
        deallocate (Shallow_microphys%cldamt      )
        deallocate (Shallow_microphys%droplet_number )
      endif

!---------------------------------------------------------------------


end subroutine cloud_spec_dealloc 



!#####################################################################

subroutine cloud_spec_end

!---------------------------------------------------------------------
!    cloud_spec_end is the destructor for cloud_spec_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    close the modules that were initialized by this module.
!--------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!-------------------------------------------------------------------
!    mgroup prescribed clouds.
!-------------------------------------------------------------------
        if (Cldrad_control%do_mgroup_prescribed) then
          call mgrp_prscr_clds_end 

!-------------------------------------------------------------------
!    rh-based diagnostic clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_rh_clouds) then
          call rh_based_clouds_end 

!-------------------------------------------------------------------
!    zonal or observed clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_zonal_clouds .or.  &
                 Cldrad_control%do_obs_clouds)  then
          call specified_clouds_W_end             

!-------------------------------------------------------------------
!    klein predicted clouds. if this option is active, donner_deep 
!    clouds may also be active. additionally, this may also have been
!    activated when running in standalone columns mode, in which case
!    standalone_clouds_end must be called.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_strat_clouds) then
          if (Cldrad_control%do_specified_strat_clouds .or. &
              Cldrad_control%do_specified_clouds ) then 
            call standalone_clouds_end
          endif

!-------------------------------------------------------------------
!    gordon diagnostic clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_diag_clouds) then
          call diag_clouds_W_end

!-------------------------------------------------------------------
!    zetac clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_zetac_clouds) then
          call zetac_clouds_W_end

!-------------------------------------------------------------------
!    standalone specified clouds.
!-------------------------------------------------------------------
        else if (Cldrad_control%do_specified_clouds) then
          call standalone_clouds_end
        endif

!------------------------------------------------------------------
!    cloud types which may coexist must be processed outside of if loop
!------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then
          call strat_clouds_W_end
        endif
        if (Cldrad_control%do_donner_deep_clouds) then
          call donner_deep_clouds_W_end
        endif
        if (Cldrad_control%do_uw_clouds) then
          call uw_clouds_W_end
        endif
      endif  ! (not do_no_clouds)

!--------------------------------------------------------------------
!    mark the module as no longer initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine cloud_spec_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#####################################################################
! <SUBROUTINE NAME="initialize_cldamts">
!  <OVERVIEW>
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!  </OVERVIEW>
!  <DESCRIPTION>
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_cldamts (ix, jx, kx, Lsc_microphys,    &
!                               Meso_microphys, Cell_microphys, Cld_spec)
!  </TEMPLATE>
!  <IN NAME="ix, jx, kx" TYPE="integer">
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!  </IN>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </INOUT>
! </SUBROUTINE>
! 
subroutine initialize_cldamts (ix, jx, kx, Lsc_microphys,    &
                               Meso_microphys, Cell_microphys,  &
                               Shallow_microphys, Cld_spec)

!---------------------------------------------------------------------
!    initialize_cldamts allocates and initializes the array components 
!    of the structures used to specify the model cloud and microphysics
!    fields.
!---------------------------------------------------------------------

integer,                      intent(in)     :: ix, jx, kx
type(microphysics_type),      intent(inout)  :: Lsc_microphys,   &
                                                Meso_microphys, &
                                                Cell_microphys, &
                                                Shallow_microphys
type(cld_specification_type), intent(inout)  :: Cld_spec

!----------------------------------------------------------------------
!    intent(in) variables:
! 
!       ix             size of i dimension of physics window
!       jx             size of j dimension of physics window
!       kx             size of k dimension of physics window
!
!    intent(inout) variables:
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
!
!            the following elements are components of these  
!            microphysics_type variables which are allocated and
!            initialized here:
!
!            %conc_ice   ice particle concentration [ g / m**3 ]
!            %conc_drop  cloud droplet concentration [ g / m**3 ]
!            %conc_rain  rain drop concentration [ g / m**3 ]
!            %conc_snow  snow concentration [ g / m**3 ]
!            %size_ice   effective ice crystal diameter [ microns ]
!            %size_drop  effective cloud drop diameter [ microns ]
!            %size_rain  effective rain drop diameter [ microns ]
!            %size_snow  effective snow flake diameter [ microns ]
!            %cldamt     total cloud fraction (crystal + droplet)
!                        [ dimensionless ]
!            %lw_stoch_conc_ice
!                        ice particle concentration as a function of
!                        lw parameterization band [ g / m**3 ]
!            %lw_stoch_conc_drop
!                        cloud droplet concentration as a function of
!                        lw parameterization band [ g / m**3 ]
!            %lw_stoch_size_ice
!                        effective ice crystal diameter as a function
!                        of lw parameterization band [ microns ]
!            %lw_stoch_size_drop
!                        effective cloud drop diameter as a function of
!                        lw parameterization band [ microns ]
!            %lw_stoch_cldamt
!                        total cloud fraction (crystal + droplet) as a
!                        function of lw parameterization band
!                        [ dimensionless ]
!            %sw_stoch_conc_ice
!                        ice particle concentration as a function of
!                        sw parameterization band [ g / m**3 ]
!            %sw_stoch_conc_drop
!                        cloud droplet concentration as a function of
!                        sw parameterization band [ g / m**3 ]
!            %sw_stoch_size_ice
!                        effective ice crystal diameter as a function
!                        of sw parameterization band [ microns ]
!            %sw_stoch_size_drop
!                        effective cloud drop diameter as a function of
!                        sw parameterization band [ microns ]
!            %sw_stoch_cldamt
!                        total cloud fraction (crystal + droplet) as a
!                        function of sw parameterization band
!                        [ dimensionless ]
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths, 
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type 
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!            the following elements are components of this  
!            cld_specification_type variable which is allocated and
!            initialized here:
!
!            %cmxolw         amount of maximally overlapped longwave 
!                            clouds [ dimensionless ]
!            %crndlw         amount of randomly overlapped longwave 
!                            clouds [ dimensionless ]
!            %nmxolw         number of maximally overlapped longwave 
!                            clouds in each grid column.
!            %nrndlw         number of maximally overlapped longwave 
!                            clouds in each grid column.
!            %camtsw         shortwave cloud amount. the sum of the max-
!                            imally overlapped and randomly overlapped 
!                            longwave cloud amounts. [ dimensionless ]
!            %ncldsw         number of shortwave clouds in each grid 
!                            column.
!            %camtsw_band    shortwave cloud amount. the sum of the max-
!                            imally overlapped and randomly overlapped
!                            longwave cloud amounts, differing with sw
!                            parameterization band. [ dimensionless ]
!            %crndlw_band    amount of randomly overlapped longwave
!                            clouds, differing with lw parameterization
!                            band [ dimensionless ]
!            %hi_cloud       logical mask for high clouds 
!            %mid_cloud      logical mask for middle clouds
!            %low_cloud      logical mask for low clouds
!            %ice_cloud      logical mask for ice clouds
!            %iwp            ice water path  [ kg / m**2 ]
!            %lwp            liquid water path [ kg / m**2 ]
!            %reff_liq       effective cloud drop radius  used with
!                            bulk cloud physics scheme [ microns ]
!            %reff_ice       effective ice crystal radius used with
!                            bulk cloud physics scheme [ microns ]
!            %reff_liq_micro effective cloud drop radius used with 
!                            microphysically based scheme [ microns ]
!            %reff_ice_micro effective ice crystal radius used with
!                            microphysically based scheme [ microns ]
!            %tau            extinction optical path  [ dimensionless ]
!            %liq_frac       fraction of cloud in a box which is liquid
!                            [ dimensionless ]
!            %cld_thickness  number of model layers contained in cloud  
!            %cloud_water    liquid cloud content [ kg liq / kg air ]
!            %cloud_ice      ice cloud content [ kg ice / kg air ]
!            %cloud_area     saturated volume fraction [ dimensionless ]
!
!---------------------------------------------------------------------

      integer  :: n

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    large-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      allocate (Lsc_microphys%conc_drop  (ix, jx, kx) )
      allocate (Lsc_microphys%conc_ice   (ix, jx, kx) )
      allocate (Lsc_microphys%conc_rain  (ix, jx, kx) )
      allocate (Lsc_microphys%conc_snow  (ix, jx, kx) )
      allocate (Lsc_microphys%size_drop  (ix, jx, kx) )
      allocate (Lsc_microphys%size_ice   (ix, jx, kx) )
      allocate (Lsc_microphys%size_rain  (ix, jx, kx) )
      allocate (Lsc_microphys%size_snow  (ix, jx, kx) )
      allocate (Lsc_microphys%cldamt     (ix, jx, kx) )
      allocate (Lsc_microphys%droplet_number (ix, jx, kx) )
       allocate (Lsc_microphys%ice_number (ix, jx, kx) )
      Lsc_microphys%conc_drop(:,:,:) = 0.
      Lsc_microphys%conc_ice(:,:,:)  = 0.
      Lsc_microphys%conc_rain(:,:,:) = 0.
      Lsc_microphys%conc_snow(:,:,:) = 0.
      Lsc_microphys%size_drop(:,:,:) = 1.0e-20 
      Lsc_microphys%size_ice(:,:,:)  = 1.0e-20 
      Lsc_microphys%size_rain(:,:,:) = 1.0e-20        
      Lsc_microphys%size_snow(:,:,:) = 1.0e-20 
      Lsc_microphys%cldamt(:,:,:)    = 0.0
      Lsc_microphys%droplet_number(:,:,:)    = 0.0
      Lsc_microphys%ice_number(:,:,:)= 0.0
      if (Cldrad_control%do_stochastic_clouds) then
        allocate (Lsc_microphys%stoch_conc_drop  &
                                  (ix, jx, kx, Cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_conc_ice &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_size_drop  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_size_ice   &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_cldamt  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_droplet_number  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
        allocate (Lsc_microphys%stoch_ice_number  &
                                  (ix, jx, kx, cldrad_control%nlwcldb + Solar_spect%nbands) )
  
        Lsc_microphys%lw_stoch_conc_drop => Lsc_microphys%stoch_conc_drop(:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_conc_drop => Lsc_microphys%stoch_conc_drop(:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_conc_ice  => Lsc_microphys%stoch_conc_ice (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_conc_ice  => Lsc_microphys%stoch_conc_ice (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_size_drop => Lsc_microphys%stoch_size_drop(:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_size_drop => Lsc_microphys%stoch_size_drop(:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_size_ice  => Lsc_microphys%stoch_size_ice (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_size_ice  => Lsc_microphys%stoch_size_ice (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_cldamt    => Lsc_microphys%stoch_cldamt   (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_cldamt    => Lsc_microphys%stoch_cldamt   (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_droplet_number    => Lsc_microphys%stoch_droplet_number   (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_droplet_number    => Lsc_microphys%stoch_droplet_number   (:, :, :, Cldrad_control%nlwcldb+1:)
        Lsc_microphys%lw_stoch_ice_number    => Lsc_microphys%stoch_ice_number   (:, :, :, 1:Cldrad_control%nlwcldb)
        Lsc_microphys%sw_stoch_ice_number    => Lsc_microphys%stoch_ice_number   (:, :, :, Cldrad_control%nlwcldb+1:)

       do n=1,Cldrad_control%nlwcldb + Solar_spect%nbands
        Lsc_microphys%stoch_conc_drop(:,:,:,n) = 0.
        Lsc_microphys%stoch_conc_ice(:,:,:,n)  = 0.
        Lsc_microphys%stoch_size_drop(:,:,:,n) = 1.0e-20
        Lsc_microphys%stoch_size_ice(:,:,:,n)  = 1.0e-20
        Lsc_microphys%stoch_cldamt(:,:,:,n)    = 0.0
        Lsc_microphys%stoch_droplet_number(:,:,:,n) = 0.0
        Lsc_microphys%stoch_ice_number(:,:,:,n)= 0.0
       end do
      endif

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    meso-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
      allocate (Meso_microphys%conc_drop  (ix, jx, kx) )
      allocate (Meso_microphys%conc_ice   (ix, jx, kx) )
      allocate (Meso_microphys%conc_rain  (ix, jx, kx) )
      allocate (Meso_microphys%conc_snow  (ix, jx, kx) )
      allocate (Meso_microphys%size_drop  (ix, jx, kx) )
      allocate (Meso_microphys%size_ice   (ix, jx, kx) )
      allocate (Meso_microphys%size_rain  (ix, jx, kx) )
      allocate (Meso_microphys%size_snow  (ix, jx, kx) )
      allocate (Meso_microphys%cldamt     (ix, jx, kx) )
      allocate (Meso_microphys%droplet_number   (ix, jx, kx) )
      Meso_microphys%conc_drop = 0.
      Meso_microphys%conc_ice  = 0.
      Meso_microphys%conc_rain = 0.
      Meso_microphys%conc_snow = 0.
      Meso_microphys%size_drop = 1.0e-20 
      Meso_microphys%size_ice  = 1.0e-20  
      Meso_microphys%size_rain = 1.0e-20                            
      Meso_microphys%size_snow = 1.0e-20 
      Meso_microphys%cldamt    = 0.0                               
      Meso_microphys%droplet_number   = 0.0                               
      endif

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    cell-scale cloud scheme, including any precipitation fields and
!    the total cloud fraction. concentrations and fractions are init-
!    ialized to 0.0, and the effective sizes are set to small numbers to
!    avoid potential divides by zero.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
      allocate (Cell_microphys%conc_drop  (ix, jx, kx) )
      allocate (Cell_microphys%conc_ice   (ix, jx, kx) )
      allocate (Cell_microphys%conc_rain  (ix, jx, kx) )
      allocate (Cell_microphys%conc_snow  (ix, jx, kx) )
      allocate (Cell_microphys%size_drop  (ix, jx, kx) )
      allocate (Cell_microphys%size_ice   (ix, jx, kx) )
      allocate (Cell_microphys%size_rain  (ix, jx, kx) )
      allocate (Cell_microphys%size_snow  (ix, jx, kx) )
      allocate (Cell_microphys%cldamt     (ix, jx, kx) )
      allocate (Cell_microphys%droplet_number  (ix, jx, kx) )
      Cell_microphys%conc_drop = 0.
      Cell_microphys%conc_ice  = 0.
      Cell_microphys%conc_rain = 0.
      Cell_microphys%conc_snow = 0.
      Cell_microphys%size_drop = 1.0e-20 
      Cell_microphys%size_ice  = 1.0e-20  
      Cell_microphys%size_rain = 1.0e-20                          
      Cell_microphys%size_snow = 1.0e-20 
      Cell_microphys%cldamt     = 0.
      Cell_microphys%droplet_number    = 0.
      endif

!---------------------------------------------------------------------
!    allocate the arrays defining the microphysical parameters of the 
!    clouds of the shallow convection scheme, including any precip-
!    itation fields and the total cloud fraction. concentrations and 
!    fractions are initialized to 0.0, and the effective sizes are set 
!    to small numbers to avoid potential divides by zero.
!---------------------------------------------------------------------
      if (Cldrad_control%do_uw_clouds) then
        allocate (Shallow_microphys%conc_drop  (ix, jx, kx) )
        allocate (Shallow_microphys%conc_ice   (ix, jx, kx) )
        allocate (Shallow_microphys%conc_rain  (ix, jx, kx) )
        allocate (Shallow_microphys%conc_snow  (ix, jx, kx) )
        allocate (Shallow_microphys%size_drop  (ix, jx, kx) )
        allocate (Shallow_microphys%size_ice   (ix, jx, kx) )
        allocate (Shallow_microphys%size_rain  (ix, jx, kx) )
        allocate (Shallow_microphys%size_snow  (ix, jx, kx) )
        allocate (Shallow_microphys%cldamt     (ix, jx, kx) )
        allocate (Shallow_microphys%droplet_number  (ix, jx, kx) )
        Shallow_microphys%conc_drop = 0.
       Shallow_microphys%conc_ice  = 0.
       Shallow_microphys%conc_rain = 0.
       Shallow_microphys%conc_snow = 0.
       Shallow_microphys%size_drop = 1.0e-20
       Shallow_microphys%size_ice  = 1.0e-20
       Shallow_microphys%size_rain = 1.0e-20
       Shallow_microphys%size_snow = 1.0e-20
       Shallow_microphys%cldamt    = 0.0
       Shallow_microphys%droplet_number    = 0.
     endif

!---------------------------------------------------------------------
!    allocate arrays to hold the cloud fractions seen by the shortwave
!    and the random and maximum overlap fractions seen by the longwave
!    radiation, and then the number of each of these types of cloud in
!    each column. initialize the cloud fractions and number of clouds
!    to zero.
!---------------------------------------------------------------------
      allocate ( Cld_spec%camtsw (ix, jx, kx ) )
      allocate ( Cld_spec%cmxolw (ix, jx, kx ) )
      allocate ( Cld_spec%crndlw (ix, jx, kx ) )
      allocate ( Cld_spec%ncldsw (ix, jx     ) )
      allocate ( Cld_spec%nmxolw (ix, jx     ) )
      allocate ( Cld_spec%nrndlw (ix, jx     ) )
      Cld_spec%cmxolw(:,:,:) = 0.0E+00
      Cld_spec%crndlw(:,:,:) = 0.0E+00
      Cld_spec%camtsw(:,:,:) = 0.0E+00
      Cld_spec%nmxolw (:,:)  = 0
      Cld_spec%nrndlw (:,:)  = 0
      Cld_spec%ncldsw (:,:)  = 0
      if (Cldrad_control%do_stochastic_clouds) then
        allocate ( Cld_spec%camtsw_band    &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate ( Cld_spec%ncldsw_band    &
                                 (ix, jx, Solar_spect%nbands) )
        allocate (Cld_spec%cld_thickness_sw_band & 
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%iwp_sw_band    &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%lwp_sw_band   &
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%reff_liq_sw_band   &  
                                 (ix, jx, kx, Solar_spect%nbands) )
        allocate (Cld_spec%reff_ice_sw_band   &   
                                 (ix, jx, kx, Solar_spect%nbands) )
        do n=1,Solar_spect%nbands
        Cld_spec%camtsw_band(:,:,:,n) = 0.0E+00
        Cld_spec%ncldsw_band(:,:,n) = 0
        Cld_spec%cld_thickness_sw_band(:,:,:,n) = 0              
        Cld_spec%lwp_sw_band(:,:,:,n)   = 0.0
        Cld_spec%iwp_sw_band(:,:,:,n)   = 0.0
        Cld_spec%reff_liq_sw_band(:,:,:,n)      = 10.0
        Cld_spec%reff_ice_sw_band(:,:,:,n)      = 30.0
        end do
        allocate ( Cld_spec%crndlw_band    &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate ( Cld_spec%nrndlw_band    &
                                 (ix, jx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%cld_thickness_lw_band & 
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%iwp_lw_band    &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%lwp_lw_band   &
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%reff_liq_lw_band   &  
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        allocate (Cld_spec%reff_ice_lw_band   &   
                                 (ix, jx, kx, Cldrad_control%nlwcldb) )
        do n=1, Cldrad_control%nlwcldb
        Cld_spec%crndlw_band(:,:,:,n) = 0.0E+00
        Cld_spec%nrndlw_band(:,:,n) = 0
        Cld_spec%cld_thickness_lw_band(:,:,:,n) = 0              
        Cld_spec%lwp_lw_band(:,:,:,n)   = 0.0
        Cld_spec%iwp_lw_band(:,:,:,n)   = 0.0
        Cld_spec%reff_liq_lw_band(:,:,:,n)      = 10.0
        Cld_spec%reff_ice_lw_band (:,:,:,n)     = 30.0
        end do
        allocate (Cld_spec%stoch_cloud_type   &
            (ix, jx, kx, Solar_spect%nbands + Cldrad_control%nlwcldb) )
        Cld_spec%stoch_cloud_type = 0
      endif

!--------------------------------------------------------------------
!    allocate and initialize various arrays that are used by one or
!    another cloud scheme to specify the cloud locations and amounts.
!    initialization provides values consistent with the absence of
!    cloud, with the exception of the particle size fields which are
!    set to small, non-zero values.
!---------------------------------------------------------------------
      allocate (Cld_spec%hi_cloud       (ix, jx, kx) )
      allocate (Cld_spec%mid_cloud      (ix, jx, kx) )
      allocate (Cld_spec%low_cloud      (ix, jx, kx) )
      allocate (Cld_spec%ice_cloud      (ix, jx, kx) )
      allocate (Cld_spec%iwp            (ix, jx, kx) )
      allocate (Cld_spec%lwp            (ix, jx, kx) )
      allocate (Cld_spec%reff_liq       (ix, jx, kx) )
      allocate (Cld_spec%reff_ice       (ix, jx, kx) )
      allocate (Cld_spec%reff_liq_lim   (ix, jx, kx) )
      allocate (Cld_spec%reff_ice_lim   (ix, jx, kx) )
      allocate (Cld_spec%reff_liq_micro (ix, jx, kx) )
      allocate (Cld_spec%reff_ice_micro (ix, jx, kx) )
      allocate (Cld_spec%tau            (ix, jx, kx, num_slingo_bands) )
      allocate (Cld_spec%liq_frac       (ix, jx, kx) )
      allocate (Cld_spec%cld_thickness  (ix, jx, kx) )
      allocate (Cld_spec%cloud_water    (ix,jx,kx) )
      allocate (Cld_spec%cloud_ice      (ix,jx,kx)   )
      allocate (Cld_spec%cloud_area     (ix,jx,kx)  )
      allocate (Cld_spec%cloud_droplet  (ix,jx,kx)  )
      allocate (Cld_spec%cloud_ice_num  (ix,jx,kx)  )
      allocate (Cld_spec%snow           (ix,jx,kx)  )
      allocate (Cld_spec%rain           (ix,jx,kx)  )
      allocate (Cld_spec%snow_size      (ix,jx,kx)  )
      allocate (Cld_spec%rain_size      (ix,jx,kx)  )

      Cld_spec%hi_cloud (:,:,:)     = .false.
      Cld_spec%mid_cloud(:,:,:)     = .false.
      Cld_spec%low_cloud(:,:,:)     = .false.
      Cld_spec%ice_cloud(:,:,:)     = .false.
      Cld_spec%lwp(:,:,:)    = 0.0
      Cld_spec%iwp(:,:,:)    = 0.0
      Cld_spec%reff_liq(:,:,:)      = 10.0
      Cld_spec%reff_ice(:,:,:)      = 30.0
      Cld_spec%reff_liq_lim(:,:,:)      = 10.0
      Cld_spec%reff_ice_lim(:,:,:)      = 30.0
      Cld_spec%reff_liq_micro(:,:,:) = 10.0
      Cld_spec%reff_ice_micro(:,:,:) = 30.0
      Cld_spec%liq_frac(:,:,:)      = 0.0
      Cld_spec%cld_thickness(:,:,:) = 0              
      Cld_spec%cloud_water(:,:,:)   = 0.
      Cld_spec%cloud_ice(:,:,:)     = 0.
      Cld_spec%cloud_area(:,:,:)    = 0.
      Cld_spec%cloud_droplet(:,:,:)    = 0.
      Cld_spec%cloud_ice_num(:,:,:) = 0.
      Cld_spec%snow(:,:,:)          = 0.
      Cld_spec%rain(:,:,:)          = 0.
      Cld_spec%snow_size(:,:,:)     = 0.
      Cld_spec%rain_size(:,:,:)     = 0.
      do n=1, num_slingo_bands
      Cld_spec%tau(:,:,:,n)  = 0.0
      end do

!---------------------------------------------------------------------



end subroutine initialize_cldamts              



!###################################################################
! <SUBROUTINE NAME="combine_cloud_properties">
!  <OVERVIEW>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </OVERVIEW>
!  <DESCRIPTION>
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call combine_cloud_properties (Lsc_microphys, Meso_microphys,  &
!                                     Cell_microphys, Cld_spec)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for 
!                        clouds associated with uw shallow convection
!  </IN>
! </SUBROUTINE>
! 
subroutine combine_cloud_properties (is, js, temp, Rad_time, &
                                     Lsc_microphys, Meso_microphys,  &
                                     Cell_microphys, Shallow_microphys,&
                                     Cld_spec)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud specification property 
!    arrays for the total cloud field in each grid box, using as input 
!    the specification of the component cloud types that may be present
!    (large-scale, donner mesoscale and cell-scale, uw shallow).
!----------------------------------------------------------------------

integer, intent(in)  :: is, js
real, dimension(:,:), intent(in) :: temp
type(time_type), intent(in) :: Rad_time
type(microphysics_type),        intent(in)    :: Lsc_microphys, &
                                                 Meso_microphys, &
                                                 Cell_microphys, &
                                                 Shallow_microphys
type(cld_specification_type), intent(inout)   :: Cld_spec

!----------------------------------------------------------------------
!   intent(in) variables:
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
!
!    intent(inout) variables:
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths, 
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type 
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    variables for folding Donner cloud properties into stochastic
!    cloud arrays
!------------------------------------------------------------------
      type(randomNumberStream),   &
                    dimension(size(Lsc_microphys%cldamt,1),   &
                              size(Lsc_microphys%cldamt,2)) :: streams
      real, &            
                    dimension(size(Lsc_microphys%cldamt,1),   &
                              size(Lsc_microphys%cldamt,2),   &       
                              size(Lsc_microphys%cldamt,3),   &       
                              size(Lsc_microphys%stoch_cldamt, 4))  :: &
                                                     randomNumbers
      real    :: seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)
      integer :: nn, nsubcols

      integer :: i, j, k, n

!---------------------------------------------------------------------
!    total-cloud specification properties need be defined only when
!    strat_cloud, donner_deep and/or uw shallow clouds are active.
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds .and.    &
          Cldrad_control%do_uw_clouds .and.   &
          Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud, donner_deep and uw shallow are all active, define 
!    the random overlap cloud fraction as the sum of the fractions of 
!    the large-scale, donner meso-scale and cell-scale, and the
!    uw shallow clouds.
!---------------------------------------------------------------------
       Cld_spec%crndlw = Lsc_microphys%cldamt +   &
                         Cell_microphys%cldamt +  &
                       Meso_microphys%cldamt + Shallow_microphys%cldamt
     else if (Cldrad_control%do_strat_clouds .and.    &
          Cldrad_control%do_donner_deep_clouds) then


!----------------------------------------------------------------------
!    if strat_cloud and donner_deep are both active, define the random
!    overlap cloud fraction as the sum of the fractions of the large-
!    scale, meso-scale and cell-scale clouds.
!---------------------------------------------------------------------
        Cld_spec%crndlw = Lsc_microphys%cldamt +    &
                          Cell_microphys%cldamt + Meso_microphys%cldamt

       else if (Cldrad_control%do_strat_clouds .and.    &
           Cldrad_control%do_uw_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud and uw_shallow  are both active, define the random
!    overlap cloud fraction as the sum of the fractions of the large-
!    scale and uw shallow clouds.
!---------------------------------------------------------------------
         Cld_spec%crndlw = Lsc_microphys%cldamt +    &
                           Shallow_microphys%cldamt

    else if (Cldrad_control%do_uw_clouds .and.    &
             Cldrad_control%do_donner_deep_clouds) then
 
!----------------------------------------------------------------------
!    if strat_cloud and donner_deep are both active, define the random
!    overlap cloud fraction as the sum of the fractions of the large-
!    scale, meso-scale and cell-scale clouds.
!---------------------------------------------------------------------
         Cld_spec%crndlw = Shallow_microphys%cldamt +    &
                          Cell_microphys%cldamt + Meso_microphys%cldamt

!------------------------------------------------------------------
!    if strat cloud is activated but donner_deep is not, define the 
!    total-cloud amount to be the large scale cloud amount.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_strat_clouds) then
        Cld_spec%crndlw = Lsc_microphys%cldamt

!---------------------------------------------------------------------
!    if donner_deep is active but strat cloud is not, then the mesoscale
!    and cell-scale cloud amounts are combined to define the total
!    random-overlap cloud fraction.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_donner_deep_clouds) then
        Cld_spec%crndlw = Cell_microphys%cldamt + Meso_microphys%cldamt
      else if (Cldrad_control%do_uw_clouds) then
        Cld_spec%crndlw = Shallow_microphys%cldamt
      endif

!---------------------------------------------------------------------
!    randomly-overlapped clouds are being assumed for donner_deep and 
!    strat cloud module clouds. set the max overlap cloud fraction to 
!    zero, be certain that the random overlap fraction is .le. 1. after
!    the summing of the component cloud fractions, and define the total
!    cloud fraction to be used by the sw code.
!---------------------------------------------------------------------
      Cld_spec%cmxolw = 0.0
      Cld_spec%crndlw = MIN (Cld_spec%crndlw, 1.00)
      Cld_spec%camtsw = Cld_spec%crndlw

!--------------------------------------------------------------------
!    if stochastic clouds are being used, define the cloud type to be 
!    seen by the radiation code in each stochastic subcolumn.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then

!--------------------------------------------------------------------
!   assign either a 1 or a 0 to each subcolumn indicating whether
!   lsc cloud is present or not. 
!--------------------------------------------------------------------
        nsubcols = Solar_spect%nbands + Cldrad_control%nlwcldb
        do n=1,nsubcols
          if ( n > Solar_spect%nbands) then
            nn = n - Solar_spect%nbands    
          else
            nn = n + Cldrad_control%nlwcldb    
          endif
          do k=1,size(Lsc_microphys%cldamt,3) ! Levels
            do j=1,size(Lsc_microphys%cldamt,2) ! Lons
              do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                if (Lsc_microphys%stoch_cldamt(i,j,k,nn) > 0.) then
 
!----------------------------------------------------------------------
!    fill it in with the large-scale cloud values.
!----------------------------------------------------------------------
                  Cld_spec%stoch_cloud_type(i,j,k,n) = 1 
               else
                  Cld_spec%stoch_cloud_type(i,j,k,n) = 0  
               endif
             end do
            end do
          end do
        end do

!----------------------------------------------------------------------
!    compare the cell and meso-scale cloud amounts to a random number, 
!    and replace the large-scale cloud and clear sky assignment in each 
!    subcolumn with an assignment of cell or meso-scale clouds when the 
!    number is less than the cloud fraction. use the maximum overlap 
!    assumption. treat the random number as the location with the PDF 
!    of total water. cells are at the top of the PDF; then meso-scale 
!    anvils, then large-scale clouds and clear sky.
!------------------------------------------------------------
        if (Cldrad_control%do_donner_deep_clouds) then     
          if (Cldrad_control%use_temp_for_seed) then
            do j=1,size(Lsc_microphys%cldamt,2)
              do i=1,size(Lsc_microphys%cldamt,1)
                streams(i,j) = &
                  initializeRandomNumberStream (                      &
                    ishftc(nint(temp(i,j)*seedwts),1))
              end do
            end do
          else
            do j=1,size(Lsc_microphys%cldamt,2)
              do i=1,size(Lsc_microphys%cldamt,1)
                streams(i,j) = &
                  initializeRandomNumberStream (                      &
                    constructSeed(nint(lons(is+i-1,js+j-1)),           &
                                  nint(lats(is+i-1,js+j-1)), Rad_time, &
                                  perm = 1))
              end do
            end do
          endif
 
!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
          do j=1,size(Lsc_microphys%cldamt,2) ! Lons
            do i=1,size(Lsc_microphys%cldamt,1) ! Lats
              call getRandomNumbers (streams(i,j),    &
                                                randomNumbers(i,j,1,:))
            end do
          end do
 
!----------------------------------------------------------------------
!    here is maximum overlap. we use a 3D arrary for the random numbers
!    for flexibility.
!----------------------------------------------------------------------
          do k=2,size(Lsc_microphys%cldamt,3)
            randomNumbers(:,:,k,:) = randomNumbers(:,:,1,:)
          end do
 
!----------------------------------------------------------------------
!    assign cloud types, band by band
!----------------------------------------------------------------------
    IF (ignore_donner_cells) then
          do n=1,nsubcols    
            do k=1,size(Lsc_microphys%cldamt,3) ! Levels
              do j=1,size(Lsc_microphys%cldamt,2) ! Lons
                do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                  if (randomNumbers(i,j,k,n) >    &
                           (1. -  Meso_microphys%cldamt(i, j, k))) then
 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                    Cld_spec%stoch_cloud_type(i,j,k,n) = 2 
                  endif
                end do
              end do
            end do
          end do
    ELSE
          do n=1,nsubcols    
            do k=1,size(Lsc_microphys%cldamt,3) ! Levels
              do j=1,size(Lsc_microphys%cldamt,2) ! Lons
                do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                  if (randomNumbers(i,j,k,n) >     &
                            (1. - Cell_microphys%cldamt(i,j,k))) then

!----------------------------------------------------------------------
!    it's a cell.
!----------------------------------------------------------------------
                    Cld_spec%stoch_cloud_type(i,j,k,n) = 3  
                  else if (randomNumbers(i,j,k,n) >    &
                           (1. - Cell_microphys%cldamt(i, j, k) - &
                                  Meso_microphys%cldamt(i, j, k))) then
 
!----------------------------------------------------------------------
!    it's a meso-scale.
!----------------------------------------------------------------------
                    Cld_spec%stoch_cloud_type(i,j,k,n) = 2 
                  endif
                end do
              end do
            end do
          end do
    ENDIF
        endif

!----------------------------------------------------------------------
!    compare the uw shallow cloud amount to a random number, and replace
!    the donner cloud, large-scale cloud or clear sky previously 
!    assigned in each subcolumn with an assignment of uw shallow cloud 
!    when the number is less than the cloud fraction. use the maximum 
!    overlap assumption. treat the random number as the location with 
!    the PDF of total water. uw shallow clouds are at the top of this 
!    PDF, then large-scale clouds and clear sky.
!------------------------------------------------------------
        if (Cldrad_control%do_uw_clouds) then     
          if (Cldrad_control%use_temp_for_seed) then
            do j=1,size(Lsc_microphys%cldamt,2)
              do i=1,size(Lsc_microphys%cldamt,1)
                streams(i,j) = &
                  initializeRandomNumberStream (                      &
                    ishftc(nint(temp(i,j)*seedwts),2))
              end do
            end do
          else
            do j=1,size(Lsc_microphys%cldamt,2)
              do i=1,size(Lsc_microphys%cldamt,1)
                streams(i,j) = &
                  initializeRandomNumberStream (                      &
                    constructSeed(nint(lons(is+i-1,js+j-1)),           &
                                  nint(lats(is+i-1,js+j-1)), Rad_time, &
                                  perm = 2))
              end do
            end do
          endif
 
!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
          do j=1,size(Lsc_microphys%cldamt,2) ! Lons
            do i=1,size(Lsc_microphys%cldamt,1) ! Lats
              call getRandomNumbers (streams(i,j),    &
                                                randomNumbers(i,j,1,:))
            end do
          end do
 
!----------------------------------------------------------------------
!    here is maximum overlap. we use a 3D arrary for the random numbers
!    for flexibility.
!----------------------------------------------------------------------
          do k=2,size(Lsc_microphys%cldamt,3)
            randomNumbers(:,:,k,:) = randomNumbers(:,:,1,:)
          end do
 
!----------------------------------------------------------------------
!    assign cloud type, band by band
!----------------------------------------------------------------------
          do n=1,nsubcols
            do k=1,size(Lsc_microphys%cldamt,3) ! Levels
              do j=1,size(Lsc_microphys%cldamt,2) ! Lons
                do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                  if (randomNumbers(i,j,k,n) >     &
                            (1. - Shallow_microphys%cldamt(i,j,k))) then

!----------------------------------------------------------------------
!    it's a uw shallow.
!----------------------------------------------------------------------
                    Cld_spec%stoch_cloud_type(i,j,k,n) = 4  
                  endif
                end do
              end do
            end do
          end do
        endif

!---------------------------------------------------------------------
!     define the cloud amount in each stochastic subcolumn to be either
!     1.0 if cloud is present, or 0.0 if no cloud exists.
!---------------------------------------------------------------------
        do n=1,Solar_spect%nbands
          do k=1,size(Lsc_microphys%cldamt,3) ! Levels
            do j=1,size(Lsc_microphys%cldamt,2) ! Lons
              do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,n) /= 0) then  
                  Cld_spec%camtsw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%camtsw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do
        
        do n=1,Cldrad_control%nlwcldb
          nn = Solar_spect%nbands + n
          do k=1,size(Lsc_microphys%cldamt,3) ! Levels
            do j=1,size(Lsc_microphys%cldamt,2) ! Lons
              do i=1,size(Lsc_microphys%cldamt,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,nn) /= 0) then  
                  Cld_spec%crndlw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%crndlw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do
      endif  ! (do_stochastic)


!#####################################################################


end subroutine combine_cloud_properties 



!###################################################################
! <SUBROUTINE NAME="microphs_presc_conc">
!  <OVERVIEW>
!   Subroutine to determine water droplet and ice crystal based on
!   prescribed microphysics model.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine uses prescribed microphysics model to determine
!   concentrations of water droplets and ice crystals. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
!                                 Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending indice of the x dimension in the physics domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting indice of the y dimension in the physics domain
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending indice of the y dimension in the physics domain 
!  </IN>
!  <IN NAME="deltaz" TYPE="real">
!   Height of each pressure layers.
!  </IN>
!  <IN NAME="temp" TYPE="real">
!   Temperatures of pressure levels
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </INOUT>
! </SUBROUTINE>
!
subroutine microphys_presc_conc (is, ie, js, je, deltaz, temp,      &
                                 Cld_spec, Lsc_microphys)

!---------------------------------------------------------------------
!    microphys_presc_conc defines microphysical properties based on the
!    assumption of specified total water paths for high, middle and low 
!    clouds.
!---------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
real, dimension(:,:,:),       intent(in)     :: deltaz, temp  
type(cld_specification_type), intent(in)     :: Cld_spec
type(microphysics_type),      intent(inout)  :: Lsc_microphys

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      deltaz         model vertical grid separation that is to be used
!                     for cloud calculations
!                     [meters]
!      temp           temperature at model levels (1:nlev) that is to
!                     be used in cloud calculations
!                     [ deg K ]
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution
!
!   intent(inout) variables:
!
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

      real,    dimension(size(temp,1), size(temp,2), size(temp,3)) :: &
                                                       conc

      integer, dimension(size(temp,1), size(temp,2)) :: &
                                                       nhi_clouds, &
                                                       nmid_clouds, &
                                                       nlow_clouds

      integer  :: i,j,k

!--------------------------------------------------------------------
!  local variables:
!
!      conc             droplet concentration  [ g / m**3 ]
!      nhi_clouds       number of layers with high clouds
!      nmid_clouds      number of layers with middle clouds
!      nlow_clouds      number of layers with low clouds
!      i,j,k            do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! RSH NOTE:
!
!    THE FOLLOWING treatment of diag_cloud_mod is here as an INITIAL 
!    IMPLEMENTATION to allow compilation and model execution, and 
!    provide "reasonable ?? " values.
! 
!    Code developed but NOT YET ADDED HERE reflects a later approach. 
!    That code is available under the fez release, and will be added to
!    the repository when upgrades to the cloud-radiation modules are 
!    completed.
!
!    obtain drop and ice size and concentrations, consistent with 
!    the diag_cloud scheme. As a test case, the following is a simple 
!    specification of constant concentration and size in all boxes 
!    defined as cloudy, attempting to come close to the prescribed 
!    values used for other cloud schemes. assume ice cld thickness 
!    = 2.0 km; then conc_ice=10.0E-03 => iwp = 20 g/m^2, similar to that
!    prescribed in microphys_presc_conc. assume water cld thickness 
!    = 3.5 km; then conc_drop = 20E-03 => lwp = 70 g / m^2, similar to 
!    that prescribed in microphys_presc_conc.  use sizes as used in 
!    microphys_presc_conc (50 and 20 microns). when done, radiative 
!    boundary fluxes are "similar" to non-microphysical results
!    for test case done here, and shows reasonable sensitivity to
!    variations in concentrations.
!    AGAIN, THIS IS AN INITIAL IMPLEMENTATION FOR TESTING ONLY !!!!
!
!!! RSH
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!---------------------------------------------------------------------

      if (Cldrad_control%do_diag_clouds) then
!---------------------------------------------------------------------
!    define concentrations and sizes of ice at those points which have
!    condensate and which were previously determined to be ice points
!    (Cld_spec%ice_cloud), and define concentrations and sizes of
!    liquid droplets at those points with condensate that had been 
!    determined to support liquid clouds.
!---------------------------------------------------------------------
        do k=1,size(Cld_spec%camtsw,3)
          do j=1,size(Cld_spec%camtsw,2)
            do i=1,size(Cld_spec%camtsw,1)
              if (Cld_spec%camtsw(i,j,k) > 0.0) then
                if (Cld_spec%ice_cloud(i,j,k)) then
                  Lsc_microphys%conc_ice(i,j,k) = 10.0E-03  
                  Lsc_microphys%size_ice(i,j,k) = 50.     
                else
                  Lsc_microphys%conc_drop(i,j,k) = 20.0E-03
                  Lsc_microphys%size_drop(i,j,k) = 20.
                endif
              endif
            end do
          end do
        end do
!----------------------------------------------------------------------
!    for the non-diag_cloud_mod cases, assume that the water path is 
!    preset at fixed values (lwpath_hi, _mid, _low) for "high", "mid", 
!    "low" clouds. the lwpath in each cloud layer within "hi", "mid" 
!    "low" pressure intervals is that lwpath_... divided by the number 
!    of clouds present in that pressure interval.
!----------------------------------------------------------------------
      else

!----------------------------------------------------------------------
!    define the number of high, middle, low clouds according to
!    Wetherald's criterion.
!----------------------------------------------------------------------
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            nhi_clouds(i,j)  = 0
            nmid_clouds(i,j) = 0
            nlow_clouds(i,j) = 0
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) &
                               nhi_clouds(i,j)  =  nhi_clouds(i,j)  + 1
              if (Cld_spec%mid_cloud(i,j,k)) &
                               nmid_clouds(i,j) =  nmid_clouds(i,j) + 1
              if (Cld_spec%low_cloud(i,j,k))  &
                               nlow_clouds(i,j) =  nlow_clouds(i,j) + 1
            end do
          end do
        end do

!----------------------------------------------------------------------
!    compute the water substance concentration in each layer 
!    (as water path / layer geometric path).
!----------------------------------------------------------------------
        conc(:,:,:) = 0.0E+00
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            do k=1,size(Cld_spec%camtsw,3)
              if (Cld_spec%hi_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_hi/   &
                              (nhi_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%mid_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_mid/    &
                              (nmid_clouds(i,j)*deltaz(i,j,k))
              endif
              if (Cld_spec%low_cloud(i,j,k)) then
                conc(i,j,k) = lwpath_low    /                   &
                              (nlow_clouds(i,j)*deltaz(i,j,k))
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    split conc into conc_ice and conc_drop, depending on temperature
!    criterion (T < 273.16). assume that rain and / or snow are not
!    present.
!----------------------------------------------------------------------
        do k=1,size(Cld_spec%camtsw,3)
          do j=1,size(Cld_spec%camtsw,2)
            do i=1,size(Cld_spec%camtsw,1)
              if (temp(i,j,k) .LT. 273.16) then
                Lsc_microphys%conc_ice(i,j,k) = conc(i,j,k)
              else
                Lsc_microphys%conc_drop(i,j,k) = conc(i,j,k)
              endif
            end do
          end do
        end do

!----------------------------------------------------------------------
!    define sizes of microphysical species, using namelist values. note
!    that namelist drop and rain sizes are radii, so multiply by 2 to 
!    produce diameter, as desired for the %size_ arrays.
!----------------------------------------------------------------------
        Lsc_microphys%size_drop(:,:,:) = 2.0*wtr_cld_reff
        Lsc_microphys%size_rain(:,:,:) = 2.0*rain_reff
        Lsc_microphys%size_ice (:,:,:) = ice_cld_reff
      endif ! (do_diag_clouds)

!--------------------------------------------------------------------


end subroutine microphys_presc_conc


!#################################################################



                       end module cloud_spec_mod


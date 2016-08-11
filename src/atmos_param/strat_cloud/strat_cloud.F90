!FDOC_TAG_GFDL
module strat_cloud_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Stephen Klein
! </CONTACT>
!/"/>q
! <OVERVIEW>
!   Code to compute time tendencies of stratiform clouds and diagnose
!   rain and snow flux with prognostic scheme.
!   
! </OVERVIEW>
! <DESCRIPTION>
!
!
!       The prognostic scheme returns the time tendencies of liquid,
!       ice, and saturated volume fraction that are suspended in 
!       stratiform clouds.  The scheme diagnoses the fluxes of rain
!       and snow in saturated and unsaturated areas.
!
!       The prognostic cloud scheme is responsible for determing
!       cloud volume fractions, condensed water tendencies, and
!       the stratiform precipitation rate.  It includes processes
!       for evaporation, condensation, deposition, and sublimation
!       of cloud water, conversion of cloud water to precipitation,
!       evaporation of falling precipitation, the bergeron-findeisan 
!       process, freezing of cloud liquid, accretion of cloud water 
!       by precipitation, and melting of falling precipitation.
!
!       This scheme is based on the experience the author had 
!       at the ECMWF in 1997. The saturated volume fraction formalism 
!       and type of solution follow directly from the scheme of Tiedtke
!       (1993): Monthly Weather Review, Volume 121, pages 3040-3061.
!       The form of most of the microphysics follows Rotstayn , 1997:
!       Quart. J. Roy. Met. Soc. vol 123, pages 1227-1282. The partial
!       precipitation area formulism follows Jakob and Klein, 2000:
!       Quart. J. Roy. Met. Soc. vol 126, pages 2525-2544. 
!
!       The statistical cloud scheme treatment, which is used as
!       a replacement for the Tiedtke cloud fraction scheme, is based
!       on a number of publications: Tompkins, A., 2002: J. Atmos. 
!       Sci., 59, 1917-1942, Klein et al., 2005: J. Geophys. Res., 
!       110, D15S06, doi:10.1029/2004JD005017. 
! </DESCRIPTION>
!

! <DATASET NAME="strat_cloud.res">
!   native format of the restart file
! </DATASET>
! <DATASET NAME="strat_cloud.res.nc">
!   netcdf format of the restart file
! </DATASET>


! <INFO>
!   <REFERENCE>           
!The saturation volume fraction formalism comes from:
!
!Tiedtke, M., 1993: Representation of clouds in large-scale models.
! Mon. Wea. Rev., 121, 3040-3061.
!
! </REFERENCE>
!   <REFERENCE>           
!The form of most of the microphysics follows:
!
!Rotstayn, L., 1997: A physically based scheme for the treatment of
! stratiform clouds and precipitation in large-scale models. I:
! Description and evaluation of microphysical processes. Quart. J.
! Roy. Met. Soc. 123, 1227-1282. 
! </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!  1. qmin should be chosen such that the range of {qmin, max(qa,ql
!     ,qi)} is resolved by the precision of the numbers used. (default =
!     1.E-10)
!   </NOTE>

!   <NOTE> 
!  2. Dmin will be MACHINE DEPENDENT and occur when
!   </NOTE>

!   <NOTE> 
!    a. 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small
!       Dmin
!   </NOTE>

!   AND

!   <NOTE> 
!    b. 1. - exp(-D) < D for all D > Dmin
!   </NOTE>
!   <FUTURE>               </FUTURE>

! </INFO>

use mpp_mod,                   only : mpp_clock_id, mpp_clock_begin, &
                                      mpp_clock_end, CLOCK_LOOP,  &
                                      input_nml_file
use fms_mod,                   only : file_exist, open_namelist_file,&
                                      open_file, error_mesg, FATAL, NOTE, &
                                      mpp_pe, mpp_root_pe, close_file, &
                                      read_data, write_data, stdlog,  &
                                      open_restart_file, open_ieee32_file,&
                                      mpp_error, check_nml_error,  &
                                      write_version_number, stdout
use fms_io_mod,                only : get_restart_io_mode, restore_state, &
                                      register_restart_field,  &
                                      restart_file_type, save_restart, &
                                      get_mosaic_tile_file
use constants_mod,             only : rdgas, rvgas, hlv, hls, hlf, cp_air, grav
use cloud_rad_mod,             only : cloud_rad_init
use diag_manager_mod,          only : register_diag_field, send_data
use time_manager_mod,          only : time_type, get_date, get_time
use rad_utilities_mod,         only : aerosol_type
use microphysics_mod,          only : microphysics_init, microphysics, &
                                      microphysics_end
use nc_cond_mod,               only:  nc_cond, nc_cond_init,  nc_cond_end
use polysvp_mod,               only : polysvp_init, polysvp_end, &
                                      compute_qs_a
use aerosol_cloud_mod,         only:  aerosol_cloud, aerosol_cloud_init, &
                                      aerosol_cloud_end
use strat_cloud_utilities_mod, only:  strat_cloud_utilities_init, &
                                      diag_id_type, diag_pt_type, &
                                      strat_nml_type, atmos_state_type, &
                                      particles_type, cloud_state_type, &
                                      strat_constants_type, &
                                      cloud_processes_type,  &
                                      precip_state_type
use strat_netcdf_mod,          only:  strat_netcdf, strat_netcdf_init, &
                                      strat_netcdf_end
use strat_cloud_legacy_mod,    only : strat_cloud_legacy_init,  &
                                      strat_cloud_legacy,   &
                                      strat_cloud_legacy_end
use check_nan_mod,             only : check_nan_init, check_nan


implicit none
private

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

!       The module contains the following  public interfaces:
!
!       strat_cloud_init       read namelist file, open logfile, initialize
!                              constants and fields, read restart file
!
!       strat_cloud            calls the legacy strat_cloud routine to 
!                              perform the calculations of the cloud scheme
!                              (strat_cloud_legacy_mod)
!
!       strat_cloud_new        rewritten version of strat_cloud  supporting
!                              double moment microphysics (rewritten by
!                              M. Salzmann and R. Hemler)
!
!       strat_cloud_end        writes out restart data to a restart file.
!
!       strat_cloud_sum        sum cloud scheme variables over time
!
!       strat_cloud_avg        return time average of summed cloud scheme 
!                              variables
!
!       do_strat_cloud         logical function, is the scheme on?
!
!       strat_cloud_restart    writes  module restart files
!
!-----------------------------------------------------------------------
!
!       The module contains the following  private subroutines:
!
!        fill_nml_variable     put all nml variables into a strat_nml_type
!                              variable for use by other modules
!
!        strat_debug           debug routine for additional output
!
!        inpose_realizability  make sure input cloud tracer fields are 
!                              physically realizable  
!
!        strat_alloc           allocate components of derived types used 
!                              in module
!
!        strat_dealloc         deallocate components of derived types used 
!                              in module
!
!------------------------------------------------------------------------

public  strat_cloud_init, strat_cloud, strat_cloud_new, strat_cloud_end,  &
        strat_cloud_sum, strat_cloud_avg, do_strat_cloud,  &
        strat_cloud_restart, strat_cloud_time_vary
private fill_nml_variable, strat_debug, impose_realizability, strat_alloc,&
        strat_dealloc


!------------------------------------------------------------------------
!---version number-------------------------------------------------------

Character(len=128) :: Version = '$Id: strat_cloud.F90,v 20.0 2013/12/13 23:22:09 fms Exp $'
Character(len=128) :: Tagname = '$Name: tikal $'

!------------------------------------------------------------------------
!---namelist-------------------------------------------------------------

# include "strat_nml.h"

!------------------------------------------------------------------------
!----module variables----------------------------------------------------

!-----------------------------------------------------------------------
!    accumulation arrays for fields that may be time averaged
!-----------------------------------------------------------------------
real,    allocatable, dimension (:,:,:) :: qlsum, qisum, cfsum
integer, allocatable, dimension (:,:)   :: nsum

!-----------------------------------------------------------------------
!    derived types present for duratioon of run.         
!-----------------------------------------------------------------------
type(strat_nml_type), save       :: Nml
type(strat_constants_type), save :: Constants

!------------------------------------------------------------------------
!     ------ constants used by the scheme -------
!-----------------------------------------------------------------------
real, parameter :: d608 = (rvgas - rdgas)/rdgas
  

!-----------------------------------------------------------------------
!    variables related to netcdf diagnostics.
!-----------------------------------------------------------------------
INTEGER            :: n_diag_4d, n_diag_4d_kp1
TYPE(diag_id_type) :: diag_id
TYPE(diag_pt_type) :: diag_pt
character(len=5)   :: mod_name = 'strat'
real               :: missing_value = -999.

!-----------------------------------------------------------------------
!    variables related to restart files.
!-----------------------------------------------------------------------
type(restart_file_type), pointer, save :: Str_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer, dimension(1)                  :: restart_versions = (/ 1 /)
integer                                :: vers

!-----------------------------------------------------------------------
!    variables related to debugging options.
!-----------------------------------------------------------------------
integer :: otun               ! file where debug output is written
logical :: debugo0 = .false.  ! small output
logical :: debugo1 = .false.  ! .true. !morr 
logical :: debugo3 = .false. 
logical :: debugo4 = .false.  ! when true, nrefuse will be output
integer :: ncall   = 1        ! timestep counter of calls to strat_cloud

!-----------------------------------------------------------------------
!    variables related to code timing.         
!-----------------------------------------------------------------------
integer   :: sc_loop, sc_pre_loop, sc_post_loop
integer   :: sc_micro, sc_init, sc_qs, sc_realiz, sc_aero, &
             sc_nccond, sc_after, sc_end

!-----------------------------------------------------------------------
!    miscellaneous variables 
!-----------------------------------------------------------------------
logical        :: do_average = .false. ! time average stratiform cloud 
                                       ! fields before computing radiation?
logical,PUBLIC :: strat_cloud_on = .false.
                                       ! is the stratiform cloud scheme
                                       ! operating? this variable used by 
                                       ! vert_turb_driver_mod; should
                                       ! change code to use function
                                       ! do_strat_cloud()
integer        :: overlap =  2         ! value of the overlap parameter
                                       ! obtained from cloud_rad_mod.
                                       ! overlap = 1 is maximum-random
                                       ! overlap = 2 is random
                                       ! RSH:specification method should be
                                       ! changed
logical :: module_is_initialized = .false.
                                       ! current module is initialized ?
logical :: running_old_code            ! we are running the legacy
                                       ! strat_cloud code (as opposed to
                                       ! the newer version compatible with
                                       ! double-moment microphysics)


CONTAINS



!#######################################################################

! <SUBROUTINE NAME="strat_cloud_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!       Initializes strat_cloud.  Reads namelist, calls cloud_rad_init,
!       reads restart (if present), initializes netcdf output.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_cloud_init(axes,Time,idim,jdim,kdim)
!                
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer">
!       Axes integer vector used for netcdf initialization.
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!       Time type variable used for netcdf.
!  </IN>
!  <IN NAME="idim" TYPE="integer">
!       Size of first array (usually longitude) dimension.
!  </IN>
!  <IN NAME="jdim" TYPE="integer">
!       Size of second array (usually latitude) dimension.
!  </IN>
!  <IN NAME="kdim" TYPE="integer">
!       Size of vertical array (usually height) dimension.
!  </IN>
! </SUBROUTINE>
!
!------------------------------------------------------------------------
subroutine strat_cloud_init (axes, Time, idim, jdim, kdim,  &
                                                     do_legacy_strat_cloud)

!------------------------------------------------------------------------
!    this subroutine reads the namelist file, opens a logfile, 
!    and initializes the physical constants of the routine.
!------------------------------------------------------------------------

integer,          intent(in)            :: idim, jdim, kdim, axes(4)
type(time_type),  intent(in)            :: Time
logical,          intent(in), optional  :: do_legacy_strat_cloud

!-----------------------------------------------------------------------
!
!    axes           integers corresponding to the
!                   x,y,z,z_half axes types
!
!    Time           time type variable
!
!    idim,jdim      number of points in first 
!                   and second dimensions
!
!    kdim           number of points in vertical
!                   dimension
!
!    do_legacy_strat_cloud
!                   activate the older version of the strat_cloud module
!                   rather than the latest ?
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!---local variables------------------------------------------------------

      integer                          :: id_restart
      integer                          :: unit,io,ierr, logunit
      integer                          :: vers2
      character(len=4)                 :: chvers
      character(len=64)                :: restart_file, fname
      character(len=64)                :: otname =  'debug_strt_cld.txt'

!-----------------------------------------------------------------------
!
!       unit           unit number for namelist and restart file
!
!       io             internal variable for reading of namelist file
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------ 
!    if module is already initialized, return.
!------------------------------------------------------------------------ 
      if (module_is_initialized) return

!----------------------------------------------------------------------- 
!    process namelist.
!-----------------------------------------------------------------------  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=strat_cloud_nml, iostat=io)
      ierr = check_nml_error(io,'strat_cloud_nml')
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=strat_cloud_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'strat_cloud_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!----------------------------------------------------------------------- 
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------  
      call write_version_number(version, tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() )  &
                    write (logunit, nml=strat_cloud_nml)

!----------------------------------------------------------------------- 
!    obtain info on form of restart file to be written.
!-----------------------------------------------------------------------  
      call get_restart_io_mode (do_netcdf_restart)

!-----------------------------------------------------------------------
!    check for acceptable namelist values:
!    qthalfwidth must be greater than 0.001
!    nsublevels must be greater than 0
!-----------------------------------------------------------------------
      if (qthalfwidth .lt. 1.e-03) then
        call error_mesg ( 'strat_cloud_mod', &
                       'qthalfwidth must be greater than 0.001', FATAL)
      endif
      if (nsublevels .lt. 1) then
        call error_mesg ( 'strat_cloud_mod', &
                            'nsublevels must be greater than 0', FATAL)
      endif

!------------------------------------------------------------------------
!    call subroutine to place nml variables in strat_nml_type variable.
!------------------------------------------------------------------------
      call fill_nml_variable

!------------------------------------------------------------------------
!    define logicals defining microphysics scheme which is active.
!------------------------------------------------------------------------
      if (trim(microphys_scheme) =='rotstayn_klein') then
        Constants%do_rk_microphys = .true.
        Constants%do_mg_microphys = .false.
        Constants%do_mg_ncar_microphys = .false.
        Constants%do_predicted_ice_number = .false.
      else if (trim(microphys_scheme) == 'morrison_gettelman') then
        Constants%do_rk_microphys = .false.
        Constants%do_mg_microphys = .true.
        Constants%do_mg_ncar_microphys = .false.
        Constants%do_predicted_ice_number = .true.
      else if (trim(microphys_scheme) == 'mg_ncar') then
        Constants%do_rk_microphys = .false.
        Constants%do_mg_microphys = .false.
        Constants%do_mg_ncar_microphys = .true.
        Constants%do_predicted_ice_number = .true.
      else
        call error_mesg ('strat_cloud_init', &
         'invalid expression supplied for nml variable microphys_scheme', &
                                                                    FATAL)
      endif
   
!------------------------------------------------------------------------
!    define logicals defining macrophysics scheme which is active.
!------------------------------------------------------------------------
     if (trim(macrophys_scheme) == 'tiedtke') then
       Constants%tiedtke_macrophysics = .true.
     else
       Constants%tiedtke_macrophysics = .false.
     endif

!------------------------------------------------------------------------
!    define logicals defining aerosol activation scheme which is active.
!------------------------------------------------------------------------
     if (trim(aerosol_activation_scheme) == 'dqa') then
       Constants%dqa_activation = .true.
       Constants%total_activation = .false.
     else if (trim(aerosol_activation_scheme) == 'total') then
       Constants%dqa_activation = .false.
       Constants%total_activation = .true.
     else
       call error_mesg ('strat_cloud_init',   &
           'invalid value for aerosol_activation_scheme specified', FATAL)
     endif

!-----------------------------------------------------------------------
!    pass values of qmin, N_land, N_ocean and as needed, do_liq_num and 
!    do_mg_microphys to cloud_rad_mod for use there. retrieve the value of
!    overlap from cloud_rad_mod for use here. this is done to assure 
!    consistency between these values in these two modules. (this process
!    should be changed to avoid potential problems --RSH)
!-----------------------------------------------------------------------
      if (do_liq_num) then
        if (Constants%do_rk_microphys) then 
          call cloud_rad_init (axes, Time, qmin_in=qmin, N_land_in=N_land,&
                               N_ocean_in=N_ocean,  &
                               prog_droplet_in=do_liq_num,  &
                               overlap_out=overlap)
        else if (Constants%do_mg_microphys .or.  &
                                     Constants%do_mg_ncar_microphys) then
          call cloud_rad_init (axes, Time, qmin_in=qmin, N_land_in=N_land,&
                               N_ocean_in=N_ocean,  &
                               prog_droplet_in=do_liq_num,  &
                               overlap_out=overlap,  &
                               qcvar_in = Nml%qcvar, &
                          prog_ice_num_in=Constants%do_mg_microphys .or.&
                                            Constants%do_mg_ncar_microphys)
        endif
      else
        call cloud_rad_init (axes, Time, qmin_in=qmin, N_land_in=N_land,&
                             N_ocean_in=N_ocean, overlap_out=overlap)
      endif

!-------------------------------------------------------------------------
!    save the returned value of overlap in a derived type variable
!-------------------------------------------------------------------------
      Constants%overlap = overlap

!-----------------------------------------------------------------------
!    set logicals indicating that strat_cloud is active and whether the
!    legacy or new version is being executed.
!-----------------------------------------------------------------------
      strat_cloud_on = .TRUE.
      if (present(do_legacy_strat_cloud)) then
        if (do_legacy_strat_cloud) then
          running_old_code = .true.
        else
          running_old_code = .false.
        endif
      else
        running_old_code = .true.
      endif

!-----------------------------------------------------------------------
!    allocate the arrays needed to accumulate cloud fields which may need
!    to be time-averaged. these fields are also saved to the restart file.
!-----------------------------------------------------------------------
      allocate (nsum(idim, jdim),      &
                qlsum(idim,jdim,kdim), &
                qisum(idim,jdim,kdim), &
                cfsum(idim,jdim,kdim)  )

!------------------------------------------------------------------------
!    register the restart fields to be written and/or read.
!------------------------------------------------------------------------- 
      if (do_netcdf_restart) then
        restart_file = 'strat_cloud.res.nc'
        call get_mosaic_tile_file (restart_file, fname, .false. ) 
        allocate(Str_restart)
        if (trim(restart_file) == trim(fname)) then
          Til_restart => Str_restart
          in_different_file = .false.
        else
          in_different_file = .true.
          allocate(Til_restart)
        endif
        id_restart = register_restart_field  &
             (Str_restart, restart_file, 'vers', vers, no_domain = .true.)
        id_restart = register_restart_field  &
                             (Til_restart, restart_file, 'nsum', nsum)
        id_restart = register_restart_field  &
                             (Til_restart, restart_file, 'qlsum', qlsum)
        id_restart = register_restart_field  &
                             (Til_restart, restart_file, 'qisum', qisum)
        id_restart = register_restart_field  &
                             (Til_restart, restart_file, 'cfsum', cfsum)
      endif

!-----------------------------------------------------------------------
!    see if restart file exists
!-----------------------------------------------------------------------
      if (file_exist('INPUT/strat_cloud.res.nc') ) then
        if (mpp_pe() == mpp_root_pe() )    &
                   call mpp_error ('strat_cloud_mod', &
             'Reading netCDF formatted restart file:  &
                                          &INPUT/strat_cloud.res.nc', NOTE)
!------------------------------------------------------------------------
!    make sure do_netcdf_restart is true.
!------------------------------------------------------------------------
        if (.not. do_netcdf_restart)   &
                     call error_mesg ('strat_cloud_mod', &
                         'netcdf format restart file &
                            &INPUT/strat_cloud.res.nc exist, but  &
                                    &do_netcdf_restart is false.', FATAL)
        call restore_state (Str_restart)
        if (in_different_file) call restore_state (Til_restart)
      else
        if (file_exist('INPUT/strat_cloud.res')) Then
          unit = open_restart_file (FILE='INPUT/strat_cloud.res', &
                                                            ACTION='read')
          if (mpp_pe() == mpp_root_pe() ) call mpp_error   &
                                 ('strat_cloud_mod', &
                           'Reading native formatted restart file.', NOTE)
          read (unit, iostat=io, err=142) vers, vers2
142       continue
          if (io == 0) then

!--------------------------------------------------------------------
!    if eor is not encountered, then the file includes radturbten.
!    that data is not needed, simply continue by reading next record.
!--------------------------------------------------------------------
            call error_mesg ('strat_cloud_mod',  &
                    'reading pre-version number strat_cloud.res file, &
                                              &ignoring radturbten', NOTE)

!--------------------------------------------------------------------
!    the file is a newer one with a version number included. read the 
!    version number. if it is not a valid version, stop execution with
!    a message.
!--------------------------------------------------------------------
          else
            if (.not. any(vers == restart_versions) ) then
              write (chvers, '(i4)') vers
              call error_mesg ('strat_cloud_mod',  &
                    'restart version ' // chvers//' cannot be read &
                             &by this version of strat_cloud_mod.', FATAL)
            endif
          endif
          call read_data (unit, nsum)
          call read_data (unit, qlsum)
          call read_data (unit, qisum)
          call read_data (unit, cfsum)
          call close_file (unit)
        else
          qlsum=0.0; qisum=0.0; cfsum=0.0; nsum=0
        endif
      endif
      vers = restart_versions(size(restart_versions(:)))

!-----------------------------------------------------------------------
!    call strat_netcdf_init to set up the netcdf diagnostic output
!    requested via the diag_table. 
!-----------------------------------------------------------------------
      call strat_netcdf_init (axes, Time, diag_id, diag_pt, n_diag_4d, &
                              n_diag_4d_kp1)

!-----------------------------------------------------------------------
!    initialize the other modules which are used by strat_cloud_mod.
!-----------------------------------------------------------------------
      call strat_cloud_utilities_init
      if (running_old_code) then
        call strat_cloud_legacy_init (do_pdf_clouds)
      else
        call microphysics_init (Nml)
        call aerosol_cloud_init (Constants)
        call nc_cond_init (do_pdf_clouds)
        call polysvp_init
        call check_nan_init
      endif

!-----------------------------------------------------------------------
!    set up clocks to time various code portions.
!-----------------------------------------------------------------------
      sc_pre_loop = mpp_clock_id ('strat_cloud: vertical loop setup', &
                                                        grain=CLOCK_LOOP)
      sc_loop = mpp_clock_id ('strat_cloud: main vertical level loop',&
                                                         grain=CLOCK_LOOP)
      sc_post_loop = mpp_clock_id ('strat_cloud: diagnostic send-data', &
                                                         grain=CLOCK_LOOP)
      sc_init = mpp_clock_id ('strat_cloud: start of routine', &
                                                          grain=CLOCK_LOOP)
      sc_qs = mpp_clock_id ('strat_cloud: qs call',  grain=CLOCK_LOOP)
      sc_realiz = mpp_clock_id ('strat_cloud: realizability',  &
                                                         grain=CLOCK_LOOP)
      sc_aero = mpp_clock_id ('strat_cloud: aerosol', grain=CLOCK_LOOP)
      sc_nccond = mpp_clock_id ('strat_cloud: nccond', grain=CLOCK_LOOP)
      sc_after = mpp_clock_id ('strat_cloud: after nccond',  &
                                                         grain=CLOCK_LOOP)
      sc_micro = mpp_clock_id ('strat_cloud: microphysics', &
                                                         grain=CLOCK_LOOP)
      sc_end = mpp_clock_id ('strat_cloud: end of routine',  &
                                                          grain=CLOCK_LOOP)

!-------------------------------------------------------------------------
!    if any debugging is desired, open a file to hold the output. set
!    debugging at level o0 to occur only on pe 0.
!-------------------------------------------------------------------------
      if ((debugo .or. debugo0 .or.  debugo1 .or. debugo3)) then
        otun = open_file (otname, threading = 'multi', action = 'append')
      else 
        otun = 0
      endif
      if ( mpp_pe() .NE. 0 ) then
        debugo0 = .FALSE.
      endif

!------------------------------------------------------------------------
!    mark the module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------


end subroutine strat_cloud_init



!#######################################################################

!-----------------------------------------------------------------------
! <SUBROUTINE NAME="strat_cloud">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!       
!  </DESCRIPTION>
!  <TEMPLATE>
! subroutine strat_cloud (Time,is,ie,js,je,dtcloud,pfull,phalf,radturbten2,
!                         T,qv,ql,qi,qa,omega,Mc,diff_t,LAND,              
!                         ST,SQ,SL,SI,SA,rain3d,snow3d,snowclr3d,surfrain, 
!                         surfsnow,qrat,ahuco,limit_conv_cloud_frac,MASK, 
!                         qn, Aerosol, SN)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!         Time
!  </IN>
!  <IN NAME="is" TYPE="integer">
!         Index of starting point in the longitude direction of the current
!         physics window               
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!         Index of ending point in the longitude direction of the current  
!         physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!         Index of starting point in the latitude direction of the current
!         physics window
!  </IN>
!  <IN NAME="je" TYPE="integer">
!         Index of ending point in the latitude direction of the current
!         physics window
!  </IN>
!  <IN NAME="dtcloud" TYPE="real">
!         Physics time step (sec)
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!         Pressure on model full levels (Pa)
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!         Pressure on model half levels (Pa)
!  </IN>
!  <IN NAME="radturbten2" TYPE="real">
!         Sum of the tendencies of temperature from turbulence and 
!         radiation schemes (K/s)
!  </IN>
!  <IN NAME="T" TYPE="real">
!         Temperature (K)         
!  </IN>
!  <IN NAME="qv" TYPE="real">
!         Water vapor specific humidity (kg vapor/kg air)
!  </IN>
!  <IN NAME="ql" TYPE="real">
!         Grid-box mean liquid water specific humidity (kg liquid/kg air)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!         Grid-box mean ice water specific humidity (kg ice/kg air)
!  </IN>
!  <IN NAME="qa" TYPE="real">
!         Cloud fraction (3d array and a prognostic variable) (fraction)
!  </IN>
!  <IN NAME="qn" TYPE="real">
!         Cloud droplet number (3d array and a prognostic variable) 
!         (#/kg air)
!  </IN>
!  <IN NAME="omega" TYPE="real">
!         Vertical pressure velocity (Pa/sec)
!  </IN>
!  <IN NAME="Mc" TYPE="real">
!         Cumulus mass flux (defined positive as upward) (kg air/m2/sec)
!  </IN>
!  <IN NAME="diff_t" TYPE="real">
!         Vertical diffusion coefficient for temperature and tracer from 
!         vertical diffusion scheme (m2/sec) 
!  </IN>
!  <IN NAME="LAND" TYPE="real">
!         Fraction of surface that contains land (fraction)
!  </IN>
!  <OUT NAME="ST" TYPE="real">
!         Change in temperature due to strat_cloud (K) 
!  </OUT>
!  <OUT NAME="SQ" TYPE="real">
!         Change in water vapor due to strat_cloud (kg vapor/kg air) 
!  </OUT>
!  <OUT NAME="SL" TYPE="real">
!         Change in cloud liquid due to strat_cloud (kg liquid/kg air)
!  </OUT>
!  <OUT NAME="SI" TYPE="real">
!         Change in cloud ice due to strat_cloud (kg ice/kg air)
!  </OUT>
!  <OUT NAME="SA" TYPE="real">
!         Change in cloud fraction due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="SN" TYPE="real">
!         Change in cloud droplet number due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="surfrain" TYPE="real">
!         Surface rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="surfsnow" TYPE="real">
!         Surface snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <OUT NAME="rain3d" TYPE="real">
!         3D rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="snow3d" TYPE="real">
!         3D snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <IN NAME="qrat" TYPE="real">
!         Ratio of large-scale specific humidity to specific humidity in 
!         environment outside of activated convective systems 
!         (donner_deep, uw) 
!  </IN>
!  <IN NAME="ahuco" TYPE="real">
!         The fraction of the grid box containing either cumulus cells or 
!         the  mesoscale circulation from donner_deep, and any uw shallow 
!         clouds.
!  </IN>
!  <IN NAME="MASK" TYPE="real">
!         Optional input real array indicating the point is above the 
!         surface if equal to 1.0 and indicating the point is below the 
!         surface if equal to 0.
!         Used only in eta vertical coordinate model.
!  </IN>
! </SUBROUTINE>
!
!-----------------------------------------------------------------------

subroutine strat_cloud    &
         (Time, is, ie, js, je, dtcloud, pfull, phalf, radturbten2,&
          T, qv, ql, qi ,qa, omega, Mc, diff_t, LAND,              &
          ST, SQ, SL, SI, SA, f_snow_berg, rain3d, snow3d, snowclr3d,   &
          surfrain, surfsnow, qrat, ahuco, limit_conv_cloud_frac, MASK,  &
          qn, Aerosol, SN)

!-------------------------------------------------------------------------
type(time_type),        intent (in)            :: Time
integer,                intent (in)            :: is, ie, js, je
real,                   intent (in)            :: dtcloud
real, dimension(:,:,:), intent (in)            :: pfull, phalf, T, qv,  &
                                                  ql, qi, qa, omega, Mc, &
                                                  diff_t,  qrat, ahuco, &
                                                  radturbten2
logical, intent(in)                            :: limit_conv_cloud_frac
real, dimension(:,:),   intent (in)            :: LAND
real, dimension(:,:,:), intent (out)           :: ST, SQ, SL, SI, SA,  &
                                                  rain3d, snow3d,  &
                                                  snowclr3d, f_snow_berg
real, dimension(:,:),   intent (out)           :: surfrain, surfsnow
real, dimension(:,:,:), intent (in),  optional :: MASK, qn
type(aerosol_type),     intent (in),  optional :: Aerosol  
real, dimension(:,:,:), intent (out), optional :: SN
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!---local variables------------------------------------------------------

      real, allocatable, dimension(:,:,:,:) :: diag_4d, diag_4d_kp1    
      real, allocatable, dimension(:,:,:)   :: diag_3d
      integer                               :: kdim

!-----------------------------------------------------------------------
      kdim = size (T,3)

!------------------------------------------------------------------------
!    allocate arrays to hold netcdf diagnostics.
!------------------------------------------------------------------------
      if (allocated(diag_3d)) deallocate (diag_3d)
      allocate(diag_3d(size(T,1),size(T,2),0:n_diag_4d))
 
      diag_3d(:,:,0:) = 0.

      if (allocated(diag_4d)) deallocate (diag_4d)
      allocate(diag_4d(size(T,1),size(T,2),size(T,3),0:n_diag_4d))

      diag_4d(:,:,:,0:) = 0.

      if (allocated(diag_4d_kp1)) deallocate (diag_4d_kp1)
      allocate(  &
            diag_4d_kp1(size(T,1),size(T,2),size(T,3)+1,0:n_diag_4d_kp1))

      diag_4d_kp1(:,:,:,0:) = 0.

!-----------------------------------------------------------------------
!    call strat_cloud_legacy to execute the pre-double-moment-capable
!    strat_cloud code.
!------------------------------------------------------------------------
      call strat_cloud_legacy       &
            (Nml, diag_id, diag_pt, n_diag_4d, n_diag_4d_kp1, diag_4d, &
             diag_4d_kp1, diag_3d, Time, is, &
             ie, js, je, dtcloud, pfull, phalf, radturbten2, T, qv, ql, &
             qi, qa, omega, Mc, diff_t, LAND, ST, SQ, SL, SI, SA,  &
             f_snow_berg, rain3d, &
             snow3d, snowclr3d, surfrain, surfsnow, qrat, ahuco,   &
             limit_conv_cloud_frac, MASK, qn, Aerosol, SN)
        
!-----------------------------------------------------------------------
!    call strat_netcdf to process diagnostics.
!-----------------------------------------------------------------------
      call strat_netcdf (diag_id, diag_pt, diag_4d, diag_4d_kp1,  &
                         diag_3d, Time, is, js, kdim, MASK)
 
!-------------------------------------------------------------------------
!    deallocate the local allocatable variables. call strat_dealloc to 
!    deallocate the derived type variable components.
!-------------------------------------------------------------------------
      deallocate ( diag_4d )
      deallocate ( diag_4d_kp1 )
      deallocate ( diag_3d )

!---------------------------------------------------------------------


end subroutine strat_cloud


subroutine strat_cloud_time_vary (dtcloud, limit_conv_cloud_frac)

real, intent(in) :: dtcloud
logical, intent(in) :: limit_conv_cloud_frac

      Constants%dtcloud = dtcloud
      Constants%inv_dtcloud = 1.0/dtcloud
      Constants%limit_conv_cloud_frac = limit_conv_cloud_frac

end subroutine strat_cloud_time_vary

!#########################################################################

!------------------------------------------------------------------------
! <SUBROUTINE NAME="strat_cloud_new">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!       
!  </DESCRIPTION>
!  <TEMPLATE>
!subroutine strat_cloud_new (Time, is, ie, js, je, dtcloud, pfull,  
!                            phalf, zhalf, zfull, radturbten2, 
!                            T_in, qv_in, ql_in, qi_in, qa_in, omega, Mc, 
!                            diff_t, LAND, ST_out, SQ_out, SL_out, SI_out, 
!                            SA_out, rain3d, snow3d, snowclr3d, surfrain, 
!                            surfsnow, qrat, ahuco, limit_conv_cloud_frac,
!                            Aerosol, MASK3d, qn_in, SN_out, qni_in,  
!                            SNi_out, lsc_snow, lsc_rain, lsc_snow_size, 
!                            lsc_rain_size )
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!         Time
!  </IN>
!  <IN NAME="is" TYPE="integer">
!         Index of starting point in the longitude direction of the current
!         physics window               
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!         Index of ending point in the longitude direction of the current  
!         physics window
!  </IN>
!  <IN NAME="js" TYPE="integer">
!         Index of starting point in the latitude direction of the current
!         physics window
!  </IN>
!  <IN NAME="je" TYPE="integer">
!         Index of ending point in the latitude direction of the current
!         physics window
!  </IN>
!  <IN NAME="dtcloud" TYPE="real">
!         Physics time step (sec)
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!         Pressure on model full levels (Pa)
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!         Pressure on model half levels (Pa)
!  </IN>
!  <IN NAME="zfull" TYPE="real">
!         Height on model full levels (m)
!  </IN>
!  <IN NAME="zhalf" TYPE="real">
!         Height on model half levels (m)
!  </IN>
!  <IN NAME="radturbten2" TYPE="real">
!         Sum of the tendencies of temperature from turbulence and 
!         radiation schemes (K/s)
!  </IN>
!  <IN NAME="T_in" TYPE="real">
!         Temperature (K)         
!  </IN>
!  <IN NAME="qv_in" TYPE="real">
!         Water vapor specific humidity (kg vapor/kg air)
!  </IN>
!  <IN NAME="ql_in" TYPE="real">
!         Grid-box mean liquid water specific humidity (kg liquid/kg air)
!  </IN>
!  <IN NAME="qi_in" TYPE="real">
!         Grid-box mean ice water specific humidity (kg ice/kg air)
!  </IN>
!  <IN NAME="qa_in" TYPE="real">
!         Cloud fraction (3d array and a prognostic variable) (fraction)
!  </IN>
!  <IN NAME="omega" TYPE="real">
!         Vertical pressure velocity (Pa/sec)
!  </IN>
!  <IN NAME="Mc" TYPE="real">
!         Cumulus mass flux (defined positive as upward) (kg air/m2/sec)
!  </IN>
!  <IN NAME="diff_t" TYPE="real">
!         Vertical diffusion coefficient for temperature and tracer from 
!         vertical diffusion scheme (m2/sec) 
!  </IN>
!  <IN NAME="LAND" TYPE="real">
!         Fraction of surface that contains land (fraction)
!  </IN>
!  <OUT NAME="ST_out" TYPE="real">
!         Change in temperature due to strat_cloud_new (K) 
!  </OUT>
!  <OUT NAME="SQ_out" TYPE="real">
!         Change in water vapor due to strat_cloud_new (kg vapor/kg air) 
!  </OUT>
!  <OUT NAME="SL_out" TYPE="real">
!         Change in cloud liquid due to strat_cloud_new (kg liquid/kg air)
!  </OUT>
!  <OUT NAME="SI_out" TYPE="real">
!         Change in cloud ice due to strat_cloud_new (kg ice/kg air)
!  </OUT>
!  <OUT NAME="SA_out" TYPE="real">
!         Change in cloud fraction due to strat_cloud_new (fraction)
!  </OUT>
!  <OUT NAME="rain3d" TYPE="real">
!         3D rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="snow3d" TYPE="real">
!         3D snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <OUT NAME="snowclr3d" TYPE="real">
!         3D snow fall outside of clouds over time step dtcloud (kg ice/m2)
!  </OUT>
!  <OUT NAME="surfrain" TYPE="real">
!         Surface rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="surfsnow" TYPE="real">
!         Surface snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <IN NAME="qrat" TYPE="real">
!         Ratio of large-scale specific humidity to specific humidity in 
!         environment outside convective systems (donner_deep, uw) 
!  </IN>
!  <IN NAME="ahuco" TYPE="real">
!         The fraction of the grid box containing either cumulus cells or 
!         the  mesoscale circulation from donner_deep, and any uw shallow 
!         clouds.
!  </IN>
!  <IN NAME="limit_conv_cloud_frac" TYPE="logical">
!         Total cloud fraction in a grid box should be limited to 100% ?
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!         Model aerosol amounts needed for nuclei activation calculations.
!  </IN>
!  <IN NAME="MASK3d" TYPE="real">
!         Optional input real array indicating the point is above the 
!         surface if equal to 1.0 and indicating the point is below the 
!         surface if equal to 0. Used only in eta vertical coordinate 
!         model.
!  </IN>
!  <IN NAME="qn_in" TYPE="real">
!         Cloud droplet number (3d array and a prognostic variable) 
!         (#/kg air)
!  </IN>
!  <OUT NAME="SN_out" TYPE="real">
!         Change in cloud droplet number due to strat_cloud_new (fraction)
!  </OUT>
!  <IN NAME="qni_in" TYPE="real">
!         Ice particle number (3d array and a prognostic variable) 
!         (#/kg air)
!  </IN>
!  <OUT NAME="SNi_out" TYPE="real">
!         Change in ice particle number due to strat_cloud_new (fraction)
!  </OUT>
!  <OUT NAME="lsc_snow" TYPE="real">
!         3D snow field (mass mixing ratio)from morrison-gettelman 
!         calculation.   (kg/kg)
!  </OUT>
!  <OUT NAME="lsc_rain" TYPE="real">
!         3D rain field (mass mixing ratio) from morrison-gettelman 
!         calculation.   (kg/kg)
!  </OUT>
!  <OUT NAME="lsc_snow_size" TYPE="real">
!         Snow flake size from morrison-gettelman calculation.   (microns)
!  </OUT>
!  <OUT NAME="lsc_rain_size" TYPE="real">
!         Rain drop size from morrison-gettelman calculation.   (microns)
!  </OUT>
! </SUBROUTINE>
!
!-----------------------------------------------------------------------


subroutine strat_cloud_new (Time, is, ie, js, je, dtcloud, pfull,  &
                            phalf, zhalf, zfull, radturbten2, &
                            T_in, qv_in, ql_in, qi_in, qa_in, omega, Mc, &
                            diff_t, LAND, ST_out, SQ_out, SL_out, SI_out, &
                            SA_out, f_snow_berg, rain3d, snow3d,    &
                            snowclr3d, surfrain, &
                            surfsnow, qrat, ahuco, limit_conv_cloud_frac, &
                            Aerosol, MASK3d, qn_in, SN_out, qni_in,  &
                            SNi_out, lsc_snow, lsc_rain, lsc_snow_size, &
                            lsc_rain_size )

!------------------------------------------------------------------------
type(time_type),        intent (in)            :: Time
integer,                intent (in)            :: is,ie,js,je
real,                   intent (in)            :: dtcloud
real, dimension(:,:,:), intent (in)            :: pfull, phalf, zhalf, &
                                                  zfull, T_in, qv_in,  &
                                                  ql_in, qi_in, qa_in, &
                                                  omega,Mc, diff_t, &
                                                  qrat, ahuco, radturbten2
logical,                intent(in)             :: limit_conv_cloud_frac
type(aerosol_type),     intent (in)            :: Aerosol  
real, dimension(:,:),   intent (in)            :: LAND
real, dimension(:,:,:), intent (out)           :: ST_out, SQ_out, SL_out, &
                                                  SI_out, SA_out, rain3d, &
                                                  snow3d, snowclr3d
real, dimension(:,:,:), intent (out)           :: f_snow_berg         
real, dimension(:,:),   intent (out)           :: surfrain,surfsnow
real, dimension(:,:,:), intent (in),  optional :: MASK3d, qn_in, qni_in
real, dimension(:,:,:), intent (out), optional :: SN_out, SNi_out,  &
                                                  lsc_snow, lsc_rain, &
                                                  lsc_snow_size, &
                                                  lsc_rain_size


!------------------------------------------------------------------------
!---local variables------------------------------------------------------

!------------------------------------------------------------------------
!  variables used in calculation of particle number diagnostics:
      real, dimension(size(T_in,1),size(T_in,2))  ::    &
              N3D_col, N3Di_col, N3D_col250, gb_N3D_col, gb_N3Di_col
      real :: dum, qa_new, qn_new, ql_new, qi_new, qni_new

!------------------------------------------------------------------------
!  variables used to hold particle number fields:
!       N3D   number of cloud drops per unit volume in liquid clouds
!             [ 1/(m*m*m) ]
!       N3Di  number of ice particles per unit volume in ice clouds
!             [ 1/(m*m*m) ]
      real, dimension(size(T_in,1),size(T_in,2),size(T_in,3)) ::     &
              N3Di, N3D

!------------------------------------------------------------------------
!  variables allocated to hold all netcdf diagnostic fields:
      real, allocatable, dimension(:,:,:,:) :: diag_4d, diag_4d_kp1
      real, allocatable, dimension(:,:,:)   :: diag_3d

!------------------------------------------------------------------------
!  derived type variables used to pass fields between modules:
      type(atmos_state_type)     :: Atmos_state
      type(cloud_state_type)     :: Cloud_state
      type(particles_type)       :: Particles 
      type(precip_state_type)    :: Precip_state
      type(cloud_processes_type) :: Cloud_processes

!------------------------------------------------------------------------
!  variables needed with column diagnostics:
      integer  :: unit, ipt, jpt, outunit

!------------------------------------------------------------------------
!  dimensions of physics window:
      integer  :: idim, jdim, kdim

!------------------------------------------------------------------------
!  do-loop indices:
      integer :: i, j ,k, nn
 
!------------------------------------------------------------------------
!  counter of columns in which mg_micro is not computed due to negative
!  water in column (activated by setting debugo4 to .true.)
      integer :: nrefuse 

      outunit = stdout()

!-----------------------------------------------------------------------
!    check for consistent arguments and options.
!-----------------------------------------------------------------------
      call mpp_clock_begin (sc_init)
      if (.not. module_is_initialized) call error_mesg  &
            ('strat_cloud_new', 'strat_cloud_new is not initialized',FATAL)

      IF (Constants%do_rk_microphys ) THEN
        IF ( present(MASK3d) ) THEN
          Constants%mask_present = .true.
          Constants%mask = MASK3d
        ELSE
          Constants%mask_present = .false.
          Constants%mask = 1.0       
        END IF
      ELSE if (Constants%do_mg_microphys .or.    &
                                     Constants%do_mg_ncar_microphys) then
        IF ( .NOT. present(SNi_out)) THEN
          call error_mesg ('strat_cloud_new_mod', &
             'morrison gettelman microp requires progn. ice num ',  FATAL) 
        END IF
        IF ( .NOT. PRESENT ( lsc_snow ) ) THEN
          call error_mesg ( 'strat_cloud_new_mod', &
            'need lsc_snow for morrrison-gettelman microphysics', FATAL)
        END IF   
        IF ( .NOT. PRESENT ( lsc_rain ) ) THEN
          call error_mesg ( 'strat_cloud_new_mod', &
            'need lsc_rain for morrrison-gettelman microphysics', FATAL)
        END IF   
        IF ( .NOT. PRESENT ( lsc_snow_size ) ) THEN
          call error_mesg ( 'strat_cloud_new_mod', &
          'need lsc_snow_size for morrrison-gettelman microphysics', FATAL)
        END IF   
        IF ( .NOT. PRESENT ( lsc_rain_size ) ) THEN
          call error_mesg ( 'strat_cloud_new_mod', &
          'need lsc_rain_size for morrrison-gettelman microphysics', FATAL)
        END IF   
        IF (PRESENT(MASK3d)) then
          call error_mesg ( 'strat_cloud_new_mod',  &
              "ERROR: mask not implementd with m-g rescaling", FATAL)
        END IF
      ENDIF  ! (do_rk_microphys)
      if (debugo) then
        IF( MAXVAL( ahuco ) .GT. 1. ) WRITE(outunit,*) "AHUCO WARNING"
      endif

!-------------------------------------------------------------------------
!    define spatial dimensions.                
!-------------------------------------------------------------------------
      idim = SIZE(T_in, 1)
      jdim = SIZE(T_in, 2)
      kdim = SIZE(T_in, 3) 

!-------------------------------------------------------------------------
!    initialize debug / diagnostic variables. be sure isamp, jsamp, ksamp
!    are valid. output initial data to debug file (otun). if debugo4 is
!    activated, output call counter to stdout to track model progress 
!    during debug, and then increment the counter.
!-------------------------------------------------------------------------
      nrefuse = 0 
      isamp = min (idim, isamp)
      jsamp = min (jdim, jsamp)
      Nml%isamp = isamp
      Nml%jsamp = jsamp
      Nml%ksamp = ksamp
      if ( ncall .eq. 1 .and.   &
                ( debugo .or. debugo0 .or.  debugo1 .or. debugo3 ) ) then
        write(otun,*) "MODIF TEST"
        write(otun,*) "is,ie,js,je ", is, ie, js, je
        if (.not. debugo1 .and. .not. debugo3 )    &
          write(otun,*) "  isamp, jsamp, ksamp ", isamp, jsamp, ksamp
      end if

      if ( debugo4 .and. mpp_pe() .EQ. 0 ) then   
        write (outunit,*) "in stratcloud mod", ncall
      endif

      if (debugo .or. debugo0 .or. debugo3 )then
        write(otun,*) "-----------------------------------------------" 
        write(otun,*) "ncall ", ncall
      endif

      ncall = ncall +1

!------------------------------------------------------------------------
!    allocate and initialize the diagnostic variables.
!------------------------------------------------------------------------
      allocate (diag_3d(idim,jdim,0:n_diag_4d))
      allocate (diag_4d(idim,jdim,kdim,0:n_diag_4d))
      allocate (diag_4d_kp1(idim,jdim,kdim+1,0:n_diag_4d_kp1))
      diag_3d(:,:,0:) = 0.
      diag_4d(:,:,:,0:) = 0.
      diag_4d_kp1(:,:,:,0:) = 0.
       
!------------------------------------------------------------------------
!    call strat_alloc to allocate and initialize the derived type variables
!    of the module.
!------------------------------------------------------------------------
      call strat_alloc (idim, jdim, kdim, pfull, phalf, zhalf, zfull,&
                        radturbten2, T_in, qv_in, ql_in, qi_in, qa_in, &
                        omega, Mc, diff_t, qrat, ahuco,  &
                        Atmos_state, Particles, Cloud_state, &
                        Precip_state, Cloud_processes, &
                        qn_in=qn_in, qni_in=qni_in)

!-----------------------------------------------------------------------
!    initialize remaining output variables. define constants to be passed
!    to other modules.
!-----------------------------------------------------------------------
      ST_out = 0.
      SQ_out = 0.
      call mpp_clock_end (sc_init)

!-----------------------------------------------------------------------
!    calculate saturation specific humidity and its temperature 
!    derivative, thermal conductivity plus vapor diffusivity factor, 
!    and relative humidity.
!-----------------------------------------------------------------------
      call mpp_clock_begin (sc_qs)
      call compute_qs_a (idim, jdim, kdim, Nml, Atmos_state, Cloud_state) 
      if (debugo ) then
        write(otun, *) " T, pfull ", T_in(isamp,jsamp,ksamp),   &
                                              pfull(isamp,jsamp,ksamp)
        write(otun, *) " aaa dqsdT, qs ",   &
                                Atmos_state%dqsdT(isamp,jsamp,ksamp),  &
                                      Atmos_state%qs(isamp,jsamp,ksamp)
      endif
      call mpp_clock_end (sc_qs)

!-----------------------------------------------------------------------
!    define cloud droplet number for non-predicted-droplet-number case.
!    for the predicted case, N3D is defined in subroutine 
!    impose_realizability, called below.
!-----------------------------------------------------------------------
      call mpp_clock_begin (sc_realiz)
      if ( .not. do_liq_num) then
        do k=1,kdim
          N3D (:,:,k) = N_land*LAND(:,:) + N_ocean*(1. - LAND(:,:))
        end do
      endif 

!-----------------------------------------------------------------------
!    call impose_realizability to account for the fact that other processes
!    may have created negative tracer or extremely small values of tracer 
!    fields. the general reason for the extremely small values of the 
!    tracer fields is due to vertical diffusion, advection of condensate or!    cumulus induced subsidence (also a form of advection) of condensate.
!
!    in this step any values of the prognostic variables which are less 
!    than qmin are reset to zero, while conserving total moisture.
!
!    note that this is done slightly different for the Tiedtke cloud 
!    fraction than it is for pdf clouds. In the former, the filling 
!    requires that cloud liquid, cloud ice, and cloud fraction are greater
!    than qmin. For PDF clouds, cloud fraction need not be considered since
!    it is diagnosed from the PDF clouds.
!-----------------------------------------------------------------------
      call impose_realizability (idim, jdim, kdim, Atmos_state, &
                                 Cloud_state, SQ_out, ST_out, N3d,   &
                                 N3Di, otun, n_diag_4d, diag_4d)
      call mpp_clock_end (sc_realiz)

!------------------------------------------------------------------------
!    call aerosol_cloud to determine available condensation nuclei.
!------------------------------------------------------------------------
      call mpp_clock_begin (sc_aero)
      CALL aerosol_cloud (idim, jdim, kdim, n_diag_4d, Nml, Constants, &
                          Atmos_state, Particles, Cloud_state%qa_upd,  &
                          Aerosol, diag_4d, diag_id, diag_pt, otun)        
      call mpp_clock_end (sc_aero)

!------------------------------------------------------------------------
!    call nc_cond to calculate non-convective condensation using either
!    Tiedtke, Tompkins et al., or Klein simple PDF. this contains the 
!    dcond_ls calc., etc. from the Tiedtke scheme and the Klein pdf scheme.
!------------------------------------------------------------------------
      call mpp_clock_begin (sc_nccond)
      CALL nc_cond (idim, jdim, kdim, Nml, Constants, Atmos_state,   &
                    Cloud_state, ST_out, SQ_out, Cloud_processes,   &
                    Particles, n_diag_4d, diag_4d, diag_id, diag_pt, otun) 
      call mpp_clock_end (sc_nccond)
      
!------------------------------------------------------------------------
!    define the mean droplet number after this timestep's activation. note
!    that if the cloud area has not increased during the timestep in r-k 
!    microphysics, then droplet number does not increase.
!    delta_cf:  A_dt * (1.-qabar)   where A_dt = A*dt , A source rate
!    Eq. 7 of Yi's 2007 paper
!------------------------------------------------------------------------
      call mpp_clock_begin (sc_after)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (diag_id%potential_droplets > 0 .and.   &
                      Cloud_processes%da_ls(i,j,k) <= 0.0)   &
                           diag_4d(i,j,k,diag_pt%potential_droplets) = 0.0
            if (diag_id%subgrid_w_variance > 0 .and.   &
                      Cloud_processes%da_ls(i,j,k) <= 0.0)   &
                           diag_4d(i,j,k,diag_pt%subgrid_w_variance) = 0.0
            if (Cloud_processes%da_ls(i,j,k) > 0.0 .or.  &
                                  Constants%do_mg_ncar_microphys .or.  &
                                          Constants%do_mg_microphys) then
              Cloud_state%qn_mean(i,j,k) =   &
                    Cloud_state%qn_upd(i,j,k) +   &
                        max(Cloud_processes%delta_cf(i,j,k),0.)*  &
                           Particles%drop1(i,j,k)*1.e6/  &
                                               Atmos_state%airdens(i,j,k)
            else
              Particles%drop1(i,j,k) = 0.                
              Cloud_state%qn_mean(i,j,k) = Cloud_state%qn_upd(i,j,k)
            end if
          end do
        end do
      end do
      call mpp_clock_end (sc_after)

!-------------------------------------------------------------------------
!    call microphysics to compute those terms.
!-------------------------------------------------------------------------
      call mpp_clock_begin (sc_micro)
      call microphysics (idim, jdim, kdim, Nml, Constants, N3D,  &
                         Atmos_state, Cloud_state, Cloud_processes,  &
                         Particles, n_diag_4d, diag_4d, diag_id, diag_pt, &
                         n_diag_4d_kp1, diag_4d_kp1, ST_out, SQ_out,  &
                         Precip_state, otun, ncall, &
                         Cloud_state%qa_upd_0, Cloud_state%SA_0, &
                         nrefuse, isamp, jsamp, ksamp, debugo, debugo0, &
                         debugo1)    
      if (Constants%do_mg_ncar_microphys .or.  &
                                          Constants%do_mg_microphys) then
        if (diag_id%qadt_limits + diag_id%qa_limits_col > 0)    &
            diag_4d(:,:,:,diag_pt%qadt_limits) =    &
                   (Cloud_state%SA_out(:,:,:) - &
                            diag_4d(:,:,:,diag_pt%qadt_limits)) *  &
                                               Constants%inv_dtcloud
      endif
      call mpp_clock_end (sc_micro)

!-----------------------------------------------------------------------
!    place some output fields in the appropriate locations.
!-----------------------------------------------------------------------
      call mpp_clock_begin(sc_end)
      SL_out = Cloud_state%SL_out
      SI_out = Cloud_state%SI_out
      SA_out = Cloud_state%SA_out
      rain3d = Precip_state%rain3d
      snow3d = Precip_state%snow3d
!RSH
!  for r-k, snow in cloud is included in cloud ice. For mg and ncar,
!  snow is not included in cloud ice, so all snow must be put into
!  snowclr3d which for those schemes is used to hold the total 
!  precipitating ice field.
      if (Constants%do_mg_ncar_microphys .or.  &
                                          Constants%do_mg_microphys) then
        snowclr3d = Precip_state%snow3d
      else
        snowclr3d = Precip_state%snowclr3d
      endif
      surfrain = Precip_state%surfrain
      surfsnow = Precip_state%surfsnow
      f_snow_berg = Cloud_processes%f_snow_berg

      if (present(lsc_snow)) lsc_snow = Precip_state%lsc_snow
      if (present(lsc_rain)) lsc_rain = Precip_state%lsc_rain
      if (present(lsc_snow_size)) lsc_snow_size =   &
                                                Precip_state%lsc_snow_size
      if (present(lsc_rain_size)) lsc_rain_size =  &
                                                Precip_state%lsc_rain_size
    
      if (present(SN_out)) SN_out(:,:,:) = Cloud_state%SN_out(:,:,:)
      if (present(SNi_out)) SNi_out(:,:,:) = Cloud_state%SNI_out(:,:,:)

!-----------------------------------------------------------------------
!    define some diagnostics.
!-----------------------------------------------------------------------
      if (diag_id%SA3d + diag_id%SA2d > 0) then
        diag_4d(:,:,:,diag_pt%SA3d) = SA_out(:,:,:)*Constants%inv_dtcloud
      endif
      if (diag_id%ST3d + diag_id%ST2d > 0) then
        diag_4d(:,:,:,diag_pt%ST3d) = ST_out(:,:,:)*Constants%inv_dtcloud
      endif
      if (diag_id%SQ3d + diag_id%SQ2d > 0) then
        diag_4d(:,:,:,diag_pt%SQ3d) = SQ_out(:,:,:)*Constants%inv_dtcloud
      endif
      if (diag_id%SL3d + diag_id%SL2d > 0) then
        diag_4d(:,:,:,diag_pt%SL3d) = SL_out(:,:,:)*Constants%inv_dtcloud
      endif
      if (diag_id%SI3d + diag_id%SI2d > 0) then
        diag_4d(:,:,:,diag_pt%SI3d) = SI_out(:,:,:)*Constants%inv_dtcloud
      endif
      if (diag_id%SN3d + diag_id%SN2d > 0) then
        diag_4d(:,:,:,diag_pt%SN3d) = Cloud_state%SN_out(:,:,:)*   &
                                                    Constants%inv_dtcloud
      endif
      if (diag_id%SNI3d + diag_id%SNI2d > 0) then
        diag_4d(:,:,:,diag_pt%SNI3d) = Cloud_state%SNI_out(:,:,:)*   &
                                                    Constants%inv_dtcloud
      endif
        
!-----------------------------------------------------------------------
!    define some diagnostics.
!-----------------------------------------------------------------------
      if (diag_id%SA_imb + diag_id%SA_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SA_imb) =    &
             diag_4d(:,:,:,diag_pt%SA3d) -   (       &
                diag_4d(:,:,:,diag_pt%qadt_lsform)   &
             +  diag_4d(:,:,:,diag_pt%qadt_lsdiss)   &
             +  diag_4d(:,:,:,diag_pt%qadt_rhred)    &
             +  diag_4d(:,:,:,diag_pt%qadt_eros)     &
             +  diag_4d(:,:,:,diag_pt%qadt_fill)     &
             +  diag_4d(:,:,:,diag_pt%qadt_super)    &
             +  diag_4d(:,:,:,diag_pt%qadt_destr)    &
             +  diag_4d(:,:,:,diag_pt%qadt_limits)   &
             +  diag_4d(:,:,:,diag_pt%qadt_ahuco)    &
                                                          )
      endif
      if (diag_id%SL_imb + diag_id%SL_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SL_imb) =  &
             diag_4d(:,:,:,diag_pt%SL3d) -   (            &
                diag_4d(:,:,:,diag_pt%qldt_cond )         &
              + diag_4d(:,:,:,diag_pt%qldt_evap  )        &
              + diag_4d(:,:,:,diag_pt%qldt_eros  )        &
              + diag_4d(:,:,:,diag_pt%qldt_berg)          &
              + diag_4d(:,:,:,diag_pt%qldt_freez )        &
              + diag_4d(:,:,:,diag_pt%liq_adj    )        &
              + diag_4d(:,:,:,diag_pt%qldt_rime  )        &
              + diag_4d(:,:,:,diag_pt%qldt_accr)          &
              + diag_4d(:,:,:,diag_pt%qldt_auto)          &
              + diag_4d(:,:,:,diag_pt%qldt_fill  )        &
              + diag_4d(:,:,:,diag_pt%qldt_destr )        &
              + diag_4d(:,:,:,diag_pt%qldt_freez2)        &
              + diag_4d(:,:,:,diag_pt%qldt_sedi  )        &
              + diag_4d(:,:,:,diag_pt%qldt_accrs)         &
              + diag_4d(:,:,:,diag_pt%qldt_bergs)         &
              + diag_4d(:,:,:,diag_pt%qldt_HM_splinter)   &
              - diag_4d(:,:,:,diag_pt%qidt_melt2 )        &
              - diag_4d(:,:,:,diag_pt%qidt_accrs)         &
              - diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid) &    
                                                             )
      endif
      if (diag_id%SI_imb + diag_id%SI_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SI_imb) =     &
             diag_4d(:,:,:,diag_pt%SI3d) -   (          &
              - diag_4d(:,:,:,diag_pt%qldt_berg)        &
              - diag_4d(:,:,:,diag_pt%qldt_freez )      &
              - diag_4d(:,:,:,diag_pt%qldt_rime  )      &
              - diag_4d(:,:,:,diag_pt%qldt_freez2)      &
              - diag_4d(:,:,:,diag_pt%qldt_HM_splinter) &
              + diag_4d(:,:,:,diag_pt%qidt_dep  )       &
              + diag_4d(:,:,:,diag_pt%qidt_subl  )      &
              + diag_4d(:,:,:,diag_pt%qidt_fall  )      &
              + diag_4d(:,:,:,diag_pt%qidt_eros  )      &
              + diag_4d(:,:,:,diag_pt%qidt_melt  )      &
              + diag_4d(:,:,:,diag_pt%qidt_melt2 )      &
              + diag_4d(:,:,:,diag_pt%qidt_fill  )      &
              + diag_4d(:,:,:,diag_pt%qidt_destr )      &
              + diag_4d(:,:,:,diag_pt%qidt_qvdep )      &
              + diag_4d(:,:,:,diag_pt%qidt_auto)        &
              + diag_4d(:,:,:,diag_pt%qidt_accr)        &
              + diag_4d(:,:,:,diag_pt%qidt_accrs)       &
              + diag_4d(:,:,:,diag_pt%ice_adj    )      &
              - diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)  &
                                                          )
      endif


      if (diag_id%SN_imb + diag_id%SN_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SN_imb) =  &
             diag_4d(:,:,:,diag_pt%SN3d) -   ( &
                diag_4d(:,:,:,diag_pt%qndt_cond  )      &
              + diag_4d(:,:,:,diag_pt%qndt_evap  )      &
              + diag_4d(:,:,:,diag_pt%qndt_fill  )      &
              + diag_4d(:,:,:,diag_pt%qndt_berg  )      &
              + diag_4d(:,:,:,diag_pt%qndt_destr )      &
              + diag_4d(:,:,:,diag_pt%qndt_super )      &
              + diag_4d(:,:,:,diag_pt%qndt_freez )      &
              + diag_4d(:,:,:,diag_pt%qndt_sacws )      &
              + diag_4d(:,:,:,diag_pt%qndt_sacws_o)     &
              + diag_4d(:,:,:,diag_pt%qndt_eros  )      &
              + diag_4d(:,:,:,diag_pt%qndt_pra   )      &
              + diag_4d(:,:,:,diag_pt%qndt_auto  )      &
              + diag_4d(:,:,:,diag_pt%qndt_nucclim)     &
              + diag_4d(:,:,:,diag_pt%qndt_sedi )       &
              + diag_4d(:,:,:,diag_pt%qndt_melt)        &
              + diag_4d(:,:,:,diag_pt%qndt_ihom)        &
              + diag_4d(:,:,:,diag_pt%qndt_size_adj)    &
              + diag_4d(:,:,:,diag_pt%qndt_fill2)       &
              + diag_4d(:,:,:,diag_pt%qndt_contact_frz) &
              + diag_4d(:,:,:,diag_pt%qndt_cleanup)     &
              + diag_4d(:,:,:,diag_pt%qndt_cleanup2)    &
                                                            )
      endif

      if (diag_id%SNi_imb + diag_id%SNi_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SNi_imb) =     &
             diag_4d(:,:,:,diag_pt%SNi3d) -   ( &
                diag_4d(:,:,:,diag_pt%qnidt_fill )     &
              + diag_4d(:,:,:,diag_pt%qnidt_nnuccd)    &
              + diag_4d(:,:,:,diag_pt%qnidt_nsubi)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nerosi)    &
              + diag_4d(:,:,:,diag_pt%qnidt_nprci)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nprai)     &
              + diag_4d(:,:,:,diag_pt%qnidt_nucclim1)  &
              + diag_4d(:,:,:,diag_pt%qnidt_nucclim2)  &
              + diag_4d(:,:,:,diag_pt%qnidt_sedi  )    &
              + diag_4d(:,:,:,diag_pt%qnidt_melt  )    &
              + diag_4d(:,:,:,diag_pt%qnidt_size_adj ) &
              + diag_4d(:,:,:,diag_pt%qnidt_fill2  )   &
              + diag_4d(:,:,:,diag_pt%qnidt_super )    &
              + diag_4d(:,:,:,diag_pt%qnidt_ihom )     &
              + diag_4d(:,:,:,diag_pt%qnidt_destr )    &
              + diag_4d(:,:,:,diag_pt%qnidt_cleanup)   &
              + diag_4d(:,:,:,diag_pt%qnidt_cleanup2)  &
              + diag_4d(:,:,:,diag_pt%qnidt_nsacwi)    &
                                                           )
      endif

      if (diag_id%SQ_imb + diag_id%SQ_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%SQ_imb) =     &
             diag_4d(:,:,:,diag_pt%SQ3d) -   (                &
              - diag_4d(:,:,:,diag_pt%qldt_cond  )            &
              - diag_4d(:,:,:,diag_pt%qldt_evap  )            &
              - diag_4d(:,:,:,diag_pt%qldt_eros  )            &
              - diag_4d(:,:,:,diag_pt%liq_adj    )            &
              - diag_4d(:,:,:,diag_pt%qldt_fill  )            &
              - diag_4d(:,:,:,diag_pt%qldt_destr )            &
              - diag_4d(:,:,:,diag_pt%qidt_dep  )             &
              - diag_4d(:,:,:,diag_pt%qidt_subl  )            &
              - diag_4d(:,:,:,diag_pt%qidt_eros  )            &
              - diag_4d(:,:,:,diag_pt%qidt_fill  )            &
              - diag_4d(:,:,:,diag_pt%qidt_destr )            &
              - diag_4d(:,:,:,diag_pt%qidt_qvdep )            &
              - diag_4d(:,:,:,diag_pt%ice_adj    )            &
              + diag_4d(:,:,:,diag_pt%rain_evap  )            &
              + diag_4d(:,:,:,diag_pt%qdt_sedi_ice2vapor)     & 
              + diag_4d(:,:,:,diag_pt%qdt_sedi_liquid2vapor)  &  
              + diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)        &
              + diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid)     &    
              + diag_4d(:,:,:,diag_pt%qdt_snow_sublim  )      &
              + diag_4d(:,:,:,diag_pt%qdt_snow2vapor    )     &
                                                                )
      endif

      if (diag_id%ST_imb + diag_id%ST_imb_col > 0) then
        diag_4d(:,:,:,diag_pt%ST_imb) =     &
             diag_4d(:,:,:,diag_pt%ST3d) -    (              &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_berg)         &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_freez )       &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_rime  )       &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_freez2)       &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_accrs)        &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_bergs)        &
              - hlf*diag_4d(:,:,:,diag_pt%qldt_HM_splinter)  &
              + hlf*diag_4d(:,:,:,diag_pt%qidt_melt  )       &
              + hlf*diag_4d(:,:,:,diag_pt%qidt_melt2 )       &
              + hlf*diag_4d(:,:,:,diag_pt%qidt_accrs)        &
              + hlf*diag_4d(:,:,:,diag_pt%rain_freeze)       &
              - hlf*diag_4d(:,:,:,diag_pt%srfrain_accrs )    &
              - hlf*diag_4d(:,:,:,diag_pt%srfrain_freez )    &
              - hlf*diag_4d(:,:,:,diag_pt%snow_melt)         &
              
              + hlv*diag_4d(:,:,:,diag_pt%qldt_cond  )           &
              + hlv*diag_4d(:,:,:,diag_pt%qldt_evap  )           &
              + hlv*diag_4d(:,:,:,diag_pt%qldt_eros  )           &
              + hlv*diag_4d(:,:,:,diag_pt%liq_adj    )           &
              + hlv*diag_4d(:,:,:,diag_pt%qldt_fill  )           &
              + hlv*diag_4d(:,:,:,diag_pt%qldt_destr )           &
              - hlv*diag_4d(:,:,:,diag_pt%rain_evap  )           &
              - hlv*diag_4d(:,:,:,diag_pt%qdt_sedi_liquid2vapor) &  
              - hlv*diag_4d(:,:,:,diag_pt%qdt_cleanup_liquid)    &      

              + hls*diag_4d(:,:,:,diag_pt%qidt_dep  )          &
              + hls*diag_4d(:,:,:,diag_pt%qidt_subl  )         &
              + hls*diag_4d(:,:,:,diag_pt%qidt_eros  )         &
              + hls*diag_4d(:,:,:,diag_pt%qidt_fill  )         &
              + hls*diag_4d(:,:,:,diag_pt%qidt_destr )         &
              + hls*diag_4d(:,:,:,diag_pt%qidt_qvdep )         &
              + hls*diag_4d(:,:,:,diag_pt%ice_adj    )         &
              - hls*diag_4d(:,:,:,diag_pt%qdt_sedi_ice2vapor)  & 
              - hls*diag_4d(:,:,:,diag_pt%qdt_cleanup_ice)     &
              - hls*diag_4d(:,:,:,diag_pt%qdt_snow_sublim  )   &
              - hls*diag_4d(:,:,:,diag_pt%qdt_snow2vapor    )  &
                                                              )/cp_air 
              endif 

              if (diag_id%rain3d > 0) then 
                 diag_4d_kp1(:,:,:,diag_pt%rain3d) = Precip_state%rain3d(:,:,:) 
              endif 
              if (diag_id%snow3d > 0) then 
                 diag_4d_kp1(:,:,:,diag_pt%snow3d) = Precip_state%snow3d(:,:,:) 
              endif 

      if ( diag_id%cf_ice_init > 0 ) then  
        diag_4d(:,:,:,diag_pt%cf_ice_init) =   &
                                MIN(diag_4d(:,:,:,diag_pt%cf_ice_init), 1.)
      end if
      if (diag_id%droplets > 0) then
        diag_4d(:,:,:,diag_pt%droplets) = N3D(:,:,:)
      end if
      if (diag_id%droplets_wtd > 0) then
        diag_4d(:,:,:,diag_pt%droplets_wtd) = N3D(:,:,:)*ql_in(:,:,:)
      end if
      if (diag_id%ql_wt > 0) then
        diag_4d(:,:,:,diag_pt%ql_wt) = ql_in(:,:,:)
      end if
      if (diag_id%nice > 0) then
        diag_4d(:,:,:,diag_pt%nice) = N3Di(:,:,:)
      end if
      if (diag_id%qrout > 0) then
        diag_4d(:,:,:,diag_pt%qrout) = Precip_state%lsc_rain  (:,:,:)
      end if
      if (diag_id%qsout > 0) then
        diag_4d(:,:,:,diag_pt%qsout) = Precip_state%lsc_snow  (:,:,:) 
      end if

!-------------------------------------------------------------------------
!    call strat_debug to output data to file otun, if requested.
!-------------------------------------------------------------------------
      if (debugo) then
        call strat_debug (otun, ST_out, SQ_out, Cloud_state,  &
                                                 Precip_state, Atmos_State)
      endif

      
      if ( debugo4 .and. nrefuse .gt. 0)then
        write(outunit,*) "WARNING: Refusing to do two moment microphysics &
           &in columns containing points with negative total water: " ,  &
                                                                   nrefuse
      end if

!-----------------------------------------------------------------------
!   INSTANTANEOUS OUTPUT DIAGNOSTICS
!-----------------------------------------------------------------------
      if (num_strat_pts > 0) then
        do nn=1,num_strat_pts
          if (strat_pts(1,nn) >= is .and. strat_pts(1,nn) <= ie .and.  &
                  strat_pts(2,nn) >= js .and. strat_pts(2,nn) <= je) then
            ipt = strat_pts(1,nn); jpt = strat_pts(2,nn)
            i = ipt - is + 1; j = jpt - js + 1
            unit = open_ieee32_file ('strat.data', action='append')
            write (unit) ipt, jpt, ql_in(i,j,:) + SL_out(i,j,:)
            write (unit) ipt, jpt, qi_in(i,j,:) + SI_out(i,j,:)
            write (unit) ipt, jpt, qa_in(i,j,:) + SA_out(i,j,:)
            write (unit) ipt, jpt, T_in(i,j,:) + ST_out(i,j,:) 
            write (unit) ipt, jpt, qv_in(i,j,:) + SQ_out(i,j,:)
            write (unit) ipt, jpt, pfull(i,j,:)
            call close_file (unit)
          endif
        enddo
      endif

!------------------------------------------------------------------------
!    define column pressure-weighting factor.
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            Atmos_state%deltpg(i,j,k) =   &
                                  (phalf(i,j,k+1) - phalf(i,j,k))/grav
            if (present(MASK3d)) then
              Atmos_state%deltpg(i,j,k) =   &
                                Atmos_state%deltpg(i,j,k)*MASK3d(i,j,k)
            endif
          end do
        end do
      end do
           
!------------------------------------------------------------------------
!    generate column integrated diagnostics.
!------------------------------------------------------------------------
      do nn=1, n_diag_4d
        do k =kdim,1, -1
          diag_3d(:,:,nn) = diag_3d(:,:,nn) &
                          + diag_4d(:,:,k,nn)*Atmos_state%deltpg(:,:,k)
        enddo
      enddo
!------------------------------------------------------------------------
!SPECIAL CASES not following general pattern of above:
!rain, snow cloud ice and cloud liquid fallout -- only column integrals 
!are valid
!------------------------------------------------------------------------
      if (Constants%do_mg_microphys .or.   &
                                     Constants%do_mg_ncar_microphys) then
        if (diag_id%cld_liq_imb + diag_id%cld_liq_imb_col > 0) then
          diag_3d(:,:,diag_pt%cld_liq_imb) =   - ( &
                diag_3d(:,:,diag_pt%qldt_sedi )  +   &
                diag_3d(:,:,diag_pt%qdt_sedi_liquid2vapor)  ) 
        endif
        if (diag_id%cld_ice_imb + diag_id%cld_ice_imb_col > 0) then
          diag_3d(:,:,diag_pt%cld_ice_imb) =   - ( &
                diag_3d(:,:,diag_pt%qidt_fall )  +   &
                diag_3d(:,:,diag_pt%qdt_sedi_ice2vapor)  ) 
        endif
      endif
      IF ( diag_id%neg_rain > 0   ) &
          diag_3d(:,:,diag_pt%neg_rain) =   &
                          diag_4d(:,:,1,diag_pt%neg_rain) 
      IF ( diag_id%neg_snow > 0   ) &
          diag_3d(:,:,diag_pt%neg_snow) =   &
                          diag_4d(:,:,1,diag_pt%neg_snow) 
      if (diag_id%rain_imb + diag_id%rain_imb_col > 0) then
        diag_3d(:,:,diag_pt%rain_imb) =    &
             Precip_state%surfrain(:,:)*Constants%inv_dtcloud -   &
                diag_3d(:,:,diag_pt%cld_liq_imb)   +  (  &
                diag_3d(:,:,diag_pt%qldt_accr)      &
              + diag_3d(:,:,diag_pt%qldt_auto )     &
              + diag_3d(:,:,diag_pt%qidt_melt  )    &
              + diag_3d(:,:,diag_pt%rain_evap)      &
              + diag_3d(:,:,diag_pt%rain_freeze)    &
              - diag_3d(:,:,diag_pt%srfrain_accrs)  &
              - diag_3d(:,:,diag_pt%srfrain_freez)  &
              + diag_3d(:,:,diag_pt%neg_rain)       &
              - diag_3d(:,:,diag_pt%snow_melt )     &
              + diag_3d(:,:,diag_pt%qdt_snow2vapor) &
                                                         )
      endif
      if (diag_id%snow_imb + diag_id%snow_imb_col > 0) then

        if (Constants%do_rk_microphys) then
          diag_3d(:,:,diag_pt%snow_imb) =    &
             Precip_state%surfsnow(:,:)*Constants%inv_dtcloud -    &
               diag_3d(:,:,diag_pt%cld_ice_imb)  + (  &
                diag_3d(:,:,diag_pt%qldt_accrs)      &
              + diag_3d(:,:,diag_pt%qldt_bergs)      &
              + diag_3d(:,:,diag_pt%qidt_fall)       &
              + diag_3d(:,:,diag_pt%qidt_auto  )     &
              + diag_3d(:,:,diag_pt%qidt_accr  )     &
              - diag_3d(:,:,diag_pt%rain_freeze)     &
              + diag_3d(:,:,diag_pt%srfrain_accrs)   &
              + diag_3d(:,:,diag_pt%srfrain_freez)   &
              + diag_3d(:,:,diag_pt%snow_melt )      &
              + diag_3d(:,:,diag_pt%neg_snow)        &
              + diag_3d(:,:,diag_pt%qdt_snow_sublim) &
                                                         )
        else
          diag_3d(:,:,diag_pt%snow_imb) =    &
             Precip_state%surfsnow(:,:)*Constants%inv_dtcloud -    &
               diag_3d(:,:,diag_pt%cld_ice_imb)  + (  &
                diag_3d(:,:,diag_pt%qldt_accrs)       &
              + diag_3d(:,:,diag_pt%qldt_bergs)       &
              + diag_3d(:,:,diag_pt%qidt_auto  )      &
              + diag_3d(:,:,diag_pt%qidt_accr  )      &
              - diag_3d(:,:,diag_pt%rain_freeze)      &
              + diag_3d(:,:,diag_pt%srfrain_accrs)    &
              + diag_3d(:,:,diag_pt%srfrain_freez)    &
              + diag_3d(:,:,diag_pt%snow_melt )       &
              + diag_3d(:,:,diag_pt%neg_snow)         &
              + diag_3d(:,:,diag_pt%qdt_snow_sublim)  &
                                                        )
        endif
      endif
!------------------------------------------------------------------------
!SPECIAL CASES not following general pattern of above:
!yim: in-cloud droplet column burden
!------------------------------------------------------------------------
      IF ( diag_id%rain_mass_conv > 0   ) &
          diag_3d(:,:,diag_pt%rain_mass_conv) =   &
                          diag_4d(:,:,1,diag_pt%rain_mass_conv) 
      IF ( diag_id%snow_mass_conv > 0   ) &
          diag_3d(:,:,diag_pt%snow_mass_conv) =   &
                          diag_4d(:,:,1,diag_pt%snow_mass_conv) 

      if (diag_id%droplets_col > 0 .or. diag_id%gb_droplets_col > 0  .or. &
                                       diag_id%droplets_col250 > 0 ) then
        if (present (qn_in)) then
          N3D_col(:,:) = 0.
          N3D_col250(:,:) = 0.
          gb_N3D_col(:,:) = 0.
          do k =1,kdim
            do j=1,jdim
              do i=1,idim
! the current code:
                qa_new = qa_in(i,j,k) + SA_out(i,j,k)
                ql_new = ql_in(i,j,k) + SL_out(i,j,k)
                qn_new = qn_in(i,j,k) + SN_out(i,j,k)
                if (ql_new > qmin .and. &
                    qa_new > qmin .and. &
                    qn_new > qmin ) then      
!RSH 12/22/11 fix as per email from yim 11/3/11:
!                 dum = qn_new*Atmos_state%airdens(i,j,k)*  &
                  dum = qn_new*                             &
!RSH 12/22/11 fix as per email from yim 11/3/11:
!                                    Atmos_state%deltpg(i,j,k)*1.e-6
                                Atmos_state%deltpg(i,j,k)*1.e-4
                  if (qa_new > 0.05) then !count only columns with qa > 5% 
                    N3D_col(i,j) = N3D_col(i,j) + dum /min(qa_new,1.)
                    if (T_in(i,j,k) + st_out(i,j,k)  .ge. 250.) then
                      N3D_col250(i,j)  = N3D_col250(i,j) +   &
                                                     dum /min (qa_new, 1.)
                    endif
                  endif
                  gb_N3D_col(i,j) = gb_N3D_col(i,j) + dum
                endif
              end do
            end do
          end do
! end of current code
! the legacy equivalent
!NOTE: in addition to the error correction, the legacy code differs from 
!  the new code in that it is appplied to input fields rather than output i
!  fields.
!               if (ql_in(i,j,k) > qmin .and. &
!                   qa_in(i,j,k) > qmin .and. &
!                   qn_in(i,j,k) > qmin ) then      
!                 dum = qn_in(i,j,k)*Atmos_state%airdens(i,j,k)*  &
!                                    Atmos_state%deltpg(i,j,k)*1.e-6
!                 if (qa_new > 0.05) then !count only columns with qa > 5% 
!                   N3D_col(i,j) = N3D_col(i,j) + dum /min(qa_in(i,j,k),1.)
!                   if (T_in(i,j,k) + st_out(i,j,k)  .ge. 250.) then
!                     N3D_col250(i,j)  = N3D_col250(i,j) +   &
!                                               dum /min (qa_in(i,j,k), 1.)
!                   endif
!                 endif
!                 gb_N3D_col(i,j) = gb_N3D_col(i,j) + dum
!               endif
!             end do
!           end do
!         end do
! end of the legacy equivalent
          diag_3d(:,:,diag_pt%droplets_col) = N3D_col
          diag_3d(:,:,diag_pt%droplets_col250) = N3D_col250
          diag_3d(:,:,diag_pt%gb_droplets_col) = gb_N3D_col
        endif
      endif

!------------------------------------------------------------------------
!SPECIAL CASES not following general pattern of above:
!yim: in-cloud ice crystal column burden
!------------------------------------------------------------------------
      if (diag_id%nice_col > 0 .or. diag_id%gb_nice_col > 0) then
        if (present (qni_in)) then
          N3Di_col(:,:) = 0.
          gb_N3Di_col(:,:) = 0.
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                qa_new = qa_in(i,j,k) + SA_out(i,j,k)
                qi_new = qi_in(i,j,k) + SI_out(i,j,k)
                qni_new = qni_in(i,j,k) + SNi_out(i,j,k)
                if (qi_new > qmin .and. &
                    qa_new > qmin .and. &
                    qni_new  > qmin ) then
                  dum =  qni_new*Atmos_state%airdens(i,j,k)*  &
                                        Atmos_state%deltpg(i,j,k)*1.e-6
                  if (qa_new > 0.05) then !count only columns with qa > 5% 
                    N3Di_col(i,j) = N3Di_col(i,j) + dum /min(qa_new,1.)
                  endif
                  gb_N3Di_col(i,j) = gb_N3Di_col(i,j) + dum
                endif
              end do
            end do
          end do
          diag_3d(:,:,diag_pt%nice_col) = N3Di_col
          diag_3d(:,:,diag_pt%gb_nice_col) = gb_N3Di_col
        endif
      endif



!-------------------------------------------------------------------------
!    call strat_netcdf to output the requested netcdf diagnostics.
!-------------------------------------------------------------------------
      call strat_netcdf (diag_id, diag_pt, diag_4d, diag_4d_kp1, diag_3d, &
                         Time, is, js, kdim, MASK3d)

!-------------------------------------------------------------------------
!    deallocate the local allocatable variables. call strat_dealloc to 
!    deallocate the derived type variable components.
!-------------------------------------------------------------------------
      deallocate ( diag_4d )
      deallocate ( diag_4d_kp1 )
      deallocate ( diag_3d )
      call strat_dealloc (Atmos_state, Particles, Cloud_state,  &
                          Precip_state, Cloud_processes)

!------------------------------------------------------------------------
!    indicate exit from module for pe 0.
!------------------------------------------------------------------------
      if ( debugo4 .and. mpp_pe() .EQ.  0) then  
        write(outunit,*) "out stratcloud mod"
      endif
      call mpp_clock_end (sc_end)

!-----------------------------------------------------------------------


end subroutine strat_cloud_new




!#######################################################################

! <SUBROUTINE NAME="strat_cloud_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   This writes out a restart (if needed).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_cloud_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine strat_cloud_end()

!-------------------------------------------------------------------------
!---local variables-------------------------------------------------------

      integer  :: unit

!------------------------------------------------------------------------
      if (.not. module_is_initialized) return
       
!------------------------------------------------------------------------
!    write out restart file.
!------------------------------------------------------------------------
      if (do_netcdf_restart) then
        call strat_cloud_restart
      else
        if (mpp_pe() == mpp_root_pe()) then
          call mpp_error ('strat_cloud_mod',    &
                           'Writing native formatted restart file.', NOTE)
        endif
        unit = open_restart_file('RESTART/strat_cloud.res', ACTION='write')
        if (mpp_pe() == mpp_root_pe()) then
          write (unit) restart_versions(size(restart_versions(:)))
        endif
        call write_data (unit, nsum)
        call write_data (unit, qlsum)
        call write_data (unit, qisum)
        call write_data (unit, cfsum)
        call close_file (unit)
      endif
 
!--------------------------------------------------------------------------
!    call destructors for modules used by strat_cloud_mod.
!--------------------------------------------------------------------------
      if (running_old_code) then
        call strat_cloud_legacy_end (do_pdf_clouds)
      else
        call microphysics_end (Nml)
        call aerosol_cloud_end
        call nc_cond_end ( do_pdf_clouds )
        call polysvp_end
      endif
      call strat_netcdf_end

!------------------------------------------------------------------------
!    mark the module as uninitialized.
!------------------------------------------------------------------------
      module_is_initialized = .false.


end subroutine strat_cloud_end


!#######################################################################

!-----------------------------------------------------------------------
! <SUBROUTINE NAME="strat_cloud_sum">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     This increments cloud variables for passing to radiation.
!     It is expected that this will become obsolete soon.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  strat_cloud_sum (is, js, ql, qi, cf)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!        Starting integer for longitude window.
!  </IN>
!  <IN NAME="js" TYPE="integer">
!        Starting integer for latitude window.
!  </IN>
!  <IN NAME="ql" TYPE="real">
!        Cloud liquid water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!        Cloud ice water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="cf" TYPE="real">
!        Cloud fraction (fraction, 0-1)
!  </IN>
! </SUBROUTINE>
!------------------------------------------------------------------------

subroutine strat_cloud_sum (is, js, ql, qi, cf)


!-----------------------------------------------------------------------
integer,                intent(in) :: is, js
real, dimension(:,:,:), intent(in) :: ql, qi, cf
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!----locla variables-----------------------------------------------------

      integer :: ie, je

!------------------------------------------------------------------------
      if (.not.module_is_initialized) then
        call error_mesg('strat_cloud_sum',     &
                              'strat_cloud_mod is not initialized',FATAL)
      endif

!-------------------------------------------------------------------------
!    define end indices of window.
!-------------------------------------------------------------------------
      ie = is + SIZE(ql,1) - 1
      je = js + SIZE(ql,2) - 1
     
!-------------------------------------------------------------------------
!    if radiation is being supplied time-averaged cloud fields, add the
!    current values to the accumulating sum.
!-------------------------------------------------------------------------
      if (do_average) then
       nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
       qlsum(is:ie,js:je,:) = qlsum(is:ie,js:je,:) + ql
       qisum(is:ie,js:je,:) = qisum(is:ie,js:je,:) + qi
       cfsum(is:ie,js:je,:) = cfsum(is:ie,js:je,:) + cf

!-------------------------------------------------------------------------
!    if radiation is being supplied instantaneous cloud fields, save the 
!    current values in the accumulating sum arrays.
!-------------------------------------------------------------------------
      else
        nsum(is:ie,js:je)   =  1
        qlsum(is:ie,js:je,:) = ql
        qisum(is:ie,js:je,:) = qi
        cfsum(is:ie,js:je,:) = cf
      endif

!-----------------------------------------------------------------------


 end subroutine strat_cloud_sum


!#######################################################################

!-----------------------------------------------------------------------
! <SUBROUTINE NAME="strat_cloud_avg">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!      Averaging routine for cloud variables to be passed to radiation.
!      Expected to be removed shortly.
!  </DESCRIPTION>
!  <TEMPLATE>
!RSH!   call  strat_cloud_new_avg (is, js, ql, qi, cf, ierr)
!   call  strat_cloud_avg (is, js, ql, qi, cf, ierr)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!      Starting integer for longitude window.
!  </IN>
!  <IN NAME="js" TYPE="integer">
!      Starting integer for latitude window.
!  </IN>
!  <OUT NAME="ql" TYPE="real">
!      Cloud liquid water specific humidity (kg/kg)
!  </OUT>
!  <OUT NAME="qi" TYPE="real">
!      Cloud ice water specific humidity (kg/kg)
!  </OUT>
!  <OUT NAME="cf" TYPE="real">
!      Cloud fraction (0-1)
!  </OUT>
!  <OUT NAME="ierr" TYPE="integer">
!      Error integer.
!  </OUT>
! </SUBROUTINE>
!-----------------------------------------------------------------------

subroutine strat_cloud_avg (is, js, ql, qi, cf, ierr)


!-----------------------------------------------------------------------
integer,                intent(in)  :: is, js
real, dimension(:,:,:), intent(out) :: ql, qi, cf
integer,                intent(out) :: ierr
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!----local variables---------------------------------------------------

      integer :: ie, je, num, k

!-----------------------------------------------------------------------
!    be sure module is initialized. check for agreement in size of  k-index
!    of argument and of module variable holding output.
!-----------------------------------------------------------------------
      if (.not.module_is_initialized) then
        call error_mesg ('strat_cloud_avg',    &
                                  'strat_cloud is not initialized',FATAL)
      endif
      if (SIZE(ql,3) /= SIZE(qlsum,3)) then
        call error_mesg ('strat_cloud_avg in strat_cloud_mod',  &
                              'input argument has the wrong SIZE',FATAL)
      endif

!-------------------------------------------------------------------------
!    define ending indices of current window.
!-------------------------------------------------------------------------
      ie = is + SIZE(ql,1) - 1
      je = js + SIZE(ql,2) - 1

!------------------------------------------------------------------------
!    determine if all values are available. if any are missing, return an
!    error code; otherwise compute average, and return the desired fields.
!------------------------------------------------------------------------
      num = count(nsum(is:ie,js:je) == 0)
      if (num > 0) then
        ierr = 1
      else
        do k=1, SIZE(ql,3)
          ql(:,:,k) = qlsum(is:ie,js:je,k)/float(nsum(is:ie,js:je))
          qi(:,:,k) = qisum(is:ie,js:je,k)/float(nsum(is:ie,js:je))
          cf(:,:,k) = cfsum(is:ie,js:je,k)/float(nsum(is:ie,js:je))
        enddo
        ierr = 0
      endif

!-------------------------------------------------------------------------
!    initialize the accumulation arrays so they are ready for the next 
!    access.
!-------------------------------------------------------------------------
      nsum (is:ie,js:je)   = 0
      qlsum(is:ie,js:je,:) = 0.0
      qisum(is:ie,js:je,:) = 0.0
      cfsum(is:ie,js:je,:) = 0.0

!-----------------------------------------------------------------------


 end subroutine strat_cloud_avg


!#######################################################################


!-------------------------------------------------------------------------
! <FUNCTION NAME="do_strat_cloud">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     Logical function to indicate whether or not strat_cloud is running.
!  </DESCRIPTION>
!  <TEMPLATE>
!   result =  do_strat_cloud ( ) result (answer)
!
!  </TEMPLATE>
! </FUNCTION>
!-------------------------------------------------------------------------

function do_strat_cloud ( ) result (answer)


logical :: answer

!------------------------------------------------------------------------
!    define output value; it will be the value of module variable
!    strat_cloud_on.
!------------------------------------------------------------------------
      answer = strat_cloud_on


end function do_strat_cloud


!##########################################################################

!-------------------------------------------------------------------------
! <SUBROUTINE NAME="strat_cloud_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents 
!                                      the model time, used for writing 
!                                      restart. timestamp will append to
!                                      any restart file name as a prefix. 
! </DESCRIPTION>
!
! </SUBROUTINE> 
!-------------------------------------------------------------------------

subroutine strat_cloud_restart(timestamp)

character(len=*), intent(in), optional :: timestamp

!-----------------------------------------------------------------------
!    write message indicating status of restart files. call routine to
!    output restart fields.
!-----------------------------------------------------------------------
      if (do_netcdf_restart) then
        if (mpp_pe() == mpp_root_pe()) then
          call mpp_error ('strat_cloud_mod',   &
            'Writing netCDF formatted restart file:  &
                                      &RESTART/strat_cloud.res.nc', NOTE)
        endif
        call save_restart(Str_restart, timestamp)
        if (in_different_file) call  save_restart(Til_restart, timestamp)
      else
        call error_mesg ('strat_cloud_mod', &
            'Native intermediate restart files are not supported.', FATAL)
      endif

!------------------------------------------------------------------------

end subroutine strat_cloud_restart




!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!
!                         PRIVATE INTERFACES
!
!
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||





!########################################################################

subroutine fill_nml_variable

!------------------------------------------------------------------------
!    place the strat_cloud_nml variables in strat_nml_type variable Nml.
!------------------------------------------------------------------------
      Nml%do_netcdf_restart = do_netcdf_restart
      Nml%U00 = U00
      Nml%U00_profile = U00_profile
      Nml%rthresh = rthresh
      Nml%use_kk_auto = use_kk_auto 
      Nml%var_limit = var_limit
      Nml%use_online_aerosol = use_online_aerosol
      Nml%sea_salt_scale = sea_salt_scale
      Nml%om_to_oc = om_to_oc
      Nml%N_land = N_land
      Nml%use_sub_seasalt = use_sub_seasalt
      Nml%N_ocean = N_ocean
      Nml%U_evap = U_evap
      Nml%eros_scale = eros_scale
      Nml%eros_choice = eros_choice
      Nml%eros_scale_t = eros_scale_t
      Nml%eros_scale_c = eros_scale_c
      Nml%mc_thresh = mc_thresh
      Nml%diff_thresh = diff_thresh
      Nml%super_choice = super_choice
      Nml%tracer_advec = tracer_advec
      Nml%qmin = qmin
      Nml%Dmin = Dmin
      Nml%num_strat_pts = num_strat_pts

      if (num_strat_pts > 0) then
        allocate (Nml%strat_pts(2, num_strat_pts))
        Nml%strat_pts = strat_pts
      endif

      Nml%efact = efact
      Nml%vfact = vfact
      Nml%cfact = cfact
      Nml%do_old_snowmelt = do_old_snowmelt
      Nml%retain_cm3_bug  = retain_cm3_bug 
      Nml%do_pdf_clouds = do_pdf_clouds
      Nml%betaP = betaP
      Nml%iwc_crit = iwc_crit
      Nml%vfall_const2 = vfall_const2
      Nml%vfall_exp2 = vfall_exp2
      Nml%qthalfwidth = qthalfwidth
      Nml%nsublevels = nsublevels
      Nml%kmap = kmap
      Nml%kord = kord
      Nml%do_liq_num = do_liq_num
      Nml%do_dust_berg = do_dust_berg
      Nml%N_min = N_min
      Nml%num_mass_ratio1 = num_mass_ratio1
      Nml%num_mass_ratio2 = num_mass_ratio2
      Nml%microphys_scheme = microphys_scheme              
      Nml%macrophys_scheme = macrophys_scheme              
      Nml%aerosol_activation_scheme = aerosol_activation_scheme        
      Nml%mass_cons = mass_cons
      Nml%super_ice_opt = super_ice_opt
      Nml%pdf_org = pdf_org
      Nml%do_ice_nucl_wpdf = do_ice_nucl_wpdf
      Nml%do_hallet_mossop = do_hallet_mossop
      Nml%activate_all_ice_always = activate_all_ice_always
      Nml%debugo = debugo
      Nml%isamp = isamp
      Nml%jsamp = jsamp
      Nml%ksamp = ksamp

      Nml%qcvar = qcvar

!----------------------------------------------------------------------

end subroutine fill_nml_variable



!#######################################################################

subroutine strat_debug (otun, ST_out, SQ_out, Cloud_state, Precip_state,  &
                        Atmos_State)

integer,                 intent(in) :: otun
real, dimension(:,:,:),  intent(in) :: ST_out, SQ_out
type(cloud_state_type),  intent(in) :: Cloud_state
type(atmos_state_type),  intent(in) :: Atmos_state
type(precip_state_type), intent(in) :: Precip_state

      integer, dimension(3) :: maxl, minl

!------------------------------------------------------------------------
!    write numerous diagnostics to file otun to aid in debugging.
!------------------------------------------------------------------------
      write(otun, *)  "eee max, min ST ", MAXVAL(ST_out), MINVAL(ST_out)
      write(otun, *)  "eee maxloc, minloc  ", MAXLOC(ST_out),   &
                                                            MINLOC(ST_out)

      call check_nan (ST_out,'ST_out')


      write(otun, *)  "eee max, min SQ ", MAXVAL(SQ_out), MINVAL(SQ_out)
      write(otun, *)  "eee maxloc, minloc  ", MAXLOC(SQ_out),   &
                                                            MINLOC(SQ_out)
      call check_nan (SQ_out,'SQ_out')

      write(otun, *)  "eee max, min SL ", MAXVAL(Cloud_state%SL_out), &
                                                 MINVAL(Cloud_state%SL_out)
      write(otun, *)  "eee maxloc, minloc  ", MAXLOC(Cloud_state%SL_out),&
                                                 MINLOC(Cloud_state%SL_out)

      write(otun, *)  "eee max, min SI ", MAXVAL(Cloud_state%SI_out), &
                                                 MINVAL(Cloud_state%SI_out)
      write(otun, *)  "eee maxloc, minloc  ", MAXLOC(Cloud_state%SI_out),&
                                                 MINLOC(Cloud_state%SI_out)
      call check_nan   (Cloud_state%SI_out,'SI_out')

      write(otun, *)  "eee max, min SA ", MAXVAL(Cloud_state%SA_out), &
                                                MINVAL(Cloud_state%SA_out)
      write(otun, *)  "eee maxloc, minloc  ", MAXLOC(Cloud_state%SA_out),&
                                                 MINLOC(Cloud_state%SA_out)
      call check_nan   (Cloud_state%SA_out,'SA_out')

      write(otun, *)  "eee max, min SN ", MAXVAL(Cloud_state%SN_out), &
                                                MINVAL(Cloud_state%SN_out)
      write(otun, *)  "eee maxloc, minloc  ", &
                     MAXLOC(Cloud_state%SN_out), MINLOC(Cloud_state%SN_out)
      call check_nan   (Cloud_state%SN_out,'SN_out')

      write(otun, *)  "eee max, min SNi ", MAXVAL(Cloud_state%SNi_out),&
                                                MINVAL(Cloud_state%SNi_out)
      write(otun, *)  "eee maxloc, minloc  ", &
                   MAXLOC(Cloud_state%SNi_out), MINLOC(Cloud_state%SNi_out)
      call check_nan   (Cloud_state%SNi_out,'SNi_out')

      write(otun, *)  "--"
      write(otun, *)  "eee max, min T+ST ",  &
                                 MAXVAL(Atmos_state%T_in + ST_out),  &  
                                         MINVAL(Atmos_state%T_in + ST_out)

      write(otun, *)  "eee max, min qv+SQ ", &
                                 MAXVAL(Atmos_state%qv_in + SQ_out), &
                                        MINVAL(Atmos_state%qv_in + SQ_out)

      write(otun, *)  "eee max, min ql+ SL ",   &
                           MAXVAL(Cloud_state%ql_in+Cloud_state%SL_out),  &
                               MINVAL(Cloud_state%ql_in+Cloud_state%SL_out)

      write(otun, *)  "eee max, min qi +SI ",  &
                           MAXVAL(Cloud_state%qi_in+ Cloud_state%SI_out), &
                             MINVAL(Cloud_state%qi_in + Cloud_state%SI_out)

      write(otun, *)  "eee max, min qa + SA ",  &
                           MAXVAL(Cloud_state%qa_in+Cloud_state%SA_out),  &
                              MINVAL(Cloud_state%qa_in+Cloud_state%SA_out)

      write(otun, *)  "eee max, min qn + SN ",  &
                         MAXVAL(Cloud_state%qn_in + Cloud_state%SN_out), &
                            MINVAL(Cloud_state%qn_in + Cloud_state%SN_out)

      write(otun, *)  "eee max, min qni SNi ",   &
                        MAXVAL(Cloud_state%qni_in+ Cloud_state%SNi_out), &
                            MINVAL(Cloud_state%qni_in+Cloud_state%SNi_out)

      write(otun, *)  "--"
      write(otun, *)  "--"

      write(otun, *)  "eee max, min rain3d ",  &
                  MAXVAL(Precip_state%rain3d), MINVAL(Precip_state%rain3d)
      write(otun, *)  "eee maxloc, minloc  ",  &
                  MAXLOC(Precip_state%rain3d), MINLOC(Precip_state%rain3d)
      call check_nan   (Precip_state%rain3d,'rain3d')

      write(otun, *)  "eee max, min snow3d ",   &
                  MAXVAL(Precip_state%snow3d), MINVAL(Precip_state%snow3d)
      write(otun, *)  "eee maxloc, minloc  ",  &
                   MAXLOC(Precip_state%snow3d), MINLOC(Precip_state%snow3d)
      call check_nan   (Precip_state%snow3d,'snow3d')

      write(otun, *)  "--"
      write(otun, *)  "eee max, min surfrain ",   &
               MAXVAL(Precip_state%surfrain), MINVAL(Precip_state%surfrain)
      write(otun, *)  "eee maxloc, minloc  ",    &
               MAXLOC(Precip_state%surfrain), MINLOC(Precip_state%surfrain)

      write(otun, *)  "eee max, min surfsnow ",   &
               MAXVAL(Precip_state%surfsnow), MINVAL(Precip_state%surfsnow)
      write(otun, *)  "eee maxloc, minloc  ",   &
               MAXLOC(Precip_state%surfsnow), MINLOC(Precip_state%surfsnow)
      write(otun, *)  "--"

      IF ( MAXVAL(SQ_out + Cloud_state%ql_in) .GT. 1.e-1 )  &
                                            write(otun, *) " MMMMM Q1 "
      IF ( MAXVAL(SQ_out + Cloud_state%ql_in) .LT. 0. )   &
                                            write(otun, *) " MMMMM Q2 "

      IF ( MAXVAL(Cloud_state%SI_out + Cloud_state%qi_in) .LT. 0. )  &
                                            write(otun, *) " MMMMM I11 "
      IF ( MAXVAL(Cloud_state%SL_out + Cloud_state%ql_in) .LT. 0. )   &
                                            write(otun, *) " MMMMM L11 "

      IF ( MAXVAL(  &
              Cloud_state%qa_in+Cloud_state%SA_out+Atmos_state%ahuco)   &
                                               .GT. 1.00000000001  ) THEN
        write(otun, *) " MMMMMA1 ahuco "
        maxl =  maxloc (Cloud_state%qa_in + Cloud_state%SA_out +  &
                                                        Atmos_state%ahuco) 
        write(otun, *) "  maxloc(qa+SA) ",  &
             maxloc (Cloud_state%qa_in + Cloud_state%SA_out +  &
                                                         Atmos_state%ahuco)
        write(otun, *) " qa+SA+ahuco ",    &
                      Cloud_state%qa_in(maxl(1),maxl(2),maxl(3)) +     &
                          Cloud_state%SA_out(maxl(1),maxl(2),maxl(3)) +&
                              Atmos_state%ahuco(maxl(1),maxl(2),maxl(3))
        write(otun, *) " qa, ahuco, SA ",   &
                     Cloud_state%qa_in(maxl(1),maxl(2),maxl(3)),   &
                     Atmos_state%ahuco(maxl(1),maxl(2),maxl(3)),  &
                     Cloud_state%SA_out(maxl(1),maxl(2),maxl(3)) 

      END IF

      IF ( MINVAL(Cloud_state%qa_in+Cloud_state%SA_out) .LT. 0.  ) THEN
        minl =  minloc(Cloud_state%qa_in+Cloud_state%SA_out)
        write(otun, *) " MMMMMA2"
        write(otun, *) "  minloc(qa+SA) ",   &
                              minloc(Cloud_state%qa_in+Cloud_state%SA_out)
        write(otun, *) " qa+SA ", &
                         Cloud_state%qa_in(minl(1),minl(2),minl(3)) +  &
                                Cloud_state%SA_out(minl(1),minl(2),minl(3))
        write(otun, *) "    qa, SA ",  &
                           Cloud_state%qa_in(minl(1),minl(2),minl(3)),  &
                                Cloud_state%SA_out(minl(1),minl(2),minl(3))
      END IF

      IF ( MINVAL(Cloud_state%qn_in+Cloud_state%SN_out) .LT. 0.  ) THEN
        write(otun, *) " MMMMMN1"
        minl =  minloc(Cloud_state%qn_in+Cloud_state%SN_out)
        write(otun, *) "  minloc(qn+SN) ",  &
                           minloc(Cloud_state%qn_in+Cloud_state%SN_out)
        write(otun, *) "    qn, SN ",   &
                          Cloud_state%qn_in(minl(1),minl(2),minl(3)),  &
                          Cloud_state%SN_out(minl(1),minl(2),minl(3)) 
      END IF

      IF ( MAXVAL(Cloud_state%qi_in +    &
                                   Cloud_state%SI_out) .GT. 1.e-2 )   then
        maxl = MAXLOC(Cloud_state%qi_in + Cloud_state%SI_out)
        write(otun, *) " MMMMMII"
        write(otun, *) "  maxloc(qi+SI) ",    &
                              maxloc(Cloud_state%qi_in+Cloud_state%SI_out)
        write(otun, *) "  qi+SI " ,    &
                           Cloud_state%qi_in(maxl(1),maxl(2),maxl(3))+  &
                                Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
        write(otun, *) "  qi, SI " ,    &
                      Cloud_state%qi_in(maxl(1),maxl(2),maxl(3)),  &
                               Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
        write(otun, *) "  T ",    &
                     Atmos_state%T_in(maxl(1),maxl(2),maxl(3)),   &
                               Cloud_state%SI_out(maxl(1),maxl(2),maxl(3))
      END IF

      IF ( MAXVAL(Cloud_state%ql_in + Cloud_state%SL_out) .GT. 1.e-2 )   &
                                               write(otun, *) " MMMMMLL"
      IF ( MINVAL(Cloud_state%qni_in+Cloud_state%SNi_out) .LT. -1.e-5) &
                                                write(otun, *) " MMMMMN2"
      IF ( MAXVAL(ST_out) .GT. 7. ) write(otun, *) " MMMMMT1 "
      IF ( MINVAL(ST_out) .LT. - 7. ) write(otun, *) " MMMMMT2 "
      IF ( MAXVAL(Atmos_state%T_in+ST_out) .GT. 330. )   &
                                            write(otun, *) " MMMMMT3 "
      IF ( MINVAL(Atmos_state%T_in+ST_out) .LT. 170. )    &
                                               write(otun, *) " MMMMMT4 "

      IF  ( MINVAL(Precip_state%rain3d) .LT. 0. )   &
                                                write(otun, *) " MMMMMR1 "
      IF  ( MINVAL(Precip_state%snow3d) .LT. 0. )   &
                                                write(otun, *) " MMMMMS1 "

      IF ( MINVAL(Precip_state%surfrain) .LT. 0. )   &
                                                write(otun, *) " MMMMMX1 "
      IF ( MINVAL(Precip_state%surfsnow) .LT. 0. )   &
                                                write(otun, *) " MMMMMX2 "

!-------------------------------------------------------------------------


end subroutine strat_debug 



!#######################################################################

subroutine impose_realizability (idim, jdim, kdim, Atmos_state,  &
                                 Cloud_state, SQ_out, ST_out, N3d, N3Di,  &
                                 otun, n_diag_4d, diag_4d)

!-----------------------------------------------------------------------
!    account for the fact that other processes may have created negative 
!    tracer or extremely small values of tracer fields. the general reason
!    for the extremely small values of the tracer fields is due to vertical
!    diffusion, advection of condensate or cumulus induced subsidence (also
!    a form of advection) of condensate.
!    in this step any values of the prognostic variables which are less 
!    than qmin are reset to zero, while conserving total moisture.
!    note that this is done slightly different for the Tiedtke cloud 
!    fraction than it is for pdf clouds. In the former, the filling 
!    requires that cloud liquid, cloud ice, and cloud fraction are greater
!    than qmin. For PDF clouds, cloud fraction need not be considered 
!    since it is diagnosed later from the PDF cloud field.
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
integer,                                     intent(in)    :: idim, jdim, &
                                                              kdim
type( atmos_state_type),                     intent(inout) :: Atmos_state
type( cloud_state_type),                     intent(inout) :: Cloud_state
real, dimension(idim,jdim,kdim),             intent(inout) :: ST_out,  &
                                                              SQ_out, &
                                                              N3D, N3Di
integer,                                     intent(in)    :: otun,  &
                                                              n_diag_4d
real, dimension(idim,jdim,kdim,0:n_diag_4d), intent(inout) :: diag_4d
!-------------------------------------------------------------------------

!------------------------------------------------------------------------
!---local variables------------------------------------------------------

      logical, dimension(idim, jdim, kdim) :: ql_too_small, qi_too_small 
      integer :: i,j,k

!-----------------------------------------------------------------------
!    for the non-pdf scheme,  assure that cloud fraction is greater than 
!    qmin.  if it is not, set it to 0.0 and save the tendency and updated
!    value. save a diagnostic if requested.
!-----------------------------------------------------------------------
      if (.not. Nml%do_pdf_clouds) then
        where (Cloud_state%qa_in .le. Nml%qmin)
          Cloud_state%SA_out= Cloud_state%SA_out - Cloud_state%qa_in
          Cloud_state%qa_upd = 0.
        elsewhere
          Cloud_state%qa_upd = Cloud_state%qa_in     
        end where
        if (diag_id%qadt_fill + diag_id%qa_fill_col > 0 ) then
          where (Cloud_state%qa_in .le. Nml%qmin)
            diag_4d(:,:,:,diag_pt%qadt_fill) =  -Cloud_state%qa_in*   &
                                                    Constants% inv_dtcloud
          endwhere
        end if

!------------------------------------------------------------------------ 
!    define the max cloud area permitted (U01), which is equal to the 
!    grid box RH, under the assumption that the cloudy air is saturated 
!    and the temperature inside and outside of the cloud are ~ the same.
!    save a diagnostic indicating the amount cloud area was reduced due
!    to this requirenment,
!------------------------------------------------------------------------ 
        Atmos_state%U01 = min(Atmos_state%U01, 1.)
        if (diag_id%qadt_rhred + diag_id%qa_rhred_col >0) then
          where (Cloud_state%qa_upd .gt. Atmos_state%U01)
            diag_4d(:,:,:,diag_pt%qadt_rhred) = -(Cloud_state%qa_upd -  &
                                Atmos_state%U01)*Constants%inv_dtcloud
          endwhere
        end if
        
        where (Cloud_state%qa_upd .gt. Atmos_state%U01)
          Cloud_state%SA_out = Cloud_state%SA_out + Atmos_state%U01 -   &
                                                      Cloud_state%qa_upd
          Cloud_state%qa_upd = Atmos_state%U01      
        end where
      endif

!-------------------------------------------------------------------------
!    define the conditions under which liquid and ice filling must be done,
!    for both the pdf and non-pdf cases.
!    for the non-pdf scheme, the filling requires that cloud liquid, cloud
!    ice, and cloud fraction are greater than qmin.
!    for pdf clouds, cloud fraction need not be considered since it is 
!    diagnosed later from the PDF cloud field.
!-------------------------------------------------------------------------
      if (.not. Nml%do_pdf_clouds) then
        if (Nml%do_liq_num) then
          ql_too_small = (Cloud_state%ql_in .le. Nml%qmin .or.   &
                          Cloud_state%qa_in .le. Nml%qmin .or.   &
                          Cloud_state%qn_in .le. Nml%qmin)
        else
          ql_too_small = (Cloud_state%ql_in .le. Nml%qmin .or.   &
                          Cloud_state%qa_in .le. Nml%qmin)
        endif
        if (Constants%do_predicted_ice_number) then
             qi_too_small = (Cloud_state%qi_in .le.  Nml%qmin .or.   &
                          Cloud_state%qa_in .le. Nml%qmin .or.   &
                          Cloud_state%qni_in .le. Nml%qmin)
        else
          qi_too_small = (Cloud_state%qi_in .le.  Nml%qmin .or.   &
                          Cloud_state%qa_in .le. Nml%qmin)
        endif
      else
        if ( Nml%do_liq_num) then
          ql_too_small = (Cloud_state%ql_in .le.  Nml%qmin  .or.   &
                          Cloud_state%qn_in .le.  Nml%qmin)
        else
          ql_too_small = (Cloud_state%ql_in .le.  Nml%qmin)
        endif
        if (Constants%do_predicted_ice_number) then
             qi_too_small = (Cloud_state%qi_in .le.  Nml%qmin .or.   &
                          Cloud_state%qni_in .le. Nml%qmin)
        else
          qi_too_small = (Cloud_state%qi_in .le.  Nml%qmin )
        endif
      endif
  
!------------------------------------------------------------------------
!    determine if ql needs filling. if so, fill. save a diagnostic if 
!    requested. define tendency terms and update cloud liquid fields. 
!    include temperature and specific humidity adjustments to assure 
!    conservation. 
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (ql_too_small(i,j,k)) then
              Cloud_state%SL_out(i,j,k) = Cloud_state%SL_out(i,j,k) -   &
                                                   Cloud_state%ql_in(i,j,k)
              SQ_out(i,j,k) = SQ_out(i,j,k) + Cloud_state%ql_in(i,j,k)
              ST_out(i,j,k) = ST_out(i,j,k) -    &
                                        hlv*Cloud_state%ql_in(i,j,k)/cp_air
              Cloud_state%ql_upd(i,j,k) = 0.
            else
              Cloud_state%ql_upd(i,j,k) = Cloud_state%ql_in(i,j,k)
            endif
          end do
        end do
      end do
      if ( diag_id%qldt_fill  + diag_id%ql_fill_col > 0 ) then
        where (ql_too_small )
          diag_4d(:,:,:,diag_pt%qldt_fill) =   &
                            -1.*Cloud_state%ql_in*Constants%inv_dtcloud
        endwhere
      end if
      if ( diag_id%qdt_liquid_init > 0 ) then
        where (ql_too_small )
          diag_4d(:,:,:,diag_pt%qdt_liquid_init) =   &
                            Cloud_state%ql_in*Constants%inv_dtcloud
        endwhere
      end if

!------------------------------------------------------------------------
!    adjust cloud droplet numbers as needed when those fields
!    are being predicted. include diagnostics if requested. if droplet
!    number not being predicted, values are as set above, dependent on
!    whether column is over land or ocean.
!------------------------------------------------------------------------
      if (Nml%do_liq_num) then
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (ql_too_small(i,j,k)) then
                Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -  &
                                                 Cloud_state%qn_in(i,j,k)
                Cloud_state%qn_upd(i,j,k) = 0.
                N3D(i,j,k)    = 0.
              else
                Cloud_state%qn_upd(i,j,k) = Cloud_state%qn_in(i,j,k)
                N3D(i,j,k) = Cloud_state%qn_in(i,j,k)*   &
                                       Atmos_state%airdens(i,j,k)*1.e-6
              endif
            end do
          end do
        end do
        if (diag_id%cf_liq_init   > 0)  then
          do k = 1,kdim
            do j=1,jdim
              do i = 1,idim
                if (ql_too_small(i,j,k)) then
                  diag_4d(i,j,k,diag_pt%cf_liq_init  ) = 0.
                else
                   diag_4d(i,j,k,diag_pt%cf_liq_init  ) =    &
                                         min(Cloud_state%qa_in(i,j,k),1.)
                endif
              end do
            end do
          end do
        endif 
        if ( diag_id%qndt_fill  + diag_id%qn_fill_col + &
                    diag_id%qldt_fill + diag_id%ql_fill_col > 0 ) then
          where (ql_too_small )
            diag_4d(:,:,:,diag_pt%qndt_fill) =    &
                           -1.*Cloud_state%qn_in*Constants%inv_dtcloud
          endwhere
        end if
      endif

!------------------------------------------------------------------------
!    determine if qi needs filling. if so, fill. save a diagnostic if 
!    requested. define tendency terms and update cloud ice fields. 
!    include temperature and specific humidity adjustments to assure 
!    conservation. adjust ice particle numbers as needed when those fields
!    are being predicted.
!------------------------------------------------------------------------
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            if (qi_too_small(i,j,k)) then
              Cloud_state%SI_out(i,j,k) = Cloud_state%SI_out(i,j,k) -   &
                                                   Cloud_state%qi_in(i,j,k)
              SQ_out(i,j,k) = SQ_out(i,j,k) + Cloud_state%qi_in(i,j,k)
              ST_out(i,j,k) = ST_out(i,j,k) -    &
                                        hls*Cloud_state%qi_in(i,j,k)/cp_air
              Cloud_state%qi_upd(i,j,k) = 0.
            else
              Cloud_state%qi_upd(i,j,k) = Cloud_state%qi_in(i,j,k)
            endif
          end do
        end do
      end do
      if (diag_id%qidt_fill  + diag_id%qi_fill_col > 0 ) then
        where (qi_too_small )
          diag_4d(:,:,:,diag_pt%qidt_fill) =     &
                              -1.*Cloud_state%qi_in *Constants%inv_dtcloud
        endwhere
      end if
      if (diag_id%qdt_ice_init > 0 ) then
        where (qi_too_small )
          diag_4d(:,:,:,diag_pt%qdt_ice_init) =     &
                             Cloud_state%qi_in *Constants%inv_dtcloud
        endwhere
      end if
      if (diag_id%cf_ice_init > 0) then
        where (qi_too_small)
          diag_4d(:,:,:,diag_pt%cf_ice_init) = 0.
        elsewhere
          diag_4d(:,:,:,diag_pt%cf_ice_init) = Cloud_state%qa_in
        end where
      endif
      if (Constants%do_mg_microphys .or. &
              Constants%do_mg_ncar_microphys) then
        if (Nml%debugo) then
          write(otun, *) " SNi 00 ",   &
                Cloud_state%SNi_out(Nml%isamp,Nml%jsamp,NMl%ksamp)*  &
                                                    Constants%inv_dtcloud
        end if
        where (qi_too_small) 
          Cloud_state%SNi_out  = Cloud_state%SNi_out - Cloud_state%qni_in
          Cloud_state%qni_upd=0.
          N3Di = 0.
        elsewhere
          Cloud_state%qni_upd = Cloud_state%qni_in
          N3Di  = Cloud_state%qni_in*Atmos_state%airdens*1.e-6
        end where

        if (diag_id%qnidt_fill  + diag_id%qni_fill_col > 0 ) then
          where (qi_too_small) 
            diag_4d(:,:,:,diag_pt%qnidt_fill) =  &
                           -1.*Cloud_state%qni_in *Constants%inv_dtcloud
          endwhere
        end if

        if (Nml%debugo) then
          write(otun, *) " SNi 01 ",   &
            Cloud_state%SNi_out(Nml%isamp,Nml%jsamp,Nml%ksamp)*   &
                                                  Constants%inv_dtcloud
          write(otun, *) "     qnidt_fill ",   &
               diag_4d(Nml%isamp,Nml%jsamp,Nml%ksamp,diag_pt%qnidt_fill) 
        end if
      else
        N3Di = 0.
      endif

!------------------------------------------------------------------------
!    save the cloud area tendency and updated area values at this point
!    for later use in the m-g microphysics.
!------------------------------------------------------------------------
      if (Constants%do_mg_microphys .or. &
              Constants%do_mg_ncar_microphys) then
        Cloud_state%qa_upd_0 = Cloud_state%qa_upd
        Cloud_state%SA_0 = Cloud_State%SA_out
      endif

!-------------------------------------------------------------------------

end subroutine impose_realizability 




!#########################################################################

subroutine strat_alloc (idim, jdim, kdim, pfull, phalf, zhalf, &
                        zfull, radturbten2, T_in, qv_in, ql_in, qi_in, &
                        qa_in, omega, Mc, diff_t, qrat, ahuco, &
                        Atmos_state, Particles, Cloud_state, Precip_state,&
                        Cloud_processes, qn_in, qni_in)

integer,                    intent(in)           :: idim, jdim, kdim
real,dimension(:,:,:),      intent(in)           ::    &
                      pfull, phalf, zfull, zhalf, radturbten2, &
                      T_in, qv_in, ql_in, qi_in, qa_in, omega, Mc,    &
                      diff_t, qrat, ahuco
type(atmos_state_type),     intent(inout)        :: Atmos_state
type(cloud_state_type),     intent(inout)        :: Cloud_state
type(precip_state_type),    intent(inout)        :: Precip_state
type(cloud_processes_type), intent(inout)        :: Cloud_processes
type(particles_type),       intent(inout)        :: Particles
real,dimension(:,:,:),      intent(in), optional :: qn_in, qni_in    


!-------------------------------------------------------------------------
!----local variables------------------------------------------------------

      integer :: k

!-----------------------------------------------------------------------
!    allocate and initialize the components of the atmos_state_type 
!    variable Atmos_state.
!-----------------------------------------------------------------------
      allocate (Atmos_state%pfull          (idim, jdim, kdim) )
      allocate (Atmos_state%phalf          (idim, jdim, kdim+1) )
      allocate (Atmos_state%zfull          (idim, jdim, kdim) )
      allocate (Atmos_state%zhalf          (idim, jdim, kdim+1) )
      allocate (Atmos_state%radturbten2    (idim, jdim, kdim) )
      allocate (Atmos_state%T_in           (idim, jdim, kdim) )
      allocate (Atmos_state%qv_in          (idim, jdim, kdim) )
      allocate (Atmos_state%omega          (idim, jdim, kdim) )
      allocate (Atmos_state%Mc             (idim, jdim, kdim) )
      allocate (Atmos_state%diff_t         (idim, jdim, kdim) )
      allocate (Atmos_state%qrat           (idim, jdim, kdim) )
      allocate (Atmos_state%ahuco          (idim, jdim, kdim) )
      allocate (Atmos_state%airdens        (idim, jdim, kdim) )
      allocate (Atmos_state%tn             (idim, jdim, kdim) )
      allocate (Atmos_state%qvn            (idim, jdim, kdim) )
      allocate (Atmos_state%qs             (idim, jdim, kdim) )
      allocate (Atmos_state%dqsdT          (idim, jdim, kdim) )
      allocate (Atmos_state%qsi            (idim, jdim, kdim) )
      allocate (Atmos_state%qsl            (idim, jdim, kdim) )
      allocate (Atmos_state%rh_crit        (idim, jdim, kdim) )
      allocate (Atmos_state%rh_crit_min    (idim, jdim, kdim) )
      allocate (Atmos_state%gamma          (idim, jdim, kdim) )
      allocate (Atmos_state%esat0          (idim, jdim, kdim) )
      allocate (Atmos_state%U_ca           (idim, jdim, kdim) )
      allocate (Atmos_state%delp           (idim, jdim, kdim) )
      allocate (Atmos_state%U01            (idim, jdim, kdim) )
      allocate (Atmos_state%deltpg         (idim, jdim, kdim) )

      Atmos_state%pfull       = pfull
      Atmos_state%phalf       = phalf
      Atmos_state%zhalf       = zhalf
      Atmos_state%zfull       = zfull
      Atmos_state%radturbten2 = radturbten2
      Atmos_state%T_in        = T_in
      Atmos_state%qv_in       = qv_in
      Atmos_state%omega       = omega
      Atmos_state%Mc          = Mc  
      Atmos_state%diff_t      = diff_t
      Atmos_state%qrat        = qrat
      Atmos_state%ahuco       = ahuco

!-----------------------------------------------------------------------
!    calculate pressure thickness of layer.
!-----------------------------------------------------------------------
      do k=1, kdim
        Atmos_state%deltpg(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/grav
      end do         

!-----------------------------------------------------------------------
!    calculate air density.   
!-----------------------------------------------------------------------
      where (qrat .gt. 0.) 
        Atmos_state%airdens = pfull/(rdgas*T_in*  &
                                (1. + (d608*qv_in/qrat) - ql_in - qi_in) )
      elsewhere
        Atmos_state%airdens = pfull/(rdgas*T_in*(1. - ql_in - qi_in) )
      end where

      Atmos_state%tn            = T_in
      Atmos_state%qvn           = 0.
      Atmos_state%qs            = 0.
      Atmos_state%dqsdT         = 0.
      Atmos_state%qsi           = 0.
      Atmos_state%qsl           = 0.
      Atmos_state%rh_crit       = 1.
      Atmos_state%rh_crit_min   = 1.
      Atmos_state%gamma         = 0.
      Atmos_state%esat0         = 0.
      Atmos_state%U_ca          = 0.
      do k = 1, kdim
        Atmos_state%delp(:,:,k) = phalf(:,:,k+1) - phalf(:,:,k)
      enddo     
      Atmos_state%U01           = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the particles_type 
!    variable Particles.
!-----------------------------------------------------------------------
      allocate (Particles%concen_dust_sub   (idim, jdim, kdim) )
      allocate (Particles%drop1             (idim, jdim, kdim) )
      allocate (Particles%drop2             (idim, jdim, kdim) )
      allocate (Particles%crystal1          (idim, jdim, kdim) )
      allocate (Particles%rbar_dust         (idim, jdim, kdim) )
      allocate (Particles%ndust             (idim, jdim, kdim) )
      allocate (Particles%hom               (idim, jdim, kdim) )

      Particles%concen_dust_sub   = 0.
      Particles%drop1    = 0.
      Particles%drop2    = 0.
      Particles%crystal1    = 0.
      Particles%rbar_dust   = 0.
      Particles%ndust   = 0.
      Particles%hom   = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the cloud_state_type 
!    variable Cloud_state.
!-----------------------------------------------------------------------
      allocate (Cloud_state%ql_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qi_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qa_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qn_upd    (idim, jdim, kdim) )
      allocate (Cloud_state%qni_upd   (idim, jdim, kdim) )
  
      allocate (Cloud_state%ql_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qi_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qa_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qn_mean   (idim, jdim, kdim) )
      allocate (Cloud_state%qni_mean  (idim, jdim, kdim) )
 
      allocate (Cloud_state%ql_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qi_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qa_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qn_in     (idim, jdim, kdim) )
      allocate (Cloud_state%qni_in    (idim, jdim, kdim) )

      allocate (Cloud_state%SL_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SI_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SA_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SN_out    (idim, jdim, kdim) )
      allocate (Cloud_state%SNi_out   (idim, jdim, kdim) )

      allocate (Cloud_state%qa_upd_0     (idim, jdim, kdim) )
      allocate (Cloud_state%SA_0         (idim, jdim, kdim) )

      Cloud_state%ql_upd    = 0.
      Cloud_state%qi_upd    = 0.
      Cloud_state%qa_upd    = 0.
      Cloud_state%qn_upd    = 0.
      Cloud_state%qni_upd   = 0.

      Cloud_state%ql_mean    = 0.      
      Cloud_state%qi_mean    = 0.      
      Cloud_state%qa_mean    = 0.      
      Cloud_state%qn_mean    = 0.      
      Cloud_state%qni_mean   = 0.      

      Cloud_state%ql_in = ql_in
      Cloud_state%qi_in = qi_in
      Cloud_state%qa_in = qa_in
      if (present(qn_in)) then
        Cloud_state%qn_in = qn_in
      else
        Cloud_state%qn_in = 0.
      end if

      if (present(qni_in)) then
        Cloud_state%qni_in  = qni_in
      else
        Cloud_state%qni_in  = 0.
      endif

      Cloud_state%SL_out  = 0.
      Cloud_state%SI_out  = 0.
      Cloud_state%SA_out  = 0.
      Cloud_state%SN_out  = 0.
      Cloud_state%SNi_out = 0.
    
      Cloud_state%qa_upd_0 = 0.
      Cloud_state%SA_0        = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the precip_state_type 
!    variable Precip_state.
!-----------------------------------------------------------------------
      allocate (Precip_state%lsc_snow      (idim, jdim, kdim) )
      allocate (Precip_state%lsc_rain      (idim, jdim, kdim) )
      allocate (Precip_state%lsc_snow_size (idim, jdim, kdim) )
      allocate (Precip_state%lsc_rain_size (idim, jdim, kdim) )
      allocate (Precip_state%qsout3d_mg    (idim, jdim, kdim) )
      allocate (Precip_state%qrout3d_mg    (idim, jdim, kdim) )
      allocate (Precip_state%rain3d        (idim, jdim, kdim+1) )
      allocate (Precip_state%snow3d        (idim, jdim, kdim+1) )
      allocate (Precip_state%snowclr3d     (idim, jdim, kdim+1) )
      allocate (Precip_state%surfrain      (idim, jdim) )
      allocate (Precip_state%surfsnow      (idim, jdim) )

      Precip_state%lsc_snow      = 0.
      Precip_state%lsc_rain      = 0.
      Precip_state%lsc_snow_size = 0.
      Precip_state%lsc_rain_size = 0.
      Precip_state%qsout3d_mg    = 0.
      Precip_state%qrout3d_mg    = 0.
      Precip_state%rain3d        = 0.
      Precip_state%snow3d        = 0.
      Precip_state%snowclr3d     = 0.
      Precip_state%surfrain      = 0.
      Precip_state%surfsnow      = 0.

!-----------------------------------------------------------------------
!    allocate and initialize the components of the cloud_processes_type
!    variable Cloud_processes.
!-----------------------------------------------------------------------
      allocate  (Cloud_processes%da_ls         (idim, jdim, kdim) )
      allocate  (Cloud_processes%D_eros        (idim, jdim, kdim) )
      allocate  (Cloud_processes%qvg           (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls      (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls_ice  (idim, jdim, kdim) )
      allocate  (Cloud_processes%dcond_ls_tot  (idim, jdim, kdim) )
      allocate  (Cloud_processes%delta_cf      (idim, jdim, kdim) )
      allocate  (Cloud_processes%f_snow_berg   (idim, jdim, kdim) )

      Cloud_processes%da_ls          = 0.
      Cloud_processes%D_eros         = 0.
      Cloud_processes%qvg            = 0.
      Cloud_processes%dcond_ls       = 0.
      Cloud_processes%dcond_ls_ice   = 0.
      Cloud_processes%dcond_ls_tot   = 0.
      Cloud_processes%delta_cf       = 0.
      Cloud_processes%f_snow_berg    = 0.

!--------------------------------------------------------------------------


end subroutine strat_alloc




!##########################################################################

subroutine strat_dealloc (Atmos_state, Particles, Cloud_State, &
                          Precip_state, Cloud_Processes)
      
!-----------------------------------------------------------------------
type(atmos_state_type),     intent(inout) :: Atmos_state
type(particles_type),       intent(inout) :: Particles   
type(cloud_state_type),     intent(inout) :: Cloud_state
type(precip_state_type),    intent(inout) :: Precip_state
type(cloud_processes_type), intent(inout) :: Cloud_processes


!-----------------------------------------------------------------------
!    deallocate the components of the atmos_state_type variable 
!    Atmos_State.
!-----------------------------------------------------------------------
      deallocate (Atmos_state%pfull       )
      deallocate (Atmos_state%phalf       )
      deallocate (Atmos_state%zfull       )
      deallocate (Atmos_state%zhalf       )
      deallocate (Atmos_state%radturbten2 )
      deallocate (Atmos_state%T_in        )
      deallocate (Atmos_state%qv_in       )
      deallocate (Atmos_state%omega       )
      deallocate (Atmos_state%Mc          )
      deallocate (Atmos_state%diff_t      )
      deallocate (Atmos_state%qrat        )
      deallocate (Atmos_state%ahuco       )
      deallocate (Atmos_state%airdens     )
      deallocate (Atmos_state%tn          )
      deallocate (Atmos_state%qvn         )
      deallocate (Atmos_state%qs          )
      deallocate (Atmos_state%dqsdT       )
      deallocate (Atmos_state%qsi         )
      deallocate (Atmos_state%qsl         )
      deallocate (Atmos_state%rh_crit     )
      deallocate (Atmos_state%rh_crit_min )
      deallocate (Atmos_state%gamma       )
      deallocate (Atmos_state%esat0       )
      deallocate (Atmos_state%U_ca        )
      deallocate (Atmos_state%delp        )
      deallocate (Atmos_state%U01         )
      deallocate (Atmos_state%deltpg      )
 
!-----------------------------------------------------------------------
!    deallocate the components of the particles_type variable Particles.
!-----------------------------------------------------------------------
      deallocate (Particles%concen_dust_sub )
      deallocate (Particles%drop1           )
      deallocate (Particles%drop2           )
      deallocate (Particles%crystal1        )
      deallocate (Particles%rbar_dust       )
      deallocate (Particles%ndust           )
      deallocate (Particles%hom             )

!-----------------------------------------------------------------------
!    deallocate the components of the cloud_state_type variable 
!    Cloud_state.
!-----------------------------------------------------------------------
      deallocate (Cloud_state%ql_upd  )
      deallocate (Cloud_state%qi_upd  )
      deallocate (Cloud_state%qa_upd  )
      deallocate (Cloud_state%qn_upd  )
      deallocate (Cloud_state%qni_upd )
      deallocate (Cloud_state%ql_mean )
      deallocate (Cloud_state%qi_mean )
      deallocate (Cloud_state%qa_mean )
      deallocate (Cloud_state%qn_mean )
      deallocate (Cloud_state%qni_mean)
      deallocate (Cloud_state%ql_in   )
      deallocate (Cloud_state%qi_in   )
      deallocate (Cloud_state%qa_in   )
      deallocate (Cloud_state%qn_in   )
      deallocate (Cloud_state%qni_in  )
      deallocate (Cloud_state%SL_out  )
      deallocate (Cloud_state%SI_out  )
      deallocate (Cloud_state%SA_out  )
      deallocate (Cloud_state%SN_out  )
      deallocate (Cloud_state%SNi_out )
      deallocate (Cloud_state%qa_upd_0)
      deallocate (Cloud_state%SA_0    )

!-----------------------------------------------------------------------
!    deallocate the components of the precip_state_type variable 
!    Precip_state.
!-----------------------------------------------------------------------
      deallocate (Precip_state%lsc_snow      )
      deallocate (Precip_state%lsc_rain      )
      deallocate (Precip_state%lsc_snow_size )
      deallocate (Precip_state%lsc_rain_size )
      deallocate (Precip_state%qsout3d_mg    )
      deallocate (Precip_state%qrout3d_mg    )
      deallocate (Precip_state%rain3d        )
      deallocate (Precip_state%snow3d        )
      deallocate (Precip_state%snowclr3d     )
      deallocate (Precip_state%surfrain      )
      deallocate (Precip_state%surfsnow      )

!-----------------------------------------------------------------------
!    deallocate the components of the cloud_processes_type variable 
!    Cloud_processes.
!-----------------------------------------------------------------------
      deallocate (Cloud_processes%da_ls       )
      deallocate (Cloud_processes%D_eros      )
      deallocate (Cloud_processes%qvg         )
      deallocate (Cloud_processes%dcond_ls    )
      deallocate (Cloud_processes%dcond_ls_ice)
      deallocate (Cloud_processes%dcond_ls_tot)
      deallocate (Cloud_processes%delta_cf    )
      deallocate (Cloud_processes%f_snow_berg )

!-------------------------------------------------------------------------

end subroutine strat_dealloc




!#######################################################################



end module strat_cloud_mod

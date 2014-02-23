                 module strat_clouds_W_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  sak
! </REVIEWER>
! <OVERVIEW>
!    strat_clouds_W_mod obtains the cloud specification variables
!    for the klein strat cloud parameterization from cloud_rad_mod
!    and makes them available to the radiation package.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!   shared modules:

use constants_mod,          only: radian
use time_manager_mod,       only: time_type, time_manager_init
use mpp_mod,                only: input_nml_file
use fms_mod,                only: open_namelist_file, mpp_pe, &
                                  mpp_root_pe, stdlog,  fms_init, &
                                  write_version_number, file_exist, &
                                  check_nml_error, error_mesg,   &
                                  FATAL, close_file

!   shared radiation package modules:

use rad_utilities_mod,      only: rad_utilities_init, &
                                  cldrad_properties_type,  &
                                  cld_specification_type, &
                                  solar_spectrum_type, &
                                  microphysics_type, Cldrad_control
use esfsw_parameters_mod,   only: Solar_spect, esfsw_parameters_init

!   cloud parameterization module:

use cloud_rad_mod,          only: cloud_rad_init, cloud_summary3, &
                                  lw_emissivity, sw_optical_properties, &
                                  snow_and_rain

!    stochastic cloud generator module
use random_numbers_mod,     only: randomNumberStream,           &
                                  initializeRandomNumberStream, &
                                  constructSeed
use cloud_generator_mod,    only: cloud_generator_init, &
                                  generate_stochastic_clouds,&
                                  cloud_generator_end
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    strat_clouds_W_mod obtains the cloud specification variables
!    for the klein strat cloud parameterization from cloud_rad_mod
!    and makes them available to the radiation package.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: strat_clouds_W.F90,v 19.0 2012/01/06 20:24:43 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          strat_clouds_W_init, strat_clouds_amt, obtain_bulk_lw_strat, &
          obtain_bulk_sw_strat, strat_clouds_W_end

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: do_stochastic_clouds = .false.
integer   :: seedperm = 0
logical   :: one_generator_call = .false.


namelist /strat_clouds_W_nml /                      &
                                do_stochastic_clouds, seedperm, &
                                one_generator_call


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


logical                               :: module_is_initialized = .false.  ! module is initialized ?
real,    dimension(:,:),    allocatable :: lats, lons  ! lats and lons in this processor window (degrees)
!----------------------------------------------------------------------
!----------------------------------------------------------------------



                              contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="strat_clouds_W_init">
!  <OVERVIEW>
!    strat_clouds_W_init is the constructor for strat_clouds_W_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    strat_clouds_W_init is the constructor for strat_clouds_W_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_clouds_W_init
!  </TEMPLATE>
! </SUBROUTINE>
!  
subroutine strat_clouds_W_init(latb, lonb)
  real, dimension(:,:), intent( in) :: latb, lonb
!---------------------------------------------------------------------
!    strat_clouds_W_init is the constructor for strat_clouds_W_mod.
!---------------------------------------------------------------------
!       lonb      2d array of model longitudes on cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]


!----------------------------------------------------------------------
!   local variables:

      integer   ::   unit, ierr, io, logunit

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    Save copies of lat and lon values for this window
!---------------------------------------------------------------------
      allocate(lats(size(latb,1),size(latb,2)), lons(size(lonb,1),size(lonb,2)))
      lats(:,:) = latb(:,:) * radian
      lons(:,:) = lonb(:,:) * radian
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call rad_utilities_init
      call esfsw_parameters_init
      call cloud_rad_init

!---------------------------------------------------------------------
!    read namelist.         
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=strat_clouds_W_nml, iostat=io)
   ierr = check_nml_error(io,'strat_clouds_W_nml')
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=strat_clouds_W_nml, iostat=io, end=10) 
        ierr = check_nml_error (io, 'strat_clouds_W_nml')
        enddo                       
10      call close_file (unit)      
      endif                         
#endif
                                    
!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                 write (logunit, nml=strat_clouds_W_nml)

!--------------------------------------------------------------------
!    save the flags indicating whether stochastic clouds are to be
!    used.
!--------------------------------------------------------------------
      Cldrad_control%do_stochastic_clouds = do_stochastic_clouds
      Cldrad_control%do_stochastic_clouds_iz = .true.
      if (do_stochastic_clouds) &
                 call cloud_generator_init

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------

end subroutine strat_clouds_W_init



!######################################################################
! <SUBROUTINE NAME="strat_clouds_amt">
!  <OVERVIEW>
!    strat_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    strat_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. if a microphysically-based cloud parameter-
!    ization is being used, particle sizes and concentrations are also
!    provided.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_clouds_amt (is, ie, js, je, Rad_time, pflux, press, 
!                          temp, qv, land, Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation calculation is to apply
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!   pressure values at flux levels (average of pressure values at
!   model grid points
!  </IN>
!  <IN NAME="press" TYPE="real">
!   pressure values at model grid points. surface 
!                   pressure is stored at index value nlev+1
!  </IN>
!  <IN NAME="temp" TYPE="real">
!    temperature at model levels (1:nlev), to be used
!                   in cloud calculations
!  </IN>
!  <IN NAME="qv" TYPE="real">
!    water vapor specific humidity at model levels (1:nlev), to be used
!                   in cloud calculations
!  </IN>
!  <IN NAME="land" TYPE="real">
!   fraction of grid box covered by land
!  </IN>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphys_type">
!   microphysics_type variable containing the size,
!                   concentration and fraction of the four condensate 
!                   types (cloud drop, cloud ice, rain, snow) in the 
!                   grid box, present when microphysically-based
!                   cloud radiation properties are desired.
!  </INOUT>
! </SUBROUTINE>
!
subroutine strat_clouds_amt (is, ie, js, je, Rad_time, pflux, &
                             press, temp, qv, &
                             land, Cld_spec, Lsc_microphys)

!---------------------------------------------------------------------
!    strat_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. if a microphysically-based cloud parameter-
!    ization is being used, particle sizes and concentrations are also
!    provided.
!----------------------------------------------------------------------

integer,                      intent(in)        :: is, ie, js, je
type(time_type),              intent(in)        :: Rad_time
real,    dimension(:,:,:),    intent(in)        :: pflux, press, temp, qv
real,    dimension(:,:),      intent(in)        :: land
type(cld_specification_type), intent(inout)     :: Cld_spec      
type(microphysics_type),      intent(inout)     :: Lsc_microphys

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Rad_time     time type variable containing radiation time
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ] 
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      temp         temperature at model levels (1:nlev), to be used
!                   in cloud calculations
!                   [ deg K ]
!      qv           water vapor specific humidity at model levels
!                   (1:nlev), to be used in cloud calculations
!      land         fraction of grid box covered by land
!                   [ non-dimensional ]
!
!   intent(inout), optional variables:
!
!      Cld_spec     cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!
!               the following elements of Cld_spec are defined here:
!
!                  %cmxolw  fraction of maximally overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %crndlw  fraction of randomly overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %camtsw  cloud fraction seen by the shortwave
!                           radiation; the sum of the maximally
!                           overlapped and randomly overlapped 
!                           longwave cloud fractions  [ dimensionless ]
!                  %nmxolw  number of maximally overlapped longwave 
!                           clouds in each grid column.
!                  %nrndlw  number of randomly overlapped longwave 
!                           clouds in each grid column.
!                  %ncldsw  number of clouds seen by the shortwave
!                           radiation in each grid column.
!                  %cloud_thickness
!                           number of model layers over which the cloud
!                           in this grid box extends
!                  %lwp     liquid water path 
!                           [ kg / m^2 ]
!                  %iwp     ice water path
!                           [ kg / m^2 ]
!                  %reff_liq
!                           effective drop radius [ microns ]
!                  %reff_ice
!                           effective ice particle size [ microns ]
!
!      Lsc_microphys
!                   microphysics_type variable containing the size,
!                   concentration and fraction of the four condensate 
!                   types (cloud drop, cloud ice, rain, snow) in the 
!                   grid box, present when microphysically-based
!                   cloud radiation properties are desired.
!
!               the following components of this variable are output 
!               from this routine when microphysically-based properties
!               are desired:
!
!                  %conc_ice  ice particle concentration [ g /m^3 ]
!                  %conc_drop cloud droplet concentration [ g /m^3 ]
!                  %size_ice  ice particle effective diameter 
!                  [ microns ]
!                  %size_drop cloud droplet effective diameter
!                  [ microns ]
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    local variables

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1) ::      cldamt

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1, Cldrad_control%nlwcldb) :: &
                         ql_stoch_lw2, qi_stoch_lw2, qa_stoch_lw2, &
                         qn_stoch_lw2, qni_stoch_lw2

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1, Solar_spect%nbands) :: &
                         ql_stoch_sw2, qi_stoch_sw2, qa_stoch_sw2, &
                         qn_stoch_sw2, qni_stoch_sw2

      real, dimension (size(pflux,1), size(pflux,2),                 &
                       size(pflux,3)-1,                              &
                       Cldrad_control%nlwcldb + Solar_spect%nbands), &
             target :: ql_stoch, qi_stoch, qa_stoch, qn_stoch, qni_stoch
 
      real, dimension(:, :, :, :), pointer :: &
                  ql_stoch_lw, qi_stoch_lw, qa_stoch_lw, qn_stoch_lw, &
                  ql_stoch_sw, qi_stoch_sw, qa_stoch_sw, qn_stoch_sw, &
                  qni_stoch_lw, qni_stoch_sw
      
!      integer, dimension(size(Cld_spec%cld_thickness_lw_band, 1), &
!                         size(Cld_spec%cld_thickness_lw_band, 2), &
!                         size(Cld_spec%cld_thickness_lw_band, 3), &
      integer, dimension(size(temp, 1), size(temp, 2), size(temp, 3), &
          Cldrad_control%nlwcldb + Solar_spect%nbands) ::         &
                        cld_thickness

      integer, dimension (size(pflux,1), size(pflux,2), &
                          size(pflux,3)-1) ::   ktop, kbtm

      integer, dimension (size(pflux,1), size(pflux,2)) :: &
                                                 ncldlvls

      type(randomNumberStream), &
                dimension(size(pflux,1), size(pflux,2)) :: streams
      integer     ::    kx 
      integer     ::    i, j, k, kc, nb
      real        ::   seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)

!-------------------------------------------------------------------
!    local variables:
!
!       cldamt          cloud fraction, in cloud-space when microphysics
!                       not being used, in model-space when microphysics
!                       is active 
!                       [ dimensionless ]
!       lwp             cloud liquid water path in cloud-space 
!                       [ kg condensate / m^2 ]
!       iwp             cloud ice path, in cloud-space 
!                       [ kg condensate / m^2 ]
!       reff_liq        effective radius for liquid clouds, 
!                       in cloud-space [ microns ]
!       reff_ice        effective particle size for ice clouds 
!                       in cloud-space [ microns ]
!       ktop            index of the model level which is cloud top
!       kbtm            index of the model level which is cloud base
!       ncldlvls        number of layers with cloud in a column
!       kx              number of model layers
!       i,j,k,kc        do-loop indices  
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kx = size (press,3) - 1

!----------------------------------------------------------------------
!    compute the cloud specification properties under the assumption
!    of random cloud overlap.
!----------------------------------------------------------------------
      if (Cldrad_control%do_random_overlap) then

!----------------------------------------------------------------------
!    if microphysically-based radiative properties are needed, call
!    cloud_summary3 with the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        if (Cldrad_control%do_pred_cld_microphys .or. &
            Cldrad_control%do_presc_cld_microphys) then

!--------------------------------------------------------------------
!    call cloud_summary3 with the full cloud field, regardless of 
!    whether or not stochastic clouds are active.
!    the full cloud field is assumed random overlap when stochastic
!    clouds are activated.
!--------------------------------------------------------------------
            where (Cld_spec%cloud_area(:,:,:) > 0.0) 
              Cld_spec%cld_thickness(:,:,:) = 1
            end where
          call cloud_summary3 (is, js, land,   &
                               Cldrad_control%using_fu2007, &
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &       
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, ncldlvls, &
                               cldamt, Cld_spec%lwp, Cld_spec%iwp,   &
                               Cld_spec%reff_liq, Cld_spec%reff_ice, &
                               conc_drop= Lsc_microphys%conc_drop, &
                               conc_ice = Lsc_microphys%conc_ice, &
                               size_drop =Lsc_microphys%size_drop,   &
                               size_ice = Lsc_microphys%size_ice, &
                         droplet_number = Lsc_microphys%droplet_number, &
                            ice_number = Lsc_microphys%ice_number )
 

          call snow_and_rain(Cld_spec%cloud_area, press(:,:,1:kx),  &
                             pflux, temp,  cldamt, Cld_spec%snow,  &
                             Cld_spec%rain, Cld_spec%snow_size,  &
                             Cld_spec%rain_size,   &
                             Lsc_microphys%conc_rain,   &
                             Lsc_microphys%conc_snow,   &
                             Lsc_microphys%size_rain,   &
                             Lsc_microphys%size_snow )


          cldamt = MIN (cldamt, 1.0)
          if (.not. Cldrad_control%do_specified_strat_clouds) then
            Cld_spec%ncldsw        = ncldlvls
            Cld_spec%nrndlw        = ncldlvls         
            Cld_spec%camtsw        = cldamt           
            Cld_spec%crndlw        = cldamt
            Lsc_microphys%cldamt   = cldamt            
          endif

!---------------------------------------------------------------------
!    if using stochastic clouds for either sw or lw, Initialize the random number streams, 
!       one per grid cell, with unique and replicable integer based 
!       on grid location and model date/time
!---------------------------------------------------------------------
          if (do_stochastic_clouds) then
            if (Cldrad_control%use_temp_for_seed) then
              do j = 1, size(Cld_spec%cloud_water, 2)
                do i = 1, size(Cld_spec%cloud_water, 1)
                  streams(i, j) =   &
                           initializeRandomNumberStream(  &
                               ishftc(nint(temp(i,j,1)*seedwts),seedperm))

                end do
              end do
             else
            do j = 1, size(Cld_spec%cloud_water, 2)
              do i = 1, size(Cld_spec%cloud_water, 1)
                streams(i, j) =   &
                         initializeRandomNumberStream(  &
                            constructSeed(nint(lons(is + i - 1, js + j - 1)), &
                                          nint(lats(is + i - 1, js + j - 1)), &
                                                    Rad_time, seedperm))
              end do
            end do
           endif

            if (one_generator_call) then
!---------------------------------------------------------------------
!    then generate all the subcolumns at once and divide them into 
!    those needed  for the sw and lw bands.
!    call routine to obtain  band-dependent values of ql, qi and qa. 
!---------------------------------------------------------------------
              call generate_stochastic_clouds (        &
                      streams,                &
                      Cld_spec%cloud_water,   &
                      Cld_spec%cloud_ice,     &
                      Cld_spec%cloud_area,    &
                      Cld_spec%cloud_droplet, &
                      Cld_spec%cloud_ice_num, &
                      pFull = press(:, :, :kx),&
                      pHalf = pflux, &
                      temperature = temp(:, :, :kx),        &
                      qv= qv(:, :, :kx), &
                      cld_thickness = cld_thickness, &
                      ql_stoch = ql_stoch, &
                      qi_stoch = qi_stoch, &
                      qa_stoch = qa_stoch, &
                      qn_stoch = qn_stoch, &
                      qni_stoch = qni_stoch )

          ql_stoch_lw => ql_stoch(:, :, :, 1:Cldrad_control%nlwcldb)
          qi_stoch_lw => qi_stoch(:, :, :, 1:Cldrad_control%nlwcldb)
          qa_stoch_lw => qa_stoch(:, :, :, 1:Cldrad_control%nlwcldb)
          qn_stoch_lw => qn_stoch(:, :, :, 1:Cldrad_control%nlwcldb)
          qni_stoch_lw => qni_stoch(:, :, :, 1:Cldrad_control%nlwcldb)
          Cld_spec%cld_thickness_lw_band = &
                       cld_thickness(:, :, :, 1:Cldrad_control%nlwcldb)

          ql_stoch_sw => ql_stoch(:, :, :, Cldrad_control%nlwcldb +1:)
          qi_stoch_sw => qi_stoch(:, :, :, Cldrad_control%nlwcldb +1:)
          qa_stoch_sw => qa_stoch(:, :, :, Cldrad_control%nlwcldb +1:)
          qn_stoch_sw => qn_stoch(:, :, :, Cldrad_control%nlwcldb +1:)
          qni_stoch_sw => qni_stoch(:, :, :, Cldrad_control%nlwcldb +1:)
          Cld_spec%cld_thickness_sw_band = &
                     cld_thickness(:, :, :, Cldrad_control%nlwcldb +1:)

!---------------------------------------------------------------------
!    call cloud_summary3 for each lw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,Cldrad_control%nlwcldb
              call cloud_summary3 (          &
                is, js, land, &
                Cldrad_control%using_fu2007, &
                ql_stoch_lw(:,:,:,nb),&
                qi_stoch_lw(:,:,:,nb), qa_stoch_lw(:,:,:,nb),&
                qn_stoch_lw(:,:,:,nb), &
                qni_stoch_lw(:,:,:,nb),&
                press(:,:,1:kx), pflux, temp, ncldlvls, &
                cldamt, Cld_spec%lwp_lw_band(:,:,:,nb),&
                Cld_spec%iwp_lw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_lw_band(:,:,:,nb), &
                Cld_spec%reff_ice_lw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%lw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%lw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%lw_stoch_droplet_number(:,:,:,nb), &
          ice_number =  Lsc_microphys%lw_stoch_ice_number(:,:,:,nb) )
              
              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column. 
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)     
              Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.                  
               
              if (.not. Cldrad_control%do_specified_strat_clouds) then
                Cld_spec%nrndlw_band(:,:,nb) = ncldlvls(:,:)
                Lsc_microphys%lw_stoch_cldamt(:,:,:,nb) = cldamt
              endif
            end do

!---------------------------------------------------------------------
!    call cloud_summary3 for each sw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
          do nb=1,size(Lsc_microphys%sw_stoch_conc_ice,4)
            call cloud_summary3 (                            &
              is, js, land,  &
              Cldrad_control%using_fu2007, &
              ql_stoch_sw(:,:,:,nb), &
              qi_stoch_sw(:,:,:,nb), qa_stoch_sw(:,:,:,nb),&
              qn_stoch_sw(:,:,:,nb), &
              qni_stoch_sw(:,:,:,nb), &
              press(:,:,1:kx), pflux, temp, ncldlvls, &
              cldamt, Cld_spec%lwp_sw_band(:,:,:,nb), &
              Cld_spec%iwp_sw_band(:,:,:,nb),   &
              Cld_spec%reff_liq_sw_band(:,:,:,nb), &
              Cld_spec%reff_ice_sw_band(:,:,:,nb), &
              conc_drop= Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb), &
              conc_ice = Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb), &
              size_drop =Lsc_microphys%sw_stoch_size_drop(:,:,:,nb), &
              size_ice = Lsc_microphys%sw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%sw_stoch_droplet_number(:,:,:,nb), &
         ice_number = Lsc_microphys%sw_stoch_ice_number(:,:,:,nb) )

         !now that the vertical cloud fraction has been used to
         !properly calculate the in-cloud particle size, rescale
         !the concentrations and cloud amounts to that the cloud
         !amount is unity in any partially cloudy sub-column. 
         !
         !This is necessary so that the radiation code will not
         !do cloud fraction weights of cloudy and clear sky fluxes.
         !
         !The rescaling of the concentrations is necessary so that         the
         !total optical depth of the layer is constant.  Note that         this
         !works because cloud extinction is linear in the concentra        tion
            Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) = &
             Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)
            Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) = &
             Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
             where (cldamt .gt. 0.) cldamt = 1.

             if (.not. Cldrad_control%do_specified_strat_clouds) then
               Cld_spec%ncldsw_band(:,:,nb) = ncldlvls(:,:)
               Lsc_microphys%sw_stoch_cldamt(:,:,:,nb) = cldamt
             endif
           end do
 
        else  ! (one_call)
!---------------------------------------------------------------------
!    call routine to obtain lw band-dependent values of ql, qi and qa. 
!---------------------------------------------------------------------
            call generate_stochastic_clouds (        &
                     streams,                &
                     Cld_spec%cloud_water,   &
                     Cld_spec%cloud_ice,     &
                     Cld_spec%cloud_area,    &
                     Cld_spec%cloud_droplet, &
                     Cld_spec%cloud_ice_num, &
                     pFull    = press(:, :, :kx),     &
                     pHalf    = pflux,&
                     temperature = temp(:, :, :kx),   &
                     qv = qv(:,:, :kx),   &
                     cld_thickness = Cld_spec%cld_thickness_lw_band, &
                     ql_stoch = ql_stoch_lw2, &
                     qi_stoch = qi_stoch_lw2, &
                     qa_stoch = qa_stoch_lw2, &
                     qn_stoch = qn_stoch_lw2, &
                     qni_stoch = qni_stoch_lw2 ) 

!---------------------------------------------------------------------
!    call routine to obtain sw band-dependent values of ql, qi and qa. 
!---------------------------------------------------------------------
            call generate_stochastic_clouds (        &
                     streams,                &
                     Cld_spec%cloud_water,   &
                     Cld_spec%cloud_ice,     &
                     Cld_spec%cloud_area,    &
                     Cld_spec%cloud_droplet, &
                     Cld_spec%cloud_ice_num, &
                     pFull    = press(:, :, :kx),     &
                     pHalf    = pflux,&
                     temperature = temp(:, :, :kx),   &
                     qv = qv(:,:, :kx),   &
                     cld_thickness = Cld_spec%cld_thickness_sw_band, &
                     ql_stoch = ql_stoch_sw2, &
                     qi_stoch = qi_stoch_sw2, &
                     qa_stoch = qa_stoch_sw2, &
                     qn_stoch = qn_stoch_sw2, &
                     qni_stoch = qni_stoch_sw2 )

!---------------------------------------------------------------------
!    call cloud_summary3 for each lw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,Cldrad_control%nlwcldb
              call cloud_summary3 (          &
                is, js, land,  &
                Cldrad_control%using_fu2007, &
                ql_stoch_lw2(:,:,:,nb),&
                qi_stoch_lw2(:,:,:,nb), qa_stoch_lw2(:,:,:,nb),&
                qn_stoch_lw2(:,:,:,nb), &
                qni_stoch_lw2(:,:,:,nb), &
                press(:,:,1:kx), pflux, temp, ncldlvls, &
                cldamt, Cld_spec%lwp_lw_band(:,:,:,nb),&
                Cld_spec%iwp_lw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_lw_band(:,:,:,nb), &
                Cld_spec%reff_ice_lw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%lw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%lw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%lw_stoch_droplet_number(:,:,:,nb), &
         ice_number =   Lsc_microphys%lw_stoch_ice_number(:,:,:,nb) )
              
              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column. 
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)     
              Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.                  
               
              if (.not. Cldrad_control%do_specified_strat_clouds) then
                Cld_spec%nrndlw_band(:,:,nb) = ncldlvls(:,:)
                Lsc_microphys%lw_stoch_cldamt(:,:,:,nb) = cldamt
              endif
            end do

!---------------------------------------------------------------------
!    call cloud_summary3 for each sw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,size(Lsc_microphys%sw_stoch_conc_ice,4)
              call cloud_summary3 (                            &
                is, js, land, &
                Cldrad_control%using_fu2007, &
                ql_stoch_sw2(:,:,:,nb), &
                qi_stoch_sw2(:,:,:,nb), qa_stoch_sw2(:,:,:,nb),&
                qn_stoch_sw2(:,:,:,nb), &
                qni_stoch_sw2(:,:,:,nb), &
                press(:,:,1:kx), pflux, temp, ncldlvls, &
                cldamt, Cld_spec%lwp_sw_band(:,:,:,nb), &
                Cld_spec%iwp_sw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_sw_band(:,:,:,nb), &
                Cld_spec%reff_ice_sw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%sw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%sw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%sw_stoch_droplet_number(:,:,:,nb), &
         ice_number =  Lsc_microphys%sw_stoch_ice_number(:,:,:,nb))
              
              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column. 
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)     
              Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.                  
               
              if (.not. Cldrad_control%do_specified_strat_clouds) then
                Cld_spec%ncldsw_band(:,:,nb) = ncldlvls(:,:)
                Lsc_microphys%sw_stoch_cldamt(:,:,:,nb) = cldamt
              endif
            end do
         endif ! (one_generator_call)
       endif  ! (do_stochastic_clouds)

!---------------------------------------------------------------------
!    if microphysically-based radiative properties are not needed, call
!    cloud_summary3 without the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        else  ! (not micro)

!--------------------------------------------------------------------
!    define the cloud thickness to be 1 at those points with cloud
!    present.
!---------------------------------------------------------------------
          where (Cld_spec%cloud_area(:,:,:) > 0.0) 
            Cld_spec%cld_thickness(:,:,:) = 1
          end where
          call cloud_summary3 (is, js, land,  &
                               Cldrad_control%using_fu2007, &      
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, ncldlvls, &
                               cldamt, Cld_spec%lwp,   &
                               Cld_spec%iwp, Cld_spec%reff_liq,   &
                               Cld_spec%reff_ice)
             cldamt = MIN (cldamt, 1.0)
          if (.not. Cldrad_control%do_specified_strat_clouds) then
            Cld_spec%ncldsw        = ncldlvls
            Cld_spec%nrndlw        = ncldlvls         
            Cld_spec%camtsw        = cldamt           
            Cld_spec%crndlw        = cldamt
            Lsc_microphys%cldamt   = cldamt            
          endif
        endif

!----------------------------------------------------------------------
!    in gcm or sa_gcm mode, all clouds are assumed to be randomly 
!    overlapped. in columns mode, no change is made to the specified 
!    input characteristics.
!----------------------------------------------------------------------
        if (.not. Cldrad_control%do_specified_strat_clouds) then
          Cld_spec%nmxolw        = 0 
          Cld_spec%cmxolw        = 0.0E+00
        endif

!---------------------------------------------------------------------
!    define cloud specification properties when max-random overlap is
!    assumed. in this case cloud in adjacent layers is assumed to be
!    part of the same cloud.
!---------------------------------------------------------------------
      else if (Cldrad_control%do_max_random_overlap) then  

!----------------------------------------------------------------------
!    microphysically-based radiative properties are not implemented
!    with the max-random overlap assumption.
!----------------------------------------------------------------------
        if (Cldrad_control%do_pred_cld_microphys .or. &
            Cldrad_control%do_presc_cld_microphys) then
          call error_mesg ('strat_clouds_W_mod', &
               'must use random overlap cloud assumption with strat '//&
              'clouds when microphysics are desired', FATAL)

!---------------------------------------------------------------------
!    if microphysically-based radiative properties are not needed, call
!    cloud_summary3 without the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        else
          call cloud_summary3 (is, js, land,  &
                               Cldrad_control%using_fu2007, &    
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, ncldlvls, &
                               Cld_spec%camtsw, Cld_spec%lwp,   &
                               Cld_spec%iwp, Cld_spec%reff_liq,   &
                               Cld_spec%reff_ice,  &
                               ktop=ktop, kbot=kbtm)

!---------------------------------------------------------------------
!    when only bulk properties are returned, they are in cloud space,
!    and must be converted to physical space before being stored in 
!    Cld_spec. random overlap and max overlap properties are assigned 
!    according to the cloud thickness - multi layer clouds are assumed 
!    to be max overlap.
!-------------------------------------------------------------------
          do j=1, size(press,2)
            do i=1, size(press,1)
              Cld_spec%ncldsw(i,j) = ncldlvls(i,j)
              do kc=1, Cld_spec%ncldsw(i,j)         
                do k=ktop(i,j,kc), kbtm(i,j,kc)
                  if (ktop(i,j,kc) == kbtm(i,j,kc)) then
                    Cld_spec%crndlw(i,j,k) = Cld_spec%camtsw(i,j,k)
                    Cld_spec%cmxolw(i,j,k) = 0.0             
                    Cld_spec%cld_thickness(i,j,k) = 1
                  else
                    Cld_spec%cmxolw(i,j,k) = Cld_spec%camtsw(i,j,k)
                    Cld_spec%crndlw(i,j,k) = 0.0
                    Cld_spec%cld_thickness(i,j,k) = kbtm(i,j,kc) -    &
                                                    ktop(i,j,kc) + 1
                  endif
                end do
                if (ktop(i,j,kc) == kbtm(i,j,kc)) then
                  Cld_spec%nrndlw(i,j) = Cld_spec%nrndlw(i,j) + 1
                else
                  Cld_spec%nmxolw(i,j) = Cld_spec%nmxolw(i,j) + 1
                endif
              end do
            end do
          end do
        endif ! (do_pred_micro or do_presc_micro)
      endif ! (do_random_overlap)

!---------------------------------------------------------------------



end subroutine strat_clouds_amt  



!#####################################################################
! <SUBROUTINE NAME="obtain_bulk_lw_strat">
!  <OVERVIEW>
!   obtain_bulk_lw_strat defines bulk longwave cloud radiative 
!    properties for the klein strat cloud scheme. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_bulk_lw_strat defines bulk longwave cloud radiative 
!    properties for the klein strat cloud scheme.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_strat (is, ie, js, je, Cld_spec, Cldrad_props)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!  </IN>
!  <INOUT NAME="cldrad_properties" TYPE="microphys_type">
!   cloud radiative properties on model grid
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_bulk_lw_strat (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_strat defines bulk longwave cloud radiative 
!    properties for the klein strat cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(in)    :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Cld_spec     cloud specification arrays defining the 
!                   location, water paths and effective particle
!                   sizes of clouds that are present, provides 
!                   input to this subroutine
!                   [ cld_specification_type ]
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
!
!               the following components of this variable are output 
!               from this routine:
!
!                    %emrndlw   longwave cloud emissivity for 
!                               randomly overlapped clouds
!                               in each of the longwave
!                               frequency bands  [ dimensionless ]
!                    %emmxolw   longwave cloud emissivity for 
!                               maximally overlapped clouds
!                               in each of the longwave 
!                               frequency bands  [ dimensionless ]
!
!---------------------------------------------------------------------
 
!-------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2),  &
                       size(Cld_spec%lwp,3)) :: emcld

      integer       :: max_cld
      integer       :: i, j, k

!-------------------------------------------------------------------
!   local variables:
!
!         emcld      longwave cloud emissivity [ dimensionless ]
!         max_cld    maximum number of clouds in any column in the
!                    window
!         i,j,k      do-loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!   find maximum number of clouds in any column in the window.
!---------------------------------------------------------------------
      max_cld = MAXVAL(Cld_spec%ncldsw(:,:))

!---------------------------------------------------------------------
!    if cloud is present in the window, call lw_emissivity to compute 
!    the longwave emissivity. otherwise, leave the emissivity arrays 
!    with their previously initialized values.
!---------------------------------------------------------------------
      if (max_cld > 0) then
!---------------------------------------------------------------------
!    call lw_emissivity to obtain the longwave cloud emissivity.
!---------------------------------------------------------------------
        call lw_emissivity (is, js, Cld_spec%lwp, Cld_spec%iwp,  &
                         Cld_spec%reff_liq_lim, Cld_spec%reff_ice_lim,&
                            Cld_spec%ncldsw, emcld)

!---------------------------------------------------------------------
!    define both the random and max overlap cloud emissivities to be
!    that value returned from lw_emissivity.
!-------------------------------------------------------------------
        do k=1,size(Cld_spec%lwp,3)
          do j=1,size(Cld_spec%lwp,2)
            do i=1,size(Cld_spec%lwp,1)
              Cldrad_props%emrndlw(i,j,k,:,1) = emcld(i,j,k)
              Cldrad_props%emmxolw(i,j,k,:,1) = emcld(i,j,k)
            end do
          end do
        end do
      endif

!--------------------------------------------------------------------


end subroutine obtain_bulk_lw_strat



!#####################################################################
! <SUBROUTINE NAME="obtain_bulk_sw_strat">
!  <OVERVIEW>
!   obtain_bulk_lw_strat defines bulk shortwave cloud radiative 
!    properties for the klein strat cloud scheme. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   obtain_bulk_lw_strat defines bulk shortwave cloud radiative 
!    properties for the klein strat cloud scheme.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_strat (is, ie, js, je, cosz, Cld_spec, Cldrad_props)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="cosz" TYPE="real">
!   cosine of the solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!  </IN>
!  <INOUT NAME="cldrad_properties" TYPE="microphys_type">
!   cloud radiative properties on model grid
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_bulk_sw_strat (is, ie, js, je, cosz, Cld_spec,   &
                                 Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_strat defines bulk shortwave cloud radiative 
!    properties for the klein strat cloud scheme.
!---------------------------------------------------------------------
 
integer,                      intent(in)    ::  is, ie, js, je
real, dimension(:,:),         intent(in)    ::  cosz
type(cld_specification_type), intent(in)    ::  Cld_spec
type(cldrad_properties_type), intent(inout) ::  Cldrad_props
 
!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle  [ dimensionless ]
!      Cld_spec     cloud specification arrays defining the 
!                   location, amount and type (hi, middle, lo)
!                   of clouds that are present, provides input 
!                   to this subroutine
!                   [ cld_specification_type ]
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
!
!               the following components of this variable are output 
!               from this routine:
!
!                    %cirabsw   absorptivity of clouds in the 
!                               infrared frequency band
!                               [ dimensionless ]
!                    %cirrfsw   reflectivity of clouds in the 
!                               infrared frequency band
!                               [ dimensionless ]
!                    %cvisrfsw  reflectivity of clouds in the 
!                               visible frequency band
!                               [ dimensionless ]
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:

      integer   :: max_cld

!---------------------------------------------------------------------  
!   local variables:
!
!          max_cld    maximum number of clouds in any column in the
!                     window
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!   find maximum number of clouds in any column in the window.
!---------------------------------------------------------------------
      max_cld  = MAXVAL(Cld_spec%ncldsw(:,:))

!---------------------------------------------------------------------
!    if cloud is present in the window, call sw_optical_properties to 
!    compute cloud optical properties and then radiative properties. 
!    otherwise, leave the arrays with their previously initialized 
!    values.
!---------------------------------------------------------------------
      if (max_cld > 0) then
        call sw_optical_properties (Cld_spec%ncldsw, Cld_spec%lwp,   &
                                  Cld_spec%iwp, Cld_spec%reff_liq_lim, &
                                    Cld_spec%reff_ice_lim, cosz,   &
                                    Cldrad_props%cvisrfsw, &
                                    Cldrad_props%cirrfsw,  &
                                    Cldrad_props%cirabsw)
      endif   

!-------------------------------------------------------------------


end subroutine obtain_bulk_sw_strat

!####################################################################
! <SUBROUTINE NAME="strat_clouds_W_end">
!  <OVERVIEW>
!    strat_clouds_W_end is the destructor for strat_clouds_W_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    strat_clouds_W_end is the destructor for strat_clouds_W_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
subroutine strat_clouds_W_end
       
!----------------------------------------------------------------------
!    strat_clouds_W_end is the destructor for strat_clouds_W_mod.
!----------------------------------------------------------------------
        
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
! close cloud_generator
!---------------------------------------------------------------------
      deallocate (lats, lons)
      if (do_stochastic_clouds) &
                call cloud_generator_end()

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------
 
 
end subroutine strat_clouds_W_end



!#################################################################




                    end module strat_clouds_W_mod


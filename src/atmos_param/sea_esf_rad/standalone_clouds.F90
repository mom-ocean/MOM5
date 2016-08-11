!FDOC_TAG_GFDL
 
                 module standalone_clouds_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! fil   
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use mpp_mod,                    only: input_nml_file
use fms_mod,                    only: fms_init, open_namelist_file, &
                                      write_version_number, mpp_pe, &
                                      mpp_root_pe, stdlog,   &
                                      file_exist, check_nml_error,   &
                                      error_mesg, FATAL, close_file
use rad_utilities_mod,          only: rad_utilities_init, &
                                      Cldrad_control, &
                                      cld_specification_type, &
                                      cldrad_properties_type,  &
                                      microphysics_type,  &
                                      microrad_properties_type, &
                                      Sw_control

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!   standalone cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: standalone_clouds.F90,v 19.0 2012/01/06 20:24:41 fms Exp $'
  character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          standalone_clouds_init,                           &
          standalone_clouds_end,                           &
          define_column_properties, &
      standalone_clouds_amt, obtain_micro_lw_sa, obtain_micro_sw_sa,  &
       obtain_bulk_lw_sa, obtain_bulk_sw_sa

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16) :: cldht_type_form='   '
character(len=16) :: cloud_data_form='   '
character(len=16) :: cloud_overlap_form='   '
character(len=16) :: lhsw_cld_prop_form='   '
character(len=16) :: lw_cld_prop_form='   '
integer           :: cloud_data_points = 0

namelist /standalone_clouds_nml /     &
                          cldht_type_form, cloud_data_form, &
                          cloud_data_points, cloud_overlap_form, &
                          lhsw_cld_prop_form,  &
                          lw_cld_prop_form


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!logical ::          do_esfsw

!------------------------------------------------------------------
!    test (singlecolumn) values for cloud amount, height index,
!    infrared emissivity and shortwave reflectivity and absorptivity
!------------------------------------------------------------------

integer                               :: ich, icm, ict, icb
real                                  :: ch, cm, cl
real, dimension(:,:), allocatable     :: cloud_amount_in
real, dimension(:,:), allocatable     :: max_cloud_amount_in
real, dimension(:,:), allocatable     :: rnd_cloud_amount_in
integer, dimension(:), allocatable    :: nmax_cloud_in
integer, dimension(:), allocatable    :: nrnd_cloud_in




 real, dimension(:,:), allocatable     :: cvis_rf_in, cir_rf_in,  &
                                          cir_abs_in

 real, dimension(:,:,:), allocatable   :: emlw_band_in
real, dimension(:,:), allocatable     :: emlw_in

!----------------------------------------------------------------------
!   define default values for shortwave cloud absorptivity and
!   reflectivity for use only in the lacis-hansen implementation
!   of shortwave radiative transfer
!      (these are the values previously used in SKYHI)
!----------------------------------------------------------------------
real                 :: lowcloud_refl_visband = 0.66E+00
real                 :: midcloud_refl_visband = 0.54E+00
real                 :: highcloud_refl_visband = 0.21E+00
real                 :: lowcloud_refl_nearirband = 0.50E+00
real                 :: midcloud_refl_nearirband = 0.46E+00
real                 :: highcloud_refl_nearirband = 0.19E+00
real                 :: lowcloud_abs_visband = 0.0E+00
real                 :: midcloud_abs_visband = 0.0E+00
real                 :: highcloud_abs_visband = 0.0E+00
real                 :: lowcloud_abs_nearirband = 0.30E+00
real                 :: midcloud_abs_nearirband = 0.20E+00
real                 :: highcloud_abs_nearirband = 0.04E+00

!----------------------------------------------------------------------
!   define default (grey) values for longwave emissivity for low, mid,
!   high clouds. these are used if no microphysics parameterizations
!   (or assumptions) are used in the longwave radiative transfer
!      (these are the values previously used in SKYHI)
!----------------------------------------------------------------------
real                 :: lowcloud_emiss = 1.00E+00
real                 :: midcloud_emiss = 1.00E+00
real                 :: highcloud_emiss = 1.00E+00

real  ::       pie 

!---------------------------------------------------------------------
!        cldhp             sigma value defining high-middle cloud boun-
!                          dary at the pole [ dimensionless ]
!        cldhe             sigma value defining high-middle cloud boun-
!                          dary at the equator [ dimensionless ]
!        cldmp             sigma value defining middle-low cloud boun-
!                          dary at the pole [ dimensionless ]
!        cldme             sigma value defining middle-low cloud boun-
!                          dary at the equator [ dimensionless ]
!---------------------------------------------------------------------
      real        :: cldhp = 0.7                   
      real        :: cldhe = 0.4            
      real        :: cldmp = 0.85                  
      real        :: cldme = 0.7                 

     logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------



                        contains 


!####################################################################

! <SUBROUTINE NAME="standalone_clouds_init">
!  <OVERVIEW>
!    subroutine standalone_clouds_init is the constructor for the
!    standalone_clouds_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine standalone_clouds_init is the constructor for the
!    standalone_clouds_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call standalone_clouds_init (pref, lonb, latb)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!       lonb      2d array of model longitudes at cell corners [ radians ]
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
!       latb      2d array of model latitudes at cell corners [radians]
! 
!  </IN>
! </SUBROUTINE>
!
subroutine standalone_clouds_init (pref, lonb, latb)

!--------------------------------------------------------------------
!    subroutine standalone_clouds_init is the constructor for the
!    standalone_clouds_mod.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in)    ::  pref        
real, dimension(:,:), intent(in)    ::  lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes at cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real               :: max_cld_calc
      integer            :: unit, ierr, io, logunit
      integer            :: ktop, kbot
      integer            :: idf, jdf, kx
      integer            :: i, k, kk

!---------------------------------------------------------------------
!   local variables:
!
!      max_cld_calc maximum cloud amount in any layer within a maximum
!                   random overlap cloud [ dimensionless ]
!      unit         io unit for reading nml file and writing logfile
!      io           error status returned from io operation  
!      ierr         error code
!      ktop         cloud top index
!      kbot         cloud bottom index
!      idf          number of longitudes assigned to the processor
!      jdf          number of latitudes assigned to the processor
!      kx           number of layers in the model
!      i,k,kk       do loop indices
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
 
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=standalone_clouds_nml, iostat=io)
   ierr = check_nml_error(io,"standalone_clouds_nml")
#else
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=standalone_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'standalone_clouds_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                       write (logunit, nml=standalone_clouds_nml)

!---------------------------------------------------------------------
!    define size of processor's domain.
!---------------------------------------------------------------------
      idf = size(lonb,1) - 1
      jdf = size(latb,2) - 1
      kx  = size(pref,1) - 1
      
!---------------------------------------------------------------------
!    define the value of pi for later use.
!---------------------------------------------------------------------
      pie = 4.0*ATAN(1.0)

!---------------------------------------------------------------------
!    save the number of requested cloud columns in a derived type 
!    variable for later use. ensure that this quantity is  at least one
!    but no more than the number of longitudes assigned to the processor
!    (idf). 
!---------------------------------------------------------------------
      Cldrad_control%cloud_data_points = cloud_data_points
      if (cloud_data_points < 1 .or.   &
          cloud_data_points > idf) then
        call error_mesg( 'standalone_clouds_mod',  &
                 ' cloud_data_points must be greater than zero but '//&
                   'no larger than the number of model longitudes', &
                                                               FATAL)
      endif

!---------------------------------------------------------------------
!    allocate the module variables which will hold the cloud specific-
!    ation information (total amount, max rnd amount, rnd amount, max
!    rnd cloud number, rnd cloud number). initialize the values to zero,
!    corresponding to the absence of clouds.
!---------------------------------------------------------------------
      allocate (cloud_amount_in     (idf, kx))
      allocate (max_cloud_amount_in (idf, kx))
      allocate (rnd_cloud_amount_in (idf, kx))
      allocate (nmax_cloud_in       (idf))
      allocate (nrnd_cloud_in       (idf))
      cloud_amount_in     = 0.0E+00
      max_cloud_amount_in = 0.0E+00
      rnd_cloud_amount_in = 0.0E+00
      nmax_cloud_in       = 0
      nrnd_cloud_in       = 0

!---------------------------------------------------------------------
!    if cloud data is to be supplied via input file, open the input
!    file, read the data for each column into cloud-amount_in, and then
!    close the file.
!---------------------------------------------------------------------
      if (trim(cloud_data_form) == 'input') then
        unit = open_namelist_file ('INPUT/cld_amt_file')
        do i=1,cloud_data_points
          read (unit,FMT = '(5e18.10)') (cloud_amount_in(i,k),k=1,kx)
        end do
        call close_file (unit)
      
!---------------------------------------------------------------------
!    if random overlap is assumed, compute the number of max_random and
!    random overlap clouds. the number of max_random is assumed to
!    be zero, while each layer with cloud present is assumed to be
!    a separate raandom overlap cloud. 
!---------------------------------------------------------------------
        if (trim(cloud_overlap_form) == 'random') then
          nrnd_cloud_in = COUNT(cloud_amount_in > 0.0E+00,DIM=2)
          nmax_cloud_in = 0
          max_cloud_amount_in = 0.0E+00
          rnd_cloud_amount_in = cloud_amount_in

!---------------------------------------------------------------------
!    if max_random overlap is assumed, compute the number of max_random
!    and random overlap clouds. 
!---------------------------------------------------------------------
        else if (trim(cloud_overlap_form) == 'max_random') then
          do i=1,cloud_data_points
            ktop = 1
            kbot = 1
            do k=1,kx-1
!----------------------------------------------------------------------
!    the cloud code below accounts for clouds in the (ktop, kbot) 
!    layers. cycle k to (kbot+1) to continue processing.
!----------------------------------------------------------------------
              if (k  > 1 .AND. k <= kbot) CYCLE

!----------------------------------------------------------------------
!    march downward in the column; find the next lower cloud top layer.
!----------------------------------------------------------------------
              if (cloud_amount_in(i,k) > 0.0E+00) then

!----------------------------------------------------------------------
!   determine the thickness of this cloud.
!----------------------------------------------------------------------
                ktop = k
                kbot = k
                do kk=ktop+1,kx   
!----------------------------------------------------------------------
!   find the base of the current cloud. at that point exit the loop.
!----------------------------------------------------------------------
                  if ( cloud_amount_in(i,kk) == 0.0E+00) EXIT
                  kbot = kk
                end do

!----------------------------------------------------------------------
!    if it is a single-layer cloud, assign it random cloud overlap prop-
!    erties.
!----------------------------------------------------------------------
                if (ktop == kbot) then
                  max_cloud_amount_in(i,k) = 0.0E+00
                  rnd_cloud_amount_in(i,k) = cloud_amount_in(i,k)
                  nrnd_cloud_in(i) = nrnd_cloud_in(i) + 1

!----------------------------------------------------------------------
!    if it is a multi-layer cloud, treat it as a max-random overlap
!    cloud. TEMPORARILY, set max cloud amount = max(cloud amounts from 
!    ktop to kbot) and set rnd cloud amount = zero.
!----------------------------------------------------------------------
                else
                  max_cld_calc = MAXVAL (cloud_amount_in(i,ktop:kbot))
                  max_cloud_amount_in(i,ktop:kbot) = max_cld_calc
                  rnd_cloud_amount_in(i,ktop:kbot) = 0.0E+00
                  nmax_cloud_in(i) = nmax_cloud_in(i) + 1
                endif

!----------------------------------------------------------------------
!    if cloud is not present at this level,  check the next level.
!----------------------------------------------------------------------
              else
                ktop = k
                kbot = k
              endif
            end do

!----------------------------------------------------------------------
!    deal with special case of cloud in lowest layer, no cloud in
!    next lowest layer (should never happen).
!----------------------------------------------------------------------
            if (cloud_amount_in(i,kx) > 0.0E+00 .AND.     &
                cloud_amount_in(i,kx-1) == 0.0E+00   ) then
              rnd_cloud_amount_in(i,kx   ) = cloud_amount_in(i,kx   )
              max_cloud_amount_in(i,kx   ) = 0.0E+00
              nrnd_cloud_in(i) = nrnd_cloud_in(i) + 1
            endif
          end do

!----------------------------------------------------------------------
!    error case: invalid specification of cloud_overlap_form.
!----------------------------------------------------------------------
        else
          call error_mesg( 'standalone_clouds_mod',  &
             ' cloud_overlap_form is not an acceptable value.', FATAL)
        endif
 
!----------------------------------------------------------------------
!    if cloud data is specified, define cloud amounts for high, middle
!    and low clouds. 
!----------------------------------------------------------------------
      else if (trim(cloud_data_form) == 'specified') then
        ch   = 0.159E+00
        cm   = 0.070E+00
        cl   = 0.269E+00

!----------------------------------------------------------------------
!    define the model k indices of the high and middle clouds (single 
!    layer) and the k indices of low cloud tops and bases (multi-layer 
!    clouds)  for levels corresponding to 3 different gcms: L40 SKYHI, 
!    the L18 NMC model, and the R30L14 supersource model. if column_type
!    is not recognized, print an error message and stop.
!----------------------------------------------------------------------
        if (trim(cldht_type_form)         == 'skyl40') then
          ich  = 29 
          icm  = 34
          ict  = 35 
          icb  = 37 
        else if (trim(cldht_type_form)         == 'nmcl18') then
          ich  =  5 
          icm  = 11 
          ict  = 12 
          icb  = 14 
        else if (trim(cldht_type_form)         == 'r30l14') then
          ich  = 6 
          icm  = 9 
          ict  = 10 
          icb  = 12 
        else if (trim(cldht_type_form)         == 'fmsl18') then
          ich  = 6   
          icm  = 9       
          ict  = 10  
          icb  = 12
        else
          call error_mesg ('standalone_clouds_mod', &
            'cloud properties have not been specified for'//&
             ' this column type', FATAL)
        endif

!---------------------------------------------------------------------
!    verify that the specified cloud levels are compatible with the
!    the model grid being used in this experiment.
!---------------------------------------------------------------------
        if (ich > kx  .or. icm > kx .or. ict > kx .or. icb > kx) then
          call error_mesg ('standalone_clouds_mod', &
            'specified cloud level index not within model grid', FATAL)
        endif
        if (ich < 1   .or. icm < 1  .or. ict < 1  .or. icb < 1) then
          call error_mesg ('standalone_clouds_mod', &
            'specified cloud level index not within model grid', FATAL)
        endif

!---------------------------------------------------------------------
!    save the cloud level indices for use elsewhere in defining the 
!    cloud radiative properties.
!---------------------------------------------------------------------
        Cldrad_control%ich = ich
        Cldrad_control%icm = icm
        Cldrad_control%ict = ict
        Cldrad_control%icb = icb

!---------------------------------------------------------------------
!    now specify the characteristics of three cloud columns; more can
!    be specified if desired, limited by the number of longitudes owned
!    by the processor.
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    the first cloud column contains 3 clouds, including a max-random 
!    overlap cloud.
!----------------------------------------------------------------------
        rnd_cloud_amount_in(1,ich) = ch
        rnd_cloud_amount_in(1,icm) = cm
        do k=ict,icb
          max_cloud_amount_in(1,k) = cl
        end do
        nmax_cloud_in(1) = 1
        nrnd_cloud_in(1) = 2

!---------------------------------------------------------------------
!    the second cloud column is cloudless.
!---------------------------------------------------------------------
        if (cloud_data_points > 1) then
          max_cloud_amount_in(2,:) = 0.0
          rnd_cloud_amount_in(2,:) = 0.0
          nmax_cloud_in(2) = 0
          nrnd_cloud_in(2) = 0
        endif

!---------------------------------------------------------------------
!    the third cloud column contains 5 random overlap clouds.
!---------------------------------------------------------------------
        if (cloud_data_points > 2) then
          rnd_cloud_amount_in(3,ich) = ch
          rnd_cloud_amount_in(3,icm) = cm
          do k=ict,icb
            rnd_cloud_amount_in(3,k) = cl
          enddo
          nmax_cloud_in(3) = 0
          nrnd_cloud_in(3) = 5
        endif

!---------------------------------------------------------------------
!    error case: invalid cloud_data_form
!---------------------------------------------------------------------
      else
        call error_mesg( 'standalone_clouds_mod',  &
             ' cloud_data_form is not an acceptable value.', FATAL)
      endif

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------


end subroutine standalone_clouds_init



!####################################################################

! <SUBROUTINE NAME="define_column_properties">
!  <OVERVIEW>
!    subroutine define_column_properties defines values for lw emiss-
!    ivity, visible and nir reflectivity and nir absorption to be used
!    with standalone clouds.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine define_column_properties defines values for lw emiss-
!    ivity, visible and nir reflectivity and nir absorption to be used
!    with standalone clouds.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_column_properties (pref, lonb, latb)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!       lonb      2d array of model longitudes at cell corners [ radians ]
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
!       latb      2d array of model latitudes at cell corners [radians]
! 
!  </IN>
! </SUBROUTINE>
!
subroutine define_column_properties (pref, lonb, latb)

!---------------------------------------------------------------------
!    subroutine define_column_properties defines values for lw emiss-
!    ivity, visible and nir reflectivity and nir absorption to be used
!    with standalone clouds.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes at cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      integer  ::   unit
      integer  ::   idf, jdf, kx
      integer  ::   ich, icm, ict, icb
      integer  ::   i, k, n

!---------------------------------------------------------------------
!    local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!      idf      x dimension of physics window
!      jdf      y dimension of physics window
!      kx       number of model layers
!      ich      model level index corresponding to level of high cloud
!      icm      model level index corresponding to level of middle cloud
!      ict      model level index corresponding to low-cloud cloud top
!      icb      model level index corresponding to low-cloud cloud base
!      i,k,n    do-loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the processor's domain dimensions. 
!---------------------------------------------------------------------
        jdf = size(latb,2) - 1
        idf = size(lonb,1) - 1
        kx  = size(pref,1) - 1

!-------------------------------------------------------------------
!    define the model levels corresponding to high cloud location, 
!    middle cloud location, low cloud top and low cloud base. 
!-------------------------------------------------------------------
        ich = Cldrad_control%ich
        icm = Cldrad_control%icm
        ict = Cldrad_control%ict
        icb = Cldrad_control%icb

!-------------------------------------------------------------------
!    allocate lw emissivity array and fill with either input or
!    specified values. default emissivity value is unity.
!-------------------------------------------------------------------
        allocate (emlw_band_in(idf, kx, Cldrad_control%nlwcldb))
        emlw_band_in(:,:,:) = 1.00E+00

!-------------------------------------------------------------------
!    if lw emissivity is to be obtained from an input file, allocate
!    an array into which it may be read, and then read the input data
!    for each desired cloud column. broadcast this data over all bands
!    of the emissivity variable.  close the input file unit. note that 
!    there must be consistency between the cloud fraction / altitude 
!    values and the emissivity values. this will be dealt with later.
!-------------------------------------------------------------------
        if (trim(lw_cld_prop_form) == 'input') then
          allocate (emlw_in (idf, kx))
          unit = open_namelist_file ('INPUT/lw_cld_prop_file')
          do i=1,Cldrad_control%cloud_data_points
            read (unit,FMT = '(5e18.10)')  (emlw_in(i,k),k=1,kx)
          end do
          do n=1,CLdrad_control%nlwcldb
            emlw_band_in(:,:,n) = emlw_in(:,:)
          end do
          call close_file (unit)

!---------------------------------------------------------------------
!    if using specified lw cloud emissivity, assign the values specified
!    in this module to the specified levels containing high, middle and
!    low clouds for each column that is to be integrated.
!---------------------------------------------------------------------
        else if (trim(lw_cld_prop_form) == 'specified') then
          if (trim(cloud_data_form) == 'specified') then
            do n=1,Cldrad_control%nlwcldb

!---------------------------------------------------------------------
!    the first specified column contains 3 clouds (max-random).
!----------------------------------------------------------------------
              emlw_band_in(1,ich,n) = highcloud_emiss
              emlw_band_in(1,icm,n) = midcloud_emiss
              do k=ict,icb
                emlw_band_in(1,k,n) = lowcloud_emiss
              end do
            end do

!---------------------------------------------------------------------
!    the second specified column contains no clouds - use default
!    value to which the array was initialized.
!----------------------------------------------------------------------
!           if (Cldrad_control%cloud_data_points > 1) then
!             emlw_band_in(2,:,:) = 1.00E+00
!           endif

!---------------------------------------------------------------------
!    the third specified column contains 5 random-overlap clouds.
!----------------------------------------------------------------------
            if (Cldrad_control%cloud_data_points > 2) then
              do n=1,CLdrad_control%nlwcldb
                emlw_band_in(3,ich,n) = highcloud_emiss
                emlw_band_in(3,icm,n) = midcloud_emiss
                do k=ict,icb
                  emlw_band_in(3,k,n) = lowcloud_emiss
                end do
              end do
            endif

!----------------------------------------------------------------------
!    if lw_cld_prop_form was specified then cloud_data_form must also be
!    specified. if this is not the case, write error message and stop
!    execution.
!---------------------------------------------------------------------
          else
            call error_mesg( 'standalone_clouds_mod',  &
                 ' if lw_cld_prop_form is specified  cloud_data_form'//&
                 ' must be specified.',                         FATAL)
          endif

!----------------------------------------------------------------------
!    if lw_cld_prop_form was neither specified nor input, write an
!    error message and stop execution.
!----------------------------------------------------------------------
        else
          call error_mesg ('standalone_clouds_mod', &
             'lw cld properties have not been specified correctly', &
                                                                FATAL)
        endif

!-------------------------------------------------------------------
!   if lhsw is active, allocate and initialize the lhsw cloud property 
!   arrays and either read an input file containing values for them or 
!   define them using the values specified in this module. the default
!   values are total reflection and zero absorption.
!-------------------------------------------------------------------
        if (Sw_control%do_lhsw) then
          allocate (cvis_rf_in (idf, kx))
          allocate (cir_rf_in  (idf, kx))
          allocate (cir_abs_in (idf, kx))
          cvis_rf_in(:,:) = 1.00E+00
          cir_rf_in (:,:) = 1.00E+00
          cir_abs_in(:,:) = 0.0

!----------------------------------------------------------------------
!    if values are supplied from an input file, open the file and read 
!    the values into the appropriate module variables.
!----------------------------------------------------------------------
          if (trim(lhsw_cld_prop_form) == 'input') then
            unit = open_namelist_file ('INPUT/lhsw_cld_prop_file')
            do i=1,Cldrad_control%cloud_data_points
              read (unit,FMT = '(5e18.10)') (cvis_rf_in(i,k),k=1,kx)
              read (unit,FMT = '(5e18.10)') (cir_rf_in(i,k),k=1,kx)
              read (unit,FMT = '(5e18.10)') (cir_abs_in(i,k),k=1,kx)
            end do
            call close_file (unit)
      
!----------------------------------------------------------------------
!    if one is to use the values provided by this module, define the
!    cldrad variables with the appropriate high, middle and low cloud
!    values. levels containing cloud are defined by ich, icm, ict and 
!    icb.
!---------------------------------------------------------------------
          else if (trim(lhsw_cld_prop_form) == 'specified') then

!---------------------------------------------------------------------
!    the first point contains 3 max-random overlap clouds. 
!---------------------------------------------------------------------
            cvis_rf_in(1,ich) = highcloud_refl_visband
            cvis_rf_in(1,icm) = midcloud_refl_visband
            cir_rf_in(1,ich)  = highcloud_refl_nearirband
            cir_rf_in(1,icm)  = midcloud_refl_nearirband
            cir_abs_in(1,ich) = highcloud_abs_nearirband
            cir_abs_in(1,icm) = midcloud_abs_nearirband
            do k=ict,icb
              cvis_rf_in(1,k) = lowcloud_refl_visband
              cir_rf_in(1,k)  = lowcloud_refl_nearirband
              cir_abs_in(1,k) = lowcloud_abs_nearirband
            end do

!---------------------------------------------------------------------
!    the second point contains no clouds -- use initialized values. 
!---------------------------------------------------------------------
!           if (Cldrad_control%cloud_data_points > 1) then
!             cvis_rf_in(2,:) = 1.00E+00
!             cir_rf_in(2,:)  = 1.00E+00
!             cir_abs_in(2,:) = 0.0
!           endif

!---------------------------------------------------------------------
!    the third point contains 5 random overlap clouds. 
!---------------------------------------------------------------------
            if (Cldrad_control%cloud_data_points > 2) then
              cvis_rf_in(3,ich) = highcloud_refl_visband
              cvis_rf_in(3,icm) = midcloud_refl_visband
              cir_rf_in(3,ich)  = highcloud_refl_nearirband
              cir_rf_in(3,icm)  = midcloud_refl_nearirband
              cir_abs_in(3,ich) = highcloud_abs_nearirband
              cir_abs_in(3,icm) = midcloud_abs_nearirband
              do k=ict,icb
                cvis_rf_in(3,k) = lowcloud_refl_visband
                cir_rf_in(3,k)  = lowcloud_refl_nearirband
                cir_abs_in(3,k) = lowcloud_abs_nearirband
              end do
            endif

!---------------------------------------------------------------------
!    if lhsw_cld_props was neither specified or input, write an error
!    message and stop execution.
!---------------------------------------------------------------------
          else
            call error_mesg ('standalone_clouds_mod', &
                    'lhsw_cld_prop_form has an improper value',  FATAL)
          endif
        endif !  (do_lhsw)

!---------------------------------------------------------------------




end subroutine define_column_properties




!######################################################################

! <SUBROUTINE NAME="standalone_clouds_end">
!  <OVERVIEW>
!    standalone_clouds_end is the destructor for standalone_clouds_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    standalone_clouds_end is the destructor for standalone_clouds_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call standalone_clouds_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine standalone_clouds_end
        
!----------------------------------------------------------------------
!    standalone_clouds_end is the destructor for standalone_clouds_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
        
!--------------------------------------------------------------------


end subroutine standalone_clouds_end



!#################################################################

! <SUBROUTINE NAME="standalone_clouds_amt">
!  <OVERVIEW>
!    standalone_clouds_amt defines the number, amount (cloud fraction),
!    and type (hi, mid, low) of clouds present on the model grid.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    standalone_clouds_amt defines the number, amount (cloud fraction),
!    and type (hi, mid, low) of clouds present on the model grid.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call standalone_clouds_amt (is, ie, js, je, lat, press_mks,  &
!                Cld_spec)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!      is,ie,js,je  starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <IN NAME="lat" TYPE="real">
!      lat          latitude of model points  [ radians ]
! 
!  </IN>
!  <IN NAME="press_mks" TYPE="real">
!      press_mks    pressure at model levels (1:nlev), surface
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
! 
!  </IN>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
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
!                  %ncldsw  number of clouds seen by he shortwave
!                           radiation in each grid column.
!                  %hi_cld  logical flag indicating the presence of
!                           high clouds in a grid box
!                 %mid_cld  logical flag indicating the presence of
!                           middle clouds in a grid box
!                 %low_cld  logical flag indicating the presence of
!                           low clouds in a grid box
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine standalone_clouds_amt (is, ie, js, je, lat, press_mks,  &
                                  Cld_spec)

!---------------------------------------------------------------------
!    standalone_clouds_amt defines the number, amount (cloud fraction), 
!    and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)     ::  is, ie, js, je
real,    dimension(:,:),      intent(in)     ::  lat  
real,    dimension(:,:,:),    intent(in)     ::  press_mks
type(cld_specification_type), intent(inout)  ::  Cld_spec

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      press_mks    pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!
!   intent(inout) variables:
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
!                  %ncldsw  number of clouds seen by he shortwave
!                           radiation in each grid column.
!                  %hi_cld  logical flag indicating the presence of 
!                           high clouds in a grid box
!                 %mid_cld  logical flag indicating the presence of 
!                           middle clouds in a grid box
!                 %low_cld  logical flag indicating the presence of 
!                           low clouds in a grid box
!                                                                  
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables

      real, dimension (size(Cld_spec%camtsw,2) ) :: cldhm, cldml

      real, dimension (size(Cld_spec%camtsw,1),                    &
                       size(Cld_spec%camtsw,2) ) ::  press_hm, press_ml

      integer     :: kx
      integer     :: i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!        cldhm             sigma value defining high-middle cloud boun-
!                          daries at window latitudes [ dimensionless ]
!        cldml             sigma value defining middle-low cloud boun-
!                          daries at window latitudes [ dimensionless ]
!        press_hm          pressure corresponding to the high-middle
!                          cloud boundary [ (kg /( m s^2) ]
!        press_ml          pressure corresponding to the middle-low
!                          cloud boundary [ (kg /( m s^2) ]
!        kx                number of model layers
!        i,j,k             do loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kx = size(Cld_spec%camtsw,3)

!---------------------------------------------------------------------
!    define the cloud fractions and number of clouds per column (random 
!    and max overlap) from the cloud specification data that was input. 
!    each cloud data point is replicated on each latitude row in the 
!    window.
!---------------------------------------------------------------------
      do j=1,size(Cld_spec%camtsw,2)
        do i=1,cloud_data_points
          Cld_spec%crndlw(i,j,:) = rnd_cloud_amount_in(i,:)
          Cld_spec%cmxolw(i,j,:) = max_cloud_amount_in(i,:)
          Cld_spec%nmxolw(i,j) = nmax_cloud_in(i)
          Cld_spec%nrndlw(i,j) = nrnd_cloud_in(i)
        end do
      end do

!!! NOTE:
!!!!  define Cld_spec%cld_thickness here (if needed)

!--------------------------------------------------------------------
!    sum up the random and max overlap cloud fractions and number of 
!    clouds per column to produce the values seen by the shortwave 
!    radiation. insure that cloud fraction at all grid points is <= 1.
!--------------------------------------------------------------------
      do j=1, size(Cld_spec%camtsw,2)
        do i=1,cloud_data_points
          Cld_spec%camtsw(i,j,:) = Cld_spec%crndlw(i,j,:) +    &
                                   Cld_spec%cmxolw(i,j,:)
          Cld_spec%ncldsw(i,j)   = Cld_spec%nmxolw(i,j) +      &
                                   Cld_spec%nrndlw(i,j)
        end do
      end do
      Cld_spec%camtsw = MIN (Cld_spec%camtsw,1.00)
      Cld_spec%cloud_area = Cld_spec%camtsw

!---------------------------------------------------------------------
!    determine the sigmas defining the hi-middle-lo cloud transitions
!    at each grid point. the sigmas vary with latitude between equator
!    and pole, and are linearly interpolated to the model latitude.
!---------------------------------------------------------------------
      do j=1,size(Cld_spec%camtsw,2)
        cldhm(j) = cldhp + ( (0.5*pie - abs(lat(1,j)))/(0.5*pie))* &
                   (cldhe - cldhp)
        cldml(j) = cldmp + ( (0.5*pie - abs(lat(1,j)))/(0.5*pie))* &
                   (cldme - cldmp)
      end do

!---------------------------------------------------------------------
!    define the pressures at the high-middle-low cloud boundaries.
!---------------------------------------------------------------------
      do j=1,size(Cld_spec%camtsw,2)
        do i=1,size(Cld_spec%camtsw,1)
          press_hm(i,j) = cldhm(j)*press_mks(i,j,kx+1)
          press_ml(i,j) = cldml(j)*press_mks(i,j,kx+1)
        end do
      end do

!--------------------------------------------------------------------
!    define flags indicating grid points at which either high, middle
!    or low clouds exist.
!---------------------------------------------------------------------
      do k=1,kx
        do j=1,size(Cld_spec%camtsw,2)
          do i=1,size(Cld_spec%camtsw,1)
            if (Cld_spec%camtsw(i,j,k) > 0.0) then
              if (press_mks(i,j,k) <= press_hm(i,j)) then
                Cld_spec%hi_cloud(i,j,k) = .true.
              else if (press_mks(i,j,k) > press_hm(i,j) .and.      &
                       press_mks(i,j,k) <= press_ml(i,j)) then
                Cld_spec%mid_cloud(i,j,k) = .true.
              else if (press_mks(i,j,k) > press_ml(i,j) ) then
                Cld_spec%low_cloud(i,j,k) = .true.
              endif
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------




end subroutine standalone_clouds_amt   



!#####################################################################

! <SUBROUTINE NAME="obtain_micro_lw_sa">
!  <OVERVIEW>
!    obtain_micro_lw_sa defines microphysically-based longwave cloud
!    radiative properties when the code is executed in standalone
!    columns mode.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_micro_lw_sa defines microphysically-based longwave cloud
!    radiative properties when the code is executed in standalone
!    columns mode.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_micro_lw_sa (is, ie, js, je, Lsc_microphys, &
!                Meso_microphys, Cell_microphys, &
!                Lscrad_props,  Mesorad_props, &
!                Cellrad_props)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!      is,ie,js,je  starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!      Lsc_microphys     microphysical specification for large-scale
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
! 
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!      Meso_microphys    microphysical specification for meso-scale
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
! 
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!      Cell_microphys    microphysical specification for cell-scale
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
! 
!  </INOUT>
!  <INOUT NAME="Lscrad_props" TYPE="microrad_properties_type">
!      Lscrad_props      cloud radiative properties on model grid,
!                        [ microrad_properties_type ]
! 
!  </INOUT>
!  <INOUT NAME="Mesorad_props" TYPE="microrad_properties_type">
!      Mesorad_props     meso-scale cloud radiative properties on
!                        model grid, [ microrad_properties_type ]
! 
!  </INOUT>
!  <INOUT NAME="Cellrad_props" TYPE="microrad_properties_type">
!      Cellrad_props     cell-scale cloud radiative properties on
!                        model grid, [ microrad_properties_type ]
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_micro_lw_sa (is, ie, js, je, Lsc_microphys, &
                               Meso_microphys, Cell_microphys, &
                               Lscrad_props,  Mesorad_props, &
                               Cellrad_props)

!---------------------------------------------------------------------
!    obtain_micro_lw_sa defines microphysically-based longwave cloud 
!    radiative properties when the code is executed in standalone 
!    columns mode.
!---------------------------------------------------------------------

integer,                        intent(in)    :: is, ie, js, je
type(microphysics_type),        intent(inout) :: Lsc_microphys, &
                                                 Meso_microphys, &
                                                 Cell_microphys
type(microrad_properties_type), intent(inout) :: Lscrad_props, &
                                                 Mesorad_props, &
                                                 Cellrad_props
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Lsc_microphys     microphysical specification for large-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Meso_microphys    microphysical specification for meso-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Cell_microphys    microphysical specification for cell-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Lscrad_props      cloud radiative properties on model grid,
!                        [ microrad_properties_type ]
!      Mesorad_props     meso-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Cellrad_props     cell-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!
!               the following component of the **_props variables is 
!               output from this routine:
!
!                    %abscoeff  absorption coefficient for  
!                               clouds in each of the longwave 
!                               frequency bands  [ km **(-1) ]
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------
!    call microphys_rad_driver2 to obtain the microphysically-based
!    sw cloud radiative quantities.
!---------------------------------------------------------------------
!     call microphys_rad_driver2 (is, ie, js, je, Lsc_microphys,   &
!     call microphys_lw_driver2 (is, ie, js, je, Lsc_microphys,   &
!                                 abscoeff=Lscrad_props%abscoeff)

!---------------------------------------------------------------------
!    if donner_deep is activated, process the cloud radiative properties
!    associated with the mesoscale and cellscale components of that
!    parameterization.
!---------------------------------------------------------------------
!     if (do_donner_deep_clouds) then
!     if (Cldrad_control%do_donner_deep_clouds) then
!        call obtain_micro_lw_donner_deep (is, ie, js, je,  &
!                                          Cell_microphys, & 
!                                          Meso_microphys,  &
!                                          Cellrad_props, Mesorad_props)
!     endif

!--------------------------------------------------------------------


end subroutine obtain_micro_lw_sa     




!#####################################################################

! <SUBROUTINE NAME="obtain_micro_sw_sa">
!  <OVERVIEW>
!    obtain_micro_sw_sa defines microphysically-based shortwave cloud
!    radiative properties for the standalone cloud scheme when run in
!    columns mode.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_micro_sw_sa defines microphysically-based shortwave cloud
!    radiative properties for the standalone cloud scheme when run in
!    columns mode.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_micro_sw_sa (is, ie, js, je, Lsc_microphys,   &
!                Meso_microphys, Cell_microphys,   &
!                Lscrad_props, Mesorad_props,   &
!                Cellrad_props)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <INOUT NAME="Lscrad_props" TYPE="microrad_properties_type">
! 
!  </INOUT>
!  <INOUT NAME="Mesorad_props" TYPE="microrad_properties_type">
! 
!  </INOUT>
!  <INOUT NAME="Cellrad_props" TYPE="microrad_properties_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_micro_sw_sa (is, ie, js, je, Lsc_microphys,   &
                               Meso_microphys, Cell_microphys,   &
                               Lscrad_props, Mesorad_props,   &
                               Cellrad_props)

!--------------------------------------------------------------------
!    obtain_micro_sw_sa defines microphysically-based shortwave cloud 
!    radiative properties for the standalone cloud scheme when run in 
!    columns mode.
!---------------------------------------------------------------------

integer,                         intent(in)    ::  is, ie, js, je
type(microphysics_type),         intent(inout) ::  Lsc_microphys, &
                                                   Meso_microphys,   &
                                                   Cell_microphys
type(microrad_properties_type),  intent(inout) ::  Lscrad_props,   &
                                                   Mesorad_props,  &
                                                   Cellrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Lsc_microphys     microphysical specification for large-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Meso_microphys    microphysical specification for meso-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Cell_microphys    microphysical specification for cell-scale 
!                        clouds, provides input to this subroutine
!                        [ microphysics_type ]
!      Lscrad_props      large-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Mesorad_props     meso-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!      Cellrad_props     cell-scale cloud radiative properties on 
!                        model grid, [ microrad_properties_type ]
!
!               the following components of the microrad_properties
!               variables are output from this routine:
!
!                   %cldext    sw extinction coefficient for  
!                              clouds in each of the shortwave 
!                              frequency bands  [ km **(-1) ]
!                   %cldsct    sw scattering coefficient for
!                              clouds in each of the shortwave
!                              frequency bands  [ km **(-1) ]
!                   %cldasymm  sw asymmetry factor for
!                              clouds in each of the shortwave 
!                              frequency bands  [ dimensionless ]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!    call microphys_rad_driver2 to obtain the microphysically-based
!    sw cloud radiative quantities.
!---------------------------------------------------------------------
!     call microphys_rad_driver2 (is, ie, js, je, Lsc_microphys, &
!     call microphys_sw_driver2 (is, ie, js, je, Lsc_microphys, &
!                                 cldext=Lscrad_props%cldext,   &
!                                 cldsct=Lscrad_props%cldsct,   &
!                                 cldasymm=Lscrad_props%cldasymm )

!---------------------------------------------------------------------
!    if donner_deep is activated, process the cloud radiative properties
!    associated with the mesoscale and cellscale components of that
!    parameterization.
!---------------------------------------------------------------------
!     if (do_donner_deep_clouds) then
!     if (CLdrad_control%do_donner_deep_clouds) then
!       call  obtain_micro_sw_donner_deep (is, ie, js, je,     &
!                                          Cell_microphys,&
!                                          Meso_microphys,     &
!                                          Cellrad_props, Mesorad_props)
!     endif

!--------------------------------------------------------------------



end subroutine obtain_micro_sw_sa     




!#####################################################################

! <SUBROUTINE NAME="obtain_bulk_lw_sa">
!  <OVERVIEW>
!    obtain_bulk_lw_sa defines bulk longwave cloud radiative properties
!    when using specified clouds in the standalone columns mode.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_lw_sa defines bulk longwave cloud radiative properties
!    when using specified clouds in the standalone columns mode.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_sa (is, ie, js, je, Cldrad_props)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!      is,ie,js,je  starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
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
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_bulk_lw_sa (is, ie, js, je, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_sa defines bulk longwave cloud radiative properties 
!    when using specified clouds in the standalone columns mode.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
type(cldrad_properties_type), intent(inout) :: Cldrad_props
 
!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
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

!----------------------------------------------------------------------
!   local variables:

      integer       :: i, j     ! do-loop indices

!--------------------------------------------------------------------
!    define the lw bulk cloud radiative properties from the values
!    previously defined in this module.
!--------------------------------------------------------------------
      do j=1,size(Cldrad_props%emrndlw,2)
        do i=1,cloud_data_points
          Cldrad_props%emmxolw(i,j,:,:,1) = emlw_band_in(i,:,:)
          Cldrad_props%emrndlw(i,j,:,:,1) = emlw_band_in(i,:,:)
        end do
      end do

!--------------------------------------------------------------------


end subroutine obtain_bulk_lw_sa     


!#####################################################################

! <SUBROUTINE NAME="obtain_bulk_sw_sa">
!  <OVERVIEW>
!    obtain_bulk_sw_sa defines bulk shortwave cloud radiative
!    properties for the specified cloud scheme when running in
!    standalone columns mode.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_sw_sa defines bulk shortwave cloud radiative
!    properties for the specified cloud scheme when running in
!    standalone columns mode.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_sa (is, ie, js, je, Cldrad_props)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine obtain_bulk_sw_sa (is, ie, js, je, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_sa defines bulk shortwave cloud radiative 
!    properties for the specified cloud scheme when running in 
!    standalone columns mode.
!---------------------------------------------------------------------

integer,                      intent(in)    ::   is, ie, js, je
type(cldrad_properties_type), intent(inout) ::   Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [ cldrad_properties_type ]
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

!----------------------------------------------------------------------
!   local variables

      integer       :: i, j    ! do-loop indices

!--------------------------------------------------------------------
!    define the sw bulk cloud radiative properties from the values
!    previously defined in this module.
!--------------------------------------------------------------------
      do j=1,size(Cldrad_props%cirabsw,2)
        do i=1,cloud_data_points
          Cldrad_props%cirabsw(i,j,:) = cir_abs_in(i,:)
          Cldrad_props%cirrfsw(i,j,:) = cir_rf_in(i,:)
          Cldrad_props%cvisrfsw(i,j,:) = cvis_rf_in(i,:)
        end do
      end do

!--------------------------------------------------------------------


end subroutine obtain_bulk_sw_sa     




!###################################################################

!subroutine find_nearest_index (latb, jindx2)

!real, dimension(:), intent(in) :: latb
!integer, dimension(:), intent(out)  :: jindx2
 

!      integer :: jd, j, jj
!     real   :: diff_low, diff_high
!      real, dimension(size(latb,1)-1) :: lat

 
!      jd = size(latb,1) - 1

!      do j = 1,jd
!        lat(j) = 0.5*(latb(j) + latb(j+1))
!      do jj=1, LATOBS
!        if (lat(j)*radians_to_degrees >= cloud_lats(jj)) then
!         diff_low = lat(j)*radians_to_degrees - cloud_lats(jj)
!          diff_high = cloud_lats(jj+1) - lat(j)*radians_to_degrees
!          if (diff_high <= diff_low) then
!            jindx2(j) = jj+1
!          else
!            jindx2(j) = jj
!          endif
!        endif
!      end do
!   end do






!end subroutine find_nearest_index

!#####################################################################




       end module standalone_clouds_mod




 

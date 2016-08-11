!FDOC_TAG_GFDL
               module bulkphys_rad_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!  
! </REVIEWER>
! <OVERVIEW>
!    bulkphys_rad_mod defines cloud radiative properties based on
!    bulk cloud physics values in contrast to microphysically-based
!    properties.
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!
 
!    shared modules:

use mpp_mod,                only: input_nml_file
use fms_mod,                only: open_namelist_file, mpp_pe, &
                                  fms_init, mpp_root_pe, stdlog,  &
                                  write_version_number, file_exist, & 
                                  check_nml_error, error_mesg,   &
                                  FATAL, close_file

!    shared radiation package modules:

use rad_utilities_mod,      only:  rad_utilities_init, &
                                  cldrad_properties_type, &
                                  cld_specification_type, &
                                  Cldrad_control

!    individual cloud modules:

use rh_based_clouds_mod,    only: rh_based_clouds_init, &
                                  obtain_bulk_sw_rh,   &
                                  obtain_bulk_lw_rh
use diag_clouds_W_mod,      only: diag_clouds_W_init, &
                                  obtain_bulk_sw_diag, &
                                  obtain_bulk_lw_diag
use strat_clouds_W_mod,     only: strat_clouds_W_init, &
                                  obtain_bulk_lw_strat, &
                                  obtain_bulk_sw_strat
use standalone_clouds_mod,  only: standalone_clouds_init,  &
                                  define_column_properties, &
                                  obtain_bulk_sw_sa, obtain_bulk_lw_sa

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    bulkphys_rad_mod defines cloud radiative properties based on
!    bulk cloud physics values in contrast to microphysically-based
!    properties.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: bulkphys_rad.F90,v 19.0 2012/01/06 20:13:07 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public                                            &
          bulkphys_rad_init, bulkphys_lw_driver,  &
          bulkphys_sw_driver, bulkphys_rad_end

!---------------------------------------------------------------------
!-------- namelist  ---------
 
integer     :: dummy = 1

namelist /bulkphys_rad_nml /  dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    visible band reflectivity, nir band reflectivity and nir absorp-
!    tivities are given for high, middle and low clouds. two separate
!    sets of parameters are present; set 1 has been used for all cloud
!    parameterizations except for mgrp_prscr_clds, which used the
!    values in data set 2.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   crfvis: visible band reflectivity. 
!--------------------------------------------------------------------
real  :: crfvis_hi_1  = 0.21
real  :: crfvis_hi_2  = 0.21
real  :: crfvis_mid_1 = 0.45
real  :: crfvis_mid_2 = 0.48
real  :: crfvis_low_1 = 0.59
real  :: crfvis_low_2 = 0.69

!--------------------------------------------------------------------
!   crfir: near-ir band reflectivity.
!--------------------------------------------------------------------
real  :: crfir_hi_1   = 0.21
real  :: crfir_hi_2   = 0.21
real  :: crfir_mid_1  = 0.45
real  :: crfir_mid_2  = 0.48
real  :: crfir_low_1  = 0.59
real  :: crfir_low_2  = 0.69

!--------------------------------------------------------------------
!   cabir: near-ir band absorptivity.
!--------------------------------------------------------------------
real  :: cabir_hi_1   = 0.005
real  :: cabir_hi_2   = 0.005
real  :: cabir_mid_1  = 0.02
real  :: cabir_mid_2  = 0.02
real  :: cabir_low_1  = 0.035
real  :: cabir_low_2  = 0.035

!--------------------------------------------------------------------
!   cldem:  infrared emissivity.
!--------------------------------------------------------------------
real  :: cldem_hi   = 1.00
real  :: cldem_mid  = 1.00
real  :: cldem_low  = 1.00

!------------------------------------------------------------------
!    these variables hold the values that are to be used for the
!    cloud radiative properties for the cloud parameterization 
!    activated.
!------------------------------------------------------------------
real  :: crfvis_hi, crfvis_mid, crfvis_low,    &
         crfir_hi,  crfir_mid,  crfir_low,  &
         cabir_hi,  cabir_mid,  cabir_low

real  :: min_cld_drop_rad, max_cld_drop_rad, &
         min_cld_ice_size, max_cld_ice_size

!-------------------------------------------------------------------
!    logical flag.
!-------------------------------------------------------------------
logical  :: module_is_initialized = .false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                           contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################

! <SUBROUTINE NAME="bulkphys_rad_init">
!  <OVERVIEW>
!    subroutine bulkphys_rad_init is the constructor for
!    bulkphys_rad_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine bulkphys_rad_init is the constructor for
!    bulkphys_rad_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_rad_init (pref, lonb, latb)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!       lonb      2d array of model longitudes on cell corners
!                 [ radians ]
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
!       latb      2d array of model latitudes at cell corners [radians]
! 
!  </IN>
! </SUBROUTINE>
!
subroutine bulkphys_rad_init (min_cld_drop_rad_in, max_cld_drop_rad_in,&
                              min_cld_ice_size_in, max_cld_ice_size_in,&
                              pref, lonb, latb)

!---------------------------------------------------------------------
!    subroutine bulkphys_rad_init is the constructor for 
!    bulkphys_rad_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,                 intent(in) :: min_cld_drop_rad_in, &
                                    max_cld_drop_rad_in,&
                                    min_cld_ice_size_in,  &
                                    max_cld_ice_size_in
real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes on cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      integer  ::   unit, ierr, io, logunit
      integer  ::   idum

!---------------------------------------------------------------------
!    local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!      idum     dummy integer argument needed to satisfy 
!               diag_clouds_W_init call
!
!--------------------------------------------------------------------

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
      if (Cldrad_control%do_diag_clouds)  call diag_clouds_W_init (idum)
      if (Cldrad_control%do_strat_clouds) call strat_clouds_W_init(latb, lonb)
      if (Cldrad_control%do_rh_clouds)    call rh_based_clouds_init 

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=bulkphys_rad_nml, iostat=io)
      ierr = check_nml_error(io,'bulkphys_rad_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=bulkphys_rad_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'bulkphys_rad_nml')
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
                      write (logunit, nml=bulkphys_rad_nml)

      min_cld_drop_rad = min_cld_drop_rad_in
      max_cld_drop_rad = max_cld_drop_rad_in
      min_cld_ice_size = min_cld_ice_size_in
      max_cld_ice_size = max_cld_ice_size_in

!--------------------------------------------------------------------
!    when executing the gcm or the standalone gcm, define the values
!    of absorptivity and reflectivity to be used for low, middle and
!    high clouds. the values used are different for different cloud 
!    parameterizations.
!---------------------------------------------------------------------
        if (Cldrad_control%do_mgroup_prescribed_iz) then
        if (Cldrad_control%do_mgroup_prescribed) then
          crfvis_hi  = crfvis_hi_2
          crfir_hi   = crfir_hi_2
          crfvis_mid = crfvis_mid_2
          crfir_mid  = crfir_mid_2
          crfvis_low = crfvis_low_2
          crfir_low  = crfir_low_2
          cabir_hi   = cabir_hi_2
          cabir_mid  = cabir_mid_2
          cabir_low  = cabir_low_2
        else
          crfvis_hi  = crfvis_hi_1
          crfir_hi   = crfir_hi_1
          crfvis_mid = crfvis_mid_1
          crfir_mid  = crfir_mid_1
          crfvis_low = crfvis_low_1
          crfir_low  = crfir_low_1
          cabir_hi   = cabir_hi_1
          cabir_mid  = cabir_mid_1
          cabir_low  = cabir_low_1
        endif
      else
        call error_mesg ('bulkphys_rad_mod', &
        ' do_mgroup_prescribed not yet defined', FATAL)
      endif

!---------------------------------------------------------------------
!    when running in standalone columns mode, call 
!    define_column_properties to define the cloud radiative properties.
!---------------------------------------------------------------------
    if (Cldrad_control%do_specified_strat_clouds_iz .and.  &
        Cldrad_control%do_specified_clouds_iz) then
       if (Cldrad_control%do_specified_strat_clouds  .or.  &
           Cldrad_control%do_specified_clouds ) then 
        call standalone_clouds_init   (pref, lonb, latb)
        call define_column_properties (pref, lonb, latb)
       endif 
    else
        call error_mesg ('bulkphys_rad_mod', &
        ' do_specified_strat not yet defined', FATAL)
   endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------



end subroutine bulkphys_rad_init



!#################################################################

! <SUBROUTINE NAME="bulkphys_sw_driver">
!  <OVERVIEW>
!    bulkphys_sw_driver obtains bulk shortwave cloud radiative
!    properties for the active cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    bulkphys_sw_driver obtains bulk shortwave cloud radiative
!    properties for the active cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_sw_driver (is, ie, js, je, cosz, Cld_spec,   &
!                Cldrad_props)
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
!  <IN NAME="cosz" TYPE="real">
!      cosz         cosine of the zenith angle [ dimensionless ]
! 
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!      Cld_spec     cloud specification arrays defining the
!                   location, amount and type (hi, middle, lo)
!                   of clouds that are present, provides input
!                   to this subroutine
!                   [ cld_specification_type ]
! 
!  </IN>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
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
!  </INOUT>
! </SUBROUTINE>
!
subroutine bulkphys_sw_driver (is, ie, js, je, cosz, Cld_spec,   &
                              Cldrad_props)

!---------------------------------------------------------------------
!    bulkphys_sw_driver obtains bulk shortwave cloud radiative 
!    properties for the active cloud scheme.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
real,    dimension(:,:),      intent(in)    :: cosz
type(cld_specification_type), intent(in)    :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle [ dimensionless ]
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

      integer     ::    i, j, k       ! do-loop indices

!---------------------------------------------------------------------
!    call obtain_bulk_sw_rh to obtain the cloud-radiative properties 
!    for rh_based_clouds.
!---------------------------------------------------------------------
      if (Cldrad_control%do_rh_clouds) then
        call obtain_bulk_sw_rh (is, ie, js, je, cosz, Cld_spec,   &
                                Cldrad_props)

!--------------------------------------------------------------------
!    values for zonal_clouds, prescribed_clouds and obs_clouds are
!    specified and constant in time. assign the proper values for cloud
!    absorptivity and reflectivity to each grid box with cloudiness, 
!    dependent on whether the cloud in that box is defined as being 
!    high, middle or low cloud. 
!----------------------------------------------------------------------
      else if(Cldrad_control%do_zonal_clouds .or. &
              Cldrad_control%do_mgroup_prescribed .or.  &
              Cldrad_control%do_obs_clouds) then
        do k=1, size(Cld_spec%hi_cloud,3)              
          do j=1,size(Cld_spec%hi_cloud,2)
            do i=1,size(Cld_spec%hi_cloud,1)
              if (Cld_spec%hi_cloud(i,j,k)) then
                Cldrad_props%cirabsw(i,j,k)  = cabir_hi 
                Cldrad_props%cirrfsw(i,j,k)  = crfir_hi
                Cldrad_props%cvisrfsw(i,j,k) = crfvis_hi
              else if (Cld_spec%mid_cloud(i,j,k)) then
                Cldrad_props%cirabsw(i,j,k)  = cabir_mid
                Cldrad_props%cirrfsw(i,j,k)  = crfir_mid
                Cldrad_props%cvisrfsw(i,j,k) = crfvis_mid
              else if (Cld_spec%low_cloud(i,j,k)) then
                Cldrad_props%cirabsw(i,j,k)  = cabir_low
                Cldrad_props%cirrfsw(i,j,k)  = crfir_low
                Cldrad_props%cvisrfsw(i,j,k) = crfvis_low
              endif
            end do
          end do
        end do

!---------------------------------------------------------------------
!    call obtain_bulk_sw_diag to define the cloud radiative properties
!    for the gordon diagnostic cloud scheme.
!---------------------------------------------------------------------
      else if (Cldrad_control%do_diag_clouds) then
        call obtain_bulk_sw_diag (is, ie, js, je, cosz, Cld_spec,  &
                                  Cldrad_props)

!---------------------------------------------------------------------
!    call obtain_bulk_sw_strat to define the cloud radiative properties
!    for the klein prognostic cloud scheme.
!---------------------------------------------------------------------
      else if (Cldrad_control%do_strat_clouds) then
        call obtain_bulk_sw_strat (is, ie, js, je, cosz, Cld_spec,   &
                                   Cldrad_props)

!-------------------------------------------------------------------
!    call obtain_bulk_sw_sa to define the cloud radiative properties
!    when specified clouds are used in standalone columns mode.
!-------------------------------------------------------------------
      else if (Cldrad_control%do_specified_clouds) then
        call obtain_bulk_sw_sa (is, ie, js, je, Cldrad_props)
      endif

!----------------------------------------------------------------------
 

end subroutine bulkphys_sw_driver



!####################################################################

! <SUBROUTINE NAME="bulkphys_lw_driver">
!  <OVERVIEW>
!    bulkphys_lw_driver defines bulk longwave cloud radiative
!    properties for the active cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    bulkphys_lw_driver defines bulk longwave cloud radiative
!    properties for the active cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_lw_driver (is, ie, js, je, Cld_spec, Cldrad_props)
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
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!      Cld_spec          cloud specification arrays defining the
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input
!                        to this subroutine
!                        [ cld_specification_type ]
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
subroutine bulkphys_lw_driver (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    bulkphys_lw_driver defines bulk longwave cloud radiative 
!    properties for the active cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(in)    :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
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

      integer    :: i, j, k      ! do-loop indices

!------------------------------------------------------------------
!    call obtain_bulk_lw_rh to define long-wave cloud emissivity for
!    rh_based_clouds_mod.
!-------------------------------------------------------------------
      if (Cldrad_control%do_rh_clouds) then
        call obtain_bulk_lw_rh (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    assign the proper values for cloud emissivity to each grid box
!    with cloudiness, dependent on whether the cloud in that box is 
!    defined as being high, middle or low cloud. high and middle clouds
!    are assumed to be random overlap, low clouds are assume to be
!    maximum overlap.
!----------------------------------------------------------------------
      else if (Cldrad_control%do_zonal_clouds .or.   &
               Cldrad_control%do_mgroup_prescribed .or.  &
               Cldrad_control%do_obs_clouds)  then
        do k=1, size(Cld_spec%hi_cloud,3)              
          do j=1,size(Cld_spec%hi_cloud,2)
            do i=1,size(Cld_spec%hi_cloud,1)
              if (Cld_spec%hi_cloud(i,j,k)) then
                Cldrad_props%emrndlw(i,j,k,:,1)  = cldem_hi
              else if (Cld_spec%mid_cloud(i,j,k)) then
                Cldrad_props%emrndlw(i,j,k,:,1)  = cldem_mid
              else if (Cld_spec%low_cloud(i,j,k)) then
                Cldrad_props%emmxolw(i,j,k,:,1)  = cldem_low
              endif
            end do
          end do
        end do

!------------------------------------------------------------------
!    call obtain_bulk_lw_diag to define long-wave cloud emissivity for
!    diag_based_clouds_mod.
!-------------------------------------------------------------------
      else if (Cldrad_control%do_diag_clouds) then
        call obtain_bulk_lw_diag (is, ie, js, je, Cld_spec,  &
                                  Cldrad_props)
 
!------------------------------------------------------------------
!    call obtain_bulk_lw_strat to define long-wave cloud emissivity for
!    strat_clouds_mod.
!-------------------------------------------------------------------
      else if (Cldrad_control%do_strat_clouds) then
        call obtain_bulk_lw_strat (is, ie, js, je, Cld_spec,   &
                                   Cldrad_props)

!--------------------------------------------------------------------
!    call obtain_bulk_lw_sa to define long-wave cloud emissivity for
!    specified clouds when running in standalone columns mode.
!-------------------------------------------------------------------
      else if (Cldrad_control%do_specified_clouds) then
        call obtain_bulk_lw_sa (is, ie, js, je, Cldrad_props)
      endif

!---------------------------------------------------------------------



end subroutine bulkphys_lw_driver



!###################################################################
 
! <SUBROUTINE NAME="bulkphys_rad_end">
!  <OVERVIEW>
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call bulkphys_rad_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine bulkphys_rad_end

!-------------------------------------------------------------------
!    bulkphys_rad_end is the destructor for bulkphys_rad_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
     module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine bulkphys_rad_end



                     end module bulkphys_rad_mod

 

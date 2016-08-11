!FDOC_TAG_GFDL

                 module specified_clouds_W_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
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

use  time_manager_mod,  only:  time_type
use   cloud_zonal_mod,  only:  getcld
use     cloud_obs_mod,  only:  cloud_obs, cloud_obs_init
use           mpp_mod,  only:  input_nml_file
use           fms_mod,  only:  open_namelist_file, file_exist, &
                               check_nml_error, &
                               close_file, &
                               mpp_pe, mpp_root_pe, &
                               write_version_number, stdlog
use rad_utilities_mod,  only:  cldrad_properties_type, &
                               Cldrad_control, &
                               cld_specification_type
use     constants_mod,  only:  pstd_mks

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!             specified clouds radiative properties module;
!             used with cloud_obs_mod and cloud_zonal_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: specified_clouds_W.F90,v 19.0 2012/01/06 20:24:09 fms Exp $'
  character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          specified_clouds_W_init, specified_clouds_amt, &
          specified_clouds_W_end

!---------------------------------------------------------------------
!-------- namelist  ---------


integer    :: dummy=0           



namelist /specified_clouds_W_nml /     &
                                         dummy   


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

integer, parameter  :: NOFCLDS_SP=3  ! total number of clouds per column
integer, parameter  :: NOFMXOLW=1    ! number of max overlap clouds
integer, parameter  :: NOFRNDLW=2    ! number of random overlap clouds
 

!--------------------------------------------------------------------
!   crfvis   :  visible band reflectivity
!--------------------------------------------------------------------
real  :: crfvis_hi  = 0.21
real  :: crfvis_mid = 0.45
real  :: crfvis_low = 0.59

!--------------------------------------------------------------------
!   crfir    :  near-ir band reflectivity
!--------------------------------------------------------------------
real  :: crfir_hi   = 0.21
real  :: crfir_mid  = 0.45
real  :: crfir_low  = 0.59

!--------------------------------------------------------------------
!   cldem    :  infrared emissivity
!--------------------------------------------------------------------
real  :: cldem_hi   = 1.00
real  :: cldem_mid  = 1.00
real  :: cldem_low  = 1.00

!--------------------------------------------------------------------
!   cabir    :  near-ir band absorptivity
!--------------------------------------------------------------------
real  :: cabir_hi   = 0.005
real  :: cabir_mid  = 0.02
real  :: cabir_low  = 0.035


logical  :: module_is_initialized = .false.

!------------------------------------------------------------------
!------------------------------------------------------------------



contains 





! <SUBROUTINE NAME="specified_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call specified_clouds_W_init (lonb, latb)
!
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine specified_clouds_W_init (lonb, latb)


real, dimension(:,:), intent(in) :: lonb, latb


      integer          :: unit, ierr, io, logunit

!---------------------------------------------------------------------
!-----  read namelist  ------
  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=specified_clouds_W_nml, iostat=io)
      ierr = check_nml_error(io,'specified_clouds_W_nml')
#else   
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=specified_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'specified_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
        call write_version_number(version, tagname)
        logunit = stdlog()
        write (logunit,nml=specified_clouds_W_nml)
      endif

!---------------------------------------------------------------------
!    if observed clouds is active, initialize that module.
!---------------------------------------------------------------------
      if (Cldrad_control%do_obs_clouds) then
        call cloud_obs_init (lonb, latb)
      endif
 
      module_is_initialized = .true.

end subroutine specified_clouds_W_init

! <SUBROUTINE NAME="specified_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call specified_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine specified_clouds_W_end
        
!----------------------------------------------------------------------
!    specified_clouds_W_end is the destructor for specified_clouds_W_mod.
!----------------------------------------------------------------------
        
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------
 

end subroutine specified_clouds_W_end


!######################################################################

! <SUBROUTINE NAME="specified_clouds_amt">
!  <OVERVIEW>
!    specified_clouds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    specified_clouds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call specified_clouds_amt (is, ie, js, je, Rad_time, lat, pflux, &
!                Cld_spec)
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
!  <IN NAME="Rad_time" TYPE="time_type">
!      Rad_time     time at which the climatologically-determined,
!                   time-varying specified cloud fields should apply
!                   [ time_type, days and seconds]
! 
!  </IN>
!  <IN NAME="lat" TYPE="real">
!      lat          latitude of model points  [ radians ]
! 
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!      pflux        average of pressure at adjacent model levels
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
subroutine specified_clouds_amt (is, ie, js, je, Rad_time, lat, pflux, &
                                 Cld_spec)

!----------------------------------------------------------------------
!    specified_clouds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------
 
!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
type(time_type),              intent(in)    :: Rad_time
real, dimension(:,:),         intent(in)    :: lat
real, dimension(:,:,:),       intent(in)    :: pflux
type(cld_specification_type), intent(inout) :: Cld_spec
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Rad_time     time at which the climatologically-determined, 
!                   time-varying specified cloud fields should apply
!                   [ time_type, days and seconds]
!      lat          latitude of model points  [ radians ]
!      pflux        average of pressure at adjacent model levels
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

!---------------------------------------------------------------------
!   local variables:

      real, dimension    (size (Cld_spec%camtsw,1),                 &
                          size (Cld_spec%camtsw,2),                 &
                          NOFCLDS_SP       ) ::      camtsw3
      integer, dimension (size (Cld_spec%camtsw,1),                 &
                          size (Cld_spec%camtsw,2),                 &
                          NOFCLDS_SP       ) ::      ktopsw3, kbtmsw3
      real, dimension    (size (Cld_spec%camtsw,1),                 &
                          size (Cld_spec%camtsw,2),                 &
                          size (Cld_spec%camtsw,3)+1)  ::   phaf      

      integer  ::     k, j, i
      integer  ::     kerad

!---------------------------------------------------------------------
!    local variables:
!
!
!         camtsw3   cloud fraction in cloud space for the specified 
!                   cloud types;  currently, k = 1 is for hi cloud,
!                   k = 2 is for mid cloud, k = 3 is for low cloud
!         ktopsw3   model k index of cloud top for the specified cloud 
!                   types (hi, mid, low)
!         kbtmsw3   model k index of the cloud base for the specified 
!                   cloud types (hi, mid, low)
!            phaf   pressure at model interface levels adjusted so that
!                   a sigma value based on the mean sea level pressure 
!                   is the same as the sigma value based on actual 
!                   surface pressure
!         i, j, k   do loop indices
!           kerad   number of model layers
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kerad = size (Cld_spec%camtsw,3)

!-----------------------------------------------------------------------
!    define the number of random and maximally overlapped clouds and the
!    total number of clouds in a column. these numbers are prescribed.
!-----------------------------------------------------------------------
      Cld_spec%nmxolw(:,:) = NOFMXOLW
      Cld_spec%nrndlw(:,:) = NOFRNDLW
      Cld_spec%ncldsw(:,:) = NOFCLDS_SP 

!---------------------------------------------------------------------
!    define the interface pressure that produces the same value of
!    sigma based on the mean sea level pressure as is obtained using
!    the actual interface pressures and surface pressure.
!---------------------------------------------------------------------
      do k=1, kerad+1        
        phaf(:,:,k) = pflux(:,:,k)*pstd_mks/pflux(:,:,kerad+1)
      end do

!------------------------------------------------------------------
!    call getcld to obtain the cloud fractions, tops and bases for
!    the specified clouds. these are returned in cloud-space.
!------------------------------------------------------------------
      call getcld (Rad_time, lat, phaf, ktopsw3, kbtmsw3, camtsw3 )
      if (Cldrad_control%do_obs_clouds) then
        call cloud_obs (is, js, Rad_time, camtsw3)
      endif


!---------------------------------------------------------------------
!    map the cloud-space arrays obtained above to model space arrays. 
!    define the logical arrays which denote the grid boxes containing
!    the various cloud types (hi, middle, low).
!-------------------------------------------------------------------
      do j=1, size(Cld_spec%hi_cloud,2)
        do i=1,size(Cld_spec%hi_cloud,1)

!--------------------------------------------------------------------
!    high clouds
!--------------------------------------------------------------------
          do k=ktopsw3(i,j,1), kbtmsw3(i,j,1)-1
            Cld_spec%camtsw(i,j,k) = camtsw3(i,j,1)
            Cld_spec%hi_cloud(i,j,k) = .true.
            Cld_spec%crndlw(i,j,k) = camtsw3(i,j,1)
          end do

!--------------------------------------------------------------------
!    middle clouds
!--------------------------------------------------------------------
          do k=ktopsw3(i,j,2), kbtmsw3(i,j,2)-1
            Cld_spec%camtsw(i,j,k) = camtsw3(i,j,2)
            Cld_spec%mid_cloud(i,j,k) = .true.
            Cld_spec%crndlw(i,j,k) = camtsw3(i,j,2)
          end do

!--------------------------------------------------------------------
!    low clouds
!--------------------------------------------------------------------
          do k=ktopsw3(i,j,3), kbtmsw3(i,j,3)-1
            Cld_spec%camtsw(i,j,k) = camtsw3(i,j,3)
            Cld_spec%low_cloud(i,j,k) = .true.
            Cld_spec%cmxolw(i,j,k) = camtsw3(i,j,3)
          end do
        end do
      end do

!---------------------------------------------------------------------


end subroutine specified_clouds_amt 


!####################################################################

! <SUBROUTINE NAME="obtain_bulk_lw_specified">
!  <OVERVIEW>
!    obtain_bulk_lw_specified defines bulk longwave cloud radiative
!    properties for the specified cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_lw_specified defines bulk longwave cloud radiative
!    properties for the specified cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_specified (is, ie, js, je, Cld_spec,   &
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
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!      Cld_spec          cloud specification arrays defining the
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input
!                        to this subroutine
!                        [ cld_specification_type ]
! 
!  </INOUT>
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
subroutine obtain_bulk_lw_specified (is, ie, js, je, Cld_spec,   &
                                     Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_specified defines bulk longwave cloud radiative 
!    properties for the specified cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(inout) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
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

!---------------------------------------------------------------------
!   local variables:

      integer     ::  i,j,k    ! do loop indices

!---------------------------------------------------------------------
!    assign the proper values for cloud emissivity to each grid box
!    with cloudiness, dependent on whether the cloud in that box is 
!    defined as being high, middle or low cloud. high and middle clouds
!    are assumed to be random overlap, low clouds are assume to be
!    maximum overlap.
!----------------------------------------------------------------------
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
          enddo
        end do
      end do

!---------------------------------------------------------------------



end subroutine obtain_bulk_lw_specified



!####################################################################

! <SUBROUTINE NAME="obtain_bulk_sw_specified">
!  <OVERVIEW>
!    obtain_bulk_sw_specified defines bulk shortwave cloud radiative
!    properties for the specified cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_sw_specified defines bulk shortwave cloud radiative
!    properties for the specified cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_specified (is, ie, js, je, Cld_spec, &
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
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!      Cld_spec          cloud specification arrays defining the
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input
!                        to this subroutine
!                        [ cld_specification_type ]
! 
!  </INOUT>
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
subroutine obtain_bulk_sw_specified (is, ie, js, je, Cld_spec, &
                                 Cldrad_props)                     

!---------------------------------------------------------------------
!    obtain_bulk_sw_specified defines bulk shortwave cloud radiative 
!    properties for the specified cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(cld_specification_type), intent(inout) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props
!-------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!
!      Cld_spec          cloud specification arrays defining the 
!                        location, amount and type (hi, middle, lo)
!                        of clouds that are present, provides input 
!                        to this subroutine
!                        [ cld_specification_type ]
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

!---------------------------------------------------------------------
!   local variables:

      integer     ::  i,j,k    ! do loop indices

!---------------------------------------------------------------------
!    assign the proper values for cloud absorptivity and reflectivity
!    to each grid box with cloudiness, dependent on whether the cloud 
!    in that box is defined as being high, middle or low cloud. 
!----------------------------------------------------------------------
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



end subroutine obtain_bulk_sw_specified 




!######################################################################



                 end module specified_clouds_W_mod




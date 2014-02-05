!FDOC_TAG_GFDL

                 module diag_clouds_W_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!    fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <OVERVIEW>
!           diag cloud radiative properties module
!            currently a wrapper until SKYHI goes away and this
!            module can be consolidated with diag_cloud_mod
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use time_manager_mod,       only: time_type
use diag_cloud_mod,         only: diag_cloud_avg2, diag_cloud_driver2
use diag_cloud_rad_mod,     only: cloud_opt_prop_tg_lw,  &
                                  cloud_opt_prop_tg_sw
use mpp_mod,                only: input_nml_file
use fms_mod,                only: open_namelist_file, file_exist, &
                                  check_nml_error, &
                                  write_version_number, &
                                  mpp_pe, mpp_root_pe, &
                                  close_file, stdlog
use rad_utilities_mod,      only: microphysics_type, &
                                  cld_specification_type, &
                                  cldrad_properties_type, &
                                  Sw_control

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!           diag cloud radiative properties module
!            currently a wrapper until SKYHI goes away and this
!            module can be consolidated with diag_cloud_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: diag_clouds_W.F90,v 19.0 2012/01/06 20:14:45 fms Exp $'
   character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          diag_clouds_W_init,    &
          diag_clouds_W_end,    &
          diag_clouds_amt,  &
          obtain_bulk_lw_diag, &
          obtain_bulk_sw_diag

!---------------------------------------------------------------------
!-------- namelist  ---------

real       :: taucrit = 1.0     !  critical optical depth at which 
                                !  solar beam is treated as diffuse
                                !  rather than direct
integer    :: num_slingo_bands = 4



namelist /diag_clouds_W_nml /     &
                               num_slingo_bands, &
                               taucrit


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





! <SUBROUTINE NAME="diag_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_clouds_W_init  (num_slingo_bands_out)
!
!  </TEMPLATE>
!  <OUT NAME="num_slingo_bands_out" TYPE="integer">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine diag_clouds_W_init  (num_slingo_bands_out)


integer, intent(out) :: num_slingo_bands_out

      integer            :: unit, ierr, io, logunit


      if (module_is_initialized) return
!---------------------------------------------------------------------
!-----  read namelist  ------
  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=diag_clouds_W_nml, iostat=io)
      ierr = check_nml_error(io,'diag_clouds_W_nml')
#else   
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=diag_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'diag_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif
#endif

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit,nml=diag_clouds_W_nml)
      endif

      num_slingo_bands_out = num_slingo_bands
! (ultimately get from diag_cloud_mod when call diag_cloud_init)
      module_is_initialized = .true.

end subroutine diag_clouds_W_init


!#####################################################################

! <SUBROUTINE NAME="diag_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine diag_clouds_W_end
 
!----------------------------------------------------------------------
!    diag_clouds_W_end is the destructor for diag_clouds_W_mod.
!----------------------------------------------------------------------
 
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine diag_clouds_W_end



!#################################################################

! <SUBROUTINE NAME="diag_clouds_amt">
!  <OVERVIEW>
!    diag_clouds_amt defines the location, amount (cloud fraction),
!    number, optical depth, thickness and liquid percentage of clouds
!    present on the model grid.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_clouds_amt defines the location, amount (cloud fraction),
!    number, optical depth, thickness and liquid percentage of clouds
!    present on the model grid.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_clouds_amt (is, ie, js, je, lat, pflux, press,   &
!                Rad_time, Cld_spec, Lsc_microphys) 
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
!  <IN NAME="lat" TYPE="real">
!      lat          latitude of model points  [ radians ]
! 
!  </IN>
!  <IN NAME="pflux" TYPE="real">
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ]
! 
!  </IN>
!  <IN NAME="press" TYPE="real">
!      press        pressure at model levels (1:nlev), surface
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
! 
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!      Rad_time     time at which the climatologically-determined,
!                   time-varying zonal cloud fields should apply
!                   [ time_type, days and seconds]
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
!                  %liq_frac
!                           percentage of cloud condensate in a grid
!                           box which is liquid  [ dimensionless ]
!                  %tau     cloud optical depth  [ dimensionless ]
!                  %cloud_thickness
!                           number of model layers over which the cloud
!                           in this grid box extends
!                  %ice_cloud
!                           logical variable, which if true, indicates
!                           that the grid box will contain ice cloud;
!                           if false, the box will contain liquid cloud
! 
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine diag_clouds_amt (is, ie, js, je, lat, pflux, press,   &
                            Rad_time, Cld_spec, Lsc_microphys) 

!----------------------------------------------------------------------
!    diag_clouds_amt defines the location, amount (cloud fraction), 
!    number, optical depth, thickness and liquid percentage of clouds 
!    present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)     ::  is, ie, js, je
real,    dimension(:,:),      intent(in)     ::  lat
real,    dimension(:,:,:),    intent(in)     ::  pflux, press
type(time_type),              intent(in)     ::  Rad_time     
type(cld_specification_type), intent(inout)  ::  Cld_spec
type(microphysics_type),      intent(inout)  ::  Lsc_microphys

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ] 
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      Rad_time     time at which the climatologically-determined, 
!                   time-varying zonal cloud fields should apply
!                   [ time_type, days and seconds]
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
!                  %liq_frac 
!                           percentage of cloud condensate in a grid 
!                           box which is liquid  [ dimensionless ]
!                  %tau     cloud optical depth  [ dimensionless ]
!                  %cloud_thickness
!                           number of model layers over which the cloud
!                           in this grid box extends
!                  %ice_cloud  
!                           logical variable, which if true, indicates 
!                           that the grid box will contain ice cloud; 
!                           if false, the box will contain liquid cloud
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables

      logical, dimension (size(Cld_spec%camtsw,1),                 &
                          size(Cld_spec%camtsw,2),                 &
                          size(Cld_spec%camtsw,3)) :: ice_cloud_cs

      real, dimension    (size(Cld_spec%camtsw,1),                 &
                          size(Cld_spec%camtsw,2),                 &
                          size(Cld_spec%camtsw,3)) ::  cldamt,     &
                                                       liq_frac

      real, dimension    (size(Cld_spec%camtsw,1),                 &
                          size(Cld_spec%camtsw,2),                 &
                          size(Cld_spec%camtsw,3),                 &
                          size(Cld_spec%tau,4))    ::  tau

      integer, dimension (size(Cld_spec%camtsw,1),                 &
                          size(Cld_spec%camtsw,2),                 &
                          size(Cld_spec%camtsw,3)) ::  ktop, kbtm

      integer     ::   i,j, k, kc           
      integer     ::   kx                     

!----------------------------------------------------------------------
!  local variables:
!
!      ice_cloud_cs  logical flag indicating whether a given cloud (in
!                    cloud-space) is made up of ice (.true.) or liquid
!                    (.false.)
!      cldamt        fractional cloudiness in a grid box 
!                    (in cloud-space)  [ dimensionless ]
!      liq_frac      the fraction of cloud in a grid box which is liquid
!                    (in cloud-space) [ dimensionless ]
!      tau           cloud optical depth, in cloud-space   
!                    [ dimensionless ]
!      ktop          model index of the cloud top
!      kbtm          model index of the cloud base
!      i, j, k, kc   do loop indices
!      kx            number of model layers
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    define number of model layers.
!--------------------------------------------------------------------
      kx = size(Cld_spec%camtsw,3)

!---------------------------------------------------------------------
!    call diag_cloud_driver2 to obtain the cloud specification variables
!    in cloud space.
!---------------------------------------------------------------------
      call diag_cloud_driver2 (is, js, press, pflux, lat, Rad_time,  &
                               Cld_spec%ncldsw, ktop, kbtm, cldamt, &
                               liq_frac, tau, ice_cloud_cs)

!---------------------------------------------------------------------
!    map the various cloud specification arrays from cloud-space to 
!    physical space so that they may be exported for use elsewhere.
!-------------------------------------------------------------------
      do j=1,size(Cld_spec%camtsw,2)
        do i=1,size(Cld_spec%camtsw,1)
          do kc=1,Cld_spec%ncldsw(i,j)  
            do k=ktop(i,j,kc), kbtm(i,j,kc)
              Cld_spec%camtsw(i,j,k) = cldamt(i,j,kc) 
              Cld_spec%liq_frac(i,j,k) = liq_frac(i,j,kc)
              Cld_spec%tau(i,j,k,:) = tau(i,j,kc,:)
              Cld_spec%cld_thickness(i,j,k) =    &
                                      kbtm(i,j,kc) - ktop(i,j,kc) + 1
              Cld_spec%ice_cloud(i,j,k) = ice_cloud_cs(i,j,kc)
      
!---------------------------------------------------------------------
!    determine if max overlap or random overlap assumption is made. if
!    cloud is more than one layer deep and using lhsw, then max overlap
!    assumed; otherwise random overlap is assumed.
!----------------------------------------------------------------------
              if (ktop(i,j,kc) == kbtm(i,j,kc)) then
                Cld_spec%crndlw(i,j,k) = cldamt(i,j,kc)
                Cld_spec%cmxolw(i,j,k) = 0.0             
              else
                if (Sw_control%do_esfsw) then
                  Cld_spec%crndlw(i,j,k) = cldamt(i,j,kc)  
                  Cld_spec%cmxolw(i,j,k) = 0.0             
                else
                  Cld_spec%cmxolw(i,j,k) = cldamt(i,j,kc)
                  Cld_spec%crndlw(i,j,k) = 0.0
                endif       
              endif
            end do
            if (ktop(i,j,kc) == kbtm(i,j,kc)) then
              Cld_spec%nrndlw(i,j) = Cld_spec%nrndlw(i,j) + 1
            else
              if (Sw_control%do_esfsw) then
                Cld_spec%nrndlw(i,j) = Cld_spec%nrndlw(i,j) + 1
              else
                Cld_spec%nmxolw(i,j) = Cld_spec%nmxolw(i,j) + 1
              endif     
            endif
          end do
        end do
      end do
 

!--------------------------------------------------------------------

end subroutine diag_clouds_amt 





!#####################################################################

! <SUBROUTINE NAME="obtain_bulk_lw_diag">
!  <OVERVIEW>
!    obtain_bulk_lw_diag defines bulk longwave cloud radiative
!    properties for the gordon diag cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_lw_diag defines bulk longwave cloud radiative
!    properties for the gordon diag cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_diag (is, ie, js, je, Cld_spec, Cldrad_props)
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
subroutine obtain_bulk_lw_diag (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_diag defines bulk longwave cloud radiative 
!    properties for the gordon diag cloud scheme.
!---------------------------------------------------------------------
 
integer,                     intent(in)     :: is, ie, js, je
type(cld_specification_type), intent(in   ) :: Cld_spec
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

!-------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cldrad_props%emrndlw,1),               &
                       size(Cldrad_props%emrndlw,2),                &
                       size(Cldrad_props%emrndlw,3)) ::  emcld

      integer    :: max_cld
      integer    :: i,j,k

!---------------------------------------------------------------------
!   local variables:
!
!            emcld    longwave cloud emissivity, assuming a single 
!                     band, returned from diag_cloud_rad_mod
!            max_cld  maximum number of clouds in any column in the 
!                     window
!            i,j,k    do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    determine if any clouds are present in the window.
!---------------------------------------------------------------------
      max_cld = maxval(Cld_spec%ncldsw)

!---------------------------------------------------------------------
!    if cloud is nowhere present, return. 
!---------------------------------------------------------------------
      if (max_cld > 0) then

!--------------------------------------------------------------------
!    initialize property array in cloud-space to zero.
!--------------------------------------------------------------------
        emcld = 0.

!---------------------------------------------------------------------
!    call cloud_opt_prop_tg_lw to obtain the cloud longwave enmissivity.
!---------------------------------------------------------------------
        call cloud_opt_prop_tg_lw (Cld_spec%tau, Cld_spec%liq_frac,  &
                                   emcld)

!---------------------------------------------------------------------
!    assign the emissivity value returned to both the random and max-
!    imum overlap case. the value will only be used when there is a
!    non-zero cloud fraction of a particular type.
!-------------------------------------------------------------------
        do k=1,size(Cldrad_props%emrndlw,3)
          do j=1,size(Cldrad_props%emrndlw,2)
            do i=1,size(Cldrad_props%emrndlw,1)
              Cldrad_props%emrndlw(i,j,k,:,1) = emcld(i,j,k) 
              Cldrad_props%emmxolw(i,j,k,:,1) = emcld(i,j,k) 
            end do
          end do
        end do
      endif  
 
!---------------------------------------------------------------------


end subroutine obtain_bulk_lw_diag




!#####################################################################

! <SUBROUTINE NAME="obtain_bulk_sw_diag">
!  <OVERVIEW>
!    obtain_bulk_sw_diag defines bulk shortwave cloud radiative
!    properties for the gordon diag cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_sw_diag defines bulk shortwave cloud radiative
!    properties for the gordon diag cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_diag (is, ie, js, je, cosz, Cld_spec,  &   
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
!      cosz         cosine of the zenith angle  [ dimensionless ]
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
subroutine obtain_bulk_sw_diag (is, ie, js, je, cosz, Cld_spec,  &   
                                Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_diag defines bulk shortwave cloud radiative 
!    properties for the gordon diag cloud scheme.
!---------------------------------------------------------------------

integer,                      intent(in)    ::  is, ie, js, je
real, dimension(:,:),         intent(in)    ::  cosz
type(cld_specification_type), intent(in   ) ::  Cld_spec
type(cldrad_properties_type), intent(inout) ::  Cldrad_props

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle  [ dimensionless ]
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

!-------------------------------------------------------------------
!   local variables:


      real, dimension (size(Cldrad_props%cirabsw,1),                 &
                       size(Cldrad_props%cirabsw,2),                 &
                       size(Cldrad_props%cirabsw,3)) ::  cuvab

      real, dimension (size(Cldrad_props%cirabsw,1),                 &
                       size(Cldrad_props%cirabsw,2)) :: qmix_kx

      logical, dimension (size(Cldrad_props%cirabsw,1),    &
                          size(Cldrad_props%cirabsw,2), &
                          size(Cldrad_props%cirabsw,3)) ::  direct

      real        ::    taucum
      integer     ::    max_cld, ierr, kcld
      integer     ::    i, j, k, kk

!-------------------------------------------------------------------
!   local variables:
!
!        cuvab      absorptivity of clouds in the visible frequency
!                   bands [ nondimensional ]
!        qmix_kx    mixing ratio at the lowest model level
!                   [ nondimensional ]
!        direct     logical variable indicating whether solar beam is
!                   treated as direct or diffuse
!        taucum     sum of cloud extinction coefficients from model
!                   top to the current level; if taucum is > taucrit
!                   then the solar beam is considered to be diffuse 
!                   at lower model levels [ dimensionless ]
!        max_cld    maximum number of clouds in any physics window 
!                   column
!        ierr       error flag
!        kcld       next model layer to check for the presence of cloud
!                   in the calculation of taucum. this will not be the
!                   next lower level when multi-layer clouds are present
!        i,j,k,kk   do-loop indices
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    define the nature of the solar beam at the top of atmosphere.
!---------------------------------------------------------------------
      direct(:,:,:) = .true.

!----------------------------------------------------------------------
!    in each column, integrate downward, defining the nature of the 
!    solar beam at each model level, and summing the cloud optical 
!    depths. once the optical depth exceeds taucrit, the solar beam is 
!    assigned diffuse properties rather than the properties of a direct 
!    solar beam. the logical variable direct is defined to indicate how
!    the beam is to be treated at each level.
!---------------------------------------------------------------------
      do j=1,size(Cld_spec%tau,2)
        do i=1,size(Cld_spec%tau,1)
          kcld = 1
          taucum = 0.
          do k=1,size(Cld_spec%tau,3)
            if (k >= kcld) then

!---------------------------------------------------------------------
!    once taucrit is exceeded, mark all lower levels as being diffuse,
!    and begin the next column.
!----------------------------------------------------------------------
              if (taucum > taucrit) then
                do kk=k,size(Cld_spec%tau,3)
                  direct(i,j,kk) = .false.
                end do
                exit
              endif

!---------------------------------------------------------------------
!    when multi-layer clouds are encountered, the optical depth for the
!    cloud must be added to taucum only once. kcld is defined as the 
!    first level below the current cloud at which to again start check-
!    ing taucum vs taucrit.
!---------------------------------------------------------------------
              if (Cld_spec%cld_thickness(i,j,k) > 0) then
                taucum = taucum + Cld_spec%tau(i,j,k,1)
                kcld = kcld + Cld_spec%cld_thickness(i,j,k)
              else
                kcld = kcld + 1
              endif
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------
!    determine if any clouds are present in the window.
!---------------------------------------------------------------------
      max_cld = maxval(Cld_spec%ncldsw)

!--------------------------------------------------------------------
!    define the absorptivity in the visible spectrum to be 0.0.
!--------------------------------------------------------------------
      cuvab = 0.

!---------------------------------------------------------------------
!    if cloud is present, define cloud radiative properties. otherwise 
!    the default values corresponding to no clouds will be used.
!---------------------------------------------------------------------
      if (max_cld > 0) then

!--------------------------------------------------------------------
!    obtain the properly averaged (or instantaneous) mixing ratio at
!    the lowest model level. it will be passed to cloud_opt_prop_tg_sw
!    and used to calculate anomalous absorption.
!---------------------------------------------------------------------
        call diag_cloud_avg2 (is, js, qmix_kx, ierr)

!---------------------------------------------------------------------
!    call cloud_opt_prop_tg_sw to obtain short-wave cloud radiative 
!    properties. 
!---------------------------------------------------------------------
        call cloud_opt_prop_tg_sw (Cld_spec%liq_frac, Cld_spec%tau, &
                                   direct,  qmix_kx, cosz,   &
                                   Cldrad_props%cvisrfsw, &
                                   Cldrad_props%cirrfsw, &
                                   cuvab, Cldrad_props%cirabsw)
      endif  ! (max_cld > 0)

!--------------------------------------------------------------------




end subroutine obtain_bulk_sw_diag



!####################################################################


       end module diag_clouds_W_mod




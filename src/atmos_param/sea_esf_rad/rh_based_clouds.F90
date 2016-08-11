!FDOC_TAG_GFDL

                 module rh_based_clouds_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!           module which defines cloud locations
!                     based on model relative humidity
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use mpp_mod,           only: input_nml_file
use fms_mod,           only: fms_init, open_namelist_file, mpp_pe, &
                             mpp_root_pe, stdlog,  &
                             write_version_number, file_exist, & 
                             check_nml_error, error_mesg,   &
                             FATAL, close_file
use rh_clouds_mod,     only: rh_clouds_avg      
use rad_utilities_mod, only: rad_utilities_init, &
                             cldrad_properties_type, &
                             cld_specification_type
use constants_mod,     only: radian
                                 

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!           module which defines cloud locations
!                     based on model relative humidity
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: rh_based_clouds.F90,v 19.0 2012/01/06 20:23:01 fms Exp $'
  character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          rh_based_clouds_init,  &
          rh_clouds_amt,  &
          obtain_bulk_lw_rh, obtain_bulk_sw_rh, &
          rh_based_clouds_end, &
          cldalb, albcld_lw, albcld_sw




!---------------------------------------------------------------------
!-------- namelist  ---------


character(len=8)             :: cirrus_cld_prop_form  = 'full'

!    logical variables derived from namelist input

logical                      :: do_part_black_cirrus=.false.
logical                      :: do_full_black_cirrus=.false.






namelist /rh_based_clouds_nml /     &
       cirrus_cld_prop_form


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!     define radiative properties for low, middle and high clouds
!     indices 1, 2 and 3, respectively).
!     cldem     : infrared emissivity
!     crfvis    : visible band reflectivity
!     crfir     : near-ir band reflectivity
!     cabir     : near-ir band absorptivity
!     crz       : cloud fraction
!--------------------------------------------------------------------
 
integer, parameter             :: NOFCLDS_SP=3  

real, dimension(NOFCLDS_SP)    ::                                  &
  crfvis,     crfir,     cabir,  &
                                  cldem
real                           :: crz


data crfvis / 0.59E+00, 0.45E+00, 0.21E+00 /
data crfir  / 0.59E+00, 0.45E+00, 0.21E+00 /
data cabir  / 0.40E+00, 0.30E+00, 0.04E+00 /
data cldem  / 1.0E+00, 1.0E+00, 1.0E+00 /  
data crz    / 1.00 /                     

!---------------------------------------------------------------------

!--------------------------------------------------------------------
!     these arrays define the cloud reflectivities as a function of 
!     zenith angle and the radiation band (visible and infrared) for 
!     high, middle and low clouds.
!     NREFL_BDS = number of radiative bands over which reflectivities
!              are provided.
!     NANGS     = number of zenith angles at which reflectivity values
!              are given.
!--------------------------------------------------------------------
 
integer, parameter                    ::  NANGS=17
integer, parameter                    ::  NREFL_BDS=2
real, dimension(NANGS,NREFL_BDS)      ::  albch, albcm, albcl 
 

!---------------------------------------------------------------------
!     albedos for high clouds at zenith angles from 0-80 deg. at 5 deg.
!     intervals for 1) visible and 2) infrared radiation.
!---------------------------------------------------------------------

data albch /      &
                 .04,.05,.05,.05,.06,.06,.07,.07,.08,.11,.13,.16,.21,  &
                 .28,.39,.48,.61,                                     &
                 .04,.05,.05,.05,.06,.06,.07,.07,.08,.10,.11,.14,.19, &
                 .26,.35,.44,.55 /

!---------------------------------------------------------------------
!     albedos for middle clouds at zenith angles from 0-80 deg. at 5 deg
!     intervals for 1) visible and 2) infrared radiation.
!----------------------------------------------------------------------

data albcm /     &
               .18,.18,.19,.20,.21,.23,.24,.26,.29,.33,.37,.42,.47,  &
               .55,.64,.71,.79,                                      &
               .14,.14,.15,.16,.17,.18,.18,.20,.23,.25,.29,.32,.37, &
               .43,.50,.55,.61 /

!-----------------------------------------------------------------------
!     albedos for low clouds at zenith angles from 0-80 deg. at 5 deg
!     intervals for 1) visible and 2) infrared radiation.
!-----------------------------------------------------------------------

data albcl /     &
                .50,.50,.51,.51,.52,.53,.54,.56,.58,.62,.65,.67,.69, &
                .73,.78,.82,.86,                                     &
                .42,.42,.43,.43,.44,.45,.46,.48,.50,.52,.55,.57,.59, &
                .63,.66,.70,.74 /


 
!-------------------------------------------------------------------
!     this array defines the zenith angle dependent albedo for each of 
!     the different cloud types (NOFCLDS_SP) for each of the radiative 
!     bands (NREFL_BDS). currently NREFL_BDS are the visible and the
!     infrared.
!-------------------------------------------------------------------
 
real, dimension (NOFCLDS_SP,NREFL_BDS)  ::  zza 
 
!--------------------------------------------------------------------
!     these variables define the boundaries (in sigma coordinates) 
!     between high and middle and middle and low clouds at the poles
!     and at the equator. 
!--------------------------------------------------------------------

real    :: cldhp = 0.7E+00
real    :: cldhe = 0.4E+00
real    :: cldmp = 0.85E+00
real    :: cldme = 0.7E+00
!-----------------------------------------------------------------

! cloud is present when relative humidity >= rh_crit, which varies liearly
!   in sigma from rh_crit_top at sigma = 0 to rh_crit_bot at sigma = 1

real :: rh_crit_bot    = 1.00
real :: rh_crit_top    = 0.90

logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!----------------------------------------------------------------------




                           contains 


! <SUBROUTINE NAME="rh_based_clouds_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rh_based_clouds_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rh_based_clouds_init 



!--------------------------------------------------------------------
     integer :: unit, ierr, io, logunit


      if (module_is_initialized) return
      call fms_init
      call rad_utilities_init
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=rh_based_clouds_nml, iostat=io)
      ierr = check_nml_error(io,"rh_based_clouds_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file (                          )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=rh_based_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'rh_based_clouds_nml')
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
           write (logunit, nml=rh_based_clouds_nml)

!--------------------------------------------------------------------
! define the "blackness" of the cirrus clouds. cirrus clouds are either 
! "part black" with emissivity of 0.6 or are "full black" with emis-
! sivity of 1.0
!--------------------------------------------------------------------
        if (trim(cirrus_cld_prop_form) == 'part') then
          do_part_black_cirrus = .true.
        else if (trim(cirrus_cld_prop_form) == 'full') then
          do_full_black_cirrus = .true.
        else
          call error_mesg( 'cloudrad_package_init',  &
                ' cirrus_cld_prop_form is not an acceptable value.', & 
                                                      FATAL)
        endif


       module_is_initialized = .true.

end subroutine rh_based_clouds_init

!####################################################################

! <SUBROUTINE NAME="rh_based_clouds_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rh_based_clouds_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rh_based_clouds_end

!----------------------------------------------------------------------
!    rh_based_clouds_end is the destructor for rh_based_cloouds_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine rh_based_clouds_end



!######################################################################

! <SUBROUTINE NAME="rh_clouds_amt">
!  <OVERVIEW>
!    rh_clouds_amt defines the location, amount (cloud fraction), number
!    and type (hi, mid, low) of clouds present on the model grid.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    rh_clouds_amt defines the location, amount (cloud fraction), number
!    and type (hi, mid, low) of clouds present on the model grid.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rh_clouds_amt (is, ie, js, je, press, lat, Cld_spec)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   is,ie,js,je  starting/ending subdomain i,j indices of data in
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
!  <IN NAME="press" TYPE="real">
!      press        pressure at model levels (1:nlev), surface
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
! 
!  </IN>
!  <IN NAME="lat" TYPE="real">
!      lat          latitude of model points  [ radians ]
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
subroutine rh_clouds_amt (is, ie, js, je, press, lat, Cld_spec)

!----------------------------------------------------------------------
!    rh_clouds_amt defines the location, amount (cloud fraction), number
!    and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer,                      intent(in)    ::  is, ie, js, je
real,    dimension(:,:,:),    intent(in)    ::  press
real,    dimension(:,:),      intent(in)    ::  lat                    
type(cld_specification_type), intent(inout) ::  Cld_spec       

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      lat          latitude of model points  [ radians ]
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
 
!-------------------------------------------------------------------
!   local variables:
 
      real, dimension (size(press,1), size(press,2),    &
                       size(press,3)-1)   ::  ccover, rh, sigma

      real, dimension (size(lat,1), size(lat,2)) ::     &
                                                 cldhm, cldml, rh_crit

      logical   ::  inside_max_ovlp_cld
      integer   ::  kmax 
      integer   ::  ierr
      integer   ::  k, j, i

!----------------------------------------------------------------------
!   local variables:
!
!           ccover       cloud fraction in a grid box [ dimensinless ]
!           rh           relative humidity [ dimensionless ]
!           sigma        ratio of pressure to surface pressure 
!                        [ dimensionless ]
!           cldhm        sigma value to use as the transition between
!                        middle and high clouds [ dimensionless ]
!           cldml        sigma value to use as the transition between
!                        low and middle clouds [ dimensionless ]
!           rh_crit      value of relative humidity at which clouds 
!                        are assumed present, as a function of sigma
!                        [ dimensionless ]
!           inside_max_ovlp_cld
!                        logical flag indicating whether the model grid
!                        box under consideration is part of a maximum
!                        overlap cloud with its base below the current
!                        level
!           kmax         number of model layers
!           ierr         error flag returned from rh_clouds_mod, if 
!                        non-zero indicates that all or some of the 
!                        relative humidity field was not retrievable
!           i,j,k        do loop indices
!
!--------------------------------------------------------------------- 


!--------------------------------------------------------------------
!    define the number of model layers.
!----------------------------------------------------------------------
      kmax = size (Cld_spec%camtsw,3)

!---------------------------------------------------------------------
!    define the sigma values marking the transitions between high,
!    middle and low clouds.  it must be calculated in this routine 
!    because model latitudes are not available during initialization
!    phase in the bgrid model core.
!---------------------------------------------------------------------
      cldhm(:,:) = cldhp + (90.0E+00-abs(lat(:,:)*   &
                   radian))*(cldhe-cldhp)/90.0E+00
      cldml(:,:) = cldmp + (90.0E+00-abs(lat(:,:)*    &
                   radian))*(cldme-cldmp)/90.0E+00

!---------------------------------------------------------------------- 
!    call rh_clouds_avg to obtain the appropriate array of relative 
!    humidity to use in defining cloud locations. this may be the inst-
!    antaneous field from the last time step or the average field in the
!    interval from the previous call to rh_clouds_avg.
!----------------------------------------------------------------------
      call rh_clouds_avg (is, js, rh, ierr)

!----------------------------------------------------------------------
!    if relative humidity data was present in rh_clouds_mod and values 
!    have been returned, determine those grid points where the relative
!    humidity exceeds the critical value for that grid box. define a
!    non-zero cloud fraction for those boxes; other boxes are given a
!    zero cloud fraction. note that cloud is not allowed to be present
!    in the topmost layer.
!----------------------------------------------------------------------
      if (ierr ==  0) then
        ccover(:,:,1) = 0.0
        do k=2,kmax 
          do j=1,size(press,2)
            do i=1,size(press,1)
              sigma(i,j,k) = press(i,j,k)/press(i,j,kmax+1)
              rh_crit(i,j) = rh_crit_top + sigma(i,j,k)*  &
                                           (rh_crit_bot - rh_crit_top)
              if (rh(i,j,k) >= rh_crit(i,j))  then 
                ccover(i,j,k) = crz
              else 
                ccover(i,j,k) = 0.0
              endif     
            end do
          end do
        end do

!---------------------------------------------------------------------
!    if relative humidity data could not be returned, set the clouds
!    to zero everywhere. this is a valid occurrence on the first
!    time step of a run, when radiation is calculated prior to the
!    moist_processes physics, and so no relative humidity data is 
!    present. if it occurs after the first step, an error is likely
!    present.
!---------------------------------------------------------------------
      else
        ccover(:,:,:) = 0.0
      endif

!---------------------------------------------------------------------
!    define the cloud specification arrays used by rh_clouds -- cloud 
!    fractions for random and maximally overlapped clouds and total 
!    clouds, number of each type of cloud in each column, and flags to
!    indicate whether a given cloud is to be given high, middle or low 
!    cloud properties.
!---------------------------------------------------------------------
      do j=1,size(press,2)
        do i=1,size(press,1)

!--------------------------------------------------------------------
!    define an indicator as to whether the current level is part of a
!    maximum-overlap cloud with a base at a lower level. the lowest
!    level obviously is not.
!---------------------------------------------------------------------
          inside_max_ovlp_cld = .false.

!---------------------------------------------------------------------
!    move upwards from the surface to level 2, since clouds are not
!    allowed at the topmost model level.
!----------------------------------------------------------------------
          do k=kmax,2,-1

!---------------------------------------------------------------------
!    if inside_max_ovlp_cld is .false., then this level is not part of 
!    a multi-layer cloud with a cloud base at a lower level, either 
!    because it does not contain cloud or because it contains a single 
!    level cloud, or because it is the lowest model level.
!---------------------------------------------------------------------
            if (.not. inside_max_ovlp_cld) then

!---------------------------------------------------------------------
!    if there is cloud at this level, determine if cloud exists at the
!    level above. if it doesn't then this is a one layer cloud and it
!    is counted as a random overlap cloud. if cloud does exist at the
!    level above, then this is the lowest level of a multi-layer cloud,
!    and it is denoted as a maximum overlap cloud. the flag 
!    inside_max_ovlp_cld is set to .true. to indicate that the next 
!    level is within a maximum overlap cloud.
!---------------------------------------------------------------------
              if (ccover(i,j,k) .NE. 0.0E+00) then
                if (ccover(i,j,k-1) .EQ. 0.0E+00 ) then
                  Cld_spec%nrndlw(i,j) = Cld_spec%nrndlw(i,j) + 1
                  Cld_spec%crndlw(i,j,k) = ccover(i,j,k)
                else
                  inside_max_ovlp_cld = .true.
                  Cld_spec%cmxolw(i,j,k) = ccover(i,j,k)
                endif
              endif
!--------------------------------------------------------------------
!    if inside_max_ovlp_cld is .true., then the current level is 
!    contained within a multi-layer cloud (max overlap) with a cloud 
!    base at a lower level. define the cloud on this level as maximum 
!    overlap. the counter of maximum overlap clouds is incremented when 
!    either the model top (k = 2, since no clouds allowed at level 1) 
!    is reached, or the level above does not contain cloudiness. if the 
!    level above does not contain clouds, then the flag 
!    inside_max_ovlp_cld is set .false. to indicate that the level above
!    is not part of a maximum overlap cloud with a cloud base at a 
!    lower level.
!--------------------------------------------------------------------
            else
              Cld_spec%cmxolw(i,j,k) = ccover(i,j,k)
              if (k .EQ. 2) then
                Cld_spec%nmxolw(i,j) = Cld_spec%nmxolw(i,j) + 1
              endif
              if (ccover(i,j,k-1) .EQ. 0.0E+00) then
                Cld_spec%nmxolw(i,j) = Cld_spec%nmxolw(i,j) + 1
                inside_max_ovlp_cld = .false.
              endif
            endif

!---------------------------------------------------------------------
!    define the cloud fraction seen by the shortwave radiation as the 
!    sum of the random and maximum overlap cloud fractions.
!---------------------------------------------------------------------
            Cld_spec%camtsw(i,j,k) = Cld_spec%cmxolw(i,j,k) +  &
                                     Cld_spec%crndlw(i,j,k)

!----------------------------------------------------------------------
!    if cloud is present, it must be designated as high, middle or low
!    so that the grid box can be given the radiation characteristics of
!    that cloud type. the arrays cldhm and cldml define the sigma
!    boundaries of the different cloud types.
!----------------------------------------------------------------------
            if (Cld_spec%camtsw(i,j,k) > 0.0) then
              if (sigma(i,j,k) <= cldhm(i,j) ) then
                Cld_spec%hi_cloud(i,j,k) = .true.
              else if (sigma(i,j,k) > cldhm(i,j)  .and.  &
                       sigma(i,j,k) < cldml(i,j) ) then
                Cld_spec%mid_cloud(i,j,k) = .true.
              else if (sigma(i,j,k) >= cldml(i,j) ) then
                Cld_spec%low_cloud(i,j,k) = .true.
              endif
            endif
          end do

!-------------------------------------------------------------------
!      define the total number of clouds present in each column.
!-------------------------------------------------------------------
          Cld_spec%ncldsw(i,j) = Cld_spec%nmxolw(i,j) +   &
                                 Cld_spec%nrndlw(i,j)
        end do
      end do

!--------------------------------------------------------------------


end subroutine rh_clouds_amt 





!####################################################################

! <SUBROUTINE NAME="obtain_bulk_lw_rh">
!  <OVERVIEW>
!    obtain_bulk_lw_rh defines bulk longwave cloud radiative
!    properties for the rh cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_lw_rh defines bulk longwave cloud radiative
!    properties for the rh cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_rh (is, ie, js, je, Cld_spec, Cldrad_props)
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
subroutine obtain_bulk_lw_rh (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_rh defines bulk longwave cloud radiative 
!    properties for the rh cloud scheme.
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


!------------------------------------------------------------------
!     call albcld_lw to define long-wave cloud emissivities.
!-------------------------------------------------------------------
        call albcld_lw (Cld_spec%hi_cloud, Cld_spec%mid_cloud,  &
                        Cld_spec%low_cloud, Cld_spec%cmxolw,  &
                        Cld_spec%crndlw,    &
                        Cldrad_props%emmxolw(:,:,:,:,1),  &
                        Cldrad_props%emrndlw(:,:,:,:,1))


end subroutine obtain_bulk_lw_rh



!######################################################################

! <SUBROUTINE NAME="obtain_bulk_sw_rh">
!  <OVERVIEW>
!    obtain_bulk_sw_rh defines bulk shortwave cloud radiative
!    properties for the rh cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_sw_rh defines bulk shortwave cloud radiative
!    properties for the rh cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_rh (is, ie, js, je, cosz, Cld_spec,   &
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
subroutine obtain_bulk_sw_rh (is, ie, js, je, cosz, Cld_spec,   &
                              Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_rh defines bulk shortwave cloud radiative 
!    properties for the rh cloud scheme.
!---------------------------------------------------------------------
 
integer,                      intent(in)    :: is, ie, js, je
real,    dimension(:,:),      intent(in)    :: cosz
type(cld_specification_type), intent(in   ) :: Cld_spec
type(cldrad_properties_type), intent(inout) :: Cldrad_props

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      cosz         cosine of the zenith angle [ dimensionless ]
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
!    local variables:

      integer      ::  j, i  ! do-loop indices
 
!-------------------------------------------------------------------
!    define the bulk sw cloud radiative properties at each grid point 
!    which contains cloud.
!-------------------------------------------------------------------
      do j=1,size(Cld_spec%ncldsw,2)          
        do i=1,size(Cld_spec%ncldsw,1)          
          if (Cld_spec%ncldsw(i,j) > 0 ) then

!-----------------------------------------------------------------------
!    call cldalb to define the zenith angle dependent cloud visible
!    and infrared reflectivities for high, middle and low clouds in
!    each model column containing cloud.
!-----------------------------------------------------------------------
            call cldalb (cosz(i,j))

!-----------------------------------------------------------------------
!    call albcld_sw to assign the zenith-angle-dependent visible and 
!    infrared reflectivities and infrared absorptivities to the model 
!    clouds in the current column.
!-----------------------------------------------------------------------
            call albcld_sw (i, j, Cld_spec%hi_cloud, &
                            Cld_spec%mid_cloud, Cld_spec%low_cloud,  &
                            Cld_spec%camtsw, Cld_spec%cmxolw,   &
                            Cld_spec%crndlw, Cldrad_props%cvisrfsw, &
                            Cldrad_props%cirrfsw, Cldrad_props%cirabsw)
          endif
        end do
      end do

!----------------------------------------------------------------------
 

 end subroutine obtain_bulk_sw_rh



!####################################################################

! <SUBROUTINE NAME="cldalb">
!  <OVERVIEW>
!     cldalb calculates a zenith angle dependency for the cloud albedos.
!     the cloud albedos are interpolated using data adapted from fritz
!     (1954).  the solar zenith angle is the only input required.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     cldalb calculates a zenith angle dependency for the cloud albedos.
!     the cloud albedos are interpolated using data adapted from fritz
!     (1954).  the solar zenith angle is the only input required.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cldalb (zenith)
!
!  </TEMPLATE>
!  <IN NAME="zenith" TYPE="real">
!   zenith angle
!  </IN>
! </SUBROUTINE>
!
subroutine cldalb (zenith)

!---------------------------------------------------------------------
!     cldalb calculates a zenith angle dependency for the cloud albedos.
!     the cloud albedos are interpolated using data adapted from fritz 
!     (1954).  the solar zenith angle is the only input required.
!-----------------------------------------------------------------------

real, intent(in)           ::  zenith

!--------------------------------------------------------------------
      real                 :: zangle, remain
      integer              :: nband, indx

!--------------------------------------------------------------------
!     define zenith angle in degrees. for original Skyhi results, use 
!     a zenith angle specified as 60.00001 and the skyhi albedo values.
!-----------------------------------------------------------------------

        zangle = ACOS(zenith)*radian

!-----------------------------------------------------------------------
!     define reflectivities for each cloud level.
!-----------------------------------------------------------------------
        if (zangle .GE. 80.0E+00) then
!-----------------------------------------------------------------------
!     if zenith angle is greater than 80 degrees, define the reflect-
!     ivities as those values in the table at 80 degrees (last entry).
!-----------------------------------------------------------------------
          do nband=1,NREFL_BDS
            zza(3,nband) = albch(NANGS,nband)
            zza(2,nband) = albcm(NANGS,nband)
            zza(1,nband) = albcl(NANGS,nband)
          end do  
        else
!-----------------------------------------------------------------------
!     if zenith angle is less than 80 degrees, interpolate albedos from 
!     tables for each cloud level.
!-----------------------------------------------------------------------
          indx   = IFIX(zangle/5.0E+00) + 1
          remain = AMOD(zangle, 5.0E+00)
          do nband=1,NREFL_BDS
            zza(3,nband) = albch(indx,nband) + (remain/5.0E+00)*  &
                           (albch(indx+1,nband) - albch(indx,nband))
            zza(2,nband) = albcm(indx,nband) + (remain/5.0E+00)*  &
                           (albcm(indx+1,nband) - albcm(indx,nband))
            zza(1,nband) = albcl(indx,nband) + (remain/5.0E+00)*   &
                           (albcl(indx+1,nband) - albcl(indx,nband))
          end do 
        endif



end subroutine cldalb




!##################################################################

! <SUBROUTINE NAME="albcld_lw">
!  <OVERVIEW>
!     albcld_lw computes the lw cloud emissivities. This calculation is
!     based on sigma and cloud thickness in the old scheme (cldht60)
!     and sigma, cloud thickness and latitude in the new scheme
!     (cldht93).
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     albcld_lw computes the lw cloud emissivities. This calculation is
!     based on sigma and cloud thickness in the old scheme (cldht60)
!     and sigma, cloud thickness and latitude in the new scheme
!     (cldht93).
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call albcld_lw(hi_cloud, mid_cloud, low_cloud,       &
!                cmxolw, crndlw, emmxolw, emrndlw)
!
!  </TEMPLATE>
!  <IN NAME="hi_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="mid_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="low_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="cmxolw" TYPE="real">
! 
!  </IN>
!  <IN NAME="crndlw" TYPE="real">
! 
!  </IN>
!  <INOUT NAME="emmxolw" TYPE="real">
! 
!  </INOUT>
!  <INOUT NAME="emrndlw" TYPE="real">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine albcld_lw(hi_cloud, mid_cloud, low_cloud,       &
             cmxolw, crndlw, emmxolw, emrndlw)

!-----------------------------------------------------------------------
!     albcld_lw computes the lw cloud emissivities. This calculation is 
!     based on sigma and cloud thickness in the old scheme (cldht60) 
!     and sigma, cloud thickness and latitude in the new scheme 
!     (cldht93).
!-----------------------------------------------------------------------

real, dimension(:,:,:),   intent(in)    :: cmxolw, crndlw
real, dimension(:,:,:,:), intent(inout) :: emmxolw, emrndlw
logical, dimension(:,:,:),intent(in)    :: hi_cloud, mid_cloud,   &
   low_cloud

!---------------------------------------------------------------------
!    local variables
!---------------------------------------------------------------------

       integer  ::  i, j, k, kk
      integer  :: israd, ierad, jsrad, jerad, ksrad, kerad
      israd= 1
      jsrad= 1
      ksrad= 1
      ierad = size (cmxolw,1)
      jerad = size (cmxolw,2)
      kerad = size (cmxolw,3)

!-----------------------------------------------------------------------
!     compute the emissivities for each cloud in the column. 
!-----------------------------------------------------------------------

         do k=KSRAD,KERAD
   do j=JSRAD,JERAD
     do i=ISRAD,IERAD
               if ((cmxolw(i,j,k) + crndlw(i,j,k) ) > 0.0) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that thick cloud properties are deter-
!     mined by the height of the base of the thick cloud, so that if 
!     there are two adjacent cirrus level clouds, they are assigned
!     cirrus cloud properties, incontrast to the cldht60 treatment.
!-----------------------------------------------------------------------
                 if (cmxolw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick high cloud.
!-----------------------------------------------------------------------
                   if (hi_cloud(i,j,k)) then
                       emmxolw(i,j,k,:) = cldem(3)*0.6
!-----------------------------------------------------------------------
!     case of a thick middle cloud.
!-----------------------------------------------------------------------
                   else if (mid_cloud(i,j,k)) then
             emmxolw(i,j,k,:) = cldem(2)
!-----------------------------------------------------------------------
!     case of a thick low cloud.
!-----------------------------------------------------------------------
                   else if (low_cloud (i,j,k)) then
             emmxolw(i,j,k,:) = cldem(1)
                   endif
                 endif
         if (crndlw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
                   if (hi_cloud(i,j,k)) then
             if (do_full_black_cirrus) then
                       emrndlw(i,j,k,:) = cldem(3)
             else if (do_part_black_cirrus) then
               emrndlw(i,j,k,:) = 0.6E+00*cldem(3)
             endif
!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
                   else if (mid_cloud(i,j,k)) then 
             emrndlw(i,j,k,:) = cldem(2)
!-----------------------------------------------------------------------
!     case of a thin low cloud.
!-----------------------------------------------------------------------
                   else if (low_cloud(i,j,k)) then
             emrndlw(i,j,k,:) = cldem(1)
           endif
                 endif
       endif
             end do     
           end do
         end do
!------------------------------------------------------------------
!! for fms formulation, set thick cloud properties based on cloud base,
!  even if some cloud layers extend out of the cloud base type region
!-------------------------------------------------------------------

           do k=KERAD,KSRAD+1,-1
     do j=JSRAD,JERAD
       do i=ISRAD,IERAD
                 if (cmxolw(i,j,k) /= 0.0) then
                   kk = k-1
                   if (cmxolw(i,j,kk) /= 0.0) then
                     emmxolw(i,j,kk,:) = emmxolw(i,j,k,:)
                   endif
                 endif
               end do
             end do
           end do
     
 
end subroutine albcld_lw


!####################################################################

! <SUBROUTINE NAME="albcld_sw">
!  <OVERVIEW>
!     albcld_sw computes the cloud albedos. This calculation is based on
!     sigma and cloud thickness in the old scheme (cldht60) and sigma,
!     cloud thickness  and latitude in the new scheme (cldht93).
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     albcld_sw computes the cloud albedos. This calculation is based on
!     sigma and cloud thickness in the old scheme (cldht60) and sigma,
!     cloud thickness  and latitude in the new scheme (cldht93).
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call albcld_sw(i,j, hi_cloud, mid_cloud, low_cloud,         &
!                camtsw, cmxolw, crndlw, cvisrfsw, cirrfsw, cirabsw)
!
!  </TEMPLATE>
!  <INOUT NAME="i" TYPE="real">
! 
!  </INOUT>
!  <IN NAME="j" TYPE="integer">
! 
!  </IN>
!  <IN NAME="hi_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="mid_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="low_cloud" TYPE="logical">
! 
!  </IN>
!  <IN NAME="camtsw" TYPE="real">
! 
!  </IN>
!  <IN NAME="cmxolw" TYPE="real">
! 
!  </IN>
!  <IN NAME="crndlw" TYPE="real">
! 
!  </IN>
!  <INOUT NAME="cvisrfsw" TYPE="real">
! 
!  </INOUT>
!  <INOUT NAME="cirrfsw" TYPE="real">
! 
!  </INOUT>
!  <INOUT NAME="cirabsw" TYPE="real">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine albcld_sw(i,j, hi_cloud, mid_cloud, low_cloud,         &
     camtsw, cmxolw, crndlw, cvisrfsw, cirrfsw, cirabsw)

!-----------------------------------------------------------------------
!     albcld_sw computes the cloud albedos. This calculation is based on
!     sigma and cloud thickness in the old scheme (cldht60) and sigma, 
!     cloud thickness  and latitude in the new scheme (cldht93).
!-----------------------------------------------------------------------

!real, dimension(:,:,:),   intent(inout) :: camtsw, cmxolw, crndlw
real, dimension(:,:,:),   intent(in) :: camtsw, cmxolw, crndlw
real, dimension(:,:,:), intent(inout) :: cvisrfsw, cirrfsw, cirabsw
logical, dimension(:,:,:),intent(in)    :: hi_cloud, mid_cloud,   &
   low_cloud
integer,                  intent(in)    :: i, j
!---------------------------------------------------------------------

       integer  ::  k,  kk

      integer  :: israd, ierad, jsrad, jerad, ksrad, kerad

      israd = 1
      jsrad = 1
      ksrad = 1
      ierad = size (camtsw,1)
      jerad = size (camtsw,2)
      kerad = size (camtsw,3)
!-----------------------------------------------------------------------
!     compute the reflectivities and absorptivities for each cloud in 
!     the column. cldhm and cldml are sigma levels which serve as 
!     boundaries between low, middle and high clouds. 
!-----------------------------------------------------------------------
         do k=KSRAD+1,KERAD
           if (camtsw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that thick cloud properties are deter-
!     mined by the height of the base of the thick cloud, so that if 
!     there are two adjacent cirrus level clouds, they are assigned
!     cirrus cloud properties.
!-----------------------------------------------------------------------
             if (cmxolw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick high cloud.
!-----------------------------------------------------------------------
               if (hi_cloud(i,j,k)) then
                 cvisrfsw(i,j,k) = zza(3,1)
                 cirrfsw(i,j,k) = zza(3,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                            cabir(3) )

!-----------------------------------------------------------------------
!     case of a thick middle cloud.
!-----------------------------------------------------------------------
               else if (mid_cloud(i,j,k)) then
                 cvisrfsw(i,j,k) = zza(2,1)
                 cirrfsw(i,j,k) = zza(2,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                            cabir(2) )

!-----------------------------------------------------------------------
!     case of a thick low cloud.
!-----------------------------------------------------------------------
               else if (low_cloud(i,j,k)) then
                 cvisrfsw(i,j,k) = zza(1,1)
                 cirrfsw(i,j,k) = zza(1,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                    cabir(1) )
               endif
             endif
     if (crndlw(i,j,k) .NE. 0.0E+00) then

!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
               if (hi_cloud(i,j,k)) then
                 cvisrfsw(i,j,k) = zza(3,1)
                 cirrfsw(i,j,k) = zza(3,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                            cabir(3) )

!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
               else if (mid_cloud(i,j,k))  then
                 cvisrfsw(i,j,k) = zza(2,1)
                 cirrfsw(i,j,k) = zza(2,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                            cabir(2) )

!-----------------------------------------------------------------------
!     case of a thin low cloud.
!-----------------------------------------------------------------------
               else if (low_cloud(i,j,k)) then
                 cvisrfsw(i,j,k) = zza(1,1)
                 cirrfsw(i,j,k) = zza(1,2)
                   cirabsw(i,j,k) = MIN(0.99-cirrfsw(i,j,k), &
                                            cabir(1) )
               endif
             endif
           endif
         end do     

!------------------------------------------------------------------
!! for fms formulation, set thick cloud properties based on cloud base,
!  even if some cloud layers extend out of the cloud base type region
!-------------------------------------------------------------------
           do k=KERAD,KSRAD+1,-1
             if (cmxolw(i,j,k) /= 0.0) then
               kk = k-1
               if (cmxolw(i,j,kk) /= 0.0) then
                 cvisrfsw(i,j,kk) = cvisrfsw(i,j,k)
                 cirrfsw(i,j,kk) = cirrfsw(i,j,k)
                 cirabsw(i,j,kk) = cirabsw(i,j,k)
               endif
             endif
           end do
     
 
end subroutine albcld_sw



       end module rh_based_clouds_mod



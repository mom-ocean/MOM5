                 module zetac_clouds_W_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  sak
! </REVIEWER>
! <OVERVIEW>
!    zetac_clouds_W_mod obtains the cloud specification variables
!    for the zetac cloud parameterization from microphys_cloud_mod
!    and makes them available to the radiation package.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!   shared modules:
use mpp_mod,                only: input_nml_file
use fms_mod,                only: open_namelist_file, mpp_pe, &
                                  mpp_root_pe, stdlog,  fms_init, &
                                  write_version_number, file_exist, &
                                  check_nml_error, error_mesg,   &
                                  FATAL, close_file
use constants_mod,          only: GRAV, constants_init

!   shared radiation package modules:

use rad_utilities_mod,      only: rad_utilities_init, &
                                  cld_specification_type, &
                                  microphysics_type

!   cloud parameterization module:

use microphys_cloud_mod,     only: microphys_cloud_init, microphys_cloud

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    zetac_clouds_W_mod obtains the cloud specification variables
!    for the zetac cloud parameterization from microphys_cloud_mod
!    and makes them available to the radiation package.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: zetac_clouds_W.F90,v 19.0 2012/01/06 20:25:17 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          zetac_clouds_W_init, zetac_clouds_amt, zetac_clouds_W_end

!---------------------------------------------------------------------
!-------- namelist  ---------

integer   :: dummy = 0


namelist /zetac_clouds_W_nml /                      &
                                 dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


logical :: module_is_initialized = .false.  ! module is initialized ?


!----------------------------------------------------------------------
!----------------------------------------------------------------------



                              contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="zetac_clouds_W_init">
!  <OVERVIEW>
!    zetac_clouds_W_init is the constructor for zetac_clouds_W_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    zetac_clouds_W_init is the constructor for zetac_clouds_W_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zetac_clouds_W_init
!  </TEMPLATE>
! </SUBROUTINE>
!  
subroutine zetac_clouds_W_init 

!---------------------------------------------------------------------
!    zetac_clouds_W_init is the constructor for zetac_clouds_W_mod.
!---------------------------------------------------------------------

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
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call constants_init
      call microphys_cloud_init

!---------------------------------------------------------------------
!    read namelist.         
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=zetac_clouds_W_nml, iostat=io)
      ierr = check_nml_error(io,"zetac_clouds_W_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read (unit, nml=zetac_clouds_W_nml, iostat=io, end=10) 
        ierr = check_nml_error (io, 'zetac_clouds_W_nml')
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
                 write (logunit, nml=zetac_clouds_W_nml)

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


end subroutine zetac_clouds_W_init



!######################################################################
! <SUBROUTINE NAME="zetac_clouds_amt">
!  <OVERVIEW>
!    zetac_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    zetac_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. if a microphysically-based cloud parameter-
!    ization is being used, particle sizes and concentrations are also
!    provided.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zetac_clouds_amt (is, ie, js, je, z_half, z_full, land, &
!                             phalf, deltaz, Cld_spec, Lsc_microphys)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
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
subroutine zetac_clouds_amt (is, ie, js, je, z_half, z_full, land, &
                             phalf, deltaz,  &
                             Cld_spec, Lsc_microphys)

!---------------------------------------------------------------------
!    zetac_clouds_amt defines the location, amount (cloud fraction), 
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop 
!    and crystal sizes. if a microphysically-based cloud parameter-
!    ization is being used, particle sizes and concentrations are also
!    provided.
!----------------------------------------------------------------------

integer,                      intent(in)        :: is, ie, js, je
real,    dimension(:,:,:),    intent(in)        :: z_half, z_full, &
                                                   phalf, deltaz
real,    dimension(:,:),      intent(in)        :: land
type(cld_specification_type), intent(inout)     :: Cld_spec      
type(microphysics_type),      intent(inout)     :: Lsc_microphys

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ] 
!      press        pressure at model levels (1:nlev), surface 
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      temp         temperature at model levels (1:nlev), to be used
!                   in cloud calculations
!                   [ deg K ]
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

      integer :: idim, jdim, kdim
      integer :: i, j
      real, dimension (size(z_full,1), size(z_full,2),   &
                       size(z_full,3)) ::    deltap, diam_ice, diam_liq

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('zetac_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif


      idim = size (z_full,1)
      jdim = size (z_full,2)
      kdim = size (z_full,3)

      where ( Cld_spec%cloud_area > 0.0 ) 
        Cld_spec%cld_thickness = 1
      elsewhere
        Cld_spec%cld_thickness = 0
      endwhere         
 
      do i=1,idim
        do j=1,jdim
          Cld_spec%ncldsw(i,j) = sum( Cld_spec%cld_thickness(i,j,:) )
        enddo
      enddo

      Cld_spec%nrndlw = Cld_spec%ncldsw
      Cld_spec%nmxolw = 0 
      Cld_spec%camtsw = Cld_spec%cloud_area  
      Cld_spec%crndlw = Cld_spec%cloud_area
      Cld_spec%cmxolw = 0.0

      deltap = phalf(:,:,2:kdim+1)  - phalf(:,:,1:kdim)

      Cld_spec%lwp = deltap/GRAV*Cld_spec%cloud_water
      Cld_spec%iwp = deltap/GRAV*Cld_spec%cloud_ice

      call microphys_cloud (z_half, z_full, diam_liq, diam_ice)

      Lsc_microphys%size_drop = diam_liq
      Lsc_microphys%size_ice  = diam_ice
      Lsc_microphys%conc_drop = Cld_spec%lwp*1.0e3/deltaz
      Lsc_microphys%conc_ice  = Cld_spec%iwp*1.0e3/deltaz
      Lsc_microphys%cldamt    = Cld_spec%cloud_area 


!---------------------------------------------------------------------



end subroutine zetac_clouds_amt  




!####################################################################
! <SUBROUTINE NAME="zetac_clouds_W_end">
!  <OVERVIEW>
!    zetac_clouds_W_end is the destructor for zetac_clouds_W_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    zetac_clouds_W_end is the destructor for zetac_clouds_W_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zetac_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
subroutine zetac_clouds_W_end
       
!----------------------------------------------------------------------
!    zetac_clouds_W_end is the destructor for zetac_clouds_W_mod.
!----------------------------------------------------------------------
        
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('zetac_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------
 
 
end subroutine zetac_clouds_W_end



!#################################################################




                    end module zetac_clouds_W_mod


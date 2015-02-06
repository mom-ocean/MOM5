                         module radiation_diag_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Module that provides a diagnostic output file of radiation-
!    related variables in user-specified atmospheric columns for the
!    sea_esf_rad radiation package.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!  shared modules:

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, fms_init, &
                              mpp_pe, mpp_root_pe, stdlog, &
                              file_exist, write_version_number, &
                              check_nml_error, error_mesg, &
                              FATAL, close_file, &
                              open_file     
use constants_mod,      only: constants_init, radcon_mks, radian

!  shared radiation package modules:

use rad_utilities_mod,  only: Lw_control,  Sw_control,&
                              rad_utilities_init, radiative_gases_type,&
                              astronomy_type, atmos_input_type, &
                              surface_type, cldrad_properties_type, &
                              cld_specification_type,  &
                              cld_space_properties_type, &
                              lw_diagnostics_type, lw_table_type, &
                              sw_output_type, lw_output_type, &
                              Rad_control

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    radiation_diag_mod provides a diagnostic output file of radiation-
!    related variables in user-specified atmospheric columns for the
!    sea_esf_rad radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: radiation_diag.F90,v 19.0 2012/01/06 20:22:27 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!------    interfaces   ------

public    &
         radiation_diag_init,  &
         radiation_diag_driver,   &
         radiation_diag_end


private   &
         radiag


!---------------------------------------------------------------------
!---------- namelist  -----

integer, parameter                 :: MAX_PTS = 20
integer, dimension (MAX_PTS)       :: iradprt_gl=0, jradprt_gl=0
real, dimension(MAX_PTS)           :: latradprt=999., lonradprt=999.
integer                            :: num_pts_ij = 0
integer                            :: num_pts_latlon = 0


namelist / radiation_diag_nml /  &
                                  iradprt_gl, jradprt_gl, &
                                  num_pts_ij, num_pts_latlon,       &
                                  latradprt, lonradprt

!----------------------------------------------------------------------
!---  public data ---


!----------------------------------------------------------------------
!---  private data ---

!----------------------------------------------------------------------
!    bandlo and bandhi are the lower and upper frequencies defining 
!    each of the nblw longwave radiation bands in the longwave spectrum
!    (0 - 3000 cm-1).
!----------------------------------------------------------------------
real,    dimension(:), allocatable  :: bandlo, bandhi

!----------------------------------------------------------------------
!    bdlocm and bdhicm are the lower and upper frequencies defining 
!    each of the nbly frequency bands used in the exact cool-to-space
!    calculation.
! --------------------------------------------------------------------
real,    dimension(:), allocatable  :: bdlocm, bdhicm

!---------------------------------------------------------------------
!    iband is frequency band index to be used in combined band 
!    calculations. 
!---------------------------------------------------------------------
integer, dimension(:), allocatable  :: iband

!---------------------------------------------------------------------
!    deglon1 and deglat1 are the longitude and latitude of the columns
!    at which diagnostics will be calculated (degrees).
!---------------------------------------------------------------------
real,    dimension(:), allocatable  :: deglon1, deglat1

!---------------------------------------------------------------------
!    iradprt and jradprt are the processor-based i and j coordinates 
!    of the desired diagnostics columns.
!---------------------------------------------------------------------
integer, dimension(:), allocatable  :: jradprt, iradprt

!---------------------------------------------------------------------
!    do_raddg is an array of logicals indicating which latitude rows
!    belonging to the processor contain diagnostics columns.
!---------------------------------------------------------------------
logical, dimension(:), allocatable  :: do_raddg


integer     :: nbly              ! number of frequency bands for exact 
                                 ! cool-to-space calculation
integer     :: n_continuum_bands ! number of bands in the h2o continuum
integer     :: radiag_unit       ! i/o unit to which output file is 
                                 ! written
integer     :: num_pts           !  total number of columns in which
                                 !  diagnostics are desired
logical     :: module_is_initialized = .false.        
                                 ! module initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------





                       contains


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                    PUBLIC SUBROUTINES
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="radiation_diag_init">
!  <OVERVIEW> 
!   Constructor of the radiation_diag_mod module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Constructor of the radiation_diag_mod module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_diag_init (latb, lonb, Lw_tables)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="Lw_tables" TYPE="lw_table_type">
!    lw_tables_type variable containing various longwave
!    table specifiers needed by radiation_diag_mod.
!  </IN>
! </SUBROUTINE>
!
subroutine radiation_diag_init (latb, lonb, Lw_tables)

!---------------------------------------------------------------------
!    radiation_diag_init is the constructor for radiation_diag_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:), intent(in)  ::  latb, lonb         
type(lw_table_type),  intent(in)  ::  Lw_tables
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      2d array of model latitudes at cell corners [radians]
!       lonb      2d array of model longitudes at cell corners [radians]
!       Lw_tables lw_tables_type variable containing various longwave
!                 table specifiers needed by radiation_diag_mod.
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      integer     :: unit, ierr, io, logunit
      integer     :: nn, j, i, nblw
      real        :: dellat, dellon

!--------------------------------------------------------------------
!  local variables
!
!     unit
!
!-------------------------------------------------------------------

!--------------------------------------------------------------------
!    if routine has already been executed, return.
!--------------------------------------------------------------------
      if (module_is_initialized)  return
 
!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call fms_init
      call constants_init
      call rad_utilities_init

!-----------------------------------------------------------------------
!    read namelist.              
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=radiation_diag_nml, iostat=io)
      ierr = check_nml_error(io,'radiation_diag_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiation_diag_nml, iostat=io, end=10) 
        ierr = check_nml_error(io,'radiation_diag_nml')
        end do                   
10      call close_file (unit)   
      endif                      
#endif
                                  
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                         write (logunit, nml=radiation_diag_nml)

!---------------------------------------------------------------------
!    allocate and initialize a flag array which indicates the latitudes
!    containing columns where radiation diagnostics are desired.
!---------------------------------------------------------------------
      allocate (do_raddg (size(latb,2)-1) )
      do_raddg(:) = .false.

!-------------------------------------------------------------------
!    define the total number of points at which diagnostics are desired.
!    points may be specified either by lat-lon pairs or by global index
!    pairs. 
!-------------------------------------------------------------------
      num_pts = num_pts_latlon + num_pts_ij

!-------------------------------------------------------------------
!    continue on only if diagnostics are desired in at least one column.
!-------------------------------------------------------------------
      if (num_pts > 0) then

!-------------------------------------------------------------------
!    if more points are desired than space has been reserved for, print 
!    a message.
!-------------------------------------------------------------------
        if (num_pts > MAX_PTS) then
          call error_mesg ( 'radiation_diag_mod', &
         'must reset MAX_PTS or reduce number of diagnostics points', &
                                                             FATAL)
        endif

!-------------------------------------------------------------------
!    allocate space for arrays which will contain the lat and lon and
!    processor-local i and j indices.
!-------------------------------------------------------------------
        allocate ( deglon1 (num_pts))
        allocate ( deglat1 (num_pts))
        allocate ( jradprt (num_pts))
        allocate ( iradprt (num_pts))

!---------------------------------------------------------------------
!    if any points for diagnostics are specified by (i,j) global 
!    indices, determine their lat-lon coordinates. assumption is made 
!    that the deltas of latitude and longitude are uniform over 
!    the globe.
!---------------------------------------------------------------------
        do nn=1,num_pts_ij
          dellat = latb(1,2) - latb(1,1)
          dellon = lonb(2,1) - lonb(1,1)
          latradprt(nn + num_pts_latlon) =     &
                      (-0.5*acos(-1.0) + (jradprt_gl(nn) - 0.5)*  &
                                           dellat) * radian
          lonradprt(nn + num_pts_latlon) =                & 
                       (iradprt_gl(nn) - 0.5)*dellon*radian
        end do

!--------------------------------------------------------------------
!    determine if the lat/lon values are within the global grid,
!    latitude between -90 and 90 degrees and longitude between 0 and
!    360 degrees.
!--------------------------------------------------------------------
        do nn=1,num_pts
          jradprt(nn) = 0
          iradprt(nn) = 0
          deglat1(nn) = 0.0
          deglon1(nn) = 0.0
          if (latradprt(nn) .ge. -90. .and. &
              latradprt(nn) .le.  90.) then
          else
            call error_mesg ('radiation_diag_mod', &
                ' invalid latitude for radiation diagnostics ', FATAL)
          endif

          if (lonradprt(nn) .ge. 0. .and. &
              lonradprt(nn) .le. 360.) then
          else
            call error_mesg ('radiation_diag_mod', &
                ' invalid longitude for radiation diagnostics ', FATAL)
          endif

!--------------------------------------------------------------------
!    determine if the diagnostics column is within the current 
!    processor's domain. if so, set a logical flag indicating the
!    presence of a diagnostic column on the particular row, define the 
!    i and j processor-coordinates and the latitude and longitude of 
!    the diagnostics column.
!--------------------------------------------------------------------
          do j=1,size(latb,2) - 1
            if (latradprt(nn) .ge. latb(1,j)*radian .and.  &
                latradprt(nn) .lt. latb(1,j+1)*radian) then
              do i=1,size(lonb,1) - 1
                if (lonradprt(nn) .ge. lonb(i,1)*radian   &
                                  .and.&
                    lonradprt(nn) .lt. lonb(i+1,1)*radian)  &
                                   then
                  do_raddg(j) = .true.
                  jradprt(nn) = j
                  iradprt(nn) = i
                  deglon1(nn) = 0.5*(lonb(i,1) + lonb(i+1,1))*  &
                                radian
                  deglat1(nn) = 0.5*(latb(1,j) + latb(1,j+1))*   &
                                radian
                  exit
                endif
              end do
              exit
            endif
          end do
        end do

!----------------------------------------------------------------------
!    open a unit for the radiation diagnostics output.
!---------------------------------------------------------------------
        radiag_unit = open_file ('radiation_diag.out', action='write', &
                                 threading='multi', form='formatted')

!----------------------------------------------------------------------
!    save the input fields from the lw_tables_type variable that will
!    be used by this module.
!----------------------------------------------------------------------
        nbly              = size(Lw_tables%bdlocm(:))
        nblw              = size(Lw_tables%bandlo(:))
        n_continuum_bands = size(Lw_tables%iband(:))

        allocate ( bdlocm (nbly) )
        allocate ( bdhicm (nbly) )
        allocate ( iband  (n_continuum_bands) )
        allocate ( bandlo (nblw) )
        allocate ( bandhi (nblw) )

        bdlocm(:) = Lw_tables%bdlocm(:)
        bdhicm(:) = Lw_tables%bdhicm(:)
        iband (:) = Lw_tables%iband(:)
        bandlo(:) = Lw_tables%bandlo(:)
        bandhi(:) = Lw_tables%bandhi(:)

      endif     ! (num_pts > 0)

!---------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!---------------------------------------------------------------------
      module_is_initialized = .true.



end subroutine radiation_diag_init





!####################################################################
! <SUBROUTINE NAME="radiation_diag_driver">
!  <OVERVIEW>
!   Subroutine to  determine if a diagnostics column is present
!    in the current physics window, and, if so, calls radiag to 
!    obtain the desired variables in that column and output them to
!    a data file.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to  determine if a diagnostics column is present
!    in the current physics window, and, if so, calls radiag to 
!    obtain the desired variables in that column and output them to
!    a data file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_diag_driver (is, ie, js, je, Atmos_input, Surface, Astro, &
!                                  Rad_gases, Cldrad_props,   &
!                                  Cld_spec, Sw_output,   &
!                                  Lw_output, Lw_diagnostics,   &
!                                  Cldspace_rad)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!  </IN>
!  <IN NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomy_type variable containing the astronomical
!     input fields needed by the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_prperties_type">
!   cldrad_prperties_type variable containing the cloud radiative
!   property input fields needed by the radiation package
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cld_specification_type variable containing 
!                   cloud information relevant to the radiation package
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   sw_output_type variable containing shortwave 
!                   radiation output data 
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </IN>
!  <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  </IN>
!  <IN NAME="Cldspace_rad" TYPE="cld_space_properties_type">
!   cld_space_properties_type variable containing infor-
!                   mation on cloud properties seen by the radiation 
!                   package in cloud-space coordinates
!  </IN>
! </SUBROUTINE>
!
subroutine radiation_diag_driver (is, ie, js, je, Atmos_input,  &
                                  Surface, Astro, Rad_gases,   &
                                  Cldrad_props, Cld_spec,  &
                                  Sw_output, Lw_output, Lw_diagnostics,&
                                  Cldspace_rad)

!---------------------------------------------------------------------
!    radiation_diag_driver determines if a diagnostics column is present
!    in the current physics window, and, if so, calls radiag to 
!    obtain the desired variables in that column and output them to
!    a data file.
!----------------------------------------------------------------------

integer,                         intent(in) :: is, ie, js, je
type(atmos_input_type),          intent(in) :: Atmos_input
type(surface_type),              intent(in) :: Surface
type(astronomy_type),            intent(in) :: Astro
type(radiative_gases_type),      intent(in) :: Rad_gases
type(cldrad_properties_type),    intent(in) :: Cldrad_props
type(cld_specification_type),    intent(in) :: Cld_spec       
type(sw_output_type), dimension(:), intent(in) :: Sw_output
type(lw_output_type), dimension(:), intent(in) :: Lw_output
type(lw_diagnostics_type),       intent(in) :: Lw_diagnostics
type(cld_space_properties_type), intent(in) :: Cldspace_rad

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      Surface      surface_type variable containing surface variables
!                   needed by the radiation package
!      Astro        astronomy_type variable containing the astronomical
!                   input fields needed by the radiation package
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!      Cldrad_props cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!      Cld_spec     cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!      Sw_output    sw_output_type variable containing shortwave 
!                   radiation output data 
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!      Cldspace_rad cld_space_properties_type variable containing infor-
!                   mation on cloud properties seen by the radiation 
!                   package in cloud-space coordinates
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      integer    :: j   ! do-loop index

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('radiation_diag_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    if this physics window includes a point at which diagnostics are 
!    to be calculated, call radiag to compute and write out the desired 
!    column data to a data file.
!---------------------------------------------------------------------
      do j=js,je
        if (do_raddg(j)) then
          call radiag (is, ie, js, je, j, Atmos_input, Surface, Astro, &
                       Rad_gases, Cldrad_props, Cld_spec, Sw_output, &
                       Lw_output, Lw_diagnostics, Cldspace_rad)
        endif
      end do

!--------------------------------------------------------------------


end subroutine radiation_diag_driver



!###################################################################
! <SUBROUTINE NAME="radiation_diag_end">
!  <OVERVIEW>
!   radiation_diag_end is the destructor for radiation_diag_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_diag_end is the destructor for radiation_diag_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_diag_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine radiation_diag_end

!--------------------------------------------------------------------
!    radiation_diag_end is the destructor for radiation_diag_mod.
!--------------------------------------------------------------------
      deallocate ( do_raddg )

      if (num_pts > 0 ) then
!--------------------------------------------------------------------
!    close the radiation_diag.out file.
!--------------------------------------------------------------------
        call close_file (radiag_unit)

!--------------------------------------------------------------------
!    deallocate module arrays.
!--------------------------------------------------------------------
        deallocate ( deglon1, deglat1, jradprt, iradprt,  &
                     bdlocm, bdhicm, iband, bandlo, bandhi )   
      endif

!--------------------------------------------------------------------
      module_is_initialized = .false.    

!---------------------------------------------------------------------


end subroutine radiation_diag_end


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                     PRIVATE SUBROUTINES
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#####################################################################
! <SUBROUTINE NAME="radiag">
!  <OVERVIEW>
!   radiag calculates and outputs radiation diagnostics in user-
!    specified columns.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiag calculates and outputs radiation diagnostics in user-
!    specified columns.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiag (is, ie, js, je, jrow, Atmos_input, Surface, Astro, &
!                   Rad_gases, Cldrad_props, Cld_spec, Sw_output, Lw_output,  &
!                   Lw_diagnostics, Cldspace_rad)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="jrow" TYPE="integer">
!   the current physics-window j index, which contains 
!                   a radiation diagnostics column
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!  </IN>
!  <IN NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomy_type variable containing the astronomical
!     input fields needed by the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_prperties_type">
!   cldrad_prperties_type variable containing the cloud radiative
!   property input fields needed by the radiation package
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cld_specification_type variable containing 
!                   cloud information relevant to the radiation package
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   sw_output_type variable containing shortwave 
!                   radiation output data 
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </IN>
!  <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  </IN>
!  <IN NAME="Cldspace_rad" TYPE="cld_space_properties_type">
!   cld_space_properties_type variable containing infor-
!                   mation on cloud properties seen by the radiation 
!                   package in cloud-space coordinates
!  </IN>
! </SUBROUTINE>
!
subroutine radiag (is, ie, js, je, jrow, Atmos_input, Surface, Astro, &
                   Rad_gases, Cldrad_props, Cld_spec, Sw_output,  &
                   Lw_output, Lw_diagnostics, Cldspace_rad)

!--------------------------------------------------------------------
!    radiag calculates and outputs radiation diagnostics in user-
!    specified columns.
!--------------------------------------------------------------------

integer,                         intent(in) :: is, ie, js, je, jrow
type(atmos_input_type),          intent(in) :: Atmos_input
type(surface_type),              intent(in) :: Surface
type(astronomy_type),            intent(in) :: Astro
type(radiative_gases_type),      intent(in) :: Rad_gases
type(cldrad_properties_type),    intent(in) :: Cldrad_props
type(cld_specification_type),    intent(in) :: Cld_spec       
type(sw_output_type), dimension(:), intent(in) :: Sw_output
type(lw_output_type), dimension(:), intent(in) :: Lw_output
type(lw_diagnostics_type),       intent(in) :: Lw_diagnostics
type(cld_space_properties_type), intent(in) :: Cldspace_rad

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      jrow         the current physics-window j index, which contains 
!                   a radiation diagnostics column
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      Surface      surface_type variable containing surface variables
!                   needed by the radiation package
!      Astro        astronomy_type variable containing the astronomical
!                   input fields needed by the radiation package
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!      Cldrad_props cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!      Cld_spec     cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!      Sw_output    sw_output_type variable containing shortwave 
!                   radiation output data 
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!      Cldspace_rad cld_space_properties_type variable containing infor-
!                   mation on cloud properties seen by the radiation 
!                   package in cloud-space coordinates
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      type(lw_output_type)   :: Lw_output_ad
      type(sw_output_type)   :: Sw_output_ad

!---------------------------------------------------------------------
!    these variables are dimensioned by the number of layers in the 
!    grid:
!
!     cts        approximate cool-to-space heating rates for 160-560
!                and 800-990, 1070-1200 cm-1 ranges. (h2o bands)
!                [ degrees K /day ]
!     ctsco2     approximate cool-to-space heating rates for 560-800
!                cm-1 range. (15 micron co2 band)
!                [ degrees K /day ]
!     ctso3      approximate cool-to-space heating rates for 990-1070
!                cm-1 range.  (9.6 micron o3 band)
!                [ degrees K /day ]
!     ctst       approximate cool-to-space heating rates for 160-1200
!                cm-1 range (sum of above 3 variables).
!                [ degrees K /day ]
!     hlwsw      total radiation heating rates.
!                [ degrees K /day ]
!     hlwswcf    total radiation heating rates in the absence of clouds.
!                [ degrees K /day ]
!     convert    factor to convert flux difference (cgs units) to 
!                heating rate in degrees/day.
!                [ (degrees K * sec**3) / ( grams * day) ]
!     cmxolw     amount of maximally overlapped longwave cloud.
!                [ dimensionless ]
!     crndlw     amount of randomly overlapped longwave cloud.
!                [ dimensionless ]
!     camtsw     shortwave cloud amount. 
!                [ dimensionless ]
!     htem       emissivity heating rate, summed over bands 0 - 160,
!                560 - 2200 cm-1.
!                [ degrees K /day ]
!     htem1      emissivity heating rate for 0-160, 1200-2200 cm-1 band,
!                when ch4, n2o are not active,  emissivity heating rate 
!                for 0-160, 1400-2200 cm-1 band when ch4, n2o are 
!                active.
!                [ degrees K /day ]
!     htem2      emissivity heating rate for 560-800 cm-1 band.
!                [ degrees K /day ]
!     htem3      emissivity heating rate for 800-900 cm-1 band
!                [ degrees K /day ]
!     htem4      emissivity heating rate for 900-990 cm-1 band
!                [ degrees K /day ]
!     htem5      emissivity heating rate for 990-1070 cm-1 band.
!                [ degrees K /day ]
!     htem6      emissivity heating rate for 1070-1200 cm-1 band
!                [ degrees K /day ]
!     htem7t     emissivity heating rate over 1200-1400 cm-1 band
!                [ degrees K /day ]
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,3) - 1 ) ::   &
                                           cts, ctsco2, ctso3, ctst, &
                                           hlwsw, hlwswcf, convert, &
                                           cmxolw, crndlw, camtsw, &
                                           htem, htem1, htem2, htem3, &
                                           htem4, htem5, htem6, htem7t


!---------------------------------------------------------------------
!    these variables are dimensioned by the number of flux levels
!    (number of layers + 1):
!
!     flx1       net emissivity flux for 0-160, 1200-2200 cm-1 band 
!                when ch4, n2o not active, net emissivity flux for 
!                0-160, 1400-2200 cm-1 band when ch4, n2o are active.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx2       flux for 560-800 cm-1 band (as one band).
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx3       flux for 800-900 cm-1 band.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx4       flux for 900-990 cm-1 band.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx5       flux for 990-1070 cm-1 band.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx6       flux for 1070-1200 cm-1 band.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx1cf     same as flx1, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx2cf     same as flx2, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx3cf     same as flx3, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx4cf     same as flx4, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx5cf     same as flx5, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx6cf     same as flx6, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     press      pressure at data levels of model.
!                [ pascals ]
!     pflux      pressure at flux levels of model.
!                [ pascals ]
!     flwsw      sum of net lw and sw fluxes.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flxem      sum of emissivity fluxes over bands 0 - 160,
!                560 - 2200 cm-1.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flxemch4n2o  
!                sum of emissivity fluxes over the nbtrge bands in the 
!                1200 - 1400 cm-1 band.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flwswcf    sum of net lw and sw fluxes in the absence of clouds.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!
!----------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,3)   ) ::   &
                                         flx1, flx2, flx3, flx4, flx5, &
                                         flx6, flx1cf, flx2cf, flx3cf, &
                                         flx4cf, flx5cf, flx6cf,    &
                                         press, pflux, flwsw, flxem, &
                                         flxemch4n2o, flwswcf      


!---------------------------------------------------------------------
!    these variables are dimensioned by the number of flux levels  
!    and the number of bands in the 1200 - 1400 cm-1 range when ch4 
!    and n2o are radiatively active gases.
!     flx7       when ch4, n2o are active, flux for 1200-1400 cm-1 
!                band (nbtrge bands).
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!     flx7cf     same as flx7, but for cloud-free conditions.
!                [ ergs / (cm**2  sec), or gram / sec**3 ]
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,3),        &
                       size (Lw_diagnostics%flx1e1f,3) ) :: flx7, flx7cf


!---------------------------------------------------------------------
!    this variable is dimensioned by the number of model layers 
!    and the number of bands in the 1200 - 1400 cm-1 range when ch4 
!    and n2o are radiatively active gases.
!     htem7      emissivity heating rate for each of the nbtrge bands
!                in the 1200-1400 cm-1 band
!                [ degrees K /day ]
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,3)-1,    &
                       size (Lw_diagnostics%flx1e1f,3))  :: htem7 


!---------------------------------------------------------------------
!    this variable is dimensioned by the number of model layers 
!    and the number of frequency bands used in the exact cool-to-space 
!    calculation.
!     exctsn     exact cool-to-space heating rates for each band.
!                [ degrees K /day ]
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,3)-1, nbly) :: exctsn


!---------------------------------------------------------------------
!    these variables are dimensioned by the number of frequency bands 
!    used in the exact cool-to-space calculation.
!     ftopn      total cool to space flux at toa for each of the nbly 
!                bands
!                [ Watts / m**2 , or kg / sec**3 ]
!     ftopac     total cool to space flux at toa summed between band 1
!                and band n
!                [ Watts / m**2 , or kg / sec**3 ]
!     vsumac     cool-to-space flux at the ground summed between bands 
!                1 and n.
!                [ Watts / m**2 , or kg / sec**3 ]
!---------------------------------------------------------------------
      real, dimension (nbly) ::  ftopac, ftopn, vsumac


!---------------------------------------------------------------------
!    various scalar variables:
!---------------------------------------------------------------------

      integer :: ncldsw    ! number of clouds at each grid point.
                           ! [ dimensionless ]
      integer :: nrndlw    ! number of randomly overlapped longwave 
                           ! clouds at each grid point.
                           ! [ dimensionless ]
      integer :: nmxolw    ! number of maximally overlapped longwave 
                           ! clouds at each grid point.
                           ! [ dimensionless ]
      real    :: fdiff     ! difference in fluxes (toa - grd)
                           ! [ Watts / m**2 , or kg / sec**3 ]
      real    ::  qsum     ! sum of cool-to-space fluxes through the 
                           ! column in a given band
                           ! [ Watts / m**2 , or kg / sec**3 ]
      real    :: ftopeft   ! flux at toa summed over nbtrge bands in 
                           ! 1200- 1400 cm-1 range.
                           ! [ Watts / m**2 , or kg / sec**3 ]
      real    :: cldext    ! exp (-cldext) = cloud sw extinction     
                           ! [ dimensionless ]
      real    :: cldssalb  ! cloud single-scattering albedo 
                           ! [ dimensionless ]
      real    :: cvisrfgd_dir  ! visible direct beam sfc albedo
                           ! [ dimensionless ]
      real    :: cvisrfgd_dif  ! visible diffuse beam sfc albedo
                           ! [ dimensionless ]
      real    :: cirrfgd_dir   ! infrared direct beam sfc albedo
                           ! [ dimensionless ]
      real    :: cirrfgd_dif   ! infrared diffuse beam sfc albedo
                           ! [ dimensionless ]
      real    :: dfsw_nir, dfsw_nir_dir, dfsw_nir_dif
!     real    :: ufsw_nir, ufsw_nir_dir, ufsw_nir_dif
      real    :: ufsw_nir,               ufsw_nir_dif
      integer :: ks=1      ! index of top level of radiation grid
      integer :: ke        ! index of lowest radiation grid level
      integer :: nlwcldb   ! number of lw cloud bands
      integer :: nbtrge    ! number of h2o bands in the 1200-1400 
                           ! cm-1 range
      integer :: iloc      ! physics window x coordinate of 
                           ! diagnostics column
      integer :: jloc      ! physics window y coordinate of 
                           ! diagnostics column

!---------------------------------------------------------------------
!  miscellaneous indices
!---------------------------------------------------------------------
      integer ::   kc, nprt, ny, nx, k, nn, m, n
      integer :: nz

      if (Rad_control%do_swaerosol_forcing) then
        Sw_output_ad = Sw_output(Rad_control%indx_swaf)
      endif
      if (Rad_control%do_lwaerosol_forcing) then
        Lw_output_ad = Lw_output(Rad_control%indx_lwaf)
      endif

!---------------------------------------------------------------------
!    define the vertical grid dimension and the number of h2o bands in 
!    the 1200 - 1400 cm-1 range.
!---------------------------------------------------------------------
      ke     = size (Atmos_input%press,3) - 1
      nbtrge = size (Lw_diagnostics%flx1e1f,3)

!--------------------------------------------------------------------
!    check for the diagnostic point(s) in the current jrow. define the 
!    physics window coordinates for those points and process the 
!    radiation diagnostic data.
!---------------------------------------------------------------------
      do nn=1,num_pts
        if (jrow == jradprt(nn)) then
          if ( (iradprt(nn) >= is) .and. (iradprt(nn) <= ie) ) then
            iloc = iradprt(nn) - is + 1
            jloc = jradprt(nn) - js + 1

!---------------------------------------------------------------------
!    write out the latitude and longitude of the model point for which
!    diagnostics will be produced.
!---------------------------------------------------------------------
            write (radiag_unit,99000) deglon1(nn), deglat1(nn)

!----------------------------------------------------------------------
!    write longwave cloud data. determine if any clouds are present
!    in the column. if there are, define the number of lw cloud bands.
!----------------------------------------------------------------------
            write (radiag_unit,9009)
            nmxolw     = Cld_spec%nmxolw(iloc, jloc)
            nrndlw     = Cld_spec%nrndlw(iloc, jloc)
            if (nmxolw > 0 .OR. nrndlw > 0) then
              write (radiag_unit,9010) nmxolw, nrndlw    
              nlwcldb = size (Cldrad_props%emmxolw,4)
!----------------------------------------------------------------------
!    write longwave cloud amounts and emissivities for each cloud band.
!    %emmxolw =  lw cloud emissivity for maximally overlapped clouds.
!                [ dimensionless ]
!    %emrndlw =  lw cloud emissivity for randomly overlapped clouds.
!                [ dimensionless ]
!----------------------------------------------------------------------
              do n=1,nlwcldb
                write (radiag_unit,9041) n
                do k = ks,ke
                  cmxolw(k) = Cld_spec%cmxolw(iloc, jloc, k)
                  crndlw(k) = Cld_spec%crndlw(iloc, jloc, k)
                  if (cmxolw(k) > 0.0 .or. crndlw(k) > 0.0)  then
                    write (radiag_unit,9030)   k,    &
                      cmxolw(k), Cldrad_props%emmxolw(iloc,jloc,k,n,1),&
                      crndlw(k), Cldrad_props%emrndlw(iloc,jloc,k,n,1)
                  endif
                end do
              end do

!--------------------------------------------------------------------
!    if no clouds are present, write a message.
!--------------------------------------------------------------------
            else
              write (radiag_unit, 9052)
            endif 

!----------------------------------------------------------------------
!    define the number of shortwave clouds.
!----------------------------------------------------------------------
            ncldsw = Cld_spec%ncldsw(iloc, jloc)
            write (radiag_unit, 9018)
            write (radiag_unit, 9019) ncldsw

!----------------------------------------------------------------------
!    if clouds exist, write shortwave cloud data.
!----------------------------------------------------------------------
            if (ncldsw /= 0) then

!---------------------------------------------------------------------
!     write out the relevant cloud-radiation variables for the lacis-
!     hansen parameterization from the data in Cldspace_rad:
!     %camtswkc   shortwave cloud amounts. their locations are specified
!                 in the ktopsw/kbtmsw indices. 
!                 [ dimensionless ]
!     %cvisrfswkc reflectivity of clouds in the visible frequency band.
!                 [ dimensionless ]
!     %cirabswkc  absorptivity of clouds in the infrared frequency band.
!                 [ dimensionless ]
!     %cirrfswkc  reflectivity of clouds in the infrared frequency band.
!                 [ dimensionless ]
!     %kbtmswkc   index of flux level pressure of cloud bottom.  
!     %ktopswkc   index of flux level pressure of cloud top. 
!---------------------------------------------------------------------
              if (Sw_control%do_lhsw) then
                write (radiag_unit,9035) 
                write (radiag_unit,9036)   (kc,    &
                         Cldspace_rad%camtswkc  (iloc, jloc,kc), &
                         Cldspace_rad%ktopswkc  (iloc,jloc,kc),   &
                         Cldspace_rad%kbtmswkc  (iloc,jloc,kc)   ,   &
                         Cldspace_rad%cvisrfswkc(iloc,jloc,kc),&
                         Cldspace_rad%cirrfswkc (iloc,jloc,kc),    &
                         Cldspace_rad%cirabswkc (iloc,jloc,kc),&
                                                      kc=ncldsw,1,-1)

!---------------------------------------------------------------------
!     write out the relevant cloud-radiation variables for the expo-
!     nential-sum-fit parameterization from the data in Cldspace_rad:
!---------------------------------------------------------------------
              else if (Sw_control%do_esfsw) then

!---------------------------------------------------------------------
!    for each shortwave cloud band, define an exinction factor,
!    single-scattering albedo and asymmetry factor and write them
!    out along with the cloud amount. use the contents of Cldrad_props,
!    Cld_spec and Atmos_input for these calculations:
!    %cldext      cloud extinction coefficient [ km -1 ]  
!    %cldsct      cloud scattering coefficient [ dimensionless ]
!    %cldsasymm   cloud asymmetry factor  [ dimensionless ]
!    %deltaz      model vertical grid interval [ meters ]
!---------------------------------------------------------------------
                do n=1,size(Cldrad_props%cldext,4)
                  write (radiag_unit,9040) n
                  do k=ks,ke
                    if (Cld_spec%camtsw(iloc,jloc,k) > 0.0) then
                      cldext = Cldrad_props%cldext(iloc,jloc,k,n,1)*   &
                               Atmos_input%clouddeltaz(iloc,jloc,k)* &
                               1.0E-03
                      if (cldext > 0.0) then           
                        cldssalb = Cldrad_props%cldsct(iloc,jloc,k,n,1)/ &
                                   Cldrad_props%cldext(iloc,jloc,k,n,1)
                      else
                        cldssalb = 0.0
                      endif
                      write (radiag_unit,9050) k,     &
                             Cld_spec%camtsw (iloc,jloc,k),   &
                             cldext, cldssalb,                      &
                             Cldrad_props%cldasymm(iloc,jloc,k,n,1)
                    endif
                  end do
                end do

!----------------------------------------------------------------------
!    if no shortwave scheme has been specified, abort.
!----------------------------------------------------------------------
              else
                call error_mesg ('radiation_diag_mod', &
                       'no shortwave clouds are activated', FATAL)
              endif

!----------------------------------------------------------------------
!    if no clouds are present in the column, write out a message. 
!----------------------------------------------------------------------
            else
              write (radiag_unit, 9053)
            endif

!--------------------------------------------------------------------
!    if microphysics has been activated, write microphysical parameters
!    at those levels where cloud is present. use variables found in 
!    Cld_spec.
!    %lwp        liquid water path                   [ kg / m**2 ]
!    %iwp        ice water path                      [ kg / m**2 ]
!    %reff_liq_micro effective cloud drop size       [ microns ]
!    %reff_ice_micro effective ice crystal size      [ microns ]
!----------------------------------------------------------------------
            if (ncldsw /= 0) then
              if (Lw_control%do_lwcldemiss .or.    &
                  Sw_control%do_esfsw) then
                write (radiag_unit, 9510)
                do k=ks,ke     
                  if (cmxolw(k) > 0.0 .or. crndlw(k) > 0.0) then
                    write (radiag_unit, 9520)   k,                &
                     1.0e03*Cld_spec%lwp      (iloc,jloc,k),   &
                     1.0e03*Cld_spec%iwp      (iloc,jloc,k),    &
                     Cld_spec%reff_liq_micro(iloc,jloc,k),  &
                     Cld_spec%reff_ice_micro(iloc,jloc,k)
                  endif
                end do
              endif
            endif

!--------------------------------------------------------------------
!    write out the visible and infrared reflectivities at the ground.
!    currently these are both given the same value (Surface%asfc), 
!    but could be different in the future.
!--------------------------------------------------------------------
!           cvisrfgd = Surface%asfc(iloc,jloc)
!           cirrfgd  = Surface%asfc(iloc,jloc)
            cvisrfgd_dir = Surface%asfc_vis_dir(iloc,jloc)
            cirrfgd_dir  = Surface%asfc_nir_dir(iloc,jloc)
            cvisrfgd_dif = Surface%asfc_vis_dif(iloc,jloc)
            cirrfgd_dif  = Surface%asfc_nir_dif(iloc,jloc)
            write (radiag_unit,9059)
!           write (radiag_unit,9060) cvisrfgd, cirrfgd 
            write (radiag_unit,9060) cvisrfgd_dir, cirrfgd_dir, &
                                     cvisrfgd_dif, cirrfgd_dif

!----------------------------------------------------------------------
!     write out the amounts of the radiative gases that the radiation
!     code sees.     
!--------------------------------------------------------------------
            write (radiag_unit,9069)
            write (radiag_unit,9070) Rad_gases%rrvco2
            write (radiag_unit,9071) Rad_gases%rrvf11
            write (radiag_unit,9072) Rad_gases%rrvf12
            write (radiag_unit,9075) Rad_gases%rrvf113
            write (radiag_unit,9076) Rad_gases%rrvf22
            write (radiag_unit,9073) Rad_gases%rrvch4
            write (radiag_unit,9074) Rad_gases%rrvn2o
 
!---------------------------------------------------------------------
!    define the shortwave parameterization being employed. define the 
!    assumption used to specify the solar zenith angle.
!---------------------------------------------------------------------
            write (radiag_unit, 9079)
            if (Sw_control%do_diurnal) then
              write (radiag_unit,99020) 
            else if (Sw_control%do_annual) then
              write (radiag_unit,99025)
            else if (Sw_control%do_daily_mean) then
              write (radiag_unit,99030)
            else ! (if all 3 are false)
              write (radiag_unit,99040)
            endif

!----------------------------------------------------------------------
!     write out the astronomical data used for this shortwave
!     calculation. use the contents in Astro:
!     %fracday         fraction of averaging period that has daylight.
!                      [ dimensionless ]
!     %cosz            mean cosine of zenith angle for all longitudes.
!                      [ dimensionless ]
!     %solar_constant  solar flux at toa at mean earth-sun radius. 
!                      [ Watts / m**2 , or kg / sec**3 ]
!     %rrsun           earth-sun distance relative to mean distance 
!                      [ dimensionless ]
!----------------------------------------------------------------------
            write (radiag_unit,9080)    &
                               Sw_control%solar_constant*Astro%rrsun, &
                               Astro%cosz(iloc,jloc), &
                               Astro%fracday(iloc,jloc)
            if (Rad_control%hires_coszen) then
              write (radiag_unit, 9081)
              do nz = 1,Rad_control%nzens
                write (radiag_unit, 9082) nz
            write (radiag_unit,9080)    &
                               Sw_control%solar_constant*Astro%rrsun, &
                               Astro%cosz_p(iloc,jloc,nz), &
                               Astro%fracday_p(iloc,jloc,nz)
              end do
             endif
        
!----------------------------------------------------------------------
!    write out atmospheric input data and longwave fluxes and heating
!    rates. Use atmospheric fields from Atmos_input and Rad_gases and 
!    longwave data from Lw_output:
!    %rh2o       mass mixing ratio of h2o at model data levels.
!                [ dimensionless ]
!    %qo3        mass mixing ratio of o3 at model data levels.
!                [ dimensionless ]
!    %temp       temperature at data levels of model.
!                [ degrees K ]
!    %heatra     lw heating rate.
!                [ degrees K / day ]
!    %flxnet     net longwave flux at model flux levels (including the 
!                ground and the top of the atmosphere).
!                [ Watts / m**2 , or kg / sec**3 ]
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    define the pressures at model and flux levels.
!----------------------------------------------------------------------
            do k=ks,ke+1
              pflux(k)  = Atmos_input%pflux(iloc,jloc,k)
              press(k)  = Atmos_input%press(iloc,jloc,k)
            end do

            write(radiag_unit,9090)
            write(radiag_unit,9100)
            write(radiag_unit,9110) (k, press (k),   &
                                     Atmos_input%temp  (iloc,jloc,k),&
                                     Atmos_input%rh2o  (iloc,jloc,k), &
                                     Rad_gases%qo3     (iloc,jloc,k),&
                                     Lw_output(1)%heatra  (iloc,jloc,k),  &
                                     Lw_output(1)%flxnet  (iloc,jloc,k), &
                                     pflux             (k), k=ks,ke)
            write(radiag_unit,9120) press (ke+1), &
                             Atmos_input%temp  (iloc,jloc,ke+1), &
                             Lw_output(1)%flxnet  (iloc,jloc,ke+1), &
                             pflux             (ke+1)

!----------------------------------------------------------------------
!    write out shortwave fluxes and heating rates. use the data in
!    Sw_output:
!    %dfsw       downward short-wave radiation
!                [ Watts / m**2 , or kg / sec**3 ]
!    %fsw        net radiation (up-down) 
!                [ Watts / m**2 , or kg / sec**3 ]
!    %ufsw       upward short-wave radiation.
!                [ Watts / m**2 , or kg / sec**3 ]
!    %hsw        sw radiation heating rates.
!                [ degrees K / day ]
!----------------------------------------------------------------------
            write (radiag_unit,9130)
            if (Sw_control%do_esfsw) then
              write (radiag_unit, 99016)
            else if (Sw_control%do_lhsw) then
              write (radiag_unit, 99018)
            endif
          do nz=1,Rad_control%nzens
                write (radiag_unit, 9082) nz
            write (radiag_unit,9140)
            write (radiag_unit,9150) (k, press(k),   &
                                     Sw_output(1)%hsw  (iloc,jloc,k,nz), &
                                     Sw_output(1)%fsw  (iloc,jloc,k,nz), &
                                     Sw_output(1)%dfsw (iloc,jloc,k,nz),    &
                                     Sw_output(1)%ufsw (iloc,jloc,k,nz),&
                                     pflux          (k), k=ks,ke)
            write (radiag_unit,6556) press(ke+1),    &
                                     Sw_output(1)%fsw  (iloc,jloc,ke+1,nz), &
                                     Sw_output(1)%dfsw (iloc,jloc,ke+1,nz), &
                                     Sw_output(1)%ufsw (iloc,jloc,ke+1,nz), &
                                     pflux          (ke+1)

            if (Sw_control%do_esfsw) then
              dfsw_nir = Sw_output(1)%dfsw(iloc,jloc,ke+1,nz) -   &
                           Sw_output(1)%dfsw_vis_sfc(iloc,jloc,nz)
              dfsw_nir_dir = Sw_output(1)%dfsw_dir_sfc(iloc,jloc,nz) -   &
                             Sw_output(1)%dfsw_vis_sfc_dir(iloc,jloc,nz)
              dfsw_nir_dif = Sw_output(1)%dfsw_dif_sfc(iloc,jloc,nz) -   &
                             Sw_output(1)%dfsw_vis_sfc_dif(iloc,jloc,nz)
              ufsw_nir = Sw_output(1)%ufsw(iloc,jloc,ke+1,nz) -   &
                         Sw_output(1)%ufsw_vis_sfc(iloc,jloc,nz)
              ufsw_nir_dif = Sw_output(1)%ufsw_dif_sfc(iloc,jloc,nz) -   &
                            Sw_output(1)%ufsw_vis_sfc_dif(iloc,jloc,nz)
              write (radiag_unit, 99026)    &
                           Sw_output(1)%dfsw_vis_sfc(iloc,jloc,nz), &
                           Sw_output(1)%ufsw_vis_sfc(iloc,jloc,nz), &
                           Sw_output(1)%dfsw_vis_sfc_dir(iloc,jloc,nz), &
                           Sw_output(1)%dfsw_vis_sfc_dif(iloc,jloc,nz), &
                           Sw_output(1)%ufsw_vis_sfc_dif(iloc,jloc,nz), &
                           dfsw_nir, ufsw_nir,  &
                           dfsw_nir_dir,               &
                           dfsw_nir_dif, ufsw_nir_dif
            endif

!----------------------------------------------------------------------
!    compute and write out total radiative heating and total fluxes 
!    (lw + sw, up-down).
!----------------------------------------------------------------------
            do k=ks,ke+1
              flwsw(k) = Lw_output(1)%flxnet (iloc,jloc,k) +    &
                         Sw_output(1)%fsw    (iloc,jloc,k,nz)
            end do
            do k=ks,ke
              hlwsw(k) = Sw_output(1)%hsw    (iloc,jloc,k,nz) +    &
                         Lw_output(1)%heatra (iloc,jloc,k)
            end do
            write (radiag_unit,9160)
            write (radiag_unit,9170)
            write (radiag_unit,9190) (k, press(k),    &
                                      hlwsw(k), flwsw(k), &
                                      pflux(k), k=ks,ke)
            write (radiag_unit,9180)  press(ke+1), flwsw(ke+1),   &
                                      pflux(ke+1)
      end do

            if (Rad_control%do_totcld_forcing) then
!----------------------------------------------------------------------
!    write out atmospheric input data and longwave fluxes and heating
!    rates for the cloud-free case. use the following data from 
!    Lw_output:
!    %heatracf   lw heating rate in the absence of clouds.
!                [ degrees K / day ]
!    %flxnetcf   net longwave flux at model flux levels (including the 
!                ground and the top of the atmosphere) in the absence of
!                cloud.
!                [ Watts / m**2 , or kg / sec**3 ]
!----------------------------------------------------------------------
              write (radiag_unit,9400)
              write (radiag_unit,9100)
              write (radiag_unit,9110) (k, press (k),  &
                                    Atmos_input%temp  (iloc,jloc,k), &
                                    Atmos_input%rh2o  (iloc,jloc,k), &
                                    Rad_gases%qo3     (iloc,jloc,k), &
                                    Lw_output(1)%heatracf(iloc,jloc,k), &
                                    Lw_output(1)%flxnetcf(iloc,jloc,k), &
                                    pflux (k), k=ks,ke)
              write (radiag_unit,9120)  press (ke+1),  &
                                    Atmos_input%temp  (iloc,jloc,ke+1),&
                                    Lw_output(1)%flxnetcf(iloc,jloc,ke+1),&
                                    pflux(ke+1)

!----------------------------------------------------------------------
!    write out shortwave fluxes and heating rates for the cloud-free 
!    case. use the data in Sw_output:
!    %dfswcf     downward short-wave radiation in absence of clouds
!                [ Watts / m**2 , or kg / sec**3 ]
!    %fswcf      net radiation (up-down) in absence of clouds   
!                [ Watts / m**2 , or kg / sec**3 ]
!    %ufswcf     upward short-wave radiation in absence of clouds
!                [ Watts / m**2 , or kg / sec**3 ]
!    %hswcf      sw radiation heating rates in the absence of clouds.
!                [ degrees K / day ]
!----------------------------------------------------------------------
              write (radiag_unit,9410)
       do nz=1,Rad_control%nzens
                write (radiag_unit, 9082) nz
            write (radiag_unit,9140)
              write (radiag_unit,9150) (k, press(k),   &
                                        Sw_output(1)%hswcf (iloc,jloc,k,nz), &
                                        Sw_output(1)%fswcf (iloc,jloc,k,nz), &
                                        Sw_output(1)%dfswcf(iloc,jloc,k,nz),&
                                        Sw_output(1)%ufswcf(iloc,jloc,k,nz), &
                                        pflux(k), k=ks,ke)
              write (radiag_unit,6556)    &
                                    press(ke+1), &
                                    Sw_output(1)%fswcf(iloc,jloc,ke+1,nz),&
                                    Sw_output(1)%dfswcf(iloc,jloc,ke+1,nz),  &
                                    Sw_output(1)%ufswcf(iloc,jloc,ke+1,nz), &
                                    pflux(ke+1)

!----------------------------------------------------------------------
!    compute and write out total radiative heating and total fluxes 
!    (lw + sw, up-down) for the cloud-free case
!----------------------------------------------------------------------
              do k=ks,ke+1
                flwswcf(k) = Lw_output(1)%flxnetcf(iloc,jloc,k) +    &
                             Sw_output(1)%fswcf   (iloc,jloc,k,nz)
              end do
              do k=ks,ke
                hlwswcf(k) = Sw_output(1)%hswcf   (iloc,jloc,k,nz) +    &
                             Lw_output(1)%heatracf(iloc,jloc,k)
              end do

              write (radiag_unit,9420)
              write (radiag_unit,9170)
              write (radiag_unit,9190) (k, press(k), hlwswcf(k),   &
                                       flwswcf(k), pflux(k), k=ks,ke)
              write (radiag_unit,9180) press(ke+1), flwswcf(ke+1),  &
                                       pflux(ke+1)
      end do
            endif

            if (Rad_control%do_lwaerosol_forcing) then
!----------------------------------------------------------------------
!    write out atmospheric input data and longwave fluxes and heating
!    rates for the aerosol forcing case. use the following data from 
!    Lw_output_ad:
!    %heatra   lw heating rate.
!                [ degrees K / day ]
!    %flxnet   net longwave flux at model flux levels (including the 
!                ground and the top of the atmosphere.
!                [ Watts / m**2 , or kg / sec**3 ]
!    %heatracf   lw heating rate in the absence of clouds.
!                [ degrees K / day ]
!    %flxnetcf   net longwave flux at model flux levels (including the 
!                ground and the top of the atmosphere) in the absence of
!                cloud.
!                [ Watts / m**2 , or kg / sec**3 ]
!----------------------------------------------------------------------
              if (Lw_control%do_lwaerosol) then
! climate includes aerosol effects, lw aerosol forcing by exclusion
                write (radiag_unit,9603)
              else
! climate includes no aerosol effects, lw aerosol forcing by inclusion
                write (radiag_unit,9601)
              endif
              write (radiag_unit,9100)
              write (radiag_unit,9110) (k, press (k),  &
                                    Atmos_input%temp  (iloc,jloc,k), &
                                    Atmos_input%rh2o  (iloc,jloc,k), &
                                    Rad_gases%qo3     (iloc,jloc,k), &
                                    Lw_output_ad%heatra(iloc,jloc,k), &
                                    Lw_output_ad%flxnet(iloc,jloc,k), &
                                    pflux (k), k=ks,ke)
              write (radiag_unit,9120)  press (ke+1),  &
                                    Atmos_input%temp  (iloc,jloc,ke+1),&
                                    Lw_output_ad%flxnet(iloc,jloc,ke+1),&
                                    pflux(ke+1)
! clear-sky results
              if (Lw_control%do_lwaerosol) then
! climate includes aerosol effects, lw aerosol forcing by exclusion
                write (radiag_unit,9604)
              else
! climate includes no aerosol effects, lw aerosol forcing by inclusion
                write (radiag_unit,9602)
              endif
              write (radiag_unit,9100)
              write (radiag_unit,9110) (k, press (k),  &
                                    Atmos_input%temp  (iloc,jloc,k), &
                                    Atmos_input%rh2o  (iloc,jloc,k), &
                                    Rad_gases%qo3     (iloc,jloc,k), &
                                    Lw_output_ad%heatracf(iloc,jloc,k), &
                                    Lw_output_ad%flxnetcf(iloc,jloc,k), &
                                    pflux (k), k=ks,ke)
              write (radiag_unit,9120)  press (ke+1),  &
                                    Atmos_input%temp  (iloc,jloc,ke+1),&
                                    Lw_output_ad%flxnetcf(iloc,jloc,ke+1),&
                                    pflux(ke+1)
            endif

            if (Rad_control%do_swaerosol_forcing) then
!----------------------------------------------------------------------
!    write out atmospheric input data and shortwave fluxes and heating
!    rates for the aerosol forcing case. use the following data from 
!    Sw_output_ad:
!    %dfsw       downward short-wave radiation
!                [ Watts / m**2 , or kg / sec**3 ]
!    %fsw        net radiation (up-down)
!                [ Watts / m**2 , or kg / sec**3 ]
!    %ufsw       upward short-wave radiation
!                [ Watts / m**2 , or kg / sec**3 ]
!    %hsw        sw radiation heating rates
!                [ degrees K / day ]
!    %dfswcf     downward short-wave radiation in absence of clouds
!                [ Watts / m**2 , or kg / sec**3 ]
!    %fswcf      net radiation (up-down) in absence of clouds   
!                [ Watts / m**2 , or kg / sec**3 ]
!    %ufswcf     upward short-wave radiation in absence of clouds
!                [ Watts / m**2 , or kg / sec**3 ]
!    %hswcf      sw radiation heating rates in the absence of clouds.
!                [ degrees K / day ]
!----------------------------------------------------------------------
              if (Sw_control%do_swaerosol) then
! climate includes aerosol effects, sw aerosol forcing by exclusion
                write (radiag_unit,9703)
              else
! climate includes no aerosol effects, sw aerosol forcing by inclusion
                write (radiag_unit,9701)
              endif
              write (radiag_unit,9140)
    do nz = 1, Rad_control%nzens
              write (radiag_unit,9150) (k, press(k),   &
                                        Sw_output_ad%hsw (iloc,jloc,k,nz), &
                                        Sw_output_ad%fsw (iloc,jloc,k,nz), &
                                        Sw_output_ad%dfsw(iloc,jloc,k,nz),&
                                        Sw_output_ad%ufsw(iloc,jloc,k,nz), &
                                        pflux(k), k=ks,ke)
              write (radiag_unit,6556)    &
                                    press(ke+1), &
                                    Sw_output_ad%fsw(iloc,jloc,ke+1,nz),&
                                    Sw_output_ad%dfsw(iloc,jloc,ke+1,nz),  &
                                    Sw_output_ad%ufsw(iloc,jloc,ke+1,nz), &
                                    pflux(ke+1)
! clear-sky results
              if (Sw_control%do_swaerosol) then
! climate includes aerosol effects, sw aerosol forcing by exclusion
                write (radiag_unit,9704)
              else
! climate includes no aerosol effects, sw aerosol forcing by inclusion
                write (radiag_unit,9702)
              endif
              write (radiag_unit,9140)
              write (radiag_unit,9150) (k, press(k),   &
                                        Sw_output_ad%hswcf (iloc,jloc,k,nz), &
                                        Sw_output_ad%fswcf (iloc,jloc,k,nz), &
                                        Sw_output_ad%dfswcf(iloc,jloc,k,nz),&
                                        Sw_output_ad%ufswcf(iloc,jloc,k,nz), &
                                        pflux(k), k=ks,ke)
              write (radiag_unit,6556)    &
                                    press(ke+1), &
                                    Sw_output_ad%fswcf(iloc,jloc,ke+1,nz),&
                                    Sw_output_ad%dfswcf(iloc,jloc,ke+1,nz),  &
                                    Sw_output_ad%ufswcf(iloc,jloc,ke+1,nz), &
                                    pflux(ke+1)
        end do
            endif

            if (Rad_control%do_lwaerosol_forcing .and.   &
                Rad_control%do_swaerosol_forcing) then
       do nz=1,Rad_control%nzens
!----------------------------------------------------------------------
!    compute and write out total radiative heating and total fluxes 
!    (lw + sw, up-down) for the total-sky and cloud-free case
!    with lw and sw aerosol forcing
!----------------------------------------------------------------------
              do k=ks,ke+1
                flwsw(k) = Lw_output_ad%flxnet(iloc,jloc,k) +    &
                             Sw_output_ad%fsw   (iloc,jloc,k,nz)
                flwswcf(k) = Lw_output_ad%flxnetcf(iloc,jloc,k) +    &
                             Sw_output_ad%fswcf   (iloc,jloc,k,nz)
              end do
              do k=ks,ke
                hlwsw(k) = Sw_output_ad%hsw   (iloc,jloc,k,nz) +    &
                             Lw_output_ad%heatra(iloc,jloc,k)
                hlwswcf(k) = Sw_output_ad%hswcf   (iloc,jloc,k,nz) +    &
                             Lw_output_ad%heatracf(iloc,jloc,k)
              end do

              write (radiag_unit,9801)
              write (radiag_unit,9170)
              write (radiag_unit,9190) (k, press(k),    &
                                        hlwsw(k), flwsw(k), &
                                        pflux(k), k=ks,ke)
              write (radiag_unit,9180)  press(ke+1), flwsw(ke+1),   &
                                        pflux(ke+1)

              write (radiag_unit,9802)
              write (radiag_unit,9170)
              write (radiag_unit,9190) (k, press(k), hlwswcf(k),   &
                                       flwswcf(k), pflux(k), k=ks,ke)
              write (radiag_unit,9180) press(ke+1), flwswcf(ke+1),  &
                                       pflux(ke+1)
     end do
            endif

!----------------------------------------------------------------------
!    define emissivity fluxes, both the standard and cloud-free case.
!    note that Lw_diagnostics%fluxn is in cgs units (ergs/(cm**2 sec).
!----------------------------------------------------------------------
            do k=ks,ke+1
              flx1(k) = Lw_diagnostics%fluxn(iloc,jloc,k,1)
              flx2(k) = Lw_diagnostics%fluxn(iloc,jloc,k,2)
              flx3(k) = Lw_diagnostics%fluxn(iloc,jloc,k,3)
              flx4(k) = Lw_diagnostics%fluxn(iloc,jloc,k,4)
              flx5(k) = Lw_diagnostics%fluxn(iloc,jloc,k,5)
              flx6(k) = Lw_diagnostics%fluxn(iloc,jloc,k,6)
              if (Rad_control%do_totcld_forcing) then
                flx1cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,1)
                flx2cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,2)
                flx3cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,3)
                flx4cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,4)
                flx5cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,5)
                flx6cf(k) = Lw_diagnostics%fluxncf(iloc, jloc, k,6)
              endif
              if (nbtrge > 0) then
                do m=1,nbtrge
                  if (Rad_control%do_totcld_forcing) then
                    flx7cf(k,m) =     &
                             Lw_diagnostics%fluxncf(iloc, jloc, k,6+m)
                  endif
                  flx7(k,m) = Lw_diagnostics%fluxn(iloc, jloc, k,6+m)
                end do
              endif
            end do

!--------------------------------------------------------------------
!    define the factor used to convert a flux divergence in cgs units 
!    to a heating rate in units of degrees per day. the 1.0e-03 
!    converts the fluxes (in cgs) to mks units [ (ergs/cm^2/s)  X 
!    1.0e-03  ---> (J/m^2/s) ]. radcon includes the conversion to 
!    degrees/day from degrees/second.
!---------------------------------------------------------------------
            do k=ks,ke
              convert(k) = 1.0e-03*radcon_mks*    &
                           (1.0/(pflux(k+1) - pflux(k)))
            end do

!----------------------------------------------------------------------
!    compute emissivity heating rates in degrees K / day.
!----------------------------------------------------------------------
            do k=ks,ke
              htem1(k) = (flx1(k+1) - flx1(k))*convert(k)
              htem2(k) = (flx2(k+1) - flx2(k))*convert(k)
              htem3(k) = (flx3(k+1) - flx3(k))*convert(k)
              htem4(k) = (flx4(k+1) - flx4(k))*convert(k)
              htem5(k) = (flx5(k+1) - flx5(k))*convert(k)
              htem6(k) = (flx6(k+1) - flx6(k))*convert(k)
            end do
            if (nbtrge > 0) then
              do m=1,nbtrge
                do k=ks,ke
                  htem7(k,m) = (flx7(k+1,m) - flx7(k,m))* convert(k)
                end do
              end do

!--------------------------------------------------------------------
!    define the sum of the heating rates in the 1200 - 1400 cm-1 band.
!--------------------------------------------------------------------
              do k=ks,ke
                htem7t(k) = 0.0E+00
                do m=1,nbtrge
                  htem7t(k) = htem7t(k) + htem7(k,m)
                end do
              end do
            endif

!----------------------------------------------------------------------
!    define the emissivity heating rate summed over all frequencies.
!----------------------------------------------------------------------
            do k=ks,ke
              htem(k) = htem1(k) + htem2(k) + htem3(k) + htem4(k) +   &
                        htem5(k) + htem6(k)
            end do
            if (nbtrge > 0) then
              do k=ks,ke
                htem(k) = htem(k) + htem7t(k)
              enddo
            endif

!----------------------------------------------------------------------
!    write approximate emissivity heating rates.
!----------------------------------------------------------------------
            if (nbtrge == 0) then
              write (radiag_unit,9200)
              write (radiag_unit,9210) (k, press(k), htem1(k),   &
                                        htem2(k), htem3(k), htem4(k), &
                                        htem5(k), htem6(k), htem(k), &
                                        k=ks,ke)
            else
              if (nbtrge .EQ. 1) then
                write (radiag_unit,9201)
                write (radiag_unit,9211) (k, press(k), htem1(k),  &
                                         htem2(k), htem3(k), htem4(k), &
                                         htem5(k), htem6(k),  &
                                         htem7(k,1), htem(k),k=ks,ke)
              else if (nbtrge .EQ. 2) then
                write (radiag_unit,9202)
                write (radiag_unit,9212) (k, press(k), htem1(k),  &
                                         htem2(k), htem3(k), htem4(k), &
                                         htem5(k), htem6(k),   &
                                         (htem7(k,n),n=1,nbtrge),  &
                                         htem7t(k), htem(k),k=ks,ke)
              else if (nbtrge .EQ. 4) then
                write (radiag_unit,9203)
                write (radiag_unit,9213) (k, press(k), htem1(k),   &
                                         htem2(k), htem3(k), htem4(k),&
                                         htem5(k), htem6(k),   &
                                         (htem7(k,n),n=1,nbtrge),   &
                                         htem7t(k), htem(k),k=ks,ke)
              else if (nbtrge .EQ. 10) then
                write (radiag_unit,9201)
                write (radiag_unit,9211) (k, press(k), htem1(k),  &
                                         htem2(k), htem3(k), htem4(k), &
                                         htem5(k), htem6(k), &
                                         htem7t(k), &
                                         htem(k),k=ks,ke)
              else if (nbtrge .EQ. 20) then
                write (radiag_unit,9201)
                write (radiag_unit,9211) (k, press(k), htem1(k),  &
                                         htem2(k), htem3(k), htem4(k), &
                                         htem5(k), htem6(k), &
                                         htem7t(k), &
                                         htem(k),k=ks,ke)
              endif
            endif

!----------------------------------------------------------------------
!    compute and write out approximate cool-to-space heating rates 
!    for the h2o, 15 micron co2 and 9.6 micron o3 bands individually, 
!    and for their sum. 
!----------------------------------------------------------------------
            do k=ks,ke
              ctsco2(k) = Lw_diagnostics%cts_out(iloc,jloc,k,2)
              ctso3(k)  = Lw_diagnostics%cts_out(iloc,jloc,k,5)
              cts (k)   = Lw_diagnostics%cts_out(iloc,jloc,k,1)
              cts (k)   = cts(k) +      &
                          Lw_diagnostics%cts_out(iloc,jloc,k,3)
              cts (k)   = cts(k) +     &
                          Lw_diagnostics%cts_out(iloc,jloc,k,4)
              cts (k)   = cts(k) +    &
                          Lw_diagnostics%cts_out(iloc,jloc,k,6)
            end do
            do k=ks,ke
              ctst(k) = ctso3(k) + ctsco2(k) + cts(k)
            end do

            write (radiag_unit,9220)
            write (radiag_unit,9230) (k, press(k), cts(k), ctsco2(k),  &
                                     ctso3(k), ctst(k), k=ks,ke)

!----------------------------------------------------------------------
!    write out exact cool-to-space heating rates, total and for each
!    individual band.
!    Lw_diagnostics%excts      exact cool-to-space heating rates for 
!                              160-1200 cm-1 range, when using ckd2.1
!                              continuum, or 560 - 1200 cm-1 range
!                              when using Roberts continuum.
!                              [ degrees K / day ]
!----------------------------------------------------------------------
            do n=1,nbly
              do k=ks,ke
                exctsn(k,n) =  Lw_diagnostics%exctsn(iloc,jloc,k,n)
              end do
            end do
            write (radiag_unit,9240)
            write (radiag_unit,9250) (k, press(k),            &
                                     Lw_diagnostics%excts(iloc,jloc,k),&
                                     (exctsn(k,n), n=1,7) , k=ks,ke)
            write (radiag_unit,9260)
            write (radiag_unit,9250) (k, press(k),   &
                                     (exctsn(k,n), n=8,15) , k=ks,ke)
            if (nbly == 48) then
              write (radiag_unit,9261)
              write (radiag_unit,9250) (k, press(k),   &
                                       (exctsn(k,n), n=16,23) , k=ks,ke)
              write (radiag_unit,9262)
              write (radiag_unit,9250) (k, press(k),   &
                                       (exctsn(k,n), n=24,31) , k=ks,ke)
              write (radiag_unit,9263)
              write (radiag_unit,9250) (k, press(k),   &
                                       (exctsn(k,n), n=32,39) , k=ks,ke)
              write (radiag_unit,9264)
              write (radiag_unit,9250) (k, press(k),   &
                                      (exctsn(k,n), n=40,47) , k=ks,ke)
            endif

!----------------------------------------------------------------------
!    compute net flux at each level summed over all bands.
!----------------------------------------------------------------------
            do k=ks,ke+1
              flxem(k) = flx1(k) + flx2(k) + flx3(k) + flx4(k) +  &
                         flx5(k) + flx6(k)
            end do
            if (nbtrge > 0) then
              flxemch4n2o(:) = 0.0E+00
              do m=1,nbtrge
                do k=ks,ke+1 
                  flxemch4n2o(k) = flxemch4n2o(k) + flx7(k,m)
                end do
              end do
              do k=ks,ke+1
                flxem(k) = flxem(k) + flxemch4n2o(k)
              end do
            endif

!----------------------------------------------------------------------
!    compute sum of flux through atmosphere for each band by converting
!    back the heating rates. the flux at toa is the difference 
!    between the surface flux and this net flux through the atmosphere.
!    %fctsg     cool-to-space flux at the ground for each band.
!                [ Watts / m**2 , or kg / sec**3 ]
!--------------------------------------------------------------------
            do n=1,nbly-1
              qsum = 0.0E+00
              do k=ks,ke
                qsum = qsum + 1.0e-03*exctsn(k,n)/convert(k)
              end do
              ftopn(n) = Lw_diagnostics%fctsg(iloc,jloc,n) - qsum
            end do
            ftopn(nbly) = 0.0E+00

!----------------------------------------------------------------------
!    compute the accumulated sum over bands from band 1 to band n for
!    the surface and toa fluxes.
!---------------------------------------------------------------------- 
            ftopac(1) = ftopn(1)
            vsumac(1) = Lw_diagnostics%fctsg(iloc,jloc,1)
            do n=2,nbly
              ftopac(n) = ftopac(n-1) + ftopn(n)
              vsumac(n) = vsumac(n-1) +           &
                                     Lw_diagnostics%fctsg(iloc,jloc,n)
            end do

!----------------------------------------------------------------------
!    write toa and surface fluxes and the differences between them.
!    %gxcts      flux at top of atmosphere for 160-1200 cm-1 range. 
!                [ Watts / m**2 , or kg / sec**3 ]
!    %flx1e1     flux at top of atmosphere for 0-160, 1200-2200
!                cm-1 range.
!                [ Watts / m**2 , or kg / sec**3 ]
!    %flx1e1f    flux at top of atmosphere for nbtrge bands in 1200-
!                1400 cm-1 range.
!                [ Watts / m**2 , or kg / sec**3 ]
!----------------------------------------------------------------------
            fdiff = Lw_diagnostics%gxcts(iloc,jloc) +   &
                    Lw_diagnostics%flx1e1(iloc,jloc) -    &
                    Lw_output(1)%flxnet(iloc,jloc,ke+1)
            write (radiag_unit,9270)      &
                                Lw_diagnostics%gxcts(iloc,jloc), &
                                Lw_diagnostics%flx1e1(iloc,jloc), &
                                Lw_diagnostics%gxcts(iloc,jloc)+ &
                                    Lw_diagnostics%flx1e1(iloc,jloc),&
                                Lw_output(1)%flxnet(iloc,jloc,ke+1), &
                                fdiff
            if (nbtrge > 0) then
              do m=1,nbtrge
                write (radiag_unit,9271) m,    &
                                   Lw_diagnostics%flx1e1f(iloc,jloc,m)
              end do
              ftopeft   = 0.0E+00
              do m=1,nbtrge
                ftopeft   = ftopeft +    &
                                    Lw_diagnostics%flx1e1f(iloc,jloc,m)
              end do
              write (radiag_unit,9272) ftopeft
            endif

!----------------------------------------------------------------------
!    write out toa and sfc fluxes for 8 combined continuum bands between
!    160-560 cm-1 when ckd2.1 is not active, toa and sfc fluxes for 40 
!    combined continuum bands when ckd2.1 is active.
!----------------------------------------------------------------------
            write(radiag_unit,9280)
            do ny=1,nbly-8   
              nprt = 1
              do nx=1,n_continuum_bands
                if (iband(nx) .EQ. ny) then
                  if (nprt .EQ. 1) then
                    write (radiag_unit,9290) ny,   &
                          bandlo(nx+16), bandhi(nx+16), &
                          ftopn(ny), ftopac(ny),    &
                          Lw_diagnostics%fctsg(iloc,jloc,ny), vsumac(ny)
                    nprt = 0
                  else
                    write (radiag_unit,9300) bandlo(nx+16),   &
                                             bandhi(nx+16)
                  endif
                endif
              end do
            end do

!----------------------------------------------------------------------
!    write out toa and sfc fluxes for remaining bands.
!----------------------------------------------------------------------
            do ny =nbly-7, nbly
              write (radiag_unit,9290) ny,     &
                        bdlocm(ny), bdhicm(ny), ftopn(ny),ftopac(ny), &
                        Lw_diagnostics%fctsg(iloc,jloc,ny), vsumac(ny)
            end do

!----------------------------------------------------------------------
!    write out emissivity fluxes.
!----------------------------------------------------------------------
            write (radiag_unit,9310)
            if (nbtrge == 0) then
              write (radiag_unit,9320) 
              write (radiag_unit,9330) (k, flx1(k),flx2(k), flx3(k),  &
                                       flx4(k), flx5(k), flx6(k), &
                                       flxem(k),  k=ks,ke+1)
            else
              if (nbtrge .EQ. 1) then
                write (radiag_unit,9321)
                write (radiag_unit,9331) (k, flx1(k),flx2(k), &
                                         flx3(k),flx4(k), &
                                         flx5(k), flx6(k), &
                                         flxemch4n2o(k), flxem(k), &
                                         k=ks,ke+1)
              else if (nbtrge .EQ. 2) then
                write (radiag_unit,9322)
                write (radiag_unit,9332) (k, flx1(k),flx2(k), &
                                         flx3(k),flx4(k), &
                                         flx5(k), flx6(k), &
                                         (flx7(k,m   ),m=1,nbtrge), &
                                         flxemch4n2o(k), flxem(k), &
                                         k=ks,ke+1)
              else if (nbtrge .EQ. 4) then
                write (radiag_unit,9323)
                write (radiag_unit,9333) (k, flx1(k),flx2(k), &
                                         flx3(k),flx4(k), &
                                         flx5(k), flx6(k),&
                                         (flx7(k,m),m=1,nbtrge),&
                                         flxemch4n2o(k), flxem(k), &
                                         k=ks,ke+1)
              endif
            endif
          endif
        endif
      end do    ! (num_pts loop)


!----------------------------------------------------------------------
!     format statements.
!----------------------------------------------------------------------

99000  format (/////' GRID POINT LOCATION (DEGREES) : LON = ', &
               F10.5, 2X, ' LAT = ', F10.5)
99016  format (/, ' THIS RUN USED THE EXPONENTIAL-SUM-FIT SW &
               &PARAMETERIZATION.',//)
99018  format (/, ' THIS RUN USED THE LACIS-HANSEN SW &
               &PARAMETERIZATION.',//)
99020  format (' SHORTWAVE CALCULATIONS BASED ON &
               &DIURNALLY VARYING ZENITH ANGLES')
99025  format (' SHORTWAVE CALCULATIONS BASED ON  &
               &ANNUAL MEAN ZENITH ANGLES')
99026  format ( '      VISIBLE SFC SW FLUXES:', //, &
                ' total downward   = ', F12.6,   &
                '   total upward   = ', F12.6,   /, &
                ' downward direct  = ', F12.6,  &
                '   upward direct  =   NONEXISTENT', /,    &
                ' downward diffuse = ', F12.6,  &
                '   upward diffuse = ', F12.6,  //,   &
                '       NIR SFC SW FLUXES:', //, &
                ' total downward   = ', F12.6, &
                '   total upward   = ', F12.6,  /,  &
                ' downward direct  = ', F12.6, &
                '   upward direct  =   NONEXISTENT', /,    &
                ' downward diffuse = ', F12.6, &
                '   upward diffuse = ', F12.6)
99030  format (' SHORTWAVE CALCULATIONS BASED ON &
               &DIURNALLY AVERAGED ZENITH ANGLES')
99040  format (' SHORTWAVE CALCULATIONS BASED ON &
               &SPECIFIED ASTRONOMICAL INPUTS')
9009   format (///, ' ************ LONGWAVE CLOUD DATA ***************')
9010   format (/,' NO. MAX OVERLAP CLOUDS= ',I2,    &
                  ' NO. RANDOM OVERLAP CLOUDS= ',I2)
9018   format (///, ' ************ SHORTWAVE CLOUD DATA **************')
9019   format (/,' NO. SW CLOUDS = ',I2)        
9030   format (I4,7X,4F14.6)
9035   format (22X,' SW CLOUD DATA '/,&
               ' CLD. NO',8X,'CLD. AMT.',2X, &
               'CLD TOP INDEX',2X,'CLD BOT INDEX',2X,'VIS. REFL',3X, &
               ' IR REFL',4X,' IR ABS.')  
9036   format (I5,7X,F12.6,I8,I15,6X,3F12.6)
9040   format (12X,      &
               ' SW CLOUD DATA, BAND = ',i2, &
               /,' MDL. LVL',7X, 'CLD. AMT.',7X,'EXT OP DEP.',3X, &
               ' SSALB.',2X,' ASYMM. PAR.')
9041   format (27X,' LW CLOUD DATA, BAND = ',i2,/,' MDL. LVL',4X, &
               'MXO CLD AMT.',2X,'MXO CLD EMIS',2X, &
               'RNDO CLD AMT.',2X,'RNDO CLD EMIS')
9050   format (I5,7X,F12.6,6X,3F12.6)
9052   format (/, ' THIS LW RADIATION CALL SEES NO CLOUDS.')
9053   format (/, ' THIS SW RADIATION CALL SEES NO CLOUDS.')
9059   format (//, ' *********** SURFACE ALBEDO DATA ****************')
!9060   format (/,10X,'VIS. SFC. ALBEDO=',F12.6,' IR SFC. ALBEDO=',  &
!              F12.6)
9060   format (/,'ALBEDO, VIS. SFC. DIRECT =',F10.6,' ALBEDO, NIR SFC. DIRECT =',  &
               F10.6,   &
              /,'ALBEDO, VIS. SFC. DIFFUSE=',F10.6,' ALBEDO, NIR SFC. DIFFUSE=',  &
               F10.6)
9069   format (//, ' *********** RADIATIVE GAS DATA ****************')
9070   format (/,' CO2 VOL.  MIXING RATIO = ', 6PF10.2,' ppmv')
9071   format (' F11 VOL.  MIXING RATIO = ',12PF10.2,' pptv')
9072   format (' F12 VOL.  MIXING RATIO = ',12PF10.2,' pptv')
9075   format (' F113 VOL. MIXING RATIO = ',12PF10.2,' pptv')
9076   format (' F22 VOL.  MIXING RATIO = ',12PF10.2,' pptv')
9073   format (' CH4 VOL.  MIXING RATIO = ', 9PF10.2,' ppbv')
9074   format (' N2O VOL.  MIXING RATIO = ', 9PF10.2,' ppbv')
9079   format (//, ' *********** ASTRONOMICAL DATA ****************')
9080   format (/,' INCOMING SOLAR FLUX =',F12.6,' W/M**2',/,  &
               ' COS(AZIMUTH)=',F12.6,10X,' FRACTION SUNUP=',F12.6)
9081   format (//, 'FOR THE HIRES ZENITH ANGLES OF THIS STEP:')
9082   format (/, ' ZENITH ANGLE NUMBER :', I4)
9090   format (//,'********* LW HEATING RATES AND FLUXES ***********',/)
9100   format ('  LVL',' PRESSURE   ',4X,' TEMP.     ','H2O MMR',5X,&
               'O3 MMR',7X,'HEAT RATE',2X,'NET FLUX',3X,'FLUX PRESS.')
9110   format (I4,E13.6,F12.4,2E12.5,2F12.6,E13.6)
9120   format (4X,E13.6,F12.4,36X,F12.6,E13.6)
9130   format (/,'*************** SW HEATING RATES AND FLUXES ******',/)
9140   format ('  LVL',' PRESSURE    ',3X,'HEAT RATE',2X,'NET FLUX',  &
               4X,'DN FLUX',6X,'UP FLUX',3X,'FLUX PRESS.') 
6556   format (4X,E13.6,12X,3F12.6,E13.6)
9150   format (I4,E13.6,4F12.6,E13.6)
9160   format (/,'*********** COMBINED HEATING RATES AND FLUXES ****',/)
9170   format ('  LVL',' PRESSURE    ',4X,'HEAT RATE',2X,'NET FLUX',  &
               3X,'FLUX PRESS.') 
9180   format (4X,E13.6,12X,F12.6,E13.6)
9190   format (I4,E13.6,2F12.6,E13.6)
9200   format (/,'****   APPROXIMATE HEATING RATES  (Q(APPROX)) ****'/ &
               '  LVL',' PRESSURE   ',5X,'  0-160,1200-2200 ', &
               '    560-800       ','     800-900      ',   &
               '      900-990     ', '    990-1070      ',  &
               '     1070-1200    ', '       TOTAL')
9210   format (I4,E13.6,7F18.6)
9201   format (/,'****   APPROXIMATE HEATING RATES  (Q(APPROX)) ****'/ &
               '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',  &
               '    560-800       ','     800-900      ',  &
               '      900-990     ', '    990-1070      ',   &
               '     1070-1200    ','     1200-1400    ', &
               '       TOTAL')
9211   format (I4,E13.6,8F18.6)
9202   format (/,'****   APPROXIMATE HEATING RATES  (Q(APPROX)) ****'/ &
               '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',   &
               '    560-800       ','     800-900      ', &
               '      900-990     ',  '    990-1070      ',  &
               '     1070-1200    ','  1200-1300 ', &
               '  1300-1400 ','  1200-1400 ','       TOTAL')
9212   format (I4,E13.6,6F18.6,3F12.6,F18.6)
9203   format (/,'****   APPROXIMATE HEATING RATES  (Q(APPROX)) ****'/ &
               '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',  &
               '    560-800       ','     800-900      ',   &
               '      900-990     ', '    990-1070      ',  &
               '     1070-1200    ','  1200-1250 ',  &
               '  1250-1300 ','  1300-1350 ','  1350-1400 ', &
               '  1200-1400 ','       TOTAL')
9213   format (I4,E13.6,6F18.6,5F12.6,F18.6)
9220   format (/,'*******APPROXIMATE CTS HEATING RATES *****'/    &
               '  LVL',' PRESSURE',&
               7X,' H2O BANDS    ',' 15 UM BAND   ',  &
               ' 9.6 UM BAND  ',' TOTAL')
9230   format (I4,E13.6,4F14.6)
9240   format (/,'********EXACT CTS HEATING RATES, BY BAND *******'/   &
               '  LVL',' PRESSURE   ','    TOTAL    ',5X,'1',11X,  &
               '2',11X,'3',  11X,'4',11X,'5',11X,'6',11X,'7',/)
9250   format (I4,E13.6,8F12.6)
9260   format ('  LVL PRESSURE   ',7X,'8',11X,'9',10X,'10',10X,'11', &
               10X,'12',10X,'13',10X,'14',10X,'15')
9261   format ('  LVL PRESSURE   ',6X,'16',10X,'17',10X,'18',10X,'19',&
               10X,'20',10X,'21',10X,'22',10X,'23')
9262   format ('  LVL PRESSURE   ',6X,'24',10X,'25',10X,'26',10X,'27',&
               10X,'28',10X,'29',10X,'30',10X,'31')
9263   format ('  LVL PRESSURE   ',6X,'32',10X,'33',10X,'34',10X,'35',&
               10X,'36',10X,'37',10X,'38',10X,'39')
9264   format ('  LVL PRESSURE   ',6X,'40',10X,'41',10X,'42',10X,'43',&
               10X,'44',10X,'45',10X,'46',10X,'47')
9270   format ( 40X,'   FLUXES'/   &
               ' FLUX AT TOP,160-1200 CM-1       =',F14.6,' W/M**2'/ &
               ' FLUX AT TOP,0-160,1200-2200 CM-1=',F14.6,' W/M**2'/&
               ' FLUX AT TOP,0-2200 CM-1         =',F14.6,' W/M**2'/&
               ' NET FLUX AT GROUND,0-2200 CM-1  =',F14.6,' W/M**2'/  &
               ' NET FLUX DIFFERENCE,0-2200 CM-1 =',F14.6,' W/M**2')
9271   format ( /,  &
               ' FLUX AT TOP, BAND ',I2,' 1200-1400 CM-1 RANGE =',  &
                F14.6, ' W/M**2')
9272   format ( /,   &
               ' FLUX AT TOP, 1200-1400 CM-1 BAND  =',F14.6, &
                 ' W/M**2')
9280   format (/,'**********CTS FLUXES **********'/   &
               1X,'BAND NO',8X,'LOFREQ',9X,'HIFREQ',9X,'F(1)',  &
               11X,'ACCUM. F(1)',4X,'CTS F(GRD)',5X,'ACCUM. CTS F(GRD)')
9290   format (I11,6F15.6)
9300   format (11X,2F15.6)
9310   format (/,'********* EMISSIVITY FLUXES ***********')
9320   format (/,2x,' lvl ',2x,   &
               ' h2o emiss ',' 560-800   ',' 800-900   ',  &
               ' 900-990   ',' 990-1070  ',' 1070-1200 ','  total    ') 
9321   format (/,2x,' lvl ',2x,   &
               ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ',&
               ' 990-1070  ',' 1070-1200 ',' 1200-1400 ',  &
               '  total    ')
9322   format (/,2x,' lvl ',2x,  &
               ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ',&
               ' 990-1070  ',' 1070-1200 ',' 1200-1300 ', &
               ' 1300-1400 ',' 1200-1400 ',' total     ')
9323   format (/,2x,' lvl ',2x,   &
               ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ',&
               ' 990-1070  ',' 1070-1200 ',' 1200-1250 ',' 1250-1300 ',&
               ' 1300-1350 ',' 1350-1400 ',' 1200-1400 ',  &
               '  total    ')
9330   format (i5,-3p,7f11.5)
9331   format (i5,-3p,8f11.5)
9332   format (i5,-3p,10f11.5)
9333   format (i5,-3p,12f11.5)
9400   format (/,'***** CLEAR-SKY LW HEATING RATES AND FLUXES ******',/)
9410   format (/,'**** CLEAR-SKY SW HEATING RATES AND FLUXES******',/)
9420   format (/,'*** COMBINED CLEAR-SKY HEATING RATES AND FLUXES **',/)
9510   format (///, '********* CLOUD MICROPHYSICAL PARAMETERS ******', &
               /, 2x,'lyr',8x, 'liq water path', 3x,    &
               'ice water path',3x, 'eff diam water', ' eff diam ice')
9520   format (I5,7X,F14.6,3x,F14.6,3x,F14.6,3x,F14.6)
9601   format (/,'***** TOTAL-SKY LW HEATING RATES AND FLUXES (AEROSOLS INCLUDED) ******',/)
9602   format (/,'***** CLEAR-SKY LW HEATING RATES AND FLUXES (AEROSOLS INCLUDED) ******',/)
9603   format (/,'***** TOTAL-SKY LW HEATING RATES AND FLUXES (AEROSOLS EXCLUDED) ******',/)
9604   format (/,'***** CLEAR-SKY LW HEATING RATES AND FLUXES (AEROSOLS EXCLUDED) ******',/)
9701   format (/,'***** TOTAL-SKY SW HEATING RATES AND FLUXES (AEROSOLS INCLUDED) ******',/)
9702   format (/,'***** CLEAR-SKY SW HEATING RATES AND FLUXES (AEROSOLS INCLUDED) ******',/)
9703   format (/,'***** TOTAL-SKY SW HEATING RATES AND FLUXES (AEROSOLS EXCLUDED) ******',/)
9704   format (/,'***** CLEAR-SKY SW HEATING RATES AND FLUXES (AEROSOLS EXCLUDED) ******',/)
9801   format (/,'**** COMBINED HEATING RATES AND FLUXES -- AEROSOL FORCING ****',/)
9802   format (/,'*** COMBINED CLEAR-SKY HEATING RATES AND FLUXES -- AEROSOL FORCING **',/)

!--------------------------------------------------------------------


end subroutine radiag         



!##################################################################



                    end module radiation_diag_mod

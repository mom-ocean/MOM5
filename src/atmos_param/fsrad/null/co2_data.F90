
      Module CO2_Data_Mod

!-----------------------------------------------------------------------

      use fs_profile_mod, ONLY:  fs_profile
      Use     co2int_mod, ONLY:  co2int, TRNS
      Use        fms_mod, ONLY:  open_namelist_file, mpp_pe,  &
                                 Error_Mesg, FATAL, close_file,  &
                                 write_version_number, mpp_root_pe

implicit none
private

!      Private 
!-----------------------------------------------------------------------
!
!   Pretabulated co2 transmission functions, evaluated using the
!   methods of Fels and Schwarzkopf (1981) and Schwarzkopf and
!   Fels (1985). 
!
!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 560-800 cm-1 band. also included are the
!   standard temperatures and the weighting function.
!   This data was formerly in COMMON /CO2BD3/.
!
!       CO251    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO258    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE) 
!                     WITH P(SFC)= ^810 MB
!       CDT51    =  FIRST TEMPERATURE DERIVATIVE OF CO251 
!       CDT58    =  FIRST TEMPERATURE DERIVATIVE OF CO258 
!       C2D51    =  SECOND TEMPERATURE DERIVATIVE OF CO251
!       C2D58    =  SECOND TEMPERATURE DERIVATIVE OF CO251
!       CO2M51   =  TRANSMISSION FCTNS FOR T0 FOR ADJACENT PRESSURE 
!                      LEVELS, WITH NO PRESSURE QUADRATURE. USED FOR
!                      NEARBY LAYER COMPUTATIONS. P(SFC)=1013.25 MB 
!       CO2M58   =  SAME AS CO2M51,WITH P(SFC)= ^810 MB 
!       CDTM51   =  FIRST TEMPERATURE DERIVATIVE OF CO2M51
!       CDTM58   =  FIRST TEMPERATURE DERIVATIVE OF CO2M58
!       C2DM51   =  SECOND TEMPERATURE DERIVATIVE OF CO2M51 
!       C2DM58   =  SECOND TEMPERATURE DERIVATIVE OF CO2M58 
!       STEMP    =  STANDARD TEMPERATURES FOR MODEL PRESSURE LEVEL
!                      STRUCTURE WITH P(SFC)=1013.25 MB 
!       GTEMP    =  WEIGHTING FUNCTION FOR MODEL PRESSURE LEVEL 
!                      STRUCTURE WITH P(SFC)=1013.25 MB.
!       B0       =  TEMP. COEFFICIENT USED FOR CO2 TRANS. FCTN. 
!                      CORRECTION FOR T(K). (SEE REF. 4 AND BD3)
!       B1       =  TEMP. COEFFICIENT, USED ALONG WITH B0 
!       B2       =  TEMP. COEFFICIENT, USED ALONG WITH B0 
!       B3       =  TEMP. COEFFICIENT, USED ALONG WITH B0 

      Real, Allocatable, Dimension(:,:) :: CO251,CO258,CDT51,CDT58
      Real, Allocatable, Dimension(:,:) :: C2D51,C2D58
      Real, Allocatable, Dimension(:)   :: CO2M51,CO2M58,CDTM51,CDTM58
      Real, Allocatable, Dimension(:)   :: C2DM51,C2DM58
      Real, Allocatable, Dimension(:)   :: STEMP,GTEMP
      Real                              :: B0,B1,B2,B3

!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 560-670 cm-1 part of the 15 um co2 band. 
!   This data was formerly in COMMON /CO2BD2/.
!
!       CO231    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO238    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB
!       CDT31    =  FIRST TEMPERATURE DERIVATIVE OF CO231
!       CDT38    =  FIRST TEMPERATURE DERIVATIVE OF CO238
!       C2D31    =  SECOND TEMPERATURE DERIVATIVE OF CO231
!       C2D38    =  SECOND TEMPERATURE DERIVATIVE OF CO231

      Real, Allocatable, Dimension(:) :: CO231,CO238,CDT31,CDT38
      Real, Allocatable, Dimension(:) :: C2D31,C2D38

!-----------------------------------------------------------------------
!
!   co2 transmission functions and temperature and pressure
!   derivatives for the 670-800 part of the 15 um co2 band.
!   This data was formerly in COMMON /CO2BD4/.
!
!       CO271    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO278    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB
!       CDT71    =  FIRST TEMPERATURE DERIVATIVE OF CO271
!       CDT78    =  FIRST TEMPERATURE DERIVATIVE OF CO278
!       C2D71    =  SECOND TEMPERATURE DERIVATIVE OF CO271
!       C2D78    =  SECOND TEMPERATURE DERIVATIVE OF CO271

      Real, Allocatable, Dimension(:) :: CO271,CO278,CDT71,CDT78
      Real, Allocatable, Dimension(:) :: C2D71,C2D78

!-----------------------------------------------------------------------
!
!   co2 transmission functions for the 2270-2380 part of the
!   4.3 um co2 band. THis data was formerly in COMMON /CO2BD5/.
!
!       CO211    =  TRANSMISSION FCTNS FOR T0 (STD. PROFILE)
!                     WITH P(SFC)=1013.25 MB
!       CO218    =  TRANSMISSION FCTNS. FOR T0 (STD. PROFILE)
!                     WITH P(SFC)= ^810 MB

      Real, Allocatable, Dimension(:) :: CO211,CO218


!-----------------------------------------------------------------------
!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: co2_data.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------

      Public   CO2_Data, Write_CO2_Data, Read_CO2_Data, &
               co2_data_init, co2_data_end

      Public   CO251,  CO258,  CDT51,  CDT58,  C2D51,  C2D58,  &
               CO2M51, CO2M58, CDTM51, CDTM58, C2DM51, C2DM58, &
               STEMP,  GTEMP,  B0,     B1,     B2,     B3,     &
               CO231,  CO238,  CDT31,  CDT38,  C2D31,  C2D38,  &
               CO271,  CO278,  CDT71,  CDT78,  C2D71,  C2D78,  &
               CO211,  CO218

      CONTAINS

!#######################################################################

      Subroutine CO2_Data (co2std, ratio, Pref)

      Implicit None
!-----------------------------------------------------------------------
      Real, Intent(IN) :: co2std, ratio
      Real, Intent(IN) :: Pref(:,:)
!-----------------------------------------------------------------------
!   CO2STD = standard co2 vol. mixing ratio (either 300 or 330 ppmv)
!   RATIO  = co2 vol. mixing ratio in units of the standard vol. 
!            mixing ratio (must lie between 0.5 and 4.0)
!   PREF   = reference pressure levels
!-----------------------------------------------------------------------

      call error_mesg('CO2_Data', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine CO2_Data

!#######################################################################

      Subroutine Write_CO2_Data


!---------------------------------------------------------------------

      call error_mesg('Write_CO2_Data', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine Write_CO2_Data

!#######################################################################

      Subroutine Read_CO2_Data (nlev)

      Implicit None
!-----------------------------------------------------------------------
!     Reads co2 transmission functions from file = INPUT/CO2.data
!-----------------------------------------------------------------------
      Integer, Intent(IN) :: nlev


!---------------------------------------------------------------------

      call error_mesg('Read_CO2_Data', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine Read_CO2_Data

!#######################################################################

      Subroutine co2_data_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('co2_data_init', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine co2_data_init

!#######################################################################

      Subroutine co2_data_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('co2_data_end', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine co2_data_end

!#######################################################################



      End Module CO2_Data_Mod


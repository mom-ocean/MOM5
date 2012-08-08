MODULE DIAG_CLOUD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       DIAGNOSTIC CLOUD PREDICTION - Gordon (1992)            
!
!       1999 Feb -> 2000 July
!       Contact persons: Bill Stern (for code structure information)
!                        Tony Gordon (for cloud scheme information)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------
!  Calculates cloud fractions diagnostically using relative humidity,
!  omega and stability 
!-------------------------------------------------------------------

 use       fms_mod, only: error_mesg, FATAL, file_exist,    &
                          check_nml_error, open_namelist_file,    &
                          mpp_pe, mpp_root_pe,  close_file,           &
                          write_version_number, stdlog
 use  Constants_Mod, only: Cp_Air, rdgas, rvgas, Kappa, HLv
 use time_manager_mod, only:  TIME_TYPE
 use  cloud_zonal_mod, only:  CLOUD_ZONAL_INIT, GETCLD
 use  diag_cloud_rad_mod, only:  CLOUD_TAU_DRIVER, diag_cloud_rad_INIT,&
                                 cloud_pres_thick_for_tau,  &
                                 cloud_opt_prop_tg_lw, &
                                 cloud_opt_prop_tg_sw, &
                                 cloud_optical_depths, &
                                 cloud_optical_depths2
 use  sat_vapor_pres_mod, ONLY: ESCOMP
 use  shallow_conv_mod, ONLY: SHALLOW_CONV_INIT,MYLCL

!-----------------------------------------------------------------------
 implicit none
!-----------------------------------------------------------------------

 private


!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: diag_cloud.F90,v 17.0 2009/07/21 02:54:10 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
!-----------------------------------------------------------------------

 logical :: module_is_initialized = .false.


 public diag_cloud_driver, diag_cloud_init, diag_cloud_end
 public diag_cloud_driver2
 public diag_cloud_sum, diag_cloud_avg, diag_cloud_avg2, do_diag_cloud
 public diag_cloud_restart

 contains

!#############################################################################      

 SUBROUTINE DIAG_CLOUD_DRIVER (is,js, &
                    temp,qmix,rhum,omega,lgscldelq,cnvcntq,convprc, &
                    pfull,phalf,psfc,coszen,lat,time, &
                    nclds,cldtop,cldbas,cldamt,r_uv,r_nir,ab_uv,ab_nir, &
                    em_lw, conc_drop, conc_ice, size_drop, size_ice, &
    kbot)

! Arguments (intent in)

 integer, intent(in)   ::  is,js
 type(time_type), intent(in)  :: time
 real, intent(in)  :: lat(:,:)
 real, intent(in), dimension (:,:,:) ::  temp,qmix,rhum,omega
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq,pfull, phalf
 real, intent(in), dimension (:,:)   ::  convprc,psfc, coszen
 
 integer, intent(in), OPTIONAL, dimension(:,:) :: kbot


!      INPUT
!      ------

!      IS,JS    starting i,j indices from the full horizontal grid
!      IX, IY   Horizontal dimensions for global storage arrays
!      TEMP     Temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      convprc Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      PFULL    Pressure at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      PHALF    Pressure at half model levels
!                   (dimensioned IDIM x JDIM x kx+1)
!      PSFC     Surface pressure field
!                   (dimensioned IDIM x JDIM)
!      COSZEN     cosine of the zenith angle
!                   (dimensioned IDIM x JDIM)
!      TIME       time of year (time_type)
!      LAT        latitudes in radians, dimensioned by (1xJDIM)   
!      KBOT      OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
!===================================================================
! Arguments (intent out)

 integer, intent(out), dimension(:,:,:) :: cldtop,cldbas
 integer, intent(out), dimension(:,:)  ::  nclds

   real, intent(out), dimension(:,:,:), optional :: r_uv,r_nir,ab_uv, &
                                                  ab_nir,em_lw, &
                                                  conc_drop, conc_ice, &
                                                  size_drop, size_ice
   real, intent(out), dimension(:,:,:) :: cldamt

!      OUTPUT
!      ------

!       NCLDS   number of (random overlapping) clouds in column and also
!                     the current # for clouds to be operating on
!                   (dimensioned IDIM x JDIM )
!      CLDTOP   index of cloud tops (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDBAS   index of cloud bottoms (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      CLDAMT   cloud amount (fraction) (at cloud levels)
!                   (dimensioned IDIM x JDIM x kx)
!      R_UV     fractional amount of ultraviolet radiation
!                     reflected by the clouds (at cloud levels)
!      R_NIRfractional amount of near inrared radiation
!                     reflected by the clouds (at cloud levels)
!      AB_UVfractional amount of ultraviolet radiation
!                     absorbed by the clouds (at cloud levels)
!      AB_NIRfractional amount of near inrared radiation
!                     absorbed by the clouds (at cloud levels)
!      EM_LWemissivity for the clouds (at cloud levels)



end SUBROUTINE DIAG_CLOUD_DRIVER

!---------------------------------------------------------------------

subroutine diag_cloud_driver2 (is, js, press, pflux, lat, time, nclds, &
                               cldtop, cldbas, cldamt, liq_frac, tau, &
                               ice_cloud, kbot) 

!--------------------------------------------------------------------- 
!    diag_cloud_driver2 returns the cloud specification arrays for the 
!    gordon diag cloud scheme. returned are the number of clouds per 
!    column, the cloud top, cloud base and fractional coverage of each 
!    cloud, the amount of the cloud which is liquid, its optical depth 
!    and an indicator as to whether it is ice or liquid (a different 
!    criterion than is used for the liquid fraction determination).
!----------------------------------------------------------------------
 
integer,                     intent(in)             ::  is,js
real,    dimension (:,:,:),  intent(in)             ::  press, pflux 
real,    dimension(:,:),     intent(in)             ::  lat
type(time_type),             intent(in)             ::  time
integer, dimension(:,:),     intent(inout)          ::  nclds
integer, dimension(:,:,:),   intent(out)            ::  cldtop,cldbas
real,    dimension(:,:,:),   intent(out)            ::  cldamt, liq_frac
real,    dimension(:,:,:,:), intent(out)            ::  tau
logical, dimension(:,:,:),   intent(out)            ::  ice_cloud
integer, dimension(:,:),     intent(in), optional   ::  kbot

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js             starting subdomain i,j indices of data 
!                        in the physics_window being integrated
!      press             pressure at model levels (1:nlev), surface    
!                        pressure is stored at index value nlev+1   
!                        [ (kg /( m s^2) ]
!      pflux             average of pressure at adjacent model levels  
!                        [ (kg /( m s^2) ]
!      lat               latitude of model points  [ radians ]
!      time              time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ]
!
!   intent(inout) variables:
!
!      nclds             total number of clouds in each grid column
!
!   intent(out) variables:
!
!      cldtop            k index of cloud top for each cloud
!      cldbas            k index of cloud base for each cloud
!      cldamt            fractional cloudiness for each cloud
!                        [ dimensionless ]
!      liq_frac          fraction of cloud which is liquid 
!                        [ dimensionless ]
!      tau               cloud optical depth  [ dimensionless ]
!      ice_cloud         logical flag indicating whether cloud is liquid
!                        or ice
!
!    intent(in), optional variables:
!
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground
!
!---------------------------------------------------------------------



end subroutine diag_cloud_driver2


!#######################################################################

  SUBROUTINE DIAG_CLOUD_INIT( ix,iy,kx, ierr )

!=======================================================================
! ***** INITIALIZE Predicted Cloud Scheme
!=======================================================================


!---------------------------------------------------------------------
! Arguments (Intent in)
!  parmameter mxband = max number of radiative bands to be considered for some
!              cloud properties (defined at top of module)
!---------------------------------------------------------------------
 integer, intent(in) :: ix, iy, kx
!      INPUT
!      ------

!      IX, IY, KX   Dimensions for global storage arrays (2- horiz, vert)
!---------------------------------------------------------------------
! Arguments (Intent out)
!---------------------------------------------------------------------
 integer, intent(out) :: ierr

!      OUTPUT
!      ------

!      IERR     Error flag

 
!=====================================================================
  end SUBROUTINE DIAG_CLOUD_INIT

!#######################################################################

  SUBROUTINE DIAG_CLOUD_END

 
!=====================================================================
  end SUBROUTINE DIAG_CLOUD_END

!#######################################################################

 function do_diag_cloud ( ) result (answer)
   logical :: answer

!  returns logical value for whether diag_cloud has been initialized
!  presumably if initialized then diag_cloud will be used

   answer = module_is_initialized

 end function do_diag_cloud

!#######################################################################

 SUBROUTINE DIAG_CLOUD_SUM (is,js, &
                    temp,qmix,rhum,omega,lgscldelq,cnvcntq,convprc,kbot)

!-----------------------------------------------------------------------
 integer, intent(in)                 :: is,js
 real, intent(in), dimension (:,:,:) ::  temp,qmix,rhum,omega
 real, intent(in), dimension (:,:,:) ::  lgscldelq,cnvcntq
 real, intent(in), dimension (:,:)   ::  convprc

 integer, intent(in), OPTIONAL, dimension(:,:) :: kbot

!      INPUT
!      ------

!      IS,JS    starting i,j indices from the full horizontal grid
!      TEMP     Temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Mixing Ratio at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      RHUM     Relative humidity fraction at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IDIM x JDIM x kx)
!      convprc Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IDIM x JDIM)
!      KBOT      OPTIONAL; lowest model level index array
!                   (dimensioned IDIM x JDIM)
! ******* kbot will be used to select only those qmix values that are really
! ******* needed (typically this will be the bottom level except for 
! ******* step mountains
!-----------------------------------------------------------------------


 end SUBROUTINE DIAG_CLOUD_SUM

!#######################################################################

 subroutine DIAG_CLOUD_AVG (is, js, temp,qmix,rhum,omega, &
                           lgscldelq,cnvcntq,convprc,      ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(inout), dimension(:,:,:) :: temp,qmix,rhum,omega
      real, intent(inout), dimension(:,:,:) :: lgscldelq,cnvcntq
      real, intent(inout), dimension(:,:)   :: convprc
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
     
!-----------------------------------------------------------------------

 end SUBROUTINE DIAG_CLOUD_AVG

!#######################################################################

 subroutine DIAG_CLOUD_AVG2 (is, js, qmix, ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(inout), dimension(:,:) :: qmix
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
     
!-----------------------------------------------------------------------

 end SUBROUTINE DIAG_CLOUD_AVG2

!#######################################################################

!#######################################################################
! <SUBROUTINE NAME="diag_cloud_restart">
!
! <DESCRIPTION>
!  Dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine diag_cloud_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

end subroutine diag_cloud_restart
! </SUBROUTINE> NAME="diag_cloud_restart"


end MODULE DIAG_CLOUD_MOD

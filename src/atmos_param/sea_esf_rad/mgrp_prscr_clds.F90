!FDOC_TAG_GFDL

                 module mgrp_prscr_clds_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil   
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!       mgroup prescribed cloud properties module
!               (this module runnable in SKYHI and FMS;
!                zonal_clouds_mod is FMS native equivalent)
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use mpp_mod,           only: input_nml_file
use fms_mod,           only: fms_init, file_exist, &
                             open_namelist_file,  &
                             check_nml_error, close_file,   &
                             write_version_number, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             error_mesg, FATAL

use constants_mod,     only: radian
use rad_utilities_mod, only: rad_utilities_init, &
                             cldrad_properties_type, &
                             cld_specification_type


!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!       mgroup prescribed cloud properties module
!               (this module runnable in SKYHI and FMS; 
!                zonal_clouds_mod is FMS native equivalent)
!
!!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: mgrp_prscr_clds.F90,v 19.0 2012/01/06 20:19:41 fms Exp $'
  character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public     mgrp_prscr_clds_init,  &
           mgrp_prscr_clds_end,  &
          prscr_clds_amt, &      
          obtain_bulk_lw_prscr, obtain_bulk_sw_prscr 


private    cldht, cldint  


!---------------------------------------------------------------------
!-------- namelist  ---------



integer    ::     dummy=0



namelist /mgrp_prscr_clds_nml /     &
                                           dummy


!----------------------------------------------------------------------
!----  public data -------

!----------------------------------------------------------------------
!----  private data -------


!--------------------------------------------------------------------
!     these variables define cloud amounts and radiative properties  
!     on a global (j,k) model grid.
!--------------------------------------------------------------------
  
logical, dimension(:,:),   allocatable :: zhi_cloud_gl, zmid_cloud_gl, &
                                         zlow_cloud_gl
real, dimension(:,:),   allocatable :: zcamtmxo, zcamtrnd

 
!--------------------------------------------------------------------
!     these variables define cloud tops and bottoms, amounts and rad-
!     iative properties on an input (LATOBS,NOFCLDS_SP) grid. the input 
!     values are the original Skyhi values. 
!--------------------------------------------------------------------
 
integer, parameter                    :: NOFCLDS_SP=3  
integer, parameter                    :: NOFMXOLW=1  
integer, parameter                    :: NOFRNDLW=2  
integer, parameter                    :: LATOBS=19

!-------------------------------------------------------------------
!   default low, middle, high cloud properties 
!   cldem    : infrared emissivity
!   crfvis   : visible band reflectivity
!   crfir    : near-ir band reflectivity
!   cabir    : near-ir band absorptivity
!-------------------------------------------------------------------
real  :: cldem_hi = 1.0     
real ::  cldem_mid= 1.0
real ::  cldem_low = 1.0          
real :: crfvis_hi = 0.21 
real :: crfvis_mid =0.48          
real :: crfvis_low =0.69          
real :: crfir_hi = 0.21
real :: crfir_mid =0.48     
real :: crfir_low = 0.69          
real :: cabir_hi =  0.005         
real :: cabir_mid =0.02           
real :: cabir_low =0.035

!-------------------------------------------------------------------
!   prescribed high, mid, low cloud amounts on 19 latitudes
!-------------------------------------------------------------------

integer                                 ::  jj, kkc
real, dimension(LATOBS)                 ::  ccd_low, ccd_mid, ccd_hi  
real, dimension(LATOBS)                 ::  cloud_lats

data cloud_lats / -90., -80., -70., -60., -50., -40., -30., -20., &
                  -10., 0.0, 10., 20., 30., 40., 50., 60., 70., 80., &
                   90. /

data ccd_low  /    &
     &  0.360E+00, 0.401E+00, 0.439E+00, 0.447E+00, 0.417E+00, &
     &  0.343E+00, 0.269E+00, 0.249E+00, 0.290E+00, 0.330E+00, &
     &  0.290E+00, 0.249E+00, 0.269E+00, 0.343E+00, 0.417E+00, &
     &  0.447E+00, 0.439E+00, 0.401E+00, 0.360E+00/               

data ccd_mid  /    &
     &  0.090E+00, 0.102E+00, 0.117E+00, 0.128E+00, 0.122E+00, &
     &  0.095E+00, 0.070E+00, 0.060E+00, 0.068E+00, 0.080E+00, &
     &  0.068E+00, 0.060E+00, 0.070E+00, 0.095E+00, 0.122E+00, &
     &  0.128E+00, 0.117E+00, 0.102E+00, 0.090E+00/              

data ccd_hi   /    &
     &  0.198E+00, 0.231E+00, 0.254E+00, 0.250E+00, 0.227E+00, &
     &  0.192E+00, 0.159E+00, 0.168E+00, 0.205E+00, 0.241E+00, &
     &  0.205E+00, 0.168E+00, 0.159E+00, 0.192E+00, 0.227E+00, &
     &  0.250E+00, 0.254E+00, 0.231E+00, 0.198E+00/ 
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!    NLWCLDB is the number of frequency bands for which lw
!    emissitivies are defined.
!----------------------------------------------------------------------

logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!----------------------------------------------------------------------




 contains 


! <SUBROUTINE NAME="mgrp_prscr_clds_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mgrp_prscr_clds_init (    pref, latb      )
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine mgrp_prscr_clds_init (    pref, latb      )

!------------------------------------------------------------------
real, dimension(:), intent(in)             ::  latb      
real, dimension(:,:), intent(in)             :: pref          
!--------------------------------------------------------------------

      integer            :: unit, ierr, io, logunit
      integer            :: j, k,                 li
      integer            :: jdf
      integer, dimension ( LATOBS, NOFCLDS_SP)  :: kkbh, kkth
         
 integer :: kerad

      integer, dimension(:), allocatable :: jindx2

      if (module_is_initialized) return

      call fms_init
      call rad_utilities_init

!---------------------------------------------------------------------
!-----  read namelist  ------
  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=mgrp_prscr_clds_nml, iostat=io)
      ierr = check_nml_error(io,"mgrp_prscr_clds_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=mgrp_prscr_clds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'mgrp_prscr_clds_nml')
        enddo
10      call close_file (unit)
      endif
#endif

      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )   &
         write (logunit, nml = mgrp_prscr_clds_nml)

!--------------------------------------------------------------------
!  retrieve module variables that come from other modules
!--------------------------------------------------------------------

       kerad = size(pref,1) - 1
      jdf = size(latb,1) - 1

         allocate (jindx2  (size(latb,1)-1)) ; jindx2 = 0.0
         call find_nearest_index (latb, jindx2)

!---------------------------------------------------------------------
!    define the number of cloud emissivity bands for use in this module.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!   allocate space to hold the cloud radiative properties on the model
!   grid (j,k)
!---------------------------------------------------------------------
      allocate( zcamtmxo    (jdf , KERAD)) ; zcamtmxo = 0.0
      allocate( zcamtrnd    (jdf , KERAD)) ; zcamtrnd = 0.0

      allocate( zhi_cloud_gl(jdf ,     kerad  ))
      allocate( zmid_cloud_gl(jdf ,    kerad   ))
      allocate( zlow_cloud_gl(jdf ,    kerad   ))
      zhi_cloud_gl  = .false.
      zmid_cloud_gl = .false.
      zlow_cloud_gl = .false.

!---------------------------------------------------------------------
!     define index arrays for cloud tops and cloud bottoms(kkth, kkbh). 
!     a program courtesy of Ron Stouffer is used. cldht specifies cloud 
!     height data in mb for each cloud type, at 10 deg intervals from 
!     90S to 90N.(ie, 19 lats).
!---------------------------------------------------------------------
 
      call cldht  (pref(:,1), kkbh, kkth)


!---------------------------------------------------------------------
!    perform latitude "interpolation" to the (JD) model latitudes
!    in default case, no interpolation is actually done; the nearest
!    latitude available (using NINT function) is used.
!---------------------------------------------------------------------

       zcamtmxo = 0.
       zcamtrnd = 0.

      do j=1,jdf
        li = jindx2(j)
do k=kkth(li,1), kkbh(li,1)
           zlow_cloud_gl(j,k) = .true.
  zcamtmxo(j,k) = ccd_low(li)        
        end do
        do k=kkth(li,2), kkbh(li,2)
           zmid_cloud_gl(j,k) = .true.
   zcamtrnd(j,k) = ccd_mid(li)
         end do
         do k=kkth(li,3), kkbh(li,3)
            zhi_cloud_gl(j,k) = .true.
    zcamtrnd(j,k) = ccd_hi(li)
        end do
      end do

!---------------------------------------------------------------------
!    if a cloud microphysics scheme is to be employed with the cloud
!    scheme, initialize the microphysics_rad module.
!--------------------------------------------------------------------

      module_is_initialized = .true.

end subroutine mgrp_prscr_clds_init

!######################################################################

! <SUBROUTINE NAME="mgrp_prscr_clds_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mgrp_prscr_clds_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine mgrp_prscr_clds_end
        
!----------------------------------------------------------------------
!    mgrp_prscr_clds_end is the destructor for mgrp_prscr_clds_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
        
!--------------------------------------------------------------------
 
 
end subroutine mgrp_prscr_clds_end




!#####################################################################

! <SUBROUTINE NAME="find_nearest_index">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call find_nearest_index (latb, jindx2)
!
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
!  <OUT NAME="jindx2" TYPE="integer">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine find_nearest_index (latb, jindx2)

real, dimension(:), intent(in) :: latb
integer, dimension(:), intent(out)  :: jindx2


      integer :: jd, j, jj
      real   :: diff_low, diff_high
      real, dimension(size(latb,1)-1) :: lat


      jd = size(latb,1) - 1

      do j = 1,jd
        lat(j) = 0.5*(latb(j) + latb(j+1))
      do jj=1, LATOBS         
         if (lat(j)*radian >= cloud_lats(jj)) then
           diff_low = lat(j)*radian - cloud_lats(jj)
           diff_high = cloud_lats(jj+1) - lat(j)*radian
          if (diff_high <= diff_low) then
            jindx2(j) = jj+1
          else
            jindx2(j) = jj
          endif
        endif
      end do
      end do






end subroutine find_nearest_index 





!######################################################################

! <SUBROUTINE NAME="prscr_clds_amt">
!  <OVERVIEW>
!    prscr_clds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    prscr_clds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call prscr_clds_amt (is, ie, js, je, Cld_spec)
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
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine prscr_clds_amt (is, ie, js, je, Cld_spec)

!---------------------------------------------------------------------
!    prscr_clds_amt defines the location, amount (cloud fraction), 
!    number and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------

integer, intent(in)                          :: is, ie, js, je
type(cld_specification_type), intent(inout)  :: Cld_spec       

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
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
 
      integer  ::  i,j,k    ! do loop indices


!---------------------------------------------------------------------
!    define the number of clouds in each column. the assumption is made
!    that all grid columns have NOFCLDS_SP clouds, with NOFMXOLW of 
!    these being maximally overlapped and the remainder (NOFRNDLW) 
!    randomly overlapped. The default case, corresponding to all 
!    previous model simulations, is 3 clouds (1 maximally overlapped, 
!    2 randomly overlapped).
!----------------------------------------------------------------------
      Cld_spec%nmxolw(:,:) = NOFMXOLW
      Cld_spec%nrndlw(:,:) = NOFRNDLW
      Cld_spec%ncldsw(:,:) = NOFCLDS_SP 

!----------------------------------------------------------------------
!    define the fractions of random and maximally overlapped clouds and
!    the total cloud fraction seen by the shortwave radiation. define
!    the hi-mid-low cloud flag arrays based on the pre-defined latitude
!    and height dependent specifications that were defined during init-
!    ialization.
!----------------------------------------------------------------------
      do k=1,size (Cld_spec%hi_cloud,3)
        do j=1,size (Cld_spec%hi_cloud,2)
          do i=1,size (Cld_spec%hi_cloud,1)
            Cld_spec%cmxolw(i,j,k) = zcamtmxo(j+js-1  ,k)
            Cld_spec%crndlw(i,j,k) = zcamtrnd(j+js-1  ,k)
            Cld_spec%camtsw(i,j,k) = Cld_spec%cmxolw(i,j,k) +  &
                                     Cld_spec%crndlw(i,j,k)
            if (Cld_spec%camtsw(i,j,k) > 0.0) then
              if (zhi_cloud_gl(js+j-1,k)) then
                Cld_spec%hi_cloud(i,j,k) = .true.
              else if (zmid_cloud_gl(js+j-1,k)) then
                Cld_spec%mid_cloud(i,j,k) = .true.
              else if (zlow_cloud_gl(js+j-1,k) ) then
                Cld_spec%low_cloud(i,j,k) = .true.
              else
                call error_mesg ('mgrp_prscr_clds_mod',  &
                    'model level is not mapped to a cloud type', FATAL)
              endif
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------



end subroutine prscr_clds_amt


!######################################################################

! <SUBROUTINE NAME="obtain_bulk_lw_prscr">
!  <OVERVIEW>
!    obtain_bulk_lw_prscr defines bulk longwave cloud radiative
!    properties for the mgrp_prscr_clds cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_lw_prscr defines bulk longwave cloud radiative
!    properties for the mgrp_prscr_clds cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_lw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)
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
subroutine obtain_bulk_lw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_lw_prscr defines bulk longwave cloud radiative 
!    properties for the mgrp_prscr_clds cloud scheme.
!---------------------------------------------------------------------

integer,                     intent(in)     :: is, ie, js, je
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
 
      integer   ::    i,j,k   ! do-loop indices

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
          end do
        end do
      end do

!---------------------------------------------------------------------



end subroutine obtain_bulk_lw_prscr 



!#####################################################################

! <SUBROUTINE NAME="obtain_bulk_sw_prscr">
!  <OVERVIEW>
!    obtain_bulk_sw_zonal defines bulk shortwave cloud radiative
!    properties for the zonal cloud scheme.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_bulk_sw_zonal defines bulk shortwave cloud radiative
!    properties for the zonal cloud scheme.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_bulk_sw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)
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
subroutine obtain_bulk_sw_prscr (is, ie, js, je, Cld_spec, Cldrad_props)

!---------------------------------------------------------------------
!    obtain_bulk_sw_zonal defines bulk shortwave cloud radiative 
!    properties for the zonal cloud scheme.
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

!----------------------------------------------------------------------



end subroutine obtain_bulk_sw_prscr 




!######################################################################




!##################################################################


! <SUBROUTINE NAME="cldht">
!  <OVERVIEW>
!  This subroutine computes the heights of the cloud tops
!  and bottoms for the fixed cloud model.  The observed data
!  are from London (1954, 1957).  This data is a function of 10 deg.
!  latitude bands (0-10, 10-20 and etc.), season and height
!  in the orginal paper and only for
!  the Northern Hemisphere for various cloud types.
!  Dick and Suki averaged the four seasons together to get annual
!  mean cloud heights for three type of clouds (hi, middle and low).
!  Somebody also interpolated the data from the 10 deg latitude
!  bands to 5 deg bands.  At the equator, this interpolation
!  was more like an extrapolation.
!  These heights were then put in pressure coordinates using
!  a Skew-T diagram which assumes a "standard atmosphere".
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!  This subroutine computes the heights of the cloud tops
!  and bottoms for the fixed cloud model.  The observed data
!  are from London (1954, 1957).  This data is a function of 10 deg.
!  latitude bands (0-10, 10-20 and etc.), season and height
!  in the orginal paper and only for
!  the Northern Hemisphere for various cloud types.
!  Dick and Suki averaged the four seasons together to get annual
!  mean cloud heights for three type of clouds (hi, middle and low).
!  Somebody also interpolated the data from the 10 deg latitude
!  bands to 5 deg bands.  At the equator, this interpolation
!  was more like an extrapolation.
!  These heights were then put in pressure coordinates using
!  a Skew-T diagram which assumes a "standard atmosphere".
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cldht (plevel, kkbh, kkth)
!
!  </TEMPLATE>
!  <IN NAME="plevel" TYPE="real">
! 
!  </IN>
!  <OUT NAME="kkbh" TYPE="integer">
! 
!  </OUT>
!  <OUT NAME="kkth" TYPE="integer">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine cldht (plevel, kkbh, kkth)

real, dimension(:), intent(in) :: plevel
integer, dimension(:,:), intent(out) :: kkbh, kkth
 
!------------------------------------------------------------------
!  This subroutine computes the heights of the cloud tops
!  and bottoms for the fixed cloud model.  The observed data
!  are from London (1954, 1957).  This data is a function of 10 deg.
!  latitude bands (0-10, 10-20 and etc.), season and height 
!  in the orginal paper and only for
!  the Northern Hemisphere for various cloud types.
!  Dick and Suki averaged the four seasons together to get annual 
!  mean cloud heights for three type of clouds (hi, middle and low).
!  Somebody also interpolated the data from the 10 deg latitude
!  bands to 5 deg bands.  At the equator, this interpolation
!  was more like an extrapolation.
!  These heights were then put in pressure coordinates using
!  a Skew-T diagram which assumes a "standard atmosphere".
!
!  The Dick and Suki (Ron's pressures) data follow:
!
! TYPE:      hi      |   middle    |            low
!     |              |             |     ctop    |   cbase
! Lat |   ht   press |  ht   press |  ht   press |  ht   press
! 0   |  9.80  272   | 4.35  590   | 3.00  700   | 1.40  855
! 5   |  9.82  271   | 4.40  586   | 3.04  698   | 1.47  848
! 10  | 10.13  259   | 4.45  581   | 3.08  693   | 1.61  836
! 15  | 10.35  250   | 4.50  578   | 3.08  693   | 1.70  828 
! 20  | 10.50  244   | 4.50  578   | 3.01  699   | 1.72  825
! 25  | 10.50  244   | 4.41  584   | 2.91  710   | 1.71  826
! 30  | 10.38  248   | 4.26  595   | 2.80  719   | 1.70  828
! 35  | 10.03  263   | 4.10  614   | 2.70  729   | 1.65  830
! 40  |  9.44  285   | 3.92  621   | 2.60  735   | 1.58  839
! 45  |  8.65  322   | 3.79  633   | 2.47  750   | 1.50  846
! 50  |  7.97  357   | 3.67  647   | 2.35  760   | 1.40  859
! 55  |  7.55  379   | 3.56  651   | 2.24  770   | 1.31  867
! 60  |  7.29  392   | 3.51  657   | 2.17  780   | 1.25  871
! 65  |  7.13  401   | 3.50  658   | 2.10  783   | 1.20  875
! 70  |  7.03  406   | 3.48  659   | 2.03  788   | 1.12  881
! 75  |  7.01  409   | 3.44  660   | 1.98  795   | 1.05  890
! 80  |  6.99  410   | 3.43  661   | 1.91  800   | 1.02  891
! 85  |  6.98  411   | 3.43  661   | 1.88  803   | 1.00  896
! 90  |  6.98  411   | 3.43  661   | 1.87  804   | 1.00  896
!
!  Note that the heights are in kilometers and the pressures
!  in millibars.
!--------------------------------------------------------------------

      real, dimension(10)   :: ciht, asht, cltop, clbase
      integer               :: n, nl, j
 
!------------------------------------------------------------------
!  The data in the following arrays run pole to equator,
!  starting at 90, 80, 70....10, 0.
!------------------------------------------------------------------
 
!  Observed high cloud heights running pole to equator (mb)
      data ciht / 411, 410, 406, 392, 357, 285, 248, 244, 259, 272 /
 
!  Observed middle cloud heights running pole to equator (mb)
      data asht / 661, 661, 659, 657, 647, 621, 595, 578, 581, 590 /
 
!  Observed low cloud top heights running pole to equator (mb)
!     data cltop / 804, 800, 788, 780, 760, 735, 719, 699, 693, 700 /
      data cltop / 804, 800, 788, 780, 766, 735, 719, 699, 693, 700 /
!  The above data statement was changed so that this code would
!  reproduce exactly the indexes used in the 9 level model.  The 6 mb
!  error should not affect the results very much....if at all.
 
!  Observed low cloud bottom heights running pole to equator (mb)
      data clbase / 896, 891, 881, 871, 859, 839, 828, 825, 836, 855 /
 
 
!--------------------------------------------------------------------
!   for (at present) unexplained reasons, the middle cloud 
!   specification does not agree with the current gcm specification.
!   a glance at the telegadas and london paper suggests differences
!   caused by 1) use of cldmid rather than clotop and cldbase heights
!   2) seeming errors in polar cloud height locations. for now,
!   the foregoing kluge gives cloud positions in agreement with
!   SKYHI, PROVIDED:
!    it is realized that these indexes are the true values, which
!   differ by 1 from the skyhiindices. in the gcm, a subtraction
!   is done to get the indices to be correct (see radmn.F).
!--------------------------------------------------------------------

      do n=8,10
asht(n) = asht(n) - 25.
      enddo
 
      clbase = clbase *1.0E02
      cltop  = cltop  *1.0E02
     asht  = asht*100.
     ciht = ciht*100.

 
!-------------------------------------------------------------------
!  First compute the cloud top indexes for the high clouds
!  nl is the index for cloud type nl=3 => high
!--------------------------------------------------------------------  

      nl = 3
      call Cldint(plevel, ciht ,kkth, nl)
 
!---------------------------------------------------------------------
!  Set cloud bottom height equal to cloud top height.
!  This assumes cirrus clouds are one level thick.
!---------------------------------------------------------------------

      do j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
      end do
 
!-------------------------------------------------------------------
!  Second compute the cloud top indexes for the middle clouds
!  nl is the index for cloud type nl=2 => middle
!-------------------------------------------------------------------
  
      nl = 2
      call Cldint(plevel, asht,kkth,  nl)
 
!-------------------------------------------------------------------
!  Set cloud bottom height equal to cloud top height.
!  This assumes middle clouds are one level thick.
!-------------------------------------------------------------------

      do j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
      end do
 
!-------------------------------------------------------------------
!  Third compute the cloud top indexes for the low clouds
!  nl is the index for cloud type nl=1 => low
!-------------------------------------------------------------------
  
      nl = 1
      call Cldint(plevel, cltop,kkth,  nl)
 
!-------------------------------------------------------------------
!  Lastly compute the cloud bottom indexes for the low clouds
!  This assumes that low clouds can be thicker than one level.
!  nl is the index for cloud type nl=1 => low
!-------------------------------------------------------------------
  
      nl = 1
      call Cldint(plevel, clbase, kkbh, nl)
 
end subroutine cldht


!##################################################################

! <SUBROUTINE NAME="cldint">
!  <OVERVIEW>
!  This subroutine computes the indexes for the heights of the cloud
!  tops and bottoms for the fixed cloud model.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!  This subroutine computes the indexes for the heights of the cloud
!  tops and bottoms for the fixed cloud model.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cldint(plevel, cldobs, kindex, nl)
!
!  </TEMPLATE>
!  <IN NAME="plevel" TYPE="real">
! 
!  </IN>
!  <IN NAME="cldobs" TYPE="real">
! 
!  </IN>
!  <OUT NAME="kindex" TYPE="integer">
! 
!  </OUT>
!  <IN NAME="nl" TYPE="integer">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine cldint(plevel, cldobs, kindex, nl)
 
!-------------------------------------------------------------------
!  This subroutine computes the indexes for the heights of the cloud 
!  tops and bottoms for the fixed cloud model. 
!-------------------------------------------------------------------
 
real, dimension(10), intent(in)                    :: cldobs
integer, dimension(LATOBS,NOFCLDS_SP), intent(out) :: kindex
integer,                               intent(in)  :: nl
real, dimension(:), intent(in) :: plevel
!------------------------------------------------------------------
 
      real, dimension(LATOBS)     :: cldlat
      real                        :: prsmid
      integer                     :: j, k
      integer                     :: kerad
 
      kerad = size(plevel(:)) - 1
!---------------------------------------------------------------------
!  Fill in Southern hemisphere cloud heights
!---------------------------------------------------------------------

      do j=1,10
       cldlat(j) = cldobs(j)
      end do  
      do j=1,9
       cldlat(j+10) = cldobs(10-j)
      end do  
  
!-------------------------------------------------------------------
!  Start latitude loop to compute index at each latitude.
!-------------------------------------------------------------------
      do j=1,LATOBS
 
!-------------------------------------------------------------------
!  Find first place where the pressure on the model level
!  is greater than the pressure of the cloud height.  Starting
!  from the top of the atm and going down the column.
!-------------------------------------------------------------------
        if (plevel(KERAD)       .lt. cldlat(j)) then       
        call error_mesg ('mgrp_prscr_clds_mod',  &
         'no level found with pressure greater than cloud pressure', & 
FATAL)
endif
        if (plevel(1)       .gt. cldlat(j)) then       
        call error_mesg ('mgrp_prscr_clds_mod',  &
                ' cloud is above highest model level', FATAL)
        endif
do k=1,KERAD
          if (plevel(k) .gt. cldlat(j)) then       
!-------------------------------------------------------------------
!  k is the index of the first model level below the cloud height.
!  compute the pressure half way between the model levels
!-------------------------------------------------------------------
            prsmid = (plevel(k)+plevel(k-1))      *0.5   
 
!-------------------------------------------------------------------
!  If prsmid is greater than cldlat (cloud height) then the
!  level above is closer to cloud height, otherwise it is the
!  level below.
!-------------------------------------------------------------------
            if (prsmid .gt. cldlat(j)) then
              kindex(j,nl) = k-1
            else
              kindex(j,nl) = k
            endif
            exit
  endif
        end do
      end do
!--------------------------------------------------------------------

 

end subroutine cldint




!####################################################################

       end module mgrp_prscr_clds_mod





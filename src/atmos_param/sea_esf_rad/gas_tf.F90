                       module gas_tf_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Module that calculates gas transmission functions
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>

!  shared modules:

use mpp_mod,             only : input_nml_file
use fms_mod,             only : open_namelist_file, fms_init, &
                                mpp_pe, mpp_root_pe, stdlog, &
                                file_exist, write_version_number, &
                                check_nml_error, error_mesg, &
                                FATAL, close_file, &
                                open_restart_file
use constants_mod,       only : constants_init, RDGAS, GRAV, pstd

!   shared radiation package modules:

use rad_utilities_mod,   only : rad_utilities_init, Lw_control, &
                                atmos_input_type, gas_tf_type
use longwave_params_mod, only : longwave_params_init, NBCO215

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    gas_tf_mod is the gas transmission functions module.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: gas_tf.F90,v 19.0 2012/01/06 20:16:25 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!------    interfaces   ------

 public    &
         gas_tf_init, co2coef, transcol,    &
         transcolrow, trans_nearby, trans_sfc,   &
         put_co2_stdtf_for_gas_tf,  &
         put_co2_nbltf_for_gas_tf, &
         put_ch4_stdtf_for_gas_tf,  &
         put_n2o_stdtf_for_gas_tf, &
         get_control_gas_tf, &
         process_co2_input_file, &
         process_ch4_input_file,  &
         process_n2o_input_file, &
         gas_tf_dealloc, gas_tf_end

private   &

!  called from gas_tf_init:
         ptz,    &

!  called from ptz:
         antemp,   &

!  called from process_co2_input_file,  
!  process_ch4_input_file and process_n2o_input_file:
         process_gas_input_file, &

!  called from co2coef:
         transfn


!---------------------------------------------------------------------
!------     namelist  -----

character(len=16) :: interp_form='log'
logical           :: do_calcstdco2tfs = .true.,   &
                     do_writestdco2tfs = .false., &
                     do_readstdco2tfs = .false.
logical           :: do_calcstdch4tfs = .true.,   &
                     do_writestdch4tfs = .false., &
                     do_readstdch4tfs = .false.
logical           :: do_calcstdn2otfs = .true.,   &
                     do_writestdn2otfs = .false., &
                     do_readstdn2otfs = .false.


namelist / gas_tf_nml /  &
                            interp_form, &
                            do_calcstdch4tfs,  &
                            do_writestdch4tfs, &
                            do_readstdch4tfs, &
                            do_calcstdn2otfs,  &
                            do_writestdn2otfs,  &
                            do_readstdn2otfs, &
                            do_calcstdco2tfs,  &
                            do_writestdco2tfs,  &
                            do_readstdco2tfs

!---------------------------------------------------------------------
!---- public data -------


!---------------------------------------------------------------------
!----   private data  -------


!--------------------------------------------------------------------- 
!    the following arrays are co2 transmission functions, temperature
!    and pressure derivatives for the 560-800 cm-1 band, and standard 
!    temperature and weighting functions.
! 
!      co251   =  transmission functions for t0 (standard profile)
!                 with p(surface)=1013.25 mb.
!
!      co258   =  transmission functions for t0 (standard profile) 
!                 with p(surface)=810 mb.
!
!      cdt51   =  first temperature derivative of co251.
! 
!      cdt58   =  first temperature derivative of co258.
! 
!      c2d51   =  second temperature derivative of co251.
!
!      c2d58   =  second temperature derivative of co258.
!
!      co2m51  =  transmission functions for t0 for adjacent pressure 
!                 levels, with no pressure quadrature.  used for nearby
!                 layer computations.  p(surface)=1013.25 mb.
! 
!      co2m58  =  transmission functions for t0 for adjacent pressure 
!                 levels, with no pressure quadrature.  used for nearby
!                 layer computations.  p(surface)=810 mb.
!
!      cdtm51  =  first temperature derivative of co2m51.
!
!      cdtm58  =  first temperature derivative of co2m58.
!
!      c2dm51  =  second temperature derivative of co2m51.
! 
!      c2dm58  =  second temperature derivative of co2m58.
!--------------------------------------------------------------------- 

real, allocatable, dimension (:,:)       ::  co251, co258,     &    
                                             cdt51, cdt58,     &    
                                             c2d51, c2d58
real, allocatable, dimension (:)         ::  co2m51, co2m58,   &    
                                             cdtm51, cdtm58,   &    
                                             c2dm51, c2dm58

!--------------------------------------------------------------------- 
!    the following arrays are co2 transmission functions for the 2270- 
!    2380 cm-1 part of the 4.3 um co2 band.
! 
!       co211    =  transmission functions for t0 (standard profile)
!                   with p(surface)=1013.25 mb.
!
!       co218    =  transmission functions for t0 (standard profile) 
!                   with p(surface)=810 mb.
!--------------------------------------------------------------------- 

real, allocatable, dimension (:)         ::  co211, co218

!--------------------------------------------------------------------- 
!    the following arrays are co2 transmission functions and temperature
!    and pressure derivatives for (NBCO215) narrow bands in the 15um
!    co2 band.
!
!        co215nbps1    =  transmission functions for USSTD profile
!                        with p(surface)=1013.25 mb.
!        co215nbps8    =  transmission functions for USSTD profile
!                        with p(surface)=810.2 mb.
!
!        co2dt15nbps1  = temperature derivative of co215nbps1.
!
!        co2dt15nbps8  = temperature derivative of co215nbps8.
!
!        co2d2t15nbps1 = second temperature derivative of co215nbps1.
!
!        co2d2t15nbps8 = second temperature derivative of co215nbps8.
!--------------------------------------------------------------------- 

real, allocatable, dimension (:,:)       ::  co215nbps1,       &    
                                             co215nbps8,       &    
                                             co2dt15nbps1,     &    
                                             co2dt15nbps8,     &    
                                             co2d2t15nbps1,    &    
                                             co2d2t15nbps8

!--------------------------------------------------------------------- 
!    the following arrays are ch4 and n2o transmission functions for
!    the 1200-1400 cm-1 band.
! 
!       ch451   =  ch4 transmission functions for t0 (standard profile)
!                  with p(surface)=1013.25 mb.
!
!       ch458   =  ch4 transmission functions for t0 (standard profile) 
!                  with p(surface)=810 mb.
!
!       ch4dt51 =  first temperature derivative of ch4 transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       ch4d2t51=  second temperature derivative of ch4 transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       ch4dt58 =  first temperature derivative of ch4 transmission
!                  functions for t0 profile with p(surface)=810 mb.
!
!       ch4d2t58=  second temperature derivative of ch4 transmission
!                  functions for t0 profile with p(surface)=810 mb.
!
!       n2o51   =  n2o transmission functions for t0 (standard profile)
!                  with p(surface)=1013.25 mb.
!
!       n2o58   =  n2o transmission functions for t0 (standard profile) 
!                  with p(surface)=810 mb.
!
!       n2odt51 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2od2t51=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2odt58 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!
!       n2od2t58=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!--------------------------------------------------------------------- 
  
real, allocatable, dimension (:,:)       ::  ch451, ch458,     &    
                                             ch4dt51, ch4dt58, &    
                                             ch4d2t51, ch4d2t58
real, allocatable, dimension (:,:)       ::  n2o51, n2o58,     &    
                                             n2odt51, n2odt58, &    
                                             n2od2t51, n2od2t58

!--------------------------------------------------------------------- 
!    the following arrays are n2o transmission functions for
!    the 560-630 cm-1 band.
! 
!       n2o71   =  n2o transmission functions for t0 (standard profile)
!                  with p(surface)=1013.25 mb.
!
!       n2o78   =  n2o transmission functions for t0 (standard profile) 
!                  with p(surface)=810 mb.
!
!       n2odt71 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2od2t71=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2odt78 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!
!       n2od2t78=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!--------------------------------------------------------------------- 
!
  
real, allocatable, dimension (:,:)       ::  n2o71, n2o78,     &    
                                             n2odt71, n2odt78, &    
                                             n2od2t71, n2od2t78

!--------------------------------------------------------------------- 
!    the following arrays are n2o transmission functions for
!    the 1070-1200 cm-1 band.
! 
!       n2o91   =  n2o transmission functions for t0 (standard profile)
!                  with p(surface)=1013.25 mb.
!
!       n2o98   =  n2o transmission functions for t0 (standard profile) 
!                  with p(surface)=810 mb.
!
!       n2odt91 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2od2t91=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=1013.25 mb.
!
!       n2odt98 =  first temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!
!       n2od2t98=  second temperature derivative of n2o transmission
!                  functions for t0 profile with p(surface)=810 mb.
!--------------------------------------------------------------------- 
  
real, allocatable, dimension (:,:)       ::  n2o91, n2o98,     &    
                                             n2odt91, n2odt98, &    
                                             n2od2t91, n2od2t98

!----------------------------------------------------------------------
!    stemp   =  standard temperatures for model pressure level
!               structure with p(surface)=1013.25 mb.
!    gtemp   =  weighting function for model pressure level 
!               structure with p(surface)=1013.25 mb.
!----------------------------------------------------------------------
real, dimension (:),       allocatable   :: stemp, gtemp 

!----------------------------------------------------------------------
!    define 4 coefficients (formerly in Id3).
!    b0, b1, b2, b3 are coefficients used to correct for the use of
!    250k in the planck function used in evaluating planck-weighted co2
!    transmission functions. (see reference(1).)
!----------------------------------------------------------------------
real     :: b0 = -0.51926410E-04
real     :: b1 = -0.18113332E-03
real     :: b2 = -0.10680132E-05
real     :: b3 = -0.67303519E-07

integer               :: ksrad, kerad
integer, parameter    :: nvalids=1
integer               :: ixprkminh2o
logical               :: do_linearlblint, do_loglblint
character(len=8)      :: co2_name_save, ch4_name_save, n2o_name_save
real                  :: co2_amount_save, ch4_amount_save, &
                         n2o_amount_save
integer               :: nstdlvls_save, kbegin_save, kend_save

real, dimension(:), allocatable      :: pa_save, pd_save, plm_save
character(len=8), dimension(nvalids) :: valid_versions= 'v1.00'
logical :: module_is_initialized = .false.  ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!###############################################################
! <SUBROUTINE NAME="gas_tf_init">
!  <OVERVIEW>
!   Initialize gas transmission function calculation from input 
!   namelist, model pressure coordinate system, etc.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initialize gas transmission function calculation from input 
!   namelist, model pressure coordinate system, etc.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call gas_tf_init(pref)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   Model pressure coordinate array
!  </IN>
! </SUBROUTINE>
!
subroutine gas_tf_init (pref)                        

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

real, dimension(:,:), intent(in) :: pref

!---------------------------------------------------------------------
!  intent(in)  variables:
!
!     pref
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables;

      real, dimension (size(pref,1)) :: plm
      real, dimension (size(pref,1)) :: pd

      real     ::  prkminh2o = 28.0
      integer  ::  kmin, kmax
      integer  ::  ks = 1
      integer  ::  unit, ierr, io, logunit
      integer  ::  k

!--------------------------------------------------------------------
!  local variables:
!
!       plm
!       pd
!       prkminh2o   pressure above which h2o-co2 overlap affects
!                   nearby layer transmissivities  [ mb ]
!       kmin
!       kmax
!       ks
!       unit
!       ierr
!       io
!       k
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
 !---------------------------------------------------------------------
 !    verify that modules used by this module that are not called later
 !    have already been initialized.
 !---------------------------------------------------------------------
       call fms_init
       call constants_init
       call rad_utilities_init
       call longwave_params_init

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=gas_tf_nml, iostat=io)
      ierr = check_nml_error(io,"gas_tf_nml")
#else
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=gas_tf_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'gas_tf_nml')
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
                       write (logunit, nml=gas_tf_nml)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      kmin = 1
      kerad = size(pref  ,1) - 1
      kmax = kerad
      ksrad = 1
      ks = 1

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      pd (:) = pref(:,1)
      plm (kmin) = 0.
      do k=kmin+1,kmax
        plm (k) = 0.5*(pd (k-1) + pd (k))
      enddo
      plm (kmax+1) = pd (kmax+1)

!---------------------------------------------------------------------
!    convert plm to mb.
!---------------------------------------------------------------------
      plm = plm   *1.0E-02
      pd = pd   *1.0E-02

!--------------------------------------------------------------------
!    check on consistency between namelist values
!--------------------------------------------------------------------
      if ( (Lw_control%do_ch4lbltmpint) .and.    &
           (.not.(Lw_control%do_ch4)) ) then
        call error_mesg ( 'gas_tf_mod', &
       'cannot have do_ch4lbltmpint active when do_ch4 is off',& 
                                                                 FATAL)
      endif
      if ( (Lw_control%do_n2olbltmpint) .and.    &
           (.not.(Lw_control%do_n2o)) ) then
        call error_mesg ( 'gas_tf_mod', &
       'cannot have do_n2olbltmpint active when do_n2o is off',& 
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(interp_form) == 'log') then
        do_loglblint    = .true.
        do_linearlblint = .false.
      else if (trim(interp_form) == 'linear') then
        do_loglblint    = .false.
        do_linearlblint = .true.
      else
        call error_mesg ( 'gas_tf_mod', &
       'improper specification for gas lbl interpolation scheme', FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Lw_control%do_co2_iz) then
        if (Lw_control%do_co2) then
          if (do_writestdco2tfs .and. do_readstdco2tfs) then
            call error_mesg ( 'gas_tf_mod', &
                ' cannot read and write std tfs in same job', FATAL)
          endif

          if (do_writestdco2tfs .and. .not. do_calcstdco2tfs) then
            call error_mesg ( 'gas_tf_mod', &
              ' cannot write std tfs without calculating them', FATAL)
          endif
        endif
      else
        call error_mesg ('gas_tf_mod', &
              'do_co2 has not yet been defined', FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_iz) then
        if (Lw_control%do_ch4) then
          if (do_writestdch4tfs .and. do_readstdch4tfs) then
            call error_mesg ( 'gas_tf_mod', &
              ' cannot read and write std tfs in same job', FATAL)
          endif
          if (do_writestdch4tfs .and. .not. do_calcstdch4tfs) then
            call error_mesg ( 'gas_tf_mod', &
              ' cannot write std tfs without calculating them', FATAL)
          endif
        endif
      else
        call error_mesg ('gas_tf_mod', &
                 'do_ch4 has not yet been defined', FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Lw_control%do_n2o_iz) then
        if (Lw_control%do_n2o) then
          if (do_writestdn2otfs .and. do_readstdn2otfs) then
            call error_mesg ( 'gas_tf_mod', &
              ' cannot read and write std tfs in same job', FATAL)
          endif
          if (do_writestdn2otfs .and. .not. do_calcstdn2otfs) then
            call error_mesg ( 'gas_tf_mod', &
               ' cannot write std tfs without calculating them', FATAL)
          endif
        endif
      else
        call error_mesg ('gas_tf_mod', &
              'do_n2o has not yet been defined', FATAL)
      endif

!--------------------------------------------------------------------
!    call ptz to compute standard temps and a pressure coefficient 
!    (gtemp) used in the radiation algorithm. 
!--------------------------------------------------------------------
      call ptz (plm, pd)

!--------------------------------------------------------------------
!    convert pressure specification for top (flux) pressure level
!    for nearby layer calculation into an index (ixprkminh2o)
!    note: minimum value of ixprkminh2o is KSRAD . (but if KSRAD = 1,
!    plm(1) is zero, so minimum value of KSRAD is at least 2).
!    if all levels used for radiative calculations are at pressures
!    less than 28 mb, nearby layer effects are going to be ignored,
!    so ixprkminh2o is set to KERAD+1 to avoid loop calculations.
!--------------------------------------------------------------------
      if (plm(ks) >= prkminh2o) then
        ixprkminh2o = 1
      else if (plm(kmax ) < prkminh2o) then
        ixprkminh2o = (kmax - ks + 1) + 1
      else
        do k=ks+1,kmax  
          if ((plm(k) - prkminh2o) .LT. 0.0) then
!! ixprkminh2o in radiation grid coordianates, not model grid coords
          else
            ixprkminh2o = k-ks + 1
            exit
          endif
        enddo
      endif
   
!----------------------------------------------------------------------
!    allocate co2 transmission function arrays to hold data which will 
!    either be read in or will be coming from lw_gases_stdtf module.
!----------------------------------------------------------------------
     if (Lw_control%do_co2) then
      allocate (cdtm51(KSRAD:KERAD) , &
                co2m51(KSRAD:KERAD) , &
                c2dm51(KSRAD:KERAD) , &
                cdtm58(KSRAD:KERAD) , &
                co2m58(KSRAD:KERAD) , &
                c2dm58(KSRAD:KERAD)   )
      allocate (co2dt15nbps1(KSRAD:KERAD+1, NBCO215) , &
                co215nbps1(KSRAD:KERAD+1, NBCO215) , &
                co2d2t15nbps1(KSRAD:KERAD+1, NBCO215) , &
                co2dt15nbps8(KSRAD:KERAD+1, NBCO215) , &
                co215nbps8(KSRAD:KERAD+1, NBCO215) , &
                co2d2t15nbps8(KSRAD:KERAD+1, NBCO215)   )
      allocate (cdt51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                co251(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                c2d51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                cdt58(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                co258(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                c2d58(KSRAD:KERAD+1, KSRAD:KERAD+1)   )
      allocate (co211(KSRAD:KERAD+1) , &
                co218(KSRAD:KERAD+1) )
     endif

!----------------------------------------------------------------------
!    allocate ch4 and n2o transmission function arrays to hold data 
!    which will either be read in or will be coming from 
!    lw_gases_stdtf module.
!----------------------------------------------------------------------
      if (Lw_control%do_ch4) then
        allocate (ch4dt51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  ch451(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  ch4d2t51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  ch4dt58(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  ch458(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  ch4d2t58(KSRAD:KERAD+1, KSRAD:KERAD+1)   )
      endif
      if (Lw_control%do_n2o) then
        allocate (n2odt51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t51(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2odt58(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o58(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t58(KSRAD:KERAD+1, KSRAD:KERAD+1)   )
        allocate (n2odt91(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o91(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t91(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2odt98(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o98(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t98(KSRAD:KERAD+1, KSRAD:KERAD+1)   )
        allocate (n2odt71(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o71(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t71(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2odt78(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2o78(KSRAD:KERAD+1, KSRAD:KERAD+1) , &
                  n2od2t78(KSRAD:KERAD+1, KSRAD:KERAD+1)   )
      endif

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.
 
!---------------------------------------------------------------------


end subroutine gas_tf_init



!###################################################################
! <SUBROUTINE NAME="co2coef">
!  <OVERVIEW>
!   Calculate CO2 absorption coefficients and transmission function
!  </OVERVIEW>
!  <DESCRIPTION>
!   Calculate CO2 absorption coefficients and transmission function
!  </DESCRIPTION>
!  <TEMPLATE>
!   call co2coef(Atmos_input, Gas_tf)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   The input data of the atmosphere structure and gas concentration
!  </IN>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   The gas transmission function table
!  </INOUT>
! </SUBROUTINE>
!
subroutine co2coef (Atmos_input, Gas_tf)
 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

type(gas_tf_type),      intent(inout) :: Gas_tf
type(atmos_input_type), intent(in)    :: Atmos_input

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     Atmos_input
!
!  intent(inout) variables:
!
!     Gas_tf
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(Atmos_input%pflux,1), &
                       size(Atmos_input%pflux,2), &
                       size(Atmos_input%pflux,3)-1) :: pdflux

      real, dimension (size(Atmos_input%pflux,1), &
                       size(Atmos_input%pflux,2), &
                       size(Atmos_input%pflux,3)  ) ::  &
                             tdif, press, temp, pflux, tflux

      real                                 ::  palog8, alogps8 
      integer                              ::  i, j, k   
      integer                              :: israd, ierad, jsrad, jerad

!---------------------------------------------------------------------
!  local variables:
!
!      pdflux
!      tdif
!      press
!      temp
!      pflux
!      tflux
!      palog8
!      alogp8
!      i,j,k
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    convert press and pflux to cgs.
!---------------------------------------------------------------------
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
      pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)
      tflux(:,:,:) = Atmos_input%tflux(:,:,:)
      temp(:,:,:) = Atmos_input%temp(:,:,:)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      israd = 1
      ierad = size(Atmos_input%press,1)
      jsrad = 1
      jerad = size(Atmos_input%press,2)
    
!---------------------------------------------------------------------
!    allocate module variables
!---------------------------------------------------------------------
      allocate (Gas_tf%a1      (ISRAD:IERAD, JSRAD:JERAD           ))
      allocate (Gas_tf%a2      (ISRAD:IERAD, JSRAD:JERAD           ))
      allocate (Gas_tf%tdav    (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%tlsqu   (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%tmpdiff (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%tstdav  (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%co2nbl  (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD  ))
      allocate (Gas_tf%n2o9c   (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%tn2o17  (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
      allocate (Gas_tf%co2spnb (ISRAD:IERAD, JSRAD:JERAD,  &
                                               KSRAD:KERAD+1,  NBCO215))
      Gas_tf%co2nbl  = 1.0                                         
      Gas_tf%co2spnb = 1.0                                           
      Gas_tf%n2o9c  = 0.                                          
      Gas_tf%tn2o17  = 0.0                                        

!--------------------------------------------------------------------
!    compute temperature difference between model profile and 
!    standard profile
!---------------------------------------------------------------------
      do k = KSRAD,KERAD+1
        Gas_tf%tmpdiff(:,:,k) = temp(:,:,k) - stemp(k)
      enddo

!-----------------------------------------------------------------------
!    compute weighted temperature difference (tdav) and pressure
!    integrals (tstdav) from level KSRAD to level KERAD. the quotient 
!    of these will be used to obtain  the difference (dift) between the
!    model temperature profile and the standard profile.
!-----------------------------------------------------------------------
      Gas_tf%tstdav(:,:,KSRAD) = 0.0E+00
      Gas_tf%tdav  (:,:,KSRAD) = 0.0E+00
      do k=KSRAD,KERAD
        pdflux(:,:,k) = pflux(:,:,k+1) - pflux(:,:,k)
        Gas_tf%tstdav(:,:,k+1) = Gas_tf%tstdav(:,:,k) + gtemp(k)* &
                                 pdflux(:,:,k)
        Gas_tf%tdav  (:,:,k+1) = Gas_tf%tdav  (:,:,k) + gtemp(k)* &
                                 pdflux(:,:,k)*Gas_tf%tmpdiff(:,:,k)
      enddo

!----------------------------------------------------------------------
!    evaluate coefficients for co2 pressure interpolation (a1, a2).
!    a linear interpolation is presently assumed, with the 2 
!    2nd pressure profile having pressures 0.8* the first, thus
!    accounting for the 0.8 and 0.2 factors.
!----------------------------------------------------------------------
      if (do_linearlblint) then
        Gas_tf%a1(:,:) = (press(:,:,KERAD+1) - pstd*0.8E+00)/ &
                         (pstd*0.2E+00)
        Gas_tf%a2(:,:) = (pstd - press(:,:,KERAD+1))/(pstd*0.2E+00)

!----------------------------------------------------------------------
!    a logarithmic interpolation is presently assumed, with the 2 
!    2nd pressure profile having pressures 0.8* the first, thus
!    accounting for the 0.8 and 0.2 factors. The denominator, which
!    is (log(pstd) - log(0.8*pstd)) is a constant (-log(0.8)) so the
!    expression can be replaced by the quantity palog8.
!----------------------------------------------------------------------
      else if (do_loglblint) then 
        alogps8 = ALOG(pstd*0.8E+00)
        palog8 = -ALOG(0.8E+00)
        Gas_tf%a1(:,:) = (ALOG(press(:,:,KERAD+1)) - alogps8)/palog8
        Gas_tf%a2(:,:) = 1.0E+00 - Gas_tf%a1(:,:)
      else 
        call error_mesg ('gas_tf_mod', &
              'neither linearlblint nor loglblint was specified.', &
                                                                FATAL)
      endif

!----------------------------------------------------------------------
!    compute temperature coefficient based on tflux. see fels and
!    schwarzkopf (1981) for details.
!----------------------------------------------------------------------
      tdif(:,:,:) = tflux(:,:,:) - 2.5E+02
      do k=KSRAD,KERAD+1
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            if (tflux(i,j,k) .LE. 2.5E+02) then
              Gas_tf%tlsqu(i,j,k) = b0 +     tdif (i,j,k)  *  &
                                   (b1 +     tdif (i,j,k)  *  &
                                   (b2 + b3* tdif (i,j,k)  )) 
            else 
              Gas_tf%tlsqu(i,j,k) = b0 
            endif
          enddo
        enddo
      enddo


!----------------------------------------------------------------------
!    call transfn to compute temperature-corrected co2 transmission 
!    functions (co2spnb and co2nbl). 
!---------------------------------------------------------------------
    if (Lw_control%do_co2) then
      call transfn (Gas_tf)
    endif

!-------------------------------------------------------------------

end subroutine co2coef





!#####################################################################
! <SUBROUTINE NAME="transcol">
!  <OVERVIEW>
!   Subroutine to compute temperature-corrected co2 transmission 
!   functions at a particular (krow).
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute temperature-corrected co2 transmission 
!   functions at a particular (krow).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call transcol (kcol, krow, kcols, kcole, co21c, Gas_tf)
!  </TEMPLATE>
!  <IN NAME="kcol" TYPE="integer">
!   Not used
!  </IN>
!  <IN NAME="krow" TYPE="integer">
!   The row index where co2 transmission is calculated
!  </IN>
!  <IN NAME="kcols" TYPE="integer">
!   The starting column index number
!  </IN>
!  <IN NAME="kcole" TYPE="integer">
!   The ending column index number
!  </IN>
!  <OUT NAME="co21c" TYPE="real">
!   The column of transmission functions
!  </OUT>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   The pre temperature-corrected co2 transmission functions
!  </IN>
! </SUBROUTINE>
!
subroutine transcol (kcol, krow, kcols, kcole, co21c, Gas_tf)        

!---------------------------------------------------------------------
!    transcol computes temperature-corrected co2 transmission 
!    functions at a particular (krow).
!    author: c. l. kerr
!    revised: 11/11/93
!    certified:  radiation version 1.0
!---------------------------------------------------------------------

integer,                intent(in)  :: krow, kcol, kcols, kcole
real, dimension(:,:,:), intent(out) :: co21c
type(gas_tf_type),      intent(in)  :: Gas_tf

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      krow
!      kcol
!      kcols
!      kcole
!      Gas_tf
!
!  intent(out) variables:
!
!      co21c    column of transmission functions.
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(Gas_tf%tdav,1),&
                       size(Gas_tf%tdav,2), &
                       size(Gas_tf%tdav,3)  ) ::   &  
                                            co2r, dift,  d2cdt2, dco2dt
      integer    ::   kp

!---------------------------------------------------------------------
!  local variables:
!
!    co2r
!    dift
!    d2cdt2
!    dco2dt
!    k,kp
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      co21c(:,:,KSRAD:KERAD+1) = 1.0E+00

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do kp = kcols,kcole
        if (kp .NE. krow) then
          dift(:,:,kp) = (Gas_tf%tdav  (:,:,kp) -   &
                          Gas_tf%tdav(:,:,krow))/  &
                         (Gas_tf%tstdav(:,:,kp) -  &
                          Gas_tf%tstdav(:,:,krow))
        else if (krow .NE. KSRAD) then
          dift(:,:,kp) = 0.5E+00*(Gas_tf%tmpdiff(:,:,kp) +  &
                                  Gas_tf%tmpdiff(:,:,kp-1))
        else
          dift(:,:,kp) = 0.0E+00
        endif
      end do

!----------------------------------------------------------------------
!    obtain transmission functions used for the flux at a fixed level
!    (krow). ie, tf's  from varying flux levels (kp) to (krow)
!       pressure interpolation
!----------------------------------------------------------------------
      do kp=kcols,kcole
        co2r  (:,:,kp) = Gas_tf%a1(:,:)*co251(kp,krow) + &
                         Gas_tf%a2(:,:)*co258(kp,krow)
        dco2dt(:,:,kp) = 1.0E-02*(Gas_tf%a1(:,:)*cdt51(kp,krow) +   &
                                  Gas_tf%a2(:,:)*cdt58(kp,krow))
        d2cdt2(:,:,kp) = 1.0E-03*(Gas_tf%a1(:,:)*c2d51(kp,krow) +  &
                                  Gas_tf%a2(:,:)*c2d58(kp,krow))
      enddo
 
!----------------------------------------------------------------------
!    temperature interpolation
!----------------------------------------------------------------------
      do kp=kcols,kcole
        co21c (:,:,kp) = co2r(:,:,kp) + dift(:,:,kp)*(dco2dt(:,:,kp) + &
                         0.5E+00*dift(:,:,kp)*d2cdt2(:,:,kp))
      enddo
 
!----------------------------------------------------------------------
!    correction for finite width of co2 bands
!    (Eqs. 7a-7c, Ref. (2))
!----------------------------------------------------------------------
      do kp=kcols,kcole
        co21c(:,:,kp) = co21c(:,:,kp)*(1.0E+00 -  &
                        Gas_tf%tlsqu(:,:,kp)) +  &
                        Gas_tf%tlsqu(:,:,kp)
      enddo

!----------------------------------------------------------------


end subroutine transcol



 
!#####################################################################
! <SUBROUTINE NAME="transcolrow">
!  <OVERVIEW>
!   Subroutine to compute temperature-corrected co2/ch4/n2o transmission 
!   functions at a particular row and particular column.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute temperature-corrected co2/ch4/n2o transmission 
!   functions at a particular row and particular column.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call transcolrow (Gas_tf, kcol, krow, kcols, kcole, krows, krowe,    &
!                     co21c, co21r, tch4n2oe)   
!  </TEMPLATE>
!  <IN NAME="kcol" TYPE="integer">
!   The column index of temperature-corrected transmission function
!  </IN>
!  <IN NAME="krow" TYPE="integer">
!   The row index of temperature-corrected transmission function
!  </IN>
!  <IN NAME="kcols" TYPE="integer">
!   The starting column index number
!  </IN>
!  <IN NAME="kcole" TYPE="integer">
!   The ending column index number
!  </IN>
!  <IN NAME="krows" TYPE="integer">
!   The starting row index number
!  </IN>
!  <IN NAME="krowe" TYPE="integer">
!   The ending row index number
!  </IN>
!  <OUT NAME="co21c" TYPE="real">
!   The column of transmission functions
!  </OUT>
!  <OUT NAME="co21r" TYPE="real">
!   The row of transmission functions
!  </OUT>
!  <INOUT NAME="tch4n2oe" TYPE="real">
!   The ch4 and n2o transmission functions
!  </INOUT>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   The pre temperature-corrected co2 transmission functions
!  </INOUT>
! </SUBROUTINE>
!
subroutine transcolrow (Gas_tf, kcol, krow, kcols, kcole, krows, krowe,&
                        co21c, co21r, tch4n2oe)                     

!----------------------------------------------------------------------
!    transcolrow computes the temperature-corrected co2 transmission 
!    functions for a particular (krow) (varying column index) and for 
!    a particular (kcol) (varying row index).
!    transcolrow also computes the pressure-interpolated ch4 and n2o
!    transmission functions for a particular (krow) and for a parti-
!    cular (kcol). By assumption, no correction for finite bandwidth
!    is performed.
!    author: c. l. kerr
!    revised: 11/11/93
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

integer,                   intent(in)    ::  kcol, krow, kcols, kcole, &
                                             krows, krowe
type(gas_tf_type),         intent(inout) ::  Gas_tf
real, dimension (:,:,:),   intent(out)   :: co21c, co21r
real, dimension (:,:,:,:), intent(inout) :: tch4n2oe

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     kcol
!     krow
!     kcols
!     kcole
!     krows
!     krowe
!
!  intent(inout) variables:
!
!     Gas_tf
!     tch4n2oe
!
!  intent(out) variables:
!     
!     co21c    column of transmission functions (fixed krow).
!     co21r    column of transmission functions (fixed kcol).
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (size(Gas_tf%tdav,1),  &
                       size(Gas_tf%tdav,2), &
                       size(Gas_tf%tdav,3)  ) ::   &  
                                           ch41c, n2o1c, n2o17c,&
                                           co2p, dift, d2cdt2, dco2dt,&
                                           ch4p, d2ch4dt2, dch4dt, &
                                           d2n2odt2, dn2odt,    &
                                           d2n2o17dt2, dn2o17dt, &
                                           d2n2o9dt2, dn2o9dt, n2op,  &
                                           n2o17p, n2o9p

      integer    :: kp

!--------------------------------------------------------------------
!   local variables:
!
!      ch41c  = column of ch4 transmission functions (fixed krow).
!      ch41r  = column of ch4 transmission functions (fixed kcol).
!      n2o1c  = column of n2o transmission functions (fixed krow).
!      n2o1r  = column of n2o transmission functions (fixed kcol).
!      n2o17c = column of n2o 17 um transmission functions (fixed krow).
!      n2o17r = column of n2o 17 um transmission functions (fixed kcol).
!      n2o9c = column of n2o 9 um transmission functions (fixed krow).
!      n2o9r = column of n2o 9 um transmission functions (fixed kcol).
!      kp
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------------
!    initialization.
!-----------------------------------------------------------------------
      co21c(:,:,kcols:kcole  ) = 1.0E+00
      co21r(:,:,KSRAD:KERAD+1) = 1.0E+00
      ch41c(:,:,kcols:kcole  ) = 1.0E+00
      n2o1c(:,:,kcols:kcole) = 1.0E+00
      n2o17c(:,:,KSRAD:KERAD+1) = 1.0E+00
      Gas_tf%n2o9c(:,:,KSRAD:KERAD+1) = 1.0E+00

!-----------------------------------------------------------------------
!    temperature difference averaged between levels k and kp
!-----------------------------------------------------------------------
      do kp = kcols,kcole
        if (kp .NE. krow) then
          dift(:,:,kp) = (Gas_tf%tdav  (:,:,kp) -  &
                          Gas_tf%tdav  (:,:,krow))/  &
                         (Gas_tf%tstdav(:,:,kp) -  &
                          Gas_tf%tstdav(:,:,krow))
        elseif (krow .NE. KSRAD) then
          dift(:,:,kp) = 0.5E+00*(Gas_tf%tmpdiff(:,:,kp) + &
                                  Gas_tf%tmpdiff(:,:,kp-1))
        else
          dift(:,:,kp) = 0.0E+00
        endif
      end do

!-----------------------------------------------------------------------
!    obtain transmission functions used for the flux at a fixed level
!    (krow). ie, tf's  from varying flux levels (kp) to (krow)
!       pressure interpolation
!-----------------------------------------------------------------------
      do kp=kcols,kcole
        if (Lw_control%do_co2) then
        co2p  (:,:,kp) = Gas_tf%a1(:,:)*co251(kp,krow) +  &
                         Gas_tf%a2(:,:)*co258(kp,krow)
        dco2dt(:,:,kp) = 1.0E-02*(Gas_tf%a1(:,:)*cdt51(kp,krow) +   &
                                  Gas_tf%a2(:,:)*cdt58(kp,krow))
        d2cdt2(:,:,kp) = 1.0E-03*(Gas_tf%a1(:,:)*c2d51(kp,krow) +   &
                                  Gas_tf%a2(:,:)*c2d58(kp,krow))
        endif
        if (Lw_control%do_ch4) then
          if (Lw_control%do_ch4lbltmpint) then
            ch4p (:,:,kp)  = Gas_tf%a1(:,:)*ch451(kp,krow) +  &
                             Gas_tf%a2(:,:)*ch458(kp,krow)
            dch4dt(:,:,kp) = 1.0E-02*(Gas_tf%a1(:,:)*ch4dt51(kp,krow) +&
                                      Gas_tf%a2(:,:)*ch4dt58(kp,krow))
            d2ch4dt2(:,:,kp) = 1.0E-03*  &
                                  (Gas_tf%a1(:,:)*ch4d2t51(kp,krow) +  &
                                   Gas_tf%a2(:,:)*ch4d2t58(kp,krow))
          else
            ch41c(:,:,kp)  = Gas_tf%a1(:,:)*ch451(kp,krow) +  &
                             Gas_tf%a2(:,:)*ch458(kp,krow)
          endif
        endif
        if (Lw_control%do_n2o) then
          if (Lw_control%do_n2olbltmpint) then
            n2op (:,:,kp)  = Gas_tf%a1(:,:)*n2o51(kp,krow) + &
                             Gas_tf%a2(:,:)*n2o58(kp,krow)
            n2o17p(:,:,kp) = Gas_tf%a1(:,:)*n2o71(kp,krow) + &
                             Gas_tf%a2(:,:)*n2o78(kp,krow)
            n2o9p(:,:,kp) = Gas_tf%a1(:,:)*n2o91(kp,krow) + &
                            Gas_tf%a2(:,:)*n2o98(kp,krow)
            dn2odt(:,:,kp) = 1.0E-02*(Gas_tf%a1(:,:)*n2odt51(kp,krow) +&
                                    Gas_tf%a2(:,:)*n2odt58(kp,krow))
            dn2o17dt(:,:,kp) = 1.0E-02* &
                                    (Gas_tf%a1(:,:)*n2odt71(kp,krow) +&
                                     Gas_tf%a2(:,:)*n2odt78(kp,krow))
            dn2o9dt(:,:,kp) = 1.0E-02*  &
                                    (Gas_tf%a1(:,:)*n2odt91(kp,krow) + &
                                     Gas_tf%a2(:,:)*n2odt98(kp,krow))
            d2n2odt2(:,:,kp) = 1.0E-03*  &
                                  (Gas_tf%a1(:,:)*n2od2t51(kp,krow) +  &
                                   Gas_tf%a2(:,:)*n2od2t58(kp,krow))
            d2n2o17dt2(:,:,kp) = 1.0E-03* &
                                  (Gas_tf%a1(:,:)*n2od2t71(kp,krow) +  &
                                   Gas_tf%a2(:,:)*n2od2t78(kp,krow))
            d2n2o9dt2(:,:,kp) = 1.0E-03*  &
                                  (Gas_tf%a1(:,:)*n2od2t91(kp,krow) +  &
                                   Gas_tf%a2(:,:)*n2od2t98(kp,krow))
          else
            n2o1c(:,:,kp)  = Gas_tf%a1(:,:)*n2o51(kp,krow) + &
                             Gas_tf%a2(:,:)*n2o58(kp,krow)
            n2o17c(:,:,kp) = Gas_tf%a1(:,:)*n2o71(kp,krow) + &
                             Gas_tf%a2(:,:)*n2o78(kp,krow)
            Gas_tf%n2o9c(:,:,kp) = Gas_tf%a1(:,:)*n2o91(kp,krow) +  &
                                   Gas_tf%a2(:,:)*n2o98(kp,krow)
          endif
        endif
      enddo

!----------------------------------------------------------------------
!    temperature interpolation
!----------------------------------------------------------------------
      do kp=kcols,kcole
        if (Lw_control%do_co2) then
        co21c (:,:,kp) = co2p(:,:,kp) + dift(:,:,kp)*(dco2dt(:,:,kp) + &
                         0.5E+00*dift(:,:,kp)*d2cdt2(:,:,kp))
        endif
        if (Lw_control%do_ch4lbltmpint) then
          ch41c (:,:,kp) = ch4p(:,:,kp) + dift(:,:,kp)* &
                           (dch4dt(:,:,kp) + 0.5E+00*dift(:,:,kp)* &
                            d2ch4dt2(:,:,kp))
        endif
        if (Lw_control%do_n2olbltmpint) then
          n2o1c (:,:,kp) = n2op(:,:,kp) + dift(:,:,kp)* &
                           (dn2odt(:,:,kp) + 0.5E+00*dift(:,:,kp)* &
                            d2n2odt2(:,:,kp))
          n2o17c(:,:,kp) = n2o17p(:,:,kp) +  &
                           dift(:,:,kp)*(dn2o17dt(:,:,kp) +  &
                           0.5E+00*dift(:,:,kp)*d2n2o17dt2(:,:,kp))
          Gas_tf%n2o9c(:,:,kp) = n2o9p(:,:,kp) +   &
                                 dift(:,:,kp)*(dn2o9dt(:,:,kp) +  &
                                 0.5E+00*dift(:,:,kp)*d2n2o9dt2(:,:,kp))
        endif
      enddo

!---------------------------------------------------------------------
!    correction for finite width of co2 bands
!    (Eqs. 7a-7c, Ref. (2))
!---------------------------------------------------------------------
      do kp=kcols,kcole
        if (Lw_control%do_co2) then
        co21c(:,:,kp) = co21c(:,:,kp)*(1.0E+00 -  &
                        Gas_tf%tlsqu(:,:,kp)) + Gas_tf%tlsqu(:,:,kp)
        endif
      enddo

!-----------------------------------------------------------------------
!    obtain transmission functions used for the flux for varying levels
!    (krow) from a fixed level (kcol). ie, tf's  from a fixed flux
!    level (kcol) to varying levels (krow).
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!    temperature difference averaged between levels k and kp. This 
!    computation is made unless krow = kcol, and range (krows,krowe) is
!    entirely within (kcols,kcole), in which case the dift computed
!    for column tfs is applicable to row tfs.
!-----------------------------------------------------------------------
      if (Lw_control%do_co2) then
        if (kcol  .NE. krow   .or. krows .LT. kcols  .or.  &
            krowe .GT. kcole)     then
          do kp = krows,krowe
            if (kp .NE. krow) then
              dift(:,:,kp) = (Gas_tf%tdav(:,:,kp) -  &
                              Gas_tf%tdav(:,:,krow))/  &
                             (Gas_tf%tstdav(:,:,kp) -  &
                              Gas_tf%tstdav(:,:,krow))
            elseif (krow .NE. KSRAD) then
              dift(:,:,kp) = 0.5E+00*(Gas_tf%tmpdiff(:,:,kp) +  &
                                      Gas_tf%tmpdiff(:,:,kp-1))
            else
              dift(:,:,kp) = 0.0E+00
            endif
          end do
        endif
      endif  ! (do_co2)

!--------------------------------------------------------------------
!    pressure interpolation
!--------------------------------------------------------------------
      do kp=krows,krowe
        if (Lw_control%do_co2) then
        co2p  (:,:,kp) = Gas_tf%a1(:,:)*co251(kcol,kp) +  &
                         Gas_tf%a2(:,:)*co258(kcol,kp)
        dco2dt(:,:,kp) = 1.0E-02*(Gas_tf%a1(:,:)*cdt51(kcol,kp) +   &
                                  Gas_tf%a2(:,:)*cdt58(kcol,kp))
        d2cdt2(:,:,kp) = 1.0E-03*(Gas_tf%a1(:,:)*c2d51(kcol,kp) +   &
                                  Gas_tf%a2(:,:)*c2d58(kcol,kp))
        endif
        if (Lw_control%do_ch4) then
          if (Lw_control%do_ch4lbltmpint) then
            ch4p (:,:,kp)  = Gas_tf%a1(:,:)*ch451(kcol,kp) +  &
                             Gas_tf%a2(:,:)*ch458(kcol,kp)
            dch4dt(:,:,kp) = 1.0E-02* &
                             (Gas_tf%a1(:,:)*ch4dt51(kcol,kp) +   &
                              Gas_tf%a2(:,:)*ch4dt58(kcol,kp))
            d2ch4dt2(:,:,kp) = 1.0E-03*  &
                               (Gas_tf%a1(:,:)*ch4d2t51(kcol,kp) +  &
                                Gas_tf%a2(:,:)*ch4d2t58(kcol,kp))
          endif
        endif
        if (Lw_control%do_n2o) then
          if (Lw_control%do_n2olbltmpint) then
            n2op (:,:,kp)  = Gas_tf%a1(:,:)*n2o51(kcol,kp) +&
                             Gas_tf%a2(:,:)*n2o58(kcol,kp)
            n2o17p(:,:,kp) = Gas_tf%a1(:,:)*n2o71(kcol,kp) + &
                             Gas_tf%a2(:,:)*n2o78(kcol,kp)
            n2o9p(:,:,kp) = Gas_tf%a1(:,:)*n2o91(kcol,kp) + &
                            Gas_tf%a2(:,:)*n2o98(kcol,kp)
            dn2odt(:,:,kp) = 1.0E-02* &
                             (Gas_tf%a1(:,:)*n2odt51(kcol,kp) +  &
                              Gas_tf%a2(:,:)*n2odt58(kcol,kp))
            dn2o17dt(:,:,kp) = 1.0E-02* &
                               (Gas_tf%a1(:,:)*n2odt71(kcol,kp) +   &
                                Gas_tf%a2(:,:)*n2odt78(kcol,kp))
            dn2o9dt(:,:,kp) = 1.0E-02* &
                              (Gas_tf%a1(:,:)*n2odt91(kcol,kp) +   &
                               Gas_tf%a2(:,:)*n2odt98(kcol,kp))
            d2n2odt2(:,:,kp) = 1.0E-03* &
                               (Gas_tf%a1(:,:)*n2od2t51(kcol,kp) +   &
                                Gas_tf%a2(:,:)*n2od2t58(kcol,kp))
            d2n2o17dt2(:,:,kp) = 1.0E-03* &
                                 (Gas_tf%a1(:,:)*n2od2t71(kcol,kp) +   &
                                  Gas_tf%a2(:,:)*n2od2t78(kcol,kp))
            d2n2o9dt2(:,:,kp) = 1.0E-03* &
                                (Gas_tf%a1(:,:)*n2od2t91(kcol,kp) +   &
                                 Gas_tf%a2(:,:)*n2od2t98(kcol,kp))
          endif
        endif
      enddo

!---------------------------------------------------------------------
!    temperature interpolation
!---------------------------------------------------------------------
      do kp=krows,krowe
        if (Lw_control%do_co2) then
        co21r (:,:,kp) = co2p(:,:,kp) + dift(:,:,kp)*(dco2dt(:,:,kp) +&
                         0.5E+00*dift(:,:,kp)*d2cdt2(:,:,kp))
        endif
      enddo

!---------------------------------------------------------------------
!    correction for finite width of co2 bands
!    (Eqs. 7a-7c, Ref. (2))
!---------------------------------------------------------------------
      do kp=krows,krowe
        if (Lw_control%do_co2) then
        co21r(:,:,kp) = co21r(:,:,kp)*(1.0E+00 -  &
                        Gas_tf%tlsqu(:,:,kcol)) +  &
                        Gas_tf%tlsqu(:,:,kcol)
        endif
      enddo

!----------------------------------------------------------------------
!    save the values which are needed elsewhere
!    tn2o17 results are for 560-630 cm-1 band. (if NBCO215=3)
!----------------------------------------------------------------------
      if (Lw_control%do_ch4 .or. Lw_control%do_n2o) then
        if (kcols == 1) then
          tch4n2oe(:,:,kcols,1) = ch41c(:,:,kcols)*n2o1c(:,:,kcols)
        endif
        tch4n2oe(:,:,kcols+1:kcole,1) = ch41c(:,:,kcols+1:kcole)*  &
                                        n2o1c(:,:,kcols+1:kcole)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      Gas_tf%tn2o17(:,:,:) = n2o17c(:,:,:)

!---------------------------------------------------------------------



end subroutine transcolrow




!#####################################################################
! <SUBROUTINE NAME="trans_nearby">
!  <OVERVIEW>
!   Compute nearby layer transmission functions at certain level in
!   the frequency band at 15 um
!  </OVERVIEW>
!  <DESCRIPTION>
!   Compute "nearby  layer" transmission functions at level k 
!  ( tau(p(k),p(k))) in the frequency band at 15 um. include all
!  gases (co2, h2o, h2o cont) used in computing fluxes in this band.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call trans_nearby (Gas_tf, Atmos_input, overod, co2diag)
!  </TEMPLATE>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   The gas transmission functions at model coordinate system
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   The atmospheric input data
!  </IN>
!  <IN NAME="overod" TYPE="real">
!   CO2 data
!  </IN>
!  <OUT NAME="co21diag" TYPE="real">
!   CO2 transmission function
!  </OUT>
! </SUBROUTINE>
subroutine trans_nearby (Gas_tf, Atmos_input, overod, co21diag)
 
!-------------------------------------------------------------------
!    compute "nearby  layer" transmission functions at level k 
!    ( tau(p(k),p(k))) in the frequency band at 15 um. include all
!    gases (co2, h2o, h2o cont) used in computing fluxes in this band.
!    the algorithm assumes that at pressures (p') near the pressure
!    at level k (p(k)), the transmission function may be written as:
!
!              tau(p',p(k)) = EXP(-alpha*SQRT(p'-p(k)))
!
!    with alpha determined by the boundary condition(s)
!    tau(p(k+1),p(k)) and tau(p(k-1),p(k)) = the values from  "normal"
!    calculations. An integration is performed over the "nearby" 
!    pressure layer to obtain tau(p(k),p(k)).
!    the computation is not done for levels from KSRAD to KMINH2O-1 (if
!    different), where it is assumed that the h2o transmissivities 
!    are near unity, and that the precomputed co2 transmissivities may
!    be used.
!     two "special case" transmissivities, viz.,
!    tau(p(KERAD),p(KERAD+1)) and tau(p(KERAD+1),p(KERAD)) are also 
!    evaluated using the above assumptions and an integration.
!-------------------------------------------------------------------

type(gas_tf_type),       intent(in)  :: Gas_tf
type(atmos_input_type),  intent(in)  :: Atmos_input
real, dimension (:,:,:), intent(in)  :: overod  
real, dimension (:,:,:), intent(out) :: co21diag

!--------------------------------------------------------------------
! intent(in) variables:
!
!    Gas_tf
!    Atmos_input
!    overod
!
! intent(out) variables:
!
!    co21diag
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:




      real, dimension (size(Atmos_input%pflux,1),  &
                       size(Atmos_input%pflux,2),  &
                       size(Atmos_input%pflux,3)-1) :: pdfinv

      real, dimension (size(Atmos_input%pflux,1), &
                       size(Atmos_input%pflux,2),  &
                       size(Atmos_input%pflux,3)) ::         &
               press, pflux, alpa, alpb, ca, cb, delpr1, delpr2, rlog

      integer   :: k, km, kmp1

!-------------------------------------------------------------------
!  local variables:
!
!    pdfl
!    pdfinv
!    press
!
!    k
!
!    kmp1
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    convert press and pflux to cgs.
!---------------------------------------------------------------------
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
      pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      pdfinv(:,:,KSRAD:KERAD) = 1.0/(pflux(:,:,KSRAD+1:KERAD+1) -   &
                                     pflux(:,:,KSRAD:KERAD) )

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      km= MAX(ixprkminh2o - 1, KSRAD)
      kmp1 = MAX(ixprkminh2o - 1, KSRAD+1)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      rlog  (:,:,km:KERAD)     = LOG(Gas_tf%co2nbl(:,:,km:KERAD)*   &
                                 overod(:,:,km+1:KERAD+1))
      delpr1(:,:,kmp1:KERAD)   = pdfinv (:,:,kmp1:KERAD)*  &
                                 (press(:,:,kmp1:KERAD) - &
                                  pflux(:,:,kmp1:KERAD)) 
      alpb  (:,:,kmp1:KERAD)   = -SQRT(delpr1(:,:,kmp1:KERAD))*  &
                                  rlog(:,:,kmp1:KERAD)
      delpr2(:,:,kmp1:KERAD+1) = pdfinv(:,:,kmp1-1:KERAD)*  &
                                 (pflux(:,:,kmp1:KERAD+1) -  &
                                  press(:,:,kmp1-1:KERAD)) 
      alpa  (:,:,km:KERAD)     = -SQRT(delpr2(:,:,km+1:KERAD+1))*  &
                                 rlog(:,:,km:KERAD)
      alpa  (:,:,KERAD+1)      = -rlog(:,:,KERAD)
      alpb  (:,:,KERAD+1)      = -rlog(:,:,KERAD)* &
                                  SQRT(pdfinv(:,:,KERAD)*&
                                 (pflux(:,:,KERAD+1) -   &
                                  press(:,:,KERAD-1)))
      ca(:,:,km:KERAD+1) = alpa(:,:,km:KERAD+1)*(-0.66667E+00 +  &
                           alpa(:,:,km:KERAD+1)*(0.25E+00 +   &
                           alpa(:,:,km:KERAD+1)*(-0.066667E+00)))
      cb(:,:,kmp1:KERAD+1) = alpb(:,:,kmp1:KERAD+1)*(-0.66667E+00 +  &
                             alpb(:,:,kmp1:KERAD+1)*(0.25E+00 +    &
                             alpb(:,:,kmp1:KERAD+1)*(-0.066667E+00)))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=ixprkminh2o,KERAD
        co21diag(:,:,k) = 1.0E+00 + 0.5E+00*(cb(:,:,k) + ca(:,:,k-1))
      enddo
      co21diag(:,:,KERAD+1) = 1.0E+00 + ca(:,:,KERAD)

!-------------------------------------------------------------------




end subroutine trans_nearby

!#####################################################################
! <SUBROUTINE NAME="trans_sfc">
!  <OVERVIEW>
!   Compute nearby layer transmission functions at certain level in
!   the frequency band at 15 um
!  </OVERVIEW>
!  <DESCRIPTION>
!   Compute "nearby  layer" transmission functions at level k 
!  ( tau(p(k),p(k))) in the frequency band at 15 um. include all
!  gases (co2, h2o, h2o cont) used in computing fluxes in this band.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call trans_sfc    (Gas_tf, Atmos_input, overod, co21c, co21r)
!  </TEMPLATE>
!  <IN NAME="Gas_tf" TYPE="gas_tf_type">
!   The gas transmission functions at model coordinate system
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   The atmospheric input data
!  </IN>
!  <IN NAME="overod" TYPE="real">
!   CO2 data
!  </IN>
!  <OUT NAME="co21c" TYPE="real">
!   CO2 transmission function
!  </OUT>
!  <OUT NAME="co21r" TYPE="real">
!   CO2 transmission function
!  </OUT>
! </SUBROUTINE>
!subroutine trans_sfc    (Gas_tf, Atmos_input, overod, co21c, co21diag, co21r)
subroutine trans_sfc (Gas_tf, Atmos_input, overod, co21c, co21r)
 
!-------------------------------------------------------------------
!    compute "nearby  layer" transmission functions at level k 
!    ( tau(p(k),p(k))) in the frequency band at 15 um. include all
!    gases (co2, h2o, h2o cont) used in computing fluxes in this band.
!    the algorithm assumes that at pressures (p') near the pressure
!    at level k (p(k)), the transmission function may be written as:
!
!              tau(p',p(k)) = EXP(-alpha*SQRT(p'-p(k)))
!
!    with alpha determined by the boundary condition(s)
!    tau(p(k+1),p(k)) and tau(p(k-1),p(k)) = the values from  "normal"
!    calculations. An integration is performed over the "nearby" 
!    pressure layer to obtain tau(p(k),p(k)).
!    the computation is not done for levels from KSRAD to KMINH2O-1 (if
!    different), where it is assumed that the h2o transmissivities 
!    are near unity, and that the precomputed co2 transmissivities may
!    be used.
!     two "special case" transmissivities, viz.,
!    tau(p(KERAD),p(KERAD+1)) and tau(p(KERAD+1),p(KERAD)) are also 
!    evaluated using the above assumptions and an integration.
!-------------------------------------------------------------------

type(gas_tf_type),       intent(in)  :: Gas_tf
type(atmos_input_type),  intent(in)  :: Atmos_input
real, dimension (:,:,:), intent(in)  :: overod  
real, dimension (:,:),   intent(out) :: co21c, co21r

!--------------------------------------------------------------------
! intent(in) variables:
!
!    Gas_tf
!    Atmos_input
!    overod
!
! intent(out) variables:
!
!    co21c
!    co21r
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(Atmos_input%pflux,1),  &
                       size(Atmos_input%pflux,2)) :: pdfl

      real, dimension (size(Atmos_input%pflux,1),  &
                       size(Atmos_input%pflux,2),  &
                       size(Atmos_input%pflux,3)-1) :: pdfinv

      real, dimension (size(Atmos_input%pflux,1), &
                       size(Atmos_input%pflux,2),  &
                       size(Atmos_input%pflux,3)) ::         &
               press, pflux, alpa, alpb, ca, cb, delpr1, delpr2, rlog

      integer   ::  km, kmp1

!-------------------------------------------------------------------
!  local variables:
!
!    pdfl
!    pdfinv
!    press
!
!    k
!
!    kmp1
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    convert press and pflux to cgs.
!---------------------------------------------------------------------
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
      pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      pdfinv(:,:,KSRAD:KERAD) = 1.0/(pflux(:,:,KSRAD+1:KERAD+1) -   &
                                     pflux(:,:,KSRAD:KERAD) )

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      km= MAX(ixprkminh2o - 1, KSRAD)
      kmp1 = MAX(ixprkminh2o - 1, KSRAD+1)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      rlog  (:,:,km:KERAD)     = LOG(Gas_tf%co2nbl(:,:,km:KERAD)*   &
      overod(:,:,km+1:KERAD+1))
      delpr1(:,:,kmp1:KERAD)   = pdfinv (:,:,kmp1:KERAD)*  &
                                 (press(:,:,kmp1:KERAD) - &
                                  pflux(:,:,kmp1:KERAD)) 
      alpb  (:,:,kmp1:KERAD)   = -SQRT(delpr1(:,:,kmp1:KERAD))*  &
                                 rlog(:,:,kmp1:KERAD)
      delpr2(:,:,kmp1:KERAD+1) = pdfinv(:,:,kmp1-1:KERAD)*  &
                                 (pflux(:,:,kmp1:KERAD+1) -  &
                                  press(:,:,kmp1-1:KERAD)) 
      alpa  (:,:,km:KERAD)     = -SQRT(delpr2(:,:,km+1:KERAD+1))*  &
                                 rlog(:,:,km:KERAD)
      alpa  (:,:,KERAD+1)      = -rlog(:,:,KERAD)
      alpb  (:,:,KERAD+1)      = -rlog(:,:,KERAD)* &
                                 SQRT(pdfinv(:,:,KERAD)*&
                                 (pflux(:,:,KERAD+1) -  &
                                  press(:,:,KERAD-1)))
      ca(:,:,km:KERAD+1) = alpa(:,:,km:KERAD+1)*(-0.66667E+00 +  &
                           alpa(:,:,km:KERAD+1)*(0.25E+00 +   &
                           alpa(:,:,km:KERAD+1)*(-0.066667E+00)))
      cb(:,:,kmp1:KERAD+1) = alpb(:,:,kmp1:KERAD+1)*(-0.66667E+00 +  &
                             alpb(:,:,kmp1:KERAD+1)*(0.25E+00 +    &
                             alpb(:,:,kmp1:KERAD+1)*(-0.066667E+00)))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      pdfl(:,:) = pflux(:,:,KERAD+1) - pflux(:,:,KERAD)
      co21c(:,:        ) = 1.0E+00 +    &
                           (pdfl  (:,:      )*ca(:,:,KERAD+1) -    &
                           (press(:,:,KERAD) - pflux(:,:,KERAD))*   &
                            cb(:,:,KERAD))/   &
                            (pflux(:,:,KERAD+1) - press(:,:,KERAD))
      co21r(:,:        )   = 1.0E+00 +    &
                             ((pflux(:,:,KERAD+1) -  &
                               press(:,:,KERAD-1))*  &
                             cb(:,:,KERAD+1) -   &
                             (pflux(:,:,KERAD+1) - press(:,:,KERAD))*  &
                             ca(:,:,KERAD))/ &
                             (press(:,:,KERAD) - press(:,:,KERAD-1))

!-------------------------------------------------------------------




end subroutine trans_sfc


!####################################################################
! <SUBROUTINE NAME="put_co2_stdtf_for_gas_tf">
!  <OVERVIEW>
!   Assign co2 transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Assign co2 transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call put_co2_stdtf_for_gas_tf  (nf,        &
!                                      co251_o, co258_o,   &
!                                      cdt51_o, cdt58_o,   &
!                                      c2d51_o, c2d58_o)
!  </TEMPLATE>
!  <IN NAME="nf" TYPE="integer">
!   index variable
!  </IN>
!  <IN NAME="co251_o co258_o  cdt51_o cdt58_o  c2d51_o c2d58_o" TYPE="real">
!   CO2 transmission functions
!  </IN>
! </SUBROUTINE>
!
subroutine put_co2_stdtf_for_gas_tf (nf, co251_o, co258_o,   &
                                     cdt51_o, cdt58_o,   &
                                     c2d51_o, c2d58_o)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

integer,              intent(in)  :: nf
real, dimension(:,:), intent(in)  :: co251_o, co258_o,   &
                                     cdt51_o, cdt58_o,   &
                                     c2d51_o, c2d58_o

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     nf
!     co251_o
!     co258_o
!     cdt51_o
!     cdt58_o
!     c2d51_o
!     c2d58_o
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (nf == 1) then
        co251 = co251_o
        co258 = co258_o
        cdt51 = cdt51_o
        cdt58 = cdt58_o
        c2d51 = c2d51_o
        c2d58 = c2d58_o
      else if (nf == 5) then
        co211(:) = co251_o(:,1)
        co218(:) = co258_o(:,1)
      endif

!--------------------------------------------------------------------


end subroutine put_co2_stdtf_for_gas_tf




!#####################################################################
! <SUBROUTINE NAME="put_co2_nbltf_for_gas_tf">
!  <OVERVIEW>
!   Assign co2 transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Assign co2 transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call put_co2_nbltf_for_gas_tf  (nf,       &
!                                      co2m51_o, cdtm51_o, c2dm51_o, &
!                                      co2m58_o, cdtm58_o, c2dm58_o, &
!                                      co215nbps1_o, co215nbps8_o,     &
!                                      co2dt15nbps1_o, co2dt15nbps8_o, &
!                                      co2d2t15nbps1_o, co2d2t15nbps8_o )
!  </TEMPLATE>
!  <IN NAME="nf" TYPE="integer">
!   index variable
!  </IN>
!  <IN NAME="co2m51_o, cdtm51_o, c2dm51_o, co2m58_o, cdtm58_o, c2dm58_o, co215nbps1_o, co215nbps8_o" TYPE="real">
!   CO2 transmission functions
!  </IN>
!  <IN NAME="co2dt15nbps1_o, co2dt15nbps8_o, co2d2t15nbps1_o, co2d2t15nbps8_o" TYPE="real">
!   CO2 transmission functions
!  </IN>
! </SUBROUTINE>
!
subroutine put_co2_nbltf_for_gas_tf  (nf,       &
                                      co2m51_o, cdtm51_o, c2dm51_o, &
                                      co2m58_o, cdtm58_o, c2dm58_o, &
                                      co215nbps1_o, co215nbps8_o,     &
                                      co2dt15nbps1_o, co2dt15nbps8_o, &
                                      co2d2t15nbps1_o, co2d2t15nbps8_o )

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
integer,              intent(in)  :: nf
real, dimension(:,:), intent(in)  :: co2m51_o, cdtm51_o, c2dm51_o, &
                                     co2m58_o, cdtm58_o, c2dm58_o
real, dimension(:),   intent(in)  :: co215nbps1_o, co215nbps8_o,     &
                                     co2dt15nbps1_o, co2dt15nbps8_o, &
                                     co2d2t15nbps1_o, co2d2t15nbps8_o

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     nf
!     co2m51_o
!     cdtm51_o
!     c2dm51_o
!     co2m58_o
!     cdtm58_o
!     c2dm58_o
!     co215nbsp1_o
!     co215nbsp8_o
!     co2dt15nbsp1_o
!     co2dt15nbsp8_o
!     co2d2t15nbps1_o
!     co2d2t15nbps8_o
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   local variables:

      integer    :: k    ! do-loop index

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
             'module has not been initialized', FATAL )
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (nf == 1) then
        do k=KSRAD,KERAD
          co2m51(k) = co2m51_o(k,k+1)
          co2m58(k) = co2m58_o(k,k+1)
          cdtm51(k) = cdtm51_o(k,k+1)
          cdtm58(k) = cdtm58_o(k,k+1)
          c2dm51(k) = c2dm51_o(k,k+1)
          c2dm58(k) = c2dm58_o(k,k+1)
        end do
      endif
      if (nf >= 2 .and. nf <= 4) then
        co215nbps1(:,nf-1) = co215nbps1_o(:)
        co215nbps8(:,nf-1) = co215nbps8_o(:)
        co2dt15nbps1(:,nf-1) = co2dt15nbps1_o(:)
        co2dt15nbps8(:,nf-1) = co2dt15nbps8_o(:)
        co2d2t15nbps1(:,nf-1) = co2d2t15nbps1_o(:)
        co2d2t15nbps8(:,nf-1) = co2d2t15nbps8_o(:)
      endif

!--------------------------------------------------------------------

end subroutine put_co2_nbltf_for_gas_tf




!#####################################################################
! <SUBROUTINE NAME="put_ch4_stdtf_for_gas_tf">
!  <OVERVIEW>
!   Assign ch4 transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Assign ch4 transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call put_ch4_stdtf_for_gas_tf  (             &
!                                      ch451_o, ch458_o,   &
!                                      ch4dt51_o, ch4dt58_o,   &
!                                      ch4d2t51_o, ch4d2t58_o)
!  </TEMPLATE>
!  <IN NAME="nf" TYPE="integer">
!   index variable
!  </IN>
!  <IN NAME="ch451_o, ch458_o, ch4dt51_o, ch4dt58_o, ch4d2t51_o, ch4d2t58_o" TYPE="real">
!   CH4 transmission functions
!  </IN>
! </SUBROUTINE>
!
subroutine put_ch4_stdtf_for_gas_tf  (             &
                                       ch451_o, ch458_o,   &
                                       ch4dt51_o, ch4dt58_o,   &
                                       ch4d2t51_o, ch4d2t58_o)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

real, dimension (:,:), intent(in) :: ch451_o, ch458_o,     &    
                                     ch4dt51_o, ch4dt58_o, &    
                                     ch4d2t51_o, ch4d2t58_o

!------------------------------------------------------------------
!  intent(in) variables:
!
!    ch451_o
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ch451 = ch451_o
      ch458 = ch458_o
      ch4dt51 = ch4dt51_o
      ch4dt58 = ch4dt58_o
      ch4d2t51 = ch4d2t51_o
      ch4d2t58 = ch4d2t58_o

!--------------------------------------------------------------------


end subroutine put_ch4_stdtf_for_gas_tf



!#####################################################################
! <SUBROUTINE NAME="put_n2o_stdtf_for_gas_tf">
!  <OVERVIEW>
!   Assign n2o transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Assign n2o transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call put_n2o_stdtf_for_gas_tf  (nf,         &
!                                      n2o1_o, n2o8_o,   &
!                                      n2odt1_o, n2odt8_o,   &
!                                      n2od2t1_o, n2od2t8_o)
!  </TEMPLATE>
!  <IN NAME="nf" TYPE="integer">
!   index variable
!  </IN>
!  <IN NAME="n2o1_o, n2o8_o, n2odt1_o, n2odt8_o, n2od2t1_o, n2od2t8_o" TYPE="real">
!   N2O transmission functions
!  </IN>
! </SUBROUTINE>
!
subroutine put_n2o_stdtf_for_gas_tf  (nf,         &
                                      n2o1_o, n2o8_o,   &
                                      n2odt1_o, n2odt8_o,   &
                                      n2od2t1_o, n2od2t8_o)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

integer,               intent(in) :: nf
real, dimension (:,:), intent(in) :: n2o1_o, n2o8_o,     &    
                                     n2odt1_o, n2odt8_o, &    
                                     n2od2t1_o, n2od2t8_o

!------------------------------------------------------------------
!  intent(in) variables:
!
!    nf       
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (nf == 1) then
        n2o51 = n2o1_o
        n2o58 = n2o8_o
        n2odt51 = n2odt1_o
        n2odt58 = n2odt8_o
        n2od2t51 = n2od2t1_o
        n2od2t58 = n2od2t8_o
      else if (nf == 3) then
        n2o71 = n2o1_o
        n2o78 = n2o8_o
        n2odt71 = n2odt1_o
        n2odt78 = n2odt8_o
        n2od2t71 = n2od2t1_o
        n2od2t78 = n2od2t8_o
      else if (nf == 2) then
        n2o91 = n2o1_o
        n2o98 = n2o8_o
        n2odt91 = n2odt1_o
        n2odt98 = n2odt8_o
        n2od2t91 = n2od2t1_o
        n2od2t98 = n2od2t8_o
      endif

!---------------------------------------------------------------------

end subroutine put_n2o_stdtf_for_gas_tf



!#####################################################################
! <SUBROUTINE NAME="get_control_gas_tf">
!  <OVERVIEW>
!   Turn on gas transmission function flag
!  </OVERVIEW>
!  <DESCRIPTION>
!   Turn on gas transmission function flag
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_control_gas_tf (calc_co2, calc_ch4, calc_n2o)
!  </TEMPLATE>
!  <OUT NAME="calc_co2, calc_ch4, calc_n2o" TYPE="logical">
!   logical variables that determine whether gas transmission functions
!   should be calculated.
!  </OUT>
! </SUBROUTINE>
subroutine get_control_gas_tf (calc_co2, calc_ch4, calc_n2o)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

logical, intent(out), optional     :: calc_co2, calc_ch4, calc_n2o

!------------------------------------------------------------------
!  intent(out),optional variables:
!
!    calc_co2 
!    calc_ch4 
!    calc_n2o 
!
!-------------------------------------------------------------------


!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
            'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      calc_co2 = do_calcstdco2tfs
      calc_ch4 = do_calcstdch4tfs
      calc_n2o = do_calcstdn2otfs

!---------------------------------------------------------------------


end subroutine get_control_gas_tf


!####################################################################
! <SUBROUTINE NAME="gas_tf_dealloc">
!  <OVERVIEW>
!    gas_tf_dealloc deallocates the array components of the gas_tf_type
!    input variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    gas_tf_dealloc deallocates the array components of the gas_tf_type
!    input variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call gas_tf_dealloc (Gas_tf)  
!  </TEMPLATE>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   gas_tf_type variable containing information needed
!   to define the gas transmission functions
!  </INOUT>
! </SUBROUTINE>
subroutine gas_tf_dealloc (Gas_tf)

!------------------------------------------------------------------
!    gas_tf_dealloc deallocates the array components of the gas_tf_type
!    input variable.
!--------------------------------------------------------------------

type(gas_tf_type), intent(inout) :: Gas_tf

!---------------------------------------------------------------------
!  intent(inout) variables:
!
!     Gas_tf     gas_tf_type variable containing information needed
!                to define the gas transmission functions
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    deallocate the array components of Gas_tf.
!---------------------------------------------------------------------
      deallocate (Gas_tf%tdav)
      deallocate (Gas_tf%tlsqu         )
      deallocate (Gas_tf%tmpdiff       )
      deallocate (Gas_tf%tstdav        )
      deallocate (Gas_tf%co2nbl        )
      deallocate (Gas_tf%n2o9c         )
      deallocate (Gas_tf%tn2o17        )
      deallocate (Gas_tf%co2spnb       )
      deallocate (Gas_tf%a1            )
      deallocate (Gas_tf%a2            )

!--------------------------------------------------------------------


end subroutine gas_tf_dealloc



!###################################################################
! <SUBROUTINE NAME="gas_tf_end">
!  <OVERVIEW>
!   End and clean up gas tranmission function calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   End and clean up gas tranmission function calculation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call gas_tf_end
!  </TEMPLATE>
! </SUBROUTINE>
!   
subroutine gas_tf_end

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variable:

      integer :: tfsunit  ! unit number for i/o

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'gas_tf_mod', &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_writestdco2tfs) then
        if (mpp_pe() == mpp_root_pe() ) then
          tfsunit = open_restart_file ('stdco2tfs', action='write')
          write (tfsunit) valid_versions(nvalids)
          write (tfsunit) co2_name_save, co2_amount_save,   &
                           nstdlvls_save, kbegin_save, kend_save
          write (tfsunit) pd_save, plm_save, pa_save
          write (tfsunit) co215nbps1, co215nbps8, co2dt15nbps1,    & 
                          co2dt15nbps8, co2d2t15nbps1, co2d2t15nbps8
          write (tfsunit) co251, co258, cdt51, cdt58, c2d51, c2d58, &
                          co2m51, co2m58, cdtm51, cdtm58, c2dm51, c2dm58
          write (tfsunit) co211, co218
          call close_file (tfsunit)
        endif
      endif

      if (do_writestdn2otfs) then
        if (mpp_pe() == mpp_root_pe()) then
          tfsunit = open_restart_file ('stdn2otfs', action='write')
          write (tfsunit) valid_versions(nvalids)
          write (tfsunit) n2o_name_save, n2o_amount_save,   &
                          nstdlvls_save, kbegin_save, kend_save
          write (tfsunit) pd_save, plm_save, pa_save
          write (tfsunit) n2o51, n2o58, n2odt51, n2odt58, n2od2t51,   &
                          n2od2t58
          write (tfsunit) n2o71, n2o78, n2odt71, n2odt78, n2od2t71,  &
                          n2od2t78
          write (tfsunit) n2o91, n2o98, n2odt91, n2odt98, n2od2t91, &
                          n2od2t98
          call close_file (tfsunit)
        endif
      endif

      if (do_writestdch4tfs) then
        if (mpp_pe() == mpp_root_pe()) then
          tfsunit = open_restart_file ('stdch4tfs', action='write')
          write (tfsunit) valid_versions(nvalids)
          write (tfsunit) ch4_name_save, ch4_amount_save,   &
                          nstdlvls_save, kbegin_save, kend_save
          write (tfsunit) pd_save, plm_save, pa_save
          write (tfsunit) ch451, ch458, ch4dt51, ch4dt58, ch4d2t51,  &
                          ch4d2t58
          call close_file (tfsunit)
        endif
      endif

!--------------------------------------------------------------------
!    mark this module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------


end subroutine gas_tf_end




!###################################################################

! <SUBROUTINE NAME="process_co2_input_file">
!  <OVERVIEW>
!   Subroutine to process co2 input file information
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to process co2 input file information
!  </DESCRIPTION>
!  <TEMPLATE>
!   call process_co2_input_file (gas_name, gas_amount, nstdlvls,  &
!                                   kbegin, kend, pd, plm, pa)
!  </TEMPLATE>
!  <IN NAME="gas_name" TYPE="character">
!   Name of the gas specy
!  </IN>
!  <IN NAME="gas_amount" TYPE="real">
!   Amount of the gas specy
!  </IN>
!  <IN NAME="nstdlvls" TYPE="integer">
!   Number of standard levels
!  </IN>
!  <IN NAME="kbegin, kend" TYPE="integer">
!   Index of the starting and ending vertical levels
!  </IN>
!  <IN NAME="pd, plm, pa" TYPE="real">
!   Pressure coordinate variables, at boundaries, mid points.
!  </IN>
! </SUBROUTINE>
! 
subroutine process_co2_input_file (gas_name, gas_amount, nstdlvls,  &
                                   kbegin, kend, pd, plm, pa)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),   intent(in)      :: gas_name
real,               intent(in)      :: gas_amount
integer,            intent(in)      :: nstdlvls, kbegin, kend
real, dimension(:), intent(in)      :: pd, plm, pa

!------------------------------------------------------------------
!  intent(in) variables:
!
!    gas_name 
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:




      integer               ::  unit



!---------------------------------------------------------------------
!  local variables:
!
!       pd_file
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('gas_tf_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_readstdco2tfs) then

!--------------------------------------------------------------------
!    read the input tf file and verify that the current model config-
!    uration matches that for which the tfs were generated.
!--------------------------------------------------------------------
        unit = open_restart_file ('INPUT/stdco2tfs', action='read')
        call process_gas_input_file (gas_name, gas_amount, nstdlvls,  &
                                     kbegin, kend, pd, plm, pa, unit)

        read  (unit) co215nbps1, co215nbps8, co2dt15nbps1,           & 
                     co2dt15nbps8, co2d2t15nbps1, co2d2t15nbps8
        read  (unit) co251, co258, cdt51, cdt58, c2d51, c2d58,       &  
                     co2m51, co2m58, cdtm51, cdtm58, c2dm51, c2dm58
        read  (unit) co211, co218
        call close_file (unit)
      else if (do_writestdco2tfs) then

!--------------------------------------------------------------------
!    save the data necessary to write a tf file at the end of this job
!    if that is desired  
!--------------------------------------------------------------------
       co2_name_save = gas_name
       co2_amount_save = gas_amount
       nstdlvls_save = nstdlvls
       kbegin_save = kbegin
       kend_save = kend
       if (.not. (allocated(pd_save) ) ) then
         allocate ( pd_save(kbegin:kend))
         allocate ( plm_save(kbegin:kend))
         allocate ( pa_save(nstdlvls))
       endif
       pd_save(:) = pd(:)
       plm_save(:) = plm(:)
       pa_save(:) = pa(:)
     endif

!--------------------------------------------------------------------

end subroutine process_co2_input_file



!####################################################################
! <SUBROUTINE NAME="process_ch4_input_file">
!  <OVERVIEW>
!   Subroutine to process ch4 input file information
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to process ch4 input file information
!  </DESCRIPTION>
!  <TEMPLATE>
!   call process_ch4_input_file (gas_name, gas_amount, nstdlvls,  &
!                                   kbegin, kend, pd, plm, pa)
!  </TEMPLATE>
!  <IN NAME="gas_name" TYPE="character">
!   Name of the gas specy
!  </IN>
!  <IN NAME="gas_amount" TYPE="real">
!   Amount of the gas specy
!  </IN>
!  <IN NAME="nstdlvls" TYPE="integer">
!   Number of standard levels
!  </IN>
!  <IN NAME="kbegin, kend" TYPE="integer">
!   Index of the starting and ending vertical levels
!  </IN>
!  <IN NAME="pd, plm, pa" TYPE="real">
!   Pressure coordinate variables, at boundaries, mid points.
!  </IN>
! </SUBROUTINE>
!
subroutine process_ch4_input_file (gas_name, gas_amount, nstdlvls,  &
                                   kbegin, kend, pd, plm, pa)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),   intent(in)      :: gas_name
real,               intent(in)      :: gas_amount
integer,            intent(in)      :: nstdlvls, kbegin, kend
real, dimension(:), intent(in)      :: pd, plm, pa

!------------------------------------------------------------------
!  intent(in) variables:
!
!    gas_name 
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:




   integer               ::  unit



!---------------------------------------------------------------------
!  local variables:
!
!       pd_file
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'gas_tf_mod', &
              'module has not been initialized', FATAL )
      endif
 
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_readstdch4tfs) then

!--------------------------------------------------------------------
!    read the input tf file and verify that the current model config-
!    uration matches that for which the tfs were generated.
!--------------------------------------------------------------------
        unit = open_restart_file ('INPUT/stdch4tfs', action='read')
        call process_gas_input_file (gas_name, gas_amount, nstdlvls,  &
                                     kbegin, kend, pd, plm, pa, unit)
        read  (unit) ch451, ch458, ch4dt51, ch4dt58, ch4d2t51, ch4d2t58
        call close_file (unit)

!--------------------------------------------------------------------
!    save the data necessary to write a tf file at the end of this job
!    if that is desired  
!--------------------------------------------------------------------
      else if (do_writestdch4tfs) then
        ch4_name_save = gas_name
        ch4_amount_save = gas_amount
        nstdlvls_save = nstdlvls
        kbegin_save = kbegin
        kend_save = kend
        if (.not. (allocated(pd_save) ) ) then
          allocate ( pd_save(kbegin:kend))
          allocate ( plm_save(kbegin:kend))
          allocate ( pa_save(nstdlvls))
        endif
        pd_save(:) = pd(:)
        plm_save(:) = plm(:)
        pa_save(:) = pa(:)
      endif

!---------------------------------------------------------------------


end subroutine process_ch4_input_file



!####################################################################
! <SUBROUTINE NAME="process_n2o_input_file">
!  <OVERVIEW>
!   Subroutine to process n2o input file information
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to process n2o input file information
!  </DESCRIPTION>
!  <TEMPLATE>
!   call process_n2o_input_file (gas_name, gas_amount, nstdlvls,  &
!                                   kbegin, kend, pd, plm, pa)
!  </TEMPLATE>
!  <IN NAME="gas_name" TYPE="character">
!   Name of the gas specy
!  </IN>
!  <IN NAME="gas_amount" TYPE="real">
!   Amount of the gas specy
!  </IN>
!  <IN NAME="nstdlvls" TYPE="integer">
!   Number of standard levels
!  </IN>
!  <IN NAME="kbegin, kend" TYPE="integer">
!   Index of the starting and ending vertical levels
!  </IN>
!  <IN NAME="pd, plm, pa" TYPE="real">
!   Pressure coordinate variables, at boundaries, mid points.
!  </IN>
! </SUBROUTINE>
!
subroutine process_n2o_input_file (gas_name, gas_amount, nstdlvls,  &
                                   kbegin, kend, pd, plm, pa)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),   intent(in)      :: gas_name
real,               intent(in)      :: gas_amount
integer,            intent(in)      :: nstdlvls, kbegin, kend
real, dimension(:), intent(in)      :: pd, plm, pa

!------------------------------------------------------------------
!  intent(in) variables:
!
!    gas_name 
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:




   integer               ::  unit



!---------------------------------------------------------------------
!  local variables:
!
!       pd_file
!
!--------------------------------------------------------------------
        
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ( 'gas_tf_mod', &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_readstdn2otfs) then

!--------------------------------------------------------------------
!    read the input tf file and verify that the current model config-
!    uration matches that for which the tfs were generated.
!--------------------------------------------------------------------
        unit = open_restart_file ('INPUT/stdn2otfs', action='read')
        call process_gas_input_file (gas_name, gas_amount, nstdlvls,  &
                                     kbegin, kend, pd, plm, pa, unit)
        read  (unit) n2o51, n2o58, n2odt51, n2odt58, n2od2t51, n2od2t58
        read  (unit) n2o71, n2o78, n2odt71, n2odt78, n2od2t71, n2od2t78
        read  (unit) n2o91, n2o98, n2odt91, n2odt98, n2od2t91, n2od2t98
        call close_file (unit)

!--------------------------------------------------------------------
!    save the data necessary to write a tf file at the end of this job
!    if that is desired  
!--------------------------------------------------------------------
      else if (do_writestdn2otfs) then
        n2o_name_save = gas_name
        n2o_amount_save = gas_amount
        nstdlvls_save = nstdlvls
        kbegin_save = kbegin
        kend_save = kend
       if (.not. (allocated(pd_save) ) ) then
         allocate ( pd_save(kbegin:kend))
         allocate ( plm_save(kbegin:kend))
         allocate ( pa_save(nstdlvls))
       endif
       pd_save(:) = pd(:)
       plm_save(:) = plm(:)
       pa_save(:) = pa(:)
     endif

!--------------------------------------------------------------------


end subroutine process_n2o_input_file




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!####################################################################
! <SUBROUTINE NAME="ptz">
!  <OVERVIEW>
!   Subroutine to calculate temperatures at up to 200 user     
!   specified pressures.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate temperatures at up to 200 user     
!   specified pressures.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call ptz(plm, pd)
!  </TEMPLATE>
!  <IN NAME="plm" TYPE="real">
!   pressure at midpoint of layer (average of adjacent
!   pd values) 
!  </IN>
!  <IN NAME="pd" TYPE="real">
!   pressures (mb) for layer boundaries. (also known
!   as flux levels).
!  </IN>
! </SUBROUTINE>
!
subroutine ptz (plm, pd)

!--------------------------------------------------------------------
!    this program calculates temperatures at up to 200 user     
!    specified pressures. it makes use of an analytical       
!    function which approximates  the us standard atm(1976).  
!    this is calculated in function 'antemp', which is called
!    by ptz.  the form of the analytical function was    
!    suggested to me (S.B. Fels) in 1971 by richard s. lindzen. 
!
!--------------------------------------------------------------------

real, dimension(:), intent(in)    :: plm, pd

!--------------------------------------------------------------------
!  intent(in) variables:
!
!     pd:  pressures (mb) for layer boundaries. (also known
!          as flux levels).
!     plm  pressure at midpoint of layer (average of adjacent
!          pd values)
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension(:), allocatable  ::  press, altquad, tempquad, & 
                                           prsint, plmcgs,  tmpint

      real       ::  delzap = 0.5
      real       ::  dlogp, znint, dz, ht, rk1,  &
                     rk2, rk3, rk4
      integer    ::  k, nlay, nint, m

!--------------------------------------------------------------------
!   local variables:
!
!            plmcgs : plm in cgs units. needed for gtemp calc.
!            prsint: same as pd, but with indices reversed,index 1 at
!                    the surface, and index (nlev+1) at the top  level.
!            press: pressures used for quadratures (4 quad. pts.)  
!                     indices run as in prsint.
!         tempquad: temperature at quad. pts. index 1 is at the sfc.
!           tmpint: temperature at quad. pts, saved over both height
!                   and quadrature index. values come from tempquad.
!            tmpout: as in temp, but with indices reversed (index 1=
!                    temp at top data level),others are layer averages
!            altquad:  height values (km) generated by antemp. lowest
!                     index = surface.
!--------------------------------------------------------------------

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
      nlay = KERAD - KSRAD + 1

!-------------------------------------------------------------------
!
!-------------------------------------------------------------------
      allocate ( gtemp    (KSRAD:KERAD+1) )
      allocate ( stemp    (KSRAD:KERAD+1) )
      allocate ( plmcgs   (KSRAD:KERAD+1) )
      allocate ( press    (1:nlay+1) )
      allocate ( altquad  (1:nlay+1) )
      allocate ( tempquad (1:nlay+1) )
      allocate ( prsint   (1:nlay+1) )
      allocate ( tmpint   (1:nlay+1) )

!--------------------------------------------------------------------
!    the gtemp code below assumes plmcgs in cgs units
!--------------------------------------------------------------------

      plmcgs(KSRAD:KERAD+1) = pd(KSRAD:KERAD+1)*1.0E+03
      do k = KSRAD,KERAD 
        gtemp(k) = plmcgs(k)**0.2*(1.+plmcgs(k)/30000.)**0.8/1013250.
      enddo
      gtemp(KERAD+1) = 1.0
 
      altquad(1)=0.0
      tempquad(1)=antemp(0.0)
 
      do k=KSRAD, KERAD+1
        prsint(k)=plm(nlay+2-k)
      enddo

!-------------------------------------------------------------------
!    obtain layer-mean quantities by quadrature over the layers.
!    the calculation is made  to find the temperature at
!    the layer-mean (plm)
!
!    calculations are done (oddly!) 1 quad. interval at a time, with
!    each going from the sfc upward to the top layer.
!-------------------------------------------------------------------- 
      press(1) = prsint(1)
      do k=2, nlay+1
        press(k) = pd(nlay+KSRAD+1-k)
      enddo

!---------------------------------------------------------------------
!    press is the pressure at the quadrature point; alt and temp
!    are the heights and pressures for each
!    such quadrature point. these are saved as tmpint and a.
!      note that press is not explicitly saved.
!--------------------------------------------------------------------
      do k=1,nlay
 
!-------------------------------------------------------------------
!    establish computational levels between user levels at   
!    intervals of approximately 'delzap' km.                
!    special care is needed for the topmost layer, which usually
!    goes to zero pressure.
!    (special definition not needed; top pressure is nonzero)
!------------------------------------------------------------------
        dlogp=7.0*ALOG(press(k)/press(k+1))
        nint=dlogp/delzap
        nint=nint+1
        znint=nint

!------------------------------------------------------------------
!    the conversion factor is used to convert dz from cm (using the
!    model's values for rgas and grav) to km (as in this program)
!------------------------------------------------------------------
        dz  = 1.0E-05*(1.0E+02*RDGAS)*dlogp/(7.0*GRAV*znint)
        ht=altquad(k)
 
!---------------------------------------------------------------------
!    calculate height at next user level by means of       
!    runge-kutta integration.                   
!-------------------------------------------------------------------
        do m=1,nint
          rk1=antemp(ht)*dz
          rk2=antemp(ht+0.5*rk1)*dz
          rk3=antemp(ht+0.5*rk2)*dz
          rk4=antemp(ht+rk3)*dz
          ht=ht+0.16666667*(rk1+rk2+rk2+rk3+rk3+rk4)
        enddo
        altquad(k+1)=ht
        tempquad(k+1)=antemp(ht)
      enddo
 
!--------------------------------------------------------------------
!    save temperature (tmpint) at quadrature
!    points for layer-mean evaluations by simpsons rule.
!--------------------------------------------------------------------
      do k=1,nlay+1
        tmpint(k)=tempquad(k)
      enddo

!--------------------------------------------------------------------
!    end of quadrature loop.
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    stemp is layer-mean temperature with index 1 
!    at the top.applies for nq=5
!---------------------------------------------------------------------
      do k=KSRAD,KERAD+1
        stemp(k)=tmpint(nlay+KSRAD+1-k)
      enddo
 
!--------------------------------------------------------------------
      deallocate ( tmpint   )
      deallocate ( prsint   )
      deallocate ( tempquad )
      deallocate ( altquad  )
      deallocate ( press    )
      deallocate ( plmcgs   )

!------------------------------------------------------------------

end subroutine ptz



!###################################################################
! <SUBROUTINE NAME="antemp">
!  <OVERVIEW>
!   Subroutine to analytically calculate temperature profiles of
!   atmosphere with arbitrary levels
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to analytically calculate temperature profiles of
!   atmosphere with arbitrary levels
!  </DESCRIPTION>
!  <TEMPLATE>
!   temp = antemp(z)
!  </TEMPLATE>
!  <IN NAME="z" TYPE="real">
!   Height
!  </IN>
! </SUBROUTINE>
!
real function antemp (z)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
real, intent(in)   ::  z

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    z
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      real, dimension (10)   :: zb, delta
      real, dimension (11)   :: c
      real                   :: tstar, temp, expo, x, y, zlog, expp,   &
                                faclog
      integer                :: n, nlast

      data zb/   &
                             11.0,  20.0,  32.0,  47.0,  51.0,      &
                             71.0,  84.8520,  90.0,  91.0,  92.0/
      data c/  &
                             -6.5,   0.0,   1.0,   2.80,  0.0,    &
                             -2.80, -2.00,  0.0,   0.0,   0.0,  0.0/
      data delta/    &
                              0.3,   1.0,   1.0,   1.0,   1.0,    &
                              1.0,   1.0,   1.0,   1.0,   1.0/
 
      data tstar/    &
                  288.15/
 
!-------------------------------------------------------------------
!  local variables:
!
!    zb
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!
!-------------------------------------------------------------------
      temp = tstar + c(1)*z

!-------------------------------------------------------------------
!
!-------------------------------------------------------------------
      nlast = 10
      do n = 1,nlast
        expo = (z - zb(n))/delta(n)
        if (ABS(expo) .LE. 60.) then
          x = EXP(expo)
          y = x + 1.0/x
          zlog = ALOG(y)
        else
          zlog = ABS(expo)
        endif
        expp = zb(n)/delta(n)
        if (ABS(expp) .LE. 60.) then
          x = EXP(expp)
          y = x + 1.0/x
          faclog = ALOG(y)
        else
          faclog = ABS(expp)
        endif
        temp = temp + (c(n+1) - c(n))*0.5*(z + delta(n)*     &
               (zlog - faclog))
      enddo
      antemp = temp

!--------------------------------------------------------------------


end function antemp


!#####################################################################
! <SUBROUTINE NAME="transfn">
!  <OVERVIEW>
!   Subroutine to compute the temperature corrected co2 nearby layer
!   transmission functions
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute the temperature corrected co2 nearby layer
!   transmission functions
!  </DESCRIPTION>
!  <TEMPLATE>
!   call transfn( Gas_tf)
!  </TEMPLATE>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   The output variable of temperature corrected co2 transmission
!   functions
!  </INOUT>
! </SUBROUTINE>
!
subroutine transfn (Gas_tf)

!----------------------------------------------------------------------
!    transfn computes the temperature-corrected co2 nearby layer
!    transmission functions.
!    author: m. d. schwarzkopf
!    revised: 1/1/93
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

type(gas_tf_type), intent(inout) :: Gas_tf

!--------------------------------------------------------------------
!  intent(inout) variables:
!    
!     Gas_tf
!       co2nbl =  co2 transmission functions (not pressure-integrated) 
!                 for adjacent levels, over the 560-800 cm-1 range.
!       co2spnb = co2 transmission functions between a flux level and
!                 space, for each of (NBCO215) frequency bands in
!                 the 15 um range. used for cool-to-space calculations.
! 
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(Gas_tf%tdav,1), &
                       size(Gas_tf%tdav,2), &
                       size(Gas_tf%tdav,3)-1) :: co2m2d, co2md, co2mr

      real, dimension (size(Gas_tf%tdav,1), &
                       size(Gas_tf%tdav,2), &
                       size(Gas_tf%tdav,3)  ) :: dift                 

      real, dimension (size(Gas_tf%tdav,1), &
                       size(Gas_tf%tdav,2), &
                       size(Gas_tf%tdav,3), NBCO215  ) ::       &
                                    co215nb, co2dt15nb, co2d2t15nb

      integer    ::  inb, k

!---------------------------------------------------------------------
!  local variables:
!
!     co2m2d
!     co2md
!     co2mr
!     dift
!     co215nb
!     co2dt15nb
!     co2d2t15nb
!     inb
!     k
! 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    perform co2 pressure interpolation on all inputted transmission 
!    functions and temperature derivatives successively computing 
!    co2r, dco2dt, and d2cdt2.
!----------------------------------------------------------------------
      do inb=1,NBCO215
        do k=KSRAD,KERAD+1
          co215nb(:,:,k,inb)    = Gas_tf%a1(:,:)*co215nbps1(k,inb) +  &
                                  Gas_tf%a2(:,:)*co215nbps8(k,inb)
          co2dt15nb(:,:,k,inb)  = 1.0E-2*(Gas_tf%a1(:,:)* &
                                  co2dt15nbps1(k,inb) +&
                                  Gas_tf%a2(:,:)*co2dt15nbps8(k,inb))
          co2d2t15nb(:,:,k,inb) = 1.0E-3*(Gas_tf%a1(:,:)* &
                                  co2d2t15nbps1(k,inb)+&
                                  Gas_tf%a2(:,:)*co2d2t15nbps8(k,inb))
        enddo
      enddo

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do k=KSRAD,KERAD
        co2mr (:,:,k) = Gas_tf%a1(:,:)*co2m51(k) +  &
                        Gas_tf%a2(:,:)*co2m58(k)
        co2md (:,:,k) = 1.0E-02*(Gas_tf%a1(:,:)*cdtm51(k) +   &
                                 Gas_tf%a2(:,:)*cdtm58(k))
        co2m2d(:,:,k) = 1.0E-03*(Gas_tf%a1(:,:)*c2dm51(k) +   &
                                 Gas_tf%a2(:,:)*c2dm58(k))
      enddo

!----------------------------------------------------------------------
!    perform the temperature interpolation for these transmissivities
!----------------------------------------------------------------------
      dift(:,:,KSRAD+1:KERAD+1) = Gas_tf%tdav(:,:,KSRAD+1:KERAD+1)/   &
                                  Gas_tf%tstdav(:,:,KSRAD+1:KERAD+1)
      do inb=1,NBCO215
        Gas_tf%co2spnb(:,:,KSRAD,inb) = 1.0E+00
        Gas_tf%co2spnb(:,:,KSRAD+1:KERAD+1,inb) =     &
                                   co215nb(:,:,KSRAD+1:KERAD+1,inb) +  &
                                   dift(:,:,KSRAD+1:KERAD+1)*    &
                                  (co2dt15nb(:,:,KSRAD+1:KERAD+1,inb) +&
                                  0.5E+00*dift(:,:,KSRAD+1:KERAD+1)*   &
                                  co2d2t15nb(:,:,KSRAD+1:KERAD+1,inb))
        do k=KSRAD,KERAD+1
          Gas_tf%co2spnb(:,:,k,inb) = Gas_tf%co2spnb(:,:,k,inb)*  &
                                      (1.0E+00 -   &
                                      Gas_tf%tlsqu(:,:,KSRAD)) +    &
                                      Gas_tf%tlsqu(:,:,KSRAD)
        enddo
      enddo

!----------------------------------------------------------------------
!    compute special nearby layer transmission functions for combined
!    band in 15 um range. the transmissivities are not layer-averaged.
!----------------------------------------------------------------------
      Gas_tf%co2nbl(:,:,KSRAD:KERAD) = co2mr(:,:,KSRAD:KERAD) +   &
                                       Gas_tf%tmpdiff(:,:,KSRAD:KERAD)*&
                                       (co2md (:,:,KSRAD:KERAD) +   &
                                       0.5E+00*  &
                                       Gas_tf%tmpdiff(:,:,KSRAD:KERAD)*&
                                       co2m2d(:,:,KSRAD:KERAD))

      Gas_tf%co2nbl(:,:,KSRAD:KERAD) = Gas_tf%co2nbl(:,:,KSRAD:KERAD)*&
                                       (1.0E+00 -   &
                                       Gas_tf%tlsqu(:,:,KSRAD:KERAD)) +&
                                       Gas_tf%tlsqu(:,:,KSRAD:KERAD) 

!--------------------------------------------------------------------


end subroutine transfn



!###################################################################
! <SUBROUTINE NAME="process_gas_input_file">
!  <OVERVIEW>
!   Subroutine to process gas input file information
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to process gas input file information
!  </DESCRIPTION>
!  <TEMPLATE>
!   call process_gas_input_file (gas_name, gas_amount, nstdlvls,  &
!                                   kbegin, kend, pd, plm, pa, unit)
!  </TEMPLATE>
!  <IN NAME="gas_name" TYPE="character">
!   Name of the gas specy
!  </IN>
!  <IN NAME="gas_amount" TYPE="real">
!   Amount of the gas specy
!  </IN>
!  <IN NAME="nstdlvls" TYPE="integer">
!   Number of standard levels
!  </IN>
!  <IN NAME="kbegin, kend" TYPE="integer">
!   Index of the starting and ending vertical levels
!  </IN>
!  <IN NAME="pd, plm, pa" TYPE="real">
!   Pressure coordinate variables, at boundaries, mid points.
!  </IN>
!  <IN NAME="unit" TYPE="integer">
!   The input file descriptor
!  </IN>
! </SUBROUTINE>
! 
subroutine process_gas_input_file (gas_name, gas_amount, nstdlvls,  &
                                   kbegin, kend, pd, plm, pa, unit)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

character(len=*),   intent(in)      :: gas_name
real,               intent(in)      :: gas_amount
integer,            intent(in)      :: nstdlvls, kbegin, kend, unit
real, dimension(:), intent(in)      :: pd, plm, pa

!---------------------------------------------------------------------
! intent(in) variables:
!
!    gas_nmae
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension(:), allocatable :: pd_file, plm_file, pa_file

      logical               :: valid=.false.
      integer               :: k, n
      character(len=8)      :: gastf_version, gas_file
      real                  :: gas_amount_file
      integer               :: nstdlvls_file, kbegin_file, kend_file

!--------------------------------------------------------------------
!  local variables:
!
!     pd_file
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      read (unit) gastf_version
      do n=1, nvalids
        if (gastf_version == valid_versions(n)) then
          valid = .true.
          exit
        endif
      end do
      if (.not.valid) then
        call error_mesg ( 'gas_tf_mod', &
           ' old gastf file -- no info on its contents ---'//&
           ' for safety, please recreate by running this code'// &
           ' with do_calc and do_write of stdtfs activated and'//&
           ' save the file thus generated to be read in future jobs',  &
                                                                FATAL)
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      read (unit) gas_file, gas_amount_file, nstdlvls_file, &
                  kbegin_file, kend_file

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      if (trim(gas_file) /= trim(gas_name) ) then
        call error_mesg ( 'gas_tf_mod', &
      'inconsistency in gas name between input file and current job', &
                                                              FATAL)
      endif
      if (gas_amount /= gas_amount_file) then
        call error_mesg ( 'gas_tf_mod', &
      'inconsistency in gas amount between input file and current job',&
                                                               FATAL)
      endif
      if (nstdlvls /= nstdlvls_file) then
        call error_mesg ( 'gas_tf_mod', &
        'inconsistency in nstdlvls between input file and current job',&
                                                                FATAL)
      endif
      if (kbegin_file /= KSRAD) then
        call error_mesg ( 'gas_tf_mod', &
         'inconsistency in KSRAD between input file and current job',&
                                                                 FATAL)
      endif
      if (kend_file /= KERAD) then
        call error_mesg ( 'gas_tf_mod', &
           'inconsistency in KERAD between input file and current job',&
                                                               FATAL)
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      allocate ( pd_file (kbegin_file:kend_file) )
      allocate ( plm_file(kbegin_file:kend_file) )
      allocate ( pa_file (nstdlvls   ) )

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      read (unit) pd_file, plm_file, pa_file

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do k=KSRAD,KERAD
        if (pd(k) /= pd_file(k)) then
          call error_mesg ( 'gas_tf_mod', &
         'inconsistency in pd  between input file and current job'//&
              '  -- have the input files been constructed using'//&
              ' current model pressure levels ?',  FATAL)
        endif
        if (plm(k) /= plm_file(k)) then
          call error_mesg ( 'gas_tf_mod', &
         'inconsistency in plm  between input file and current job'//&
          ' -- have the input files been constructed using'//&
          ' current model pressure levels ?',  FATAL)
        endif
      end do     

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do k=1,nstdlvls
        if (pa(k) /= pa_file(k)) then
          call error_mesg ( 'gas_tf_mod', &
         'inconsistency in pa between input file and current job',&
                                                           FATAL)
        endif
      end do

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      deallocate ( pa_file  )
      deallocate ( plm_file )
      deallocate ( pd_file  )

!--------------------------------------------------------------------

end subroutine process_gas_input_file


!#####################################################################

 
      end module gas_tf_mod




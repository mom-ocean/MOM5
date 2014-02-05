                 module isccp_clouds_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!    isccp_clouds partitions the model cloud fields into the isccp
!    cloud categories, by cld top height and cld optical thickness
!    and provides netcdf output.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>


! shared modules:

use mpp_mod,                 only: input_nml_file
use fms_mod,                 only: fms_init, open_namelist_file, &
                                   write_version_number, mpp_pe, &
                                   mpp_root_pe, stdlog, file_exist,  &
                                   check_nml_error, error_mesg,   &
                                   FATAL, close_file
use time_manager_mod,        only: time_type, time_manager_init
use diag_manager_mod,        only: register_diag_field, send_data, &
                                   diag_manager_init

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    isccp_clouds partitions the model cloud fields into the isccp
!    cloud categories, by cld top height and cld optical thickness
!    and provides netcdf output.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: isccp_clouds.F90,v 19.0 2012/01/06 20:16:57 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public   isccp_clouds_init, isccp_output,               &
         isccp_cloudtypes, isccp_cloudtypes_stochastic, &
         isccp_clouds_end

private          &
!   called from isccp_clouds_init:
         diag_field_init, &
!   called from isccp_cloudtypes:
         ran0


!---------------------------------------------------------------------
!-------- namelist  ---------
!
!                         ISCCP CLOUD PROCESSING
!    
!
!       top_height
!
!                      integer variable indicating whether 
!                      or not to simulate 10.5 micron brightness
!                      temperatures to adjust top heights according
!                      to the emissivity of the cloud. 
!                     
!                      1 = adjust top height using both a computed
!                          infrared brightness temperature and the
!                          visible optical depth to adjust cloud top
!                          pressure. Note that this calculation is
!                          most appropriate to compare to ISCCP data
!                          during sunlit hours.
!       
!                      2 = do not adjust top height, that is cloud top
!                          pressure is the actual cloud top pressure
!                          in the model
!       
!                      3 = adjust top height using only the computed
!                          infrared brightness temperature. Note that
!                          this calculation is most appropriate to
!                          compare to ISCCP IR only algortihm (i.e.
!                          you can compare to nighttime ISCCP data
!                          with this option)
!                      
!       ncol           number of columns used in ISCCP cloud type
!                      simulations
!                      NOTE: This parameter is ignored when using 
!                      stochastic clouds. 
! 
!       isccp_taumin   minimum optical depth ISCCP can see
!                      [ dimensionless ]
!
!       emsfclw        assumed constant fraction of longwave emissivity
!                      of the surface [ dimensionless ]
! 
!       do_sunlit_only should ISCCP diagnostics be done during sunlit
!                      hours only? [ logical ]
!
!       overlap        variable indicating which overlap assumption to 
!                      use in ISCCP processing.
!
!                      NOTE THIS HAS NO IMPACT ON ANYTHING ELSE BUT
!                      ISCCP DIAGNOSIS FROM THIS SUBROUTINE
!
!                      NOTE: This parameter is ignored when using 
!                      stochastic clouds. 
!
!                      overlap = 1. means condensate in all levels 
!                                   is treated as part of the same cloud
!                                   i.e. maximum overlap
!
!                      overlap = 2. means condensate in adjacent levels 
!                                   is treated as different clouds
!                                   i.e. random overlap
!       
!                      overlap = 3. means condensate in adjacent levels 
!                                   is treated as part of the same cloud
!                                   i.e. maximum-random overlap
!
!      minColsInhomo   Minimum number of cloudy subcolumns required
!                      to do calculation of the inhomogeneity parameter
!                      epsilon = 1 - tau**/tau_bar
!
!                      where tau_bar = linear average of tau
!
!                            tau**  = exponential of the
!                                     linear average of logarithm of tau
!
!
!----------------------------------------------------------------------

integer  ::  ncol = 50
integer  ::  top_height = 1
real     ::  isccp_taumin = 0.3
real     ::  emsfclw = 0.94
logical  ::  do_sunlit_only = .false.
integer  ::  overlap = 2
integer  ::  minColsInhomo = 3

namelist /isccp_clouds_nml /  ncol, top_height, isccp_taumin, emsfclw,&
                              do_sunlit_only, overlap, minColsInhomo

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

integer, parameter :: numIsccpPressureIntervals     = 7, &
                      numIsccpOpticalDepthIntervals = 7

real, parameter     :: taumin = 1.E-06  ! minimum value allowed for 
                                        ! optical depth 
                                        ! [ dimensionless ]
real                :: qmin = 1.E-06    ! small number used in a couple
                                        ! of places in the code

!----------------------------------------------------------------------
!    diagnostics variables.     
!----------------------------------------------------------------------
character(len=5)    :: mod_name = 'isccp'
real                :: missing_value = -999.

integer :: id_deep,         id_cirrostratus,  id_cirrus,           &
           id_nimbostratus, id_altostratus,   id_altocumulus,      &
           id_stratus,      id_stratocumulus, id_cumulus,          &
           id_hithin,       id_midthin,       id_lowthin,          &
           id_high,         id_mid,           id_low,              &
           id_total,        id_allclouds,                          &          
           id_pc1tau0,id_pc1tau1,id_pc1tau2,id_pc1tau3,id_pc1tau4, &
           id_pc1tau5,id_pc1tau6, &
           id_pc2tau0,id_pc2tau1,id_pc2tau2,id_pc2tau3,id_pc2tau4, &
           id_pc2tau5,id_pc2tau6, &
           id_pc3tau0,id_pc3tau1,id_pc3tau2,id_pc3tau3,id_pc3tau4, &
           id_pc3tau5,id_pc3tau6, &
           id_pc4tau0,id_pc4tau1,id_pc4tau2,id_pc4tau3,id_pc4tau4, &
           id_pc4tau5,id_pc4tau6, &
           id_pc5tau0,id_pc5tau1,id_pc5tau2,id_pc5tau3,id_pc5tau4, &
           id_pc5tau5,id_pc5tau6, &
           id_pc6tau0,id_pc6tau1,id_pc6tau2,id_pc6tau3,id_pc6tau4, &
           id_pc6tau5,id_pc6tau6, &
           id_pc7tau0,id_pc7tau1,id_pc7tau2,id_pc7tau3,id_pc7tau4, &
           id_pc7tau5,id_pc7tau6, &
           id_nisccp, id_ninhomog, id_inhomogeneity

logical :: module_is_initialized =   .false.    ! module  initialized ?


!----------------------------------------------------------------------
!----------------------------------------------------------------------


  !
  ! Overloaded procedures
  !
  interface fluxToTb
    module procedure fluxToTb_1D, fluxToTb_2D, fluxToTb_3d
  end interface ! fluxToTb
  
  interface TbToFlux
    module procedure TbToFlux_1D, TbToFlux_2D, TbToFlux_3d
  end interface ! TbToFlux
  
  interface computeRadiance
    module procedure computeRadiance_1D, computeRadiance_2D
  end interface ! computeRadiance

                        contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="isccp_clouds_init">
!  <OVERVIEW>
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_clouds_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_clouds_init (axes, Time)

!---------------------------------------------------------------------
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!------------------------------------------------------------------

integer, dimension(4),   intent(in)              :: axes
type(time_type),         intent(in)              :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       axes             diagnostic variable axes
!       Time             current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: unit, io, ierr, logunit

!---------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      io       error status returned from io operation  
!      ierr     error code
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
      call time_manager_init
      call diag_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=isccp_clouds_nml, iostat=io)
      ierr = check_nml_error(io,'isccp_clouds_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=isccp_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'isccp_clouds_nml')
        enddo
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                       write (logunit, nml=isccp_clouds_nml)
 

!-------------------------------------------------------------------
!    initialize the netcdf diagnostics provided with this module.
!-------------------------------------------------------------------
      call diag_field_init (Time, axes)

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine isccp_clouds_init



!######################################################################
! <SUBROUTINE NAME="isccp_output">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_output (is, js, fq_isccp, npoints, 
!                      inhomogeneity_parameter, ninhomog, Time)
!  </TEMPLATE>
!  <IN NAME="is, js" TYPE="integer">
!   starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                      diagnostic output 
!  </IN>
!  <IN NAME="fq_isccp" TYPE="real">
!   matrix of fractional area covered by cloud
!                      types of a given optical depth and cloud
!                      top pressure range.  The matrix is 7x7 for
!                      7 cloud optical depths and 7 cloud top 
!                      pressure ranges
!  </IN>
!  <IN NAME="npoints" TYPE="real">
!   flag indicating whether isccp cloud is present
!                      in column (cloud + daylight needed)
!  </IN>
!  <OUT NAME="inhomogeneity_parameter" TYPE="real">
!   Cloud inhomogeneity parameter (between 0 and 1 if valid
!   point, -1. if not computed at this point [ dimensionless ]
!  </OUT>
!  <IN NAME="ninhomog" TYPE="real">
!   flag indicating cloud inhomogeneity calculations have been
!      performed [1.=True, 0.=False]
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_output (is, js, fq_isccp, npoints, &
                         inhomogeneity_parameter, ninhomog, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                      intent(in)   :: is,js
real, dimension(:,:,:,:),     intent(in)   :: fq_isccp
real, dimension(:,:),         intent(in)   :: npoints, ninhomog, &
                                              inhomogeneity_parameter
type(time_type),              intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      fq_isccp        matrix of fractional area covered by cloud
!                      types of a given optical depth and cloud
!                      top pressure range.  The matrix is 7x7 for
!                      7 cloud optical depths and 7 cloud top 
!                      pressure ranges
!      npoints         flag indicating whether isccp cloud is present
!                      in column (cloud + daylight needed)
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!local variable:
     
     real, dimension(size(npoints,1),size(npoints,2)) :: tmpmat
     
     logical :: used     !  flag returned from send_data indicating
                         !  whether diag_manager_mod has received 
                         !  data that was sent

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------

      !------------------------------------------------------------
      ! do individual types
      
      used = send_data (id_pc1tau0, fq_isccp(:,:,1,1), Time, &
                        is, js )
      used = send_data (id_pc1tau1, fq_isccp(:,:,2,1), Time, &
                        is, js )
      used = send_data (id_pc1tau2, fq_isccp(:,:,3,1), Time, &
                        is, js )
      used = send_data (id_pc1tau3, fq_isccp(:,:,4,1), Time, &
                        is, js )
      used = send_data (id_pc1tau4, fq_isccp(:,:,5,1), Time, &
                        is, js )
      used = send_data (id_pc1tau5, fq_isccp(:,:,6,1), Time, &
                        is, js )
      used = send_data (id_pc1tau6, fq_isccp(:,:,7,1), Time, &
                        is, js )
      used = send_data (id_pc2tau0, fq_isccp(:,:,1,2), Time, &
                        is, js )
      used = send_data (id_pc2tau1, fq_isccp(:,:,2,2), Time, &
                        is, js )
      used = send_data (id_pc2tau2, fq_isccp(:,:,3,2), Time, &
                        is, js )
      used = send_data (id_pc2tau3, fq_isccp(:,:,4,2), Time, &
                        is, js )
      used = send_data (id_pc2tau4, fq_isccp(:,:,5,2), Time, &
                        is, js )
      used = send_data (id_pc2tau5, fq_isccp(:,:,6,2), Time, &
                        is, js )
      used = send_data (id_pc2tau6, fq_isccp(:,:,7,2), Time, &
                        is, js )
      used = send_data (id_pc3tau0, fq_isccp(:,:,1,3), Time, &
                        is, js )
      used = send_data (id_pc3tau1, fq_isccp(:,:,2,3), Time, &
                        is, js )
      used = send_data (id_pc3tau2, fq_isccp(:,:,3,3), Time, &
                        is, js )
      used = send_data (id_pc3tau3, fq_isccp(:,:,4,3), Time, &
                        is, js )
      used = send_data (id_pc3tau4, fq_isccp(:,:,5,3), Time, &
                        is, js )
      used = send_data (id_pc3tau5, fq_isccp(:,:,6,3), Time, &
                        is, js )
      used = send_data (id_pc3tau6, fq_isccp(:,:,7,3), Time, &
                        is, js )
      used = send_data (id_pc4tau0, fq_isccp(:,:,1,4), Time, &
                        is, js )
      used = send_data (id_pc4tau1, fq_isccp(:,:,2,4), Time, &
                        is, js )
      used = send_data (id_pc4tau2, fq_isccp(:,:,3,4), Time, &
                        is, js )
      used = send_data (id_pc4tau3, fq_isccp(:,:,4,4), Time, &
                        is, js )
      used = send_data (id_pc4tau4, fq_isccp(:,:,5,4), Time, &
                        is, js )
      used = send_data (id_pc4tau5, fq_isccp(:,:,6,4), Time, &
                        is, js )
      used = send_data (id_pc4tau6, fq_isccp(:,:,7,4), Time, &
                        is, js )
      used = send_data (id_pc5tau0, fq_isccp(:,:,1,5), Time, &
                        is, js )
      used = send_data (id_pc5tau1, fq_isccp(:,:,2,5), Time, &
                        is, js )
      used = send_data (id_pc5tau2, fq_isccp(:,:,3,5), Time, &
                        is, js )
      used = send_data (id_pc5tau3, fq_isccp(:,:,4,5), Time, &
                        is, js )
      used = send_data (id_pc5tau4, fq_isccp(:,:,5,5), Time, &
                        is, js )
      used = send_data (id_pc5tau5, fq_isccp(:,:,6,5), Time, &
                        is, js )
      used = send_data (id_pc5tau6, fq_isccp(:,:,7,5), Time, &
                        is, js )
      used = send_data (id_pc6tau0, fq_isccp(:,:,1,6), Time, &
                        is, js )
      used = send_data (id_pc6tau1, fq_isccp(:,:,2,6), Time, &
                        is, js )
      used = send_data (id_pc6tau2, fq_isccp(:,:,3,6), Time, &
                        is, js )
      used = send_data (id_pc6tau3, fq_isccp(:,:,4,6), Time, &
                        is, js )
      used = send_data (id_pc6tau4, fq_isccp(:,:,5,6), Time, &
                        is, js )
      used = send_data (id_pc6tau5, fq_isccp(:,:,6,6), Time, &
                        is, js )
      used = send_data (id_pc6tau6, fq_isccp(:,:,7,6), Time, &
                        is, js )
      used = send_data (id_pc7tau0, fq_isccp(:,:,1,7), Time, &
                        is, js )
      used = send_data (id_pc7tau1, fq_isccp(:,:,2,7), Time, &
                        is, js )
      used = send_data (id_pc7tau2, fq_isccp(:,:,3,7), Time, &
                        is, js )
      used = send_data (id_pc7tau3, fq_isccp(:,:,4,7), Time, &
                        is, js )
      used = send_data (id_pc7tau4, fq_isccp(:,:,5,7), Time, &
                        is, js )
      used = send_data (id_pc7tau5, fq_isccp(:,:,6,7), Time, &
                        is, js )
      used = send_data (id_pc7tau6, fq_isccp(:,:,7,7), Time, &
                        is, js )
      used = send_data (id_nisccp, npoints, Time, is, js )

      used = send_data (id_ninhomog, ninhomog, Time, is, js )
     
      used = send_data (id_inhomogeneity, inhomogeneity_parameter, &
                        Time, is, js )

      !------------------------------------------------------------
      ! do summed fields

      !do hithin
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,1:1,1:3), dim = 4), dim = 3)
      used = send_data (id_hithin, tmpmat(:,:), Time, is, js )
      !do cirrus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:3,1:3), dim = 4), dim = 3)
      used = send_data (id_cirrus, tmpmat(:,:), Time, is, js )
      
      !do cirrostratus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,4:5,1:3), dim = 4), dim = 3)
      used = send_data (id_cirrostratus, tmpmat(:,:), Time, is, js )
      
      !do deep
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,6:7,1:3), dim = 4), dim = 3)
      used = send_data (id_deep, tmpmat(:,:), Time, is, js )
      
      !do high
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:7,1:3), dim = 4), dim = 3)
      used = send_data (id_high, tmpmat(:,:), Time, is, js )
      
      !do midthin
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,1:1,4:5), dim = 4), dim = 3)
      used = send_data (id_midthin, tmpmat(:,:), Time, is, js )
      
      !do altocumulus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:3,4:5), dim = 4), dim = 3)
      used = send_data (id_altocumulus, tmpmat(:,:), Time, is, js )
      
      !do altostratus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,4:5,4:5), dim = 4), dim = 3)
      used = send_data (id_altostratus, tmpmat(:,:), Time, is, js )
      
      !do nimbostratus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,6:7,4:5), dim = 4), dim = 3)
      used = send_data (id_nimbostratus, tmpmat(:,:), Time, is, js )
      
      !do mid
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:7,4:5), dim = 4), dim = 3)
      used = send_data (id_mid, tmpmat(:,:), Time, is, js )
      
      !do lowthin
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,1:1,6:7), dim = 4), dim = 3)
      used = send_data (id_lowthin, tmpmat(:,:), Time, is, js )
      
      !do cumulus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:3,6:7), dim = 4), dim = 3)
      used = send_data (id_cumulus, tmpmat(:,:), Time, is, js )
      
      !do stratocumulus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,4:5,6:7), dim = 4), dim = 3)
      used = send_data (id_stratocumulus, tmpmat(:,:), Time, is, js )
      
      !do stratus
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,6:7,6:7), dim = 4), dim = 3)
      used = send_data (id_stratus, tmpmat(:,:), Time, is, js )
      
      !do low
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:7,6:7), dim = 4), dim = 3)
      used = send_data (id_low, tmpmat(:,:), Time, is, js )
            
      !do total
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,2:7,1:7), dim = 4), dim = 3)
      used = send_data (id_total, tmpmat(:,:), Time, is, js )
             
      !do all clouds
      tmpmat(:,:) = sum(sum(fq_isccp(:,:,1:7,1:7), dim = 4), dim = 3)
      used = send_data (id_allclouds, tmpmat(:,:), Time, is, js )
       
!--------------------------------------------------------------------
         

end subroutine isccp_output


!######################################################################
! <SUBROUTINE NAME="isccp_cloudtypes">
!  <OVERVIEW>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_cloudtypes (sunlit, pfull, phalf, qv, at, skt, cc, &
!                          dtau_s, dem_s, fq_isccp, nisccp,&
!                          inhomogeneity_parameter, ninhomog)
!  </TEMPLATE>
!  <IN NAME="sunlit" TYPE="integer">
!   integer flag indicating whether or not a given point is sunlit
!                        [1 = True, 0 = False ]
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!   pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!  </IN>
!  <IN NAME="qv" TYPE="real">
!   water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="at" TYPE="real">
!   temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="skt" TYPE="real">
!   skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <INOUT NAME="cc" TYPE="real">
!   cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!  </INOUT>
!  <INOUT NAME="dtau_s" TYPE="real">
!   mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!  </INOUT>
!  <INOUT NAME="dem_s" TYPE="real">
!   10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!  </INOUT>
!  <OUT NAME="fq_isccp" TYPE="real">
!   matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!  </OUT>
!  <OUT NAME="nisccp" TYPE="real">
!   real flag indicating whether or not isccp_cloudtypes produced
!                       valid output [ 1.=True, 0.=False ]
!  </OUT>
!  <OUT NAME="inhomogeneity_parameter" TYPE="real">
!   Cloud inhomogeneity parameter (between 0 and 1 if valid
!   point, -1. if not computed at this point [ dimensionless ]
!  </OUT>
!  <IN NAME="ninhomog" TYPE="real">
!   flag indicating cloud inhomogeneity calculations have been
!      performed [1.=True, 0.=False]
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_cloudtypes (sunlit, pfull, phalf, qv, at, skt, cc, &
                             dtau_s, dem_s, fq_isccp, nisccp, &
                             inhomogeneity_parameter, ninhomog)

!---------------------------------------------------------------------
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!
!---------------------------------------------------------------------

integer, dimension(:,:),  intent(in)      :: sunlit
real,  dimension(:,:,:),  intent(in)      :: pfull, phalf, qv, at
real,    dimension(:,:),  intent(in)      :: skt
real,  dimension(:,:,:),  intent(in)      :: cc, dtau_s, dem_s
real,dimension(:,:,:,:),  intent(out)     :: fq_isccp
!  Last two dimensions should be numIsccpOpticalDepthIntervals, numIsccpPressureIntervals (7, 7)
real,    dimension(:,:),  intent(out)     :: nisccp
real,  dimension(:,:),    intent(out)     :: inhomogeneity_parameter, &
                                             ninhomog

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       sunlit           integer indicating whether or not a given
!                        point is sunlit
!       pfull            pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model [ Pa ]
!       phalf            pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!                        [ Pa ]
!       qv               water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!                        [ kg vapor / kg air ]
!       at               temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!       skt              skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!
!   intent(inout) variables:
!
!       cc               cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!       dtau_s           mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!       dem_s            10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!
!       NOTE :  OPTION TO RUN WITH CONVECTIVE CLOUDS IS NOT
!               IMPLEMENTED YET
!
!       conv             convective cloud cover in each model 
!                        level (fraction) this includes convective 
!                        clouds if any
!  
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!                         
!       dtau_c           mean 0.67 micron optical depth of convective
!                        clouds in each model level
!
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!
!   intent(out) variable:
!
!       fq_isccp        matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!
!       nisccp          integer indicating whether or not isccp diagnostics
!                       were calculated for this point
!
!       inhomogeneity_parameter
!
!                         = 1 - tau**/tau_bar
!
!                      where tau_bar = linear average of tau
!
!                            tau**  = exponential of the
!                                     linear average of logarithm of tau
!
!       ninhomog       flag indicating that the inhomogeneity 
!                      parameter was calculated
!
!---------------------------------------------------------------------
! Local variables
!   The ISCCP simulator code is vectorized over one set of input GCM grid cells
!   The input variables here are over an x-y window, so we'll loop 
!   over y. 
    integer :: nPoints, nYPoints, nLev, i, j
    real,    dimension(size(pFull, 1), nCol, size(pFull, 3)) :: frac_out, dtau, dem
    real,    dimension(size(pFull, 1),       size(pFull, 3)) :: conv, strat
    real,    dimension(size(pFull, 1),       size(pFull, 3)) :: dem_wv
    integer, dimension(size(pFull, 1))                       :: seed
    real,    dimension(size(pFull, 1), nCol)                 :: boxPtop, boxTau
    
!---------------------------------------------------------------------
    nPoints  = size(pFull, 1)
    nYPoints = size(pFull, 2) 
    nLev     = size(pFull, 3)
    
    if(.not. module_is_initialized) &
      call error_mesg("isccp_clouds_mod", "module has not been initialized", FATAL)
    
    !
    ! Don't compute statistics if do_sunlit_only flag is set and this point is in darkness
    !   (sunlit = 0).  
    !
    do j = 1, nYPoints
      if(do_sunlit_only .and. .not. any(sunlit(:, j) == 1)) then
        fq_isccp(:, j, :, :) = 0.
        nisccp(:, j) = 0
        inhomogeneity_parameter(:,j) = 0.
        ninhomog(:,j) = 0.
      else
        strat(:, :) = cc(:, j, :)
        conv (:, :) = 0.
        ! Make sure all the values are resonable
        where (strat(:, :) < 0.) strat(:, :) = 0
        where (strat(:, :) > 1.) strat(:, :) = 1.
        !
        ! Generate sub-cloud structures
        !
        seed(:) = (pfull(:, j, nlev) - int(pfull(:, j, nlev))) * 100 + 1
        call scops(strat, conv, seed, frac_out)
        
        !
        ! Take scops predictions of cloud fraction and fill in emmissivity and optical 
        !   depth arrays
        !
        where(nint(frac_out(:, :, :)) == 0)
          dem (:, :, :) = 0. 
          dtau(:, :, :) = 0. 
        elsewhere(nint(frac_out(:, :, :)) == 1)
          dem (:, :, :) = spread(dem_s (:, j, :), dim = 2, nCopies = nCol)
          dtau(:, :, :) = spread(dtau_s(:, j, :), dim = 2, nCopies = nCol) 
        end where 
        ! Make sure all the values are resonable
        where (dtau(:, :, :) < 0.) dtau(:, :, :) = 0
        where (dem (:, :, :) < 0.) dem (:, :, :) = 0
        where (dem (:, :, :) > 1.) dem (:, :, :) = 1.
        
        if(top_height == 1 .or. top_height == 3) then 
          ! 
          ! We're looking for adjusted cloud tops. Compute water vapor emissivity
          !
          call computeWaterVaporEmissivity(pfull(:, j, :), phalf(:, j, :), &
                                           qv(:, j, :),    at(:, j, :), dem_wv)
          
          ! Call Icarus...
          call icarus(dtau, pFull(:, j, :),                  & 
                      dem, dem_wv,  at(:, j, :),  skt(:, j), &
                       (/ (emsfclw, j = 1, nPoints) /),      &
                      boxtau = boxtau, boxptop = boxptop)
        else 
          ! We're asking for the real cloud tops. 
          !   We don't correct very optically thin clouds either. 
          call icarus(dtau, pFull(:, j, :),                  & 
                      boxtau = boxtau, boxptop = boxptop)
        end if 
        
        ! Compute histograms
        fq_isccp(:, j, :, :) = computeIsccpJointHistograms(boxtau, boxptop, sunlit(:, j))
        nisccp(:, j) = 1

        ! Compute inhomogeneity parameter
        inhomogeneity_parameter(:,j) = computeInhomogeneityParameter(boxtau,boxptop,sunlit(:,j))

        where(inhomogeneity_parameter(:,j)<-0.5)
              inhomogeneity_parameter(:,j) = 0.
              ninhomog(:,j) = 0.
        elsewhere
              ninhomog(:,j) = 1.
        endwhere   
        
        ! Zero out frequency histograms if the point is dark and we want statistics
        !    only for sunlit points. 
        if(do_sunlit_only) then
          where(sunlit(:, j) == 0) 
            nisccp(:, j) = 0
            inhomogeneity_parameter(:,j) = 0.
            ninhomog(:,j) = 0.
          endwhere
          do i = 1, nPoints
            if(sunlit(i, j) == 0) fq_isccp(i, j, :, :) = 0. 
          end do 
        end if 
        
      end if
    end do

end subroutine isccp_cloudtypes

!######################################################################
! <SUBROUTINE NAME="isccp_cloudtypes_stochastic">
!  <OVERVIEW>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure). 
!    This version uses the columns generated for the McICA treatment
!    of radiation.  
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure). 
!    For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!    This version uses the columns generated for the McICA treatment
!    of radiation; overlap is imposed in the "cloud generator" that 
!    takes the place of SCOPS, and internal inhomogeneity can be 
!    added too.  
!
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_cloudtypes (sunlit, pfull, phalf, qv, at, skt, cc, &
!                          dtau_s, dem_s, fq_isccp, nisccp,&
!                          inhomogeneity_parameter, ninhomog)
!  </TEMPLATE>
!  <IN NAME="sunlit" TYPE="integer">
!   integer flag indicating whether or not a given point is sunlit
!                        [1 = True, 0 = False ]
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!   pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!  </IN>
!  <IN NAME="qv" TYPE="real">
!   water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="at" TYPE="real">
!   temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="skt" TYPE="real">
!   skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <INOUT NAME="cc" TYPE="real">
!   cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!  </INOUT>
!  <INOUT NAME="dtau_s" TYPE="real">
!   mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!  </INOUT>
!  <INOUT NAME="dem_s" TYPE="real">
!   10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!  </INOUT>
!  <OUT NAME="fq_isccp" TYPE="real">
!   matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!  </OUT>
!  <OUT NAME="nisccp" TYPE="real">
!   real flag indicating whether or not isccp_cloudtypes produced
!                       valid output [ 1.=True, 0.=False ]
!  </OUT>
!  <OUT NAME="inhomogeneity_parameter" TYPE="real">
!   Cloud inhomogeneity parameter (between 0 and 1 if valid
!   point, -1. if not computed at this point [ dimensionless ]
!  </OUT>
!  <IN NAME="ninhomog" TYPE="real">
!   flag indicating cloud inhomogeneity calculations have been
!      performed [1.=True, 0.=False]
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_cloudtypes_stochastic (sunlit, pfull, phalf, qv, at, skt, cc, &
                                        dtau_s, dem_s, fq_isccp, nisccp, &
                                        inhomogeneity_parameter, ninhomog)

!---------------------------------------------------------------------
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure).
!    For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!    This version uses the columns generated for the McICA treatment
!    of radiation; overlap is imposed in the "cloud generator" that 
!    takes the place of SCOPS, and internal inhomogeneity can be 
!    added too.  
!
!---------------------------------------------------------------------

integer, dimension(:,:),     intent(in)  :: sunlit
real,    dimension(:,:,:),   intent(in)  :: pfull, phalf, qv, at
real,    dimension(:,:),     intent(in)  :: skt
real,    dimension(:,:,:,:), intent(in)  :: cc, dtau_s, dem_s
real,    dimension(:,:,:,:), intent(out) :: fq_isccp
!  Last two dimensions should be numIsccpOpticalDepthIntervals, numIsccpPressureIntervals (7, 7)
real,    dimension(:,:),     intent(out) :: nisccp
real,  dimension(:,:),    intent(out)    :: inhomogeneity_parameter, &
                                            ninhomog      

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       sunlit           integer indicating whether or not a given
!                        point is sunlit
!       pfull            pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model [ Pa ]
!       phalf            pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!                        [ Pa ]
!       qv               water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!                        [ kg vapor / kg air ]
!       at               temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!       skt              skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!
!   intent(inout) variables:
!
!       cc               cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!       dtau_s           mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!       dem_s            10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!
!
!   intent(out) variable:
!
!       fq_isccp        matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!
!       nisccp          integer indicating whether or not isccp diagnostics
!                       were calculated for this point
!
!       inhomogeneity_parameter
!
!                         = 1 - tau**/tau_bar
!
!                      where tau_bar = linear average of tau
!
!                            tau**  = exponential of the
!                                     linear average of logarithm of tau
!
!
!       ninhomog       flag indicating that the inhomogeneity 
!                      parameter was calculated
!
!
!---------------------------------------------------------------------
! Local variables
!   The ISCCP simulator code is vectorized over one set of input GCM grid cells
!   The input variables here are over an x-y window, so we'll loop 
!   over y. 
    integer :: nPoints, nYPoints, nLev, nColumns, i, j, col
    real,    dimension(size(cc, 1), size(cc, 4), size(cc, 3)) :: dtau, dem
    real,    dimension(size(cc, 1), size(cc, 3))              :: dem_wv
    real,    dimension(size(cc, 1), size(cc, 4))              :: boxPtop, boxTau
    
!---------------------------------------------------------------------
    nPoints  = size(cc, 1)
    nYPoints = size(cc, 2) 
    nLev     = size(cc, 3)
    nColumns = size(cc, 4)
    
    if(.not. module_is_initialized) &
      call error_mesg("isccp_clouds_mod", "module has not been initialized", FATAL)
    
    !
    ! Don't compute statistics if do_sunlit_only flag is set and this point is in darkness
    !   (sunlit = 0).  
    !
    do j = 1, nYPoints
      if(do_sunlit_only .and. .not. any(sunlit(:, j) == 1)) then
        fq_isccp(:, j, :, :) = 0.
        nisccp(:, j) = 0
        inhomogeneity_parameter(:,j) = 0.
        ninhomog(:,j) = 0.
      else
        !
        ! Reorder optical depth and emissivity arrays for this set 
        !   of Y points. In calling routines the order is x, y, lev, col; 
        !   icarus expects x, col, lev.  
        !
        do col = 1, nColumns
          dtau(:, col, :) = dtau_s(:, j, :, col)
          dem (:, col, :) = dem_s (:, j, :, col)
        end do
        !
        ! Make sure all the values are resonable
        !
        where (dtau(:, :, :) < 0.) dtau(:, :, :) = 0
        where (dem (:, :, :) < 0.) dem (:, :, :) = 0
        where (dem (:, :, :) > 1.) dem (:, :, :) = 1.
        
        if(top_height == 1 .or. top_height == 3) then 
          ! 
          ! We're looking for adjusted cloud tops. Compute water vapor emissivity
          !
          call computeWaterVaporEmissivity(pfull(:, j, :), phalf(:, j, :), &
                                           qv(:, j, :),    at(:, j, :), dem_wv)
          
          ! Call Icarus...
          call icarus(dtau, pFull(:, j, :),                  & 
                      dem, dem_wv,  at(:, j, :),  skt(:, j), &
                       (/ (emsfclw, j = 1, nPoints) /),       &
                      boxtau = boxtau, boxptop = boxptop)
        else 
          ! We're asking for the real cloud tops. 
          !   We don't correct very optically thin clouds either. 
          call icarus(dtau, pFull(:, j, :),                 & 
                      boxtau = boxtau, boxptop = boxptop)
        end if 
        
        ! Compute histograms
        fq_isccp(:, j, :, :) = computeIsccpJointHistograms(boxtau, boxptop, sunlit(:, j))
        nisccp(:, j) = 1

        ! Compute inhomogeneity parameter
        inhomogeneity_parameter(:,j) = computeInhomogeneityParameter(boxtau,boxptop,sunlit(:,j))

        where(inhomogeneity_parameter(:,j)<-0.5)
              inhomogeneity_parameter(:,j) = 0.
              ninhomog(:,j) = 0.
        elsewhere
              ninhomog(:,j) = 1.
        endwhere   
        
        ! Zero out frequency histograms if the point is dark and we want statistics
        !    only for sunlit points. 
        if(do_sunlit_only) then
          where(sunlit(:, j) == 0) 
            nisccp(:, j) = 0
            inhomogeneity_parameter(:,j) = 0.
            ninhomog(:,j) = 0.
          endwhere
          do i = 1, nPoints
            if(sunlit(i, j) == 0) fq_isccp(i, j, :, :) = 0. 
          end do 
        end if 
        
      end if
    end do

end subroutine isccp_cloudtypes_stochastic
!###################################################################

! <SUBROUTINE NAME="isccp_clouds_end">
!  <OVERVIEW>
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_clouds_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine isccp_clouds_end

!-------------------------------------------------------------------
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine isccp_clouds_end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (Time, axes )
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   initialization time for the netcdf output fields
!  </IN>
! </SUBROUTINE>
subroutine diag_field_init (Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       Time      initialization time for the netcdf output fields
!       axes      diagnostic variable axes
!
!---------------------------------------------------------------------


!-----------------------------------------------------------------------
!    register the isccp diagnostic fields with diag_manager_mod.
!-----------------------------------------------------------------------
       
       id_allclouds     = register_diag_field ( mod_name, &
                        'all_clouds_is', axes(1:2), Time, &
       'Cloud fraction as ISCCP would see it', 'fraction' )       
       id_total         = register_diag_field ( mod_name, &
                       'tot_cld_amt_is', axes(1:2), Time, &
       'Cloud fraction as ISCCP would see it', 'fraction' )       
       id_high          = register_diag_field ( mod_name, &
                      'high_cld_amt_is', axes(1:2), Time, &
               '    pc<440                  ', 'fraction' )       
       id_mid           = register_diag_field ( mod_name, &
                       'mid_cld_amt_is', axes(1:2), Time, &
               '440<pc<680                  ', 'fraction' )
       id_low           = register_diag_field ( mod_name, &
                       'low_cld_amt_is', axes(1:2), Time, &
               '680<pc                      ', 'fraction' )       
       id_deep          = register_diag_field ( mod_name, &
                                 'deep', axes(1:2), Time, &
               '    pc<440;    23<tau       ', 'fraction' )       
       id_cirrostratus  = register_diag_field ( mod_name, &
                                   'cs', axes(1:2), Time, &
               '    pc<440;   3.6<tau<23    ', 'fraction' )
       id_cirrus        = register_diag_field ( mod_name, &
                                   'ci', axes(1:2), Time, &
               '    pc<440;taumin<tau<3.6   ', 'fraction' )
       id_hithin        = register_diag_field ( mod_name, &
                               'hithin', axes(1:2), Time, &
               '    pc<440;     0<tau<taumin', 'fraction' )
       id_nimbostratus  = register_diag_field ( mod_name, &
                                   'ns', axes(1:2), Time, &
               '440<pc<680;    23<tau       ', 'fraction' )
       id_altostratus   = register_diag_field ( mod_name, &
                                   'as', axes(1:2), Time, &
               '440<pc<680;   3.6<tau<23    ', 'fraction' )
       id_altocumulus   = register_diag_field ( mod_name, &
                                   'ac', axes(1:2), Time, &
               '440<pc<680;taumin<tau<3.6   ', 'fraction' )
       id_midthin       = register_diag_field ( mod_name, &
                              'midthin', axes(1:2), Time, &
               '440<pc<680;     0<tau<taumin', 'fraction' )
       id_stratus       = register_diag_field ( mod_name, &
                                   'st', axes(1:2), Time, &
               '680<pc    ;    23<tau       ', 'fraction' )
       id_stratocumulus = register_diag_field ( mod_name, &
                                   'sc', axes(1:2), Time, &
               '680<pc    ;   3.6<tau<23    ', 'fraction' )
       id_cumulus       = register_diag_field ( mod_name, &
                                   'cu', axes(1:2), Time, &
               '680<pc    ;taumin<tau<3.6   ', 'fraction' )
       id_lowthin       = register_diag_field ( mod_name, &
                              'lowthin', axes(1:2), Time, &
               '680<pc    ;     0<tau<taumin', 'fraction' )
       
        id_pc1tau0 = register_diag_field ( mod_name, &
                             'pc1tau0', axes(1:2), Time, &
                      '    pc<180;     0<tau<taumin', 'fraction' )
        id_pc1tau1 = register_diag_field ( mod_name, &
                             'pc1tau1', axes(1:2), Time, &
                      '    pc<180;taumin<tau<1.3   ', 'fraction' )
        id_pc1tau2 = register_diag_field ( mod_name, &
                             'pc1tau2', axes(1:2), Time, &
                      '    pc<180;   1.3<tau<3.6   ', 'fraction' )
        id_pc1tau3 = register_diag_field ( mod_name, &
                             'pc1tau3', axes(1:2), Time, &
                     '    pc<180;   3.6<tau<9.4   ', 'fraction' )
        id_pc1tau4 = register_diag_field ( mod_name, &
                             'pc1tau4', axes(1:2), Time, &
                      '    pc<180;   9.4<tau<23    ', 'fraction' )
        id_pc1tau5 = register_diag_field ( mod_name, &
                             'pc1tau5', axes(1:2), Time, &
                      '    pc<180;    23<tau<60    ', 'fraction' )
        id_pc1tau6 = register_diag_field ( mod_name, &
                             'pc1tau6', axes(1:2), Time, &
                      '    pc<180;    60<tau       ', 'fraction' )
        id_pc2tau0 = register_diag_field ( mod_name, &
                             'pc2tau0', axes(1:2), Time, &
                      '180<pc<310;     0<tau<taumin', 'fraction' )
        id_pc2tau1 = register_diag_field ( mod_name, &
                             'pc2tau1', axes(1:2), Time, &
                      '180<pc<310;taumin<tau<1.3   ', 'fraction' )
        id_pc2tau2 = register_diag_field ( mod_name, &
                             'pc2tau2', axes(1:2), Time, &
                      '180<pc<310;   1.3<tau<3.6   ', 'fraction' )
        id_pc2tau3 = register_diag_field ( mod_name, &
                             'pc2tau3', axes(1:2), Time, &
                      '180<pc<310;   3.6<tau<9.4   ', 'fraction' )
        id_pc2tau4 = register_diag_field ( mod_name, &
                             'pc2tau4', axes(1:2), Time, &
                      '180<pc<310;   9.4<tau<23    ', 'fraction' )
        id_pc2tau5 = register_diag_field ( mod_name, &
                             'pc2tau5', axes(1:2), Time, &
                      '180<pc<310;    23<tau<60    ', 'fraction' )
        id_pc2tau6 = register_diag_field ( mod_name, &
                             'pc2tau6', axes(1:2), Time, &
                      '180<pc<310;    60<tau       ', 'fraction' )
        id_pc3tau0 = register_diag_field ( mod_name, &
                             'pc3tau0', axes(1:2), Time, &
                      '310<pc<440;     0<tau<taumin', 'fraction' )
        id_pc3tau1 = register_diag_field ( mod_name, &
                             'pc3tau1', axes(1:2), Time, &
                      '310<pc<440;taumin<tau<1.3   ', 'fraction' )
        id_pc3tau2 = register_diag_field ( mod_name, &
                             'pc3tau2', axes(1:2), Time, &
                      '310<pc<440;   1.3<tau<3.6   ', 'fraction' )
        id_pc3tau3 = register_diag_field ( mod_name, &
                             'pc3tau3', axes(1:2), Time, &
                      '310<pc<440;   3.6<tau<9.4   ', 'fraction' )
        id_pc3tau4 = register_diag_field ( mod_name, &
                             'pc3tau4', axes(1:2), Time, &
                      '310<pc<440;   9.4<tau<23    ', 'fraction' )
        id_pc3tau5 = register_diag_field ( mod_name, &
                             'pc3tau5', axes(1:2), Time, &
                      '310<pc<440;    23<tau<60    ', 'fraction' )
        id_pc3tau6 = register_diag_field ( mod_name, &
                             'pc3tau6', axes(1:2), Time, &
                      '310<pc<440;    60<tau       ', 'fraction' )
        id_pc4tau0 = register_diag_field ( mod_name, &
                             'pc4tau0', axes(1:2), Time, &
                      '440<pc<560;     0<tau<taumin', 'fraction' )
        id_pc4tau1 = register_diag_field ( mod_name, &
                             'pc4tau1', axes(1:2), Time, &
                      '440<pc<560;taumin<tau<1.3   ', 'fraction' )
        id_pc4tau2 = register_diag_field ( mod_name, &
                             'pc4tau2', axes(1:2), Time, &
                      '440<pc<560;   1.3<tau<3.6   ', 'fraction' )
        id_pc4tau3 = register_diag_field ( mod_name, &
                             'pc4tau3', axes(1:2), Time, &
                      '440<pc<560;   3.6<tau<9.4   ', 'fraction' )
        id_pc4tau4 = register_diag_field ( mod_name, &
                             'pc4tau4', axes(1:2), Time, &
                      '440<pc<560;   9.4<tau<23    ', 'fraction' )
        id_pc4tau5 = register_diag_field ( mod_name, &
                             'pc4tau5', axes(1:2), Time, &
                      '440<pc<560;    23<tau<60    ', 'fraction' )
        id_pc4tau6 = register_diag_field ( mod_name, &
                             'pc4tau6', axes(1:2), Time, &
                      '440<pc<560;    60<tau       ', 'fraction' )
        id_pc5tau0 = register_diag_field ( mod_name, &
                             'pc5tau0', axes(1:2), Time, &
                      '560<pc<680;     0<tau<taumin', 'fraction' )
        id_pc5tau1 = register_diag_field ( mod_name, &
                             'pc5tau1', axes(1:2), Time, &
                      '560<pc<680;taumin<tau<1.3   ', 'fraction' )
        id_pc5tau2 = register_diag_field ( mod_name, &
                             'pc5tau2', axes(1:2), Time, &
                      '560<pc<680;   1.3<tau<3.6   ', 'fraction' )
        id_pc5tau3 = register_diag_field ( mod_name, &
                             'pc5tau3', axes(1:2), Time, &
                      '560<pc<680;   3.6<tau<9.4   ', 'fraction' )
        id_pc5tau4 = register_diag_field ( mod_name, &
                             'pc5tau4', axes(1:2), Time, &
                      '560<pc<680;   9.4<tau<23    ', 'fraction' )
        id_pc5tau5 = register_diag_field ( mod_name, &
                             'pc5tau5', axes(1:2), Time, &
                      '560<pc<680;    23<tau<60    ', 'fraction' )
        id_pc5tau6 = register_diag_field ( mod_name, &
                             'pc5tau6', axes(1:2), Time, &
                      '560<pc<680;    60<tau       ', 'fraction' )
        id_pc6tau0 = register_diag_field ( mod_name, &
                             'pc6tau0', axes(1:2), Time, &
                      '680<pc<800;     0<tau<taumin', 'fraction' )
        id_pc6tau1 = register_diag_field ( mod_name, &
                             'pc6tau1', axes(1:2), Time, &
                     '680<pc<800;taumin<tau<1.3   ', 'fraction' )
        id_pc6tau2 = register_diag_field ( mod_name, &
                             'pc6tau2', axes(1:2), Time, &
                      '680<pc<800;   1.3<tau<3.6   ', 'fraction' )
        id_pc6tau3 = register_diag_field ( mod_name, &
                             'pc6tau3', axes(1:2), Time, &
                      '680<pc<800;   3.6<tau<9.4   ', 'fraction' )
        id_pc6tau4 = register_diag_field ( mod_name, &
                             'pc6tau4', axes(1:2), Time, &
                      '680<pc<800;   9.4<tau<23    ', 'fraction' )
        id_pc6tau5 = register_diag_field ( mod_name, &
                             'pc6tau5', axes(1:2), Time, &
                      '680<pc<800;    23<tau<60    ', 'fraction' )
        id_pc6tau6 = register_diag_field ( mod_name, &
                             'pc6tau6', axes(1:2), Time, &
                      '680<pc<800;    60<tau       ', 'fraction' )
        id_pc7tau0 = register_diag_field ( mod_name, &
                             'pc7tau0', axes(1:2), Time, &
                      '800<pc    ;     0<tau<taumin', 'fraction' )
        id_pc7tau1 = register_diag_field ( mod_name, &
                             'pc7tau1', axes(1:2), Time, &
                      '800<pc    ;taumin<tau<1.3   ', 'fraction' )
        id_pc7tau2 = register_diag_field ( mod_name, &
                             'pc7tau2', axes(1:2), Time, &
                      '800<pc    ;   1.3<tau<3.6   ', 'fraction' )
        id_pc7tau3 = register_diag_field ( mod_name, &
                             'pc7tau3', axes(1:2), Time, &
                      '800<pc    ;   3.6<tau<9.4   ', 'fraction' )
        id_pc7tau4 = register_diag_field ( mod_name, &
                             'pc7tau4', axes(1:2), Time, &
                      '800<pc    ;   9.4<tau<23    ', 'fraction' )
        id_pc7tau5 = register_diag_field ( mod_name, &
                             'pc7tau5', axes(1:2), Time, &
                      '800<pc    ;    23<tau<60    ', 'fraction' )
        id_pc7tau6 = register_diag_field ( mod_name, &
                             'pc7tau6', axes(1:2), Time, &
                      '800<pc    ;    60<tau       ', 'fraction' )
        id_nisccp = register_diag_field ( mod_name, &
                             'nisccp', axes(1:2), Time, &
                      'frequency of ISCCP calculations', 'fraction' )
       
        id_ninhomog = register_diag_field ( mod_name, &
                             'ninhomog', axes(1:2), Time, &
                      'frequency of cloud inhomogeneity calculations', &
                       'fraction' )
        
        id_inhomogeneity = register_diag_field ( mod_name, &
                      'inhomog_param', axes(1:2), Time, &
                      'Ratio of logarthmic to linear mean optical thickness', &
                       'fraction' )
       
        
!---------------------------------------------------------------------

 
end subroutine diag_field_init




!######################################################################
! <FUNCTION NAME="ran0">
!  <OVERVIEW>
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!  </OVERVIEW>
!  <DESCRIPTION>
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!  </DESCRIPTION>
!  <TEMPLATE>
!   x = ran0 (idum)
!  </TEMPLATE>
!  <IN NAME="idum" TYPE="real">
!   seed for random number generator
!  </IN>
! </FUNCTION>
function ran0 (idum)

!--------------------------------------------------------------------- 
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!---------------------------------------------------------------------
       
real                    :: ran0               
integer, intent(inout)  :: idum
                                 
!--------------------------------------------------------------------
!  intent(out) variable:
!
!      ran0           random number generated by this function
!
!  intent(inout) variable:
!
!      idum           seed for random number generator
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer :: ia = 16807        ! constant in random number generator
      integer :: im = 2147483647   ! constant in random number generator
      integer :: iq = 127773       ! constant in random number generator
      integer :: ir = 2836         ! constant in random number generator
      real    :: am                ! constant in random number generator
      integer :: k                 ! work variable              

!---------------------------------------------------------------------
!    define a needed  variable.
!---------------------------------------------------------------------
      am = 1.0/im

!---------------------------------------------------------------------
!    verify that the seed is valid.
!---------------------------------------------------------------------
      if (idum == 0) then
        call error_mesg ('isccp_clouds_mod', &
         'ZERO seed not allowed in ran0', FATAL)
      endif
 
!---------------------------------------------------------------------
!    compute next random number in sequence, using input value of
!    idum. return a new value of idum for use on next call.
!---------------------------------------------------------------------
      k = idum/iq
      idum = ia*(idum - k*iq) - ir*k
      if (idum < 0) idum = idum + im
      ran0 = am*idum

!---------------------------------------------------------------------


end function ran0

!######################################################################
! Pincus additions start here
! 
!######################################################################
  subroutine scops(cc, conv, seed, frac_out)
    real,    dimension(:, :), &  ! Dimensions nPoints, nLev
      intent( in) :: cc, &       !  Cloud cover in each model level (fraction) 
                                 !    NOTE:  This is the HORIZONTAL area of each grid box covered by clouds
                     conv        !  Convective cloud cover in each model level (fraction) 
                                 !    NOTE:  This is the HORIZONTAL area of each grid box covered by convective clouds
    integer, dimension(:), & ! Dimensions nPoints
      intent(inout) :: seed        !  seed value for random number generator ( see Numerical Recipes Chapter 7)
                                 !  It is recommended that the seed is set to a different value for each model
                                 !  gridbox it is called on, as it is possible that the choice of the same 
                                 !  seed value every time may introduce some statistical bias in the results, 
                                 !  particularly for low values of NCOL.

    real, dimension(:, :, :), &  ! Dimensions nPoints, nCol, nLev
      intent(out) :: frac_out    ! boxes gridbox divided up into
                                 ! Equivalent of BOX in original version, but indexed by column then row, 
                                 ! rather than by row then column
!     -----
!     Internal variables 
!     -----
    integer :: nPoints, nCol, nLev
!     -----
    integer :: ilev, ibox
    real, dimension(size(frac_out, 1), &
                  0:size(frac_out, 3)) :: tca ! total cloud cover in each model level (fraction)
                                              ! with extra layer of zeroes on top
                                              ! in this version this just contains the values input
                                              ! from cc but with an extra level
    real, dimension(size(frac_out, 1), &
                    size(frac_out, 3)) :: cca ! convective cloud cover in each model level (fraction) from conv 

    real, dimension(size(frac_out, 1), &
                    size(frac_out, 2)) :: threshold, &   ! pointer to position in gridbox
                                           maxocc,   &   ! Flag for max overlapped conv cld
                                           maxosc,   &   ! Flag for max overlapped strat cld
                                           boxpos,   &   ! ordered pointer to position in gridbox
                                           threshold_min ! minimum value to define range in with new threshold is chosen
    real, dimension(size(frac_out, 1), &
                    size(frac_out, 2)) :: ran            ! vector of random numbers
      
! ---------------------------------------------------------------------------
      nPoints = size(frac_out, 1)
      nCol    = size(frac_out, 2)
      nLev    = size(frac_out, 3)

!     -----------------------------------------------------!
!  Error checking - not needed in this version because we reset these variables in the calling routine.
!     ---------------------------------------------------!
!      call check_bounds(cc,   "cloud fraction",             0., 1.)
!      call check_bounds(conv, "convective cloud fraction",  0., 1.)

!     ---------------------------------------------------!

!     assign 2d tca array using 1d input array cc
      tca(:, 0 )   = 0.
      tca(:, 1:)   = cc(:, :)
 !     assign 2d cca array using 1d input array conv
      cca(:, :)    = conv(:, :)
      boxpos(:, :) = spread((/ (ibox - 0.5, ibox = 1, nCol) /) , &
                             dim = 1, nCopies = nPoints) / real(nCol)

!     ---------------------------------------------------!
!     Initialise output variable
!     ---------------------------------------------------!

!   Initialised frac_out to zero
      frac_out(:, :, :) = 0.
   
!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO ilev = 1,nlev
!     Initialise threshold
        IF (ilev == 1) then
          ! If max overlap 
          IF (overlap == 1) then
            ! select pixels spread evenly across the gridbox
            threshold(:, :) = boxpos(:, :)
          ELSE
            DO ibox=1, ncol
              call ran0_vec(seed, ran(:, ibox))
              ! select random pixels from the non-convective
              ! part the gridbox ( some will be converted into
              ! convective pixels below )
              threshold(:,ibox) = cca(:,ilev) + (1 - cca(:,ilev)) * ran(:, ibox)
            end do
          ENDIF
        end if

        ! All versions of overlap
        where(boxpos(:, :) <= spread(cca(:, ilev), dim = 2, nCopies = nCol))
          maxocc(:, :) = 1
        elsewhere
          maxocc(:, :) = 0
        end where
        
        !
        ! Apply the overlap assumption
        !
        select case(overlap)
          case(1) ! Max overlap
            threshold_min(:, :) = spread(cca(:,ilev), dim = 2, nCopies = nCol)
            maxosc(:, :) = 1
          case(2)! Random overlap
            threshold_min(:, :) = spread(cca(:,ilev), dim = 2, nCopies = nCol)
            maxosc(:, :) = 0
          case(3)  ! Max/Random overlap
            threshold_min(:, :) = spread(max(cca(:, ilev), min(tca(:, ilev-1), tca(:, ilev))), &
                                         dim = 2, nCopies = nCol)
            where (threshold(:, :) < spread(min(tca(:, ilev-1), tca(:, ilev)), dim = 2, nCopies = nCol) .and. &
                   threshold(:, :) > spread(cca(:, ilev),                      dim = 2, nCopies = nCol)) 
              maxosc(:, :)= 1
            elsewhere
              maxosc(:, :)= 0
            end where
        end select
  
        ! Reset threshold 
        DO ibox=1,ncol
          call ran0_vec(seed, ran(:, ibox))
        end do
        threshold(:, :) = &
            !if max overlapped conv cloud     
            (    maxocc(:, :)) * ( boxpos(:, :)                                                ) +   &
            ! else
                !if max overlapped strat cloud; threshold=boxpos                       
            (1 - maxocc(:, :)) * ( (    maxosc(:, :)) * (threshold(:, :)                       ) +   &
                !else threshold_min=random[thrmin,1] 
                                   (1 - maxosc(:, :)) * (threshold_min(:, :) + &
                                                        (1 - threshold_min(:, :)) * ran(:, :) )    &
                                 )

        ! Fill frac_out with 1's where tca is greater than the threshold
        where (spread(tca(:,ilev), dim = 2, nCopies = nCol) > threshold(:, :)) 
          frac_out(:, :, ilev) = 1
        elsewhere
          frac_out(:, :,ilev) = 0
        end where               

       ! Partition boxes into stratiform and convective parts 
        where (spread(cca(:,ilev), dim = 2, nCopies = nCol) > threshold(:, :)) 
          frac_out(:, :, ilev) = 2
        end where               

    end do    !loop over nlev
      
  end subroutine scops
! -------------------------------------------------------------------
  subroutine icarus(dtau, pfull,                              & ! Required 
                    dem, dem_wv, at, skt, emsfc_lw, iTrop,    & ! Optional
                    boxtau, boxptop)
    !
    ! Required input arguments
    !
    real,    dimension(:, :, :), & ! Dimensions nPoints, nCol, nLev
      intent( in)  :: dtau         !  mean 0.67 micron optical depth of clouds in each model level
                                   !  NOTE:  this the cloud optical depth of only the cloudy part of the grid box, 
                                   !  it is not weighted with the 0 cloud optical depth of the clear part of the grid box
    real,    dimension(:, :),    &
      intent( in) :: pfull         !  pressure of full model levels (Pascals)
                                   !  pfull(npoints,1)    is    top level of model
                                   !  pfull(npoints,nlev) is bottom level of model
    !
    ! Optional input arguments - for computing radiative cloud-top pressure 
    !   All variables except itrop are needed if top_height == 1 .or. top_height == 3
    !
    real,    dimension(:, :, :), optional, & ! Dimensions nPoints, nCol, nLev
     intent( in) ::  dem         !  10.5 micron longwave emissivity of clouds in each model level.  
                                 !  Same note applies as in dtau.
   real,   dimension(:, :),      optional, & ! Dimensions nPoints, nLev
     intent( in) :: dem_wv, &    ! Water vapor emissivity.  
                    at           ! Air temperature in each model level (K)
   real,   dimension(:),         optional, &
     intent( in) :: skt, &       !  skin Temperature (K)
                    emsfc_lw     !  10.5 micron emissivity of surface (fraction)
   ! Users may supply their own tropopause height indicator
   integer, dimension(:),        optional, & ! Dimension nPoints
     intent( in) :: itrop        ! Index of the tropopause location in each column
   
    !     ------
    !     Output
    !     ------
    real, dimension(:, :), &   !  dimension nPoints, nCol
      intent(out) :: boxtau,         &   !  optical thickness in each column
                     boxptop             !  cloud top pressure (mb) in each sub-column
                  
    !     --------------------------------------------
    !     Local variables and parameters
    integer :: nPoints, nCol, nLev 
    !     ------
    integer                           ::  i, j, ilev, ibox, icycle
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: tau, pTop, tb
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2), &
                       size(dtau, 3)) :: dem_local
    
    ! Variables for adjusting cloud top pressure based on TOA IR radiance
    real,    dimension(size(dtau, 1)) :: fluxtop_clrsky
    real,    dimension(size(dtau, 1), &
                       size(dtau, 2)) :: emcld, fluxtop, tauir, taumin, &
                                         fluxtopinit, transmax, btcmin
    logical, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: intermedTrans
  
    ! Historgram quanitities 
    integer, dimension(size(dtau, 1), &
                       size(dtau, 2)) :: levMatch
    
    ! Tropopause 
    integer, dimension(size(dtau, 1)) :: itrop_local    ! index to tropopause level
    real,    dimension(size(dtau, 1)) :: pTrop, atTrop  ! tropopause pressure, temperature
                                   
    
    !     ------
    !     Local constants
    !     ------
    real,    parameter ::  VisTauToIrTauIce = 1./2.13, VisTauToIrTauLiquid = 1./2.56, &
                           assumedFreezingTb = 260.

    ! -----------------------------------------------------------------------------------
    nPoints = size(dtau, 1); nCol    = size(dtau, 2); nLev    = size(dtau, 3)
    
    ! Error checking - not needed in this version because we ensured these values
    !   made sense in the calling routine.
!   call check_Lbound(dtau, "cloud optical depth", 0.)
!   if(present(dem)) &
!     call check_bounds(dem, "cloud emissivity",  0., 1.)
    
    ! -----------------------------------------------------------------------------------
    ! The easy part - what's the total cloud optical depth in each column? 
    ! 
    tau(:, :) = sum(dtau, dim = 3)

    ! -----------------------------------------------------------------------------------
    ! If we're adjusting cloud-top heights using the VIS and/or IR radiances 
    !    we need to know 1) the TOA radiance, and 2) the position of the tropopause 
    !    (the latter is the pressure level to which really cold clouds are assigned)
    !
    if(top_height == 1 .or. top_height == 3) then
      dem_local(:, :, :) = spread(dem_wv(:, :), dim = 2, nCopies = nCol)

      !
      ! Clear sky calculation - only need one col per GCM grid point
      !
      call computeRadiance(dem_local(:, 1, :), at, skt, emsfc_lw, fluxtop_clrsky(:))
      !
      ! Add contribution to emissivity from cloud
      !
      where(dem(:, :, :) > tiny(dem)) 
        dem_local(:, :, :) = 1. - ( (1. - dem_local(:, :, :)) * (1. -  dem(:, :, :)) )
      end where
      !
      ! And now the all-sky radiance
      !
      call computeRadiance(dem_local, at, skt, emsfc_lw, fluxtop)
    end if
    !     ---------------------------------------------------!
    ! Account for ISCCP procedures to determine cloud top temperature
  
    ! Account for partially transmitting cloud and recompute flux 
    !    ISCCP would see assuming a single layer cloud
    ! Note choice here of VisTauToIrTauIce = 1/2.13, as it is primarily ice
    !    clouds which have partial emissivity and need the adjustment 
    !
    ! If the cloud brightness temperature is greater than 260K,  the liquid cloud conversion
    !   factor is used.
    !
    ! This is discussed on pages 85-87 of the ISCCP D level documentation (Rossow et al. 1996)

    if (top_height == 1 .or. top_height == 3) then
      ! Tropoause height needed if cloud tops are to be adjusted
      !
      if(present(itrop)) then
        itrop_local(:) = itrop(:)
      else 
        call diagnoseTropPressure(pfull, at, itrop_local)
      end if
      do j = 1, nPoints
        ptrop(j)  = pfull(j, itrop_local(j)) 
        attrop(j) =    at(j, itrop_local(j))
      end do

      if (top_height == 1) then
        !compute minimum brightness temperature and optical depth
        btcmin(:, :) = spread(TbToFlux(attrop(:) - 5.), dim = 2, nCopies = nCol)
        
        !note that the initial setting of tauir(:) is needed so that
        !tauir(:) has a realistic value 
        tauir(:, :) = tau(:, :) * VisTauToIrTauIce
        transmax(:, :) = (fluxtop(:,:) - btcmin(:, :)) / &
                         (spread(fluxtop_clrsky(:), dim = 2, nCopies = nCol) - btcmin(:, :))
        taumin(:, :) = -log(max(min(transmax(:, :), 1. - spacing(1.)), tiny(transmax)))

        intermedTrans(:, :) = transmax(:, :) > tiny(transmax) .and. &
                              transmax(:, :) < 1. - spacing(1.)
        where (intermedTrans) 
          fluxtopinit(:, :) = fluxtop(:,:)
          tauir(:, :) = tau(:,:) * VisTauToIrTauIce
        end where
        
!       do icycle=1,2
!         where (tau(:,:) > tiny(tau) .and. intermedTrans) 
!            emcld(:,:) = 1. - exp(-tauir(:, :))
!            fluxtop(:,:) = fluxtopinit(:, :) - ( (1. - emcld(:,:)) * &
!                                                 spread(fluxtop_clrsky(:), dim = 2, nCopies = nCol) )
!            fluxtop(:,:) = max(tiny(fluxtop), fluxtop(:,:)/emcld(:,:))
!            tb(:, :) = fluxToTb(fluxtop(:, :))
!            where (tb(:, :) .gt. assumedFreezingTb) tauir(:, :) = tau(:,:) * VisTauToIrTauLiquid
!         end where
!       enddo

        do icycle=1,2
          do j= 1, size(tau,2)
          do i= 1, size(tau,1)
            if (tau(i,j) > tiny(tau)  .and. intermedTrans(i,j)) then
             emcld(i,j) = 1. - exp(-tauir(i, j))
             fluxtop(i,j) = fluxtopinit(i,j) - ( (1. - emcld(i,j)) * &
                      fluxtop_clrsky(i))
             fluxtop(i,j) = max(tiny(fluxtop), fluxtop(i,j)/emcld(i,j))
!            tb(i,j) = fluxToTb(fluxtop(i ,j))
              tb(i,j) = 1307.27/ log(1.+(1./fluxtop(i,j)))
             if (tb(i,j) .gt. assumedFreezingTb)  then
               tauir(i,j) = tau(i,j) * VisTauToIrTauLiquid
             endif
            endif
           end do
           end do
        enddo
      end if  
      
      where(tau(:, :) > tiny(tau)) 
        tb(:, :) = fluxToTb(fluxtop(:, :))
      elsewhere
        ! Clear sky brightness temperature
        tb(:,:) = spread(fluxToTb(fluxtop_clrsky(:)), dim = 2, nCopies = nCol)
      end where 
      
      if(top_height == 1) then
        ! Adjust the brightness temperature and optical depth of very thin clouds
        where (tau(:, :) > tiny(tau) .and. tauir(:, :) < taumin(:, :)) 
           tb(:, :) = spread(attrop(:) - 5., dim= 2, nCopies = nCol)
          tau(:, :) = taumin(:, :) / VisTauToIrTauIce
        end where
      end if
    end if
    ! -----------------------------------------------------------------------------------
    ! Determine cloud-top pressure. Three choices:
    !     Radiatively determined cloud top pressure (top_height = 1 for VIS/R;  3 for IR)
    !     Physical cloud top pressure (top_height = 2)
    
    if (top_height .eq. 1 .or. top_height .eq. 3) then  
      !segregate according to optical thickness (??? - RP)
      levmatch(:, :) = 0
      do ibox=1,ncol
        ! find lowest (?) level whose temperature
        ! most closely matches brightness temperature
        !
        do ilev = 1, nlev-1
          where((at(:,ilev) >= tb(:,ibox) .and. at(:,ilev+1) < tb(:,ibox)) .or. &
                (at(:,ilev) <= tb(:,ibox) .and. at(:,ilev+1) > tb(:,ibox)) )
            where(abs(at(:,ilev) - tb(:,ibox)) < abs(at(:,ilev+1) - tb(:,ibox)))
              levmatch(:, ibox) = ilev
            elsewhere 
              levmatch(:, ibox) = ilev + 1
            end where
          end where 
        end do

        ! If we've found a matching level use it, otherwise use the boundary value
        !  
        do j = 1, nPoints
          if(levmatch(j, ibox) >= 1) then
            ptop(j, ibox) = pfull(j, levmatch(j, ibox))
          else if (tb(j, ibox) < minval(at(j, :))) then
            levmatch(j, ibox) = itrop_local(j)
            ptop(j, ibox) = ptrop(j)
          else if (tb(j, ibox) > maxval(at(j, :))) then
            levmatch(j, ibox) = nLev
            ptop(j, ibox) = pFull(j, nLev)
          end if
        end do
      end do
      
    else !  When top_height .eq. 2
      ! The highest cloud top (clouds being where tau > 0). 
      ptop(:, :) = 0.
      do ibox = 1, nCol
        do ilev = 1, nlev
          where(ptop(:, ibox) <= 0. .and. dtau(:, ibox, ilev) > 0.)
            ptop(:, ibox) = pfull(:, ilev)
            levmatch(:,ibox) = ilev
          end where
        end do
      end do
    end if                            
          
    ! No matter how cloud top pressure is determined, 
    !   pTop and levmatch are 0 when the column is clear.   
    ! 
    where (tau(:,:) <= tiny(tau)) 
      ptop(:,:) = 0.
      levmatch(:,:) = 0      
    end where
    !     ---------------------------------------------------!
    ! Output
    !    
    boxptop(:, :) = ptop(:, :) / 100. 
    boxtau(:, :)  = tau(:, :)
    
  end subroutine icarus
  ! -------------------------------------------------------------------
  ! ------------------------------------------------------ 
  function computeIsccpJointHistograms(tau, ptop, sunlit) result(isccpJointHistogram)
    ! Dimensions nGridCell, nSubColumn 
    real, dimension(:, :),    intent( in) :: tau, ptop
    ! Dimensions nGridCell 
    integer, dimension(:),    intent( in) :: sunlit
    ! Dimensions nGridCells, nTauLevels (7), nPressureLevels(7)
    real, dimension(size(tau, 1), numIsccpOpticalDepthIntervals, numIsccpPressureIntervals) &
                                          :: isccpJointHistogram
    
    ! Local parameters
    real, dimension(numIsccpPressureIntervals + 1),      parameter ::       &
             isccpPressureBinEdges = (/ tiny(isccp_taumin),                 &
                                        180., 310., 440., 560., 680., 800., &
                                        huge(isccp_taumin) /) 
    real, dimension(numIsccpOpticalDepthIntervals + 1)             ::            &
            isccpOpticalDepthBinEdges ! Set at runtime, since isccp_taumin is variable.

    ! Local variables
    integer                     :: i, j
    logical, dimension(size(tau, 1), size(tau, 2)) &
                                :: box_cloudy
    ! --------------------
    isccpOpticalDepthBinEdges = (/ tiny(isccp_taumin),                    &
                                   isccp_taumin, 1.3, 3.6, 9.4, 23., 60., & 
                                   huge(isccp_taumin) /) 
    box_cloudy(:, :) = tau(:, :) > tiny(tau) .and. ptop(:, :) > 0 
    
    !
    ! Construct the histogram
    !
    do i = 1, numIsccpPressureIntervals 
      do j = 1, numIsccpOpticalDepthIntervals 
        isccpJointHistogram(:, i, j) = count(box_cloudy(:, :) .and.                               &
                                             tau(:, :)  >= isccpOpticalDepthBinEdges(i)     .and. &
                                             tau(:, :)  <  isccpOpticalDepthBinEdges(i + 1) .and. &
                                             pTop(:, :) >= isccpPressureBinEdges(j)         .and. &
                                             pTop(:, :) <  isccpPressureBinEdges(j + 1), dim = 2)
      end do
    end do
    isccpJointHistogram(:, :, :)  = isccpJointHistogram(:, :, :)/size(tau, 2)
  end function computeIsccpJointHistograms
  ! ------------------------------------------------------ 
  ! ------------------------------------------------------ 
  function computeInhomogeneityParameter (tau,  ptop, sunlit)  result(inhomog_number)
    ! Dimensions nGridCell, nSubColumn 
    real, dimension(:, :),    intent( in) :: tau, ptop
    ! Dimensions nGridCell 
    integer, dimension(:),    intent( in) :: sunlit
    ! Dimensions nGridCell
    real, dimension(size(tau, 1))         :: inhomog_number
    
    ! Local variables
    real,   dimension(size(tau, 1))               :: logAve, linearAve, tmp
    logical, dimension(size(tau, 1), size(tau, 2)) :: isCloudy
    integer :: i,j
    !
    ! compute linear and logarithmic averages
    
    isCloudy(:, :) = tau(:,:) > isccp_taumin .and. ptop(:, :) > 0
    
    !
    ! compute inhomogeneity parameter
    
!    where(count(isCloudy(:, :), dim = 2) > minColsInhomo)
!      logAve(:)    = sum(log(tau(:, :)), dim = 2, mask = isCloudy(:, :)) / &
!                     count(isCloudy(:, :), dim = 2)
!      linearAve(:) = sum(tau(:, :),      dim = 2, mask = isCloudy(:, :)) / &
!                     count(isCloudy(:, :), dim = 2)
!      inhomog_number(:) = 1. - ( exp(logAve(:))/linearAve(:) )        
!    elsewhere
!      inhomog_number(:) = -1.
!    endwhere
    tmp = count(isCloudy(:, :), dim = 2)
    inhomog_number(:) = -1.
    do i= 1,size(logave,1)
      if ( tmp(i) > minColsInhomo) then
      logAve(i)    = 0.0
      linearAve(i) = 0.0
        do j = 1, size(tau,2)
          if (isCloudy(i,j) ) then
            logAve(i)    = logAve(i) + log(tau(i,j))
            linearAve(i) = linearAve(i) +  tau(i,j)
          endif
        enddo
      logAve(i)    = logAve(i)   /tmp(i)
      linearAve(i) = linearAve(i)/tmp(i)
      endif
    enddo  
    do i= 1,size(logave,1)
      if ( tmp(i) > minColsInhomo) &
        inhomog_number(i) = 1. - ( exp(logAve(i))/linearAve(i) )
    enddo  
    
  end function computeInhomogeneityParameter
  ! ------------------------------------------------------   
  ! -------------------------------------------------------------------      
  subroutine computeWaterVaporEmissivity(pfull, phalf, qv, at, dem_wv)
     real, dimension(:, :), & ! nPoints, nLev
       intent( in) :: pFull, pHalf, qv, at
     real, dimension(:, :), & ! nPoints, nLev
       intent(out) :: dem_wv
     
     ! Local variables
     integer :: nPoints, nLev
     integer :: iLev
     real, dimension(size(dem_wv, 1)) :: press, dpress, atmden, rvh20, wk, rhoave, rh20s, &
                                         rfrgn, tmpexp, tauwv
    
     real, parameter :: wtmair = 28.9644, wtmh20 = 18.01534, Navo = 6.023E+23, grav = 9.806650E+02, &
                        pstd = 1.013250E+06, t0 = 296.
    ! -------------------------------------------
    nPoints = size(pFull, 1); nLev    = size(pFull, 2)

    !compute water vapor continuum emissivity
    !this treatment follows Schwarkzopf and Ramasamy
    !JGR 1999,vol 104, pages 9467-9499.
    !the emissivity is calculated at a wavenumber of 955 cm-1, 
    !or 10.47 microns 
    do ilev=1,nlev
      ! press and dpress are dyne/cm2 = Pascals *10
      press(:) = pfull(:,ilev)*10.
      dpress(:) = (phalf(:,ilev+1)-phalf(:,ilev))*10
      !atmden = g/cm2 = kg/m2 / 10 
      atmden(:) = dpress(:)/grav
      rvh20(:) = qv(:,ilev)*wtmair/wtmh20
      wk(:) = rvh20(:)*Navo*atmden(:)/wtmair
      rhoave(:) = (press(:)/pstd)*(t0/at(:,ilev))
      rh20s(:) = rvh20(:)*rhoave(:)
      rfrgn(:) = rhoave(:)-rh20s(:)
      tmpexp(:) = exp(-0.02*(at(:,ilev)-t0))
      tauwv(:) = wk(:)*1.e-20*((0.0224697*rh20s(:)*tmpexp(:)) + (3.41817e-7*rfrgn(:)) )*0.98
      dem_wv(:,ilev) = 1. - exp( -1. * tauwv(:))
    end do
  end subroutine computeWaterVaporEmissivity 
  ! -------------------------------------------------------------------
  subroutine diagnoseTropPressure(pfull, at, itrop)
    real,    dimension(:, :), intent( in) :: pFull, at
    integer, dimension(:),    intent(out) :: itrop
    
    integer                         :: nPoints, nLev
    real, dimension(size(pFull, 1)) :: attropmin
    integer                         :: ilev

    nPoints = size(pFull, 1); nLev = size(pFull, 2)
    attropmin(:) = 400.
    itrop(:)     = 1
  
    do  ilev=1,nlev
      where(pfull(:, ilev) < 40000. .and. pfull(:, ilev) > 5000. .and. &
            at(:, ilev) < attropmin(:)) 
        attropmin(:) = at(:, ilev)
        itrop(:)=ilev
      end where
    end do
  end subroutine diagnoseTropPressure 
  ! -------------------------------------------------------------------      
  ! -------------------------------------------------------------------      
  subroutine check_bounds(array, name, minAllowed, maxAllowed)
    implicit none
    ! Input variables
    real, dimension(:, :), intent( in) :: array
    character(len = *),    intent( in) :: name
    real,                  intent( in) :: minAllowed, maxAllowed
    
    ! ---------------------
    if(any(array(:, :) < minAllowed .or. array(:, :) > maxAllowed)) then
      call error_mesg ('isccp_clouds_mod', &
                       'Values in array out of bounds', FATAL)
    end if
  end subroutine check_bounds
  ! -------------------------------------------------------------------      
  subroutine check_Lbound(array, name, minAllowed)
    implicit none
    ! Input variables
    real, dimension(:, :, :), intent( in) :: array
    character(len = *),       intent( in) :: name
    real,                     intent( in) :: minAllowed
    
    ! ---------------------
    if(any(array(:, :, :) < minAllowed )) then
      call error_mesg ('isccp_clouds_mod', &
                       'Array values smaller than lower bound', FATAL)
    end if
  end subroutine check_lBound
  ! -------------------------------------------------------------------      
  ! Functions to compute brightness temperature given a 11 micron radiance
  !
  ! -------------------------------------------------------------------      
  function fluxToTb_1D(flux) result(tb)
    real, dimension(:),       intent( in) :: flux
    real, dimension(size(flux))           :: tb
    
    tb(:) = 1307.27 / log(1. + (1./flux(:)))
  end function fluxToTb_1D
  !     ---------------------------------------------------!
  function fluxToTb_2D(flux) result(tb)
    real, dimension(:, :),    intent( in) :: flux
    real, dimension(size(flux, 1), &
                    size(flux, 2))        :: tb
    
    tb(:, :) = 1307.27 / log(1. + (1./flux(:, :)))
  end function fluxToTb_2D
  !     ---------------------------------------------------!
  function fluxToTb_3D(flux) result(tb)
    real, dimension(:, :, :), intent( in) :: flux
    real, dimension(size(flux, 1), &
                    size(flux, 2), &
                    size(flux, 3))        :: tb
    
    tb(:, :, :) = 1307.27 / log(1. + (1./flux(:, :, :)))
  end function fluxToTb_3D
  ! -------------------------------------------------------------------      
  ! Functions to compute 11 micron radiance given a brightness temperature
  !
  ! -------------------------------------------------------------------      
  function TbToFlux_1D(tb) result(flux)
    real, dimension(:),       intent( in) :: tb
    real, dimension(size(tb))             :: flux
    
    flux(:) = 1. / ( exp(1307.27/tb(:)) - 1. )
  end function TbToFlux_1D
  !     ---------------------------------------------------!
  function TbToFlux_2D(tb) result(flux)
    real, dimension(:, :),    intent( in) :: tb
    real, dimension(size(tb, 1), &  
                    size(tb, 2))          :: flux
    
    flux(:, :) = 1. / ( exp(1307.27/tb(:, :)) - 1. )
  end function TbToFlux_2D
  ! ---------------------------------------------------!
  function TbToFlux_3D(tb) result(flux)
    real, dimension(:, :, :), intent( in) :: tb
    real, dimension(size(tb, 1), &
                    size(tb, 2), &
                    size(tb, 3))          :: flux
    
    flux(:, :, :) = 1. / ( exp(1307.27/tb(:, :, :)) - 1. )
  end function TbToFlux_3D
  ! -------------------------------------------------------------------      
  ! -------------------------------------------------------------------      
  subroutine computeRadiance_1D(dem, at, skt, emsfc_lw, TOAradiance)
    real,    dimension(:, :), &  ! Dimensions nPoints, nLev
      intent( in) :: dem,     &  !   10.5 micron emissivity of water vapor
                     at          !   air temperature
    real,   dimension(:),     &  ! Dimension nPoints
      intent( in) :: skt,     &  !   skin Temperature (K)
                     emsfc_lw    !   10.5 micron emissivity of surface (fraction) 
    real, dimension(:),       &  ! Dimension nPoint, nCol
      intent(out) :: TOAradiance !   10.5 micron nadir radiance at TOA
    
  !     ------
  !     Local variables and parameters
    integer :: nPoints, nLev
    !     ------
    integer ::  ilev
    real,    dimension(size(dem, 1)) :: trans_layers_above
    real,    dimension(size(dem, 1)) :: bb
                                                          
    !----------------------------------------------------------------------
    ! Computes radiance at TOA from an emitting/absorbing atmosphere
    !     TOAradiance is the 10.5 micron radiance at the top of the
    !              atmosphere
    !     trans_layers_above is the total transmissivity in the layers
    !             above the current layer
    !----------------------------------------------------------------------
  
    !initialize variables
    nPoints = size(dem, 1); nLev    = size(dem, 2)
    TOAradiance(:) = 0.; trans_layers_above(:) = 1.
  
    do ilev=1,nlev
      ! Black body emission at temperature of the layer
      bb(:) = TbToFlux(at(:,ilev))
  
      ! increase TOA flux by flux emitted from layer
      ! times total transmittance in layers above
      TOAradiance(:) = TOAradiance(:) + dem(:, ilev) * trans_layers_above(:) * bb(:)
        
      ! update trans_layers_above with transmissivity
      ! from this layer for next time around loop
      trans_layers_above(:) = trans_layers_above(:) * (1. - dem(:, ilev))
    enddo ! ilev
  
    !surface emission
    bb(:) = TbToFlux(skt(:)) 
  
    !add in surface emission
    TOAradiance(:) = TOAradiance(:) + trans_layers_above(:) * emsfc_lw(:) * bb(:)
                       
  end subroutine computeRadiance_1D
 ! -------------------------------------------------------------------      
  subroutine computeRadiance_2D(dem, at, skt, emsfc_lw, TOAradiance)
    real,    dimension(:, :, :), &  ! Dimensions nPoints, nCol, nLev
      intent( in) :: dem            !   10.5 micron emissivity of water vapor
   real,    dimension(:, :),     &  ! Dimensions nPoints, nLev
      intent( in) :: at             !   air temperature
    real,   dimension(:),        &  ! Dimension nPoints
      intent( in) :: skt,        &  !   skin Temperature (K)
                     emsfc_lw       !   10.5 micron emissivity of surface (fraction) 
    real, dimension(:, :),       &  ! Dimension nPoint, nCol
      intent(out) :: TOAradiance    !   10.5 micron nadir radiance at TOA
    
  !     ------
  !     Local variables and parameters
    integer :: nPoints, nCol, nLev 
    !     ------
    integer ::  ilev
    real,    dimension(size(dem, 1), &
                       size(dem, 2)) :: trans_layers_above
    real,    dimension(size(dem, 1)) :: bb
 
    !----------------------------------------------------------------------
    ! Computes radiance at TOA from an emitting/absorbing atmosphere
    !     TOAradiance is the 10.5 micron radiance at the top of the
    !              atmosphere
    !     trans_layers_above is the total transmissivity in the layers
    !             above the current layer
    !----------------------------------------------------------------------
  
    !initialize variables
    nPoints = size(dem, 1); nCol    = size(dem, 2); nLev    = size(dem, 3)
    TOAradiance(:, :) = 0.; trans_layers_above(:, :) = 1.
  
    do ilev=1,nlev
      ! Black body emission at temperature of the layer
      bb(:) = TbToFlux(at(:,ilev))
  
      ! increase TOA flux by flux emitted from layer
      ! times total transmittance in layers above
      TOAradiance(:,:) = TOAradiance(:,:) + &
                         dem(:,:, ilev) * trans_layers_above(:,:) * spread(bb(:), dim = 2, nCopies = nCol) 
        
      ! update trans_layers_above with transmissivity
      ! from this layer for next time around loop
      trans_layers_above(:, :) = trans_layers_above(:, :) * (1. - dem(:, :, ilev))
    enddo ! ilev
  
    !surface emission
    bb(:) = TbToFlux(skt(:)) 
  
    !add in surface emission
    TOAradiance(:,:) = TOAradiance(:,:) +  &
                       trans_layers_above(:,:) * spread(emsfc_lw(:) * bb(:), dim = 2, nCopies = nCol)
  
  end subroutine computeRadiance_2D
 ! -------------------------------------------------------------------      
  subroutine ran0_vec(idum, ran0)
    integer, dimension(:), intent(inout) :: idum
    real,    dimension(:), intent(  out) :: ran0
!     $Id: isccp_clouds.F90,v 19.0 2012/01/06 20:16:57 fms Exp $
!     Platform independent random number generator from
!     Numerical Recipies
!     Mark Webb July 1999
!     Fortran 90 implementation Robert Pincus, Dec 2003
    
    integer, parameter :: IA = 16807, IM = 2147483647,  IQ = 127773, IR = 2836
    real,    parameter :: AM = 1.0/IM
    
    integer, dimension(size(idum)) :: k

    if(any(idum(:) == 0))                  &
      call error_mesg ('isccp_clouds_mod', &
                       'idum: ZERO seed not allowed', FATAL)

    k(:) = idum(:) / IQ
    idum (:) = IA * (idum(:) - k(:) * IQ) - IR * k(:)
    where (idum(:) < 0) idum(:) = idum(:) + IM
    ran0(:) = AM * idum(:)

  end subroutine ran0_vec
!#####################################################################

                end module isccp_clouds_mod



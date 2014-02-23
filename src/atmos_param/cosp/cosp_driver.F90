#include "cosp_defs.H"

module cosp_driver_mod

! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: open_namelist_file, open_file,  &
                                    close_file, error_mesg, FATAL, &
                                    file_exist, mpp_pe, mpp_root_pe,   &
                                    check_nml_error, write_version_number, &
                                    stdlog
use sat_vapor_pres_mod,       only: compute_qs
use time_manager_mod,         only: set_date, time_type, operator (+), &
                                    operator(-), operator(<),    &
                                    operator(>), operator(<=), &
                                    operator(>=),  get_date, print_date, &
                                    get_calendar_type, NOLEAP, &
                                    assignment(=), set_time
USE MOD_COSP_TYPES,           only: cosp_config, cosp_gridbox,    &
                                    cosp_subgrid, cosp_sgradar,   &
                                    cosp_sglidar, cosp_isccp, &
#ifdef RTTOV
                                    cosp_rttov, &
#endif
                                    cosp_vgrid, cosp_radarstats,  &
                                    cosp_lidarstats, &
                                    cosp_sghydro,  cosp_misr, &
                                    construct_cosp_gridbox,  &
                                    construct_cosp_misr,  &
                                    construct_cosp_vgrid,  &
                                    construct_cosp_subgrid, &
                                    construct_cosp_sghydro, &
                                    construct_cosp_sgradar, &
                                    construct_cosp_radarstats, &
                                    construct_cosp_sglidar, &
                                    construct_cosp_lidarstats, &
                                    construct_cosp_isccp, &           
#ifdef RTTOV
                                    construct_cosp_rttov, &           
                                    free_cosp_rttov, &           
#endif
                                    free_cosp_gridbox,  &
                                    free_cosp_misr,  &
                                    free_cosp_vgrid,  &
                                    free_cosp_subgrid, &
                                    free_cosp_sghydro, &
                                    free_cosp_sgradar, &
                                    free_cosp_radarstats, &
                                    free_cosp_sglidar, &
                                    free_cosp_lidarstats, &
                                    free_cosp_isccp
USE MOD_COSP,                 only: cosp
USE MOD_COSP_IO,              only: read_cosp_output_nl,  &
                                    map_ll_to_point, map_point_to_ll
                       
use MOD_COSP_CONSTANTS,       only: PARASOL_NREFL, I_LSCLIQ, I_LSCICE, &
                                    I_CVCLIQ, I_CVCICE, I_LSGRPL, &
                                    I_LSRAIN, I_LSSNOW, I_CVRAIN, I_CVSNOW,&
                                    N_HYDRO, RTTOV_MAX_CHANNELS
use MOD_COSP_Modis_Simulator, only: COSP_MODIS, FREE_COSP_MODIS,  &
                                    CONSTRUCT_COSP_MODIS
use cosp_diagnostics_mod,     only: cosp_diagnostics_init,   &
                                    cosp_diagnostics_time_vary, &
                                    output_cosp_fields, &
                                    cosp_diagnostics_endts, &
                                    cosp_diagnostics_end
use mod_cosp_utils,           only: flip_vert_index

IMPLICIT NONE

public cosp_driver, cosp_driver_init, cosp_driver_end, cosp_driver_endts, &
       cosp_driver_time_vary

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cosp_driver.F90,v 20.0 2013/12/13 23:15:41 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

!---------------------------------------------------------------------
!namelist variables

character(len=512) :: finput ! Input file name, not used in FMS
character(len=512) :: cmor_nl = '   '  ! not used in FMS
integer :: overlap = 3   !  overlap type: 1=max, 2=rand, 3=max/rand
integer :: isccp_topheight = 1  
                    ! 1 = adjust top height using both a computed
                    !     infrared brightness temperature and the visible
                    !     optical depth to adjust cloud top pressure. Note
                    !     that this calculation is most appropriate to 
                    !     compare to ISCCP data during sunlit hours.
                    ! 2 = do not adjust top height, that is cloud top
                    !     pressure is the actual cloud top pressure
                    !     in the model
                    ! 3 = adjust top height using only the computed
                    !     infrared brightness temperature. Note that this
                    !     calculation is most appropriate to compare to 
                    !     ISCCP IR only algortihm (i.e. you can compare to 
                    !     nighttime ISCCP data with this option)

integer :: isccp_topheight_direction = 2   
                    !     direction for finding atmosphere pressure level
                    !     with interpolated temperature equal to the 
                    !     radiance determined cloud-top temperature
                    ! 1 = find the *lowest* altitude (highest pressure) 
                    !     level with interpolated temperature equal to the 
                    !     radiance determined cloud-top temperature
                    ! 2 = find the *highest* altitude (lowest pressure) 
                    !     level with interpolated temperature equal to the 
                    !     radiance determined cloud-top temperature
                    !     ONLY APPLICABLE IF top_height EQUALS 1 or 3

integer :: Nlr = 40 !     Number of levels in statistical outputs
integer :: Npoints_it = 20000  
                    !     Max number of gridpoints to be processed in 
                    !     one iteration
real :: radar_freq = 94.       ! CloudSat radar frequency (GHz)
real :: k2= -1.                ! |K|^2, -1=use frequency dependent default
integer :: surface_radar = 0   !  surface=1, spaceborne=0 
integer :: use_mie_tables = 0  ! use a precomputed lookup table? yes=1,no=0
integer :: use_gas_abs = 1     ! include gaseous absorption? yes=1,no=0
integer :: do_ray = 0          ! calculate/output Rayleigh refl=1, not=0
integer :: melt_lay = 0        ! melting layer model off=0, on=1
integer :: Nprmts_max_hydro = 12 
                               ! Max number of parameters for hydrometeor 
                               ! size distributions
integer :: Naero = 1           ! Number of aerosol species (Not used)
integer :: Nprmts_max_aero = 1 ! Max number of parameters for aerosol 
                               ! size distributions (Not used)
integer :: lidar_ice_type = 0  ! Ice particle shape in lidar calculations 
                               ! (0=ice-spheres ; 1=ice-non-spherical)
logical :: use_vgrid = .true.  ! Use fixed vertical grid for outputs? 
                               ! (if .true. then you need to define number 
                               ! of levels with Nlr)
logical :: csat_vgrid =.true.  ! CloudSat vertical? 
                               ! (if .true. then the CloudSat standard grid
                               ! is used for the outputs.
logical :: use_precipitation_fluxes =.true. 
                               ! True if precipitation fluxes are input 
                               ! to the algorithm
logical :: use_reff =.true.    ! True if you want effective radius to be 
                               ! used by radar simulator (always used 
                               ! by lidar)
logical :: use_input_file = .false.
logical :: produce_cmor_output_fields = .false.
real    :: emsfc_lw_nml=0.94
logical :: use_rh_wrt_liq = .true.
!-------------------------------------------------------------------------
!-------------- RTTOV inputs
!-------------------------------------------------------------------------
integer :: platform = 1    ! satellite platform
integer :: satellite = 15  ! satellite
integer :: Instrument = 0  ! instrument
integer :: Nchannels = 8   ! Number of channels to be computed
real :: ZenAng = 50.       ! Satellite Zenith Angle
real :: co2 = 5.241e-04    ! mixing ratio of trace gas
real :: ch4 = 9.139e-07    ! mixing ratio of trace gas
real :: n2o = 4.665e-07    ! mixing ratio of trace gas
real :: co = 2.098e-07     ! mixing ratio of trace gas
integer,dimension(RTTOV_MAX_CHANNELS) :: Channels = 0
                           ! Channel numbers (please be sure that 
                           ! you supply Nchannels)
real,dimension(RTTOV_MAX_CHANNELS) :: Surfem = 0.0 
                           ! Surface emissivity (please be sure that 
                           ! you supply Nchannels)

namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight, &
                    isccp_topheight_direction, &
                    use_vgrid,nlr,csat_vgrid,  &
                    npoints_it,finput, &
                    radar_freq,surface_radar,use_mie_tables, &
                    use_input_file, produce_cmor_output_fields, &
                    emsfc_lw_nml, use_rh_wrt_liq, &
                    use_gas_abs,do_ray,melt_lay,k2,Nprmts_max_hydro,  &
                    Naero,Nprmts_max_aero,lidar_ice_type, &
                    use_precipitation_fluxes,use_reff, &
                    platform,satellite,Instrument,Nchannels, &
                    Channels,Surfem,ZenAng,co2,ch4,n2o,co

! Local variables

character(len=64)  :: cosp_output_nl='cosp_output_nl.txt'
integer            :: Ncolumns ! Number of subcolumns in SCOPS
integer            :: Nlevels  ! Number of levels
type(cosp_config)  :: cfg   ! Configuration options
integer            :: geomode
real               :: emsfc_lw
double precision   :: time=1.D0
double precision   :: time_bnds(2)= (0.5D0, 1.5D0)

!---------------- End of declaration of variables --------------



contains

!######################################################################

subroutine cosp_driver_init (lonb, latb, Time_diag, axes, kd_in, ncol_in)

real, dimension(:,:),  intent(in) :: lonb, latb
type(time_type),       intent(in) :: Time_diag
integer, dimension(4), intent(in) :: axes
integer,               intent(in) :: kd_in, ncol_in

   integer :: io, unit, ierr, logunit
   integer :: imax, jmax

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=cosp_input, iostat=io)
    ierr = check_nml_error(io,"cosp_input")
#else
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
    if ( file_exist('input.nml')) then
       unit =  open_namelist_file ()
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=cosp_input, iostat=io, end=10)
      ierr = check_nml_error(io,'cosp_input')
      enddo
10    call close_file (unit)
    endif
#endif
        
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
    call write_version_number (version, tagname)
    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe() )    &
                        write (logunit, nml=cosp_input)

    if (use_mie_tables /= 0) then
      call error_mesg ('cosp_driver', &
            'use_mie_tables must be set to 0 currently', FATAL)
    endif

    nlevels = kd_in
    ncolumns = ncol_in 
    imax = size(lonb,1) - 1
    jmax = size(lonb,2) - 1

    call read_cosp_output_nl(cosp_output_nl,cfg)

    call cosp_diagnostics_init      &
            (imax, jmax, Time_diag, axes, nlevels, ncolumns, cfg, &
             use_vgrid, csat_vgrid, nlr)     

!---------------------------------------------------------------------
!   COSP takes a single, spacially independent value for surface
!   emissivity. it may be supplied via namelist.
!---------------------------------------------------------------------
    emsfc_lw = emsfc_lw_nml
 
!--------------------------------------------------------------------
!   variable geomode indicates that the grid (i,j) => (lon,lat)
!--------------------------------------------------------------------
    geomode = 2
 

end subroutine cosp_driver_init


!######################################################################

subroutine cosp_driver_time_vary (Time_diag)

type(time_type), intent(in)  :: Time_diag

    call cosp_diagnostics_time_vary (Time_diag)

end subroutine cosp_driver_time_vary


!######################################################################

subroutine cosp_driver_endts

    call cosp_diagnostics_endts

end subroutine cosp_driver_endts


!####################################################################

subroutine cosp_driver   &
        (lat_in, lon_in, daytime_in, phalf_plus, p_full_in, zhalf_plus,&
         z_full_in, u_wind_in, v_wind_in, mr_ozone_in, &
         T_in, sh_in, tca_in, cca_in, lsliq_in, lsice_in, ccliq_in, &
         ccice_in, fl_lsrain_in, fl_lssnow_in, fl_lsgrpl_in, &
         fl_ccrain_in,  &
         fl_ccsnow_in, reff_lsclliq_in, reff_lsclice_in,   &
         reff_lsprliq_in, reff_lsprice_in, reff_ccclliq_in,  &
         reff_ccclice_in, reff_ccprliq_in, reff_ccprice_in,  &
         skt_in, land_in, Time_diag, is, js, stoch_mr_liq_in, &
         stoch_mr_ice_in, stoch_size_liq_in, stoch_size_frz_in, &
         tau_stoch_in, lwem_stoch_in, stoch_cloud_type_in)
!--------------------------------------------------------------------
!    subroutine cosp_driver is the interface between the cosp simulator 
!    code and the AM model.
!--------------------------------------------------------------------
real, dimension(:,:),   intent(in) :: lat_in, lon_in, skt_in, land_in, &
                                      u_wind_in, v_wind_in
real, dimension(:,:), intent(in) :: daytime_in
real, dimension(:,:,:), intent(in) :: phalf_plus, p_full_in, &
        zhalf_plus, z_full_in, T_in, sh_in, &
        tca_in, cca_in, lsliq_in, lsice_in, ccliq_in, ccice_in, &
        fl_lsrain_in, fl_lssnow_in, fl_lsgrpl_in, fl_ccrain_in, &
        fl_ccsnow_in, mr_ozone_in, &
        reff_lsclliq_in, reff_lsclice_in, reff_lsprliq_in, &
        reff_lsprice_in, reff_ccclliq_in, reff_ccclice_in, &
        reff_ccprliq_in, reff_ccprice_in
real, dimension(:,:,:,:), intent(in), optional ::  &
               tau_stoch_in, lwem_stoch_in, stoch_cloud_type_in, &
               stoch_mr_liq_in, stoch_mr_ice_in, stoch_size_liq_in, &
               stoch_size_frz_in
type(time_type), intent(in) :: Time_diag

!local variables:

  integer, intent(in) :: is, js
  real, dimension(size(T_in,1)*size(T_in,2), size(T_in,3), &
                              ncolumns)  :: y3, y3a, y4, y5, y6, y7, y8
  integer :: i, j, n, l
  integer :: k

  type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
  type(cosp_subgrid) :: sgx     ! Subgrid outputs
  type(cosp_sghydro) :: sghydro ! Subgrid condensate
  type(cosp_sgradar) :: sgradar ! Output from radar simulator
  type(cosp_sglidar) :: sglidar ! Output from lidar simulator
  type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
  type(cosp_modis)   :: modis   ! Output from MODIS simulator
  type(cosp_misr)    :: misr    ! Output from MISR simulator
#ifdef RTTOV 
  type(cosp_rttov)   :: rttov   ! Output from RTTOV 
#endif 
  type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  real, dimension ( size(T_in,1)*size(T_in,2), size(T_in,3)) :: &
                    fl_lsrain, &
                    fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
  real, dimension ( size(T_in,1)*size(T_in,2),    &
                                  Ncolumns, size(T_in,3)) :: &
                                                          cloud_type
  real, dimension ( size(T_in,1),size(T_in,2), size(T_in,3)) :: &
                                                p_half_in, z_half_in
  integer :: nlon,nlat,npoints

  !---------------- End of declaration of variables --------------
   
      nlon = size(T_in,1)
      nlat = size(T_in,2)
      npoints = nlon*nlat
 
 
      p_half_in(:,:,1:size(T_in,3)) = phalf_plus(:,:,2:size(T_in,3)+1)
      z_half_in(:,:,1:size(T_in,3)) = zhalf_plus(:,:,2:size(T_in,3)+1)
      if (present (tau_stoch_in)         .and. &
          present (lwem_stoch_in)        .and. &
          present (stoch_cloud_type_in)  .and. &
          present (stoch_mr_liq_in)      .and. &
          present (stoch_mr_ice_in)      .and. &
          present (stoch_size_liq_in)    .and. &
          present (stoch_size_frz_in) ) then
        sgx%cols_input_from_model = .true.
      else
        sgx%cols_input_from_model = .false.
      endif


      call construct_cosp_gridbox    &
            (time, time_bnds, radar_freq, surface_radar, use_mie_tables,  &
             use_gas_abs, do_ray, melt_lay, k2, Npoints, Nlevels,   &
             Ncolumns, N_HYDRO, Nprmts_max_hydro, Naero, Nprmts_max_aero, &
             Npoints_it, lidar_ice_type, isccp_topheight,  &
             isccp_topheight_direction,overlap, emsfc_lw,  &
             use_precipitation_fluxes, use_reff, &
             Platform, Satellite, Instrument, Nchannels, ZenAng, &
             channels(1:Nchannels), surfem(1:Nchannels), co2, ch4, n2o, co,&
             gbx)
  
      call produce_cosp_input_fields ( Npoints, Nlevels, N_hydro,  &
              lon_in, lat_in, daytime_in, p_half_in, p_full_in, z_half_in, &
              z_full_in, u_wind_in, v_wind_in, mr_ozone_in, T_in, &
              sh_in, tca_in, &
              cca_in, lsliq_in, lsice_in, ccliq_in, ccice_in,  &
              fl_lsrain_in,  &
              fl_lssnow_in, fl_lsgrpl_in, fl_ccrain_in, fl_ccsnow_in, &
              reff_lsclliq_in, reff_lsclice_in, reff_lsprliq_in, &
              reff_lsprice_in, reff_ccclliq_in, reff_ccclice_in,  &
              reff_ccprliq_in, reff_ccprice_in, tau_stoch_in,  &
              lwem_stoch_in, stoch_cloud_type_in, skt_in, land_in, &
              cloud_type, gbx)  
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Define new vertical grid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)
        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Allocate memory for other types
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)

      call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,N_hydro,sghydro)

      if (sgx%cols_input_from_model) then
!---------------------------------------------------------------------
!    convert the stochastic column inputs from lon-lat to npoints, then
!    save the column values, reversing the vertical indices, into the
!    cosp_subgrid_type variables.
!---------------------------------------------------------------------
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= stoch_mr_liq_in, y3= y3  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= stoch_mr_ice_in   (:,:,:,:), y3= y4  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= stoch_size_liq_in, y3= y5  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= stoch_size_frz_in, y3= y6  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= tau_stoch_in, y3= y7  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= lwem_stoch_in, y3= y8  )
        call map_ll_to_point(nlon,nlat,npoints, &
                        x4= stoch_cloud_type_in(:,:,:,:), y3= y3a  )
        do l=1, NCOLUMNS
          do k=1,nlevels
            do n=1,npoints
              if (y3a(n,k,l) == 1.0) then
                sghydro%mr_hydro(n,l,nlevels-k+1,I_LSCLIQ) = y3(n,k,l)
                sghydro%mr_hydro(n,l,nlevels+1-k,I_LSCICE) = y4(n,k,l)
                if ( sghydro%mr_hydro(n,l,nlevels-k+1,I_LSCLIQ) > 0.0) then
                  sghydro%Reff(n,l,nlevels-k+1,I_LSCLIQ) = y5(n,k,l)
                else
                  sghydro%Reff(n,l,nlevels+1-k,I_LSCLIQ) = 0.0          
                endif
                if (sghydro%mr_hydro(n,l,nlevels+1-k,I_LSCICE) > 0.0) then
                  sghydro%Reff(n,l,nlevels+1-k,I_LSCICE) = y6(n,k,l)
                else
                  sghydro%Reff(n,l,nlevels+1-k,I_LSCICE) = 0.0          
                endif
              else
                sghydro%mr_hydro(n,l,nlevels+1-k,I_LSCLIQ) = 0.0          
                sghydro%mr_hydro(n,l,nlevels+1-k,I_LSCICE) = 0.0          
                sghydro%Reff(n,l,nlevels+1-k,I_LSCLIQ) = 0.0          
                sghydro%Reff(n,l,nlevels+1-k,I_LSCICE) = 0.0          
              endif
              if (y3a(n,k,l) == 2.0) then
                sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCLIQ) = y3(n,k,l)
                sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCICE) = y4(n,k,l)
                if (sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCLIQ) > 0.0) then
                  sghydro%Reff(n,l,nlevels+1-k,I_CVCLIQ) = y5(n,k,l)
                else
                  sghydro%Reff(n,l,nlevels+1-k,I_CVCLIQ) = 0.0          
                endif
                if (sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCICE) > 0.0) then
                  sghydro%Reff(n,l,nlevels+1-k,I_CVCICE) = y6(n,k,l)
                else
                  sghydro%Reff(n,l,nlevels+1-k,I_CVCICE) = 0.0          
                endif
              else
                sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCLIQ) = 0.0          
                sghydro%mr_hydro(n,l,nlevels+1-k,I_CVCICE) = 0.0          
                sghydro%Reff(n,l,nlevels+1-k,I_CVCLIQ) = 0.0          
                sghydro%Reff(n,l,nlevels+1-k,I_CVCICE) = 0.0          
              endif
              sgx%dtau_col(n,l,k) = y7(n,k,l)
              sgx%dem_col(n,l,k) = y8(n,k,l)
            end do
          end do
        end do
      endif
      call construct_cosp_sgradar  &
                              (cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
      call construct_cosp_radarstats   &
                        (cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)
      call construct_cosp_sglidar   &
                (cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
      call construct_cosp_lidarstats  &
          (cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
      call construct_cosp_isccp  &
                               (cfg,Npoints,Ncolumns,Nlevels,isccp)
      call construct_cosp_modis    &
                                          (cfg,Npoints,Ncolumns,modis)
      call construct_cosp_misr   &
                               (cfg,Npoints,misr)
#ifdef RTTOV 
      call construct_cosp_rttov   &
                                 (Npoints,Nchannels,rttov) 
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call simulator
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (Ncolumns == 1) then
        if (gbx%use_precipitation_fluxes) then
          call error_mesg ('cosp_driver:cosp_driver', &
              'Use of precipitation fluxes not supported in&
                                & CRM mode (Ncolumns=1)', FATAL)
        endif
        if ((maxval(gbx%dtau_c) > 0.0).or.(maxval(gbx%dem_c) > 0.0)) then
          call error_mesg ('cosp_driver:cosp_driver', &
              ' dtau_c > 0.0 or dem_c > 0.0. In CRM mode (Ncolumns=1) &
             &the optical depth (emmisivity) of all clouds must be &
             &passed through dtau_s (dem_s)', FATAL)
        endif
      endif

!-------------------------------------------------------------------------
! save the grid-box mean of these fields for diagnostic output. The gbx%
! arrays will be changed to sub-column values in subroutine cosp.
!------------------------------------------------------------------------
      fl_lsrain = gbx%rain_ls
      fl_lssnow = gbx%snow_ls
      fl_lsgrpl = gbx%grpl_ls
      fl_ccrain = gbx%rain_cv
      fl_ccsnow = gbx%snow_cv
#ifdef RTTOV 
      call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,  &
                isccp,misr,modis,rttov,stradar,stlidar, sghydro, cloud_type)
#else
      call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,  &
                      isccp,misr,modis,stradar,stlidar, sghydro, cloud_type)
#endif

!-------------------------------------------------------------------------
!    define the effective size of precipitation particles in the gridbox. 
!    each subcolumn with precip type will have the same size.
!-------------------------------------------------------------------------
      gbx%Reff(:,:,I_LSRAIN) = 0.
      gbx%Reff(:,:,I_LSSNOW) = 0.
      gbx%Reff(:,:,I_LSGRPL) = 0.
      gbx%Reff(:,:,I_CVRAIN) = 0.
      gbx%Reff(:,:,I_CVSNOW) = 0.
      do l=1,Ncolumns
        gbx%Reff(:,:,I_LSRAIN) =    &
                  Max(sghydro%reff(:,l,:,I_LSRAIN), gbx%Reff(:,:,I_LSRAIN))
        gbx%Reff(:,:,I_LSSNOW) =    &
                  Max(sghydro%reff(:,l,:,I_LSSNOW), gbx%Reff(:,:,I_LSSNOW))
        gbx%Reff(:,:,I_LSGRPL) =    &
                  Max(sghydro%reff(:,l,:,I_LSGRPL), gbx%Reff(:,:,I_LSGRPL))
        gbx%Reff(:,:,I_CVRAIN) =    &
                  Max(sghydro%reff(:,l,:,I_CVRAIN), gbx%Reff(:,:,I_CVRAIN))
        gbx%Reff(:,:,I_CVSNOW) =    &
                   Max(sghydro%reff(:,l,:,I_CVSNOW), gbx%Reff(:,:,I_CVSNOW))
      end do

!-------------------------------------------------------------------------
!  return these fields to grid-box mean values for diagnostic 
!  output purposes.
!-------------------------------------------------------------------------
      gbx%rain_ls = fl_lsrain
      gbx%snow_ls = fl_lssnow
      gbx%grpl_ls = fl_lsgrpl
      gbx%rain_cv = fl_ccrain
      gbx%snow_cv = fl_ccsnow

!-------------------------------------------------------------------------
!  call output_cosp_fields to produce netcdf output of desired COSP fields.
!-------------------------------------------------------------------------
      call output_cosp_fields (nlon, nlat, npoints, geomode, stlidar, &
                               stradar, isccp, modis, misr, sgradar, &
                               sglidar, sgx, Time_diag, is, js, &
                               cloud_type, gbx, cfg, phalf_plus, zhalf_plus)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Deallocate memory in derived types
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call free_cosp_gridbox(gbx)
      call free_cosp_subgrid(sgx)
      call free_cosp_sghydro(sghydro)
      call free_cosp_sgradar(sgradar)
      call free_cosp_radarstats(stradar)
      call free_cosp_sglidar(sglidar)
      call free_cosp_lidarstats(stlidar)
      call free_cosp_isccp(isccp)
      call free_cosp_misr(misr)
      call free_cosp_modis(modis)
#ifdef RTTOV 
      call free_cosp_rttov(rttov) 
#endif
      call free_cosp_vgrid(vgrid)  

 
end subroutine cosp_driver

!#####################################################################

subroutine cosp_driver_end 

      call cosp_diagnostics_end 


end subroutine cosp_driver_end

!#####################################################################

subroutine produce_cosp_input_fields   &
        (Npnts, Nl, N_hydro, lon_in, lat_in, daytime_in, p_half_in,   &
         p_full_in, z_half_in, z_full_in, u_wind_in, v_wind_in,    &
         mr_ozone_in, T_in, sh_in, tca_in, cca_in, lsliq_in, &
         lsice_in, ccliq_in, ccice_in, fl_lsrain_in, fl_lssnow_in, &
         fl_lsgrpl_in, fl_ccrain_in, fl_ccsnow_in, reff_lsclliq_in,   &
         reff_lsclice_in, reff_lsprliq_in, reff_lsprice_in,    &
         reff_ccclliq_in, reff_ccclice_in, reff_ccprliq_in,    &
         reff_ccprice_in, tau_stoch_in, lwem_stoch_in,    &
         stoch_cloud_type_in, skt_in, land_in, cloud_type, gbx) 

!--------------------------------------------------------------------
!    subroutine produce_cosp_input_fields converts inputs from AM3 
!    to the form needed by COSP.
!--------------------------------------------------------------------

integer,                  intent(in) :: Npnts, Nl, N_hydro
real,dimension(:,:),      intent(in) :: lon_in,lat_in, skt_in, land_in,&
                                        u_wind_in, v_wind_in
real, dimension(:,:),  intent(in) :: daytime_in
real,dimension(:,:,:),    intent(in) ::    &
            p_half_in, p_full_in, z_half_in, z_full_in, T_in, sh_in,  &
            tca_in, cca_in, lsliq_in, lsice_in, ccliq_in, ccice_in, &
            fl_lsrain_in, fl_lssnow_in, fl_lsgrpl_in, fl_ccrain_in, &
            fl_ccsnow_in, mr_ozone_in, &
            reff_lsclliq_in, reff_lsclice_in, reff_lsprliq_in, &
            reff_lsprice_in, reff_ccclliq_in, reff_ccclice_in, &
            reff_ccprliq_in, reff_ccprice_in
real,dimension(:,:,:,:),intent(in)   :: tau_stoch_in, lwem_stoch_in, &
                                        stoch_cloud_type_in
real,dimension(Npnts,Ncolumns,Nl),  intent(out) :: cloud_type
type(cosp_gridbox), intent(out) :: gbx ! Gridbox information. Input for COSP

!  local variables:

     real,dimension(Npnts,Nl)             :: y2                         
     real,dimension(Npnts,Nl,Ncolumns)    :: y3,y3a                   
     real, dimension(Npnts,Nl) :: qs
     real, dimension(size(T_in,1), size(T_in,2), size(T_in,3)) :: &
              qs_in, tau_stoch_mean, lwem_stoch_mean, tau_s_in, &
              tau_c_in, lwem_s_in, lwem_c_in
     integer :: nxdir, nydir, npts
     integer :: n, i, j, k
     real :: sum_s1, sum_s2, sum_c1, sum_c2
     integer :: ctr_s, ctr_c 


!--------------------------------------------------------------------
!   define array dimensions; verify consistency.
!--------------------------------------------------------------------
     nxdir = size(lat_in,1)
     nydir = size(lat_in, 2)
     npts = nxdir*nydir
     if (npts /= Npnts) then
       call error_mesg ('cosp_driver/produce_cosp_input_fields', &
                                     'ERROR -- i*j /= npts', FATAL)
     endif
   
!---------------------------------------------------------------------
!   map the 2d lon-lat arrays to 1D (npoints).
!---------------------------------------------------------------------
   call map_ll_to_point(nxdir,nydir,npts,x2=daytime_in, y1=gbx%sunlit)
   call map_ll_to_point(nxdir,nydir,npts,x2=lat_in, y1=gbx%latitude)
   call map_ll_to_point(nxdir,nydir,npts,x2=lon_in, y1=gbx%longitude)
   call map_ll_to_point(nxdir,nydir,npts,x2=skt_in, y1=gbx%skt)
   call map_ll_to_point(nxdir,nydir,npts,x2=land_in, y1=gbx%land)
   call map_ll_to_point(nxdir,nydir,npts,x2=u_wind_in, y1=gbx%u_wind)
   call map_ll_to_point(nxdir,nydir,npts,x2=v_wind_in, y1=gbx%v_wind)

!---------------------------------------------------------------------
!   map the 3d lon-lat-k arrays to 2D (npoints,k), and flip their
!   vertical indices (index 1 nearest ground in COSP).
!---------------------------------------------------------------------
   call map_ll_to_point(nxdir,nydir,npts,x3=p_full_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%p  )

   call map_ll_to_point(nxdir,nydir,npts,x3=p_half_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%ph )

   gbx%psfc = gbx%ph(:,1)

   call map_ll_to_point(nxdir,nydir,npts,x3=z_full_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%zlev  )

   call map_ll_to_point(nxdir,nydir,npts,x3=z_half_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%zlev_half )

   call map_ll_to_point(nxdir,nydir,npts,x3=mr_ozone_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%mr_ozone )

   call map_ll_to_point(nxdir,nydir,npts,x3=T_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%T  )

   call map_ll_to_point(nxdir,nydir,npts,x3=sh_in, y2=y2)
   call flip_vert_index    (y2, nl,gbx%sh )

!---------------------------------------------------------------------
!   define surface height
!---------------------------------------------------------------------
   gbx%sfc_height(:) = gbx%zlev_half(:,1)

!--------------------------------------------------------------------
!   compute qs and then the relative humidity.
!   a limit my be imposed (nml control) to account for slightly inexact
!   values near saturation.
!--------------------------------------------------------------------
   if (use_rh_wrt_liq) then
     call compute_qs (T_in, p_full_in, qs_in, q=sh_in,  &
                                          es_over_liq = use_rh_wrt_liq)
   else
     call compute_qs (T_in, p_full_in, qs_in, q=sh_in)
   endif
   call map_ll_to_point(nxdir,nydir,npts,x3=qs_in, y2=y2)
   call flip_vert_index    (y2, nl,qs )
   gbx%q = gbx%sh/qs

   call map_ll_to_point(nxdir,nydir,npts,x3=tca_in, y2=y2 )
   call flip_vert_index    (y2, nl,gbx%tca)

   call map_ll_to_point(nxdir,nydir,npts,x3=cca_in, y2=y2 )
   call flip_vert_index    (y2, nl,gbx%cca)

   call map_ll_to_point(nxdir,nydir,npts,x3=lsliq_in, y2=y2      )
   call flip_vert_index    (y2, nl,gbx%mr_hydro(:,:,I_LSCLIQ) )

   call map_ll_to_point(nxdir,nydir,npts,x3=lsice_in, y2=y2      )
   call flip_vert_index    (y2, nl,gbx%mr_hydro(:,:,I_LSCICE) )

   call map_ll_to_point(nxdir,nydir,npts,x3=ccliq_in, y2=y2      )
   call flip_vert_index    (y2, nl,gbx%mr_hydro(:,:,I_CVCLIQ) )

   call map_ll_to_point(nxdir,nydir,npts,x3=ccice_in, y2=y2      )
   call flip_vert_index    (y2, nl,gbx%mr_hydro(:,:,I_CVCICE) )

   call map_ll_to_point(nxdir,nydir,npts,x3=fl_lsrain_in, y2=y2       )
   call flip_vert_index    (y2, nl,gbx%rain_ls)

   call map_ll_to_point(nxdir,nydir,npts,x3=fl_lssnow_in, y2=y2       )
   call flip_vert_index    (y2, nl,gbx%snow_ls)

   call map_ll_to_point(nxdir,nydir,npts,x3=fl_lsgrpl_in, y2=y2       )
   call flip_vert_index    (y2, nl,gbx%grpl_ls)

   call map_ll_to_point(nxdir,nydir,npts,x3=fl_ccrain_in, y2=y2       )
   call flip_vert_index    (y2, nl,gbx%rain_cv)

   call map_ll_to_point(nxdir,nydir,npts,x3=fl_ccsnow_in, y2=y2       )
   call flip_vert_index    (y2, nl,gbx%snow_cv)

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_lsclliq_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_lscliq ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_lsclice_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_lscice ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_lsprliq_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_lsrain ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_lsprice_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_lssnow ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_ccclliq_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_cvcliq ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_ccclice_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_cvcice ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_ccprliq_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_cvrain ))

   call map_ll_to_point(nxdir,nydir,npts,x3=reff_ccprice_in, y2=y2   )
   call flip_vert_index    (y2, nl,gbx%reff(:,:,i_cvsnow ))

   gbx%reff(:,:,i_lsgrpl) = 0.0

!---------------------------------------------------------------------
!   the values of tau and lwem are passed in for each stochastic column.
!   here grid box mean values are obtained for the convective and
!   large-scale components
!---------------------------------------------------------------------

   do k=1, size(tau_stoch_in,3)
     do j=1, size(tau_stoch_in,2)
       do i=1, size(tau_stoch_in,1)
         ctr_s = 0
         ctr_c = 0
         sum_s1 = 0.
         sum_c1 = 0.
         sum_s2 = 0.
         sum_c2 = 0.
         do n=1, size(tau_stoch_in,4)
           if (stoch_cloud_type_in(i,j,k,n) == 1. ) then 
             ctr_s = ctr_s + 1
             sum_s1 = sum_s1 +  tau_stoch_in(i,j,k,n)
             sum_s2 = sum_s2 +  lwem_stoch_in(i,j,k,n)
           else if(stoch_cloud_type_in(i,j,k,n) == 2. ) then 
             ctr_c = ctr_c + 1
             sum_c1 = sum_c1 +  tau_stoch_in(i,j,k,n)
             sum_c2 = sum_c2 +  lwem_stoch_in(i,j,k,n)
           endif
         end do
         if (ctr_s > 0) then
           tau_s_in(i,j,k) = sum_s1/ctr_s
           lwem_s_in(i,j,k) = sum_s2/ctr_s
         else
           tau_s_in(i,j,k) = 0.             
           lwem_s_in(i,j,k) = 0.               
         endif
         if (ctr_c > 0) then
           tau_c_in(i,j,k) = sum_c1/ctr_c
           lwem_c_in(i,j,k) = sum_c2/ctr_c
         else
           tau_c_in(i,j,k) = 0.             
           lwem_c_in(i,j,k) = 0.               
         endif
       end do
     end do
   end do
       
   call map_ll_to_point(nxdir,nydir,npts,x3=tau_s_in(:,:,:), y2=y2)
   call flip_vert_index    (y2, nl,gbx%dtau_s )

   call map_ll_to_point(nxdir,nydir,npts,x3=tau_c_in(:,:,:), y2=y2)
   call flip_vert_index    (y2, nl,gbx%dtau_c )

   call map_ll_to_point(nxdir,nydir,npts,x3=lwem_s_in(:,:,:), y2=y2)
   call flip_vert_index    (y2, nl,gbx%dem_s)

   call map_ll_to_point(nxdir,nydir,npts,x3=lwem_c_in(:,:,:), y2=y2)
   call flip_vert_index    (y2, nl,gbx%dem_c)

!----------------------------------------------------------------------
!    stoch_cloud_type is not flipped here; it will be used in subroutine
!    cosp where it is needed with index 1 being TOA. however the
!    column and vertical indices do need to be reversed.
!----------------------------------------------------------------------
   call map_ll_to_point(nxdir,nydir,npts,  &
                        x4=stoch_cloud_type_in(:,:,:,:), y3=y3        )
   do j=1,nl
     do i=1,Ncolumns
       cloud_type(:,i,j) = y3 (:,j,i)
     end do
   end do
   
!########################################################################

end subroutine produce_cosp_input_fields



    
!#####################################################################



end module cosp_driver_mod



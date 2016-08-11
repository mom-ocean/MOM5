#include "cosp_defs.H"

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



module cosp_diagnostics_mod

use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: open_namelist_file, open_file,  &
                                    close_file, error_mesg, FATAL, &
                                    file_exist, mpp_pe, mpp_root_pe,   &
                                    check_nml_error, write_version_number,&
                                    stdlog
use time_manager_mod,         only: set_date, time_type, operator (+), &
                                    operator(-), operator(<),    &
                                    operator(>), operator(<=), &
                                    operator(>=),  get_date, print_date, &
                                    get_calendar_type, NOLEAP, &
                                    assignment(=), set_time
use diag_grid_mod,            only: get_local_indexes2
use diag_manager_mod,         only: register_diag_field, send_data,  &
                                    diag_axis_init, register_static_field
USE MOD_COSP_TYPES,           only: cosp_config, cosp_gridbox,   &
                                    cosp_subgrid, cosp_sgradar,  &
                                    cosp_sglidar, cosp_isccp, &
#ifdef RTTOV
                                    cosp_rttov, &
#endif
                                    cosp_vgrid, cosp_radarstats,  &
                                    cosp_lidarstats, &
                                    cosp_sghydro,  cosp_misr, &
                                    construct_cosp_vgrid,  &
                                    free_cosp_vgrid
USE MOD_COSP_IO,              only: map_point_to_ll            
                       
use MOD_COSP_CONSTANTS,       only: DBZE_BINS,SR_BINS, PARASOL_NREFL,  &
                                    PARASOL_SZA, CFAD_ZE_MIN,    &
                                    CFAD_ZE_WIDTH, &
                                    LIDAR_UNDEF, ISCCP_PC_BNDS, ISCCP_TAU,&
                                    I_LSCLIQ, I_LSCICE, I_CVCLIQ,   &
                                    I_CVCICE, I_LSGRPL, &
                                    I_LSRAIN, I_LSSNOW, I_CVRAIN,   &
                                    I_CVSNOW, &
                                    N_HYDRO, ISCCP_TAU_BNDS,&
                                    RTTOV_MAX_CHANNELS, MISR_N_CTH,  &
                                    MISR_CTH_BNDS
use MOD_LMD_IPSL_STATS,       only: define_srbval
use MOD_COSP_Modis_Simulator, only: COSP_MODIS
use mod_modis_sim,            only: numTauHistogramBins,   &
                                    numPressureHistogramBins, &
                                    tauHistogramBoundaries, &
                                    nominalTauHistogramBoundaries, &
                                    nominalTauHistogramCenters, &
                                    nominalPressureHistogramBoundaries
use mod_cosp_utils,           only: flip_vert_index

IMPLICIT NONE

public cosp_diagnostics_init, output_cosp_fields, cosp_diagnostics_end, &
       cosp_diagnostics_time_vary, cosp_diagnostics_endts

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id $'
character(len=128)  :: tagname =  '$Name $'

!---------------------------------------------------------------------
!namelist variables

logical :: output_p_and_z_by_index = .false.
logical :: generate_orbital_output = .false.
character (len = 128) :: orbital_filename =  '  '
integer, dimension(6) :: sat_begin_time = (/0,0,0,0,0,0/)
integer :: sat_period     = 0  ! [seconds]
integer :: num_sat_periods = 0
integer :: max_sdgs_per_sat_period = 3500

namelist/cosp_diagnostics_nml/ output_p_and_z_by_index, &
                    generate_orbital_output, orbital_filename, &
                    sat_begin_time, sat_period, num_sat_periods, &
                    max_sdgs_per_sat_period

! Local variables

character(len=16)       :: mod_name = 'cosp'

integer, dimension(14)  :: cosp_axes

integer :: id_lat, id_lon, id_p, id_ph, id_z, id_zh, id_T, id_sh, &
           id_u_wind, id_v_wind, id_mr_ozone, &
           id_tot_h2o, &
           id_rh, id_tca, id_cca, id_lsliq, id_lsice, id_ccliq, &
           id_ccice, id_fl_lsrain, id_fl_lssnow, id_fl_lsgrpl, &
           id_fl_ccrain, id_fl_ccsnow, &
           id_reff_lsclliq, id_reff_lsclice, &
           id_reff_lsprliq, id_reff_lsprice, &
           id_reff_ccclliq, id_reff_ccclice, &
           id_reff_ccprliq, id_reff_ccprice, &
           id_reff_lsclliq_cmip, id_reff_ccclliq_cmip, &
           id_lsca_cmip, id_cca_cmip, &
           id_dtau_s, id_dtau_c, id_dem_s, id_dem_c, id_skt, id_land, &
           id_sfcht, id_sunlit
integer :: id_cltcalipso_sat, id_cllcalipso_sat, id_clmcalipso_sat,  &
           id_clhcalipso_sat
integer :: id_cltcalipso, id_cllcalipso, id_clmcalipso, id_clhcalipso, &
           id_cltlidarradar, id_tclisccp, id_ctpisccp, id_tauisccp, &
           id_tbisccp, id_tbclrisccp, &
           id_betamol532, &
           id_albisccp, id_clcalipso, id_clcalipso2, &
           id_clcalipso_sat, id_clcalipso2_sat, &
           id_clcalipso_mdl, id_clcalipso2_mdl, &
           id_boxtauisccp, id_boxptopisccp, id_parasolrefl, &
           id_parasolrefl_sat, &
           id_sampling_sat, id_location_sat, id_lat_sat, id_lon_sat
integer :: id_tclmodis, id_lclmodis, id_iclmodis, id_ttaumodis, &
           id_ltaumodis, id_itaumodis, id_tlogtaumodis, &
           id_llogtaumodis, id_ilogtaumodis, id_lremodis, &
           id_badlremodis, id_badiremodis, &
           id_locldmodis, id_mdcldmodis, id_hicldmodis, &
           id_iremodis, id_ctpmodis, id_lwpmodis, id_iwpmodis
integer, allocatable, dimension(:) :: id_dbze94, id_cloudsatcfad, &
                                      id_cloudsatcfad_sat, &
                                      id_atb532, id_calipsosrcfad, &
                                      id_calipsosrcfad_sat, &
                                      id_cloud_type, id_boxtauisccp_n, &
                                      id_boxptopisccp_n, &
                                      id_taumodis_n, id_ptopmodis_n, &
                                      id_badsizemodis_n, &
                                      id_sizemodis_n, id_phasemodis_n
integer, allocatable, dimension(:) :: id_cloudsatcfad_mdl, &
                                      id_calipsosrcfad_mdl
integer , dimension(7)            :: id_clisccp
integer , dimension(7,7)          :: id_clisccp_n
integer , dimension(MISR_N_CTH)   :: id_misr    
integer , dimension(7,MISR_N_CTH) :: id_misr_n
integer , dimension(numTauHistogramBins, numPressureHistogramBins) ::  &
                                                         id_tauctpmodis_n
integer , dimension(numPressureHistogramBins) :: id_tauctpmodis

real  :: missing_value = -1.0E30

real, dimension(:,:,:), allocatable        :: location   
logical, dimension(:,:,:), allocatable     :: lflag_array
logical, dimension(:,:,:,:), allocatable   :: lflag_array_temp, &
                                              lflag_array_parasol
real, dimension(:,:,:), allocatable        :: flag_array
type(time_type), dimension(:), allocatable :: Time_start, Time_end
integer   :: imax, jmax, nlr, nlevels, ncolumns
integer   :: nsat_time_prev
integer   :: nsat_time
logical   :: use_vgrid, csat_vgrid

!---------------- End of declaration of variables --------------

include 'netcdf.inc'

contains

!######################################################################

subroutine cosp_diagnostics_init     &
             (imax_in, jmax_in, Time, axes, nlevels_in, ncolumns_in, cfg, &
              use_vgrid_in, csat_vgrid_in, nlr_in)       

type(time_type), intent(in) :: Time
integer, dimension(4), intent(in) :: axes
integer, intent(in) :: imax_in, jmax_in           
integer, intent(in) :: nlevels_in, ncolumns_in
logical, intent(in) :: use_vgrid_in
type(cosp_config), intent(in) :: cfg   ! Configuration options
logical, intent(in) :: csat_vgrid_in
integer, intent (in) :: nlr_in

   integer :: io, unit, ierr, logunit

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=cosp_diagnostics_nml, iostat=io)
    ierr = check_nml_error(io,"cosp_diagnostics_nml")
#else
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
    if ( file_exist('input.nml')) then
       unit =  open_namelist_file ()
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=cosp_diagnostics_nml, iostat=io, end=10)
      ierr = check_nml_error(io,'cosp_diagnostics_nml')
      enddo
10    call close_file (unit)
    endif
#endif
        
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
    call write_version_number(version, tagname)
    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe() )    &
                        write (logunit, nml=cosp_diagnostics_nml)

!----------------------------------------------------------------------
!    save i and j dimensions.
!----------------------------------------------------------------------
    imax = imax_in
    jmax = jmax_in
    nlr = nlr_in
    nlevels = nlevels_in
    ncolumns = ncolumns_in
    use_vgrid = use_vgrid_in
    csat_vgrid = csat_vgrid_in

    if (generate_orbital_output) then
      if (sat_begin_time(1) == 0 .or. sat_begin_time(2) == 0 .or. &
          sat_begin_time(3) ==0) then
        call error_mesg ('cosp_diagnostics_init', &
           'requesting orbital output but not supplying &
                                               &valid start time', FATAL)
      endif
      if (sat_period == 0) then
        call error_mesg ('cosp_diagnostics_init', &
           'satellite sampling period [seconds] must be non-zero', FATAL)
      endif
      if (num_sat_periods == 0) then
        call error_mesg ('cosp_diagnostics_init', &
         'must define number of satellite periods to be processed', FATAL)
      endif
      if (trim(orbital_filename) == '') then
        call error_mesg ('cosp_diagnostics_init', &
              'filename for orbital specification not provided', FATAL)
      endif
    endif

    call diag_field_init (Time, axes, cfg)

    if (generate_orbital_output) then
      allocate (location    (imax,jmax, 1:num_sat_periods))
      allocate (lflag_array (imax,jmax, 0:num_sat_periods))
      allocate (lflag_array_temp (imax,jmax, nlr, 0:num_sat_periods))
      allocate (lflag_array_parasol   &
                            (imax,jmax, PARASOL_NREFL, 0:num_sat_periods))
      allocate (flag_array(imax,jmax,12))
      allocate (Time_start(num_sat_periods))
      allocate (Time_end  (num_sat_periods))
      call read_cloudsat_orbit  
      nsat_time_prev = 1
    endif

end subroutine cosp_diagnostics_init 

!#####################################################################

subroutine diag_field_init (Time, axes, cfg)

type(time_type), intent(in) :: Time
integer, dimension(4), intent(in) :: axes
type(cosp_config), intent(in) :: cfg   ! Configuration options

   real :: column_ax(Ncolumns)
   real :: level_ax(Nlevels )
   real :: isccp_ax(7)           
   real :: modis_ax(numTauHistogramBins)
   real :: dbze_ax(DBZE_BINS)
   real :: lidar_ax(SR_BINS)
   real :: sratio_bounds(2, SR_BINS)
   real :: srbval(SR_BINS)
   real :: csat_ax(NLR)
   real :: month_ax(12)
   real :: hr_ax(num_sat_periods)
   integer :: parasol_ax(PARASOL_NREFL)
   integer, dimension(3) :: halfindx = (/1,2,4/)
   integer, dimension(3) :: columnindx = (/1,2,5/)
   integer, dimension(3) :: levelindx = (/1,2,11/)
   integer, dimension(3) :: parasolindx = (/1,2,6/)
   integer, dimension(3) :: dbzeindx = (/1,2,7/)
   integer, dimension(3) :: lidarindx = (/1,2,8/)
   integer, dimension(3) :: tauindx = (/1,2,9/)
   integer, dimension(3) :: modistauindx = (/1,2,12/)
   integer, dimension(3) :: csatindx = (/1,2,10/)
   integer, dimension(3) :: samplingindx = (/1,2,13/)
   integer, dimension(3) :: samplingindx2 = (/1,2,14/)
   integer :: i, n, m
   integer :: id_columnindx, id_parasolindx, id_dbzeindx, id_lidarindx
   integer :: id_levelindx
   integer :: id_tauindx
   integer :: id_modistauindx
   integer :: id_csatindx
   integer :: id_monindx
   integer :: id_3hrindx
   character(len=2) :: chvers, chvers4
   character(len=8) :: chvers2, chvers3, chvers5, chvers6
   type(cosp_gridbox) :: gbx_t ! Gridbox information. Input for COSP
   type(cosp_vgrid)   :: vgrid_t   ! Information on vertical grid of stats


!--------------------------------------------------------------------
!    define the varisous axes needed for this data.
!--------------------------------------------------------------------
   cosp_axes(1:4) = axes(1:4)

!--------------------------------------------------------------------
! a level counter:
!--------------------------------------------------------------------
   do i=1,Nlevels 
     level_ax(i) = float(i)
   end do
   id_levelindx = diag_axis_init  ('levelindx', level_ax, &
          'levels', 'n', 'level number', & 
           set_name =  mod_name)
   cosp_axes(11) = id_levelindx

!--------------------------------------------------------------------
! a stochastic column counter:
!--------------------------------------------------------------------
   do i=1,Ncolumns
     column_ax(i) = float(i)
   end do
   id_columnindx = diag_axis_init  ('columnindx', column_ax, &
          'subcol', 'n', 'subcolumn number', & 
           set_name =  mod_name)
   cosp_axes(5) = id_columnindx

!--------------------------------------------------------------------
!  a PARASOL index counter:
!--------------------------------------------------------------------
   id_parasolindx = diag_axis_init  ('parasolindx', PARASOL_SZA, &
          'parasolindx', 'n', 'parasol reflectivity index', & 
           set_name =  mod_name)
   cosp_axes(6) = id_parasolindx

!--------------------------------------------------------------------
!  a radar bin counter:
!--------------------------------------------------------------------
   do i=1,DBZE_BINS
      dbze_ax(i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i-0.5)
   end do
   id_dbzeindx = diag_axis_init  ('dbzeindx', dbze_ax, &
          'dbzeindx', 'n', 'dbze', & 
           set_name =  mod_name)
   cosp_axes(7) = id_dbzeindx

!--------------------------------------------------------------------
!  a lidar bin counter:
!--------------------------------------------------------------------

   call define_srbval (srbval)

   sratio_bounds(1,:) = srbval(:)
   sratio_bounds(2,1:SR_BINS-1) = srbval(2:SR_BINS)
   sratio_bounds(2,SR_BINS) = srbval(SR_BINS) +10.0
   lidar_ax(1:SR_BINS) = (sratio_bounds(1,1:SR_BINS) +    &
                                           sratio_bounds(2,1:SR_BINS))/2.0
   id_lidarindx = diag_axis_init  ('lidarindx', lidar_ax, &
          'lidarindx', 'n', 'lidar scattering', & 
           set_name =  mod_name)
   cosp_axes(8) = id_lidarindx

!--------------------------------------------------------------------
!  an isccp tau bin counter:
!--------------------------------------------------------------------
   isccp_ax = isccp_tau
   id_tauindx = diag_axis_init  ('tauindx', isccp_ax, &
          'tauindx', 'n', 'isccp tau category', & 
           set_name =  mod_name)
   cosp_axes(9) = id_tauindx

!--------------------------------------------------------------------
!  a modis tau bin counter:
!--------------------------------------------------------------------
   modis_ax = nominalTauHistogramCenters
   id_modistauindx = diag_axis_init  ('modistauindx', modis_ax, &
          'modistauindx', 'n', 'modis tau category', &
           set_name =  mod_name)
   cosp_axes(12) = id_modistauindx

!--------------------------------------------------------------------
!  a specified vertical index needed when use_vgrid = .true. 
!--------------------------------------------------------------------
   gbx_t%Npoints = 256       
   gbx_t%Ncolumns = ncolumns    
   gbx_t%Nlevels = Nlevels
   allocate(gbx_t%zlev(256    , nlevels))
   allocate(gbx_t%zlev_half(256    , nlevels))
   gbx_t%zlev = 0.0
   gbx_t%zlev_half = 0.0
   call construct_cosp_vgrid(gbx_t,Nlr,use_vgrid,csat_vgrid,vgrid_t)
   csat_ax = vgrid_t%z
   id_csatindx = diag_axis_init  ('csatindx', csat_ax, &
          'csatindx', 'z', 'csat vert index', & 
           set_name =  mod_name)
   cosp_axes(10) = id_csatindx
   deallocate (gbx_t%zlev, gbx_t%zlev_half) 
   call free_cosp_vgrid (vgrid_t)
   do i=1,12
     month_ax(i) = i
   end do
   id_monindx = diag_axis_init  ('samplingindx', month_ax, &
          'samplingindx', 'n', 'month index', & 
           set_name =  mod_name)
   cosp_axes(13) = id_monindx
   
   do i=1,num_sat_periods
     hr_ax(i) = i
   end do
   id_3hrindx = diag_axis_init  ('samplingindx2', hr_ax, &
          'samplingindx2', 'n', '3hr index', & 
           set_name =  mod_name)
   cosp_axes(14) = id_3hrindx
   
!--------------------------------------------------------------------
!    register input fields with diag_manager.
!--------------------------------------------------------------------
   id_lat        = register_diag_field &
      (mod_name, 'lat', axes(1:2), Time, 'Latitude',  'degrees N')

   id_lon        = register_diag_field &
      (mod_name, 'lon', axes(1:2), Time, 'Longitude',  'degrees E')

   id_u_wind     = register_diag_field &
      (mod_name, 'u_wind', axes(1:2), Time, 'sfc u wind',  'm / s')

   id_v_wind     = register_diag_field &
      (mod_name, 'v_wind', axes(1:2), Time, 'sfc v wind',  'm / s')

   if (output_p_and_z_by_index) then
     id_p          = register_diag_field &
       (mod_name, 'p', cosp_axes(levelindx), Time,  &
                                        'P at full levels',  'Pa  ')
     id_ph         = register_diag_field &
       (mod_name, 'ph', cosp_axes(levelindx), Time, &
                                        'p at half levels',  'Pa')
     id_z        = register_diag_field &
       (mod_name, 'z', cosp_axes(levelindx), Time, 'height  ', 'meters')
     id_zh        = register_diag_field &
       (mod_name, 'zh', cosp_axes(levelindx), Time, &
                                      'height at half levs',  'meters')
   else
     id_p          = register_diag_field &
      (mod_name, 'p', axes(1:3), Time, 'P at full levels',  'Pa  ')
     id_ph         = register_diag_field &
      (mod_name, 'ph', axes(halfindx), Time, 'p at half levels',  'Pa')
     id_z        = register_diag_field &
      (mod_name, 'z', axes(1:3), Time,  'height  ',  'meters  ')
     id_zh        = register_diag_field &
      (mod_name, 'zh', axes(halfindx), Time, 'height at half levs', &
                                                              'meters')
   endif

   id_mr_ozone   = register_diag_field &
      (mod_name, 'ozone', axes(1:3), Time, 'Ozone mass mixing ratio', &
                                                   'kg (o3) / kg (air)')

   id_T          = register_diag_field &
      (mod_name, 'T', axes(1:3), Time, 'Temp at full levels',  'deg K ')

   id_sh         = register_diag_field &
      (mod_name, 'sh', axes(1:3), Time, &
        'vapor specific humidity at full levels',  'kg(h2o) / kg(air) ')

   id_rh         = register_diag_field &
      (mod_name, 'relhum', axes(1:3), Time, &
                      'relative humidity at full levels',  'fraction ')

   id_tot_h2o   = register_diag_field &
      (mod_name, 'tot_h2o', axes(1:3), Time, &
                                  'total water substance',  &
                            'kg(h2o) / kg(air) ' )

   id_lsca_cmip       = register_diag_field &
      (mod_name, 'lsca_cmip', axes(1:3), Time, &
                'ls liq cld fraction',  'fraction ', &
                mask_variant = .true., &
                   missing_value = missing_value)

   id_cca_cmip   = register_diag_field &
      (mod_name, 'cca_cmip', axes(1:3), Time, &
                 'convective liq cld fraction',  'fraction ', &
                mask_variant = .true., &
                   missing_value = missing_value)

   id_tca       = register_diag_field &
      (mod_name, 'tca', axes(1:3), Time, &
                                  'total cld fraction',  'fraction ')

   id_cca        = register_diag_field &
      (mod_name, 'cca', axes(1:3), Time, &
                           'convective cld fraction',  'fraction ')

   id_lsliq      = register_diag_field &
      (mod_name, 'lsliq', axes(1:3), Time, &
                                  'large scale cld liq',  'kg / kg  ')

   id_lsice      = register_diag_field &
      (mod_name, 'lsice', axes(1:3), Time, &
                                   'large scale cld ice',  'kg / kg  ')

   id_ccliq      = register_diag_field &
      (mod_name, 'ccliq', axes(1:3), Time, &
                                   'convective  cld liq',  'kg / kg  ')

   id_ccice      = register_diag_field &
      (mod_name, 'ccice', axes(1:3), Time, &
                                   'convective  cld ice',  'kg / kg  ')

   id_fl_lsrain  = register_diag_field &
      (mod_name, 'fl_lsrain', axes(1:3), Time, &
                             'large scale rain flx',  'kg / (m**2 s)')

   id_fl_lssnow  = register_diag_field &
      (mod_name, 'fl_lssnow', axes(1:3), Time, &
                             'large scale snow flx',  'kg / (m**2 s)')

   id_fl_lsgrpl  = register_diag_field &
      (mod_name, 'fl_lsgrpl', axes(1:3), Time, &
                           'large scale graupel flx',  'kg / (m**2 s)')

   id_fl_ccrain  = register_diag_field &
      (mod_name, 'fl_ccrain', axes(1:3), Time, &
                            'cnvctv scale rain flx',  'kg / (m**2 s)')

   id_fl_ccsnow  = register_diag_field &
      (mod_name, 'fl_ccsnow', axes(1:3), Time, &
                            'cnvctv scale snow flx',  'kg / (m**2 s)')

   id_reff_lsclliq_cmip  = register_diag_field &
      (mod_name, 'reff_lsclliq_cmip', axes(1:3), Time, &
           'ls liq cld drop size*cfrac ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_ccclliq_cmip  = register_diag_field &
      (mod_name, 'reff_ccclliq_cmip', axes(1:3), Time, &
         'cv liq cld drop size*cfrac ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_lsclliq  = register_diag_field &
      (mod_name, 'reff_lsclliq', axes(1:3), Time, &
               'ls liq cld drop size ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_lsclice  = register_diag_field &
      (mod_name, 'reff_lsclice', axes(1:3), Time, &
                'ls ice cld drop size ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_lsprliq  = register_diag_field &
      (mod_name, 'reff_lsprliq', axes(1:3), Time, &
                                       'ls liq prcp drop size ',  'm')

   id_reff_lsprice  = register_diag_field &
      (mod_name, 'reff_lsprice', axes(1:3), Time, &
                                        'ls ice prcp drop size ',  'm')

   id_reff_ccclliq  = register_diag_field &
      (mod_name, 'reff_ccclliq', axes(1:3), Time, &
             'cv liq cld drop size ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_ccclice  = register_diag_field &
      (mod_name, 'reff_ccclice', axes(1:3), Time, &
          'cv ice cld drop size ',  'm', mask_variant = .true., &
                   missing_value = missing_value)

   id_reff_ccprliq  = register_diag_field &
      (mod_name, 'reff_ccprliq', axes(1:3), Time, &
                                        'cv liq prcp drop size ',  'm')

   id_reff_ccprice  = register_diag_field &
      (mod_name, 'reff_ccprice', axes(1:3), Time, &
                                        'cv ice prcp drop size ',  'm')

   id_dtau_s  = register_diag_field &
      (mod_name, 'dtau_s', axes(1:3), Time, &
                   'ls cloud optical depth ',  'dimensionless')

   id_dtau_c  = register_diag_field &
      (mod_name, 'dtau_c', axes(1:3), Time, &
                    'cv cloud optical depth ',  'dimensionless')

   id_dem_s  = register_diag_field &
      (mod_name, 'dem_s', axes(1:3), Time, &
                             'ls cloud emissivity ',  'dimensionless')

   id_dem_c  = register_diag_field &
      (mod_name, 'dem_c', axes(1:3), Time, &
                             'cv cloud emissivity  ',  'dimensionless')

   id_skt        = register_diag_field &
      (mod_name, 'skt', axes(1:2), Time, 'skin temp',  'deg K')

   id_sunlit     = register_diag_field &
      (mod_name, 'sunlit', axes(1:2), Time, 'sun is shining?',  'none')

   id_land       = register_diag_field &
      (mod_name, 'land', axes(1:2), Time, 'land frac',  'fraction')

   id_sfcht      = register_diag_field &
      (mod_name, 'sfc_ht', axes(1:2), Time, 'height of surface',   &
                                                             'meters')

!---------------------------------------------------------------------
!    COSP output fields:
!---------------------------------------------------------------------
   allocate (id_dbze94(Ncolumns))
   if (use_vgrid) then
     allocate (id_cloudsatcfad(DBZE_BINS))
     allocate (id_calipsosrcfad(SR_BINS ))
     allocate (id_cloudsatcfad_sat(DBZE_BINS))
     allocate (id_calipsosrcfad_sat(SR_BINS ))
   else
     allocate (id_cloudsatcfad_mdl(DBZE_BINS))
     allocate (id_calipsosrcfad_mdl(SR_BINS ))
   endif
   allocate (id_cloud_type     (Ncolumns ))
   do n=1, size(id_cloud_type,1)
     if (n <= 9) then
       write (chvers, '(i1)') n
     else if (n <=99) then
       write (chvers, '(i2)') n
     else
       call error_mesg ('cosp_driver', &      
        'can not process over 99 levels', FATAL)
     endif
     id_cloud_type(n) = register_diag_field &
         (mod_name, 'cloud_type_' // trim(chvers), axes(1:3), Time, &
           'Cloud type present in column ' // trim(chvers), 'none')
   end do

   if (cfg%Llidar_sim) then
     id_cltcalipso = register_diag_field &
      (mod_name, 'cltcalipso', axes(1:2), Time, &
          'Lidar Total Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

     id_cllcalipso = register_diag_field &
      (mod_name, 'cllcalipso', axes(1:2), Time, &
          'Lidar Low-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

     id_clmcalipso = register_diag_field &
      (mod_name, 'clmcalipso', axes(1:2), Time, &
          'Lidar Mid-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

     id_clhcalipso = register_diag_field &
      (mod_name, 'clhcalipso', axes(1:2), Time, &
          'Lidar High-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

     if (generate_orbital_output) then

       id_cltcalipso_sat = register_diag_field &
      (mod_name, 'cltcalipso_sat', axes(1:2), Time, &
          'Lidar Total Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

       id_cllcalipso_sat = register_diag_field &
      (mod_name, 'cllcalipso_sat', axes(1:2), Time, &
          'Lidar Low-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

       id_clmcalipso_sat = register_diag_field &
      (mod_name, 'clmcalipso_sat', axes(1:2), Time, &
          'Lidar Mid-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

       id_clhcalipso_sat = register_diag_field &
      (mod_name, 'clhcalipso_sat', axes(1:2), Time, &
          'Lidar High-level Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)

       id_clcalipso_sat = register_diag_field &
      (mod_name, 'clcalipso_sat', cosp_axes(csatindx), Time, &
       'Lidar Cloud Fraction (532 nm)', 'percent', &
          mask_variant = .true., missing_value=missing_value)

       id_sampling_sat = register_static_field &
      (mod_name, 'sampling_sat', cosp_axes(samplingindx),       &
       'Times sampled by Cloudsat', 'number', &
           missing_value=missing_value)

       id_location_sat = register_static_field &
      (mod_name, 'location_sat', cosp_axes(samplingindx2),       &
       'Satellite location index', 'counter', &
           missing_value=missing_value)

       id_lon_sat = register_diag_field &
      (mod_name, 'lon_sat', axes(1:2),  Time,      &
       'Satellite longitude', 'degrees E', &
          mask_variant = .true.,  missing_value=missing_value)

       id_lat_sat = register_diag_field &
      (mod_name, 'lat_sat', axes(1:2), Time,      &
       'Satellite latitude', 'degrees N', &
      mask_variant = .true.,     missing_value=missing_value)

       id_parasolrefl_sat = register_diag_field &
      (mod_name, 'parasol_refl_sat', cosp_axes(parasolindx), Time, &
      'PARASOL-like mono-directional reflectance', 'fraction', &
          mask_variant = .true., missing_value=missing_value)

     endif

     id_clcalipso = register_diag_field &
      (mod_name, 'clcalipso', cosp_axes(csatindx), Time, &
       'Lidar Cloud Fraction (532 nm)', 'percent', &
          mask_variant = .true., missing_value=missing_value)

     id_clcalipso_mdl = register_diag_field &
      (mod_name, 'clcalipso_mdl', axes(1:3), Time, &
       'Lidar Cloud Fraction (532 nm)', 'percent', &
          mask_variant = .true., missing_value=missing_value)
     id_parasolrefl = register_diag_field &
      (mod_name, 'parasol_refl', cosp_axes(parasolindx), Time, &
      'PARASOL-like mono-directional reflectance', 'fraction', &
          mask_variant = .true., missing_value=missing_value)
     id_betamol532 = register_diag_field &
        (mod_name, 'betamol532', axes(1:3       ), Time, &
           'Lidar Molecular Backscatter (532 nm)', &
           '(m sr)**(-1)', &
          mask_variant = .true., missing_value=missing_value)
     allocate (id_atb532(Ncolumns))
     do n=1, size(id_atb532,1)
       if (n <= 9) then
         write (chvers, '(i1)') n
       else if (n <=99) then
         write (chvers, '(i2)') n
       else
         call error_mesg ('cosp_driver', &      
          'can not process over 99 columns', FATAL)
       endif
       id_atb532(n) = register_diag_field &
        (mod_name, 'atb532_' // trim(chvers), axes(1:3       ), Time, &
           'Lidar Attenuated Total Backscatter (532 nm) column# ' // &
          & trim(chvers), '(m sr)**(-1)', &
          mask_variant = .true., missing_value=missing_value)
     end do

     do n=1, SR_BINS                       
       if (n <= 9) then
         write (chvers, '(i1)') n
       else if (n <=99) then
         write (chvers, '(i2)') n
       else
         call error_mesg ('cosp_driver', &      
          'can not process over 99 levels', FATAL)
       endif
       if (n == 1) then
         write (chvers2, '(f8.2)') -100.0            
       else
         write (chvers2, '(f8.2)') srbval(n-1)
       endif
       write (chvers3, '(f8.2)') srbval(n)
       if (use_vgrid) then
         id_calipsosrcfad(n) = register_diag_field &
          (mod_name, 'calipsosrcfad_' // trim(chvers),  &
            cosp_axes(csatindx ), Time, &
              'Fractional area with Lidar 532 nm Scattering Ratio  &
              &between' // trim(chvers2) // ' and' // trim(chvers3) // &
                    ' -- bin' // trim(chvers),  'fraction', &
                    mask_variant = .true., missing_value=missing_value)
         if (generate_orbital_output) then
           id_calipsosrcfad_sat(n) = register_diag_field &
             (mod_name, 'calipsosrcfad_sat_' // trim(chvers),  &
            cosp_axes(csatindx ), Time, &
              'Fractional area with Lidar 532 nm Scattering Ratio  &
              &between' // trim(chvers2) // ' and' // trim(chvers3) // &
                    ' -- bin' // trim(chvers),  'fraction', &
                    mask_variant = .true., missing_value=missing_value)
         endif
       else
         id_calipsosrcfad_mdl(n) = register_diag_field &
           (mod_name, 'calipsosrcfad_mdl_' // trim(chvers), axes(1:3), &
          Time, 'Fractional area with Lidar 532 nm Scattering Ratio  &
           &between' // trim(chvers2) // ' and' // trim(chvers3) // &
                ' -- bin' // trim(chvers), 'fraction', &
                mask_variant = .true., missing_value=missing_value)
       endif
     end do
   endif  !(Llidar_sim)

 if (cfg%Lradar_sim) then
   do n=1, size(id_dbze94,1)
     if (n <= 9) then
       write (chvers, '(i1)') n
     else if (n <=99) then
       write (chvers, '(i2)') n
     else
       call error_mesg ('cosp_driver', &      
        'can not process over 99 levels', FATAL)
     endif
     id_dbze94(n) = register_diag_field &
       (mod_name, 'dbze94_' // trim(chvers), axes(1:3), Time, &
      'Radar Effective Reflectivity Factor in dBZe (94 GHz) column# ' &
            // trim(chvers), 'dBZe')
   end do

   do n=1, DBZE_BINS              
     if (n <= 9) then
       write (chvers, '(i1)') n
     else if (n <=99) then
       write (chvers, '(i2)') n
     else
       call error_mesg ('cosp_driver', &      
        'can not process over 99 levels', FATAL)
     endif
     write (chvers2, '(i6)') INT(cfad_ze_min + float(n-1)*cfad_ze_width)
     write (chvers3, '(i6)') INT(cfad_ze_min + float(n)*cfad_ze_width)
     if (use_vgrid) then
       id_cloudsatcfad(n) = register_diag_field &
          (mod_name, 'cloudsatcfad_' // trim(chvers),   &
           cosp_axes(csatindx), Time, &
           'Fractional area with radar reflectivity (94 GHz) between ' &
              // trim(chvers2) //  ' and' // trim(chvers3) //  &
               ' dbZe -- bin # '  //  trim(chvers),   'fraction', &
                mask_variant = .true., missing_value=missing_value)
       if (generate_orbital_output) then
         id_cloudsatcfad_sat(n) = register_diag_field &
          (mod_name, 'cloudsatcfad_sat_' // trim(chvers),   &
           cosp_axes(csatindx), Time, &
           'Fractional area with radar reflectivity (94 GHz) between ' &
              // trim(chvers2) //  ' and' // trim(chvers3) //  &
               ' dbZe -- bin # '  //  trim(chvers),   'fraction', &
                mask_variant = .true., missing_value=missing_value)
       endif
     else
       id_cloudsatcfad_mdl(n) = register_diag_field &
           (mod_name, 'cloudsatcfad_mdl_' // trim(chvers), axes(1:3), &
              Time, 'Fractional area with radar reflectivity &
             &(94 GHz) between ' // trim(chvers2) //  ' and' // &
             & trim(chvers3) //  ' dbZe -- bin # '  //  trim(chvers),  &
             'fraction', &
             mask_variant = .true., missing_value=missing_value)
     endif
   end do
 endif ! (Lradar_sim)


 if (cfg%Lradar_sim .and. cfg%Llidar_sim) then
   id_cltlidarradar = register_diag_field &
      (mod_name, 'cltlidarradar', axes(1:2), Time, &
          'Lidar and Radar Total Cloud Fraction',  'percent', &
          mask_variant = .true., missing_value=missing_value)
   id_clcalipso2 = register_diag_field &
      (mod_name, 'clcalipso2', cosp_axes(csatindx), Time, &
'Cloud frequency of occurrence as seen by CALIPSO but not CloudSat', &
         'percent', &
          mask_variant = .true., missing_value=missing_value)

   if (generate_orbital_output) then
     id_clcalipso2_sat = register_diag_field &
      (mod_name, 'clcalipso2_sat', cosp_axes(csatindx), Time, &
'Cloud frequency of occurrence as seen by CALIPSO but not CloudSat', &
         'percent', &
          mask_variant = .true., missing_value=missing_value)
   endif

   id_clcalipso2_mdl = register_diag_field &
      (mod_name, 'clcalipso2_mdl', axes(1:3), Time, &
'Cloud frequency of occurrence as seen by CALIPSO but not CloudSat', &
         'percent', &
          mask_variant = .true., missing_value=missing_value)
 endif !(cfg%Lradar_sim .and. cfg%Llidar_sim) 

 if (cfg%Lisccp_sim) then
   id_tclisccp = register_diag_field &
      (mod_name, 'tclisccp', axes(1:2), Time, &
          'Total Cloud Fraction as Calculated by the ISCCP Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_ctpisccp = register_diag_field &
      (mod_name, 'ctpisccp', axes(1:2), Time, &
       'Mean Cloud Top Pressure *CPCT as Calculated by the ISCCP Simulator', &
         'Pa', mask_variant = .true., missing_value=missing_value)

   id_tbisccp = register_diag_field &
      (mod_name, 'tbisccp', axes(1:2), Time, &
       'Mean All-sky 10.5 micron brightness temp -- ISCCP Simulator', &
         'deg K', mask_variant = .true., missing_value=missing_value)

   id_tbclrisccp = register_diag_field &
      (mod_name, 'tbclrisccp', axes(1:2), Time, &
       'Mean Clr-sky 10.5 micron brightness temp -- ISCCP Simulator', &
         'deg K', mask_variant = .true., missing_value=missing_value)

   id_tauisccp = register_diag_field &
      (mod_name, 'tauisccp', axes(1:2), Time, &
       'Mean Optical Depth *CPCT as Calculated by the ISCCP Simulator', &
         'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_albisccp = register_diag_field &
      (mod_name, 'albisccp', axes(1:2), Time, &
       'Mean Cloud Albedo *CPCT as Calculated by the ISCCP Simulator', &
         'fraction', &
          mask_variant = .true., missing_value=missing_value)
   id_boxtauisccp = register_diag_field &
      (mod_name, 'boxtauisccp', cosp_axes(columnindx), Time, &
         'Optical Depth  from the ISCCP Simulator', 'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_boxptopisccp = register_diag_field &
      (mod_name, 'boxptopisccp', cosp_axes(columnindx), Time, &
          'Cloud Top Pressure from the ISCCP Simulator', 'Pa')
   allocate (id_boxtauisccp_n(Ncolumns))
   allocate (id_boxptopisccp_n(Ncolumns))
   do n=1,Ncolumns
     if (n <= 9) then
       write (chvers, '(i1)') n
     else if (n <=99) then
       write (chvers, '(i2)') n
     else
       call error_mesg ('cosp_driver', &      
                   'can not process over 99 levels', FATAL)
     endif

     id_boxtauisccp_n(n) = register_diag_field &
        (mod_name, 'boxtauisccp_' // trim(chvers), axes(1:2), Time, &
          'Optical Depth in stochastic Column' // trim(chvers) //  &
            ' from the ISCCP Simulator', 'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

     id_boxptopisccp_n(n) = register_diag_field &
       (mod_name, 'boxptopisccp_' // trim(chvers), axes(1:2), Time, &
          'Cloud Top Pressure in stochastic column' // trim(chvers)  &
             //' from the ISCCP Simulator', 'Pa', &
          mask_variant = .true., missing_value=missing_value)
   end do
   do n=1,7
     write (chvers, '(i1)') n
     write (chvers2, '(i6)') INT(isccp_pc_bnds(1,n)*1.0e-02)
     write (chvers3, '(i6)') INT(isccp_pc_bnds(2,n)*1.0e-02)
     id_clisccp(n) = register_diag_field &
       (mod_name, 'clisccp_'// trim(chvers), cosp_axes(tauindx), &
          Time, 'ISCP Cld Frac for clouds between ' // trim(chvers2) &
             // ' and' // trim(chvers3) // ' hPa', 'percent', &
                  mask_variant = .true., missing_value=missing_value)
   end do

   do m=1,7
     write (chvers4, '(i1)') m
     write (chvers5, '(f4.1)') isccp_tau_bnds(1,m)
     write (chvers6, '(f8.1)') isccp_tau_bnds(2,m)
     do n=1,7
       write (chvers, '(i1)') n
       write (chvers2, '(i5)') INT(isccp_pc_bnds(1,n)*1.0e-02)
       write (chvers3, '(i5)') INT(isccp_pc_bnds(2,n)*1.0e-02)
       id_clisccp_n(m,n) = register_diag_field &
         (mod_name, 'clisccp_'// trim(chvers4)//'_' // trim(chvers), &
          axes(1:2), Time, 'ISCCP CldFrac - tau between ' // &
           trim(chvers5) // ' and ' // trim(chvers6) //  &
           ' , pr between ' // trim(chvers2) // ' and' // &
             trim(chvers3) // ' hPa',  'percent', &
                mask_variant = .true., missing_value=missing_value)
     end do
   end do
 endif !(Lisccp_sim)

  if (cfg%Lmisr_sim) then
   do n=1,MISR_N_CTH
       if (n <=9) then
       write (chvers, '(i1)') n
       else
       write (chvers, '(i2)') n
       endif
     write (chvers2, '(f6.1)') 1.0e-03*MISR_CTH_BNDS(1,n)
     write (chvers3, '(f6.1)') 1.0E-03*MISR_CTH_BNDS(2,n)
     id_misr(n) = register_diag_field &
       (mod_name, 'misr_'// trim(chvers), cosp_axes(tauindx), &
          Time, 'MISR Cld Frac for clouds with top between ' // trim(chvers2) &
             // ' and' // trim(chvers3) // ' km', 'percent', &
                  mask_variant = .true., missing_value=missing_value)
   end do

   do m=1,7
     write (chvers4, '(i1)') m
     write (chvers5, '(f4.1)') isccp_tau_bnds(1,m)
     write (chvers6, '(f8.1)') isccp_tau_bnds(2,m)
     do n=1,MISR_N_CTH
       if (n <=9) then
       write (chvers, '(i1)') n
       else
       write (chvers, '(i2)') n
       endif
       write (chvers2, '(f6.1)') 1.0e-03*MISR_CTH_BNDS(1,n)
       write (chvers3, '(f6.1)') 1.0e-03*MISR_CTH_BNDS(2,n)
       id_misr_n(m,n) = register_diag_field &
         (mod_name, 'misr_'// trim(chvers4)//'_' // trim(chvers), &
          axes(1:2), Time, 'MISR CldFrac - tau between ' // &
           trim(chvers5) // ' and ' // trim(chvers6) //  &
           ' , top between ' // trim(chvers2) // ' and' // &
             trim(chvers3) // ' km', 'percent', &
                mask_variant = .true., missing_value=missing_value)
     end do
   end do
 endif !(Lmisr_sim)

  if (cfg%Lmodis_sim) then

   id_tclmodis = register_diag_field &
      (mod_name, 'tclmodis', axes(1:2), Time, &
          'Total Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_locldmodis = register_diag_field &
      (mod_name, 'locldmodis', axes(1:2), Time, &
          'Low Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_mdcldmodis = register_diag_field &
      (mod_name, 'mdcldmodis', axes(1:2), Time, &
          'Middle Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_hicldmodis = register_diag_field &
      (mod_name, 'hicldmodis', axes(1:2), Time, &
          'High Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_lclmodis = register_diag_field &
      (mod_name, 'lclmodis', axes(1:2), Time, &
          'Total Liquid Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_iclmodis = register_diag_field &
      (mod_name, 'iclmodis', axes(1:2), Time, &
          'Total Ice Cloud Fraction as Calculated by the MODIS Simulator', &
          'percent', &
          mask_variant = .true., missing_value=missing_value)

   id_ttaumodis = register_diag_field &
      (mod_name, 'ttaumodis', axes(1:2), Time, &
          'Total Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_ltaumodis = register_diag_field &
      (mod_name, 'ltaumodis', axes(1:2), Time, &
          'Total Liquid Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_itaumodis = register_diag_field &
      (mod_name, 'itaumodis', axes(1:2), Time, &
          'Total Ice Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_tlogtaumodis = register_diag_field &
      (mod_name, 'tlogtaumodis', axes(1:2), Time, &
          'Total Log Mean Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_llogtaumodis = register_diag_field &
      (mod_name, 'llogtaumodis', axes(1:2), Time, &
          'Total Log Mean Liquid Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_ilogtaumodis = register_diag_field &
      (mod_name, 'ilogtaumodis', axes(1:2), Time, &
          'Total Log Mean Ice Optical Thickness*CPCT as Calculated by the MODIS Simulator', &
          'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

   id_lremodis = register_diag_field &
      (mod_name, 'lremodis', axes(1:2), Time, &
          ' Liquid Water particle Size*CPCT as Calculated by the MODIS Simulator', &
          'm', &
          mask_variant = .true., missing_value=missing_value)

   id_badlremodis = register_diag_field &
      (mod_name, 'badlsizemodis', axes(1:2), Time, &
          ' Flag for liquid size retrieval failure in the MODIS Simulator', &
          '1', &
          mask_variant = .true., missing_value=missing_value)

   id_badiremodis = register_diag_field &
      (mod_name, 'badisizemodis', axes(1:2), Time, &
          ' Flag for ice size retrieval failure in the MODIS Simulator', &
          '1', &
          mask_variant = .true., missing_value=missing_value)

   id_iremodis = register_diag_field &
      (mod_name, 'iremodis', axes(1:2), Time, &
          ' Ice Water particle Size*CPCT as Calculated by the MODIS Simulator', &
          'm', &
          mask_variant = .true., missing_value=missing_value)

   id_ctpmodis = register_diag_field &
      (mod_name, 'ctpmodis', axes(1:2), Time, &
          ' Mean Cloud Top Pressure*CPCT as Calculated by the MODIS Simulator', &
          'Pa', &
          mask_variant = .true., missing_value=missing_value)

   id_lwpmodis = register_diag_field &
      (mod_name, 'lwpmodis', axes(1:2), Time, &
          ' Mean Liquid Water Path*CPCT as Calculated by the MODIS Simulator', &
          'kg / ( m**2)',   &
          mask_variant = .true., missing_value=missing_value)

   id_iwpmodis = register_diag_field &
      (mod_name, 'iwpmodis', axes(1:2), Time, &
          ' Mean Ice Water Path*CPCT as Calculated by the MODIS Simulator', &
          'kg / ( m**2)',  &
          mask_variant = .true., missing_value=missing_value)

   allocate (id_taumodis_n(Ncolumns))
   allocate (id_ptopmodis_n(Ncolumns))
   allocate (id_sizemodis_n(Ncolumns))
   allocate (id_badsizemodis_n(Ncolumns))
   allocate (id_phasemodis_n(Ncolumns))
   do n=1,Ncolumns
     if (n <= 9) then
       write (chvers, '(i1)') n
     else if (n <=99) then
       write (chvers, '(i2)') n
     else
       call error_mesg ('cosp_driver', &      
                   'can not process over 99 levels', FATAL)
     endif

     id_taumodis_n(n) = register_diag_field &
        (mod_name, 'taumodis_' // trim(chvers), axes(1:2), Time, &
          'Optical Depth in stochastic Column' // trim(chvers) //  &
            ' from the MODIS Simulator', 'dimensionless', &
          mask_variant = .true., missing_value=missing_value)

     id_ptopmodis_n(n) = register_diag_field &
       (mod_name, 'ptopmodis_' // trim(chvers), axes(1:2), Time, &
          'Cloud Top Pressure in stochastic column' // trim(chvers)  &
             //' from the MODIS Simulator', 'hPa', &
          mask_variant = .true., missing_value=missing_value)

     id_sizemodis_n(n) = register_diag_field &
        (mod_name, 'sizemodis_' // trim(chvers), axes(1:2), Time, &
          'Particle Size in stochastic Column' // trim(chvers) //  &
            ' from the MODIS Simulator', 'meters', &
          mask_variant = .true., missing_value=missing_value)

     id_badsizemodis_n(n) = register_diag_field &
        (mod_name, 'badsizemodis_' // trim(chvers), axes(1:2), Time, &
          'Particle Size failures in stochastic Column' // trim(chvers) //  &
            ' from the MODIS Simulator', 'meters', &
          mask_variant = .true., missing_value=missing_value)

     id_phasemodis_n(n) = register_diag_field &
        (mod_name, 'phasemodis_' // trim(chvers), axes(1:2), Time, &
          'Phase in stochastic Column' // trim(chvers) //  &
            ' from the MODIS Simulator', 'unitless', &
          mask_variant = .true., missing_value=missing_value)

   end do
   do n=numPressureHistogramBins,1,-1
       if (n <=9) then
       write (chvers, '(i1)') n
       else
       write (chvers, '(i2)') n
       endif
     write (chvers2, '(f8.1)') nominalPressureHistogramBoundaries(1,n)
     write (chvers3, '(f8.1)') nominalPressureHistogramBoundaries(2,n)
     id_tauctpmodis(n) = register_diag_field &
       (mod_name, 'tauctpmodis_'// trim(chvers), cosp_axes(modistauindx), &
          Time, 'MODIS Cld Frac for clouds with top between ' // trim(chvers2) &
             // ' and' // trim(chvers3) // ' Pa', 'percent', &
                  mask_variant = .true., missing_value=missing_value)
   end do

   do m=1,numTauHistogramBins
     write (chvers4, '(i1)') m + 1
     write (chvers5, '(f6.1)') nominalTauHistogramBoundaries(1,m)
     write (chvers6, '(f6.1)') nominalTauHistogramBoundaries(2,m)
     do n=numPressureHistogramBins,1,-1
       if (n <=9) then
       write (chvers, '(i1)') n
       else
       write (chvers, '(i2)') n
       endif
       write (chvers2, '(f8.1)') nominalPressureHistogramBoundaries(1,n)
       write (chvers3, '(f8.1)') nominalPressureHistogramBoundaries(2,n)
       id_tauctpmodis_n(m,n) = register_diag_field &
         (mod_name, 'tauctpmodis_'// trim(chvers4)//'_' // trim(chvers), &
          axes(1:2), Time, 'MODIS CldFrac - tau between ' // &
           trim(chvers5) // ' and ' // trim(chvers6) //  &
           ' , top between ' // trim(chvers2) // ' and' // &
             trim(chvers3) // ' Pa', 'percent', &
                mask_variant = .true., missing_value=missing_value)
     end do
   end do
 endif !(Lmodis_sim)




  end subroutine diag_field_init 

!####################################################################

subroutine cosp_diagnostics_time_vary (Time_diag)

type(time_type), intent(in) :: Time_diag

    integer :: n

      if (generate_orbital_output) then
!----------------------------------------------------------------------
!    determine the time index of the current time in the satellite 
!    orbit data.
!----------------------------------------------------------------------
        do n= nsat_time_prev, num_sat_periods  
          if (Time_diag >= Time_start(n) .and.   &
                                           Time_diag < Time_end(n)) then
            nsat_time = n
            nsat_time_prev = nsat_time
            exit
          else
!   set nsat_time to 0 if current time not within sampling region
            nsat_time = 0
          endif
        end do
      endif



end subroutine cosp_diagnostics_time_vary



!####################################################################

subroutine cosp_diagnostics_endts      

    return

end subroutine cosp_diagnostics_endts

!####################################################################

subroutine output_cosp_fields   &
        (nlon,nlat,npoints, geomode, stlidar, stradar, isccp, modis,   &
         misr, sgradar, sglidar, sg, Time_diag, is, js, cloud_type,   &
         gbx, cfg, phalf_plus, zhalf_plus)

!---------------------------------------------------------------------
!     subroutine output_cosp_fields outputs fields relevant to the
!     cosp ismulator, both input and output.
!---------------------------------------------------------------------

integer,                            intent(in) :: nlon,nlat,npoints
integer,                            intent(in) :: geomode
integer,                            intent(in) :: is, js
type(cosp_config),                  intent(in) :: cfg   
real, dimension(npoints, ncolumns, nlevels),  intent(in) :: cloud_type
type(cosp_lidarstats), intent(in) :: stlidar
type(cosp_radarstats), intent(in) :: stradar
type(cosp_isccp     ), intent(in) :: isccp  
type(cosp_modis     ), intent(in) :: modis
type(cosp_misr      ), intent(in) :: misr   
type(cosp_sgradar   ), intent(in) :: sgradar
type(cosp_sglidar   ), intent(in) :: sglidar
type(cosp_subgrid   ), intent(in) :: sg
type(time_type)      , intent(in) :: Time_diag
type(cosp_gridbox)   , intent(in) :: gbx
real, dimension(nlon,nlat, nlevels+1), intent(in) :: phalf_plus, zhalf_plus

!   local variables:

      logical :: used
      integer :: n, m
      real, dimension(Nlon,Nlat) :: y2, y2save, alpha, y2sunlit 
      real, dimension(Nlon,Nlat) :: y2lsave, y2isave
      real, dimension(Nlon,Nlat,Nlevels) :: y3 
      real, dimension(Nlon,Nlat,Nlevels) :: y31,y32, y33,y34, y35, y36,y37 
      real, dimension(Nlon,Nlat,Nlevels) :: y3a
      real, dimension(Nlon,Nlat,Nlr    ) :: z3 
      real, dimension(Nlon,Nlat,Nlr    ) :: z3a
      real, dimension(Nlon,Nlat,Ncolumns) :: y4 
      real, dimension(Nlon,Nlat,PARASOL_NREFL) :: y5 
      real, dimension(Nlon,Nlat,Ncolumns,Nlevels) :: y6,y6a 
      real, dimension(Nlon,Nlat,Ncolumns,Nlr    ) :: z6,z6a 
      real, dimension(Nlon,Nlat,DBZE_BINS,Nlevels) :: y7,y7a 
      real, dimension(Nlon,Nlat,DBZE_BINS,Nlr    ) :: z7,z7a 
      real, dimension(Nlon,Nlat,SR_BINS,Nlevels) :: y8, y8a
      real, dimension(Nlon,Nlat,SR_BINS,Nlr    ) :: z8, z8a
      real, dimension(Nlon,Nlat,7,7            ) :: y9 
      real, dimension(Nlon,Nlat,numTauHistogramBins,  &
                                      numPressureHistogramBins  ) :: y13
      real, dimension(Nlon,Nlat,numTauHistogramBins+1,  &
                                      numPressureHistogramBins  ) :: y12
      real, dimension(Nlon,Nlat,7,MISR_N_CTH   ) :: y10
      logical, dimension (Nlon,Nlat,Nlevels) :: mask_y3a
      logical, dimension (Nlon,Nlat) :: lmsk
      integer  :: ie, je

     ie = is + nlon -1
     je = js + nlat -1
!----------------------------------------------------------------------
!    output the input fields to COSP. fields must be converted from
!    2d arrays (i,j)
!----------------------------------------------------------------------

!   2D fields:
   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%latitude, y2 = y2)
   used = send_data (id_lat       , y2, Time_diag, is, js )
  if (generate_orbital_output) then
   used = send_data (id_lat_sat   , y2, Time_diag, is, js,  mask =  &
                                           lflag_array(is:ie,js:je,nsat_time))
endif

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%longitude, y2 = y2)
   used = send_data (id_lon       , y2, Time_diag, is, js )
  if (generate_orbital_output) then
   used = send_data (id_lon_sat   , y2, Time_diag, is, js,  mask =  &
                                           lflag_array(is:ie,js:je,nsat_time))
endif

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%sunlit, y2 = y2sunlit)
   used = send_data (id_sunlit    , y2sunlit, Time_diag, is, js )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%skt, y2 = y2)
   used = send_data (id_skt       , y2, Time_diag, is, js )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%land, y2 = y2)
   used = send_data (id_land      , y2, Time_diag, is, js )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%u_wind, y2 = y2)
   used = send_data (id_u_wind    , y2, Time_diag, is, js )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%v_wind, y2 = y2)
   used = send_data (id_v_wind    , y2, Time_diag, is, js )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=gbx%sfc_height, y2 = y2)
   used = send_data (id_sfcht     , y2, Time_diag, is, js )

!   3D fields:

       used = send_data (id_ph    , phalf_plus, Time_diag, is, js, 1 )
       used = send_data (id_zh     , zhalf_plus, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%p,  y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_p         , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%zlev, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_z         , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%mr_ozone, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_mr_ozone  , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%T, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_T         , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%sh, y3 = y3)
   call flip_vert_index    (y3, nlevels,y37   )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%q, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_rh        , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%tca, y3 = y3)
   call flip_vert_index    (y3, nlevels,y35   )
   used = send_data (id_tca       , y35, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%cca, y3 = y3)
   call flip_vert_index    (y3, nlevels,y36   )
   used = send_data (id_cca       , y36, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%mr_hydro(:,:,I_lscliq), y3 = y3)
   call flip_vert_index    (y3, nlevels,y31   )
   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%mr_hydro(:,:,i_lscice), y3 = y3)
   call flip_vert_index    (y3, nlevels,y32   )
   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%mr_hydro(:,:,i_cvcliq), y3 = y3)
   call flip_vert_index    (y3, nlevels,y33   )
   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%mr_hydro(:,:,I_cvcice), y3 = y3)
   call flip_vert_index    (y3, nlevels,y34   )

   used = send_data (id_lsca_cmip  , y35-y36, Time_diag, is, js, 1, &
                                           mask =  y31 > 0)
   used = send_data (id_cca_cmip  , y36, Time_diag, is, js, 1, &
                                           mask =  y33 > 0)

   used = send_data (id_lsliq     , (y35-y36)*y31/((1.0+y36*(y33+y34))*&
                                       (1+y31)), Time_diag, is, js, 1 )

   used = send_data (id_lsice     , (y35-y36)*y32/((1.0+y36*  &
                            (y33+y34))*(1+y32)), Time_diag, is, js, 1 )

   used = send_data (id_ccliq     , y36*y33/((1.0+y36*(y33+y34))*  &
                                       (1+y33)), Time_diag, is, js, 1 )

   used = send_data (id_ccice     , y36*y34/((1.0+y36*(y33+y34))* &
                                       (1+y34)), Time_diag, is, js, 1 )

  used = send_data (id_sh        , y37/(1.+y36*(y33+y34)),  &
                                                 Time_diag, is, js, 1 )
   used = send_data (id_tot_h2o   , (y37 + (y35-y36)*y31/(1.+y31)+ &
          (y35-y36)*y32/(1.+y32)+y36*(y33/(1.+y33)+y34/(1.+y34)))/ &
                     ((1.0+y36*(y33+y34) )), Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%rain_ls, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_fl_lsrain , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%snow_ls, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_fl_lssnow , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%grpl_ls, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_fl_lsgrpl , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%rain_cv, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_fl_ccrain , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%snow_cv, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_fl_ccsnow , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_lscliq),&
                                                             y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a )
   used = send_data (id_reff_lsclliq , 0.5*y3a, Time_diag, is, js, 1, &
                   mask = y31 > 0.0 )
   used = send_data (id_reff_lsclliq_cmip , 0.5*y3a*(y35-y36) , Time_diag, is, js, 1, &
                   mask = y31 > 0.0 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_lscice),&
                                                               y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_lsclice , 0.5*y3a, Time_diag, is, js, 1 , &
                   mask = y32 > 0.0)

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_lsrain),&
                                                               y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_lsprliq , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_lssnow),&
                                                               y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_lsprice , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_cvcliq),&
                                                               y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_ccclliq , 0.5*y3a, Time_diag, is, js, 1 , &
                   mask = y33 > 0.0)
   used = send_data (id_reff_ccclliq_cmip , 0.5*y3a*y36 , Time_diag, is, js, 1 , &
                   mask = y33 > 0.0)

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_cvcice),&
                                                                y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_ccclice , 0.5*y3a, Time_diag, is, js, 1 , &
                   mask = y34 > 0.0)

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_cvrain),&
                                                              y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_ccprliq , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%reff(:,:,i_cvsnow),&
                                                              y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_reff_ccprice , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%dtau_s, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_dtau_s       , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%dtau_c, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_dtau_c       , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%dem_s, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_dem_s       , y3a, Time_diag, is, js, 1 )

   call map_point_to_ll (Nlon, Nlat, geomode, x2=gbx%dem_c, y3 = y3)
   call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_dem_c       , y3a, Time_diag, is, js, 1 )

!---------------------------------------------------------------------
!    process COSP output variables
!---------------------------------------------------------------------

 if (cfg%Llidar_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode, x1=stlidar%cldlayer(:,4),&
                                                               y2 = y2)
   used = send_data (id_cltcalipso,      y2, Time_diag, is, js , &
                                          mask = y2 /= missing_value )

   if (generate_orbital_output) then
     used = send_data (id_cltcalipso_sat,      y2, Time_diag, is, js , &
                                     mask = y2 /= missing_value  .and. &
                                          lflag_array(is:ie,js:je,nsat_time))
   endif

   call map_point_to_ll (Nlon, Nlat, geomode, x1=stlidar%cldlayer(:,1),&
                                                               y2 = y2)
   used = send_data (id_cllcalipso,      y2, Time_diag, is, js , &
                           mask = y2 /= missing_value )

   if (generate_orbital_output) then
     used = send_data (id_cllcalipso_sat,      y2, Time_diag, is, js , &
                                     mask = y2 /= missing_value  .and. &
                                          lflag_array(is:ie,js:je,nsat_time))
   endif

   call map_point_to_ll (Nlon, Nlat, geomode, x1=stlidar%cldlayer(:,2),&
                                                               y2 = y2)
   used = send_data (id_clmcalipso,      y2, Time_diag, is, js , &
                           mask = y2 /= missing_value )

   if (generate_orbital_output) then
     used = send_data (id_clmcalipso_sat,      y2, Time_diag, is, js , &
                                     mask = y2 /= missing_value  .and. &
                                          lflag_array(is:ie,js:je,nsat_time))
   endif

   call map_point_to_ll (Nlon, Nlat, geomode, x1=stlidar%cldlayer(:,3),&
                                                               y2 = y2)
   used = send_data (id_clhcalipso,      y2, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   if (generate_orbital_output) then
     used = send_data (id_clhcalipso_sat,      y2, Time_diag, is, js , &
                                     mask = y2 /= missing_value  .and. &
                                          lflag_array(is:ie,js:je,nsat_time))
   endif
 endif

 if(cfg%Lradar_sim .and.cfg%Llidar_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,  &
                                    x1=stradar%radar_lidar_tcc,y2 = y2)
   used = send_data (id_cltlidarradar, y2, Time_diag, is, js , &
                                           mask = y2 /= missing_value )
 endif

 if (cfg%Lisccp_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%totalcldarea,&
                                                       y2 = y2save)
   used = send_data (id_tclisccp,      y2save, Time_diag, is, js , &
                                           mask = y2save /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%meanptop,&
                                                              y2 = y2)
   where (y2save== 0.0 .and. y2sunlit == 1.0)
     alpha = 0.0
   elsewhere
     alpha =     y2*y2save
   endwhere

   used = send_data (id_ctpisccp , alpha     , Time_diag, is, js , &
                                           mask = y2sunlit == 1.0    )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%meantb,&
                                                              y2 = y2)

   used = send_data (id_tbisccp  , y2, Time_diag, is, js , &
                                           mask = y2sunlit == 1.0 )
!---

   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%meantbclr,&
                                                              y2 = y2)
   used = send_data (id_tbclrisccp  , y2, Time_diag, is, js , &
                                           mask = y2sunlit == 1.0 )

!----
   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%meantaucld,&
                                                             y2 = y2)
   where (y2save== 0.0 .and. y2sunlit == 1.0)
     alpha = 0.0
   elsewhere
     alpha = y2*y2save
   endwhere

   used = send_data (id_tauisccp  , alpha    , Time_diag, is, js , &
                                           mask = y2sunlit == 1.0 )

   call map_point_to_ll (Nlon, Nlat, geomode, x1=isccp%meanalbedocld,&
                                                              y2 = y2)
   where (y2save== 0.0 .and. y2sunlit == 1.0)
     alpha = 0.0
   elsewhere
     alpha = y2*y2save
   endwhere

   used = send_data (id_albisccp  , alpha, Time_diag, is, js , &
                                           mask = y2sunlit == 1.0 )
 endif

 if (cfg%Lmodis_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_Total_Mean,   &
                                                           y2 = y2save)
   used = send_data (id_tclmodis  , y2save, Time_diag, is, js , &
                                           mask = y2save /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_High_Mean,   &
                                                           y2 = y2)
   used = send_data (id_hicldmodis  , y2, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_Mid_Mean,   &
                                                           y2 = y2)
   used = send_data (id_mdcldmodis  , y2, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_Low_Mean,   &
                                                           y2 = y2)
   used = send_data (id_locldmodis  , y2, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_Water_Mean,   &
                                                              y2 = y2lsave)
   used = send_data (id_lclmodis  , y2lsave, Time_diag, is, js , &
                                           mask = y2lsave /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Fraction_Ice_Mean,   &
                                                              y2 = y2isave)
   used = send_data (id_iclmodis  , y2isave, Time_diag, is, js , &
                                           mask = y2isave /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                           x1=modis%Optical_Thickness_Total_Mean,   &
                                                              y2 = y2)

   where (y2save == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2save
   endwhere

   used = send_data (id_ttaumodis  , alpha, Time_diag, is, js , &
                                         mask = y2sunlit == 1.0 )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                           x1=modis%Optical_Thickness_Water_Mean,   &
                                                              y2 = y2)
   where (y2lsave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2lsave
   endwhere

   used = send_data (id_ltaumodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                           x1=modis%Optical_Thickness_Ice_Mean,   &
                                                              y2 = y2)
   where (y2isave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2save
     alpha = y2*y2isave
   endwhere

   used = send_data (id_itaumodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                        x1=modis%Optical_Thickness_Total_LogMean,   &
                                                              y2 = y2)
   where (y2save == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2save
   endwhere

   used = send_data (id_tlogtaumodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                           x1=modis%Optical_Thickness_Water_LogMean,   &
                                                              y2 = y2)
   where (y2lsave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2lsave
   endwhere

   used = send_data (id_llogtaumodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                           x1=modis%Optical_Thickness_Ice_LogMean,   &
                                                              y2 = y2)
   where (y2isave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2isave
   endwhere

   used = send_data (id_ilogtaumodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Particle_Size_Water_Mean,   &
                                                              y2 = y2)
   where (y2lsave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2lsave
   endwhere

   used = send_data (id_lremodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )
     lmsk(:,:) = (y2(:,:) < 0.0) .and. (y2(:,:) > -1.0)
   used = send_data (id_badlremodis  , y2, Time_diag, is, js , &
                                           mask = lmsk                )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Particle_Size_Ice_Mean,   &
                                                              y2 = y2)
  where (y2isave == 0.0 .and. y2sunlit == 1.0) 
    alpha = 0.
  elsewhere
     alpha = y2*y2isave
   endwhere

   used = send_data (id_iremodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )
     lmsk(:,:) = (y2(:,:) < 0.0) .and. (y2(:,:) > -1.0)
   used = send_data (id_badiremodis  , y2, Time_diag, is, js , &
                                           mask = lmsk                )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Cloud_Top_Pressure_Total_Mean,   &
                                                              y2 = y2)
   where (y2save == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2save
   endwhere

   used = send_data (id_ctpmodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Liquid_Water_Path_Mean,   &
                                                              y2 = y2)
   where (y2lsave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2lsave
   endwhere

   used = send_data (id_lwpmodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )

   call map_point_to_ll (Nlon, Nlat, geomode,  &
                               x1=modis%Ice_Water_Path_Mean,   &
                                                              y2 = y2)
   where (y2isave == 0.0 .and. y2sunlit == 1.0) 
     alpha = 0.
   elsewhere
     alpha = y2*y2isave
   endwhere

   used = send_data (id_iwpmodis  , alpha, Time_diag, is, js , &
                                           mask = y2 /= missing_value )


   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=modis%Column_Optical_Thickness, y3 = y4)
   do n=1,ncolumns
     used = send_data (id_taumodis_n(n), y4(:,:,n), Time_diag,  &
                       is, js, mask = y4(:,:,n) /= missing_value )
   end do

   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=modis%Column_Cloud_Top_Pressure, y3 = y4)
   do n=1,ncolumns
     used = send_data (id_ptopmodis_n(n), 0.01*y4(:,:,n), Time_diag, &
                       is, js, mask = y4(:,:,n) /= missing_value )
   end do
   
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=modis%Column_Particle_Size, y3 = y4)
   do n=1,ncolumns
     used = send_data (id_sizemodis_n(n), y4(:,:,n), Time_diag, &
                       is, js, mask = y4(:,:,n) > 0.0 )
    
     lmsk(:,:) = (y4(:,:,n) < 0.0) .and. (y4(:,:,n) > -1.0)
     used = send_data (id_badsizemodis_n(n), y4(:,:,n), Time_diag, &
                       is, js, mask =  lmsk       )
   end do

   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=modis%retrievedPhase      , y3 = y4)
   do n=1,ncolumns
     used = send_data (id_phasemodis_n(n), y4(:,:,n), Time_diag, &
                       is, js, mask = y4(:,:,n) /= missing_value )
   end do

 endif 

 if (use_vgrid) then
   if (cfg%Llidar_sim) then
     
     call map_point_to_ll (Nlon, Nlat, geomode, x2=stlidar%lidarcld,&
                                                             y3 = z3)
     used = send_data (id_clcalipso,      z3 , Time_diag, is, js, 1,  &
                                  mask = z3 (:,:,:) /= missing_value )
     if (generate_orbital_output) then
       used = send_data (id_clcalipso_sat,   z3 , Time_diag, is, js, 1,  &
                               mask = (z3 (:,:,:) /= missing_value) .and.& 
                                         lflag_array_temp(is:ie,js:je,:,nsat_time))
     endif
   endif
   if(cfg%Lradar_sim .and. cfg%Llidar_sim) then
     call map_point_to_ll (Nlon, Nlat, geomode,   &
                            x2=stradar%lidar_only_freq_cloud, y3 = z3)
     used = send_data (id_clcalipso2,      z3 , Time_diag, is, js, 1 , &
                                 mask = z3 (:,:,:) /= missing_value )
     if (generate_orbital_output) then
       used = send_data (id_clcalipso2_sat,  z3 , Time_diag, is, js, 1,  &
                               mask = (z3 (:,:,:) /= missing_value) .and.& 
                                         lflag_array_temp(is:ie,js:je,:,nsat_time))
     endif
   endif
 else
   if (cfg%Llidar_sim) then
     call map_point_to_ll (Nlon, Nlat, geomode, x2=stlidar%lidarcld,&
                                                              y3 = y3)
     call flip_vert_index    (y3, nlevels,y3a   )
     used = send_data (id_clcalipso_mdl, y3a, Time_diag, is, js, 1,  &
                                   mask = y3a(:,:,:) /= missing_value )
   endif
   if(cfg%Lradar_sim .and. cfg%Llidar_sim) then
     call map_point_to_ll (Nlon, Nlat, geomode,   &
                             x2=stradar%lidar_only_freq_cloud, y3 = y3)
     call flip_vert_index    (y3, nlevels,y3a   )
     used = send_data (id_clcalipso2_mdl, y3a, Time_diag, is, js, 1 , &
                                  mask = y3a(:,:,:) /= missing_value )
   endif
 endif

!3d arrays (i,j,columns):
 if (cfg%Lisccp_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=isccp%boxtau, y3 = y4)
   used = send_data (id_boxtauisccp, y4, Time_diag, is, js,  &
                           mask = y4 /= missing_value )
   do n=1,ncolumns
     used = send_data (id_boxtauisccp_n(n), y4(:,:,n), Time_diag,  &
                       is, js, mask = y4(:,:,n) /= missing_value )
   end do

   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x2=isccp%boxptop, y3 = y4)
   used = send_data (id_boxptopisccp, y4, Time_diag, is, js )
   do n=1,ncolumns
     used = send_data (id_boxptopisccp_n(n),      y4(:,:,n), Time_diag, &
                       is, js, mask = y4(:,:,n) /= missing_value )
   end do
 endif

!3d arrays (i,j,parasol_nrefl):
 if (cfg%Llidar_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                                       x2=stlidar%parasolrefl, y3 = y5)
   used = send_data (id_parasolrefl, y5, Time_diag, is, js, 1 , &
                                          mask = y5 /= missing_value )
   if (generate_orbital_output) then
     used = send_data (id_parasolrefl_sat, y5, Time_diag, is, js, 1 , &
                                     mask = y5 /= missing_value  .and. &
                                  lflag_array_parasol(is:ie,js:je,:,nsat_time))
   endif
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                                       x2=sglidar%beta_mol, y3 = y3)
     call flip_vert_index    (y3, nlevels,y3a   )
   used = send_data (id_betamol532, y3a, Time_diag, is, js, 1 , &
                                          mask = y3 /= missing_value )
 endif

!4d array (i,j,columns, levels):
   call map_point_to_ll (Nlon, Nlat, geomode, x3=sg%frac_out, y4 = y6)
     call flip_vert_index    (y6, nlevels,y6a   )
   do n=1, size(id_cloud_type,1)
     used = send_data (id_cloud_type(n), y6a(:,:,n,:),  &
                                                 Time_diag, is, js,1 )
   end do

!4d array (i,j,columns, levels):
 if(cfg%Lradar_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x3=sgradar%Ze_tot, y4 = y6)
   call flip_vert_index    (y6, nlevels,y6a   )
   do n=1, size(id_dbze94,1)
     used = send_data (id_dbze94(n), y6a(:,:,n,:), Time_diag, is, js,1 )
   end do

!4d array (i,j, dbze_bins, levels):
   if (use_vgrid) then
     call map_point_to_ll (Nlon, Nlat, geomode, x3=stradar%cfad_ze, &
                                                              y4 = z7)
     do n=1, size(id_cloudsatcfad,1)
       used = send_data (id_cloudsatcfad(n), z7(:,:,n,:), Time_diag, &
                        is, js, 1, mask = z7(:,:,n,:) /= missing_value )
       if (generate_orbital_output) then
         used = send_data (id_cloudsatcfad_sat(n), z7(:,:,n,:), Time_diag,&
                  is, js, 1, mask = (z7(:,:,n,:) /= missing_value) .and. & 
                       lflag_array_temp(is:ie,js:je,:,nsat_time))
       endif
     end do
   else
     call map_point_to_ll (Nlon, Nlat, geomode,   &
                                         x3=stradar%cfad_ze, y4 = y7)
     call flip_vert_index    (y7, nlevels,y7a   )
     do n=1, size(id_cloudsatcfad_mdl,1)
       used = send_data (id_cloudsatcfad_mdl(n), y7a(:,:,n,:),  &
                             Time_diag, is, js,1 , &
                                mask = y7a(:,:,n,:) /= missing_value )
     end do
   endif
endif

!4d array (i,j,columns, levels   ):
 if (cfg%Llidar_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                                          x3=sglidar%beta_tot, y4 = y6)
   call flip_vert_index    (y6, nlevels,y6a   )
   do n=1, size(id_atb532,1)
     used = send_data (id_atb532(n), y6a(:,:,n,:), Time_diag, is,  &
                        js, 1, mask = y6a(:,:,n,:) /= missing_value )
   end do

!4d array (i,j, sr_bins,levels):
   if (use_vgrid) then
     call map_point_to_ll (Nlon, Nlat, geomode,   &
                          x3=stlidar%cfad_sr, y4 = z8)
     do n=1, size(id_calipsosrcfad,1)
       used = send_data (id_calipsosrcfad(n), z8(:,:,n,:),    &
                          Time_diag, is, js,1 , &
                                 mask = z8 (:,:,n,:) /= missing_value )
       if (generate_orbital_output) then
         used = send_data (id_calipsosrcfad_sat(n), z8(:,:,n,:),    &
                          Time_diag, is, js,1 , &
                            mask = (z8 (:,:,n,:) /= missing_value) .and. & 
                                       lflag_array_temp(is:ie,js:je,:,nsat_time))
       endif
     end do
   else
     call map_point_to_ll (Nlon, Nlat, geomode,   &
                                          x3=stlidar%cfad_sr, y4 = y8)
     call flip_vert_index    (y8, nlevels,y8a   )
     do n=1, size(id_calipsosrcfad_mdl,1)
       used = send_data (id_calipsosrcfad_mdl(n), y8a (:,:,n,:),    &
                          Time_diag, is, js,1 , &
                                 mask = y8a(:,:,n,:) /= missing_value )
     end do
   endif
 endif

!4d array (i,j, isccp_tau,isccp_press):
 if (cfg%Lisccp_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                                            x3=isccp%fq_isccp, y4 = y9)
   do n=1, 7                           
     used = send_data (id_clisccp(n),      y9(:,:,:,n), Time_diag, is, &
                           js, 1, mask = y9(:,:,:,n) /= missing_value )
   end do

   do m=1,7
     do n=1, 7                           
       used = send_data (id_clisccp_n(m,n), y9(:,:,m,n), Time_diag, &
                           is, js, mask = y9(:,:,m,n) /= missing_value )
     end do
   end do
 endif


!4d array (i,j, modis_tau,modis_press):
 if (cfg%Lmodis_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
             x3=modis%Optical_Thickness_vs_Cloud_Top_Pressure, y4 = y12)
   y13(:,:,1:numTauHistogramBins,:) = y12(:,:,2:numTauHistogramBins+1,:)
   do n=1, numPressureHistogramBins   
     used = send_data (id_tauctpmodis(n), y13(:,:,:,n), Time_diag, is, &
                           js, 1, mask = y13(:,:,:,n) /= missing_value )
   end do

   do m=1,numTauHistogramBins
     do n=1, numPressureHistogramBins   
       used = send_data (id_tauctpmodis_n(m,n), y13(:,:,m,n), Time_diag, &
                           is, js, mask = y13(:,:,m,n) /= missing_value )
     end do
   end do
 endif

!4d array (i,j, isccp_tau,MISR_N_CTH ):
 if (cfg%Lmisr_sim) then
   call map_point_to_ll (Nlon, Nlat, geomode,   &
                                            x3=misr%fq_misr, y4 = y10)
   do n=1, MISR_N_CTH                  
     used = send_data (id_misr(n), y10(:,:,:,n), Time_diag, is, &
                           js, 1, mask = y10(:,:,:,n) /= missing_value )
   end do

   do m=1,7
     do n=1, MISR_N_CTH                  
       used = send_data (id_misr_n(m,n), y10(:,:,m,n), Time_diag, &
                          is, js, mask = y10(:,:,m,n) /= missing_value )
     end do
   end do
 endif

!-------------------------------------------------------------------
 
 
end subroutine output_cosp_fields



!#####################################################################

subroutine cosp_diagnostics_end 

!-------------------------------------------------------------------
!    deallocate the module arrays.
!-------------------------------------------------------------------
    deallocate (id_dbze94)
    deallocate (id_cloud_type)
    if (allocated(id_atb532)) deallocate (id_atb532)
    if (use_vgrid) then
      deallocate (id_cloudsatcfad)
      deallocate (id_calipsosrcfad)
      deallocate (id_cloudsatcfad_sat)
      deallocate (id_calipsosrcfad_sat)
    else
      deallocate (id_cloudsatcfad_mdl)
      deallocate (id_calipsosrcfad_mdl)
    endif
    
    if (generate_orbital_output) then
      deallocate (location, lflag_array, flag_array, lflag_array_temp, &
                  lflag_array_parasol, Time_start, Time_end)
    endif

end subroutine cosp_diagnostics_end

!#######################################################################

subroutine read_cloudsat_orbit 


!------------------------------------------------------------------------
!    subroutine read_cloudsat_orbit reads a netcdf file containing the
!    orbital position of the satellites as a function of time.
!------------------------------------------------------------------------

      real*4, dimension(:), allocatable    :: lat_in, lon_in
      integer*2, dimension(:), allocatable :: year_in
      byte, dimension(:), allocatable      ::  mon_in
      byte, dimension(:), allocatable      :: day_in, hour_in
      byte, dimension(:), allocatable      :: min_in
      real*4, dimension(:), allocatable    :: sec_in
      integer, dimension(:), allocatable   :: int_year_in
      integer, dimension(:), allocatable   ::  int_mon_in
      integer, dimension(:), allocatable   :: int_day_in, int_hour_in
      integer, dimension(:), allocatable   :: int_min_in
      real*8, dimension(:,:), allocatable  :: lat_out, lon_out

      character (len = *), parameter :: LAT_NAME  = "lat"
      character (len = *), parameter :: LON_NAME  = "lon"
      character (len = *), parameter :: YEAR_NAME = "year"
      character (len = *), parameter ::  MON_NAME = "month"
      character (len = *), parameter ::  DAY_NAME = "day"
      character (len = *), parameter :: HOUR_NAME = "hour"
      character (len = *), parameter ::  MIN_NAME = "minute"
      character (len = *), parameter ::  SEC_NAME = "second"

      integer          :: lat_varid, lon_varid, year_varid, day_varid,  &
                          mon_varid, hour_varid, min_varid, sec_varid
      integer          :: ncid
      integer          :: nlocs
      integer (kind=4) :: rcode, recdim
      type(time_type)  :: Time
      integer          :: k, mm, ptctr, n, ll, j, i
      integer          :: yeara, montha, daya, houra, minutea, seconda
      integer          :: yearb, monthb, dayb, hourb, minuteb, secondb
      integer          :: is, ie, js, je
      real             :: UNSET = -500.
      integer          :: calendar, nstart
      logical          :: used
      integer          :: ndims, nvars, ngatts
      integer          :: ndsize
      character*31     :: dummy
   
!------------------------------------------------------------------------
!    open the netcdf file. 
!------------------------------------------------------------------------
      ncid = ncopn (orbital_filename,   0, rcode)

!------------------------------------------------------------------------
!    determine number of dimensions (ndims); current file has 
!    only 1 ("location")
!------------------------------------------------------------------------
      call ncinq (ncid, ndims, nvars, ngatts, recdim, rcode)

!------------------------------------------------------------------------
!    determine value of the location dimension (nlocs) to use to dimension
!    arrays allocated below.
!------------------------------------------------------------------------
      do n=1,ndims
        call ncdinq(ncid, n, dummy, ndsize, rcode)
        if (trim(dummy) == 'location') then
          nlocs = ndsize
        endif
      end do

!------------------------------------------------------------------------
!    allocate arrays to hold the data read from the file.
!------------------------------------------------------------------------
      allocate (lat_in(nlocs), lon_in(nlocs), year_in(nlocs),  &
                mon_in(nlocs), day_in(nlocs), hour_in(nlocs),  &
                min_in(nlocs), sec_in(nlocs), int_year_in(nlocs), &
                int_mon_in(nlocs), int_day_in(nlocs), int_hour_in(nlocs), &
                int_min_in(nlocs) )
      allocate (lat_out(num_sat_periods, max_sdgs_per_sat_period), &
                lon_out(num_sat_periods, max_sdgs_per_sat_period) )
 
!------------------------------------------------------------------------
!    obtain the var_ids for the needed variables.
!------------------------------------------------------------------------

      lat_varid = ncvid(ncid, LAT_NAME , rcode)
      lon_varid = ncvid(ncid, LON_NAME , rcode)
      year_varid = ncvid(ncid, YEAR_NAME , rcode)
      mon_varid = ncvid(ncid, MON_NAME , rcode)
      day_varid = ncvid(ncid, DAY_NAME , rcode)
      hour_varid = ncvid(ncid, HOUR_NAME , rcode)
      min_varid = ncvid(ncid, MIN_NAME , rcode)
      sec_varid = ncvid(ncid, SEC_NAME , rcode)

!------------------------------------------------------------------------
!    read the netcdf data.
!------------------------------------------------------------------------
      call ncvgt (ncid, lat_varid, 1, nlocs, lat_in, rcode)
      call ncvgt (ncid, lon_varid, 1, nlocs, lon_in, rcode)
      call ncvgt (ncid, year_varid, 1, nlocs, year_in, rcode)
      call ncvgt (ncid, mon_varid, 1, nlocs, mon_in, rcode)
      call ncvgt (ncid, day_varid, 1, nlocs, day_in, rcode)
      call ncvgt (ncid, hour_varid, 1, nlocs, hour_in, rcode)
      call ncvgt (ncid, min_varid, 1, nlocs, min_in, rcode)
      call ncvgt (ncid, sec_varid, 1, nlocs, sec_in, rcode)

      call ncclos (ncid, rcode)

!------------------------------------------------------------------------
!    convert non-integer fields to integers.
!------------------------------------------------------------------------
      int_year_in = year_in
      int_mon_in = mon_in
      int_day_in = day_in
      int_hour_in = hour_in
      int_min_in = min_in

!------------------------------------------------------------------------
!    convert longitude to lie between 0 --> 360, rather than -180 --> 180.
!------------------------------------------------------------------------
      do  mm=1, size(lon_in)
        if (lon_in(mm) < 0.) then
          lon_in(mm) = lon_in(mm) + 360.
        endif
      end do

!------------------------------------------------------------------------
!    define the start and end of each time period for which the satellite 
!    orbital curtain data is desired. it is centered on sat_begin_time from
!    the cosp_input namelist.
!------------------------------------------------------------------------
      Time_start(1) = set_date (sat_begin_time(1), sat_begin_time(2),  &
                                sat_begin_time(3), sat_begin_time(4),  &
                                sat_begin_time(5), sat_begin_time(6))  - &
                                                   set_time(sat_period/2,0)
      Time_end(1) = Time_start(1) + set_time(sat_period, 0)

      do mm = 2,num_sat_periods 
        Time_start(mm) = Time_start(mm-1) + set_time(sat_period, 0)      
        Time_end  (mm) = Time_end  (mm-1) + set_time(sat_period, 0)      
      end do

!------------------------------------------------------------------------
!    initialize output variables.
!------------------------------------------------------------------------
      lat_out = UNSET
      lon_out = UNSET
      flag_array = 0.
      lflag_array = .false.
      location = 0.

!------------------------------------------------------------------------
!    define the latitudes/longitudes coordinates over which the satellite 
!    passes during each of the requested model sampling periods.
!------------------------------------------------------------------------
      calendar = get_calendar_type()

      nstart = 1
      do k=1,num_sat_periods      
        ptctr = 0
        do n=nstart, nlocs
          if (calendar == NOLEAP) then
!------------------------------------------------------------------------
!    ignore 2/29 when using the noleap calendar
!------------------------------------------------------------------------
            if (int_mon_in(n) == 2 .and. int_day_in(n) == 29) cycle
          endif

!-------------------------------------------------------------------------
!    determine if satellite observation time n is in any of the requested 
!    sampling periods. if it is before the first sampling period, cycle. 
!    if it is within sampling period k, increment the counter of obser-
!    vation times ptctr and enter the satellite location in the output 
!    arrays as the ptctr occurrence for sampling period k. if the sampling 
!    period has ended, exit the loop.
!-------------------------------------------------------------------------
          Time = set_date(int_year_in(n), int_mon_in(n), int_day_in(n), &
                          int_hour_in(n), int_min_in(n), INT(sec_in(n)))
          if (Time < Time_start(k)) then
            cycle
          else if (Time >= Time_start(k) .and. Time < Time_end(k)) then
            ptctr = ptctr + 1
            if (ptctr >= max_sdgs_per_sat_period) then
              call error_mesg ('cosp_driver:read_cloudsat_orbit', &
                    ' Need to increase &cosp_input variable &
                                       &max_sdgs_per_sat_period', FATAL)
            endif
            lat_out(k, ptctr) = lat_in(n)
            lon_out(k,ptctr) = lon_in(n)
          else if (Time >= Time_end(k))  then

!-------------------------------------------------------------------------
!    reset starting index into observations for next sampling period.
!-------------------------------------------------------------------------
            nstart = n - 1
            exit
          endif
        end do  ! n

!-------------------------------------------------------------------------
!    reset counter for next sampling period.
!-------------------------------------------------------------------------
        ptctr = 0
      end do   ! k

!-------------------------------------------------------------------------
!    call get_local_indexes2 to map the latitudes/longitudes seen by the 
!    satellite during sampling period k to the closest model grid  point 
!    (is,js). set a logical to indicate that grid point (is,js) is seen 
!    during time period k. 
!-------------------------------------------------------------------------
      do k=1,num_sat_periods   
        do ll = 1,max_sdgs_per_sat_period
          if (lat_out(k,ll) == UNSET .and. lon_out(k,ll) == UNSET) exit
          call get_local_indexes2(lat_out(k,ll),lon_out(k,ll), is,js)
          if (is /= 0 .and. js /= 0 .and. is <= imax .and. js <= jmax) then
            lflag_array(is,js,k) = .true.
            location(is,js,k) = ll
          endif
        end do

!-------------------------------------------------------------------------
!     collect sampling frequency diagnostic, if desired.
!-------------------------------------------------------------------------
        if (id_sampling_sat > 0) then
          call get_date(Time_end(k), yearb, monthb, dayb, hourb,    &
                                                         minuteb, secondb)
          do j=1,jmax
            do i=1,imax
              if (lflag_array(i,j,k)) then
                flag_array(i,j,monthb) = flag_array(i,j,monthb) + 1.
              endif
            end do
          end do
        endif
      end do

!-------------------------------------------------------------------------
!    define additional flag arrays for other diagnostics.
!-------------------------------------------------------------------------
      do k=1,PARASOL_NREFL
        lflag_array_parasol(:,:,k,:) = lflag_array(:,:,:)
      end do
      do k=1,nlr
        lflag_array_temp(:,:,k,:) = lflag_array(:,:,:)
      end do
   
!-------------------------------------------------------------------------
!    output the satellite sampling frequency at each point for each
!    month of the year for which data is requested.
!-------------------------------------------------------------------------
      used = send_data (id_sampling_sat, flag_array,   &
                                             is_in=1, js_in=1, ks_in=1) 

!-------------------------------------------------------------------------
!    output the satellite location index for each sampling period
!    for which data is requested.
!-------------------------------------------------------------------------
      used = send_data (id_location_sat, location,   &
                          is_in=1, js_in=1, ks_in=1, mask = location > 0.) 

!-----------------------------------------------------------------------
!    deallocate local variables.
!-----------------------------------------------------------------------
      deallocate (lat_in, lon_in, year_in, mon_in,&
                  day_in, hour_in, min_in, sec_in,&
                  int_year_in, int_mon_in, int_day_in, &
                  int_hour_in, int_min_in, lat_out, lon_out )

end subroutine read_cloudsat_orbit

!#####################################################################



end module cosp_diagnostics_mod


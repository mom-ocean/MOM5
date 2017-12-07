module ocean_bbc_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matthew Harrison
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Hyun-Chul Lee
!</CONTACT>
!
!<OVERVIEW>
! Set bottom boundary conditions 
!</OVERVIEW>
!
!<DESCRIPTION>
! Set bottom boundary conditions 
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, 2012: Elements of MOM
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bbc_nml">
!
!  <DATA NAME="cdbot" UNITS="dimensionless" TYPE="real">
! Dimensionless coefficient for quadratic bottom drag. 
!  </DATA> 
!
!  <DATA NAME="bmf_implicit" TYPE="logical">
!  For incorporating the bottom momentum drag implicitly in time.
!  Default is bmf_implicit=.false. 
!  </DATA> 
!
!  <DATA NAME="cdbot_law_of_wall" TYPE="logical">
!  For determining bottom drag coefficient using a constant roughness length.
!  Will take maximum between cdbot and the computed value using law of
!  wall log-profile.  This option of use when have very very 
!  refined vertical resolution (say on order of meters) near the bottom.
!  Terrain following coordinates should use this option since they generally 
!  have very refined vertical grid spacing on topography. 
!  Default is cdbot_law_of_wall=.false. 
!  </DATA> 
!  <DATA NAME="law_of_wall_rough_length" UNITS="metre" TYPE="real">
!  Bottom roughness length.  Default is law_of_wall_rough_length=0.01m, following
!  the default used in the Princeton Ocean Model (POM). This value 
!  corresponds to "Law of Wall" physics.  
!  </DATA> 
!
!  <DATA NAME="cdbot_roughness_length" TYPE="logical">
!  For determining bottom drag coefficient using a map of the roughness length.
!  This approach is more relevant for coarse models
!  than the constant roughness length used in the cdbot_law_of_wall option. 
!  Default is cdbot_roughness_length=.false. 
!  </DATA> 
!
!  <DATA NAME="cdbot_roughness_uamp" TYPE="logical">
!  For determining bottom drag coefficient using a map of the roughness length
!  and tidal velocity amplitude.
!  This approach is more relevant for coarse models
!  than the constant roughness length used in the cdbot_law_of_wall option.
!  cdbot_lo <= cdbot(i,j) <= cdbot_hi.
!  Default is cdbot_roughness_length=.false.
!  </DATA>
!
!  <DATA NAME="cdbot_HH" UNITS="m" TYPE="real">
!  H0 in a parameterization of cdbot_roughness_uamp.
!  Default is cdbot_HH=1100.0.
!  </DATA>
!
!  <DATA NAME="cdbot_UU" UNITS="m/s" TYPE="real">
!  U0 in a parameterization of cdbot_roughness_uamp.
!  Default is cdbot_UU=1.0.
!  </DATA>
!
!  <DATA NAME="cdbot_wave" TYPE="logical">
!  For determining bottom drag coefficient using a map of the roughness length
!  and the surface wind wave field. The modified drag coefficient is calculated
!  following Grant and Mattsen. Likewise this method can be improved using
!  more sophisticated wave models including swell.
!  Default is cdbot_wave=.false. 
!  </DATA> 
!
!  <DATA NAME="uresidual" UNITS="m/s" TYPE="real">
!  Residual bottom velocity due to unresolved fluctuations (e.g., waves and tides)
!  that contribute to bottom dissipation.  Should be set to zero when running 
!  with explicit representation of tidal forcing and when waves are well resolved.
!  Default is uresidual=.05.
!  </DATA> 
!
!  <DATA NAME="uvmag_max" UNITS="m/s" TYPE="real">
!  Maximum magnitude of the bottom velocity used to compute the bottom 
!  momentum drag.  Default is uvmag_max=10.0.
!  </DATA> 
!
!  <DATA NAME="bmf_max" UNITS="N/m2" TYPE="real">
!  Maximum magnitude of the bottom momentum drag.
!  Default is bmf_max=1.0.
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>
!
use constants_mod,        only: epsln, pi
use diag_manager_mod,     only: register_diag_field, register_static_field, send_data
use fms_mod,              only: open_namelist_file, check_nml_error, write_version_number
use fms_mod,              only: close_file, read_data
use mpp_mod,              only: input_nml_file, mpp_error, NOTE, FATAL, WARNING, stdout, stdlog
use mpp_domains_mod,      only: mpp_update_domains, BGRID_NE

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, TERRAIN_FOLLOWING
use ocean_parameters_mod, only: onefourth, rho0, rho0r, cp_ocean, von_karman
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID
use ocean_types_mod,      only: ocean_velocity_type, ocean_domain_type
use ocean_types_mod,      only: ocean_grid_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_density_type, ocean_options_type
use ocean_types_mod,      only: ocean_external_mode_type
use ocean_workspace_mod,  only: wrk1_2d, wrk2_2d
use wave_types_mod,       only: ocean_wave_type
use ocean_wave_mod,       only: wave_model_is_initialised
use ocean_util_mod,       only: diagnose_2d, diagnose_2d_u, diagnose_sum

implicit none

private

#include <ocean_memory.h>

integer :: num_prog_tracers

! for Bgrid or Cgrid
integer :: horz_grid

! diagnostics 
integer :: id_eta_tend_geoheat      =-1
integer :: id_eta_tend_geoheat_glob =-1
integer :: id_geo_heat              =-1
integer :: id_gamma_bmf             =-1
integer :: id_drag_coeff            =-1
integer :: id_cdbot                 =-1
integer :: id_bmf_u                 =-1
integer :: id_bmf_v                 =-1
integer :: id_cur_wav_dr            =-1
integer :: id_wave_s                =-1
integer :: id_wave_u                =-1
integer :: id_iter                  =-1
logical :: used

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

character(len=128) :: version=&
     '$Id: ocean_bbc.F90,v 20.0 2013/12/14 00:10:36 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

public  :: ocean_bbc_init
public  :: get_ocean_bbc
private :: current_wave_drag_diag 
private :: wave_u_diag

logical :: module_is_initialized = .false.
logical :: first_call            = .true.

real    :: vonkarman2          
real    :: rho0_cdbot
real    :: law_of_wall_rough_length_r

real    :: cellarea_r
real    :: cp_r
real    :: uresidual2
real    :: small=1.e-20

integer :: index_temp=-1
integer :: index_salt=-1
integer :: index_redist_heat=-1

integer :: bbc_geothermal 
integer :: stdoutunit,stdlogunit 

real,    allocatable, dimension(:,:) :: cdbot_array       ! dimensionless drag coefficient
real,    allocatable, dimension(:,:) :: cdbot_lowall      ! dimensionless drag coefficient from law of wall
real,    allocatable, dimension(:,:) :: drag_coeff        ! rho0 * dimensionless drag coefficient
real,    allocatable, dimension(:,:) :: roughness_length  ! for computing bottom drag coefficient
real,    allocatable, dimension(:,:) :: tidal_uamp        ! for computing bottom drag coefficient
real,    allocatable, dimension(:,:) :: geo_heat          ! geothermal heating
real,    allocatable, dimension(:,:) :: wave_u            ! wave action at the sea floor, wave boundary layer       
real,    allocatable, dimension(:,:) :: wave_s            ! wave action at the sea floor, skin layer       
real,    allocatable, dimension(:,:) :: data       
integer, allocatable, dimension(:,:) :: grd_kbot          ! bottom k-level 

real, allocatable, dimension(:, :) :: uresidual2_2d
real, allocatable, dimension(:, :) :: tide_speed_t

! nml variables 
logical :: bmf_implicit             = .false. 
logical :: debug_this_module        = .false. 
logical :: cdbot_law_of_wall        = .false. 
logical :: cdbot_roughness_length   = .false.  
logical :: cdbot_wave               = .false.
logical :: cdbot_roughness_uamp     = .false.
logical :: use_geothermal_heating   = .false. 
real    :: convert_geothermal       = .001     ! for converting units from mW to W
real    :: law_of_wall_rough_length = 0.01     ! metre
real    :: cdbot                    = 2.5e-3   ! dimensionless
real    :: uresidual                = .05      ! m/sec
real    :: cdbot_hi                 = 3.0e-3   ! hi-end of cdbot for cdbot_roughness_length 
real    :: cdbot_lo                 = 1.0e-3   ! lo-end of cdbot for cdbot_roughness_length 
real    :: cdbot_gamma              = 40.0     ! for setting exp decay for cdbot w/ cdbot_roughness_length 
real    :: uvmag_max                = 10.0     ! max velocity scale (m/s) for computing bmf
real    :: bmf_max                  = 1.0      ! max bmf (N/m^2)
real    :: cdbot_HH                 = 1100.0   ! H0 of cdbot_roughness_uamp (m)
real    :: cdbot_UU                 = 1.0      ! U0 of cdbot_roughness_uamp (m/s)

logical :: read_tide_speed = .false.
real    :: uresidual2_max = 0.05

 namelist /ocean_bbc_nml/ bmf_implicit, cdbot, uresidual, cdbot_law_of_wall, law_of_wall_rough_length, &
                         cdbot_roughness_length, use_geothermal_heating, convert_geothermal,           &
                         cdbot_hi, cdbot_lo, cdbot_gamma, uvmag_max, bmf_max, debug_this_module,       &
                         cdbot_roughness_uamp, cdbot_HH, cdbot_UU, cdbot_wave

 namelist /ocean_bbc_ofam_nml/ read_tide_speed, uresidual2_max

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bbc_init">
!
! <DESCRIPTION>
! Initialize the bottom boundary condition module
! </DESCRIPTION>
!
subroutine ocean_bbc_init(Grid, Domain, Time, T_prog, Velocity, &
                          Ocean_options, vert_coordinate_type, hor_grid)

type(ocean_grid_type),   target, intent(in)    :: Grid
type(ocean_domain_type), target, intent(in)    :: Domain
type(ocean_time_type),           intent(in)    :: Time
type(ocean_prog_tracer_type),    intent(inout) :: T_prog(:)
type(ocean_velocity_type),       intent(inout) :: Velocity
type(ocean_options_type),        intent(inout) :: Ocean_options
integer,                         intent(in)    :: vert_coordinate_type
integer,                         intent(in)    :: hor_grid

logical :: read_roughness=.false. 
integer :: i,j,n, ioun, io_status, ierr
integer :: grd_axes(3)
real    :: depth, cdbot_term1, cdbot_term2

stdoutunit=stdout()
stdlogunit=stdlog() 

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_bbc_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bbc_nml')
read (input_nml_file, nml=ocean_bbc_ofam_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bbc_ofam_nml')
#else
ioun = open_namelist_file()
read(ioun, ocean_bbc_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bbc_nml')
rewind(ioun)
read(ioun, ocean_bbc_ofam_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bbc_ofam_nml')
call close_file(ioun)
#endif
write (stdoutunit,'(/)')
write (stdoutunit, ocean_bbc_nml)
write (stdlogunit, ocean_bbc_nml)
write (stdoutunit,'(/)')
write (stdoutunit, ocean_bbc_ofam_nml)
write (stdlogunit, ocean_bbc_ofam_nml)

#ifndef MOM_STATIC_ARRAYS
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid
horz_grid = hor_grid

num_prog_tracers = size(T_prog(:))
index_temp=-1;index_salt=-1;index_redist_heat=-1
do n=1,num_prog_tracers
   if (T_prog(n)%name == 'temp')        index_temp = n
   if (T_prog(n)%name == 'salt')        index_salt = n
   if (T_prog(n)%name == 'redist_heat') index_redist_heat = n
enddo

do n=1, num_prog_tracers
#ifndef MOM_STATIC_ARRAYS
   allocate(T_prog(n)%btf(isd:ied,jsd:jed))
#endif
   T_prog(n)%btf = 0.0
enddo

allocate (cdbot_array(isd:ied,jsd:jed))
allocate (cdbot_lowall(isd:ied,jsd:jed))
allocate (drag_coeff(isd:ied,jsd:jed))
allocate (roughness_length(isd:ied,jsd:jed))
allocate (tidal_uamp(isd:ied,jsd:jed))
allocate (geo_heat(isd:ied,jsd:jed))
allocate (data(isd:ied,jsd:jed))
allocate (wave_u(isd:ied,jsd:jed))
allocate (wave_s(isd:ied,jsd:jed))
allocate (grd_kbot(isd:ied,jsd:jed))

rho0_cdbot            = rho0*cdbot
vonkarman2            = von_karman*von_karman
geo_heat(:,:)         = 0.0
data(:,:)             = 0.0
roughness_length(:,:) = 0.0
tidal_uamp(:,:)       = 0.0
wave_u(:,:)           = 0.0
wave_s(:,:)           = 0.0

if(horz_grid == MOM_BGRID) then 
   grd_kbot(:,:) = Grd%kmu(:,:)
   grd_axes(:)   = Grd%vel_axes_uv(:)
else
   grd_kbot(:,:) = Grd%kmt(:,:)
   grd_axes(:)   = Grd%tracer_axes(:)
endif 

if(cdbot_wave .or. cdbot_roughness_length .or. cdbot_roughness_uamp) then 
  read_roughness = .true.
endif

drag_coeff(:,:)   = rho0*cdbot*Grd%mask(:,:,1) 
cdbot_array(:,:)  = cdbot*Grd%mask(:,:,1) 
cdbot_lowall(:,:) = cdbot*Grd%mask(:,:,1)

cellarea_r                 = 1.0/(epsln + Grd%tcellsurf)
cp_r                       = 1.0/cp_ocean
uresidual2                 = uresidual**2

if (read_tide_speed) then
    allocate(uresidual2_2d(isd:ied, jsd:jed))
    allocate(tide_speed_t(isd:ied, jsd:jed))
    call read_data('INPUT/tideamp.nc', 'tideamp', tide_speed_t, Domain%domain2d)
    write (stdout(),*) '==>ocean_bbc_mod: Completed read of tide_speed.'
    uresidual2_2d(:,:) = min(tide_speed_t(:,:)**2, uresidual2_max)
    call mpp_update_domains(uresidual2_2d(:,:), Dom%domain2d)
    deallocate(tide_speed_t)
else
   write(stdout(),'(a)') &
      '==>Note: NOT reading tide_speed for ocean_vert_tidal_mod.'
   call mpp_error(NOTE, &
      '==>ocean_vert_tidal_mod: Setting tide_speed to default value.')
endif

law_of_wall_rough_length_r = 1.0/law_of_wall_rough_length

if(cdbot_law_of_wall) then 
   write(stdoutunit,'(a)') &
   '==>Note from ocean_bbc_mod: Computing bottom drag coefficient from roughness length.'
   Ocean_options%bottom_roughness = 'Used bottom drag coefficient from Law of Wall roughness length.'
elseif (cdbot_roughness_length) then 
   Ocean_options%bottom_roughness = 'Used bottom drag from map of bottom roughness using law of wall.'
elseif (cdbot_roughness_uamp) then
   Ocean_options%bottom_roughness = 'Used bottom drag from map of bottom roughness w/ tidal amplitude from law of wall.'
elseif (cdbot_wave) then 
   Ocean_options%bottom_roughness = 'Used bottom drag from map of bottom roughness using law of wall.'&
                                    //' Also considered influence of the wind wave field.'
else 
   Ocean_options%bottom_roughness = 'Used bottom drag from specified drag coefficient.'
endif 

if(bmf_implicit) then 
   write(stdoutunit,'(a)') &
   '==>Note from ocean_bbc_mod: Computing bottom drag implicitly in time.'
   Ocean_options%bmf_implicit = 'Used bottom drag computed implicitly in time.'
   Velocity%bmf_implicit = .true.
else
   Ocean_options%bmf_implicit = 'Used bottom drag computed explicitly in time.'
   Velocity%bmf_implicit = .false.
endif 

! compute drag coefficient using law of wall with a roughness length 
! read in from a map.  This approach is more relevant for large-scale 
! coarse models than the constant roughness length used in the 
! cdbot_law_of_wall option. 


! read in topographic roughness_length. 
! assume roughness_length on B-grid velocity point if using MOM_BGRID.
! assume roughness_length on tracer point if using MOM_CGRID.

data = 0.0
if(cdbot_wave) then 
    data = law_of_wall_rough_length
endif 

if(read_roughness) then 
    call read_data('INPUT/roughness_cdbot.nc','roughness_cdbot', data, Domain%domain2d)
    if(horz_grid == MOM_BGRID) then 
       write (stdoutunit,*) '==>ocean_bbc_mod: read in topographic roughness length, assumed to be on B-grid velocity point.'
    else 
       write (stdoutunit,*) '==>ocean_bbc_mod: read in topographic roughness length, assumed to be on tracer point.'
    endif 
    do j=jsc,jec
       do i=isc,iec
          roughness_length(i,j) = Grd%mask(i,j,1)*max(0.0,data(i,j))
       enddo
    enddo
    call mpp_update_domains(roughness_length(:,:), Dom%domain2d) 
endif 

! bottom drag from law of wall using roughness_length
if(cdbot_roughness_length) then

    cdbot_lowall=0.0
    do j=jsd,jed
       do i=isd,ied
          depth = Grd%ht(i,j)                 
          if(depth > 0.0) then 
              cdbot_lowall(i,j)= (von_karman/( alog(depth/(roughness_length(i,j)+small)) + small ) )**2
          endif
       enddo
    enddo

    ! bottom drag coefficient; the range of cdbot is
    ! cd_beta <= cdbot(i,j) <= cd_alpha
    do j=jsd,jed
       do i=isd,ied
          cdbot_array(i,j) = Grd%mask(i,j,1)*cdbot_hi*(1.0-exp(-cdbot_gamma*cdbot_lowall(i,j)))
          cdbot_array(i,j) = max(cdbot_lo,cdbot_array(i,j))
       enddo
    enddo

endif

! scheme coded by H.C. Lee (May 2012) 
if(cdbot_roughness_uamp) then

    ! read in tidal velocity amplitude; assumed on U-grid
    data=0.0
    call read_data('INPUT/tideamp.nc','TIDEAMP', data, Domain%domain2d)
    if(horz_grid == MOM_BGRID) then 
       write (stdoutunit,*) '==>ocean_bbc_mod: completed read of tidal velocity amplitude, assumed to be on B-grid velocity point.'
    else 
       write (stdoutunit,*) '==>ocean_bbc_mod: completed read of tidal velocity amplitude, assumed to be on tracer point.'
    endif 
    do j=jsc,jec
       do i=isc,iec
          tidal_uamp(i,j) = Grd%mask(i,j,1)*max(0.0,data(i,j))
       enddo
    enddo
    call mpp_update_domains(tidal_uamp(:,:), Dom%domain2d)


    ! bottom drag from law of wall using roughness_length and tidal velocity amplitude
    cdbot_lowall=0.0
    do j=jsd,jed
       do i=isd,ied
          depth = Grd%ht(i,j)
          cdbot_term1= cdbot_HH/(roughness_length(i,j)+small)
          cdbot_term2= cdbot_UU/(tidal_uamp(i,j)+small)
          if(depth > 0.0) then
              cdbot_array(i,j)= (von_karman/( alog(cdbot_term1*cdbot_term2) + small ) )**2
          endif
       enddo
    enddo

    ! bottom drag coefficient; the range of cdbot is
    ! cdbot_lo <= cdbot(i,j) <= cdbot_hi
    do j=jsd,jed
       do i=isd,ied
          cdbot_array(i,j) = min(cdbot_hi,cdbot_array(i,j))
          cdbot_array(i,j) = max(cdbot_lo,cdbot_array(i,j))
       enddo
    enddo

endif

    
! save cdbot into Velocity type 
do j=jsd,jed
   do i=isd,ied
      Velocity%cdbot_array(i,j) = cdbot_array(i,j)
   enddo
enddo


if(use_geothermal_heating) then 
    write(stdoutunit,'(a)') &
         '==>Note from ocean_bbc_mod: Geothermal heating introduced at ocean bottom.'
    Ocean_options%geothermal_heating = 'Used geothermal heating introduced at ocean bottom.'

    ! fill geothermal heating array 
    geo_heat = 0.0
    call read_data('INPUT/geothermal_heating.nc', 'geo_heat', data(:,:), Dom%domain2d, timelevel=1)

    ! make sure there are no spurious negative values.
    ! also convert to units appropriate for temperature btf. 
    do j=jsc,jec
       do i=isc,iec
          geo_heat(i,j) = Grd%tmask(i,j,1)*max(0.0,data(i,j))*cp_r*convert_geothermal
       enddo
    enddo

else
    Ocean_options%geothermal_heating = 'Did NOT introduce geothermal heating.'
endif


if(vert_coordinate_type == TERRAIN_FOLLOWING .and. .not. cdbot_law_of_wall) then 
   call mpp_error(WARNING,&
   '==>Warning from ocean_bbc_init: recommend law_of_wall=.true. for TERRAIN_FOLLOWING coordinates.')
endif 


id_drag_coeff = register_diag_field ('ocean_model', 'drag_coeff',         &
                grd_axes(1:2), Time%model_time,                           &
                'Dimensionless bottom drag coefficient', 'dimensionless', &
                missing_value=missing_value, range=(/-1.0,1.e3/))
id_gamma_bmf  = register_diag_field ('ocean_model', 'gamma_bmf',       &
                grd_axes(1:2), Time%model_time,                        &
                'Bottom drag factor rho0*cdbot*uvmag', 'kg/(m^2 sec)', &
                missing_value=missing_value, range=(/-1.0,1.e12/))
id_bmf_u      = register_diag_field ('ocean_model', 'bmf_u',  &
                grd_axes(1:2), Time%model_time,               &
                'Bottom u-stress via bottom drag', 'N/m^2',   &
                missing_value=missing_value, range=(/-1.0,1.e2/))
id_bmf_v      = register_diag_field ('ocean_model', 'bmf_v',  &
                grd_axes(1:2), Time%model_time,               &
                'Bottom v-stress via bottom drag', 'N/m^2',   &
                missing_value=missing_value, range=(/-1.0,1.e2/))

id_cdbot      = register_static_field ('ocean_model', 'cdbot',            &
                grd_axes(1:2), 'Bottom drag coefficient', 'dimensionless',&     
                missing_value=missing_value, range=(/-1.0,1e6/))
if (id_cdbot > 0) then 
    used = send_data (id_cdbot, cdbot_array(:,:), Time%model_time, &
           rmask=Grd%mask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif

id_geo_heat = register_static_field ('ocean_model', 'geo_heat',   &
              Grd%tracer_axes(1:2), 'Geothermal heating', 'W/m^2',&     
              missing_value=missing_value, range=(/-10.0,1e6/),   &
              standard_name='upward_geothermal_heat_flux_at_sea_floor')
call diagnose_2d(Time, Grd, id_geo_heat, geo_heat(:,:)*cp_ocean)

id_eta_tend_geoheat = register_diag_field('ocean_model',             &
       'eta_tend_geoheat', Grd%tracer_axes(1:2),Time%model_time   ,  &
       'non-Bouss steric sea level tendency from geothermal heating',&
       'm/s',missing_value=missing_value,range=(/-1.e10,1.e10/))   

id_eta_tend_geoheat_glob = register_diag_field('ocean_model',        &
       'eta_tend_geoheat_glob', Time%model_time,                     &
       'non-Bouss steric sea level tendency from geothermal heating',&
       'm/s',missing_value=missing_value,range=(/-1.e10,1.e10/))   

id_cur_wav_dr = register_diag_field ('ocean_model', 'current_wave_stress',&
                Grd%tracer_axes(1:2), Time%model_time,                    &
                'combined current wave bottom stress', 'N/m^2',           &
                missing_value=missing_value, range=(/-1.0,1.e3/))

id_wave_s = register_diag_field ('ocean_model', 'wave_s', &
            Grd%tracer_axes(1:2), Time%model_time,        &
            'wave skin friction velocity', 'm/s',         &
            missing_value=missing_value, range=(/-1.0,1.e3/))

id_wave_u = register_diag_field ('ocean_model', 'wave_u',&
            Grd%tracer_axes(1:2), Time%model_time,       &
            'wave friction velocity', 'm/s',             &
            missing_value=missing_value, range=(/-1.0,1.e3/))
       
id_iter   = register_diag_field ('ocean_model', 'iter',&
            grd_axes(1:2), Time%model_time,            &
            'number of ustar iterations', 'm/s',       &
            missing_value=missing_value, range=(/0.0,1.e3/))


end subroutine ocean_bbc_init
! </SUBROUTINE> NAME="ocean_bbc_init"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_bbc">
!
! <DESCRIPTION>
! Set bottom boundary conditions for velocity and tracer.
!
! Dimensions of bottom momentum flux are 
! N/m^2 = (kg/m^3)*(m^2/s^2).
!
! Note the use of rho0 for the conversion from m^2/s^2 to 
! (kg/m^3)*(m^2/s^2).  We do not know the precise value 
! of cdbot, so the rho0 approximate value is well within 
! our level of uncertainty.  No reason therefore to 
! use in situ rho for this conversion, even when using 
! non-Boussinesq version of MOM.   
!
!
! Note that bmf needs to be computed on the data domain since the 
! halo values are accessed in ocean_vert_gotm.F90.
!
! </DESCRIPTION>
!
subroutine get_ocean_bbc(Time, Thickness, Dens, Velocity, T_prog, Waves)

type(ocean_time_type),          intent(in)    :: Time
type(ocean_thickness_type),     intent(in)    :: Thickness
type(ocean_density_type),       intent(in)    :: Dens 
type(ocean_velocity_type),      intent(inout) :: Velocity
type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
type(ocean_wave_type),          intent(inout) :: Waves

integer :: i, j, kbot, n
integer :: taum1, tau, tstep
real    :: uvmag, argument, distance 
real    :: umag, vmag
real    :: rhobot_inv, alphabot 

if (.not. module_is_initialized) then 
   call mpp_error(FATAL,'==>Error from ocean_bbc_mod (get_ocean_bbc): module must be initialized')
endif 

taum1 = Time%taum1
tau   = Time%tau

if(bmf_implicit) then 
  tstep=tau
else
  tstep=taum1
endif 

! Bottom tracer flux defaulted to zero.
do n=1,num_prog_tracers
   do j=jsd,jed
      do i=isd,ied
         T_prog(n)%btf(i,j) = 0.0
      enddo
   enddo
enddo

! T_prog(index_temp)%btf<0 means heat enters ocean; hence the minus sign. 
if (use_geothermal_heating) then
    do j=jsc,jec
       do i=isc,iec
          T_prog(index_temp)%btf(i,j) = -geo_heat(i,j)
       enddo
    enddo
    if(index_redist_heat > 0) then 
      do j=jsc,jec
         do i=isc,iec
            T_prog(index_redist_heat)%btf(i,j) = -geo_heat(i,j)
         enddo
      enddo
    endif   
endif


! set bottom drag coefficient
drag_coeff(:,:) = 0.0
if(horz_grid == MOM_BGRID) then 

   if(cdbot_law_of_wall) then 
       do j=jsd,jed
          do i=isd,ied
             kbot = grd_kbot(i,j)   
             if(kbot > 0) then 
                 distance        = Thickness%depth_zwu(i,j,kbot) - Thickness%depth_zu(i,j,kbot)
                 argument        = law_of_wall_rough_length_r*(distance + law_of_wall_rough_length)
                 drag_coeff(i,j) = rho0*max(cdbot,vonkarman2/(log(argument))**2)
             endif
          enddo
       enddo
   else 
       do j=jsd,jed
          do i=isd,ied
             drag_coeff(i,j) = rho0*cdbot_array(i,j)*Grd%umask(i,j,1)     
          enddo
       enddo
   endif

else 

   if(cdbot_law_of_wall) then 
       do j=jsd,jed
          do i=isd,ied
             kbot = grd_kbot(i,j)   
             if(kbot > 0) then 
                 distance        = Thickness%depth_zwt(i,j,kbot) - Thickness%depth_zt(i,j,kbot)
                 argument        = law_of_wall_rough_length_r*(distance + law_of_wall_rough_length)
                 drag_coeff(i,j) = rho0*max(cdbot,vonkarman2/(log(argument))**2)
             endif
          enddo
       enddo
   else 
       do j=jsd,jed
          do i=isd,ied
             drag_coeff(i,j) = rho0*cdbot_array(i,j)*Grd%tmask(i,j,1)     
          enddo
       enddo
   endif

endif ! horz_grid iftest


! modified stress if use idealized wave model 
if(cdbot_wave) then
    call current_wave_drag_diag(Thickness, Velocity, drag_coeff, Waves, Time, tstep)
endif


! set quadratic bottom momentum drag
! units of bmf are N/m^2
! introduce ceilings for uvmag and bmf
! to minimize potential for instability. 
Velocity%bmf(:,:,:) = 0.0

! if cgrid, ignore offset in velocity components for uvmag calculation 
do j=jsd,jed
   do i=isd,ied
       kbot = grd_kbot(i,j)
       Velocity%gamma(i,j) = 0.0
       if (kbot > 0) then
           if (read_tide_speed) then
               uvmag = sqrt(uresidual2_2d(i, j) &
                            + Velocity%u(i, j, kbot, 1, taum1)**2 &
                            + Velocity%u(i, j, kbot, 2, taum1)**2)
           else
               uvmag = sqrt(uresidual2 + Velocity%u(i, j, kbot, 1, tstep)**2 &
                            + Velocity%u(i, j, kbot, 2, tstep)**2)
           end if
           uvmag = min(uvmag,uvmag_max)
           Velocity%gamma(i,j)  = drag_coeff(i,j)*uvmag
           umag                 = abs(Velocity%u(i,j,kbot,1,tstep))
           vmag                 = abs(Velocity%u(i,j,kbot,2,tstep))
           Velocity%bmf(i,j,1)  = min(bmf_max, Velocity%gamma(i,j)*umag)*sign(1.0,Velocity%u(i,j,kbot,1,tstep))
           Velocity%bmf(i,j,2)  = min(bmf_max, Velocity%gamma(i,j)*vmag)*sign(1.0,Velocity%u(i,j,kbot,2,tstep))
       endif
   enddo
enddo

if (read_tide_speed) then
    deallocate(uresidual2_2d)
end if

! send diagnostics 

if (id_drag_coeff > 0) then 
  used = send_data(id_drag_coeff, rho0r*drag_coeff(:,:), Time%model_time, &
         rmask=Grd%mask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_gamma_bmf > 0) then 
  used = send_data(id_gamma_bmf, Velocity%gamma(:,:), Time%model_time, &
         rmask=Grd%mask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_bmf_u > 0) then 
  used = send_data(id_bmf_u, Velocity%bmf(:,:,1), Time%model_time, &
         rmask=Grd%mask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_bmf_v > 0) then 
  used = send_data(id_bmf_v, Velocity%bmf(:,:,2), Time%model_time, &
         rmask=Grd%mask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 


! some sea level diagnostics from geothermal heating 
if (id_eta_tend_geoheat > 0 .or. id_eta_tend_geoheat_glob > 0) then 

    wrk1_2d = 0.0
    do j=jsc,jec
       do i=isc,iec
          kbot=Grd%kmt(i,j)
          if(kbot > 0) then 
              rhobot_inv   = 1.0/(epsln+Dens%rho(i,j,kbot,tau))
              alphabot     = -rhobot_inv*Dens%drhodT(i,j,kbot)
              wrk1_2d(i,j) = -alphabot*rhobot_inv*T_prog(index_temp)%btf(i,j)
          endif
       enddo
    enddo

    call diagnose_2d(Time, Grd, id_eta_tend_geoheat, wrk1_2d(:,:))
    call diagnose_sum(Time, Grd, Dom, id_eta_tend_geoheat_glob, wrk1_2d, cellarea_r)

endif

end subroutine get_ocean_bbc
! </SUBROUTINE> NAME="get_ocean_bbc"

!#######################################################################
! <SUBROUTINE NAME="current_wave_drag_diag">
!
! <DESCRIPTION>
! calculates wave-current bottom shear stress 
! using model of Grand and Madsen(1979) J. Geophys. Res. 84, C4, 1797
! see Signell et. al (1990)  J. Geophys. Res. 95, C6, 9671
!
! input bot_vel: current velocity at u points
! Note: assumed that this is the velocity just above the bottom boundary layer.
!
! A relation between grain size, ripples steepness and and roughness is 
! assumed. More general relations are possible but not needed, since 
! the output is used only to parameterise resuspension of organic matter 
! in the ecosystem model ocean_shared/generic_tracers/generic_ERGOM.F90.
!
! output: effective drag coefficient drag_coeff. 
! It is valid for momentum flux from currents into the bottom, but
! from combined waves+current action.
! Velocity%current_wave_stress is the stress from waves and currents 
! to the sediment.
! 
! Note:
!
! 1/ drag_coeff in this routine arises from both (waves+current); 
! That is, ustar**2/uref**2 1 meter above the bottom.
!
! 2/ bottom velocity is taken just above the bottom boundary layer, and 
! assumed here to be at lowest u point.
!
! 3/ to understand the calculation of shear stress acting on grains 
! at sediment surface ("wrk1_2d(i,j) = (ustar2/ucomb)*0.3152"), 
! see Kuhrts et al. (2004) Eqs. (4,5,6,7).
! Assume a thin skin friction layer according to Smith, McLean (1977)
! thickness of the skin friction layer zskin scales with roughness length ruff
! grained sediments are characterised by median diameter d50[m]
! ripples at sea bottom have spacing lambda and height eta with steepness eta/lambda=0.1
! common approx. for grain roughness length = d50/12
! Nielsen (1983) form drag roughness length = (8/30)*eta*(eta/lambda) ==> (8/3000)*lambda
! Yalin (1977) lambda = 1000*d50 ==> form drag = (8/3)*d50
! ruff = grain + form drag = (1/12+8/3)*d50 = (33/12)*d50 = 33*grain
! Smith, McLean (1977) zskin = 0.09*grain*(lambda/grain)**0.8 = 165*grain = 5*ruff
! matching at zskin leads to log(zskin/ruff)/log(zskin/grain)) = log(5)/log(165) = 0.3152
! the current induced skin friction velocity is derived from matching skin friction 
! to wave boundary layer.
!
! 4/ Algorithm has yet to be updated for Cgrid. 
!
! April 2012
! martin.schmidt@io-warnemuende.de 
!
! </DESCRIPTION>
!
subroutine current_wave_drag_diag(Thickness, Velocity, drag_coeff, Waves, Time, tstep)
type(ocean_thickness_type),     intent(in)    :: Thickness
type(ocean_velocity_type),      intent(inout) :: Velocity
real,  dimension(isd:,jsd:),    intent(out)   :: drag_coeff    
type(ocean_wave_type),          intent(inout) :: Waves
type(ocean_time_type),          intent(in)    :: Time
integer,                        intent(in)    :: tstep

real,parameter:: twopi=2.*pi

  real    :: omega, cdrag, cold, cnew, ruff, rauh, dist, wave_p_u, wave_u_u
  real    :: bot_vel, u_vel, v_vel   ! Velocity in the bottommost layer. 
  real    :: ustar, ustar2, ucomb, ucskin  ! friction velocities in wave & skin boundary layers
  integer:: i, j, icount, kmu
  integer, parameter::  itmax=10

  if (first_call) then
     if ( .not. wave_model_is_initialised() ) then
       call mpp_error(FATAL,&
       '==>Error from ocean_bbc_mod (get_ocean_bbc): enable a wave model for option cdbot_wave ')
     endif
     first_call = .false.
  endif
  
  ! get the wave induced bottom friction velocities 
  ! wave_u applies inside the wave-boundary layer, wave_s inside the skin friction layer
  if (id_iter > 0) wrk2_2d = 0.0
  wrk1_2d    = 0.0
  drag_coeff = 0.0
  call wave_u_diag(wave_u, wave_s, Waves) 

  do j=jsc,jec 
     do i=isc,iec 

        kmu=Grd%kmu(i,j) 
        if (kmu >= 1) then

            ruff  = roughness_length(i,j) ! bottom roughness length (grain+form drag) [m]
            dist  = Thickness%depth_zwu(i,j,kmu) - Thickness%depth_zu(i,j,kmu)
            cdrag = (von_karman/log(dist/ruff))**2  ! drag coefficient from currents at height=dist above sea bed
            u_vel = Velocity%u(i,j,kmu,1,tstep)
            v_vel = Velocity%u(i,j,kmu,2,tstep)
            bot_vel=sqrt(u_vel*u_vel + v_vel*v_vel)

            ! find the effective bottom roughness rauh seen by bot_vel above the wave-boundary layer
            ! determined by matching log-profiles of outer and wave-boundary layer, see Kuhrts et al. (2004) Eq. (1,2)
            ! Waves%-quantities are defined at the t-grid. 
                
            ! put wave_u onto the u-grid
            wave_u_u = max(wave_u(i,j),   wave_u(i+1,j), &
                           wave_u(i,j+1), wave_u(i+1,j+1))
            
            if (wave_u_u > epsln) then

                ! put wave_p onto the u-grid
                wave_p_u = max(Waves%wave_p(i,j),   Waves%wave_p(i+1,j), &
                               Waves%wave_p(i,j+1), Waves%wave_p(i+1,j+1))
                ! wave orbital frequency derived from peak frequency [1/s]
                omega=twopi*wave_p_u 
 
                ! iterate now the bottom drag above the wave-boundary layer (Grand and Madsen)

                cnew = cdrag ! use cdrag as first guess                                    ! [dimensionless]
                do icount=1,itmax
                   cold  = cnew                                                             ! [dimensionless]

                   ! calculate current friction velocity with effective bottom drag 
                   ustar  = bot_vel*sqrt(cold)                                              ! [m/s]

                   ! add wave friction to current friction velocity assuming same direction
                   ucomb = sqrt(ustar*ustar + wave_u_u*wave_u_u)                            ! [m/s]

                   ! determine effective bottom roughness rauh from physical bottom roughness ruff
                   rauh  = ruff*(2*von_karman*ucomb/(omega*ruff+epsln))**&                   ! [m]
                        (1.0-ustar/(ucomb+epsln))

                   ! calculate bottom drag above the wave-boundary layer
                   cnew = (von_karman/log(dist/rauh))**2                                     ! [dimensionless]
                   if (id_iter > 0) wrk2_2d(i,j) = icount
                   if (abs(cold-cnew) .le. 0.01*abs(cnew)) exit

                     !TS iteration cnew(iter) better than cnew=cdrag
                     !!! in case of bad convergence, but this must never happen
                     !!            if (icount.eq.itmax) cnew = cdrag

                enddo


            else ! no waves 
                cnew=cdrag
            endif

            ! the effective drag coefficient seen by model velocities above the wave  boundary layer 
            drag_coeff(i,j) = cnew * rho0

            ! the combined current-wave induced stress in wave  boundary layer 
            ustar2 = bot_vel*bot_vel*cnew                                                ! [m/s]
            ucomb = sqrt(ustar2 + wave_u_u*wave_u_u)                                     ! [m/s]


            ! calculate shear stress acting on grains at sediment surface   
            wrk1_2d(i,j) = (ustar2/ucomb)*0.3152

            ! the wave induced friction velocity wave_s is derived applying Nielsen (1992) for grain roughness
            ! store the total wave-current bottom stress in skin friction layer for later use in sediment dynamics

        endif ! kmu.ge.1

     enddo
  enddo
  call mpp_update_domains(wrk1_2d(:,:), drag_coeff(:,:), Dom%domain2d) 

  ! place onto the T-grid
  do j=jsc,jec 
     do i=isc,iec 
        ucskin = max(wrk1_2d(i,j  ), wrk1_2d(i-1,j  ),  &
                     wrk1_2d(i,j-1), wrk1_2d(i-1,j-1))
        Velocity%current_wave_stress(i,j) = rho0*(ucskin*ucskin + wave_s(i,j)*wave_s(i,j))
     enddo
  enddo

  ! diagnostics 
  call diagnose_2d(Time, Grd, id_cur_wav_dr, Velocity%current_wave_stress(:,:))
  call diagnose_2d(Time, Grd, id_wave_s, wave_s(:,:))
  call diagnose_2d(Time, Grd, id_wave_u, wave_u(:,:))
  call diagnose_2d_u(Time, Grd, id_iter, wrk2_2d(:,:))

  if(debug_this_module) then 
      write(stdoutunit,*) 'ocean_wave_model: end wave_drag_diag'
  endif 


end subroutine current_wave_drag_diag
! </SUBROUTINE> NAME="current_wave_drag_diag"


!#######################################################################
! <SUBROUTINE NAME="wave_u_diag">
!
! <DESCRIPTION>
! calculates wave bottom shear stress velocity 
! wave friction factor is parametrized by approximation of 
! Nielsen (1992), Coastal bottom boundary layers and sediment transport!
!
! April 2012
! martin.schmidt@io-warnemuende.de 
!
! </DESCRIPTION>
!
subroutine wave_u_diag(wave_u, wave_s, Waves)
  real, dimension(isd:,jsd:), intent(inout) :: wave_u
  real, dimension(isd:,jsd:), intent(inout) :: wave_s
  type(ocean_wave_type),      intent(inout) :: Waves

  real    :: omega, U_m, ampli, f_w, wk_times_depth, ruff, grain
  integer :: i, j
  real,parameter:: twopi=2.*pi

  ! these fields have units m/s
  wave_u=0.0
  wave_s=0.0

  do j=jsd,jed 
     do i=isd,ied 

        ! water_depth*wave number [dimensionless]
        wk_times_depth = Grd%ht(i,j) * Waves%wave_k(i,j) 

        !TS uncleasr why wk_times_depth .gt. 0.3 ??? possibly to exclude land ht=0
        !!      if (wk_times_depth .gt. 0.3 .and. wk_times_depth.le.50.) then

        if (wk_times_depth.le.50.) then
            ruff = roughness_length(i,j)                                 ! bottom roughness length (grain+form drag) [m]
            omega=twopi*Waves%wave_p(i,j)                                ! wave orbital frequency [1/s]
            ampli=0.5*Waves%height(i,j)/(sinh(wk_times_depth)+epsln)     ! wave amplitude at sea bottom [m]
            U_m=ampli*omega                                              ! wave orbital velocity at sea bottom [m/s]

            !Kuhrts et al. (2004) Eq. (3) fw >= 0.3 if ampli >= 43.9315 * ruff 
            if (ampli.ge.(43.9315*ruff)) then
                f_w= exp(5.5*(30.*ruff/(ampli+epsln))**0.2 - 6.3)          ! [dimensionless]
            else
                f_w= 0.3                                                   ! [dimensionless]
            endif
            wave_u(i,j)=sqrt(0.5*f_w)*U_m                          ! [m/s]

            !Kuhrts et al. (2004) Eq. (7) applied for grain roughness
            grain = ruff/33.0  
            if (ampli.ge.(43.9315*grain)) then
                f_w= exp(5.5*(30.*grain/(ampli+epsln))**0.2 - 6.3)         ! [dimensionless]
            else
                f_w= 0.3                                                   ! [dimensionless]
            endif
            wave_s(i,j)=sqrt(0.5*f_w)*U_m                                ! [m/s]
        endif
     enddo
  enddo

  if(debug_this_module) then 
      write(stdoutunit,*) 'ocean_wave_model: end wave_u_diag'
  endif


end subroutine wave_u_diag
! </SUBROUTINE> NAME="wave_u_diag"

end module ocean_bbc_mod



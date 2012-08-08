! ============================================================================
! snow model module
! ============================================================================
module snow_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : error_mesg, file_exist, check_nml_error, &
     stdlog, write_version_number, close_file, mpp_pe, mpp_root_pe, FATAL, NOTE
use time_manager_mod,   only: time_type, increment_time, time_type_to_real
use constants_mod,      only: tfreeze, hlv, hlf, PI

use land_constants_mod, only : NBANDS
use snow_tile_mod, only : &
     snow_tile_type, snow_prog_type, read_snow_data_namelist, &
     snow_data_thermodynamics, snow_data_area, snow_data_radiation, snow_data_diffusion, &
     snow_data_hydraulics, max_lev, cpw, clw, csw

use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=)
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use land_data_mod,      only : land_state_type, lnd
use land_tile_io_mod, only : create_tile_out_file, read_tile_data_r1d_fptr, &
     write_tile_data_r1d_fptr, print_netcdf_error, get_input_restart_name, &
     sync_nc_files
use nf_utils_mod, only : nfu_def_dim, nfu_put_att
use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: read_snow_namelist
public :: snow_init
public :: snow_end
public :: save_snow_restart
public :: snow_get_sfc_temp
public :: snow_get_depth_area
public :: snow_radiation
public :: snow_diffusion
public :: snow_step_1
public :: snow_step_2
! =====end of public interfaces ==============================================


! ==== module variables ======================================================
character(len=*), parameter, private   :: &
       module_name = 'snow_mod' ,&
       version     = '$Id: snow.F90,v 19.0 2012/01/06 20:42:42 fms Exp $' ,&
       tagname     = '$Name: siena_201207 $'

! ==== module variables ======================================================

!---- namelist ---------------------------------------------------------------
logical :: retro_heat_capacity  = .false.
logical :: lm2  = .false.
logical :: steal = .false.
character(len=16):: albedo_to_use = ''  ! or 'brdf-params'
real    :: max_snow             = 1000.
real    :: wet_max              = 0.0  ! TEMP, move to snow_data
real    :: snow_density         = 300. ! TEMP, move to snow_data and generalize
   real :: init_temp = 260.   ! cold-start snow T
   real :: init_pack_ws   =   0.  
   real :: init_pack_wl   =   0.
   real :: min_snow_mass = 0.

namelist /snow_nml/ retro_heat_capacity, lm2, steal, albedo_to_use, &
                    max_snow, wet_max, snow_density, &
                    init_temp, init_pack_ws, init_pack_wl, &
                    min_snow_mass
!---- end of namelist --------------------------------------------------------

logical         :: module_is_initialized =.FALSE.
logical         :: use_brdf
type(time_type) :: time
real            :: delta_time
integer         :: num_l    ! # of snow layers
! next three 'z' variables are all normalized by total snow pack depth
real            :: dz (max_lev) ! relative thicknesses of layers
real            :: z  (max_lev) ! relative depths of layer bounds
real            :: zz (max_lev) ! relative depths of layer centers
real            :: heat_capacity_retro = 1.6e6
real            :: mc_fict

! ==== end of module variables ===============================================

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains

! ============================================================================
subroutine read_snow_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call read_snow_data_namelist(num_l,dz,mc_fict)

  call write_version_number(version, tagname)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=snow_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snow_nml)
  endif

  ! -------- set up vertical discretization --------
  zz(1) = 0
  do l = 1, num_l
     zz(l+1) = zz(l) + dz(l)
     z(l)    = 0.5*(zz(l+1) + zz(l))
  enddo

end subroutine read_snow_namelist


! ============================================================================
! initialize snow model
subroutine snow_init ( id_lon, id_lat )
  integer, intent(in)               :: id_lon  ! ID of land longitude (X) axis  
  integer, intent(in)               :: id_lat  ! ID of land latitude (Y) axis

  ! ---- local vars ----------------------------------------------------------
  integer :: unit         ! unit for various i/o
  type(land_tile_enum_type)     :: te,ce ! tail and current tile list elements
  type(land_tile_type), pointer :: tile  ! pointer to current tile
  character(len=256) :: restart_file_name
  logical :: restart_exists

  module_is_initialized = .TRUE.
  time       = lnd%time
  delta_time = time_type_to_real(lnd%dt_fast)

  ! -------- initialize snow state --------
  call get_input_restart_name('INPUT/snow.res.nc',restart_exists,restart_file_name)
  if (restart_exists) then
     call error_mesg('snow_init',&
          'reading NetCDF restart "'//trim(restart_file_name)//'"',&
          NOTE)
     __NF_ASRT__(nf_open(restart_file_name,NF_NOWRITE,unit))
     call read_tile_data_r1d_fptr(unit, 'temp', snow_temp_ptr  )
     call read_tile_data_r1d_fptr(unit, 'wl'  , snow_wl_ptr )
     call read_tile_data_r1d_fptr(unit, 'ws'  , snow_ws_ptr )
     __NF_ASRT__(nf_close(unit))     
  else
     call error_mesg('snow_init',&
          'cold-starting snow',&
          NOTE)
     te = tail_elmt (lnd%tile_map)
     ce = first_elmt(lnd%tile_map)
     do while(ce /= te)
        tile=>current_tile(ce)  ! get pointer to current tile
        ce=next_elmt(ce)       ! advance position to the next tile
        
        if (.not.associated(tile%snow)) cycle
     
        tile%snow%prog(1:num_l)%wl = init_pack_wl * dz(1:num_l)
        tile%snow%prog(1:num_l)%ws = init_pack_wl * dz(1:num_l)
        tile%snow%prog(1:num_l)%T  = init_temp
     enddo
  endif

  if (trim(albedo_to_use)=='') then
     use_brdf = .false.
  elseif (trim(albedo_to_use)=='brdf-params') then
     use_brdf = .true.
  else
     call error_mesg('snow_init',&
          'option albedo_to_use="'//&
          trim(albedo_to_use)//'" is invalid, use "" or "brdf-params"',&
          FATAL)
  endif

end subroutine snow_init


! ============================================================================
subroutine snow_end ()

  module_is_initialized =.FALSE.

end subroutine snow_end


! ============================================================================
subroutine save_snow_restart (tile_dim_length, timestamp)
  integer, intent(in) :: tile_dim_length ! length of tile dim. in the output file
  character(*), intent(in) :: timestamp ! timestamp to add to the file name

  ! ---- local vars ----------------------------------------------------------
  integer :: unit            ! restart file i/o unit

  call error_mesg('snow_end','writing NetCDF restart',NOTE)
  ! create output file, including internal structure necessary for tile output
  call create_tile_out_file(unit,'RESTART/'//trim(timestamp)//'snow.res.nc', &
          lnd%coord_glon, lnd%coord_glat, snow_tile_exists, tile_dim_length )

  ! additionally, define vertical coordinate
  if (mpp_pe()==lnd%io_pelist(1)) then
     __NF_ASRT__(nfu_def_dim(unit,'zfull',zz(1:num_l),'depth of level centers'))
     __NF_ASRT__(nfu_put_att(unit,'zfull','positive','down'))
  endif
  call sync_nc_files(unit)

  ! write fields
  call write_tile_data_r1d_fptr(unit,'temp',snow_temp_ptr,'zfull','snow temperature','degrees_K')
  call write_tile_data_r1d_fptr(unit,'wl'  ,snow_wl_ptr,  'zfull','snow liquid water content','kg/m2')
  call write_tile_data_r1d_fptr(unit,'ws'  ,snow_ws_ptr,  'zfull','snow solid water content','kg/m2')
  ! close output file
  if (mpp_pe()==lnd%io_pelist(1)) &
       __NF_ASRT__(nf_close(unit))

end subroutine save_snow_restart

! ============================================================================
subroutine snow_get_sfc_temp(snow, snow_T)
  type(snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_T
  
  snow_T = snow%prog(1)%T
end subroutine


! ============================================================================
subroutine snow_get_depth_area(snow, snow_depth, snow_area)
  type(snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_depth, snow_area

  integer :: l

  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%prog(l)%ws
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )
end subroutine


! ============================================================================
! compute snow properties needed to do soil-canopy-atmos energy balance
subroutine snow_radiation ( snow_T, cosz, &
     snow_refl_dir, snow_refl_dif, snow_refl_lw, snow_emis )
  real, intent(in) :: snow_T  ! snow temperature, deg K
  real, intent(in) :: cosz ! cosine of zenith angle
  real, intent(out) :: snow_refl_dir(NBANDS), snow_refl_dif(NBANDS), snow_refl_lw, snow_emis

  call snow_data_radiation (snow_T, snow_refl_dir, snow_refl_dif, &
                                snow_emis, cosz, use_brdf )
  snow_refl_lw = 1 - snow_emis
end subroutine 


! ============================================================================
! compute snow properties needed to do soil-canopy-atmos energy balance
subroutine snow_diffusion ( snow, snow_z0s, snow_z0m )
  type(snow_tile_type), intent(in) :: snow
  real, intent(out) :: snow_z0s, snow_z0m

  call snow_data_diffusion ( snow_z0s, snow_z0m )
end subroutine


! ============================================================================
! update snow properties explicitly for time step.
! integrate snow-heat conduction equation upward from bottom of snow
! to surface, delivering linearization of surface ground heat flux.
subroutine snow_step_1 ( snow, snow_G_Z, snow_G_TZ, &
                         snow_active, snow_T, snow_rh, snow_liq, snow_ice, &
                         snow_subl, snow_area, snow_G0, snow_DGDT )
  type(snow_tile_type), intent(inout) :: snow
  real,                 intent(in) :: snow_G_Z
  real,                 intent(in) :: snow_G_TZ
  logical,              intent(out):: snow_active
  real,                 intent(out):: &
       snow_T, snow_rh, snow_liq, snow_ice, &
       snow_subl, snow_area, snow_G0, snow_DGDT

  ! ---- local vars
  real :: snow_depth, bbb, denom, dt_e
  real, dimension(num_l):: aaa, ccc, thermal_cond, dz_phys, heat_capacity
  integer :: l

! ----------------------------------------------------------------------------
! in preparation for implicit energy balance, determine various measures
! of water availability, so that vapor fluxes will not exceed mass limits
! ----------------------------------------------------------------------------

  snow_T = tfreeze
  snow_T = snow%prog(1)%T

  call snow_data_thermodynamics ( snow_rh, thermal_cond )
  snow_depth= 0.0
  do l = 1, num_l
     snow_depth = snow_depth + snow%prog(l)%ws
  enddo
  snow_depth = snow_depth / snow_density
  call snow_data_area (snow_depth, snow_area )
  ! ---- only liquid in the top snow layer is available to freeze implicitly
  snow_liq =     snow%prog(1)%wl
  ! ---- snow in any layer can be melted implicitly
  snow_ice = sum(snow%prog(:)%ws)

! ---- fractionate evaporation/sublimation according to sfc phase ratios
!  where (max(snow%prog(1)%ws,0.)+max(snow%prog(1)%wl,0.)>0)
!      snow_subl = max(snow%prog(1)%ws,0.) &
!       /(max(snow%prog(1)%ws,0.)+max(snow%prog(1)%wl,0.))
!    elsewhere
!      snow_subl = 0
!    endwhere
!  snow_active = snow_subl>0.
  if (snow_depth>0) then
     snow_subl = 1.
  else
     snow_subl = 0
  endif
  snow_active = snow_subl>0.

  do l = 1, num_l
     dz_phys(l) = dz(l)*snow_depth
  enddo

  if (retro_heat_capacity) then
     do l = 1, num_l
        heat_capacity(l) = heat_capacity_retro*dz_phys(l)
     enddo
  else
     do l = 1, num_l
        heat_capacity(l) = mc_fict*dz(l) + &
             clw*snow%prog(l)%wl + csw*snow%prog(l)%ws
     enddo
  endif

!  if(num_l > 1) then
  if (snow_depth > 0) then
     do l = 1, num_l-1
        dt_e = 2 / ( dz_phys(l+1)/thermal_cond(l+1) &
                     + dz_phys(l)/thermal_cond(l)   )
        aaa(l+1) = - dt_e * delta_time / heat_capacity(l+1)
        ccc(l)   = - dt_e * delta_time / heat_capacity(l)
     enddo

     bbb = 1.0 - aaa(num_l) + delta_time*snow_G_TZ/heat_capacity(num_l)
     denom = bbb
     dt_e = aaa(num_l)*(snow%prog(num_l)%T - snow%prog(num_l-1)%T) &
          - delta_time*snow_G_Z/heat_capacity(num_l)
     snow%e(num_l-1) = -aaa(num_l)/denom
     snow%f(num_l-1) = dt_e/denom

     do l = num_l-1, 2, -1
        bbb = 1.0 - aaa(l) - ccc(l)
        denom = bbb + ccc(l)*snow%e(l)
        dt_e = - ( ccc(l)*(snow%prog(l+1)%T - snow%prog(l)%T  ) &
                  -aaa(l)*(snow%prog(l)%T   - snow%prog(l-1)%T) )
        snow%e(l-1) = -aaa(l)/denom
        snow%f(l-1) = (dt_e - ccc(l)*snow%f(l))/denom
     enddo

     denom = delta_time/heat_capacity(1)
     snow_G0    = ccc(1)*(snow%prog(2)%T- snow%prog(1)%T &
          + snow%f(1)) / denom
     snow_DGDT  = (1 - ccc(1)*(1-snow%e(1))) / denom    
  endif

!    else  ! one-level case
!      denom = delta_time/heat_capacity(1)
!      snow_G0    = 0.
!      snow_DGDT  = 1. / denom
!    end if

  if (snow_depth <= 0) then
     snow_G0   = snow_G_Z
     snow_DGDT = snow_G_TZ
  endif

  if(is_watch_point()) then
     write(*,*) 'snow_depth', snow_depth 
     write(*,*) '############ snow_step_1 output'
     write(*,*) 'mask      ', .true.
     write(*,*) 'snow_T    ', snow_T     
     write(*,*) 'snow_rh   ', snow_rh    
     write(*,*) 'snow_liq  ', snow_liq   
     write(*,*) 'snow_ice  ', snow_ice   
     write(*,*) 'snow_subl ', snow_subl  
     write(*,*) 'snow_area ', snow_area  
     write(*,*) 'snow_G_Z  ', snow_G_Z   
     write(*,*) 'snow_G_TZ ', snow_G_TZ  
     write(*,*) 'snow_G0   ', snow_G0    
     write(*,*) 'snow_DGDT ', snow_DGDT  
     write(*,*) '############ end of snow_step_1 output'
  endif

end subroutine snow_step_1



! ============================================================================
! apply boundary flows to snow water and move snow water vertically.
  subroutine snow_step_2 ( snow, snow_subl,                     &
                           vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec, &
                           DTg,  Mg_imp,  evapg,  fswg,  flwg,  sensg,  &
                           use_tfreeze_in_grnd_latent, subs_DT, &
                           subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens,  &
                           snow_fsw, snow_flw, snow_sens, &
                           snow_levap, snow_fevap, snow_melt, &
                           snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
                           snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, &
                           snow_avrg_T )
  type(snow_tile_type), intent(inout) :: snow
  real, intent(in) :: &
     snow_subl, vegn_lprec, vegn_fprec, vegn_hlprec, vegn_hfprec
  real, intent(in) :: &
     DTg, Mg_imp, evapg, fswg, flwg, sensg
  logical, intent(in) :: use_tfreeze_in_grnd_latent
  real, intent(out) :: &
         subs_DT, subs_M_imp, subs_evap, subs_fsw, subs_flw, subs_sens, &
         snow_fsw, snow_flw, snow_sens, &
         snow_levap, snow_fevap, snow_melt, &
         snow_lprec, snow_hlprec, snow_lrunf, snow_frunf, &
         snow_hlrunf, snow_hfrunf, snow_Tbot, snow_Cbot, snow_C, snow_avrg_T

  ! ---- local vars
  real, dimension(num_l) :: del_t, M_layer
  real :: depth, &
         cap0, dW_l, dW_s, dcap, dheat,&
         melt, melt_per_deg, drain,       &
         snow_mass, sum_liq, &
         sum_heat, sum_sno, &
         snow_transfer, frac,&
         liq_rate, hliq_rate,&
         sno_rate, hsno_rate, fict_heat, &
         evapg_lm2, vegn_fprec_lm2, &
         snow_LMASS, snow_FMASS, snow_HEAT
  integer :: l, l_old
  type(snow_prog_type) :: new_prog(num_l)
  ! --------------------------------------------------------------------------

  depth= 0.
  do l = 1, num_l
    depth = depth + snow%prog(l)%ws
  enddo
  depth = depth / snow_density

  if(is_watch_point()) then
     write(*,*) '############ snow_step_2 input'
     write(*,*) 'mask       ', .TRUE.       
     write(*,*) 'snow_subl  ', snow_subl  
     write(*,*) 'vegn_lprec ', vegn_lprec 
     write(*,*) 'vegn_fprec ', vegn_fprec 
     write(*,*) 'vegn_hlprec', vegn_hlprec
     write(*,*) 'vegn_hfprec', vegn_hfprec
     write(*,*) 'DTg        ', DTg        
     write(*,*) 'Mg_imp     ', Mg_imp     
     write(*,*) 'evapg      ', evapg      
     write(*,*) 'fswg       ', fswg       
     write(*,*) 'flwg       ', flwg       
     write(*,*) 'sensg      ', sensg      
     write(*,*) '############ end of snow_step_2 input'

     write(*,*) 'depth   ', depth
     do l = 1, num_l
        write(*,'(i2,3(x,a,g))') l,&
             ' wl=', snow%prog(l)%wl,&
             ' ws=', snow%prog(l)%ws,&
             ' T =', snow%prog(l)%T
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l;
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 1.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- record fluxes -------------------------------------------------------
  if (lm2.and.steal) then
    if (snow_FMASS-Mg_imp > 0.) then
      if (evapg <= (snow_FMASS-Mg_imp)/delta_time) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec
      else if (evapg <= (snow_FMASS-Mg_imp)/delta_time+vegn_fprec) then
          evapg_lm2 = evapg
          vegn_fprec_lm2 = vegn_fprec - evapg + (snow_FMASS-Mg_imp)/delta_time
      else
          evapg_lm2 = (snow_FMASS-Mg_imp)/delta_time+vegn_fprec
          vegn_fprec_lm2 = 0.
      endif
    else
      evapg_lm2 = 0.
      vegn_fprec_lm2 = vegn_fprec
    endif
  else
     evapg_lm2 = evapg
     vegn_fprec_lm2 = vegn_fprec
  endif
  vegn_fprec_lm2 = vegn_fprec
  if (depth>0) then
        snow_fsw   = fswg
        snow_flw   = flwg
        snow_sens  = sensg
        snow_levap = evapg_lm2*(1-snow_subl)
        snow_fevap = evapg_lm2*   snow_subl
  else
        snow_fsw    = 0
        snow_flw    = 0
        snow_sens   = 0
        snow_levap  = 0
        snow_fevap  = 0
  endif
  subs_fsw = fswg - snow_fsw
  subs_flw = flwg - snow_flw
  subs_evap = evapg - snow_levap - snow_fevap
  subs_sens = sensg - snow_sens 

  ! ---- load surface temp change and perform back substitution --------------
  if (depth>0) then
      del_t(1) = DTg
      snow%prog(1)%T  = snow%prog(1)%T + del_t(1)
  endif
  if ( num_l > 1) then
     do l = 1, num_l-1
        if (depth>0) then
            del_t(l+1) = snow%e(l) * del_t(l) + snow%f(l)
            snow%prog(l+1)%T = snow%prog(l+1)%T + del_t(l+1)
        endif
     enddo
  endif
  if (depth>0) then
    subs_DT = del_t(num_l)
  else
    subs_DT = DTg
  endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2 ***** '
     do l = 1, num_l
        write(*,'(i2,a,g)') l,' T =', snow%prog(l)%T
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- evaporation and sublimation -----------------------------------------
  if (depth>0) then
      snow%prog(1)%wl = snow%prog(1)%wl - snow_levap*delta_time
      snow%prog(1)%ws = snow%prog(1)%ws - snow_fevap*delta_time
      cap0 = mc_fict*dz(1) + clw*snow%prog(1)%wl + csw*snow%prog(1)%ws
      ! T adjustment for nonlinear terms (del_T)*(del_W)
      dheat = delta_time*(clw*snow_levap+csw*snow_fevap)*del_T(1)
      ! take out extra heat not claimed in advance for evaporation
      if (use_tfreeze_in_grnd_latent) dheat = dheat &
            - delta_time*((cpw-clw)*snow_levap+(cpw-csw)*snow_fevap) &
                               *(snow%prog(1)%T-del_T(1)-tfreeze)
      snow%prog(1)%T  = snow%prog(1)%T  + dheat/cap0
    endif

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.5 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g))')l,&
             ' wl=', snow%prog(l)%wl,&
             ' ws=', snow%prog(l)%ws,&
             ' T =', snow%prog(l)%T
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 2.51 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- distribute implicit phase change downward through snow layers -------
  if (depth>0) then 
      snow_melt = Mg_imp/delta_time
  else
      snow_melt = 0
  endif
  M_layer = 0.
  subs_M_imp = Mg_imp
  do l = 1, num_l
    if (depth>0 .and. subs_M_imp.gt.0) then
        M_layer(l) =  min( subs_M_imp, max(0.,snow%prog(l)%ws) )
        subs_M_imp = subs_M_imp - M_layer(l)
    endif
  enddo
  if (depth>0) then
      M_layer(1) = M_layer(1) + subs_M_imp
      subs_M_imp = 0.
  endif
  do l = 1, num_l
    if (depth>0) then
          cap0 = mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws
          snow%prog(l)%wl = snow%prog(l)%wl + M_layer(l)
          snow%prog(l)%ws = snow%prog(l)%ws - M_layer(l)
          snow%prog(l)%T  = tfreeze + (cap0*(snow%prog(l)%T-tfreeze) ) &
                                                          / ( cap0 + (clw-csw)*M_layer(l) )
    endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 3 ***** '
     do l = 1, num_l
        write(*,'(i2,3(a,g))') l,&
             ' wl=', snow%prog(l)%wl,&
             ' ws=', snow%prog(l)%ws,&
             '  T=', snow%prog(l)%T
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
    snow_LMASS = snow_LMASS + snow%prog(l)%wl
    snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 3.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

! ----------------------------------------------------------------------------
!  call snow_data_hydraulics (pars, snow%prog%wl, psi, hyd_cond )

! ---- remainder of mass fluxes and associated sensible heat fluxes ----------
  liq_rate = vegn_lprec
  sno_rate = vegn_fprec_lm2
  hliq_rate = vegn_hlprec
  if (vegn_fprec.ne.0.) then
          hsno_rate = vegn_hfprec*(vegn_fprec_lm2/vegn_fprec)
  else
          hsno_rate = 0.
  endif

  do l = 1, num_l
    if(depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- mix inflow with existing snow and water ---------------------------
          cap0 = mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws
          dW_l = liq_rate*delta_time
          dW_s = sno_rate*delta_time
          dcap = clw*dW_l + csw*dW_s
          snow%prog(l)%ws = snow%prog(l)%ws + dW_s
          snow%prog(l)%wl = snow%prog(l)%wl + dW_l
          snow%prog(l)%T  = tfreeze &
                                        + (cap0*(snow%prog(l)%T-tfreeze) &
                                        + (hsno_rate+hliq_rate)*delta_time) /(cap0 + dcap)
    endif

    if(is_watch_point()) then
       write(*,*) ' ***** snow_step_2 checkpoint 4a ***** '
       write(*,'(i2,3(a,g))') l,&
            ' wl=', snow%prog(l)%wl,&
            ' ws=', snow%prog(l)%ws,&
            '  T=', snow%prog(l)%T
    endif

    if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute explicit melt/freeze --------------------------------------
          melt_per_deg = (cap0+dcap)/hlf
          if (snow%prog(l)%ws>0 .and. snow%prog(l)%T>tfreeze) then
                  melt =  min(snow%prog(l)%ws, (snow%prog(l)%T-tfreeze)*melt_per_deg)
      elseif (snow%prog(l)%wl>0 .and. snow%prog(l)%T<tfreeze) then
                  melt = -min(snow%prog(l)%wl, (tfreeze-snow%prog(l)%T)*melt_per_deg)
      else
                  melt = 0
          endif
          snow_melt = snow_melt + melt/delta_time
          snow%prog(l)%wl = snow%prog(l)%wl + melt
          snow%prog(l)%ws = snow%prog(l)%ws - melt
!        where (cap0+dcap.ne.0.) &
!        snow%prog(l)%T  = snow%prog(l)%T  - melt/melt_per_deg
          snow%prog(l)%T = tfreeze &
                 + ((cap0+dcap)*(snow%prog(l)%T-tfreeze) - hlf*melt) &
                                                          / ( cap0+dcap + (clw-csw)*melt )
    endif

   if(is_watch_point()) then
      write(*,*) ' ***** snow_step_2 checkpoint 4b ***** '
      write(*,'(i2,3(a,g))')l,&
            ' wl=', snow%prog(l)%wl,&
            ' ws=', snow%prog(l)%ws,&
            '  T=', snow%prog(l)%T
   endif

   if (depth>0 .or. vegn_fprec_lm2>0) then
    ! ---- compute drainage from this layer to next --------------------------
        drain = max (0., snow%prog(l)%wl - wet_max*snow%prog(l)%ws)
        snow%prog(l)%wl = snow%prog(l)%wl - drain
        liq_rate = drain / delta_time
        hliq_rate = clw*liq_rate*(snow%prog(l)%T-tfreeze)
        sno_rate = 0
        hsno_rate = 0
    endif
  enddo

  snow_lprec  = liq_rate
  snow_hlprec = hliq_rate

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo


  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4c ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

! ---- conceptually remove fictitious mass/heat for the moment ---------------
  fict_heat = 0.
  do l = 1, num_l
    fict_heat = fict_heat + dz(l)*snow%prog(l)%T     ! (*mc_fict)
  enddo

  snow_mass  = sum(snow%prog%ws)
  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 4d ***** '
     write(*,*) 'max_snow    ', max_snow
     write(*,*) 'snow_mass   ', snow_mass
  endif

! ---- remove any isolated snow molecules (!) or sweep any excess snow from top of pack ----

  snow_lrunf  = 0.
  snow_frunf  = 0.
  snow_hlrunf = 0.
  snow_hfrunf = 0.
  if (0. < snow_mass .and. snow_mass < min_snow_mass ) then
        do l = 1, num_l
          snow_hlrunf = snow_hlrunf  &
            + clw*snow%prog(l)%wl*(snow%prog(l)%T-tfreeze)
          snow_hfrunf = snow_hfrunf  &
            + csw*snow%prog(l)%ws*(snow%prog(l)%T-tfreeze)
          enddo
        snow_lrunf  = sum(snow%prog%wl)
        snow_frunf  = snow_mass
        snow_mass   = 0.
        snow%prog%ws = 0.
        snow%prog%wl = 0.
    else if (max_snow < snow_mass) then
        snow_frunf  = snow_mass - max_snow
        snow_mass  = max_snow
        sum_sno  = 0
        snow_transfer = 0
        do l = 1, num_l
          if (sum_sno + snow%prog(l)%ws > snow_frunf) then
              snow_transfer = snow_frunf - sum_sno
            else
              snow_transfer = snow%prog(l)%ws
            endif
          if (snow%prog(l)%ws > 0) then
              frac = snow_transfer / snow%prog(l)%ws
            else
              frac = 1.
            endif
          sum_sno  = sum_sno  + snow_transfer
          snow_lrunf  = snow_lrunf  +     frac*snow%prog(l)%wl
          snow_hlrunf = snow_hlrunf + clw*frac*snow%prog(l)%wl*(snow%prog(l)%T-tfreeze)
          snow_hfrunf = snow_hfrunf + csw*frac*snow%prog(l)%ws*(snow%prog(l)%T-tfreeze)
          snow%prog(l)%ws = (1-frac)*snow%prog(l)%ws
          snow%prog(l)%wl = (1-frac)*snow%prog(l)%wl
          enddo
    endif
  snow_lrunf  = snow_lrunf  / delta_time
  snow_frunf  = snow_frunf  / delta_time
  snow_hlrunf = snow_hlrunf / delta_time
  snow_hfrunf = snow_hfrunf / delta_time

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5 ***** '
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g))')l,&
             ' wl=', snow%prog(l)%wl,&
             ' ws=', snow%prog(l)%ws,&
             ' T =', snow%prog(l)%T 
     enddo
  endif

  snow_LMASS = 0; snow_FMASS = 0; snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                                                * (snow%prog(l)%T-tfreeze)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.01 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  depth= 0.
  new_prog%ws=0
  new_prog%wl=0
  new_prog%T=0
  do l = 1, num_l
    depth = depth + snow%prog(l)%ws
  enddo
  depth = depth / snow_density

!************************** fudge to avoid T=NaN from too-small mass **
!   if(depth*snow_density < min_snow_mass .and. depth>0.) then
!       depth = 0
!       snow%prog%ws = 0
!       snow%prog%wl = 0
!     endif

! ---- re-layer the snowpack ------------------------------------------------
  do l = 1, num_l
     if (depth > 0) then
        new_prog(l)%ws = snow_mass*dz(l)
        sum_sno = 0
        sum_liq = 0
        sum_heat = 0
     endif
     do l_old = 1, num_l
        if (depth > 0) then
           if (sum_sno + snow%prog(l_old)%ws > new_prog(l)%ws) then
              snow_transfer = new_prog(l)%ws - sum_sno
           else
              snow_transfer = snow%prog(l_old)%ws
           endif
           if (snow%prog(l_old)%ws .ne. 0.) then
              frac = snow_transfer / snow%prog(l_old)%ws
           else
              frac = 1
           endif
           sum_sno  = sum_sno  + snow_transfer
           sum_liq  = sum_liq  + frac*     snow%prog(l_old)%wl
           sum_heat = sum_heat + frac*&
                (clw*snow%prog(l_old)%wl + csw*snow%prog(l_old)%ws)&
                *snow%prog(l_old)%T
           snow%prog(l_old)%ws = (1.-frac)*snow%prog(l_old)%ws
           snow%prog(l_old)%wl = (1.-frac)*snow%prog(l_old)%wl
           if(is_watch_point()) then
              write(*,*) 'l=',l, ' l_old=',l_old,snow_transfer,frac,&
                   sum_sno,sum_liq, sum_heat
           endif
        endif
        
     enddo
     if (depth > 0) then
        new_prog(l)%wl = sum_liq
        new_prog(l)%T  = sum_heat &
             / (clw*new_prog(l)%wl + csw*new_prog(l)%ws)
     endif
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 5.1 ***** '
     write(*,*) 'depth             ', depth
     write(*,*) 'fict_heat         ', fict_heat
     do l = 1, num_l
        write(*,'(i2,3(a,g))')l,&
             ' new_wl=', new_prog(l)%wl,&
             ' new_ws=', new_prog(l)%ws,&
             ' new_T =', new_prog(l)%T 
     enddo
  endif

! add back fictional mass/heat
  do l = 1, num_l
    if (depth > 0) &
    new_prog(l)%T = ( &
    (clw*new_prog(l)%wl + csw*new_prog(l)%ws)*new_prog(l)%T  &
      + mc_fict*dz(l)*fict_heat ) &
      / (clw*new_prog(l)%wl + csw*new_prog(l)%ws + dz(l)*mc_fict)
  enddo
    
  do l = 1, num_l
!!    where (mask .and. snow_mass > 0) snow%prog(l) = new_prog(l)
    if (depth > 0) snow%prog(l) = new_prog(l)
  enddo

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 6 ***** '
     write(*,*) 'evap         ', subs_evap
     write(*,*) 'snow_lprec', snow_lprec
     write(*,*) 'depth        ', depth
     do l = 1, num_l
        write(*,'(i2,3(a,g))')l,&
             ' wl=', snow%prog(l)%wl,&
             ' ws=', snow%prog(l)%ws,&
             ' T =', snow%prog(l)%T 
     enddo
  endif

  snow_LMASS = 0
  snow_FMASS = 0
  snow_HEAT = 0
  do l = 1, num_l
        snow_LMASS = snow_LMASS + snow%prog(l)%wl
        snow_FMASS = snow_FMASS + snow%prog(l)%ws
        snow_HEAT = snow_HEAT + &
          (mc_fict*dz(l) + clw*snow%prog(l)%wl + csw*snow%prog(l)%ws)  &
                            * (snow%prog(l)%T-tfreeze)
  enddo
  snow_Tbot = snow%prog(num_l)%T
  snow_Cbot = mc_fict*dz(num_l) &
        + clw*snow%prog(num_l)%wl + csw*snow%prog(num_l)%ws
  snow_C = sum(mc_fict*dz(1:num_l) &
        + clw*snow%prog(1:num_l)%wl + csw*snow%prog(1:num_l)%ws)
  snow_avrg_T = snow_HEAT/snow_C+tfreeze

  if(is_watch_point()) then
     write(*,*) ' ***** snow_step_2 checkpoint 7 ***** '
     write(*,*) 'LMASS         ', snow_LMASS
     write(*,*) 'FMASS         ', snow_FMASS
     write(*,*) 'HEAT          ', snow_HEAT
  endif

  ! ---- increment time and do diagnostics -----------------------------------
  time = increment_time(time, int(delta_time), 0)

end subroutine snow_step_2

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function snow_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   snow_tile_exists = associated(tile%snow)
end function snow_tile_exists

! ============================================================================
! accessor functions: given a pointer to a land tile, they return pointer
! to the desired member of the land tile, of NULL if this member does not
! exist.
subroutine snow_temp_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr=>tile%snow%prog%T
   endif
end subroutine snow_temp_ptr

subroutine snow_wl_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr=>tile%snow%prog%wl
   endif
end subroutine snow_wl_ptr

subroutine snow_ws_ptr(tile, ptr)
   type(land_tile_type), pointer :: tile
   real                , pointer :: ptr(:)
   ptr=>NULL()
   if(associated(tile)) then
      if(associated(tile%snow)) ptr=>tile%snow%prog%ws
   endif
end subroutine snow_ws_ptr

end module snow_mod




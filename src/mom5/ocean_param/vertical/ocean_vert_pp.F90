module ocean_vert_pp_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Vertical viscosity and diffusivity according Pacanowski and Philander (1981)
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes vertical viscosity and diffusivity according to 
! Pacanowski and Philander (1981).  This scheme is most effective for 
! studies of the tropical circulation.  It computes the vertical mixing
! coefficient based on the Richardson number.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and G. Philander 
! Parametrization of vertical mixing in numerical models of the tropical ocean
! Journal of Physical Oceanography (1981) vol 11, pages 1442--1451
! </REFERENCE>
!
! <NOTE>
! This parameterization was designed for equatorial models
! and may not do a good job in mid or high latitudes. Simulations
! in these regions (where vertical shear is small) are improved with
! the addition of solar short wave penetration into the ocean which 
! reduces buoyancy and enhances vertical mixing.
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_pp_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false. 
!  </DATA> 
!
!  <DATA NAME="wndmix" UNITS="m^2/sec" TYPE="real">
!  Minimum viscosity at bottom of 1st level to simulate 
!  missing high frequency windstress components.
!  </DATA> 
!  <DATA NAME="fricmx" UNITS="m^2/sec" TYPE="real">
!  Maximum mixing
!  </DATA> 
!  <DATA NAME="diff_cbt_back_pp" UNITS="m^2/sec" TYPE="real">
!  Space-time independent background vertical diffusivity 
!  thought to be that arising from internal waves. Note that 
!  if using Bryan-Lewis background diffusivity, then should 
!  set diff_cbt_back_pp=0.0. 
!  </DATA> 
!  <DATA NAME="visc_cbu_back_pp" UNITS="m^2/sec" TYPE="real">
!  Background vertical viscosity
!  </DATA> 
!</NAMELIST>

use constants_mod,       only: pi, epsln
use diag_manager_mod,    only: register_diag_field
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_mod,             only: input_nml_file, mpp_error

use ocean_density_mod,    only: density_delta_z
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, rho0r, grav 
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,      only: ocean_time_steps_type, ocean_time_type, ocean_thickness_type
use ocean_vert_util_mod,  only: ri_for_cgrid 
use ocean_workspace_mod,  only: wrk1, wrk1_v
use ocean_util_mod,       only: diagnose_3d, diagnose_3d_u

implicit none

public ocean_vert_pp_init
public vert_mix_pp
private ri_for_pp

private

real, private, dimension(:,:,:), allocatable :: riu  ! Richardson number at base of U-cells
real, private, dimension(:,:,:), allocatable :: rit  ! Richardson number at base of T-cells

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

integer :: index_temp=-1
integer :: index_salt=-1
integer :: num_prog_tracers=-1

! for Bgrid or Cgrid
integer :: horz_grid

character(len=128) :: version = &
     '$Id: ocean_vert_pp.F90,v 20.0 2013/12/14 00:16:50 fms Exp $'

character (len=128) :: tagname = &
     '$Name: tikal $'

! for diagnostics 
integer :: id_diff_cbt_pp= -1
integer :: id_visc_cbt_pp= -1
integer :: id_visc_cbu_pp= -1

logical :: module_is_initialized = .FALSE.

logical :: use_this_module = .false.
real :: fricmx           = 50.0e-4 ! (m^2/s) max vertical mixing coefficient
real :: wndmix           = 10.0e-4 ! (m^2/s) min vertical mixing in level 1 to simulate wind mixing
real :: diff_cbt_back_pp = 1.0e-5  ! (m^2/s) background "diff_cbt"
real :: visc_cbu_back_pp = 1.0e-4  ! (m^2/s) background "visc_cbu"
real :: visc_cbu_limit   = 1.0e2   ! (m^2/s) largest allowable "visc_cbu" (reset below)
real :: diff_cbt_limit   = 1.0e2   ! (m^2/s) largest allowable "diff_cbt" (reset below)

namelist /ocean_vert_pp_nml/ use_this_module, wndmix, fricmx,    &
                             diff_cbt_back_pp, visc_cbu_back_pp, &
                             visc_cbu_limit, diff_cbt_limit                     

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_pp_init">
!
! <DESCRIPTION>
! Initialization for the Pacanowski/Philander vertical mixing scheme
!
! input:
!
!   dzt    = thickness of vertical levels (m)
!
!   nk     = number of vertical levels
!
!   yt     = latitude of grid points (deg)
!
!   nj     = number of latitudes
!
!   error  = logical to signal problems
!
! output:
!
!   wndmix = min value for mixing at surface to simulate high freq
!
!            wind mixing (if absent in forcing). (m^2/sec)
!
!   fricmx = maximum mixing (m^2/sec)
!
!   diff_cbt_back_pp = background "diff_cbt" (m^2/sec)
!
!   visc_cbu_back_pp = background "visc_cbu" (m^2/sec)
!
!   diff_cbt_limit = largest "diff_cbt" (m^2/sec)
!
!   visc_cbu_limit = largest "visc_cbu" (m^2/sec)
!
!   error  = true if some inconsistency was found
!
! </DESCRIPTION>
!
subroutine ocean_vert_pp_init (Grid, Domain, Time, Time_steps, T_prog, hor_grid)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_time_steps_type),  intent(in)         :: Time_steps
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  integer,                      intent(in)         :: hor_grid

  real :: dzmin, extlat
  integer :: k, j, ioun, io_status, ierr, n
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_vert_pp_mod: module is already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_pp_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_pp_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_vert_pp_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_vert_pp_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_pp_nml)  
  write (stdlogunit, ocean_vert_pp_nml)

  write(stdoutunit,'(/a/)')'==>USING Pacanowski-Philander (pp) vertical mixing scheme.'
  write(stdoutunit,'(/a,f10.2)') &
  '==>Note from ocean_vert_pp_mod: using forward time step for vert-frict of (secs)', Time_steps%dtime_u 
  write(stdoutunit,'(/a,f10.2)') &
  '==>Note from ocean_vert_pp_mod: using forward time step for vert-diff  of (secs)', Time_steps%dtime_t 

  num_prog_tracers = size(T_prog)
  horz_grid        = hor_grid 

  do n = 1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_vert_pp_mod: temp and/or salt not present in tracer array')
  endif 

  Dom => Domain
  Grd => Grid

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(use_this_module) then 
    call mpp_error(NOTE, '==>Note: USING ocean_vert_pp_mod')
  else
    call mpp_error(NOTE, '==>Note: NOT using ocean_vert_pp_mod')
    return 
  endif


  allocate (riu(isd:ied,jsd:jed,nk))
  riu(:,:,:)    = 0.0

  allocate (rit(isd:ied,jsd:jed,nk))
  rit(:,:,:)    = 0.0

  dzmin  = 1.e10  ! meters
  do k=1,nk
    dzmin = min(dzmin,Grid%dzt(k))
  enddo
  if (dzmin >= 25.0) then
    write (stdoutunit,'(/,(1x,a))') &
    '==> Warning: "ppmix" may not work well with coarse vertical resolution'
  endif

  extlat = 0.0
  do j=jsc,jec
    extlat = max(abs(Grid%yt(isc,j)),extlat)
  enddo
  if (extlat > 10.0) then
    write (stdoutunit,'(/,(1x,a))')'==> Warning: "ppmix" may not work well outside tropics where '&
  ,'                             vertical shear is small, unless use solar shortwave penetration. '
  endif

  if (Time_steps%aidif < 1.0) then ! Explicit mixing 
    do k=1,nk
      if ((Time_steps%dtime_t*fricmx)/Grid%dzt(k)**2 >= 0.5) then
        write (stdoutunit,'(/,(1x,a))')&
        '==> Error: vertical diffusive criteria exceeded for "fricmx".  use a smaller'&
       ,'            "dtts" and/or  "fricmx" .... or use "aidif=1.0"'
        write (stdoutunit,'(a48,i3)') ' at level =',k
        call mpp_error(FATAL, &
        '==>Error in ocean_vert_pp_mod: vertical diffusive criteria exceeded for "fricmx" ')
      endif
      if ((Time_steps%dtime_t*diff_cbt_limit)/Grid%dzt(k)**2 >= 0.5) then
        write (stdoutunit,'(/,(1x,a))') &
        '==> Error: vertical diffusive criteria exceeded for "diff_cbt_limit". use a smaller'&
        ,'            "dtts" and/or  "diff_cbt_limit" ...or use "aidif=1.0"'
        write (stdoutunit,'(a48,i3)') ' at level =',k
        call mpp_error(FATAL, &
        '==> Error in ocean_vert_pp_mod: vertical diffusive criteria exceeded for "diff_cbt_limit" ')
      endif
    enddo

    if ((Time_steps%dtime_u*fricmx)/dzmin**2 >= 0.5) then
      write (stdoutunit,'(/,(1x,a))') &
       '==> Error: vertical diffusive criteria exceeded for "fricmx". use a smaller'&
      ,'            "dtuv" and/or "fricmx"  ...or use "aidif=1.0"'
      call mpp_error(FATAL, &
      '==> Error in ocean_vert_pp_mod: vertical diffusive criteria exceeded for "fricmx" ')
    endif

    if ((Time_steps%dtime_u*visc_cbu_limit)/dzmin**2 >= 0.5) then
      write (stdoutunit,'(/,(1x,a))') &
      '==> Error: vertical diffusive criteria exceeded for "visc_cbu_limit". use a smaller'&
      ,'            "dtuv" or "visc_cbu_limit"  ...or use "aidif=1.0"'
      call mpp_error(FATAL, &
      '==> Error in ocean_vert_pp_mod: vertical diffusive criteria exceeded for "visc_cbu_limit" ')
    endif
  endif

  id_diff_cbt_pp = register_diag_field('ocean_model','diff_cbt_pp',Grd%tracer_axes_wt(1:3),&
       Time%model_time, 'vert diffusivity from pp', 'm^2/sec',                             &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbt_pp = register_diag_field('ocean_model','visc_cbt_pp',Grd%tracer_axes_wt(1:3),&
       Time%model_time, 'vert viscosity from pp', 'm^2/sec',                               &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_visc_cbu_pp = register_diag_field('ocean_model','visc_cbu_pp',Grd%vel_axes_wu(1:3),&
       Time%model_time, 'vert viscosity from pp', 'm^2/sec',                            &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

end subroutine ocean_vert_pp_init
! </SUBROUTINE>  NAME="ocean_vert_pp_init"


!#######################################################################
! <SUBROUTINE NAME="vert_mix_pp">
!
! <DESCRIPTION>
! This subroutine computes the vertical diffusivity and viscosity
! according to the Pacanowski and Philander scheme. Mixing coefficients  
! are space and time dependent. 
!
! inputs:
!
!  nk              = number of vertical levels                             
!  grav            = gravity (m/sec^2)                                     
!  fricmx          = max viscosity (m^2/sec)                               
!  wndmix          = min viscosity at bottom of 1st level to simulate      
!                    missing high frequency windstress components (m^2/sec)
!  visc_cbu_back_pp = background "visc_cbu" (m^2/sec)                      
!  diff_cbt_back_pp = background "diff_cbt" (m^2/sec)                      
!  visc_cbu_limit  = largest "visc_cbu" in regions of gravitational        
!                    instability (m^2/sec)                                 
!  diff_cbt_limit  = largest "diff_cbt" in regions of gravitational        
!                    instability (m^2/sec)                                 
!  riu             = richardson number at bottom of U cells                
!  rit             = richardson number at bottom of T cells                
!
! outputs:
!
!  visc_cbu = viscosity at bottom of U cells (m^2/s)  
!  visc_cbt = viscosity at bottom of T cells (m^2/s)  
!  diff_cbt = diffusion at bottom of T cells (m^2/s) 
!
! </DESCRIPTION>
!
  subroutine vert_mix_pp(Time, Thickness, Velocity, T_prog, Dens, visc_cbu, visc_cbt, diff_cbt)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_prog_tracer_type),   intent(in)    :: T_prog(:)
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbu
  real, dimension(isd:,jsd:,:),   intent(inout) :: visc_cbt
  real, dimension(isd:,jsd:,:,:), intent(inout) :: diff_cbt

  integer :: i, j, k
  integer :: taum1
  real    :: t1, t2, t3

  if(.not. use_this_module) then 
    return 
  endif 

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_vert_pp (vert_mix_pp): module must be initialized')
  endif 

  taum1 = Time%taum1

  ! compute gradient richardson number at base of T cells and U cells

  if(horz_grid == MOM_BGRID) then
     call ri_for_pp (Time, Thickness, Velocity, T_prog(index_temp)%field(:,:,:,taum1), &
                     Dens%rho_salinity(:,:,:,taum1), Dens%pressure_at_depth, Dens%rho(:,:,:,taum1))
  else 
      ! average horizontal C-grid velocity onto T-cell center
      wrk1_v(:,:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = Grd%tmask(i,j,k)                                                      &
                    *(Velocity%u(i,j,k,taum1,1)*Grd%dtw(i,j)+Velocity%u(i-1,j,k,taum1,1)*Grd%dte(i,j)) &
                    /(epsln+Grd%dxt(i,j))
               wrk1_v(i,j,k,2) = Grd%tmask(i,j,k)                                                      &
                    *(Velocity%u(i,j,k,taum1,2)*Grd%dts(i,j)+Velocity%u(i,j-1,k,taum1,2)*Grd%dtn(i,j)) &
                    /(epsln+Grd%dyt(i,j))
            enddo
         enddo
      enddo
      call ri_for_cgrid(Time, Thickness%dzwt, Dens%drhodT, Dens%drhodS,          &
           T_prog(index_temp)%field(:,:,:,taum1), Dens%rho_salinity(:,:,:,taum1),&
           wrk1_v(:,:,:,1), wrk1_v(:,:,:,2), rit, riu)
  endif 

  ! viscosity and diffusivity are on bottom of T and U cells.

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        t1                = 1.0/(1.0 + 5.0*riu(i,j,k))
        visc_cbu(i,j,k)   = fricmx*t1**2 + visc_cbu_back_pp
        t2                = 1.0/(1.0 + 5.0*rit(i,j,k))
        diff_cbt(i,j,k,:) = fricmx*t2**3 + diff_cbt_back_pp 
        t3                = 1.0/(1.0 + 5.0*rit(i,j,k))
        visc_cbt(i,j,k)   = fricmx*t3**2 + visc_cbu_back_pp
      enddo
    enddo
  enddo

  ! limit mixing coeffs on bottom of cells in unstable regions

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        if (riu(i,j,k) < 0.0) visc_cbu(i,j,k)   = visc_cbu_limit
        if (rit(i,j,k) < 0.0) visc_cbt(i,j,k)   = visc_cbu_limit
        if (rit(i,j,k) < 0.0) diff_cbt(i,j,k,:) = diff_cbt_limit
      enddo
    enddo
  enddo

  ! approximation for high freq wind mixing near the surface
  ! set no flux through bottom of bottom level "nk"

  k = 1
  where (diff_cbt(isc:iec,jsc:jec,k,:) < wndmix) diff_cbt(isc:iec,jsc:jec,k,:) = wndmix
  diff_cbt(isc:iec,jsc:jec,nk,:) = 0.0
  where (visc_cbu(isc:iec,jsc:jec,k) < wndmix) visc_cbu(isc:iec,jsc:jec,k) = wndmix
  visc_cbu(isc:iec,jsc:jec,nk) = 0.0
  where (visc_cbt(isc:iec,jsc:jec,k) < wndmix) visc_cbt(isc:iec,jsc:jec,k) = wndmix
  visc_cbt(isc:iec,jsc:jec,nk) = 0.0

  diff_cbt(isc:iec,jsc:jec,:,2) = diff_cbt(isc:iec,jsc:jec,:,1)

  call diagnose_3d(Time, Grd, id_diff_cbt_pp, diff_cbt(:,:,:,1))
  call diagnose_3d(Time, Grd, id_visc_cbt_pp, visc_cbt(:,:,:))
  call diagnose_3d_u(Time, Grd, id_visc_cbu_pp, visc_cbu(:,:,:))

end subroutine vert_mix_pp
! </SUBROUTINE> NAME="vert_mix_pp"


!#######################################################################
! <SUBROUTINE NAME="ri_for_pp">
!
! <DESCRIPTION>
! Compute richardson number for the pp scheme
! </DESCRIPTION>
!
subroutine ri_for_pp (Time, Thickness, Velocity, theta, salinity, pressure_at_depth, rho_taum1)

  type(ocean_time_type),         intent(in) :: Time
  type(ocean_thickness_type),    intent(in) :: Thickness
  type(ocean_velocity_type),     intent(in) :: Velocity
  real, dimension(isd:,jsd:,:),  intent(in) :: theta
  real, dimension(isd:,jsd:,:),  intent(in) :: salinity
  real, dimension(isd:,jsd:,:),  intent(in) :: pressure_at_depth
  real, dimension(isd:,jsd:,:),  intent(in) :: rho_taum1

  integer :: i, j, k, taum1
  real    :: active_cells, fx, t1

  fx   = -0.25*grav*rho0r

  taum1 = Time%taum1

  ! compute density difference across bottom of T cells at tau-1

  wrk1(:,:,:) = density_delta_z(rho_taum1(:,:,:), salinity(:,:,:), &
                                theta(:,:,:), pressure_at_depth(:,:,:))

  ! compute richardson numbers on bottom of U cells

  do k=1,nk-1
    do j=jsd,jed-1
      do i=isd,ied-1
        t1 = fx*Thickness%dzwt(i,j,k)
        riu(i,j,k) = t1*Grd%umask(i,j,k+1)*(wrk1(i,j+1,k) + wrk1(i+1,j+1,k) + wrk1(i,j,k)   + wrk1(i+1,j,k)) /&
                      ((Velocity%u(i,j,k,1,taum1) - Velocity%u(i,j,k+1,1,taum1))**2 + &
                       (Velocity%u(i,j,k,2,taum1) - Velocity%u(i,j,k+1,2,taum1))**2 + epsln)
      enddo
    enddo
  enddo

  ! compute richardson numbers on bottom of T cells as average
  ! of four nearest richardson numbers on bottom of U cells.
  ! (do not consider land cells in the average... only active ones)

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        active_cells = Grd%umask(i,j,k+1)   + Grd%umask(i-1,j,k+1) &
                     + Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln
        rit(i,j,k)   = (riu(i,j,k) + riu(i-1,j,k) + riu(i,j-1,k) + riu(i-1,j-1,k))/active_cells

        ! make sure no static instability exists (one that is not seen
        ! by the Richardson number).  This may happen due to
        ! horizontal averaging used in calculating the Richardson
        ! number.

        if (rit(i,j,k) > 0.0 .and. wrk1(i,j,k) > 0.0) then
          rit(i,j,k) = -10.0
        endif
      enddo
    enddo
  enddo

end subroutine ri_for_pp
! </SUBROUTINE> NAME="ri_for_pp"


end module ocean_vert_pp_mod

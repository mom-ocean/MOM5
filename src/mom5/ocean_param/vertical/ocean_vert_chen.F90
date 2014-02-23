module ocean_vert_chen_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Vertical viscosity and diffusivity according Chen scheme 
!</OVERVIEW>
!
!<DESCRIPTION>
! This scheme was originally developed by researchers 
! at the CSIRO Marine and Atmospheric Research and 
! Bureau of Meteorology, both in Australia.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Chen, D., L.M. Rothstein and A.J. Busalacchi, 1994:
! A hybrid vertical mixing scheme and its application 
! to tropical ocean models, 
! J. Phys. Oceanogr., 24, 2156-2179
! </REFERENCE>
!
! <NOTE>
! Surface fresh water contributes to surface buoyancy via conversion to 
! a locally implied salt flux. 
! </NOTE>
!
! <NOTE>
! This module typically runs with the ocean_shortwave_csiro scheme.
! </NOTE>
!
! <NOTE>
! Use smf_bgrid since this array contains the primary smf array read in from 
! from the coupler in ocean_core/ocean_sbc.F90, when using the FMS coupler.
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_chen_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false. 
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  Default debug_this_module=.false. 
!  </DATA> 
!
!  <DATA NAME="diff_cbt_iw" UNITS="m^2/sec" TYPE="real">
!  Background vertical diffusivity.  Note that if using Bryan-Lewis as a 
!  background diffusivity, then should set diff_cbt_iw=0.0. 
!  </DATA> 
!  <DATA NAME="visc_cbu_iw" UNITS="m^2/sec" TYPE="real">
!  Background vertical viscosity
!  </DATA> 
!  <DATA NAME="visc_cbu_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical viscosity due to shear instability 
!  </DATA> 
!  <DATA NAME="diff_cbt_limit" UNITS="m^2/sec" TYPE="real">
!  Enhanced vertical diffusivity due to shear instability 
!  </DATA> 
!  <DATA NAME="bulk_tn" UNITS="" TYPE="real">
!  Bulk turblence parameter n_0
!  </DATA> 
!  <DATA NAME="bulk_tm" UNITS="" TYPE="real">
!  Bulk turblence parameter m_0
!  </DATA> 
!  <DATA NAME="hbl_growth_max" UNITS="m/sec" TYPE="real">
!  Maximum growth rate of kraus mixed layer 
!  </DATA> 
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: FATAL, NOTE, WARNING, stdout, stdlog, file_exist
use fms_mod,          only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type
use mpp_domains_mod,  only: mpp_update_domains, NUPDATE, EUPDATE
use mpp_mod,          only: input_nml_file, mpp_error

use ocean_density_mod,         only: density, density_delta_z, density_delta_sfc
use ocean_domains_mod,         only: get_local_indices
use ocean_parameters_mod,      only: missing_value, rho0r, rho_cp, cp_ocean, grav
use ocean_parameters_mod,      only: MOM_BGRID, MOM_CGRID
use ocean_shortwave_csiro_mod, only: F_vis, ssw_atten_depth
use ocean_types_mod,           only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,           only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,           only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_util_mod,            only: write_timestamp, diagnose_2d, write_chksum_2d
use ocean_vert_util_mod,       only: ri_for_cgrid
use ocean_workspace_mod,       only: wrk1, wrk1_v

implicit none

private

public ocean_vert_chen_init
public ocean_vert_chen_end
public vert_mix_chen
public ocean_vert_chen_restart

private kraus_turner
private ri_for_chen
private nr_for_chen

#ifdef MOM_STATIC_ARRAYS

integer, parameter :: isc = 1, iec = NI_LOCAL_ , jsc = 1, jec = NJ_LOCAL_
integer, parameter :: isd = 0, ied = NI_LOCAL_ + 1 , jsd = 0, jed = NJ_LOCAL_ + 1
integer, parameter :: isg = 1, ieg = NI_ , jsg = 1, jeg = NJ_
integer, parameter :: nk = NK_


real, dimension(isd:ied,jsd:jed)                :: ustar      ! surface friction velocity       (m/s)
real                :: BoT        ! thermal turb buoy. forcing  (m^2/s^3)
real, dimension(isd:ied,jsd:jed)                :: Bosol      ! radiative buoy forcing      (m^2/s^3)
real                :: BoS        ! salinity buoy forcing      (m^2/s^3)
real                :: BoTot      ! total buoy forcing      (m^2/s^3)
real, dimension(isd:ied,jsd:jed)             :: dbloc      ! local delta buoy at interfaces(m/s^2)

integer, dimension(isd:ied,jsd:jed)             :: kbl        ! index of first grid level below hbl

real, private, dimension(isd:ied,jsd:jed,nk)   :: riu      ! Richardson number at base of U-cells
real, private, dimension(isd:ied,jsd:jed,nk)   :: rit      ! Richardson number at base of T-cells
real, private, dimension(isd:ied,jsd:jed)      :: hbl      ! boundary layer depth

#else

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk


real, dimension(:,:), allocatable        :: rhosfc     ! potential density of sfc layer(kg/m^3)
real, dimension(:,:), allocatable        :: ustar      ! surface friction velocity       (m/s)
real        :: BoT        ! 
real, dimension(:,:), allocatable        :: Bosol      ! 
real        :: BoS        !
real        :: BoTot      !
real, dimension(:,:), allocatable      :: dbloc      ! local delta buoy at interfaces(m/s^2)
integer, dimension(:,:), allocatable     :: kbl        ! index of first grid level below hbl

real, private, dimension(:,:,:), allocatable   :: riu      ! Richardson number at base of U-cells
real, private, dimension(:,:,:), allocatable   :: rit      ! Richardson number at base of T-cells
real, private, dimension(:,:),   allocatable   :: hbl      ! boundary layer depth

#endif

real, parameter :: epsilon = 0.1

! for diagnostics 
integer :: id_diff_cbt_chen_t, id_diff_cbt_chen_s,id_hbl
integer :: id_wmix,id_bte,id_dbloc,id_tke
logical :: used

! for restart
type(restart_file_type), save :: Che_restart

real :: tracer_timestep=0.0
real :: aidif 

integer :: num_prog_tracers=0, index_temp, index_salt


! Namelist 

real :: visc_cbu_limit     = 50.0e-4 ! max visc due to shear instability
real :: diff_cbt_limit     = 50.0e-4 ! max diff due to shear instability
real :: visc_con_limit     = 0.1     ! m^2/s. visc due to convective instability
real :: diff_con_limit     = 0.1     ! m^2/s. diff due to convective instability
real :: visc_cbu_iw        = 1.0e-4  ! m^2/s. visc background due to internal waves
real :: diff_cbt_iw        = 0.1e-4  ! m^2/s. diffusivity background due to internal waves
 

! bulk turbulence parameters
real :: bulk_tm        = 0.4
real :: bulk_tn        = 0.18
real :: hbl_growth_max = 100.    ! (m/s) set large by default

! for Bgrid or Cgrid
integer :: horz_grid

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer :: Grd =>NULL()

character(len=256) :: version=&
     '$Id: ocean_vert_chen.F90,v 20.0 2013/12/14 00:16:36 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: module_is_initialized = .FALSE.
logical :: use_this_module       = .false.
logical :: debug_this_module     = .false. 

namelist /ocean_vert_chen_nml/ use_this_module, debug_this_module, &
                               diff_cbt_iw, visc_cbu_iw,           &
                               visc_cbu_limit, diff_cbt_limit,     &
                               bulk_tn, bulk_tm, hbl_growth_max,   &
                               visc_con_limit, diff_con_limit

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_chen_init">
!
! <DESCRIPTION>
! Initialization for the Chen vertical mixing scheme
!
!     input:
!       dzt    = thickness of vertical levels (m)                        
!       km     = number of vertical levels                               
!       yt     = latitude of grid points (deg)                           
!       dtts   = density time step (sec)                                 
!       dtuv   = internal mode time step (sec)                           
!       error  = logical to signal problems                              
!       vmixset= logical to determine if a vertical mixing scheme was    
!                chosen
!
!     output:
!       visc_cbu_limit = visc max due to shear instability  (m**2/sec)   
!       diff_cbt_limit = diffusivity ..                     (m**2/sec)   
!       visc_cbu_iw  = visc background due to internal waves(m**2/sec)   
!       diff_cbt_iw  = diffusivity ..                       (m**2/sec)   
!       error  = true if some inconsistancy was found                    
!
! </DESCRIPTION>
!
subroutine ocean_vert_chen_init (Grid, Domain, Time, Time_steps, T_prog, hor_grid, debug)
  
  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
  integer,                      intent(in)           :: hor_grid 
  logical,                      intent(in), optional :: debug 

  integer :: i, j, k, ioun, io_status, ierr, n
  integer :: id_restart
  real    :: hbl_init=25.0

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
       call mpp_error(FATAL,  '==>Error: ocean_vert_chen_init attempted reinitialization')
  endif 
  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_vert_chen_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_vert_chen_nml')
#else
  ioun =  open_namelist_file ()
  read  (ioun, ocean_vert_chen_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_chen_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_chen_nml)  
  write (stdlogunit, ocean_vert_chen_nml)

  if(use_this_module) then
      write(stdoutunit,'(/1x,a)')'==> NOTE: USING Chen vertical mixing scheme.'
      write(stdoutunit,'(1x,a/)')'==> NOTE: Chen is to be run with shortwave flux from ocean_shortwave_csiro_mod'
  else
      write(stdoutunit,'(/1x,a)')'==> NOTE: NOT USING Chen vertical mixing scheme.'
      return
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_vert_chen_mod with debug_this_module=.true.'  
  endif 

  tracer_timestep = Time_steps%dtts
  aidif           = Time_steps%aidif 
  horz_grid       = hor_grid
  
  index_temp=-1;index_salt=-1
  num_prog_tracers = size(T_prog)
  do n= 1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,&
    '==>Error: ocean_vert_chen_mod: temp and/or salt not present in tracer array')
  endif 
  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  allocate (riu(isd:ied,jsd:jed,nk))
  allocate (rit(isd:ied,jsd:jed,nk))
  allocate (hbl(isd:ied,jsd:jed))

  allocate (ustar(isd:ied,jsd:jed))     ! surface friction velocity       (m/s)
  allocate (dbloc(isd:ied,jsd:jed))     ! local delta buoy at interfaces(m/s^2)
  allocate (kbl(isd:ied,jsd:jed))       ! index of first grid level below hbl
  allocate(Bosol(isd:ied,jsd:jed))
#endif

  riu(:,:,:)    = 0.0
  rit(:,:,:)    = 0.0
  id_restart = register_restart_field(Che_restart, 'kraus.res.nc', 'kbl',kbl, &
               domain=Dom%domain2d)
  id_restart = register_restart_field(Che_restart, 'kraus.res.nc', 'hbl',hbl, &
               domain=Dom%domain2d)
  if(file_exist('INPUT/kraus.res.nc')) then
      call restore_state(Che_restart)
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_vert_chen: reading hbl from kraus.res.nc'
      call write_timestamp(Time%model_time)
      call write_chksum_2d('hbl', hbl(COMP)*Grd%tmask(COMP,1))
  else
      hbl(:,:) = hbl_init*Grd%tmask(:,:,1)
      kbl(:,:)=1
      do k=nk,1,-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%zt(k) > hbl(i,j)) then
                  kbl(i,j)=k
               endif
            enddo
         enddo
      enddo
  endif

  call mpp_update_domains(kbl, Dom%domain2d)
  call mpp_update_domains(hbl, Dom%domain2d)

!-----------------------------------------------------------------------
! error checks and diagnostics setup 
!-----------------------------------------------------------------------

  id_wmix = register_diag_field('ocean_model','wmix',Grd%tracer_axes_wt(1:2),&
       Time%model_time, 'wind mixing', 'm^3/s^3',                            &
       missing_value = missing_value, range=(/-1.e5,1.e5/))
  id_bte = register_diag_field('ocean_model','bte',Grd%tracer_axes_wt(1:2),&
       Time%model_time, 'bulk_Tke term', 'm^2/s^3',                        &
       missing_value = missing_value, range=(/-1.e5,1.e5/))
  id_dbloc = register_diag_field('ocean_model','dbloc',Grd%tracer_axes_wt(1:2),&
       Time%model_time, 'delta density ', 'm/s^2',                             &
       missing_value = missing_value, range=(/-1.e5,1.e5/))
  id_tke = register_diag_field('ocean_model','tke',Grd%tracer_axes_wt(1:2),&
       Time%model_time, 'Rhs of tke', 'm^3/s^3',                           &
       missing_value = missing_value, range=(/-1.e5,1.e5/))
  id_hbl = register_diag_field('ocean_model','hbl',Grd%tracer_axes_wt(1:2),&
       Time%model_time, 'Kraus mixed layer depth', 'm',                    &
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_chen_t = register_diag_field('ocean_model','diff_cbt_chen_t',            &
       Grd%tracer_axes_wt(1:3), Time%model_time, 'vert diffusivity from chen for temp',&
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_chen_s = register_diag_field('ocean_model','diff_cbt_chen_s',                 &
       Grd%tracer_axes_wt(1:3), Time%model_time, 'vert diffusivity from chen for salinity', &
       'm^2/sec', missing_value = missing_value, range=(/-1.e5,1.e5/))

end subroutine ocean_vert_chen_init
! </SUBROUTINE> NAME="ocean_vert_mix_coeff_init">


!#######################################################################
! <SUBROUTINE NAME="vert_mix_chen">
!
! <DESCRIPTION>
!
! --Compute interior mixing everywhere:                               
!   interior mixing gets computed at all cell interfaces due to constant
!   internal wave background activity ("visc_cbu_iw" and "diff_cbt_iw").
!   Additionally, mixing can be enhanced by contribution from shear 
!   instability which is a function of the local Ri.
!
! --Boundary layer:
!
!   (A) Boundary layer depth:
!       at every gridpoint the depth of the Kraus boundary layer 
!       ("hbl") gets computed.
!
!   (B) Boundary layer mixing:
!       within the boundary layer, above hbl, vertical mixing is 
!       set to a maximum
!
! inputs
!
! outputs
!
!  hbl      = boundary layer depth (meters)
!  visc_cbu = viscosity coefficient at bottom of U cells (m^2/s) <BR/> 
!  diff_cbt = diffusion coefficient at bottom of T cells (m^2/s) <BR/>
! </DESCRIPTION>
!
  subroutine vert_mix_chen(Time, Thickness, Velocity, T_prog, Dens, &
                           swflx,  pme, visc_cbu, visc_cbt, diff_cbt)

  type(ocean_time_type),                 intent(in)    :: Time
  type(ocean_thickness_type),            intent(in)    :: Thickness
  type(ocean_velocity_type),             intent(in)    :: Velocity
  type(ocean_prog_tracer_type),          intent(in)    :: T_prog(:)
  type(ocean_density_type),              intent(in)    :: Dens
  real, dimension(isd:,jsd:),            intent(in)    :: swflx
  real, dimension(isd:,jsd:),            intent(in)    :: pme
  real, dimension(isd:ied,jsd:jed,nk),   intent(inout) :: visc_cbu
  real, dimension(isd:ied,jsd:jed,nk),   intent(inout) :: visc_cbt
  real, dimension(isd:ied,jsd:jed,nk,2), intent(inout) :: diff_cbt

!
! Coefficients for richardson no. mixing.
!
  real, parameter  :: a_pt_8=6.7490227E-8
  real, parameter  :: b_pt_2=8.6118424E-4
  real, parameter  :: d_pt_10=2.0269152E-9
  real, parameter  :: e_pt_3=8.4795331E-4    
  real :: t1, t2, t3, rriu_8, rrit_10, rrit_8

  real, dimension(isd:ied,jsd:jed,nk)  :: mixmask
  integer :: i, j, k
  integer :: tau, taum1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) call &
       mpp_error(FATAL, '==>Error: Module ocean_vert_mix_coeff_mod must be initialized')

  tau   = Time%tau
  taum1 = Time%taum1

  visc_cbu(:,:,:) = 0.0
  visc_cbt(:,:,:) = 0.0

  ! compute mixed layer according to Kraus-Turner algorithm
  call kraus_turner(Time, Velocity, T_prog, Dens, swflx, pme, mixmask)

!-----------------------------------------------------------------------
!     compute gradient Ri 
!-----------------------------------------------------------------------

  if(horz_grid == MOM_BGRID) then 
     call ri_for_chen(Time, Thickness, Velocity, T_prog, Dens%pressure_at_depth, &
                      Dens%rho(:,:,:,tau), Dens%rho(:,:,:,taum1), Dens%rho_salinity)
  else
      ! average horizontal C-grid velocity onto T-cell center
      wrk1_v(:,:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = Grd%tmask(i,j,k)                                                  &
                    *(Velocity%u(i,j,k,tau,1)*Grd%dtw(i,j)+Velocity%u(i-1,j,k,tau,1)*Grd%dte(i,j)) &
                    /(epsln+Grd%dxt(i,j))
               wrk1_v(i,j,k,2) = Grd%tmask(i,j,k)                                                  &
                    *(Velocity%u(i,j,k,tau,2)*Grd%dts(i,j)+Velocity%u(i,j-1,k,tau,2)*Grd%dtn(i,j)) &
                    /(epsln+Grd%dyt(i,j))
            enddo
         enddo
      enddo
      call ri_for_cgrid(Time, Thickness%dzwt, Dens%drhodT, Dens%drhodS,          &
           T_prog(index_temp)%field(:,:,:,taum1), Dens%rho_salinity(:,:,:,taum1),&
           wrk1_v(:,:,:,1), wrk1_v(:,:,:,2), rit, riu)
  endif 

! Do Richardson number mixing according to whatever law you like.
! Perhaps should make this a distinct user supplied/customisable routine.
!
! Law used here is
! visc= a/Ri^8  + b/(1+5Ri)^2
! diff= d/Ri^10 + e/(1+Ri)^3
!
! compute C-grid viscosity on T-points and 
! assume east value same as north value.
!
  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        t1                   = 1.0/(1.0 + 5.0*riu(i,j,k))
        t3                   = 1.0/(1.0 + 5.0*rit(i,j,k))
        rriu_8               = 1.0/(riu(i,j,k)**8+epsln)
        visc_cbu(i,j,k)      = mixmask(i,j,k)*visc_cbu_limit + &
                               (a_pt_8*rriu_8+b_pt_2*t1**2)+visc_cbu_iw
        rrit_8               = 1.0/(rit(i,j,k)**8+epsln)
        visc_cbt(i,j,k)      = mixmask(i,j,k)*visc_cbu_limit + &
                              (a_pt_8*rrit_8+b_pt_2*t3**2)+visc_cbu_iw
        t2                   = 1.0/(1.0 + 5.0*rit(i,j,k))
        rrit_10              = 1.0/(rit(i,j,k)**10+epsln)
        diff_cbt(i,j,k,1)    = mixmask(i,j,k)*diff_cbt_limit + &
                               (d_pt_10*rrit_10+e_pt_3*t2**3)+diff_cbt_iw
      enddo
    enddo
  enddo

! Limit
  
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           diff_cbt(i,j,k,1)   = min(diff_cbt(i,j,k,1),  diff_cbt_limit)*Grd%tmask(i,j,k+1)
           visc_cbt(i,j,k)     = min(visc_cbt(i,j,k),    visc_cbu_limit)*Grd%tmask(i,j,k+1)
           visc_cbu(i,j,k)     = min(visc_cbu(i,j,k),    visc_cbu_limit)*Grd%umask(i,j,k+1)
        enddo
     enddo
  enddo

  ! Convection
  if(aidif == 1.0) then
     do k=1,nk-1
        do j=jsc,jec
           do i=isc,iec
              if(rit(i,j,k) < 0.0) then
                 diff_cbt(i,j,k,1) = diff_con_limit*Grd%tmask(i,j,k+1)
                 visc_cbu(i,j,k)   = visc_con_limit*Grd%umask(i,j,k+1)  
                 visc_cbt(i,j,k)   = visc_con_limit*Grd%tmask(i,j,k+1)  
              endif
           enddo
        enddo
     enddo
  endif


  ! default salinity coefficient equal to thermal (assume no double diffusion) 
  diff_cbt(:,:,:,2) = diff_cbt(:,:,:,1)

  if(debug_this_module) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_vert_chen: intermediate chksum'
      call write_timestamp(Time%model_time)
      call write_chksum_2d('hbl', hbl(COMP)*Grd%tmask(COMP,1))
      call write_chksum_2d('diff_cbt', diff_cbt(COMP,:,1)*Grd%tmask(COMP,:))
  endif
  
end subroutine vert_mix_chen
! </SUBROUTINE> NAME="vert_mix_chen"


!#######################################################################
! <SUBROUTINE NAME="kraus_turner">
!
! <DESCRIPTION>
!
! Calculate the Kraus mixed layer depth
!
! Note: This formulation assumes a single exponential decay in the solar
! shortwave penetration.
!
! Use smf_bgrid since this array contains the primary smf array read in from 
! from the coupler in ocean_core/ocean_sbc.F90, when using the FMS coupler.
!
!      real dbloc(ij_bounds)  = local delta buoyancy         (m/s^2)        
!      real ustar(ij_bounds)  = surface friction velocity     (m/s)      
!      real Bo(ij_bounds)     = surface turbulent buoyancy forcing(m^2/s^3)
!      real Bosol(ij_bounds)  = radiative buoyancy forcing (m^2/s^3)       
!
!  output
!      real hbl(ij_bounds)        ! boundary layer depth              (m)      
!      real mixmask(ij_bounds,nk) ! fraction of cell which resides in mixed layer
!      integer kbl(ij_bounds)     ! index of first grid level below hbl        
!
! </DESCRIPTION>
subroutine kraus_turner(Time, Velocity, T_prog, Dens, swflx, pme, mixmask)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_velocity_type),    intent(in)  :: Velocity
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)
  type(ocean_density_type),     intent(in)  :: Dens
  real, dimension(isd:,jsd:),   intent(in)  :: swflx
  real, dimension(isd:,jsd:),   intent(in)  :: pme
  real, dimension(isd:,jsd:,:), intent(out) :: mixmask

  real, dimension(isd:ied,jsd:jed) :: talpha !-d(rho)/ d(pot.temperature)  (kg/m^3/C)
  real, dimension(isd:ied,jsd:jed) :: sbeta  !d(rho)/ d(salinity)          (kg/m^3/PSU)
  real, dimension(isd:ied,jsd:jed) :: tempbl,saltbl,pressbl  ! Values across boundary layer 
  real, dimension(isd:ied,jsd:jed) :: density_mixed, density_mixed_1
  real, dimension(isd:ied,jsd:jed) :: wind_mixing,tke,bulk_tke

  integer, dimension(isd:ied,jsd:jed) :: kblm1  ! level immediately above  mixed layer
  integer, dimension(isd:ied,jsd:jed) :: kblp1  ! next level immediately below  mixed layer

  integer :: i,j,k, kbl_m1
  integer :: tau, taum1
  real    :: shallow
  real    :: smftu, smftv, active_cells
  real    :: pen_ke,exdepth,hbltemp  

  tau   = Time%tau
  taum1 = Time%taum1   
!
! Layers above and below mixed layer interface
!
  kblm1=max(kbl-1,1)
  kblp1=max(min(kbl+1,Grd%kmt),1)

! result from density_derivs is d(rho)/d(theta), yet need minus for this routine. 
  do j=jsd,jed
     do i=isd,ied
        talpha(i,j) = -Dens%drhodT(i,j,1)
        sbeta(i,j)  =  Dens%drhodS(i,j,1)
     enddo
  enddo

  ! use smf_bgrid since this is directly from the FMS coupler. 
  do j=jsc,jec
     do i=isc,iec
        active_cells = Grd%umask(i,j,1) + Grd%umask(i-1,j,1) +  Grd%umask(i,j-1,1) + Grd%umask(i-1,j-1,1) + epsln
        smftu = (Velocity%smf_bgrid(i,j,1)   + Velocity%smf_bgrid(i-1,j,1) &
               + Velocity%smf_bgrid(i,j-1,1) + Velocity%smf_bgrid(i-1,j-1,1) )/active_cells
        smftv = (Velocity%smf_bgrid(i,j,2)   + Velocity%smf_bgrid(i-1,j,2) &
               + Velocity%smf_bgrid(i,j-1,2) + Velocity%smf_bgrid(i-1,j-1,2) )/active_cells
        ustar(i,j) = sqrt(rho0r*sqrt( smftu**2 + smftv**2 ))
     enddo
  enddo


! Surface Buoyancy fluxes
! Note stf is now ( heat_flux/cp_ocean)
! swflux is proper flux.
! Heat flux term : -(alpha*g/rho_cp)*(Heat flux into ocean)
! SSW term       : -(alpha*g/rho_cp)*(penetrating component of SSW)
! Salt term      : (beta*g)*salinty(surface)*(freshwater flux into ocean)
!

  do j=jsc,jec
    do i=isc,iec
!
! Calculate buoyancy fluxes in mixed layer
!
      BoT=-(grav*rho0r)*(talpha(i,j)/ (Dens%rho(i,j,1,tau)+epsln)) * T_prog(index_temp)%stf(i,j)
      Bosol(i,j) =-(grav*rho0r)* (talpha(i,j)/(Dens%rho(i,j,1,tau)+epsln)) * F_vis * swflx(i,j)/cp_ocean
      BoS=-(grav*rho0r) * Dens%rho_salinity(i,j,1,tau) *pme(i,j) * sbeta(i,j) / (Dens%rho(i,j,1,tau)+epsln)
!
      BoTot = BoT + BoS - Bosol(i,j)

! Penetrative SSW contribution
      exdepth=exp(-hbl(i,j)/ssw_atten_depth(i,j))
      pen_ke= Bosol(i,j)*(hbl(i,j)*(1+exdepth)-2.0*ssw_atten_depth(i,j)*(1.0-exdepth))

      wind_mixing(i,j)= 2.0 * bulk_tm * ustar(i,j)**3
      bulk_tke(i,j)=((1+bulk_tn)*BoTot-(1.0-bulk_tn)*abs(BoTot))/2.0
      tke(i,j)=Grd%tmask(i,j,1)*(wind_mixing(i,j)+hbl(i,j)*bulk_tke(i,j)+pen_ke)
    enddo
  enddo
  do j=jsd,jed
    do i=isd,ied
      if(kbl(i,j) > 0) then
         saltbl(i,j) = Dens%rho_salinity(i,j,kbl(i,j),taum1)
         tempbl(i,j) = T_prog(index_temp)%field(i,j,kbl(i,j),taum1)
         pressbl(i,j)= Dens%pressure_at_depth(i,j,kbl(i,j))
      end if
    enddo
  enddo
!
! Calculate densities across mixed layer
!
  density_mixed=density(saltbl,tempbl,pressbl)

  do j=jsd,jed
    do i=isd,ied
      saltbl(i,j) = Dens%rho_salinity(i,j,kblm1(i,j),taum1)
      tempbl(i,j) = T_prog(index_temp)%field(i,j,kblm1(i,j),taum1)
    enddo
  enddo
  density_mixed_1=density(saltbl,tempbl,pressbl)
!
! Calculate and limit buoyancy differential
!
  dbloc=max(2.5e-4,grav*(density_mixed-density_mixed_1)*rho0r)


!
! Entrain mixed layer where necessary and limit by depth of next layer
! and maximum growth rate
!
  do j=jsc,jec
    do i=isc,iec
      hbltemp=hbl(i,j)+tracer_timestep*min(max(tke(i,j),0.0)/(dbloc(i,j)*hbl(i,j)+epsln),hbl_growth_max)
!
! Alternate update which will be better when winds dominate
!
! hbltemp=min(sqrt(hbl(i,j)**2+2.0*tracer_timestep*(max(tke(i,j),0.0)/dbloc(i,j)), &
!              hbl(1,j)+tracer_timestep*hbl_growth_max)
      hbl(i,j)=min(hbltemp,Grd%zt(kblp1(i,j)))*Grd%tmask(i,j,1)
    enddo
  enddo


  kbl(:,:)=1
  do j=jsc,jec
!
! Newton Raphson scheme for detrainment.
!
! Convergence is checked on a row by row basis in order to vectorise.
!
    call nr_for_chen(wind_mixing(isc:iec,j),bulk_tke(isc:iec,j),Bosol(isc:iec,j),ssw_atten_depth(isc:iec,j), &
                tke(isc:iec,j),hbl(isc:iec,j))
    mixmask(:,j,:)=0.0
    do k=1,nk
      do i=isc,iec
        if(hbl(i,j) >= Grd%zt(k)) then 
          kbl(i,j)=min(Grd%kmt(i,j),kbl(i,j)+1)
          mixmask(i,j,k)=1.0
        endif
      enddo
    enddo
    kbl=max(kbl,1)
    do i=isc,iec
      kbl_m1=max(1,kbl(i,j)-1)
      shallow=max(hbl(i,j),Grd%zt(1))-Grd%zt(kbl_m1)
      mixmask(i,j,kbl(i,j))=shallow/(Grd%zt(kbl(i,j))-Grd%zt(kbl_m1)+ epsln)
    enddo
  enddo

  call diagnose_2d(Time, Grd, id_wmix, wind_mixing(:,:))
  call diagnose_2d(Time, Grd, id_bte, bulk_tke(:,:))
  call diagnose_2d(Time, Grd, id_dbloc, dbloc(:,:))
  call diagnose_2d(Time, Grd, id_tke, tke(:,:))
  call diagnose_2d(Time, Grd, id_hbl, hbl(:,:))

end subroutine kraus_turner
! </SUBROUTINE> NAME="kraus_turner"


subroutine nr_for_chen(a,b,c,h,tkeinit,z)
!
! Hybrid Newton Raphson, half interval search
!
! Equation is of the form
!
! f(z)=a+b*z+c*(z*(1+exp(z/h)-2*h*(1-exp(z/h))
!
! Where a,b,c,h are constant for each point.
!
! We wish to solve f(z)=0 for those points where detrainment takes place. That
! is when f(z_init)=tkinit < 0
!
! Note that a>=0 so we must have a solution in the interval [0,z_init)
!
! A global convergence criteria of 0.01m is used
! No changes are made if the local change is less than crit. This should enable reproduciblity.
! 

  real, dimension(isc:iec),intent(in)    :: a,b,c,h,tkeinit
  real, dimension(isc:iec),intent(inout) :: z


  real, dimension(isc:iec)   :: z_left,z_half
  real, dimension(isc:iec)   :: dz_old,dz
  real,dimension(isc:iec)    :: f,df
  real                       :: texp
  logical:: converge
  integer i,its
  real :: crit
  crit = 0.01

  do i=isc,iec
    z_half(i)=0.0
    z_left(i)=z(i)
    if(tkeinit(i) < 0) then
      dz_old(i)=z_left(i)
      dz(i)=z_left(i)
    else
      dz_old(i)=0.0
      dz(i)=0.0
    endif
    z(i)=z_left(i)
    f(i)=tkeinit(i)
    df(i)=b(i)+c(i)*(1.0-exp(-z(i)/h(i)))*(1.0+z(i)/h(i))
  enddo
  converge=.false.
  its=0
  do while(its .lt. 20 .and. (.not. converge))
    its=its+1
    do i=isc,iec
      if(tkeinit(i) < 0.0 .and. abs(dz(i)) > crit) then
        if(((z(i)-z_half(i))*df(i)-f(i))*((z(i)-z_left(i))*df(i)-f(i)) >= 0.0 .or. &
           abs(2.0*f(i)) > abs(dz_old(i)*df(i))) then

! If proposed solution is outside range or slow convergence use half interval
! search
           
          dz_old(i)=dz(i)
          dz(i)=0.5*(z_half(i)-z_left(i))
          z(i)=z_left(i)+dz(i)
        else
          dz_old(i)=dz(i)
          dz(i)=f(i)/df(i)
          z(i)=z(i)-dz(i)
        endif
        texp=exp(-z(i)/h(i))
        f(i)=a(i)+z(i)*b(i)+c(i)*(z(i)*(1+texp)-2.0*h(i)*(1.0-texp))
        df(i)=b(i)+c(i)*(1-texp*(1.0+z(i)/h(i)))
        if(f(i) < 0.0) then
          z_left(i)=z(i)
        else
          z_half(i)=z(i)
        endif
      endif
    enddo
    converge=maxval(abs(dz)) <= crit
  enddo
  end subroutine nr_for_chen


!#######################################################################
! <SUBROUTINE NAME="ri_for_chen">
!
! <DESCRIPTION>
! Compute Richardson number on tracer and velocity cell bottoms. 
!  rit = richardson number at bottom of T cells <BR/>
!  riu = richardson number at bottom of U cells
! </DESCRIPTION>
!
subroutine ri_for_chen (Time, Thickness, Velocity, T_prog, pressure_at_depth, rho, rho_taum1, rho_salinity)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_thickness_type),     intent(in) :: Thickness
  type(ocean_velocity_type),      intent(in) :: Velocity
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:)
  real, dimension(isd:,jsd:,:),   intent(in) :: pressure_at_depth
  real, dimension(isd:,jsd:,:),   intent(in) :: rho
  real, dimension(isd:,jsd:,:),   intent(in) :: rho_taum1
  real, dimension(isd:,jsd:,:,:), intent(in) :: rho_salinity

  integer :: i, j, k, tlev,  mr, tau, taum1
  real :: active_cells, fx, t1 

  real  :: ri_temp1(isd:ied),ri_temp2(isd:ied)

  logical :: smooth_richardson_number = .true.
  integer :: num_smoothings = 1 ! for vertical smoothing of Richardson number

  real,parameter :: epsln_chen=1.e-10


  fx   = -0.25*grav*rho0r

  tau = Time%tau
  taum1 = Time%taum1

  ! compute density difference across bottom of T cells at tau-1

  if (aidif == 1.0) then
    tlev = tau
    wrk1(:,:,:) = density_delta_z(rho(:,:,:), rho_salinity(:,:,:,tau),             &
                                              T_prog(index_temp)%field(:,:,:,tau), & 
                                              pressure_at_depth(:,:,:))
  else
    tlev = taum1
    wrk1(:,:,:) = density_delta_z(rho_taum1(:,:,:), rho_salinity(:,:,:,taum1), &
                                 T_prog(index_temp)%field(:,:,:,taum1), pressure_at_depth(:,:,:))
  endif

  ! compute richardson numbers on bottom of U cells
  do k=1,nk-1
    do j=jsd,jed-1
      do i=isd,ied-1
        t1 = fx*Thickness%dzwt(i,j,k)
        riu(i,j,k) = t1*Grd%umask(i,j,k+1)*(wrk1(i,j+1,k) + wrk1(i+1,j+1,k) + wrk1(i,j,k)   + wrk1(i+1,j,k)) /&
                      ((Velocity%u(i,j,k,1,tlev) - Velocity%u(i,j,k+1,1,tlev))**2 + &
                       (Velocity%u(i,j,k,2,tlev) - Velocity%u(i,j,k+1,2,tlev))**2 + epsln_chen)
      enddo
    enddo
  enddo

  if (smooth_richardson_number) then

    ! smooth Richardson number in the vertical using a 1-2-1 filter:

    do mr = 1,num_smoothings
      do j=jsd,jed-1
         do i=isd,ied-1
            ri_temp2(i)    =  riu(i,j,1)
         enddo
         do k=2,nk-2
!CDIR NODEP
            do i=isd,ied-1
              if(k < Grd%kmu(i,j) -1) then
              ri_temp1(i)        =  ri_temp2(i)
              ri_temp2(i)        =  riu(i,j,k)
              riu(i,j,k) = 0.25*ri_temp1(i) + 0.5 * riu(i,j,k) + 0.25 * riu(i,j,k+1)
              endif
            enddo
         enddo
      enddo
    enddo
  endif

  riu=max(min(riu,10.0),-10.0)
  ! compute richardson numbers on bottom of T cells as average
  ! of four nearest richardson numbers on bottom of U cells.
  ! (do not consider land cells in the average... only active ones)

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        active_cells = Grd%umask(i,j,k+1) + Grd%umask(i-1,j,k+1) + Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln_chen
        rit(i,j,k)   = (riu(i,j,k) + riu(i-1,j,k) + riu(i,j-1,k) + riu(i-1,j-1,k))/active_cells

        ! make sure no static instability exists (one that is not seen
        ! by the Richardson number).  This may happen due to
        ! horizontal averaging used in calculating the Richardson
        ! number.

        if (rit(i,j,k) > 0.0 .and. wrk1(i,j,k) > 0.0) then
          rit(i,j,k) = -0.01
        endif
      enddo
    enddo
  enddo
  rit=max(min(rit,10.0),-10.0)

end subroutine ri_for_chen
! </SUBROUTINE> NAME="ri_for_chen"


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_chen_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_vert_chen_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  call save_restart(Che_restart, time_stamp)

end subroutine ocean_vert_chen_restart
! </SUBROUTINE> NAME="ocean_vert_chen_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_chen_end">
!
! <DESCRIPTION>
! Save the Kraus boundary layer depth to start next time step. 
! </DESCRIPTION>
!
subroutine ocean_vert_chen_end(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: stdoutunit 
  stdoutunit=stdout() 

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_vert_chen: ending hbl checksum for kraus.res.nc'
  call write_timestamp(Time%model_time)
  call write_chksum_2d('hbl', hbl(COMP)*Grd%tmask(COMP,1))

  return 

end subroutine ocean_vert_chen_end  
! </SUBROUTINE> NAME="ocean_vert_chen_end"


end module ocean_vert_chen_mod


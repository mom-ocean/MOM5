module ocean_lapcst_friction_mod
#define COMP isc:iec,jsc:jec
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! This module computes the thickness weighted and density weighted
! acceleration for horizontal velocity arising from horizontal 
! Laplacian friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness weighted and density weighted
! time tendency for horizontal velocity arising from horizontal
! Laplacian friction. 
!
! The viscosity used to determine the strength of the tendency 
! can be a general function of space yet it is constant in time.  
! A namelist option exists that determines this viscosity 
! as a local function of the grid spacing. 
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
!  R. J. Murray and C. J. C. Reason,
!  A curvilinear version of the Bryan-Cox ocean model
!  Journal of Computational Physics (2002), vol 171, pages 1--46
! </REFERENCE>
!
! <REFERENCE>
!  S.M. Griffies, Elements of MOM (2012)
! </REFERENCE>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! <NOTE>
! This scheme has been found to be faster than the Smagorinsky viscosity 
! scheme.  However, the algorithm here is less robust since it 
! contains null modes in the terms associated with the sphericity 
! of the earth.  Hence, there may be flow configurations that are 
! not dissipated.
! </NOTE>
!
! <NOTE>
! The sink from drag next to partial cells has been dropped 
! from MOM4p1.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lapcst_friction_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false. 
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging by printing checksums.  
!  </DATA> 
!
!  <DATA NAME="alap" UNITS="m^2/sec" TYPE="real">
!  This is the value for the space-time constant Laplacian viscosity. 
!  </DATA> 
!  <DATA NAME="velocity_mix_micom" TYPE="logical">
!  If .true., then the viscosity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field
use fms_mod,             only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog, read_data
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: input_nml_file, mpp_sum, mpp_pe, mpp_max, mpp_error

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: FAX, FAY, BAX, BAY, FMX, FMY
use ocean_operators_mod,  only: BDX_EU, BDY_NU, FDX_U, FDY_U, FDX_NT, FDY_ET
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type, ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d_u, diagnose_2d_u, write_chksum_3d
use ocean_workspace_mod,  only: wrk1_v, wrk1  

implicit none

public ocean_lapcst_friction_init
public lapcst_friction
public lapcst_viscosity_check
public lapcst_reynolds_check

private

! length of forward time step  
real :: dtime

! for diagnostics 
integer :: id_visc      =-1
integer :: id_lap_fric_u=-1
integer :: id_lap_fric_v=-1
logical :: used

real :: alap                  = 0e5     ! constant horz viscosity for momentum (m^2/sec)
real :: vel_micom             = 0.0     ! background scaling velocity for micom viscosity (m/sec) 
real :: alap_taper_start_lat  = 60.     ! for tapering viscosity at high latitudes
real :: alap_taper_end_lat    = 90.     ! for tapering viscosity at high latitudes
real :: alap_taper_min_fac    = 30.     ! for tapering viscosity at high latitudes
logical :: alap_taper_hi_lats = .false. ! for tapering viscosity at high latitudes
logical :: velocity_mix_micom =.false.  ! if true, viscosity made a function of the grid spacing

real, dimension(:,:,:), allocatable :: visc_ceu  ! viscosity on east face of U cells (m^2/s)
real, dimension(:,:,:), allocatable :: visc_cnu  ! viscosity on north face of U cells (m^2/s)
real, dimension(:,:,:), allocatable :: visc_cu   ! variable viscosity averaged to U cell (m^2/s)

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: ocean_lapcst_friction.f90 ($Id: ocean_lapcst_friction.F90,v 20.0 2013/12/14 00:14:24 fms Exp $)'
character (len=128) :: tagname = &
     '$Name: tikal $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.
logical :: read_alap              = .FALSE.

namelist /ocean_lapcst_friction_nml/ use_this_module, debug_this_module, alap,                     &
                                     velocity_mix_micom, vel_micom,                                &
                                     alap_taper_hi_lats, alap_taper_start_lat, alap_taper_end_lat, &
                                     alap_taper_min_fac, read_alap

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_lapcst_friction_init">
!
! <DESCRIPTION>
! Initialize the horizontal Laplacian friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_lapcst_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                      obc, use_lapcst_friction, debug)

  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_lapcst_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k, num

  integer, save :: unit=6     

  real    :: cbeta, alap_munk, dau_max
  real    :: alap_taper_min            
  real, dimension(Domain%isd:Domain%ied,Domain%jsd:Domain%jed) :: alap_const

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcst_friction_mod (ocean_lapcst_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_lapcst_friction_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_lapcst_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_lapcst_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_lapcst_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_lapcst_friction_nml)  
  write (stdlogunit,ocean_lapcst_friction_nml)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING ocean_lapcst_friction_mod')
      Ocean_options%horz_lap_friction = 'Used constant horizontal Laplacian friction.'
      use_lapcst_friction = .true. 
  else
      call mpp_error(NOTE, '==>Note: NOT using ocean_lapcst_friction_mod')
      Ocean_options%horz_lap_friction = 'Did NOT use horizontal Laplacian friction.'
      use_lapcst_friction = .false.
      return 
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  Grd => Grid
  Dom => Domain
  dtime = d_time

  write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_lapcst_friction_mod: using forward time step of (secs)', dtime 

  have_obc       = obc
  if (have_obc) write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_lapcst_friction_mod: considering obc' 

  alap_const=alap

  if (read_alap) then
      call read_data('INPUT/alap.nc','alap',alap_const,Domain%domain2d)
      call mpp_update_domains (alap_const, Dom%domain2d)      
  endif
  
  if (velocity_mix_micom) then
      alap_const(:,:) = vel_micom*(2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/&
           (epsln + Grd%dxu(:,:) + Grd%dyu(:,:))
  else if (.not. read_alap) then
      ! use alternate velocity scaling based on cell area
      dau_max = 0.0
      do j=jsc,jec
         do i=isc,iec
            if (Grd%dau(i,j) > dau_max) dau_max = Grd%dau(i,j)
         enddo
      enddo
      call mpp_max(dau_max)
      alap_const(:,:) = alap*(Grd%dau(:,:)/dau_max)
  endif

  ! enhance background viscosity near open boundaries if needed
  if (have_obc) call ocean_obc_enhance_visc_back(alap_const)
  
  if (alap_taper_hi_lats) then
    alap_const = alap
    alap_taper_min = alap/alap_taper_min_fac
    do j=jsc,jec
      do i=isc,iec
        if (abs(Grid%yu(i,j)) > alap_taper_start_lat) then
           alap_const(i,j) = alap - &
            (alap-alap_taper_min)*(alap_taper_start_lat-abs(Grid%yu(i,j)))**2 &
            /(alap_taper_start_lat-alap_taper_end_lat)**2
        endif
        if (alap_const(i,j) < alap_taper_min) alap_const(i,j) = alap_taper_min
      enddo
    enddo
    call mpp_update_domains (alap_const, Dom%domain2d)
    do j=jsc,jec
      write (unit,'(a,i5,f10.3,f14.4)') ' lateral_friction_init: j, lat, alap_const =', &
                      j+Dom%joff, Grid%yu(isc,j), alap_const(isc,j)
    enddo
  endif

  allocate (visc_ceu(isd:ied,jsd:jed,nk)) ! east face of U cell
  allocate (visc_cnu(isd:ied,jsd:jed,nk)) ! north face of U cell
  allocate (visc_cu(isd:ied,jsd:jed,nk))  ! for metric term
 
  do k=1,nk
    visc_ceu(:,:,k) = FAX(alap_const(:,:))
    visc_cnu(:,:,k) = FAY(alap_const(:,:))
    visc_cu(:,:,k)  = alap_const(:,:)
  enddo
  call mpp_update_domains (visc_ceu, Dom%domain2d)
  call mpp_update_domains (visc_cnu, Dom%domain2d)
  call mpp_update_domains (visc_cu,  Dom%domain2d)


  ! Munk boundary check
  if (dtime /= 0.0) then
    num = 0
    write (stdoutunit,'(/,(1x,a))') &
    '==> Warning: locations (if any) where the Munk boundary layer is unresolved.'
    do j=jsc,jec
      do i=isc,iec
        cbeta = 2.28e-11
        alap_munk = cbeta*(Grd%h1u(i,j)/radius)*(Grd%dxu(i,j)*sqrt(3.0)/pi)**3
        if (visc_cu(i,j,1) <= alap_munk .and. num <= 10 .and. i == isc) then
          num = num + 1
          write (stdoutunit,'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') ' at (i,j)= (',i,',',j,'),(lon,lat)= (', &
                           Grd%xt(i,j),',',Grd%yt(i,j),&
                           '),  "alap" = ', visc_cu(i,j,1),&
                           ' m^2/s. the critical value =',alap_munk,' m^2/s'
        endif
      enddo
    enddo
    call mpp_sum (num)
    if (num > 10) write (stdoutunit,*) ' (a max of 10 violations per PE are listed)'
  endif

  ! register diagnostics and send static viscosity 
  id_visc = register_static_field ('ocean_model', 'lap_visc', Grd%vel_axes_uv(1:2),&
            'laplacian viscosity', 'm^2/sec',missing_value=missing_value,          &
            range=(/-10.0,1.e10/),                                                 &
            standard_name='ocean_momentum_xy_laplacian_diffusivity')
  call diagnose_2d_u(Time, Grd, id_visc, visc_cu(:,:,1))

  id_lap_fric_u = register_diag_field('ocean_model','lap_fric_u',Grd%vel_axes_uv(1:3), &
                  Time%model_time,'Thickness and rho wghtd horz lap frict on u',       &
                  '(kg/m^2)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_lap_fric_v = register_diag_field('ocean_model','lap_fric_v',Grd%vel_axes_uv(1:3), &
                  Time%model_time,'Thickness and rho wghtd horz lap frict on v',       &
                  '(kg/m^2)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))


end subroutine ocean_lapcst_friction_init
! </SUBROUTINE>  NAME="ocean_lapcst_friction_init"


!#######################################################################
! <SUBROUTINE NAME="lapcst_friction">
!
! <DESCRIPTION>
! This routine computes the rho*thickness weighted time tendency for
! horizontal velocity from horizontal Laplacian friction.
! </DESCRIPTION>
!
subroutine lapcst_friction(Time, Thickness, Velocity, lap_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
  real, dimension(isd:,jsd:), intent(inout) :: lap_viscosity

  real, dimension(isd:ied,jsd:jed)   :: metric
  real, dimension(isd:ied,jsd:jed,2) :: uk
  real, dimension(isd:ied,jsd:jed)   :: tmp, tmp1, tmp2, m1, m2, n1, n2
  integer :: i, j, k, n, taum1, tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) then 

      if(energy_analysis_step) then 
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Velocity%wrkv(i,j,k,n) = 0.0
                   enddo
                enddo
             enddo
          enddo
      endif

      return 
  endif

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcst_friction_mod (lapcst_friction): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  metric(:,:)     = 0.0
  wrk1_v(:,:,:,:) = 0.0

  ! laplacian friction operator 
  do k=1,nk

     m1(:,:)   = 2.0*visc_cu(:,:,k)*Grd%dh1dy(:,:) + BAY(FDY_U(visc_cu(:,:,k)))
     m2(:,:)   = 2.0*visc_cu(:,:,k)*Grd%dh2dx(:,:) + BAX(FDX_U(visc_cu(:,:,k)))
     tmp1(:,:) = visc_ceu(:,:,k)*Grd%dyue(:,:)
     tmp2(:,:) = visc_cnu(:,:,k)*Grd%dxun(:,:)
     n1(:,:)   = -Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh2dx(:,:))) &
                 -Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh1dy(:,:)))
     n2(:,:)   = +Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh1dy(:,:))) &
                 -Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh2dx(:,:)))

     do n=1,2
        uk(:,:,:)   = Velocity%u(:,:,k,:,taum1)
        metric(:,:) = (3-2*n)*(m1(:,:)*FDX_NT(BAX(uk(:,:,3-n)))-m2(:,:)*FDY_ET(BAY(uk(:,:,3-n)))  )&
                             + n1(:,:)*uk(:,:,n) + (3-2*n)*n2(:,:)*uk(:,:,3-n)
        tmp(:,:) = ( BDX_EU(visc_ceu(:,:,k)*FMX(Thickness%rho_dzu(:,:,k,tau))*FDX_U(uk(:,:,n)))&
                    +BDY_NU(visc_cnu(:,:,k)*FMY(Thickness%rho_dzu(:,:,k,tau))*FDY_U(uk(:,:,n))) ) &
                    + metric(:,:)*Thickness%rho_dzu(:,:,k,tau)
        do j=jsc,jec
           do i=isc,iec
              wrk1_v(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)
           enddo
        enddo

     enddo

  enddo

  if(energy_analysis_step) then 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

  else 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      
      ! vertically averaged viscosity at U-cell centre
      lap_viscosity(:,:) = 0.0
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               lap_viscosity(i,j) = lap_viscosity(i,j) &
                + Grd%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)*visc_cu(i,j,k)
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            lap_viscosity(i,j) = Grd%umask(i,j,1)*lap_viscosity(i,j)/(epsln+Thickness%mass_u(i,j,tau))
         enddo
      enddo


      ! diagnostics
      call diagnose_3d_u(Time, Grd, id_lap_fric_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_lap_fric_v, wrk1_v(:,:,:,2))
  endif

  if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_lapcst_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('lapcst friction(1)', wrk1_v(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('lapcst friction(2)', wrk1_v(COMP,:,2)*Grd%umask(COMP,:))
  endif 


end subroutine lapcst_friction
! </SUBROUTINE> NAME="lapcst_friction"


!#######################################################################
! <SUBROUTINE NAME="lapcst_viscosity_check">
!
! <DESCRIPTION>
! Perform linear stability check for the Laplacian operator 
! given a value for the horizontal laplacian viscosity.
! </DESCRIPTION>

subroutine lapcst_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: crit, dsmin

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcst_friction_mod (lapcst_viscosity_check): needs initialization')
  endif 

  write (stdoutunit,'(/60x,a/)') &
  ' Excessive horizontal Laplacian friction summary:'
  write (stdoutunit,'(1x,a/)')   &
  'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then 
          dsmin = min(Grd%dxu(i,j),Grd%dxu(i+1,j),Grd%dyu(i,j),Grd%dyu(i,j+1))
          crit = 0.5*dsmin**2/max(dtime,epsln)
          if (visc_cu(i,j,k) > crit .and. num < max_num) then
            num = num + 1
            write (unit,9600) &
            'visc_cu(',visc_cu(i,j,k), crit, i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo

  call mpp_sum(num)
  if (num > max_num) then 
    write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'
  endif 

9600  format(/' Warning: ',a,es10.3,' m^2/s) exceeds max value (',es10.3,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')

end subroutine lapcst_viscosity_check
! </SUBROUTINE> NAME="lapcst_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="lapcst_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the Laplacian grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine lapcst_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lapcst_friction_mod (lapcst_reynolds_check): needs initialization')
  endif 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(visc_cu(i,j,k) + epsln)
        ramn = 1.0/(visc_cu(i,j,k) + epsln)
        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j))*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j))*ramn
        if (reyy > reyny) then
          ireyny = i
          jreyny = j
          kreyny = k
          reyny  = reyy
          reynv  = Velocity%u(i,j,k,2,tau)
          reynmv = 1.0/ramn
        endif
      enddo
    enddo
  enddo
  write (stdoutunit,'(/60x,a/)') ' Horizontal Laplacian Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, Grd%xu(ireynx,jreynx), &
                       Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, Grd%xu(ireyny,jreyny), &
                       Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine lapcst_reynolds_check
! </SUBROUTINE> NAME="lapcst_reynolds_check"


end module ocean_lapcst_friction_mod
      
      

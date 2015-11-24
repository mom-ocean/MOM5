module ocean_bihcst_friction_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted and density weighted time tendency for velocity
! from horizontal biharmonic friction 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the thickness and density weighted time tendency
! for horizontal velocity arising from horizontal biharmonic friction. 
! The viscosity used to determine the strength of the tendency 
! can be a general function of space yet it is constant in time.  
! A namelist option exists that determines this viscosity 
! as a local function of the grid spacing. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
!  R. J. Murray and C.J.C. Reason,
!  A curvilinear version of the Bryan-Cox ocean model
!  Journal Computational Physics (2002), vol 171, pages 1-46
! </REFERENCE>
!
! <REFERENCE>
!  Elements of MOM (2012), S.M. Griffies 
! </REFERENCE>
!
! <NOTE>
! The numerical implementation requires one call to mpp_update_domains if 
! running the model with halo=1.  
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
! The model can generally run with both Laplacian and biharmonic friction
! enabled at the same time. Such has been found useful for some eddying 
! ocean simulations. 
! </NOTE>
!
! <NOTE>
! The numerical implementation in mom4p1 and mom5 does not include the sink due 
! to partial bottom cells. This sink was implemented in mom2 and mom3
! and mom4.  It is not known how relevant this sink is.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bihcst_friction_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module.  Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging by printing checksums.  
!  </DATA> 
!
!  <DATA NAME="abih" UNITS="m^4/sec" TYPE="real">
!  This is the value for the space-time constant biharmonic viscosity. 
!  </DATA> 
!  <DATA NAME="velocity_mix_micom" TYPE="logical">
!  If .true., then the viscosity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity within a latitude
!  band surrounding the equator.  This is useful for some models of enhanced equatorial
!  resolution that can maintain numerical integrity in this region with less friction 
!  than outside the tropical band.  
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity is set according 
!  to eq_vel_micom.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: pi, radius, epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,          only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,  only: mpp_update_domains
use mpp_mod,          only: input_nml_file, mpp_error, mpp_pe, mpp_max, mpp_sum

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: FAX, FAY, BAX, BAY, FMX, FMY
use ocean_operators_mod,  only: BDX_EU, BDY_NU, FDX_U, FDY_U, FDX_NT, FDY_ET
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type, ocean_options_type
use ocean_util_mod,       only: write_timestamp, diagnose_3d_u, diagnose_2d_u, write_chksum_3d
use ocean_workspace_mod,  only: wrk1_v

implicit none

public ocean_bihcst_friction_init
public bihcst_friction
public bihcst_viscosity_check
public bihcst_reynolds_check
private delsq_velocity

private

! time step 
real ::  dtime = 0.0  ! forward time step for friction (2*dtuv if threelevel, dtuv if twolevel) 

! for diagnostics 
integer :: id_aiso      =-1
integer :: id_bih_fric_u=-1
integer :: id_bih_fric_v=-1
logical :: used

real, dimension(:,:,:,:), allocatable :: del2_vel  ! del^2 of horizontal velocity
real, dimension(:,:,:),   allocatable :: visc_ceu  ! viscosity on east face of U cells (m^4/s)
real, dimension(:,:,:),   allocatable :: visc_cnu  ! viscosity on north face of U cells (m^4/s)
real, dimension(:,:,:),   allocatable :: visc_cu   ! variable viscosity averaged to U cell (m^4/s)

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
  '=>Using: ocean_bihcst_friction.f90 ($Id: ocean_bihcst_friction.F90,v 20.0 2013/12/14 00:14:14 fms Exp $)'

character (len=128) :: tagname = &
  '$Name: tikal $'

logical :: module_is_initialized = .FALSE.
logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: have_obc              = .FALSE.
real ::    abih                  = 0.0    ! constant horz biharmonic viscosity (m^4/s)
real ::    vel_micom             = 0.0    ! background scaling velocity for micom viscosity (m/sec) 
real ::    eq_vel_micom          = 0.0    ! background scaling velocity (m/sec) w/i equatorial band for iso-visc
real ::    eq_lat_micom          = 0.0    ! equatorial latitude band for micom (degrees)
logical :: velocity_mix_micom    =.false. ! if true, diffusivity made a function of the grid spacing

namelist /ocean_bihcst_friction_nml/ use_this_module, debug_this_module, abih, &
                                     velocity_mix_micom, vel_micom, eq_vel_micom, eq_lat_micom

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bihcst_friction_init">
!
! <DESCRIPTION>
! Initialize the horizontal biharmonic friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bihcst_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                      obc, use_bihcst_friction, debug)
   
  type(ocean_grid_type),   target, intent(in)    :: Grid
  type(ocean_domain_type), target, intent(in)    :: Domain
  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_options_type),        intent(inout) :: Ocean_options
  real,                            intent(in)    :: d_time
  logical,                         intent(in)    :: obc
  logical,                         intent(inout) :: use_bihcst_friction
  logical,               optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real, dimension(Domain%isd:Domain%ied,Domain%jsd:Domain%jed) :: Biso

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bihcst_friction_init: module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_bihcst_friction_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bihcst_friction_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_bihcst_friction_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_bihcst_friction_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bihcst_friction_nml)  
  write (stdlogunit,ocean_bihcst_friction_nml)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING ocean_bihcst_friction_mod')
      Ocean_options%horz_bih_friction = 'Used constant horizontal biharmonic friction.'
      use_bihcst_friction=.true.
  else
      call mpp_error(NOTE, '==>Note: NOT using ocean_bihcst_friction_mod')
      Ocean_options%horz_bih_friction = 'Did NOT use horizontal biharmonic friction.'
      use_bihcst_friction=.false.
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

  allocate (visc_ceu(isd:ied,jsd:jed,nk))   ! east face of U cell
  allocate (visc_cnu(isd:ied,jsd:jed,nk))   ! north face of U cell
  allocate (visc_cu(isd:ied,jsd:jed,nk))    ! for metric term
  allocate (del2_vel(isd:ied,jsd:jed,nk,2)) ! laplacian of velocity on U cell

  write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_bihcst_friction_mod: using forward time step of (secs)', dtime 

  have_obc = obc
  if (have_obc) then 
    write(stdoutunit,'(/a,f10.2)')'==> Note from ocean_bihcst_friction_mod: considering obc' 
  endif 

  Biso(:,:) = 0.0
  if(velocity_mix_micom) then 
     do j=jsd,jed
        do i=isd,ied
           Biso(i,j) = vel_micom*((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           if(abs(Grd%yu(i,j)) < eq_lat_micom) then
             Biso(i,j) = eq_vel_micom*((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           endif
       enddo  
     enddo
  else
    Biso(:,:) = abs(abih)
  endif 

  ! enhance background viscosity near open boundaries if needed
  if (have_obc) call ocean_obc_enhance_visc_back(Biso, 'bih')

  do k=1,nk
    visc_ceu(:,:,k) = FAX(Biso(:,:))  
    visc_cnu(:,:,k) = FAY(Biso(:,:))  
    visc_cu(:,:,k)  = Biso(:,:)          
  enddo
  call mpp_update_domains (visc_ceu, Dom%domain2d)
  call mpp_update_domains (visc_cnu, Dom%domain2d)
  call mpp_update_domains (visc_cu,  Dom%domain2d)


  ! diagnostics
  id_aiso = register_static_field ('ocean_model', 'aiso_bih', Grd%vel_axes_uv(1:2),      &
            'bih viscosity','m^4/sec',missing_value=missing_value, range=(/-10.0,1.e20/),&
            standard_name='ocean_momentum_xy_biharmonic_diffusivity')
  call diagnose_2d_u(Time, Grd, id_aiso, visc_cu(:,:,1))

  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%vel_axes_uv(1:3), &
                  Time%model_time,'Thickness and rho wghtd horz bih-frict on u',          &
                  '(kg/m^3)*(m^2/s^2)',missing_value=missing_value, range=(/-10.0,1.e10/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time,'Thickness and rho wghtd horz bih-frict on v',          &
                  '(kg/m^3)*(m^2/s^2)',missing_value=missing_value, range=(/-10.0,1.e10/))


end subroutine ocean_bihcst_friction_init
! </SUBROUTINE>  NAME="ocean_bihcst_friction_init"


!#######################################################################
! <SUBROUTINE NAME="bihcst_friction">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted acceleration on
! horizontal velocity arising from horizontal biharmonic friction.
! </DESCRIPTION>
!
subroutine bihcst_friction(Time, Thickness, Velocity, bih_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  real, dimension(isd:,jsd:), intent(inout) :: bih_viscosity
  logical,                    intent(in)    :: energy_analysis_step

  real, dimension(isd:ied,jsd:jed)   :: metric
  real, dimension(isd:ied,jsd:jed,2) :: uk
  real, dimension(isd:ied,jsd:jed)   :: tmp, tmp1, tmp2, m1, m2, n1, n2
  integer :: i, j, k, n, tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcst_friction_mod (horz_bih_friction): module needs initialization')
  endif 

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

  tau = Time%tau 

  metric(:,:)     = 0.0
  wrk1_v(:,:,:,:) = 0.0

  call delsq_velocity (Time, Thickness, Velocity)
  if(Dom%xhalo==1 .or. Dom%yhalo==1 ) call mpp_update_domains (del2_vel(:,:,:,:), Dom%domain2d)

  ! Murray and Reason Eqn 29

  do k=1,nk

    ! metric terms for generalized coordinate system
    m1(:,:) = 2.0*visc_cu(:,:,k)*Grd%dh1dy(:,:) + BAY(FDY_U(visc_cu(:,:,k)))
    m2(:,:) = 2.0*visc_cu(:,:,k)*Grd%dh2dx(:,:) + BAX(FDX_U(visc_cu(:,:,k)))
    tmp1(:,:) = visc_ceu(:,:,k)*Grd%dyue(:,:)
    tmp2(:,:) = visc_cnu(:,:,k)*Grd%dxun(:,:)
    n1(:,:) = - Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh2dx(:,:))) &
              - Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh1dy(:,:)))
    n2(:,:) = + Grd%dyur(:,:)*BDX_EU(tmp1(:,:)*FAX(Grd%dh1dy(:,:))) &
              - Grd%dxur(:,:)*BDY_NU(tmp2(:,:)*FAY(Grd%dh2dx(:,:)))

    do n=1,2
       uk(:,:,:) = -del2_vel(:,:,k,:)
       metric(:,:) = (3-2*n)*(m1(:,:)*FDX_NT(BAX(uk(:,:,3-n)))  -&
                              m2(:,:)*FDY_ET(BAY(uk(:,:,3-n)))  )&
            +  n1(:,:)*uk(:,:,n) + (3-2*n)*n2(:,:)*uk(:,:,3-n)
       tmp(:,:) = ( BDX_EU(visc_ceu(:,:,k)*FMX(Thickness%rho_dzu(:,:,k,tau))*FDX_U(uk(:,:,n)))   &
                   +BDY_NU(visc_cnu(:,:,k)*FMY(Thickness%rho_dzu(:,:,k,tau))*FDY_U(uk(:,:,n))) ) &
            +metric(:,:)*Thickness%rho_dzu(:,:,k,tau)

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
      bih_viscosity(:,:) = 0.0
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               bih_viscosity(i,j) = bih_viscosity(i,j) &
                + Grd%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)*visc_cu(i,j,k)
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            bih_viscosity(i,j) = Grd%umask(i,j,1)*bih_viscosity(i,j)/(epsln+Thickness%mass_u(i,j,tau))
         enddo
      enddo

      ! diagnostics 
      call diagnose_3d_u(Time, Grd, id_bih_fric_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_bih_fric_v, wrk1_v(:,:,:,2))
  endif

  if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_bihcst_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_dzu(tau)', Thickness%rho_dzu(COMP,:,tau)*Grd%umask(COMP,:))
      call write_chksum_3d('bihcst friction(1)', wrk1_v(COMP,:,1)*Grd%umask(COMP,:))
      call write_chksum_3d('bihcst friction(2)', wrk1_v(COMP,:,2)*Grd%umask(COMP,:))
  endif 


end subroutine bihcst_friction
! </SUBROUTINE> NAME="bihcst_friction"


!#######################################################################
! <SUBROUTINE NAME="delsq_velocity">
!
! <DESCRIPTION>
! Subroutine computes the laplacian operator acting on velocity.
! </DESCRIPTION>
!
subroutine delsq_velocity (Time, Thickness, Velocity)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type),  intent(in) :: Velocity

  real, dimension(isd:ied,jsd:jed)       :: fe, fn, metric, n1, n2, tmp
  integer                                :: i, j, k, n, taum1, tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihcst_friction_mod (delsq_velocity): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  metric(:,:)  = 0.0

  ! metric terms for generalized system
  n1(:,:) = - Grd%dyur(:,:)*BDX_EU(Grd%dyue(:,:)*FAX(Grd%dh2dx(:,:)))      &
            - Grd%dxur(:,:)*BDY_NU(Grd%dxun(:,:)*FAY(Grd%dh1dy(:,:)))
  n2(:,:) = + Grd%dyur(:,:)*BDX_EU(Grd%dyue(:,:)*FAX(Grd%dh1dy(:,:)))      &
            - Grd%dxur(:,:)*BDY_NU(Grd%dxun(:,:)*FAY(Grd%dh2dx(:,:)))
  do k=1,nk

    do n=1,2

      ! Viscous flux across east and north face of U cells
      fe(:,:) = FDX_U(Velocity%u(:,:,k,n,taum1))*FMX(Thickness%rho_dzu(:,:,k,tau))
      fn(:,:) = FDY_U(Velocity%u(:,:,k,n,taum1))*FMY(Thickness%rho_dzu(:,:,k,tau))

      metric(:,:) = (3-2*n)*2.0*(Grd%dh1dy(:,:)*FDX_NT(BAX(Velocity%u(:,:,k,3-n,taum1)))-&
                                 Grd%dh2dx(:,:)*FDY_ET(BAY(Velocity%u(:,:,k,3-n,taum1)))  )& 
       + n1(:,:)*Velocity%u(:,:,k,n,taum1) + (3-2*n)*n2(:,:)*Velocity%u(:,:,k,3-n,taum1)

      ! note: do not use rho_dzur since rho_dzur is at time taup1
      tmp(:,:) =  BDX_EU(fe(:,:)) + BDY_NU(fn(:,:))
      do j=jsc,jec
         do i=isc,iec  
            del2_vel(i,j,k,n) = ( tmp(i,j)/Thickness%rho_dzu(i,j,k,tau) + metric(i,j) )*Grd%umask(i,j,k) 
         enddo
      enddo

    enddo
  enddo 

end subroutine delsq_velocity
! </SUBROUTINE> NAME="delsq_tracer"


!#######################################################################
! <SUBROUTINE NAME="bihcst_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the biharmonic
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine bihcst_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: dsmin, crit

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 
 
  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bihcst_friction_mod (bih_viscosity_check): needs initialization')
  endif 

  write (stdoutunit,'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdoutunit,'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then            
          dsmin = min(Grd%dxu(i,j),Grd%dxu(i+1,j),Grd%dyu(i,j),Grd%dyu(i,j+1))
          crit = 0.0625*dsmin**4/max(dtime,epsln)
          if (visc_cu(i,j,k) > crit .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'visc_cu(',visc_cu(i,j,k), crit, i, j, k, Grd%xu(i,j), Grd%yu(i,j), Grd%zt(k), mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es10.3,' m^4/s) exceeds max value (',es10.3,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine bihcst_viscosity_check
! </SUBROUTINE> NAME="bihcst_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="bihcst_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the biharmonic grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine bihcst_reynolds_check(Time, Velocity)

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
    call mpp_error(FATAL, '==>Error in ocean_bihcst_friction_mod (bih_reynolds_check): needs initialization')
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

        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j)**3)*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j)**3)*ramn
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
  write (stdoutunit,'(/60x,a/)') ' Horizontal biharmonic Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, &
                       Grd%xu(ireynx,jreynx), Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, &
                       Grd%xu(ireyny,jreyny), Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine bihcst_reynolds_check
! </SUBROUTINE> NAME="bihcst_reynolds_check"


end module ocean_bihcst_friction_mod
      
      

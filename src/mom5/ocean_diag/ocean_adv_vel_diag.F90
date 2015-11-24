module ocean_adv_vel_diag_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Numerical diagnostics for advection velocity related quantities. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Numerical diagnostics for advection velocity related quantities.
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_adv_vel_diag_nml">
!  <DATA NAME="max_cfl_value" UNITS="dimensionless" TYPE="real">
!  Critical value for Courant number, above which the model will be brought down.
!  </DATA> 
!  <DATA NAME="large_cfl_value" UNITS="dimensionless" TYPE="real">
!  Large value for Courant number, above which will write some diagnostics. 
!  </DATA> 
!  <DATA NAME="verbose_cfl" TYPE="logical">
!  For printing out lots of information about regions of large Courant numbers.
!  </DATA> 
!
!  <DATA NAME="diag_step" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, send_data, need_data
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, stdout, stdlog
use mpp_mod,             only: input_nml_file, mpp_error, mpp_max, mpp_pe
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,    only: time_type, increment_time

use ocean_domains_mod,    only: get_local_indices
use ocean_operators_mod,  only: BAX, BAY, BDX_ET, BDY_NT, BDX_EU, BDY_NU
use ocean_operators_mod,  only: REMAP_ET_TO_EU, REMAP_NT_TO_NU, REMAP_BT_TO_BU
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID
use ocean_parameters_mod, only: missing_value, rho0r 
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,      only: ocean_adv_vel_type, ocean_density_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type
use ocean_util_mod,       only: matrix, diagnose_2d, diagnose_3d, diagnose_3d_u
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_2d

implicit none

private

#include <ocean_memory.h>

real :: dtts
real :: dtuv
real :: dtime_t
real :: dtime_u

integer :: index_temp
integer :: diag_step=-1

! for output
integer :: unit=6

! for diagnostics clocks 
integer :: id_adv_vel_numerics
integer :: id_transport_on_s
integer :: id_transport_on_rho
integer :: id_transport_on_nrho
integer :: id_transport_on_theta

! for CFL checks 
real    :: max_cfl_value   = 100.0           
real    :: large_cfl_value = 10.0 
logical :: verbose_cfl     =.false.

! for transport_on_s 
integer :: id_tx_trans       =-1
integer :: id_ty_trans       =-1
integer :: id_tz_trans       =-1
integer :: id_tx_trans_int_z =-1
integer :: id_ty_trans_int_z =-1
integer :: id_tz_trans_int_z =-1
integer :: id_tz_trans_sq    =-1
integer :: id_ux_trans       =-1
integer :: id_uy_trans       =-1
integer :: id_uz_trans       =-1

! for transport_on_rho
integer :: id_tx_trans_rho=-1
integer :: id_ty_trans_rho=-1

! for transport_on_nrho
integer :: id_tx_trans_nrho=-1
integer :: id_ty_trans_nrho=-1

! for transport_on_theta
integer :: id_tx_trans_theta=-1
integer :: id_ty_trans_theta=-1

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

logical :: used

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

! for Bgrid or Cgrid
integer :: horz_grid

#ifdef MOM_STATIC_ARRAYS
 real, dimension(isd:ied,jsd:jed)  :: tmp 
#else
 real, dimension(:,:), allocatable :: tmp 
#endif

logical :: module_is_initialized = .FALSE.

character(len=128) :: version=&
     '$Id: ocean_adv_vel_diag.F90,v 20.0 2013/12/14 00:12:49 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

public ocean_adv_vel_diag_init
public ocean_adv_vel_diagnostics 

private remapping_check 
private cfl_check1
private cfl_check2
private maximum_bottom_w
private vertical_reynolds_check
public max_continuity_error

private transport_on_s
private transport_on_rho
private transport_on_nrho
private transport_on_theta

namelist /ocean_adv_vel_diag_nml/ max_cfl_value, large_cfl_value, verbose_cfl, diag_step


contains

!#######################################################################
! <SUBROUTINE NAME="ocean_adv_vel_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_adv_vel_diag module containing subroutines
! diagnosing advection velocity related properties of the simulation.
! </DESCRIPTION>
!
subroutine ocean_adv_vel_diag_init(Grid, Domain, Time, Time_steps, T_prog, Dens, hor_grid, cmip_units)

type(ocean_grid_type),    target, intent(in) :: Grid
type(ocean_domain_type),  target, intent(in) :: Domain
type(ocean_time_type),            intent(in) :: Time
type(ocean_time_steps_type),      intent(in) :: Time_steps
type(ocean_prog_tracer_type),     intent(in) :: T_prog(:)
type(ocean_density_type),         intent(in) :: Dens
integer,                          intent(in) :: hor_grid
logical,                          intent(in) :: cmip_units

integer :: n, ioun, io_status, ierr

integer :: stdoutunit,stdlogunit 
stdoutunit=stdout();stdlogunit=stdlog() 

if (module_is_initialized) return

module_is_initialized = .TRUE.

call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_adv_vel_diag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_adv_vel_diag_nml')
#else
ioun = open_namelist_file()
read(ioun, ocean_adv_vel_diag_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_adv_vel_diag_nml')
call close_file(ioun)
#endif
write (stdoutunit,'(/)')
write (stdoutunit, ocean_adv_vel_diag_nml)
write (stdlogunit, ocean_adv_vel_diag_nml)

if (diag_step == 0) diag_step = 1

#ifndef MOM_STATIC_ARRAYS
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
allocate ( tmp(isd:ied,jsd:jed) )
#endif

tmp = 0.0
Dom => Domain
Grd => Grid

dtuv    = Time_steps%dtuv
dtts    = Time_steps%dtts
dtime_t = Time_steps%dtime_t
dtime_u = Time_steps%dtime_u
horz_grid = hor_grid

! set index for potential temperature 
do n=1,size(T_prog,1)
  if(T_prog(n)%name == 'temp' ) then 
    index_temp = n
  endif 
enddo 

! for units placed on diagnosed transport 
if(cmip_units) then
  transport_convert=1.0
  transport_dims   = 'kg/s'
else
  transport_convert=1.0e-9 
  transport_dims   = 'Sv (10^9 kg/s)'
endif 



! register fields for diagnostic output

id_tx_trans = register_diag_field ('ocean_model','tx_trans', Grd%tracer_axes_flux_x(1:3),&
              Time%model_time, 'T-cell i-mass transport',trim(transport_dims),           &
              missing_value=missing_value, range=(/-1e20,1e20/),                         &
              standard_name='ocean_x_mass_transport')
id_ty_trans = register_diag_field ('ocean_model','ty_trans', Grd%tracer_axes_flux_y(1:3), &
              Time%model_time, 'T-cell j-mass transport',trim(transport_dims),            &
              missing_value=missing_value, range=(/-1e20,1e20/),                          &
              standard_name='ocean_y_mass_transport')
id_tx_trans_int_z = register_diag_field ('ocean_model','tx_trans_int_z',        &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                     &
              'T-cell i-mass transport vertically summed',trim(transport_dims), &
              missing_value=missing_value, range=(/-1e20,1e20/))
id_ty_trans_int_z = register_diag_field ('ocean_model','ty_trans_int_z',        &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                     &
              'T-cell j-mass transport vertically summed',trim(transport_dims), &
              missing_value=missing_value, range=(/-1e20,1e20/))

id_tz_trans = register_diag_field ('ocean_model','tz_trans', Grd%tracer_axes_wt(1:3), &
              Time%model_time, 'T-cell k-mass transport',trim(transport_dims),        &
              missing_value=missing_value, range=(/-1e20,1e20/),                      &
              standard_name='upward_ocean_mass_transport')
id_tz_trans_int_z = register_diag_field ('ocean_model','tz_trans_int_z',        &
              Grd%tracer_axes_wt(1:2), Time%model_time,                         &
              'T-cell k-mass transport vertically summed',trim(transport_dims), &
              missing_value=missing_value, range=(/-1e20,1e20/))
id_tz_trans_sq = register_diag_field ('ocean_model','tz_trans_sq', Grd%tracer_axes_wt(1:3),&
              Time%model_time, 'Square of T-cell k-mass transport (positive upwards)'      &
              ,trim(transport_dims)//'^2',                                                 &
              missing_value=missing_value, range=(/-1e20,1e20/),                           &
              standard_name='square_of_upward_ocean_mass_transport')

id_ux_trans = register_diag_field ('ocean_model','ux_trans', Grd%vel_axes_flux_x(1:3), &
              Time%model_time, 'U-cell i-mass transport',trim(transport_dims),         &
              missing_value=missing_value, range=(/-1e20,1e20/))
id_uy_trans = register_diag_field ('ocean_model','uy_trans', Grd%vel_axes_flux_y(1:3), &
              Time%model_time, 'U-cell j-mass transport',trim(transport_dims),         &
              missing_value=missing_value, range=(/-1e20,1e20/))
id_uz_trans = register_diag_field ('ocean_model','uz_trans', Grd%vel_axes_wu(1:3), &
              Time%model_time, 'U-cell k-mass transport',trim(transport_dims),     &
              missing_value=missing_value, range=(/-1e20,1e20/))

id_tx_trans_rho = register_diag_field ('ocean_model','tx_trans_rho', Dens%potrho_axes_flux_x(1:3),&
                  Time%model_time, 'T-cell i-mass transport on pot_rho',trim(transport_dims),     &
                  missing_value=missing_value, range=(/-1e20,1e20/))
id_ty_trans_rho = register_diag_field ('ocean_model','ty_trans_rho', Dens%potrho_axes_flux_y(1:3),&
                  Time%model_time, 'T-cell j-mass transport on pot_rho',trim(transport_dims),     &
                  missing_value=missing_value, range=(/-1e20,1e20/))

id_tx_trans_nrho = register_diag_field ('ocean_model','tx_trans_nrho', Dens%neutralrho_axes_flux_x(1:3),&
                   Time%model_time, 'T-cell i-mass transport on neutral rho',trim(transport_dims),      &
                   missing_value=missing_value, range=(/-1e20,1e20/))
id_ty_trans_nrho = register_diag_field ('ocean_model','ty_trans_nrho', Dens%neutralrho_axes_flux_y(1:3),&
                   Time%model_time, 'T-cell j-mass transport on neutral rho',trim(transport_dims),      &
                   missing_value=missing_value, range=(/-1e20,1e20/))

id_tx_trans_theta = register_diag_field ('ocean_model','tx_trans_theta', Dens%theta_axes_flux_x(1:3),&
                    Time%model_time, 'T-cell i-mass transport on theta',trim(transport_dims),        &
                    missing_value=missing_value, range=(/-1e20,1e20/))
id_ty_trans_theta = register_diag_field ('ocean_model','ty_trans_theta', Dens%theta_axes_flux_y(1:3),&
                    Time%model_time, 'T-cell j-mass transport on theta',trim(transport_dims),        &
                    missing_value=missing_value, range=(/-1e20,1e20/))

! set ids for clocks
id_adv_vel_numerics   = mpp_clock_id('(Ocean adv_vel_diag: numerics)'    ,grain=CLOCK_ROUTINE)
id_transport_on_s     = mpp_clock_id('(Ocean adv_vel_diag: s-trans)'     ,grain=CLOCK_ROUTINE)
id_transport_on_rho   = mpp_clock_id('(Ocean adv_vel_diag: rho-trans)'   ,grain=CLOCK_ROUTINE)
id_transport_on_nrho  = mpp_clock_id('(Ocean adv_vel_diag: nrho-trans)'  ,grain=CLOCK_ROUTINE)
id_transport_on_theta = mpp_clock_id('(Ocean adv_vel_diag: theta-trans)' ,grain=CLOCK_ROUTINE)

if(horz_grid == MOM_BGRID) then
  call remapping_check
endif

end subroutine ocean_adv_vel_diag_init
! </SUBROUTINE>  NAME="ocean_adv_vel_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_adv_vel_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the velocity. 
! </DESCRIPTION>

subroutine ocean_adv_vel_diagnostics(Time, Thickness, Adv_vel, T_prog, Dens, visc_cbt)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: visc_cbt
  
  call mpp_clock_begin(id_adv_vel_numerics)
  if (diag_step > 0) then
      if (mod(Time%itt,diag_step) == 0) then
          call cfl_check1(Time, Thickness, Adv_vel)
          call cfl_check2(Time, Thickness, Adv_vel)
          call maximum_bottom_w(Adv_vel)
          call max_continuity_error(Adv_vel, Thickness)
          call vertical_reynolds_check(Adv_vel, Thickness, visc_cbt)
      endif
  endif
  call mpp_clock_end(id_adv_vel_numerics)
  
  call mpp_clock_begin(id_transport_on_s)
    call transport_on_s(Time, Adv_vel)
  call mpp_clock_end(id_transport_on_s)

  call mpp_clock_begin(id_transport_on_rho)
    call transport_on_rho(Time, Dens, Adv_vel)
  call mpp_clock_end(id_transport_on_rho)
  
  call mpp_clock_begin(id_transport_on_nrho)
    call transport_on_nrho(Time, Dens, Adv_vel)
  call mpp_clock_end(id_transport_on_nrho)
  
  call mpp_clock_begin(id_transport_on_theta)
    call transport_on_theta(Time, Dens, T_prog(index_temp), Adv_vel)
  call mpp_clock_end(id_transport_on_theta)


end subroutine ocean_adv_vel_diagnostics
! </SUBROUTINE>  NAME="ocean_adv_vel_diagnostics"



!#######################################################################
! <SUBROUTINE NAME="remapping_check">
!
! <DESCRIPTION>
! Compute remapping error for mapping from T to U cell fields.
!
! This error is only relevant for Bgrid calculations, where
! we need to remap some of the T-grid transports to U-grid.
! The C-grid does not have such remapping operations.
!
! The remapping error will be roundoff only for model
! grids where the tracer and velocity grid cell distances are 
! linearly related.  The spherical version of MOM satisfies the
! appropriate relation, and so should maintain roundoff for the 
! remapping error.  The tripolar version of MOM does not have 
! tracer and velocity grids related linearly, and so the 
! "remapping error" is nontrivial.  The significance of this error
! is unclear.  No adverse effects have been identified.
! </DESCRIPTION>
!
subroutine remapping_check

  real, dimension(isd:ied,jsd:jed) :: uet, vnt, ueu, vnu, u, v, dw, err
  integer                          :: i, j, ierr, jerr
  real                             :: remap_err, remap_err0, fudge

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_adv_vel_diag_mod (remapping_check) : module needs initialization ')
  endif 

  u(:,:)   = 1.0*Grd%umask(:,:,1)
  v(:,:)   = 1.0*Grd%umask(:,:,1)
  uet(:,:) = BAY(u(:,:)*Grd%dyu(:,:))/Grd%dyte(:,:)
  vnt(:,:) = BAX(v(:,:)*Grd%dxu(:,:))/Grd%dxtn(:,:)
  dw(:,:)  = REMAP_BT_TO_BU(BDX_ET(uet(:,:)) + BDY_NT(vnt(:,:)))
  
  ueu(:,:) = REMAP_ET_TO_EU(uet(:,:))  
  vnu(:,:) = REMAP_NT_TO_NU(vnt(:,:))
  err(:,:) = abs(dw(:,:)-(BDX_EU(ueu(:,:))   + BDY_NU(vnu(:,:))))

  ! used to distinguish processors when result is independent of processor
  fudge = 1 + 1.e-12*mpp_pe() 

  remap_err=0.0; ierr=isc; jerr=jsc
  do j=jsc,jec
    do i=isc,iec
      if (abs(err(i,j)) > abs(remap_err)) then 
        remap_err  = err(i,j)
        ierr = i
        jerr = j
      endif
    enddo
  enddo
  remap_err  = remap_err*fudge
  remap_err0 = remap_err
  remap_err  = abs(remap_err)
  call mpp_max(remap_err)

  if (abs(remap_err0) == remap_err) then
    remap_err = remap_err0
    write (unit,9000)  remap_err/fudge, ierr+Dom%ioff, jerr+Dom%joff, Grd%xu(ierr,jerr), Grd%yu(ierr,jerr)
    if(Grd%tripolar) write (unit,'(a)') &
    '==>Note: T-->U remapping error will be small (i.e., order 1e-20) only for spherical grids.'
  endif

9000  format(//'==>Maximum T-->U remapping error = ',es10.3,' m/s  at (i,j) = ','(',i4,',',i4,'),',' (lon,lat) = (',f7.2,',',f7.2,')')

end subroutine remapping_check
! </SUBROUTINE> NAME="remapping_check"


!#######################################################################
! <SUBROUTINE NAME="cfl_check1">
!
! <DESCRIPTION>
! Perform the first of two CFL checks for vertical velocity component. 
!
! Vectorized version from Russell.Fiedler@csiro.au computes cfl
! values at a single latitude. The location of the maximum at this
! latitude is calculated via the maxloc() intrinsic. The maximum 
! value for this processor is then updated if necessary.
! 
! </DESCRIPTION>
subroutine cfl_check1(Time, Thickness, Adv_vel)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_adv_vel_type),   intent(in) :: Adv_vel

  real     :: cflwtp                ! percent of cfl criteria reached by vertical velocity component
  real     :: cflwtm                ! vertical velocity component which comes closest to its cfl criteria
  integer  :: icflwt,jcflwt,kcflwt  ! "i","j","k" coordinate of "cflwtm"
  real     :: dtmax
  real     :: cflwtp0
  real     :: fudge
  integer  :: i, j, k, tau

  ! work array containing cfl values 
  real , dimension(isc:iec) :: tempcfl

  ! array required for maxloc() intrinsic function
  integer, dimension(1) :: itemp

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_adv_vel_diag_mod (cfl_check1): module needs initialization ')
  endif 

  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Vertical CFL summary I:  '
  write (stdoutunit,'(1x,a/)')'Location of largest vertical CFL.'

  ! to distinguish processors when tracer is independent of processor
  fudge = 1 + 1.e-12*mpp_pe() 

  icflwt=isc; jcflwt=jsc; kcflwt=1; cflwtp=epsln; cflwtm=0.0
  dtmax = max(dtime_t, dtime_u)  

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tempcfl(i)=abs(Grd%tmask(i,j,k)*Adv_vel%wrho_bt(i,j,k)/Thickness%dzwt(i,j,k))
        enddo
        itemp=maxloc(tempcfl)+isc-1
        if(tempcfl(itemp(1))>cflwtp) then
            cflwtp=tempcfl(itemp(1))
            icflwt=itemp(1)
            jcflwt=j
            kcflwt=k
        endif
     enddo
  enddo

  ! multiply by 100.0 to convert to percentages 
  ! multiply by dtmax to convert to dimensionless CFL number 
  ! multipby cflwtp by rho0r for dimensional reasons as well
  cflwtp = 100.0*cflwtp*rho0r*dtmax
  cflwtm = Adv_vel%wrho_bt(icflwt,jcflwt,kcflwt)*rho0r
  cflwtp  = cflwtp*fudge
  cflwtp0 = cflwtp
  call mpp_max(cflwtp)

  if (cflwtp == cflwtp0) then
    write (unit,'(a,g10.3,a,f7.2,a,g10.3,a,i4,a,i4,a,i3,a,a,f12.3,a,f12.3,a,f10.3,a)')&
   ' w_bt (',cflwtm,' m/s) is ',cflwtp/fudge,' % of CFL (',abs(100.0*cflwtm/(cflwtp/fudge)),&
   ' m/s) at (i,j,k) = (',icflwt+Dom%ioff,',',jcflwt+Dom%joff,',',kcflwt,'),',&
   ' (lon,lat,dpt) = (',Grd%xt(icflwt,jcflwt),',',Grd%yt(icflwt,jcflwt),',',  &
   Thickness%depth_zwt(icflwt,jcflwt,kcflwt),' m)'
    write (unit,'(a,e12.3,a,e12.3,a,f10.3)')&
   ' where grid is (m) (dxte,dytn,dzwu) = (',&
   Grd%dxte(icflwt,jcflwt),',',Grd%dytn(icflwt,jcflwt),',',Thickness%dzwu(icflwt,jcflwt,kcflwt)
   write(unit,'(a,i6)') ' The number of cells in the column are kmt = ', Grd%kmt(icflwt,jcflwt)
   write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
         Thickness%dst(icflwt,jcflwt,kcflwt),', ',Thickness%rho_dzt(icflwt,jcflwt,kcflwt,tau)
  endif

end subroutine cfl_check1
! </SUBROUTINE> NAME="cfl_check1"


!#######################################################################
! <SUBROUTINE NAME="cfl_check2">
!
! <DESCRIPTION>
! Perform the second of two vertical CFL checks.  
!
! Bring the model down if too many large Courant numbers detected.  
! </DESCRIPTION>
subroutine cfl_check2(Time, Thickness, Adv_vel)

 type(ocean_time_type),      intent(in) :: Time
 type(ocean_thickness_type), intent(in) :: Thickness
 type(ocean_adv_vel_type),   intent(in) :: Adv_vel

  real :: dtmax, cflw
  real :: wmax, pcflw, scl
  real, dimension(isd:ied,jsd:jed) :: w
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error ocean_adv_vel_diag_mod (cfl_check2): module needs initialization ')
  endif 

  tau = Time%tau   

  write (stdoutunit,'(/60x,a/)') ' Vertical CFL summary II:  '
  write (stdoutunit,'(1x,a,f5.2/)') &
  'Locations (if any) where vertical Courant number exceeds ', large_cfl_value

  dtmax = max(dtime_t, dtime_u)  
  do k=1,nk
    w(:,:) = rho0r*Adv_vel%wrho_bt(:,:,k)
    do j=jsc,jec
      do i=isc,iec
        cflw  = abs((dtmax/Thickness%dzwt(i,j,k))*w(i,j))*Grd%tmask(i,j,k)

        if (cflw >= large_cfl_value) then

          write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')      &
           ' Note: CFL velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',large_cfl_value 

          wmax  = Thickness%dzwt(i,j,k)/dtmax
          pcflw = abs(100.0*w(i,j)/wmax)
          write (unit,'(a,f8.2,a,g15.8,a)') ' w_bt reached', pcflw,' % of the local CFL limit (',wmax,' m/s)'
          write (unit,'(a,e12.3,a,e12.3,a,f10.3)') &
          ' where the grid cell has dimensions (metres) (dxt,dyt,dzwt) = (',&
          Grd%dxt(i,j),',',Grd%dyt(i,j),',',Thickness%dzwt(i,j,k)
          write(unit,'(a,i6)') ' The number of cells in the column are kmt = ', Grd%kmt(i,j)
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
          Thickness%dst(i,j,k),', ',Thickness%rho_dzt(i,j,k,tau)
        endif 

        if (cflw >= max_cfl_value) then

           write (unit,'(/,a,i4,a1,i3,a,i3,a,f6.3)')     &
           ' Note: CFL vertical velocity limit exceeded at (i,j,k) = (',&
           i+Dom%ioff, ',',j+Dom%joff,',',k,') by factor =',max_cfl_value 

           write(*,'(/a)') '==>Error in ocean_adv_vel_diag: max_cfl_value exceeded.'
           write(*,'(a)')  '   The model is bringing itself down in order for the user to determine '
           write(*,'(a)')  '   whether the simulation is going to produce Inf or NaN, or whether it '
           write(*,'(a)')  '   will successfully run through the present phase with large Courant numbers.'
           write(*,'(a)')  '   It is often the case that early spin-ups produce very large Courant'
           write(*,'(a)')  '   numbers, with the model later adjusting to stable behaviour with modest'
           write(*,'(a)')  '   Courant numbers. To test if the model will reach stable behaviour,'
           write(*,'(a)')  '   try significally increasing "max_cfl_value" in ocean_adv_vel_diag_nml.'

          if(verbose_cfl) then 

              write(*,'(/a,f6.3/)') &
              ' Information about region where vertical Courant number is larger than',max_cfl_value
              is = max(isc,i-3)
              ie = min(iec,i+3)
              js = max(jsc,j-3)
              je = min(jec,j+3)
              scl = 0.0
              write (*,9100)  'w_bt (m/s)', Time%itt, k, Thickness%depth_zwt(is,js,k), &
                   Grd%xt(is,j), Grd%xt(ie,j), Grd%yt(i,js), Grd%yt(i,je), scl
              call matrix (w(is:ie,js:je), is, ie, js, je, scl)

          endif

          call mpp_error(FATAL,'==>Error in ocean_adv_vel_diag_mod: max_cfl_value exceeded.')

        endif

      enddo
    enddo
  enddo

9100 format(1x,a12,1x,'ts=',i10,1x,',k=',i3,', z(m)=',f8.1,', lon(deg):',f6.2,' --> ',f6.2,&
            ', lat(deg):',f6.2,' --> ',f6.2,', scaling=',1pg10.3)

end subroutine cfl_check2
! </SUBROUTINE> NAME="cfl_check2"


!#######################################################################
! <SUBROUTINE NAME="maximum_bottom_w">
!
! <DESCRIPTION>
! Compute maximum vertical velocity on the bottom of tracer and velocity cells.
! The vertical velocity at bottom of a column of tracer cells should be roundoff.  
! For flat bottom simulations, the vertical velocity on the bottom of the 
! velocity cell column should also be roundoff.  For simulations with topography,
! the vertical velocity on the bottom of a velocity cell column will not vanish
! due to the effects of topography.  
! </DESCRIPTION>
!
subroutine maximum_bottom_w(Adv_vel)

  type(ocean_adv_vel_type), intent(in) :: Adv_vel

  real    :: wtbot, wubot, wtbot0, wubot0, fudge
  integer :: i, j, k, iwtbot, jwtbot, kwtbot, iwubot, jwubot, kwubot

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_adv_vel_diag_mod (maximum_bottom_w): module needs initialization ')
  endif 

  ! Find Max error in "w_bt" at bottom

  wtbot=0.0; iwtbot=isc; jwtbot=jsc; kwtbot=1
  do j=jsc,jec
    do i=isc,iec
      k = Grd%kmt(i,j)
      if (k /= 0 .and. (abs(Adv_vel%wrho_bt(i,j,k)) > abs(wtbot))) then
        wtbot  = Adv_vel%wrho_bt(i,j,k)
        iwtbot = i
        jwtbot = j
        kwtbot = k
      endif
    enddo
  enddo
  wtbot = rho0r*abs(wtbot)

  ! Find Max "wrho_bu" at bottom
  ! Relevant only for Bgrid.
  ! Note that wrho_bu is nonzero at the bottom of a U-column
  ! since it is part of the slope velocity.

  wubot=0.0; iwubot=isc; jwubot=jsc; kwubot=1
  do j=jsc,jec
    do i=isc,iec
      k = Grd%kmu(i,j)
      if (k /= 0 .and. (abs(Adv_vel%wrho_bu(i,j,k)) > abs(wubot))) then 
        wubot  = Adv_vel%wrho_bu(i,j,k)
        iwubot = i
        jwubot = j
        kwubot = k
      endif
    enddo
  enddo
  wubot = rho0r*abs(wubot)

  write (stdoutunit,'(//60x,a/)') ' Bottom vertical velocity summary:'
  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  wtbot  = wtbot*fudge
  wubot  = wubot*fudge
  wtbot0 = wtbot
  wubot0 = wubot
  wtbot  = abs(wtbot)
  wubot  = abs(wubot)
  call mpp_max(wtbot)
  call mpp_max(wubot)

  if (abs(wtbot0) == wtbot) then
    wtbot = wtbot0/fudge
    write (unit,9112) wtbot, iwtbot+Dom%ioff, jwtbot+Dom%joff, kwtbot, &
                      Grd%xt(iwtbot,jwtbot), Grd%yt(iwtbot,jwtbot), Grd%zt(kwtbot)
  endif
  
  if(horz_grid == MOM_BGRID) then
     if (abs(wubot0) == wubot) then
       wubot = wubot0/fudge
       write (unit,9113) wubot, iwubot+Dom%ioff, jwubot+Dom%joff, kwubot, &
                         Grd%xu(iwubot,jwubot), Grd%yu(iwubot,jwubot), Grd%zt(kwubot)
     endif
  endif

9112  format(/' Maximum T-cell bottom velocity (',es10.3,' m/s){error}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)')
9113  format(' Maximum U-cell bottom velocity (',es10.3,' m/s){slope}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)'/)

end subroutine maximum_bottom_w
! </SUBROUTINE> NAME="maximum_bottom_w"


!#######################################################################
! <SUBROUTINE NAME="max_continuity_error">
!
! <DESCRIPTION>
! Compute continuity error. Should be roundoff if all is working well.  
! </DESCRIPTION>
!
subroutine max_continuity_error(Adv_vel, Thickness)

  type(ocean_adv_vel_type),   intent(in) :: Adv_vel
  type(ocean_thickness_type), intent(in) :: Thickness
  
  real    :: bigt, bigu, bigt0, bigu0, fudge
  integer :: i, j, k, iu, ju, ku, it, jt, kt

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_adv_vel_diag_mod (max_continuity_error): module needs initialization ')
  endif 

  do k=1,nk
       
     tmp(:,:)    = Thickness%rho_dzt_tendency(:,:,k) - Thickness%mass_source(:,:,k)
     wrk1(:,:,k) =  (tmp(:,:)                        &
                   + BDX_ET(Adv_vel%uhrho_et(:,:,k)) &
                   + BDY_NT(Adv_vel%vhrho_nt(:,:,k)) &
                   + Adv_vel%wrho_bt(:,:,k-1)-Adv_vel%wrho_bt(:,:,k))*Grd%tmask(:,:,k)

     ! relevant only for Bgrid
     tmp(:,:)    = Grd%umask(:,:,k)*REMAP_BT_TO_BU(tmp(:,:))
     wrk2(:,:,k) =  (tmp(:,:)                        &
                   + BDX_EU(Adv_vel%uhrho_eu(:,:,k)) &
                   + BDY_NU(Adv_vel%vhrho_nu(:,:,k)) &
                   + Adv_vel%wrho_bu(:,:,k-1)-Adv_vel%wrho_bu(:,:,k))*Grd%umask(:,:,k)

  enddo

  bigt = epsln !find location of largest continuity error
  bigu = epsln 
  it = isc; jt = jsc; kt = 1
  iu = isc; ju = jsc; ku = 1
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (abs(wrk1(i,j,k)) > abs(bigt)) then
          bigt = wrk1(i,j,k)
          it = i
          jt = j
          kt = k
        endif 
        if (abs(wrk2(i,j,k)) > abs(bigu)) then
          bigu = wrk2(i,j,k)
          iu = i
          ju = j
          ku = k
        endif 
      enddo
    enddo
  enddo
  
  ! to distinguish processors when tracer is independent of processor
  fudge = 1 + 1.e-12*mpp_pe() 

  bigt  = bigt*fudge*rho0r
  bigu  = bigu*fudge*rho0r
  bigt0 = bigt
  bigu0 = bigu
  bigt  = abs(bigt)
  bigu  = abs(bigu)
  call mpp_max(bigt)
  call mpp_max(bigu)

  write (stdoutunit,'(/60x,a/)') ' Continuity error summary:'

  if (abs(bigt0) == bigt) then
    bigt = bigt0/fudge
    write (unit,'(a,es10.3,a,i4,a1,i3,a1,i3,a,f9.2,a,f9.2,a,f9.2,a/)') ' Maximum T-cell Continuity Error (',bigt,&
    ' m/s) is at (i,j,k) = (',it+Dom%ioff,',',jt+Dom%joff,',',kt,'),  (lon,lat,dpt) = ('&
     ,Grd%xt(it,jt),',',Grd%yt(it,jt),',',  Grd%zt(kt),' m)'
  endif

  if(horz_grid==MOM_BGRID) then
     if (abs(bigu0) == bigu) then
       bigu = bigu0/fudge
       write (unit,'(/,a,es10.3,a,i4,a1,i3,a1,i3,a,f9.2,a,f9.2,a,f9.2,a)') ' Maximum U-cell Continuity Error (',bigu,&
       ' m/s) is at (i,j,k) = (',iu+Dom%ioff,',',ju+Dom%joff,',',ku,'),  (lon,lat,dpt) = ('&
        ,Grd%xu(iu,ju),',',Grd%yu(iu,ju),',',  Grd%zt(ku),' m)'
     endif
  endif

  
end subroutine max_continuity_error
! </SUBROUTINE> NAME="max_continuity_error"


!#######################################################################
! <SUBROUTINE NAME="transport_on_s">
!
! <DESCRIPTION>
! Compute transports on s-levels (defined by same k-level)
! and send to diag_manager.
! </DESCRIPTION>
!
subroutine transport_on_s(Time, Adv_vel)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  integer :: i, j, k
  
  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_adv_vel_diag_mod (transport_on_s): module needs initialization ')
  endif 

  ! mass transports leaving faces of grid cells
  wrk1=0.0

  if (id_tx_trans > 0) then 
    do k=1,nk
         do j=jsc,jec 
            do i=isc,iec
               wrk1(i,j,k) = Adv_vel%uhrho_et(i,j,k)*Grd%dyte(i,j)*transport_convert
            enddo
         enddo
    enddo 
    call diagnose_3d(Time, Grd, id_tx_trans, wrk1(:,:,:))
  endif 
  if (id_tx_trans_int_z > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec 
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + Adv_vel%uhrho_et(i,j,k)*Grd%dyte(i,j)*transport_convert
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_tx_trans_int_z, wrk1_2d(:,:))
  endif
 
  if (id_ty_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%vhrho_nt(i,j,k)*Grd%dxtn(i,j)*transport_convert 
         enddo 
      enddo
    enddo 
    call diagnose_3d(Time, Grd, id_ty_trans, wrk1(:,:,:))
  endif 
  if (id_ty_trans_int_z > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec 
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + Adv_vel%vhrho_nt(i,j,k)*Grd%dxtn(i,j)*transport_convert 
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_ty_trans_int_z, wrk1_2d(:,:))
  endif

  if (id_tz_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%wrho_bt(i,j,k)*Grd%dat(i,j)*transport_convert 
         enddo 
      enddo
    enddo 
    call diagnose_3d(Time, Grd, id_tz_trans, wrk1(:,:,:))
  endif 
  if (id_tz_trans_int_z > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec 
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + Adv_vel%wrho_bt(i,j,k)*Grd%dat(i,j)*transport_convert
            enddo
         enddo
      enddo
      call diagnose_2d(TIme, Grd, id_tz_trans_int_z, wrk1_2d(:,:))
  endif
  if (id_tz_trans_sq > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = (Adv_vel%wrho_bt(i,j,k)*Grd%dat(i,j)*transport_convert)**2 
         enddo 
      enddo
    enddo 
    call diagnose_3d(Time, Grd, id_tz_trans_sq, wrk1(:,:,:))
  endif 

  ! relevant only for Bgrid 
  if (id_ux_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%uhrho_eu(i,j,k)*Grd%dyue(i,j)*transport_convert
         enddo 
      enddo
    enddo 
    call diagnose_3d_u(Time, Grd, id_ux_trans, wrk1(:,:,:))
  endif 

  ! relevant only for Bgrid 
  if (id_uy_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%vhrho_nu(i,j,k)*Grd%dxun(i,j)*transport_convert 
         enddo 
      enddo
    enddo 
    call diagnose_3d_u(Time, Grd, id_uy_trans, wrk1(:,:,:))
  endif 

  ! relevant only for Bgrid 
  if (id_uz_trans > 0) then 
    do k=1,nk
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,k) = Adv_vel%wrho_bu(i,j,k)*Grd%dxu(i,j)*Grd%dyu(i,j)*transport_convert
         enddo 
      enddo 
    enddo 
    call diagnose_3d_u(Time, Grd, id_uz_trans, wrk1(:,:,:))
  endif 

end subroutine transport_on_s
! </SUBROUTINE> NAME="transport_on_s"


!#######################################################################
! <SUBROUTINE NAME="transport_on_nrho">
!
! <DESCRIPTION>
! Classify horizontal transport according to neutral density classes. 
!
! NOTE: This diagnostic works with transport through cell faces.
! To get transport_on_nrho, a binning must be done, rather than 
! a remapping (as done for trans_rho_gm).  
!
! code history:
!
! 2002: Stephen.Griffies
! Zhi.Liang
! Alexander.Pletzer
!
! 2007: updated algorithm with weighting as done in the paper
! Lee, Nurser, Coward, and de Cuevas, 2007:
! "Eddy advective and diffusive transports of heat and salt 
! in the Southern Ocean" JPO, vol 37, pages 1376-1393
!
! 2010: updated by Rusty.Benson for optimization.
!
! 2010: Corrected weighting by Stephen.Griffies
!
! 2010: Removed the weighting in favor of a strait rebinning
!       approach.  The weighting approach was 
!       unnecessary, and added more cost to the scheme.  
! 
! </DESCRIPTION>
!
subroutine transport_on_nrho (Time, Dens, Adv_vel)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  type(time_type)                      :: next_time 

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: work1(isc:iec,jsc:jec,size(Dens%neutralrho_ref))
  real    :: work2(isc:iec,jsc:jec,size(Dens%neutralrho_ref))
  real    :: tmp(2,isc:iec,jsc:jec)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_adv_vel_diag_mod (transport_on_nrho): module needs initialization ')
  endif 

  neutralrho_nk = size(Dens%neutralrho_ref(:))
  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_tx_trans_nrho,next_time) .or. need_data(id_ty_trans_nrho,next_time)) then

      work1(:,:,:) = 0.0
      work2(:,:,:) = 0.0

      do k_rho=1,neutralrho_nk

         tmp(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if (k_rho == 1) then
                     if(Dens%neutralrho(i,j,k) < Dens%neutralrho_bounds(k_rho)) then 
                         tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                         tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  elseif(k_rho < neutralrho_nk) then
                     if( (Dens%neutralrho_bounds(k_rho) <= Dens%neutralrho(i,j,k)) .and.  &
                         (Dens%neutralrho(i,j,k)        <  Dens%neutralrho_bounds(k_rho+1)) ) then 
                           tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                           tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  else    ! if (k_rho == neutralrho_nk) then
                     if(Dens%neutralrho_bounds(k_rho) <= Dens%neutralrho(i,j,k)) then 
                         tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)             
                         tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  endif
               enddo
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               work1(i,j,k_rho) = (tmp(1,i,j)+work1(i,j,k_rho))*Grd%dyte(i,j)*transport_convert*Grd%tmask(i,j,1)
               work2(i,j,k_rho) = (tmp(2,i,j)+work2(i,j,k_rho))*Grd%dxtn(i,j)*transport_convert*Grd%tmask(i,j,1)
            enddo
         enddo

      enddo

      if (id_tx_trans_nrho > 0) then 
          used = send_data (id_tx_trans_nrho, work1, Time%model_time, &
                            ks_in=1, ke_in=neutralrho_nk)
      endif
      if (id_ty_trans_nrho > 0) then 
          used = send_data (id_ty_trans_nrho, work2, Time%model_time, &
                            ks_in=1, ke_in=neutralrho_nk)
      endif

  endif

end subroutine transport_on_nrho
! </SUBROUTINE> NAME="transport_on_nrho"


!#######################################################################
! <SUBROUTINE NAME="transport_on_rho">
!
! <DESCRIPTION>
! Classify horizontal transport according to potential density classes. 
!
! Diagnostic makes sense when potrho is monotonically increasing with 
! depth, although the algorithm does not explicitly make this assumption.  
!
! NOTE: This diagnostic works with transport through cell faces.
! To get transport_on_rho, a binning must be done, rather than 
! a remapping (as done for trans_rho_gm).  
!
! Code history 
!
! 2002: Stephen.Griffies
! Zhi.Liang
! Alexander.Pletzer
!
! 2007: updated algorithm with weighting as done in the paper
! Lee, Nurser, Coward, and de Cuevas, 2007:
! "Eddy advective and diffusive transports of heat and salt 
! in the Southern Ocean" JPO, vol 37, pages 1376-1393
!
! 2010: updated by Rusty.Benson for optimization.
!       corrected weigthing by Stephen.Griffies
!
! 2010: Removed the weighting in favor of a strait rebinning
!       approach.  The weighting approach was 
!       unnecessary, and added more cost to the scheme.  
!
! </DESCRIPTION>
!
subroutine transport_on_rho (Time, Dens, Adv_vel)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  type(ocean_adv_vel_type), intent(in) :: Adv_vel
  type(time_type)                      :: next_time 

  integer :: i, j, k, k_rho, potrho_nk
  real    :: work1(isc:iec,jsc:jec,size(Dens%potrho_ref))
  real    :: work2(isc:iec,jsc:jec,size(Dens%potrho_ref))
  real    :: tmp(2,isc:iec,jsc:jec)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_adv_vel_diag_mod (transport_on_rho): module needs initialization ')
  endif 

  potrho_nk = size(Dens%potrho_ref(:))
  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_tx_trans_rho,next_time) .or. need_data(id_ty_trans_rho,next_time)) then

      work1(:,:,:) = 0.0
      work2(:,:,:) = 0.0

      do k_rho=1,potrho_nk
         tmp(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if (k_rho == 1) then 
                     if(Dens%potrho(i,j,k) < Dens%potrho_bounds(k_rho)) then 
                        tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                        tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  elseif(k_rho < potrho_nk) then 
                     if( (Dens%potrho_bounds(k_rho) <= Dens%potrho(i,j,k)) .and. &
                         (Dens%potrho(i,j,k)        <  Dens%potrho_bounds(k_rho+1)) ) then 
                           tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                           tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  else  ! if (k_rho == potrho_nk) then
                     if(Dens%potrho_bounds(k_rho) <= Dens%potrho(i,j,k)) then 
                        tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)             
                        tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  endif
               enddo
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               work1(i,j,k_rho) = (tmp(1,i,j)+work1(i,j,k_rho))*Grd%dyte(i,j)*transport_convert*Grd%tmask(i,j,1)
               work2(i,j,k_rho) = (tmp(2,i,j)+work2(i,j,k_rho))*Grd%dxtn(i,j)*transport_convert*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_rho > 0) then 
          used = send_data (id_tx_trans_rho, work1, Time%model_time, &
                            ks_in=1, ke_in=potrho_nk)
      endif
      if (id_ty_trans_rho > 0) then 
          used = send_data (id_ty_trans_rho, work2, Time%model_time, &
                            ks_in=1, ke_in=potrho_nk)
      endif

  endif

end subroutine transport_on_rho
! </SUBROUTINE> NAME="transport_on_rho"


!#######################################################################
! <SUBROUTINE NAME="transport_on_theta">
!
! <DESCRIPTION>
! Classify horizontal transport of mass according to potential 
! temperature classes.  This diagnostic is useful to deduce the 
! heat that is transported between potential temperature classes. 
!
! Diagnostic makes sense when theta is monotonically decreasing
! with depth, although the algorithm does not explicitly make
! this assumption.  
!
! NOTE: This diagnostic works with transport through cell faces.
! To get transport_on_theta, a binning must be done, rather than 
! a remapping (as done for trans_rho_gm).  
!
! Code history 
!
! 2002: Stephen.Griffies
! Zhi.Liang
! Alexander.Pletzer
!
! 2007: updated algorithm with weighting as done in the paper
! Lee, Nurser, Coward, and de Cuevas, 2007:
! "Eddy advective and diffusive transports of heat and salt 
! in the Southern Ocean" JPO, vol 37, pages 1376-1393
!
! 2010: updated by Rusty.Benson for optimization.
!       corrected weigthing by Stephen.Griffies
!
! 2010: Removed the weighting in favor of a strait rebinning
!       approach.  The weighting approach was 
!       unnecessary, and added more cost to the scheme.  
!
! </DESCRIPTION>
!
subroutine transport_on_theta (Time, Dens, Theta, Adv_vel)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: Theta
  type(ocean_adv_vel_type),     intent(in) :: Adv_vel

  type(time_type)                          :: next_time 

  integer :: i, j, k, k_theta, theta_nk, tau
  real    :: work1(isc:iec,jsc:jec,size(Dens%theta_ref))
  real    :: work2(isc:iec,jsc:jec,size(Dens%theta_ref))
  real    :: tmp(2,isc:iec,jsc:jec)

  tau = Time%tau

  theta_nk = size(Dens%theta_ref(:))
  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_tx_trans_theta,next_time) .or. need_data(id_ty_trans_theta,next_time)) then

      work1(:,:,:)  = 0.0
      work2(:,:,:)  = 0.0

      do k_theta=1,theta_nk
         tmp(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if (k_theta == 1) then
                     if(Theta%field(i,j,k,tau) < Dens%theta_bounds(k_theta)) then 
                         tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                         tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  elseif(k_theta < theta_nk) then
                     if( (Dens%theta_bounds(k_theta) <= Theta%field(i,j,k,tau)) .and. &
                         (Theta%field(i,j,k,tau)     <  Dens%theta_bounds(k_theta+1)) ) then 
                          tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                          tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  else ! if (k_theta == theta_nk) then
                     if(Dens%theta_bounds(k_theta) <= Theta%field(i,j,k,tau)) then 
                         tmp(1,i,j) = tmp(1,i,j) + Adv_vel%uhrho_et(i,j,k)
                         tmp(2,i,j) = tmp(2,i,j) + Adv_vel%vhrho_nt(i,j,k)
                     endif
                  endif
               enddo
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               work1(i,j,k_theta) = (tmp(1,i,j)+work1(i,j,k_theta))*Grd%dyte(i,j)*transport_convert*Grd%tmask(i,j,1)
               work2(i,j,k_theta) = (tmp(2,i,j)+work2(i,j,k_theta))*Grd%dxtn(i,j)*transport_convert*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_theta > 0) then 
          used = send_data (id_tx_trans_theta, work1, Time%model_time, &
                            ks_in=1, ke_in=theta_nk)
      endif
      if (id_ty_trans_theta > 0) then 
          used = send_data (id_ty_trans_theta, work2, Time%model_time, &
                            ks_in=1, ke_in=theta_nk)
      endif

  endif

end subroutine transport_on_theta
! </SUBROUTINE> NAME="transport_on_theta"


!#######################################################################
! <SUBROUTINE NAME="vertical_reynolds_check">
!
! <DESCRIPTION>
! This subroutine computes the Reynolds number associated with vertical
! friction visc_cbt and vertical velocity wrho_bt.  This check is
! appropriate for either Bgrid or Cgrid. 
! </DESCRIPTION>
!
subroutine vertical_reynolds_check(Adv_vel, Thickness, visc_cbt)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: visc_cbt

  real    :: ramb, reyz
  real    :: fudge
  real    :: reynz0, reynz, reynw, reynmw
  integer :: ireynz, jreynz, kreynz
  integer :: i, j, k, kk

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynz=isc; jreynz=jsc; kreynz=1
  reynz=0.0 ; reynw=0.0 ; reynmw=0.0

  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        ramb = 1.0/(visc_cbt(i,j,k) + epsln)

        kk = min(k+1,nk)
        if (k >= Grd%kmt(i,j)) then
          reyz = 0.0
        else 
          reyz = Grd%tmask(i,j,kk)*abs(Adv_vel%wrho_bt(i,j,k)*Thickness%dzwt(i,j,k))*ramb
        endif
        if (reyz > reynz) then
          ireynz = i
          jreynz = j
          kreynz = k
          reynz  = reyz
          reynw  = Adv_vel%wrho_bt(i,j,k)
          reynmw = 1.0/ramb
        endif
      enddo
    enddo
  enddo
  write (stdoutunit,'(/60x,a/)') ' Vertical Laplacian Grid Reynolds number summary'


  fudge = 1 + 1.e-12*mpp_pe()  ! used to distinguish processors when result is independent of processor
  reynz = reynz*fudge 
  if(reynz == 0.0) then 
    reynz = reynz + 1.e-12*mpp_pe() 
  endif 

  reynz0 = reynz
  call mpp_max(reynz)
  
  if (reynz == reynz0) then
     write(unit,'(/,a,es9.2,a,i4,a1,i4,a1,i3,a,f7.2,a,f7.2,a,f8.2,a,/,1x,a,es9.2,a,es9.2,/)') &
           ' Maximum vertical grid Re # = ', reynz,&
           ' at (i,j,k) = (',ireynz+Dom%ioff,',',jreynz+Dom%joff,',',kreynz,&
           '),  (lon,lat,dpt) = (',Grd%xt(ireynz,jreynz),',', &
           Grd%yt(ireynz,jreynz),',', Grd%dzt(kreynz),' m). ', 'Vertical velocity "w_bt" is (m/s) ',&
            rho0r*reynw, ' and vertical viscosity (m^2/s) "visc_cbt" is ', reynmw
  endif

end subroutine vertical_reynolds_check
! </SUBROUTINE> NAME="vertical_reynolds_check"


end module ocean_adv_vel_diag_mod



module ocean_velocity_advect_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Velocity advective transport. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes advection of velocity using one of the 
! following advection schemes:
! 1/ second order centered 
! 2/ first order upwind 
!</DESCRIPTION>
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
! A Technical Guide to MOM4 (2004)
! </REFERENCE>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! </REFERENCE>
!
!<REFERENCE>
! Hundsdorder and Trompert (1994), "Method of lines and 
! direct discretization: a comparison for linear 
! advection", Applied Numerical Mathematics,
! pages 469--490.
!</REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_velocity_advect_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging
!  </DATA> 
!  <DATA NAME="zero_velocity_advect_horz" TYPE="logical">
!  For debugging, it is often useful to remove horizontal advection of velocity. 
!  </DATA> 
!  <DATA NAME="zero_velocity_advect_vert" TYPE="logical">
!  For debugging, it is often useful to remove vertical advection of velocity. 
!  </DATA> 
!
!  <DATA NAME="velocity_advect_centered" TYPE="logical">
!  For using the standard second order centered method for 
!  computing the advection of linear momentum. This is the 
!  default: velocity_advect_centered=.true.
!  </DATA> 
!  <DATA NAME="velocity_advect_upwind" TYPE="logical">
!  For using the first order upwind method for 
!  computing the advection of linear momentum. 
!  Default: velocity_advect_upwind=.false.
!  </DATA> 
!
!</NAMELIST>

use diag_manager_mod, only: register_diag_field
use fms_mod,          only: mpp_error, FATAL, NOTE, stdout, stdlog, write_version_number
use fms_mod,          only: read_data, open_namelist_file, check_nml_error, close_file
use mpp_mod,          only: input_nml_file

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_update_boundary, ocean_obc_zero_boundary
use ocean_operators_mod,  only: FAX, FAY, BDX_EU, BDY_NU, REMAP_BT_TO_BU
use ocean_parameters_mod, only: missing_value, onefourth
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID
use ocean_parameters_mod, only: VEL_ADVECT_UPWIND, VEL_ADVECT_2ND_ORDER
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_adv_vel_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type
use ocean_workspace_mod,  only: wrk1_v, wrk1_2d, wrk2_2d   
use ocean_util_mod,       only: write_timestamp, diagnose_2d_u, diagnose_3d_u, diagnose_3d_en, write_chksum_3d

implicit none

private

! for Bgrid or Cgrid
integer :: horz_grid

public horz_advection_of_velocity
public vert_advection_of_velocity
public ocean_velocity_advect_init

private horz_advection_centered
private vert_advection_centered
private horz_advection_upwind
private vert_advection_upwind 


! for diagnostics 
logical :: used 
integer :: id_hadv_u=-1
integer :: id_hadv_v=-1
integer :: id_vadv_u=-1
integer :: id_vadv_v=-1
integer :: id_surf_accel(2)=-1
integer :: id_pme_u  =-1
integer :: id_river_u=-1

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS
real, dimension(isd:ied,jsd:jed)  :: tmp, tmp1, tmp2
#else
real, dimension(:,:), allocatable :: tmp, tmp1, tmp2
#endif

real, dimension(:),   allocatable :: kmask

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

character(len=128) :: version=&
  '$Id: ocean_velocity_advect.F90,v 20.0 2013/12/14 00:12:43 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: advection_scheme           = 2
logical :: module_is_initialized      = .false.
logical :: debug_this_module          = .false.
logical :: zero_velocity_advect_horz  = .false. 
logical :: zero_velocity_advect_vert  = .false. 
logical :: have_obc                   = .false.
logical :: velocity_advect_centered   = .true.
logical :: velocity_advect_upwind     = .false.

namelist /ocean_velocity_advect_nml/ debug_this_module,         &
          zero_velocity_advect_horz, zero_velocity_advect_vert, &
          velocity_advect_centered, velocity_advect_upwind 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_advect_init">
!
! <DESCRIPTION>
! Initialize the velocity advection module.
! </DESCRIPTION>
!
subroutine ocean_velocity_advect_init(Grid, Domain, Time, obc, hor_grid, debug)

  type(ocean_grid_type),   intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target   :: Domain
  type(ocean_time_type),   intent(in), target   :: Time
  logical,                 intent(in)           :: obc
  integer,                 intent(in)           :: hor_grid
  logical,                 intent(in), optional :: debug

  integer :: ioun, io_status, ierr
  integer :: num_methods=0 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_velocity_advect_mod (ocean_velocity_advect_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  Grd => Grid
  Dom => Domain
  horz_grid = hor_grid
  have_obc  = obc

  ! provide for namelist over-ride
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_velocity_advect_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_velocity_advect_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_velocity_advect_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_velocity_advect_nml')
  call close_file (ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_velocity_advect_nml)  
  write (stdlogunit, ocean_velocity_advect_nml)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_velocity_advect with debug_this_module=.true.'  
  endif 
  
#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grid%nk
  allocate(tmp(isd:ied,jsd:jed),tmp1(isd:ied,jsd:jed),tmp2(isd:ied,jsd:jed))
#endif
  tmp =0.0
  tmp1=0.0
  tmp2=0.0

  ! for help with masking flux at ocean bottom 
  allocate(kmask(nk))
  kmask(1:nk-1) = 1.0
  kmask(nk)     = 0.0

  if(zero_velocity_advect_horz) then 
    write(stdoutunit,*) &
    '==>Warning: running MOM with zero horizontal advection of velocity. Unrealistic simulation.'
  endif 
  if(zero_velocity_advect_vert) then 
    write(stdoutunit,*) &
    '==>Warning: running MOM with zero vertical advection of velocity. Unrealistic simulation.'
  endif 

  if(velocity_advect_centered) then 
    write(stdoutunit,*) &
    '==>Note: running MOM with traditional second order centred advection of linear momentum. '
    num_methods = num_methods+1
    advection_scheme = VEL_ADVECT_2ND_ORDER
  endif 

  if(velocity_advect_upwind) then 
    if(horz_grid == MOM_CGRID) then 
       write(stdoutunit,'(/a)') &
       '==>Error from ocean_velocity_advect_mod: Cgrid has yet to code upwind velocity advection.'
       call mpp_error(FATAL,&
       '==>Error from ocean_velocity_advect_mod: Cgrid has yet to code upwind velocity advection.')
    endif 
    write(stdoutunit,'(/a)') &
    '==>Note: running MOM with first order upwind advection of linear momentum. '
    write(stdoutunit,'(a)') &
    '   This method is for testing purposes.  It is not generally recommened.'
    write(stdoutunit,'(a)') &
    '   Diagnosed energy conversion errors will be large, as energy analysis assumes 2nd order centred advection.'
    num_methods = num_methods+1
    advection_scheme = VEL_ADVECT_UPWIND
  endif 

  if(num_methods==0) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_velocity_advect_mod: No momentum advection scheme chosen. Must choose a scheme.')
  endif 
  if(num_methods>1) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_velocity_advect_mod: More than one momentum advect scheme chosen. Choose only one.')
  endif 


  ! register for diagnostics manager 

  id_hadv_u = register_diag_field ('ocean_model', 'hadv_u', Grd%vel_axes_u(1:3), &
              Time%model_time, 'Thickness and rho wghtd horz advection of u',    &
              '(kg/m^3)*(m^2/s^2)',missing_value=-1e10, range=(/-1e10,1e10/))
  id_hadv_v = register_diag_field ('ocean_model', 'hadv_v', Grd%vel_axes_v(1:3), &
              Time%model_time, 'Thickness and rho wghtd horz advection of v',    &
              '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))

  id_vadv_u = register_diag_field ('ocean_model', 'vadv_u', Grd%vel_axes_u(1:3), &
              Time%model_time, 'Thickness and rho  wghtd vert advection of u',   &
              '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_vadv_v = register_diag_field ('ocean_model', 'vadv_v', Grd%vel_axes_v(1:3), &
              Time%model_time, 'Thickness and rho  wghtd vert advection of v',   &
              '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))

  id_surf_accel(1) = register_diag_field ('ocean_model', 'surf_accel_u', Grd%vel_axes_u(1:2), &
                     Time%model_time, 'uh-forcing from fresh water',                          &
                     '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_surf_accel(2) = register_diag_field ('ocean_model', 'surf_accel_v', Grd%vel_axes_v(1:2), &
                     Time%model_time, 'vh-forcing from fresh water',                          &
                     '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e10,1e10/))

  id_pme_u         = register_diag_field ('ocean_model', 'pme_u', Grd%vel_axes_uv(1:2), &
                     Time%model_time, 'pme on u-cell', '(kg/m^3)*(m/s)',                &
                     missing_value=missing_value, range=(/-1e10,1e10/))
  id_river_u       = register_diag_field ('ocean_model', 'river_u', Grd%vel_axes_uv(1:2), &
                     Time%model_time, 'river on u-cell', '(kg/m^3)*(m/s)',                &
                     missing_value=missing_value, range=(/-1e10,1e10/))


end subroutine ocean_velocity_advect_init
! </SUBROUTINE> NAME="ocean_velocity_advect_init"


!#######################################################################
! <SUBROUTINE NAME="horz_advection_of_velocity">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to horizontal (constant k-level)
! advection of velocity.
!
! </DESCRIPTION>
!
subroutine horz_advection_of_velocity(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
 
  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL,&
      '==>Error from ocean_velocity_advect_mod (horz_advection_of_velocity): module not yet initialized')
  endif 

  if(zero_velocity_advect_horz) then 
      Velocity%wrkv(:,:,:,:) = 0.0
      return 
  endif

  if(advection_scheme==VEL_ADVECT_2ND_ORDER) then 
    call horz_advection_centered(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)  
  elseif(advection_scheme==VEL_ADVECT_UPWIND) then
    call horz_advection_upwind(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)  
  endif 

end subroutine horz_advection_of_velocity
! </SUBROUTINE> NAME="horz_advection_of_velocity"


!#######################################################################
! <SUBROUTINE NAME="horz_advection_centered">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to horizontal (constant k-level)
! advection of velocity.
!
! Use second order centered method for either Bgrid or Cgrid. 
!
! </DESCRIPTION>
!
subroutine horz_advection_centered(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step

  real    :: metric_force
  integer :: i, j, k, n
  integer :: tau, tau_m0
  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau    = Time%tau
  tau_m0 = Time%tau_m0
  wrk1_v = 0.0

  if(horz_grid == MOM_BGRID) then 

     do k=1,nk
        do n=1,2

           tmp1(:,:) = Adv_vel%uhrho_eu(:,:,k)*FAX(Velocity%u(:,:,k,n,tau))
           tmp2(:,:) = Adv_vel%vhrho_nu(:,:,k)*FAY(Velocity%u(:,:,k,n,tau))

           ! compute horizontal divergence of flux 
           tmp(:,:)  = BDX_EU(tmp1) + BDY_NU(tmp2) 

           ! remove advection of momentum but leave the metric term unchanged
           if(have_obc) then
             call ocean_obc_zero_boundary(tmp(:,:), 'C')
           endif
           do j=jsc,jec
             do i=isc,iec
               tmp(i,j) = tmp(i,j) &
                           + (3-2*n)*Thickness%rho_dzu(i,j,k,tau)*Velocity%u(i,j,k,3-n,tau) &
                            *(Grd%dh1dy(i,j)*Velocity%u(i,j,k,1,tau)-Grd%dh2dx(i,j)*Velocity%u(i,j,k,2,tau))
               wrk1_v(i,j,k,n) = Grd%umask(i,j,k)*tmp(i,j)
             enddo
           enddo
        enddo
     enddo

  else ! cgrid 

     do k=1,nk

        ! minus the advection metric frequency
        wrk1_2d(:,:) = Grd%dh1dy(:,:)*FAY(Velocity%u(:,:,k,1,tau)) - Grd%dh2dx(:,:)*FAX(Velocity%u(:,:,k,2,tau))

        ! zonal velocity 
        n=1
        tmp1(:,:) = FAX(Adv_vel%uhrho_et(:,:,k))*FAX(Velocity%u(:,:,k,n,tau))
        tmp2(:,:) = FAX(Adv_vel%vhrho_nt(:,:,k))*FAY(Velocity%u(:,:,k,n,tau))
          
        do j=jsc,jec
           do i=isc,iec
              metric_force = onefourth*( (Adv_vel%vhrho_nt(i,j,k)   + Adv_vel%vhrho_nt(i+1,j,k))  *wrk1_2d(i,j) &
                                        +(Adv_vel%vhrho_nt(i,j-1,k) + Adv_vel%vhrho_nt(i+1,j-1,k))*wrk1_2d(i,j-1) ) 
              tmp(i,j) = (tmp1(i,j)-tmp1(i-1,j))*Grd%dxter(i,j) + (tmp2(i,j)-tmp2(i,j-1))*Grd%dyter(i,j) &
                         + metric_force 
              wrk1_v(i,j,k,n) = Grd%tmasken(i,j,k,n)*tmp(i,j)
           enddo
        enddo  

        ! meridional velocity 
        n=2
        tmp1(:,:) = FAY(Adv_vel%uhrho_et(:,:,k))*FAX(Velocity%u(:,:,k,n,tau))
        tmp2(:,:) = FAY(Adv_vel%vhrho_nt(:,:,k))*FAY(Velocity%u(:,:,k,n,tau))
          
        do j=jsc,jec
           do i=isc,iec
              metric_force = -onefourth*((Adv_vel%uhrho_et(i-1,j,k) + Adv_vel%uhrho_et(i-1,j+1,k))*wrk1_2d(i-1,j) &
                                        +(Adv_vel%uhrho_et(i,j,k)   + Adv_vel%uhrho_et(i,j+1,k))  *wrk1_2d(i,j-1) )
              tmp(i,j) =  (tmp1(i,j)-tmp1(i-1,j))*Grd%dxtnr(i,j) + (tmp2(i,j)-tmp2(i,j-1))*Grd%dytnr(i,j) &
                          + metric_force 
              wrk1_v(i,j,k,n) = Grd%tmasken(i,j,k,n)*tmp(i,j)
           enddo
        enddo  

     enddo   ! k-loop


  endif  ! for horz_grid 


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
                  Velocity%advection(i,j,k,n,tau_m0) = Velocity%advection(i,j,k,n,tau_m0) + wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_en(Time, Grd, id_hadv_u, id_hadv_v, wrk1_v(:,:,:,:))

      if(debug_this_module) then 
          write(stdoutunit,*) ' '
          write(stdoutunit,*) 'From ocean_velocity_advect_mod: horz_advection_centered chksums'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('horz_advection_of_velocity(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('horz_advection_of_velocity(2)', wrk1_v(COMP,:,2))
      endif

  endif

end subroutine horz_advection_centered
! </SUBROUTINE> NAME="horz_advection_centered"


!#######################################################################
! <SUBROUTINE NAME="horz_advection_upwind">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to horizontal (constant k-level)
! advection of velocity.
!
! Use first order upwind method. Not coded for Cgrid, since
! first order upwind method is rarely used.  
!
! </DESCRIPTION>
!
subroutine horz_advection_upwind(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step
 
  real, dimension(isd:ied,jsd:jed) :: fe, fn 
  real                             :: vel, upos, uneg 
  integer                          :: i, j, k, n
  integer                          :: tau, taum1, tau_m0

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau    = Time%tau
  taum1  = Time%taum1
  tau_m0 = Time%tau_m0
  wrk1_v = 0.0

  do k=1,nk
     do n=1,2

        ! i-flux 
        fe = 0.0
        do j=jsd,jed
           do i=isd,iec
              vel     = 0.5*Adv_vel%uhrho_eu(i,j,k)
              upos    = vel + abs(vel)
              uneg    = vel - abs(vel)
              fe(i,j) = (upos*Velocity%u(i,j,k,n,taum1) + uneg*Velocity%u(i+1,j,k,n,taum1)) &
                         *Grd%umask(i,j,k)*Grd%umask(i+1,j,k)
           enddo
        enddo

        ! j-flux
        fn = 0.0
        do j=jsd,jec
           do i=isd,ied
              vel     = 0.5*Adv_vel%vhrho_nu(i,j,k)
              upos    = vel + abs(vel)
              uneg    = vel - abs(vel)
              fn(i,j) = (upos*Velocity%u(i,j,k,n,taum1) + uneg*Velocity%u(i,j+1,k,n,taum1)) &
                        *Grd%umask(i,j,k)*Grd%umask(i,j+1,k)
           enddo
        enddo

        ! compute horizontal divergence of flux 
        tmp(:,:) = BDX_EU(fe) + BDY_NU(fn) 


        ! remove advection of momentum but leave the metric term unchanged
        if(have_obc) then
          call ocean_obc_zero_boundary(tmp(:,:), 'C')
        endif

        do j=jsc,jec
          do i=isc,iec
            tmp(i,j) = tmp(i,j) &
                        + (3-2*n)*Thickness%rho_dzu(i,j,k,tau)*Velocity%u(i,j,k,3-n,tau) &
                         *(Grd%dh1dy(i,j)*Velocity%u(i,j,k,1,tau)-Grd%dh2dx(i,j)*Velocity%u(i,j,k,2,tau))
            wrk1_v(i,j,k,n) = Grd%umask(i,j,k)*tmp(i,j)
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
                  Velocity%advection(i,j,k,n,tau_m0) = Velocity%advection(i,j,k,n,tau_m0) + wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_hadv_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_hadv_v, wrk1_v(:,:,:,2))

      if(debug_this_module) then 
          write(stdoutunit,*) ' '
          write(stdoutunit,*) 'From ocean_velocity_advect_mod: horz_advection_upwind chksums'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('horz_advection_of_velocity(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('horz_advection_of_velocity(2)', wrk1_v(COMP,:,2))
      endif

  endif


end subroutine horz_advection_upwind
! </SUBROUTINE> NAME="horz_advection_upwind"



!#######################################################################
! <SUBROUTINE NAME="vert_advection_of_velocity">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to vertical advection of velocity.  
!
! Include vertical advection due to fresh water entering surface cells.  
!
! </DESCRIPTION>
! 
subroutine vert_advection_of_velocity(Time, Adv_vel, Velocity, pme, river, &
                                      upme, uriver, energy_analysis_step)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_velocity_type),    intent(inout) :: Velocity
  real, dimension(isd:,jsd:),   intent(in)    :: pme
  real, dimension(isd:,jsd:),   intent(in)    :: river
  real, dimension(isd:,jsd:,:), intent(in)    :: upme
  real, dimension(isd:,jsd:,:), intent(in)    :: uriver 
  logical,                      intent(in)    :: energy_analysis_step

  real, dimension(isd:ied,jsd:jed)   :: pme_u
  real, dimension(isd:ied,jsd:jed)   :: river_u
  real, dimension(isd:ied,jsd:jed,2) :: surf_accel 
  integer :: i, j, n


  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL,&
      '==>Error from ocean_velocity_advect_mod (vert_advection_of_velocity): module not yet initialized')
  endif 

  if(zero_velocity_advect_vert) then 
      Velocity%wrkv(:,:,:,:) = 0.0
      return 
  endif

  ! fresh water on U-cell
  pme_u    = 0.0
  river_u  = 0.0
  pme_u    = REMAP_BT_TO_BU(pme(:,:))
  river_u  = REMAP_BT_TO_BU(river(:,:))

  if(advection_scheme==VEL_ADVECT_2ND_ORDER) then 
    call vert_advection_centered(Time, Adv_vel, Velocity, pme, pme_u, river, river_u, &
                                 upme, uriver, energy_analysis_step)
  elseif(advection_scheme==VEL_ADVECT_UPWIND) then
    call vert_advection_upwind(Time, Adv_vel, Velocity, pme_u, river_u,   &
                               upme, uriver, energy_analysis_step)
  endif 
  

  ! diagnostics
  surf_accel = 0.0
  if(id_surf_accel(1) > 0 .or. id_surf_accel(2) > 0) then   
      do n=1,2
         do j=jsc,jec
            do i=isc,iec
               surf_accel(i,j,n) = Grd%umask(i,j,1)*(pme_u(i,j)*upme(i,j,n) + river_u(i,j)*uriver(i,j,n))
            enddo
         enddo
      enddo
  endif

  call diagnose_2d_u(Time, Grd, id_surf_accel(1), surf_accel(:,:,1))
  call diagnose_2d_u(Time, Grd, id_surf_accel(2), surf_accel(:,:,2))
  call diagnose_2d_u(Time, Grd, id_pme_u, pme_u(:,:))
  call diagnose_2d_u(Time, Grd, id_river_u, river_u(:,:))

end subroutine vert_advection_of_velocity
! </SUBROUTINE> NAME="vert_advection_of_velocity"


!#######################################################################
! <SUBROUTINE NAME="vert_advection_centered">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to vertical advection of velocity.  
!
! Include vertical advection due to fresh water entering surface cells.  
! Assume upme and uriver are the most fundamental fields, passed through
! the FMS coupler on the Bgrid velocity point. So use these to average 
! onto the C-grid velocity points.  Likewise, pme and river are fundamental,
! and sit on the T-grid.  
!
! Use second order centered method here.
!
! </DESCRIPTION>
! 
subroutine vert_advection_centered(Time, Adv_vel, Velocity, pme, pme_u, river, river_u, &
                                   upme, uriver, energy_analysis_step)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_velocity_type),    intent(inout) :: Velocity
  real, dimension(isd:,jsd:),   intent(in)    :: pme
  real, dimension(isd:,jsd:),   intent(in)    :: pme_u
  real, dimension(isd:,jsd:),   intent(in)    :: river
  real, dimension(isd:,jsd:),   intent(in)    :: river_u
  real, dimension(isd:,jsd:,:), intent(in)    :: upme
  real, dimension(isd:,jsd:,:), intent(in)    :: uriver 
  logical,                      intent(in)    :: energy_analysis_step

  real, dimension(isd:ied,jsd:jed) :: ft1
  real, dimension(isd:ied,jsd:jed) :: ft2

  integer :: i, j, k, kp1, n
  integer :: tau, tau_m0

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau    = Time%tau
  tau_m0 = Time%tau_m0
  wrk1_v = 0.0

  if(horz_grid == MOM_BGRID) then 

     do n=1,2

        ! fresh water contribution at the surface 
        do j=jsc,jec
           do i=isc,iec
              ft1(i,j) = -pme_u(i,j)*upme(i,j,n) -river_u(i,j)*uriver(i,j,n) 
              ft2(i,j) = 0.0
           enddo
        enddo

        do k=1,nk
           kp1 = min(k+1,nk)
           do j=jsc,jec
              do i=isc,iec
                 ft2(i,j) = Adv_vel%wrho_bu(i,j,k)*0.5*(Velocity%u(i,j,k,n,tau) + kmask(k)*Velocity%u(i,j,kp1,n,tau))
                 wrk1_v(i,j,k,n) = Grd%umask(i,j,k)*(ft1(i,j)-ft2(i,j))
                 ft1(i,j) = ft2(i,j)
              enddo
           enddo

           if(have_obc) then
             call ocean_obc_zero_boundary(wrk1_v(:,:,k,n), 'C')
           endif

        enddo
     enddo

  else  ! cgrid  

     n=1 ! zonal velocity 

     ! fresh water contribution at the surface 
     do j=jsc,jec
        do i=isc,iec
           ft1(i,j) = -onefourth*Grd%tmasken(i,j,1,n)*(pme(i,j)  +pme(i+1,j))  *(upme(i,j,n)  +upme(i,j-1,n)  ) &
                      -onefourth*Grd%tmasken(i,j,1,n)*(river(i,j)+river(i+1,j))*(uriver(i,j,n)+uriver(i,j-1,n)) 
           ft2(i,j) = 0.0
        enddo
     enddo
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsc,jec
           do i=isc,iec
              ft2(i,j) = onefourth*(Adv_vel%wrho_bt(i,j,k)+Adv_vel%wrho_bt(i+1,j,k)) &
                                  *(Velocity%u(i,j,k,n,tau) + kmask(k)*Velocity%u(i,j,kp1,n,tau))
              wrk1_v(i,j,k,n) = Grd%tmasken(i,j,k,n)*(ft1(i,j)-ft2(i,j))
              ft1(i,j) = ft2(i,j)
           enddo
        enddo
     enddo


     n=2 ! meridional velocity 

     ! fresh water contribution at the surface 
     do j=jsc,jec
        do i=isc,iec
           ft1(i,j) = -onefourth*Grd%tmasken(i,j,1,n)*(pme(i,j)  +pme(i,j+1))  *(upme(i,j,n)  +upme(i-1,j,n)  ) &
                      -onefourth*Grd%tmasken(i,j,1,n)*(river(i,j)+river(i,j+1))*(uriver(i,j,n)+uriver(i-1,j,n)) 
           ft2(i,j) = 0.0
        enddo
     enddo
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsc,jec
           do i=isc,iec
              ft2(i,j) = onefourth*(Adv_vel%wrho_bt(i,j,k)+Adv_vel%wrho_bt(i,j+1,k)) &
                                  *(Velocity%u(i,j,k,n,tau) + kmask(k)*Velocity%u(i,j,kp1,n,tau))
              wrk1_v(i,j,k,n) = Grd%tmasken(i,j,k,n)*(ft1(i,j)-ft2(i,j))
              ft1(i,j) = ft2(i,j)
           enddo
        enddo
     enddo


  endif  ! horz_grid if-test


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
                  Velocity%advection(i,j,k,n,tau_m0) = Velocity%advection(i,j,k,n,tau_m0) + wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_vadv_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_vadv_v, wrk1_v(:,:,:,2))

      if(debug_this_module) then
          write(stdoutunit,*) ' '
          write(stdoutunit,*) 'From ocean_velocity_advect_mod: vert_advection_centered chksums'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('vert_advection_of_velocity(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('vert_advection_of_velocity(2)', wrk1_v(COMP,:,2))
      endif

  endif

end subroutine vert_advection_centered
! </SUBROUTINE> NAME="vert_advection_centered"


!#######################################################################
! <SUBROUTINE NAME="vert_advection_upwind">
!
! <DESCRIPTION>
!
! Compute thickness weighted and density weighted acceleration 
! (kg/m^3)*(m^2/s^2) due to vertical advection of velocity.  
!
! Include vertical advection due to fresh water entering surface cells.  
!
! Use first order upwind method here. Not coded for Cgrid, since
! first order upwind method is rarely used. 
!
! </DESCRIPTION>
! 
subroutine vert_advection_upwind(Time, Adv_vel, Velocity, pme_u, river_u, &
                                 upme, uriver, energy_analysis_step)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_velocity_type),    intent(inout) :: Velocity
  real, dimension(isd:,jsd:),   intent(in)    :: pme_u
  real, dimension(isd:,jsd:),   intent(in)    :: river_u
  real, dimension(isd:,jsd:,:), intent(in)    :: upme
  real, dimension(isd:,jsd:,:), intent(in)    :: uriver 
  logical,                      intent(in)    :: energy_analysis_step

  real, dimension(isd:ied,jsd:jed) :: ft1
  real, dimension(isd:ied,jsd:jed) :: ft2

  integer :: i, j, k, kp1, n
  integer :: taum1, tau_m0
  real    :: vel, wpos, wneg

  integer :: stdoutunit 
  stdoutunit=stdout() 

  taum1  = Time%taum1
  tau_m0 = Time%tau_m0
  wrk1_v = 0.0

  do n=1,2

     ! fresh water contribution at the surface 
     do j=jsc,jec
        do i=isc,iec
           ft1(i,j) = -pme_u(i,j)*upme(i,j,n) -river_u(i,j)*uriver(i,j,n)
           ft2(i,j) = 0.0
        enddo
     enddo

     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsc,jec
           do i=isc,iec
              vel      = 0.5*Adv_vel%wrho_bu(i,j,k)
              wpos     = vel + abs(vel) 
              wneg     = vel - abs(vel) 
              ft2(i,j) = (wneg*Velocity%u(i,j,k,n,taum1) + kmask(k)*wpos*Velocity%u(i,j,kp1,n,taum1)) &
                         *Grd%umask(i,j,k)*Grd%umask(i,j,kp1) 
              wrk1_v(i,j,k,n) = Grd%umask(i,j,k)*(ft1(i,j)-ft2(i,j))
              ft1(i,j) = ft2(i,j)
           enddo
        enddo

        if(have_obc) then
          call ocean_obc_zero_boundary(wrk1_v(:,:,k,n), 'C')
        endif

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
                  Velocity%advection(i,j,k,n,tau_m0) = Velocity%advection(i,j,k,n,tau_m0) + wrk1_v(i,j,k,n) 
               enddo
            enddo
         enddo
      enddo

      call diagnose_3d_u(Time, Grd, id_vadv_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_vadv_v, wrk1_v(:,:,:,2))

      if(debug_this_module) then
          write(stdoutunit,*) ' '
          write(stdoutunit,*) 'From ocean_velocity_advect_mod: vert_advection_upwind chksums'
          call write_timestamp(Time%model_time)
          call write_chksum_3d('vert_advection_of_velocity(1)', wrk1_v(COMP,:,1))
          call write_chksum_3d('vert_advection_of_velocity(2)', wrk1_v(COMP,:,2))
      endif

  endif


end subroutine vert_advection_upwind
! </SUBROUTINE> NAME="vert_advection_upwind"





end module ocean_velocity_advect_mod

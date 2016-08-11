
module bgrid_conserve_energy_mod

!------------------------------------------------------------------
!  Enforces energy conservation by computing a correction term
!------------------------------------------------------------------
! This module computes a correction to the temperature tendency to
! conserve the globally-averaged total energy, TE = cp*T + KE.
! The assumption is made that the input tendencies are only
! due to the dynamical core.
!
! The global average of TE is weighted by pressure and area,
! and therefore the time rate of change of pressure must also be
! considered. The correction is applied as a global constant value.
!
! Diagnostics fields that are available (through the diag_manager)
! have few spatial patterns because of the uniform global correction.
!------------------------------------------------------------------

use       bgrid_horiz_mod, only: horiz_grid_type, update_np, VELGRID => VGRID
use        bgrid_vert_mod, only: vert_grid_type, compute_pres_depth
use       bgrid_masks_mod, only: grid_mask_type
use    bgrid_prog_var_mod, only: prog_var_type
use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID

use      time_manager_mod, only: time_type
use      diag_manager_mod, only: register_diag_field, send_data
use               fms_mod, only: error_mesg, FATAL, stdlog, &
                                 write_version_number,      &
                                 mpp_pe, mpp_root_pe,       &
                                 mpp_clock_id,              &
                                 mpp_clock_begin,           &
                                 mpp_clock_end,             &
                                 MPP_CLOCK_SYNC, CLOCK_MODULE
use         constants_mod, only: CP_AIR, GRAV
use       mpp_domains_mod, only: mpp_global_sum, BITWISE_EXACT_SUM

 implicit none
 private

 public :: bgrid_conserve_energy_init, bgrid_conserve_energy, &
           bgrid_conserve_energy_end

!------------------------------------------------------------------
! private module data
 character(len=128) :: version = '$Id: bgrid_conserve_energy.F90,v 19.0 2012/01/06 19:53:25 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'

! saves the total correction for the entire run
 real :: total_correction

! for diagnostics
 character(len=8) :: mod_name = 'dynamics'
 integer :: id_tdt_diss, id_diss_heat

 integer :: id_clock
 logical :: do_init = .true.

contains
!################################################################

 subroutine bgrid_conserve_energy ( dt, Time, Hgrid, Vgrid, Masks, &
                                    Var, Var_dt )
 real,                  intent(in)    :: dt
 type      (time_type), intent(in)    :: Time
 type(horiz_grid_type), intent(in)    :: Hgrid
 type (vert_grid_type), intent(in)    :: Vgrid
 type (grid_mask_type), intent(in)    :: Masks
 type  (prog_var_type), intent(in)    :: Var
 type  (prog_var_type), intent(inout) :: Var_dt

 real, dimension(Hgrid%ilb:Hgrid%iub, &
                 Hgrid%jlb:Hgrid%jub, Vgrid%nlev) :: tcor, vcor, pwt, dpde, dpde_vel, &
                                                     dpde_old, dpde_old_vel
 real, dimension(Hgrid%ilb:Hgrid%iub, Hgrid%jlb:Hgrid%jub) :: diss_heat, pssl_new
 real    :: sum_tcor, sum_vcor, sum_pwt, correction
 integer :: j, k
 logical :: used

    if (do_init) call error_mesg ('bgrid_conserve_energy', &
                                  'initializtion not called', FATAL)
    call mpp_clock_begin (id_clock)

  ! compute mass weights (at tau and tau+1)
    pssl_new = Var%pssl + dt*Var_dt%pssl
    call compute_pres_depth ( Vgrid, Var%pssl, dpde_old )
    call compute_pres_depth ( Vgrid, pssl_new, dpde     )
    call change_grid ( Hgrid, TEMP_GRID, WIND_GRID, dpde_old,     dpde,    &
                                                    dpde_old_vel, dpde_vel )

  ! compute local correction terms
  ! also need to take into account time rate of change of pressure
    vcor =  ((Var%u+.5*dt*Var_dt%u)*Var_dt%u + &
             (Var%v+.5*dt*Var_dt%v)*Var_dt%v) * dpde_vel + &
            0.5*(Var%u**2+Var%v**2)*(dpde_vel-dpde_old_vel)/dt
    tcor = Var_dt%t*dpde + Var%t*(dpde-dpde_old)/dt
    do k = 1, Vgrid%nlev
    do j = Hgrid%jlb, Hgrid%jub
      vcor(:,j,k) = Hgrid%Vel%area(j) * vcor (:,j,k)
      tcor(:,j,k) = Hgrid%Tmp%area(j) * tcor (:,j,k)
      pwt (:,j,k) = Hgrid%Tmp%area(j) * dpde (:,j,k)
    enddo
    enddo

  ! compute bit-reproducible global sums
  ! pass only data domain for velocity grid
    sum_vcor = mpp_global_sum ( Hgrid%Vel%Domain, vcor(:,Hgrid%Vel%jsd:Hgrid%Vel%jed,:), &
                                                        flags=BITWISE_EXACT_SUM )
    sum_tcor = mpp_global_sum ( Hgrid%Tmp%Domain, tcor, flags=BITWISE_EXACT_SUM )
    sum_pwt  = mpp_global_sum ( Hgrid%Tmp%Domain, pwt,  flags=BITWISE_EXACT_SUM )

  ! global correction
    correction = -( sum_tcor + sum_vcor/CP_AIR ) / sum_pwt

  ! add on tendency
    Var_dt%t = Var_dt%t + correction
  ! save total correction (as a potential diagnostic)
    total_correction = total_correction + correction

  !------ diagnostics section ------

    if ( id_tdt_diss > 0 .or. id_diss_heat > 0 ) then
         tcor = correction
         if (.not.Masks%sigma) tcor = tcor*Masks%Tmp%mask ! zero-out below surface
    endif

    if ( id_tdt_diss > 0 ) then
     ! this will be a very boring 3d diagnostic
       used = send_data ( id_tdt_diss, tcor(Hgrid%Tmp%is:Hgrid%Tmp%ie,Hgrid%Tmp%js:Hgrid%Tmp%je,:), &
        Time, mask=Masks%Tmp%mask(Hgrid%Tmp%is:Hgrid%Tmp%ie,Hgrid%Tmp%js:Hgrid%Tmp%je,:) > 0.5 )
    endif

    if ( id_diss_heat > 0 ) then
      ! vertical integral of ke dissipation
        diss_heat = CP_AIR/GRAV * sum( tcor*dpde, 3 ) 
        used = send_data ( id_diss_heat, &
                           diss_heat(Hgrid%Tmp%is:Hgrid%Tmp%ie,Hgrid%Tmp%js:Hgrid%Tmp%je), Time )
    endif

    call mpp_clock_end (id_clock)

 end subroutine bgrid_conserve_energy

!################################################################

 subroutine bgrid_conserve_energy_init ( Time, axes )
 type(time_type), intent(in) :: Time
 integer,         intent(in) :: axes(3)

   call write_version_number(version, tagname)

   id_clock = mpp_clock_id ('BGRID: conserve_energy', &
                  flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE)

! initialize diagnostics related to the energy conservation

   id_tdt_diss = register_diag_field ( mod_name, 'tdt_diss_dynam', axes(1:3), Time, &       
                           'Dissipative heating from dynamical core', & 
                                 'deg_k/s', missing_value=-1.e10   )
       
   id_diss_heat = register_diag_field &
        ( mod_name, 'diss_heat_dynam', axes(1:2), Time,      &
         'Integrated dissipative heating for dynamical core', & 
         'W/m2' )

   total_correction = 0.

   do_init = .false.

 end subroutine bgrid_conserve_energy_init

!################################################################

 subroutine bgrid_conserve_energy_end

 integer :: logunit 
 ! rewrite module name and version info to log file
   call write_version_number(version, tagname)
 ! write diagnostic to log file
   logunit = stdlog()
   if (mpp_pe() == mpp_root_pe()) write (logunit,*)  &
    'Total dissipative heating from dynamical core (deg/s) = ', total_correction

 end subroutine bgrid_conserve_energy_end

!################################################################

end module bgrid_conserve_energy_mod


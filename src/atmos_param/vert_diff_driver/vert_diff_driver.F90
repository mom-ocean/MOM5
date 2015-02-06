
module vert_diff_driver_mod

!-----------------------------------------------------------------------
!   module performs vertical diffusion of atmospheric variables
!-----------------------------------------------------------------------

use    vert_diff_mod, only:  surf_diff_type,     &
                             vert_diff_init, &
                             gcm_vert_diff_down, &
                             gcm_vert_diff_up,   &
                             vert_diff_end


use diag_manager_mod, only:  register_diag_field, send_data

use time_manager_mod, only:  time_type

use          mpp_mod, only:  input_nml_file
use          fms_mod, only:  file_exist, open_namelist_file, error_mesg,  &
                             check_nml_error, FATAL, mpp_pe, mpp_root_pe, &
                             close_file, write_version_number, stdlog

use    constants_mod, only:  CP_AIR, GRAV

use   field_manager_mod, only: MODEL_ATMOS
use  tracer_manager_mod, only: get_number_tracers, get_tracer_names, &
                               NO_TRACER

!-----------------------------------------------------------------------

implicit none
private

public :: vert_diff_driver_down, vert_diff_driver_up,  &
          vert_diff_driver_init, vert_diff_driver_end
public :: surf_diff_type


!-----------------------------------------------------------------------
!---- namelist ----

logical :: do_conserve_energy         = .false.
logical :: do_mcm_no_neg_q            = .false.
logical :: use_virtual_temp_vert_diff = .true.
logical :: do_mcm_plev                = .false.
logical :: do_mcm_vert_diff_tq        = .false.

namelist /vert_diff_driver_nml/ do_conserve_energy,         &
                                do_mcm_no_neg_q,            &
                                use_virtual_temp_vert_diff, &
                                do_mcm_plev, do_mcm_vert_diff_tq

!-----------------------------------------------------------------------
! tracer storage is used 
type :: tracer_storage_type
   integer :: id_tr_dt     = 0 ! diag id of the tracer tendency due 
                               ! to vert diff
   integer :: id_tr_dt_int = 0 ! diag id of the vertically-integrated 
                               ! tracer tendency
   real, pointer :: &
        buffer(:,:,:) => NULL() ! buffer for tendency calculations
end type
type(tracer_storage_type), allocatable :: tr_store(:)

real, allocatable, dimension(:,:,:) :: dt_t_save, dt_q_save

!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_vdif, id_qdt_vdif, id_udt_vdif, id_vdt_vdif,  &
           id_sens_vdif, id_evap_vdif,                          &
           id_tdt_diss_vdif, id_diss_heat_vdif

real :: missing_value = -999.

character(len=9), parameter :: mod_name = 'vert_diff'

!-----------------------------------------------------------------------
!---- version number ----

character(len=128) :: version = '$Id: vert_diff_driver.F90,v 20.0 2013/12/13 23:22:29 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

logical :: module_is_initialized = .false.


contains

!#######################################################################

 subroutine vert_diff_driver_down (is, js, Time, delt, p_half, p_full, &
                                   z_full, diff_mom, diff_heat,        &
                                   u, v, t, q, trs,                    &
                                   dtau_du, dtau_dv, tau_x, tau_y,     &
                                   dt_u, dt_v, dt_t, dt_q, dt_trs,     &
                                   Surf_diff, diff_t_clubb, mask, kbot ) !cjg

integer, intent(in)                     :: is, js
type(time_type),   intent(in)           :: Time
real, intent(in)                        :: delt
real, intent(in)   , dimension(:,:,:)   :: p_half, p_full, z_full,  &
                                           diff_mom, diff_heat
real, intent(in),    dimension(:,:,:)   :: u, v, t, q
real, intent(in),    dimension(:,:,:,:) :: trs
real, intent(in),    dimension(:,:)     :: dtau_du, dtau_dv

real, intent(inout), dimension(:,:)     :: tau_x, tau_y
real, intent(inout), dimension(:,:,:)   :: dt_u, dt_v, dt_t, dt_q
real, intent(inout), dimension(:,:,:,:) :: dt_trs

type(surf_diff_type), intent(inout)     :: Surf_diff

!-->cjg
real, intent(in)   , dimension(:,:,:), optional :: diff_t_clubb
!<--cjg
real   , intent(in), dimension(:,:,:), optional :: mask
integer, intent(in), dimension(:,:),   optional :: kbot

real, dimension(size(t,1),size(t,2),size(t,3)) :: tt, dpg, q_2
real, dimension(size(t,1),size(t,2),size(t,3)) :: dissipative_heat
integer :: k, ntp, tr
logical :: used
real, dimension(size(t,1),size(t,2)) :: diag2
integer :: ie, je

!-----------------------------------------------------------------------

  if (.not. module_is_initialized) call error_mesg       &
                  ('vert_diff_driver_mod',  &
                   'vert_diff_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------

  ntp = size(dt_trs,4) ! number of prognostic tracers
  if (size(trs,4) < ntp) call error_mesg ('vert_diff_driver', &
             'Number of tracers .lt. number of tracer tendencies',FATAL)

  ie = is + size(t,1) -1
  je = js + size(t,2) -1


    if(do_mcm_vert_diff_tq) then
      dt_t_save(is:ie,js:je,:) = dt_t
      dt_q_save(is:ie,js:je,:) = dt_q
      dt_t = 0.0
      dt_q = 0.0
    endif

!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t, dt_q, dt_u, and dt_v at this point add 
!-----in the negative value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------




!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, -2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, -2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif


!------- diagnostics for dt/dt_diff -------
    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, -3.*dt_t, Time, is, js, 1, &
                           rmask=mask )
    endif
    
!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, -2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

    ! store values of tracer tendencies before the vertical diffusion, if 
    ! requested -- availability of storage serves as an indicatior that 
    ! storing is necessary
    do tr = 1,size(tr_store(:))
       if( associated(tr_store(tr)%buffer) ) &
            tr_store(tr)%buffer(is:ie,js:je,:) = dt_trs(:,:,:,tr)
    enddo

!-----------------------------------------------------------------------
!---- local temperature ----
!     (avoids divid by zero when computing virtual temp)

   tt = t
   if (present(mask)) where (mask < 0.5) tt = 200.

!-----------------------------------------------------------------------
!---- momentum diffusion ----
!---- heat/moisture diffusion (down only) ----
!---- tracer diffusion (no surface flux) ----

 q_2 = q
 if (do_mcm_no_neg_q) then
   where (q_2 < 0.0)  q_2 = 0.0
 endif

!--> cjg
 call gcm_vert_diff_down (is, js, delt, u, v, tt, q_2, trs(:,:,:,1:ntp), &
                          diff_mom, diff_heat,                           &
                          p_half, p_full, z_full,                        &
                          tau_x, tau_y, dtau_du, dtau_dv,                &
                          dt_u, dt_v, dt_t, dt_q, dt_trs(:,:,:,1:ntp),   &
                          dissipative_heat, Surf_diff,  diff_t_clubb, kbot ) !cjg
!<--cjg

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- to do diagnostics on dt_u, and dt_v at this point add 
!-----in the value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 .or. id_diss_heat_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/GRAV
            enddo
            if (present(mask)) dpg = dpg*mask
    endif
    

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, -2.*CP_AIR*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, -2.*diag2, Time, is, js )
    endif
    
    

!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, 2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, 2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif
    
!------- diagnostics for dissipative heating -------
    if ( id_tdt_diss_vdif > 0 ) then
       used = send_data ( id_tdt_diss_vdif, dissipative_heat, Time, &
                          is, js, 1, &
                          rmask=mask)
    endif

!------- diagnostics for dt/dt_diff -------
    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, -3.*dissipative_heat, Time, is, js, 1, &
                           rmask=mask )
    endif

!------- diagnostics for vertically integrated dissipative heating -------
    if ( id_diss_heat_vdif > 0 ) then
          diag2 = sum( CP_AIR*dissipative_heat*dpg, 3 )
          used = send_data ( id_diss_heat_vdif, diag2, Time, is, js )
    endif

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_down

!#######################################################################

 subroutine vert_diff_driver_up (is, js, Time, delt, p_half, &
                                 Surf_diff, dt_t, dt_q, dt_tr, mask, kbot)

 integer,           intent(in)            :: is, js
 type(time_type),   intent(in)            :: Time
 real,    intent(in)                      :: delt
 real,    intent(in),    dimension(:,:,:) :: p_half
 type(surf_diff_type),   intent(in)       :: Surf_diff
 real,    intent(inout), dimension(:,:,:) :: dt_t, dt_q
 real,    intent(inout), dimension(:,:,:,:) :: dt_tr
 real   , intent(in), dimension(:,:,:), optional :: mask
 integer, intent(in),    dimension(:,:), optional :: kbot

 integer :: k, tr
 logical :: used
 real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: dpg
 real, dimension(size(p_half,1),size(p_half,2)) :: diag2
 integer :: ie, je

!-----------------------------------------------------------------------
    ie = is + size(p_half,1) -1
    je = js + size(p_half,2) -1
!-----------------------------------------------------------------------

    call gcm_vert_diff_up (is, js, delt, Surf_diff, dt_t, dt_q, dt_tr, kbot)


!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t and dt_q at this point add in the
!---- the postive value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------
!------- diagnostics for dt/dt_diff -------

    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, 3.*dt_t, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, 2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/GRAV
            enddo
            if (present(mask)) dpg = dpg*mask
    endif

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, 2.*CP_AIR*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, 2.*diag2, Time, is, js )
    endif

!------- diagnostics of tracer tendencies ------- 
    do tr = 1, size(tr_store(:))
       if(tr_store(tr)%id_tr_dt > 0) then
          used = send_data(tr_store(tr)%id_tr_dt, &
               dt_tr(:,:,:,tr)-tr_store(tr)%buffer(is:ie,js:je,:),Time,is,js)
       endif
       
       if(tr_store(tr)%id_tr_dt_int > 0) then
          diag2 = sum((dt_tr(:,:,:,tr)-tr_store(tr)%buffer(is:ie,js:je,:))*dpg,3)
          used = send_data(tr_store(tr)%id_tr_dt_int, diag2, Time, is, js)
       endif
    enddo

    if(do_mcm_vert_diff_tq) then
      dt_t = dt_t + dt_t_save(is:ie,js:je,:)
      dt_q = dt_q + dt_q_save(is:ie,js:je,:)
    endif

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_up

!#######################################################################

 subroutine vert_diff_driver_init ( Surf_diff, idim, jdim, kdim,  &
                                    axes, Time, do_clubb )  !cjg

 type(surf_diff_type), intent(inout) :: Surf_diff
 integer             , intent(in)    :: idim, jdim, kdim, axes(4)
 type(time_type)     , intent(in)    :: Time
!-->cjg
 integer, intent(in)                 :: do_clubb
!<--cjg

 integer :: unit, io, ierr, tr, logunit
 integer :: ntprog ! number of prognostic tracers in the atmosphere
 character(len=32)  :: name, units ! name of the tracer
 character(len=128) :: longname    ! long name of the tracer

!-----------------------------------------------------------------------
!------ read namelist ------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=vert_diff_driver_nml, iostat=io)
   ierr = check_nml_error(io,'vert_diff_driver_nml')
#else   
   if ( file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=vert_diff_driver_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'vert_diff_driver_nml')
      enddo
 10   call close_file (unit)
   endif
#endif

!--------- write version number and namelist ------------------

   call write_version_number ( version, tagname )
   logunit = stdlog()
   if(mpp_pe() == mpp_root_pe() ) write(logunit,nml=vert_diff_driver_nml)

!-------- initialize gcm vertical diffusion ------

   call vert_diff_init (Surf_diff, idim, jdim, kdim, do_conserve_energy, &
                        use_virtual_temp_vert_diff, do_mcm_plev, do_clubb) !cjg

!-----------------------------------------------------------------------

   if(do_mcm_vert_diff_tq) then
     allocate(dt_t_save(idim,jdim,kdim)) ; dt_t_save = 0.0
     allocate(dt_q_save(idim,jdim,kdim)) ; dt_q_save = 0.0
   endif

!--------------- initialize diagnostic fields --------------------

   id_tdt_vdif = &
   register_diag_field ( mod_name, 'tdt_vdif', axes(1:3), Time, &
                        'Temperature tendency from vert diff',  &
                        'deg_K/s', missing_value=missing_value  )

   id_qdt_vdif = &
   register_diag_field ( mod_name, 'qdt_vdif', axes(1:3), Time, &
                        'Spec humidity tendency from vert diff',&
                        'kg/kg/s', missing_value=missing_value  )

   id_udt_vdif = &
   register_diag_field ( mod_name, 'udt_vdif', axes(1:3), Time, &
                        'Zonal wind tendency from vert diff',   &
                        'm/s2', missing_value=missing_value     )

   id_vdt_vdif = &
   register_diag_field ( mod_name, 'vdt_vdif', axes(1:3), Time,    &
                        'Meridional wind tendency from vert diff', &
                        'm/s2', missing_value=missing_value        )

   id_sens_vdif = &
   register_diag_field ( mod_name, 'sens_vdif', axes(1:2), Time,  &
                        'Integrated heat flux from vert diff',    &
                        'W/m2' )

   id_evap_vdif = &
   register_diag_field ( mod_name, 'evap_vdif', axes(1:2), Time,    &
                        'Integrated moisture flux from vert diff',  &
                        'kg/m2/s' )

   id_tdt_diss_vdif = &
   register_diag_field ( mod_name, 'tdt_diss_vdif', axes(1:3), Time,  &
                        'Dissipative heating from vert_diff', 'deg_K/s', &
                         missing_value=missing_value  ) 

   id_diss_heat_vdif = &
   register_diag_field ( mod_name, 'diss_heat_vdif', axes(1:2), Time,  &
                        'Integrated dissipative heating from vert diff',  &
                        'W/m2' )

   ! initialize diagnostics tracers
   call get_number_tracers(MODEL_ATMOS, num_prog=ntprog)
   allocate(tr_store(ntprog))
   do tr = 1,ntprog
      call get_tracer_names( MODEL_ATMOS, tr, name, longname, units )
      tr_store(tr)%id_tr_dt = &
        register_diag_field ( mod_name, trim(name)//'dt_vdif', axes(1:3), Time, &
           'Tendency of '//trim(longname)//' from vert diff', trim(units)//'/s')
      tr_store(tr)%id_tr_dt_int = &
        register_diag_field ( mod_name, trim(name)//'dt_vint_vdif', axes(1:2), Time, &
           'Integrated tendency of '//trim(longname)//' from vert diff',&
           trim(units)//' kg/(m2 s)')

      if(tr_store(tr)%id_tr_dt>0 .or.tr_store(tr)%id_tr_dt_int>0 ) &
           allocate(tr_store(tr)%buffer(idim,jdim,kdim))
   enddo

!-----------------------------------------------------------------------

   module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_init

!#######################################################################

 subroutine vert_diff_driver_end

   integer :: tr ! tracer index

   call vert_diff_end
   if(do_mcm_vert_diff_tq) deallocate(dt_t_save, dt_q_save)
   ! deallocate tracer diagnostics storage
   do tr = 1,size(tr_store(:))
      if(associated(tr_store(tr)%buffer)) &
           deallocate(tr_store(tr)%buffer)
   enddo
   deallocate(tr_store)

!-----------------------------------------------------------------------

   module_is_initialized = .false.

!-----------------------------------------------------------------------
 end subroutine vert_diff_driver_end

!#######################################################################

end module vert_diff_driver_mod


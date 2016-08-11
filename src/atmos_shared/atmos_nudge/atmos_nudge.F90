module atmos_nudge_mod

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: check_nml_error, close_file, &
                   stdlog, mpp_pe, mpp_root_pe, write_version_number, &
                   error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type, set_time, get_date, &
                            operator( + ), operator( < )
use data_override_mod,only: data_override
use diag_manager_mod, only: register_diag_field, send_data
use mpp_mod, only: mpp_min,mpp_max

implicit none
private

public :: atmos_nudge_init, get_atmos_nudge, atmos_nudge_end, do_ps

character(len=128), parameter :: version = '$Id: atmos_nudge.F90,v 19.0 2012/01/06 20:28:34 fms Exp $'
character(len=128), parameter :: tagname = '$Name: tikal $'

logical :: module_is_initialized = .false.

integer :: freq = 0   ! frequency in seconds
real ::  u_tau = -1.  ! relaxation time in seconds (no insertion if < 0)
real ::  v_tau = -1.
real ::  t_tau = -1.
real ::  q_tau = -1.
real :: ps_tau = -1.
integer :: skip_top_v = 2            ! momentum
integer :: skip_bot_v = 0
integer :: skip_top_t = 0            ! temperature
integer :: skip_bot_t = 21
integer :: skip_top_q = 8            ! specific humidity
integer :: skip_bot_q = 0

namelist /atmos_nudge_nml/ freq, u_tau, v_tau, t_tau, q_tau, ps_tau, &
                           skip_top_v, skip_bot_v,               &
                           skip_top_t, skip_bot_t,               &
                           skip_top_q, skip_bot_q

type(time_type) :: Time_next
integer :: id_udt, id_vdt, id_tdt, id_qdt, id_psdt
logical :: do_u, do_v, do_t, do_q, do_ps

contains

!-----------------------------------------------------------------------

subroutine get_atmos_nudge(Time, dt, beglon, endlon, beglat, endlat, nlev,  &
                           ng, ps, u, v, t, q, psdt, udt, vdt, tdt, qdt )
type (time_type),       intent(in)    :: Time
real,                   intent(in)    :: dt

integer, intent(in):: beglon, endlon, beglat, endlat, nlev, ng

real, intent(inout):: ps(beglon:endlon, beglat:endlat)
real, intent(inout):: psdt(beglon:endlon, beglat:endlat)

real, intent(inout):: u(beglon:endlon, beglat:endlat, nlev)
real, intent(inout):: v(beglon:endlon, beglat:endlat, nlev)
real, intent(inout):: t(beglon:endlon, beglat-ng:endlat+ng, nlev)
real, intent(inout):: q(beglon:endlon, beglat-ng:endlat+ng, nlev)

real, intent(inout):: udt(beglon:endlon, beglat:endlat, nlev)
real, intent(inout):: vdt(beglon:endlon, beglat:endlat, nlev)
real, intent(inout):: tdt(beglon:endlon, beglat:endlat, nlev)
real, intent(inout):: qdt(beglon:endlon, beglat:endlat, nlev)

real ::  obs(beglon:endlon, beglat:endlat, nlev)
real :: tend(beglon:endlon, beglat:endlat, nlev)
real ::  obs2(beglon:endlon, beglat:endlat)
real :: tend2(beglon:endlon, beglat:endlat)
real :: factor(nlev,3)
logical :: sent, done
integer :: i,j,k

   if (.not.module_is_initialized) then
       call error_mesg ('atmos_nudge_mod', 'module not initialized', FATAL)
   endif

 ! no data forcing override
   if (freq <= 0) then
       return
   endif

 ! is it time for data forcing
   if (Time < Time_next) then
       return
   endif
   Time_next = Time_next + set_time(freq)

! vertically dependent Tau
! very crude - zero at top/bottom + linear increase with level downward/upward

   factor = 1.

!------------------------------------------------------------------
! Momentum:
   if (skip_top_v > 0) then
      factor(1,1) = 0.
      do k = 2, skip_top_v
         factor(k,1) = factor(k-1,1) + 1./real(skip_top_v)
      enddo
   endif
   if (skip_bot_v > 0) then
      factor(nlev,1) = 0.
      do k = nlev-1, nlev-skip_bot_v+1, -1
         factor(k,1) = factor(k+1,1) + 1./real(skip_bot_v)
      enddo
   endif

! temperature
   if (skip_top_t > 0) then
      factor(1,2) = 0.
      do k = 2, skip_top_t
         factor(k,2) = factor(k-1,2) + 1./real(skip_top_t)
      enddo
   endif
   if (skip_bot_t > 0) then
         factor(nlev-skip_bot_t-1,2) = 0.5
         factor(nlev-skip_bot_t,  2) = 0.25
      do k=nlev-skip_bot_t+1,nlev
         factor(k,2) = 0.
      enddo
   endif
   
! Specific humidity
   if (skip_top_q > 0) then
      do k = 1, skip_top_q
         factor(k,3) = 0.
      enddo
         factor(skip_top_q+1,3) = 0.25
         factor(skip_top_q+2,3) = 0.5
   endif
   if (skip_bot_q > 0) then
      factor(nlev,3) = 0.
      do k = nlev-1, nlev-skip_bot_q+1, -1
         factor(k,3) = factor(k+1,3) + 1./real(skip_bot_q)
      enddo
   endif
!------------------------------------------------------------------

! zonal wind component
   if (do_u .or. id_udt>0) then
       call data_override ('ATM', 'u_obs', obs, Time, override=done)
       if (.not.done) call override_error (Time,'zonal wind')
   endif
   if (do_u) then
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = (obs(i,j,k) - u(i,j,k)) / (u_tau + dt) * factor(k,1)
                   u(i,j,k) = u(i,j,k) + dt*tend(i,j,k)
                 udt(i,j,k) = udt(i,j,k) + tend(i,j,k)
             enddo
          enddo
       enddo
   else
       if (id_udt > 0) then
!          tend = 0.
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
! Report the error if not nudging
                tend(i,j,k) = u(i,j,k) - obs(i,j,k)
             enddo
          enddo
       enddo
       endif
   endif
   if (id_udt > 0) sent = send_data (id_udt, tend, Time) ! masking?

! meridional wind component
   if (do_v .or. id_vdt>0) then
       call data_override ('ATM', 'v_obs', obs, Time, override=done)
       if (.not.done) call override_error (Time,'meridional wind')
   endif
   if (do_v) then
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = (obs(i,j,k) - v(i,j,k)) / (v_tau + dt) * factor(k,1)
                   v(i,j,k) = v(i,j,k) + dt*tend(i,j,k)
                 vdt(i,j,k) = vdt(i,j,k) + tend(i,j,k)
             enddo
          enddo
       enddo
   else
       if (id_vdt > 0) then
!          tend = 0.
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = v(i,j,k) - obs(i,j,k)
             enddo
          enddo
       enddo
       endif
   endif
   if (id_vdt > 0) sent = send_data (id_vdt, tend, Time) ! masking?

! temperature
   if (do_t .or. id_tdt>0) then
       call data_override ('ATM', 't_obs', obs, Time, override=done)
       if (.not.done) call override_error (Time,'temperature')
   endif
   if (do_t) then
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = (obs(i,j,k) - t(i,j,k)) / (t_tau + dt) * factor(k,2)
                   t(i,j,k) = t(i,j,k) + dt*tend(i,j,k)
                 tdt(i,j,k) = tdt(i,j,k) + tend(i,j,k)
             enddo
          enddo
       enddo
   else
       if (id_tdt > 0) then
!          tend = 0.
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = t(i,j,k) - obs(i,j,k)
             enddo
          enddo
       enddo
       endif
   endif
   if (id_tdt > 0) sent = send_data (id_tdt, tend, Time) ! masking?

! specific humidity
   if (do_q .or. id_qdt>0) then
       call data_override ('ATM', 'q_obs', obs, Time, override=done)
       if (.not.done) call override_error (Time,'specific humidity')
   endif
   if (do_q) then
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = (max(1.e-8,obs(i,j,k)) - q(i,j,k)) / (q_tau + dt) * factor(k,3)
                   q(i,j,k) = q(i,j,k) + dt*tend(i,j,k)
                 qdt(i,j,k) = qdt(i,j,k) + tend(i,j,k)
             enddo
          enddo
       enddo
   else
       if (id_qdt > 0) then
!          tend = 0.
       do k=1,nlev
          do j=beglat,endlat
             do i=beglon,endlon
                tend(i,j,k) = q(i,j,k) - obs(i,j,k)
             enddo
          enddo
       enddo
       endif
   endif
   if (id_qdt > 0) sent = send_data (id_qdt, tend, Time) ! masking?

! surface pressure
   if (do_ps .or. id_psdt>0) then
       call data_override ('ATM', 'ps_obs', obs2, Time, override=done)
       if (.not.done) call override_error (Time,'surface pressure')
   endif
   if (do_ps) then
       do j=beglat,endlat
          do i=beglon,endlon
             tend2(i,j) = (obs2(i,j) - ps(i,j)) / (ps_tau + dt)
                ps(i,j) = ps(i,j) + dt*tend2(i,j)
              psdt(i,j) = psdt(i,j) + tend2(i,j)
          enddo
       enddo
   else
       if (id_psdt > 0) then
!          tend2 = 0.
          do j=beglat,endlat
             do i=beglon,endlon
                tend2(i,j) = ps(i,j) - obs2(i,j)
             enddo
          enddo
       endif
   endif
   if (id_psdt > 0) sent = send_data (id_psdt, tend2, Time) ! masking?
!if (mpp_pe()==mpp_root_pe()) print *, 'Leaving atmos_nudge'

end subroutine get_atmos_nudge

!-----------------------------------------------------------------------

subroutine atmos_nudge_init ( Time, axes, flag )
type (time_type),      intent(in)  :: Time
integer, dimension(3), intent(in)  :: axes
logical, optional,     intent(out) :: flag
integer :: ierr, io, unit, logunit
real :: eps
character(len=64) :: desc
real :: missing_value = -1.e10

 ! read namelist
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=atmos_nudge_nml, iostat=io)
   ierr = check_nml_error(io, 'atmos_nudge_nml')
#else
   unit = open_namelist_file()
   ierr=1  
   do while (ierr /= 0)
     read (unit, nml=atmos_nudge_nml, iostat=io, end=10) 
     ierr = check_nml_error (io, 'atmos_nudge_nml')
   enddo   
10 call close_file (unit)
#endif
   call write_version_number(version, tagname)
   logunit=stdlog()
   if (mpp_pe() == mpp_root_pe()) write (logunit, nml=atmos_nudge_nml)

 ! initialize flags
   eps = 1.e-10
   do_u  = .false.; if ( u_tau > -eps) do_u  = .true.
   do_v  = .false.; if ( v_tau > -eps) do_v  = .true.
   do_t  = .false.; if ( t_tau > -eps) do_t  = .true.
   do_q  = .false.; if ( q_tau > -eps) do_q  = .true.
   do_ps = .false.; if (ps_tau > -eps) do_ps = .true.

 ! namelist dummy checks
 ! if no overrides turned on then set freq = 0
   if (freq > 0) then
       if ( .not.do_u .and. .not.do_v .and. .not.do_t .and. &
            .not.do_q .and. .not.do_ps ) then
!           call error_mesg ('atmos_nudge_mod', 'no variables specified '//&
!                            'for override, resetting freq = 0', WARNING)
!           freq = 0
            call error_mesg ('atmos_nudge_mod', 'no variables specified '//&
                             'for override', WARNING)
       endif
   else
       if ( do_u .or. do_v .or. do_t .or.  do_q .or. do_ps ) then
            call error_mesg ('atmos_nudge_mod', 'variables specified '//&
                             'for override when freq = 0', FATAL)
       endif
       freq = 0
   endif

 ! return flag = true when override is needed
   if (present(flag)) then
       flag = freq .gt. 0
   endif

 ! what is the next time for data insertion

   Time_next = Time + set_time(freq)

 ! initialize diagnostics

   desc = ' tendency due to data override/insertion'

   id_udt = register_diag_field ('atmos_nudge', 'udt_nudge', axes, Time, &
                                 'zonal wind'//trim(desc), 'm/s2', missing_value=missing_value)
   id_vdt = register_diag_field ('atmos_nudge', 'vdt_nudge', axes, Time, &
                                 'meridional wind'//trim(desc), 'm/s2',missing_value=missing_value)
   id_tdt = register_diag_field ('atmos_nudge', 'tdt_nudge', axes, Time, &
                                 'temperature'//trim(desc), 'degK/s',missing_value=missing_value)
   id_qdt = register_diag_field ('atmos_nudge', 'qdt_nudge', axes, Time, &
                                 'specific humidity'//trim(desc), 'kg/kg/s',missing_value=missing_value)
   id_psdt = register_diag_field ('atmos_nudge', 'psdt_nudge', axes(1:2), Time, &
                                 'surface pressure'//trim(desc), 'Pa/s',missing_value=missing_value)

   module_is_initialized = .true.

end subroutine atmos_nudge_init

!-----------------------------------------------------------------------

subroutine atmos_nudge_end

    u_tau = -1.
    v_tau = -1.
    t_tau = -1.
    q_tau = -1.
   ps_tau = -1.
   module_is_initialized = .false.

end subroutine atmos_nudge_end

!-----------------------------------------------------------------------

subroutine override_error ( Time, field )
type (time_type), intent(in) :: Time
character(len=*), intent(in) :: field
integer :: date(6)
character(len=19) :: cdate

! private routine for handling data override errors
! prints out field name and time of error

   call get_date (Time,date(1),date(2),date(3),date(4),date(5),date(6))
   write (cdate,'(i4,5(a1,i2.2))') date(1),'-',date(2),'-',date(3),' ', &
                                  date(4),':',date(5),':',date(6)
   call error_mesg ('atmos_nudge_mod', &
     'data override not done for '//trim(field)//', date = '//cdate, FATAL)

end subroutine override_error

!#######################################################################

end module atmos_nudge_mod


module damping_driver_mod

!-----------------------------------------------------------------------
!
!       This module controls four functions:
!
!   (1) rayleigh friction applied to momentum fields at levels
!       1 to kbot (i.e., momentum is damped toward zero).
!
!   (2) mountain gravity wave drag module may be called
!
!   (3) Alexander-Dunkerton gravity wave drag may be called
!
!   (4) Garner topo_drag module may be called
!
!-----------------------------------------------------------------------

 use      mg_drag_mod, only:  mg_drag, mg_drag_init, mg_drag_end, &
                              mg_drag_restart
 use      cg_drag_mod, only:  cg_drag_init, cg_drag_calc, cg_drag_end, &
                              cg_drag_time_vary, cg_drag_endts, &
                              cg_drag_restart
 use    topo_drag_mod, only:  topo_drag_init, topo_drag, topo_drag_end, &
                              topo_drag_restart
 use          mpp_mod, only:  input_nml_file
 use          fms_mod, only:  file_exist, mpp_pe, mpp_root_pe, stdlog, &
                              write_version_number, &
                              open_namelist_file, error_mesg, &
                              check_nml_error,                   &
                              FATAL, close_file
 use diag_manager_mod, only:  register_diag_field,  &
                              register_static_field, send_data
 use time_manager_mod, only:  time_type
 use    constants_mod, only:  cp_air, grav

 implicit none
 private

 public   damping_driver, damping_driver_init, damping_driver_end
 public   damping_driver_time_vary, damping_driver_endts
 public   damping_driver_restart

!-----------------------------------------------------------------------
!---------------------- namelist ---------------------------------------

   real     :: trayfric = 0.
   integer  :: nlev_rayfric = 1
   logical  :: do_mg_drag = .false.
   logical  :: do_cg_drag = .false.
   logical  :: do_topo_drag = .false.
   logical  :: do_conserve_energy = .false.

   namelist /damping_driver_nml/  trayfric,   nlev_rayfric,  &
                                  do_cg_drag, do_topo_drag, &
                                  do_mg_drag, do_conserve_energy

!
!   trayfric = damping time in seconds for rayleigh damping momentum
!              in the top nlev_rayfric layers (if trayfric < 0 then time
!              in days)
!                 [real, default: trayfric=0.]
!
!   nlev_rayfric = number of levels at the top of the model where
!                  rayleigh friction of momentum is performed, if
!                  trayfric=0. then nlev_rayfric has no effect
!                    [integer, default: nlev_rayfric=1]
!
!-----------------------------------------------------------------------
!----- id numbers for diagnostic fields -----

integer :: id_udt_rdamp,  id_vdt_rdamp,   &
           id_udt_gwd,    id_vdt_gwd,     &
                          id_sgsmtn,      &
           id_udt_cgwd,   id_taus
integer    id_vdt_cgwd

integer :: id_tdt_diss_rdamp,  id_diss_heat_rdamp, &
           id_tdt_diss_gwd,    id_diss_heat_gwd,   &
           id_tdt_diss_topo,   id_diss_heat_topo

integer :: id_udt_topo,   id_vdt_topo,   id_taubx,  id_tauby

!----- missing value for all fields ------

real :: missing_value = -999.

character(len=7) :: mod_name = 'damping'

!-----------------------------------------------------------------------

 logical :: do_rayleigh

 real, parameter ::  daypsec=1./86400.
 logical :: module_is_initialized =.false.

 real :: rfactr

!   note:  
!     rfactr = coeff. for damping momentum at the top level

 character(len=128) :: version = '$Id: damping_driver.F90,v 19.0 2012/01/06 20:05:04 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine damping_driver (is, js, lat, Time, delt, pfull, phalf, zfull, zhalf, &
                            u, v, t, q, r,  udt, vdt, tdt, qdt, rdt,  &
!                                   mask, kbot)
                            z_pbl,  mask, kbot)
 
!-----------------------------------------------------------------------
 integer,         intent(in)                :: is, js
 real, dimension(:,:), intent(in)           :: lat
 type(time_type), intent(in)                :: Time
 real,            intent(in)                :: delt
 real,    intent(in),    dimension(:,:,:)   :: pfull, phalf, &
                                               zfull, zhalf, &
                                               u, v, t, q
 real,    intent(in),    dimension(:,:,:,:) :: r
 real,    intent(inout), dimension(:,:,:)   :: udt,vdt,tdt,qdt
 real,    intent(inout), dimension(:,:,:,:) :: rdt
 real, dimension(:,:), intent(in)           :: z_pbl
 real,    intent(in),    dimension(:,:,:), optional :: mask
 integer, intent(in),    dimension(:,:),   optional :: kbot

!-----------------------------------------------------------------------
 real, dimension(size(udt,1),size(udt,2))             :: diag2
 real, dimension(size(udt,1),size(udt,2))             :: taubx, tauby
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: taus
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: utnd, vtnd, &
                                                         ttnd, pmass, &
                                                         p2
 integer :: k
 logical :: used

!-----------------------------------------------------------------------

   if (.not.module_is_initialized) call error_mesg ('damping_driver',  &
                     'damping_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------- r a y l e i g h   d a m p i n g ---------------------
!-----------------------------------------------------------------------
   if (do_rayleigh) then

       p2 = pfull * pfull
       call rayleigh (delt, p2, u, v, utnd, vtnd, ttnd)
       udt = udt + utnd
       vdt = vdt + vtnd
       tdt = tdt + ttnd

!----- diagnostics -----

       if ( id_udt_rdamp > 0 ) then
            used = send_data ( id_udt_rdamp, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_rdamp > 0 ) then
            used = send_data ( id_vdt_rdamp, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_tdt_diss_rdamp > 0 ) then
            used = send_data ( id_tdt_diss_rdamp, ttnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_diss_heat_rdamp > 0 ) then
            do k = 1,size(u,3)
              pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
            enddo
            diag2 = cp_air/grav * sum(ttnd*pmass,3)
            used = send_data ( id_diss_heat_rdamp, diag2, Time, is, js )
       endif

   endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------- m t n   g r a v i t y   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
   if (do_mg_drag) then

       call mg_drag (is, js, delt, u, v, t, pfull, phalf, zfull, zhalf,  &
                     utnd, vtnd, ttnd, taubx,tauby,taus,        kbot)
       udt = udt + utnd
       vdt = vdt + vtnd
       tdt = tdt + ttnd

!----- diagnostics -----

       if ( id_udt_gwd > 0 ) then
            used = send_data ( id_udt_gwd, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_gwd > 0 ) then
            used = send_data ( id_vdt_gwd, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_taubx > 0 ) then
            used = send_data ( id_taubx, taubx, Time, is, js )
       endif

       if ( id_tauby > 0 ) then
            used = send_data ( id_tauby, tauby, Time, is, js )
       endif

       if ( id_taus > 0 ) then
           used = send_data ( id_taus, taus, Time, is, js, 1, &
                              rmask=mask )
       endif

       if ( id_tdt_diss_gwd > 0 ) then
            used = send_data ( id_tdt_diss_gwd, ttnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_diss_heat_gwd > 0 ) then
            do k = 1,size(u,3)
              pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
            enddo
            diag2 = cp_air/grav * sum(ttnd*pmass,3)
            used = send_data ( id_diss_heat_gwd, diag2, Time, is, js )
       endif

   endif

!   Alexander-Dunkerton gravity wave drag

   if (do_cg_drag) then

     call cg_drag_calc (is, js, lat, pfull, zfull, t, u, v, Time,    &
                        delt, utnd, vtnd)

     udt =  udt + utnd
     vdt =  vdt + vtnd

!----- diagnostics -----

     if ( id_udt_cgwd > 0 ) then
        used = send_data ( id_udt_cgwd, utnd, Time, is, js, 1, &
                          rmask=mask )
     endif
      if ( id_vdt_cgwd > 0 ) then
        used = send_data ( id_vdt_cgwd, vtnd, Time, is, js, 1, &
                          rmask=mask )
     endif


   endif

!-----------------------------------------------------------------------
!---------topographic   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
   if (do_topo_drag) then

     call topo_drag ( is, js, delt, u, v, t, pfull, phalf, zfull, zhalf,  &
                      utnd, vtnd, ttnd, taubx, tauby, taus, kbot )

     udt = udt + utnd
     vdt = vdt + vtnd

!----- diagnostics -----

     if ( id_udt_topo > 0 ) then
        used = send_data ( id_udt_topo, utnd, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_vdt_topo > 0 ) then
        used = send_data ( id_vdt_topo, vtnd, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_taubx > 0 ) then
       used = send_data ( id_taubx, taubx, Time, is, js )
     endif

     if ( id_tauby > 0 ) then
        used = send_data ( id_tauby, tauby, Time, is, js )
     endif

     if ( id_taus > 0 ) then
        used = send_data ( id_taus, taus, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_tdt_diss_topo > 0 ) then
        used = send_data ( id_tdt_diss_topo, ttnd, Time, is, js, 1, &
                               rmask=mask )
     endif

     if ( id_diss_heat_topo > 0 ) then
          do k = 1,size(u,3)
             pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
          enddo
          diag2 = cp_air/grav * sum(ttnd*pmass,3)
          used = send_data ( id_diss_heat_topo, diag2, Time, is, js )
     endif

 endif

!-----------------------------------------------------------------------

 end subroutine damping_driver

!#######################################################################

 subroutine damping_driver_init ( lonb, latb, pref, axes, Time, sgsmtn)

 real,            intent(in) :: lonb(:,:), latb(:,:), pref(:)
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 real, dimension(:,:), intent(out) :: sgsmtn
!-----------------------------------------------------------------------
!     lonb  = longitude in radians of the grid box corners
!     latb  = latitude  in radians of the grid box corners
!     axes  = axis indices, (/x,y,pf,ph/)
!               (returned from diag axis manager)
!     Time  = current time (time_type)
!     sgsmtn = subgrid scale topography variance
!-----------------------------------------------------------------------
 integer :: unit, ierr, io, logunit
 logical :: used

!-----------------------------------------------------------------------
!----------------- namelist (read & write) -----------------------------

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=damping_driver_nml, iostat=io)
   ierr = check_nml_error(io,"damping_driver_nml")
#else
   if (file_exist('input.nml')) then
      unit = open_namelist_file ()
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=damping_driver_nml, iostat=io, end=10)
         ierr = check_nml_error (io, 'damping_driver_nml')
      enddo
 10   call close_file (unit)
   endif
#endif

   call write_version_number(version, tagname)
   logunit = stdlog()
   if(mpp_pe() == mpp_root_pe() ) then
        write (logunit,nml=damping_driver_nml)
   endif

!-----------------------------------------------------------------------
!--------- rayleigh friction ----------

   do_rayleigh=.false.

   if (abs(trayfric) > 0.0001 .and. nlev_rayfric > 0) then
      if (trayfric > 0.0) then
         rfactr=(1./trayfric)
      else
         rfactr=(1./abs(trayfric))*daypsec
      endif
         do_rayleigh=.true.
   else
         rfactr=0.0
   endif

!-----------------------------------------------------------------------
!----- mountain gravity wave drag -----

   if (do_mg_drag) call mg_drag_init (lonb, latb, sgsmtn)

!--------------------------------------------------------------------
!----- Alexander-Dunkerton gravity wave drag -----
 
   if (do_cg_drag)  then
     call cg_drag_init (lonb, latb, pref, Time=Time, axes=axes)
   endif

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

if (do_rayleigh) then

   id_udt_rdamp = &
   register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time,       &
                       'u wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_rdamp = &
   register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time,       &
                       'v wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )

   id_tdt_diss_rdamp = &
   register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), Time,  &
                      'Dissipative heating from Rayleigh damping',&
                             'deg_k/s', missing_value=missing_value   )
       
   id_diss_heat_rdamp = &
   register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), Time,   &
                'Integrated dissipative heating from Rayleigh damping',&
                  'W/m2' )
endif

if (do_mg_drag) then

 ! register and send static field
   id_sgsmtn = &
   register_static_field ( mod_name, 'sgsmtn', axes(1:2), &
               'sub-grid scale topography for gravity wave drag', 'm')
   if (id_sgsmtn > 0) used = send_data (id_sgsmtn, sgsmtn, Time)

 ! register non-static field
   id_udt_gwd = &
   register_diag_field ( mod_name, 'udt_gwd', axes(1:3), Time,        &
                     'u wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_gwd = &
   register_diag_field ( mod_name, 'vdt_gwd', axes(1:3), Time,        &
                     'v wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_taubx = &
   register_diag_field ( mod_name, 'taubx', axes(1:2), Time,        &
                         'x base flux for grav wave drag', 'kg/m/s2', &
                         missing_value=missing_value               )

   id_tauby = &
   register_diag_field ( mod_name, 'tauby', axes(1:2), Time,        &
                         'y base flux for grav wave drag', 'kg/m/s2', &
                         missing_value=missing_value )

   id_taus = &
   register_diag_field ( mod_name, 'taus', axes(1:3), Time,        &
                       'saturation flux for gravity wave drag', 'kg/m/s2', &
                      missing_value=missing_value               )

   id_tdt_diss_gwd = &
   register_diag_field ( mod_name, 'tdt_diss_gwd', axes(1:3), Time,    &
                          'Dissipative heating from gravity wave drag',&
                              'deg_k/s', missing_value=missing_value   )
       
   id_diss_heat_gwd = &
   register_diag_field ( mod_name, 'diss_heat_gwd', axes(1:2), Time,      &
                'Integrated dissipative heating from gravity wave drag',&
                                 'W/m2' )
endif

   if (do_cg_drag) then

    id_udt_cgwd = &
    register_diag_field ( mod_name, 'udt_cgwd', axes(1:3), Time,        &
                 'u wind tendency for cg gravity wave drag', 'm/s2', &
                      missing_value=missing_value               )


    id_vdt_cgwd = &
    register_diag_field ( mod_name, 'vdt_cgwd', axes(1:3), Time,        &
                 'v wind tendency for cg gravity wave drag', 'm/s2', &
                      missing_value=missing_value               )


   endif

!-----------------------------------------------------------------------
!----- topo wave drag -----



  if (do_topo_drag) then
          call topo_drag_init (lonb, latb)
          sgsmtn(:,:) = -99999.
  endif



  if (do_topo_drag) then

   id_udt_topo = &
   register_diag_field ( mod_name, 'udt_topo', axes(1:3), Time,        &
                       'u wind tendency for topo wave drag', 'm/s2', &
                        missing_value=missing_value               )

  id_vdt_topo = &
   register_diag_field ( mod_name, 'vdt_topo', axes(1:3), Time,        &
                       'v wind tendency for topo wave drag', 'm/s2', &
                         missing_value=missing_value               )

   id_taubx = &
   register_diag_field ( mod_name, 'taubx', axes(1:2), Time,        &
                     'x base flux for topo wave drag', 'kg/m/s2', &
                        missing_value=missing_value               )

    id_tauby = &
    register_diag_field ( mod_name, 'tauby', axes(1:2), Time,        &
                    'y base flux for topo wave drag', 'kg/m/s2', &
                      missing_value=missing_value )

    id_taus = &
   register_diag_field ( mod_name, 'taus', axes(1:3), Time,        &
                  'saturation flux for topo wave drag', 'kg/m/s2', &
                     missing_value=missing_value               )

   id_tdt_diss_topo = &
   register_diag_field ( mod_name, 'tdt_diss_topo', axes(1:3), Time,    &
                          'Dissipative heating from topo wave drag',&
                              'deg_k/s', missing_value=missing_value   )
       
   id_diss_heat_topo = &
   register_diag_field ( mod_name, 'diss_heat_topo', axes(1:2), Time,      &
                'Integrated dissipative heating from topo wave drag',&
                                 'W/m2' )
 endif


!-----------------------------------------------------------------------

   module_is_initialized =.true.

!******************** end of initialization ****************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 end subroutine damping_driver_init


!#####################################################################

subroutine damping_driver_time_vary (delt)

real, intent(in) :: delt


       call cg_drag_time_vary (delt)

end subroutine damping_driver_time_vary



!#####################################################################

subroutine damping_driver_endts


     call cg_drag_endts

end subroutine damping_driver_endts



!######################################################################
!#######################################################################

 subroutine damping_driver_end

     if (do_mg_drag)   call mg_drag_end
     if (do_cg_drag)   call cg_drag_end
     if (do_topo_drag) call topo_drag_end

     module_is_initialized =.false.


 end subroutine damping_driver_end

!#######################################################################
! <SUBROUTINE NAME="damping_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine damping_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

     if (do_mg_drag)   call mg_drag_restart(timestamp)
     if (do_cg_drag)   call cg_drag_restart(timestamp)
     if (do_topo_drag) call topo_drag_restart(timestamp)

end subroutine damping_driver_restart
! </SUBROUTINE> NAME="damping_driver_restart" 

!#######################################################################

 subroutine rayleigh (dt, p2, u, v, udt, vdt, tdt)

  real,    intent(in)                      :: dt
  real,    intent(in),  dimension(:,:,:)   :: p2, u, v
  real,    intent(out), dimension(:,:,:)   :: udt, vdt, tdt

  real, dimension(size(u,1),size(u,2)) :: fact
  integer :: k
!-----------------------------------------------------------------------
!--------------rayleigh damping of momentum (to zero)-------------------

   do k = 1, nlev_rayfric
     fact(:,:) = rfactr*(1.+(p2(:,:,1)-p2(:,:,k))/(p2(:,:,1)+p2(:,:,k)))
     udt(:,:,k) = -u(:,:,k)*fact(:,:)
     vdt(:,:,k) = -v(:,:,k)*fact(:,:)
   enddo

   do k = nlev_rayfric+1, size(u,3)
     udt(:,:,k) = 0.0
     vdt(:,:,k) = 0.0
   enddo

!  total energy conservation
!  compute temperature change loss due to ke dissipation

   if (do_conserve_energy) then
       do k = 1, nlev_rayfric
          tdt(:,:,k) = -((u(:,:,k)+.5*dt*udt(:,:,k))*udt(:,:,k) +  &
                         (v(:,:,k)+.5*dt*vdt(:,:,k))*vdt(:,:,k)) / cp_air
       enddo
       do k = nlev_rayfric+1, size(u,3)
          tdt(:,:,k) = 0.0
       enddo
   else
       tdt = 0.0
   endif

!-----------------------------------------------------------------------

 end subroutine rayleigh

!#######################################################################

end module damping_driver_mod


module ocean_momentum_source_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Internal momentum sources, such as sink from Rayleigh drag.   
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted and density weighted tendency
! [velocity*(kg/m^3)*meter/sec] for velocity associated with 
! momentun sources in the interior.  Primary application is to 
! specify Rayleigh drag (a momentum sink).  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, Elements of MOM (2012)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_momentum_source_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Needs to be true in order to use this scheme. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose initialization information.
!  </DATA> 
!  <DATA NAME="use_rayleigh_damp_table" TYPE="logical">
!  For reading in Rayleigh damping times from a table.
!  </DATA> 
!  <DATA NAME="rayleigh_damp_exp_from_bottom" TYPE="logical">
!  For computing a Rayleigh damping time with largest at bottom and 
!  decaying towards surface.
!  </DATA> 
!  <DATA NAME="rayleigh_damp_exp_scale" TYPE="real"  UNITS="metre">
!  Exponential decay scale from bottom upwards for computing 
!  the Rayleigh damping time. 
!  </DATA> 
!  <DATA NAME="rayleigh_damp_exp_max" TYPE="real"  UNITS="sec">
!  Damping time at the bottom for rayleigh_damp_exp_from_bottom=.true.
!  </DATA> 
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field, register_static_field
use field_manager_mod, only: MODEL_OCEAN, parse, find_field_index
use field_manager_mod, only: get_field_methods, method_type, get_field_info
use fms_mod,           only: stdout, stdlog, FATAL, NOTE, WARNING
use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
use mpp_domains_mod,   only: mpp_update_domains 
use mpp_mod,           only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_parameters_mod, only: missing_value, sec_in_day_r
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_time_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_options_type
use ocean_workspace_mod,  only: wrk1_v, wrk2
use ocean_util_mod,       only: diagnose_3d_u

implicit none

private 

public ocean_momentum_source_init
public momentum_source

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

real, dimension(:,:,:), allocatable  :: rayleigh_damp   ! inverse Rayleigh damping time (1/sec)
                                                        ! set to zero for places with no damping  

! for reading in points where wish to set a static Rayleigh drag
integer, dimension(:),   allocatable :: itable
integer, dimension(:),   allocatable :: jtable     
integer, dimension(:,:), allocatable :: ktable     
real,    dimension(:),   allocatable :: rayleigh_damp_table

! for diagnostics 
logical :: used
integer :: id_rayleigh_drag_u    =-1
integer :: id_rayleigh_drag_v    =-1
integer :: id_rayleigh_damp      =-1
integer :: id_rayleigh_damp_table=-1
integer :: id_rayleigh_drag_power=-1

character(len=128) :: version=&
       '$Id: ocean_momentum_source.F90,v 20.0 2013/12/14 00:16:02 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk 

real    :: rayleigh_damp_exp_scale_r
real    :: rayleigh_damp_exp_time_r

logical :: use_rayleigh_damp_table       =.false.                     
logical :: rayleigh_damp_exp_from_bottom =.false.
real    :: rayleigh_damp_exp_scale       =100.0
real    :: rayleigh_damp_exp_time        =8.64e5

logical :: verbose_init                  = .true.
logical :: use_this_module               = .false.
logical :: debug_this_module             = .false.
logical :: module_is_initialized         = .FALSE.

namelist /ocean_momentum_source_nml/ verbose_init, use_this_module, debug_this_module,       &
                                     use_rayleigh_damp_table, rayleigh_damp_exp_from_bottom, &
                                     rayleigh_damp_exp_scale, rayleigh_damp_exp_time

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_momentum_source_init">
!
! <DESCRIPTION>
!    Initial set up for momentum source/sink 
! </DESCRIPTION>
!
subroutine ocean_momentum_source_init(Grid, Domain, Time, Ocean_options, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  logical,                      intent(in), optional :: debug

  integer :: io_status, ioun, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if (module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_momentum_source_mod (ocean_momentum_source_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_momentum_source_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_momentum_source_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_momentum_source_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_momentum_source_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_momentum_source_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_momentum_source_nml)

  if(use_this_module) then 
      write(stdoutunit,'(a)') '==>Note from ocean_momentum_source_mod: USING this module' 
      Ocean_options%momentum_source = 'Used momentum sources.'
  else 
      call mpp_error(NOTE, '==>Note from ocean_momentum_source_mod: NOT USING this module')
      Ocean_options%momentum_source = 'Did NOT use momentum sources.'
      return
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
  endif 
  if(debug_this_module) then 
     write(stdoutunit,'(a)') '==>Note: running ocean_momentum_source_mod with debug_this_module=.true.'  
  endif 

  allocate(rayleigh_damp(isd:ied,jsd:jed,nk))
  rayleigh_damp(:,:,:) = 0.0

  rayleigh_damp_exp_scale_r = 1.0/(epsln+rayleigh_damp_exp_scale)
  rayleigh_damp_exp_time_r  = 1.0/(epsln+rayleigh_damp_exp_time)

  ! read in static rayleigh drag damping times from a table
  call rayleigh_damp_table_init(Time)


end subroutine ocean_momentum_source_init
! </SUBROUTINE> NAME="ocean_momentum_source_init"


!#######################################################################
! <SUBROUTINE NAME="momentum_source">
!
! <DESCRIPTION>
!    Compute the momentum source/sink (N/m2).
! </DESCRIPTION>
!
subroutine momentum_source(Time, Thickness, Velocity, energy_analysis_step)

  type(ocean_time_type),       intent(in)    :: Time
  type(ocean_thickness_type),  intent(in)    :: Thickness
  type(ocean_velocity_type),   intent(inout) :: Velocity
  logical,                     intent(in)    :: energy_analysis_step 

  integer :: i,j,k,n
  integer :: tau, taum1
  integer :: kmu
  real    :: argument, rayleigh_damp_scalar

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
    '==>Error in ocean_momentum_source_mod (momentum_source): module needs to be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  wrk1_v(:,:,:,:) = 0.0

  ! Inverse Rayleigh damping times, with largest damping near the bottom. 
  if(rayleigh_damp_exp_from_bottom) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               kmu = Grd%kmu(i,j)
               if(k <= kmu) then 
                   argument = -(Thickness%depth_zu(i,j,kmu)-Thickness%depth_zu(i,j,k))*rayleigh_damp_exp_scale_r 
                   rayleigh_damp_scalar = rayleigh_damp_exp_time_r*exp(argument)
                   rayleigh_damp(i,j,k) = Grd%umask(i,j,k)*max(rayleigh_damp(i,j,k), rayleigh_damp_scalar)
               endif
            enddo
         enddo
      enddo
  endif

  ! drag from Rayleigh damping 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              wrk1_v(i,j,k,n) = -Thickness%rho_dzu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1)*rayleigh_damp(i,j,k)
           enddo
        enddo
     enddo
  enddo

  ! energetic analysis 
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

  ! add to acceleration 
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
      call diagnose_3d_u(Time, Grd, id_rayleigh_drag_u, wrk1_v(:,:,:,1))
      call diagnose_3d_u(Time, Grd, id_rayleigh_drag_v, wrk1_v(:,:,:,2))
      if (id_rayleigh_drag_power > 0) then 
          wrk2(:,:,:) = 0.0
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      wrk2(i,j,k) = wrk2(i,j,k) + Velocity%u(i,j,k,n,tau)*wrk1_v(i,j,k,n)*Grd%dau(i,j)
                   enddo
                enddo
             enddo
          enddo
          call diagnose_3d_u(Time, Grd, id_rayleigh_drag_power, wrk2(:,:,:))
      endif

  endif 


end subroutine momentum_source
! </SUBROUTINE> NAME="momentum_source"


!#######################################################################
! <SUBROUTINE NAME="rayleigh_damp_table_init">
!
! <DESCRIPTION>
! Read in static Rayleigh drag dissipation times entered to the 
! table "rayleigh_damp_table".  The dissipation times should be 
! entered in units of seconds.  
! </DESCRIPTION>
!
subroutine rayleigh_damp_table_init(Time)

  type(ocean_time_type), intent(in) :: Time
  type(method_type), allocatable, dimension(:) :: rayleigh_damp_methods

  character(len=32) :: fld_type, fld_name
  integer :: ntable, ntable_points, model, parse_ok
  integer :: itbl, jtbl, ntbl
  integer :: k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ntable = find_field_index(MODEL_OCEAN,'rayleigh_damp_table')
  if (ntable < 1) then 
    call mpp_error(NOTE, &
    '==>Warning: ocean_momentum_source_init NO table for Rayleigh damp time. Nothing read.')  
  endif
  if(ntable > 1 .and. .not. use_rayleigh_damp_table) then 
    call mpp_error(NOTE, &
    '==>Warning: ocean_momentum_source_init found rayleigh_damp_table, yet use_rayleigh_damp_table=.false.')  
  endif 
  if(ntable > 1 .and. use_rayleigh_damp_table) then 
    write(stdoutunit,*) &
    '==>Note: ocean_momentum_source_init will read Rayleigh damping coefficients from rayleigh_damp_table.'  
  endif 

  ! read in Rayleigh damping times specified from "rayleigh_damp_table"
  wrk2(:,:,:) = 0.0 
  if(use_rayleigh_damp_table .and. ntable > 1) then 

      call get_field_info(ntable,fld_type,fld_name,model,ntable_points)
      allocate(rayleigh_damp_methods(ntable_points))

      allocate (itable(ntable_points))
      itable(:) = 0
      allocate (jtable(ntable_points))
      jtable(:) = 0
      allocate (ktable(ntable_points,2))
      ktable(:,:) = 0
      allocate (rayleigh_damp_table(ntable_points))
      rayleigh_damp_table(:) = 0.0

      call get_field_methods(ntable,rayleigh_damp_methods)
      do ntbl=1,ntable_points
         parse_ok = parse(rayleigh_damp_methods(ntbl)%method_control,'itable',itable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_momentum_source_init: rayleigh_damp_table entry "itable" error')
         parse_ok = parse(rayleigh_damp_methods(ntbl)%method_control,'jtable',jtable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_momentum_source_init: rayleigh_damp_table entry "jtable" error')
         parse_ok = parse(rayleigh_damp_methods(ntbl)%method_control,'ktable_1',ktable(ntbl,1))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_momentum_source_init: rayleigh_damp_table entry "ktable_1" error')
         parse_ok = parse(rayleigh_damp_methods(ntbl)%method_control,'ktable_2',ktable(ntbl,2))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_momentum_source_init: rayleigh_damp_table entry "ktable_2" error')
         parse_ok = parse(rayleigh_damp_methods(ntbl)%method_control,'rayleigh_damp_table',rayleigh_damp_table(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,&
         '==>Error ocean_momentum_source_init: rayleigh_damp_table entry "rayleigh_damp_table" error')
         if(ktable(ntbl,2) < ktable(ntbl,1) ) then 
             call mpp_error(FATAL,&
             '==>Error ocean_momentum_source_init: rayleigh_damp_table entry needs ktable_2 >= ktable_1')
         endif
      enddo

      do ntbl=1,ntable_points

         if(on_comp_domain(ntbl)) then

             itbl=itable(ntbl)-Dom%ioff
             jtbl=jtable(ntbl)-Dom%joff

             do k=ktable(ntbl,1),ktable(ntbl,2)

                if(Grd%umask(itbl,jtbl,k) == 0.0) then 
                    write(*,'(a,i4,a,i4,a,i4,a)')                                                  &
                    '==>Warning ocean_momentum_source_init: ignored nonzero rayleigh_damp_table(', &
                    itable(ntbl),',',jtable(ntbl),',',k,') set over land. Is your table correct?'
                endif

                if(rayleigh_damp_table(ntbl) > 0.0) then 
                   wrk2(itbl,jtbl,k) = Grd%umask(itbl,jtbl,k)/rayleigh_damp_table(ntbl)
                else 
                   wrk2(itbl,jtbl,k) = 0.0
                endif  
                rayleigh_damp(itbl,jtbl,k) = rayleigh_damp(itbl,jtbl,k) + wrk2(itbl,jtbl,k)

             enddo

         endif

      enddo

  endif

  id_rayleigh_damp_table = register_static_field ('ocean_model', 'rayleigh_damp_table',            &
                           Grd%vel_axes_uv(1:3), 'Inverse Rayleigh damping time from table', '1/s',&
                           missing_value=missing_value, range=(/-1.0,1e6/))

  id_rayleigh_damp = register_static_field ('ocean_model', 'rayleigh_damp',       &
                     Grd%vel_axes_uv(1:3), 'Inverse Rayleigh damping time', '1/s',&
                     missing_value=missing_value, range=(/-1.0,1e6/))

  call diagnose_3d_u(Time, Grd, id_rayleigh_damp_table, wrk2(:,:,:))
  call diagnose_3d_u(TIme, Grd, id_rayleigh_damp, rayleigh_damp(:,:,:))

  id_rayleigh_drag_u = register_diag_field ('ocean_model', 'rayleigh_drag_u',  &
                       Grd%vel_axes_uv(1:3), Time%model_time,                  &
                       'Rayleigh drag on i-component of velocity', 'N/m2',     &
                       missing_value=missing_value, range=(/-1e6,1e6/))

  id_rayleigh_drag_v = register_diag_field ('ocean_model', 'rayleigh_drag_v',  &
                       Grd%vel_axes_uv(1:3), Time%model_time,                  &
                       'Rayleigh drag on j-component of velocity', 'N/m2',     &
                       missing_value=missing_value, range=(/-1e6,1e6/))

  id_rayleigh_drag_power = register_diag_field ('ocean_model', 'rayleigh_drag_power',  &
                       Grd%vel_axes_uv(1:3), Time%model_time,                          &
                       'Energy sink from Rayleigh drag', 'Watt',                       &
                       missing_value=missing_value, range=(/-1e16,1e16/))


end subroutine rayleigh_damp_table_init
! </SUBROUTINE> NAME="rayleigh_damp_table_init"


!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor.
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular point in a table.
! </IN>
function on_comp_domain(ntable)

  integer, intent(in) :: ntable
  logical             :: on_comp_domain

  if(isc+Dom%ioff <= itable(ntable) .and. itable(ntable) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jtable(ntable) .and. jtable(ntable) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"



end module ocean_momentum_source_mod
      
      





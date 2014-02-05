
module bgrid_core_driver_mod

!-----------------------------------------------------------------------
!
!      Driver module for running the FMS B-grid dynamical core.
!
!        * reads namelist
!        * sets up the B-grid core
!        * packages the B-grid core with diagnostic routines
!
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
use bgrid_core_mod           , only: bgrid_dynam_type,  &
                                     bgrid_core_init, update_bgrid_core,  &
                                     bgrid_core_end
use bgrid_horiz_mod          , only: horiz_grid_type, horiz_grid_init
use bgrid_vert_mod           , only: vert_grid_type, vert_grid_init
use bgrid_prog_var_mod       , only: prog_var_type, prog_var_init, &
                                     prog_var_time_diff, var_init, &
                                     open_prog_var_file,           &
                                     read_prog_var, write_prog_var
use bgrid_diagnostics_mod    , only: bgrid_diagnostics,      &
                                     bgrid_diagnostics_tend, &
                                     bgrid_diagnostics_init
use bgrid_integrals_mod      , only: bgrid_integrals, bgrid_integrals_init, &
                                     bgrid_integrals_end
use bgrid_conserve_energy_mod, only: bgrid_conserve_energy_init, &
                                     bgrid_conserve_energy,      &
                                     bgrid_conserve_energy_end

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers
use   time_manager_mod, only: time_type, get_time
use            fms_mod, only: error_mesg, FATAL, file_exist, open_namelist_file,  &
                              check_nml_error, write_version_number,     &
                              mpp_pe, mpp_root_pe, close_file, stdlog

use    mpp_domains_mod, only: domain2d
!-----------------------------------------------------------------------

implicit none
private

public  bgrid_dynam_type, bgrid_core_driver_init, &
        bgrid_core_driver, bgrid_core_driver_end, &
        bgrid_core_time_diff, get_bottom_data, put_bottom_data, &
        atmosphere_domain

!-----------------------------------------------------------------------
character(len=128) :: version =  '$Id: bgrid_core_driver.F90,v 19.0 2012/01/06 19:53:59 fms Exp $'
character(len=128) :: tag =  '$Name: tikal $'
!-----------------------------------------------------------------------
!
!             NAMELIST INPUT: bgrid_core_driver_nml
!
!           This namelist is read from file input.nml.
!           See the on-line documentation for more details.
!
!-----------------------------------------------------------------------
!
!  num_adjust_dt       The number of adjustment time steps for each advection
!                      time step, where num_adjust_dt >= 1.
!
!  num_advec_dt        The number of advection time steps for each
!                      atmospheric/physics time step, where num_advec_dt >= 1.

   integer :: num_adjust_dt = 3
   integer :: num_advec_dt  = 3

!  layout              The domain decomposition, where layout(1) = x-axis
!                      decomposition, layout(2) = y-axis decomposition.
!                      * If layout(1)*layout(2) does not equal the number
!                        of processors the model will fail.
!                      * If layout(1)=layout(2)=0 then the decomposition is
!                        determined by MPP_DEFINE_LAYOUT.

   integer, dimension(2) :: layout = (/0,0/)

!  filter_option       Determines how polar filtering is performed.
!                      Possible values are :
!                        filter_option = 0, no polar filtering (decrease time step)
!                        filter_option = 1, obsolete scheme (NO NOT USE)
!                        filter_option = 2, default scheme (refer to technical doc)
! 
!  filter_weight       Weight applied to the polar filter that will
!                      increase (or decrease) the strength of the standard
!                      polar filter response function.
!
!  ref_lat_filter      The reference latitude at which polar filtering
!                      (in each hemisphere) will begin to be applied.

   integer   ::     filter_option = 2
   integer   ::     filter_weight = 1
   real      ::     ref_lat_filter = 60.

!  do_conserve_energy  If TRUE the temperature tendency will be updated to
!                      guarantee that the dynamical core conserves total energy.
!                      The correction is applied as a global uniform value.

   logical   :: do_conserve_energy = .false.

!  pgf_scheme          The scheme used to compute the pressure gradient.
!                      Specify one of the following: 'default', 'finite_volume'.
!                      The default scheme is that of Simmons and Burridge.

   character(len=24) :: pgf_scheme = 'default'   ! default, finite_volume

!  restart_output_format   Format used for the output restart file.
!                          The only possible values are: 'native' or 'netcdf'.

   character(len=24) :: restart_output_format = 'netcdf'   ! native, netcdf

!  do_average_omega    If TRUE the omega diagostic returned by the dynamical core
!                      is averaged over all adjustment time steps. If FALSE then
!                      omega for the last adjustment step is returned.

   logical   :: do_average_omega = .false.

!  ddamp_coeff         Damping coefficient for divergence damping.

   real :: ddamp_coeff = 0.0

!  verbose             Flag that control additional printed output.
!                      Currently, this option is not being used.

   integer :: verbose = 0


   namelist /bgrid_core_driver_nml/  num_adjust_dt, num_advec_dt,   &
                                     layout, filter_option,         & 
                                     filter_weight, ref_lat_filter, &
                                     do_conserve_energy,            &
                                     pgf_scheme,                    &
                                     restart_output_format,         &
                                     do_average_omega,              &
                                     ddamp_coeff, verbose

!-----------------------------------------------------------------------
!------ private data ------

real, dimension(:,:),   pointer :: fis, res        ! topography data
real, dimension(:,:,:), pointer :: div, mfew, mfns ! diagnostic fields
real, dimension(:), allocatable :: eta, peta       ! vertical grid


! derived type data containing horizonal and vertical grid constants
! other pointer data may "point at" these data at any time
type  (horiz_grid_type), target, save :: Hgrid
type   (vert_grid_type), target, save :: Vgrid

! axis indices returned by diagnostics manager
integer, dimension(4) :: mass_axes, vel_axes

real :: dt_atmos ! atmospheric time step

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine bgrid_core_driver_init ( Time_init, Time, Time_step, &
                                     Var, Var_dt, Dynam, phys_axes )

!-----------------------------------------------------------------------
! Time_init = initial time
! Time      = current time
! Time_step = atmospheric model time step
! Var       = prognostic variables
! Var_dt    = prognostic variable tendencies
! Dynam     = data type for dynamical core constants and data
! phys_axes = axis indices for the grid used by the atmospheric physics
!-----------------------------------------------------------------------

 type       (time_type), intent(in)    :: Time_init, Time, Time_step
 type   (prog_var_type), intent(inout) :: Var, Var_dt
 type(bgrid_dynam_type), intent(inout) :: Dynam
 integer,                intent(out)   :: phys_axes(4)

 integer :: unit, io, ierr, logunit
 integer :: ix, jx, kx
 integer :: sec, ntrace, ntprog, ntdiag
!-----------------------------------------------------------------------
!  ----- read namelist -----

    if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=bgrid_core_driver_nml, iostat=io)
        ierr = check_nml_error(io,'bgrid_core_driver_nml')
#else
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=bgrid_core_driver_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'bgrid_core_driver_nml')
        enddo
 10     call close_file (unit)
#endif
    endif
    logunit = stdlog()

!-----------------------------------------------------------------------
!  ----- read restart header records and set up grid resolution -----

   call open_prog_var_file (ix, jx, kx)

!  ---- horizontal grid initialization ----

   call horiz_grid_init ( Hgrid, ix, jx, layout=layout )

! how many tracers have been registered?
   call get_number_tracers ( MODEL_ATMOS, num_tracers=ntrace, num_prog=ntprog, num_diag=ntdiag )

!  ----- write version, namelist and tracer info to log file -----

    call write_version_number (version, tag)
    if (mpp_pe() == mpp_root_pe()) then
        write (logunit, nml=bgrid_core_driver_nml)
        write (logunit, '(a,i3)') 'Number of tracers =', ntrace
        write (logunit, '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (logunit, '(a,i3)') 'Number of diagnostic tracers =', ntdiag
    endif

!  ---- prognostic variable initialization -----

   call prog_var_init (Hgrid, kx, ntrace, Var)     ! prognostic+diagnostic tracers
   call prog_var_init (Hgrid, kx, ntprog, Var_dt)

!----- read data -----

   fis   => var_init (Hgrid)
   res   => var_init (Hgrid)

   allocate (eta(kx+1), peta(kx+1))

   call read_prog_var ( Hgrid, Var, eta, peta, fis, res )


!---- vertical grid initialization ----

   call vert_grid_init ( Vgrid, eta, peta )

   deallocate (eta, peta)

!---- diagnostic fields ----

   div  => var_init (Hgrid,kx)
   mfew => var_init (Hgrid,kx)
   mfns => var_init (Hgrid,kx)

!---- compute time step in seconds ----

   call get_time (Time_step, sec)
   dt_atmos = real(sec)

!-----------------------------------------------------------------------
!  ----- initialize dynamical core -----

   call bgrid_core_init ( Dynam, Hgrid, Vgrid, fis, res, dt_atmos,     &
                          num_adjust_dt, num_advec_dt, pgf_scheme,     &
                          filter_option, filter_weight, ref_lat_filter,&
                          ddamp_coeff, do_average_omega, verbose )

!-----------------------------------------------------------------------
!---- initialize history (netcdf) file and integrals -------

   call bgrid_diagnostics_init ( Time, Hgrid, Vgrid, Var,      &
                                 fis, res, mass_axes, vel_axes )
   phys_axes = mass_axes

   call bgrid_integrals_init (Time_init, Time)

!---- initialize integrals ----

   call bgrid_integrals (Time, Hgrid, Vgrid, Var, Dynam%Masks)

!---- initialize energy conservation module ----

   if (do_conserve_energy) call bgrid_conserve_energy_init (Time, mass_axes)

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver_init

!#######################################################################

 subroutine bgrid_core_driver (Time_diag, Var, Var_dt, Dynam, omega)

!-----------------------------------------------------------------------
! Time_diag = time used for diagnostic output
!              (typically time at end of time step)
! Var       = prognostic variables
! Var_dt    = prognostic variable tendencies
! Dynam     = data type for dynamical core constants and data
! omega     = omega (vertical velocity) diagnostic (Pa/s)
!-----------------------------------------------------------------------
   type       (time_type), intent(in)    :: Time_diag
   type   (prog_var_type), intent(in)    :: Var
   type   (prog_var_type), intent(inout) :: Var_dt
   type(bgrid_dynam_type), intent(inout) :: Dynam
   real,                   intent(out)   :: omega(:,:,:)
!-----------------------------------------------------------------------
!  dynamics

   call update_bgrid_core (Var, Var_dt, Dynam, omega, div, mfew, mfns )

!  energy conservation

   if (do_conserve_energy) then
       call bgrid_conserve_energy ( dt_atmos, Time_diag, Hgrid, Vgrid, &
                                    Dynam%Masks, Var, Var_dt )
   endif

!  diagnostics for dynamics tendencies

   call bgrid_diagnostics_tend ( Hgrid, Var_dt, Dynam%Masks, Time_diag )

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver

!#######################################################################

 subroutine bgrid_core_time_diff ( omega, Time_diag, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------
! omega     = omega (vertical velocity) diagnostic (Pa/s)
! Time_diag = time used for diagnostic output
!              (typically time at end of time step)
! Dynam     = data type for dynamical core constants and data
! Var       = prognostic variables
! Var_dt    = prognostic variable tendencies
!-----------------------------------------------------------------------
   real,                   intent(in)    :: omega(:,:,:)
   type       (time_type), intent(in)    :: Time_diag
   type(bgrid_dynam_type), intent(in)    :: Dynam
   type   (prog_var_type), intent(inout) :: Var
   type   (prog_var_type), intent(inout) :: Var_dt
!-----------------------------------------------------------------------

!  time differencing

   call prog_var_time_diff ( dt_atmos, Dynam%Masks, Var_dt, Var )

!  global integrals and diagnostics

   call bgrid_integrals ( Time_diag, Hgrid, Vgrid, Var, Dynam%Masks )

   call bgrid_diagnostics ( Hgrid, Vgrid, Var, Dynam%Masks,   &
                            Time_diag, omega, div, mfew, mfns ) 

!-----------------------------------------------------------------------

 end subroutine bgrid_core_time_diff

!#######################################################################

 subroutine bgrid_core_driver_end ( Var, Dynam )

!-----------------------------------------------------------------------
! Var       = prognostic variables
! Dynam     = data type for dynamical core constants and data
!-----------------------------------------------------------------------
   type   (prog_var_type), intent(in)    :: Var
   type(bgrid_dynam_type), intent(inout) :: Dynam
!-----------------------------------------------------------------------
!  terminate dynamics

   call bgrid_core_end ( Dynam )

!  write restart for prognostic variables

   call write_prog_var ( Var, Hgrid, Vgrid, fis, res, &
                         format=restart_output_format )

!  terminate integrals
   call bgrid_integrals_end

!  only prints diagnostics
   call bgrid_conserve_energy_end

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver_end

!#######################################################################
! The following routines do not really belong in this module but since
! they are not used within the core itself they will reside here.
!#######################################################################

 subroutine get_bottom_data ( a, b, a_bot, b_bot, k_bot ) 

!-------- Extract data from the lowest model level ---------
! a, b     = 3-D data fields
! a_bot,
!    b_bot = 2-D data fields containing data at lowest level
! k_bot    = index of lowest model level (optional)
!-----------------------------------------------------------

  real   , intent(in) , dimension(:,:,:) :: a    , b
  real   , intent(out), dimension(:,:)   :: a_bot, b_bot
  integer, intent(in) , dimension(:,:), optional :: k_bot

! returns the lowest level data (a_bot,b_bot) from 3d fields (a,b)

  integer :: i, j, kb, kb_min
    
                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot) 

  if ( kb_min == size(a,3) ) then
          a_bot = a(:,:,kb_min)
          b_bot = b(:,:,kb_min)
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a_bot(i,j) = a(i,j,kb)
          b_bot(i,j) = b(i,j,kb)
       enddo   
       enddo   
  endif

 end subroutine get_bottom_data

!#######################################################################

 subroutine put_bottom_data ( a_bot, b_bot, a, b, k_bot ) 

!-------- Insert data into the lowest model level ---------
! a_bot,
!    b_bot = 2-D data fields containing data for the lowest level
! a, b     = 3-D data fields with data inserted into lowest level
! k_bot    = index of lowest model level (optional)
!-----------------------------------------------------------

  real   , intent(in)   , dimension(:,:)   :: a_bot, b_bot
  real   , intent(inout), dimension(:,:,:) :: a    , b
  integer, intent(in)   , dimension(:,:), optional :: k_bot

! inserts the lowest level data (a_bot,b_bot) into 3d fields (a,b)

  integer :: i, j, kb, kb_min

                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot)

  if ( kb_min == size(a,3) ) then
          a(:,:,kb_min) = a_bot
          b(:,:,kb_min) = b_bot
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a(i,j,kb) = a_bot(i,j)
          b(i,j,kb) = b_bot(i,j)
       enddo
       enddo
  endif

 end subroutine put_bottom_data

!#######################################################################

 subroutine atmosphere_domain(Domain)
 type(domain2d), intent(inout) :: Domain

 Domain = Hgrid%Tmp%Domain

 end subroutine atmosphere_domain
!#######################################################################

end module bgrid_core_driver_mod


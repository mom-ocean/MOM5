#include <fms_platform.h>

module transport_matrix_mod  !{
!
!<CONTACT EMAIL="spk@ldeo.columbia.edu"> Samar Khatiwala
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Rick Slater
!</REVIEWER>
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies 
!</REVIEWER>
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Jennifer Simeon
!</REVIEWER>
!
!<OVERVIEW>
! Transport Matrix Method
!</OVERVIEW>
!
!<DESCRIPTION>
! Transport Matrix Method, for use in finding an approximate steady 
! state of tracers.  
!
! Ported to MOM4p0d by Samar Khatiwala spk@ldeo.columbia.edu 
! June-July 2007
!
! Ported to MOM4p1 by Stephen.Griffies 
! July 2007
!
! Some code clean-up; increased conformity to tracer module standard 
! jes.13JUN08
!
! further code cleanup by Richard.Slater Jan2009.
!
! **Preliminary testing only.
!
! This module saves out the components of the explicit and implicit
! transport matrix.  spk uses the output from the transport matrix module  
! to assemble the actual transport matrix with his matlab code.
! 
! In conjunction with the Newton-Krylov (NK) solver, the transport matrix 
! is used to do an accelerated forward model in an iterative manner
! which allows the NK solver to reach solution convergence with greater 
! efficiency.
!
! To create the transport matrix building blocks, passive, prognostic
! "tracers" are simulated.  These "tracers" are, in practice, dye
! tracers initialized to 1 at one or more grid points and initialized 
! to 0 elsewhere.  Since the "area of influence" of
! a tracer at a particular grid point is determined by the advection scheme,
! it becomes more practical to introduce a number of dye tracers at other
! grid points, such that the "tracers" are staggered and will have 
! non-overlapping areas of influence. These non-overlapping regions are
! called by spk, "tiles".  The total number of dye tracers 
! needed is independent of the horizontal grid resolution, but
! rather dependent upon the advection scheme's area of influence, i.e.
! the total number of grid points contained in a "tile". For a
! simple linear advection scheme, the number of tracers needed is about
! 10 x (number of vertical levels).
! spk has some special matlab code that generates the initial condition 
! for the "tracers", given the grid_spec of the model.
!
! The behavior of the tracers' advection and diffusion is averaged over 
! time (which the user may specify) and is saved out as a component 
! of the later-to-be-assembled transport matrix.
!
! To run the transport_matrix module, an initial condition
! must be created by spk.
! 
! Include the transport_matrix field_table and diag_table entries
! for the xml. Examples of these follow.
!
! Sample field table.
! -------------------------------------------------
! "tracer_packages","ocean_mod","transport_matrix"
!
! names = '01', '02', '03'
! horizontal-advection-scheme = mdfl_sweby
! vertical-advection-scheme = mdfl_sweby
! /
!
! --------------------------------------------------
!
! Sample diag table entry: enter as many "tracers" as you need
!                          with the naming convention exp_tm_# and imp_tm_#
!                          where # is a string as given abaove in "names"
! --------------------------------------------------
!"transport_matrix","exp_tm_01",  "exp_tm_01" ,"ocean_transport_matrix","all",.false.,"none",1
!"transport_matrix","imp_tm_01",  "imp_tm_01" ,"ocean_transport_matrix","all",.false.,"none",1
!"transport_matrix","exp_tm_02",  "exp_tm_02" ,"ocean_transport_matrix","all",.false.,"none",1
!"transport_matrix","imp_tm_02",  "imp_tm_02" ,"ocean_transport_matrix","all",.false.,"none",1
!"transport_matrix","exp_tm_03",  "exp_tm_03" ,"ocean_transport_matrix","all",.false.,"none",1
!"transport_matrix","imp_tm_03",  "imp_tm_03" ,"ocean_transport_matrix","all",.false.,"none",1
! --------------------------------------------------
!
!  
! SPK NOTES:
! 1) The calling sequence is as follows:
!  Top level driver (e.g., ocean_solo)
!    -> S/R ocean_model_init(Ocean, Time_init, Time_in, Time_step_ocean, ensemble_ocean)
!		Time%init       = Time_in .eq. Time_init
!		Time%Time_init  = Time_init
!		Time%Time_step  = Time_step_ocean
!		Time%model_time = Time_in
!		Time%itt        = 0
!		-> S/R ocean_prog_tracer_init 
!		   -> S/R ocean_tpm_init
!			  -> S/R transport_matrix_init
!		-> S/R ocean_tpm_start
!		   -> S/R transport_matrix_start
!	 Start time stepping loop
!	 do nc=1, num_cpld_calls     
!	   do no=1, num_ocean_calls
!		  ocean_seg_start = ( no .eq. 1 )     
!		  ocean_seg_end   = ( no .eq. num_ocean_calls )
!		  -> S/R update_ocean_model(Ice_ocean_boundary, Ocean_sfc, &
!                                           ocean_seg_start, ocean_seg_end, num_ocean_calls)
!			 Time%model_time = Time%model_time + Time%Time_step
!			 Time%itt        = Time%itt+1
!			 -> S/R ocean_tracer, S/R update_ocean_tracer              
!			 	   do explicit transport
!				   -> S/R transport_matrix_store_explicit 
!                                   (accumulate explicit matrix and reset tracer field to initial condition)
!  	 			   do implicit transport
!			 -> S/R ocean_tpm_tracer
!			    -> S/R transport_matrix_store_implicit 
!                           (accumulate implicit matrix and reset tracer field to initial condition)
!				   -> S/R transport_matrix_write(.FALSE.) (time average and write matrices)
!		  Time = Time + Time_step_ocean
!	   enddo
!	 enddo
!    -> S/R ocean_model_end
!       -> S/R ocean_tpm_end
!          -> S/R transport_matrix_end
!             -> S/R transport_matrix_write(.TRUE.) 
!                (time average and write matrices for multi year runs)
!
! 2) Time counters are incremented BEFORE calling S/R update_ocean_tracer, so the first time 
!    transport_matrix_store_explicit is called, itt (and hence myIter) will be 1.
!    The namelist parameter matrixStoreStartIter indicating the iteration number to begin 
!    accumulating matrices at should be RELATIVE to the current model run start (unlike 
!    in the MIT GCM where it refers to an absolute counter.
!    I am not entirely certain this is handled correctly below. Things might be off by 
!    1 time step.
!
!
!

!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
!      Khatiwala, S., M. Visbeck, M.A. Cane, 2005.
!  Accelerated simulation of passive tracers in ocean circulation models. 
!  Ocean Modelling, 9, 51-69.
! </REFERENCE>
! </INFO>
!
! $Id: transport_matrix.F90,v 20.0 2013/12/14 00:17:24 fms Exp $
!

use field_manager_mod,  only: fm_string_len, fm_path_name_len, fm_field_name_len
use field_manager_mod,  only: fm_get_length, fm_get_value
use mpp_mod,            only: stdout, mpp_error, FATAL
use diag_manager_mod,   only: send_data, register_diag_field

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use ocean_types_mod,    only: ocean_time_type, ocean_prog_tracer_type

!
!       force all variables to be "typed"
!

implicit none

!
!       Set all variables to be private by default

private

!
!       Public routines
!

public  :: transport_matrix_init
public  :: transport_matrix_start
public  :: transport_matrix_store_explicit
public  :: transport_matrix_store_implicit

!
!       Private routines
!

!
!       Public parameters
!

!
!       Private parameters
!

character(len=fm_field_name_len), parameter     :: package_name = 'transport_matrix'
character(len=48), parameter                    :: diag_name = 'transport_matrix'
character(len=48), parameter                    :: mod_name = 'transport_matrix_mod'
character(len=fm_string_len), parameter         :: default_restart_file = 'transport_matrix.res.nc'

!
!       Public variables
!

logical, public :: do_transport_matrix

!
!       Private types
!

type instance_type
  character(len=fm_field_name_len)                   :: name
  real, dimension(:,:,:),  _ALLOCATABLE              :: initial_value _NULL
  integer                                            :: index
  integer                                            :: id_extm
  integer                                            :: id_imtm
end type instance_type

!
!       Private variables
!

type(instance_type), dimension(:), allocatable          :: instance
integer                                                 :: instances
integer                                                 :: package_index

contains

!#######################################################################
! <SUBROUTINE NAME="transport_matrix_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by tracer package manager
! </DESCRIPTION>
!

subroutine transport_matrix_init  !{

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'transport_matrix_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_field_name_len+3)                      :: long_suffix

  integer :: stdoutunit 
  stdoutunit=stdout() 

!
!       Initialize the ocean transport matrix package
!

package_index = otpm_set_tracer_package(package_name,               &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')',       &
     restart_file = default_restart_file)

!
!       Check whether to use this package
!

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!
!       Check some things
!

write (stdoutunit,*)
if (instances .eq. 0) then  !{
  write (stdoutunit,*)                                            &
       trim(note_header), ' No instances'
  do_transport_matrix = .false.
else  !}{
  write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
  do_transport_matrix = .true.
endif  !}

!
!       Return if we don't want to use this package
!

if (.not. do_transport_matrix) then  !{
  return
endif  !}

!
!       allocate the instance structure
!

allocate ( instance(instances) )

!
!       loop over the names, saving them into the instance array

do n = 1, instances  !{

  if (fm_get_value(path_to_names, name, index = n)) then  !{
    instance(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         ' Bad field name for index ' // trim(name))
  endif  !}
   
enddo  !} n

do n = 1, instances !{

!----------------------------------------------------------------------------------------
! Register the prognostic tracers
!----------------------------------------------------------------------------------------

  ! Set up the instance name from field table
  name = instance(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

  ! determine the tracer names for this instance and assign indices

   instance(n)%index = otpm_set_prog_tracer('ptr' // trim(suffix),                              &
        package_name, longname = 'Transport matrix passive tracer' // trim(long_suffix),        &
        units = 'unit/kg', flux_units = 'unit/m^2/s',                                           &
        caller = trim(mod_name)//'('//trim(sub_name)//')')

enddo  !} n

return

end subroutine transport_matrix_init  !}
! </SUBROUTINE> NAME="transort_matrix_init"

!#######################################################################
! <SUBROUTINE NAME="transport_matrix_start">
!
! <DESCRIPTION>
!
! </DESCRIPTION>
!
subroutine transport_matrix_start(Time, T_prog, isd, ied, jsd, jed, nk, &
                                  grid_tracer_axes)  !{

!-----------------------------------------------------------------------
!       modules (have to come first)
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

  type(ocean_time_type),                      intent(in) :: Time
  type(ocean_prog_tracer_type), dimension(:), intent(in) :: T_prog
  integer,                                    intent(in) :: isd
  integer,                                    intent(in) :: ied
  integer,                                    intent(in) :: jsd
  integer,                                    intent(in) :: jed
  integer,                                    intent(in) :: nk
  integer, dimension(3),                      intent(in) :: grid_tracer_axes

!
!       local parameters
!

!
!       local variables
!

integer                             :: n
character(len=fm_field_name_len)    :: name
character(len=fm_field_name_len+3)  :: long_suffix
character(len=fm_field_name_len+1)  :: suffix

!
!       save initial values
!

do n = 1, instances  !{
  allocate( instance(n)%initial_value(isd:ied,jsd:jed,nk) )
  instance(n)%initial_value(:,:,:) = T_prog(instance(n)%index)%field(:,:,:,Time%tau)
enddo  !} n


!-----------------------------------------------------------------------
!     Set up analyses
!-----------------------------------------------------------------------
!

!
!       register the diag fields, i.e. the TMs
!

do n = 1, instances  !{

  ! Set up the instance name from field table
  name = instance(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

  ! determine the tracer names for this instance and register
   
  instance(n)%id_extm = register_diag_field(diag_name,                                          &
                 'exp_tm' // trim(suffix), grid_tracer_axes(1:3), Time%model_time,              &
                 'Explicit transport matrix for passive tracer ' // trim(long_suffix), ' ',     &
                 missing_value = -1.0e+10)
  instance(n)%id_imtm = register_diag_field(diag_name,                                          &
                 'imp_tm' // trim(suffix), grid_tracer_axes(1:3), Time%model_time,              &
                 'Implicit transport matrix for passive tracer ' // trim(long_suffix), ' ',     &
                 missing_value = -1.0e+10)

enddo  !} n

return

end subroutine transport_matrix_start  !}
! </SUBROUTINE> NAME="transport_matrix_start"


!#######################################################################
! <SUBROUTINE NAME="transport_matrix_store_explicit">
!
! <DESCRIPTION>
! For the time explicit tendencies. 
! </DESCRIPTION>
!

subroutine transport_matrix_store_explicit(Time, T_prog, isd, ied, jsd, jed, nk,        &
     isc, iec, jsc, jec, dtts, grid_tmask) 

  type(ocean_time_type),        intent(in)                  :: Time
  type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog
  integer, intent(in)                                       :: isd
  integer, intent(in)                                       :: ied
  integer, intent(in)                                       :: jsd
  integer, intent(in)                                       :: jed
  integer, intent(in)                                       :: nk
  integer, intent(in)                                       :: isc
  integer, intent(in)                                       :: iec
  integer, intent(in)                                       :: jsc
  integer, intent(in)                                       :: jec
  real, intent(in)                                          :: dtts
  real, dimension(isd:ied,jsd:jed,nk), intent(in)           :: grid_tmask

!
!       local parameters
!

!
!       local variables
!

integer :: n
real, dimension(isd:ied,jsd:jed,nk)     :: temp
logical                                 :: used

if(.not. do_transport_matrix) return

do n = 1, instances  !{

  if (instance(n)%id_extm .gt. 0) then !{
    temp(:,:,:) = (T_prog(instance(n)%index)%field(:,:,:,Time%taup1) - instance(n)%initial_value(:,:,:)) / dtts
    used = send_data(instance(n)%id_extm, temp, Time%model_time, rmask = grid_tmask,    &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif !}

  T_prog(instance(n)%index)%field(:,:,:,Time%taup1) = instance(n)%initial_value(:,:,:)

enddo  !} n

return

end subroutine transport_matrix_store_explicit  !}
! </SUBROUTINE> NAME="transport_matrix_store_explicit"

!#######################################################################
! <SUBROUTINE NAME="transport_matrix_store_implicit">
!
! <DESCRIPTION>
! For the time implicit tendencies. 
! </DESCRIPTION>
!

subroutine transport_matrix_store_implicit(Time, T_prog, isd, ied, jsd, jed, nk,        &
     isc, iec, jsc, jec, grid_tmask) !{

!---------------------------------------------------------------
! Modules first
!---------------------------------------------------------------

!---------------------------------------------------------------
! Arguments
!---------------------------------------------------------------


  type(ocean_time_type), intent(in)                         :: Time
  type(ocean_prog_tracer_type), dimension(:), intent(inout) :: T_prog
  integer, intent(in)                                       :: isd
  integer, intent(in)                                       :: ied
  integer, intent(in)                                       :: jsd
  integer, intent(in)                                       :: jed
  integer, intent(in)                                       :: nk
  integer, intent(in)                                       :: isc
  integer, intent(in)                                       :: iec
  integer, intent(in)                                       :: jsc
  integer, intent(in)                                       :: jec
  real, dimension(isd:ied,jsd:jed,nk), intent(in)           :: grid_tmask

!
!       local parameters
!

!
!       local variables
!

integer :: n
logical :: used

do n = 1, instances  !{

  if (instance(n)%id_imtm .gt. 0) then !{
      used = send_data(instance(n)%id_imtm, T_prog(instance(n)%index)%field(:,:,:,Time%taup1),  &
           Time%model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif !}

  T_prog(instance(n)%index)%field(:,:,:,Time%taup1) = instance(n)%initial_value(:,:,:)

enddo  !} n

return

end subroutine transport_matrix_store_implicit  !}
! </SUBROUTINE> NAME="transport_matrix_store_implicit"

end module transport_matrix_mod  !}

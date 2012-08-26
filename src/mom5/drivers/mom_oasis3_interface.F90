module mom_oasis3_interface_mod
!
!<CONTACT EMAIL="Dave.Bi@csiro.au"> Dave Bi (for OASIS3 hooks)
!</CONTACT>
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russ Fiedler (for OASIS3 hooks)
!</CONTACT>
!
!<OVERVIEW>
! Interface to OASIS3/PRISM2-5 coupling.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Interface to the OASIS3 

! This module serves as a hook between PRISM2-5/OASIS3 and mom4p1
! Only serial coupling has been tested at this stage (Oct 2009) although there are
! hints for the user for implementing different strategies.
!
! Standard mom surface quantities In Ocean_Ice_Boundary and Ocean_sfc may be passed
!
! The  namelist ocean_oasis3_interface_nml contains the quantities to be passed to and from the ocean 
! </DESCRIPTION>
!
! <NAMELIST NAME="mom_oasis3_interface_nml">
!
!   <DATA NAME="num_fields_in" TYPE="integer" DEFAULT="0">
!    The number of fields to be passed to mom from external sources
!   </DATA>
!   <DATA NAME="num_fields_out" TYPE="integer" DEFAULT="0">
!    The number of fields to be passed from mom to external progams
!   </DATA>
!   <DATA NAME="fields_in" TYPE="character(maxlen=8)" DEFAULT="''">
!    The fields to be passed to mom from external progams. Currently
!    valid names correspond to variables in the Ice_ocean_boundary structure.
!    These names must agree with names in the OASIS namcouple file
!    WARNING! Note truncation of names of ssw flux components and salt_flx.
!   </DATA>
!   <DATA NAME="fields_out" TYPE="character(maxlen=8)" DEFAULT="''">
!    The fields to be passed from mom to external progams. Currently
!    valid names correspond to variables in the Ocean_sfc structure.
!    These names must agree with the 'interpolated' names in the OASIS namcouple file
!   </DATA>
!   <DATA NAME="send_before_ocean_update" TYPE="logical" DEFAULT=".FALSE.">
!    TRUE if coupling strategy requires we send data to coupler BEFORE updating the ocean
!   </DATA>
!   <DATA NAME="send_after_ocean_update" TYPE="logical" DEFAULT=".FALSE.">
!    TRUE if coupling strategy requires we send data to coupler AFTER updating the ocean
!   </DATA>
!
! </NAMELIST>
!
!
#ifdef OASIS3

!prism stuff
use mod_prism_proto 
use mod_prism_def_partition_proto
use mod_prism_put_proto
use mod_prism_get_proto
use mod_comprism_proto

!MOM4 modules: 
use fms_mod,         only: file_exist
use fms_mod,         only: write_version_number, open_namelist_file, close_file, check_nml_error
use mpp_domains_mod, only: mpp_update_domains, domain2d
use mpp_mod,         only: input_nml_file, mpp_broadcast, mpp_pe, mpp_npes, mpp_root_pe
use mpp_mod,         only: mpp_error, FATAL, WARNING, NOTE
use mpp_mod,         only: stdlog, stdout
use mpp_domains_mod, only: mpp_get_compute_domain, &
                           mpp_get_data_domain, &
                           mpp_get_global_domain, &
                           mpp_global_field
use ocean_types_mod, only: ice_ocean_boundary_type, &
                           ocean_public_type, &
                           ocean_domain_type
use time_manager_mod, only: time_type

implicit none

public :: mom_prism_init, mom_prism_terminate, coupler_init, from_coupler, into_coupler, &
          write_coupler_restart

private

include "mpif.h" 

!
! OASIS3 variables:
!
integer :: id_component     ! Component ID
integer :: id_partition     ! Local partition ID
integer :: num_total_proc   ! Total number of processes
integer :: num_coupling_proc   ! Number of processes involved in the coupling

logical :: parallel_coupling

integer :: jf

!
! Coupling associated variables
!
integer :: ierr
integer :: limt, ljmt, i, j
integer :: step

character(len=6), parameter :: cp_modnam='mom4p1'  ! Component model name same as in namcouple

integer :: imt_global, jmt_global                  ! 2D global layout 
integer :: imt_local, jmt_local                    ! 2D global layout 
integer iisc,iiec,jjsc,jjec
integer iisd,iied,jjsd,jjed

integer, parameter :: max_fields_in=16
integer, parameter :: max_fields_out=6

integer, dimension(max_fields_in)  :: id_var_in  ! ID for fields to be rcvd
integer, dimension(max_fields_out) :: id_var_out ! ID for fields to be sent

character(len=8), dimension(max_fields_in)  :: mom_name_read  ! Standard mom Ice_Ocean_boundary names
character(len=8), dimension(max_fields_out) :: mom_name_write ! Standard mom Ocean_Sfc names


! Namelist variables
integer :: num_fields_in    ! Number of fields to be sent
integer :: num_fields_out   ! Number of fields to be rcvd
character(len=8),dimension(max_fields_out) :: fields_out
character(len=8),dimension(max_fields_in) :: fields_in
logical :: send_before_ocean_update=.FALSE.,       &
           send_after_ocean_update=.FALSE.  ! Control when to pass fields to coupler.

! Work array
real, allocatable,dimension(:,:)  :: vwork

contains

!-----------------------------------------------------------------------------------
subroutine mom_prism_init(mom4_local_comm)

integer :: mom4_local_comm
!--------------------!

!
! Initialize PRISM.
!

! Initialise MPI

call MPI_INIT(ierr)

call prism_init_comp_proto (id_component, cp_modnam, ierr)

if (ierr /= PRISM_Ok) then 
  call prism_abort_proto(id_component, 'mom4 prism_init','STOP 1')
endif

!
! PRISM attribution of local communicator. 
! 

call prism_get_localcomm_proto(mom4_local_comm, ierr) 

if (ierr /= PRISM_Ok) then
  call prism_abort_proto(id_component, 'mom4 prism_init','STOP 2')
endif

end subroutine mom_prism_init

!-----------------------------------------------------------------------------------
subroutine coupler_init(Dom, dt_cpld, Time, Time_step_coupled, Run_len)

! In this routine we set up all our arrays and determine which fields are to be passed to and fro.
! Determine the style of coupling
! 
! Time, Time_step_coupled & run_len not used at the moment. They are provided for convenience if the user wishes to
! use them for initialising.

type(domain2d)  :: Dom  
integer         :: dt_cpld
type(time_type), optional :: Time, Time_step_coupled, Run_len

integer, dimension(5) :: il_paral
integer, dimension(2) :: var_num_dims ! see below
integer, dimension(4) :: var_shape  ! see below

integer isg,ieg,jsg,jeg

integer ioun, io_status

integer ifield

logical fmatch 

integer :: stdoutunit,stdlogunit

namelist /mom_oasis3_interface_nml/ num_fields_in, num_fields_out,fields_in,fields_out, &
          send_before_ocean_update, send_after_ocean_update

! all processors read the namelist--

stdoutunit=stdout();stdlogunit=stdlog()

#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=mom_oasis3_interface_nml, iostat=io_status)
ierr = check_nml_error(io_status,'mom_oasis3_interface_nml')
#else
ioun = open_namelist_file()
read  (ioun, mom_oasis3_interface_nml,iostat=io_status)
ierr = check_nml_error(io_status,'mom_oasis3_interface_nml')
call close_file (ioun)
#endif
write (stdoutunit,'(/)')
write (stdoutunit, mom_oasis3_interface_nml)
write (stdlogunit, mom_oasis3_interface_nml)

call mpp_get_compute_domain(Dom,iisc,iiec,jjsc,jjec)
call mpp_get_data_domain(Dom,iisd,iied,jjsd,jjed)
call mpp_get_global_domain(Dom,isg,ieg,jsg,jeg)   

!  Global domain extent

imt_global=ieg-isg+1
jmt_global=jeg-jsg+1

!B: for parallel coupling, the above global layout does NOT make sense...
! Use local domains
!
! Original did all this on the data domain using temp arrays
!limt=iied-iisd+1
!ljmt=jjed-jjsd+1
! Use computational domain directly 

imt_local=iiec-iisc+1
jmt_local=jjec-jjsc+1
  
! Can't we get these numbers from elsewhere? Shouldn't be hardwired.  RASF
! It would be better to read in num_coupling_proc and parallel_coupling and perform sanity checking.
! Can we change some of the variable names? It's really hard to figure out what they mean.
! (OASIS3 naming convention problem?).

num_total_proc = mpp_npes()
!num_coupling_proc = num_total_proc     !multi-process coupling (real parallel cpl)!
num_coupling_proc = 1             !mono-process coupling	


! Type of coupling

if (num_coupling_proc == num_total_proc .and. num_total_proc /= 1) then
  parallel_coupling = .TRUE.
else
  parallel_coupling = .FALSE.
endif


!
!! vwork is the 'coupling' array that directly communicates with oasis3.
!! -- global (if parallel_coupling=.f.) domain or 
!! -- local  (if parallel_coupling=.t.) domain.
!! 
!
if (parallel_coupling) then
  allocate(vwork(iisc:iiec,jjsc:jjec))
else
  allocate(vwork(isg:ieg,jsg:jeg))
endif
  vwork=0.0

  ! Define name (as in namcouple) and declare each field sent/received by oce:
  !ice ==> ocn
  mom_name_read(:)=''

  mom_name_read(1)='u_flux'
  mom_name_read(2)='v_flux'
  mom_name_read(3)='t_flux'
  mom_name_read(4)='q_flux'
  mom_name_read(5)='salt_flx'
  mom_name_read(6)='lw_flux'
  mom_name_read(7)='sw_vdir'
  mom_name_read(8)='sw_vdif'
  mom_name_read(9)='sw_irdir'
  mom_name_read(10)='sw_irdif'
  mom_name_read(11)='lprec'
  mom_name_read(12)='fprec'
  mom_name_read(13)='runoff'
  mom_name_read(14)='calving'
  mom_name_read(15)='p'
  mom_name_read(16)='sw_flux'   ! For partitioning shortwave flux

  !ocn ==> ice
  mom_name_write(:)=''

  mom_name_write(1)='t_surf'
  mom_name_write(2)='s_surf'
  mom_name_write(3)='u_surf'
  mom_name_write(4)='v_surf'
  mom_name_write(5)='sea_lev'
  mom_name_write(6)='frazil'

  fmatch = .false.
  do jf = 1,num_fields_in
     do ifield=1,max_fields_in
       if ( trim(mom_name_read(ifield)) ==  trim(fields_in(jf)) ) then
         fmatch = .true.
         exit
       endif
     enddo
     if (.not. fmatch) then
       call mpp_error(FATAL,'coupler_init: Illegal input name' )
     endif
     fmatch = .false.
  enddo

  fmatch = .false.
  do jf = 1,num_fields_out
     do ifield=1,max_fields_in
       if ( trim(mom_name_write(ifield))  ==  trim(fields_out(jf)) ) then
         fmatch = .true.
         exit
       endif
     enddo
     if (.not. fmatch) then 
       call mpp_error(FATAL,'coupler_init: Illegal output name' )
     endif
     fmatch = .false.
  enddo

  il_paral (:) = 0
  if( parallel_coupling ) then
! ???
  !il_paral ( clim_strategy ) = clim_Box
  !il_paral ( clim_offset   ) = (ieg-isg+1)*(jjsc-1)+(iisc-1)
  !il_paral ( clim_SizeX    ) = iiec-iisc+1
  !il_paral ( clim_SizeY    ) = jjec-jjsc+1
  !il_paral ( clim_LdX      ) = ieg-isg+1
  ! send_before_ocean_update = .true. ! ?????
  ! send_after_ocean_update = .false. ! ?????
  
  else

    il_paral ( clim_strategy ) = clim_serial
    il_paral ( clim_offset   ) = 0
    il_paral ( clim_length   ) = imt_global * jmt_global
    send_before_ocean_update = .false.
    send_after_ocean_update = .true.
  endif

  !
  !  PRISM coupling fields declaration
  !
  var_num_dims(1)= 2      ! rank of coupling field
  var_num_dims(2)= 1      ! number of bundles in coupling field (always 1)
  var_shape(1)= lbound(vwork,1)       ! min index for the coupling field local dimension
  var_shape(2)= ubound(vwork,1)     ! max index for the coupling field local dim
  var_shape(3)= lbound(vwork,2)       ! min index for the coupling field local dim
  var_shape(4)= ubound(vwork,2)     ! max index for the coupling field local dim


if (mpp_pe() == mpp_root_pe() .or. parallel_coupling) then
  !
  ! The following steps need to be done:
  ! -> by the process if mom is monoprocess;
  ! -> only by the master process, if mom is parallel and only 
  !    master process is involved in the coupling;
  ! -> by all processes, if mom is parallel and all processes 
  ! are involved in the coupling.
  !
  ! - Define parallel partitions and allocate coupling fields accordingly: 
  ! (Refer to oasis/psmile/prism/modules/mod_prism_proto.F90 for integer value
  !  of clim_xxxx parameters)
  ! For each coupling field, association of a port to its symbolic name
  !       -Define the parallel decomposition associated to the port of each
  !        field; here no decomposition for all ports.
  !



  call prism_def_partition_proto (id_partition, il_paral, ierr)

  !

  do jf=1, num_fields_out
    call prism_def_var_proto (id_var_out(jf),fields_out(jf), id_partition, &
         var_num_dims, PRISM_Out, var_shape, PRISM_Real, ierr)
  enddo


  do jf=1, num_fields_in
    call prism_def_var_proto (id_var_in(jf), fields_in(jf), id_partition, &
         var_num_dims, PRISM_In, var_shape, PRISM_Real, ierr)
  enddo

  !
  !  PRISM end of declaration phase
  !
  call prism_enddef_proto (ierr)
  if (ierr /= PRISM_Ok) then
    call mpp_error(FATAL,'coupler_init: *** problem with prism_enddef_proto! ***' )
  endif

endif   ! 

end subroutine  coupler_init

!=======================================================================
subroutine into_coupler(step, Ocean_sfc, before_ocean_update, Time)
!------------------------------------------!

implicit none

type (ocean_public_type) :: Ocean_sfc
type (time_type), optional :: Time

integer, intent(in) :: step
logical, intent(in) :: before_ocean_update ! Flag to indicate whether
                                           ! we are calling before or after updating the ocean
integer:: jf, ierr, i, j

real, pointer, dimension(:,:) :: vwork_local

! Check if  we want to couple at this call.

if ( before_ocean_update .and. (.not. send_before_ocean_update) ) return
if ( (.not. before_ocean_update) .and. (.not. send_after_ocean_update) ) return


do jf = 1,num_fields_out

! Just point to array. Don't need to copy.

  coupling_fields_out: select case( trim(fields_out(jf)))
  case('t_surf')
    vwork_local => Ocean_sfc%t_surf
  case('s_surf')
    vwork_local => Ocean_sfc%s_surf
  case('u_surf')
    vwork_local => Ocean_sfc%u_surf
  case('v_surf')
    vwork_local => Ocean_sfc%v_surf
  case('sea_lev')
    vwork_local => Ocean_sfc%sea_lev
  case('frazil')
    vwork_local => Ocean_sfc%frazil
   case DEFAULT
   call mpp_error(FATAL,&
      '==>Error from into_coupler: Unknown quantity.')
   end select coupling_fields_out


  if (.not. parallel_coupling) then 

    call mpp_global_field(Ocean_sfc%domain, vwork_local(iisc:iiec,jjsc:jjec), vwork)

  endif

  if (mpp_pe() == mpp_root_pe() .or. parallel_coupling) then

    print *, 'MOM4: into_coupler putting field at step = ', trim(fields_out(jf)), ' ', step

    if (parallel_coupling) then 
      call prism_put_proto(id_var_out(jf),step,vwork_local,ierr)
    else
      call prism_put_proto(id_var_out(jf),step,vwork,ierr)
    endif

    if (ierr /= prism_ok.and. ierr < prism_sent) then
      call prism_abort_proto(id_component, 'MOM4 into_cpl','stop 1')
    endif

  
  endif 

enddo


return
end subroutine into_coupler

!-----------------------------------------------------------------------------------
subroutine from_coupler(step,Ice_ocean_boundary, Time)

! This is all highly user dependent. 

use constants_mod, only: hlv    ! 2.500e6 J/kg
implicit none

type (ice_ocean_boundary_type) :: Ice_ocean_boundary
type (time_type), optional :: Time

integer, intent(in) :: step

real :: frac_vis_dir=0.5*0.43, frac_vis_dif=0.5*0.43,             &
        frac_nir_dir=0.5*0.57, frac_nir_dif=0.5*0.57 ! shortwave partitioning

do jf = 1, num_fields_in

  if (mpp_pe() == mpp_root_pe() .or. parallel_coupling) then

    print *, 'MOM4: from_coupler getting fields at step = ', trim(fields_in(jf)), ' ', step
    call prism_get_proto(id_var_in(jf),step,vwork,ierr)

    if (ierr /= PRISM_Ok.and. ierr < prism_recvd) then
      call prism_abort_proto(id_component, 'MOM4 _get_ ','stop 1')
    endif

  endif

  if ( .not. parallel_coupling) then
    call mpp_broadcast(vwork, imt_global*jmt_global, mpp_root_pe())
  endif

! Handle conversions from raw field to what mom requires
! Sign conventions  & scaling  factors. Partioning of SSW etc.
! Would prefer to select on a generic name for clarity. but OASIS has 8 character limita on fields to be passed

  coupling_fields_in: select case ( trim(fields_in(jf) ))
  case('u_flux')
     Ice_ocean_boundary%u_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('v_flux')
     Ice_ocean_boundary%v_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('lprec')
     Ice_ocean_boundary%lprec(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('salt_flx')
     Ice_ocean_boundary%salt_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('calving') ! 
     Ice_ocean_boundary%calving(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('sw_vdir')
     Ice_ocean_boundary%sw_flux_vis_dir(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('sw_vdif')
     Ice_ocean_boundary%sw_flux_vis_dif(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('sw_irdir')
     Ice_ocean_boundary%sw_flux_nir_dir(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('sw_irdif')
     Ice_ocean_boundary%sw_flux_nir_dif(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('sw_flux')
     Ice_ocean_boundary%sw_flux_vis_dir(iisc:iiec,jjsc:jjec) =  frac_vis_dir*vwork(iisc:iiec,jjsc:jjec)
     Ice_ocean_boundary%sw_flux_vis_dif(iisc:iiec,jjsc:jjec) =  frac_vis_dif*vwork(iisc:iiec,jjsc:jjec)
     Ice_ocean_boundary%sw_flux_nir_dir(iisc:iiec,jjsc:jjec) =  frac_nir_dir*vwork(iisc:iiec,jjsc:jjec)
     Ice_ocean_boundary%sw_flux_nir_dif(iisc:iiec,jjsc:jjec) =  frac_nir_dif*vwork(iisc:iiec,jjsc:jjec)
  case('q_flux')
     Ice_ocean_boundary%q_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)/hlv
  case('t_flux')
     Ice_ocean_boundary%t_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('lw_flux')
     Ice_ocean_boundary%lw_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('runoff')
     Ice_ocean_boundary%runoff(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('p')
     Ice_ocean_boundary%p(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('fprec')
     Ice_ocean_boundary%fprec(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case DEFAULT
! Probable error. Leave as warning for the moment. RASF
   call mpp_error(WARNING,&
      '==>Warning from from_coupler: Unknown quantity.')
  end select coupling_fields_in

enddo    !jf

end subroutine from_coupler

!-----------------------------------------------------------------------------------
subroutine write_coupler_restart(step,write_restart)

logical, intent(in) :: write_restart
integer, intent(in) :: step

if ( write_restart ) then
   if (mpp_pe() == mpp_root_pe() .or. parallel_coupling) then
      do jf = 1,num_fields_out
         print *, 'MOM4: (write_coupler_restart) calling _put_restart at ', step, ' ', fields_out(jf)
         call prism_put_restart_proto(id_var_out(jf),step,ierr)

         if (ierr == PRISM_ToRest ) then  
            call mpp_error(NOTE,&
                 '==>Note from into_coupler: Written field to o2i file.')
         else
            call mpp_error(FATAL,&
                 '==>Error from into_coupler: writing field into o2i file failed.')
         endif
      enddo
   endif
endif
end subroutine write_coupler_restart

!-----------------------------------------------------------------------------------
subroutine mom_prism_terminate

! PRISM termination
!
! deallocate all coupling associated arrays here ......
! Note: prism won't terminate MPI
!
call prism_terminate_proto (ierr)

if (ierr .ne. PRISM_Ok) then
   call mpp_error(FATAL,&
      '==>Error from coupler_termination: Failure in prism_terminate_proto.')

endif

end subroutine mom_prism_terminate

!=====
#endif

end module mom_oasis3_interface_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!   <DATA NAME="frac_vis_dir" TYPE="real" DEFAULT="0.215">
!    Fraction of ssw apportion to direct visible radiation
!   </DATA>
!   <DATA NAME="frac_vis_dif" TYPE="real" DEFAULT="0.215">
!    Fraction of ssw apportioned to diffuse visible radiation
!   </DATA>
!   <DATA NAME="frac_nir_dir" TYPE="real" DEFAULT="0.285">
!    Fraction of ssw apportioned to direct near infrared radiation
!   </DATA>
!   <DATA NAME="frac_nir_dif" TYPE="real" DEFAULT="0.285">
!    Fraction of ssw apportioned to diffuse near infrared radiation
!   </DATA>
!
! </NAMELIST>
!
!

use mod_prism

!MOM4 modules: 
use fms_mod,         only: file_exist
use fms_mod,         only: write_version_number, open_namelist_file, close_file, check_nml_error
use mpp_domains_mod, only: mpp_update_domains, domain2d
use mpp_mod,         only: mpp_broadcast, mpp_pe, mpp_npes, mpp_root_pe
use mpp_mod,         only: mpp_error, FATAL, WARNING, NOTE
use mpp_mod,         only: stdlog, stdout
use mpp_domains_mod, only: mpp_get_compute_domain, &
                           mpp_get_data_domain, &
                           mpp_get_global_domain, &
                           mpp_global_field
use mpp_parameter_mod, only: GLOBAL_ROOT_ONLY, XUPDATE, YUPDATE
use ocean_types_mod, only: ice_ocean_boundary_type, &
                           ocean_public_type, &
                           ocean_domain_type
use time_manager_mod, only: time_type, get_time

! Timing

  use mpp_mod,                  only: mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_MODULE

  use auscom_ice_mod, only : il_out 
  use cpl_netcdf_setup_mod

implicit none

public :: mom_prism_init, mom_prism_terminate, coupler_init, from_coupler, into_coupler, &
          write_coupler_restart

public :: iisd, iied, jjsd, jjed, iisc, iiec, jjsc, jjec

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
character(len=12) :: choceout
character(len=6) :: chout

!
! Coupling associated variables
!
integer :: ierr
integer :: limt, ljmt, i, j
integer :: step

character(len=6), parameter :: cp_modnam='mom5xx'  ! Component model name same as in namcouple

integer :: imt_global, jmt_global                  ! 2D global layout 
integer :: imt_local, jmt_local                    ! 2D global layout 
integer iisc,iiec,jjsc,jjec
integer iisd,iied,jjsd,jjed

integer, parameter :: max_fields_in=25

integer, parameter :: max_fields_out=10

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

! Timing
integer :: id_oasis_send, id_oasis_recv
integer :: id_oasis_send1, id_oasis_recv1


!global domain
integer isg,ieg,jsg,jeg

!for paralell coupling 
real, allocatable,dimension(:,:)  :: vwork_2
integer, dimension(2) :: pe_layout
integer iiscpl,iiecpl,jjscpl,jjecpl
integer :: sendsubarray, recvsubarray , resizedrecvsubarray
integer, dimension(:), allocatable :: counts, disps
integer, dimension(2) :: starts,sizes,subsizes
integer(kind=mpi_address_kind) :: start, extent
real :: realvalue
integer :: col_comm
integer :: mom4_comm

real :: frac_vis_dir=0.5*0.43, frac_vis_dif=0.5*0.43,             &
        frac_nir_dir=0.5*0.57, frac_nir_dif=0.5*0.57 ! shortwave partitioning

contains

!-----------------------------------------------------------------------------------
subroutine mom_prism_init(mom4_local_comm, accessom2_config_dir, atm_intercomm)

integer, intent(out) :: mom4_local_comm
character(len=*), intent(in) :: accessom2_config_dir
integer, optional, intent(out) :: atm_intercomm
!--------------------!

logical initialized

!
! Initialize PRISM.
!

! Initialise MPI

call MPI_Initialized(initialized, ierr)
if (.not. initialized) then
    call MPI_INIT(ierr)
endif

call prism_init_comp_proto(id_component, cp_modnam, ierr, &
                           config_dir=accessom2_config_dir)

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

  if (present(atm_intercomm)) then
      call prism_get_intercomm(atm_intercomm, 'matmxx', ierr)
  endif

mom4_comm = mom4_local_comm

end subroutine mom_prism_init

!-----------------------------------------------------------------------------------
subroutine coupler_init(Dom, Time, Time_step_coupled, dt_cpld, Run_len, &
                        coupling_field_timesteps)

! In this routine we set up all our arrays and determine which fields are to be passed to and fro.
! Determine the style of coupling
! 
! Time, Time_step_coupled & run_len not used at the moment. They are provided for convenience if the user wishes to
! use them for initialising.

type(domain2d)  :: Dom  
integer,optional         :: dt_cpld
type(time_type),optional :: Time, Time_step_coupled, Run_len
integer, dimension(:), intent(in), optional :: coupling_field_timesteps

integer, dimension(5) :: il_paral
integer, dimension(2) :: var_num_dims ! see below
integer, dimension(4) :: var_shape  ! see below
integer :: run_len_seconds

!integer isg,ieg,jsg,jeg
integer :: jcol, irow

integer ioun, io_status

integer ifield

logical fmatch 

namelist /mom_oasis3_interface_nml/ num_fields_in, num_fields_out,fields_in,fields_out, &
          send_before_ocean_update, send_after_ocean_update, frac_vis_dir, frac_vis_dif, &
          frac_nir_dir, frac_nir_dif


! all processors read the namelist--

ioun = open_namelist_file()
read  (ioun, mom_oasis3_interface_nml,iostat=io_status)
write (stdout(),'(/)')
write (stdout(), mom_oasis3_interface_nml)
write (stdlog(), mom_oasis3_interface_nml)
ierr = check_nml_error(io_status,'mom_oasis3_interface_nml')
call close_file (ioun)

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
  
  pe_layout(1)=imt_global/imt_local
  pe_layout(2)=jmt_global/jmt_local
  


! Can't we get these numbers from elsewhere? Shouldn't be hardwired.  RASF
! It would be better to read in num_coupling_proc and parallel_coupling and perform sanity checking.
! Can we change some of the variable names? It's really hard to figure out what they mean.
! (OASIS3 naming convention problem?).

num_total_proc = mpp_npes()
num_coupling_proc = num_total_proc     !multi-process coupling (real parallel cpl)!

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
  mom_name_read(13)='runof'
  mom_name_read(14)='calving'
  mom_name_read(15)='p'
  mom_name_read(16)='sw_flux'   ! For partitioning shortwave flux
  mom_name_read(17)='aice'   ! Ice fraction
  mom_name_read(18)='mh_flux'   ! Heat flux due to melting
  mom_name_read(19)='wfimelt'  !Water flux due to ice melting
  mom_name_read(20)='wfiform'  !Water flux due to ice forming 
  mom_name_read(21)='licefw'  ! Water flux from land ice
  mom_name_read(22)='liceht'  ! Heat flux from land ice
  mom_name_read(23)='wnd_io'  !
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  mom_name_read(24)='iof_nit'  !
  mom_name_read(25)='iof_alg'  !
#endif

  !ocn ==> ice
  mom_name_write(:)=''

  mom_name_write(1)='t_surf'
  mom_name_write(2)='s_surf'
  mom_name_write(3)='u_surf'
  mom_name_write(4)='v_surf'
  mom_name_write(5)='sea_lev'
  mom_name_write(6)='frazil'
  mom_name_write(7)='dssldx'
  mom_name_write(8)='dssldy'
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  mom_name_write(9)='n_surf'
  mom_name_write(10)='alg_surf'
#endif

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

!DHB
#if defined(DEBUG)
    if (mpp_pe() == mpp_root_pe() .or. parallel_coupling) then
      il_out = 85 + mpp_pe()
      write(chout,'(I6.6)'), il_out
      choceout='oceout'//trim(chout)
      open(il_out,file=choceout,form='formatted')
    endif
#endif


#if defined(DEBUG)
    write(il_out, *) "compute domain:",mpp_pe(), iisc, iiec, jjsc, jjec 
    write(il_out, *) "data domain:",mpp_pe(), iisd, iied, jjsd, jjed 
    write(il_out, *) "global domain:",mpp_pe(), isg, ieg, jsg, jeg 
    write(il_out, *) "global layout nx x ny:",pe_layout(1), pe_layout(2) 
    flush(il_out)
#endif

  il_paral (:) = 0
  if (.not. parallel_coupling) then 
    il_paral ( clim_strategy ) = clim_serial
    il_paral ( clim_offset   ) = 0
    il_paral ( clim_length   ) = imt_global * jmt_global
    send_before_ocean_update = .false.
    send_after_ocean_update = .true.
  else !if( parallel_coupling .and. mpp_pe() < pe_layout(1) ) then

  il_paral ( clim_strategy ) = clim_Box
!every cpu gets coupling
  il_paral ( clim_offset   ) = (ieg-isg+1)*(jjsc-jsg)+(iisc-isg)
  il_paral ( clim_SizeX    ) = iiec-iisc+1
  il_paral ( clim_SizeY    ) = jjec-jjsc+1
  il_paral ( clim_LdX      ) = ieg-isg+1
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


if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling )) then
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

  call prism_def_partition_proto (id_partition, il_paral, ierr, imt_global*jmt_global)

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
  if (present(run_len)) then
    call get_time(run_len, run_len_seconds)
    call prism_enddef_proto (ierr, runtime=run_len_seconds, &
                             coupling_field_timesteps=coupling_field_timesteps)
  else
    call prism_enddef_proto (ierr)
  endif

  if (ierr /= PRISM_Ok) then
    call mpp_error(FATAL,'coupler_init: *** problem with prism_enddef_proto! ***' )
  endif

endif   ! 
  id_oasis_recv = mpp_clock_id('oasis_recv', flags=MPP_CLOCK_SYNC,grain=CLOCK_MODULE)
  id_oasis_send = mpp_clock_id('oasis_send', flags=MPP_CLOCK_SYNC,grain=CLOCK_MODULE)
  id_oasis_recv1 = mpp_clock_id('oasis_recv1', grain=CLOCK_MODULE)
  id_oasis_send1 = mpp_clock_id('oasis_send1', grain=CLOCK_MODULE)

end subroutine  coupler_init

!=======================================================================
subroutine into_coupler(step, Ocean_sfc, Time, before_ocean_update)
!------------------------------------------!

use ocean_operators_mod, only : GRAD_BAROTROPIC_P       !GRAD_SURF_sealev
use auscom_ice_mod, only      : auscom_ice_heatflux_new
use auscom_ice_mod, only      : chk_o2i_fields, chk_fields_period, chk_fields_start_time

implicit none

type (ocean_public_type) :: Ocean_sfc
type (time_type),optional         :: Time

integer, intent(in) :: step
logical, intent(in),optional :: before_ocean_update ! Flag to indicate whether
                                           ! we are calling before or after updating the ocean
integer:: jf, ierr, i, j

real, pointer, dimension(:,:) :: vwork_local
real, dimension(iisd:iied,jjsd:jjed) :: vtmp
real, dimension(isg:ieg,jsg:jeg) :: gtmp

  character*80 :: fname='fields_o2i_in_ocn.nc'
  integer :: ncid,currstep,ll,ilout
  data currstep/0/
  save currstep

! Check if  we want to couple at this call.

if ( before_ocean_update .and. (.not. send_before_ocean_update) ) return
if ( (.not. before_ocean_update) .and. (.not. send_after_ocean_update) ) return

  
  if (mpp_pe() == mpp_root_pe()) then
    if (chk_o2i_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time)) then
        currstep=currstep+1
      if (currstep == 1) then
        call create_ncfile(trim(fname),ncid,imt_global,jmt_global,ll=1,ilout=il_out)
      endif
#if defined(DEBUG)
      write(il_out,*) 'opening file at nstep = ', trim(fname), step
#endif
      call ncheck( nf_open(trim(fname),nf_write,ncid) )
      call write_nc_1Dtime(real(step),currstep,'time',ncid)
    endif
  endif

!20110329: gradient is now calculated in initialize_ocean_sfc to avoid the domain boundary shift!
!!!Ocean_sfc%gradient(:,:,:) = GRAD_BAROTROPIC_P(Ocean_sfc%sea_lev(:,:))
!frazil is actually updated 
call auscom_ice_heatflux_new(Ocean_sfc)

     call mpp_clock_begin(id_oasis_send)
do jf = 1,num_fields_out
     if(jf .ne. 1) call mpp_clock_begin(id_oasis_send1)

! Just point to array. Don't need to copy.
  vtmp(:,:) = 0.0
  coupling_fields_out: select case( trim(fields_out(jf)))
  case('t_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%t_surf(iisd:iied,jjsd:jjed)
  case('s_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%s_surf(iisd:iied,jjsd:jjed)
  case('u_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%u_surf(iisd:iied,jjsd:jjed)
  case('v_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%v_surf(iisd:iied,jjsd:jjed)
  case('sea_lev')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%sea_lev(iisd:iied,jjsd:jjed)
  case('frazil')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%frazil(iisd:iied,jjsd:jjed)
  case('dssldx') 
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%gradient(iisd:iied,jjsd:jjed,1)
  case('dssldy') 
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%gradient(iisd:iied,jjsd:jjed,2)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  case('n_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%n_surf(iisd:iied,jjsd:jjed)
  case('alg_surf')
    vtmp(iisd:iied,jjsd:jjed) = Ocean_sfc%alg_surf(iisd:iied,jjsd:jjed)
#endif
  case DEFAULT
  call mpp_error(FATAL,&
      '==>Error from into_coupler: Unknown quantity.')
  end select coupling_fields_out

  if(.not. parallel_coupling) then  
      call mpp_global_field(Ocean_sfc%domain, vtmp(iisc:iiec,jjsc:jjec), vwork)
  else
      vwork = vtmp
  endif

  if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling )) then

      call prism_put_proto(id_var_out(jf),step,vwork,ierr) 
  
      if (chk_o2i_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time)) then
        if (parallel_coupling) then 
          call mpp_global_field(Ocean_sfc%domain, vtmp(iisc:iiec,jjsc:jjec), gtmp)
          if( mpp_pe() == mpp_root_pe() )  &
            call write_nc2D(ncid, trim(fields_out(jf)), gtmp, 1, imt_global,jmt_global, &
                        currstep,ilout=il_out)
        else
          if( mpp_pe() == mpp_root_pe() )  &
            call write_nc2D(ncid, trim(fields_out(jf)), vwork, 1, imt_global,jmt_global, &
                        currstep,ilout=il_out)
        endif
    endif

    if (ierr /= prism_ok.and. ierr < prism_sent) then
      call prism_abort_proto(id_component, 'MOM4 into_cpl','stop 1')
    endif

   endif 

     if(jf .ne. 1) call mpp_clock_end(id_oasis_send1)

enddo
     call mpp_clock_end(id_oasis_send)

  if (chk_o2i_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time) .and. (mpp_pe() == mpp_root_pe())) then
    call ncheck(nf_close(ncid))
  endif

return
end subroutine into_coupler

!-----------------------------------------------------------------------------------
subroutine from_coupler(step,Ocean_sfc,Ice_ocean_boundary, Time)

! This is all highly user dependent. 

use constants_mod, only  : hlv    ! 2.500e6 J/kg
use auscom_ice_mod, only : chk_i2o_fields, chk_fields_period, chk_fields_start_time
implicit none

type (ocean_public_type) :: Ocean_sfc
type (ice_ocean_boundary_type) :: Ice_ocean_boundary
type (time_type),optional         :: Time

real, dimension(isg:ieg,jsg:jeg) :: gtmp

integer, intent(in) :: step

  character*80 :: fname = 'fields_i2o_in_ocn.nc'
  integer :: ncid,currstep,ll,ilout
  data currstep/0/
  save currstep

  if (mpp_pe() == mpp_root_pe()) then
    if (chk_i2o_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time)) then 
        currstep=currstep+1
      if (currstep == 1) then
        call create_ncfile(trim(fname),ncid,imt_global,jmt_global,ll=1,ilout=il_out)
      endif
#if defined(DEBUG)
      write(il_out,*) 'opening file at nstep = ', trim(fname), step
#endif
      call ncheck( nf_open(trim(fname),nf_write,ncid) )
      call write_nc_1Dtime(real(step),currstep,'time',ncid)
    endif
  endif

     call mpp_clock_begin(id_oasis_recv)
do jf =  1, num_fields_in
     if(jf .ne. 1) call mpp_clock_begin(id_oasis_recv1)

  !if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling .and. mpp_pe() < pe_layout(1))) then
  if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling )) then
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
     ! Sign conventions for Ice_ocean_boundary%salt_flux according to component model:
     !  - For MOM, if Ice_ocean_boundary%salt_flux < 0 then there is a transfer
     !    of salt from ice to liquid.
     !  - For CICE, if Ice_ocean_boundary%salt_flux < 0 then there is a
     !    transfer of salt from liquid to ice.
     ! We change the sign on vwork to accord with the MOM sign convention when
     ! filling Ice_ocean_boundary%salt_flux, since Ice_ocean_boundary%salt_flux is
     ! used in MOM to fill its local salt flux array.
     Ice_ocean_boundary%salt_flux(iisc:iiec,jjsc:jjec) =  -vwork(iisc:iiec,jjsc:jjec)
  case('calving')
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
  case('runof')
     Ice_ocean_boundary%runoff(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('p')
     Ice_ocean_boundary%p(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('fprec')
     Ice_ocean_boundary%fprec(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('aice')
     Ice_ocean_boundary%aice(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
  case('iof_nit')
     Ice_ocean_boundary%iof_nit(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('iof_alg')
     Ice_ocean_boundary%iof_alg(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
#endif
  case('mh_flux')
     Ice_ocean_boundary%mh_flux(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('wfimelt')
     Ice_ocean_boundary%wfimelt(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('wfiform')
     Ice_ocean_boundary%wfiform(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('licefw')
     Ice_ocean_boundary%licefw(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('liceht')
     Ice_ocean_boundary%liceht(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case('wnd_io')
     Ice_ocean_boundary%wnd(iisc:iiec,jjsc:jjec) =  vwork(iisc:iiec,jjsc:jjec)
  case DEFAULT
! Probable error. Leave as warning for the moment. RASF
   call mpp_error(WARNING,&
      '==>Warning from from_coupler: Unknown quantity.')
  end select coupling_fields_in

      if (chk_i2o_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time)) then
        if (parallel_coupling ) then
          call mpp_global_field(Ocean_sfc%domain, vwork(iisc:iiec,jjsc:jjec), gtmp) 
          if( mpp_pe() == mpp_root_pe() )  &
            call write_nc2D(ncid, trim(fields_in(jf)), gtmp, 1, imt_global,jmt_global, &
                          currstep,ilout=il_out)
        else
          if( mpp_pe() == mpp_root_pe() )  &
            call write_nc2D(ncid, trim(fields_in(jf)), vwork, 1, imt_global,jmt_global, &
                          currstep,ilout=il_out)
        endif
      endif

     if(jf .ne. 1) call mpp_clock_end(id_oasis_recv1)
enddo    !jf
     call mpp_clock_end(id_oasis_recv)

  if (chk_i2o_fields .and. (mod(step, chk_fields_period) == 0) .and. (step >= chk_fields_start_time) .and. (mpp_pe() == mpp_root_pe())) then
    call ncheck(nf_close(ncid))
  endif

end subroutine from_coupler

!-----------------------------------------------------------------------------------
subroutine write_coupler_restart(step,Ocean_sfc,write_restart)

use auscom_ice_mod, only      : auscom_ice_heatflux_new

logical, intent(in) :: write_restart
integer, intent(in) :: step
type (ocean_public_type) :: Ocean_sfc

integer :: ncid,ll,ilout
real, dimension(iisd:iied,jjsd:jjed) :: vtmp
real, dimension(isg:ieg,jsg:jeg) :: gtmp
character(len=8) :: fld_ice

if ( write_restart ) then
   if ( mpp_pe() == mpp_root_pe() ) then
     call create_ncfile('o2i.nc', ncid, imt_global,jmt_global, ll=1, ilout=il_out)
     call write_nc_1Dtime(real(step), 1, 'time', ncid)
   endif
   !update frazil
   call auscom_ice_heatflux_new(Ocean_sfc)
   !if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling .and.  mpp_pe() < pe_layout(1)) ) then
   if (mpp_pe() == mpp_root_pe() .or. (parallel_coupling) ) then
      do jf = 1,num_fields_out
        vtmp(:,:) = 0.0
        select case (trim(fields_out(jf)))
          case('t_surf'); vtmp = Ocean_sfc%t_surf(iisd:iied,jjsd:jjed); fld_ice='sst_i'
          case('s_surf'); vtmp = Ocean_sfc%s_surf(iisd:iied,jjsd:jjed); fld_ice='sss_i'
          case('u_surf'); vtmp = Ocean_sfc%u_surf(iisd:iied,jjsd:jjed); fld_ice='ssu_i'
          case('v_surf'); vtmp = Ocean_sfc%v_surf(iisd:iied,jjsd:jjed); fld_ice='ssv_i'
          case('dssldx'); vtmp = Ocean_sfc%gradient(iisd:iied,jjsd:jjed,1); fld_ice='sslx_i'
          case('dssldy'); vtmp = Ocean_sfc%gradient(iisd:iied,jjsd:jjed,2); fld_ice='ssly_i'
          case('frazil'); vtmp = Ocean_sfc%frazil(iisd:iied,jjsd:jjed); fld_ice='pfmice_i'
#if defined(ACCESS_OM) && defined(CSIRO_BGC)
          case('n_surf'); vtmp = Ocean_sfc%n_surf(iisd:iied,jjsd:jjed); fld_ice='ssn_i'
          case('alg_surf'); vtmp = Ocean_sfc%alg_surf(iisd:iied,jjsd:jjed); fld_ice='ssalg_i'
#endif
        end select

        if (parallel_coupling) then
          call mpp_global_field(Ocean_sfc%domain, vtmp(iisc:iiec,jjsc:jjec), &
                                gtmp, flags=GLOBAL_ROOT_ONLY+XUPDATE+YUPDATE)
        else
          call mpp_global_field(Ocean_sfc%domain, vtmp(iisc:iiec,jjsc:jjec), &
                                vwork, flags=GLOBAL_ROOT_ONLY+XUPDATE+YUPDATE)
        end if

        if (mpp_pe() == mpp_root_pe()) then
          if (parallel_coupling) then
            call write_nc2D(ncid, trim(fld_ice), gtmp, 2, imt_global,jmt_global, &
                            1,ilout=il_out)
          else
            call write_nc2D(ncid, fld_ice, vwork, 2, imt_global,jmt_global, 1, ilout=il_out)
          end if
        end if
      end do
   end if

   if (mpp_pe() == mpp_root_pe()) call ncheck( nf_close(ncid) )

endif
end subroutine write_coupler_restart

!-----------------------------------------------------------------------------------
subroutine mom_prism_terminate

! PRISM termination
!
! deallocate all coupling associated arrays here ......
! Note: prism won't terminate MPI
!

    call MPI_Barrier(mom4_comm, ierr)


call prism_terminate_proto (ierr)

if (ierr .ne. PRISM_Ok) then
   call mpp_error(FATAL,&
      '==>Error from coupler_termination: Failure in prism_terminate_proto.')

endif

end subroutine mom_prism_terminate

end module mom_oasis3_interface_mod

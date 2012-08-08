module ocean_blob_static_bottom_mod
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Contains static bottom blob parameterisations, namely the Campin and 
! Goossee (1999) scheme and some variations.
!</OVERVIEW>
!
!<DESCRIPTION>
! Presently, there is only one static bottom blob scheme implemented 
! (although, there are three potential variants), and there are no plans
! to implement any additional schemes.
!
! The scheme that is implemented emulates the Campin and Goossee (1999) 
! scheme, as well as having two additional variations for this scheme, 
! in which the "plumbing" is change.
!
! Details of the variations can be found in Bates et al. (2010) and are 
! controlled by the namelist options overflow_no_return and
! overflow_one_return.
!</DESCRIPTION>
!<INFO>
!
! <REFERENCE>
! Bates, M.L., Griffies, S.M., England, M.H., Adcroft, A.J. (2009)
! Lagrangian blobs of buoyancy embedded in Eulerian models: a framework
! to parameterise vertical and downslope motion of gravitationally unstable
! water parcels.  Unpublished Notes.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Elements of mom4p1 (2009)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! Campin, J.-M., Goossee, H., 1999, Parameterization of density-driven
! downsloping flow for a coarse-resolution ocean model in z-coordinate.
! Tellus 51A (3), 412-430
! </REFERENCE>
!
!</INFO>
!
!<NAMELIST NAME="ocean_blob_static_bottom_nml">
!
!  <DATA NAME="blob_overflow" TYPE="logical">
!  if true, will use a Campin and Goosse (1999) style overflow 
!  formulation.
!  Default is blob_overflow=.false.
!  </DATA>
!
!  <DATA NAME="blob_overflow_mu" UNITS="1/s" TYPE="real">
!  Frictional dissipation used in blob_overflow scheme
!  Default is blob_overflow_mu=1.0e-4
!  </DATA>
!
!  <DATA NAME="blob_overflow_delta" TYPE="real">
!  Fraction of grid cell participating in overflow
!  Valid values are 0<=delta<=1
!  Default is blob_overflow_delta=1/3
!  </DATA>
!
!  <DATA NAME="blob_overflow_umax" UNITS="m/s" TYPE="real">
!  Maximum downslope speed allowed for overflow
!  Default is blob_overflow_umax=0.01
!  </DATA>
!  
!  <DATA NAME="overflow_no_return" TYPE="logical">
!  When .false. creates return blobs to replicate the
!  original Campin and Goosse scheme.  When .true. only
!  creates blobs that sink.  See further overflow_one_return
!  </DATA>
!
!  <DATA NAME="overflow_one_return" TYPE="logical">
!  Creates a single return blob when .true.  Cannot be .true.
!  when overflow_no_return is also .true.
!  </DATA>
!
!</NAMELIST>
!
use fms_mod,         only: open_namelist_file, check_nml_error, stdout, stdlog
use fms_mod,         only: stderr, FATAL, error_mesg, mpp_error, close_file
use mpp_domains_mod, only: mpp_update_domains
use mpp_mod,         only: NULL_PE, mpp_sync_self, mpp_send, mpp_recv

use ocean_blob_util_mod,  only: insert_blob, count_blob, free_blob_memory, blob_delete
use ocean_density_mod,    only: density
use ocean_parameters_mod, only: onehalf, rho0r, grav
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_density_type, ocean_blob_type, blob_grid_type

implicit none

private

type(ocean_grid_type),   pointer :: Grd  => NULL()
type(ocean_domain_type), pointer :: Dom  => NULL()
type(blob_grid_type),    pointer :: Info => NULL()

logical :: module_is_initialized=.false.
integer :: num_prog_tracers
integer :: index_temp
integer :: index_salt

integer :: isd,ied,jsd,jed,isc,iec,jsc,jec,nk
integer :: rea_buff_size

integer, dimension(4) :: ip, jq, ish, jsh, ieh, jeh

real, dimension(:,:),   allocatable :: datdtime_r   ! 1./(Grd%dat*dtime)
real, dimension(:,:,:), allocatable :: topog_step, topog_slope
real, dimension(:,:),   allocatable :: slope_x, slope_y
real, dimension(:,:),   allocatable :: max_overflow_flux
real :: dtime
real :: overflow_factor

type, private :: blob_buffer_type
   integer :: numblobs
   integer :: size
   integer, allocatable, dimension(:,:) :: integer_buff
   real,    allocatable, dimension(:,:) :: real_buff
   logical, allocatable, dimension(:,:) :: logical_buff
end type blob_buffer_type

type(blob_buffer_type), pointer :: Ebuffer_out,  Ebuffer_in
type(blob_buffer_type), pointer :: Wbuffer_out,  Wbuffer_in
type(blob_buffer_type), pointer :: Nbuffer_out,  Nbuffer_in
type(blob_buffer_type), pointer :: Sbuffer_out,  Sbuffer_in

! buffers and buffer variables to send blobs to adjacent PE's
integer, parameter :: delta_buffer = 25   ! Size by which to increment the buffer
integer, parameter :: int_buff_params = 8 ! number of integer blob attributes
integer, parameter :: log_buff_params = 1 ! number of logical blob attributes

public blob_static_bottom_init
public blob_overflow_like
public blob_static_bottom_end

! namelist defaults
logical :: use_this_module     =.false.  ! whether to use this module or not
logical :: overflow_no_return  =.false.  ! for the overflow-like scheme
logical :: overflow_one_return =.false.  ! for the overflow-like scheme
real    :: blob_overflow_mu    = 1.0e-4  ! frictional dissipation rate (sec^-1) at bottom  
real    :: blob_overflow_delta = 0.3333  ! fraction of a grid cell participating in overflow 
real    :: blob_overflow_umax  = 0.01    ! maximum downslope speed (m/s)

namelist /ocean_blob_static_bottom_nml/ use_this_module, overflow_no_return, &
                                        overflow_one_return, blob_overflow_mu, &
                                        blob_overflow_delta, blob_overflow_umax


contains

!######################################################################
! <SUBROUTINE NAME="blob_static_bottom_init">
!
! <DESCRIPTION>
! Initialises the static bottom blob module.
! </DESCRIPTION>
!
subroutine blob_static_bottom_init(Grid, Domain, PE_info, num_tracers, itemp, isalt, dtimein)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(blob_grid_type),    intent(in), target :: PE_info
  integer, intent(in) :: num_tracers
  integer, intent(in) :: itemp
  integer, intent(in) :: isalt
  real,    intent(in) :: dtimein

  integer :: i, j, m
  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit 

  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
         '==>Error in ocean_blob_static_bottom_mod (ocean_blob_static_bottom_init):' &
         //' module already initialized')
  endif 

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_static_bottom_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_static_bottom_nml)  
  write (stdlogunit,ocean_blob_static_bottom_nml)
  ierr = check_nml_error(io_status,'ocean_blob_static_bottom_nml')
  call close_file (ioun)

  module_is_initialized = .true.

  if (.not. use_this_module) return

  write(stdoutunit,'(a)') & 
       ' ==>Note: Using the Lagrangian overflow scheme.'
  if (overflow_no_return .and. overflow_one_return) then
     call mpp_error(FATAL, '  ==>Error in ocean_blobs_mod '   &
          //'(ocean_blob_init): overflow_no_return and '      &
          //'overflow_one_return are selected.  Only one may '&
          //'be selected')
  elseif(overflow_one_return) then
     write(stdoutunit,'(a)') '  ==>Note: using a single return blob'
  elseif(overflow_no_return) then
     write(stdoutunit,'(a)') '  ==>Note: not using any return blobs'
  else
     write(stdoutunit,'(a)') '  ==>Note: using the full C&G scheme'
  endif

  Grd  => Grid
  Dom  => Domain
  Info => PE_info

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  nk = Grd%nk

  dtime             = dtimein
  index_temp        = itemp
  index_salt        = isalt
  num_prog_tracers  = num_tracers

  allocate(datdtime_r(isd:ied,jsd:jed))
  datdtime_r(:,:) = Grd%datr(:,:)/dtime

  allocate(max_overflow_flux(isd:ied,jsd:jed))
  max_overflow_flux(:,:) = onehalf*Grd%dat(:,:)/dtime

  overflow_factor = grav*blob_overflow_delta*rho0r/blob_overflow_mu

  ! compute topographic slope arrays for the i-slope and j-slope
  ! slopes are centered on the i-face and j-face of tracer cells 
  allocate(slope_x(isd:ied,jsd:jed))
  allocate(slope_y(isd:ied,jsd:jed))
  slope_x(:,:) = 0.0 
  slope_y(:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
     enddo
  enddo

  call mpp_update_domains(slope_x(:,:),Dom%domain2d)
  call mpp_update_domains(slope_y(:,:),Dom%domain2d)
  
  ! labels for the four tracer cells surrounding a 
  ! central tracer cell (moving counter-clockwise)
  m=1 ; ip(m)=1  ; jq(m)=0
  m=2 ; ip(m)=0  ; jq(m)=1
  m=3 ; ip(m)=-1 ; jq(m)=0
  m=4 ; ip(m)=0  ; jq(m)=-1

  
  ! a useful index for accessing only certain parts of the halo
  ish(:) = isc; ieh(:) = iec; jsh(:) = jsc; jeh(:) = jec
  if (Info%pe_W .ne. NULL_PE) ish(1) = isc - 1
  if (Info%pe_S .ne. NULL_PE) jsh(2) = jsc - 1
  if (Info%pe_E .ne. NULL_PE) ieh(3) = iec + 1
  if (Info%pe_N .ne. NULL_PE) jeh(4) = jec + 1


  ! a useful index for accessing only certain parts of the halo
!  ish(:) = isc; ieh(:) = iec; jsh(:) = jsc; jeh(:) = jec
!  if (Info%pe_W .ne. NULL_PE) ish(1) = isc - 1
!  if (Info%pe_S .ne. NULL_PE) jsh(2) = jsc - 1
!  if (Info%pe_E .ne. NULL_PE) ieh(3) = iec + 1
!  if (Info%pe_N .ne. NULL_PE) jeh(4) = jec + 1
  
  ! compute directions from an (i,j) point where topography deepens.
  ! these directions may potentially have downslope flow. 
  ! insist that downslope flow occurs only when there are more kmt 
  ! cells in the adjacent column. 
  ! also insist that downslope flow does not involve k=1 cells. 
  allocate (topog_step(isd:ied,jsd:jed,4))
  allocate (topog_slope(isd:ied,jsd:jed,4))
  topog_step(:,:,:)  = 0.0
  topog_slope(:,:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        if(Grd%kmt(i,j) > 1) then 
           if(Grd%kmt(i+1,j) > Grd%kmt(i,j)) topog_step(i,j,1)=1.0
           if(Grd%kmt(i,j+1) > Grd%kmt(i,j)) topog_step(i,j,2)=1.0
           if(Grd%kmt(i-1,j) > Grd%kmt(i,j)) topog_step(i,j,3)=1.0
           if(Grd%kmt(i,j-1) > Grd%kmt(i,j)) topog_step(i,j,4)=1.0
        endif
     enddo
  enddo

  ! Block out the bipolar fold in order to ensure tracer 
  ! conservation.
  if(jec+Dom%joff==Dom%jeg) topog_step(:,jec,:) = 0.0

  do j=jsc,jec
     do i=isc,iec
        m=1 ; topog_slope(i,j,m) = slope_x(i,j)  
        m=2 ; topog_slope(i,j,m) = slope_y(i,j)  
        m=3 ; topog_slope(i,j,m) = slope_x(i-1,j)
        m=4 ; topog_slope(i,j,m) = slope_y(i,j-1)
     enddo
  enddo

  call mpp_update_domains(topog_step(:,:,:), Dom%domain2d)
  call mpp_update_domains(topog_slope(:,:,:),Dom%domain2d)

  ! Allocate the buffers
  rea_buff_size = num_prog_tracers+1
  call allocate_buffer(Ebuffer_out, rea_buff_size)
  call allocate_buffer(Ebuffer_in, rea_buff_size)
  
  call allocate_buffer(Wbuffer_out, rea_buff_size)
  call allocate_buffer(Wbuffer_in, rea_buff_size)

  call allocate_buffer(Nbuffer_out, rea_buff_size)
  call allocate_buffer(Nbuffer_in, rea_buff_size)

  call allocate_buffer(Sbuffer_out, rea_buff_size)
  call allocate_buffer(Sbuffer_in, rea_buff_size)

end subroutine blob_static_bottom_init
! </SUBROUTINE>  NAME="ocean_blob_static_bottom_init"


!######################################################################
! <SUBROUTINE NAME="blob_overflow_like">
!
! <DESCRIPTION>
! Run the Lagrangian blob model for the Overflow schemes.  These 
! schemes are static schemes and uses Lagrangian blobs to transport 
! tracer and mass down slopes where an onshelf/offshelf instability
! exists.
!
! This scheme is activated by setting the namelist variable 
! blob_overflow=.true.
!
! There are three flavours to this scheme, which are detailed in
! sections 2.4 and 2.5 of Bates et al. for further details.  The three
! flavours are:
! 1/ A Campin and Goosse (1999) scheme in which the full "plumbing" of
!    water going off shelf to the deep ocean and deep ocean waters
!    returning on shelf is specified (described in full in section 2.4
!    of Bates et al.)
! 2/ Only the lateral part of the plumbing is specified and all vertical
!    movement of water within the deep ocean column is taken care of 
!    (described in full in section 2.5 of Bates et al.).  This option is
!    activated by setting the namelist variable overflow_one_return=.true.
! 3/ Only the movement of shelf water to the deep ocean column is
!    explicitly dealt with (described in full in section 2.5 of Bates
!    et al.).  This option is activated by setting the namelist variable
!    overflow_no_return=.true.
! 
! Other namelist variables associated with this scheme are:
!  blob_overflow_mu (real), which is the coefficient of friction used
!    used to calculate the overflow velocity,
!  blob_overflow_delta (real), which is the fraction of a grid cell that
!    participates in any one overflow event,
!  blob_overflow_umax (real), which is the maximum overflow velocity.
!    Depending on the flavour of scheme chosen and the impact that you want
!    it to have, typical values should be O(0.01) to O(1).
! </DESCRIPTION>
!
subroutine blob_overflow_like(Time, Thickness, T_prog, Dens, head, blob_counter)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_blob_type),        pointer       :: head 
  integer, dimension(isc:iec,jsc:jec,nk), intent(inout) :: blob_counter

  !local variables
  real    :: temp_so, salt_so, press, so_dens, blob_mass
  real    :: overflow_speed, overflow_thickness, overflow_flux
  integer :: i,j,k,kdw,kup,m,n,tau,iip,jjq
  integer :: nE, nW, nN, nS
  integer :: stderrunit
  logical :: sink
  type(ocean_blob_type), pointer :: prev => NULL()
  type(ocean_blob_type), pointer :: this => NULL()
  type(ocean_blob_type), pointer :: next => NULL()

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
         '==>Error in ocean_blob_static_bottom_mod (Lagrangian blob model): '&
         //'module needs to be initialized')
  endif 

  if (.not. use_this_module) return

  stderrunit = stderr()
  tau   = Time%tau
  kup   = 0
  kdw   = 0

  ! search for horizontal instability, including in certain parts of the halo
  do m=1,4
     do j=jsh(m),jeh(m)
        do i=ish(m),ieh(m)

           ! check whether downslope flow is possible
           if (topog_step(i,j,m) == 1.0) then

              ! some convenient variables
              iip = i+ip(m)
              jjq = j+jq(m)
              kup = Grd%kmt(i,j)

              ! overflow _flux and _speed are used to calculate the mass that 
              ! is overflowing.
              overflow_flux  = 0.0 !overflow flux kg/s
              overflow_speed = 0.0 !overflow speed m/s

              ! check whether density of shelf water is more than deep water
              if(Dens%rho(i,j,kup,tau) > Dens%rho(iip,jjq,kup,tau)) then
                 
                 ! calculate the overflow velocity
                 overflow_speed = overflow_factor*topog_slope(i,j,m) &
                                  *( Dens%rho(i,j,kup,tau)-Dens%rho(iip,jjq,kup,tau) )
                 overflow_speed = min(overflow_speed, blob_overflow_umax)

                 ! find the neutral level, or the lowest ocean grid box
                 salt_so = T_prog(index_salt)%field(i,j,kup,tau)
                 temp_so = T_prog(index_temp)%field(i,j,kup,tau)
                 kdw = kup
                 do k=kup+1,Grd%kmt(iip,jjq)
                    press = Dens%pressure_at_depth(iip,jjq,k)
                    so_dens = density(salt_so,temp_so,press)
                    ! compare the density of the shelf water to the 
                    ! deep water density referenced to the pressure of
                    ! depth of the deep water
                    if (so_dens > Dens%rho(iip,jjq,k,tau)) then
                       kdw = k
                    else
                       exit
                    endif
                 enddo

                 ! compute m^2/sec flux leaving vertical face of cell  
                 if(m==1 .or. m==3) then 
                    overflow_flux = overflow_speed*Grd%dyt(i,j)
                 elseif(m==2 .or. m==4) then 
                    overflow_flux = overflow_speed*Grd%dxt(i,j)
                 endif

                 ! based on 1D diffusive stability criteria, no more than 
                 ! one-half the mass of a grid cell can be mixed per time
                 ! step, which means:
                 ! overflow_flux(m^2/sec)*rho_dzt*dtime < 0.5*rho_dzt*dat,
                 ! or equivalently overflow_flux(m^2/sec) < 0.5*dat/dtime
                 overflow_flux = min(overflow_flux, max_overflow_flux(i,j))
       
                 ! now convert overflow_flux from m^2/sec to kg/sec
                 ! using rho_dzt 
                 overflow_thickness = Thickness%rho_dzt(i,j,kup,tau)
    
                 do k=kup,kdw
                    overflow_thickness = min(overflow_thickness,&
                                         Thickness%rho_dzt(iip,jjq,k,tau))
                 enddo

                 ! overflow flux is (m^2/s).  In order to get mass (kg), we
                 ! need to multiply by overflow_thickness(kg/m^2)*dtime(s)
                 blob_mass = overflow_flux*overflow_thickness*dtime

                 ! Form the onshelf blobs.  Make sure they are in the 
                 ! compute domain.
                 if (i>=isc .and. i<=iec .and. j>=jsc .and. j<=jec) then
                    call formoverflowblob(i,j,kup,.true.)
                 endif

                 ! Form the deep ocean blobs.  Make sure they are in the compute
                 ! domain
                 if (iip>=isc .and. iip<=iec .and. jjq>=jsc .and. jjq<=jec) then
                    if(overflow_one_return) then
                       call formoverflowblob(iip,jjq,kup,.false.)
                    elseif(.not. overflow_no_return) then
                       do k=kup,kdw
                          call formoverflowblob(iip,jjq,k,.false.)
                       enddo
                    endif
                 endif
              endif !rho_so>rho_do?
           endif !topog_step==1.0?
        enddo !i
     enddo !j
  enddo !m

  ! The next portion of code performs the following operations:
  ! 1/ move the blobs
  ! 2/ any blobs that have gone outside the compute domain need to:
  ! 2a/ put into a buffer and deleted from this PE's blob list
  ! 2b/ send the buffer to the adjacent compute domain, N,S,E or W
  ! 2c/ receive buffers from adjacent compute domains
  ! 2d/ add the blobs contained in the received bufferto this PE's blob
  !     list
  ! 3/ do some more diagnostics 
  ! 4/ transfer the properties from the Lagrangian model to the Eulerian
  !    model
  ! MPI does not handle passing of objects from one PE to another.  So,
  ! we need to "pack" the blob into a "buffer", which is an array of
  ! properties.  We then pass the array to another PE, where the array is
  ! read, and objects are created on that PE.

  ! clear the outward buffers & accummulators
  call clear_buffer(Ebuffer_out);  nE  = 0
  call clear_buffer(Wbuffer_out);  nW  = 0
  call clear_buffer(Nbuffer_out);  nN  = 0
  call clear_buffer(Ebuffer_out);  nS  = 0

  if(associated(head)) then
     prev=>NULL()
     this=>head
     next=>head%next

     blobmove: do
        ! for convenience
        i    = this%i - Dom%ioff
        j    = this%j - Dom%joff
        k    = this%k
        m    = this%m
        kdw  = this%kdw
        kup  = this%kup
        sink = this%sink

        if (sink) then
           ! Move the blobs:
           ! if sink is .true. then the on shelf blob moves to kdw in the
           ! deep ocean column
           i = i + ip(m)
           j = j + jq(m)
           k = kdw
        else
           ! if sink is not .true. the the blob must be in the deep
           ! ocean column
           if (k==kup) then
              ! if the blob is in the deep ocean at kup, then move it on
              ! to the shelf, staying at kup
              i = i - ip(m)
              j = j - jq(m)
           else
              ! if the blob is not at kup, then keep it in the deep ocean
              ! column and move it up one grid box
              k = k - 1
           endif
        endif ! sink?

        this%i = i + Dom%ioff
        this%j = j + Dom%joff
        this%k = k

        ! test if this has gone out to the E,W,N or S of the 
        ! compute domain.  If it has, then pack it into an outward
        ! buffer and remove it from this PE's list of blobs
        if (i>iec) then
           nE = nE + 1
           call addtobuffer(Ebuffer_out, nE, Info%pe_E)
        elseif (i<isc) then
           nW = nW + 1
           call addtobuffer(Wbuffer_out, nW, Info%pe_W)
        elseif (j>jec) then
           nN = nN + 1
           call addtobuffer(Nbuffer_out, nN, Info%pe_N)
        elseif (j<jsc) then
           nS = nS + 1
           call addtobuffer(Sbuffer_out, nS, Info%pe_S)
        endif

        ! move down the list.  If we're at the end of the list, exit.
        if(associated(this)) prev=>this
        this=>next
        if(associated(next)) next=>next%next

        if(.not.associated(this)) then
           nullify(prev)
           nullify(this)
           nullify(next)
           exit blobmove
        endif

     enddo blobmove
  endif !head associated

  ! Tell the buffer how many blobs it contains from the accummulators
  Ebuffer_out%numblobs  = nE
  Wbuffer_out%numblobs  = nW
  Nbuffer_out%numblobs  = nN
  Sbuffer_out%numblobs  = nS
  
  ! Now, send blobs to neighbouring PE's.  Note: since blobs can only move
  ! N,S,E or W, we do not need to worry about diagonal PE's 
  ! (i.e. NE,SE,SW,NW)
  if (Info%pe_E .ne. NULL_PE) call send_buffer(Ebuffer_out, Info%pe_E, rea_buff_size)
  if (Info%pe_W .ne. NULL_PE) call send_buffer(Wbuffer_out, Info%pe_W, rea_buff_size)
  if (Info%pe_N .ne. NULL_PE) call send_buffer(Nbuffer_out, Info%pe_N, rea_buff_size)
  if (Info%pe_S .ne. NULL_PE) call send_buffer(Sbuffer_out, Info%pe_S, rea_buff_size)
  
  ! Clear the incoming buffers
  call clear_buffer(Wbuffer_in)
  call clear_buffer(Ebuffer_in)
  call clear_buffer(Sbuffer_in)
  call clear_buffer(Nbuffer_in)
  
  ! Now receive blobs from neighbouring PE's
  if (Info%pe_W .ne. NULL_PE) call receive_buffer(Wbuffer_in, Info%pe_W, rea_buff_size)
  if (Info%pe_E .ne. NULL_PE) call receive_buffer(Ebuffer_in, Info%pe_E, rea_buff_size)
  if (Info%pe_S .ne. NULL_PE) call receive_buffer(Sbuffer_in, Info%pe_S, rea_buff_size)
  if (Info%pe_N .ne. NULL_PE) call receive_buffer(Nbuffer_in, Info%pe_N, rea_buff_size)

  ! Now unpack the buffer and add the received blobs to this PE's blob list
  if (Info%pe_W .ne. NULL_PE) call unpackbuffer(Wbuffer_in)
  if (Info%pe_E .ne. NULL_PE) call unpackbuffer(Ebuffer_in)
  if (Info%pe_S .ne. NULL_PE) call unpackbuffer(Sbuffer_in)
  if (Info%pe_N .ne. NULL_PE) call unpackbuffer(Nbuffer_in)

  call mpp_sync_self()
  
  if (associated(head)) then
     prev=>NULL()
     this=>head
     next=>this%next

     blob_diag: do
        ! for convenience
        i         = this%i - Dom%ioff
        j         = this%j - Dom%joff
        k         = this%k
        m         = this%m
        kup       = this%kup
        kdw       = this%kdw
        sink      = this%sink
        blob_mass = this%mass

        if (i>=isc .and. i<=iec .and. j>=jsc .and. j<= jec) then
           ! now transfer mass from the Lagrangian model to the Eulerian model
           Thickness%mass_source(i,j,k) = Thickness%mass_source(i,j,k) &
                                          + this%mass*datdtime_r(i,j)
           this%mass = this%mass - blob_mass

           ! now transfer tracer from the Lagrangian system to the Eulerian System
           do n=1,num_prog_tracers
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) &
                                             + this%tracer(n)*datdtime_r(i,j)
              this%tracer(n) = this%tracer(n) - this%tracer(n)
           enddo

           ! move down the list
           prev=>this
           this=>next
           if(associated(next)) next=>next%next

           ! if we are at the end of the list, then exit the loop
           if(.not.associated(this)) exit blob_diag
        else
           write(stderrunit,*) 'blob not on PE when it should be! blob coordinates i=',i,'j=',j,'pe=',Info%pe_this
           call error_mesg( 'ocean_blob_static_bottom_mod, '&
                //'subroutine blob_overflow_like', 'blob outside domain',FATAL)
        endif
     enddo blob_diag
  endif !head associated?

  ! now delete all approximately zero mass blobs
  ! (which should be all of them in the overflow_like_scheme)
  call blob_delete(Time, Thickness, T_prog(:), head)

  ! While there should be identically zero nett change in the mass source
  ! we need to update the mass_source in order to get bitwise agreement
  ! across different numbers of processors.
  call mpp_update_domains(Thickness%mass_source(:,:,:),Dom%domain2d)

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that forms the overflow blobs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine formoverflowblob(this_i, this_j, this_k, this_sink)

    integer :: this_i, this_j, this_k
    logical :: this_sink

    !local variables
    type(ocean_blob_type), pointer :: new_blob => NULL()

    allocate(new_blob)
    allocate(new_blob%tracer(num_prog_tracers))

    ! allocate important position integers to the blob
    new_blob%i   = this_i + Dom%ioff
    new_blob%j   = this_j + Dom%joff
    new_blob%k   = this_k
    new_blob%m   = m
    new_blob%kup = kup
    new_blob%kdw = kdw

    call count_blob(new_blob, blob_counter)

    ! .true. indicates the onshelf blob 
    ! .false. indicates a deep ocean blob
    new_blob%sink = this_sink

    ! now transfer the properties from the Eulerian model to
    ! the Lagrangian model.  mass_source is in units 
    ! (kg/m^3)*(m/s)=kg/(m^2*s)
    Thickness%mass_source(this_i,this_j,this_k) =    &
         Thickness%mass_source(this_i,this_j,this_k) &
         - blob_mass*datdtime_r(this_i,this_j)
    new_blob%mass      = blob_mass

    do n=1,num_prog_tracers
       ! transfer the tracer from the Eulerian system to the Lagrangian 
       ! system. field is tracer concentration. th_tendency is 
       ! kg*[concentration]/(m^2 sec).  
       new_blob%tracer(n) = new_blob%mass&
                            *T_prog(n)%field(this_i,this_j,this_k,tau)
       T_prog(n)%th_tendency(this_i,this_j,this_k) =    &
            T_prog(n)%th_tendency(this_i,this_j,this_k) &
            - new_blob%tracer(n)*datdtime_r(this_i,this_j)
    enddo

    ! insert the new blob into the list
    call insert_blob(new_blob,head)

    !now disassociate the new_blob pointer from the new_blob object
    nullify(new_blob)
 
  end subroutine formoverflowblob

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that packs a blob into a buffer          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine addtobuffer(buffer, nblob, pe)

    type(blob_buffer_type), pointer :: buffer
    integer, intent(in) :: nblob, pe

    ! local variables
    integer :: pt !counter for prognostic tracers

    if (pe .eq. NULL_PE) then
       write(stderrunit,*) 'this pe =',Info%pe_this, 'other pe =',pe
       call error_mesg('ocean_blobs_mod, subroutine addtobuffer', &
            'trying to pack a buffer to be sent to a NULL_PE',FATAL)
    endif

    if (nblob>buffer%size) call increase_buffer(buffer, rea_buff_size)

    buffer%numblobs = nblob

    ! Pack all the integer data
    buffer%integer_buff(1, nblob) = this%i
    buffer%integer_buff(2, nblob) = this%j
    buffer%integer_buff(3, nblob) = this%k
    buffer%integer_buff(4, nblob) = this%m
    buffer%integer_buff(5, nblob) = this%kdw
    buffer%integer_buff(6, nblob) = this%kup
    buffer%integer_buff(7, nblob) = this%hash
    buffer%integer_buff(8, nblob) = this%number
    
    ! Pack all the real data
    do pt=1,num_prog_tracers
       buffer%real_buff(pt,nblob) = this%tracer(pt)
    enddo
    buffer%real_buff(num_prog_tracers+1,nblob) = this%mass
       
    ! Pack all the logical data
    buffer%logical_buff(1,nblob) = this%sink

    ! Now take the blob out of the list and erase from memory
    if(associated(prev)) then
       prev%next=>next
    else
       head=>next
    endif
    if(associated(next)) next%prev=>prev

!!$    deallocate(this)
!!$    nullify(this)
    call free_blob_memory(this)

  end subroutine addtobuffer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that unpacks an entire buffer            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine unpackbuffer(buffer)

    type(blob_buffer_type), pointer :: buffer

    ! local variables
    type(ocean_blob_type),pointer :: new_blob=>NULL()
    integer :: pt !counter for prognostic tracers
    integer :: nb !counter for blobs

    if (buffer%numblobs>0) then
       do nb=1,buffer%numblobs

          allocate(new_blob)
          allocate(new_blob%tracer(num_prog_tracers))

          ! Unpack all the integer data
          new_blob%i      = buffer%integer_buff(1 ,nb)
          new_blob%j      = buffer%integer_buff(2 ,nb)
          new_blob%k      = buffer%integer_buff(3 ,nb)
          new_blob%m      = buffer%integer_buff(4 ,nb)
          new_blob%kdw    = buffer%integer_buff(5 ,nb)
          new_blob%kup    = buffer%integer_buff(6 ,nb)
          new_blob%hash   = buffer%integer_buff(7 ,nb)
          new_blob%number = buffer%integer_buff(8 ,nb)

          ! Unpack all the real data
          do pt=1,num_prog_tracers
             new_blob%tracer(pt)      = buffer%real_buff(pt,nb)
          enddo
          new_blob%mass      = buffer%real_buff(num_prog_tracers+1,nb)

          ! Unpack all the logical data
          new_blob%sink = buffer%logical_buff(1,nb)

          ! Add the blob to the list
          call insert_blob(new_blob, head)

          nullify(new_blob)
          
       enddo
    endif

  end subroutine unpackbuffer

end subroutine blob_overflow_like
! </SUBROUTINE>  NAME="blob_overflow_like"

!######################################################################
! <SUBROUTINE NAME="blob_static_bottom_end">
!
! <DESCRIPTION>
! Does what is necessary to finish the run.
! </DESCRIPTION>
!
subroutine blob_static_bottom_end()

  if (.not. use_this_module) return

  call deallocate_buffer(Ebuffer_out);  call deallocate_buffer(Ebuffer_in)
  call deallocate_buffer(Wbuffer_out);  call deallocate_buffer(Wbuffer_in)
  call deallocate_buffer(Nbuffer_out);  call deallocate_buffer(Nbuffer_in)
  call deallocate_buffer(Sbuffer_out);  call deallocate_buffer(Sbuffer_in)

  nullify(Dom)
  nullify(Grd)

end subroutine blob_static_bottom_end
! </SUBROUTINE>  NAME="blob_static_bottom_end"

!######################################################################
! <SUBROUTINE NAME="allocate_buffer">
!
! <DESCRIPTION>
! Increases the buffer size for sending blobs from one PE to another.
! </DESCRIPTION>
!
subroutine allocate_buffer(buffer, rea_buff_params)
  
  type(blob_buffer_type), pointer :: buffer
  integer, intent(in) :: rea_buff_params
  
  allocate(buffer)
  allocate(buffer%integer_buff(int_buff_params, delta_buffer))
  allocate(buffer%real_buff(rea_buff_params,    delta_buffer))
  allocate(buffer%logical_buff(log_buff_params, delta_buffer))
  
  buffer%size = delta_buffer
  
end subroutine allocate_buffer

! </SUBROUTINE>  NAME="allocate_buffer"

!######################################################################
! <SUBROUTINE NAME="increase_buffer">
!
! <DESCRIPTION>
! Increases the buffer size for sending blobs from one PE to another.
! </DESCRIPTION>
!
subroutine increase_buffer(buffer, rea_buff_params)

  type(blob_buffer_type), pointer :: buffer
  integer, intent(in) :: rea_buff_params

  ! local variables
  type(blob_buffer_type), pointer :: new_buffer
  integer :: newbuffsize

  allocate(new_buffer)
  newbuffsize = buffer%size + delta_buffer

  allocate(new_buffer%integer_buff(1:int_buff_params, 1:newbuffsize))
  allocate(new_buffer%real_buff(   1:rea_buff_params, 1:newbuffsize))
  allocate(new_buffer%logical_buff(1:log_buff_params, 1:newbuffsize))
  
  new_buffer%integer_buff(:,1:buffer%size) = buffer%integer_buff(:,1:buffer%size)
  new_buffer%real_buff(   :,1:buffer%size) = buffer%real_buff(   :,1:buffer%size)
  new_buffer%logical_buff(:,1:buffer%size) = buffer%logical_buff(:,1:buffer%size)

  new_buffer%size = newbuffsize

  deallocate(buffer%integer_buff)
  deallocate(buffer%real_buff)
  deallocate(buffer%logical_buff)
  deallocate(buffer)

  nullify(buffer)
  buffer=>new_buffer

end subroutine increase_buffer
! </SUBROUTINE>  NAME="increase_buffer"

!######################################################################
! <SUBROUTINE NAME="send_buffer">
!
! <DESCRIPTION>
! Sends a buffer to an adjoining PE
! </DESCRIPTION>
!
subroutine send_buffer(buffer, pe, rea_buff_params)

  type(blob_buffer_type), pointer :: buffer
  integer, intent(in) :: pe
  integer, intent(in) :: rea_buff_params

  call mpp_send(buffer%numblobs, plen=1, to_pe=pe)
  if (buffer%numblobs>0) then
     call mpp_send(buffer%integer_buff, buffer%numblobs*int_buff_params, pe)
     call mpp_send(buffer%real_buff,    buffer%numblobs*rea_buff_params, pe)
     call mpp_send(buffer%logical_buff, buffer%numblobs*log_buff_params, pe)
  endif

end subroutine send_buffer
! </SUBROUTINE>  NAME="send_buffer"



!######################################################################
! <SUBROUTINE NAME="receive_buffer">
!
! <DESCRIPTION>
! Receives a buffer from an adjoining PE
! </DESCRIPTION>
!
subroutine receive_buffer(buffer, pe, rea_buff_params)

  type(blob_buffer_type), pointer :: buffer
  integer, intent(in) :: pe
  integer, intent(in) :: rea_buff_params

  !local variables
  integer :: incoming

  call mpp_recv(incoming, glen=1, from_pe=pe)
  if (incoming>0) then
     if (incoming>buffer%size) call increase_buffer(buffer, rea_buff_params)
     call mpp_recv(buffer%integer_buff, incoming*int_buff_params, pe)
     call mpp_recv(buffer%real_buff,    incoming*rea_buff_params, pe)
     call mpp_recv(buffer%logical_buff, incoming*log_buff_params, pe)
  endif
  buffer%numblobs = incoming

end subroutine receive_buffer
! </SUBROUTINE>  NAME="receive_buffer"


!######################################################################
! <SUBROUTINE NAME="clear_buffer">
!
! <DESCRIPTION>
! Clears the contents of a buffer
! </DESCRIPTION>
!
subroutine clear_buffer(buffer)
    
  type(blob_buffer_type), pointer :: buffer

  buffer%numblobs          = 0
  buffer%integer_buff(:,:) = 0
  buffer%real_buff(:,:)    = 0.0
  buffer%logical_buff(:,:) = .true.
  
  !note, we do not clear buffer%size
  
end subroutine clear_buffer
! </SUBROUTINE>  NAME="clear_buffer"

!######################################################################
! <SUBROUTINE NAME="deallocate_buffer">
!
! <DESCRIPTION>
! Deallocates memory from a buffer (usually at the end of a run).
! </DESCRIPTION>
!
subroutine deallocate_buffer(buffer)
  
  type(blob_buffer_type), pointer :: buffer
  
  deallocate(buffer%logical_buff)
  deallocate(buffer%real_buff)
  deallocate(buffer%integer_buff)
  deallocate(buffer)
  
end subroutine deallocate_buffer
! </SUBROUTINE>  NAME="deallocate_buffer"

end module ocean_blob_static_bottom_mod

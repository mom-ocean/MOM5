module ocean_blob_static_free_mod
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module controls and runs the free static Lagrangian blob 
! parameterisations.
!</OVERVIEW>
!
!<DESCRIPTION>
! There are three available static blob schemes.  None are sanctioned
! for their physical integrity, but they have been important in the 
! development and testing of the Lagrangian framework.
!
! The first free static scheme emulates the NCON scheme of Cox (1984).
! The second is a scheme that acts in a diffusive manner, while the third
! swaps the properties of adjacent grid cells.  Details of all three schemes
! can be found in Bates et al. (2010).
!</DESCRIPTION>
!
!<INFO>
!
! <REFERENCE>
! S.M. Griffies, Elements of mom4p1 (2009)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! Cox, M.D., 1984, A Primitive Equation, 3-Dimensional Model of the
! Ocean.  NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
!</INFO>
!<NAMELIST NAME="ocean_blob_static_free_nml">
!
!  <DATA NAME="blob_ncon_like" TYPE="logical">
!  If true, will use NCon-like formulation. 
!  Default blob_ncon_like = .true. 
!  </DATA> 
!
!  <DATA NAME="blob_diff_like" TYPE="logical">
!  If true, will use the diffusion-like formulation. 
!  Default is "blob_diff_like=.false.
!  </DATA>
!
!  <DATA NAME="blob_swap_like" TYPE="logical">
!  If true, will use the swap-like formulation.  
!  Default is blob_switch_like=.false.
!  </DATA>
!
!  <DATA NAME="ncon_blob" TYPE="integer">
!  The number of times that the water column is checked and adjusted
!  for instability.
!  </DATA>
!</NAMELIST>
!
use fms_mod, only: open_namelist_file, check_nml_error, stdout, stdlog, FATAL, mpp_error
use fms_mod, only: close_file, mpp_error

use ocean_blob_util_mod, only: blob_delete
use ocean_density_mod,   only: density
use ocean_types_mod,     only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_density_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_blob_type
use ocean_workspace_mod, only: wrk1

implicit none

private

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

logical :: module_is_initialized=.false.

real, dimension(:,:,:), allocatable :: ocean_mask
integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1
integer :: isd,ied,jsd,jed,isc,iec,jsc,jec,nk

public blob_static_free_init
public blob_static_free
public blob_static_free_end

!namelist defaults
logical :: use_this_module = .false.
logical :: blob_ncon_like  = .false. ! for NCON-like scheme.
logical :: blob_diff_like  = .false. ! for Diffusion-like scheme.
logical :: blob_swap_like  = .false. ! for Swap-like scheme.
integer :: ncon_blob       = 7       ! equivalent to traditional ncon

namelist /ocean_blob_static_free_nml/ use_this_module, blob_ncon_like, blob_diff_like, &
                                      blob_swap_like,ncon_blob

contains

!#######################################################################
! <SUBROUTINE NAME="blob_static_free_init">
!
! <DESCRIPTION>
! Initialises the free static schemes.
! </DESCRIPTION>
!
subroutine blob_static_free_init(Grid, Domain, num_tracers, itemp, isalt)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  integer, intent(in) :: num_tracers
  integer, intent(in) :: itemp
  integer, intent(in) :: isalt

  integer :: k
  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit 

  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
         '==>Error in ocean_blob_static_free_mod (ocean_blob_static_free_init):' &
         //' module already initialized')
  endif 

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_static_free_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_static_free_nml)  
  write (stdlogunit,ocean_blob_static_free_nml)
  ierr = check_nml_error(io_status,'ocean_blob_static_free_nml')
  call close_file (ioun)

  module_is_initialized=.true.

  if(.not. use_this_module) return

  if(blob_ncon_like) then 
     write(stdoutunit,'(a,i2)') & 
     ' ==>Note: Using the static Lagrangian blob NCon-like convection scheme '&
     //'equivalent to ncon=',ncon_blob
  endif

  if(blob_diff_like) then
          write(stdoutunit,'(a)') &
     ' ==>Note: Using the static Lagrangian Diffusion like convection scheme '
  endif

  if(blob_swap_like) then
          write(stdoutunit,'(a)') &
     ' ==>Note: Using the static Lagrangian Swap convection scheme '
  endif

  index_temp       = itemp
  index_salt       = isalt
  num_prog_tracers = num_tracers

  Grd => Grid
  Dom => Domain

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  nk = Grd%nk

  ! Makes an ocean/land mask where land values are k*1e30 and
  ! ocean values are 0.0
  allocate(ocean_mask(isd:ied,jsd:jed,nk))
  do k=1,nk
     ocean_mask(:,:,k) = -real(k)*1.e30*(Grd%tmask(:,:,k) - 1)
  enddo

end subroutine blob_static_free_init
! </SUBROUTINE>  NAME="blob_static_free_init"


!#######################################################################
! <SUBROUTINE NAME="blob_static_free">
!
! <DESCRIPTION>
! A subroutine that is called by the main driver blob module, 
! ocean_blob_mod.  The present module in turn calls the three available
! static schemes.  They are not mutually exclusive and so may be run
! in any combination.
!
! Note that only the NCON-like scheme has been implemented.
! </DESCRIPTION>
!
subroutine blob_static_free(Time, Thickness, T_prog, Dens, head)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_blob_type),        pointer       :: head

  if (.not. use_this_module) return

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
         '==>Error in ocean_blob_static_free_mod (Static free blob model): module '&
         //' needs to be initialized')
  endif 

  if(blob_ncon_like) call blob_ncon_like_scheme(Time, Thickness, T_prog(:), Dens, head)

  if(blob_diff_like) then
     call mpp_error(FATAL, &
      '==>Error in ocean_blobs_mod (Lagrangian blob model): '&
      //'blob_diff_like=.true. not yet implemented.') 
  endif

  if(blob_swap_like) then
     call mpp_error(FATAL, &
      '==>Error in ocean_blobs_mod (Lagrangian blob model): '&
      //'blob_swap_like=.true. not yet implemented.') 
  endif

end subroutine blob_static_free
! </SUBROUTINE>  NAME="blob_static_free"


!######################################################################
! <SUBROUTINE NAME="blob_ncon_like_scheme">
!
! <DESCRIPTION>
! Run the Lagrangian blob model for the NCon-like scheme.  This scheme
! is a static scheme and uses Lagrangian blobs to homogenise adjacent 
! vertical grid cells by by transferring mass and tracer between them.
!
! This scheme is activated by setting the namelist variable 
! blob_ncon_like=.true.  The number of times that the domain is checked
! for instability is set by the integer namelist variable ncon_blob, of 
! which the default value is 7.
!
! It is based on the original convective adjustment scheme of 
! Cox (1984).
!
! See section 2.1 of Bates et al. for further details.
! </DESCRIPTION>
!
subroutine blob_ncon_like_scheme(Time, Thickness, T_prog, Dens, head)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_blob_type),        pointer       :: head

  real    :: rhodztk, rhodztkp1, dat, rho_dzt_dat_r
  real    :: mass_blob
  integer :: i, j, k, n, nn, taup1, l, ks, kp1
  type(ocean_blob_type), pointer :: this => NULL()

  taup1 = Time%taup1

  ! ks=1: compare levels 1 to 2; 3 to 4; etc.
  ! ks=2: compare levels 2 to 3; 4 to 5; etc.
  do nn = 1,ncon_blob
     do ks = 1,2

        ! find density with reference pressure determined by vertical pairs 
        do j=jsc,jec
           do l=1,nk,2
              if (ks .eq. 1) then
                 k = min(l+1,nk)
              else
                 k = l
              endif
              do i=isc,iec
                 wrk1(i,j,l)  =density(T_prog(index_salt)%field(i,j,l,taup1),&
                                       T_prog(index_temp)%field(i,j,l,taup1),&
                                       Dens%pressure_at_depth(i,j,k))
              enddo
           enddo
        enddo
        
        do j=jsc,jec
           do l=2,nk,2
              if (ks .eq. 1) then
                 k = l
              else
                 k = min(l+1,nk)
              endif
              do i=isc,iec
                 wrk1(i,j,l)=density(T_prog(index_salt)%field(i,j,l,taup1),&
                                     T_prog(index_temp)%field(i,j,l,taup1),&
                                     Dens%pressure_at_depth(i,j,k))
              enddo
           enddo
        enddo

        ! set "heavy water" in land to stop convection
        wrk1(:,:,:) = wrk1(:,:,:) + ocean_mask(:,:,:)

        do k=ks,nk-1,2
           kp1 = k+1
           do j=jsc,jec
              do i=isc,iec
                 ! test for instability
                 if(wrk1(i,j,k)>wrk1(i,j,kp1)) then
                    ! some shorthand variables for convenience
                    rhodztk   = Thickness%rho_dzt(i,j,k,  taup1)
                    rhodztkp1 = Thickness%rho_dzt(i,j,kp1,taup1)
                    dat       = Grd%dat(i,j)
                    
                    ! compute the mass per unit area of the two blobs.  Note that 
                    ! dat is the same for both k and kp1;
                    ! mass=(dat*rhodztk * dat*rhodztkp1)/(dat(rhodztk + rhodztkp1))
                    !     =dat*rhodztk*rhodztkp1/(rhodztk + rhodztkp1)
                    mass_blob = dat*(rhodztk*rhodztkp1)/(rhodztk + rhodztkp1)

                    ! form the sinking blob
                    call nconblobform(k,  .true.)
                    
                    ! form the rising blob
                    call nconblobform(kp1,.false.)
                 endif
              enddo
           enddo
        enddo
  

        ! now that the blobs have been formed and their properties
        ! transferred from the Eulerian system to the Lagrangian system:
        ! 1/ move the blobs
        ! 2/ complete the diagnostics
        ! Note: see comments in blob_form for rationale of units and
        ! variables used for transferring properties between Eulerian and
        ! Lagrangian system.
        if(associated(head)) then
           this=>head

           blobmove: do
              ! for convenience and readability
              i = this%i
              j = this%j
              k = this%k
           
              ! now move the blobs
              if (this%sink) then
                 ! it is a sinking blob, so move it down one grid cell
                 k = k + 1
              else
                 ! it is a rising blob, so move it up one grid cell
                 k = k - 1
              endif
              this%k = k
              
              ! some shorthand variables for convenience
              rho_dzt_dat_r = Grd%datr(i,j)*Thickness%rho_dztr(i,j,k)
           
              ! transfer mass from the Lagrangian system to the Eulerian System
              Thickness%rho_dzt(i,j,k,taup1) = Thickness%rho_dzt(i,j,k,taup1)&
                                               + this%mass*Grd%datr(i,j)
              this%mass = 0.0
        
              ! transfer tracer from the Lagrangian system to the Eulerian System
              do n=1,num_prog_tracers
                 T_prog(n)%field(i,j,k,taup1) = T_prog(n)%field(i,j,k,taup1) &
                                                + this%tracer(n)*rho_dzt_dat_r
                 this%tracer(n) = 0.0
              enddo

              ! set pointer to the next blob
              this=>this%next

              !if there is no next blob, then exit the loop
              if(.not.associated(this)) exit blobmove
           enddo blobmove
        
           nullify(this)
        endif
        ! delete all approximately zero mass blobs (which should be all of 
        ! them in the ncon_like_scheme).                                   
        call blob_delete(Time, Thickness, T_prog(:), head)
     enddo !ks=1,2
  enddo !ncon_blob

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that forms the blobs for the             !
! ncon_like_scheme                                                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine nconblobform(this_k,this_sink)

    integer :: this_k
    logical :: this_sink
    
    ! local variables
    type(ocean_blob_type), pointer :: new_blob => NULL()

    ! allocate memory to the sinking blob object
    allocate(new_blob)
    allocate(new_blob%tracer(num_prog_tracers))

    ! allocate properties to the sinking blob
    new_blob%sink = this_sink

    new_blob%i = i
    new_blob%j = j
    new_blob%k = this_k

    ! transfer mass from the Eulerian system to the Lagrangian system
    ! blob%mass is in units (kg) & rho_dzt is in units (kg/m^3)*(m)
    Thickness%rho_dzt(i,j,this_k,taup1) = Thickness%rho_dzt(i,j,this_k,taup1)&
                                           - mass_blob*Grd%datr(i,j)
    new_blob%mass = mass_blob

    ! compute the tracer content of the blob and transfer 
    ! that tracer from the Eulerian system to the Lagrangian system
    rho_dzt_dat_r = Grd%datr(i,j)*Thickness%rho_dztr(i,j,this_k)
    do n=1,num_prog_tracers
       ! new_blob%tracer is blob_mass*[concentration]
       new_blob%tracer(n) = new_blob%mass*T_prog(n)%field(i,j,this_k,taup1)
       ! this calculates the reduction in concentration of a grid cell
       T_prog(n)%field(i,j,this_k,taup1) = T_prog(n)%field(i,j,this_k,taup1) &
                            - new_blob%tracer(n)*rho_dzt_dat_r
    enddo

    ! Add the new blob to the top (head) of the list.  Only one blob can exist in a grid cell
    ! at a time.  So, there is no requirement to order the blobs within the list (see the 
    ! subroutine insert_blob for more details about ordering)
    new_blob%next => head
    new_blob%prev => NULL()
    if(associated(head)) head%prev => new_blob
    head      => new_blob

    !now disassociate the new_blob pointer from the new blob object
    nullify(new_blob)
              
  end subroutine nconblobform


end subroutine blob_ncon_like_scheme
! </SUBROUTINE>  NAME="blob_ncon_like_scheme"


!######################################################################
! <SUBROUTINE NAME="blob_static_free_end">
!
! <DESCRIPTION>
! Ends the free static module.
! </DESCRIPTION>
!
subroutine blob_static_free_end()

  if (.not. use_this_module) return

  nullify(Dom)
  nullify(Grd)

end subroutine blob_static_free_end
! </SUBROUTINE>  NAME="blob_static_free_end"

end module ocean_blob_static_free_mod

module ocean_lap_tracer_mod
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted and density weighted time tendency 
! for tracer from lateral laplacian diffusion.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes lateral laplacian diffusion of a tracer.
! There are two main options.
!
! (1) The lateral tracer fluxes can be aligned with the z-coordinate surfaces,
! in which case the fluxes must be approximated if (i) we use non-geopotential
! vertical coordinates, (ii) next to partial bottom step topography.
! This form of the diffusion is not recommended since it can lead to 
! the creation of spurious extrema.  
!
! (2) The lateral tracer fluxes can be aligned surfaces of constant vertical 
! coordinate. In this case the fluxes are no longer strictly "horizontal."
! However, the operator is simpler and it ensures that no suprious
! extrema are created. It is for this reason that the simpler operator
! is preferred.   
!  
! The diffusivity used to determine the strength of the tendency can be 
! a general function of space yet it is constant in time.  A namelist 
! option exists that determines this diffusivity as a local function 
! of the grid spacing. 
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lap_tracer_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!
!  <DATA NAME="horz_z_diffuse" TYPE="logical">
!  To compute diffusion along surfaces of constant depth.
!  This operation must necessarily be approximate for the two 
!  cases (i) non-geopotential vertical coordinates, (2)
!  next to partial bottom step topography.  There are cases where
!  use of this operator can lead to spurious creation of extrema
!  due to truncation errors associated with the "slope" term.
!  For most cases where lateral diffusion is required, we 
!  will want it to be "diffusive" in the sense that no extrema are
!  created.  So the default is horz_z_diffuse=.false.
!  The option to use horz_z_diffuse=.true. is maintained for 
!  legacy purposes alone.  
!  </DATA> 
!  <DATA NAME="horz_s_diffuse" TYPE="logical">
!  To compute diffusion along surfaces of constant vertical s-coordinate. 
!  This operation is ensured of obtaining a smoothing operator
!  that does not create extrema. It is the default for this 
!  reason. 
!  </DATA> 
!
!  <DATA NAME="alap" UNITS="m^2/sec" TYPE="real">
!  This is the value for the space-time constant Laplacian diffusivity. 
!  </DATA> 
!
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the diffusivity is set according to a velocity scale times
!  the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model (MICOM).  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!
!  <DATA NAME="read_diffusivity_mask" TYPE="logical">
!  Allows for reading of a mask that to apply diffusivity
!  only in selected regions.
!  Default read_diffusivity_mask=.false.
!  </DATA>
!
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose writes during initialization 
!  </DATA> 
!</NAMELIST>
!
use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, register_static_field
use fms_mod,             only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,             only: FATAL, NOTE, stdout, stdlog, read_data
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_mod,             only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_obc_mod,        only: store_ocean_obc_tracer_flux, ocean_obc_enhance_diff_back
use ocean_operators_mod,  only: FMX, FMY, FDX_ZT, FDY_ZT, BDX_ET, BDY_NT
use ocean_operators_mod,  only: FDX_T, FDY_T
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, ocean_time_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_prog_tracer_type, ocean_options_type
use ocean_workspace_mod,  only: wrk1
use ocean_util_mod,       only: diagnose_3d

implicit none

private

public ocean_lap_tracer_init
public lap_tracer

real    :: alap              = 0.0e3  ! tracer diffusivity (m^2/sec)
real    :: vel_micom         = 0.0    ! constant velocity scale (m/s) for setting micom diffusivity  
logical :: tracer_mix_micom  =.false. ! if true, diffusivity made a function of the grid spacing
logical :: verbose_init      =.true.  ! for verbose writes during initialization 

! for diagnostics 
logical :: used
integer :: id_ah_laplacian=-1
integer, dimension(:), allocatable  :: id_xflux_diff
integer, dimension(:), allocatable  :: id_yflux_diff
integer, dimension(:), allocatable  :: id_h_diffuse

real, dimension(:,:,:), allocatable :: diff_cet ! diffusivity for eastern face of T cell (m^2/s)
real, dimension(:,:,:), allocatable :: diff_cnt ! diffusivity for northern face of T cell (m^2/s)
real, dimension(:,:,:), allocatable :: diffusivity_mask ! 3d mask to selectively apply diffusion

character(len=128) :: version=&
     '$Id: ocean_lap_tracer.F90,v 20.0 2013/12/14 00:14:20 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

real, dimension(:,:,:), allocatable :: fx, fy
integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: num_prog_tracers      = 0

logical :: use_this_module       = .false.
logical :: horz_z_diffuse        = .false.
logical :: horz_s_diffuse        = .true.
logical :: read_diffusivity_mask = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc=.false.

namelist /ocean_lap_tracer_nml/ use_this_module, alap, tracer_mix_micom, vel_micom, verbose_init, &
                                read_diffusivity_mask, horz_z_diffuse, horz_s_diffuse

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_lap_tracer_init">
!
! <DESCRIPTION>
! Initialize the laplacian diffusion module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that diffusivity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_lap_tracer_init(Grid, Domain, Time, T_prog, Ocean_options, dtime, obc, tmask_sigma)

  type(ocean_grid_type),     target, intent(in)    :: Grid
  type(ocean_domain_type),   target, intent(in)    :: Domain
  type(ocean_time_type),             intent(in)    :: Time
  type(ocean_prog_tracer_type),      intent(inout) :: T_prog(:)
  type(ocean_options_type),          intent(inout) :: Ocean_options
  real,                              intent(in)    :: dtime
  logical,                           intent(in)    :: obc
  real,                    optional, intent(in)    :: tmask_sigma(Domain%isc:Domain%iec,Domain%jsc:Domain%jec)

  integer :: ioun, io_status, i, j, k, n, ierr
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_lap_tracer_init: module already initialized')
  endif 

  module_is_initialized = .TRUE.

  have_obc = obc

  num_prog_tracers = size(T_prog(:))

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_lap_tracer_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_lap_tracer_nml')
#else
  ioun  = open_namelist_file()
  read (ioun,ocean_lap_tracer_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_lap_tracer_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_lap_tracer_nml)  
  write (stdlogunit,ocean_lap_tracer_nml)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
 
  if(use_this_module) then 
      write(stdoutunit,*)'==>Note from ocean_lap_tracer_mod: Using this module.'
      Ocean_options%horz_lap_tracer = 'Used horizontal laplacian tracer mixing option.'
  else
      write(stdoutunit,*)'==>Note from ocean_lap_tracer_mod: NOT using this module.'
      Ocean_options%horz_lap_tracer = 'Did NOT use horizontal laplacian tracer mixing.'
      return 
  endif

  if(alap > 0 .or. (tracer_mix_micom .and. vel_micom > 0)) then 
      write (stdoutunit,'(/1x,a)') ' ==> Warning: laplacian diffusion generally contributes to unphysically'
      write (stdoutunit,'(7x,a)') 'large diapycnal mixing. Such may be damaging for longer simulations.'
      write (stdoutunit,'(7x,a)') 'It is recommended that the tracer diffusivity be set to zero, and one'
      write (stdoutunit,'(7x,a)') 'instead uses a neutral physics.'
  endif

  if(horz_z_diffuse) then 
      write(stdoutunit,*)'==>Warning: Computing laplacian tracer diffusion along constant geoptoential surfaces.'
      write(stdoutunit,*)'   This operator can create spurious extrema where s-surfaces deviate from z-surfaces.'
      write(stdoutunit,*)'   For this reason, the alternative horz_s_diffuse=.true. is recommended.'
  endif

  if(horz_s_diffuse) then 
      write(stdoutunit,*)'==>Note: Computing lateral tracer diffusion along constant s-surfaces.'
      write(stdoutunit,*)'   This operator will always smooth tracer fields and will not produce extrema.'
  endif

  if(horz_s_diffuse .and. horz_z_diffuse) then 
      write(stdoutunit,'(a)')'==>Error from ocean_lap_tracer_mod: use one of horz_s_diffuse and horz_z_diffuse.'
      write(stdoutunit,'(a)')'   Presently, both are set true and this is not allowed. '
      call mpp_error(FATAL,'==>Error from ocean_lap_tracer_mod: use one of horz_s_diffuse and horz_z_diffuse.')
  endif

  if(.not. horz_s_diffuse .and. .not. horz_z_diffuse) then 
      write(stdoutunit,*)'==>Warning from ocean_lap_tracer_mod: use_this_module=.true.'
      write(stdoutunit,*)'   However, both horz_s_diffuse and horz_z_diffuse are .false.'
      write(stdoutunit,*)'   So will not produce any lateral tracer diffusion.'
  endif


  write(stdoutunit,'(/a,f10.2)')'==> Note from ocean_lap_tracer_mod: using forward time step (secs)', dtime 

  Dom => Domain
  Grd => Grid

  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux', maskmap=Dom%maskmap)

  allocate (diff_cet(isd:ied,jsd:jed,nk))
  allocate (diff_cnt(isd:ied,jsd:jed,nk))
  allocate (fx(isd:ied,jsd:jed,nk))
  allocate (fy(isd:ied,jsd:jed,nk))  

  diff_cet(:,:,:) = alap
  diff_cnt(:,:,:) = alap
  fx = 0.0
  fy = 0.0
  
  ! Micom background diffusivity
  ! Space scale set by grid size
  ! Velocity scale input via namelist 
  if(tracer_mix_micom) then
    do k=1,nk
      diff_cet(:,:,k) = vel_micom*(2.0*Grd%dxt(:,:)*Grd%dyte(:,:)/(Grd%dxt(:,:)+Grd%dyte(:,:)))
      diff_cnt(:,:,k) = vel_micom*(2.0*Grd%dxtn(:,:)*Grd%dyt(:,:)/(Grd%dxtn(:,:)+Grd%dyt(:,:)))
    enddo
    if(verbose_init) then 
      do j=jsc,jec
        write (stdoutunit,'(a,i4,a,e14.7,a)') &
        ' MICOM Laplacian diffusivity at (isc,',j,',1) = ',diff_cet(isc,j,1),' m^2/s'
      enddo
    endif 
  endif 

  ! enhance background diffusion near open boundaries if needed
  if (have_obc) call ocean_obc_enhance_diff_back(diff_cet, diff_cnt)
  
  if(PRESENT(tmask_sigma)) then 
     do j=jsc,jec
        do i=isc,iec
           if(tmask_sigma(i,j) > 0.0) then   
              k = Grd%kmt(i,j)
              diff_cet(i,j,k) = diff_cet(i,j,k)*(1.0-tmask_sigma(i,j))
              diff_cnt(i,j,k) = diff_cnt(i,j,k)*(1.0-tmask_sigma(i,j))
           endif
        enddo
     enddo
  endif

  ! apply diffusivity in selected regions
  if(read_diffusivity_mask) then
      allocate (diffusivity_mask(isd:ied,jsd:jed,nk))
      call read_data('INPUT/tracer_diffusivity_mask.nc','diffusivity_mask', diffusivity_mask, Domain%domain2d)
      write (stdoutunit,*) '==>ocean_lap_tracer_mod: Completed read of 3d tracer diffusivity mask.'
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               diff_cet(i,j,k) = diff_cet(i,j,k)*diffusivity_mask(i,j,k)
               diff_cnt(i,j,k) = diff_cnt(i,j,k)*diffusivity_mask(i,j,k)
            enddo
         enddo
      enddo
  endif

  call mpp_update_domains (diff_cet(:,:,:), Dom%domain2d)
  call mpp_update_domains (diff_cnt(:,:,:), Dom%domain2d)

  id_ah_laplacian = register_static_field ('ocean_model', 'ah_laplacian',&
                    Grd%tracer_axes(1:3),'static laplacian diffusivity', &
                    'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e20/))
  call diagnose_3d(Time, Grd, id_ah_laplacian, diff_cet(:,:,:))

  ! register for diagnostics manager 
  allocate (id_xflux_diff(num_prog_tracers))
  allocate (id_yflux_diff(num_prog_tracers))
  allocate (id_h_diffuse(num_prog_tracers))
  id_xflux_diff = -1
  id_yflux_diff = -1
  id_h_diffuse  = -1

  do n=1,num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_xflux_diff(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'i-diffusive heat flux', 'Watt',          &
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_yflux_diff(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'j-diffusive heat flux', 'Watt',          &
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'horz-diffusion of heat', 'Watt/m^2',    &
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
     else
         id_xflux_diff(n) = register_diag_field ('ocean_model',                       &
         trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3),                   &
              Time%model_time, 'i-diffusive flux of '//trim(T_prog(n)%name), 'kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))
         id_yflux_diff(n) = register_diag_field ('ocean_model',                       &
         trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3),                   &
              Time%model_time, 'j-diffusive flux of '//trim(T_prog(n)%name), 'kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model',                            &
         trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3),                        &
              Time%model_time, 'horz-diffusion of '//trim(T_prog(n)%name), 'kg/(sec*m^2)',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
     endif

  enddo

  
end subroutine ocean_lap_tracer_init
! </SUBROUTINE>  NAME="ocean_lap_tracer_init"


!#######################################################################
! <SUBROUTINE NAME="lap_tracer">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from lateral laplacian diffusion. 
! </DESCRIPTION>
!
subroutine lap_tracer (Time, Thickness, Tracer, ntracer, diag_flag)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer,                      intent(in)    :: ntracer
  logical,            optional, intent(in)    :: diag_flag 
  logical                                     :: send_diagnostics 

  real, dimension(isd:ied,jsd:jed,2)  :: tracr
  integer                             :: i, j, k
  integer                             :: tau, taum1
 
  if(.not. use_this_module) return 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_lap_tracer_mod (lap_tracer): must be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! fx = flux component through "eastern"  face of T-cells at level k
  ! fy = flux component through "northern" face of T-cells at level k

  if(horz_z_diffuse) then 

      ! tracr(:,:,1) = Tracer%field at level k-1
      ! tracr(:,:,2) = Tracer%field at level k

      tracr(:,:,1) = Tracer%field(:,:,1,taum1)
      do k=1,nk
         tracr(:,:,2) = Tracer%field(:,:,k,taum1)
         fx(:,:,k)    = diff_cet(:,:,k)*FDX_ZT(tracr(:,:,1:2),k)*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         fy(:,:,k)    = diff_cnt(:,:,k)*FDY_ZT(tracr(:,:,1:2),k)*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         tracr(:,:,1) = tracr(:,:,2)
      enddo

  elseif(horz_s_diffuse) then 

      do k=1,nk
         tracr(:,:,1) = Tracer%field(:,:,k,taum1)
         fx(:,:,k)    = diff_cet(:,:,k)*FDX_T(tracr(:,:,1))*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         fy(:,:,k)    = diff_cnt(:,:,k)*FDY_T(tracr(:,:,1))*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
      enddo

  endif

  if(Grd%tripolar) call mpp_update_domains(fx, fy, Dom_flux%domain2d, gridtype=CGRID_NE)

  do k=1,nk
     wrk1(:,:,k) = (BDX_ET(fx(:,:,k)) + BDY_NT(fy(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Grd%tmask(i,j,k)*wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! send fluxes to diag_manager 
  ! minus sign is due to MOM-convention for physics fluxes. 
  if(send_diagnostics) then 

      if (id_xflux_diff(ntracer) > 0) then
          do k=1,nk
            fx(:,:,k) = fx(:,:,k)*Grd%dyte(:,:)
          enddo
          call diagnose_3d(Time, Grd, id_xflux_diff(ntracer), -1.0*fx(:,:,:)*Tracer%conversion)
      endif
      if (id_yflux_diff(ntracer) > 0) then
          do k=1,nk
            fy(:,:,k) = fy(:,:,k)*Grd%dxtn(:,:)
          enddo
          call diagnose_3d(Time, Grd, id_yflux_diff(ntracer), -1.0*fy(:,:,:)*Tracer%conversion)
      endif
      if (id_h_diffuse(ntracer) > 0) then
         call diagnose_3d(Time, Grd, id_h_diffuse(ntracer), wrk1(:,:,:)*Tracer%conversion)
      endif

  endif

  if (have_obc) then 
     call store_ocean_obc_tracer_flux(Time, Tracer, -1.0*fx, ntracer, 'z', 'dif')
     call store_ocean_obc_tracer_flux(Time, Tracer, -1.0*fy, ntracer, 'm', 'dif')
  endif
  
end subroutine lap_tracer
! </SUBROUTINE> NAME="lap_tracer"

end module ocean_lap_tracer_mod
      
      

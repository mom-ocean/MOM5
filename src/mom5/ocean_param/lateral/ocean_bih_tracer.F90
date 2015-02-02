module ocean_bih_tracer_mod
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted and density weighted time tendency
! for tracer from biharmonic tracer mixing. 
!</OVERVIEW>
!
!<DESCRIPTION>
! There are two main options for computing the fluxes.
!
! (1) The lateral fluxes can be aligned with the z-coordinate surfaces,
! in which case the fluxes must be approximated if (i) we use non-geopotential
! vertical coordinates, (ii) next to partial bottom step topography.
! This form of the diffusion is not recommended since it can lead to 
! the creation of spurious extrema.  
!
! (2) The lateral fluxes can be aligned surfaces of constant vertical 
! coordinate. In this case the fluxes are no longer strictly "horizontal."
! Howerver, the operator is simpler and it ensures that no suprious
! extrema are created. It is for this reason that the simpler operator
! is preferred.   
!  
! The diffusivity used to determine the strength of the tendency can be 
! a general function of space yet it is constant in time.  A namelist 
! option exists that determines this diffusivity as a local function 
! of the grid spacing. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_bih_tracer_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module
!  </DATA> 
!  <DATA NAME="horz_z_diffuse" TYPE="logical">
!  To compute fluxes along surfaces of constant depth.
!  This operation must necessarily be approximate for the two 
!  cases (i) non-geopotential vertical coordinates, (2)
!  next to partial bottom step topography.  There are cases where
!  use of this operator can lead to spurious creation of extrema
!  due to truncation errors associated with the "slope" term.
!  The option to use horz_z_diffuse=.true. is maintained for 
!  legacy purposes alone.  
!  </DATA> 
!  <DATA NAME="horz_s_diffuse" TYPE="logical">
!  To compute diffusion along surfaces of constant vertical s-coordinate. 
!  </DATA> 
!  <DATA NAME="abih" UNITS="m^4/sec" TYPE="real">
!  This is the value for the space-time constant biharmonic diffusivity. 
!  </DATA> 
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the diffusivity is set according to a velocity scale times
!  the cube of the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model.  
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
!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, register_static_field
use fms_mod,          only: FATAL, NOTE, stdout, stdlog, read_data
use fms_mod,          only: write_version_number, open_namelist_file, close_file, check_nml_error
use mpp_domains_mod,  only: mpp_update_domains, CGRID_NE
use mpp_mod,          only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_operators_mod,  only: FMX, FMY, BDX_ET, BDY_NT
use ocean_operators_mod,  only: FDX_ZT, FDY_ZT, FDX_T, FDY_T
use ocean_obc_mod,        only: store_ocean_obc_tracer_flux, ocean_obc_enhance_diff_back
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, ocean_time_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_prog_tracer_type, ocean_options_type
use ocean_workspace_mod,  only: wrk1 
use ocean_util_mod,       only: diagnose_3d


implicit none

private

public ocean_bih_tracer_init
public bih_tracer
private del2_tracer

real ::    abih      = 0.0e10        ! constant horz biharmonic mixing coefficient for tracers
real ::    vel_micom = 0.0           ! constant velocity scale (m/s) for setting micom diffusivity  
logical :: tracer_mix_micom =.false. ! for spatial dependent diffusivity 

real, dimension(:,:,:), allocatable :: diff_cet         ! diffusivity for eastern face of T cell (m^4/s)
real, dimension(:,:,:), allocatable :: diff_cnt         ! diffusivity for northern face of T cell (m^4/s)
real, dimension(:,:,:), allocatable :: del2_tracer      ! Laplacian operator acting on tracers
real, dimension(:,:,:), allocatable :: diffusivity_mask ! 3d mask to selectively apply diffusion 

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xflux_diff
integer, dimension(:), allocatable  :: id_yflux_diff
integer, dimension(:), allocatable  :: id_h_diffuse
integer                             :: id_ah_biharmonic

character(len=128) :: version=&
     '=>Using: /bih/ocean_bih_tracer.F90 ($Id: ocean_bih_tracer.F90,v 20.0 2013/12/14 00:14:10 fms Exp $)'
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
logical :: module_is_initialized = .FALSE.
logical :: read_diffusivity_mask = .false.
logical :: have_obc=.false.

namelist /ocean_bih_tracer_nml/ use_this_module, abih, tracer_mix_micom, vel_micom, &
                                read_diffusivity_mask, horz_z_diffuse, horz_s_diffuse

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_tracer_init">
!
! <DESCRIPTION>
! Initialize the biharmonic tracer mixing module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that diffusivity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bih_tracer_init(Grid, Domain, Time, T_prog, Ocean_options, dtime, obc, tmask_sigma)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  real,                         intent(in)           :: dtime
  logical,                      intent(in)           :: obc
  real,                         intent(in), optional :: tmask_sigma(Domain%isc:Domain%iec,Domain%jsc:Domain%jec)

  real    :: dxdymn, ah_crit
  integer :: ioun, io_status, ierr
  integer :: i, j, k, n, num
 
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_tracer_mod (ocean_bih_tracer_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  have_obc = obc

  num_prog_tracers = size(T_prog(:))

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_bih_tracer_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_bih_tracer_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_bih_tracer_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_bih_tracer_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bih_tracer_nml)  
  write (stdlogunit,ocean_bih_tracer_nml)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(use_this_module) then 
      write(stdoutunit,*)'==>Note from ocean_bih_tracer_mod: Using this module.'
      Ocean_options%horz_bih_tracer = 'Used horizontal biharmonic tracer mixing option.'
 else
      write(stdoutunit,*)'==>Note from ocean_bih_tracer_mod: NOT using this module.'
      Ocean_options%horz_bih_tracer = 'Did NOT use horizontal biharmonic tracer mixing.'
      return 
  endif

  if(abih > 0 .or. (tracer_mix_micom .and. vel_micom > 0)) then 
      write (stdoutunit,'(1x,a)') ' ==>WARNING: biharmonic tracer mixing contributes to unphysically'
      write (stdoutunit,'(7x,a)') 'large diapycnal mixing. Such may be damaging for longer simulations.'
      write (stdoutunit,'(7x,a)') 'It is recommended that one instead uses neutral physics.'
  endif

  if(horz_z_diffuse) then 
      write(stdoutunit,*)'==>Warning: Computing biharmonic tracer mixing along constant geoptoential surfaces.'
      write(stdoutunit,*)'   This operator can create spurious extrema where s-surfaces deviate from z-surfaces.'
      write(stdoutunit,*)'   For this reason, the alternative horz_s_diffuse=.true. is recommended.'
  endif

  if(horz_s_diffuse) then 
      write(stdoutunit,*)'==>Note: Computing biharmonic tracer mixing along constant s-surfaces.'
  endif

  if(horz_s_diffuse .and. horz_z_diffuse) then 
      write(stdoutunit,*)'==>Error from ocean_bih_tracer_mod: use only one of horz_s_diffuse and horz_z_diffuse.'
      write(stdoutunit,*)'   Presently, both are set true and this is not allowed. '
      call mpp_error(FATAL,'==>Error ocean_bih_tracer_mod: use only one of horz_s_diffuse and horz_z_diffuse.')
  endif

  if(.not. horz_s_diffuse .and. .not. horz_z_diffuse) then 
      write(stdoutunit,*)'==>Warning from ocean_bih_tracer_mod: use_this_module=.true.'
      write(stdoutunit,*)'   However, both horz_s_diffuse and horz_z_diffuse are .false.'
      write(stdoutunit,*)'   So will not produce any lateral biharmonic tracer mixing.'
  endif


  write(stdoutunit,'(/a,f10.2)')'==> Note from ocean_bih_tracer_mod: using forward time step (secs)', dtime 

  Dom => Domain
  Grd => Grid
  
  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux', maskmap=Dom%maskmap)

  allocate (diff_cet(isd:ied,jsd:jed,nk))
  allocate (diff_cnt(isd:ied,jsd:jed,nk))
  allocate (del2_tracer(isd:ied,jsd:jed,nk))
  allocate (fx(isd:ied,jsd:jed,nk))
  allocate (fy(isd:ied,jsd:jed,nk))  
  
  diff_cet    = abih
  diff_cnt    = abih
  del2_tracer = 0.0
  fx          = 0.0
  fy          = 0.0
  
  ! Micom background diffusivity
  ! Space scale set by grid size
  ! Velocity scale input via namelist 
  if(tracer_mix_micom) then
    do k=1,nk
      diff_cet(:,:,k) = vel_micom*(2.0*Grd%dxt(:,:)*Grd%dyte(:,:)/(Grd%dxt(:,:)+Grd%dyte(:,:)))**3
      diff_cnt(:,:,k) = vel_micom*(2.0*Grd%dxtn(:,:)*Grd%dyt(:,:)/(Grd%dxtn(:,:)+Grd%dyt(:,:)))**3
    enddo
    do j=jsc,jec
      write (stdoutunit,'(a,i4,a,e14.7,a)') ' MICOM biharmonic diffusivity at (isc,',j,',1) = ',diff_cet(isc,j,1),' m^4/s'
    enddo
  endif 

  ! enhance background diffusion near open boundaries if needed
  if (have_obc) call ocean_obc_enhance_diff_back(diff_cet, diff_cnt, 'bih')
  
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
      write (stdoutunit,*) '==>ocean_bih_tracer_mod: Completed read of 3d tracer diffusivity mask.'
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

  if (dtime /= 0.0) then
    num = 0
    do j=jsc,jec
      do i=isc,iec
        dxdymn = (2.0/(1.0/(Grd%dxt(i,j))**2 + 1.0/Grd%dyt(i,j)**2))**2
        ah_crit = 0.0625*dxdymn/dtime
        if (diff_cet(i,j,1)*Grd%tmask(i,j,1) > ah_crit .and. num <= 10 .and. i == isc) then
          num = num + 1
          if (num == 1) write (stdoutunit,'(/,(1x,a))')&
          '==> Warning: horz diffusive criteria exceeded for "ah". use a smaller "dtts", "ah". Show first 10 violations' 
          write (stdoutunit,'(a,i4,a,i4,a,f6.2,a,f6.2,a,2(e14.7,a))') ' at (i,j)= (',i,',',j,'), (lon,lat)= (', &
                Grd%xt(i,j),',',Grd%yt(i,j), '),  "ah" = ', diff_cet(i,j,1),' m^4/s. the critical value =',ah_crit,' m^4/s'
        endif
      enddo
    enddo
  endif

  ! register for diagnostics manager 
  allocate (id_xflux_diff(num_prog_tracers))
  allocate (id_yflux_diff(num_prog_tracers))
  allocate (id_h_diffuse(num_prog_tracers))
  id_xflux_diff=-1
  id_yflux_diff=-1
  id_h_diffuse =-1

  do n=1, num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_xflux_diff(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'i-diffusive heat flux', 'Watts',         &
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_yflux_diff(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'j-diffusive heat flux', 'Watts',         &
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model',          &
              trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'horz-diffusion heating', 'Watts/m^2',   &
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
              Time%model_time, 'horz-diffusion of '//trim(T_prog(n)%name), 'kg/(m^2*sec)',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
     endif

  enddo

  id_ah_biharmonic = register_static_field ('ocean_model', 'ah_biharmonic',&
                     Grd%tracer_axes(1:3),'static biharmonic diffusivity', &
                     'm^4/sec', missing_value=missing_value, range=(/-10.0,1.e20/))
  call diagnose_3d(Time, Grd, id_ah_biharmonic, diff_cet(:,:,:))


end subroutine ocean_bih_tracer_init
! </SUBROUTINE>  NAME="ocean_bih_tracer_init"


!#######################################################################
! <SUBROUTINE NAME="bih_tracer">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from biharmonic mixing. 
! </DESCRIPTION>
!
subroutine bih_tracer (Time, Thickness, Tracer, ntracer, diag_flag)

  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_prog_tracer_type), intent(inout)        :: Tracer
  integer,                      intent(in)           :: ntracer
  logical,                      intent(in), optional :: diag_flag 
  logical :: send_diagnostics 

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

  ! fx = flux component through "eastern"  face of T-cells at level k
  ! fy = flux component through "northern" face of T-cells at level k
  ! tracr(:,:,1) = tracer at level k-1
  ! tracr(:,:,2) = tracer at level k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_tracer_mod (bih_tracer): needs initialization')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  call delsq_tracer (Thickness, Tracer%field(:,:,:,taum1), tau)
  if(Dom%xhalo==1 .and. Dom%yhalo==1) call mpp_update_domains (del2_tracer, Dom%domain2d)

  if(horz_z_diffuse) then 
      tracr(:,:,1) = -del2_tracer(:,:,1)
      do k=1,nk
         tracr(:,:,2) = -del2_tracer(:,:,k)
         fx(:,:,k)    = diff_cet(:,:,k)*FDX_ZT(tracr(:,:,1:2),k)*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         fy(:,:,k)    = diff_cnt(:,:,k)*FDY_ZT(tracr(:,:,1:2),k)*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         tracr(:,:,1) = tracr(:,:,2)
      enddo
  elseif(horz_s_diffuse) then 
      do k=1,nk
         fx(:,:,k) = -diff_cet(:,:,k)*FDX_T(del2_tracer(:,:,k))*FMX(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
         fy(:,:,k) = -diff_cnt(:,:,k)*FDY_T(del2_tracer(:,:,k))*FMY(Thickness%rho_dzt(:,:,k,tau)*Grd%tmask(:,:,k))
      enddo
  endif

  !  set redundancies for tripolar
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
  ! absence of minus sign relative to laplacian due to minus sign used for biharmonic.
  if(send_diagnostics) then 

      if (id_xflux_diff(ntracer) > 0) then
          do k=1,nk
            fx(:,:,k) = fx(:,:,k)*Grd%dyte(:,:)
          enddo
          call diagnose_3d(Time, Grd, id_xflux_diff(ntracer), fx(:,:,:)*Tracer%conversion)
      endif
      if (id_yflux_diff(ntracer) > 0) then
          do k=1,nk
            fy(:,:,k) = fy(:,:,k)*Grd%dxtn(:,:)
          enddo
          call diagnose_3d(Time, Grd, id_yflux_diff(ntracer), fy(:,:,:)*Tracer%conversion)
      endif
      if (id_h_diffuse(ntracer) > 0) then
         call diagnose_3d(Time, Grd, id_h_diffuse(ntracer), wrk1(:,:,:)*Tracer%conversion)
      endif

  endif

  if (have_obc) then 
     call store_ocean_obc_tracer_flux(Time, Tracer, -1.0*fx, ntracer, 'z', 'dif')
     call store_ocean_obc_tracer_flux(Time, Tracer, -1.0*fy, ntracer, 'm', 'dif')
  endif
  
end subroutine bih_tracer
! </SUBROUTINE> NAME="bih_tracer"


!#######################################################################
! <SUBROUTINE NAME="delsq_tracer">
!
! <DESCRIPTION>
! Subroutine computes the laplacian operator acting on tracer with unit 
! diffusivity. Units of del2_tracer are tracer/length^2
! </DESCRIPTION>
!
subroutine delsq_tracer (Thickness, tracer, index)

  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:,:), intent(in) :: tracer
  integer,                      intent(in) :: index

  real, dimension(isd:ied,jsd:jed)   :: fe, fn
  real, dimension(isd:ied,jsd:jed,2) :: tracr
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_bih_tracer_mod (delsq_tracer): module needs to be initialized')
  endif 


  ! fe = Diffusive flux across east face of T cells. Use unit diffusivity
  ! fn = Diffusive flux across north face of T cells. Use unit diffusivity

  if(horz_z_diffuse) then 
      tracr(:,:,1) = tracer(:,:,1)
      do k=1,nk
         tracr(:,:,2) = tracer(:,:,k)
         fe(:,:) = FDX_ZT(tracr(:,:,1:2),k)*FMX(Thickness%rho_dzt(:,:,k,index)*Grd%tmask(:,:,k))
         fn(:,:) = FDY_ZT(tracr(:,:,1:2),k)*FMY(Thickness%rho_dzt(:,:,k,index)*Grd%tmask(:,:,k))
         del2_tracer(:,:,k) = Grd%tmask(:,:,k)*(BDX_ET(fe(:,:)) + BDY_NT(fn(:,:)))/Thickness%rho_dzt(:,:,k,index)
         tracr(:,:,1) = tracr(:,:,2)
      enddo
  elseif(horz_s_diffuse) then
      do k=1,nk
         fe(:,:) = FDX_T(tracer(:,:,k))*FMX(Thickness%rho_dzt(:,:,k,index)*Grd%tmask(:,:,k))
         fn(:,:) = FDY_T(tracer(:,:,k))*FMY(Thickness%rho_dzt(:,:,k,index)*Grd%tmask(:,:,k))
         del2_tracer(:,:,k) = Grd%tmask(:,:,k)*(BDX_ET(fe(:,:)) + BDY_NT(fn(:,:)))/Thickness%rho_dzt(:,:,k,index)
      enddo
  endif


end subroutine delsq_tracer
! </SUBROUTINE> NAME="delsq_tracer"


end module ocean_bih_tracer_mod
      
      

module ocean_mixdownslope_mod
#define COMP isc:iec,jsc:jec
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Mixing of tracer between dense shallow parcel and 
! deeper parcels downslope. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Mixing of tracer properties as dense shallow parcel is discharged
! into deeper water to approach the parcel's depth of neutral buoyancy.
! This module can be characterized as a mixture of the approach from
! Campin and Goosse (1999) and slope convection. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Campin and Goosse (1999): Parameterization of density-driven downsloping flow 
! for a coarse-resolution model in z-coordinate", Tellus 51A, pages 412-430
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_mixdownslope_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  For using this module.  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!
!  <DATA NAME="mixdownslope_npts" TYPE="integer">
!  Number of horizontally distant points used in search downslope. 
!  Note: it is not possible to have 
!  mixdownslope_npts greater than or equal to the computational domain
!  extents, as this would require updates across multiple processors. 
!  Default mixdownslope_npts=1. 
!  </DATA> 
!  <DATA NAME="mixdownslope_frac_central" TYPE="real">
!  Fraction of the central cell that participates in downslope mixing
!  in any particular direction.  Default mixdownslope_frac_central=0.25
!  </DATA> 
!
!  <DATA NAME="mixdownslope_weight_far" TYPE="logical">
!  To place more weight on points further from central point.  This may 
!  be done to enhance properties getting downslope.  Default is 
!  mixdownslope_weight_far=.false.
!  </DATA> 
!  <DATA NAME="mixdownslope_width" TYPE="integer">
!  Width of the re-weighting function used to emphasize points further 
!  along in the search for exchange points.  Default mixdownslope_width=1.
!  </DATA> 
!
!  <DATA NAME="read_mixdownslope_mask" TYPE="logical">
!  For reading in a mask that selects regions of the domain 
!  where mixdownslope is allowed to function (mask=1) or not 
!  to function (mask=0). Default read_mixdownslope_mask=.false., 
!  whereby mixdownslope_mask is set to tmask(k=1). 
!  </DATA> 
!  <DATA NAME="mixdownslope_mask_gfdl" TYPE="logical">
!  For modifying the mixdownslope mask based on reading in 
!  the GFDL regional mask. Default mixdownslope_mask_gfdl=.false.
!  </DATA> 
!
!</NAMELIST>
!

use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, register_static_field
use fms_mod,             only: write_version_number, error_mesg, read_data, FATAL, NOTE
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_domains_mod,     only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,             only: input_nml_file, mpp_error

use ocean_density_mod,     only: density
use ocean_domains_mod,     only: get_local_indices, set_ocean_domain
use ocean_parameters_mod,  only: missing_value, rho0, rho0r
use ocean_parameters_mod,  only: TERRAIN_FOLLOWING
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_time_type, ocean_options_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_density_type, ocean_thickness_type
use ocean_util_mod,        only: write_timestamp, diagnose_2d, diagnose_3d, diagnose_sum
use ocean_util_mod,        only: write_chksum_2d, write_chksum_3d, write_chksum_2d_int
use ocean_tracer_util_mod, only: diagnose_3d_rho
use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_v

implicit none

private 

public ocean_mixdownslope_init
public mixdownslope

private watermass_diag_init
private watermass_diag

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

#include <ocean_memory.h>

integer, dimension(:,:), allocatable :: ip    ! shifting in i-direction
integer, dimension(:,:), allocatable :: jq    ! shifting in j-direction

real, dimension(:,:),     allocatable :: data              ! for reading in mask information 
real, dimension(:,:),     allocatable :: area_k            ! area on k-level, including tmask 
real, dimension(:,:),     allocatable :: slope_x           ! topog slope (dimensionless) in the i-direction
real, dimension(:,:),     allocatable :: slope_y           ! topog slope (dimensionless) in the j-direction 
real, dimension(:,:),     allocatable :: mixdownslope_mask ! mask=1 when apply mixdownslope, 0 otherwise
real, dimension(:,:,:),   allocatable :: topog_step        ! where downslope flow is topgraphically possible 
real, dimension(:,:,:),   allocatable :: tend_mix          ! (rho*dzt*tracer) tendency
real, dimension(:,:,:,:), allocatable :: mixdownslope_frac ! fraction of mixing occurring between two cells 

! fields with extended halos 
integer, dimension(:,:),     allocatable :: kmt_ex
integer, dimension(:,:,:),   allocatable :: kup_ex
integer, dimension(:,:,:,:), allocatable :: kdw_ex
real, dimension(:,:,:),   allocatable :: topog_slope_ex    
real, dimension(:,:,:),   allocatable :: topog_step_ex
real, dimension(:,:,:),   allocatable :: tracer_ex
real, dimension(:,:,:),   allocatable :: temp_ex
real, dimension(:,:,:),   allocatable :: salt_ex
real, dimension(:,:,:),   allocatable :: press_ex
real, dimension(:,:,:),   allocatable :: mass_ex
real, dimension(:,:,:),   allocatable :: rho_ex

! domain for mixdownslope 
type(ocean_domain_type), save :: Mixdownslope_domain 

! work array for eta_tend 
real, dimension(:,:,:),   allocatable :: theta_tend
real, dimension(:,:,:),   allocatable :: salt_tend


! internally set for computing watermass diagnstics
logical :: compute_watermass_diag = .false. 

! for global area normalization
real    :: cellarea_r

! for diagnostic manager 
logical :: used
integer :: id_slope_x=-1
integer :: id_slope_y=-1
integer :: id_topog_step_1       =-1
integer :: id_topog_step_2       =-1
integer :: id_topog_step_3       =-1
integer :: id_topog_step_4       =-1
integer :: id_mixdownslope_mask  =-1

integer :: id_neut_rho_mixdown          =-1
integer :: id_wdian_rho_mixdown         =-1
integer :: id_tform_rho_mixdown         =-1
integer :: id_neut_rho_mixdown_on_nrho  =-1
integer :: id_wdian_rho_mixdown_on_nrho =-1
integer :: id_tform_rho_mixdown_on_nrho =-1

integer :: id_eta_tend_mixdown     =-1
integer :: id_eta_tend_mixdown_glob=-1

integer :: id_neut_temp_mixdown          =-1
integer :: id_wdian_temp_mixdown         =-1
integer :: id_tform_temp_mixdown         =-1
integer :: id_neut_temp_mixdown_on_nrho  =-1
integer :: id_wdian_temp_mixdown_on_nrho =-1
integer :: id_tform_temp_mixdown_on_nrho =-1

integer :: id_neut_salt_mixdown          =-1
integer :: id_wdian_salt_mixdown         =-1
integer :: id_tform_salt_mixdown         =-1
integer :: id_neut_salt_mixdown_on_nrho  =-1
integer :: id_wdian_salt_mixdown_on_nrho =-1
integer :: id_tform_salt_mixdown_on_nrho =-1


integer, dimension(:), allocatable :: id_mixdownslope 


character(len=128) :: version=&
       '=>Using: ocean_mixdownslope.f90 ($Id: ocean_mixdownslope.F90,v 20.0 2013/12/14 00:14:28 fms Exp $)'
character (len=128) :: tagname=&
     '$Name: tikal $'

! number of prognostic tracers
integer :: num_prog_tracers=0

! processor zero writes to unit 6
integer :: unit=6  

! initialization flag 
logical :: module_is_initialized=.false.

! flag for mpp_global_sum
integer :: global_sum_flag      

! halo 
integer :: ijhalo=1 

! time step 
real :: dtime 
real :: dtimer

! set from nml
logical :: use_this_module        = .false. ! must be set .true. in nml to enable this scheme
logical :: debug_this_module      = .false. ! for debugging
logical :: do_bitwise_exact_sum   = .false. ! set true to get slower sum that is same for all PEs. 
logical :: mixdownslope_weight_far= .false. ! to place more weight on points further from central point
logical :: read_mixdownslope_mask = .false. ! for applying a mask where do not apply mixdownslope
logical :: mixdownslope_mask_gfdl = .false. ! for case when applying GFDL regional mask
integer :: mixdownslope_npts      = 1       ! number horizontal points searched for the downslope mixing
integer :: mixdownslope_width     = 1       ! width for exponential that damps points near to central point
real    :: mixdownslope_frac_central =0.25  ! fraction of central cell participating in downslope mixing


namelist /ocean_mixdownslope_nml/ use_this_module, debug_this_module, mixdownslope_npts, &
                                  mixdownslope_width, mixdownslope_weight_far,           &
                                  mixdownslope_frac_central, do_bitwise_exact_sum,       &
                                  read_mixdownslope_mask, mixdownslope_mask_gfdl
                                 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_mixdownslope_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers into the abyss next to topography.
! </DESCRIPTION>
!
subroutine ocean_mixdownslope_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, &
                                     vert_coordinate_type, dtim, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in), target   :: Time
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  integer,                      intent(in)           :: vert_coordinate_type 
  real,                         intent(in)           :: dtim
  logical,                      intent(in), optional :: debug

  integer :: io_status, ioun, ierr
  integer :: i,j,m,n,kbot

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
           '==>Error from ocean_mixdownslope_mod (ocean_mixdownslope_init): module already initialized')
  endif

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)
#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif
  Dom => Domain
  Grd => Grid

  dtime      = dtim 
  dtimer     = 1.0/dtime 
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

  num_prog_tracers = size(T_prog(:))

! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_mixdownslope_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_mixdownslope_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_mixdownslope_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_mixdownslope_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_mixdownslope_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_mixdownslope_nml)

  if(do_bitwise_exact_sum) then
      global_sum_flag = BITWISE_EXACT_SUM
  else
      global_sum_flag = NON_BITWISE_EXACT_SUM
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
  endif
  if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_mixdownslope_mod with debug_this_module=.true.'  
  endif

  if(.not. use_this_module) then 
      call mpp_error(NOTE,&
           '==>From ocean_mixdownslope_mod: NOT using downslope mixing scheme.')
      Ocean_options%mixdownslope = 'Did NOT use downslope mixing scheme.'
      return
  else 
      if(vert_coordinate_type == TERRAIN_FOLLOWING) then 
          call mpp_error(FATAL, &
               '==>ocean_mixdownslope_mod: this module is NOT for use with TERRAIN_FOLLOWING vert coodinates.')
      endif
      Ocean_options%mixdownslope = 'Used the downslope mixing scheme. This scheme is experimental!'
      call mpp_error(NOTE,&
           '==>From ocean_mixdownslope_mod: USING downslope mixing scheme.')
  endif
  if(debug_this_module) then 
      call mpp_error(NOTE,'==>From ocean_mixdownslope_mod: USING debug_this_module')
  endif

  ! allocate for tracer tendency array 
  allocate(tend_mix(isd:ied,jsd:jed,nk))
  tend_mix(:,:,:) = 0.0

  if(mixdownslope_npts<1) then
      call mpp_error(FATAL, &
           '==>ocean_mixdownslope_mod: mixdownslope_npts < 1 not allowed. In nml, set mixdownslope_npts>=1.')
  endif
  write(stdoutunit,'(a,i3)')' In ocean_mixdownslope_mod: mixdownslope_npts = ',mixdownslope_npts  
  write(stdoutunit,'(a)')   ' Be sure this number is smaller than dimensions of computational domain.'


  ! set mixdownslope_mask=0.0 in those regions where mixdownslope is NOT applied
  ! and mixdownslope_mask=1.0 in regions where mixdownslope is applied.  Default is 
  ! mixdownslope everywhere (mixdownslope_mask(:,:) = Grd%tmask(:,:,1)).  
  allocate(mixdownslope_mask(isd:ied,jsd:jed))
  mixdownslope_mask(:,:) = Grd%tmask(:,:,1)
  if(read_mixdownslope_mask) then 
      allocate(data(isd:ied,jsd:jed))
      call read_data('INPUT/mixdownslope_mask','mixdownslope_mask',data,Domain%domain2d)
      do j=jsc,jec
         do i=isc,iec
            mixdownslope_mask(i,j) = data(i,j)
         enddo
      enddo

      ! the following is specific to the mask used at GFDL for
      ! the OM3 ocean model configuration. 
      ! here, remove the Black Sea (labelled with 7.0), as this
      ! region contains some odd water masses which cause instabilities
      ! with the mixdownslope scheme. 
      if(mixdownslope_mask_gfdl) then 
          do j=jsc,jec
             do i=isc,iec
                if(mixdownslope_mask(i,j)==0.0 .or. mixdownslope_mask(i,j)==7.0) then 
                    mixdownslope_mask(i,j)=0.0
                else 
                    mixdownslope_mask(i,j)=1.0
                endif
             enddo
          enddo
      endif
      call mpp_update_domains(mixdownslope_mask(:,:), Dom%domain2d)
  endif


  ! halo size for extended domain and larger arrays 
  ijhalo=mixdownslope_npts

  ! k-level of central point in search 
  allocate(kup_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
  kup_ex(:,:,:) = 0

  ! k-level of deeper points in horizontally adjacent columns 
  allocate(kdw_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,ijhalo,4))
  kdw_ex(:,:,:,:) = 0

  ! fraction of a cell participating in mixing 
  allocate(mixdownslope_frac(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,ijhalo,4))
  mixdownslope_frac(:,:,:,:) = 0.0

  ! define extended domain 
  call set_ocean_domain(Mixdownslope_domain,Grd,xhalo=ijhalo,yhalo=ijhalo,&
                        name='mixdownslope',maskmap=Dom%maskmap)

  ! time independent arrays defined over Mixdownslope_domain 
  allocate(kmt_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo))
  kmt_ex(:,:) = 0
  kmt_ex (isc:iec,jsc:jec) = Grd%kmt (isc:iec,jsc:jec)
  call mpp_update_domains (kmt_ex(:,:), Mixdownslope_domain%domain2d)

  ! time dependent arrays defined over Mixdownslope_domain
  allocate(temp_ex   (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  allocate(salt_ex   (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  allocate(tracer_ex (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  allocate(press_ex  (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  allocate(rho_ex    (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  allocate(mass_ex   (isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,nk))
  temp_ex(:,:,:)   = 0.0
  salt_ex(:,:,:)   = 0.0
  tracer_ex(:,:,:) = 0.0
  press_ex(:,:,:)  = 0.0
  rho_ex(:,:,:)    = 0.0
  mass_ex(:,:,:)   = 0.0

  ! area of cells on k-level, including tmask 
  allocate (area_k(isd:ied,jsd:jed))
  area_k(:,:) = 0.0

  ! compute topographic slope arrays for the i-slope and j-slope
  ! slopes are centered on the i-face and j-face of tracer cells 
  allocate (slope_x(isd:ied,jsd:jed))
  allocate (slope_y(isd:ied,jsd:jed))
  slope_x = 0.0 
  slope_y = 0.0
  do j=jsc,jec
     do i=isc,iec
        slope_x(i,j) = (Grd%ht(i+1,j)-Grd%ht(i,j))*Grd%dxter(i,j)*Grd%tmask(i,j,1)*Grd%tmask(i+1,j,1)
        slope_y(i,j) = (Grd%ht(i,j+1)-Grd%ht(i,j))*Grd%dytnr(i,j)*Grd%tmask(i,j,1)*Grd%tmask(i,j+1,1)
        slope_x(i,j) = abs(slope_x(i,j))
        slope_y(i,j) = abs(slope_y(i,j))
     enddo
  enddo
  call mpp_update_domains(slope_x(:,:),Dom%domain2d)
  call mpp_update_domains(slope_y(:,:),Dom%domain2d)

  ! topographic slope for the four surrounding directions 
  allocate (topog_slope_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
  topog_slope_ex(:,:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        m=1 ; topog_slope_ex(i,j,m) = slope_x(i,j)  
        m=2 ; topog_slope_ex(i,j,m) = slope_y(i,j)  
        m=3 ; topog_slope_ex(i,j,m) = slope_x(i-1,j)
        m=4 ; topog_slope_ex(i,j,m) = slope_y(i,j-1)
     enddo
  enddo
  call mpp_update_domains (topog_slope_ex(:,:,:), Mixdownslope_domain%domain2d)

  ! compute directions from an (i,j) point where topography deepens.
  ! these directions may potentially have downslope mixing. 
  ! insist that downslope mixing occurs only when there are more 
  ! kmt cells in the adjacent column. also insist that downslope 
  ! mixing does not involve k=1 cells. 
  allocate (topog_step(isd:ied,jsd:jed,4))
  topog_step(:,:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        kbot = Grd%kmt(i,j)
        if(kbot > 1) then 
            if(Grd%kmt(i+1,j) > kbot) topog_step(i,j,1)=1.0
            if(Grd%kmt(i,j+1) > kbot) topog_step(i,j,2)=1.0
            if(Grd%kmt(i-1,j) > kbot) topog_step(i,j,3)=1.0
            if(Grd%kmt(i,j-1) > kbot) topog_step(i,j,4)=1.0
        endif
     enddo
  enddo
  ! block out the bipolar fold in order to ensure tracer conservation.
  ! The reason we do so is related to how the algorithm reaches between  
  ! adjacent columns of tracer points.  When the column straddles 
  ! the bipolar fold, the code logic is not general and so it actually 
  ! attempts to reach to a non-adjacent column.  Shy of generalizing 
  ! the logic, we simply do not consider mixdownslope for points along fold. 
  if(jec+Dom%joff==Dom%jeg) then 
      topog_step(:,jec,:) = 0.0
  endif

  ! fill the larger array 
  allocate (topog_step_ex(isc-ijhalo:iec+ijhalo,jsc-ijhalo:jec+ijhalo,4))
  topog_step_ex(:,:,:) = 0.0
  topog_step_ex(isc:iec,jsc:jec,:) = topog_step(isc:iec,jsc:jec,:)
  call mpp_update_domains (topog_step_ex(:,:,:), Mixdownslope_domain%domain2d)

  ! labels for the four quadranats surrounding a 
  ! central tracer cell (moving counter-clockwise)
  allocate(ip(0:ijhalo,4))
  allocate(jq(0:ijhalo,4))

  do n=0,ijhalo
     m=1 
     ip(n,m)   = 1+(n-1)
     jq(n,m)   = 0  

     m=2 
     ip(n,m)   = 0  
     jq(n,m)   = 1+(n-1)  

     m=3 
     ip(n,m)   = -1-(n-1) 
     jq(n,m)   = 0  

     m=4 
     ip(n,m)   =  0  
     jq(n,m)   = -1-(n-1)
  enddo

  ! register/send diagnostics 
  id_mixdownslope_mask = register_static_field ('ocean_model', 'mixdownslope_mask', &
       Grd%tracer_axes(1:2), &
       'mixdownslope mask', 'dimensionless', missing_value=missing_value, range=(/-1.0,1e4/))
  call diagnose_2d(Time, Grd, id_mixdownslope_mask, mixdownslope_mask(:,:))

  id_slope_x = register_static_field ('ocean_model', 'slope_x', Grd%tracer_axes_flux_x(1:2), &
       '|d(ht)/dx| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
  call diagnose_2d(Time, Grd, id_slope_x, slope_x(:,:))

  id_slope_y = register_static_field ('ocean_model', 'slope_y', Grd%tracer_axes_flux_y(1:2), &
       '|d(ht)/dy| on T-cell face', 'm/m', missing_value=missing_value, range=(/-1.e9,1.e9/))
  call diagnose_2d(Time, Grd, id_slope_y, slope_y(:,:))

  id_topog_step_1 = register_static_field ('ocean_model', 'topog_step_1', Grd%tracer_axes(1:2), &
       'topog_step_1', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
  call diagnose_2d(Time, Grd, id_topog_step_1, topog_step(:,:,1))

  id_topog_step_2 = register_static_field ('ocean_model', 'topog_step_2', Grd%tracer_axes(1:2), &
       'topog_step_2', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
  call diagnose_2d(Time, Grd, id_topog_step_2, topog_step(:,:,2))

  id_topog_step_3 = register_static_field ('ocean_model', 'topog_step_3', Grd%tracer_axes(1:2), &
       'topog_step_3', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
  call diagnose_2d(Time, Grd, id_topog_step_3, topog_step(:,:,3))

  id_topog_step_4 = register_static_field ('ocean_model', 'topog_step_4', Grd%tracer_axes(1:2), &
       'topog_step_4', 'dimensionless', missing_value=missing_value, range=(/-1.0,1.0/))
  call diagnose_2d(Time, Grd, id_topog_step_4, topog_step(:,:,4))


  allocate (id_mixdownslope(num_prog_tracers))
  id_mixdownslope = -1

  do n=1,num_prog_tracers
     if(T_prog(n)%name == 'temp') then 
         id_mixdownslope(n) = register_diag_field ('ocean_model', &
              'mixdownslope_'//trim(T_prog(n)%name),              &
              Grd%tracer_axes(1:3), Time%model_time,              &
              'cp*mixdownslope*rho*dzt*temp',                     &
              'Watt/m^2', missing_value=missing_value, range=(/-1.e9,1.e9/))
     else
         id_mixdownslope(n) = register_diag_field ('ocean_model', 'mixdownslope_'//trim(T_prog(n)%name), &
              Grd%tracer_axes(1:3), Time%model_time,                                                     &
              'mixdownslope*rho*dzt*tracer for '//trim(T_prog(n)%name),                                  &
              trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e9,1.e9/))
     endif
  enddo

  call watermass_diag_init(Time,Dens)

  if (debug_this_module) then
      write(stdoutunit,*) ' '
      write(stdoutunit,*) '==Global sums from initialization of ocean_mixdownslope_mod== '
      call write_timestamp(Time%model_time)
      write(stdoutunit,'(a,es24.17)') &
     'slope_x          = ',mpp_global_sum(Dom%domain2d,slope_x(:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'slope_y          = ',mpp_global_sum(Dom%domain2d,slope_y(:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_slope_ex(1)= ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_slope_ex(:,:,1), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_slope_ex(2)= ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_slope_ex(:,:,2), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_slope_ex(3)= ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_slope_ex(:,:,3), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_slope_ex(4)= ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_slope_ex(:,:,4), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step(1)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,1), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step(2)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,2), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step(3)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,3), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step(4)    = ',mpp_global_sum(Dom%domain2d,topog_step(:,:,4), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step_ex(1) = ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_step_ex(:,:,1), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step_ex(2) = ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_step_ex(:,:,2), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step_ex(3) = ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_step_ex(:,:,3), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
     'topog_step_ex(4) = ',mpp_global_sum(Mixdownslope_domain%domain2d,topog_step_ex(:,:,4), global_sum_flag)
  endif

end subroutine ocean_mixdownslope_init
! </SUBROUTINE> NAME="ocean_mixdownslope_init"


!#######################################################################
! <SUBROUTINE NAME="mixdownslope">
!
! <DESCRIPTION>
! Compute thickness and density weighted tracer tendency [tracer*rho*m/s]
! due to exchange of tracer properties in regions where density-driven
! downslope transport is favorable. 
!
! Allow for exchanges to occur over horizontally 
! distant points, so long as the dense shallow parcel finds that it 
! will sit on the bottom of the horizontally adjacent columns.  Doing
! so requires a search algorithm, which requires some if-test logic
! as well as extended halos.  Note that the halos cannot be extended
! to larger than the size of the computational domain on a processor.
! This restriction limits the extent that we can search horizontally.
!
! The rates for the exchange are functions of the topographic slope
! and the density differences between parcels.  
!
! This scheme can be characterized as a slope convection based on 
! logic incorporated into the overflow and overexchange schemes. 
!
! </DESCRIPTION>
!
subroutine mixdownslope (Time, Thickness, T_prog, Dens, index_temp, index_salt)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens
  integer,                      intent(in)    :: index_temp
  integer,                      intent(in)    :: index_salt

  integer :: tau, taum1
  integer :: i, j, k, m, n, nt, npts 
  integer :: kmtij, ku, kd
  integer :: iip0, iip1, jjq0, jjq1
  integer :: iip1r, jjq1r
  integer :: nm1

  real    :: temp_so, salt_so, press, density_check
  real    :: weight, arg
  real    :: delta, delta_rho(ijhalo)
  real    :: mass_sum, tmix
  real    :: mixdownslope_total, mixdownslope_total_r
  real    :: tendency

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_mixdownslope_mod (mixdownslope): module must be initialized')
  endif 

  tau     = Time%tau
  taum1   = Time%taum1

  delta_rho(:) = 0.0

  ! extend some fields to extended domain 
  rho_ex   = 0.0
  mass_ex  = 0.0
  press_ex = 0.0
  temp_ex  = 0.0
  salt_ex  = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           rho_ex(i,j,k)   = Dens%rho(i,j,k,tau) 
           mass_ex(i,j,k)  = Thickness%rho_dzt(i,j,k,tau)*Grd%dat(i,j) 
           press_ex(i,j,k) = Dens%pressure_at_depth(i,j,k)           
           temp_ex(i,j,k)  = T_prog(index_temp)%field(i,j,k,tau)           
           salt_ex(i,j,k)  = T_prog(index_salt)%field(i,j,k,tau)           
        enddo
     enddo
  enddo

  call mpp_update_domains(rho_ex(:,:,:),   Mixdownslope_domain%domain2d, complete=.false.)
  call mpp_update_domains(mass_ex(:,:,:),  Mixdownslope_domain%domain2d, complete=.false.)
  call mpp_update_domains(press_ex(:,:,:), Mixdownslope_domain%domain2d, complete=.false.)
  call mpp_update_domains(temp_ex(:,:,:),  Mixdownslope_domain%domain2d, complete=.false.)
  call mpp_update_domains(salt_ex(:,:,:),  Mixdownslope_domain%domain2d, complete=.true.)

  ! compute details for cells participating in downslope mixing
  ! note:  "so" = "shallow ocean" cell...the central point  
  ! note:  "do" = cells within "deep ocean" columns 
  ! this part of the compuatation is independent of the tracer 

  kup_ex(:,:,:)              = 0
  kdw_ex(:,:,:,:)            = 0
  mixdownslope_frac(:,:,:,:) = 0.0

  do j=jsc,jec
     do i=isc,iec

        kmtij   = max(1,kmt_ex(i,j))
        temp_so = temp_ex(i,j,kmtij)
        salt_so = Dens%rho_salinity(i,j,kmtij,tau)

        ! 4-directions surrounding each cell 
        do m=1,4

           ! search for density favorable mixing 
           npts=0 
           nloop1_npts: do n=1,mixdownslope_npts

              iip0 = i+ip(n-1,m)
              jjq0 = j+jq(n-1,m)
              iip1 = i+ip(n,m)
              jjq1 = j+jq(n,m)

              ! check if downslope mixing is topographically possible
              if(topog_step_ex(iip0,jjq0,m)==1.0) then  

                  ! check if density of shallow ocean cell > density of deep ocean cell  
                  if(rho_ex(i,j,kmtij) > rho_ex(iip1,jjq1,kmtij)) then   

                      ! kdw = k-level in do-column where central cell is neutrally buoyant (or at bottom) 
                      kdw_ex(i,j,n,m) = 0
                      kloop1 : do k=kmtij+1,kmt_ex(iip1,jjq1)
                         press         = press_ex(iip1,jjq1,k)
                         density_check = density(salt_so,temp_so,press)
                         if(density_check > rho_ex(iip1,jjq1,k)) then 
                             kdw_ex(i,j,n,m) = k
                             delta_rho(n)    = density_check-rho_ex(iip1,jjq1,k) 
                         else 
                             exit kloop1
                         endif
                      enddo kloop1

                      ! check that the n-parcel has a k-level greater than the n-1 parcel.
                      ! if not, then do not mix with the n-column.
                      nm1 = max(n-1,1)
                      if(n > 1) then 
                          if(kdw_ex(i,j,n,m) <= kdw_ex(i,j,nm1,m)) then  
                              kdw_ex(i,j,n,m) = 0
                          endif  
                      endif

                      ! set strength of downslope mixing between central cell and n-parcel 
                      if(kdw_ex(i,j,n,m) > 0) then 

                          ! add this cell to the number of cells participating in mixing 
                          npts=npts+1

                          if(n==1) then 
                             kup_ex(i,j,m) = kmtij
                             delta = rho_ex(i,j,kmtij)-rho_ex(iip1,jjq1,kmtij)
                          else 
                             delta = delta_rho(n)
                          endif 

                          ! compute mixing weight as product of slope and density difference 
                          mixdownslope_frac(i,j,n,m) = topog_slope_ex(iip0,jjq0,m)*delta 

                      endif  ! kdw_ex(i,j,n,m) > 0 if-test 

                  endif      ! rho_ex(i,j,k) > rho_ex(iip1,jjq1,k) if-test 
              endif          ! topog_step_ex(iip0,jjq0,m)==1.0 if-test


              ! if kdw is not on the bottom, then exit n-loop since finished with search 
              if(kdw_ex(i,j,n,m) < kmt_ex(iip1,jjq1)) then 
                  exit nloop1_npts
              endif

           enddo nloop1_npts ! n-loop for horizontal search reaching out from central point

           
           ! place more weight on the farther points to 
           ! encourage tracer exchange further downslope.  
           if(mixdownslope_weight_far .and. npts > 0) then 
               weight = mixdownslope_frac(i,j,1,m)
               do n=1,npts
                  arg = float((n-npts)/mixdownslope_width)
                  mixdownslope_frac(i,j,n,m) = weight*exp(arg)
               enddo
           endif 

           ! normalize mixdownslope_frac to set fraction of
           ! a cell that is being mixed with central cell. 
           if(npts > 0) then 
               mixdownslope_total=0.0
               do n=1,npts
                  mixdownslope_total = mixdownslope_total + mixdownslope_frac(i,j,n,m)   
               enddo
               mixdownslope_total_r = 1.0/(mixdownslope_total+epsln)
               do n=1,npts
                  mixdownslope_frac(i,j,n,m) = mixdownslope_frac(i,j,n,m)*mixdownslope_total_r
               enddo
           endif

           ! if no exchange points 
           if(npts==0) then 
               mixdownslope_frac(i,j,:,m) = 0.0
           endif

        enddo  ! m-loop
     enddo     ! i-loop
  enddo        ! j-loop


  ! extend arrays to wide halos 
  call mpp_update_domains(kup_ex(:,:,:),              Mixdownslope_domain%domain2d)
  call mpp_update_domains(kdw_ex(:,:,:,:),            Mixdownslope_domain%domain2d)
  call mpp_update_domains(mixdownslope_frac(:,:,:,:), Mixdownslope_domain%domain2d)


  ! compute tracer tendency for cells participating in downslope mixing 
  do nt=1,num_prog_tracers 

     ! place tracer concentration in the wider array tracer_ex 
     if(nt==index_temp) then 
         do k=1,nk
            do j=jsc-ijhalo,jec+ijhalo
               do i=isc-ijhalo,iec+ijhalo
                  tracer_ex(i,j,k) = temp_ex(i,j,k)
               enddo
            enddo
         enddo
     elseif(nt==index_salt) then 
         do k=1,nk
            do j=jsc-ijhalo,jec+ijhalo
               do i=isc-ijhalo,iec+ijhalo
                  tracer_ex(i,j,k) = salt_ex(i,j,k)
               enddo
            enddo
         enddo
     else
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tracer_ex(i,j,k) = T_prog(nt)%field(i,j,k,tau)
               enddo
            enddo
         enddo
         call mpp_update_domains (tracer_ex(:,:,:), Mixdownslope_domain%domain2d)
     endif

     ! compute tendency, noting that each (i,j,k) cell can generally be part
     ! of mixing as either a central shallow cell (n=0), or as one of the deep cells (n>0). 
     ! for the tendency, we compute a mixing tracer concentration between two cells, 
     ! and then back-out a tendency which is placed in tend_mix.

     tend_mix(:,:,:) = 0.0
     do j=jsc,jec
        do i=isc,iec
           do m=1,4
              do n=1,mixdownslope_npts

                 ! i,j is central cell at k=ku
                 iip1 = i+ip(n,m)
                 jjq1 = j+jq(n,m)
                 ku   = kup_ex(i,j,m) 
                 kd   = kdw_ex(i,j,n,m) 
                 if(ku > 0 .and. kd > ku) then 
                     mass_sum = mixdownslope_frac_central*mass_ex(i,j,ku) &
                               +mixdownslope_frac(i,j,n,m)*mass_ex(iip1,jjq1,kd)
                     tmix = (mixdownslope_frac_central*mass_ex(i,j,ku)*tracer_ex(i,j,ku)               &
                            +mixdownslope_frac(i,j,n,m)*mass_ex(iip1,jjq1,kd)*tracer_ex(iip1,jjq1,kd)) &
                             /mass_sum
                     tendency = dtimer*mixdownslope_frac_central*mass_ex(i,j,ku)/Grd%dat(i,j) &
                                *(tmix-tracer_ex(i,j,ku))
                     tend_mix(i,j,ku) = tend_mix(i,j,ku) + tendency 
                 endif 

                 ! i,j is deep cell at k=kd
                 iip1r = i-ip(n,m)
                 jjq1r = j-jq(n,m)
                 ku    = kup_ex(iip1r,jjq1r,m) 
                 kd    = kdw_ex(iip1r,jjq1r,n,m) 
                 if(ku > 0 .and. kd > ku) then 
                     mass_sum = mixdownslope_frac_central*mass_ex(iip1r,jjq1r,ku) &
                               +mixdownslope_frac(iip1r,jjq1r,n,m)*mass_ex(i,j,kd)
                     tmix = (mixdownslope_frac_central                          &
                             *mass_ex(iip1r,jjq1r,ku)*tracer_ex(iip1r,jjq1r,ku) &
                            +mixdownslope_frac(iip1r,jjq1r,n,m)                 &
                             *mass_ex(i,j,kd)*tracer_ex(i,j,kd))                & 
                            /mass_sum
                     tendency = dtimer*mixdownslope_frac(iip1r,jjq1r,n,m)*mass_ex(i,j,kd)/Grd%dat(i,j) &
                                *(tmix-tracer_ex(i,j,kd))
                     tend_mix(i,j,kd) = tend_mix(i,j,kd) + tendency 
                 endif 

              enddo ! n-loop 
           enddo    ! m-loop
        enddo       ! i-loop
     enddo          ! j-loop


     ! fill tracer tendency array 
     do k=1,nk
        do j=jsc,jec   
           do i=isc,iec
              T_prog(nt)%th_tendency(i,j,k) = T_prog(nt)%th_tendency(i,j,k) &
                                             +tend_mix(i,j,k)*mixdownslope_mask(i,j)
           enddo
        enddo
     enddo

     if(id_mixdownslope(nt) > 0) then
         wrk1(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  wrk1(i,j,k) = tend_mix(i,j,k)*mixdownslope_mask(i,j)*T_prog(nt)%conversion
               enddo
            enddo
         enddo
         call diagnose_3d(Time, Grd, id_mixdownslope(nt), wrk1(:,:,:))
     endif

     if(nt==index_temp) then 
         theta_tend(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  theta_tend(i,j,k) = tend_mix(i,j,k)*mixdownslope_mask(i,j)
               enddo
            enddo
         enddo
     endif
     if(nt==index_salt) then 
         salt_tend(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  salt_tend(i,j,k) = tend_mix(i,j,k)*mixdownslope_mask(i,j)
               enddo
            enddo
         enddo
     endif

     if(debug_this_module) then 
         write(stdoutunit,*) ' '
         write(stdoutunit,*) '==Global sums for tendency from ocean_mixdownslope_mod== '
         call write_timestamp(Time%model_time)
         do k=1,nk
            area_k(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,k)  
            tend_mix(:,:,k) = dtime*tend_mix(:,:,k)*area_k(:,:)*T_prog(nt)%conversion
            write(stdoutunit,'(a,i2,a,i2,a,es24.17)') 'tend_mix(',nt,',',k,') = ',&
                 mpp_global_sum(Dom%domain2d,tend_mix(:,:,k))
         enddo
         write(stdoutunit,'(a,i2,a,es24.17)') &
              'tend_mix(',nt,')  = ',&
              mpp_global_sum(Dom%domain2d,tend_mix(:,:,:), global_sum_flag)
     endif


  enddo  ! nt-end for num_prog_tracers 

  call watermass_diag(Time, Dens)

  if (debug_this_module) then
      write(stdoutunit,*) ' '
      write(stdoutunit,*) '==Global sums from ocean_mixdownslope_mod== '
      call write_timestamp(Time%model_time)

      write(stdoutunit,'(a,es24.17)') &
           'rho_ex            = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,rho_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'mass_ex        = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,mass_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'press_ex          = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,press_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'temp_ex           = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,temp_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,es24.17)') &
           'salt_ex           = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,salt_ex(:,:,:), global_sum_flag)
      write(stdoutunit,'(a,i24)') &
           'kmt_ex            = ',&
           mpp_global_sum(Mixdownslope_domain%domain2d,kmt_ex(:,:), global_sum_flag)

      do m=1,4
         write(stdoutunit,'(a,i24)') &
              'kup_ex            = ',  &
              mpp_global_sum(Mixdownslope_domain%domain2d,kup_ex(:,:,m), global_sum_flag)
         do n=1,mixdownslope_npts
            write(stdoutunit,'(a,es24.17)') &
                 'mixdownslope_frac  = ',  &
                 mpp_global_sum(Mixdownslope_domain%domain2d,mixdownslope_frac(:,:,n,m), global_sum_flag)
            write(stdoutunit,'(a,i24)') &
                 'kdw_ex            = ',  &
                 mpp_global_sum(Mixdownslope_domain%domain2d,kdw_ex(:,:,n,m), global_sum_flag)
         enddo
      enddo

      write(stdoutunit,*) ' '
      write(stdoutunit,*) '==Global check sums from ocean_mixdownslope_mod== '
      call write_timestamp(Time%model_time)
      call write_chksum_3d('rho_ex', rho_ex(COMP,:))
      call write_chksum_3d('mass_ex', mass_ex(COMP,:))
      call write_chksum_3d('press_ex', press_ex(COMP,:))
      call write_chksum_3d('temp_ex', temp_ex(COMP,:))
      call write_chksum_3d('salt_ex', salt_ex(COMP,:))
      call write_chksum_2d_int('kmt_ex', kmt_ex(COMP))
      do m=1,4
         call write_chksum_2d_int('kup_ex', kup_ex(COMP,m))
         do n=1,mixdownslope_npts
            call write_chksum_2d_int('kdw_ex', kdw_ex(COMP,n,m))
            call write_chksum_2d('mixdownslope_frac', mixdownslope_frac(COMP,n,m))
         enddo
      enddo

  endif

end subroutine mixdownslope
! </SUBROUTINE> NAME="mixdownslope"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 

  id_neut_rho_mixdown = register_diag_field ('ocean_model', 'neut_rho_mixdown',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
    'update of locally referenced potential density from mixdowslope',         &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_mixdown > 0) compute_watermass_diag = .true.

  id_wdian_rho_mixdown = register_diag_field ('ocean_model', 'wdian_rho_mixdown',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'dianeutral mass transport due to mixdowslope',                              &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_mixdown > 0) compute_watermass_diag = .true.

  id_tform_rho_mixdown = register_diag_field ('ocean_model', 'tform_rho_mixdown',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'watermass transform due to mixdowslope on levels (pre-layer binning)',      &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_mixdown > 0) compute_watermass_diag = .true.

  id_neut_rho_mixdown_on_nrho = register_diag_field ('ocean_model',                    &
    'neut_rho_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
    'update of local ref potrho from mixdownslope as binned to neutral density layers',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_rho_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_rho_mixdown_on_nrho = register_diag_field ('ocean_model',                  &
   'wdian_rho_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,           &
   'dianeutral mass transport due to mixdowslope as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_rho_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_rho_mixdown_on_nrho = register_diag_field ('ocean_model',            &
   'tform_rho_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,     &
   'watermass transform due to mixdowslope as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_rho_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_eta_tend_mixdown= -1          
  id_eta_tend_mixdown= register_diag_field ('ocean_model','eta_tend_mixdown',&
       Grd%tracer_axes(1:2), Time%model_time,                                &
       'non-Bouss steric sea level tendency from mixdown tendency', 'm/s',   &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_mixdown > 0) compute_watermass_diag=.true.

  id_eta_tend_mixdown_glob= -1          
  id_eta_tend_mixdown_glob= register_diag_field ('ocean_model', 'eta_tend_mixdown_glob',&
       Time%model_time,                                                                 &
       'global mean non-bouss steric sea level tendency from mixdown tendency',         &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_mixdown_glob > 0) compute_watermass_diag=.true.


  ! temp contributions 
  id_neut_temp_mixdown = register_diag_field ('ocean_model', 'neut_temp_mixdown',  &
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'temp related update of locally referenced potential density from mixdowslope',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_mixdown > 0) compute_watermass_diag = .true.

  id_wdian_temp_mixdown = register_diag_field ('ocean_model', 'wdian_temp_mixdown',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'temp related dianeutral mass transport due to mixdowslope',                   &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_mixdown > 0) compute_watermass_diag = .true.

  id_tform_temp_mixdown = register_diag_field ('ocean_model', 'tform_temp_mixdown',     &
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'temp related watermass transform due to mixdowslope on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_mixdown > 0) compute_watermass_diag = .true.

  id_neut_temp_mixdown_on_nrho = register_diag_field ('ocean_model',                               &
   'neut_temp_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                        &
   'temp related update of local ref potrho from mixdownslope as binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_temp_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_temp_mixdown_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_temp_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'temp related dianeutral mass transport due to mixdowslope as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_temp_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_temp_mixdown_on_nrho = register_diag_field ('ocean_model',                          &
     'tform_temp_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
     'temp related watermass transform due to mixdowslope as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_temp_mixdown_on_nrho > 0) compute_watermass_diag = .true.


  ! salt contributions 
  id_neut_salt_mixdown = register_diag_field ('ocean_model', 'neut_salt_mixdown',  &
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'salt related update of locally referenced potential density from mixdowslope',&
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_mixdown > 0) compute_watermass_diag = .true.

  id_wdian_salt_mixdown = register_diag_field ('ocean_model', 'wdian_salt_mixdown',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'salt related dianeutral mass transport due to mixdowslope',                   &
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_mixdown > 0) compute_watermass_diag = .true.

  id_tform_salt_mixdown = register_diag_field ('ocean_model', 'tform_salt_mixdown',     &
    Grd%tracer_axes(1:3), Time%model_time,                                              &
    'salt related watermass transform due to mixdowslope on levels (pre-layer binning)',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_mixdown > 0) compute_watermass_diag = .true.

  id_neut_salt_mixdown_on_nrho = register_diag_field ('ocean_model',                               &
   'neut_salt_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                        &
   'salt related update of local ref potrho from mixdownslope as binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_neut_salt_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_wdian_salt_mixdown_on_nrho = register_diag_field ('ocean_model',                               &
    'wdian_salt_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                       &
    'salt related dianeutral mass transport due to mixdowslope as binned to neutral density layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_wdian_salt_mixdown_on_nrho > 0) compute_watermass_diag = .true.

  id_tform_salt_mixdown_on_nrho = register_diag_field ('ocean_model',                          &
     'tform_salt_mixdown_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                 &
     'salt related watermass transform due to mixdowslope as binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  if(id_tform_salt_mixdown_on_nrho > 0) compute_watermass_diag = .true.


  allocate(theta_tend(isd:ied,jsd:jed,nk))
  allocate(salt_tend(isd:ied,jsd:jed,nk))
  theta_tend(:,:,:) = 0.0
  salt_tend(:,:,:)  = 0.0

  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_mixdownslope_mod w/ compute_watermass_diag=.true. to compute some watermass diagnostics.'  
  endif 


end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from mixdownslope on the watermass transformation.
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens)

  type(ocean_time_type),         intent(in) :: Time
  type(ocean_density_type),      intent(in) :: Dens

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_mixdownslope (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 

  tau = Time%tau 

  ! rho diagnostics = sum of temp + salt contributions 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)* &
             (Dens%drhodT(i,j,k)*theta_tend(i,j,k)+Dens%drhodS(i,j,k)*salt_tend(i,j,k))
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend 
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_mixdown, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_mixdown, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_mixdown, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_mixdown_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_mixdown_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_mixdown_on_nrho, wrk4)

  if(id_eta_tend_mixdown > 0 .or. id_eta_tend_mixdown_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_mixdown, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_mixdown_glob, eta_tend, cellarea_r)
  endif


  ! temp contributions 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*theta_tend(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_mixdown, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_mixdown, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_mixdown, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_mixdown_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_mixdown_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_mixdown_on_nrho, wrk4)

  ! salinity contributions 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodS(i,j,k)*salt_tend(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)  
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_mixdown, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_mixdown, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_mixdown, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_mixdown_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_mixdown_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_mixdown_on_nrho, wrk4)


end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_mixdownslope_mod

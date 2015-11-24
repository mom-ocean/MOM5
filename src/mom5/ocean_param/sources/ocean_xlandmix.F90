module ocean_xlandmix_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison 
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> K. Dixon 
!</REVIEWER>
!
!<OVERVIEW>
! Tracer and mass source from cross-land mixing.
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted and density weighted tendency
! [tracer*(kg/m^3)*meter/sec] of tracer and mass associated with 
! cross-land mixing of points. 
!
! In particular, this module provides interaction between bodies of
! water separated by land in coarse models (e.g., Med-Atlantic).  
! Interaction consists of mixing tracer and mass between water 
! columns found on each side of the land barrier.
!
! To conserve total tracer when using xlandmix with k=1, 
! also need to connect mass between remote bodies of water
! via xlandmix of eta.  
!
! When connecting cells with k>1, we do not mix mass. This 
! is ensured by making the time tendency only a function of the 
! difference in tracer concentration. 
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2004)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <NOTE>
! Algorithm ensures both total tracer and total mass is conserved.
! </NOTE>
!
! <NOTE>
! 2D domain decomposition is implemented according to following
! notions.  If the two points in an xlandmix pair live within halo 
! of each other (i.e., within same local_domain), 
! then no added mpp communication required. However, more generally 
! the two points live further away.  In this case, xland_domain
! is defined so that its halos incorporate the maximum separation 
! of xlandmix points.  New tracer, eta, and grid arrays
! are defined over this extended xland_domain.  This added domain
! size will come at some computational cost, so it is advantageous
! to limit the separation between points within an xlandmix pair. 
! </NOTE>
!
! <NOTE>
! The current implementation of xlandmix has not been generalized 
! to allow for communication between points separated by the tripolar fold.  
! The problem is in the logic used to compute xland_domain.  
! There is nothing fundamental limiting a generalization of the
! logic used to generate xland_domain.
! </NOTE>
!
! <NOTE>
! Many of the user specified values given in USER INPUT are
! model dependent since unresolved straits can become resolved 
! in finer mesh models. 
! </NOTE>
!
! <NOTE>
! Algorithm originally developed by Keith Dixon (Keith.Dixon) 
! for rigid lid and full cells and one processor (i.e., MOM1).  
! Algorithm extended to free surface and partial cells and 2D 
! domain decomposition by S.M.Griffies (Stephen.Griffies).
! Further extensions to generalized vertical coordinates by 
! Stephen.Griffies
! </NOTE>
!
! <NOTE>
! BUG FIX 30 May 2011. 
! Trying to mix across the ends of the global domain results in incorrect indexing.
! Local indices need to be incremented/decremented by the global extent.
! Russell.Fiedler@csiro.au
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_xlandmix_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Needs to be true in order to use this scheme. Default is false.
!  </DATA> 
!  <DATA NAME="xlandmix_kmt" TYPE="logical">
!  To allow xlandmixing to occur at k=kmt cell.  
!  Default is xlandmix_kmt=.false. 
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose initialization information.
!  </DATA> 
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field
use field_manager_mod, only: MODEL_OCEAN, parse, find_field_index
use field_manager_mod, only: get_field_methods, method_type, get_field_info
use fms_mod,           only: stdout, stdlog, FATAL, NOTE, WARNING
use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
use mpp_domains_mod,   only: mpp_update_domains
use mpp_mod,           only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices, get_global_indices, set_ocean_domain
use ocean_parameters_mod, only: missing_value, rho0, rho0r, onehalf
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_time_type, ocean_density_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_options_type
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_v
use ocean_util_mod,       only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_tracer_util_mod,only: diagnose_3d_rho


implicit none

private 

public ocean_xlandmix_init
public xlandmix
private xlandmix_mass
private watermass_diag_init
private watermass_diag

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()


! Variables set inside subroutine ocean_xlandmix_init
! See comments in ocean_xlandmix_init for description
integer, dimension(:,:), allocatable :: ixland     
integer, dimension(:,:), allocatable :: jxland     
integer, dimension(:,:), allocatable :: kxland     
real, dimension(:), allocatable      :: vxland     

! number of xland mixing pairs 
integer :: nxland=0 

! Variables calculated from user input
real, dimension(:,:), allocatable    :: depth_static   ! time independent depth of columns (m)
real, dimension(:,:), allocatable    :: gamma_xland    ! inverse damping time (s^-1)

integer                              :: xland_halo     ! halo size needed to perform xlandmix
integer                              :: kxland_min     ! min of the k-levels for the nxland pairs
integer                              :: kxland_max     ! max of the k-levels for the nxland pairs
type(ocean_domain_type), save        :: Xland_domain   ! domain for xlandmix 

! variables defined with halos given by xland_halo
real, dimension(:,:), allocatable    :: xland_mass      ! mass source associated with xlandmix  kg/(m^2*sec)
real, dimension(:,:), allocatable    :: xland_mass_pair ! mass source on space of xlandmix pairs kg/(m^2*sec)
real, dimension(:,:), allocatable    :: xtx             ! zonal position of tracer points (degrees)
real, dimension(:,:), allocatable    :: ytx             ! meridional position of tracer points (degrees)
real, dimension(:,:), allocatable    :: datx            ! area of tracer cell (m^2)
real, dimension(:,:,:), allocatable  :: rho_dzt_x       ! mass per area of tracer cell (kg/m^2)

integer                              :: unit=6         !processor zero writes to unit 6

logical, dimension(:), allocatable   :: error_xland    !for checking that all xland points are OK.  


! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 

! for global area normalization
real    :: cellarea_r

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xland
integer :: id_xland_mass=-1

integer :: id_neut_rho_xmix          =-1
integer :: id_neut_rho_xmix_on_nrho  =-1
integer :: id_wdian_rho_xmix         =-1
integer :: id_wdian_rho_xmix_on_nrho =-1
integer :: id_tform_rho_xmix         =-1
integer :: id_tform_rho_xmix_on_nrho =-1

integer :: id_eta_tend_xmix     =-1
integer :: id_eta_tend_xmix_glob=-1

integer :: id_neut_temp_xmix          =-1
integer :: id_neut_temp_xmix_on_nrho  =-1
integer :: id_wdian_temp_xmix         =-1
integer :: id_wdian_temp_xmix_on_nrho =-1
integer :: id_tform_temp_xmix         =-1
integer :: id_tform_temp_xmix_on_nrho =-1

integer :: id_neut_salt_xmix          =-1
integer :: id_neut_salt_xmix_on_nrho  =-1
integer :: id_wdian_salt_xmix         =-1
integer :: id_wdian_salt_xmix_on_nrho =-1
integer :: id_tform_salt_xmix         =-1
integer :: id_tform_salt_xmix_on_nrho =-1

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1

character(len=128) :: version=&
       '$Id: ocean_xlandmix.F90,v 20.0 2013/12/14 00:16:30 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: isg, ieg, jsg, jeg
integer :: ixl, jxl, ixl2, jxl2
integer :: kmt_offset=1

logical :: verbose_init          = .true.
logical :: use_this_module       = .false.
logical :: module_is_initialized = .FALSE.
logical :: xlandmix_kmt          = .false.

namelist /ocean_xlandmix_nml/ verbose_init, use_this_module, xlandmix_kmt

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_xlandmix_init">
!
! <DESCRIPTION>
!    Initial set up for crossland mixing of tracers 
!
!    (i,j,k) locations of crossland mixing points are set via 
!    field table entries.
!
!    Time invariant crossland mixing rates (in units of m**3/sec) are 
!    set via field table entries. 
!
!    Checks are performed to ensure that the crossland mixing
!    grid locations are valid according to model configuration.
!
!    A summary of the locations of crossland mixing points is written out.
!
!    User specified inputs in "USER INPUT" section:
!
!    ixland and jxland = user specified nxland pairs of i,j grid locations
!
!    kxland = user specified uppermost (ktop=kxland(n,1)) and deepest 
!             (kbot=kxland(n,2)) levels over which crossland mixing will 
!             be done for each pair of crossland mixing points.
!
!    vxland = user specified time invariant rates of crossland mixing (m3/sec)
!             Equivalent to "the flow to the east = the flow to the west"
!             Dynamic vxland is not available in MOM. 
!
!    NOTE: for ixland, jxland, and kxland pairs, values of the
!    (n,1) element should be < the corresponding (n,2) value.
! </DESCRIPTION>
!
subroutine ocean_xlandmix_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, dtime)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_density_type),     intent(in)         :: Dens
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  type(ocean_options_type),     intent(inout)      :: Ocean_options
  real,                         intent(in)         :: dtime

  type(method_type), allocatable, dimension(:) :: xland_methods
  character(len=32) :: fld_type, fld_name
  logical           :: error_check=.false.

  integer :: i, n, nt, model, imax, imin
  integer :: nxl, lx, xland_halo_test
  integer :: parse_ok
  integer :: kmt_max
  integer :: io_status, ioun, ierr
  real    :: ztop 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if (module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandmix_mod (ocean_xlandmix_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  call get_global_indices(Domain,isg,ieg,jsg,jeg)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid
  
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

  num_prog_tracers = size(T_prog(:))
  do nt=1,num_prog_tracers
     if (T_prog(nt)%name == 'temp') index_temp = nt
     if (T_prog(nt)%name == 'salt') index_salt = nt
  enddo

  n=find_field_index(MODEL_OCEAN,'xland_mix')

  if (n < 1) then 
     Ocean_options%xlandmix = 'Did NOT use xlandmix.'
     return
  endif 

  call get_field_info(n,fld_type,fld_name,model,nxland)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_xlandmix_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_xlandmix_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_xlandmix_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_xlandmix_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_xlandmix_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_xlandmix_nml)

  if(use_this_module) then 
      if(nxland == 0) then 
          write(stdoutunit, *) &
          '==>Warning: use_this_module=.true. but no mixing points chosen. Is this what you wish?'
          call mpp_error(WARNING, &
          '==>NOTE from ocean_xlandmix_mod: No cross-land mixing points chosen')
          Ocean_options%xlandmix = 'Used xlandmix, but set zero cross-land points, so module does nothing.'
          return
      else
          call mpp_error(NOTE, &
          '==>Note from ocean_xlandmix_mod (ocean_xlandmix_init): USING this module')
          write(stdoutunit, *) &
          '==>Note: Using xlandmix to connect tracers and mass between non-local ocean cells.'
          Ocean_options%xlandmix = 'Used xlandmix to transfer tracer between inland seas and open ocean.'
      endif
  else 
      call mpp_error(NOTE, '==>Note from ocean_xlandmix_mod: NOT USING this module')
      if(nxland > 0) then 
          call mpp_error(WARNING, &
          '==>Warning: nxland > 0, but use_this_module=.false. NOT using xlandmix. Is this what you wish?')
      endif
      Ocean_options%xlandmix = 'Did NOT use xlandmix.'
      return
  endif

  write(stdoutunit,'(/a,f10.2)') &
  '==> Note from ocean_xlandmix_mod: using forward time step of (secs)', dtime 

  if(xlandmix_kmt) then 
    kmt_offset=0
    write(stdoutunit,'(/a)') &
    '==> Note from ocean_xlandmix_mod: allowing xlandmix to occur at k=kmt.'
  else 
    kmt_offset=1
  endif 

  allocate(xland_methods(nxland))
  allocate (ixland(nxland,2))
  ixland(:,:) = 0
  allocate (jxland(nxland,2))
  jxland(:,:) = 0
  allocate (kxland(nxland,2))
  kxland(:,:) = 0
  allocate (vxland(nxland))
  vxland(:)   = 0.0
  allocate (depth_static(nxland,2))
  depth_static(:,:) = 0.0
  allocate (gamma_xland(nxland,2))
  gamma_xland(:,:) = 0.0
  allocate (error_xland(nxland))
  error_xland(:) = .false.

  call get_field_methods(n,xland_methods)
  do i=1, nxland
     parse_ok = parse(xland_methods(i)%method_control,'ixland_1',ixland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "ixland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'ixland_2',ixland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "ixland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_1',jxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "jxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_2',jxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "jxland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_1',kxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "kxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_2',kxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "kxland_2" error')
     parse_ok   = parse(xland_methods(i)%method_control,'vxland',vxland(i))
     if (parse_ok == 0) call mpp_error(FATAL,&
     '==>Error ocean_xlandmix_mod: xland table entry "vxland" error')
  enddo

  if(Grid%tripolar) then 
     write(stdoutunit,'(a)') &
     '==>Warning: xlandmix has not been implemented to connect points across the tripolar fold.'
  endif

  xland_halo = min(isc-isd,ied-iec,jsc-jsd,jed-jec)
  kxland_min = nk
  kxland_max = 0

  do nxl=1,nxland

     ! check for invalid crossland mixing grid locations

     if (kxland(nxl,1) > kxland(nxl,2)) then
        write (unit,994) nxl, kxland(nxl,1), nxl, kxland(nxl,2)
        error_xland(nxl) = .true.
     endif
     if (kxland(nxl,1) < kxland_min ) kxland_min = kxland(nxl,1)  
     if (kxland(nxl,2) > kxland_max ) kxland_max = kxland(nxl,2)  

     do lx = 1,2
        if (ixland(nxl,lx) < 1 .or. ixland(nxl,lx) > Grd%ni) then
           write(unit,991) nxl, lx, ixland(nxl,lx)
           error_xland(nxl) = .true.
        endif
        if (jxland(nxl,lx) < 1 .or. jxland(nxl,lx) > Grd%nj ) then
           write(unit,992) nxl, lx, jxland(nxl,lx)
           error_xland(nxl) = .true.
        endif
        if (kxland(nxl,lx) < 1 .or. kxland(nxl,lx) > Grd%nk ) then
           write(unit,993) nxl, lx, kxland(nxl,lx)
           error_xland(nxl) = .true.
        endif
     enddo

     ! determine size for xland_halo
     if(Grd%cyclic_x) then
       imax = max(ixland(nxl,1),ixland(nxl,2))
       imin = min(ixland(nxl,1),ixland(nxl,2))
       xland_halo_test = min( abs(imax-imin), abs(imax-Grd%ni-imin))
     else
       xland_halo_test = abs(ixland(nxl,1)-ixland(nxl,2))
     endif
     if (xland_halo_test > xland_halo ) xland_halo = xland_halo_test

     xland_halo_test = abs(jxland(nxl,1)-jxland(nxl,2))
     if (xland_halo_test > xland_halo ) xland_halo = xland_halo_test

  enddo  ! end of do-loop nxl=1,nxland

  if(xland_halo == 0) then 
     call mpp_error(FATAL,&
     '==>Error in ocean_xlandmix_mod (ocean_xlandmix_init): xland_halo=0 => problems defining xlandmix points.')
  endif
 
  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
    write(stdoutunit,*) &
    'Defining extra tracer arrays for xland_domain over k-levels ',kxland_min,' through ',kxland_max 
  endif 

  allocate (xland_mass_pair(nxland,2))
  xland_mass_pair(:,:)  = 0.0


  ! define arrays and xland_domain
  if(xland_halo <= Dom%xhalo .and. xland_halo <= Dom%yhalo) then 

     write(stdoutunit,*)'The model local computational domain has a x,y halo = ',Dom%xhalo,Dom%yhalo
     write(stdoutunit,*)'This is large enough for chosen points in xlandmix,'
     write(stdoutunit,*)'where we need to have a halo of size = ', xland_halo

     ! time independent arrays 
     allocate (datx(isd:ied,jsd:jed))
     allocate (xtx(isd:ied,jsd:jed))
     allocate (ytx(isd:ied,jsd:jed))
     xtx(:,:)         = Grd%xt(:,:)
     ytx(:,:)         = Grd%yt(:,:)
     datx(:,:)        = Grd%dat(:,:)

     ! time dependent arrays 
     allocate (rho_dzt_x(isd:ied,jsd:jed,nk))
     allocate (xland_mass(isd:ied,jsd:jed))
     rho_dzt_x(:,:,:) = 0.0
     xland_mass(:,:)  = 0.0
  endif 

  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
     write(stdoutunit,*)'The model local computational domain has a x,y halo = ',Dom%xhalo, Dom%yhalo
     write(stdoutunit,*)'This is smaller than the halo required for xlandmix.'
     write(stdoutunit,*)'For xlandmix, will define a new domain type with halo = ', xland_halo

     if(xland_halo > (Dom%xhalo+6)) then 
       write(stdoutunit,*)'=>Warning: ocean_xlandmix_init determined xland_halo = ', xland_halo
       write(stdoutunit,*)'           Are you sure you wish to connect such distant locations?'
       write(stdoutunit,*)'           Are you instead trying to connect points across a tripolar fold?'
       write(stdoutunit,*)'           If so, then presently this has yet to be coded into xlandmix.'
       write(stdoutunit,*)'           Alternative xlandmix points must be taken, or new code implemented.'
     endif

     ! define xlandmix domain 
     call set_ocean_domain(Xland_domain,Grd,xhalo=xland_halo,yhalo=xland_halo,&
                           name='xlandmix',maskmap=Dom%maskmap)

     ! time independent arrays defined over the xland_domain region 
     allocate (datx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     allocate (xtx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     allocate (ytx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     datx(:,:) = 0.0
     xtx(:,:)  = 0.0
     ytx(:,:)  = 0.0

     datx(isc:iec,jsc:jec) = Grd%dat(isc:iec,jsc:jec)
     xtx(isc:iec,jsc:jec)  = Grd%xt(isc:iec,jsc:jec)
     ytx(isc:iec,jsc:jec)  = Grd%yt(isc:iec,jsc:jec)
     call mpp_update_domains (datx(:,:), Xland_domain%domain2d)
     call mpp_update_domains (xtx(:,:), Xland_domain%domain2d)
     call mpp_update_domains (ytx(:,:), Xland_domain%domain2d)

     ! time dependent arrays defined over the xland_domain region 
     allocate (rho_dzt_x(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo,nk))
     allocate (xland_mass(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     rho_dzt_x(:,:,:) = 0.0
     xland_mass(:,:)  = 0.0
  endif

! Mixing across Straits of Gibraltar in CSIRO MK3.5 model resulted in 
! exceeding bounds with the Dec2009 version of MOM4p1. 
! ixl needs to be changed to reflect cyclic BCs.
! Earlier checks were made to see if ixl was negative or too 
! large wrt the global domain. However we want ixl in the current 
! expanded data domain of xtx.
  if(Grd%cyclic_x) then
     do nxl=1,nxland
        if( at_least_one_in_comp_domain(nxl) ) then
            do lx = 1,2
               if(.not. on_processor(nxl, lx)) then               ! Not in domain. This is WRONG!
                  ixl = ixland(nxl,lx)-Dom%ioff
                  if(ixl < lbound(xtx,1) ) then                   ! Gone outside lower index
                     ixl2 = ixland(nxl,lx) + ieg
                     write(unit,*) 'xlandmix init: ixland outside local domain i=',lbound(xtx,1),':',ubound(xtx,1),   &
                                   ' Changing index from ', ixland(nxl,lx) ,' to ', ixl2
                     ixland(nxl,lx) = ixl2
                  else if ( ixl > ubound(xtx,1) ) then             ! Gone outside upper index
                     ixl2 = ixland(nxl,lx) - ieg
                     write(unit,*) 'xlandmix init: ixland outside local domain i=',lbound(xtx,1),':',ubound(xtx,1),   &
                                   ' Changing index from ', ixland(nxl,lx) ,' to ', ixl2
                     ixland(nxl,lx) = ixl2
                  endif
                  if(.not. on_processor(nxl, lx)) then
                    write(unit,*) 'xlandmix init FAILED?: ixl =', ixl, ' Bounds ',isc+Dom%ioff-xland_halo,':',iec+Dom%ioff+xland_halo
                    error_xland(nxl)=.true.
                  endif
               endif
             enddo
         endif
      enddo
  endif

  do nxl=1,nxland

     ! check for attempts to mix land rather than sea.
     ! Restrict xland mixing to be above the bottom-most level in a column.      
     do lx = 1,2
        if(on_comp_domain(nxl, lx)) then
           ixl = ixland(nxl,lx)-Dom%ioff
           jxl = jxland(nxl,lx)-Dom%joff
           kmt_max = max(0,Grid%kmt(ixl,jxl)-kmt_offset)
           if (kxland(nxl,2) > kmt_max) then
              write (unit,192) nxl, kxland(nxl,2), ixland(nxl,lx),jxland(nxl,lx),Grid%kmt(ixl,jxl)
              error_xland(nxl)=.true.
           endif
        endif
     enddo

     ! write out summary information for this pair of crossland mixing points
     if (kxland(nxl,1) == 1) then
        ztop = 0.0
     else
        ztop = Grid%zw(kxland(nxl,1)-1)
     endif

     if(at_least_one_in_comp_domain(nxl)) then
         ixl = ixland(nxl,1)-Dom%ioff
         jxl = jxland(nxl,1)-Dom%joff
         ixl2 = ixland(nxl,2)-Dom%ioff
         jxl2 = jxland(nxl,2)-Dom%joff
         if(verbose_init) then 
             write(unit,191) nxl, ixland(nxl,1), jxland(nxl,1) &
              ,xtx(ixl,jxl), ytx(ixl,jxl) &
              ,ixland(nxl,2), jxland(nxl,2) &
              ,xtx(ixl2,jxl2), ytx(ixl2,jxl2) &
              ,kxland(nxl,1), kxland(nxl,2), ztop, Grid%zw(kxland(nxl,2))
         endif 
         call xland_check (nxl, dtime, error_check)
         if(error_check) then 
           error_xland(nxl) = .true.
         endif 
     endif

  enddo

  ! Bring model down here if there are problems with xlandmix points.  
  do nxl=1,nxland
    if(error_xland(nxl)) then
      call mpp_error(FATAL, &
      '==>Error detected in ocean_xlandmix_mod: Use "grep -i error" to find ALL errors.')
    endif 
  enddo 

  ! register for diag_manager 
  allocate (id_xland(num_prog_tracers))
  id_xland = -1

  do n=1,num_prog_tracers
     if(T_prog(n)%name == 'temp') then 
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xland', &
                       Grd%tracer_axes(1:3), Time%model_time,                              &
                       'xlandmix*cp*rho*dzt*temp', 'Watt/m^2',                             &
                       missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xland', &
                       Grd%tracer_axes(1:3), Time%model_time,                              &
                       'xlandmix*rho*dzt*tracer for '//trim(T_prog(n)%name),               &
                       trim(T_prog(n)%flux_units),                                         &
                       missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
  enddo

  id_xland_mass  = register_diag_field ('ocean_model', 'mass_xland', &
                   Grd%tracer_axes(1:2), Time%model_time,            &
                   'mass source from xland', 'kg/(m^2*sec)',         &
                   missing_value=missing_value, range=(/-1.e8,1.e8/))

  call watermass_diag_init(Time, Dens)


191 format(/' ===== from ocean_xland_init ====='/                                  &
       ' for crossland sea connection pair number',i3,/                            &
       ' mix  (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      &
       ' with (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      & 
       ' from level',i3,' to',i3,' [depths of ',f10.3,' to ',f10.3,'m]' ) 
192 format(/'=>Error: Problem with crossland mixing: improper k-levels requested.'/&
       'kxland(',i2,',2) =',i4, ' and kmt(',i4,',',i4,') = ',i4)
991 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds i grid location requested'/  &
       '      ixland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
992 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds j-row grid location requested'/ &
       '      jxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
993 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds k-level grid location requested'/ &
       '      kxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
994 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' improper k-values requested'/' kxland(',i3,',1)=',i5,                     &
       ' and is not less than or equal to kxland(',i3,',2)=',i5,/)

end subroutine ocean_xlandmix_init
! </SUBROUTINE> NAME="ocean_xlandmix_init"


!#######################################################################
! <SUBROUTINE NAME="xlandmix">
!
! <DESCRIPTION>
! Compute thickness and density weighted time tendency from 
! xlandmix. Units of the tendency are tracer*(kg/m^3)*meter/sec.
!
! Logic is a bit tricky in order to allow each (i,j,k) point
! to participate in an arbitrary number of xlandmix pairs.  
!
! </DESCRIPTION>
!
subroutine xlandmix (Time, Ext_mode, Dens, Thickness, T_prog)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
 
  real, dimension(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo,kxland_min:kxland_max,num_prog_tracers) :: &
                   tracerx    ! tracer concentration

  real, dimension(2) :: top_rho_dzt       
  real               :: temporary
  integer            :: i, j, k, nxl, lx, n
  integer            :: taum1, tau

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandmix_mod (xlandmix): module must be initialized')
  endif 

  if (nxland < 1) return
  if (.not. use_this_module) return

  taum1 = Time%taum1
  tau   = Time%tau

  xland_mass(:,:)      = 0.0
  xland_mass_pair(:,:) = 0.0
 
  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then

         ! for surface cell, wish to have only difference be due to eta_t differences
         ! to ensure that xlandmix at surface will readily equilize eta_t across straits. 
         ! if use rho_dzt_x = rho0*Thickness%dzt, then with zstar and pstar, there would
         ! be a far smaller difference in rho_dzt_x than for the following specification,
         ! which is based on the geopotential case.  

     do k=kxland_min,kxland_max
        rho_dzt_x(isc:iec,jsc:jec,k) = Thickness%rho_dzt(isc:iec,jsc:jec,k,taum1)
     enddo
     rho_dzt_x(isc:iec,jsc:jec,1) = rho0*(Grd%dzt(1)+Ext_mode%eta_t(isc:iec,jsc:jec,taum1))
     call mpp_update_domains (rho_dzt_x(:,:,:), Xland_domain%domain2d)
     tracerx = 0.0
     do n = 1, num_prog_tracers
        do k=kxland_min,kxland_max
            tracerx(isc:iec,jsc:jec,k,n)   = T_prog(n)%field(isc:iec,jsc:jec,k,taum1)
        enddo
     enddo
     call mpp_update_domains( tracerx(:,:,:,:), Xland_domain%domain2d)
  endif  

  do n=1,num_prog_tracers

     T_prog(n)%wrk1(:,:,:) = 0.0

     if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
         do nxl=1,nxland 

            if(at_least_one_in_comp_domain(nxl)) then  

                do lx=1,2
                   ixl = ixland(nxl,lx)-Dom%ioff; jxl = jxland(nxl,lx)-Dom%joff
                   gamma_xland(nxl,lx) = (2.0*vxland(nxl))/  &
                        (datx(ixl,jxl)*(depth_static(nxl,1)+depth_static(nxl,2)))
                enddo

                ! tracer exchange for k>1 
                do lx=1,2 
                   if(on_comp_domain(nxl,lx)) then
                       ixl  = ixland(nxl,lx)-Dom%ioff
                       jxl  = jxland(nxl,lx)-Dom%joff
                       ixl2 = ixland(nxl,3-lx)-Dom%ioff
                       jxl2 = jxland(nxl,3-lx)-Dom%joff
                       do k=max(2,kxland(nxl,1)),kxland(nxl,2)  
                          T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k) + gamma_xland(nxl,lx) &
                               *onehalf*(rho_dzt_x(ixl2,jxl2,k)+rho_dzt_x(ixl,jxl,k))                 &
                               *(tracerx(ixl2,jxl2,k,n)-tracerx(ixl,jxl,k,n))
                       enddo
                   endif
                enddo

                ! mass and tracer exchange for k=1
                if(kxland(nxl,1)==1) then 

                    k=1
                    do lx=1,2
                       ixl = ixland(nxl,lx)-Dom%ioff
                       jxl = jxland(nxl,lx)-Dom%joff                 
                       top_rho_dzt(lx) = rho_dzt_x(ixl,jxl,k)
                    enddo

                    if(n==1) then 
                        do lx=1,2
                           ixl = ixland(nxl,lx)-Dom%ioff
                           jxl = jxland(nxl,lx)-Dom%joff                   
                           xland_mass_pair(nxl,lx) = gamma_xland(nxl,lx)*( top_rho_dzt(3-lx)-top_rho_dzt(lx) )
                           xland_mass(ixl,jxl)     = xland_mass(ixl,jxl) + xland_mass_pair(nxl,lx) 
                        enddo
                    endif

                    do lx=1,2  
                       if(on_comp_domain(nxl,lx)) then
                           ixl  = ixland(nxl,lx)-Dom%ioff
                           jxl  = jxland(nxl,lx)-Dom%joff
                           ixl2 = ixland(nxl,3-lx)-Dom%ioff
                           jxl2 = jxland(nxl,3-lx)-Dom%joff
                           T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k)        &
                                +gamma_xland(nxl,lx)*(top_rho_dzt(3-lx)*tracerx(ixl2,jxl2,k,n) &
                                -top_rho_dzt(lx)*tracerx(ixl,jxl,k,n)) 
                       endif
                    enddo

                endif

            endif ! end of at_least_one_in_comp_domain(nxl) if-test 

         enddo  ! end of the nxland do-loop 

     endif ! end of xland_halo > halo if-test 

     if(xland_halo <= Dom%xhalo .and. xland_halo<= Dom%yhalo) then 

         do nxl=1,nxland 

            if(at_least_one_in_comp_domain(nxl)) then  

                do lx=1,2
                   ixl = ixland(nxl,lx)-Dom%ioff
                   jxl = jxland(nxl,lx)-Dom%joff 
                   gamma_xland(nxl,lx) = (2.0*vxland(nxl))/  &
                        (Grd%dat(ixl,jxl)*(depth_static(nxl,1)+depth_static(nxl,2)))
                enddo

                ! tracer exchange for k>1 
                do lx=1,2 
                   if(on_comp_domain(nxl,lx)) then
                       ixl = ixland(nxl,lx)-Dom%ioff
                       jxl = jxland(nxl,lx)-Dom%joff
                       ixl2 = ixland(nxl,3-lx)-Dom%ioff
                       jxl2 = jxland(nxl,3-lx)-Dom%joff
                       do k=max(2,kxland(nxl,1)),kxland(nxl,2)  
                          T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k) + gamma_xland(nxl,lx)       &
                               *onehalf*(Thickness%rho_dzt(ixl2,jxl2,k,taum1)+Thickness%rho_dzt(ixl,jxl,k,taum1)) &
                               *(T_prog(n)%field(ixl2,jxl2,k,taum1)-T_prog(n)%field(ixl,jxl,k,taum1))
                       enddo
                   endif
                enddo

                ! mass and tracer exchange for k=1
                if(kxland(nxl,1)==1) then 

                    k=1
                    do lx=1,2
                       ixl = ixland(nxl,lx)-Dom%ioff
                       jxl = jxland(nxl,lx)-Dom%joff
                       top_rho_dzt(lx) = rho0*(Grd%dzt(1)+Ext_mode%eta_t(ixl,jxl,taum1))
                    enddo

                    if(n==1) then  
                        do lx=1,2
                           ixl = ixland(nxl,lx)-Dom%ioff
                           jxl = jxland(nxl,lx)-Dom%joff                   
                           xland_mass_pair(nxl,lx) = gamma_xland(nxl,lx)*( top_rho_dzt(3-lx)-top_rho_dzt(lx) )
                           xland_mass(ixl,jxl)     = xland_mass(ixl,jxl) + xland_mass_pair(nxl,lx)
                        enddo
                    endif

                    do lx=1,2  
                       if(on_comp_domain(nxl,lx)) then
                           ixl  = ixland(nxl,lx)-Dom%ioff
                           jxl  = jxland(nxl,lx)-Dom%joff
                           ixl2 = ixland(nxl,3-lx)-Dom%ioff
                           jxl2 = jxland(nxl,3-lx)-Dom%joff
                           T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k)                      &
                                +gamma_xland(nxl,lx)*(top_rho_dzt(3-lx)*T_prog(n)%field(ixl2,jxl2,k,taum1) &
                                -top_rho_dzt(lx)  *T_prog(n)%field(ixl,jxl,k,taum1))  
                       endif
                    enddo

                endif

            endif ! end of at_least_one_in_comp_domain(nxl) if-test 

         enddo  ! end of the nxland do-loop 

     endif ! end of xland_halo > halo if-test 


     ! fill source array 
     T_prog(n)%th_tendency(isc:iec,jsc:jec,:) = T_prog(n)%th_tendency(isc:iec,jsc:jec,:) &
          + T_prog(n)%wrk1(isc:iec,jsc:jec,:)

     ! send diagnostics 
     if(id_xland(n) > 0) then 
        call diagnose_3d(Time, Grd, id_xland(n), T_prog(n)%conversion*T_prog(n)%wrk1(:,:,:))
     endif

  enddo ! end of do loop n=1,num_prog_tracers 


  ! compute effects on mass 
  call xlandmix_mass(Time, Thickness)

  ! call some diagnostics 
  call watermass_diag(Time, Dens, &
  T_prog(index_temp)%wrk1(:,:,:), T_prog(index_salt)%wrk1(:,:,:))
  

end subroutine xlandmix
! </SUBROUTINE> NAME="xlandmix"


!#######################################################################
! <SUBROUTINE NAME="xlandmix_mass">
!
! <DESCRIPTION>
! Compute the mass source kg/(m^2*sec).
! Note that xlandmix has already been called, so xland_mass has been 
! filled.  
! </DESCRIPTION>
!
subroutine xlandmix_mass (Time, Thickness)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(inout) :: Thickness 
 
  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandmix_mod (xlandmix_mass): module must be initialized')
  endif 

  Thickness%mass_source(isc:iec,jsc:jec,1) = Thickness%mass_source(isc:iec,jsc:jec,1) &
                                           + xland_mass(isc:iec,jsc:jec)

  call mpp_update_domains (Thickness%mass_source(:,:,1), Dom%domain2d)

  call diagnose_2d(Time, Grd, id_xland_mass, xland_mass(:,:))

end subroutine xlandmix_mass
! </SUBROUTINE> NAME="xlandmix_mass"


!#######################################################################
! <SUBROUTINE NAME="xland_check">
!
! <DESCRIPTION>
! Check if prescribed too much mixing
!
! In this routine the crossland mixing rate vxland(nxl) for
! the pair of points associated with index number nxl is
! converted into the fraction of the model grid boxes to be mixed
! per second, and checked.  These checks ensure that the rate of
! crossland mixing requested is valid in that it can be realized
! given the timestep length and column volumes involved.
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
! <INOUT NAME="error" TYPE="logical">
! Error flag indicates whether initialization was performed successfully.  
! </INOUT>
!
subroutine xland_check(nxl, dtime, error)

  logical, intent(inout) :: error
  integer, intent(in)    :: nxl
  real,    intent(in)    :: dtime

  integer                   :: k, lx, ll
  real, dimension(2)        :: colvol
  real                      :: factor, tslim
  real, dimension(nxland,2) :: bxland 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandmix_mod (xland_check): module must be initialized')
  endif 

  ! Calculate static depth range over which crossland mixing will be applied. 
  ! Note that do not include partial bottom cell in kxland range.  
  do lx=1,2
    depth_static(nxl,lx) = 0.0
    do k=max(1,kxland(nxl,1)),kxland(nxl,2)
       ixl = ixland(nxl,lx)-Dom%ioff
       jxl = jxland(nxl,lx)-Dom%joff                            
       depth_static(nxl,lx) = depth_static(nxl,lx) + Grd%dzt(k)
    enddo
  enddo

  do lx=1,2

! For each of the two (i,j) points to be mixed, convert from volume
! to be mixed per second vxland(nxl) to model box fraction to be
! mixed per second bxland(nxl,l)
     ixl = ixland(nxl,lx)-Dom%ioff
     jxl = jxland(nxl,lx)-Dom%joff                          
     colvol(lx)     = datx(ixl,jxl)*depth_static(nxl,lx)
     bxland(nxl,lx) = vxland(nxl) / (epsln + colvol(lx))

    ! check for attempts to mix at a rate so fast that it can not
    ! be stably achieved in a single forward timestep
    ! (for the two boxes to be thoroughly mixed, no more than
    ! one-half of the volume of a given grid box can be
    ! stably transported in effect into the other gridbox)
    factor = bxland(nxl,lx)*dtime
    if (factor > 0.5) then
      tslim = 0.5/bxland(nxl,lx)
      write (unit,997) nxl, factor, ixland(nxl,lx), jxland(nxl,lx), bxland(nxl,lx), tslim, &
             vxland(nxl), kxland(nxl,1), kxland(nxl,2), colvol(lx)
      error=.true.
      return
    endif

  enddo

  if(verbose_init) then 
    write (unit,197) nxl, ixland(nxl,1), jxland(nxl,1), ixland(nxl,2), jxland(nxl,2), kxland(nxl,1), &
         kxland(nxl,2), depth_static(nxl,1), (colvol(ll),ll=1,2), vxland(nxl), (bxland(nxl,ll),ll=1,2)
  endif 

197   format(/' ===== from xlandvchk ====='/                             &
     ' for crossland sea communication pair number',i3,/                 &
     ' mix I,J gridpoints (',i4,',',i4,') and (',i4,',',i4,')'/          &
     ' from level',i3,' to',i3,' (a depth range of ',e12.6,' m)'/        &
     ' column volumes =',e12.6,' and ',e12.6,' m^3'/                     & 
     ' simulated flow in = flow out = ',e12.6,' m^3/sec,',/,             &
     ' so mix ',e12.6,' fraction of 1st column with ',e12.6,             &
     ' of 2nd column per sec.'/)
997   format(/' Error => problem with crossland tracer mixing',          &
     ' pair #',i3,/,' Error => attempting to mix at too fast a rate:'/   &
     ' asking for more than half (',g12.6,') the volume of a grid box',  &
     ' at i, j, k: ',2i5,/,' to be mixed into its neighbor on',          &
     ' leapfrog timesteps'/                                              &  
     ' As specified now, fraction of column mixed per sec = ',e12.6,/,   &
     ' Potential Fixes: shorten dtime to be less than',e13.6,/           &
     t18,' reduce mixing rate vxland (now ',e12.6,' m^3/sec)'/           &
     t18,' or mix over a larger depth range (now from k=',i3,' to ',i3   &
     ,' => water column volume = ',e12.6,' m^3)'/)

end subroutine xland_check
! </SUBROUTINE> NAME="xland_check"


!#######################################################################
! <FUNCTION NAME="at_least_one_in_comp_domain">
!
! <DESCRIPTION>
! Function to see if at least one of the two points in a crossland pair
! is within the computational domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
!
function at_least_one_in_comp_domain(nxl)

  integer, intent(in) :: nxl
  integer             :: lx
  logical             :: in_domain(2), at_least_one_in_comp_domain

  do lx=1,2
    if(isc+Dom%ioff <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff .and. &
       jsc+Dom%joff <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jec+Dom%joff) then
       in_domain(lx) = .true.
    else 
       in_domain(lx) = .false.  
    endif
  enddo

  at_least_one_in_comp_domain = .false.
  if(in_domain(1) .and. in_domain(2)     ) at_least_one_in_comp_domain = .true.
  if(in_domain(1) .and. .not.in_domain(2)) at_least_one_in_comp_domain = .true.
  if(.not.in_domain(1) .and. in_domain(2)) at_least_one_in_comp_domain = .true.


end function at_least_one_in_comp_domain 
! </FUNCTION> NAME="at_least_one_in_comp_domain"


!#######################################################################
! <FUNCTION NAME="both_in_local_domain">
!
! <DESCRIPTION>
! Determine if both points in a crossland pair are within the
! local domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
!
function both_in_local_domain(nxl)

  integer, intent(in) :: nxl
  integer             :: lx, ilo, ihi, jlo, jhi
  logical             :: in_domain(2), both_in_local_domain

  ilo = isd+Dom%ioff
  ihi = ied+Dom%ioff
  jlo = jsd+Dom%joff
  jhi = jed+Dom%joff

  do lx=1,2
    if(ilo <= ixland(nxl,lx) .and. ixland(nxl,lx) <= ihi .and. &
       jlo <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jhi) then
       in_domain(lx) = .true.
    else 
       in_domain(lx) = .false.  
    endif
  enddo
  if(in_domain(1) .and. in_domain(2)) then
    both_in_local_domain = .true.
  else 
    both_in_local_domain = .false.
  endif

end function both_in_local_domain 
! </FUNCTION> NAME="both_in_local_domain"


!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
! <IN NAME="lx" TYPE="integer">
! lx=1,2 labels the point within an xlandmix pair
! </IN>
!
function on_comp_domain(nxl, lx)

  integer, intent(in) :: nxl, lx
  logical             :: on_comp_domain

  if(isc+Dom%ioff <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"


!#######################################################################
! <FUNCTION NAME="on_processor">
!
! <DESCRIPTION>
! Determine if the point is on processor 
! </DESCRIPTION>
!
function on_processor(nxl, lx)

  integer, intent(in) :: nxl, lx
  logical             :: on_processor

  if(isc+Dom%ioff-xland_halo <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff+xland_halo .and. &
     jsd+Dom%joff-xland_halo <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jed+Dom%joff+xland_halo) then
     on_processor = .true.
  else 
     on_processor = .false.  
  endif

end function on_processor
! </FUNCTION> NAME="on_processor"


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

  id_neut_rho_xmix = register_diag_field ('ocean_model', 'neut_rho_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                               &
    'update of locally ref potrho from xlandmix',                        &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_xmix > 0) compute_watermass_diag = .true. 

  id_neut_rho_xmix_on_nrho = register_diag_field ('ocean_model',                 &
   'neut_rho_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,          &
   'update of locally ref potrho from xlandmix binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_xmix = register_diag_field ('ocean_model', 'wdian_rho_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                                 &
    'dianeutral mass transport due to xlandmix',                           &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_xmix > 0) compute_watermass_diag = .true. 

  id_wdian_rho_xmix_on_nrho = register_diag_field ('ocean_model',             &
   'wdian_rho_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,      &
   'dianeutral mass transport from xlandmix binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_xmix = register_diag_field ('ocean_model',    &
     'tform_rho_xmix', Grd%tracer_axes(1:3), Time%model_time,&
     'watermass transform from xlandmix',                    &
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_xmix_on_nrho = register_diag_field ('ocean_model',          &
     'tform_rho_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time, &
     'watermass transform from xlandmix binned to neutral density layers', &
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_eta_tend_xmix= -1          
  id_eta_tend_xmix= register_diag_field ('ocean_model','eta_tend_xmix',    &
       Grd%tracer_axes(1:2), Time%model_time,                              &
       'non-Bouss steric sea level tendency from xlandmix tendency', 'm/s',&
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_xmix > 0) compute_watermass_diag=.true.

  id_eta_tend_xmix_glob= -1          
  id_eta_tend_xmix_glob= register_diag_field ('ocean_model', 'eta_tend_xmix_glob',&
       Time%model_time,                                                           &
       'global mean non-bouss steric sea level tendency from xlandmix tendency',  &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_xmix_glob > 0) compute_watermass_diag=.true.


  ! temperature effects 
  id_neut_temp_xmix = register_diag_field ('ocean_model', 'neut_temp_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                                 &
    'temp related update of locally ref potrho from xlandmix',             &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_xmix > 0) compute_watermass_diag = .true. 

  id_neut_temp_xmix_on_nrho = register_diag_field ('ocean_model',                             &
   'neut_temp_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                      &
   'temp related update of locally ref potrho from xlandmix binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_xmix = register_diag_field ('ocean_model', 'wdian_temp_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                                   &
    'temp related dianeutral mass transport due to xlandmix',                &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_xmix > 0) compute_watermass_diag = .true. 

  id_wdian_temp_xmix_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_temp_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'temp related dianeutral mass transport from xlandmix binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_xmix = register_diag_field ('ocean_model',    &
     'tform_temp_xmix', Grd%tracer_axes(1:3), Time%model_time,&
     'temp related watermass transform from xlandmix',        &
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_xmix_on_nrho = register_diag_field ('ocean_model',                     &
     'tform_temp_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
     'temp related watermass transform from xlandmix binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_xmix_on_nrho > 0) compute_watermass_diag = .true. 


  ! salinity effects 
  id_neut_salt_xmix = register_diag_field ('ocean_model', 'neut_salt_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                                 &
    'salt related update of locally ref potrho from xlandmix',             &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_xmix > 0) compute_watermass_diag = .true. 

  id_neut_salt_xmix_on_nrho = register_diag_field ('ocean_model',                             &
   'neut_salt_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                      &
   'salt related update of locally ref potrho from xlandmix binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_xmix = register_diag_field ('ocean_model', 'wdian_salt_xmix',&
    Grd%tracer_axes(1:3), Time%model_time,                                   &
    'salt related dianeutral mass transport due to xlandmix',                &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_xmix > 0) compute_watermass_diag = .true. 

  id_wdian_salt_xmix_on_nrho = register_diag_field ('ocean_model',                         &
   'wdian_salt_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                  &
   'salt related dianeutral mass transport from xlandmix binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_xmix = register_diag_field ('ocean_model',    &
     'tform_salt_xmix', Grd%tracer_axes(1:3), Time%model_time,&
     'salt related watermass transform from xlandmix',        &
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_xmix_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_xmix_on_nrho = register_diag_field ('ocean_model',                     &
     'tform_salt_xmix_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,            &
     'salt related watermass transform from xlandmix binned to neutral density layers',&
     'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_xmix_on_nrho > 0) compute_watermass_diag = .true. 


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_xlandmix_mod w/ compute_watermass_diag=.true.'  
  endif 


end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from xlandmix on watermass transformation.  
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens, theta_tend, salt_tend)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_density_type),       intent(in) :: Dens
  real, dimension(isd:,jsd:,:),   intent(in) :: theta_tend
  real, dimension(isd:,jsd:,:),   intent(in) :: salt_tend

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend
  real    :: eta_tend_glob

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_xlandmix (watermass_diag): module needs initialization ')
  endif 

  if (.not. compute_watermass_diag) return

  tau = Time%tau

  ! effects from both temperature and salinity 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*theta_tend(i,j,k) &
                        +Dens%drhodS(i,j,k)*salt_tend(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2)  ! for eta_tend 
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_xmix, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_xmix, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_xmix, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_xmix_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_xmix_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_xmix_on_nrho, wrk4)

  if(id_eta_tend_xmix > 0 .or. id_eta_tend_xmix_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_xmix, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_xmix_glob, eta_tend, cellarea_r)
  endif

 
  ! effects from temperature
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodT(i,j,k)*theta_tend(i,j,k) 
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_temp_xmix, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_xmix, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_xmix, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_xmix_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_xmix_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_xmix_on_nrho, wrk4)
 
  ! effects from salinity 
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Dens%drhodS(i,j,k)*salt_tend(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_salt_xmix, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_xmix, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_xmix, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_xmix_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_xmix_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_xmix_on_nrho, wrk4)

 
end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_xlandmix_mod

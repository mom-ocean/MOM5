module ocean_xlandinsert_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> K.W. Dixon
!</CONTACT>  
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison
!</REVIEWER>   
!
!<OVERVIEW>
! Tracer and mass source from cross-land insertion.
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness and density weighted tendency [rho*tracer*meter/sec]
! of tracer associated with cross-land insertion.  Also compute 
! mass source from cross-land insertion.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <NOTE>
! Algorithm ensures both total tracer and total mass is conserved.
! Algorithm sets up the insertion points as in xlandmix, and 
! transports between the columns as in rivermix. 
! </NOTE>
!
! <NOTE>
! 2D domain decomposition is implemented according to following
! notions.  If the two points in an xlandinsert pair live within halo 
! of each other (i.e., within same local_domain), 
! then no added mpp communication required. However, nore generally 
! the two points live further away.  In this case, xland_domain
! is defined so that its halos incorporate the maximum separation 
! of xlandinsert points.  New tracer, eta, and grid arrays
! are defined over this extended xland_domain.  This added domain
! size will come at some computational cost, so it is advantageous
! to limit the separation between points within an xlandinsert pair. 
! </NOTE>
!
! <NOTE>
! The current implementation of xlandinsert has not been generalized 
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
! BUG FIX 30 May 2011. 
! Trying to mix across the ends of the global domain results in incorrect indexing.
! Local indices need to be incremented/decremented by the global extent.
! Russell.Fiedler@csiro.au
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_xlandinsert_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Needs to be true in order to use this scheme. Default is false.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose initialization information.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.
!  </DATA> 
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field
use field_manager_mod, only: MODEL_OCEAN, parse, find_field_index
use field_manager_mod, only: get_field_methods, method_type, get_field_info
use fms_mod,           only: stdout, stdlog, FATAL, NOTE, WARNING
use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
use mpp_domains_mod,   only: mpp_update_domains 
use mpp_domains_mod,   only: BITWISE_EXACT_SUM
use mpp_mod,           only: input_nml_file, mpp_error

use ocean_domains_mod,    only: get_local_indices, get_global_indices, set_ocean_domain
use ocean_parameters_mod, only: rho0, rho0r, missing_value
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_external_mode_type 
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_options_type, ocean_density_type 
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk5, wrk1_v
use ocean_util_mod,       only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_tracer_util_mod,only: diagnose_3d_rho

implicit none

private 

public  ocean_xlandinsert_init
public  xlandinsert
private xlandinsert_mass
private watermass_diag_init
private watermass_diag

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()


! Variables set inside subroutine ocean_xlandinsert_init
! See comments in ocean_xlandinsert_init for description
integer, dimension(:,:), allocatable :: ixland     
integer, dimension(:,:), allocatable :: jxland     
integer, dimension(:,:), allocatable :: kxland     
real, dimension(:), allocatable      :: tauxland     
real, dimension(:), allocatable      :: tauxlandr     
real, dimension(:,:), allocatable    :: water_rate ! rate (kg/m^3)*(m/s) of water transferred between columns

! number of xland mixing pairs 
integer :: nxland=0 

integer                              :: xland_halo     ! halo size needed to perform xlandinsert
integer                              :: kxland_min     ! min of the k-levels for the nxland pairs
integer                              :: kxland_max     ! max of the k-levels for the nxland pairs
type(ocean_domain_type), save        :: Xland_domain   ! domain for xlandinsert

! variables defined with halos given by xland_halo
real, dimension(:,:), allocatable    :: xland_rho_dzt  ! mass source on (x,y) space from xlandinsert (kg/m^3)*(m/s)
real, dimension(:,:), allocatable    :: rho_dzt_x      ! surface value for rho_dzt (kg/m^3)*(m)
real, dimension(:,:), allocatable    :: xtx            ! zonal position of tracer points (degrees)
real, dimension(:,:), allocatable    :: ytx            ! meridional position of tracer points (degrees)
real, dimension(:,:), allocatable    :: datx           ! area of tracer cell (m^2)
real                                 :: dtime          ! time step

integer                              :: unit=6         ! processor zero writes to unit 6
  
logical, dimension(:), allocatable   :: error_xland    ! for checking that all xland points are OK.  

! work array for eta_tend 
real, dimension(:,:,:),   allocatable :: theta_tend
real, dimension(:,:,:),   allocatable :: salt_tend

! for global area normalization
real    :: cellarea_r

! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xland
integer :: id_xland_rho_dzt=-1

integer :: id_neut_rho_xinsert          =-1
integer :: id_neut_rho_xinsert_on_nrho  =-1
integer :: id_wdian_rho_xinsert         =-1
integer :: id_wdian_rho_xinsert_on_nrho =-1
integer :: id_tform_rho_xinsert         =-1
integer :: id_tform_rho_xinsert_on_nrho =-1

integer :: id_eta_tend_xinsert     =-1
integer :: id_eta_tend_xinsert_glob=-1

integer :: id_neut_temp_xinsert          =-1
integer :: id_neut_temp_xinsert_on_nrho  =-1
integer :: id_wdian_temp_xinsert         =-1
integer :: id_wdian_temp_xinsert_on_nrho =-1
integer :: id_tform_temp_xinsert         =-1
integer :: id_tform_temp_xinsert_on_nrho =-1

integer :: id_neut_salt_xinsert          =-1
integer :: id_neut_salt_xinsert_on_nrho  =-1
integer :: id_wdian_salt_xinsert         =-1
integer :: id_wdian_salt_xinsert_on_nrho =-1
integer :: id_tform_salt_xinsert         =-1
integer :: id_tform_salt_xinsert_on_nrho =-1

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1

character(len=128) :: version=&
       '$Id: ocean_xlandinsert.F90,v 20.0 2013/12/14 00:16:28 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: isg, ieg, jsg, jeg
integer :: ixl, jxl, ixl2, jxl2

logical :: verbose_init          = .true.
logical :: use_this_module       = .false.
logical :: module_is_initialized = .FALSE.
logical :: debug_this_module     = .false.

namelist /ocean_xlandinsert_nml/ verbose_init, use_this_module, debug_this_module

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_xlandinsert_init">
!
! <DESCRIPTION>
!    Initial set up for crossland insertion of tracers and eta 
!
!    Checks are performed to ensure that the crossland mixing
!    grid locations are valid according to model configuration.
!
!    A summary of the locations of crossland points is written out.
!
!    User specified inputs in "USER INPUT" section:
!
!    ixland and jxland = user specified nxland pairs of i,j grid locations
!
!    kxland = user specified uppermost (ktop=kxland(n,1)) and  
!             deepest (kbot=kxland(n,2)) levels over which crossland 
!             will be done for each crossland pair. Note the for 
!             xlandinsert, kxland(nxl,1)=1 is required, since the 
!             aim is to move excess volume from one column to another. 
!
!    tauxland = user specified time constant (seconds) setting the rate 
!               of transport via upwind advection.  
!
!    NOTE: for ixland, jxland, and kxland pairs, values of the
!    (n,1) element should be < the corresponding (n,2) value.
!
! </DESCRIPTION>

subroutine ocean_xlandinsert_init(Grid, Domain, Time, Dens, T_prog, Ocean_options, d_time)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_density_type),     intent(in)         :: Dens
  type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
  type(ocean_options_type),     intent(inout)      :: Ocean_options
  real,                         intent(in)         :: d_time

  type(method_type), allocatable, dimension(:) :: xland_methods
  character(len=32) :: fld_type, fld_name
  integer :: n, nt, model, imax, imin
  integer :: nxl, lx, xland_halo_test,i
  integer :: parse_ok
  integer :: io_status, ioun, ierr
  real    :: ztop 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if (module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandinsert_mod (ocean_xlandinsert_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  call get_global_indices(Domain,isg,ieg,jsg,jeg)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid
  
  dtime = d_time 
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)
  
  num_prog_tracers = size(T_prog(:))
  do nt=1,num_prog_tracers
     if (T_prog(nt)%name == 'temp') index_temp = nt
     if (T_prog(nt)%name == 'salt') index_salt = nt
  enddo

  n = find_field_index(MODEL_OCEAN,'xland_insert')

  if (n < 1) then 
    write(stdoutunit,*)' '
    write(stdoutunit,*) &
    '==>Warning: ocean_xlandinsert_init found n < 1 for xland_insert table.  Will NOT use ocean_xlandinsert.'  
    write(stdoutunit,*)' '
    Ocean_options%xlandinsert = 'Did NOT use xlandinsert.'
    return
  endif 

  call get_field_info(n,fld_type,fld_name,model,nxland)

  ! provide for namelist over-ride of default values
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_xlandinsert_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_xlandinsert_nml')
#else
  ioun =  open_namelist_file()
  read (ioun,ocean_xlandinsert_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status,'ocean_xlandinsert_nml')
  call close_file (ioun)
#endif
  write (stdlogunit,ocean_xlandinsert_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_xlandinsert_nml)

  if(use_this_module) then 
      if(nxland == 0) then 
          write(stdoutunit, *) &
          '==>Warning: use_this_module=.true. but no mixing points chosen. Is this what you wish?'
          call mpp_error(WARNING,&
          '==>NOTE from ocean_xlandinsert_mod: No cross-land mixing points chosen')
          Ocean_options%xlandinsert = 'Used xlandinsert, but set zero cross-land points, so module does nothing.'
          return
      else
          call mpp_error(NOTE, &
          '==>Note from ocean_xlandinsert_mod (ocean_xlandinsert_init): USING this module')
          write(stdoutunit, *) &
          '==>NOte: Using xlandinsert to connect tracers and mass between non-local ocean cells.'
          Ocean_options%xlandinsert = 'Used xlandinsert to transfer tracer between inland seas and open ocean.'
      endif
  else 
      call mpp_error(NOTE, '==>Note from ocean_xlandinsert_mod: NOT USING this module')
      if(nxland > 0) then 
          call mpp_error(WARNING, &
          '==>Warning: nxland > 0, but use_this_module=.false. NOT using xlandinsert. Is this what you wish?')
      endif
      Ocean_options%xlandinsert = 'Did NOT use xlandinsert.'
      return
  endif

  write(stdoutunit,'(/a,f10.2)')'==> Note from ocean_xlandinsert_mod: using time step of (secs)', dtime 

  allocate(xland_methods(nxland))

  allocate (ixland(nxland,2))
  ixland(:,:) = 0
  allocate (jxland(nxland,2))
  jxland(:,:) = 0
  allocate (kxland(nxland,2))
  kxland(:,:) = 0
  allocate (tauxland(nxland))
  tauxland(:)   = 0.0
  allocate (tauxlandr(nxland))
  tauxlandr(:)   = 0.0
  allocate (error_xland(nxland))
  error_xland(:) = .false.
  allocate (water_rate(nxland,2))
  water_rate(:,:) = 0.0

  call get_field_methods(n,xland_methods)
  do i=1, nxland
     parse_ok = parse(xland_methods(i)%method_control,'ixland_1',ixland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "ixland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'ixland_2',ixland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "ixland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_1',jxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "jxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_2',jxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "jxland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_1',kxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "kxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_2',kxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "kxland_2" error')
     parse_ok   = parse(xland_methods(i)%method_control,'tauxland',tauxland(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandinsert_mod: xland table entry "tauxland" error')
  enddo

  if(Grid%tripolar) then 
     write(stdoutunit,'(a)') &
     '==>Warning: xlandinsert has not been implemented to connect points across the tripolar fold.'
  endif

  write(stdoutunit, *) &
  'Using xlandinsert to connect tracers and surface height between non-local ocean cells.'


  xland_halo = min(isc-isd,ied-iec,jsc-jsd,jed-jec)
  kxland_min = nk
  kxland_max = 1

  do nxl=1,nxland

     ! check for invalid crossland insertion grid locations

     if(kxland(nxl,1) > 1) then 
       call mpp_error(FATAL,&
       '==>Error in ocean_xlandinsert_init: kxland set wrong--must start from top cell w/ kxland(nxl,1)=1.')
     endif 

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

     ! set tauxlandr as inverse of tauxland, except when tauxland=0.
     if(tauxland(nxl) > 0.0) tauxlandr(nxl) = 1.0/tauxland(nxl)

  enddo  ! end of do-loop nxl=1,nxland

  if(xland_halo == 0) then 
     call mpp_error(FATAL,&
     '==>Error in ocean_xlandinsert_init: xland_halo=0 => problems defining xlandinsert points.')
  endif
 
  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
    write(stdoutunit,*) &
    'Defining extra tracer arrays for xland_domain over k-levels ',kxland_min,' through ',kxland_max 
  endif 

! define arrays and xland_domain
  
  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
      write(stdoutunit,*)'The model local computational domain has a x,y halo = ',Dom%xhalo, Dom%yhalo
      write(stdoutunit,*)'This is smaller than the halo required for xlandinsert.'
      write(stdoutunit,*)'For xlandinsert, define a new domain type with halo = ', xland_halo

      if(xland_halo > (Dom%xhalo+6)) then 
          write(stdoutunit,*)'=>Warning: ocean_xlandinsert_init determined xland_halo = ', xland_halo
          write(stdoutunit,*)'           Are you sure you wish to connect such distant locations?'
          write(stdoutunit,*)'           Are you instead trying to connect points across a tripolar fold?'
          write(stdoutunit,*)'           If so, then presently this has yet to be coded into xlandinsert.'
          write(stdoutunit,*)'           Alternative xlandinsert points must be taken, or new code implemented.'
      endif

  endif

  ! define xlandinsert domain 
  call set_ocean_domain(Xland_domain,Grd,xhalo=xland_halo,yhalo=xland_halo,&
                        name='xlandinsert',maskmap=Dom%maskmap)

  ! time independent arrays defined over the xland_domain
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

  ! time dependent arrays defined over the xland_domain
  allocate (xland_rho_dzt(isc:iec,jsc:jec))
  allocate (rho_dzt_x(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
  xland_rho_dzt(:,:) = 0.0
  rho_dzt_x(:,:)     = 0.0

! Mixing across Straits of Gibraltar in CSIRO MK3.5 model resulted in 
! exceeding bounds with the Dec2009 verison of MOM4p1. 
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
                     write(unit,*) 'xlandinsert init: ixland outside local domain i=',lbound(xtx,1),':',ubound(xtx,1),   &
                                   ' Changing index from ', ixland(nxl,lx) ,' to ', ixl2
                     ixland(nxl,lx) = ixl2
                  else if ( ixl > ubound(xtx,1) ) then             ! Gone outside upper index
                     ixl2 = ixland(nxl,lx) - ieg
                     write(unit,*) 'xlandinsert init: ixland outside local domain i=',lbound(xtx,1),':',ubound(xtx,1),   &
                                   ' Changing index from ', ixland(nxl,lx) ,' to ', ixl2
                     ixland(nxl,lx) = ixl2
                  endif
                  if(.not. on_processor(nxl, lx)) then
                    write(unit,*) 'xlandinsert init?: ixl =', ixl, ' Bounds ',isc+Dom%ioff-xland_halo,':',iec+Dom%ioff+xland_halo
                    error_xland(nxl)=.true.
                  endif
               endif
             enddo
         endif
      enddo
  endif

  do nxl=1,nxland

     ! check for attempts to connect land rather than sea
     do lx = 1,2
        if(on_comp_domain(nxl, lx)) then
            ixl = ixland(nxl,lx)-Dom%ioff
            jxl = jxland(nxl,lx)-Dom%joff
            if (kxland(nxl,2) > Grid%kmt(ixl,jxl)) then
                write (unit,192) nxl, kxland(nxl,2), ixland(nxl,lx),jxland(nxl,lx),Grid%kmt(ixl,jxl)
                error_xland(nxl)=.true.
            endif
        endif
     enddo

     ! write out summary information for this pair of crossland points
     if (kxland(nxl,1) == 1) then
        ztop = 0.0
     else
        ztop = Grid%zw(kxland(nxl,1)-1)
     endif

     if(at_least_one_in_comp_domain(nxl)) then
         ixl  = ixland(nxl,1)-Dom%ioff
         jxl  = jxland(nxl,1)-Dom%joff
         ixl2 = ixland(nxl,2)-Dom%ioff
         jxl2 = jxland(nxl,2)-Dom%joff
         if(verbose_init) then 
             write(unit,191) nxl, ixland(nxl,1), jxland(nxl,1) &
              ,xtx(ixl,jxl), ytx(ixl,jxl) &
              ,ixland(nxl,2), jxland(nxl,2) &
              ,xtx(ixl2,jxl2), ytx(ixl2,jxl2) &
              ,kxland(nxl,1), kxland(nxl,2), ztop, Grid%zw(kxland(nxl,2)), &
               tauxland(nxl) 
         endif 
     endif

  enddo

  ! Bring model down here if there are problems with xlandinsert points.  
  do nxl=1,nxland
    if(error_xland(nxl)) then
      call mpp_error(FATAL, &
      '==>Error detected in ocean_xlandinsert_mod: Use "grep -i error" to find ALL errors.')
    endif 
  enddo 

  ! register for diag_manager 
  allocate (id_xland(num_prog_tracers))
  id_xland = -1

  do n=1,num_prog_tracers
     if(T_prog(n)%name == 'temp') then 
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xlandinsert', &
                       Grd%tracer_axes(1:3), Time%model_time, 'xlandinsert*cp*rho*dzt*temp',     &
                      'Watt/m^2', missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xlandinsert', &
                       Grd%tracer_axes(1:3), Time%model_time,                                    &
                      'xlandinsert*rho*dzt*tracer for '//trim(T_prog(n)%name),                   &
                      trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
  enddo

  call watermass_diag_init(Time, Dens)

191 format(/' ===== from ocean_xlandinsert_init ====='/                            &
       ' for crossland sea connection pair number',i3,/                            &
       ' mix  (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      &
       ' with (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      & 
       ' from level',i3,' to',i3,' [depths of ',f10.3,' to ',f10.3,'m]', /         &
       ' Time scale for insertion via upwind advection is (sec) ',f12.6 /) 
192 format(/'=>Error: Problem with crossland insertion: improper k-levels requested.'/&
       'kxland(',i2,',2) =',i4, ' exceeds kmt(',i4,',',i4,') = ',i4)
991 format(/'=>Error: problem with crossland tracer insertion:',/,                 &
       ' out of bounds i grid location requested'/  &
       '      ixland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
992 format(/'=>Error: problem with crossland tracer insertion:',/,                 &
       ' out of bounds j-row grid location requested'/ &
       '      jxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
993 format(/'=>Error: problem with crossland tracer insertion:',/,                 &
       ' out of bounds k-level grid location requested'/ &
       '      kxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
994 format(/'=>Error: problem with crossland tracer insertion:',/,                 &
       ' improper k-values requested'/' kxland(',i3,',1)=',i5,                     &
       ' and is not less than or equal to kxland(',i3,',2)=',i5,/)

end subroutine ocean_xlandinsert_init
! </SUBROUTINE> NAME="ocean_xlandinsert_init"



!#######################################################################
! <SUBROUTINE NAME="xlandinsert">
!
! <DESCRIPTION>
! Compute thickness and density weighted tracer source 
! [(kg/m^3)*tracer*m/s] and mass source (kg/m^3)*(m/s)
! associated with discharge of tracer from surface of a 
! thick column into a set of boxes within a thin column. 
!
! NOTE:
! Compute rho_dzt_x as if using a geopotential model, since the  
! xlandinsert rates are typically first tuned with geopotential models.
! also, we wish to have the transport proportional to differences in 
! eta_t, which for zstar and pstar are generally must larger than 
! differences in dzt.  
!
! Modified by Keith.Dixon (January 2004) 
!
! Modified by Stephen.Griffies for general 
! vertical coordinates January 2005 and August 2006.
!
! </DESCRIPTION>
!
subroutine xlandinsert (Time, Ext_mode, Dens, Thickness, T_prog)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)

  ! tracer concentration on stack
  real, dimension(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo,kxland_min:kxland_max,num_prog_tracers) :: tracerx    

  real, dimension(nxland,0:2) :: top_rho_dzt   
  real, dimension(nxland,2)   :: area       
  real, dimension(nxland,2)   :: mass     
  real, dimension(nxland)     :: rho_dzt_new_large
  real, dimension(nxland)     :: rho_dzt_new_small

  real  :: thkocean    ! thickness of column over which we insert seawater 
  real  :: delta(nk)   ! fraction inserted into a cell 

  real  :: tracernew(nk)
  real  :: blanda, tracerextra, zextra, zinsert
  real  :: tottracer_thin, tottracer_thick 
  real  :: totrate_thin, totrate_thick 

  integer :: i, j, k, nxl, nz
  integer :: tau
  integer :: n
  integer :: lx
  integer :: lx_thin(max(1,nxland))
  integer :: lx_thick(max(1,nxland))

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error in ocean_xlandinsert_mod (xlandinsert): module must be initialized')
  endif

  if (nxland < 1) return
  if (.not. use_this_module) return

  tau   = Time%tau

! tracer independent computations

  xland_rho_dzt(:,:) = 0.0
  rho_dzt_x(:,:)     = 0.0
  tracerx(:,:,:,:)   = 0.0

  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
     rho_dzt_x(isc:iec,jsc:jec) = rho0*(Grd%dzt(1)+ext_mode%eta_t(isc:iec,jsc:jec,tau))
     call mpp_update_domains (rho_dzt_x(:,:), Xland_domain%domain2d)
     do n=1,num_prog_tracers
        do k=1,kxland_max
           do j=jsc,jec
              do i=isc,iec
                 tracerx(i,j,k,n) = T_prog(n)%field(i,j,k,tau)
              enddo
           enddo
        enddo
     enddo
     call mpp_update_domains (tracerx(:,:,:,:), Xland_domain%domain2d)
  else
     rho_dzt_x(:,:) =  rho0*(Grd%dzt(1)+ext_mode%eta_t(:,:,tau))
     do n=1,num_prog_tracers
        do k=1,kxland_max
           tracerx(:,:,k,n) = T_prog(n)%field(:,:,k,tau)
        enddo
     enddo
  endif

  do nxl=1,nxland 

     ! calculate only if at least one point in the pair is on the comp-domain  
     if(at_least_one_in_comp_domain(nxl)) then  

         ! geometric factors      
         do lx=1,2
            ixl                 = ixland(nxl,lx)-Dom%ioff
            jxl                 = jxland(nxl,lx)-Dom%joff
            top_rho_dzt(nxl,lx) = rho_dzt_x(ixl,jxl) 
            area(nxl,lx)        = datx(ixl,jxl) 
            mass(nxl,lx)        = top_rho_dzt(nxl,lx)*area(nxl,lx) 
         enddo
         if(top_rho_dzt(nxl,1) > top_rho_dzt(nxl,2)) then 
             lx_thick(nxl) = 1
             lx_thin(nxl)  = 2
         else 
             lx_thick(nxl) = 2
             lx_thin(nxl)  = 1
         endif

         ! transfer rate
         top_rho_dzt(nxl,0)            = (mass(nxl,1) + mass(nxl,2))/(area(nxl,1) + area(nxl,2))
         water_rate(nxl,lx_thick(nxl)) = (top_rho_dzt(nxl,lx_thick(nxl))-top_rho_dzt(nxl,0))*tauxlandr(nxl)
         water_rate(nxl,lx_thin(nxl))  = (top_rho_dzt(nxl,0)-top_rho_dzt(nxl,lx_thin(nxl))) *tauxlandr(nxl)

         if(debug_this_module) then 
           write(*,*)' ' 
           write(*,*)'For xlandinsert pair                    = ',nxl 
           write(*,*)'top_rho_dzt(',lx_thick(nxl),') (kg/m^2) = ',top_rho_dzt(nxl,lx_thick(nxl)),'   (nxl=',nxl,')'
           write(*,*)'top_rho_dzt(',lx_thin(nxl),') (kg/m^2)  = ',top_rho_dzt(nxl,lx_thin(nxl)),'   (nxl=',nxl,')'
           write(*,*)'top_rho_dzt( 0 ) (kg/m^2)               = ',top_rho_dzt(nxl,0),'   (nxl=',nxl,')'
           write(*,*)'water_rate_thin   (kg/m^3)*(m/s)        = ',water_rate(nxl,lx_thin(nxl)),'   (nxl=',nxl,')'
           write(*,*)'water_rate_thick  (kg/m^3)*(m/s)        = ',water_rate(nxl,lx_thick(nxl)),'   (nxl=',nxl,')'
           write(*,*)'area_thin         (m^2)          = ',area(nxl,lx_thin(nxl)),'   (nxl=',nxl,')'
           write(*,*)'area_thick        (m^2)          = ',area(nxl,lx_thick(nxl)),'   (nxl=',nxl,')'
           write(*,*) &
            '{[area*water_rate](thick) - [area*water_rate](thin)}/(area(thin)+area(thick)) (kg/m^3)*(m/s) = ', &
             ( water_rate(nxl,lx_thin(nxl))*area(nxl,lx_thin(nxl))    &
              -water_rate(nxl,lx_thick(nxl))*area(nxl,lx_thick(nxl))) &
               /(area(nxl,lx_thin(nxl)) + area(nxl,lx_thick(nxl))),'   (nxl=',nxl,')'
             write(*,*) &
             ' >> check mass rates (thin)  ', water_rate(nxl,lx_thin(nxl))*area(nxl,lx_thin(nxl)),' (nxl=',nxl,')'
             write(*,*) &
             ' >> check mass rates (thick) ', water_rate(nxl,lx_thick(nxl))*area(nxl,lx_thick(nxl)),' (nxl=',nxl,')'
             write(*,*)' Insert ',water_rate(nxl,lx_thin(nxl))*dtime*area(nxl,lx_thin(nxl)),'kg', &
                       ' of water over ',dtime,' sec timestep   (nxl=',nxl,')'

             if(on_comp_domain(nxl,lx_thick(nxl))) then
                write(*,*)' the thick point is on the computational domain (nxl=',nxl,')'
             else
                write(*,*)' the thick point is NOT on the computational domain (nxl=',nxl,')'
             endif
             if(on_comp_domain(nxl,lx_thin(nxl))) then
                write(*,*)' the thin  point is on the computational domain (nxl=',nxl,')'
             else
                write(*,*)' the thin  point is NOT on the computational domain (nxl=',nxl,')'
             endif

         endif

! compute tendency to remove water from region of more mass to 
! region with less mass.
!
! xland_rho_dzt is added to Thickness%mass_source
!
! rho_dzt_new_large is the estimated tau+1 value of rho_dzt if
! the xland_insert adjustment was the only process contributing
! to the tendency
!
! note that the more massive side can be associated with more that one xland pair,
! so the tendency is accumulated.  Also, this is only calculated if the massive 
! side is a computational point.
         
         lx=lx_thick(nxl) 
         if(on_comp_domain(nxl,lx)) then
             ixl = ixland(nxl,lx)-Dom%ioff
             jxl = jxland(nxl,lx)-Dom%joff
             xland_rho_dzt(ixl,jxl) = xland_rho_dzt(ixl,jxl) - water_rate(nxl,lx)
             rho_dzt_new_large(nxl) = rho_dzt_x(ixl,jxl) - dtime*water_rate(nxl,lx)

             if(debug_this_module) then
               write(*,*) ' BEFORE: (thick) rho_dzt=',rho_dzt_x(ixl,jxl),' (nxl=',nxl,')'
               write(*,*) ' BEFORE: (thick) sfc mass=',(area(nxl,lx)*top_rho_dzt(nxl,lx)),' (nxl=',nxl,')'
               write(*,*) ' >> check sfc mass above =', mass(nxl,lx),' (nxl=',nxl,')'
               write(*,*) ' PREDICT: (thick) mass shift per timestep=',(-water_rate(nxl,lx)*area(nxl,lx)*dtime), &
                            ' (nxl=',nxl,')'
               write(*,*) ' PREDICT: rho_dzt_new_large=', rho_dzt_new_large(nxl), &
                            ' (nxl=',nxl,')'
               write(*,*) ' PREDICT: (thick) sfc mass(tau+1)=',( area(nxl,lx)* &
                                     (top_rho_dzt(nxl,lx)-(water_rate(nxl,lx)*dtime)) ), ' (nxl=',nxl,')'
             endif

         endif

! compute tendency to add water to region of less mass from region of more mass.
!
! xland_rho_dzt is added to Thickness%mass_source
!
! rho_dzt_new_small is the estimated tau+1 value of rho_dzt if
! the xland_insert adjustment was the only process contributing
! to the tendency
!
! note that the less massive side can be associated with more that one xland pair,
! so the tendency is accumulated.  Also, this is only calculated if the less massive 
! side is a computational point.
         
         lx = lx_thin(nxl)
         if(on_comp_domain(nxl,lx)) then
             ixl = ixland(nxl,lx)-Dom%ioff
             jxl = jxland(nxl,lx)-Dom%joff
             xland_rho_dzt(ixl,jxl) = xland_rho_dzt(ixl,jxl) + water_rate(nxl,lx)
             rho_dzt_new_small(nxl) = rho_dzt_x(ixl,jxl) + dtime*water_rate(nxl,lx_thin(nxl))

             if(debug_this_module) then
               write(*,*) ' BEFORE: (thin)  rho_dzt=',rho_dzt_x(ixl,jxl),' (nxl=',nxl,')'
               write(*,*) ' BEFORE: (thin)  sfc mass=',(area(nxl,lx)*top_rho_dzt(nxl,lx)),' (nxl=',nxl,')'
               write(*,*) ' PREDICT: (thin)  mass shift per timestep= ',(water_rate(nxl,lx)*area(nxl,lx)*dtime), &
                          ' (nxl=',nxl,')'
               write(*,*) ' PREDICT: rho_dzt_new_small =', rho_dzt_new_small(nxl), &
                            ' (nxl=',nxl,')'
               write(*,*) ' PREDICT: (thin)  sfc mass(tau+1)=',(area(nxl,lx)* &
                                     (top_rho_dzt(nxl,lx)+(water_rate(nxl,lx)*dtime))), ' (nxl=',nxl,')'
             endif
         endif

     endif ! end of at_least_one_in_comp_domain(nxl) if-test 

  enddo  ! end of the nxland do-loop 


  ! do the tracers 
  do n=1,num_prog_tracers

     wrk1(:,:,:)    = 0.0

     ! loop over nxland pairs 
     do nxl=1,nxland 

        ! perform calculations only if at least one point in the pair is on comp-domain  
        if(at_least_one_in_comp_domain(nxl)) then  

            ! remove tracer from top of the thick column via upwind advection  
            lx=lx_thick(nxl) 
            k =1                                ! kxland(nxl,1)=1 for xlandinsert  
            blanda = 0.0
            if(on_comp_domain(nxl,lx)) then
                ixl = ixland(nxl,lx)-Dom%ioff
                jxl = jxland(nxl,lx)-Dom%joff
                blanda = -1.0*water_rate(nxl,lx)*tracerx(ixl,jxl,k,n)
                wrk1(ixl,jxl,k) = wrk1(ixl,jxl,k) + blanda

                if(debug_this_module) then
                  ! let totrate_thick keep track of the mass weighted tendency of the tracer content and
                  ! tottracer_thick keep track of the total tracer content
                  tottracer_thick = tracerx(ixl,jxl,k,n)*top_rho_dzt(nxl,lx)*area(nxl,lx)
                  totrate_thick   = blanda*area(nxl,lx)
                  write(*,*)' ' 
                  write(*,*) 'BEFORE: tracer value = ',tracerx(ixl,jxl,k,n), &
                             '(nxl=',nxl,' tracer#=',n,' k=',k,' thick)'
                  write(*,*) 'BEFORE: thick mass  = ',(top_rho_dzt(nxl,lx)*area(nxl,lx)),' = ',&
                              top_rho_dzt(nxl,lx),'*',area(nxl,lx), &
                             '(nxl=',nxl,' tracer#=',n,' k=',k,' thick)'
                  write(*,*) 'BEFORE: tottracer_thick = ',tottracer_thick, &
                             '(nxl=',nxl,' tracer#=',n,' k=',k,')'
                  write(*,*) 'xland_insert totrate_thick = ',totrate_thick, &
                             '(nxl=',nxl,' tracer#=',n,' k=',k,' thick)'
                  write(*,*) 'PREDICT AFTER: thick mass   = ', area(nxl,lx) * &
                              ( top_rho_dzt(nxl,lx) - (water_rate(nxl,lx)*dtime) ),' (nxl=',nxl,' tracer#=',n,')'
                  write(*,*) 'PREDICT AFTER: tottracer_thick = ',( tracerx(ixl,jxl,k,n) * area(nxl,lx) * &
                              ( top_rho_dzt(nxl,lx) - (water_rate(nxl,lx)*dtime) ) ),' (nxl=',nxl,' tracer#=',n,')'
                endif

            endif

            lx = lx_thin(nxl)
            if(on_comp_domain(nxl,lx)) then

                nz   = kxland(nxl,2)             ! kxland(nxl,2) is the deepest insertion level
                ixl  = ixland(nxl,lx)  -Dom%ioff
                jxl  = jxland(nxl,lx)  -Dom%joff
                ixl2 = ixland(nxl,3-lx)-Dom%ioff
                jxl2 = jxland(nxl,3-lx)-Dom%joff

                ! compute fraction inserted into each thin grid cell                
                thkocean = 0.0
                do k=1,nz
                   thkocean = thkocean + Thickness%rho_dzt(ixl,jxl,k,tau)
                enddo
                do k=1,nz
                   delta(k) = Thickness%rho_dzt(ixl,jxl,k,tau)/(epsln+thkocean)
                enddo

                if(debug_this_module) then 
                    tottracer_thin = 0.0
                    do k=1,nz
                       if (k.eq. 1) then
                         write(*,*) k, tracerx(ixl,jxl,k,n), tracerx(ixl2,jxl2,k,n), ' (nxl=',nxl,' tracer#=',n,')' &
                             ,' thin Thickness%rho_dzt =',Thickness%rho_dzt(ixl,jxl,k,tau)
                       else
                         write(*,*) k, tracerx(ixl,jxl,k,n), '       (nxl=',nxl,' tracer#=',n,')' &
                             ,' thin Thickness%rho_dzt =',Thickness%rho_dzt(ixl,jxl,k,tau)
                       endif
                       ! for tottracer_thin, add in tracer content of destination (thin) column k-levels
                       tottracer_thin  =  tottracer_thin + &
                         ( tracerx(ixl,jxl,k,n)*Thickness%rho_dzt(ixl,jxl,k,tau)*area(nxl,lx) )
                    enddo
                    write(*,*)' Does Thickness%rho_dzt(ixl,jxl,1,tau) = ',Thickness%rho_dzt(ixl,jxl,1,tau), &
                              ' = ',top_rho_dzt(nxl,lx),' = top_rho_dzt(xnl,lx) ?', &
                           '(nxl=',nxl,' tracer#=',n,')'
                    write(*,*)' BEFORE xland_insert: tottracer_thin = ',tottracer_thin, &
                              ' (nxl=',nxl,' tracer#=',n,')'
                    write(*,*)' ' 
                    write(*,*)' tracernew(k) = ((tracerextra*zextra) +', &
                              '  (tracerx(ixl,jxl,k,n)*Thickness%rho_dzt(ixl,jxl,k,tau)) +', &
                              '  (tracerx(ixl2,jxl2,1,n)*zinsert) )', &
                              ' / (zextra +  Thickness%rho_dzt(ixl,jxl,k,tau) + zinsert)'
                endif

! Insert tracer into thin column by doing an xland_insert adjustment process.
! and then diagnose the effective tendencies from the tau+1 - tau timelevel values.
! Start at the bottom (k=nz) of the destination (thin) column and compute what the 
! tracer value would be at the end of the timestep if xland_insert was the only 
! process contributing to the tracer tendency. Then move up the thin column, one 
! k-level at a time, and at each k-level from k=nz-1 up to k=1 compute what the updated
! tracer value would be accounting for both the insertion of the water from the top of
! the thick column +and+ the bubbling up of the excess mass of water from the level
! below. When the process is done, compare the updated tracer values with the current
! tracer values and convert back into a tendency.

! zextra is the vertical thickness of the excess water bubbling up from the
!        layer below due to the accumulated effects of xland_insert
! tracerextra is the tracer value at the tau+1 time step of the water bubbling
!        up from the layer below
! zinsert is the surface height change that will occur due to xland_insert
!        at the thin (destination) point over the course of this tracer timestep.
! tracernew would be the tau+1 time level tracer value value if the xland_insert
!        adjustment was the only process contributing to the tendencies

                zextra = 0.0

                do k=nz,1,-1

                  tracernew(k) = 0.0 

                  if (k .EQ. nz) then
                    tracerextra = 0.0
                  else
                    tracerextra = tracernew(k+1)
                  endif
                  zinsert = water_rate(nxl,lx)*dtime*delta(k)
                  tracernew(k) = ((tracerextra*zextra) + &
                                  (tracerx(ixl,jxl,k,n)*Thickness%rho_dzt(ixl,jxl,k,tau)) + &
                                  (tracerx(ixl2,jxl2,1,n)*zinsert) ) &
                               / (zextra +  Thickness%rho_dzt(ixl,jxl,k,tau) + zinsert)

                  if(debug_this_module) then 
                      write(*,*)' k=',k, tracernew(k),'   (nxl=',nxl,' tracer#=',n,')'
                      write(*,*)'[nxl=',nxl,'] ',tracernew(k),' = ((',tracerextra,'*',zextra,') +'
                      write(*,*)'[nxl=',nxl,']  (',tracerx(ixl,jxl,k,n),'*',Thickness%rho_dzt(ixl,jxl,k,tau),') +'
                      write(*,*)'[nxl=',nxl,']  (',tracerx(ixl2,jxl2,1,n),'*',zinsert,') )'
                      write(*,*)'[nxl=',nxl,'] / (',zextra,' + ',Thickness%rho_dzt(ixl,jxl,k,tau),' + ',zinsert,')' 
                  endif

                  zextra = zextra + zinsert                  

                enddo

                if(debug_this_module) then 
                  write(*,*)' ' 
                  write(*,*) ' CHECK:   (thin)  sfc mass(tau+1)=',area(nxl,lx)*rho_dzt_new_small(nxl), &
                             ' (nxl=',nxl,')'
                  write(*,*) ' CHECK:   (thin) Does zextra =',zextra,' = ', &
                             '(rho_dzt_new_small(nxl) - rho_dzt_x(ixl,jxl)) = ',(rho_dzt_new_small(nxl) - rho_dzt_x(ixl,jxl)), &
                             '(dtime*water_rate(nxl,lx_thin(nxl)))=',(dtime*water_rate(nxl,lx_thin(nxl)))
 
                  ! accumulate tottracer_thin at tau+1 
                  tottracer_thin = 0.0
                  tottracer_thin =  tracernew(1)*area(nxl,lx)* &
                                    (Thickness%rho_dzt(ixl,jxl,1,tau)-rho_dzt_x(ixl,jxl)+rho_dzt_new_small(nxl))
                  if(debug_this_module) then 
                    write(*,*) ' >nxl=',nxl,' tracer#=',n,') tottracer_thin = ', tottracer_thin,' = '
                    write(*,*) ' >nxl=',nxl,'  ',tracernew(1),'*',area(nxl,lx),'*'
                    write(*,*) ' >nxl=',nxl,' (',Thickness%rho_dzt(ixl,jxl,1,tau),'-',rho_dzt_x(ixl,jxl),&
                               '+',rho_dzt_new_small(nxl),')'
                    write(*,*) ' >nxl=',nxl,' which is [',(Thickness%rho_dzt(ixl,jxl,1,tau)-rho_dzt_x(ixl,jxl)&
                               +rho_dzt_new_small(nxl)),']'
                  endif
                  do k=1,nz
                     write(*,*) k, tracernew(k), tracerx(ixl2,jxl2,k,n)
                     if (k .NE. 1) then
                       tottracer_thin   = tottracer_thin   + &
                         (tracernew(k)*area(nxl,lx)*Thickness%rho_dzt(ixl,jxl,k,tau))
                       if(debug_this_module) then 
                         write(*,*) ' >nxl=',nxl,' tottracer_thin = ', tottracer_thin,' = tottracer_thin'
                         write(*,*) ' >nxl=',nxl,' (',tracernew(k),'*',area(nxl,lx),'*',Thickness%rho_dzt(ixl,jxl,k,tau),')'
                       endif
                     endif
                  enddo
                  write(*,*)'  AFTER xland_insert: tottracer_thin = ',tottracer_thin, &
                              ' (nxl=',nxl,' tracer#=',n,')'
                endif

                ! convert change in tracer amounts between tau a tau+1 timesteps
                !  into tendencies and store into the wrk1 array
                blanda = 0.0

                do k=1,nz

                  if (k .eq. 1) then
                     blanda = ( ( tracernew(1)*(Thickness%rho_dzt(ixl,jxl,1,tau)-rho_dzt_x(ixl,jxl)+rho_dzt_new_small(nxl)) )  &
                              - ( tracerx(ixl,jxl,1,n)*Thickness%rho_dzt(ixl,jxl,1,tau) ) )/dtime
                   else
                     blanda = Thickness%rho_dzt(ixl,jxl,k,tau)*(tracernew(k) - tracerx(ixl,jxl,k,n))/dtime
                   endif
                     wrk1(ixl,jxl,k) = wrk1(ixl,jxl,k) + blanda

                   if (debug_this_module) then 
                     if (k .eq. 1) then
                       totrate_thin = 0.0
                     endif
                     totrate_thin = totrate_thin + (blanda*area(nxl,lx))
                     write(*,*) 'xland_insert tendency = ',blanda,'(nxl=',nxl,' tracer#=',n,' k=',k,' thin)'
                   endif

                end do

                if (debug_this_module) then 
                  write(*,*) ' '
                  write(*,*) 'xland_insert totrate_thin  = ',totrate_thin , &
                             '(nxl=',nxl,' tracer#=',n,' thin )'
                endif

            endif ! end of on_comp_domain

        endif ! end of at_least_one_in_comp_domain(nxl)

     enddo  ! end of nxland 


     ! fill tendency array 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + wrk1(i,j,k)
           enddo
        enddo
     enddo

     ! send diagnostics 
     if(id_xland(n) > 0) then 
        call diagnose_3d(Time, Grd, id_xland(n), T_prog(n)%conversion*wrk1(:,:,:))
     endif

     ! for watermass diagnostics 
     if(n==index_temp) then 
         theta_tend(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  theta_tend(i,j,k) = wrk1(i,j,k)
               enddo
            enddo
         enddo
     endif
     if(n==index_salt) then 
         salt_tend(:,:,:) = 0.0
         do k=1,nk
            do j=jsc,jec   
               do i=isc,iec
                  salt_tend(i,j,k) = wrk1(i,j,k)
               enddo
            enddo
         enddo
     endif


  enddo ! end of num_prog_tracers do-loop 

  ! fill mass source array 
  call xlandinsert_mass(Time, Thickness)

  ! call some diagnostics 
  call watermass_diag(Time, Dens)


end subroutine xlandinsert
! </SUBROUTINE> NAME="xlandinsert"


!#######################################################################
! <SUBROUTINE NAME="xlandinsert_mass">
!
! <DESCRIPTION>
! Compute the mass source (kg/m^3)*meter/sec. 
!
! Note that xlandinsert has already been called, so xland_rho_dzt 
! has been filled.  
! </DESCRIPTION>
!
subroutine xlandinsert_mass (Time, Thickness)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(inout) :: Thickness
  integer                                   :: i,j
 
  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_xlandinsert_mod (xlandinsert_mass): module must be initialized')
  endif 

  do j=jsc,jec
     do i=isc,iec     
        Thickness%mass_source(i,j,1) = Thickness%mass_source(i,j,1) + xland_rho_dzt(i,j)
     enddo
  enddo

  call mpp_update_domains (Thickness%mass_source(:,:,1), Dom%domain2d)

  call diagnose_2d(Time, Grd, id_xland_rho_dzt, xland_rho_dzt(:,:))

end subroutine xlandinsert_mass
! </SUBROUTINE> NAME="xlandinsert_mass"


!#######################################################################
! <FUNCTION NAME="at_least_one_in_comp_domain">
!
! <DESCRIPTION>
! Function to see if at least one of the two points in a crossland pair
! is within the computational domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandinsert pair
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
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandinsert pair
! </IN>
! <IN NAME="lx" TYPE="integer">
! lx=1,2 labels the point within an xlandinsert pair
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

  id_neut_rho_xinsert = register_diag_field ('ocean_model', 'neut_rho_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                     &
   'update of locally ref potrho from xlandinsert',                            &
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_xinsert > 0) compute_watermass_diag = .true. 

  id_neut_rho_xinsert_on_nrho = register_diag_field ('ocean_model',                    &
   'neut_rho_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,             &
   'update of locally ref potrho from xlandinsert as binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_xinsert = register_diag_field ('ocean_model', 'wdian_rho_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
    'dianeutral mass transport due to xlandinsert',                              &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_xinsert > 0) compute_watermass_diag = .true. 

  id_wdian_rho_xinsert_on_nrho = register_diag_field ('ocean_model',                &
   'wdian_rho_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,         &
   'dianeutral mass transport from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_xinsert = register_diag_field ('ocean_model',  &
   'tform_rho_xinsert', Grd%tracer_axes(1:3), Time%model_time,&
   'watermass transform from xlandinsert',                    &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_xinsert_on_nrho = register_diag_field ('ocean_model',          &
   'tform_rho_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,   &
   'watermass transform from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_eta_tend_xinsert= -1          
  id_eta_tend_xinsert= register_diag_field ('ocean_model','eta_tend_xinsert',&
       Grd%tracer_axes(1:2), Time%model_time,                                &
       'non-Bouss steric sea level tendency from xinsert tendency', 'm/s',   &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_xinsert > 0) compute_watermass_diag=.true.

  id_eta_tend_xinsert_glob= -1          
  id_eta_tend_xinsert_glob= register_diag_field ('ocean_model', 'eta_tend_xinsert_glob',&
       Time%model_time,                                                                 &
       'global mean non-bouss steric sea level tendency from xinsert tendency',         &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_xinsert_glob > 0) compute_watermass_diag=.true.

 
  ! temperature effects
  id_neut_temp_xinsert = register_diag_field ('ocean_model', 'neut_temp_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
   'temp related update of locally ref potrho from xlandinsert',                 &
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_xinsert > 0) compute_watermass_diag = .true. 

  id_neut_temp_xinsert_on_nrho = register_diag_field ('ocean_model',                                &
   'neut_temp_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                         &
   'temp related update of locally ref potrho from xlandinsert as binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_temp_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_temp_xinsert = register_diag_field ('ocean_model', 'wdian_temp_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'temp related dianeutral mass transport due to xlandinsert',                   &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_xinsert > 0) compute_watermass_diag = .true. 

  id_wdian_temp_xinsert_on_nrho = register_diag_field ('ocean_model',                            &
   'wdian_temp_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
   'temp related dianeutral mass transport from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_temp_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_xinsert = register_diag_field ('ocean_model',  &
   'tform_temp_xinsert', Grd%tracer_axes(1:3), Time%model_time,&
   'temp related watermass transform from xlandinsert',        &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_temp_xinsert_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_temp_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'temp related watermass transform from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_temp_xinsert_on_nrho > 0) compute_watermass_diag = .true. 


  ! salinity effects
  id_neut_salt_xinsert = register_diag_field ('ocean_model', 'neut_salt_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                       &
   'salt related update of locally ref potrho from xlandinsert',                 &
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_xinsert > 0) compute_watermass_diag = .true. 

  id_neut_salt_xinsert_on_nrho = register_diag_field ('ocean_model',                                &
   'neut_salt_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                         &
   'salt related update of locally ref potrho from xlandinsert as binned to neutral density layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_salt_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_salt_xinsert = register_diag_field ('ocean_model', 'wdian_salt_xinsert',&
    Grd%tracer_axes(1:3), Time%model_time,                                         &
    'salt related dianeutral mass transport due to xlandinsert',                   &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_xinsert > 0) compute_watermass_diag = .true. 

  id_wdian_salt_xinsert_on_nrho = register_diag_field ('ocean_model',                            &
   'wdian_salt_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
   'salt related dianeutral mass transport from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_salt_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_xinsert = register_diag_field ('ocean_model',  &
   'tform_salt_xinsert', Grd%tracer_axes(1:3), Time%model_time,&
   'salt related watermass transform from xlandinsert',        &
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_salt_xinsert_on_nrho = register_diag_field ('ocean_model',                      &
   'tform_salt_xinsert_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,               &
   'salt related watermass transform from xlandinsert as binned to neutral density layers',&
   'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_salt_xinsert_on_nrho > 0) compute_watermass_diag = .true. 

  allocate(theta_tend(isd:ied,jsd:jed,nk))
  allocate(salt_tend(isd:ied,jsd:jed,nk))
  theta_tend(:,:,:) = 0.0
  salt_tend(:,:,:)  = 0.0


  if(compute_watermass_diag) then 
    write(stdoutunit,'(/a/)') &
    '==>Note: running ocean_xlandinsert_mod w/ compute_watermass_diag=.true.'  
  endif 


end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from xlandinsert on watermass transformation.  
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Dens)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_density_type),       intent(in) :: Dens

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_xlandinsert (watermass_diag): module needs initialization ')
  endif 

  if(.not. compute_watermass_diag) return 

  tau = Time%tau

  ! full effects from temperature and salinity 
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

  call diagnose_3d(Time, Grd, id_neut_rho_xinsert, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_xinsert, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_xinsert, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_xinsert_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_xinsert_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_xinsert_on_nrho, wrk4)

  if(id_eta_tend_xinsert > 0 .or. id_eta_tend_xinsert_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_xinsert, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_xinsert_glob, eta_tend, cellarea_r)
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

  call diagnose_3d(Time, Grd, id_neut_temp_xinsert, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_temp_xinsert, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_temp_xinsert, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_temp_xinsert_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_temp_xinsert_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_temp_xinsert_on_nrho, wrk4)

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

  call diagnose_3d(Time, Grd, id_neut_salt_xinsert, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_salt_xinsert, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_salt_xinsert, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_salt_xinsert_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_salt_xinsert_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_salt_xinsert_on_nrho, wrk4)
 
end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_xlandinsert_mod

#include <fms_platform.h>
module fv_arrays_mod
!! module for declaring main arrays of FV model. these arrays were
!! declared in fv_pack, now they will be declared here. If pointer
!! sharing is turned on, only the pointers will be published, or by
!! default the arrays.
  use mpp_mod, only: mpp_error, mpp_pe, mpp_npes, mpp_sync, &
       stderr, mpp_root_pe, FATAL, WARNING, &
       mpp_chksum, mpp_send, mpp_recv, &
       mpp_get_current_pelist, mpp_set_current_pelist, mpp_declare_pelist
  use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_domains, &
       mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_pset_mod, only: mpp_pset_type, mpp_pset_create, mpp_pset_delete, &
       mpp_pset_sync, mpp_pset_broadcast, mpp_pset_broadcast_ptr, &
       mpp_pset_check_ptr, mpp_pset_segment_array, mpp_pset_stack_push, &
       mpp_pset_stack_reset, mpp_pset_print_chksum, mpp_pset_root, &
       mpp_pset_get_root_pelist, mpp_pset_print_stack_chksum

  use mod_comm, only: mp_init

  use ensemble_manager_mod, only: get_ensemble_size, get_ensemble_pelist, &
                                  get_ensemble_filter_pelist

  implicit none
  private

  logical :: debug=.TRUE.
!! declarations from fv_pack.F90

!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.

   real, allocatable :: u(:,:,:)    ! D grid zonal wind (m/s)
   real, allocatable :: v(:,:,:)    ! D grid meridional wind (m/s)
   real, allocatable :: pt(:,:,:)   ! temperature (K)
   real, allocatable :: delp(:,:,:) ! pressure thickness (pascal)
   real, allocatable :: q(:,:,:,:)  ! specific humidity and constituents

! phys_state: use these for physics packages that require tendencies as input
! The following 4 arrays are allocated only when use_tendency = .true.
   real, allocatable :: u_phys(:,:,:)     ! u-wind on the physics grid (m/s)
   real, allocatable :: v_phys(:,:,:)     ! v-wind on the physics grid (m/s)
   real, allocatable :: t_phys(:,:,:)     ! temperature (K)
   real, allocatable :: q_phys(:,:,:,:)   ! tracer (moist) mixing ratios

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
   real, allocatable :: ps (:,:)      ! Surface pressure (pascal)
   real, allocatable :: pe (:,:,: )   ! edge pressure (pascal)
   real, allocatable :: pk  (:,:,:)   ! pe**cappa
   real, allocatable :: peln(:,:,:)   ! ln(pe)
   real, allocatable :: pkz (:,:,:)   ! finite-volume mean pk

   real, allocatable :: ps_bp (:,:)   ! ps BEFORE physics update
                                          ! for computing omega
   real, allocatable :: phis(:,:)     ! Surface geopotential (g*Z_surf)
   real, allocatable :: omga(:,:,:)   ! Vertical pressure velocity (pa/s)
   real, allocatable :: pesouth(:,:)  ! "pe at the south" for SPMD

   real, allocatable:: ua(:,:,:)      ! (ua, va) are mostly used as the A grid winds
   real, allocatable:: va(:,:,:)

! Surface winds for coupler
   real, allocatable:: u_srf(:,:)     ! A grid u-wind at surface
   real, allocatable:: v_srf(:,:)     ! A grid v-wind at surface

   real, allocatable, dimension(:,:,:) :: u_dt, v_dt, t_dt
   real, allocatable :: q_dt(:,:,:,:)

#ifdef MARS_GCM
   real, allocatable :: delp_dt(:,:,:)
   real, allocatable :: mars_sfc_budg(:,:,:)
#endif MARS_GCM

!! declarations from fv_pack.F90.1D2D
#ifdef MARS_GCM
  integer, parameter :: NUMPTRS=28
#else
  integer, parameter :: NUMPTRS=26
#endif MARS_GCM

  integer(POINTER_KIND) :: ptrArray(NUMPTRS) !this array contains addresses to all the variables to be shared
                              !between physics and dynamics.
  integer(POINTER_KIND) :: ptr_u, &
                           ptr_v, &
                           ptr_delp, &
                           ptr_pt, &
                           ptr_q, &
                           ptr_u_phys, &
                           ptr_v_phys, &
                           ptr_t_phys, &
                           ptr_q_phys, &
                           ptr_phis, &
                           ptr_ps, &
                           ptr_omga, &
                           ptr_pkz, &
                           ptr_pk, &
                           ptr_pe, &
                           ptr_peln, &
                           ptr_pesouth, &
                           ptr_ua, &
                           ptr_va, &
                           ptr_ps_bp, &
                           ptr_u_srf, &
                           ptr_v_srf, &
                           ptr_u_dt, &
                           ptr_v_dt, &
                           ptr_t_dt, &
#ifdef MARS_GCM
                           ptr_q_dt,  ptr_delp_dt, ptr_mars_sfc_budg  
#else
                           ptr_q_dt
#endif

!note that npes/layout vars are module global: implies that the same pelist
!must call the init and other methods... should be ok since all atmos PEs
!call it
  integer :: npe, npes, npex, npey
  integer :: npel, nper !thread on the left and right
  integer :: pos !pos = modulo(npe,npex) !0.LE.pos.LE.npex-1

  integer, allocatable, dimension(:) :: pelist, ypelist
  logical :: pset_root !TRUE for PEs that allocate shared arrays (ypelist)
!compute/data/global domain limits
  integer, public :: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg
!the ..p integers are the loop indices corresponding to a parallel loop
  integer, public :: isp, iep, jsp, jep, ksp, kep
!number of levels and tracers; loop decomposition in East-West direction
  integer, public :: nlev, ncnst, n_zonal=0
  integer, public :: layout(2)=(/0,0/) !layout, public to fv_core_nml
  integer, public :: fv_pset_stack=1000000 !internal workspace for PSETs

!MPI communicator for ypelist
  integer :: commID
!PSET
  type(mpp_pset_type), public, save :: pset

  public :: fv_arrays_init, fv_arrays_exit, fv_stack_push, &
       fv_print_chksum, fv_print_chksums, fv_thread_bcast, &
       fv_array_check, fv_array_limits, fv_stack_reset, fv_array_sync
#ifdef use_shared_pointers
!arrays pointers are public
  public :: ptr_u, ptr_v, ptr_delp, ptr_pt, ptr_q, ptr_u_phys, ptr_v_phys, &
            ptr_t_phys, ptr_q_phys, ptr_phis, ptr_ps, ptr_omga, ptr_pkz, &
            ptr_pk, ptr_pe, ptr_peln, ptr_pesouth, ptr_ua, ptr_va, ptr_ps_bp, &
            ptr_u_srf, ptr_v_srf, ptr_u_dt, ptr_v_dt, ptr_t_dt, ptr_q_dt
#   ifdef MARS_GCM
  public ::  ptr_delp_dt, ptr_mars_sfc_budg
#   endif MARS_GCM

#else
!arrays themselves are public
  public :: u, v, delp, pt, q, u_phys, v_phys, &
            t_phys, q_phys, phis, ps, omga, pkz, &
            pk, pe, peln, pesouth, ua, va, ps_bp, &
            u_srf, v_srf, u_dt, v_dt, t_dt, q_dt
#   ifdef MARS_GCM
  public ::  delp_dt, mars_sfc_budg
#   endif MARS_GCM

#endif use_shared_pointers

contains

  subroutine fv_arrays_init( nx, ny, nz, nt, shalo, nhalo, ighost, domain )
    integer, intent(in) ::   nx, ny, nz, nt, shalo, nhalo, ighost 
    type(domain2D), intent(out) :: domain

    integer               :: ensemble_siz(6), num_ensemble, atmos_npes_pm, n
    integer, allocatable :: ensemble_pelist(:,:), all_pelist(:)

!convert to domain2D specification
!this requires ng_s=ng-d for now
    if( shalo.NE.nhalo )call mpp_error( FATAL, &
         'FV_ARRAYS_INIT: now requires even halos (ng_s==ng_d)' )
    npe = mpp_pe()
    npes = mpp_npes()
    if( layout(1).EQ.0 .AND. layout(2).EQ.0 )then
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    else if( layout(1).EQ.0 )then
        layout(1) = npes/layout(2)
    else if( layout(2).EQ.0 )then
        layout(2) = npes/layout(1)
    end if
    npex = layout(1); npey = layout(2)
    if( npex*npey.NE.npes )call mpp_error( FATAL, &
         'FV_ARRAYS_INIT: npex*npey.NE.npes: check if layout in fv_core_nml matches atmos_npes.' )
!latitudes must be evenly divided and more than 3 rows per PE
    if( modulo(ny,npey).NE.0 )call mpp_error( FATAL, &
         'FV_ARRAYS_INIT: mlat must be evenly decomposed.' )
    if( ny/npey.LT.3 )call mpp_error( FATAL, &
         'FV_ARRAYS_INIT: FV core requires at least 3 rows per PE.' )
    if( ighost.NE.0 )call mpp_error( FATAL, &
         'FV_ARRAYS_INIT: currently assumes phis has no halo (ighost==0)' )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain, &
         xhalo=0, yhalo=shalo, name='FV Core Arrays.' )
    call mpp_get_compute_domain( domain, is, ie, js, je  )
    call mpp_get_data_domain   ( domain, isd,ied,jsd,jed )
    call mpp_get_global_domain ( domain, isg,ieg,jsg,jeg )
    nlev = nz; ncnst = nt
!#ifndef use_shared_pointers
!    if ( n_zonal == 0 ) then
!        if ( (ny/npey) / npex >= 6 ) then
!            n_zonal = 1
!        else
!            if ( 4*(nx/4) == nx ) Then
!                n_zonal = 4
!            else
!                n_zonal = 2
!            endif
!        endif
!    endif
!#else
!in our case n_zonal MUST EQUAL npex!
!fv_arrays_init is called after reading fv_core_nml so it overrides namelist
!value... but better to remove it from namelist later
    if( n_zonal.NE.0 )then
        call mpp_error( WARNING, &
             'FV_ARRAYS_INIT: n_zonal was set in fv_core_nml '// &
             'but overridden here.' )
    end if
    n_zonal = npex
!#endif

    allocate( pelist(0:npes-1) )
    call mpp_get_current_pelist(pelist)
!set up PSETs
!z1l: mpp_pset_create need to be called for each ensemble because of mpp_declare_pelist
    ensemble_siz = get_ensemble_size()
    num_ensemble = ensemble_siz(1)
    if(num_ensemble == 1) then
       call mpp_pset_create( npex, pset, fv_pset_stack)
    else
       atmos_npes_pm = ensemble_siz(4)
       allocate(ensemble_pelist(num_ensemble, atmos_npes_pm) )
       allocate(all_pelist(num_ensemble*atmos_npes_pm) )
       call get_ensemble_pelist(ensemble_pelist, 'atmos')
       call get_ensemble_filter_pelist(all_pelist, 'atmos')
       call mpp_set_current_pelist(all_pelist)
       call mpp_get_current_pelist(all_pelist, commID = commID)
       do n = 1, num_ensemble
             call mpp_pset_create( npex, pset, fv_pset_stack, ensemble_pelist(n,:), commID ) !initialized npey PSETs of npex each
       enddo
       call mpp_set_current_pelist(pelist)
    endif

!the y pelist consists of 1 PE from each row(y) of the decomposition
    allocate( ypelist(0:npey-1) )
    call mpp_pset_get_root_pelist( pset, ypelist, commID=commID )

!all PEs allocate by default; but only root PEs allocate when use_shared_pointers
    pset_root = .TRUE.
#ifdef use_shared_pointers
    pset_root = mpp_pset_root(pset)
#endif
!only the root PEs should do these allocate()s!
    if( pset_root )then
! State variables:
        allocate (   u(isg:ieg, jsd:jed, nz))
        allocate (   v(isg:ieg, jsd:jed, nz))
        allocate (delp(isg:ieg, js:je,   nz))
        allocate (  pt(isg:ieg, jsd:jed, nz))
        allocate (   q(isg:ieg, jsd:jed, nz, nt))

!      if( use_tendency ) then
        allocate ( u_phys(isg:ieg, js:je, nz) )
        allocate ( v_phys(isg:ieg, js:je, nz) )
        allocate ( t_phys(isg:ieg, js:je, nz) )
        allocate ( q_phys(isg:ieg, js:je, nz, nt) )
!      endif

        allocate (phis(isg:ieg, js:je))
        allocate (ps  (isg:ieg, js:je))
        allocate (omga(isg:ieg, js:je,nz))    

        allocate (pkz (isg:ieg, js:je, nz))
        allocate (pk  (isg:ieg, js:je, nz+1))
        allocate (pe  (isg:ieg, nz+1, js:je)) 
        allocate (peln(isg:ieg, nz+1, js:je)) 

        allocate (pesouth(isg:ieg, nz+1))

        allocate ( ua(isg:ieg, js:je, nz) )
        allocate ( va(isg:ieg, js:je, nz) )

! For GFDL physics:
        allocate (ps_bp (isg:ieg, js:je))
        allocate ( u_srf(isg:ieg, js:je) )
        allocate ( v_srf(isg:ieg, js:je) )
        allocate( u_dt(isg:ieg, js:je, nz) )
        allocate( v_dt(isg:ieg, js:je, nz) )
        allocate( t_dt(isg:ieg, js:je, nz) )
        allocate( q_dt(isg:ieg, js:je, nz, nt) )
#ifdef MARS_GCM
        allocate( delp_dt(isg:ieg, js:je, nz) )
        allocate( mars_sfc_budg(isg:ieg, js:je, 5) )
#endif MARS_GCM
    end if
    isp = isg; iep = ieg
    jsp = js; jep = je
    ksp = 1; kep = nz
#ifdef use_shared_pointers
    if( pset_root )then
    ptr_u = LOC(u)
    ptr_v = LOC(v)
    ptr_delp = LOC(delp)
    ptr_pt = LOC(pt)
    ptr_q = LOC(q)
    ptr_u_phys = LOC(u_phys)
    ptr_v_phys = LOC(v_phys)
    ptr_t_phys = LOC(t_phys)
    ptr_q_phys = LOC(q_phys)
    ptr_phis = LOC(phis)
    ptr_ps = LOC(ps)
    ptr_omga = LOC(omga)
    ptr_pkz = LOC(pkz)
    ptr_pk = LOC(pk)
    ptr_pe = LOC(pe)
    ptr_peln = LOC(peln)
    ptr_pesouth = LOC(pesouth)
    ptr_ua = LOC(ua)
    ptr_va = LOC(va)
    ptr_ps_bp = LOC(ps_bp)
    ptr_u_srf = LOC(u_srf)
    ptr_v_srf = LOC(v_srf)
    ptr_u_dt = LOC(u_dt)
    ptr_v_dt = LOC(v_dt)
    ptr_t_dt = LOC(t_dt)
    ptr_q_dt = LOC(q_dt)
#ifdef MARS_GCM
    ptr_delp_dt = LOC(delp_dt)
    ptr_mars_sfc_budg = LOC(mars_sfc_budg)
#endif MARS_GCM
!the array ptrArray is now filled by equivalence()
    ptrArray( 1) = ptr_u
    ptrArray( 2) = ptr_v
    ptrArray( 3) = ptr_delp
    ptrArray( 4) = ptr_pt
    ptrArray( 5) = ptr_q
    ptrArray( 6) = ptr_u_phys
    ptrArray( 7) = ptr_v_phys
    ptrArray( 8) = ptr_t_phys
    ptrArray( 9) = ptr_q_phys
    ptrArray(10) = ptr_phis
    ptrArray(11) = ptr_ps
    ptrArray(12) = ptr_omga
    ptrArray(13) = ptr_pkz
    ptrArray(14) = ptr_pk
    ptrArray(15) = ptr_pe
    ptrArray(16) = ptr_peln
    ptrArray(17) = ptr_pesouth
    ptrArray(18) = ptr_ua
    ptrArray(19) = ptr_va
    ptrArray(20) = ptr_ps_bp
    ptrArray(21) = ptr_u_srf
    ptrArray(22) = ptr_v_srf
    ptrArray(23) = ptr_u_dt
    ptrArray(24) = ptr_v_dt
    ptrArray(25) = ptr_t_dt
    ptrArray(26) = ptr_q_dt
#ifdef MARS_GCM
    ptrArray(27) = ptr_delp_dt
    ptrArray(28) = ptr_mars_sfc_budg
#endif MARS_GCM
    endif
    call mpp_pset_broadcast_ptr( pset, ptrArray )
    ptr_u = ptrArray( 1)
    ptr_v = ptrArray( 2)
    ptr_delp = ptrArray( 3)
    ptr_pt = ptrArray( 4)
    ptr_q = ptrArray( 5)
    ptr_u_phys = ptrArray( 6)
    ptr_v_phys = ptrArray( 7)
    ptr_t_phys = ptrArray( 8)
    ptr_q_phys = ptrArray( 9)
    ptr_phis = ptrArray(10)
    ptr_ps = ptrArray(11)
    ptr_omga = ptrArray(12)
    ptr_pkz = ptrArray(13)
    ptr_pk = ptrArray(14)
    ptr_pe = ptrArray(15)
    ptr_peln = ptrArray(16)
    ptr_pesouth = ptrArray(17)
    ptr_ua = ptrArray(18)
    ptr_va = ptrArray(19)
    ptr_ps_bp = ptrArray(20)
    ptr_u_srf = ptrArray(21)
    ptr_v_srf = ptrArray(22)
    ptr_u_dt = ptrArray(23)
    ptr_v_dt = ptrArray(24)
    ptr_t_dt = ptrArray(25)
    ptr_q_dt = ptrArray(26)
#ifdef MARS_GCM
    ptr_delp_dt = ptrArray(27)
    ptr_mars_sfc_budg = ptrArray(28)
#endif MARS_GCM
!isp,... index limits when you divide the computational domain by npex
    isp = is
    iep = ie
    call fv_array_limits( js, je, jsp, jep )
    call fv_array_limits(  1, nz, ksp, kep )
#endif
!initialize mod_comm on the ypelist
    if( pset_root )then
        call mpp_set_current_pelist(ypelist)
        call mpp_get_current_pelist(ypelist, commID=commID)
        call mp_init( nx, ny, nz, commID )
    end if
    call mpp_set_current_pelist(pelist)
    return
  end subroutine fv_arrays_init

  subroutine fv_thread_bcast(a)
!broadcast value on the pset_root to its sub-threads
!wrapper, eliminate
    real, intent(inout) :: a

    call mpp_pset_broadcast( pset, a )
  end subroutine fv_thread_bcast

  subroutine fv_array_sync
!this is a replacement for mpp_sync, doing syncs across
!shared arrays without calling mpp_sync
!wrapper, eliminate
    call mpp_pset_sync(pset)
  end subroutine fv_array_sync

  subroutine fv_array_check(ptr)
!check if shared array
!wrapper, eliminate
#ifdef use_CRI_pointers
    real :: a
    pointer( ptr, a )
#else
    integer(POINTER_KIND), intent(in) :: ptr
#endif
    call mpp_pset_check_ptr( pset, ptr )
  end subroutine fv_array_check

  subroutine fv_array_limits( ls, le, lsp, lep )
!get parallel limits
!wrapper: eliminate
    integer, intent(in) :: ls, le
    integer, intent(out) :: lsp, lep

    call mpp_pset_segment_array( pset, ls, le, lsp, lep )
  end subroutine fv_array_limits

  subroutine fv_stack_push( ptr, len )
!mpp_malloc specialized for FV shared arrays
!len is the length of the required array
!lstack is the stack already in play
!user should zero lstack (call fv_stack_reset) when the stack is to be cleared
!wrapper: eliminate
    integer, intent(in) :: len
#ifdef use_CRI_pointers
    real :: a
    pointer( ptr, a )
#else
    integer(POINTER_KIND), intent(out) :: ptr
#endif
    call mpp_pset_stack_push( pset, ptr, len )
  end subroutine fv_stack_push

  subroutine fv_stack_reset
!reset stack... will reuse any temporary arrays! USE WITH CARE
    call mpp_pset_stack_reset(pset)
  end subroutine fv_stack_reset

  subroutine fv_print_chksum(caller, array)
!wrapper: eliminate
    character(len=*), intent(in) :: caller
    real, intent(in) :: array(:,:,:)

    call mpp_pset_print_chksum( pset, caller, array )
    return
  end subroutine fv_print_chksum

  subroutine fv_print_chksums(caller)
    character(len=*), intent(in) :: caller
#ifdef PSET_DEBUG
#include "fv_point.inc"
!checksums may not match across changes in Y decomposition!

    call mpp_pset_print_chksum( pset, trim(caller)//' u=', &
         u(isg:ieg,jsd:jed,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' v=', &
         v(isg:ieg,jsd:jed,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' delp=', &
         delp(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' pt=', &
         pt(isg:ieg,jsd:jed,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' q=', &
         q(isg:ieg,jsd:jed,:,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' u_phys=', &
         u_phys(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' v_phys=', &
         v_phys(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' t_phys=', &
         t_phys(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' q_phys=', &
         q_phys(isg:ieg,js:je,:,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' phis=', &
         phis(isg:ieg,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' ps=', &
         ps  (isg:ieg,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' omga=', &
         omga(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' pkz=', &
         pkz (isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' pk=', &
         pk  (isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' pe=', &
         pe  (isg:ieg,:,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' peln=', &
         peln(isg:ieg,:,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' pesouth=', &
         pesouth(isg:ieg,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' ua=', &
         ua(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' va=', &
         va(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' ps_bp=', &
         ps_bp (isg:ieg,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' u_srf=', &
         u_srf(isg:ieg,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' v_srf=', &
         v_srf(isg:ieg,js:je) )
    call mpp_pset_print_chksum( pset, trim(caller)//' u_dt=', &
         u_dt(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' v_dt=', &
         v_dt(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' t_dt=', &
         t_dt(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' q_dt=', &
         q_dt(isg:ieg,js:je,:,:) )
#ifdef MARS_GCM
    call mpp_pset_print_chksum( pset, trim(caller)//' delp_dt=', &
         delp_dt(isg:ieg,js:je,:) )
    call mpp_pset_print_chksum( pset, trim(caller)//' mars_sfc_budg=', &
         mars_sfc_budg(isg:ieg,js:je,:) )
#endif MARS_GCM
    call mpp_pset_print_stack_chksum(pset, caller)
#endif
    return
  end subroutine fv_print_chksums

  subroutine fv_arrays_exit
    if( pset_root )then
        deallocate ( u )
        deallocate ( v )
        deallocate ( delp )
        deallocate ( pt )
        deallocate ( q )

        deallocate ( phis )
        deallocate ( ps )
        deallocate ( omga )

        deallocate ( pkz )
        deallocate ( pk )
        deallocate ( pe )
        deallocate ( peln )

        deallocate ( pesouth )

        deallocate ( ua )
        deallocate ( va )


        deallocate( ps_bp )
        deallocate( u_srf )
        deallocate( v_srf )

!      if( use_tendency ) then
        deallocate ( u_phys )
        deallocate ( v_phys )
        deallocate ( t_phys )
        deallocate ( q_phys )
!      endif

#ifdef MARS_GCM
        deallocate ( delp_dt )
        deallocate ( mars_sfc_budg )
#endif MARS_GCM


    end if
    call mpp_pset_delete(pset)
  end subroutine fv_arrays_exit
end module fv_arrays_mod

#ifdef test_fv_arrays
  program test
    use fv_arrays_mod, only: fv_arrays_init, fv_array_limits
    use mpp_mod, only: mpp_init, mpp_exit, mpp_pe

    integer :: nx=360, ny=200, halo=2, nzonal=0
    integer :: isp, iep

    call mpp_init()

    call fv_arrays_init( nx, ny, halo, halo, nzonal )
!test fv_array_limits
!run tests at PE counts that divide le exactly, as well as counts that do not
    ls = 1; le = 1000
    call fv_array_limits( ls, le, lsp, lep )
    write( 6,'(a,4i4)' )'ls,le,lsp,lep=', ls,le,lsp,lep
    call mpp_exit()
  end program test
#endif

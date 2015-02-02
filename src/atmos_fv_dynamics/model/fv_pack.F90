! $Id: fv_pack.F90,v 19.0 2012/01/06 20:00:50 fms Exp $

module fv_pack

   use timingModule
#ifdef USE_LIMA
   use shr_kind_mod, only: r8 => shr_kind_r8
#endif



!-----------------------------------------------------------------------
! "SPMD" stands for Single Program Multiple Data. It activates
! domain decomposition with inter-process communication done
! by MPI-1, MPI-2, or SHMEM -- this is the first level of the
! parallelism. The second level parallelism uses the standard
! "OpenMP". To activate OpenMP, add "-mp" to the compiler and the load
! flags on the SGI.
!-----------------------------------------------------------------------

#ifdef SPMD
   use mod_comm,     only: gid, mp_init, mp_exit, y_decomp,  &
                           mp_send3d, mp_recv3d,                             &
                           mp_scatter3d, mp_send4d_ns, mp_recv4d_ns
#endif

! --- FMS ----
   use fms_mod,      only: file_exist,  open_namelist_file,    &
                           error_mesg,  FATAL,                 &
                           close_file,  check_nml_error,       &
                           write_version_number, stdlog, stdout, stderr
   use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, MPP_VERBOSE, &
        mpp_get_current_pelist, mpp_sync, mpp_error, input_nml_file
   use mpp_domains_mod, only: mpp_define_domains, mpp_domains_init, domain2D, &
        mpp_get_compute_domain
   use constants_mod,    only: radius, omega, grav, kappa, rdgas, cp_air, rvgas

! --- tracer manager ---
   use tracer_manager_mod, only : tm_get_number_tracers => get_number_tracers, &
                                  tm_get_tracer_index   => get_tracer_index,   &
                                  tm_get_tracer_indices => get_tracer_indices, &
                                  tm_set_tracer_profile => set_tracer_profile, &
                                  tm_get_tracer_names   => get_tracer_names,   &
                                  tm_check_if_prognostic=> check_if_prognostic,&
                                  tm_register_tracers   => register_tracers
   use field_manager_mod, only  : MODEL_ATMOS

   use fv_arrays_mod, only: fv_arrays_init, fv_arrays_exit, fv_array_check, &
        fv_stack_push, fv_array_sync, isp, iep, jsp, jep, ksp, kep, layout, &
        nlev, ncnst, n_zonal, fv_print_chksums, fv_pset_stack
#ifndef USE_LIMA
   use pft_module, only: pft_init, pft_end
#else
   use pft_module, only: pft_init
#endif


! --- FMS ----

   implicit none
   private

! --- FMS ----



!==================
! fv_core namelist:
!==================

!----------------
! Dimensionality:
!----------------
   integer nlon              ! E-W (total number of longitudes)
   integer mlat              ! N-S (total number of latitudes)
!   integer nlev              ! Vertical dimension
!   integer ncnst             ! Total number of tracers
   integer :: nt_phys = 4    ! Total number of tracers needed by physics
   integer :: nt_prog = 0    ! Number of prognostic tracers
   integer :: pnats = 0      ! Number of non-advected tracers; default to 0
   logical :: do_fms_tracer_manager = .FALSE.  ! get list of tracers from field_table
! name of fms tracers restart file
   character(len=*), parameter :: fms_tracers_file='atmos_tracers.res'

!------------------------------------------
! Parameters that control the FV algorithms:
!------------------------------------------

! Momentum transport options:
   integer :: iord_mt = 4    ! E-W
   integer :: jord_mt = 4    ! N-S
   integer :: kord_mt = 4    ! vertical mapping option

! Heat transport options:
   integer :: iord_tm = 4    ! E-W
   integer :: jord_tm = 4    ! N-S
   integer :: kord_tm = 7    ! vertical mapping option

! Tracer transport options:
   integer :: iord_tr = 4    ! E-W
   integer :: jord_tr = 4    ! N-S
   integer :: kord_tr = 7    ! vertical mapping option for tracers

   integer :: n_split = 0    ! Number of time splits for the lagrangian dynamics
                             ! Default = 0 (automatic computation of best value)
   integer :: n_spong = 1    ! Number of sponge layers from the top
   integer :: n_diffu = 0    ! Number of diffusive layers from the bottom
!   integer :: n_zonal = 0    ! loop decomposition in East-West direction
#ifndef USE_LIMA
   integer :: dyn_step= 0    ! Multi-Dimensional transport algorithms for dynamics
   integer :: tra_step= 1    ! Multi-Dimensional transport algorithms for tracers
                             ! 0: NO directional split
                             ! 1: start with meridional transport
                             ! 2: start with zonal transport
#endif

   logical :: age_tracer= .false.  ! Carry the age tracer as the last tracer
   logical :: adiabatic = .false.  ! Run without physics (full or idealized).
                                   ! this is mainly to turn off the Held-Suarez
                                   ! forcing in the "solo" mode
   logical :: full_phys = .true.   ! Run with full physics with virtual effects
                                   ! for R (gas constant) and Cp
   real:: consv_te  = 0.       ! Total energy fixer. Range[0,1.]. For example 0.75
                                   ! would provide 75% correction
! Divergence damping parameters: D = a_div + b_div * cos(lat)
   real:: a_div     = 0.6      !
   real:: b_div     = 0.6      !

#ifndef USE_LIMA
!!!1   real:: dry_mass = 98290.
#endif

#ifdef MARS_GCM
   real    :: p_ref = 700.
   real    :: fv_sponge_damp=   10.0
   real    :: fv_sponge_lev= 0.5      !   level (mb) for enhanced divergence damping 
   real    :: dry_mass = 7.01E2
#else
   real    :: p_ref = 1.E5
   real    :: dry_mass = 98290.
#endif


   logical :: cold_start      = .false.       ! read in initial condition
   logical :: change_time     = .false.       ! ignore time stamp in fv_rst.res
   logical :: adjust_dry_mass = .false.       ! adjust initial total dry air mass
   logical :: use_set_eta     = .false.       ! use set_eta(), instead of restart, for (ak,bk)
   logical :: pft_phys        = .true.        ! polar filter physical tendencies
   logical :: use_tendency    = .false.       ! use the tendency approach
   logical :: do_ch4_chem     = .false.       ! Relax H2o to Haloe data between 1 and 10 mb
   integer :: print_freq      = 0             ! Print max/min of selected fields
                                              ! 0: off
                                              ! positive n: every n hours
                                              ! negative n: every time step

#ifdef SW_DYN
   logical :: fill = .false.
#else
   logical :: fill = .true.  ! Activate a simple filling algorithm for tracers
#endif

   integer :: map_dt         ! Remapping time step (seconds); default to model
                             ! time step

! Default case for the shallow water test: Rossby wave-4
   integer :: icase = 6      ! Test case number if "SW_DYN" is defined
                             ! see tools/nit_sw_ic.f90 for details

   character(len=24) :: restart_format = 'native'   ! native or netcdf

   character(len=128) :: version = '$Id: fv_pack.F90 1.1.2.15.2.2'
   character(len=128) :: tagname = '$Name: tikal $'


#ifndef USE_LIMA
   namelist /fv_core_nml/nlon,      mlat,     nlev,                           &
                         ncnst,     nt_phys,  pnats,                          &
                         iord_mt,   iord_tm,  iord_tr,                        &
                         jord_mt,   jord_tm,  jord_tr,                        &
                         kord_mt,   kord_tm,  kord_tr,                        &
                         n_split,   n_zonal,  map_dt,                         &
                         consv_te,  restart_format,  adiabatic, full_phys,    &
                         fill,      adjust_dry_mass, n_spong, n_diffu,        &
                         change_time, a_div, b_div,  use_set_eta,             &
                         use_tendency, do_ch4_chem,  pft_phys, age_tracer,    &
                         tra_step, dyn_step, dry_mass, print_freq,  icase,    &
#ifdef MARS_GCM
                        p_ref, fv_sponge_lev, fv_sponge_damp,                      &
#endif MARS_GCM 
                         layout, fv_pset_stack
#else
   namelist /fv_core_nml/nlon,      mlat,     nlev,                           &
                         ncnst,     nt_phys,  pnats,                          &
                         iord_mt,   iord_tm,  iord_tr,                        &
                         jord_mt,   jord_tm,  jord_tr,                        &
                         kord_mt,   kord_tm,  kord_tr,                        &
                         n_split,   n_zonal,  map_dt,                         &
                         consv_te,  restart_format,  adiabatic, full_phys,    &
                         fill,      adjust_dry_mass, n_spong, n_diffu,        &
                         change_time, a_div, b_div,  use_set_eta,             &
                         use_tendency, do_ch4_chem,  pft_phys, age_tracer,    &
#ifdef MARS_GCM
                         p_ref, fv_sponge_lev, fv_sponge_damp,                      &
#endif MARS_GCM 
                         print_freq,  icase, layout

#endif

   logical mountain          ! If terrain is present
   integer nqrst             ! Total number of tracers in fv_rst file

!-----------------------------------------------------------------
! Domain decomposition parameters (not to be changed by the user):
!-----------------------------------------------------------------
   integer beglon             ! beginning lon index
   integer endlon             ! ending    lon index
   integer beglat             ! beginning lat index
   integer endlat             ! ending    lat index
   integer ng_d               ! algorithm dependent ghost zones
   integer ng_s               ! algorithm dependent ghost zones for u-wind
   logical master             ! If "master process"? used mostly for print
                              ! and read/write restart
! FMS specific:
   type(domain2D), save :: fv_domain


!-----------------------------------------------------------------------
! Eulerian vertical coordinate:
! pe(k) = ak(k) + bk(k) * ps
!-----------------------------------------------------------------------
   real, allocatable :: ak(:)
   real, allocatable :: bk(:)
   integer  ks                            ! number of layers in pure p region
   real ptop                          ! pressure at model top (pascal)

! Horizontal grid geometry
#ifndef USE_LIMA
   real, allocatable:: f_d(:,:)       !  Coriolis parm at cell center
#else
   real, allocatable:: f_d(:)       !  Coriolis parm at cell center
#endif
   real, allocatable:: dx(:)          !  grid distance (m) along latitude circles
   real, allocatable:: one_by_dx(:)   !  1 / dx
   real one_by_dy                     !  1 / dy
   real acap                          ! scaled polar area = nlon * cosp(1)
   real rcap                          ! 1 / acap

!-----------------------------------
! Geometric arrays
!-----------------------------------
   real, allocatable :: lon(:), lonb(:), rlonb(:)
   real, allocatable :: lat(:), latb(:), rlatb(:)
   real, allocatable :: rlat(:,:)
   real, allocatable :: rlon(:,:)
   real, allocatable :: area(:,:)         ! grid cell area: m**2
   real, allocatable :: sine(:)
   real, allocatable :: cose(:)
   real, allocatable :: sinp(:)
   real, allocatable :: cosp(:)
   real, allocatable :: acosp(:)
   real, allocatable :: acosu(:)          ! for u-wind averaging
   real, allocatable :: grid_weight(:)

   real, allocatable :: coslon(:)
   real, allocatable :: sinlon(:)
   real, allocatable :: cosl5(:)
   real, allocatable :: sinl5(:)

   real, parameter::cv_wv = 1463.  ! specific heat of water vapor at constant
                                       ! volume [J/(K KG)]
   real, parameter::cp_wv = 1952.  ! specific heat of water vapor at constant
                                       ! pressure [J/(K KG)]
   real cp_vir
   real pi

#ifdef USE_LIMA
! Polar filter related
   integer  ifax(13)                      !ECMWF fft
   real(r8), allocatable :: trigs(:)
   real(r8), allocatable :: dc(:,:), de(:,:), sc(:), se(:)
   integer:: ipft1, ipft2
#endif

!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
!   real, allocatable :: pesouth(:,:)  ! "pe at the south" for SPMD
   logical :: u_ghosted = .false.         ! u-wind is ghosted at the first lat to the N
   real, parameter::ptop_min = 1.E-6  ! minimum pressure (pa) at model top to avoid
                                          ! floating point exception; this is not needed
                                          ! if model top is not at zero
   real, parameter :: too_big = 1.E30      ! use this for missing values
   real :: age_time                   ! Time (sec) since init of age_tracer


#ifdef USE_LIMA
   integer:: js2g0, jn1g1
#endif
   integer:: ighost
!PUBLIC VARIABLES
   public :: nlon, mlat, nlev, beglon, endlon, beglat, endlat, &
        age_tracer, age_time, do_ch4_chem, u_ghosted, ng_d, ng_s, &
        use_tendency, pft_phys, ak, bk, get_eta_level, ptop, full_phys, &
        coslon, sinlon, ighost, cosp, sinp, cose, sine, master, &
        nt_phys, nt_prog, rlat, rlon, area, rlonb, rlatb, rcap, &
        one_by_dy, one_by_dx, f_d, acosp, ks, lon, lat, lonb, latb, &
        fv_domain, print_freq, do_fms_tracer_manager, nqrst, mountain, &
        adjust_dry_mass, use_set_eta, fms_tracers_file, restart_format, &
        cold_start, gid, acap, acosu, ptop_min, pi, n_spong, n_diffu, &
        a_div, b_div, cosl5, sinl5, dx, n_split, iord_mt, iord_tm, &
        iord_tr, jord_mt, jord_tm, jord_tr, kord_mt, kord_tm, kord_tr, &
        fill, n_zonal, cp_vir, icase, ncnst, pnats, consv_te, &
        change_time, grid_weight
#ifdef USE_LIMA
   public :: sc, se, dc, de, ipft1, ipft2, ifax, trigs, too_big
#else
   public :: dry_mass, tra_step, dyn_step, p_ref
#endif

#ifdef MARS_GCM
   public ::   fv_sponge_lev,  fv_sponge_damp
#endif

!Balaji: needed by solo driver
   public :: map_dt, adiabatic
!PUBLIC ROUTINES
   public fv_init, fv_end, d2a3d, p_var, polavg, drymadj, tracer_mass


! !REVISION HISTORY:
!
!   03.10.02     SJ Lin        Modifications from fvGCM
!   03.11.12     SJ Lin        Inline documentation
!   03.11.21     SJ Lin        Change pt from potential temperature to temperature
!   04.11.29     SJ Lin        Pre-2D decomposition modifications
!-----------------------------------------------------------------------

contains

 subroutine fv_init(ndt)
   use fms_io_mod, only : get_instance_filename

#ifdef USE_LIMA
   use pft_module, only: pft_init, fftfax
#endif
!#include "fv_arrays.h"

! SJL note:
! This routine initializes various parameters for FV and touch the arrays
! that should promote better memory placement on share memory machines (e.g., SGI)

   integer, intent(in)    :: ndt             ! model time step (seconds)

! LOCAL VARIABLES:
#ifdef USE_LIMA
   real rat, ycrit
#endif

   real:: dim0 = 180.           ! base dimension
   real:: dt0  = 1800.          ! base time step
   real:: ns0  = 4.             ! base nsplit for base dimension
   real dimx

   integer :: i, j, k, m, n
   integer :: unit, ierr, io
   integer :: kfirst, klast  ! useless (but required) input to y_decomp at this point
   integer, allocatable:: ysize(:)

! tracers
   integer :: num_family, index ! output of register_tracers
   integer, allocatable :: tracer_indices(:)
   logical :: is_in_nonadv_tracer_section
   character(len=128) :: msg, tname, restart_file_name
   integer :: log_unit, out_unit

   map_dt = ndt        ! mapping time step default to model time step

   call get_instance_filename('INPUT/fv_rst.res',restart_file_name )
   cold_start = (.not.file_exist('INPUT/fv_rst.res') &
         & .and. .not.file_exist('INPUT/fv_rst.res.nc')) &
         & .and. .not.file_exist(restart_file_name) &
         & .and. .not.file_exist('INPUT/'//fms_tracers_file)

!----- read namelist -----

#ifdef INTERNAL_FILE_NML
      read (input_nml_file,fv_core_nml,iostat=io)
      ierr = check_nml_error (io, 'fv_core_nml')
#else

    if (file_exist('input.nml')) then
          unit = open_namelist_file ( )
          ierr=1; do while (ierr /= 0)

! Every PE reading the same file? Ok on the SGI
             read (unit, nml=fv_core_nml, iostat=io, end=10)
                   ierr = check_nml_error (io, 'fv_core_nml')
          enddo
10    call close_file (unit)
    endif
#endif
!

! override number of tracers by reading field_table
    if (file_exist('field_table')) do_fms_tracer_manager = .TRUE.
    if (do_fms_tracer_manager) then !{


#ifndef MARS_GCM

#ifndef USE_LIMA
       if(full_phys .and. nt_phys /= 4) then !{
          call error_mesg('FV_INIT',  &
               & 'Wrong number of physics tracers nt_phys, should be 4! ', &
               FATAL)
       endif !}
#else
       if(nt_phys /= 4) then !{
          call error_mesg('FV_INIT',  &
               & 'Wrong number of physics tracers nt_phys, should be 4! ', &
               FATAL)
       endif !}
#endif

#endif 

       call tm_register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
       if(mpp_pe() == mpp_root_pe()) then
          out_unit = stdout()
          write(out_unit, *)'ncnst=', ncnst,' num_prog=',nt_prog,' pnats=',pnats,' num_family=',num_family
       endif


#ifndef MARS_GCM 

#ifndef USE_LIMA
       if(full_phys .and. nt_prog < nt_phys) then !{
          call error_mesg('FV_INIT',  &
               & 'Number of prognostic tracers in field_table should be >= number of physics tracers', &
               FATAL)
       endif !}
#else
       if(nt_prog < nt_phys) then !{
          call error_mesg('FV_INIT',  &
               & 'Number of prognostic tracers in field_table should be >= number of physics tracers', &
               FATAL)
       endif !}
#endif

#endif

       allocate(tracer_indices(ncnst))

       call tm_get_tracer_indices(MODEL_ATMOS, tracer_indices)

#ifndef USE_LIMA
       if(full_phys) then !{
#endif
 


#ifndef MARS_GCM
      ! Check tracer positions in field_table
       !
       ! sphum should the first entry
       ! followed by liq_wat, ice_wat, cld_amt (in any order)
       ! followed by all other advected tracers (if any)
       ! followed by non-advected tracers (if any)

          index = tm_get_tracer_index(MODEL_ATMOS, 'sphum', tracer_indices)
          if(index /= 1) then !{
             call error_mesg('FV_INIT',  &
                  & 'sphum should be the first tracer entry in field table ', &
                  FATAL)
          endif !}
          index = tm_get_tracer_index(MODEL_ATMOS, 'liq_wat', tracer_indices)
          if(index < 2 .or. index > 4) then !{
             call error_mesg('FV_INIT',  &
               & 'liq_wat should be tracer entry no 2 to 4 in field table ', &
               FATAL)
          endif !}
          index = tm_get_tracer_index(MODEL_ATMOS, 'ice_wat', tracer_indices)
          if(index < 2 .or. index > 4) then !{
             call error_mesg('FV_INIT',  &
               & 'ice_wat should be tracer entry no 2 to 4 in field table ', &
               FATAL)
          endif !}
          index = tm_get_tracer_index(MODEL_ATMOS, 'cld_amt', tracer_indices)
          if(index < 2 .or. index > 4) then !{
             call error_mesg('FV_INIT',  &
               & 'cld_amt should be tracer entry no 2 to 4 in field table ', &
               FATAL)
          endif
#endif


#ifndef USE_LIMA
       endif
#endif

       if(age_tracer) then !{
          ! age of air tracer should be present in field_table
          index = tm_get_tracer_index(MODEL_ATMOS, 'aoa', tracer_indices)
          if(index < 1 .or. index > ncnst) then !{
             call error_mesg('FV_INIT',  &
                  & 'age_tracer=.TRUE.: entry "aoa" is required in field table but could not find it', &
                  FATAL)
          endif !}
       endif !}

       ! check that non-advected tracers are defined last in field table

       is_in_nonadv_tracer_section = .FALSE.
       do i = 1, ncnst
          call tm_get_tracer_names(MODEL_ATMOS, tracer_indices(i), tname)
          if( is_in_nonadv_tracer_section .AND. &
               tm_check_if_prognostic(MODEL_ATMOS, tracer_indices(i)) ) then
             call tm_get_tracer_names(MODEL_ATMOS, tracer_indices(i), tname)
             write(msg, '(a,a,a)') 'Advected tracer ', tname,                    &
                   ' entry found after non-advected tracer entry in field table.'
             call error_mesg('FV_INIT', msg, FATAL)
          endif
          if( .NOT. tm_check_if_prognostic(MODEL_ATMOS, tracer_indices(i)) ) then
               is_in_nonadv_tracer_section = .TRUE.
          endif
       enddo

       deallocate(tracer_indices)

     endif !}

    call write_version_number ( version, tagname )
    if ( mpp_pe() == mpp_root_pe() ) then
       log_unit = stdlog()
       write (log_unit, nml=fv_core_nml)
    endif

#ifdef SPMD
      ng_d = max( 2, min(3, max(jord_mt, jord_tm, jord_tr)) )
      ng_s   = 3
      ighost = 0
      gid = -999 !the PEs not participating in mod_comm will show this value
      call fv_arrays_init( nlon, mlat, nlev, ncnst, ng_d, ng_s, ighost, &
           fv_domain )
      master = (gid == 0)
!set some variables from SPMD correctly...
      call mpp_get_compute_domain( fv_domain, beglon, endlon, beglat, endlat )
!gid.NE.mpp_pe when npex.NE.1!
      call y_decomp( mlat, nlev, beglat, endlat, kfirst, klast, gid )
#else
      beglat = 1; endlat = mlat
      beglon = 1; endlon = nlon
      ng_d   = 0;   ng_s = 0
      ighost = 0
      master = .true.
      if( n_zonal == 0 ) n_zonal = 1
#endif

      if(master) then
         write(6,*) 'SMP pseudo domain decomposition in E-W =', n_zonal
         write(6,*) 'Cold_start=', cold_start
      endif

      call timing_init()
      call timing_on('Total')
!     call timing_on('FV_init')

!----------------------------------
! Determine nsplit for the dynamics
!----------------------------------

      if ( n_split == 0 ) then
           dimx   = max ( nlon, 2*(mlat-1) )
          n_split = nint ( ns0*abs(ndt)*dimx/(dt0*dim0) + 0.49 )
          n_split = max ( 1, n_split )
          if(master) write(6,*) 'n_split is set to ', n_split, ' for (im,jm,dt)=',nlon,mlat,ndt
      else
          if(master) write(6,*) 'Using n_split from the namelist: ', n_split
      endif

     if(master) write(6,*) 'Total number of sponge layers= ', n_spong

     if ( adiabatic ) full_phys = .false.

! Set up some constants:
     if ( full_phys ) then
          cp_vir = (cv_wv + rvgas) / cp_air - 1.
     else
          cp_vir = 0.
     endif


#ifdef MARS_GCM
          cp_vir = 0.0
#endif 

!--------------------------------------------
! Set up horizontal grid geometrical factors:
!--------------------------------------------

      call set_fv_geom

#ifndef USE_LIMA
      call pft_init(nlon, mlat, beglat, endlat, ighost, cosp, cose, master)
#else
      ipft1 = max(2, beglat-ighost)
      js2g0 = max(2, beglat)

      ipft2 = min(mlat-1, endlat+ighost)
      jn1g1 = min(mlat,   endlat+1)

      allocate( sc(ipft1:ipft2),            se(js2g0:jn1g1)  )
      allocate( dc(nlon,ipft1:ipft2),  de(nlon,js2g0:jn1g1) )
      allocate( trigs(3*nlon/2+1) )

      call fftfax(nlon, ifax, trigs)

!---------------------------------------------
! Determine ycrit such that effective DX >= DY
!---------------------------------------------
      rat = float(nlon)/float(2*(mlat-1))
      ycrit = acos( min(0.81, rat) ) * (180./pi)

      if ( master ) write(6,*) 'Critical latitude for pft=', ycrit
      write (6,*) 'nlon, mlat, js2g0, jn1g1 ', nlon, mlat, js2g0, jn1g1
      write (6,*)  ' size(sc,1), size(dc,1), size(dc,2) ', size(sc,1), size(dc,1), size(dc,2)

      call pft_init(nlon, mlat, js2g0, jn1g1,      &
                    sc,   se,   dc,    de,         &
                    cosp, cose, ycrit, ipft1, ipft2)

#endif

!
! Vertical reference Eulerian coordinate for remapping:
!
      ks   = -1000    ! something large enough to cause trouble if not initialized
      ptop = too_big

      allocate ( ak(nlev+1) )
      allocate ( bk(nlev+1) )

      call initialize_arrays
      call fv_print_chksums( 'Exiting  fv_init' )

 end subroutine fv_init

 subroutine initialize_arrays
#include "fv_arrays.h"

!   real, parameter :: too_big=huge(1.)

   integer :: i, j, k, m
#include "fv_point.inc"
!---------------------------
! Initialize prognostic vars:
!---------------------------
! Note: the use of omp parallel is for better memory placement
!       (the "first touch" principle)

!$omp parallel do private (i, j, k)
   do k=ksp,kep !1,nlev
      do j=beglat-ng_d,endlat+ng_s
         do i=1,nlon
            u(i,j,k) = too_big
         enddo
      enddo

      do j=beglat-ng_d,endlat+ng_d
         do i=1,nlon
            v(i,j,k) = too_big
         enddo
      enddo

      do j=beglat-ng_d,endlat+ng_d
         do i=1,nlon
            pt(i,j,k) = too_big
         enddo
      enddo

      do j=beglat,endlat
         do i=1,nlon
           delp(i,j,k) = too_big
         enddo
      enddo

      do j=beglat,endlat
         do i=1,nlon
            pkz(i,j,k) = too_big
         enddo
      enddo
  enddo

!------------------
! Initialize tracers
!------------------

  if(do_fms_tracer_manager) then !{

     ! use tracer_manager profile spec
     do m = 1, ncnst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! should use 1D decomp here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call tm_set_tracer_profile (MODEL_ATMOS, m, q(:,:,:,m))
     enddo

  else !}{

     ! default

     do m=1,ncnst
        !$omp parallel do private (i, j, k)
!        do k=1,nlev
        do k=ksp,kep
           do j=beglat-ng_d, endlat+ng_d
              do i=1,nlon
                 q(i,j,k,m) = too_big
              enddo
           enddo
        enddo
     enddo

  endif !}

!$omp parallel do private (i, j, k)
!  do k=1,nlev+1
  do k=ksp,kep+1 !NOTE that kep+1 overlaps between PEs!
     do j=beglat,endlat
        do i=1,nlon
           pk(i,j,k) = too_big
        enddo
     enddo
     do i = 1,nlon
        pesouth(i,k) = too_big
     end do
  enddo
  call fv_array_sync()

! The follwoing var are mostly used by physics
!$omp parallel do private (i, j, k)
!  do j=beglat, endlat
   do j=jsp,jep
      do i=1,nlon
         ps(i,j) = too_big
      enddo

      do k=1,nlev
         do i=1,nlon
            omga(i,j,k) = too_big
         enddo
      enddo

      do k=1, nlev+1
         do i=1,nlon
            pe  (i,k,j) = too_big
            peln(i,k,j) = too_big
         enddo
      enddo
! required to solve a denorm problem for HS (pletzer)
      do i = 1,nlon
         ps_bp(i,j) = too_big
         phis (i,j) = too_big
         u_srf(i,j) = too_big
         v_srf(i,j) = too_big
      end do
   enddo
   call fv_array_sync()

 end subroutine initialize_arrays

!-----------------------------------------------------------------------

 subroutine get_eta_level(km, p_s, pf, ph, pscale)

  integer, intent(in) :: km
  real, intent(in)  :: p_s            ! unit: pascal
  real, intent(out) :: pf(km)
  real, intent(out) :: ph(km+1)
  real, intent(in), optional :: pscale

  integer k
  real    ak1

  if(nlev /= km)  &
  call error_mesg('get_eta_level:','dimensionally inconsistent', FATAL)

     ph(1) = ak (1)
  do k=2,nlev+1
     ph(k) = ak (k) + bk(k)*p_s
  enddo

  if ( present(pscale) ) then
      do k=1,nlev+1
         ph(k) = pscale*ph(k)
      enddo
  endif

  if( ak(1) > ptop_min ) then
     pf(1) = (ph(2) - ph(1)) / log(ph(2)/ph(1))
  else
     ak1 = (kappa + 1.) / kappa
     pf(1) = (ph(2) - ph(1)) / ak1
!    if(master) write(*,*) 'Modified p_full(1)=',pf(1), 0.5*(ph(1)+ph(2))
  endif

  do k=2,nlev
     pf(k) = (ph(k+1) - ph(k)) / log(ph(k+1)/ph(k))
  enddo

 end subroutine get_eta_level

!-----------------------------------------------------------------------

 subroutine set_fv_geom
! Sets the geometric factors used in FV core for the lat-lon grid
! This routine should only be called once per run

   integer i, j, imh
   real dp, dl
   real ph5, zamda, zam5
   real dg2r, asum, gmean

   allocate (  lon(nlon) )
   allocate (  lonb(nlon+1) )
   allocate ( rlonb(nlon+1) )

   allocate (   lat(mlat) )
   allocate (  latb(mlat+1) )
   allocate ( rlatb(mlat+1) )

   allocate ( rlat(nlon,beglat:endlat) )
   allocate ( rlon(nlon,beglat:endlat) )
   allocate ( area(nlon,beglat:endlat) )

   allocate ( cose(mlat) )
   allocate ( sine(mlat) )
   allocate ( cosp(mlat) )
   allocate ( sinp(mlat) )
   allocate (acosp(mlat) )
   allocate (acosu(mlat) )
   allocate ( grid_weight(mlat) )

   allocate ( coslon(nlon) )
   allocate ( sinlon(nlon) )

   allocate ( cosl5(nlon) )
   allocate ( sinl5(nlon) )

   allocate ( dx(mlat) )
   allocate ( one_by_dx(mlat) )
#ifndef USE_LIMA
   allocate ( f_d(nlon,beglat-ng_d:endlat+ng_d) )
#else
   allocate ( f_d(beglat-ng_d:endlat+ng_d) )
#endif


   pi = 4.d0 * atan(1.d0)
   dl = (pi+pi) / nlon
   dp = pi/(mlat-1)

         sine(1) = -1.
      do j=2,mlat
            ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(mlat-1))
         sine(j) = sin(ph5)
      enddo

      do j=1,mlat-1
         grid_weight(j) = sine(j+1) - sine(j)
      enddo
         grid_weight(mlat) = 1. - sine(mlat)

! Define area-conservative cosine at "cell center"
      do j=1,mlat
                cosp(j) = grid_weight(j) / dp
               acosp(j) = dp / grid_weight(j)
         grid_weight(j) = grid_weight(j) / (2*nlon)
      enddo

! Define cosine at edges..

         cose(1) = 0.
         cose(2) = 0.5 * cosp(2)
      do j=3,mlat-1
         cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      enddo
         cose(mlat) = cose(2)

      do j=2,mlat-1
         acosu(j) = 2. / (cose(j) + cose(j+1))
      enddo

      do j=1,mlat-1
         sinp(j) = 0.5 * (sine(j) + sine(j+1))
      enddo
         sinp(mlat) = 0.5 * (sine(mlat) + 1.)

      acap = nlon*(1.+sine(2))/dp               ! acap = nlon * cosp(1)
      rcap = 1.d0 / acap

! Coriolis parameter at cell center

#ifndef USE_LIMA
      do j = max(1,beglat-ng_d), min(mlat,endlat+ng_d)
         do i=1,nlon
            f_d(i,j) = (omega+omega) * sinp(j)
         enddo
      enddo
#else
      do j = max(1,beglat-ng_d), min(mlat,endlat+ng_d)
         f_d(j) = (omega+omega) * sinp(j)
      enddo
#endif

!-------------
! Bounding box
!-------------

      do j = 2, mlat
         latb(j) = -90. + (real(j-1)-0.5) * 180. / real(mlat-1)
      enddo

      latb( 1 ) = -90.
      latb( mlat+1 ) =  90.

      do j=1, mlat
         lat(j) = 0.5 * ( latb(j) + latb(j+1) )
      enddo

      one_by_dy = 1. / ( radius*dp )

      do j=2, mlat-1
                dx(j) = dl*radius*cosp(j)
         one_by_dx(j) = 1. / dx(j)
      enddo

            dg2r = pi/180.

      do j=1,mlat+1
         rlatb(j) = dg2r * latb(j)
      enddo

!-------------------------------
! Set longitude dependent param
!-------------------------------

      do i=1,nlon+1
#ifdef OLD_WAY
         lonb(i) = (real(i-1)-0.5) * 360. / real(nlon)
#else
!--------------------------
! FMS B grid way:
!--------------------------
         lonb(i) = real(i-1)*360./real(nlon)
#endif
      enddo

      do i=1,nlon
         lon(i) = 0.5 * ( lonb(i) + lonb(i+1) )
      enddo

      do i=1,nlon+1
         rlonb(i) = dg2r * lonb(i)
      enddo


        imh = nlon/2
      if(nlon /= 4*(imh/2)) then
         call error_mesg('FV_INIT','nlon must be dividable by 4', FATAL)
      endif

      do i=1,imh
#ifdef OLD_WAY
         zam5          = ((i-1)-0.5d0) * dl
         zamda         =  (i-1)*dl
#else
         zam5          = (i-1) * dl
         zamda         = (0.5d0+(i-1))*dl
#endif
         cosl5(i)      =  cos(zam5)
         cosl5(i+imh)  = -cosl5(i)
         sinl5(i)      =  sin(zam5)
         sinl5(i+imh)  = -sinl5(i)
         coslon(i)     =  cos(zamda)
         coslon(i+imh) = -coslon(i)
         sinlon(i)     =  sin(zamda)
         sinlon(i+imh) = -sinlon(i)
      enddo

      do j=beglat, endlat
         do i=1,nlon
            rlat(i,j) = dg2r * lat(j)
            rlon(i,j) = dg2r * lon(i)          ! Cubed sphere ready
            area(i,j) = (radius*dl*cosp(j)) * (dp*radius)
         enddo
      enddo

         asum = gmean(nlon, mlat, beglat, endlat, area)
      if(master) then
         write(6,*) 'Mean cell width (km)=', 1.e-3*sqrt(asum)
      endif

 end subroutine set_fv_geom

 subroutine fv_end(days, seconds)

integer, intent(in):: days
integer, intent(in):: seconds

      deallocate ( lon )
      deallocate ( lat )
      deallocate ( rlat )
      deallocate ( rlon )
      deallocate ( area )

      deallocate ( lonb )
      deallocate ( latb )

      deallocate ( rlonb )
      deallocate ( rlatb )

      deallocate ( cose )
      deallocate ( sine )
      deallocate ( cosp )
      deallocate ( sinp )
      deallocate ( grid_weight )

      deallocate ( acosp )
      deallocate ( acosu )

      deallocate ( coslon )
      deallocate ( sinlon )

      deallocate ( cosl5 )
      deallocate ( sinl5 )

      deallocate ( ak )
      deallocate ( bk )

      deallocate ( dx )
      deallocate ( one_by_dx )

      call fv_arrays_exit
! Clean up polar filter module:

#ifndef USE_LIMA
      call pft_end
#endif


      call timing_off('Total')

#ifdef SPMD
      call timing_prt(gid)
      call mp_exit( .false. )         !Call mp_exit to cleanly exit mod_comm
#else
      call timing_prt(0)
#endif
 end subroutine fv_end

 subroutine d2a3d(u,   v,   ua,   va,   im,   jm,   km,     &
                  jfirst, jlast, ng_d, ng_s,  coslon, sinlon)
! !ROUTINE: d2a3d -- Second order D-to-A grid transformation (3D)

! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km           ! Dimensions
      integer, intent(in):: jfirst, jlast        ! Latitude strip
      integer, intent(in):: ng_d, ng_s           ! as defined in fvgcm.F

      real, intent(in):: coslon(im),sinlon(im)    ! Sine and cosine in longitude
      real, intent(in):: v(im,jfirst-ng_d:jlast+ng_d,km)

      real, intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,km)

! !OUTPUT PARAMETERS:
      real, intent(out):: ua(im,jfirst:jlast,km)   ! U-Wind on A Grid
      real, intent(out):: va(im,jfirst:jlast,km)   ! V-Wind on A Grid

! local
      integer i, j, k
      integer imh, js2g0
      real un, vn, us, vs

      integer :: ks, ke

      imh = im/2
      js2g0 = max(jfirst, 2)

      call fv_array_check( LOC(u ) )
      call fv_array_check( LOC(v ) )
      call fv_array_check( LOC(ua) )
      call fv_array_check( LOC(va) )

#ifdef SPMD
      if ( .not. u_ghosted ) then
      call mp_send3d(gid-1, gid+1, im, jm, km,              &
                     1, im, jfirst-ng_d, jlast+ng_s, 1, km, &
                     1, im, jfirst, jfirst, 1, km, u)
      endif
#endif
      if( km.EQ.1 )then
!this means processing the last level km=nlev, surface
          ks=1; ke=1
      else if( km.EQ.nlev )then
          ks=ksp; ke=kep !shared arrays
      else
          call mpp_error( FATAL, 'D2A3D: shared pointers requires km.EQ.1 or km.EQ.nlev' )
      end if

!$omp parallel do private(i,j,k,us,vs,un,vn)
      do k=ks,ke !1,km
         do j=js2g0,jlast-1
            do i=1,im
               ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
            enddo
         enddo

         do j=js2g0,min(jlast,jm-1)
            do i=1,im-1
               va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(im,j,k) = 0.5*(v(im,j,k) + v(1,j,k))
         enddo

         if( jfirst == 1 ) then
! Projection at SP
             us = 0.
             vs = 0.
             do i=1,imh
                us = us + (ua(i+imh,2,k)-ua(i,2,k))*sinlon(i)      &
                     + (va(i,2,k)-va(i+imh,2,k))*coslon(i)
                vs = vs + (ua(i+imh,2,k)-ua(i,2,k))*coslon(i)      &
                     + (va(i+imh,2,k)-va(i,2,k))*sinlon(i)
             enddo
             us = us/im
             vs = vs/im
             do i=1,imh
                ua(i,1,k)   = -us*sinlon(i) - vs*coslon(i)
                va(i,1,k)   =  us*coslon(i) - vs*sinlon(i)
                ua(i+imh,1,k)   = -ua(i,1,k)
                va(i+imh,1,k)   = -va(i,1,k)
             enddo
         endif

         if ( jlast == jm ) then
! Projection at NP
             un = 0.
             vn = 0.
             do i=1,imh
                un = un + (ua(i+imh,jm-1,k)-ua(i,jm-1,k))*sinlon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*coslon(i)
                vn = vn + (ua(i,jm-1,k)-ua(i+imh,jm-1,k))*coslon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*sinlon(i)
             enddo

             un = un/im
             vn = vn/im
             do i=1,imh
                ua(i,jm,k) = -un*sinlon(i) + vn*coslon(i)
                va(i,jm,k) = -un*coslon(i) - vn*sinlon(i)
                ua(i+imh,jm,k) = -ua(i,jm,k)
                va(i+imh,jm,k) = -va(i,jm,k)
             enddo
         endif
      enddo

#ifdef SPMD
      if ( .not. u_ghosted ) then
      call mp_recv3d(gid+1, im, jm, km, &
                         1, im, jfirst-ng_d, jlast+ng_s, 1, km, &
                         1, im, jlast+1, jlast+1, 1, km, u)
           u_ghosted = .true.
      endif
      call fv_array_sync()

      if ( jlast < jm ) then
!balaji, switched parallelization to i, since km might be 1
!$omp parallel do private(i, k)
        do k=1,km
           do i=isp,iep
              ua(i,jlast,k) = 0.5 * ( u(i,jlast,k) + u(i,jlast+1,k) )
           enddo
        enddo
      endif
#endif
      call fv_array_sync()
    end subroutine d2a3d

 subroutine p_var(im, jm, km, jfirst, jlast, ptop,    &
                  delp,  ps,  pe, peln, pk, pkz,      &
                  cappa, q,   ng, nq, adjust_dry_mass )

! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
!

! Input:
   integer,  intent(in):: im, jm, km               ! Dimensions
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   integer,  intent(in):: nq
   integer,  intent(in):: ng
   logical,  intent(in):: adjust_dry_mass

   real, intent(in):: ptop
   real, intent(in):: cappa
   real, intent(inout):: delp(im,jfirst:jlast, km)
   real, intent(inout):: q(im,jfirst-ng:jlast+ng, km, nq)
! Output:
   real, intent(out) ::   ps(im, jfirst:jlast)
   real, intent(out) ::   pk(im, jfirst:jlast, km+1)
   real, intent(out) ::   pe(im, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) :: peln(im, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(im, jfirst:jlast, km)

! Local
   real pek
   real lnp
   real ak1
   integer i, j, k

   pek = ptop ** cappa
   ak1 = (cappa + 1.) / cappa

   call fv_array_check( LOC(delp) )
   call fv_array_check( LOC(q) )
   call fv_array_check( LOC(ps) )
   call fv_array_check( LOC(pk) )
   call fv_array_check( LOC(pe) )
   call fv_array_check( LOC(peln) )
   call fv_array_check( LOC(pkz) )

!$omp parallel do private(i,j,k)
   do j=jsp,jep !jfirst,jlast
      do i=1,im
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      do k=2,km+1
         do i=1,im
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = pe(i,k,j)**cappa
         enddo
      enddo

      do i=1,im
         ps(i,j) = pe(i,km+1,j)
      enddo

!---- GFDL modification
      if( ptop < ptop_min ) then
          do i=1,im
             peln(i,1,j) = peln(i,2,j) - ak1
          enddo
      else
             lnp = log( ptop )
          do i=1,im
             peln(i,1,j) = lnp
          enddo
      endif
!---- GFDL modification

      do k=1,km
        do i=1,im
           pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(cappa*(peln(i,k+1,j)-peln(i,k,j)))
        enddo
      enddo
   enddo
   call fv_array_sync()

! Check dry air mass:
   call drymadj( im,  jm,  km,  jfirst,  jlast, ng, mountain, cappa,   &
                 ptop,   ps, delp, pe, pk, peln, pkz, nq, q, adjust_dry_mass )

 end subroutine p_var


 subroutine polavg(p, im, jm, jfirst, jlast)
 integer, intent(in):: im, jm, jfirst, jlast
 real, intent(inout):: p(im,jfirst:jlast)

! local:
   integer i
   real sum1

   if ( jfirst == 1 ) then
           sum1 = 0.
        do i=1,im
           sum1 = sum1 + p(i,1)
        enddo
           sum1 = sum1/im

        do i=1,im
           p(i,1) = sum1
        enddo
   endif

   if ( jlast == jm ) then
           sum1 = 0.
        do i=1,im
           sum1 = sum1 + p(i,jm)
        enddo
           sum1 = sum1/im

        do i=1,im
           p(i,jm) = sum1
        enddo
   endif

 end subroutine polavg


 subroutine drymadj( im,  jm,  km,  jfirst,  jlast,  ng, &
                     moun, cappa,   ptop, ps, delp, pe,  &
                     pk, peln, pkz, nq, q, adjust_dry_mass )

! !INPUT PARAMETERS:
      integer im, jm, km     ! Dimensions
      integer jfirst, jlast  ! Latitude strip
      integer nq             ! Number of tracers
      integer ng
      logical moun
      logical, intent(in):: adjust_dry_mass
      real, intent(in):: ptop
      real, intent(in):: cappa

! !INPUT/OUTPUT PARAMETERS:
      real  delp(im,jfirst:jlast,km)     !
      real  pe(im,km+1,jfirst:jlast)     !
      real   q(im,jfirst-ng:jlast+ng,km,nq)
      real  ps(im,jfirst:jlast)        ! surface pressure

      real, intent(inout)::   pk(im, jfirst:jlast, km+1)
      real, intent(inout):: peln(im, km+1, jfirst:jlast)    ! Edge pressure
      real, intent(inout)::  pkz(im, jfirst:jlast, km)

#ifdef USE_LIMA
      real drym

!     parameter ( drym = 98222. )       ! setting for US NAVY 10 min data
      parameter ( drym = 98288. )       ! setting for USGS
!     parameter ( drym = 98290. )       ! New setting for USGS
#endif

      integer   i, j, k
      real   psmo
      real   psdry, qtot
      real   gmean
      real   dpd
      integer ic, ip
! Local
      real  psq(im,jfirst:jlast)
      real  psd(im,jfirst:jlast)     ! surface pressure  due to dry air mass
#ifdef use_shared_pointers
      pointer( p_psq, psq )
      pointer( p_psd, psd )
      call fv_stack_push( p_psq, im*(jlast-jfirst+1) )
      call fv_stack_push( p_psd, im*(jlast-jfirst+1) )
#endif
      call fv_array_check( LOC(delp) )
      call fv_array_check( LOC(pe) )
      call fv_array_check( LOC(q) )
      call fv_array_check( LOC(ps) )
      call fv_array_check( LOC(pk) )
      call fv_array_check( LOC(peln) )
      call fv_array_check( LOC(pkz) )

      ip = min(3,size(q,4))

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)
      do 1000 j=jsp,jep !jfirst,jlast
         do i=1,im
            psd(i,j) = ptop
            psq(i,j) = 0.
         enddo

#ifdef MARS_GCM
         do k=1,km
            do i=1,im
               psd(i,j) = psd(i,j) +  delp(i,j,k)
            enddo
         enddo
#else

         if( full_phys ) then
             do k=1,km
                do i=1,im
                   psd(i,j) = psd(i,j) + delp(i,j,k)*(1.-q(i,j,k,1))
!             psq(i,j) = psq(i,j) + delp(i,j,k)*(q(i,j,k,1)+q(i,j,k,2)+q(i,j,k,3))
                   psq(i,j) = psq(i,j) + delp(i,j,k)*sum( q(i,j,k,1:ip) )
                enddo
             enddo
         else
             do k=1,km
                do i=1,im
                   psd(i,j) = psd(i,j) +  delp(i,j,k)
                enddo
             enddo
         endif
#endif MARS_GCM

1000  continue
      call fv_array_sync() !probably needed before gmean

! Check global maximum/minimum
      psmo  = gmean(im, jm, jfirst, jlast, ps (1,jfirst))
      psdry = gmean(im, jm, jfirst, jlast, psd(1,jfirst))
       qtot = gmean(im, jm, jfirst, jlast, psq(1,jfirst))

      if(master) then
         write(6,*) 'Total surface pressure (mb) = ', 0.01*psmo
         if ( full_phys ) then
         write(6,*) 'mean dry surface pressure = ', 0.01*psdry
         write(6,*) 'TPW-vapor (kg/m**2) =', (psmo-psdry)/GRAV
         write(6,*) 'TPW-total (kg/m**2) =', qtot/GRAV
         endif
      endif

      if( .not. adjust_dry_mass ) return

#ifndef USE_LIMA
      if(moun) then
         dpd = dry_mass - psdry
      else
         dpd = 1000.*100. - psdry
      endif
#else
      if(moun) then
         dpd = drym - psdry
      else
         dpd = 1000.*100. - psdry
      endif
#endif

      if(master) write(6,*) 'dry mass to be added (pascals) =', dpd

!$omp  parallel do            &
!$omp  default(shared)        &
!$omp  private(i,j, ic)

      do 2000 j=jsp,jep !jfirst,jlast
         do ic=1,nq
            do i=1,im
               q(i,j,km,ic) = q(i,j,km,ic)*delp(i,j,km) / (delp(i,j,km)+dpd)
            enddo
         enddo

! Adjust the lowest Lagrangian layer
         do i=1,im
            delp(i,j,km) = delp(i,j,km) + dpd
            pe(i,km+1,j) = pe(i,km,j) + delp(i,j,km)
            ps(i,j) = pe(i,km+1,j)
! adjust pk, peln, pkz
            pk(i,j,km+1) = pe(i,km+1,j) ** cappa
            peln(i,km+1,j) = log(pe(i,km+1,j))
            pkz(i,j,km) = (pk(i,j,km+1)-pk(i,j,km)) /      &
                 (cappa*(peln(i,km+1,j)-peln(i,km,j)))
         enddo
2000  continue
      call fv_array_sync() !probably needed before gmean

      psmo = gmean(im, jm, jfirst, jlast, ps(1,jfirst) )
      if( master ) write(6,*)                                     &
              'Total (moist) surface pressure after adjustment = ', &
               0.01*psmo

 end subroutine drymadj


 subroutine tracer_mass(im, jm, km, jfirst, jlast, ng, nq, q, delp)


! !INPUT PARAMETERS:
      integer im, jm, km     ! Dimensions
      integer jfirst, jlast  ! Latitude strip
      integer nq             ! Number of tracers
      integer ng

      real  delp(im,jfirst:jlast,km)     !
      real   q(im,jfirst-ng:jlast+ng,km,nq)

! Local
      real   gmean
      real  qz(im,jfirst:jlast,nq)
      real  qmass
      integer   i, j, k
      integer iq
#ifdef use_shared_pointers
      pointer( p_qz, qz )
      call fv_stack_push( p_qz, im*(jlast-jfirst+1)*nq )
#endif

      call fv_array_check( LOC(delp) )
      call fv_array_check( LOC(q) )
      do iq=1,nq

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)
      do 1000 j=jsp,jep
        do i=1,im
           qz(i,j,iq) = 0.
        enddo

        do k=1,km
           do i=1,im
              qz(i,j,iq) = qz(i,j,iq) + delp(i,j,k)*q(i,j,k,iq)
           enddo
        enddo
1000    continue
        call fv_array_sync()

      qmass = gmean(im, jm, jfirst, jlast, qz(1,jfirst,iq))

      if(master) write(6,*) iq, 'tracer mass = ', qmass

      enddo

 end subroutine tracer_mass
end module fv_pack

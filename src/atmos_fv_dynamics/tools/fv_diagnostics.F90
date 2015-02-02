module fv_diagnostics

! Programmer: sjl, Oct 30, 2003

 use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                             register_static_field, send_data
 use constants_mod,    only: GRAV, cp_air, rdgas, kappa, WTMAIR, WTMCO2
 use time_manager_mod, only: time_type, get_date, get_time
 use fms_mod,          only: error_mesg, FATAL, stdlog, write_version_number
 use mpp_domains_mod,  only: mpp_define_domains,  mpp_domains_init, domain2D

! use fv_pack,          only: nlon, mlat, nlev, beglat, endlat, u, v, pt, delp, q,    &
!                             omga, peln, pe, ua, va, phis, ptop, ak, bk, ks, ncnst,  &
!                             ng_d, ng_s, ps, pk, pkz, lon, lat, lonb, latb,          &
!                             age_tracer, age_time, fv_domain,  print_freq,           &
!                             drymadj, tracer_mass,                              &
!                             master, get_eta_level, nt_phys, do_fms_tracer_manager
 use fv_pack,          only: nlon, mlat, beglon, endlon, beglat, endlat, &
      ptop, ak, bk, ks, ng_d, ng_s, lon, lat, lonb, latb,          &
      age_tracer, age_time, fv_domain,  print_freq,  drymadj, tracer_mass, &
      master, get_eta_level, nt_phys, do_fms_tracer_manager, area
#ifdef MARS_GCM
 use fv_pack,       only:   p_ref
 use fv_phys_mod,   only:   mars_mass_budget
#endif MARS_GCM

! --- tracer manager ---
 use tracer_manager_mod, only : get_tracer_names, get_number_tracers
 use field_manager_mod, only  : MODEL_ATMOS      
 use mpp_mod, only: mpp_sync
 use fv_arrays_mod, only: fv_stack_push, fv_array_sync, fv_print_chksums
 use pmaxmin_mod, only: pmaxmin, pmaxming


 implicit none

 integer, save ::  id_ps,   id_ua,   id_va,   id_wa,  id_ta,   id_pt,   &
                   id_divg, id_vort, id_hght, id_pv,  id_usus, id_tsts, &
                   id_vsts, id_usvs, id_h500, id_slp, id_p5km, id_aoa,  &
                   id_dnflux, id_upflux, id_area
 integer, allocatable :: id_tracer(:), id_tracer_tend(:)
! 2008/04/09  jgj: add tracer dry mass and volume mixing ratios
 integer, allocatable :: id_tracer_dmmr(:)
 integer, allocatable :: id_tracer_dvmr(:)


 real       :: missing_value = -1.e10
 real, save :: ginv
 logical used
 type(time_type) :: fv_time
 real, allocatable, save :: phalf(:)
 save fv_time

 !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: fv_diagnostics.F90,v 17.0 2009/07/21 02:53:23 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

contains

  subroutine fv_diag_init(axes, Time)
#include "fv_arrays.h" 
! FV dynamics specific diagnostics

    type(time_type), intent(in) :: Time
    integer :: id_lonb, id_lon, id_latb, id_lat, id_phalf, id_pfull
    integer :: id_bk, id_pk, id_zsurf
    character(len=8) :: mod_name = 'dynamics'
    integer, intent(out) :: axes(4)
    logical used
    integer i, j
    integer :: log_unit
    real zsurf(nlon,beglat:endlat)
    real pfull(nlev)
    real qmin, qmax

    real :: vrange(2), trange(2), arange(2), slprange(2)

! tracers
    character(len=128)   :: tname
    character(len=256)   :: tlongname, tunits
    integer              :: ntprog
#include "fv_point.inc"

    log_unit = stdlog()

    call fv_print_chksums( 'Entering  fv_diag_init' )
! valid range for some fields

    vrange = (/ -330.,  350. /)  ! winds
#ifdef MARS_GCM
    trange = (/  50.,  350. /)  ! temperature
#else
    trange = (/  100.,  350. /)  ! temperature
#endif 
    arange = (/    0.,   20. /)  ! age-of-air
    slprange = (/800.,  1200./)  ! sea-level-pressure

    ginv = 1./GRAV
    fv_time = Time

    allocate ( phalf(nlev+1) )

#ifdef MARS_GCM
    call get_eta_level(nlev, p_ref, pfull, phalf, 0.01) 
#else
    call get_eta_level(nlev, 1.E5, pfull, phalf, 0.01) 
#endif 

    id_lonb=diag_axis_init('lonb', lonb, 'degrees_E', 'x', 'longitude edges', &
         set_name=mod_name , Domain2=fv_domain )
    id_lon =diag_axis_init('lon',  lon,  'degrees_E', 'x', 'longitude',       &
         set_name=mod_name, edges=id_lonb, Domain2=fv_domain )

    id_latb=diag_axis_init('latb', latb, 'degrees_N', 'y', 'latitude edges',  &
         set_name=mod_name, Domain2=fv_domain )
    id_lat =diag_axis_init('lat',  lat,  'degrees_N', 'y', 'latitude',        &
         set_name=mod_name, edges=id_latb, Domain2=fv_domain)

    id_phalf = diag_axis_init('phalf', phalf, 'mb', 'z', &
         'ref half pressure level', direction=-1, set_name=mod_name )
    id_pfull = diag_axis_init('pfull', pfull, 'mb', 'z', &
         'ref full pressure level', direction=-1, set_name=mod_name, &
         edges=id_phalf)

    axes(1) = id_lon
    axes(2) = id_lat
    axes(3) = id_pfull
    axes(4) = id_phalf

!---- register static fields -------

    id_bk    = register_static_field ( mod_name, 'bk', (/id_phalf/), &
         'vertical coordinate sigma value', 'none' )

    id_pk    = register_static_field ( mod_name, 'pk', (/id_phalf/), &
         'pressure part of the hybrid coordinate', 'pascal' )

    id_area =  register_static_field ( mod_name, 'area', axes(1:2), &
         & 'cell area elements', 'm^2')
    if( id_area > 0 ) &
         &  used = send_data( id_area, area(beglon:endlon,beglat:endlat), Time )

    do j=beglat,endlat
       do i=1,nlon
          zsurf(i,j) = ginv * phis(i,j) 
       enddo
    enddo

    id_zsurf = register_static_field ( mod_name, 'zsurf', axes(1:2),   &
         'surface height', 'm' )

!--- Send static data

    if ( id_bk > 0 )    used = send_data ( id_bk,    bk, Time )
    if ( id_pk > 0 )    used = send_data ( id_pk,    ak, Time )
    if ( id_zsurf > 0 ) &
         used = send_data( id_zsurf, zsurf(beglon:endlon,beglat:endlat), Time )

!     if ( master ) then
!     do i=1,nlev+1
!        write(6,*) i, ak(i), bk(i)
!     enddo
!     endif
!--------------------------------------------------------------
! Register main prognostic fields: ps, (u,v), t, omega (dp/dt)
!--------------------------------------------------------------

    id_ps = register_diag_field ( mod_name, 'ps', axes(1:2), Time,        &
         'surface pressure', 'Pa', missing_value=missing_value )
    id_ua = register_diag_field ( mod_name, 'ucomp', axes(1:3), Time,        &
         'zonal wind', 'm/sec', missing_value=missing_value, range=vrange )
    id_va = register_diag_field ( mod_name, 'vcomp', axes(1:3), Time,        &
         'meridional wind', 'm/sec', missing_value=missing_value, range=vrange)
    id_wa = register_diag_field ( mod_name, 'omega', axes(1:3), Time,        &
         'omega', 'pa/sec', missing_value=missing_value )
    id_dnflux = register_diag_field (mod_name, 'dnflux', axes(1:3), Time, &
         'downward mass flux', 'pa/sec', missing_value=missing_value )
    id_upflux = register_diag_field (mod_name, 'upflux', axes(1:3), Time, &
         'upward mass flux', 'pa/sec', missing_value=missing_value )

    id_ta = register_diag_field ( mod_name, 'temp', axes(1:3), Time,        &
         'temperature', 'deg_K', missing_value=missing_value, range=trange )

    id_pt = register_diag_field ( mod_name, 'pt_dry', axes(1:3), Time,        &
         'dry poten temperature', 'deg_K', missing_value=missing_value )

!------------------
! Register tracers
!------------------

    allocate(id_tracer(max(ncnst, nt_phys)))
    id_tracer(:) = 0

! 2008/04/09  jgj: add tracer dry mass and volume mixing ratios
    allocate(id_tracer_dmmr(max(ncnst, nt_phys)))
    id_tracer_dmmr(:) = 0
    allocate(id_tracer_dvmr(max(ncnst, nt_phys)))
    id_tracer_dvmr(:) = 0

    if(.NOT. do_fms_tracer_manager) then !{

! no field table 

        id_tracer(1) = register_diag_field ( mod_name, 'sphum', axes(1:3), &
             Time, 'sphum', 'mass/air_mass', &
             missing_value=missing_value )
        id_tracer(2) = register_diag_field ( mod_name, 'liq_wat', axes(1:3), &
             Time, 'cloud liquid water', 'mass/air_mass', &
             missing_value=missing_value )
        id_tracer(3) = register_diag_field ( mod_name, 'ice_wat', axes(1:3), &
             Time, 'cloud ice water', 'mass/air_mass',&
             missing_value=missing_value )
        id_tracer(4) = register_diag_field ( mod_name, 'cld_amt', axes(1:3), &
             Time, 'cloud craction', '%', &
             missing_value=missing_value )

    else !}{

! do_fms_traces == .TRUE.

#ifdef MARS_GCM
        do i= 1, ncnst
#else
        do i = 1, nt_phys
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           id_tracer(i) = register_diag_field ( mod_name, trim(tname),  &
                axes(1:3), Time, trim(tlongname), &
                trim(tunits), missing_value=missing_value)
           if (master) then
               if (id_tracer(i) > 0) then
                   write(log_unit,'(a,a,a,a)') &
                        & 'Diagnostics available for tracer ',tname, &
                        ' in module ', mod_name
               end if
           endif
        enddo
! for higher tracer number we take the default missing values from FMS
        do i = nt_phys+1, ncnst
#endif MARS_GCM
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           id_tracer(i) = register_diag_field ( mod_name, trim(tname),  &
                axes(1:3), Time, trim(tlongname),           &
                trim(tunits), missing_value=missing_value)
           if (master) then
               if (id_tracer(i) > 0) then
                   write(log_unit,'(a,a,a,a)') &
                        & 'Diagnostics available for tracer ',trim(tname), &
                        ' in module ', trim(mod_name)
               end if
           endif
        enddo

!-------
#ifdef MARS_GCM

#else
! jgj: add co2 tracer dry mass and volume mixing ratios
        do i = nt_phys+1, ncnst
           call get_tracer_names ( MODEL_ATMOS, i, tname, tlongname, tunits )
           if (trim(tname).eq.'co2') then
             id_tracer_dmmr(i) = register_diag_field ( mod_name, trim(tname)//'_dmmr',  &
                axes(1:3), Time, trim(tlongname)//" (dry mmr)",           &
                trim(tunits), missing_value=missing_value)
             id_tracer_dvmr(i) = register_diag_field ( mod_name, trim(tname)//'_dvmr',  &
                axes(1:3), Time, trim(tlongname)//" (dry vmr)",           &
                'mol/mol', missing_value=missing_value)
             if (master) then
                if (id_tracer_dmmr(i) > 0) then
                   write(log_unit,'(a,a,a,a)') &
                        & 'Diagnostics available for co2 dry mmr ',trim(tname)//'_dmmr', &
                        ' in module ', trim(mod_name)
                end if
                if (id_tracer_dvmr(i) > 0) then
                   write(log_unit,'(a,a,a,a)') &
                        & 'Diagnostics available for co2 dry vmr ',trim(tname)//'_dvmr', &
                        ' in module ', trim(mod_name)
                end if
             endif
           endif
        enddo
#endif 
!-------

    endif !}

!----------------------
! Divergence on C grid
!----------------------
! Computed inside fvcore/trac2d()
    id_divg = register_diag_field ( mod_name, 'divg', axes(1:3), Time,        &
         'divergence', '1/sec', missing_value=missing_value )

!----------------------
! Vorticity on D grid
!----------------------
    id_vort = register_diag_field ( mod_name, 'vort', axes(1:3), Time,        &
         'vorticity', '1/sec', missing_value=missing_value )

!--------------------------
! Ertel Potential Vorticity
!--------------------------
    id_pv = register_diag_field ( mod_name, 'pv', axes(1:3), Time,        &
         'potential vorticity', '1/sec', missing_value=missing_value )

!-------------------
! Age_of_air tracer
!-------------------
    id_aoa = register_diag_field ( mod_name, 'aoa', axes(1:3), Time,    &
         'age_of_air tracer', 'sec', missing_value=missing_value,   &
         range=arange )

!---------------------------
! Height at model half level
!---------------------------
    id_hght = register_diag_field ( mod_name, 'hght', &
         (/id_lon,id_lat,id_phalf/), &
         Time, 'height', 'm', missing_value=missing_value )

!--------------
! 500 mb Height
!--------------
    id_h500 = register_diag_field ( mod_name, 'h500', axes(1:2),  Time,   &
         '500-mb hght', 'm', missing_value=missing_value )
!-------------------
! Sea-level-pressure
!-------------------
    id_slp = register_diag_field ( mod_name, 'slp', axes(1:2),  Time,   &
         'sea-level pressure', 'mb', missing_value=missing_value,  &
         range=slprange )

!-------------------
! 5-km pressure
!-------------------
! This is analogous to 500-mb height. This shows the actual highs and lows
! at the 5-km physical height.
    id_p5km = register_diag_field ( mod_name, 'p5km', axes(1:2),  Time,   &
         'pressure at 5 km', 'mb', missing_value=missing_value )
!------------------
! 2nd moment stats:
!------------------

! u-wind variance
    id_usus = register_diag_field ( mod_name, 'usus', axes(1:3), Time, &
         'u-wind variance', '(m/sec)**2', missing_value=missing_value )

! temperature variance
    id_tsts = register_diag_field ( mod_name, 'tsts', axes(1:3), Time, &
         'temperature variance', 'K**2', missing_value=missing_value )

! eddy momentum transport
    id_usvs = register_diag_field ( mod_name, 'usvs', axes(1:3), Time, &
         'eddy momentum transport', '(m/sec)**2', missing_value=missing_value )

! eddy heat transport
    id_vsts = register_diag_field ( mod_name, 'vsts', axes(1:3), Time, &
         'eddy heat transport', '(m/sec)*K', missing_value=missing_value )

    ! make sure all send_data completed before one thread exits the routine
    call fv_array_sync()

    call fv_print_chksums( 'Exiting  fv_diag_init' )
  end subroutine fv_diag_init

!=================================================


 subroutine fv_diag( Time, im, jm, km, jfirst, jlast, nq, zvir, dt_atmos, hs_phys)

   use pv_module,   only: vort_d, pv_entropy
#include "fv_arrays.h"
!INPUT PARAMETERS:
      type(time_type), intent(in) :: Time
      integer, intent(in):: im         ! dimension in east-west
      integer, intent(in):: jm         ! dimension in North-South
      integer, intent(in):: km         ! number of Lagrangian layers
      integer, intent(in):: jfirst     ! starting latitude index for MPI
      integer, intent(in):: jlast      ! ending latitude index for MPI
      integer, intent(in):: nq         ! number of tracers
      real, intent(in):: zvir
      real, intent(in):: dt_atmos     ! model time step in seconds
      logical, intent(in):: hs_phys

      integer, parameter:: kmax = 100
      real height(kmax)
      real  log_p(kmax)
      real qmax, qmin, fac

      real gg
      real gmean
      real ztop
      integer  i, j, k, ic
      integer  yr, mon, dd, hr, mn, days, seconds
      logical::  prt_minmax
      real p0
! Local:
      real :: a2(im,jfirst:jlast)     ! work array for x-y diag
      real :: ax(jfirst:jlast,km)     ! work array for zonal means
!      real, allocatable:: xx(:,:,:)   ! work array for 2nd moments
!      real, allocatable:: us(:,:,:)
!      real, allocatable:: vs(:,:,:)
!      real, allocatable:: wz(:,:,:)
      real :: xx(im,jfirst:jlast,km)
      real :: us(im,jfirst:jlast,km)
      real :: vs(im,jfirst:jlast,km)
      real :: wz(im,jfirst:jlast,km+1)

! 2008/04/09  jgj: add tracer dry vmr to compare to obs
      real :: co2_dmmr(im,jfirst:jlast,km)
      real :: co2_dvmr(im,jfirst:jlast,km)
! tracers
    character(len=128)   :: tname
    character(len=256)   :: tlongname, tunits

#ifdef use_shared_pointers
      pointer( p_a2, a2 )
      pointer( p_ax, ax )
      pointer( p_xx, xx )
      pointer( p_us, us )
      pointer( p_vs, vs )
      pointer( p_wz, wz )
!jgj
      pointer( p_co2_dmmr, co2_dmmr )
      pointer( p_co2_dvmr, co2_dvmr )
#include "fv_point.inc"

      call fv_stack_push( p_a2, im*(jlast-jfirst+1)        )
      call fv_stack_push( p_ax, (jlast-jfirst+1)*km        )
      call fv_stack_push( p_xx, im*(jlast-jfirst+1)*km     )
      call fv_stack_push( p_us, im*(jlast-jfirst+1)*km     )
      call fv_stack_push( p_vs, im*(jlast-jfirst+1)*km     )
      call fv_stack_push( p_wz, im*(jlast-jfirst+1)*(km+1) )
!jgj
      call fv_stack_push( p_co2_dmmr, im*(jlast-jfirst+1)*km     )
      call fv_stack_push( p_co2_dvmr, im*(jlast-jfirst+1)*km     )
#endif
      fv_time = Time

! fv_time used within fvcore for diagnostics purpose
      if ( hs_phys ) then
         call get_time (fv_time, seconds,  days)
         if( print_freq == 0 ) then
                 prt_minmax = .false.
         elseif( print_freq < 0 ) then
                 prt_minmax = .true.
         else
                 prt_minmax = mod(seconds, 3600*print_freq) == 0
         endif
      else
         call get_date(fv_time, yr, mon, dd, hr, mn, seconds)
         if( print_freq == 0 ) then
                 prt_minmax = .false.
         elseif( print_freq < 0 ) then
                 prt_minmax = .true.
         else
                 prt_minmax = mod(hr, print_freq) == 0 .and. mn==0 .and. seconds==0
         endif
      endif


      if( prt_minmax ) then
        if ( master ) then
             if ( hs_phys ) then
             write(*,*) Days, seconds
             else
             write(*,*) yr, mon, dd, hr, mn, seconds
             endif
        endif
             call pmaxmin('PS', ps, qmin, qmax, im*(jlast-jfirst+1),  1, 0.01)

#ifdef MARS_GCM
       call mars_mass_budget( im, jm, km, jfirst, jlast, ng_d, ps, delp, nq, q, mars_sfc_budg ) 

#else
! Check dry air mass:
             call drymadj(im, jm, km, jfirst, jlast, ng_d, .true., kappa,   &
                          ptop,  ps, delp, pe, pk, peln, pkz, nq,           &
                          q, .false.)  
#endif MARS_GCM
             call pmaxmin('U', ua, qmin, qmax, im*(jlast-jfirst+1), km, 1.)
             call pmaxmin('V', va, qmin, qmax, im*(jlast-jfirst+1), km, 1.)
    endif

    if (id_ps > 0) used = send_data ( id_ps, ps(beglon:endlon,beglat:endlat), Time )

    if ( id_ua > 0 .or. id_va > 0 ) then

         if (id_ua > 0) used = send_data ( id_ua, ua(beglon:endlon,beglat:endlat,:), Time )
         if (id_va > 0) used = send_data ( id_va, va(beglon:endlon,beglat:endlat,:), Time )

    endif

! Send vertical velocity (pa/s)
    if (id_wa > 0)  used = send_data ( id_wa, omga(beglon:endlon,beglat:endlat,:), Time )
    ! The following sync call is required to ensure that send_data completes before 
    ! array omga is touched by another thread. 
    if( id_vort>0 .or. id_pv>0 .or. id_ta>0 .or. id_pt>0 ) call fv_array_sync()
    if( prt_minmax .and. km>1 ) then
        call pmaxmin('W', omga, qmin, qmax, im*(jlast-jfirst+1), km, 1.)
    endif

!    if ( id_usus>0 .or. id_usvs>0 .or. id_vsts>0 .or. id_tsts>0 ) then
!           allocate ( xx(im,jfirst:jlast,km) )
!           allocate ( us(im,jfirst:jlast,km) )
!           allocate ( vs(im,jfirst:jlast,km) )
!    endif

!----------------------------------------------------------------
! Note: wind variances computed from A grid winds are weaker than 
!       native D grid due to two-grid averaging (from D to A)
!----------------------------------------------------------------

    if ( id_usus > 0 .or. id_usvs > 0 ) then
!$omp parallel do default(shared)                     &
!$omp private(i, j, k)
      do k=ksp,kep
         call zsmean( ua(1,jfirst,k), im, jfirst, jlast,        &
                      ax(jfirst,k),   us(1,jfirst,k), xx(1,jfirst,k) )
      enddo
      call fv_array_sync()
      if ( id_usus > 0 ) used = send_data ( id_usus, xx(beglon:endlon,beglat:endlat,:), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array xx is touched by another thread. 
      if( id_usvs>0 .or. id_tsts>0 .or. id_vsts>0 ) call fv_array_sync()
    endif

    if ( id_usvs > 0 .or. id_vsts > 0 ) then
!$omp parallel do default(shared)                     &
!$omp private(i, j, k)
      do k=ksp,kep
         call zsmean( va(1,jfirst,k), im, jfirst, jlast,        &
                      ax(jfirst,k),   vs(1,jfirst,k)   )
         if ( id_usvs > 0 ) then
         do j=jfirst,jlast
            do i=1,im
               xx(i,j,k) = us(i,j,k) * vs(i,j,k)
            enddo
         enddo
         endif
      enddo
      call fv_array_sync()
      if ( id_usvs > 0 ) used = send_data ( id_usvs, xx(beglon:endlon,beglat:endlat,:), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array xx is touched by another thread. 
      if( id_tsts>0 .or. id_vsts>0 ) call fv_array_sync()
    endif


    if ( id_vort > 0 .or. id_pv > 0 ) then
         call vort_d(im, jm, km, jfirst, jlast, u, v, omga, ng_s, ng_d)
         if ( id_vort > 0 ) used = send_data ( id_vort, omga(beglon:endlon,beglat:endlat,:), Time )             
         ! The following sync call is required to ensure that send_data completes before 
         ! array omga is touched by another thread. 
         if( id_pv>0 .or. id_ta>0 .or. id_pt>0 ) call fv_array_sync()
    endif


! Potential Vorticity

    if ( id_pv > 0 ) then
         ! Note: this is expensive computation. don't do it too often
         call pv_entropy(im, jm, km, jfirst, jlast, omga, pt, pkz, delp, ng_d, grav)
         used = send_data ( id_pv, omga(beglon:endlon,beglat:endlat,:), Time )
         ! The following sync call is required to ensure that send_data completes before 
         ! array omga is touched by another thread. 
         if( id_ta>0 .or. id_pt>0 ) call fv_array_sync()
    endif

    if ( id_ta > 0 .or. id_tsts > 0 .or. id_vsts > 0 ) then
!$omp parallel do default(shared)                     &
!$omp private(i, j, k)
      do k=ksp,kep
          do j=jfirst,jlast
            do i=1,im
              omga(i,j,k) = pt(i,j,k)
            enddo
          enddo
          if ( id_tsts > 0 ) then
               call zsmean( omga(1,jfirst,k), im, jfirst, jlast,        &
                            ax(jfirst,k),   us(1,jfirst,k), xx(1,jfirst,k) )
          endif
      enddo    ! end parallel k-loop
      call fv_array_sync()
      if ( id_ta > 0)   used = send_data ( id_ta,   omga(beglon:endlon,beglat:endlat,:), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array omga is touched by another thread. 
      if( id_pt>0 ) call fv_array_sync()
      if ( id_tsts > 0) used = send_data ( id_tsts, xx(beglon:endlon,beglat:endlat,:), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array xx is touched by another thread. 
      if( id_vsts>0 ) call fv_array_sync()
      if( prt_minmax ) then
         call pmaxmin('T', omga, qmin, qmax, im*(jlast-jfirst+1), km, 1.)
         if(qmax > 360.) call error_mesg('fv_diagnostics:','too warm', FATAL)
#ifdef MARS_GCM
         if(qmin < 50.) call error_mesg('fv_diagnostics:','too cold', FATAL)
#else
         if(qmin < 100.) call error_mesg('fv_diagnostics:','too cold', FATAL)
#endif 
      endif

    endif

    if ( id_pt > 0 ) then
       p0 = (1.E5) ** kappa  
!$omp parallel do default(shared)                     &
!$omp private(i, j, k)
      do k=ksp,kep
          do j=jfirst,jlast
            do i=1,im
              omga(i,j,k) = pt(i,j,k)/pkz(i,j,k) * p0
            enddo
          enddo
      enddo    ! end parallel k-loop
      call fv_array_sync()
      used = send_data ( id_pt, omga(beglon:endlon,beglat:endlat,:), Time )
      call fv_array_sync()
    endif


    if ( id_vsts > 0 ) then

!$omp parallel do default(shared)                     &
!$omp private(i, j, k)
      do k=ksp,kep
         do j=jfirst,jlast
            do i=1,im
               xx(i,j,k) = vs(i,j,k) * us(i,j,k)
            enddo
         enddo
      enddo
      call fv_array_sync()
      used = send_data ( id_vsts, xx(beglon:endlon,beglat:endlat,:), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array xx is touched by another thread. 
      call fv_array_sync()
    endif

    if ( id_hght > 0 .or. id_h500 > 0  .or. id_slp > 0 .or.    &
         id_p5km > 0  )   then

#ifdef SW_DYN
      gg  = ginv
#else
      gg  = rdgas * ginv
#endif

!      allocate ( wz(im,jfirst:jlast,km+1) )
      call fv_array_sync() !next loops are parallel over j

!$omp parallel do default(shared) private(i,j,k)
      do j=jsp,jep

! Compute height at layer edges
         do i=1,im
            wz(i,j,km+1) = phis(i,j) * ginv
         enddo

         do k=km,1,-1
            do i=1,im
#ifdef SW_DYN
               wz(i,j,k) = wz(i,j,k+1) + gg*(pk(i,j,k+1)-pk(i,j,k))
#else
               wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))   &
                    *(peln(i,k+1,j)-peln(i,k,j))
#endif
            enddo
         enddo
      end do
      call fv_array_sync()

      if(id_hght > 0) used = send_data ( id_hght, wz(beglon:endlon,beglat:endlat,:), Time )

!     call pmaxmin('Ztop', wz(1,jfirst,1), qmin, qmax, im, (jlast-jfirst+1), 0.001)
!     ztop = gmean(im, jm, jfirst, jlast, wz(1,jfirst,1)) * 0.001
!     if( master ) write(6,*) 'Mean model top height (km) =', ztop 

! Compute 5-km pressure
      if(id_p5km > 0) then
         height(1) = 5.E3
         call get_pressure_given_height(im, jfirst, jlast, km, wz, 1, height(1),   &
                                        pt(1,beglat, nlev), a2, 0.01)
         used = send_data ( id_p5km, a2(beglon:endlon,beglat:endlat), Time )
         ! The following sync call is required to ensure that send_data completes before 
         ! array a2 is touched by another thread. 
         if( id_slp>0 .or. id_h500>0 ) call fv_array_sync()
      endif

! Cumpute SLP (pressure at height=0)
      if(id_slp > 0) then
         height(1) = 0.
         call get_pressure_given_height(im, jfirst, jlast, km, wz, 1, height(1),   &
                                        pt(1,beglat, nlev), a2, 0.01)
         if( prt_minmax )   &
             call pmaxmin('SLP', a2, qmin, qmax, im*(jlast-jfirst+1),  1, 1.)
         used = send_data ( id_slp, a2(beglon:endlon,beglat:endlat), Time )
         ! The following sync call is required to ensure that send_data completes before 
         ! array a2 is touched by another thread. 
         if( id_h500>0 ) call fv_array_sync()
      endif

! Compute H500
      if(id_h500 > 0) then
         log_p(1) = log( 50000. )
         call get_height_given_pressure(im, jfirst, jlast, km, wz, 1, log_p, a2)
         used = send_data ( id_h500, a2(beglon:endlon,beglat:endlat), Time )
      endif

!      deallocate ( wz )
    endif

!----------------
! Output tracers:
!----------------

    if( nq /= 0 ) then
        do ic=1, nq
          if (id_tracer(ic) > 0) &
               & used = send_data (id_tracer(ic), q(beglon:endlon,beglat:endlat,1:km,ic), Time )
          if( prt_minmax ) then 
             call pmaxming('Q', q(1,jfirst-ng_d,1,ic), im, jm, km,  &
                  jfirst, jlast, ng_d, ng_d, 1.0)
          endif
        enddo
#ifdef MARS_GCM


#else
!-------
! jgj: per SJ email (jul 17 2008): q_dry = q_moist/(1-sphum)
! mass mixing ratio: q_dry = mass_tracer/mass_dryair = mass_tracer/(mass_air - mass_water) ~ q_moist/(1-sphum)
! co2_mmr = (wco2/wair) * co2_vmr
! Note: There is a check in fv_pack.F90 that ensures that tracer number one is sphum

        do ic = nt_phys+1, ncnst
          call get_tracer_names ( MODEL_ATMOS, ic, tname, tlongname, tunits )
          if (id_tracer_dmmr(ic) > 0 .and. trim(tname).eq.'co2') then
          call fv_array_sync() !next loops are parallel over j
!$omp parallel do default(shared) private(i,j,k)
            do k=1,km     
              do j=jsp,jep
                do i=1,im
                  co2_dmmr(i,j,k) = q(i,j,k,ic)/(1.0-q(i,j,k,1))
                enddo
              enddo
            enddo
            call fv_array_sync()
            used = send_data (id_tracer_dmmr(ic), co2_dmmr(beglon:endlon,beglat:endlat,1:km), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array co2_dmmr is touched by another thread. 
            call fv_array_sync()
            if( prt_minmax ) then 
               call pmaxming('co2_dmmr', co2_dmmr(1,jfirst-ng_d,1), im, jm, km,  &
                    jfirst, jlast, ng_d, ng_d, 1.0)
            endif
          endif

!! use dry vmr to compare to obs
          if (id_tracer_dvmr(ic) > 0 .and. trim(tname).eq.'co2') then
          call fv_array_sync() !next loops are parallel over j
!$omp parallel do default(shared) private(i,j,k)
            do k=1,km     
              do j=jsp,jep
                do i=1,im
                  co2_dvmr(i,j,k) = (q(i,j,k,ic)/(1.0-q(i,j,k,1))) * WTMAIR/WTMCO2   ! to compare to obs
                enddo
              enddo
            enddo
            call fv_array_sync()
            used = send_data (id_tracer_dvmr(ic), co2_dvmr(beglon:endlon,beglat:endlat,1:km), Time )
      ! The following sync call is required to ensure that send_data completes before 
      ! array co2_dmmr is touched by another thread. 
            call fv_array_sync()
            if( prt_minmax ) then 
               call pmaxming('co2_dvmr', co2_dvmr(1,jfirst-ng_d,1), im, jm, km,  &
                    jfirst, jlast, ng_d, ng_d, 1.0)
            endif
          endif
        enddo
#endif 
!-------

! Check tracer mass. Unit:
!       if( prt_minmax )        &
!           call tracer_mass(im, jm, km, jfirst, jlast, ng_d, nq, q, delp)

        if ( age_tracer ) then
! Assuming age tracer is the last tracer
             call age_of_air(im, jm, km, jfirst, jlast, ng_d, age_time,   &
                             dt_atmos, peln, delp, q(1,jfirst-ng_d,1,nq), &
                             prt_minmax)
        endif

    endif

!    if ( id_usus>0 .or. id_usvs>0 .or. id_vsts>0 .or. id_tsts>0 ) then
!         deallocate ( xx )
!         deallocate ( us )
!         deallocate ( vs )
!    endif
 
    ! make sure all send_data completed before one thread exits the routine
    call fv_array_sync()

 end subroutine fv_diag

 subroutine get_pressure_given_height(im, jfirst, jlast, km, wz, kd, height, &
      ts, a2, fac)
#include "fv_arrays.h"

 integer,  intent(in):: im, jfirst, jlast, km
 integer,  intent(in):: kd           ! vertical dimension of the ouput height
 real, intent(in):: wz(im,jfirst:jlast,km+1)
 real, intent(in):: ts(im,jfirst:jlast)
 real, intent(in):: height(kd)   ! must be monotonically decreasing with increasing k
 real, optional, intent(in):: fac
 real, intent(out):: a2(im,jfirst:jlast,kd)      ! pressure (pa)

! local:
 integer n, i,j,k
 real ptmp, tm
 integer ii, mcount
#include "fv_point.inc"

 call fv_array_check( LOC(wz) )
 call fv_array_check( LOC(ts) )
 call fv_array_check( LOC(a2) )


 do n=1,kd

!$omp parallel do private(ii,i,j,k, ptmp, tm, mcount)
    do j=jsp,jep

       do 1000 i=1,im

         if ( height(n) >= wz(i,j,km+1) ) then
!---------------------
! Search from top down
!---------------------
          do k=1,km
             if( height(n) < wz(i,j,k) .and. height(n) >= wz(i,j,k+1) ) then
! Found it!
                 ptmp = peln(i,k,j) + (peln(i,k+1,j)-peln(i,k,j)) *   &
                       (wz(i,j,k)-height(n)) / (wz(i,j,k)-wz(i,j,k+1))
                 a2(i,j,n) = exp(ptmp)
                 go to 500
             endif
          enddo

         else

! Method: use a mean tm based on zonal mean with height < local height
             tm = 0.
             mcount = 0
          do ii=1,im
             if (wz(ii,j,km+1) < wz(i,j,km+1) .and. wz(ii,j,km+1) > (height(n)+0.1)) then
! adding 0.1 to exclude most ocean points
! This algorithm works very well except the S pole.
             mcount = mcount + 1
             tm = tm + ts(ii,j)
             endif
          enddo
             if ( mcount > im/10 + 1 ) then
                tm = rdgas*ginv * tm / real(mcount)
             else
! Assuming 6.5 deg/km lapse rate
                tm = rdgas*ginv*(ts(i,j) + 3.25E-3*(wz(i,j,km)-height(n)))
             endif
          a2(i,j,n) = exp( peln(i,km+1,j) + (wz(i,j,km+1) - height(n))/tm )
         endif
500      if ( present(fac) ) a2(i,j,n) = fac * a2(i,j,n)
1000   continue
    enddo
 enddo
 call fv_array_sync()

 end subroutine get_pressure_given_height


 subroutine get_height_given_pressure(im, jfirst, jlast, km, wz, kd, log_p, a2)

#include "fv_arrays.h"
 integer,  intent(in):: im, jfirst, jlast, km
 integer,  intent(in):: kd           ! vertical dimension of the ouput height
 real, intent(in):: log_p(kd)    ! must be monotonically decreasing with increasing k
                                     ! log (p)
 real, intent(in):: wz(im,jfirst:jlast,km+1)

 real, intent(out):: a2(im,jfirst:jlast,kd)      ! height (m)

! local:
 integer n, i,j,k
#include "fv_point.inc"

 call fv_array_check( LOC(wz) )
 call fv_array_check( LOC(a2) )

 do n=1,kd

!$omp parallel do private(i,j,k)
    do j=jsp,jep
       do 1000 i=1,im
          do k=1,km
             if( log_p(n) <= peln(i,k+1,j) .and. log_p(n) >= peln(i,k,j) ) then
! Found it!
                 a2(i,j,n) = wz(i,j,k)  +  (wz(i,j,k+1) - wz(i,j,k)) *   &
                       (log_p(n)-peln(i,k,j)) / (peln(i,k+1,j)-peln(i,k,j) )
                 go to 1000
             endif
          enddo
                 a2(i,j,n) = missing_value
1000   continue
    enddo
 enddo
 call fv_array_sync()

 end subroutine get_height_given_pressure


 subroutine zsmean(a2, im, js, je, ux, us, usus)
 integer im, js, je
 real, intent(in) :: a2(im,js:je)
 real, intent(out):: ux(js:je)          ! zonal mean
 real, intent(out):: us(im,js:je)       ! departure from ux
 real, optional, intent(out):: usus(im,js:je)       ! us ** 2

! local
 integer i, j
 real  rim

      rim = 1. / im

      do j=js,je
            ux(j) = a2(1,j)
         do i=2,im
            ux(j) = ux(j) + a2(i,j)
         enddo
            ux(j) = ux(j) * rim

        do i=1,im
            us(i,j) = a2(i,j) - ux(j)
        enddo
        if ( present(usus) ) then
           do i=1,im
              usus(i,j) = us(i,j) ** 2
           enddo
        endif
      enddo
 end subroutine zsmean

 subroutine age_of_air(im, jm, km, jfirst, jlast, ng, atime, dt,   &
      peln, delp, q, print_max)
   use fv_arrays_mod, only: fv_array_check, isp, iep, jsp, jep, ksp, kep


 integer im
 integer jm
 integer km
 integer jfirst
 integer jlast  
 integer ng
 logical, intent(in):: print_max

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

  real, intent(in):: delp(im,jfirst:jlast,km)
  real, intent(in):: peln(im,km+1,jfirst:jlast)
  real, intent(inout):: atime        ! accumulated time since init
  real, intent(in):: dt
  real, intent(inout):: q(im,jfirst-ng:jlast+ng,km)

! Local
  integer i, j, k
  real p_source      ! source level (pa)
  real ascale, ra, ryear
  real tiny, pfull
  real rim, qsum, qmin, qmax
  real age(im,jfirst:jlast,km)
  parameter ( tiny = 1.e-6 )
  parameter ( p_source = 75000. )
  parameter ( ascale = 5.e-6 / 60. )
  parameter ( ra = 1./ascale )
  parameter (ryear = 1./(365.*24.*3600) )

  rim = 1./float(im)
  atime = atime + dt

  call fv_array_check( LOC(delp) )
  call fv_array_check( LOC(peln) )
  call fv_array_check( LOC(q) )

!$omp parallel do private(i, j, k, qsum, pfull)
  do k=ksp,kep
     do j=jfirst, jlast
        do i=1,im
           pfull = delp(i,j,k) / (peln(i,k+1,j)-peln(i,k,j))
           if( atime < tiny ) then
               q(i,j,k) = 0.
           elseif( pfull >= p_source ) then
               q(i,j,k) = ascale * atime
           endif
        enddo
        if(id_aoa > 0 ) then
            do i=1,im
               age(i,j,k) =  max(0., (atime - ra*q(i,j,k))) * ryear
            enddo
        endif
     enddo
  enddo
  call fv_array_sync()

  if(id_aoa > 0 ) then
     used = send_data (id_aoa, age(beglon:endlon,beglat:endlat,:), fv_time)
     if(print_max )         &
        call pmaxmin('Age-of-Air (yr)', age, qmin, qmax, im*(jlast-jfirst+1), km, 1.)
  endif

  ! make sure all send_data completed before one thread exits the routine
  call fv_array_sync()

 end subroutine age_of_air

end module fv_diagnostics

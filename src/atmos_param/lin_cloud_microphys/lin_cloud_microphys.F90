!
! Cloud micro-physics package for GFDL global cloud resolving model
! The algorithms are originally based on Lin et al 1983. Many key 
! elements have been changed/improved based on several other publications
! Developer: Shian-Jiann Lin
!
module lin_cld_microphys_mod
 use mpp_mod,           only: stdlog, mpp_pe, mpp_root_pe, mpp_clock_id, &
                              mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE, &
                              input_nml_file
 use diag_manager_mod,  only: register_diag_field, send_data
 use time_manager_mod,  only: time_type, get_date, get_time
 use constants_mod,     only: grav, rdgas, rvgas, cp_air, hlv, hlf, kappa
 use fms_mod,           only: write_version_number, open_namelist_file, &
                              check_nml_error, file_exist, close_file,  &
                              error_mesg, FATAL 

 implicit none
 private

 public  lin_cld_microphys_driver, lin_cld_microphys_init, lin_cld_microphys_end, sg_conv
 public  qsmith_init, qsmith, es2_table1d, es3_table1d, esw_table1d
 real             :: missing_value = -1.e10
 logical          :: module_is_initialized = .false.
 character(len=17) :: mod_name = 'lin_cld_microphys'

!==== fms constants ====================
!real :: rdgas = 287.04
!real :: rvgas = 461.50
 real, parameter :: cp    = cp_air          ! heat capacity at constant pressure (j/kg/k)
 real, parameter :: eps   = rdgas/rvgas     ! = 0.621971831
 real, parameter :: zvir  = rvgas/rdgas-1.  ! = 0.607789855
 real, parameter :: latv  = hlv             ! = 2.500e6
 real, parameter :: lati  = hlf             ! = 3.34e5
 real, parameter :: lats  = hlv+hlf         ! = 2.834E6
!==== fms constants ====================

 real, parameter :: qrmin  = 1.e-9
 real, parameter :: qvmin  = 1.e-20      ! min value for water vapor (treated as zero)
 real, parameter :: qcmin  = 1.e-12      ! min value for cloud condensates
 real, parameter :: sfcrho = 1.20        ! surface air density
 real, parameter :: vmin   = 1.e-2       ! minimum fall speed for rain/graupel
 real, parameter :: tice   = 273.16  ! melting  starts above tice
 real, parameter :: tice0  = 273.15  ! freezing starts below tice0
 real, parameter :: rhor   = 1.0e3  ! LFO83
 real, parameter :: f_l2s  = 50.
 real, parameter:: dz_min = 1.e-2

 real :: cracs, csacr, cgacr, cgacs, acco(3,4), csacw,          &
         craci, csaci, cgacw, cgaci, cracw, cssub(5), cgsub(5), &
         crevp(5), cgfr(2), csmlt(5), cgmlt(5)
 real :: rmi50, es0, ces0, c1brg, c2brg


 real :: dts, rdts, pie  ! these variables have been left unchanged
 real :: lcp, icp, tcp, rgrav
 real :: fac_rc
 real :: mp_count = 0.

 logical :: do_setup=.true.
 logical :: master 
 logical :: g_sum_initialized
 real, allocatable, dimension(:,:) :: l_area

 real, allocatable:: vt_r(:,:,:), vt_s(:,:,:), vt_g(:,:,:), vt_i(:,:,:)
 real, allocatable:: prec0(:,:), rain0(:,:), snow0(:,:), ice0(:,:), graupel0(:,:)
 real, allocatable:: prec1(:,:), prec_mp(:,:), cond(:,:), w_var(:,:)
 real, allocatable:: table(:), table2(:), table3(:), tablew(:), des(:), des2(:), des3(:), desw(:)

 integer:: isc, iec, jsc, jec
 integer:: id_vtr, id_vts,  id_vtg, id_vti, id_rain, id_snow, id_graupel, &
           id_ice, id_prec, id_cond, id_var

 real, parameter :: dt_fr = 6.       ! homogeneous freezing of all cloud water at t_wfr - dt_fr
 real, parameter :: t_wfr = tice-42. ! supercooled water can exist down to -48 C, which is the "absolute"
                                     ! minimum temperature water can exist (Moore & Molinero Nov. 2011, Nature)
                                     ! dt_fr can be considered as the error bar
 real, parameter :: t_00 =  t_wfr - dt_fr   ! This is the absolute freezing point for super-cooled cloud water
 integer, parameter:: ng    = 0     ! NO ghost zones required as "area" is passed from the phys driver
 integer :: lin_cld_mp_clock   ! clock for timing of driver routine

 real :: t_snow_melt = 10.      ! snow melt tempearture scale factor
 real :: q00     = 1.0e-3
!----------------------
! namelist  parameters:
!----------------------
 real :: qc_crt  = 1.0e-7  ! minimum condensate mixing ratio to allow partial cloudiness
 real :: t_min   = 165.  ! Min temperature for ice-phase micro phys
 real :: mp_time = 120.  ! maximum micro-physics time step (sec)

! The following 3 time scales are for terminal falls
 real :: tau_s  = 120.   ! snow melt
 real :: tau_g  = 150.   ! graupel melt

! Ice:
 real :: tau_frz = 600.   ! cloud water freezing time-scale (mixed phase)
 real :: tau_mlt = 10.    ! ice melting time-scale
 real :: tau_i2v = 30.    ! ice   ---> vapor
 real :: tau_v2i = 150.   ! vapor ---> ice
! cloud water
 real :: tau_l2v = 30.   ! cloud water --> vapor (evaporation)  time scale
 real :: tau_v2l = 150.  ! vapor --> cloud water (condensation) time scale
! Snow
 real :: tau_s2v = 600.   ! snow to vapor (after liquid/ice sat adj)
 real :: tau_v2s = 600.   ! vapor to snow
! Graupel
 real :: tau_g2v = 1200.  ! Grapuel sublimation time scale
 real :: tau_v2g = 1200.  ! Grapuel deposition -- make it a slow process 

 real :: dw_land  = 0.18  ! base value for subgrid deviation/variability over land 
 real :: dw_ocean = 0.14  ! base value for ocean
 real :: rh_inc = 0.05    ! parameter to control instant evap of condensates
 real :: ccn_o =  70.    
 real :: ccn_l = 200.    
 real :: rthresh = 8.0e-6     ! critical cloud drop radius (micro m)

!-------------------------------------------------------------
 real :: qi0_crt = 1.0e-4    ! ice  --> snow autocon mixing ratio threshold
 real :: qr0_crt = 2.0e-4    ! rain --> snow or graupel/hail threshold
                             ! LFO used *mixing ratio* = 1.E-4 (hail in LFO)
 real :: c_psaut = 1.0e-3   ! autoconversion rate: cloud_ice -> snow
 real :: c_psaci = 0.1      ! accretion: cloud ice --> snow (was 0.1 in Zetac)
 real :: c_piacr = 1.       ! accretion: rain --> ice:
 real :: c_cracw = 1.0      ! rain accretion efficiency

! Decreasing  clin to reduce csacw (so as to reduce cloud water ---> snow)
 real:: alin = 842.0
 real:: clin = 4.8      ! 4.8 --> 2.4?

!-----------------
! Graupel control:
!-----------------
 real :: qs0_crt = 1.0e-3   ! snow --> graupel density threshold (6.0e-4 in Purdue Lin scheme)
 real :: c_pgacs = 0.01     ! snow --> graupel "accretion" eff. (was 0.1 in Zetac)

! fall velocity tuning constants:
 real :: den_ref = sfcrho   ! Reference (surface) density for fall speed
                            ! Larger value produce larger fall speed
 real :: vr_fac = 1.
 real :: vs_fac = 1.
 real :: vg_fac = 1.
 real :: vi_fac = 1.

 logical :: z_slope  = .false.          !  use linear mono slope for autocconversions
 logical :: use_deng_mace = .true.       ! Helmfield-Donner ice speed
 logical :: do_subgrid_z = .false.       ! 2X resolution sub-grid saturation/cloud scheme
 logical :: use_ccn      = .true.
 logical :: use_ppm      = .true.
 logical :: ppm_rain_fall  = .true.
 logical :: mono_prof = .true.          ! perform terminal fall with mono ppm scheme
 logical :: mp_debug = .false.
 logical :: mp_print = .true.

 real :: p_crt   = 200.E2   ! 
 integer :: k_moist = 20 
 real:: rh_adj
 real:: fac_l2v, fac_v2l, fac_mlt, fac_sno, fac_i2v, fac_v2i
 real:: fac_s2v, fac_v2s, fac_g2v, fac_v2g

 namelist /lin_cld_microphys_nml/mp_time, t_min, tau_s, tau_g, dw_land, dw_ocean,  &
                      tau_frz, vr_fac, vs_fac, vg_fac, vi_fac,       &
                      qs0_crt, qi0_crt, qr0_crt,    &
                      rh_inc, den_ref, use_deng_mace, use_ccn, do_subgrid_z,  &
                      rthresh, ccn_l, ccn_o, qc_crt,  &
                      c_piacr, tau_mlt, tau_i2v, tau_v2i,  tau_l2v, tau_v2l,     &
                      c_psaut, c_psaci, c_pgacs, z_slope,  &
                      c_cracw, alin, clin,    &
                      use_ppm, ppm_rain_fall, mono_prof, mp_debug, mp_print

!---- version number -----
 character(len=128) :: version = '$Id: lin_cloud_microphys.F90,v 19.0.2.1 2012/06/10 04:47:37 Rusty.Benson Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'

 contains
 

  subroutine lin_cld_microphys_driver(qv,    ql,    qr,    qi,    qs,    qg,    qa,  &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      & 
                               pt_dt, pt, p3, dz,  delp, area, dt_in,                &
                               land,  rain, snow, ice, graupel,                      &
                               hydrostatic, phys_hydrostatic,                        &
                               iis,iie, jjs,jje, kks,kke, ktop, kbot, time)

  type(time_type), intent(in):: time
  logical,         intent(in):: hydrostatic, phys_hydrostatic
  integer,         intent(in):: iis,iie, jjs,jje  ! physics window
  integer,         intent(in):: kks,kke           ! vertical dimension
  integer,         intent(in):: ktop, kbot        ! vertical compute domain
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(:,:)  :: area
  real, intent(in   ), dimension(:,:)  :: land  !land fraction
  real, intent(out  ), dimension(:,:)  :: rain, snow, ice, graupel
  real, intent(in   ), dimension(:,:,:):: p3, delp, dz    ! p3 not used
  real, intent(in   ), dimension(:,:,:):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(:,:,:):: pt_dt,  qa_dt
  real, intent(inout), dimension(:,:,:):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                          qs_dt, qg_dt


! local:
  logical used
  real    :: mpdt, rdt, convt, tot_prec
  integer :: i,j,k
  integer :: is,ie, js,je  ! physics window
  integer :: ks,ke         ! vertical dimension
  integer :: seconds, days, ntimes

  is = 1
  js = 1
  ks = 1
  ie = iie-iis+1
  je = jje-jjs+1
  ke = kke-kks+1

  call mpp_clock_begin (lin_cld_mp_clock)

! tendency zero out for am moist processes should be done outside the driver

     mpdt = min(dt_in, mp_time)
      rdt = 1. / dt_in
   ntimes = nint( dt_in/mpdt )
! small time step:
      dts = dt_in / real(ntimes)
     rdts = 1./dts

  fac_l2v = 1. - exp( -dts/tau_l2v )        ! exact-in-time integration
  fac_v2l = 1. - exp( -dts/tau_v2l )        ! exact-in-time integration

  fac_mlt = 1. - exp( -dts/tau_mlt )        ! 
  fac_i2v = 1. - exp( -dts/tau_i2v )        ! 
  fac_v2i = 1. - exp( -dts/tau_v2i )        ! 

  fac_sno = 1. - exp( -dts/tau_s   )        ! 

  fac_s2v = 1. - exp( -dts/tau_s2v )
  fac_v2s = 1. - exp( -dts/tau_v2s )

 fac_g2v = 1. - exp( -dts/tau_g2v )
 fac_v2g = 1. - exp( -dts/tau_v2g )

  call get_time (time, seconds, days)


  do j=js, je
     do i=is, ie
        graupel(i,j) = 0.
           rain(i,j) = 0.
           snow(i,j) = 0.
            ice(i,j) = 0.
           cond(i,j) = 0.
     enddo
  enddo
  do j=js,je
     call mpdrv( delp, pt, qv, ql, qr, qi, qs, qg, qa, dz,  &
                 is, ie, js, je, ks, ke, ktop, kbot, j, dt_in,  & 
                 ntimes, rain(:,j), snow(:,j), graupel(:,j), &
                 ice(:,j), cond(:,j), area(:,j), land(:,j),  &
                 pt_dt, qv_dt, ql_dt, qr_dt, qi_dt,    &
                 qs_dt, qg_dt, qa_dt )
  enddo

! no clouds allowed above ktop
   if ( ks < ktop ) then
      do k=ks, ktop
         do j=js,je
            do i=is,ie
!              qa(i,j,k) = 0.
               qa_dt(i,j,k) = -qa(i,j,k) * rdt
            enddo
         enddo
      enddo
   endif

#ifdef SIM_PHYS
   if ( id_vtr> 0 ) used=send_data(id_vtr, vt_r, time)
   if ( id_vts> 0 ) used=send_data(id_vts, vt_s, time)
   if ( id_vtg> 0 ) used=send_data(id_vtg, vt_g, time)
   if ( id_vti> 0 ) used=send_data(id_vti, vt_i, time)
   if ( id_var> 0 ) used=send_data(id_var, w_var,time)
#else
   if ( id_vtr> 0 ) used=send_data(id_vtr, vt_r, time, iis, jjs)
   if ( id_vts> 0 ) used=send_data(id_vts, vt_s, time, iis, jjs)
   if ( id_vtg> 0 ) used=send_data(id_vtg, vt_g, time, iis, jjs)
   if ( id_vti> 0 ) used=send_data(id_vti, vt_i, time, iis, jjs)
   if ( id_var> 0 ) used=send_data(id_var, w_var,time, iis, jjs)
#endif

! Convert to mm/day
   convt = 86400.*rdt*rgrav
   do j=js,je
      do i=is,ie
            rain(i,j) =    rain(i,j) * convt
            snow(i,j) =    snow(i,j) * convt
             ice(i,j) =     ice(i,j) * convt
         graupel(i,j) = graupel(i,j) * convt
         prec_mp(i,j) =    rain(i,j) + snow(i,j) + ice(i,j) + graupel(i,j)
      enddo
   enddo

   if ( id_cond>0 ) then
        do j=js,je
           do i=is,ie
              cond(i,j) = cond(i,j)*rgrav
           enddo
        enddo
#ifdef SIM_PHYS
        used=send_data(id_cond, cond, time)
#else
        used=send_data(id_cond, cond, time, iis, jjs)
#endif
   endif

   if ( id_snow>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_snow,    snow,    time)
#else
        used=send_data(id_snow,    snow,    time, iis, jjs)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(snow, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean snow=', tot_prec
        endif
        snow0(:,:) = snow0(:,:) + snow(:,:)
   endif

   if ( id_graupel>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_graupel, graupel, time)
#else
        used=send_data(id_graupel, graupel, time, iis, jjs)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(graupel, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean graupel=', tot_prec
        endif
        graupel0(:,:) = graupel0(:,:) + graupel(:,:)
   endif

   if ( id_ice>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_ice, ice, time)
#else
        used=send_data(id_ice, ice, time, iis, jjs)
#endif
        if ( seconds==0 ) then
             tot_prec = g_sum(ice, is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'mean ice_mp=', tot_prec
        endif
        ice0(:,:) = ice0(:,:) + ice(:,:)
   endif

   if ( id_rain>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_rain,    rain,    time)
#else
        used=send_data(id_rain,    rain,    time, iis, jjs)
#endif
        if ( seconds==0 ) then
!            tot_prec = g_sum(rain, is, ie, js, je, ng, area, 1) 
!            if(master) write(*,*) 'mean rain=', tot_prec
        endif
        rain0(:,:) = rain0(:,:) + rain(:,:)
   endif
   

   if ( id_prec>0 ) then
#ifdef SIM_PHYS
        used=send_data(id_prec, prec_mp, time)
#else
        used=send_data(id_prec, prec_mp, time, iis, jjs)
#endif
   endif

!----------------------------------------------------------------------------

        prec0(:,:) = prec0(:,:) + prec_mp(:,:)
        prec1(:,:) = prec1(:,:) + prec_mp(:,:)
        mp_count = mp_count + 1.

        if ( seconds==0 .and. mp_print ) then
             tot_prec = g_sum(prec1*dt_in/86400., is, ie, js, je, ng, area, 1) 
             if(master) write(*,*) 'Daily prec_mp=', tot_prec
!            call prt_maxmin('prec_mp', prec1*dt_in/86400., is, ie, js, je, 0, 1, 1., master)
             prec1(:,:) = 0.
        endif
!----------------------------------------------------------------------------


!rab  if ( mp_debug ) then
!rab       call prt_maxmin('T_a_mp',    pt, is, ie, js, je, 0, kbot, 1., master)
!rab       call prt_maxmin('qg_dt_a_mp',  qg_dt, is, ie, js, je, 0, kbot, 1., master)
!rab       call prt_maxmin('prec', prec_mp, is, ie, js, je, 0,    1, 1., master)
!rab  endif

   call mpp_clock_end (lin_cld_mp_clock)

 end subroutine lin_cld_microphys_driver



 subroutine mpdrv( delp, pt, qv, ql, qr, qi, qs, qg, qa, dz,     &
                   is, ie, js, je, ks, ke, ktop, kbot, j, dt_in, ntimes,  & 
                   rain, snow, graupel, ice, &
                   cond, area1, land, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt,    &
                   qs_dt, qg_dt, qa_dt )

!-------------------------------------------------------------------
!  lin et al., 1983, jam, 1065-1092, and
!  rutledge and hobbs, 1984, jas, 2949-2972
!-------------------------------------------------------------------
! terminal fall is handled lagrangianly by conservative fv algorithm
!
! pt: temperature (k)
! 6 water species:
! 1) qv: water vapor (kg/kg)
! 2) ql: cloud water (kg/kg)
! 3) qr: rain        (kg/kg)
! 4) qi: cloud ice   (kg/kg)
! 5) qs: snow        (kg/kg)
! 6) qg: graupel     (kg/kg)

  integer,         intent(in):: j, is,ie, js,je, ks,ke
  integer,         intent(in):: ntimes, ktop, kbot
  real,            intent(in):: dt_in

  real, intent(in), dimension(is:ie,js:je,ks:ke) :: delp
  real, intent(in), dimension(is:ie):: area1, land
  real, intent(in   ), dimension(is:ie,js:je,ks:ke):: pt, qv, ql, qr, qi, qs, qg, qa, dz
  real, intent(inout), dimension(is:ie,js:je,ks:ke):: pt_dt,  qa_dt,  &
                                            qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  real, intent(out), dimension(is:ie):: rain, snow, ice, graupel, cond
!----------
! local var
!----------
  real, dimension(ktop:kbot):: qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
                               vtiz, vtsz, vtgz, vtrz, &
                               dp1, qv0, ql0, qr0, qi0, qs0, qg0, qa0, t0, den, &
                               den0, tz, p1, dz0, dz1, denfac

  real :: r1, s1, i1, g1, rdt, omq
  real :: cpaut, ccn, c_praut
  real :: dt_rain
  real :: h_var, s_leng, t_land, t_ocean
  integer :: i,k,n
! real:: x, pexp
! pexp(x) = 1.+x*(1.+x*(0.5+x/6.*(1.+x*(0.25+0.05*x))))

   dt_rain = dts * 0.5

   rdt = 1. / dt_in

   cpaut = 0.55*0.104*grav/1.717e-5

   do 2000 i=is, ie

   do k=ktop, kbot
       t0(k) = pt(i,j,k)
       tz(k) = t0(k) 
!-----------------------------------
      qvz(k) = max(qvmin, qv(i,j,k))
      qlz(k) = max(qvmin, ql(i,j,k))
      qrz(k) = max(qvmin, qr(i,j,k))
      qiz(k) = max(qvmin, qi(i,j,k))
      qsz(k) = max(qvmin, qs(i,j,k))
      qgz(k) = max(qvmin, qg(i,j,k))
!-----------------------------------
      qa0(k) = qa(i,j,k)
      qaz(k) = 0.  
      dz0(k) = dz(i,j,k)
!--------------------------
         omq = 1. - (qvz(k)+qlz(k)+qrz(k)+qiz(k)+qsz(k)+qgz(k))
      dp1(k) = delp(i,j,k) * omq         ! dry air mass * grav
     den0(k) = -dp1(k)/(grav*dz0(k))     ! density of dry air
       p1(k) = den0(k)*rdgas*t0(k)       ! dry pressure
!------------------------------
! convert to dry mixing ratios:
!------------------------------
         omq = 1. / omq
      qvz(k) = qvz(k)*omq
      qv0(k) = qvz(k)
      qlz(k) = qlz(k)*omq
      ql0(k) = qlz(k)
      qrz(k) = qrz(k)*omq
      qr0(k) = qrz(k)
      qiz(k) = qiz(k)*omq
      qi0(k) = qiz(k)
      qsz(k) = qsz(k)*omq
      qs0(k) = qsz(k)
      qgz(k) = qgz(k)*omq
      qg0(k) = qgz(k)
   enddo

! Compute dry pressure for non-hydrostatic case
!-----------------------------------------------
!  if ( .not. phys_hydrostatic ) then
!      do k=ktop, kbot
!         p1(k) = den0(k)*rdgas*t0(k)
!      enddo
!  endif
!-----------------------------------------------

! Based on Klein Eq. 15
   ccn = (ccn_l*land(i) + ccn_o*(1.-land(i))) * 1.e6
   if ( use_ccn ) then
!  CCN is formulted as CCN = CCN_surface * (den/den_surface)
       ccn = ccn * rdgas*tz(kbot)/p1(kbot)
   endif
   c_praut = cpaut * (ccn*rhor)**(-1./3.)

!--------------------------------------------------------
! Total water subgrid deviation in horizontal direction
!--------------------------------------------------------
!       default area dependent form: use dx ~ 100 km as the base 
   s_leng  = sqrt(sqrt(area1(i)/1.E10))
   t_land  = dw_land  * s_leng
   t_ocean = dw_ocean * s_leng
   h_var = t_land*land(i) + t_ocean*(1.-land(i))
   h_var = min(0.22, max(0.001, h_var))            ! cap:
   if ( id_var>0 ) w_var(i,j) = h_var

   rh_adj = 1. - h_var - rh_inc

!-------------------------
! * fix all negatives
!-------------------------

 call neg_adj(ktop, kbot, p1, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

 do 1000 n=1,ntimes

   do k=ktop, kbot
         dz1(k) = dz0(k)*tz(k)/t0(k) 
         den(k) = den0(k)*dz0(k)/dz1(k)
      denfac(k) = sqrt(sfcrho/den(k))
   enddo

!-------------------------------------------
! Time-split warm rain processes: first pass
!-------------------------------------------
!                                       call timing_on (" warm_rain")
   call warm_rain(dt_rain, ktop, kbot, p1, dp1, dz1, tz, qvz, qlz, qrz, p1, den, denfac, ccn, c_praut, h_var, vtrz, r1)
!                                       call timing_off(" warm_rain")
   rain(i) = rain(i) + r1

!------------------------------------------------
! * sedimentation of cloud ice, snow, and graupel
!------------------------------------------------
!                                       call timing_on (" terminal_fall")
   call fall_speed(ktop, kbot, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)

   call terminal_fall ( dts, ktop, kbot, tz, qvz, qlz, qrz, qgz, qsz, qiz, p1, &
                        dz1, dp1, den, vtgz, vtsz, vtiz,    &
                        r1, g1, s1, i1 )
!                                       call timing_off(" terminal_fall")

      rain(i) = rain(i)    + r1  ! from melted snow & ice that reached the ground
      snow(i) = snow(i)    + s1
   graupel(i) = graupel(i) + g1
       ice(i) = ice(i)     + i1

!-------------------------------------------
! Time-split warm rain processes: 2nd pass
!-------------------------------------------
!                                       call timing_on (" warm_rain")
   call warm_rain(dt_rain, ktop, kbot, p1, dp1, dz1, tz, qvz, qlz, qrz, p1, den, denfac, ccn, c_praut, h_var, vtrz, r1)
!                                       call timing_off(" warm_rain")
   rain(i) = rain(i) + r1

!-------------------------
! * ice-phase microphysics
!-------------------------

!                                       call timing_on (" icloud")
   call icloud( ktop, kbot, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz,  &
                dp1, den, denfac, vtsz, vtgz, vtrz, qaz, h_var )
!                                       call timing_off(" icloud")
1000  continue  ! sub-cycle

   do k = ktop, kbot
               omq = dp1(k) / delp(i,j,k)
      pt_dt(i,j,k) = pt_dt(i,j,k) + rdt*(tz(k)- t0(k)) *omq
      qv_dt(i,j,k) = qv_dt(i,j,k) + rdt*(qvz(k)-qv0(k))*omq
      ql_dt(i,j,k) = ql_dt(i,j,k) + rdt*(qlz(k)-ql0(k))*omq
      qr_dt(i,j,k) = qr_dt(i,j,k) + rdt*(qrz(k)-qr0(k))*omq
      qi_dt(i,j,k) = qi_dt(i,j,k) + rdt*(qiz(k)-qi0(k))*omq
      qs_dt(i,j,k) = qs_dt(i,j,k) + rdt*(qsz(k)-qs0(k))*omq
      qg_dt(i,j,k) = qg_dt(i,j,k) + rdt*(qgz(k)-qg0(k))*omq
      qa_dt(i,j,k) = qa_dt(i,j,k) + rdt*( qaz(k)/real(ntimes)-qa0(k))
!        dz(i,j,k) = dz1(k)
   enddo

!-----------------
! fms diagnostics:
!-----------------

   if ( id_cond>0 ) then
     do k=ktop,kbot                   ! total condensate
        cond(i) = cond(i) + dp1(k)*(qlz(k)+qrz(k)+qsz(k)+qiz(k)+qgz(k))
     enddo
   endif

   if ( id_vtr> 0 ) then
        do k=ktop, kbot
           vt_r(i,j,k) = vtrz(k)
        enddo
   endif
   if ( id_vts> 0 ) then
        do k=ktop, kbot
           vt_s(i,j,k) = vtsz(k)
        enddo
   endif
   if ( id_vtg> 0 ) then
        do k=ktop, kbot
           vt_g(i,j,k) = vtgz(k)
        enddo
   endif
   if ( id_vts> 0 ) then
        do k=ktop, kbot
           vt_i(i,j,k) = vtiz(k)
        enddo
   endif

2000  continue

 end subroutine mpdrv



 subroutine warm_rain( dt, ktop, kbot, p1, dp, dz, tz, qv, ql, qr, pm,  &
                       den, denfac, ccn, c_praut, h_var, vtr, r1)

 integer, intent(in):: ktop, kbot
 real,    intent(in):: dt                    ! time step (s)
 real,    intent(in), dimension(ktop:kbot):: p1, dp, dz, pm, den, denfac
 real,    intent(in):: ccn, c_praut, h_var
 real, intent(inout), dimension(ktop:kbot):: tz, qv, ql, qr, vtr
 real, intent(out):: r1
 
! local:
 real, parameter:: so3 = 7./3.
 real, dimension(ktop:kbot):: dl
 real, dimension(ktop:kbot+1):: ze, zt
 real:: sink, dq, qc0, qc, q_plus, q_minus
 real:: rho0, qden
 real:: zs = 0.
 real:: dt5
 integer k
!-----------------------------------------------------------------------
! fall velocity constants:
!-----------------------------------------------------------------------
 real, parameter :: vconr = 2503.23638966667
 real, parameter :: normr = 25132741228.7183
 real, parameter :: thr=1.e-10
 logical no_fall

!---------------------
! warm-rain processes:
!---------------------

  dt5 = 0.5*dt

!------------------------
! Terminal speed of rain:
!------------------------

  call check_column(ktop, kbot, qr, no_fall)
  if ( no_fall ) then
       vtr(:) = vmin
       r1 = 0.
       go to 999   ! jump to auto-conversion
  endif

  if ( den_ref < 0. ) then
       rho0 = -den_ref*den(kbot) 
  else
       rho0 = den_ref   ! default=1.2
  endif

  do k=ktop, kbot
     qden = qr(k)*den(k)
     if ( qr(k) < thr ) then
         vtr(k) = vmin
     else
         vtr(k) = max(vmin, vr_fac*vconr*sqrt(min(10., rho0/den(k)))*exp(0.2*log(qden/normr)))
     endif
  enddo

  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo
  zt(ktop) = ze(ktop)


 do k=ktop+1,kbot
    zt(k) = ze(k) - dt5*(vtr(k-1)+vtr(k))
 enddo
 zt(kbot+1) = zs - dt*vtr(kbot)

 do k=ktop,kbot
    if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
 enddo

! Evap_acc of rain for 1/2 time step
  call revap_racc( ktop, kbot, dt5, tz, qv, ql, qr, pm, den, denfac, h_var )

  if ( ppm_rain_fall ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qr, r1, mono_prof)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qr, r1)
  endif

! Finish the remaing 1/2 time step
  call revap_racc( ktop, kbot, dt5, tz, qv, ql, qr, pm, den, denfac, h_var )

999  continue

!-------------------
! * auto-conversion
!-------------------
! Assuming linear subgrid vertical distribution of cloud water
! following Lin et al. 1994, MWR

  call linear_prof( kbot-ktop+1, p1(ktop), ql(ktop), dl(ktop), z_slope, h_var )

  qc0 = fac_rc*ccn

! * Auto conversion

  do k=ktop,kbot
    if ( tz(k) > t_wfr ) then
!----------------------------------------------------------------
!    As in Klein's GFDL AM2 stratiform scheme.
!----------------------------------------------------------------
       dl(k) = max( qrmin, dl(k) )
      q_plus = ql(k) + dl(k)
      if ( use_ccn ) then
!  CCN is formulted as CCN = CCN_surface * (den/den_surface)
           qc = qc0
      else
           qc = qc0/den(k)
      endif
      if ( q_plus > qc ) then
              sink =  dt*c_praut*den(k)
           q_minus = ql(k) - dl(k)
           if ( qc > q_minus ) then
                dq = 0.25*(q_plus-qc)**2 / dl(k)
! autoconversion rate computed using average of qc and q_plus
               sink = min(dq, sink*(q_plus-qc)/(2.*dl(k))*(0.5*(qc+q_plus))**so3)
           else                                         ! qc < q_minus
               sink = min(ql(k)-qc, sink*ql(k)**so3)
           endif
           ql(k) = ql(k) - sink
           qr(k) = qr(k) + sink
      endif
    endif
  enddo


 end subroutine warm_rain


 subroutine revap_racc( ktop, kbot, dt, tz, qv, ql, qr, pm, den, denfac, h_var )
 integer, intent(in):: ktop, kbot
 real,    intent(in):: dt                 ! time step (s)
 real,    intent(in), dimension(ktop:kbot):: pm, den, denfac
 real,    intent(in)                      :: h_var
 real, intent(inout), dimension(ktop:kbot):: tz, qv, qr, ql
! local:
 real:: qsat, dqsdt, evap, tsq, qden, q_plus, q_minus, sink
 real:: qpz, dq, dqh, tin
 integer k

  do k=ktop,kbot
   if ( tz(k) > t_wfr ) then
     if ( qr(k) > qrmin ) then
            qden = qr(k)*den(k)
             tin = tz(k) - lcp*ql(k) ! presence of clouds suppresses the rain evap
            qsat = ws1d(tin, pm(k), dqsdt)
             qpz = qv(k) + ql(k)
             dqh = h_var*max(qpz, qvmin)
         q_minus = qpz - dqh
         q_plus  = qpz + dqh

! qsat must be > q_minus to activate evaporation
! qsat must be < q_plus  to activate accretion

!-------------------
! * Rain evaporation
!-------------------
         if ( qsat > q_minus ) then
              if ( qsat > q_plus ) then
                   dq = qsat - qpz
              else
! q_minus < qsat < q_plus
! dq == dqh if qsat == q_minus
                  dq = 0.25*(q_minus-qsat)**2 / dqh
              endif
               tsq = tin*tin
              evap =  crevp(1)*tsq*dq*(crevp(2)*sqrt(qden)+crevp(3)*exp(0.725*log(qden)))   & 
                   / (crevp(4)*tsq + crevp(5)*qsat*den(k))
              evap = min( qr(k), dt*evap, dq/(1.+lcp*dqsdt) )
             qr(k) = qr(k) - evap
             qv(k) = qv(k) + evap
             tz(k) = tz(k) - evap*lcp
         endif

!-------------------
! * Accretion: pracc
!-------------------
!        if ( ql(k)>1.E-8  .and.  qsat<q_plus .and. tz(k) > tice ) then
         if ( qr(k)>qrmin .and. ql(k)>1.E-8  .and.  qsat<q_plus ) then
               sink = dt*denfac(k)*cracw * exp(0.95*log(qden))
               sink = sink/(1.+sink)*ql(k)
              ql(k) = ql(k) - sink
              qr(k) = qr(k) + sink
         endif

     endif   ! rain existed
   endif   ! warm region
  enddo

 end subroutine revap_racc


 subroutine linear_prof(km, p1,  q, dm, z_var, h_var)
! Used for cloud ice and cloud water autoconversion
! qi --> ql  & ql --> qr
! Edges: qE == qbar +/- dm
 integer, intent(in):: km
 real, intent(in ):: p1(km),  q(km)
 real, intent(out):: dm(km)
 logical, intent(in):: z_var
 real, intent(in):: h_var
!
 real:: dq(km)
 integer:: k

 if ( z_var ) then
    do k=2,km
       dq(k) = 0.5*(q(k) - q(k-1))
    enddo
    dm(1) = 0.
    do k=2, km-1
! Use twice the strength of the  positive definiteness limiter (Lin et al 1994)
       dm(k) = 0.5*min(abs(dq(k) + dq(k+1)), 0.5*q(k))
       if ( dq(k)*dq(k+1) <= 0. ) then
            if ( dq(k) > 0. ) then   ! Local max
                 dm(k) = min( dm(k), dq(k), -dq(k+1) )
            else
                 dm(k) = 0.
            endif
       endif
    enddo
    dm(km) = 0.
! impose a presumed background horizontal variability that is proportional to the value itself
    do k=1, km
       dm(k) = max( dm(k), qvmin, h_var*q(k) )
    enddo
 else
    do k=1, km
       dm(k) = max( qvmin, h_var*q(k) )
    enddo
 endif

 end subroutine linear_prof


 subroutine icloud(ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
                   den, denfac, vts, vtg, vtr, qak, h_var)

!----------------------------------------------------
! Bulk cloud micro-physics; processes splitting
! with some un-split sub-grouping
! Time implicit (when possible) accretion and autoconversion
! Author: Shian-Jiann Lin, GFDL
!-------------------------------------------------------

 integer, intent(in) :: ktop, kbot
 real, intent(in),    dimension(ktop:kbot):: p1, dp1, den, denfac, vts, vtg, vtr
 real, intent(inout), dimension(ktop:kbot):: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak
 real, intent(in) :: h_var
! local:
 real, parameter:: rhos = 0.1e3    ! snow density (1/10 of water)
 real, dimension(2*(kbot-ktop-1)):: p2, den2, tz2, qv2, ql2, qr2, qs2, qi2, qg2, qa2 
 real, dimension(ktop:kbot) :: lcpk, icpk, tcpk, di
 real :: tz, qv, ql, qr, qi, qs, qg, melt
 real :: praut, pracw, pracs, psacw, pgacw, pgmlt,   &
         psmlt, prevp, psacr, pgacr, pgfr,  pgacs,   &
         pgaut, pgaci, praci, psaut, psaci, piacr
 real :: tc, tsq, dqs0, qden, qim, qsm
 real :: factor, sink
 real :: tmp1, qsw, qsi, dqsdt, dq
 real :: n0s, lamda
 real :: qc, q_plus, q_minus
 integer :: km, kn
 integer :: i, j, k, k1

 do k=ktop,kbot
!--------------------------------------
!      tmp = cp - rdgas*ptop/p1(k)
!   lcpk(k) =  latv / tmp
!   icpk(k) =  lati / tmp
!   tcpk(k) = lcpk(k) + icpk(k)
!--------------------------------------
    lcpk(k) = lcp
    icpk(k) = icp
    tcpk(k) = tcp
 enddo

! Sources of cloud ice: pihom, cold rain, and the sat_adj
! (initiation plus deposition)

! Sources of snow: cold rain, auto conversion + accretion (from cloud ice)
! sat_adj (deposition; requires pre-existing snow); initial snow comes from auto conversion

 do k=ktop, kbot
!--------------------------------------
! * pimlt: instant melting of cloud ice
!--------------------------------------
    if( tzk(k) > tice .and. qik(k) > qcmin ) then
        melt = min( qik(k), (tzk(k)-tice)/icpk(k) )
!           qim = qi0_crt / den(k)
            qim = qi0_crt
! min rain due to melted snow autoconversion
           tmp1 = min( melt, dim(qik(k),qim) )
! limit max ql amount to no greater than snow autocon threshold
           tmp1 = min( melt-tmp1, dim(qim, qlk(k)) )
         qlk(k) = qlk(k) + tmp1
         qrk(k) = qrk(k) + melt - tmp1
         qik(k) = qik(k) - melt
         tzk(k) = tzk(k) - melt*icpk(k)
    endif
 enddo

 call linear_prof( kbot-ktop+1, p1(ktop), qik(ktop), di(ktop), .false., h_var )

 do 3000 k=ktop, kbot

   if( tzk(k) < t_min ) goto 3000

   tz = tzk(k)
   qv = qvk(k)
   ql = qlk(k)
   qi = qik(k)
   qr = qrk(k)
   qs = qsk(k)
   qg = qgk(k)

!--------------------------------------
! *** Split-micro_physics_processes ***
!--------------------------------------
! Zetac: excluded (from LFO83) term: psdep
! pgwet and realidw  removed by SJL

   pgacr = 0.
   pgacw = 0.

   tc = tz-tice
if ( tc > 0. ) then

!-----------------------------
!* Melting of snow and graupel
!-----------------------------
     dqs0 = ces0/p1(k) - qv

#ifdef MORE_SNOW_MELT
     if ( qs>qvmin ) then
! Melting of snow into rain
          factor = min( 1., tc/t_snow_melt )
            sink = min( fac_sno*qs, factor*tc/icpk(k) )
          qs = qs - sink
          qr = qr + sink
          tz = tz - sink*icpk(k)    ! cooling due to snow melting
          tc = tz-tice
     endif
#endif

     if( qs>qcmin ) then 

! * accretion: cloud water --> snow
! only rate is used (for snow melt) since tc > 0.
        if( ql>qrmin ) then
            factor = denfac(k)*csacw*exp(0.8125*log(qs*den(k)))
             psacw = factor/(1.+dts*factor)*ql     ! rate
        else
             psacw = 0.
        endif

        if ( qr>qrmin ) then
! * accretion: melted snow --> rain:
             psacr = min(acr3d(vts(k), vtr(k), qr, qs, csacr, acco(1,2), den(k)), qr*rdts)
! * accretion: snow --> rain
             pracs = acr3d(vtr(k), vts(k), qs, qr, cracs, acco(1,1), den(k))
        else
             psacr = 0.
             pracs = 0.
        endif

! Total snow sink:
! * Snow melt (due to rain accretion): snow --> rain
        psmlt = max(0., smlt(tc, dqs0, qs*den(k), psacw, psacr, csmlt, den(k), denfac(k)))
        sink = min(qs, dts*(psmlt+pracs), tc/icpk(k))

        qs = qs - sink
        qr = qr + sink
        tz = tz - sink*icpk(k)    ! cooling due to snow melting
        tc = tz-tice
     endif

     if ( qg>qcmin .and. tc>0. ) then
         if ( qr>qrmin ) then
! * accretion: rain --> graupel
              pgacr = min(acr3d(vtg(k), vtr(k), qr, qg, cgacr, acco(1,3), den(k)), rdts*qr)
         endif

         qden = qg*den(k)
         if( ql>qrmin ) then
! * accretion: cloud water --> graupel
!            factor = cgacw/sqrt(den(k))*(qg*den(k))**0.875
             factor = cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
              pgacw = factor/(1.+dts*factor) * ql  ! rate
         endif

! * melting: graupel --> rain
         pgmlt = dts*gmlt(tc, dqs0, qden, pgacw, pgacr, cgmlt, den(k))
         pgmlt = min( max(0., pgmlt), qg, tc/icpk(k) )
            qg = qg - pgmlt 
            qr = qr + pgmlt 
            tz = tz - pgmlt*icpk(k)
     endif   ! graupel existed

elseif( tc < 0.0 ) then 

!------------------
! Cloud ice proc:
!------------------
  if ( qi>1.E-8 ) then

!----------------------------------------
! * accretion (pacr): cloud ice --> snow
!----------------------------------------
     if ( qs>1.E-8 )  then
! The following is originally from the "Lin Micro-physics" in Zetac
! SJL added (following Lin Eq. 23) the temperature dependency
! To reduce accretion, use Esi = exp(0.05*tc) as in Hong et al 2004
! To increase ice/reduce snow: exp(0.025*tc)
          factor = dts*denfac(k)*csaci*exp(0.05*tc + 0.8125*log(qs*den(k)))
          psaci = factor/(1.+factor) * qi
     else
          psaci = 0.
     endif

!-------------------------------------
! * autoconversion: cloud ice --> snow
!-------------------------------------
! Similar to LFO 1983: Eq. 21 solved implicitly
! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~0.8E-4
!   qim = qi0_crt / den(k)
    qim = qi0_crt

! Assuming linear subgrid vertical distribution of cloud ice
! The mismatch computation following Lin et al. 1994, MWR
    di(k) = max( di(k), qrmin )
    q_plus = qi + di(k)
    if ( q_plus > (qim+qrmin) ) then
         if ( qim > (qi - di(k)) ) then
              dq = 0.25*(q_plus-qim)**2 / di(k)
         else
              dq = qi - qim
         endif
         factor = dts*c_psaut*exp(0.025*tc)
         psaut  = factor/(1.+factor) * dq
    else
         psaut = 0.
    endif

    sink = min( qi, psaci+psaut )
      qi = qi - sink
      qs = qs + sink

!-----------------------------------
! * accretion: cloud ice --> graupel
!-----------------------------------
    if ( qg>qrmin .and. qi>1.E-7 ) then
!        factor = dts*cgaci/sqrt(den(k))*(qg*den(k))**0.875
! SJL added exp(0.025*tc) efficiency factor
         factor = dts*cgaci/sqrt(den(k))*exp(0.025*tc + 0.875*log(qg*den(k)))
          pgaci = factor/(1.+factor)*qi
             qi = qi - pgaci
             qg = qg + pgaci
    endif

  endif  ! cloud ice existed

!----------------
! Cold-Rain proc:
!----------------
! rain to ice, snow, graupel processes:

  tc = tz-tice

  if ( qr>qrmin .and. tc < 0. ) then

! * accretion: accretion of cloud ice by rain to produce snow or graupel
! (LFO: produces snow or graupel; cloud ice sink.. via psacr & pgfr)
! ice --> snow OR graupel (due to falling rain)
! No change to qr and  tz
         if ( qi > qrmin ) then
! SJL added exp(0.025*tc) efficiency factor as follows:
!           factor = dts*denfac(k)*craci*exp(0.025*tc + 0.95*log(qr*den(k)))
            factor = dts*denfac(k)*craci*exp(0.95*log(qr*den(k)))
             praci = factor/(1.+factor)*qi
!            if ( qr > qr0_crt/den(k) ) then
             if ( qr > qr0_crt ) then
                  qg = qg + praci
             else
                  qs = qs + praci
             endif
             qi = qi - praci
         endif

! *sink* terms to qr: psacr + piacr + pgfr
! source terms to qs: psacr
! source terms to qi: piacr
! source terms to qg: pgfr

! * accretion of rain by snow
      if ( qs > 1.E-8 ) then   ! if snow exists
           psacr = dts*acr3d(vts(k), vtr(k), qr, qs, csacr, acco(1,2), den(k))
      else
           psacr = 0.
      endif

! The following added by SJL (missing from Zetac)
! * piacr: accretion of rain by cloud ice [simplified from lfo 26]
! The value of c_piacr needs to be near order(1) to have significant effect
!-------------------------------------------------------------------
! rain --> ice 
      if ( qi > qrmin ) then
         factor = dts*denfac(k)*qi * c_piacr
          piacr = factor/(1.+factor)*qr
      else
          piacr = 0.
      endif

!-------------------------------------------------------------------
! * rain freezing --> graupel
!-----------------------------------------------------------------------------------
       pgfr = dts*cgfr(1)*(exp(-cgfr(2)*tc)-1.)*(qr*den(k))**1.75/den(k)
       qden = qr*den(k)
       pgfr = dts*cgfr(1)*(exp(-cgfr(2)*tc)-1.)*qden*qden/(sqrt(sqrt(qden))*den(k))
!-----------------------------------------------------------------------------------

!--- Total sink to qr
       sink = psacr + piacr + pgfr
     factor = min( sink, qr, -0.5*tc/icpk(k) ) / max( sink, qrmin )

      psacr = factor * psacr
      piacr = factor * piacr
      pgfr  = factor * pgfr

      sink = psacr + piacr + pgfr
        tz = tz + sink*icpk(k)
        qr = qr - sink
        qs = qs + psacr
        qi = qi + piacr
        qg = qg + pgfr
        tc = tz - tice
  endif  ! qr existed

!------------------
! Cloud water sink:
!------------------
  if( ql>qcmin ) then

! * pihom * homogeneous Freezing of cloud water into cloud ice:
! This is the 1st occurance of liquid water freezing in the split MP process
! done here to prevents excessive snow production
#ifdef USE_T_WFR
      if( tz < t_00 + dt_fr ) then
        factor = min( 1., (t_00 + dt_fr - tz)/dt_fr )
          sink = min( ql*factor, (t_00+dt_fr-tz)/icpk(k) )
            ql = ql - sink
            qi = qi + sink
            tz = tz + sink*icpk(k)
      endif
#else
      if( tz < t_00 ) then
          sink = min( ql, (t_00-tz)/icpk(k) )
            ql = ql - sink
            qi = qi + sink
            tz = tz + sink*icpk(k)
      endif
#endif

! * cloud water --> Snow (requires pre-existing cloud ice)

     if( qs>1.E-8 .and. ql>1.E-8 .and. qi>1.E-8 ) then
! The following originally from Zetac: PSACW
        factor = dts*denfac(k)*csacw*exp(0.8125*log(qs*den(k)))
        psacw = min( factor/(1.+factor)*ql, -tc/icpk(k) )
        qs = qs + psacw
        ql = ql - psacw
        tz = tz + psacw*icpk(k)
     endif

  endif  ! (significant) cloud water existed

!--------------------------
! Graupel production terms:
!--------------------------

  if( qs > qrmin ) then
! * accretion: snow --> graupel
      if ( qg > qrmin ) then
           sink = dts*acr3d(vtg(k), vts(k), qs, qg, cgacs, acco(1,4), den(k))
      else
           sink = 0.
      endif

      qsm = qs0_crt / den(k)
      if ( qs > qsm ) then
! * Autoconversion Snow --> graupel
           factor = dts*1.e-3*exp(0.09*(tz-tice))
             sink = sink + factor/(1.+factor)*(qs-qsm)
      endif
      sink = min( qs, sink )
        qs = qs - sink
        qg = qg + sink

  endif   ! snow existed

  if ( qg>qrmin .and. tz < tice ) then

#ifndef PGACW_OFF
! * accretion: cloud water --> graupel
     if( ql>1.E-8 ) then
!        factor = dts*cgacw/sqrt(den(k))*(qg*den(k))**0.875
! Optimized form:
           qden = qg*den(k)
         factor = dts*cgacw*qden/sqrt(den(k)*sqrt(sqrt(qden)))
          pgacw = factor/(1.+factor)*ql
     else
          pgacw = 0.
     endif
#else
          pgacw = 0.
#endif

! * accretion: rain --> graupel
     if ( qr>qrmin ) then 
          pgacr = min(dts*acr3d(vtg(k), vtr(k), qr, qg, cgacr, acco(1,3), den(k)), qr)
     else
          pgacr = 0.
     endif

       sink = pgacr + pgacw
     factor = min( sink, 0.5*(tice-tz)/icpk(k) ) / max( sink, qrmin )
      pgacr = factor * pgacr
      pgacw = factor * pgacw

     sink = pgacr + pgacw
       tz = tz + sink*icpk(k)
       qg = qg + sink
       qr = qr - pgacr
       ql = ql - pgacw

  endif    ! graupel existed

endif   ! end ice-physics 

     tzk(k) = tz
     qvk(k) = qv
     qlk(k) = ql
     qik(k) = qi
     qrk(k) = qr
     qsk(k) = qs
     qgk(k) = qg

3000 continue   ! k-loop

 if ( do_subgrid_z ) then

! Except top 2 and bottom 2 layers (4 layers total), using subgrid PPM distribution
! to perform saturation adjustment at 2X the vertical resolution

   kn = kbot - ktop + 1
   km = 2*(kbot-ktop-1)

   p2(1) =  p1(ktop  )
   p2(2) =  p1(ktop+1)
   do k=3,km-3,2
           k1 = ktop+1 + k/2
      p2(k  ) = p1(k1) - 0.25*dp1(k1) 
      p2(k+1) = p1(k1) + 0.25*dp1(k1) 
   enddo

   if ( mp_debug ) then
     if (k1 /= (kbot-2))  then
         write(*,*) 'FATAL: k1=', k1
         call error_mesg ('LIN_CLD_MICROPHYS:', 'DO_MAP2_SAT', FATAL) 
     endif
   endif

   p2(km-1) = p1(kbot-1)
   p2(km  ) = p1(kbot)

   call remap2(ktop, kbot, kn, km, dp1, tzk, tz2, 1)
   call remap2(ktop, kbot, kn, km, dp1, qvk, qv2, 1)
   call remap2(ktop, kbot, kn, km, dp1, qlk, ql2, 1)
   call remap2(ktop, kbot, kn, km, dp1, qik, qi2, 1)
   call remap2(ktop, kbot, kn, km, dp1, qsk, qs2, 1)
   call remap2(ktop, kbot, kn, km, dp1, qgk, qg2, 1)

   do k=1,km
      den2(k) = p2(k)/(rdgas*tz2(k))
       qa2(k) = 0.
   enddo

   call subgrid_z_proc(1, km, p2, den2, h_var, tz2, qv2, ql2, qi2, qs2, qg2, qa2)
 
! Remap back to original larger volumes:
   qak(ktop  ) = qak(ktop  ) + qa2(1)
   qak(ktop+1) = qak(ktop+1) + qa2(2)
  
   tzk(ktop  ) = tz2(1)
   tzk(ktop+1) = tz2(2)

   qvk(ktop  ) = qv2(1)
   qvk(ktop+1) = qv2(2)

   qlk(ktop  ) = ql2(1)
   qlk(ktop+1) = ql2(2)

   qik(ktop  ) = qi2(1)
   qik(ktop+1) = qi2(2)

   qsk(ktop  ) = qs2(1)
   qsk(ktop+1) = qs2(2)

   qgk(ktop  ) = qg2(1)
   qgk(ktop+1) = qg2(2)

   do k=3,km-3,2
          k1  = ktop+1 + k/2
      qak(k1) = qak(k1) + max(qa2(k), qa2(k+1))  ! Maximum only
! Subgrid overlap schemes: max and random parts weighted by subgrid horizontal deviation
!-------------------------------------------------------------------------------------
! Random cloud fraction = 1 - (1-a1)*(1-a2) = a1 + a2 - a1*a2
! RAND_CLOUD
!     qak(k1) = qak(k1) + (1.-h_var)*max(qa2(k), qa2(k+1))     &  ! Maximum fraction
!                       + h_var*(qa2(k)+qa2(k+1)-qa2(k)*qa2(k+1)) ! Random  fraction
!-------------------------------------------------------------------------------------
      tzk(k1) = 0.5*(tz2(k) + tz2(k+1))
      qvk(k1) = 0.5*(qv2(k) + qv2(k+1))
      qlk(k1) = 0.5*(ql2(k) + ql2(k+1))
      qik(k1) = 0.5*(qi2(k) + qi2(k+1))
      qsk(k1) = 0.5*(qs2(k) + qs2(k+1))
      qgk(k1) = 0.5*(qg2(k) + qg2(k+1))
   enddo

   qak(kbot-1) = qak(kbot-1) + qa2(km-1)
   qak(kbot  ) = qak(kbot  ) + qa2(km  )

   tzk(kbot-1) = tz2(km-1)
   tzk(kbot  ) = tz2(km  )

   qvk(kbot-1) = qv2(km-1)
   qvk(kbot  ) = qv2(km  )

   qlk(kbot-1) = ql2(km-1)
   qlk(kbot  ) = ql2(km  )

   qik(kbot-1) = qi2(km-1)
   qik(kbot  ) = qi2(km  )

   qsk(kbot-1) = qs2(km-1)
   qsk(kbot  ) = qs2(km  )

   qgk(kbot-1) = qg2(km-1)
   qgk(kbot  ) = qg2(km  )
 else
   call subgrid_z_proc(ktop, kbot, p1, den, h_var, tzk, qvk, qlk, qik, qsk, qgk, qak)
 endif

 end subroutine icloud


 subroutine remap2(ktop, kbot, kn, km, dp, q1, q2, id)
 integer, intent(in):: ktop, kbot, kn, km , id
! constant distribution if id ==0
 real, intent(in), dimension(ktop:kbot):: q1, dp
 real, intent(out):: q2(km)
! local
 real:: a4(4,ktop:kbot)
 real:: tmp
 integer:: k, k1

  q2(1) = q1(ktop  )
  q2(2) = q1(ktop+1)

  if ( id==1 ) then

      do k=ktop,kbot
         a4(1,k) = q1(k)
      enddo
      call cs_profile( a4(1,ktop), dp(ktop), kn, mono_prof )  ! non-monotonic

      do k=3,km-3,2
              k1 = ktop+1 + k/2
         q2(k  ) = min( 2.*q1(k1), max( qvmin, a4(1,k1) + 0.25*(a4(2,k1)-a4(3,k1)) ) )
         q2(k+1) = 2.*q1(k1) - q2(k)
      enddo

  else
      do k=3,km-3,2
              k1 = ktop+1 + k/2
         q2(k  ) = q1(k1)
         q2(k+1) = q1(k1)
      enddo
  endif

  q2(km-1) = q1(kbot-1)
  q2(km  ) = q1(kbot)

 end subroutine remap2



 subroutine subgrid_z_proc(ktop, kbot, p1, den, h_var, tz, qv, ql, qi, qs, qg, qa)

! Temperature sentive high vertical resolution processes:

 integer, intent(in):: ktop, kbot
 real, intent(in),    dimension(ktop:kbot):: p1, den
 real, intent(in)                         :: h_var
 real, intent(inout), dimension(ktop:kbot):: tz, qv, ql, qi, qs, qg, qa
! local:
! qstar over water may be accurate only down to -80 C with ~10% uncertainty
! must not be too large to allow PSC
 real:: denf, rh, clouds,  rqi, tin, qsw, qsi, qpz, qstar
 real:: dqsdt, dwsdt, dq, qimin, pidep, factor, tmp, pgsub
 real:: qlv, q_plus, q_minus, qi_crt, dqh, qlt
 real:: evap, sink, qden, tc, tsq, pisub, pssub, iwt, q_adj, qq, f_frz
 integer :: k

 do 4000 k=ktop,kbot

! Quick pass check
!-----------------
   if ( tz(k) < t_min ) goto 4000

! Instant evaporation/sublimation of all clouds if RH<rh_adj --> cloud free
! This segment is the only true "saturation adjustment" in this code, and it
! operates only for low RH (set by rh_adj)

!     iwt = qi(k) + qs(k)
! Do not include snow
      iwt = qi(k)
   clouds = ql(k) + iwt

   tin = tz(k) - ( lcp*clouds + icp*iwt )  ! minimum  temperature
   qpz = qv(k) + clouds                    ! conserved within subgrid_z_proc
    rh = qpz*p1(k)/(eps*es2_table(tin))    ! 2-phase (pure ice & water)

    if ( rh<rh_adj ) then  ! qpz / rh_adj < qs
         tz(k) = tin
         qv(k) = qpz
         ql(k) = 0.
         qi(k) = 0.
!        qs(k) = 0.
         goto 4000            ! cloud free
    endif

!---------------------------------------------------------------------------------------------------
! Evaporation of cloud water from EQ A12, Fowler, Randall, and Rutledge 1996, solved exactly in time
! since dq/qsw < 1 --> evap < ql; therefore ql will never be completely evaporated within this step
! Complete evaporation happens in the pre-conditioner if conditions are met
!---------------------------------------------------------------------------------------------------
        qsw = ws1d(tz(k), p1(k), dwsdt)
       evap = 0.
         dq = qsw -  qv(k)
      q_adj = dq / (1.+lcp*dwsdt) ! maximum possible change amount to qv

      if( dq > 0.0 .and. ql(k)>qcmin ) then
!         Evaporation of ql:
          evap = min( fac_l2v*(dq/qsw)*ql(k), q_adj )
      endif
      if( dq < 0.0 ) then
!         Condensation:
          evap = max( fac_v2l*dq, q_adj )
      endif

      if ( abs(evap) > qvmin ) then
           qv(k) = qv(k) + evap
           ql(k) = ql(k) - evap
           tz(k) = tz(k) - evap*lcp
      endif

   tc = tz(k) - tice

  if ( tc < 0. ) then

! *********** freezing of cloud water ********
! -- pihom --

     if( ql(k) > qcmin ) then
! SJL, May 21, 2012
! Complete freezing below t_00 (-48 C)
! tmp = estimated maximum warming (deg K) due to freezing
           tmp = ql(k)*icp*min( 1., tc/(t_00 - tice))
           tmp = min( tmp, -tc )
        factor = tc/(t_00 - tice - tmp)
        if ( factor > 0.999 ) then
             sink = min( ql(k), (t_wfr-tz(k))/icp )
        else
! freezing time scale = tau_frz*(1-factor)/factor 
! Biggs form  ~ qq**2 * density * [exp(-0.66tc)-1]
!               qq = ql(k)/q00 
!            f_frz = 1. - exp( -dts*qq*qq*den(k)*factor/(tau_frz*(1.-factor)) )
!-------------------------------------------------------------------
             f_frz = 1. - exp( -dts*factor/(tau_frz*(1.-factor)) )
              sink = min( f_frz*ql(k), -tc/icp )
        endif
        ql(k) = ql(k) - sink
        qi(k) = qi(k) + sink
        tz(k) = tz(k) + sink*icp
     endif ! significant ql existed

!--------------------------------------------------------------
! Simplified ice <--> vapor exchange scheme
!--------------------------------------------------------------
! Fowler, Randall, and Rutledge 1996, solved exactly in time
       tc = tz(k) - tice
      qsi = qs1d(tz(k), p1(k), dqsdt)
       dq = qv(k) - qsi

    if ( dq > 0. ) then
! vapor ---> ice (does not require pre-existing ice)
         pisub = min( fac_v2i*dq, dq/(1.+tcp*dqsdt), -tc/tcp )
         qi(k) = qi(k) + pisub
         qv(k) = qv(k) - pisub
         tz(k) = tz(k) + pisub*tcp
    elseif ( qi(k) > qcmin ) then      ! qi --> qv
          pisub = -max( fac_i2v*dq*qi(k)/qsi, dq/(1.+tcp*dqsdt) )
          qv(k) = qv(k) + pisub
          qi(k) = qi(k) - pisub
          tz(k) = tz(k) - pisub*tcp
    endif

!------------------------------
! * Snow sublimation-deposition
!------------------------------
    tc = tz(k) - tice
    if ( qs(k) > qcmin  .and. tc < 0. ) then
        qsi = qs1d(tz(k), p1(k), dqsdt)
         dq = qv(k) - qsi
       sink = dq/(1.+tcp*dqsdt)
      if ( dq > 0. ) then
! vapor --> snow
         pssub = min( fac_v2s*dq*qs(k)/qsi, sink, -tc/tcp )
         qs(k) = qs(k) + pssub 
         qv(k) = qv(k) - pssub 
         tz(k) = tz(k) + pssub*tcp
      else
! snow --> vapor
! This is not sufficient to completely evaporate snow
         pssub = -max( fac_s2v*dq*qs(k)/qsi, sink )
         qv(k) = qv(k) + pssub
         qs(k) = qs(k) - pssub
         tz(k) = tz(k) - pssub*tcp
      endif
    endif       ! snow existed

!---------------------------------
! * Grapuel sublimation-deposition
!---------------------------------
    if ( qg(k) > qcmin ) then
         qsi = qs1d(tz(k), p1(k), dqsdt)
          dq = qv(k) - qsi
        sink = dq/(1.+tcp*dqsdt)
       pgsub = (qv(k)/qsi-1.) * qg(k)

       if ( pgsub > 0. ) then        ! deposition
            pgsub = min( fac_v2g*pgsub, sink, (tice-tz(k))/tcp )
       else                          ! submilation
            pgsub = max( fac_g2v*pgsub, sink )
       endif

       qg(k) = qg(k) + pgsub
       qv(k) = qv(k) - pgsub 
       tz(k) = tz(k) + pgsub*tcp
    endif    ! graupel existed

  endif   ! sub-freezing check

! use cloud condensates at true temperature to determine the water/ice partition
!  iwt = qi(k) + qs(k) 

   iwt = qi(k)
   clouds = ql(k) + iwt
   qpz = qv(k) + clouds                   
   tin = tz(k) - ( lcp*clouds + icp*iwt )  ! minimum  temperature

!--------------------
! * determine qstar 
!--------------------
! Using the "liquid-frozen water temperature": tin
   if( tin <= t_wfr ) then
       qstar = iqsat(tin, p1(k))
   elseif ( tin >= tice ) then
       qstar = wqsat(tin, p1(k))
   else
! mixed phase:
       qsi = iqsat(tin, p1(k))
       qsw = wqsat(tin, p1(k))
       if( clouds > 3.E-6 ) then
           rqi = iwt / clouds
       else
! Mostly liquid water clouds at initial cloud development stage
           rqi = (tice-tin)/(tice-t_wfr)
!          rqi = rqi ** 2    ! biased towards water phase when little condensates exist
       endif
       qstar = rqi*qsi + (1.-rqi)*qsw
   endif

!-------------------------
! * cloud fraction
!-------------------------
! Assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
! binary cloud scheme

   if ( qpz > qrmin ) then
! Partial cloudiness by PDF:
            dq = max(qcmin, h_var*qpz)
       q_plus  = qpz + dq        ! cloud free if qstar > q_plus
       q_minus = qpz - dq
       if ( qstar < q_minus ) then
            qa(k) = qa(k) + 1.       ! Air fully saturated; 100 % cloud cover
       elseif ( qstar<q_plus .and. clouds>qc_crt ) then
            qa(k) = qa(k) + (q_plus-qstar)/(dq+dq)
       endif
   endif

4000 continue

 end subroutine subgrid_z_proc



 subroutine terminal_fall(dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, pm, dz, dp,  &
                          den, vtg, vts, vti, r1, g1, s1, i1)

! lagrangian control-volume method:

 real,    intent(in):: dtm                    ! time step (s)
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp, vtg, vts, vti, pm, den
 real,    intent(inout), dimension(ktop:kbot):: dz, qv, ql, qr, qg, qs, qi, tz
 real,    intent(out):: r1, g1, s1, i1
! local:
 real, dimension(ktop:kbot+1):: ze, zt
 real:: qsat, dqsdt, dt5, melt, evap, dtime
 real:: factor, frac
 real:: tmp1, qim, precip, tc, sink
 real, dimension(ktop:kbot):: lcpk, icpk
 real:: zs = 0.
 integer k, k0, m
 logical no_fall

  do k=ktop,kbot
!       tmp1 = cp - rdgas*ptop/pm(k)
!    lcpk(k) = latv / tmp1
!    icpk(k) = lati / tmp1
     lcpk(k) = lcp
     icpk(k) = icp
  enddo

  dt5 = 0.5*dtm


! Melting of cloud_ice and snow (before fall):

! find significant melting level
  k0 = kbot
  do k=ktop, kbot-1
     if ( tz(k) > tice ) then
          k0 = k
          go to 11
     endif
  enddo
11  continue

 do k=k0, kbot
!------
! * ice
!------
    tc = tz(k) - tice
    if( qi(k) > qcmin .and. tc>0. ) then
        melt = min( fac_mlt*qi(k), tc/icpk(k) )
!       qim = qi0_crt / den(k)
        qim = qi0_crt
! min rain due to melted snow autoconversion
        tmp1 = min( melt, dim(qi(k),qim) )
! limit max ql amount to no greater than ice-->snow autocon threshold
        tmp1 = min( melt-tmp1, dim(qim, ql(k)) )
        ql(k) = ql(k) + tmp1
        qr(k) = qr(k) + melt - tmp1
        qi(k) = qi(k) - melt
        tz(k) = tz(k) - melt*icpk(k)
           tc = tz(k) - tice
    endif
!------
! Snow
!------
    if ( qs(k)>qvmin .and. tc>0. ) then
         factor = min( 1., tc/t_snow_melt)
           sink = min(qs(k)*fac_sno, factor*tc/icpk(k))
         qs(k) = qs(k) - sink
         qr(k) = qr(k) + sink
         tz(k) = tz(k) - sink*icpk(k)    ! cooling due to snow melting
    endif

 enddo

  if ( dts < 75. ) k0 = kbot

!-----
! ice:
!-----

  ze(kbot+1) = zs
  do k=kbot, ktop, -1
     ze(k) = ze(k+1) - dz(k)  ! dz<0
  enddo

  zt(ktop) = ze(ktop)

  call check_column(ktop, kbot, qi, no_fall)

  if ( vi_fac < 1.e-5 .or. no_fall ) then
     i1 = 0.
  else

  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vti(k-1)+vti(k))
  enddo
  zt(kbot+1) = zs - dtm*vti(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qi(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min( 1.0, (ze(m)-ze(m+1))/(max(vmin,vti(k))*tau_mlt) )
                   melt = min( qi(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m) )
!                   qim = qi0_crt / den(m)
                    qim = qi0_crt
                   tmp1 = min( melt, dim(qi(k), qim) )      ! min rain (snow autoconversion)
                   tmp1 = min( melt-tmp1, dim(qim, ql(m)) ) ! limit max ql amount
!
                  ql(m) = ql(m) + tmp1
                  qr(m) = qr(m) - tmp1 + melt
                  tz(m) = tz(m) - melt*icpk(m)
                  qi(k) = qi(k) - melt*dp(m)/dp(k)
             endif
          enddo
     endif
  enddo
  endif

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qi, i1, mono_prof)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qi, i1)
  endif

  endif

!--------------------------------------------
! melting of falling snow (qs) into rain(qr)
!--------------------------------------------
  r1 = 0.

  call check_column(ktop, kbot, qs, no_fall)

  if ( no_fall ) then
       s1 = 0.
  else

  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vts(k-1)+vts(k))
  enddo
  zt(kbot+1) = zs - dtm*vts(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qs(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
                  dtime = min( dtm, (ze(m)-ze(m+1))/(vmin+vts(k)) )
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min(1., dtime/tau_s)
                   melt = min(qs(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m))
                  tz(m) = tz(m) - melt*icpk(m)
                  qs(k) = qs(k) - melt*dp(m)/dp(k)
                  if ( zt(k)<zs ) then
                       r1 = r1 + melt*dp(m)   ! precip as rain
                  else
!                      qr source here will fall next time step (therefore, can evap)
                       qr(m) = qr(m) + melt
                  endif
             endif
             if ( qs(k) < qrmin ) exit
          enddo
     endif
  enddo
  endif

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qs, s1, mono_prof)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qs, s1)
  endif
  endif

!----------------------------------------------
! melting of falling graupel (qg) into rain(qr)
!----------------------------------------------
  call check_column(ktop, kbot, qg, no_fall)

  if ( no_fall ) then
       g1 = 0.
  else
  do k=ktop+1,kbot
     zt(k) = ze(k) - dt5*(vtg(k-1)+vtg(k))
  enddo
  zt(kbot+1) = zs - dtm*vtg(kbot)

  do k=ktop,kbot
     if( zt(k+1)>=zt(k) ) zt(k+1) = zt(k) - dz_min
  enddo

  if ( k0 < kbot ) then
  do k=kbot-1,k0,-1
     if ( qg(k) > qrmin ) then
          do m=k+1, kbot
             if ( zt(k+1)>=ze(m) ) exit
             dtime = min( dtm, (ze(m)-ze(m+1))/vtg(k) )
             if ( zt(k)<ze(m+1) .and. tz(m)>tice ) then
                  dtime = min(1., dtime/tau_g)
                   melt = min(qg(k)*dp(k)/dp(m), dtime*(tz(m)-tice)/icpk(m))
                  tz(m) = tz(m) - melt*icpk(m)
                  qg(k) = qg(k) -  melt*dp(m)/dp(k)
                  if ( zt(k)<zs ) then
                       r1 = r1 + melt*dp(m)
                  else
                       qr(m) = qr(m) + melt
                  endif
             endif
             if ( qg(k) < qrmin ) exit
           enddo
     endif
  enddo
  endif

  if ( use_ppm ) then
       call lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, qg, g1, mono_prof)
  else
       call lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, qg, g1)
  endif
  endif


 end subroutine terminal_fall


 subroutine check_column(ktop, kbot, q, no_fall)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: q(ktop:kbot)
 logical, intent(out):: no_fall
! local:
 integer k

 no_fall = .true.
 do k=ktop, kbot
    if ( q(k) > qrmin ) then
         no_fall = .false.
         exit
    endif
 enddo

 end subroutine check_column


 subroutine lagrangian_fall_pcm(ktop, kbot, zs, ze, zt, dp, q, precip)
 real,    intent(in):: zs
 integer, intent(in):: ktop, kbot
 real,    intent(in), dimension(ktop:kbot):: dp
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(inout), dimension(ktop:kbot):: q
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: qm1, qm2
 integer k, k0, n, m

! density:
  do k=ktop,kbot
     qm1(k) = q(k)*dp(k) / (zt(k)-zt(k+1))
     qm2(k) = 0.
  enddo

   k0 = ktop
   do k=ktop,kbot
      do n=k0,kbot
      if(ze(k) <= zt(n) .and. ze(k) >= zt(n+1)) then
         if(ze(k+1) >= zt(n+1)) then
!                          entire new grid is within the original grid
            qm2(k) = qm1(n)*(ze(k)-ze(k+1))
            k0 = n
            goto 555
         else
            qm2(k) = qm1(n)*(ze(k)-zt(n+1))    ! fractional area
            do m=n+1,kbot
!                                        locate the bottom edge: ze(k+1)
               if(ze(k+1) < zt(m+1) ) then
                  qm2(k) = qm2(k) + q(m)*dp(m)
               else
                  qm2(k) = qm2(k) + qm1(m)*(zt(m)-ze(k+1))
                  k0 = m
                  goto 555
               endif
            enddo
            goto 555
         endif
      endif
      enddo
555 continue
   enddo

     precip = 0.
! direct algorithm (prevent small negatives)
     do k=ktop,kbot
        if ( zt(k+1) < zs ) then
             precip = qm1(k)*(zs-zt(k+1)) 
             if ( (k+1) > kbot ) goto 777
                  do m=k+1,kbot
                     precip = precip + q(m)*dp(m)
                  enddo
             goto 777
        endif
     enddo
777  continue

   do k=ktop,kbot
      q(k) = qm2(k) / dp(k)
   enddo

 end subroutine lagrangian_fall_pcm



 subroutine lagrangian_fall_ppm(ktop, kbot, zs, ze, zt, dp, q, precip, mono)
 integer, intent(in):: ktop, kbot
 real,    intent(in):: zs
 logical, intent(in):: mono
 real,    intent(in), dimension(ktop:kbot):: dp
 real,    intent(in), dimension(ktop:kbot+1):: ze, zt
 real,    intent(inout), dimension(ktop:kbot):: q
 real,    intent(out):: precip
! local:
 real, dimension(ktop:kbot):: qm0, qm1, qm2, dz
 real a4(4,ktop:kbot)
 real pl, pr, delz, esl
 integer k, k0, n, m
 real, parameter:: r3 = 1./3., r23 = 2./3.

! density:
  do k=ktop,kbot
      dz(k) = zt(k) - zt(k+1)      ! note: dz is positive
     qm0(k) = q(k)*dp(k)
     qm1(k) = qm0(k) / dz(k)
     qm2(k) = 0.
     a4(1,k) = qm1(k)
  enddo

! Construct qm1 profile with zt as coordinate

   call cs_profile(a4(1,ktop), dz(ktop), kbot-ktop+1, mono)

   k0 = ktop
   do k=ktop,kbot
      do n=k0,kbot
      if(ze(k) <= zt(n) .and. ze(k) >= zt(n+1)) then
         pl = (zt(n)-ze(k)) / dz(n)
         if( zt(n+1) <= ze(k+1) ) then
!                          entire new grid is within the original grid
                pr = (zt(n)-ze(k+1)) / dz(n)
            qm2(k) = a4(2,n) + 0.5*(a4(4,n)+a4(3,n)-a4(2,n))*(pr+pl) -  &
                     a4(4,n)*r3*(pr*(pr+pl)+pl**2)
            qm2(k) = qm2(k)*(ze(k)-ze(k+1))
            k0 = n
            goto 555
         else
            qm2(k) = (ze(k)-zt(n+1)) * (a4(2,n)+0.5*(a4(4,n)+   &
                      a4(3,n)-a4(2,n))*(1.+pl) - a4(4,n)*( r3*(1.+pl*(1.+pl))) )
            if ( n<kbot ) then
               do m=n+1,kbot
!                                        locate the bottom edge: ze(k+1)
                  if( ze(k+1) < zt(m+1) ) then
                     qm2(k) = qm2(k) + q(m)*dp(m)
                  else
                     delz = zt(m) - ze(k+1)
                      esl = delz / dz(m)
                     qm2(k) = qm2(k) + delz*( a4(2,m) + 0.5*esl*        &
                             (a4(3,m)-a4(2,m)+a4(4,m)*(1.-r23*esl)) )
                     k0 = m
                     goto 555
                  endif
               enddo
            endif
            goto 555
         endif
      endif
      enddo
555 continue
   enddo

   precip = 0.

   do k=ktop,kbot
      precip = precip + qm0(k) - qm2(k)
   enddo
!  precip = max(0., precip)

   do k=ktop,kbot
      q(k) = qm2(k) / dp(k)
   enddo

 end subroutine lagrangian_fall_ppm


 subroutine cs_profile(a4, del, km, do_mono)
 integer, intent(in):: km      ! vertical dimension
 real   , intent(in):: del(km)
 logical, intent(in):: do_mono
 real , intent(inout):: a4(4,km)
!-----------------------------------------------------------------------
 real  gam(km)
 real  q(km+1)
 real   d4, bet, a_bot, grat, pmp, lac
 real   pmp_1, lac_1, pmp_2, lac_2
 real  da1, da2, a6da
 integer k
 logical extm(km)

     grat = del(2) / del(1)   ! grid ratio
      bet = grat*(grat+0.5)
     q(1) = (2.*grat*(grat+1.)*a4(1,1)+a4(1,2)) / bet
   gam(1) = ( 1. + grat*(grat+1.5) ) / bet

  do k=2,km
      d4 = del(k-1) / del(k)
     bet =  2. + 2.*d4 - gam(k-1)
     q(k) = (3.*(a4(1,k-1)+d4*a4(1,k))-q(k-1))/bet
     gam(k) = d4 / bet
  enddo
 
       a_bot = 1. + d4*(d4+1.5)
     q(km+1) = (2.*d4*(d4+1.)*a4(1,km)+a4(1,km-1)-a_bot*q(km))  &
             / ( d4*(d4+0.5) - a_bot*gam(km) )

  do k=km,1,-1
     q(k) = q(k) - gam(k)*q(k+1)
  enddo

!------------------
! Apply constraints
!------------------
  do k=2,km
     gam(k) = a4(1,k) - a4(1,k-1)
  enddo

! Apply large-scale constraints to ALL fields if not local max/min

! Top:
  q(1) = max( q(1), 0. )
  q(2) = min( q(2), max(a4(1,1), a4(1,2)) )
  q(2) = max( q(2), min(a4(1,1), a4(1,2)), 0. )

! Interior:
  do k=3,km-1
     if ( gam(k-1)*gam(k+1)>0. ) then
          q(k) = min( q(k), max(a4(1,k-1),a4(1,k)) )
          q(k) = max( q(k), min(a4(1,k-1),a4(1,k)) )
     else
          if ( gam(k-1) > 0. ) then
! There exists a local max                                                                             
               q(k) = max( q(k), min(a4(1,k-1),a4(1,k)) )     
          else
! There exists a local min
               q(k) = min( q(k), max(a4(1,k-1),a4(1,k)) )
               q(k) = max( q(k), 0.0 )
          endif
     endif
  enddo

  q(km  ) = min( q(km), max(a4(1,km-1), a4(1,km)) )
  q(km  ) = max( q(km), min(a4(1,km-1), a4(1,km)), 0. )
! q(km+1) = max( q(km+1), 0.)

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
  do k=1,km-1
     a4(2,k) = q(k  )
     a4(3,k) = q(k+1)
  enddo

  do k=2,km-1
     if ( gam(k)*gam(k+1) > 0.0 ) then
          extm(k) = .false.
     else
          extm(k) = .true.
     endif
  enddo

  if ( do_mono ) then
     do k=3,km-2
        if ( extm(k) ) then
! positive definite constraint ONLY if true local extrema
           if ( extm(k-1)  .or.  extm(k+1) ) then
               a4(2,k) = a4(1,k)
               a4(3,k) = a4(1,k)
           endif
        else
           a4(4,k) = 6.*a4(1,k) - 3.*(a4(2,k)+a4(3,k))
           if( abs(a4(4,k)) > abs(a4(2,k)-a4(3,k)) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                pmp_1 = a4(1,k) - 2.0*gam(k+1)
                lac_1 = pmp_1   + 1.5*gam(k+2)
              a4(2,k) = min( max(a4(2,k), min(a4(1,k), pmp_1, lac_1)),  &
                                          max(a4(1,k), pmp_1, lac_1) )
                pmp_2 = a4(1,k) + 2.0*gam(k)
                lac_2 = pmp_2   - 1.5*gam(k-1)
              a4(3,k) = min( max(a4(3,k), min(a4(1,k), pmp_2, lac_2)),  &
                                          max(a4(1,k), pmp_2, lac_2) )
           endif
        endif
     enddo
  else
     do k=3,km-2
        if ( extm(k) .and. (extm(k-1) .or. extm(k+1)) ) then
             a4(2,k) = a4(1,k)
             a4(3,k) = a4(1,k)
        endif
     enddo
  endif

  do k=1,km-1
     a4(4,k) = 6.*a4(1,k) - 3.*(a4(2,k)+a4(3,k))
  enddo

  k = km-1
  if( extm(k) ) then
      a4(2,k) = a4(1,k)
      a4(3,k) = a4(1,k)
      a4(4,k) = 0.
  else
      da1  = a4(3,k) - a4(2,k)
      da2  = da1**2
      a6da = a4(4,k)*da1
      if(a6da < -da2) then
         a4(4,k) = 3.*(a4(2,k)-a4(1,k))
         a4(3,k) = a4(2,k) - a4(4,k)
      elseif(a6da > da2) then
         a4(4,k) = 3.*(a4(3,k)-a4(1,k))
         a4(2,k) = a4(3,k) - a4(4,k)
      endif
  endif

  call cs_limiters(km-1, a4)

! Bottom layer:
  a4(2,km) = a4(1,km)
  a4(3,km) = a4(1,km)
  a4(4,km) = 0.

 end subroutine cs_profile



 subroutine cs_limiters(km, a4)
 integer, intent(in) :: km
 real, intent(inout) :: a4(4,km)   ! PPM array
! !LOCAL VARIABLES:
 real, parameter:: r12 = 1./12.
 integer k

! Positive definite constraint

 do k=1,km
 if( abs(a4(3,k)-a4(2,k)) < -a4(4,k) ) then
     if( (a4(1,k)+0.25*(a4(3,k)-a4(2,k))**2/a4(4,k)+a4(4,k)*r12) < 0. ) then
         if( a4(1,k)<a4(3,k) .and. a4(1,k)<a4(2,k) ) then
             a4(3,k) = a4(1,k)
             a4(2,k) = a4(1,k)
             a4(4,k) = 0.
         elseif( a4(3,k) > a4(2,k) ) then
             a4(4,k) = 3.*(a4(2,k)-a4(1,k))
             a4(3,k) = a4(2,k) - a4(4,k)
         else
             a4(4,k) = 3.*(a4(3,k)-a4(1,k))
             a4(2,k) = a4(3,k) - a4(4,k)
         endif
     endif
 endif
 enddo

 end subroutine cs_limiters



 subroutine fall_speed(ktop, kbot, den, qs, qi, qg, ql, tk, vts, vti, vtg)
 integer, intent(in)                     :: ktop, kbot
 real, intent(in ), dimension(ktop:kbot) :: den, qs, qi, qg, ql, tk
 real, intent(out), dimension(ktop:kbot) :: vts, vti, vtg
! fall velocity constants:
 real, parameter :: thi = 1.0e-9   ! cloud ice threshold for terminal fall
 real, parameter :: thg = 1.0e-9
 real, parameter :: ths = 1.0e-9
 real, parameter :: vf_min = 1.0E-6
 real, parameter :: vs_max = 7.        ! max fall speed for snow
!-----------------------------------------------------------------------
! marshall-palmer constants
!-----------------------------------------------------------------------
 real :: vcons = 6.6280504, vcong = 87.2382675, vconi = 3.29
 real :: norms = 942477796.076938, &
         normg =  5026548245.74367
 real, dimension(ktop:kbot) :: ri, qden, tc
!real :: aa = -1.70704e-5, bb = -0.00319109, cc = -0.0169876, dd = 0.00410839, ee = 1.93644
 real :: aa = -4.14122e-5, bb = -0.00538922, cc = -0.0516344, dd = 0.00216078, ee = 1.9714 

 real :: rhof, rho0
 integer:: k
!-----------------------------------------------------------------------
! marshall-palmer formula
!-----------------------------------------------------------------------

! try the local air density -- for global model; the true value could be
! much smaller than sfcrho over high mountains

  if ( den_ref < 0. ) then
       rho0 = -den_ref*den(kbot) 
  else
       rho0 = den_ref   ! default=1.2
  endif

   do k=ktop, kbot
        rhof = sqrt( min(100., rho0/den(k)) )
! snow:
      if ( qs(k) < ths ) then
           vts(k) = vf_min
      else
           vts(k) = max(vf_min, vcons*rhof*exp(0.0625*log(qs(k)*den(k)/norms)))
!--------------------------------------------------------------------------------------
! What if ql == 0  (ri---> 0?)
!           ri(k) = 1./(1. + 6.e-5/(max(qcmin,ql(k)) * den(k)**1.235 * qs(k)**0.235))  !--- riming intensity
!          vts(k) = max(vf_min, vconi*rhof*exp( 0.16*log((1.0-ri(k))*qs(k)*den(k)) ) +   &
!                                19.3*rhof*exp( 0.37*log(     ri(k) *qs(k)*den(k)) ) )
!--------------------------------------------------------------------------------------
           vts(k) = min( vs_max, vs_fac*vts(k) )
      endif 

! graupel:
      if ( qg(k) < thg ) then
           vtg(k) = vf_min
      else
           vtg(k) = max(vf_min, max(vmin, vg_fac*vcong*rhof*sqrt(sqrt(sqrt(qg(k)*den(k)/normg)))))
      endif
   enddo

! ice:
   if ( use_deng_mace ) then
! ice use Deng and Mace (2008, GRL), which gives smaller fall speed than HD90 formula
       do k=ktop, kbot
          if ( qi(k) < thi ) then
               vti(k) = vf_min
          else
           qden(k) = log10( 1000.*qi(k)*den(k) )   !--- used in DM formula, in g/m^-3
             tc(k) = tk(k) - tice
            vti(k) = qden(k)*( tc(k)*(aa*tc(k) + bb) + cc ) + dd*tc(k) + ee
            vti(k) = max( vf_min, vi_fac*0.01*10.**vti(k) )
          endif
       enddo
   else
! HD90 ice speed:
       do k=ktop, kbot
          if ( qi(k) < thi ) then
               vti(k) = vf_min
          else
                 rhof = sqrt( min(100., rho0/den(k)) )
               vti(k) = max( vf_min, vconi*rhof*exp(0.16*log(qi(k)*den(k))) )
          endif
       enddo
   endif

 end subroutine fall_speed


 subroutine setupm

 real :: gcon, cd, scm3, pisq, act(8), acc(3)
 real :: vdifu, tcond
 real :: visk
 real :: ch2o, hltf
 real ::  hlts, hltc, ri50

 real :: gam263, gam275, gam290,                                &
         gam325, gam350, gam380,                                &
         gam425, gam450, gam480,                                &
         gam625, gam680

 data  gam263/1.456943/,   gam275/1.608355/,  gam290/1.827363/  &
       gam325/2.54925/,    gam350/3.323363/,  gam380/4.694155/  &
       gam425/8.285063/,   gam450/11.631769/, gam480/17.837789/ &
       gam625/184.860962/, gam680/496.604067/
!
!     physical constants (mks)
!     lin's constants(mks) except rmi50,rmi40 (cgs)
!
 real :: rnzr, rnzs, rnzg, rhos, rhog
!data alin, clin  /842.0, 4.80/
 data rnzr /8.0e6/  ! lin83
 data rnzs /3.0e6/  ! lin83
 data rnzg /4.0e6/  ! rh84
 data rhos /0.1e3/  ! lin83    (snow density; 1/10 of water)
 data rhog /0.4e3/  ! rh84     (graupel density)
 data acc/5.0,2.0,0.5/

 real den_rc
 integer :: k, i

      pie = 4.*atan(1.0)

! S. Klein's formular (EQ 16) from AM2
      fac_rc = (4./3.)*pie*rhor*rthresh**3
      den_rc = fac_rc * ccn_o*1.e6
      if(master) write(*,*) 'MP: rthresh=', rthresh, 'vi_fac=', vi_fac
      if(master) write(*,*) 'MP: for ccn_o=', ccn_o, 'ql_rc=', den_rc
      den_rc = fac_rc * ccn_l*1.e6
      if(master) write(*,*) 'MP: for ccn_l=', ccn_l, 'ql_rc=', den_rc

      vdifu=2.11e-5
      tcond=2.36e-2

      visk=1.259e-5
      hlts=2.8336e6
      hltc=2.5e6
      hltf=3.336e5

      ch2o=4.1855e3
      rmi50=3.84e-6      ! Purdue Lin scheme 4.8e-7 [g]
!     rmi40=2.46e-7
      ri50=1.e-4

      pisq = pie*pie
      scm3 = (visk/vdifu)**(1./3.)
!
      cracs = pisq*rnzr*rnzs*rhos
      csacr = pisq*rnzr*rnzs*rhor
      cgacr = pisq*rnzr*rnzg*rhor
      cgacs = pisq*rnzg*rnzs*rhos
      cgacs = cgacs*c_pgacs
!
!     act:  1-2:racs(s-r); 3-4:sacr(r-s);
!           5-6:gacr(r-g); 7-8:gacs(s-g)
!
      act(1) = pie * rnzs * rhos
      act(2) = pie * rnzr * rhor
      act(6) = pie * rnzg * rhog
      act(3) = act(2)
      act(4) = act(1)
      act(5) = act(2)
      act(7) = act(1)
      act(8) = act(6)

      do i=1,3
         do k=1,4
            acco(i,k) = acc(i)/(act(2*k-1)**((7-i)*0.25)*act(2*k)**(i*0.25))
         enddo
      enddo
!
      gcon  = 40.74 * sqrt( sfcrho )   ! 44.628
!
      csacw = pie*rnzs*clin*gam325/(4.*act(1)**0.8125)
! Decreasing  csacw to reduce cloud water ---> snow

      craci = pie*rnzr*alin*gam380/(4.*act(2)**0.95)
      csaci = csacw * c_psaci
!
      cgacw = pie*rnzg*gam350*gcon/(4.*act(6)**0.875)
      cgaci = cgacw*0.1
!
      cracw = craci            ! cracw= 3.27206196043822
      cracw = c_cracw * cracw
!
!     subl and revp:  five constants for three separate processes
!
      cssub(1) = 2.*pie*vdifu*tcond*rvgas*rnzs
      cgsub(1) = 2.*pie*vdifu*tcond*rvgas*rnzg
      crevp(1) = 2.*pie*vdifu*tcond*rvgas*rnzr
      cssub(2) = 0.78/sqrt(act(1))
      cgsub(2) = 0.78/sqrt(act(6))
      crevp(2) = 0.78/sqrt(act(2))
      cssub(3) = 0.31*scm3*gam263*sqrt(clin/visk)/act(1)**0.65625
      cgsub(3) = 0.31*scm3*gam275*sqrt(gcon/visk)/act(6)**0.6875
      crevp(3) = 0.31*scm3*gam290*sqrt(alin/visk)/act(2)**0.725
      cssub(4) = tcond*rvgas
      cssub(5) = hlts**2*vdifu
      cgsub(4) = cssub(4)
      crevp(4) = cssub(4)
      cgsub(5) = cssub(5)
      crevp(5) = hltc**2*vdifu
!
      cgfr(1) = 20.e2*pisq*rnzr*rhor/act(2)**1.75
      cgfr(2) = 0.66
!
!sk ********************************************************************
!sk   smlt:  five constants ( lin et al. 1983 )
      csmlt(1) = 2.*pie*tcond*rnzs/hltf
      csmlt(2) = 2.*pie*vdifu*rnzs*hltc/hltf
      csmlt(3) = cssub(2)
      csmlt(4) = cssub(3)
      csmlt(5) = ch2o/hltf
!sk ********************************************************************
!     gmlt:  five constants
      cgmlt(1) = 2.*pie*tcond*rnzg/hltf
      cgmlt(2) = 2.*pie*vdifu*rnzg*hltc/hltf
      cgmlt(3) = cgsub(2)
      cgmlt(4) = cgsub(3)
      cgmlt(5) = ch2o/hltf
!sk ********************************************************************
      es0 = 6.107799961e2   ! ~6.1 mb
      ces0 = eps*es0
!
!     c2brg has conversion factor of 10**3
      c1brg = dts/rmi50
!lin  c2brg = ri50**2*1.e3 ! error
      c2brg = pie*ri50**2*1.e3

 end subroutine setupm


 subroutine lin_cld_microphys_init(id, jd, kd, axes, time)
 
    integer,         intent(in) :: id, jd, kd
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time
    
    integer   :: unit, io, ierr, k, logunit
    integer   :: is, ie, js, je, ks, ke
    logical   :: flag
    real :: tmp, q1, q2

    master = (mpp_pe().eq.mpp_root_pe())

#ifdef INTERNAL_FILE_NML
    read( input_nml_file, nml = lin_cld_microphys_nml, iostat = io )
    ierr = check_nml_error(io,'lin_cloud_microphys_nml')
#else
    if( file_exist( 'input.nml' ) ) then
       unit = open_namelist_file ()
       io = 1
       do while ( io .ne. 0 )
          read( unit, nml = lin_cld_microphys_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'lin_cloud_microphys_nml')
       end do
10     call close_file ( unit )
    end if
#endif
    call write_version_number (version, tagname)
    logunit = stdlog()
    
    if ( do_setup ) then
      is = 1
      js = 1
      ks = 1
      ie = id
      je = jd
      ke = kd

      call setup_con (is, ie, js, je, ks, ke)
      call setupm
      do_setup = .false.
    endif

    if (master) write( logunit, nml = lin_cld_microphys_nml )
 
    id_vtr = register_diag_field ( mod_name, 'vt_r', axes(1:3), time,        &
         'rain fall speed', 'm/sec', missing_value=missing_value )
    id_vts = register_diag_field ( mod_name, 'vt_s', axes(1:3), time,        &
         'snow fall speed', 'm/sec', missing_value=missing_value )
    id_vtg = register_diag_field ( mod_name, 'vt_g', axes(1:3), time,        &
         'graupel fall speed', 'm/sec', missing_value=missing_value )
    id_vti = register_diag_field ( mod_name, 'vt_i', axes(1:3), time,        &
         'ice fall speed', 'm/sec', missing_value=missing_value )

    id_rain = register_diag_field ( mod_name, 'rain_lin', axes(1:2), time,        &
         'rain_lin', 'mm/day', missing_value=missing_value )
    id_snow = register_diag_field ( mod_name, 'snow_lin', axes(1:2), time,        &
         'snow_lin', 'mm/day', missing_value=missing_value )
    id_graupel = register_diag_field ( mod_name, 'graupel_lin', axes(1:2), time,  &
         'graupel_lin', 'mm/day', missing_value=missing_value )
    id_ice = register_diag_field ( mod_name, 'ice_lin', axes(1:2), time,        &
         'ice_lin', 'mm/day', missing_value=missing_value )
    id_prec = register_diag_field ( mod_name, 'prec_lin', axes(1:2), time,     &
         'prec_lin', 'mm/day', missing_value=missing_value )
!   if ( master ) write(*,*) 'prec_lin diagnostics initialized.', id_prec

    id_cond = register_diag_field ( mod_name, 'cond_lin', axes(1:2), time,     &
         'total condensate', 'kg/m**2', missing_value=missing_value )

    id_var = register_diag_field ( mod_name, 'var_lin', axes(1:2), time,     &
         'subgrid variance', 'n/a',  missing_value=missing_value )

!------------------------
! fall speed diagnostics:
!------------------------
      if ( id_vtr> 0 ) then
           allocate ( vt_r(is:ie, js:je, ks:ke) )
           vt_r = 0.
      endif
      if ( id_vts> 0 ) then
           allocate ( vt_s(is:ie, js:je, ks:ke) )
           vt_s = 0.
      endif
      if ( id_vtg> 0 ) then
           allocate ( vt_g(is:ie, js:je, ks:ke) )
           vt_g = 0.
      endif
      if ( id_vti> 0 ) then
           allocate ( vt_i(is:ie, js:je, ks:ke) )
           vt_i = 0.
      endif
      if ( id_var>0 ) then
           allocate ( w_var(is:ie, js:je) )
           w_var = 0.
      endif

      allocate (     cond(is:ie, js:je) )
      allocate (  prec_mp(is:ie, js:je) )
      allocate (    prec0(is:ie, js:je) )
      allocate (    prec1(is:ie, js:je) )
      allocate (    rain0(is:ie, js:je) )
      allocate (    snow0(is:ie, js:je) )
      allocate (     ice0(is:ie, js:je) )
      allocate ( graupel0(is:ie, js:je) )

       cond = 0.
      prec0 = 0.
      prec1 = 0.
      rain0 = 0.
      snow0 = 0.
       ice0 = 0.
   graupel0 = 0.

!   call qsmith_init

! TESTING the water vapor tables
   if ( mp_debug .and. master ) then
        write(*,*) 'TESTING water vapor tables in lin_cld_microphys'
        tmp = tice - 90.
   do k=1,25
      q1 = wqsat(tmp, 1.E5)
      q2 = iqsat(tmp, 1.E5)
      write(*,*) NINT(tmp-tice), q1, q2, 'dq=', q1-q2
      tmp = tmp + 5.
   enddo
   endif

   if ( master ) write(*,*) 'lin_cld_micrphys diagnostics initialized.'

   lin_cld_mp_clock = mpp_clock_id('Lin_cld_microphys', grain=CLOCK_ROUTINE)
   g_sum_initialized = .false.
   module_is_initialized = .true.
    
 end subroutine lin_cld_microphys_init



 subroutine lin_cld_microphys_end
   real gmp

  if ( mp_print ) then
! the g_sum call does not work if physics window is used *****
   if ( id_ice> 0 ) then
        gmp = g_sum(ice0, isc, iec, jsc, jec, ng, l_area, 1) 
        if(master) write(*,*) 'total ice=', gmp/mp_count
   endif
   if ( id_graupel> 0 ) then
        gmp = g_sum(graupel0, isc, iec, jsc, jec, ng, l_area, 1) 
        if(master) write(*,*) 'total graupel=', gmp/mp_count
   endif
   if ( id_snow> 0 ) then
        gmp = g_sum(snow0, isc, iec, jsc, jec, ng, l_area, 1) 
        if(master) write(*,*) 'total snow=', gmp/mp_count
   endif
   if ( id_rain> 0 ) then
        gmp = g_sum(rain0, isc, iec, jsc, jec, ng, l_area, 1) 
        if(master) write(*,*) 'total rain=', gmp/mp_count
   endif
!  if ( id_prec> 0 ) then
        gmp = g_sum(prec0, isc, iec, jsc, jec, ng, l_area, 1) 
        if(master) write(*,*) 'total prec=', gmp/mp_count
!  endif
  endif

   if ( id_vtr> 0 ) deallocate ( vt_r )
   if ( id_vts> 0 ) deallocate ( vt_s )
   if ( id_vti> 0 ) deallocate ( vt_i )
   if ( id_vtg> 0 ) deallocate ( vt_g )
   if ( id_var> 0 ) deallocate ( w_var )

   deallocate (  prec_mp  )
   deallocate (  prec0    )
   deallocate (  prec1    )
   deallocate (  rain0    )
   deallocate (  snow0    )
   deallocate (  ice0     )
   deallocate (  graupel0 )
   deallocate (  cond )

   deallocate ( table  )
   deallocate ( table2 )
   deallocate ( table3 )
   deallocate ( tablew )
   deallocate ( des )
   deallocate ( des2 )
   deallocate ( des3 )
   deallocate ( desw )
    
 end subroutine lin_cld_microphys_end



 subroutine setup_con( is, ie, js, je, ks, ke )
 integer, intent(in) :: is,ie, js,je, ks, ke

  master = (mpp_pe().eq.mpp_root_pe())

  isc = is;   iec = ie
  jsc = js;   jec = je

  lcp = latv / cp
  icp = lati / cp
  tcp = (latv+lati) / cp

  rgrav = 1./ grav

  call qsmith_init


 end subroutine setup_con



 real function acr3d(v1, v2, q1, q2, c, cac, rho)
 real, intent(in) :: v1, v2, c, rho
 real, intent(in) :: q1, q2    ! mixing ratio!!!
 real, intent(in) :: cac(3)
 real :: t1, s1, s2
!integer :: k
! real:: a
!     a=0.0
!     do k=1,3
!        a = a + cac(k)*( (q1*rho)**((7-k)*0.25) * (q2*rho)**(k*0.25) )
!     enddo
!     acr3d = c * abs(v1-v2) * a/rho
!----------
! Optimized
!----------
      t1 = sqrt(q1*rho)
      s1 = sqrt(q2*rho)
      s2 = sqrt(s1)       ! s1 = s2**2
      acr3d = c*abs(v1-v2)*q1*s2*(cac(1)*t1 + cac(2)*sqrt(t1)*s2 + cac(3)*s1)

 end function acr3d




 real function smlt(tc, dqs, qsrho,psacw,psacr,c,rho, rhofac)
 real, intent(in):: tc,dqs,qsrho,psacw,psacr,c(5),rho, rhofac
     
 smlt = (c(1)*tc/rho-c(2)*dqs) * (c(3)*sqrt(qsrho)+ &
         c(4)*qsrho**0.65625*sqrt(rhofac)) + c(5)*tc*(psacw+psacr)

 end function smlt
 

 real function gmlt(tc, dqs,qgrho,pgacw,pgacr,c, rho)
 real, intent(in)::  tc,dqs,qgrho,pgacw,pgacr,c(5),rho
     
!     note:  pgacw and pgacr must be calc before gmlt is called
!
 gmlt = (c(1)*tc/rho-c(2)*dqs) * (c(3)*sqrt(qgrho)+ &
         c(4)*qgrho**0.6875/rho**0.25) + c(5)*tc*(pgacw+pgacr)
 end function gmlt


 subroutine qsmith_init
  integer, parameter:: length=2621 
  integer i

  if( .not. allocated(table) ) then
!                            generate es table (dt = 0.1 deg. c)
       allocate ( table( length) )
       allocate ( table2(length) )
       allocate ( table3(length) )
       allocate ( tablew(length) )
       allocate (   des (length) )
       allocate (   des2(length) )
       allocate (   des3(length) )
       allocate (   desw(length) )

       call qs_table (length )
       call qs_table2(length )
       call qs_table3(length )
       call qs_tablew(length )

       do i=1,length-1
           des(i) = max(0.,  table(i+1) -  table(i))
          des2(i) = max(0., table2(i+1) - table2(i))
          des3(i) = max(0., table3(i+1) - table3(i))
          desw(i) = max(0., tablew(i+1) - tablew(i))
       enddo
        des(length) =  des(length-1)
       des2(length) = des2(length-1)
       des3(length) = des3(length-1)
       desw(length) = desw(length-1)
  endif
 
 end subroutine qsmith_init
 

 real function qs1d(ta, pa, dqdt)
! 2-phase tabel
  real, intent(in):: ta, pa
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
      qs1d = eps*es/pa
        it = ap1 - 0.5
      dqdt = eps10*(des2(it) + (ap1-it)*(des2(it+1)-des2(it)))/pa

 end function qs1d


 real function ws1d(ta, pa, dqdt)
! Pure water phase
  real, intent(in):: ta, pa
  real, intent(out):: dqdt
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  real, parameter:: eps10 = 10.*eps
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
      ws1d = eps*es/pa
        it = ap1 - 0.5
      dqdt = eps10*(desw(it) + (ap1-it)*(desw(it+1)-desw(it)))/pa

 end function ws1d


 real function wqsat(ta, pa)
! Pure water phase
  real, intent(in):: ta, pa
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = tablew(it) + (ap1-it)*desw(it)
     wqsat = eps*es/pa

 end function wqsat

 real function iqsat(ta, pa)
  real, intent(in):: ta, pa
! local:
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table2(it) + (ap1-it)*des2(it)
     iqsat = eps*es/pa

 end function iqsat

 real function d_sat(ta)
! Computes the difference in saturation vapor *density* between water and ice
  real, intent(in):: ta
  real, parameter:: tmin=tice - 160.
  real es_w, es_i, ap1
  integer it

       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
! over Water:
       es_w = tablew(it) + (ap1-it)*desw(it)
! over Ice:
       es_i = table2(it) + (ap1-it)*des2(it)
      d_sat = dim(es_w, es_i)/(rvgas*ta)  ! Take positive difference

 end function d_sat


 real function esw_table(ta)
! pure water phase table
  real, intent(in):: ta
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      esw_table = tablew(it) + (ap1-it)*desw(it)
 end function esw_table


 real function es2_table(ta)
! two-phase table
  real, intent(in):: ta
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer it
       ap1 = 10.*dim(ta, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
      es2_table = table2(it) + (ap1-it)*des2(it)
 end function es2_table


 subroutine esw_table1d(ta, es, n)
  integer, intent(in):: n
! For waterphase only
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = tablew(it) + (ap1-it)*desw(it)
  enddo
 end subroutine esw_table1d



 subroutine es2_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
! For sea ice model
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table2(it) + (ap1-it)*des2(it)
  enddo
 end subroutine es2_table1d


 subroutine es3_table1d(ta, es, n)
  integer, intent(in):: n
! two-phase table with -2C as the transition point for ice-water phase
  real, intent(in)::  ta(n)
  real, intent(out):: es(n)
  real, parameter:: tmin=tice - 160.
  real  ap1
  integer i, it

  do i=1, n
       ap1 = 10.*dim(ta(i), tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
     es(i) = table3(it) + (ap1-it)*des3(it)
  enddo
 end subroutine es3_table1d



 subroutine qs_tablew(n)
! 2-phase table
      integer, intent(in):: n
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
        tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!  compute es over water
!  see smithsonian meteorological tables page 350.
        aa  = -7.90298*(tbasw/tem-1.)
        b   =  5.02808*alog10(tbasw/tem)
        c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
        d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
        e   = alog10(esbasw)
        tablew(i) = 0.1 * 10**(aa+b+c+d+e)
     enddo

 end subroutine qs_tablew


 subroutine qs_table2(n)
! 2-phase table
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
        if ( i<= 1600 ) then
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table2(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between 0c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table2(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1600;  i1 = 1601
      tem0 = 0.25*(table2(i0-1) + 2.*table(i0) + table2(i0+1))
      tem1 = 0.25*(table2(i1-1) + 2.*table(i1) + table2(i1+1))
      table2(i0) = tem0
      table2(i1) = tem1

 end subroutine qs_table2



 subroutine qs_table3(n)
! 2-phase table with "-2 C" as the transition point
  integer, intent(in):: n
  real:: delt=0.1
  real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e
  integer :: i0, i1
  real :: tem0, tem1
  integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
      tmin = tbasi - 160.

     do i=1,n
        tem = tmin+delt*real(i-1)
!       if ( i<= 1600 ) then
        if ( i<= 1580 ) then  ! to -2 C
!  compute es over ice between -160c and 0 c.
!  see smithsonian meteorological tables page 350.
              aa  = -9.09718 *(tbasi/tem-1.)
              b   = -3.56654 *alog10(tbasi/tem)
              c   =  0.876793*(1.-tem/tbasi)
              e   = alog10(esbasi)
             table3(i) = 0.1 * 10**(aa+b+c+e)
        else
!  compute es over water between -2c and 102c.
!  see smithsonian meteorological tables page 350.
             aa  = -7.90298*(tbasw/tem-1.)
             b   =  5.02808*alog10(tbasw/tem)
             c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
             d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
             e   = alog10(esbasw)
             table3(i) = 0.1 * 10**(aa+b+c+d+e)
        endif
     enddo

!----------
! smoother
!----------
      i0 = 1580
      tem0 = 0.25*(table3(i0-1) + 2.*table(i0) + table3(i0+1))
      i1 = 1581
      tem1 = 0.25*(table3(i1-1) + 2.*table(i1) + table3(i1+1))
      table3(i0) = tem0
      table3(i1) = tem1

 end subroutine qs_table3


 real function qs1d_blend(t, p, q)
! Note: this routine is based on "moist" mixing ratio
! Blended mixed phase table
  real, intent(in):: t, p, q
  real es, ap1
  real, parameter:: tmin=tice - 160.
  integer it

       ap1 = 10.*dim(t, tmin) + 1.
       ap1 = min(2621., ap1)
        it = ap1
        es = table(it) + (ap1-it)*des(it)
      qs1d_blend = eps*es*(1.+zvir*q)/p

 end function qs1d_blend

 subroutine qs_table(n)
      integer, intent(in):: n
      real esupc(200)
      real:: delt=0.1
      real esbasw, tbasw, esbasi, tbasi, tmin, tem, aa, b, c, d, e, esh20 
      real wice, wh2o
      integer i

! constants
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16

!  compute es over ice between -160c and 0 c.
      tmin = tbasi - 160.
!  see smithsonian meteorological tables page 350.
      do i=1,1600
         tem = tmin+delt*real(i-1)
         aa  = -9.09718 *(tbasi/tem-1.)
         b   = -3.56654 *alog10(tbasi/tem)
         c   =  0.876793*(1.-tem/tbasi)
         e   = alog10(esbasi)
         table(i)=10**(aa+b+c+e)
      enddo

!  compute es over water between -20c and 102c.
!  see smithsonian meteorological tables page 350.
      do  i=1,1221
          tem = 253.16+delt*real(i-1)
          aa  = -7.90298*(tbasw/tem-1.)
          b   =  5.02808*alog10(tbasw/tem)
          c   = -1.3816e-07*(10**((1.-tem/tbasw)*11.344)-1.)
          d   =  8.1328e-03*(10**((tbasw/tem-1.)*(-3.49149))-1.)
          e   = alog10(esbasw)
          esh20  = 10**(aa+b+c+d+e)
          if (i <= 200) then
              esupc(i) = esh20
          else
              table(i+1400) = esh20
          endif
      enddo

!  derive blended es over ice and supercooled water between -20c and 0c
      do i=1,200
         tem  = 253.16+delt*real(i-1)
         wice = 0.05*(273.16-tem)
         wh2o = 0.05*(tem-253.16)
         table(i+1400) = wice*table(i+1400)+wh2o*esupc(i)
      enddo

      do i=1,n
         table(i) = table(i)*0.1
      enddo

 end subroutine qs_table


 subroutine qsmith(im, km, ks, t, p, q, qs, dqdt)
! input t in deg k; p (pa) : moist pressure
  integer, intent(in):: im, km, ks
  real, intent(in),dimension(im,km):: t, p, q
  real, intent(out),dimension(im,km):: qs
  real, intent(out), optional:: dqdt(im,km)
! local:
  real, parameter:: eps10 = 10.*eps
  real es(im,km)
  real ap1
  real tmin
  integer i, k, it

  tmin = tice-160.

  if( .not. allocated(table) ) then
       call  qsmith_init
  endif
 
      do k=ks,km
         do i=1,im
            ap1 = 10.*dim(t(i,k), tmin) + 1.
            ap1 = min(2621., ap1)
            it = ap1
            es(i,k) = table(it) + (ap1-it)*des(it)
            qs(i,k) = eps*es(i,k)*(1.+zvir*q(i,k))/p(i,k)
         enddo
      enddo

      if ( present(dqdt) ) then
      do k=ks,km
           do i=1,im
              ap1 = 10.*dim(t(i,k), tmin) + 1.
              ap1 = min(2621., ap1) - 0.5
              it  = ap1
              dqdt(i,k) = eps10*(des(it)+(ap1-it)*(des(it+1)-des(it)))*(1.+zvir*q(i,k))/p(i,k)
           enddo
      enddo
      endif
 
 end subroutine qsmith


 subroutine neg_adj(ktop, kbot, p1, pt, dp, qv, ql, qr, qi, qs, qg)
! 1d version:
! this is designed for 6-class micro-physics schemes
 integer, intent(in):: ktop, kbot
 real, intent(in):: dp(ktop:kbot), p1(ktop:kbot)
 real, intent(inout), dimension(ktop:kbot)::    &
                                pt, qv, ql, qr, qi, qs, qg
! local:
 real lcpk(ktop:kbot), icpk(ktop:kbot)
 real dq, tmp1
 integer k

 do k=ktop,kbot
!      tmp1 = cp - rdgas*ptop/p1(k)
!   lcpk(k) = latv / tmp1
!   icpk(k) = lati / tmp1
    lcpk(k) = latv / cp
    icpk(k) = lati / cp
 enddo

 do k=ktop, kbot
!-----------
! ice-phase:
!-----------
! if ice<0 borrow from snow
          if( qi(k) < 0. ) then
              qs(k) = qs(k) + qi(k)
              qi(k) = 0.
          endif
! if snow<0 borrow from graupel
          if( qs(k) < 0. ) then
              qg(k) = qg(k) + qs(k)
              qs(k) = 0.
          endif
! if graupel < 0 then borrow from rain
          if ( qg(k) < 0. ) then
               qr(k) = qr(k) + qg(k)
               pt(k) = pt(k) - qg(k)*icpk(k)   ! heating
               qg(k) = 0.
          endif

! liquid phase:
! fix negative rain by borrowing from cloud water
          if ( qr(k) < 0. ) then
               ql(k) = ql(k) + qr(k)
               qr(k) = 0.
          endif
! fix negative cloud water with vapor
          if ( ql(k) < 0. ) then
               qv(k) = qv(k) + ql(k)
               pt(k) = pt(k) - ql(k)*lcpk(k)
               ql(k) = 0.
          endif
 enddo

!-----------------------------------
! fix water vapor; borrow from below
!-----------------------------------
 do k=ktop,kbot-1
    if( qv(k) < 0. ) then
        qv(k+1) = qv(k+1) + qv(k)*dp(k)/dp(k+1)
        qv(k  ) = 0.
    endif
 enddo
 
! bottom layer; borrow from above
 if( qv(kbot) < 0. .and. qv(kbot-1)>0.) then
             dq = min(-qv(kbot)*dp(kbot), qv(kbot-1)*dp(kbot-1))
     qv(kbot-1) = qv(kbot-1) - dq/dp(kbot-1) 
     qv(kbot  ) = qv(kbot  ) + dq/dp(kbot  ) 
 endif
! if qv is still < 0

 end subroutine neg_adj




 subroutine sg_conv(is, ie, js, je, isd, ied, jsd, jed,               &
                    isc, iec, jsc, jec,  km, nq, dt, tau,             &
                    delp, phalf, pm, zfull, zhalf, ta, qa, ua, va, w, &
                    u_dt, v_dt, t_dt, q_dt, mcond, nqv, nql, nqi, &
                    hydrostatic, phys_hydrostatic)
! Non-precipitating sub-grid scale convective adjustment-mixing
!-------------------------------------------
      logical, intent(in):: hydrostatic, phys_hydrostatic
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: mcond
      integer, intent(in):: isc, iec, jsc, jec
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau            ! Relaxation time scale
      integer, intent(in):: nqv, nql, nqi  ! vapor, liquid, ice
      real, intent(in):: dt             ! model time step
      real, intent(in):: phalf(is:ie,js:je,km+1) 
      real, intent(in):: pm(is:ie,js:je,km)
      real, intent(in):: zfull(is:ie,js:je,km)
      real, intent(in):: zhalf(is:ie,js:je,km+1)
      real, intent(in):: delp(isd:ied,jsd:jed,km)      ! Delta p at each model level
      real, intent(in)::   ta(isd:ied,jsd:jed,km)      ! Temperature
      real, intent(in)::   qa(isd:ied,jsd:jed,km,nq)   ! Specific humidity & tracers
      real, intent(in)::   ua(isd:ied,jsd:jed,km)
      real, intent(in)::   va(isd:ied,jsd:jed,km)
      real, intent(inout):: w(isd:ied,jsd:jed,km)
! Output:
! Updated fields:
      real, intent(out):: u_dt(isd:ied,jsd:jed,km)   ! updated u-wind field
      real, intent(out):: v_dt(isd:ied,jsd:jed,km)   !         v-wind
      real, intent(out):: t_dt(isc:iec,jsc:jec,km)   !         temperature
      real, intent(out):: q_dt(isc:iec,jsc:jec,km,nq) !
!---------------------------Local variables-----------------------------
      real, dimension(is:ie,km):: tvm, u0, v0, w0, t0, gz, hd, pkz
      real, dimension(is:ie,km+1):: pk, peln
      real q0(is:ie,km,nq) 
      real gzh(is:ie)
      real pbot, ri, pt1, pt2, lf, ratio
      real rdt, dh, dh0, dhs, dq, tv, h0, mc, mx,  fra, rk, rz, rcp
      real qs1, detn
      real clouds, rqi
      integer kcond
      integer i, j, k, n, m, iq, kk, ik
      real, parameter:: ustar2 = 1.E-8
      real, parameter:: dh_min = 1.E-4

      if ( nqv /= 1 .or. nql/=2 ) then
           call error_mesg ('sg_conv', 'Tracer indexing error', FATAL) 
      endif



    rz = rvgas - rdgas          ! rz = zvir * rdgas
    rk = cp_air/rdgas + 1.
   rcp = 1./cp_air

    m = 4
    rdt = 1. / dt
    fra = dt/real(tau)

!------------
! Compute gz: center 
!------------
  do 1000 j=js,je       ! this main loop can be OpneMPed in j

    do k=mcond,km+1
       do i=is,ie
          peln(i,k) = log(phalf(i,j,k))
!           pk(i,k) = phalf(i,j,k)**kappa
            pk(i,k) = exp(kappa*peln(i,k))
       enddo
    enddo

    do k=mcond,km
       do i=is,ie
          u0(i,k) = ua(i,j,k)
          v0(i,k) = va(i,j,k)
          t0(i,k) = ta(i,j,k)
         pkz(i,k) = (pk(i,k+1)-pk(i,k))/(kappa*(peln(i,k+1)-peln(i,k)))
       enddo
    enddo

    if ( .not.hydrostatic ) then
       do k=mcond,km
          do i=is,ie
             w0(i,k) = w(i,j,k)
          enddo
       enddo
    endif

    do iq=1,nq
       do k=mcond,km
          do i=is,ie
             q0(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
    enddo


!-----------------
! K-H instability:
!-----------------
   kcond = mcond

   do n=1,m
      ratio = real(n)/real(m)

    if( phys_hydrostatic ) then
       do i=is,ie
          gzh(i) = 0.
       enddo
       do k=km, mcond,-1
          do i=is,ie
           tvm(i,k) = t0(i,k)*(1.+zvir*q0(i,k,nqv))
                tv  = rdgas*tvm(i,k)
            gz(i,k) = gzh(i) + tv*(1.-phalf(i,j,k)/pm(i,j,k))
            hd(i,k) = cp_air*tvm(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2)
             gzh(i) = gzh(i) + tv*(peln(i,k+1)-peln(i,k))
          enddo
       enddo
       do i=is,ie
          gzh(i) = 0.
       enddo
    else
       do k=mcond,km
          do i=is,ie
             gz(i,k) = grav*zfull(i,j,k)
             hd(i,k) = cp_air*t0(i,k)+gz(i,k)+0.5*(u0(i,k)**2+v0(i,k)**2+w0(i,k)**2)
          enddo
       enddo
    endif

      do k=km,kcond+1,-1
         do i=is,ie
! Richardson number at interface: g*delz * (del_theta/theta) / (del_u**2 + del_v**2)
            pt1 = t0(i,k-1)/pkz(i,k-1)
            pt2 = t0(i,k  )/pkz(i,k  )
             ri = (gz(i,k-1)-gz(i,k))*(pt1-pt2)/( 0.5*(pt1+pt2)*        &
                 ((u0(i,k-1)-u0(i,k))**2+(v0(i,k-1)-v0(i,k))**2+ustar2) )
! Dry convective mixing for K-H instability & CAT (Clear Air Turbulence):
! Compute equivalent mass flux: mc
#ifndef USE_RIP1  
            if ( ri < 0.25 ) then
                 mc = ratio * (1.-max(0.0, 4.*ri)) ** 2
#else
            if ( ri < 1. ) then
                 mc = ratio * (1.-max(0.0, ri)) ** 2
#endif 
                 mc = mc*delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                 do iq=1,nq
                              h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                    q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                    q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                 enddo
! u:
                        h0 = mc*(u0(i,k)-u0(i,k-1))
                 u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                 u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                        h0 = mc*(v0(i,k)-v0(i,k-1))
                 v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                 v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
! h:
                          h0 = mc*(hd(i,k)-hd(i,k-1))
                   hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                   hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
                if ( .not.hydrostatic ) then
                           h0 = mc*(w0(i,k)-w0(i,k-1))
                    w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                    w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
                endif
            endif
         enddo
!--------------
! Retrive Temp:
!--------------
      if ( phys_hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - phalf(i,j,kk)/pm(i,j,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1)-peln(i,kk))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,nqv) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-phalf(i,j,kk)/pm(i,j,kk))*(rdgas+rz*q0(i,kk,nqv)))
         enddo
      else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               t0(i,kk) = rcp*(hd(i,kk)-gz(i,kk)-0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2))
            enddo
         enddo
      endif
      enddo
   enddo       ! n-loop


!-------------------------
! Moist adjustment/mixing:
!-------------------------
 m = 3

 if( km>k_moist+1 ) then
   do n=1,m

    ratio = real(n)/real(m)

    if ( phys_hydrostatic ) then
       do i=is,ie
          gzh(i) = 0.
       enddo
    endif

    do k=km,max(kcond,k_moist)+1,-1
       do i=is,ie
          if ( phalf(i,j,k) > p_crt ) then
!--------------------------------------------------------------------
!           qs1 = qs1d_blend(t0(i,k-1), pm(i,j,k-1), q0(i,k-1,nqv))
!            lf = hlv + hlf*min(1.0, max(0.0, (tice-t0(i,k-1))/30.))
!--------------------------------------------------------------------
            clouds = q0(i,k-1,nql) + q0(i,k-1,nqi)
            if( clouds > 1.e-5 ) then
                rqi = q0(i,k-1,nqi) / clouds
            else
                rqi = max(0., min(1., (tice-t0(i,k-1))/30.))
            end if
            qs1 = rqi*es2_table(t0(i,k-1)) + (1.-rqi)*esw_table(t0(i,k-1))
            qs1 = eps*qs1*(1.+zvir*q0(i,k-1,nqv))/pm(i,j,k-1)
             lf = hlv + rqi*hlf

              dh0 = hd(i,k) - hd(i,k-1)
              dhs = dh0 + lf*(q0(i,k,nqv)-qs1        )
              dh  = dh0 + lf*(q0(i,k,nqv)-q0(i,k-1,nqv))
!             if ( dhs>0.0 .and. dh>dh_min ) then
              if ( dhs>0.0 .and. dh>dh_min .and. q0(i,k,nqv)>q0(i,k-1,nqv) ) then
                   mc = ratio*min(1.0, 0.5*dhs/dh)*    &
                        delp(i,j,k-1)*delp(i,j,k)/(delp(i,j,k-1)+delp(i,j,k))
                          h0 = mc*dh0
                   hd(i,k-1) = hd(i,k-1) + h0/delp(i,j,k-1)
                   hd(i,k  ) = hd(i,k  ) - h0/delp(i,j,k  )
! Perform local mixing of all advected tracers:
#ifdef DET_CON
                 if ( zhalf(i,j,k) > (1.E3+zhalf(i,j,km+1)) ) then
                      detn = min(1., zhalf(i,j,k)/7.e3)
! specific humidity:
                              h0 = mc*(q0(i,k,nqv)-q0(i,k-1,nqv))
                              dq = h0/delp(i,j,k-1)
                   q0(i,k-1,nqv) = q0(i,k-1,nqv) + dq*(1.-detn)
                   q0(i,k  ,nqv) = q0(i,k  ,nqv) - h0/delp(i,j,k  )
                   do iq=2,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
!--------------
! Condensation:
!--------------
                   dq = dq * detn
                   q0(i,k-1,nql) = q0(i,k-1,nql) + dq*(1.-rqi)
                   q0(i,k-1,nqi) = q0(i,k-1,nqi) + dq*rqi
                   hd(i,k-1) = hd(i,k-1) + dq*lf

                 else
                   do iq=1,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
                 endif
#else
                   do iq=1,nq
                                h0 = mc*(q0(i,k,iq)-q0(i,k-1,iq))
                      q0(i,k-1,iq) = q0(i,k-1,iq) + h0/delp(i,j,k-1)
                      q0(i,k  ,iq) = q0(i,k  ,iq) - h0/delp(i,j,k  )
                   enddo
#endif
! u:
                          h0 = mc*(u0(i,k)-u0(i,k-1))
                   u0(i,k-1) = u0(i,k-1) + h0/delp(i,j,k-1)
                   u0(i,k  ) = u0(i,k  ) - h0/delp(i,j,k  )
! v:
                          h0 = mc*(v0(i,k)-v0(i,k-1))
                   v0(i,k-1) = v0(i,k-1) + h0/delp(i,j,k-1)
                   v0(i,k  ) = v0(i,k  ) - h0/delp(i,j,k  )
! *** Non-hydrostatic:
                  if ( .not.hydrostatic ) then
                          h0 = mc*(w0(i,k)-w0(i,k-1))
                   w0(i,k-1) = w0(i,k-1) + h0/delp(i,j,k-1)
                   w0(i,k  ) = w0(i,k  ) - h0/delp(i,j,k  )
                  endif
! ***
              endif  ! dh check
            endif    ! p_crt check
         enddo
!--------------
! Retrive Temp:
!--------------
       if ( phys_hydrostatic ) then
         kk = k
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ( rk - phalf(i,j,kk)/pm(i,j,kk) )
              gzh(i) = gzh(i) + t0(i,kk)*(peln(i,kk+1)-peln(i,kk))
            t0(i,kk) = t0(i,kk) / ( rdgas + rz*q0(i,kk,nqv) )
         enddo
         kk = k-1
         do i=is,ie
            t0(i,kk) = (hd(i,kk)-gzh(i)-0.5*(u0(i,kk)**2+v0(i,kk)**2))  &
                     / ((rk-phalf(i,j,kk)/pm(i,j,kk))*(rdgas+rz*q0(i,kk,nqv)))
         enddo
       else
! Non-hydrostatic under constant volume heating/cooling
         do kk=k-1,k
            do i=is,ie
               t0(i,kk) = rcp*(hd(i,kk)-gz(i,kk)-0.5*(u0(i,kk)**2+v0(i,kk)**2+w0(i,kk)**2))
            enddo
         enddo
       endif
      enddo
   enddo       ! n-loop
 endif      ! k_moist check

   if ( fra < 1. ) then
      do k=mcond,km
         do i=is,ie
            t0(i,k) = ta(i,j,k) + (t0(i,k) - ta(i,j,k))*fra
            u0(i,k) = ua(i,j,k) + (u0(i,k) - ua(i,j,k))*fra
            v0(i,k) = va(i,j,k) + (v0(i,k) - va(i,j,k))*fra
         enddo
      enddo

      if ( .not.hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w0(i,k) = w(i,j,k) + (w0(i,k) - w(i,j,k))*fra
         enddo
      enddo
      endif

      do iq=1,nq
         do k=mcond,km
            do i=is,ie
               q0(i,k,iq) = qa(i,j,k,iq) + (q0(i,k,iq) - qa(i,j,k,iq))*fra
            enddo
         enddo
      enddo
   endif

!--------------------
! Update fields:
!--------------------
   do k=1,mcond-1
      do i=is,ie
         u_dt(i,j,k) = ua(i,j,k)
         v_dt(i,j,k) = va(i,j,k)
         t_dt(i,j,k) = ta(i,j,k)
      enddo
   enddo
   do k=mcond,km
      do i=is,ie
         u_dt(i,j,k) = u0(i,k)
         v_dt(i,j,k) = v0(i,k)
         t_dt(i,j,k) = t0(i,k)
      enddo
   enddo

   if ( .not.hydrostatic ) then
      do k=mcond,km
         do i=is,ie
            w(i,j,k) = w0(i,k)
         enddo
      enddo
   endif

   do iq=1,nq
      do k=1,mcond-1
         do i=is,ie
            q_dt(i,j,k,iq) = qa(i,j,k,iq)
         enddo
      enddo
      do k=mcond,km
         do i=is,ie
            q_dt(i,j,k,iq) = q0(i,k,iq)
         enddo
      enddo
   enddo

1000 continue


 end subroutine sg_conv


 real function g_sum(p, ifirst, ilast, jfirst, jlast, ngc, area, mode)
      use mpp_mod,           only: mpp_sum
      real, save :: global_area

! Fast version of globalsum
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast, ngc
      integer, intent(IN) :: mode  ! if ==1 divided by area
      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
      integer :: i,j
      real gsum

!-------------------------
! Quick local sum algorithm
!-------------------------
      if ( .not. g_sum_initialized ) then
         allocate (l_area(ifirst:ilast,jfirst:jlast))
         global_area = 0.
         do j=jfirst,jlast
           do i=ifirst,ilast
             global_area = global_area + area(i,j)
             l_area(i,j) = area(i,j)
           enddo
         enddo
         call mpp_sum(global_area)
!        if ( mpp_pe().eq.mpp_root_pe() ) write(*,*) 'Global Area=',global_area
         g_sum_initialized = .true.
      end if

      gsum = 0.
      do j=jfirst,jlast
        do i=ifirst,ilast
          gsum = gsum + p(i,j)*l_area(i,j)
        enddo
      enddo
      call mpp_sum(gsum)

      if ( mode==1 ) then
        g_sum = gsum / global_area
      else
        g_sum = gsum
      endif

 end function g_sum

end module lin_cld_microphys_mod

!

subroutine bio_v3(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, Grid, Time, dtts, Thickness, Dens, swflx, sw_frac_zt)


! Based on the pt_npzd.aos version  - which is the optimized version

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc, iec
integer, intent(in)                                             :: jsc, jec
integer, intent(in)                                             :: isd, ied
integer, intent(in)                                             :: jsd, jed
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, intent(in)                                                :: dtts
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in), dimension(isd:ied,jsd:jed)                    :: swflx        ! short wave radiation flux (W/m^2)
real, intent(in), dimension(isd:,jsd:,:)                    :: sw_frac_zt        ! short wave fraction on T grid (none)

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

integer :: i
integer :: j
integer :: k
integer :: n

logical :: used
integer :: index_temp, index_salt

! BGC parameters now read from bgc_param.nc at run time. mac, aug11.  
! include "rjm_param_201005.inc"

      integer :: ts_npzd ! number of time steps within NPZD model
      integer :: kmeuph = 20 ! deepest level of  euphotic zone
      integer :: tn, trn

      real :: pi    = 3.14159265358979  !yes, this is pi
      real :: biotr(isc:iec,grid%nk,ntr_bmax), bioma(isc:iec,grid%nk,ntr_bmax), pprod(isc:iec,jsc:jec,grid%nk)
      real :: bion,biop,bioz,biod,bioo,biocaco3,bioi
      real :: u_npz,g_npz
      real :: fx1,fx2,fx3,fx4,fu1,fu2,radbio, sw_zt, sw_zt1, sw_zt2
      real :: vpbio(isc:iec,grid%nk)
      real :: avej(isc:iec,grid%nk)
!chd auxiliary variables to prevent unnecessary computation
      real :: fbc
      real :: f11,f21,f22,f23,f31,f32,f41,f51 
      real :: epsi      = 1e-15
      real :: rdtts      !1/dtts
      real :: dtsb
      real :: adv_fb(isc:iec,1:grid%nk+1)
      real, dimension(isd:ied,jsd:jed) :: mld
      real :: caco3_bgc_change, no3_bgc_change

!
! =====================================================================
!     begin executable code
! =====================================================================

! read biotic parameters

  call time_interp_external(alphabio_id, time%model_time, alphabio)
  call time_interp_external(parbio_id, time%model_time, parbio)
  call time_interp_external(kwbio_id, time%model_time, kwbio)
  call time_interp_external(kcbio_id, time%model_time, kcbio)
  call time_interp_external(abio_id, time%model_time, abio)
  call time_interp_external(bbio_id, time%model_time, bbio)
  call time_interp_external(cbio_id, time%model_time, cbio)
  call time_interp_external(k1bio_id, time%model_time, k1bio)
  call time_interp_external(muepbio_id, time%model_time, muepbio)
  call time_interp_external(muepsbio_id, time%model_time, muepsbio)
  call time_interp_external(gam1bio_id, time%model_time, gam1bio)
  call time_interp_external(gbio_id, time%model_time, gbio)
  call time_interp_external(epsbio_id, time%model_time, epsbio)
  call time_interp_external(muezbio_id, time%model_time, muezbio)
  call time_interp_external(gam2bio_id, time%model_time, gam2bio)
  call time_interp_external(muedbio_id, time%model_time, muedbio)
  call time_interp_external(wdetbio_id, time%model_time, wdetbio)
  call time_interp_external(muecaco3_id, time%model_time, muecaco3)
  call time_interp_external(wcaco3_id, time%model_time, wcaco3)
  call time_interp_external(tscav_fe_id, time%model_time, tscav_fe)
  call time_interp_external(fe_bkgnd_id, time%model_time, fe_bkgnd)
  call time_interp_external(f_inorg_id, time%model_time, f_inorg)





! set the maximum index for euphotic depth
      do k=1,grid%nk
         if (grid%zw(k).le. 400) kmeuph=k
      enddo
!        print*,'rjm euphotic ', kmeuph

!
!
!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------
!

  
  index_temp = fm_get_index('/ocean_mod/prog_tracers/temp')
  index_salt = fm_get_index('/ocean_mod/prog_tracers/salt')


  ts_npzd= max(1,nint(dtts / 900.))
  rdtts  = 1/dtts

!  write (stdout(),*) ' AO-NPZD model will do ',ts_npzd,' time steps'
!  write (stdout(),*) ' time step in NPZD model will be ',      &
!                dtts/ts_npzd,'sec.'
       
!chd   time step within NPZD model
!chd
      dtsb=dtts/float(ts_npzd)

! Calculate the mixed layer depth.  mac, aug11.  
       call calc_mixed_layer_depth(Thickness, &
        T_prog(index_salt)%field(isd:ied,jsd:jed,:,Time%tau), &
        T_prog(index_temp)%field(isd:ied,jsd:jed,:,Time%tau), &
        Dens%rho(isd:ied,jsd:jed,:,Time%tau),                 &
        Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
        mld)

!
!       Loop over multiple instances
!



do n = 1, instances  !{

  pprod_gross(:,:,:) = 0.0
  zprod_gross(:,:,:) = 0.0
  light_limit(:,:) = 0.0
  radbio3d(:,:,:) = 0.0
  npp3d(:,:,:) = 0.0
  adic_intmld(:,:) = 0.0
  dic_intmld(:,:) = 0.0
  o2_intmld(:,:) = 0.0
  no3_intmld(:,:) = 0.0
  fe_intmld(:,:) = 0.0
  phy_intmld(:,:) = 0.0
  det_intmld(:,:) = 0.0
  pprod_gross_intmld(:,:) = 0.0
  npp_intmld(:,:) = 0.0
  radbio_intmld(:,:) = 0.0
  wdet100avg(:,:) = 0.0


! possible tracers in model
  ind_dic = biotic(n)%ind_bgc(id_dic)
  ind_adic = biotic(n)%ind_bgc(id_adic)
  ind_alk = biotic(n)%ind_bgc(id_alk)
  ind_po4 = biotic(n)%ind_bgc(id_po4)
  ind_o2  = biotic(n)%ind_bgc(id_o2)

  ind_no3  = biotic(n)%ind_bgc(id_no3)

  ind_zoo  = biotic(n)%ind_bgc(id_zoo)
  ind_phy  = biotic(n)%ind_bgc(id_phy)
  ind_det  = biotic(n)%ind_bgc(id_det)

  ind_caco3= biotic(n)%ind_bgc(id_caco3)
  ind_fe   = biotic(n)%ind_bgc(id_fe)

!chd calculate biotic source terms using Euler
!chd forward timestepping

  do j = jsc, jec  !{ 


!chd    transfer tracers at tau-1 to temporary arrays
        do k=1,grid%nk
          do i=isc,iec
            biotr(i,k,id_no3) = max(0.0,t_prog(ind_no3)%field(i,j,k,Time%taum1)) 
            biotr(i,k,id_phy) = max(0.0,t_prog(ind_phy)%field(i,j,k,Time%taum1))
            biotr(i,k,id_zoo) = max(0.0,t_prog(ind_zoo)%field(i,j,k,Time%taum1))
            biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
            biotr(i,k,id_o2) = max(0.0,t_prog(ind_o2)%field(i,j,k,Time%taum1)) 
           if (id_caco3.ne.0) &
              biotr(i,k,id_caco3) = max(0.0,t_prog(ind_caco3)%field(i,j,k,Time%taum1)) 
           if (id_fe.ne.0) &
              biotr(i,k,id_fe) = max(0.0,t_prog(ind_fe)%field(i,j,k,Time%taum1)) 
          enddo
        enddo

!chd   create mask bioma, which is 0 for biotr < epsi and 1
!chd   otherwise
       do trn=1,ntr_bmax   
        do k=1,grid%nk
          do i=isc,iec
           bioma(i,k,trn) = 0.
           if (biotr(i,k,trn) .gt. epsi) then
            bioma(i,k,trn) = 1.
           endif 
          enddo
        enddo
       enddo


       avej(:,:)= 0.0
       vpbio(:,:)= 0.0

       do k=1,kmeuph  
        do i=isc,iec

! Shortwave intensity averaged over layer
	   sw_zt2 = (sw_frac_zt(i,j,k)-sw_frac_zt(i,j,k+1))/(grid%zt(k)-grid%zt(k+1) ) * (grid%zw(k)-grid%zt(k)) + sw_frac_zt(i,j,k)   
	  if(k.eq.1) then
            sw_zt1=1.
	   else
           sw_zt1 = (sw_frac_zt(i,j,k-1)-sw_frac_zt(i,j,k))/(grid%zt(k-1)-grid%zt(k) ) * (grid%zw(k-1)-grid%zt(k-1)) + sw_frac_zt(i,j,k-1)  
	   endif
	   sw_zt = (sw_zt1 + sw_zt2)*.5 * swflx(i,j)
!

! add ice feedback on growth ......
! Use sw_thru_ice = .true. to indicate ice term already included in swflx when running with ice model
!    otherwise shortwave flux from the ocean model is corrected here.  mac, sep11.  
          if (sw_thru_ice) then
            radbio = max(0.0,parbio(i,j)*sw_zt)
          else
            radbio = max(0.0,parbio(i,j)*(1.0-fice_t(i,j))*sw_zt)
          endif


! Assume now that shortwave flux from ocean model has been modified for ice.  mac, aug12
!          radbio = max(0.0,parbio(i,j)*sw_zt)
! well, now doing an experiment without ice, so assumption no longer true.  mac, oct12. 

          radbio3d(i,j,k) = radbio

          vpbio(i,k) = abio(i,j)*bbio(i,j)**                                  &
                  (cbio(i,j)*t_prog(index_temp)%field(i,j,k,time%tau))

!          vpbio(i,k)= 0.6* (1.006)**T_prog(index_temp)%field(i,j,k,Time%tau) 

!          avej(i,k) = vpbio(i,k)*alphabio*radbio/(vpbio(i,k)+(alphabio*radbio)**2)
!    Growth term from Evans and Parslow 1985, mac Apr10
!         avej(i,k) = vpbio(i,k)*alphabio*radbio/(vpbio(i,k)**2+(alphabio*radbio)**2)**0.5

!    Growth term from Brian Griffiths
        avej(i,k) = vpbio(i,k)*(1.0-exp(-1.0*(alphabio(i,j)*radbio)/vpbio(i,k)))
        
        if (Grid%zw(k) .le. mld(i,j)) &
          light_limit(i,j)=light_limit(i,j) +   &
            (1.0-exp(-1.0*(alphabio(i,j)*radbio)/vpbio(i,k))) * &
            thickness%dzt(i,j,k)

!	avej(i,k)=0.   ! zero phyto growth
        enddo
       enddo


!chd  begin time stepping over dtts in NPZD model with Euler forward
!chd     This is the NPZD model:
!chd      (P: phytoplankton, Z: Zooplankton, N: Nitrate and D: Detritus)
!chd
!chd      dP/dt = u(N,Temp.,Light) P - p_P P - g(P) P Z
!chd
!chd      dZ/dt = a g(P) P Z - d Z - p_Z Z^2
!chd
!chd      dN/dt = r D + d Z - u(N,Temp.,Light) P  [ + r_d DOC ]
!chd
!chd      dD/dt = (1-s)[ (1-a) g(P) P Z + p_P P + p_Z Z^2] -r D + w_D dD/dz

       do tn = 1,ts_npzd  !{
        do k=1,grid%nk  !{
         do i=isc,iec  !{
          bion = max(0.0,biotr(i,k,id_no3))
          biop = max(0.0,biotr(i,k,id_phy))
          bioz = max(0.0,biotr(i,k,id_zoo))
          biod = max(0.0,biotr(i,k,id_det))
          bioo = max(0.0,biotr(i,k,id_o2))
         if (id_caco3.ne.0) biocaco3 = max(0.0,biotr(i,k,id_caco3))
         if (id_fe.ne.0) bioi = max(0.0,biotr(i,k,id_fe))

!chd   -- phytoplankton equation
!chdc
!chdc      use Liebigs Law of the Minimum (Liebig, 1845) for growth rate
!chdc      (minimum of light-limited and nutrient limited growth rates;
!chdc      although chlorophyll is not explicitly considered, this will
!chdc      later allow for a diagnostic determination of a Chl:N ratio
!chdc      depending on light- or nutrient-limited growth.
!chdc      --> Hurtt and Armstrong, 1996)
!chdc      saturation growth rate (infinite light, infinite nutrients)
!chd          vpbio = abio*bbio**(cbio*t_prog(index_temp)%field(i,j,k,tau))
!chd      growth rate
          u_npz = min(avej(i,k),vpbio(i,k)*bion/(k1bio(i,j)+bion))
! rjm - iron limitation
         if(id_fe.ne.0) u_npz = min(u_npz,vpbio(i,k)*bioi/(0.1+bioi))

!chd          if (k .eq. 1) then
!chd           write(stdout(),*) 'ave'
!chd           write(stdout(),*) avej(i,k)
!chd          endif
!chd      grazing function
          g_npz = gbio(i,j)*epsbio(i,j)*biop*biop/(gbio(i,j)+epsbio(i,j)*biop*biop)
!chd      temperature dependance of growth rates
          fbc = bbio(i,j)**(cbio(i,j)*t_prog(index_temp)%field(i,j,k,time%tau))

        f11 = u_npz*biop         * bioma(i,k,id_no3)
        f21 = g_npz*bioz         * bioma(i,k,id_phy)
        f22 = muepbio(i,j)*fbc*biop   * bioma(i,k,id_phy)
        f23 = muepsbio(i,j)*biop*biop * bioma(i,k,id_phy)
        f31 = gam2bio(i,j)*fbc*bioz   * bioma(i,k,id_zoo)
        f32 = muezbio(i,j)*bioz*bioz  * bioma(i,k,id_zoo)
        f41 = muedbio(i,j)*fbc*biod   * bioma(i,k,id_det)
! if (grid%zw(k) .ge. 180) f41 = f41*.2   ! reduce decay below 300m
! change the ratio of remineralisation of det in upper to lower ocean from 5 to 2,
!  but keep the remin rate the same at depth, so change the value in rjm_param_2010.... as well
!  mac, jun10.
        if (grid%zw(k) .ge. 180) f41 = f41*.5   ! reduce decay below 300m

        if (id_caco3.ne.0) then 
          f51 = muecaco3(i,j)*biocaco3 *bioma(i,k,id_caco3)
        endif

!chd   -- nutrient equation
          biotr(i,k,id_no3)  = biotr(i,k,id_no3)  + dtsb * (                     &
               f41 + f31 + f22 - f11) 

!chd   -- phyto plankton equation
          biotr(i,k,id_phy)  = biotr(i,k,id_phy)  + dtsb * (                     &
               f11 - f21 - f22 - f23)

! Estimate primary productivity from phytoplankton growth
          pprod_gross(i,j,k) = pprod_gross(i,j,k) + dtsb*f11
! Net primary productivity (gross PP minus linear mortality)
	  npp3d(i,j,k) = npp3d(i,j,k) + dtsb*(f11 - f21)

!chd   -- zooplankton equation
          biotr(i,k,id_zoo)  = biotr(i,k,id_zoo)  + dtsb * (                     &
               gam1bio(i,j)*f21 - f31 - f32)

! Estimate secondary productivity from zooplankton growth
          zprod_gross(i,j,k) = zprod_gross(i,j,k) + dtsb*f21

!chd   -- detritus equation
          biotr(i,k,id_det)  = biotr(i,k,id_det)  + dtsb * (                     &
               (1-gam1bio(i,j))*f21 + f23 + f32 - f41)    

!chd   -- oxygen equation
          biotr(i,k,id_o2)  = biotr(i,k,id_o2)  - bioma(i,k,id_o2) * 172 / 16 * dtsb * ( &
               f41 + f31 + f22 - f11)

! rjm -- extra equation for caco3 - alkalinity
        if (id_caco3.ne.0) &
         biotr(i,k,id_caco3) = biotr(i,k,id_caco3)  + dtsb * ( (   &
               (1-gam1bio(i,j))*f21 + f23 + f32)*f_inorg(i,j)*106/16 - f51)  ! 8% of POC 106/16*.08   
 
! rjm -- extra equation for iron
!        if (id_fe.ne.0) &
!          biotr(i,k,id_fe)  = biotr(i,k,id_fe)  +dtsb * 200./16.* (           &
!               f41 + f31 + f22 - f11)

! mac -- extra equation for iron, molar Fe:N = 1.98e-5:1.0 (Christian et al. 2002), and Fe units are micro mole/m3 cf milli mole/m3 for others.
        if (id_fe.ne.0) &
          biotr(i,k,id_fe)  = biotr(i,k,id_fe)  +dtsb * 2.0e-2* (           &
               f41 + f31 + f22 - f11)


         enddo  !} i
        enddo  !} k
       enddo  !} tn 


!chd  add biotically induced tendency to biotracers,
        do k = 1, grid%nk  !{
         do i = isc, iec  !{

        no3_bgc_change = grid%tmask(i,j,k)* &
            (biotr(i,k,id_no3) - max(0.0,t_prog(ind_no3)%field(i,j,k,time%taum1)))
        caco3_bgc_change = grid%tmask(i,j,k)* &
           (biotr(i,k,id_caco3) - max(0.0,t_prog(ind_caco3)%field(i,j,k,time%taum1)))

        t_prog(ind_no3)%field(i,j,k,time%taum1) = biotr(i,k,id_no3)
        t_prog(ind_phy)%field(i,j,k,time%taum1) = biotr(i,k,id_phy)
        t_prog(ind_zoo)%field(i,j,k,time%taum1) = biotr(i,k,id_zoo)
        t_prog(ind_det)%field(i,j,k,time%taum1) = biotr(i,k,id_det)
	t_prog(ind_o2)%field(i,j,k,time%taum1) = biotr(i,k,id_o2)
        if ( id_caco3 .ne. 0 ) &
          t_prog(ind_caco3)%field(i,j,k,time%taum1) = biotr(i,k,id_caco3)
        if ( id_fe .ne. 0 ) &
          t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe) - &
           dtts * tscav_fe(i,j) * max(0.0,(biotr(i,k,id_fe) - fe_bkgnd(i,j)) ) 

        if (id_dic.ne.0) &
          t_prog(ind_dic)%field(i,j,k,time%taum1) = &
           t_prog(ind_dic)%field(i,j,k,time%taum1) + &
           106./16. * no3_bgc_change - caco3_bgc_change

        if (id_adic.ne.0) &
          t_prog(ind_adic)%field(i,j,k,time%taum1) = &
           t_prog(ind_adic)%field(i,j,k,time%taum1) + &
           106./16. * no3_bgc_change - caco3_bgc_change 

        if (id_alk.ne.0) &
          t_prog(ind_alk)%field(i,j,k,time%taum1) =  &
           t_prog(ind_alk)%field(i,j,k,time%taum1) +  &  
           ( - 2.0 * caco3_bgc_change - no3_bgc_change)

        pprod_gross(i,j,k)=rdtts*pprod_gross(i,j,k)*grid%tmask(i,j,k)
        zprod_gross(i,j,k)=rdtts*zprod_gross(i,j,k)*grid%tmask(i,j,k)
        npp3d(i,j,k)=rdtts*npp3d(i,j,k)*grid%tmask(i,j,k)

          if (Grid%zw(k) .le. mld(i,j)) then
            adic_intmld(i,j) = adic_intmld(i,j) + biotr(i,k,id_adic) * thickness%dzt(i,j,k)
            dic_intmld(i,j) = dic_intmld(i,j) + biotr(i,k,id_dic) * thickness%dzt(i,j,k)
            o2_intmld(i,j) = o2_intmld(i,j) + biotr(i,k,id_o2) * thickness%dzt(i,j,k)
            no3_intmld(i,j) = no3_intmld(i,j) + biotr(i,k,id_no3) * thickness%dzt(i,j,k)
            fe_intmld(i,j) = fe_intmld(i,j) + biotr(i,k,id_fe) * thickness%dzt(i,j,k)
            phy_intmld(i,j) = phy_intmld(i,j) + biotr(i,k,id_phy) * thickness%dzt(i,j,k)
            det_intmld(i,j) = det_intmld(i,j) + biotr(i,k,id_det) * thickness%dzt(i,j,k)
            pprod_gross_intmld(i,j) = pprod_gross_intmld(i,j) + pprod_gross(i,j,k) * thickness%dzt(i,j,k)
            npp_intmld(i,j) = npp_intmld(i,j) + npp3d(i,j,k) * thickness%dzt(i,j,k)
            radbio_intmld(i,j) = radbio_intmld(i,j) + radbio3d(i,j,k) * thickness%dzt(i,j,k)
          endif

         enddo  !} i
        enddo  !} k

        wdet100avg(:,:) = wdetbio(:,:)*t_prog(ind_det)%field(:,:,minloc(grid%zt(:)-100,dim=1),time%taum1)


!RASF upstream sinking of detritus
! change field(..,k,..) to field(..,k-1,..), mac apr10
     do k=2,grid%nk+1
      do i=isc,iec
!      adv_fb(i,k)=wdetbio(i,j)*t_prog(ind_det)%field(i,j,k-1,time%taum1) 
       adv_fb(i,k)=wdetbio(i,j)*biotr(i,k-1,id_det)
      enddo
     enddo


!RASF no flux boundary conditions
     do i=isc,iec
      adv_fb(i,1)               = 0.0
     enddo 

! Deposit tracer to sediment as tracer sinks through base of column.  mac, nov12
     do i = isc, iec 
       k = grid%kmt(i,j)
       if (k .gt. 0) then 

         biotic(n)%det_sed_depst(i,j) = adv_fb(i,k+1)

       endif ! k .gt. 0
     enddo  ! i


!!mac remineralise tracer sinking through the base of the column.
!     do i = isc, iec 
!       k = grid%kmt(i,j)
!       if (k .gt. 0) then 
!
!        t_prog(ind_no3)%th_tendency(i,j,k) = t_prog(ind_no3)%th_tendency(i,j,k) + &
!         grid%tmask(i,j,k)  * Thickness%rho_dzt(i,j,k,time%tau)  &
!         * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!       
!       t_prog(ind_o2)%th_tendency(i,j,k) = t_prog(ind_o2)%th_tendency(i,j,k) - &
!         172.0 / 16.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!         * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!      
!       if (id_fe .ne. 0) then
!        t_prog(ind_fe)%th_tendency(i,j,k) = t_prog(ind_fe)%th_tendency(i,j,k) + &
!         2.0e-2 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!         * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!       endif
!      
!       if (id_dic .ne. 0) then
!        t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
!         106.0 / 16.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!         * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!       endif
!      
!       if (id_adic .ne. 0) then
!        t_prog(ind_adic)%th_tendency(i,j,k) =t_prog(ind_adic)%th_tendency(i,j,k) + &
!         106.0 / 16.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!         * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!       endif
!      
!       if (id_alk .ne. 0) then
!        t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) - &
!         grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!         *adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!       endif
!
!       endif ! k .gt. 0
!
!     enddo  ! i


     do k =1,grid%nk  !{
      do i =isc,iec  !{
        t_prog(ind_det)%field(i,j,k,time%taum1) = t_prog(ind_det)%field(i,j,k,time%taum1) + &
            grid%tmask(i,j,k) * dtts * &
            (-adv_fb(i,k+1) + adv_fb(i,k))/Thickness%dzt(i,j,k)
      enddo  !} i
     enddo  !} k

!RASF upstream sinking of caco3
if (id_caco3.ne.0) then
     do k=2,grid%nk+1
      do i=isc,iec
!      adv_fb(i,k)=wcaco3(i,j)*t_prog(ind_caco3)%field(i,j,k-1,time%taum1)
       adv_fb(i,k)=wcaco3(i,j)*biotr(i,k-1,id_caco3)
      enddo
     enddo
! no flux boundary conditions
     do i=isc,iec
      adv_fb(i,1)               = 0.0
     enddo 


! Deposit tracer to sediment as tracer sinks through base of column.  mac, nov12
     do i = isc, iec 
       k = grid%kmt(i,j)
       if (k .gt. 0) then 
       
        biotic(n)%caco3_sed_depst(i,j) = adv_fb(i,k+1)
       
       endif ! k .gt. 0
     enddo ! i


!!mac remineralise tracer sinking through the base of the column.
!     do i = isc, iec 
!       k = grid%kmt(i,j)
!       if (k .gt. 0) then 
!       
!        if (id_dic .ne. 0) then
!         t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
!          grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!          * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!        endif
!      
!        if (id_adic .ne. 0) then
!         t_prog(ind_adic)%th_tendency(i,j,k) =t_prog(ind_adic)%th_tendency(i,j,k) + &
!          grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!          * adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!        endif
!      
!        if (id_alk .ne. 0) then
!         t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) + &
!          2.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!          *adv_fb(i,k+1)/Thickness%dzt(i,j,k)
!        endif
!       endif ! k .gt. 0
!     enddo ! i
!


     do k =1,grid%nk  !{
      do i =isc,iec  !{
        t_prog(ind_caco3)%field(i,j,k,time%taum1)=t_prog(ind_caco3)%field(i,j,k,time%taum1) +&
            grid%tmask(i,j,k) * dtts * &
            (-adv_fb(i,k+1) + adv_fb(i,k))/Thickness%dzt(i,j,k)
      enddo  !} i
     enddo  !} k
endif                   ! end loop for caco3




 enddo  !} j

enddo  !} n

return
end subroutine bio_v3

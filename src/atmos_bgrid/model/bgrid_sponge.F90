
module bgrid_sponge_mod

!-----------------------------------------------------------------------
!
!   Eddy damping of prognostic fields at the top level of the model
!
!   Damping is done using a 5-point Shapiro filter.
!   For temperature, tracers, and zonal wind the zonal mean is
!   removed before applying the filter.  For meridional wind,
!   the entire field is damped.
!
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
 use bgrid_horiz_mod,       only: horiz_grid_type
 use bgrid_masks_mod,       only: grid_mask_type
 use bgrid_prog_var_mod,    only: prog_var_type
 use bgrid_change_grid_mod, only: change_grid, TEMP_GRID, WIND_GRID
 use bgrid_halo_mod,        only: update_halo, vel_flux_boundary, &
                                  TEMP, UWND, VWND, &
                                  NORTH, EAST, NOPOLE

 use         fms_mod, only: error_mesg, FATAL, write_version_number, &
                            file_exist, open_namelist_file, stdlog,  &
                            check_nml_error, close_file, mpp_pe,     &
                            mpp_npes, mpp_root_pe, mpp_clock_id,     &
                            mpp_clock_begin, mpp_clock_end,          &       
                            MPP_CLOCK_SYNC, CLOCK_MODULE

 use         mpp_mod, only: mpp_transmit, mpp_sync_self

 use mpp_domains_mod, only: domain2d, domain1d,         &
                            mpp_update_domains,         &
                            mpp_get_layout,             &
                            mpp_get_pelist,             &
                            mpp_get_global_domain,      &
                            mpp_get_compute_domain,     &
                            mpp_get_compute_domains,    &
                            mpp_get_domain_components,  &
                            WUPDATE, EUPDATE,           &
                            CYCLIC_GLOBAL_DOMAIN

 implicit none
 private

 public   sponge_driver, sponge_init

!-----------------------------------------------------------------------
! namelist

!  num_sponge_levels   The number of uppermost model levels where 
!                      the sponge damping is applied.  Currently,
!                      this cannot exceed one level.   
!                         [integer, default = 0]
!
!  sponge_coeff_wind   Normalized [0,1] sponge damping coefficients
!  sponge_coeff_temp   for the top model level.
!  sponge_coeff_tracer    [real, default = 0.]
!

   integer   ::     num_sponge_levels = 0
   real      ::     sponge_coeff_wind   = 0.0
   real      ::     sponge_coeff_temp   = 0.0
   real      ::     sponge_coeff_tracer = 0.0

   namelist /bgrid_sponge_nml/ num_sponge_levels, sponge_coeff_wind, &
                               sponge_coeff_temp, sponge_coeff_tracer

!-----------------------------------------------------------------------

 character(len=128) :: version='$Id: bgrid_sponge.F90,v 19.0 2012/01/06 19:54:03 fms Exp $'
 character(len=128) :: tagname='$Name: tikal $'
 logical :: do_log  = .true.
 logical :: do_init = .true.
 integer :: id_clock

 real :: small = .000001  !  damping coefficients larger than this
                          !  activate the sponge

!--- module storage for computing exact/reproducible zonal means ---
  type zsum_type
     integer :: is , ie , js , je, isize, jsize, &
                isg, ieg, nlon, imaxsize
     integer, pointer ::  pelist(:) => NULL(), &
                        sizelist(:) => NULL()
  end type


!--- module storage for sponge control parameters ---
  type sponge_control_type
     real :: coeff_vel, coeff_tmp, coeff_trs
     integer :: numlev
     type(zsum_type) :: Zdomain_tmp, Zdomain_vel
  end type sponge_control_type

  type(sponge_control_type),save :: Control

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine sponge_driver ( Hgrid, nplev, dt, dpde, Var, Var_dt )
 
!-----------------------------------------------------------------------
! Hgrid  = horizontal grid constants
! nplev  = number of "pure" pressure levels at the top of the model
! dt     = adjustment time step
! dpde   = pressure weight for model layers
! Var    = prognostic variables at the last updated time level
! Var_dt = tendency of prog variables since the last updated time level
!-----------------------------------------------------------------------
 type (horiz_grid_type), intent(inout) :: Hgrid
 integer,                intent(in)    :: nplev
 real,                   intent(in)    :: dt
 real,                   intent(in)    :: dpde(Hgrid%ilb:,Hgrid%jlb:,:)
 type  (prog_var_type),  intent(in)    :: Var
 type  (prog_var_type),  intent(inout) :: Var_dt
!-----------------------------------------------------------------------

 real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,Control%numlev) :: dpxy
 integer :: n, np, nt

!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('bgrid_sponge_mod',  &
                                 'initialization not called', FATAL)

   if ( Control%numlev == 0 ) return

   call mpp_clock_begin (id_clock)
!-----------------------------------------------------------------------

   n  = Control%numlev   ! number of sponge levels
   np = nplev    ! number of pressure levels at top of model

!---- temperature and tracers -----

   if ( Control%coeff_tmp > small .or. Control%coeff_trs > small ) then
       nt = Var_dt%ntrace
       call local_filter_mass ( Hgrid, np, dt, dpde(:,:,1:n),              &
                                Var   %t(:,:,1:n), Var   %r(:,:,1:n,1:nt), &
                                Var_dt%t(:,:,1:n), Var_dt%r(:,:,1:n,1:nt)  )
   endif

!---- momentum components -----

   if ( Control%coeff_vel > small ) then
      ! compute pressure weights at velocity points
        dpxy(:,:,:) = dpde(:,:,1:n)
        if (np < n) then
          call change_grid (Hgrid, TEMP_GRID, WIND_GRID, &
                            dpxy(:,:,np+1:n), dpxy(:,:,np+1:n))
          call update_halo (Hgrid, UWND, dpxy(:,:,np+1:n), &
                            halos=EAST+NORTH, flags=NOPOLE)
        endif

        call local_filter_vel ( Hgrid, np, dt, dpxy,                  &
                                Var   %u(:,:,1:n), Var   %v(:,:,1:n), &
                                Var_dt%u(:,:,1:n), Var_dt%v(:,:,1:n)  )
   endif

   call mpp_clock_end (id_clock)
!-----------------------------------------------------------------------

 end subroutine sponge_driver

!#######################################################################

 subroutine sponge_init (Hgrid )
 type (horiz_grid_type), intent(in) :: Hgrid ! horizontal grid constants

 integer :: unit, ierr, io, logunit
!-----------------------------------------------------------------------
! read namelist
   if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
       read (input_nml_file, nml=bgrid_sponge_nml, iostat=io)
       ierr = check_nml_error(io,'bgrid_sponge_nml')
#else
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read (unit, nml=bgrid_sponge_nml, iostat=io, end=5)
          ierr = check_nml_error (io, 'bgrid_sponge_nml')
       enddo
 5     call close_file (unit)
#endif
   endif
! write version and namelist to log 
   if (do_log) then
      call write_version_number (version,tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe()) write (logunit, nml=bgrid_sponge_nml)
      do_log = .false.
   endif

! timing routine initialization
   id_clock = mpp_clock_id ('BGRID: sponge', flags=MPP_CLOCK_SYNC, &
                                             grain=CLOCK_MODULE)

! set values for optional arguments
   Control%coeff_vel = min(max(sponge_coeff_wind  ,0.),1.)
   Control%coeff_tmp = min(max(sponge_coeff_temp  ,0.),1.)
   Control%coeff_trs = min(max(sponge_coeff_tracer,0.),1.)
   Control%numlev    = max(num_sponge_levels,0)

!  do not allow more than one sponge layer in this version
   if (Control%numlev > 1) call error_mesg ('bgrid_sponge_mod',  &
                                            'numlev > 1 ', FATAL)

!  set up domain2d types for computing bit-reproducible zonal means

   if ( Control%coeff_tmp > small .or. Control%coeff_trs > small ) then
           call zsum_init ( Hgrid%Tmp%Domain, Control%Zdomain_tmp )
   endif
   if ( Control%coeff_vel > small ) then
           call zsum_init ( Hgrid%Vel%Domain, Control%Zdomain_vel )
   endif

   do_init=.false.

!-----------------------------------------------------------------------

 end subroutine sponge_init

!#######################################################################
! sponge/filter for temperature and tracers fields

 subroutine local_filter_mass ( Hgrid, nplev, dt, dp, t, tr, tdt, trdt )

   type (horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)                   :: nplev
   real,    intent(in)                   :: dt
   real,    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, t
   real,    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: tdt
   real,    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: tr
   real,    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:,:) :: trdt

   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,size(t,3)) :: &
                                             tmp, akew, akns, akdp
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: tew, tns
   real, dimension(Hgrid%jlb:Hgrid%jub) :: akew2, akns2
   integer :: i, j, k, is, ie, js, je, n

   is  = Hgrid%Tmp%is;   ie  = Hgrid%Tmp%ie
   js  = Hgrid%Tmp%js;   je  = Hgrid%Tmp%je

 ! 2D geometric constants
   do j = js-1, je
      akew2(j) = 0.0625 * (Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j))
      akns2(j) = 0.0625 * (Hgrid%Tmp%area(j)+Hgrid%Tmp%area(j+1))
   enddo

 ! 3D mass weighted constants
   do k = 1, size(t,3)
   do j = js-1, je
      akew(is-1:ie,j,k) = akew2(j)
      akns(is-1:ie,j,k) = akns2(j)
      akdp(:,j,k) = Hgrid%Tmp%rarea(j)/dt
   enddo
   enddo
   do k = nplev+1, size(t,3)
      do j = js-1, je
      do i = is-1, ie
         akew(i,j,k) = akew(i,j,k) * (dp(i,j,k)+dp(i+1,j,k))
         akns(i,j,k) = akns(i,j,k) * (dp(i,j,k)+dp(i,j+1,k))
         akdp(i,j,k) = akdp(i,j,k) / (2.*dp(i,j,k))
      enddo
      enddo
   enddo


 ! temperature

   if ( Control%coeff_tmp > small ) then
      tmp = t + dt*tdt
     !---- remove zonal mean ----
      call remove_mean ( Control%Zdomain_tmp, tmp(is:ie,js:je,:) )
      call update_halo ( Hgrid, TEMP, tmp )
      do k = 1, size(t,3)
         do j = js-1, je
         do i = is-1, ie
            tew(i,j) = akew(i,j,k) * (tmp(i+1,j,k)-tmp(i,j,k))
            tns(i,j) = akns(i,j,k) * (tmp(i,j+1,k)-tmp(i,j,k))
         enddo
         enddo
         do j = js, je
         do i = is, ie
            tdt(i,j,k) = tdt(i,j,k) + &
                  Control%coeff_tmp*(tew(i,j)-tew(i-1,j)+tns(i,j)-tns(i,j-1))*akdp(i,j,k)
         enddo
         enddo
      enddo
   endif

 ! tracers

   if ( Control%coeff_trs > small ) then
      do n = 1, size(tr,4)
         tmp = tr(:,:,:,n) + dt*trdt(:,:,:,n)
        !---- remove zonal mean ----
         call remove_mean ( Control%Zdomain_tmp, tmp(is:ie,js:je,:) )
         call update_halo ( Hgrid, TEMP, tmp )
         do k = 1, size(tr,3)
            do j = js-1, je
            do i = is-1, ie
               tew(i,j) = akew(i,j,k) * (tmp(i+1,j,k)-tmp(i,j,k))
               tns(i,j) = akns(i,j,k) * (tmp(i,j+1,k)-tmp(i,j,k))
            enddo
            enddo
            do j = js, je
            do i = is, ie
               trdt(i,j,k,n) = trdt(i,j,k,n) +  &
                    Control%coeff_trs*(tew(i,j)-tew(i-1,j)+tns(i,j)-tns(i,j-1))*akdp(i,j,k)
            enddo
            enddo
         enddo
      enddo
   endif

 end subroutine local_filter_mass

!#######################################################################
! sponge/filter for momentum components

 subroutine local_filter_vel ( Hgrid, nplev, dt, dp, u, v, udt, vdt )

   type (horiz_grid_type), intent(inout) :: Hgrid
   integer, intent(in)                   :: nplev
   real,    intent(in)                   :: dt
   real,    intent(in),    dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: dp, u, v
   real,    intent(inout), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: udt, vdt

   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub,size(u,3)) :: &
                                                akew, akns, akdp, uu, vv
   real, dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::  uew, vew, uns, vns
   real, dimension(Hgrid%jlb:Hgrid%jub) ::  akew2, akns2
   integer :: i, j, k, is, ie, js, je

   is  = Hgrid%Tmp%is;   ie  = Hgrid%Tmp%ie
   js  = Hgrid%Vel%js;   je  = Hgrid%Vel%je

 ! 2D geometric constants
   do j = js, je+1
      akew2(j) = 0.0625 * Control%coeff_vel * (Hgrid%Vel%area(j)+Hgrid%Vel%area(j))
      akns2(j) = 0.0625 * Control%coeff_vel * (Hgrid%Vel%area(j-1)+Hgrid%Vel%area(j))
   enddo

 ! 3D mass weighted constants
   do k = 1, size(u,3)
   do j = js, je+1
      akew(is:ie+1,j,k) = akew2(j)
      akns(is:ie+1,j,k) = akns2(j)
      akdp(:,j,k) = Hgrid%Vel%rarea(j)/dt
   enddo
   enddo
   do k = nplev+1, size(u,3)
      do j = js, je+1
      do i = is, ie+1
         akew(i,j,k) = akew(i,j,k) * (dp(i,j,k)+dp(i-1,j,k))
         akns(i,j,k) = akns(i,j,k) * (dp(i,j,k)+dp(i,j-1,k))
         akdp(i,j,k) = akdp(i,j,k) / (2.*dp(i,j,k))
      enddo
      enddo
   enddo

   uu = u + dt*udt
   vv = v + dt*vdt
  !---- remove zonal mean from u comp ----
   call remove_mean ( Control%Zdomain_vel, uu(is:ie,js:je,:) )
   call update_halo ( Hgrid, UWND, uu )
   call update_halo ( Hgrid, VWND, vv )

   do k = 1, size(u,3)
      do j = js, je+1
      do i = is, ie+1
         uew(i,j) = akew(i,j,k) * (uu(i,j,k)-uu(i-1,j,k))
         vew(i,j) = akew(i,j,k) * (vv(i,j,k)-vv(i-1,j,k))
         uns(i,j) = akns(i,j,k) * (uu(i,j,k)-uu(i,j-1,k))
         vns(i,j) = akns(i,j,k) * (vv(i,j,k)-vv(i,j-1,k))
      enddo
      enddo
     !---- remove meridional gradients adjacent to poles ----
      call vel_flux_boundary (Hgrid, uns)
      call vel_flux_boundary (Hgrid, vns)

      do j = js, je
      do i = is, ie
         udt(i,j,k) = udt(i,j,k) + (uew(i+1,j)-uew(i,j)+uns(i,j+1)-uns(i,j))*akdp(i,j,k)
         vdt(i,j,k) = vdt(i,j,k) + (vew(i+1,j)-vew(i,j)+vns(i,j+1)-vns(i,j))*akdp(i,j,k)
      enddo
      enddo
   enddo

 end subroutine local_filter_vel

!#######################################################################
! initializes domain2d type for summation in zonal direction

 subroutine zsum_init ( Domain, Zonal )
 type(domain2d),  intent(in)   :: Domain
 type(zsum_type), intent(out)  :: Zonal
 integer :: layout(2)
 type(domain1D) :: Domx, Domy

! create new domain2d type with large global halo along x-axis 
  call mpp_get_layout ( Domain, layout )
  allocate ( Zonal%pelist(layout(1)), Zonal%sizelist(layout(1)) )

  call mpp_get_domain_components ( Domain, Domx, Domy )
  call mpp_get_pelist            ( Domx, Zonal%pelist )
  call mpp_get_compute_domains   ( Domx, size=Zonal%sizelist )

! get compute domain 
  call mpp_get_global_domain  ( Domain, Zonal%isg, Zonal%ieg  )
  call mpp_get_compute_domain ( Domain, Zonal%is , Zonal%ie , &
                                        Zonal%js , Zonal%je   )
  Zonal%imaxsize = maxval(Zonal%sizelist)
  Zonal%isize = Zonal%ie-Zonal%is+1
  Zonal%jsize = Zonal%je-Zonal%js+1
  Zonal%nlon  = Zonal%ieg-Zonal%isg+1

 end subroutine zsum_init

!#######################################################################

 subroutine remove_mean ( Zonal, local )
 type(zsum_type), intent(in)   :: Zonal
 real,            intent(inout):: local(Zonal%is:,Zonal%js:,:) ! compute domain only
 real, dimension(Zonal%js:Zonal%je,size(local,3)) :: zsum
 real, dimension(Zonal%imaxsize*Zonal%jsize*size(local,3)) :: data_get
 real, dimension(Zonal%isize*Zonal%jsize*size(local,3)) :: data_put
 integer :: i, j, k, m, n, is, ie, npts_put, npts_get

   zsum = 0.
   ie = Zonal%isg - 1
 ! loop thru PEs in zonal direction
   do n = 1, size(Zonal%pelist(:))
      is = ie + 1
      ie = is + Zonal%sizelist(n) - 1
      if (Zonal%pelist(n) == mpp_pe()) then
        ! data is local to current PE
          do k = 1, size(local,3)
          do j = Zonal%js, Zonal%je
          do i = is, ie
             zsum(j,k) = zsum(j,k) + local(i,j,k)
          enddo
          enddo
          enddo
      else
        ! reshape input array into 1d array (data must be contiguous)
          m = 0
          do k = 1, size(local,3)
          do j = Zonal%js, Zonal%je
          do i = Zonal%is, Zonal%ie
             m = m+1
             data_put(m) = local(i,j,k)
          enddo
          enddo
          enddo
          npts_put = size(data_put(:))
          npts_get = Zonal%sizelist(n)*size(local,2)*size(local,3)
        ! data is not local (communication required)
          call mpp_transmit ( put_data=data_put, put_len=npts_put,   to_pe=Zonal%pelist(n), &
                              get_data=data_get, get_len=npts_get, from_pe=Zonal%pelist(n)  )
        ! reshape (with summation)
          m = 0
          do k = 1, size(local,3)
          do j = Zonal%js, Zonal%je
          do i = is, ie
             m = m+1
             zsum(j,k) = zsum(j,k) + data_get(m)
          enddo
          enddo
          enddo
      endif
   enddo

 ! remove zonal mean
   do k = 1, size(local,3)
   do j = Zonal%js, Zonal%je
   do i = Zonal%is, Zonal%ie
      local(i,j,k) = local(i,j,k) - zsum(j,k)/real(Zonal%nlon)
   enddo
   enddo
   enddo

   ! required after a call to mpp_transmit
   call mpp_sync_self()
      
 end subroutine remove_mean

!#######################################################################

end module bgrid_sponge_mod


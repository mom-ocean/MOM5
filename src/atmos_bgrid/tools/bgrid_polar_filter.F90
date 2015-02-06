
module bgrid_polar_filter_mod
!-----------------------------------------------------------------------

use bgrid_horiz_mod, only: horiz_grid_type, bgrid_type, &
                           TGRID, VGRID

use         fft_mod, only: fft_init, fft_grid_to_fourier,  &
                                     fft_fourier_to_grid

use         fms_mod, only: error_mesg, FATAL,    &
                           write_version_number, &
                           mpp_npes, mpp_pe, mpp_root_pe, &
                           mpp_clock_id, mpp_clock_begin, &
                           mpp_clock_end, MPP_CLOCK_SYNC, &
                           CLOCK_ROUTINE

use mpp_domains_mod, only: domain2d, domain1d,         &
                           mpp_redistribute,           &
                           mpp_define_layout,          &
                           mpp_define_domains,         &
                           mpp_get_domain_components,  &
                           mpp_get_layout,             &
                           mpp_get_global_domain,      &
                           mpp_get_data_domain,        &
                           mpp_get_compute_domains,    &
                           CYCLIC_GLOBAL_DOMAIN

implicit none
private

 public :: polar_filter, polar_filter_wind, polar_filter_init
 public :: TGRID, VGRID

 interface polar_filter
    module procedure polar_filter_3d, polar_filter_two_3d, &
                     polar_filter_2d, polar_filter_two_2d
 end interface polar_filter

 interface polar_filter_wind
    module procedure polar_filter_wind_3d, polar_filter_wind_2d
 end interface polar_filter_wind

!--------------------------------
! derived types
 public :: pfilt_domain_type, pfilt_index_type, pfilt_control_type

 type pfilt_domain_type
    private
    type (domain2d) :: Domain
    integer :: isd, ied, jsd, jed
 end type pfilt_domain_type

 type pfilt_index_type
    private
    integer    :: leng, lenc
    integer    :: is, ie, js, je
    integer    :: isd, ied, jsd, jed
    integer    :: weight
    logical    :: sigma, rowend
    real       :: cph0
    real,    pointer :: sklm(:) =>NULL(), cph(:) =>NULL()
    integer ,pointer :: jpf(:) =>NULL()
    type (pfilt_domain_type) :: Local, Zonal, Local2, Zonal2
 end type pfilt_index_type

 type pfilt_control_type
    private
    integer    :: nlpf, nlev
    integer    :: isize, jsize
    real, pointer, dimension(:) :: slm =>NULL(), clm =>NULL()
    type (pfilt_index_type) :: Tmp, Vel
 end type pfilt_control_type
   
!--------------------------------
! version id info
 character(len=128) :: version='$Id: bgrid_polar_filter.F90,v 10.0 2003/10/24 22:00:19 fms Exp $'
 character(len=128) :: tagname='$Name: tikal $'
 logical :: do_log = .true.
!--------------------------------
! private data
 real   , allocatable :: xsum(:)
 integer, allocatable :: xknt(:)

! performance timing data
 integer :: id_pfiltr
 logical :: do_clock_init = .true.

contains

!#######################################################################
! public routines for polaring filtering
!#######################################################################

 subroutine polar_filter_3d ( Control, u, grid, mask )

!-----------------------------------------------------
!     Polar filter a single 3D field
!
! Control = control parameters for polar filter
! u       = 3D data
! grid    = grid identifier, use: TGRID, VGRID
! mask    = grid box mask for eta coordinate topography
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)               :: Control
 real,    intent(inout), dimension(:,:,:)           :: u
 integer, intent(in)                                :: grid
 real,    intent(in),    dimension(:,:,:), optional :: mask

 ! error checks
   if (size(u,1) /= Control%isize .or. size(u,2) /= Control%jsize) &
   call error_mesg ('polar_filter_mod', 'incorrect horizontal dimension', FATAL)
   if (size(u,3) /= Control%nlev) &
   call error_mesg ('polar_filter_mod', 'incorrect number of levels', FATAL )

 ! select grid
   select case (grid)
     case (TGRID)
        call filter_field ( Control%Tmp, u, mask=mask )
     case (VGRID)
        call filter_field ( Control%Vel, u, mask=mask )
     case default
        call error_mesg ('polar_filter_mod', 'invalid grid argument', FATAL)
   end select

 end subroutine polar_filter_3d

!=======================================================================

 subroutine polar_filter_two_3d ( Control, u, v, grid, mask )

!-----------------------------------------------------
!   Polar filter two 3D fields on the same grid
!
! Control = control parameters for polar filter
! u, v    = 3D data field on the same grid
! grid    = grid identifier, use: TGRID, VGRID
! mask    = grid box mask for eta coordinate topography
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)               :: Control
 real,    intent(inout), dimension(:,:,:)           :: u, v
 integer, intent(in)                                :: grid
 real,    intent(in),    dimension(:,:,:), optional :: mask

 ! error checks
   if (size(u,1) /= Control%isize .or. size(u,2) /= Control%jsize) &
   call error_mesg ('polar_filter_mod', 'incorrect horizontal dimension', FATAL)
   if (size(u,3) /= Control%nlev .or. size(v,3) /= Control%nlev) &
   call error_mesg ('polar_filter_mod', 'incorrect number of levels', FATAL )

 ! select grid
   select case (grid)
     case (TGRID)
        call filter_two_fields ( Control%Tmp, u, v, mask=mask )
     case (VGRID)
        call filter_two_fields ( Control%Vel, u, v, mask=mask )
     case default
        call error_mesg ('polar_filter_mod', 'invalid grid argument', FATAL)
   end select

 end subroutine polar_filter_two_3d

!=======================================================================

 subroutine polar_filter_wind_3d ( Control, u, v, mask )

!-----------------------------------------------------
!      Polar filter 3D momentum fields
!
! Control = control parameters for polar filter
! u, v    = 3D data field on the same grid
! mask    = grid box mask for eta coordinate topography
!
! NOTE: a polar stereographic transformation of the
!       wind components done when the filter is applied
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)               :: Control
 real,    intent(inout), dimension(:,:,:)           :: u, v
 real,    intent(in),    dimension(:,:,:), optional :: mask

 ! error checks
   if (size(u,1) /= Control%isize .or. size(u,2) /= Control%jsize) &
   call error_mesg ('polar_filter_mod', 'incorrect horizontal dimension', FATAL)
   if (size(u,3) /= Control%nlev .or. size(v,3) /= Control%nlev) &
   call error_mesg ('polar_filter_mod', 'incorrect number of levels', FATAL )

   call filter_two_fields ( Control%Vel, u, v, &
                            slm=Control%slm, clm=Control%clm, mask=mask )

 end subroutine polar_filter_wind_3d

!#######################################################################
!=============== overloaded 2D interfaces ====================

 subroutine polar_filter_2d ( Control, u, grid, mask )

!-----------------------------------------------------
!     Polar filter a single 2D field
!
! Control = control parameters for polar filter
! u       = 2D data
! grid    = grid identifier, use: TGRID, VGRID
! mask    = grid box mask for eta coordinate topography
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)             :: Control
 real,    intent(inout), dimension(:,:)           :: u
 integer, intent(in)                              :: grid
 real,    intent(in),    dimension(:,:), optional :: mask

 real, dimension(size(u,1),size(u,2),1) :: u3, m3

   u3(:,:,1) = u
   if (present(mask)) then
       m3(:,:,1) = mask
       call polar_filter_3d ( Control, u3, grid, mask=m3 )
   else
       call polar_filter_3d ( Control, u3, grid )
   endif
   u = u3(:,:,1)

 end subroutine polar_filter_2d

!=======================================================================

 subroutine polar_filter_two_2d ( Control, u, v, grid, mask )

!-----------------------------------------------------
!   Polar filter two 2D fields on the same grid
!
! Control = control parameters for polar filter
! u, v    = 2D data field on the same grid
! grid    = grid identifier, use: TGRID, VGRID
! mask    = grid box mask for eta coordinate topography
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)             :: Control
 real,    intent(inout), dimension(:,:)           :: u, v
 integer, intent(in)                              :: grid
 real,    intent(in),    dimension(:,:), optional :: mask

 real, dimension(size(u,1),size(u,2),1) :: u3, v3, m3

   u3(:,:,1) = u
   v3(:,:,1) = v
   if (present(mask)) then
       m3(:,:,1) = mask
       call polar_filter_two_3d ( Control, u3, v3, grid, mask=m3 )
   else
       call polar_filter_two_3d ( Control, u3, v3, grid )
   endif
   u = u3(:,:,1)
   v = v3(:,:,1)

 end subroutine polar_filter_two_2d

!=======================================================================

 subroutine polar_filter_wind_2d ( Control, u, v, mask )

!-----------------------------------------------------
!      Polar filter 2D momentum fields
!
! Control = control parameters for polar filter
! u, v    = 2D data field on the same grid
! mask    = grid box mask for eta coordinate topography
!
! NOTE: a polar stereographic transformation of the
!       wind components done when the filter is applied
!-------------------------------------------------------

 type(pfilt_control_type), intent(in)             :: Control
 real,    intent(inout), dimension(:,:)           :: u, v
 real,    intent(in),    dimension(:,:), optional :: mask

 real, dimension(size(u,1),size(u,2),1) :: u3, v3, m3

 ! error checks
   if (size(u,1) /= Control%isize .or. size(u,2) /= Control%jsize) &
   call error_mesg ('polar_filter_mod', 'incorrect horizontal dimension', FATAL)

   u3(:,:,1) = u
   v3(:,:,1) = v
   if (present(mask)) then
       m3(:,:,1) = mask
       call filter_two_fields ( Control%Vel, u3, v3, &
                            slm=Control%slm, clm=Control%clm, mask=m3 )
   else
       call filter_two_fields ( Control%Vel, u3, v3, &
                            slm=Control%slm, clm=Control%clm )
   endif
   u = u3(:,:,1)
   v = v3(:,:,1)

 end subroutine polar_filter_wind_2d

!#######################################################################
! filter a single field

 subroutine filter_field ( Index, u, mask )

 type(pfilt_index_type), intent(in)     :: Index
 real, intent(inout), dimension(Index%isd:,Index%jsd:,:) :: u
 real, intent(in),    dimension(Index%isd:,Index%jsd:,:), &
                                                optional :: mask

 real, dimension( Index%Local%isd:Index%Local%ied, &
                  Index%Local%jsd:Index%Local%jed ) :: g, gm
 real, dimension( Index%Zonal%isd:Index%Zonal%ied, &
                  Index%Zonal%jsd:Index%Zonal%jed ) :: z, zm

 real,    dimension( Index%lenc, Index%Zonal%jsd:Index%Zonal%jed ) :: ss
 complex, dimension( Index%lenc, Index%Zonal%jsd:Index%Zonal%jed ) :: c

 integer :: i, j, k, n, nlev, is, ie, isg, ieg
 logical :: use_mask
 real    :: zm_min

  call mpp_clock_begin (id_pfiltr)

! scalar constants
  nlev = size(u,3)
  is  = Index%is;        ie  = Index%ie
  isg = Index%Zonal%isd; ieg = Index%Zonal%ied-1
  use_mask = .not.Index%sigma .and. present(mask)
  zm_min = 1.0

!--------------------------------------------
! reorder indexing i*k*j

     n = Index%Local%jsd - 1
  do j = Index%js, Index%je
  do k = 1, nlev
     if (Index%jpf(j) /= 0) then
        n = n + 1
        do i = is, ie
          g(i,n) = u(i,j,k)
        enddo
        ! add cosine latitude to end of rows
        if (Index%rowend) g(ie+1,n) = Index%cph(j)
        ! step-mountain mask
        if (use_mask) then
          do i = is, ie
            gm(i,n) = mask(i,j,k)
          enddo
        endif
     endif
  enddo
  enddo
!--------------------------------------------
! distribute the data across all processors
  call mpp_redistribute ( Index%Local%Domain, g, Index%Zonal%Domain, z )
! mask ?
  if (use_mask) then
    call mpp_redistribute ( Index%Local%Domain, gm, Index%Zonal%Domain, zm )
    zm_min = minval(zm(isg:ieg,:))
    if (zm_min < .01) call fill_missing (zm(isg:ieg,:), z(isg:ieg,:))
  endif

! compute filter response
  call set_filter_response ( Index, z(ieg+1,:), ss )
! transform to fourier coefficients
  c = fft_grid_to_fourier (z)
! filter
  c = c * ss
! transform back to grid spce
  z = fft_fourier_to_grid (c)
! restore zonal mean if mask used
  if (zm_min < .01) call fix_missing (zm(isg:ieg,:), z(isg:ieg,:))

! distribute the data back to original processor
  call mpp_redistribute ( Index%Zonal%Domain, z, Index%Local%Domain, g )

!--------------------------------------------
! place data back in original arrays

     n = Index%Local%jsd - 1
  do j = Index%js, Index%je
  do k = 1, nlev
     ! only replace filtered data
     if (Index%jpf(j) /= 0) then
        n = n + 1
        do i = is, ie
          u(i,j,k) = g(i,n)
        enddo
     endif
     ! apply mask for step-mountain coord
     if (use_mask) then
       do i = is, ie
         u(i,j,k) = u(i,j,k) * mask(i,j,k)
       enddo
     endif
  enddo
  enddo

  call mpp_clock_end (id_pfiltr)
!--------------------------------------------

end subroutine filter_field

!#######################################################################
! filter two fields on the same grid

 subroutine filter_two_fields ( Index, u, v, slm, clm, mask )

 type(pfilt_index_type), intent(in)     :: Index
 real, intent(inout), dimension(Index%isd:,Index%jsd:,:) :: u, v
 real, intent(in),    dimension(Index%isd:), optional    :: slm, clm
 real, intent(in),    dimension(Index%isd:,Index%jsd:,:), &
                                                optional :: mask

 real, dimension( Index%Local2%isd:Index%Local2%ied, &
                  Index%Local2%jsd:Index%Local2%jed ) :: g, gm
 real, dimension( Index%Zonal2%isd:Index%Zonal2%ied, &
                  Index%Zonal2%jsd:Index%Zonal2%jed ) :: z, zm

 real,    dimension( Index%lenc, Index%Zonal2%jsd:Index%Zonal2%jed ) :: ss
 complex, dimension( Index%lenc, Index%Zonal2%jsd:Index%Zonal2%jed ) :: c

 integer :: i, j, k, n1, n2, nlev, nfld, is, ie, isg, ieg
 logical :: vectors, use_mask
 real    :: zm_min

  call mpp_clock_begin (id_pfiltr)

! initialize scalars
  nlev = size(u,3)
  nfld = 2
  is  = Index%is;         ie  = Index%ie
  isg = Index%Zonal2%isd; ieg = Index%Zonal2%ied-1
  vectors = present(slm) .and. present(clm)
  use_mask = .not.Index%sigma .and. present(mask)
  zm_min = 1.0

!--------------------------------------------
! reorder indexing i*k*j
! add cosine latitude to end of rows

     n1 = Index%Local2%jsd - 2
  do j = Index%js, Index%je
  do k = 1, nlev
     if (Index%jpf(j) /= 0) then
        n1 = n1 + 2
        n2 = n1 + 1
        if (vectors) then
           ! convert u/v components using streographic projection
           select case (Index%jpf(j))
             case (-1) ! so.hemis.
               do i = is, ie
                 g(i,n1) = -u(i,j,k)*slm(i) + v(i,j,k)*clm(i)
                 g(i,n2) = +u(i,j,k)*clm(i) + v(i,j,k)*slm(i)
               enddo
             case (+1) ! no.hemis.
               do i = is, ie
                 g(i,n1) = -u(i,j,k)*slm(i) - v(i,j,k)*clm(i)
                 g(i,n2) = +u(i,j,k)*clm(i) - v(i,j,k)*slm(i)
               enddo
           end select
        else
           ! two distinct fields - no conversion necessary
           do i = is, ie
             g(i,n1) = u(i,j,k)
             g(i,n2) = v(i,j,k)
           enddo
        endif
        ! insert cosine latitude at end of global latitude rows
        if (Index%rowend) then
            g(ie+1,n1) = Index%cph(j)
            g(ie+1,n2) = Index%cph(j)
        endif
        ! step-mountain mask (need two copies)
        if (use_mask) then
            do i = is, ie
              gm(i,n1) = mask(i,j,k)
              gm(i,n2) = mask(i,j,k)
            enddo
        endif
     endif
  enddo
  enddo
!--------------------------------------------
! distribute the data across all processors
  call mpp_redistribute ( Index%Local2%Domain, g, Index%Zonal2%Domain, z )
! mask ?
  if (use_mask) then
    call mpp_redistribute ( Index%Local2%Domain, gm, Index%Zonal2%Domain, zm )
    zm_min = minval(zm(isg:ieg,:))
    if (zm_min < .01) call fill_missing (zm(isg:ieg,:), z(isg:ieg,:))
  endif

! compute filter response
  call set_filter_response ( Index, z(ieg+1,:), ss )
! transform to fourier coefficients
  c = fft_grid_to_fourier (z)
! filter
  c = c * ss
! transform back to grid spce
  z = fft_fourier_to_grid (c)
! restore zonal mean if mask used
  if (zm_min < .01) call fix_missing (zm(isg:ieg,:), z(isg:ieg,:))

! distribute the data back to original processor
  call mpp_redistribute ( Index%Zonal2%Domain, z, Index%Local2%Domain, g )

!--------------------------------------------
! place data back in original arrays

     n1 = Index%Local2%jsd - 2
  do j = Index%js, Index%je
  do k = 1, nlev
     if (Index%jpf(j) /= 0) then
        n1 = n1 + 2
        n2 = n1 + 1
        if (vectors) then
           ! convert stereographic projection back to u/v components
           select case (Index%jpf(j))
             case (-1) ! so. hemis.
               do i = is, ie
                 u(i,j,k) = -g(i,n1)*slm(i) + g(i,n2)*clm(i)
                 v(i,j,k) = +g(i,n1)*clm(i) + g(i,n2)*slm(i)
               enddo
             case (+1) ! no. hemis.
               do i = is, ie
                 u(i,j,k) = -g(i,n1)*slm(i) + g(i,n2)*clm(i)
                 v(i,j,k) = -g(i,n1)*clm(i) - g(i,n2)*slm(i)
               enddo
           end select
        else
           ! two distinct fields
           do i = is, ie
             u(i,j,k) = g(i,n1)
             v(i,j,k) = g(i,n2)
           enddo
        endif
     endif
     ! apply mask for step-mountain coord
     if (use_mask) then
       do i = is, ie
         u(i,j,k) = u(i,j,k) * mask(i,j,k)
         v(i,j,k) = v(i,j,k) * mask(i,j,k)
       enddo
     endif
  enddo
  enddo

  call mpp_clock_end (id_pfiltr)
!--------------------------------------------

end subroutine filter_two_fields

!#######################################################################
! computes standard filter response (after arakawa & lamb)

 subroutine set_filter_response ( Index, cph, ss )
 type(pfilt_index_type), intent(in)  :: Index
 real                  , intent(in)  :: cph (:)
 real                  , intent(out) :: ss (:,:)

 real, dimension(Index%lenc) :: cph0_sklm
 integer :: j, k

  ! longitude dependent values
    do k = 2, Index%lenc
      cph0_sklm (k) = Index%cph0 * Index%sklm(k)
    enddo   

  ! mean always one
    ss (1,:) = 1.0

  ! compute response function (range 0. to 1.)
    do j = 1, size(ss,2)
       do k = 2, Index%lenc
          ss (k,j) = max(0.0, min( 1.0, cph(j)/cph0_sklm(k) ))
       enddo   
    enddo   

 end subroutine set_filter_response

!#######################################################################
!======== routines for step-mountain coordinate =========
!#######################################################################

subroutine fill_missing ( mask, dat )
real, intent(in),    dimension(:,:) :: mask
real, intent(inout), dimension(:,:) :: dat
integer :: i, j, nlon, nlat

! computes zonal mean at points where mask /= 0
! interpolates values into points where mask = 0

  nlon = size(dat,1);  nlat = size(dat,2)
  allocate ( xknt(nlat), xsum(nlat) )
  xknt = 0; xsum = 0.

  do j = 1, nlat
    do i = 1, nlon
      if (mask(i,j) < .01) cycle
      xknt(j) = xknt(j) + 1
      xsum(j) = xsum(j) + dat(i,j)
    enddo
    if (xknt(j) == nlon .or. xknt(j) == 0) cycle
    call intrp ( mask(:,j), dat(:,j) )
   !where (mask(:,j)<.01) dat(:,j)=xsum(j)/real(xknt(j))
  enddo

end subroutine fill_missing

!=======================================================================

subroutine fix_missing ( mask, dat )
real, intent(in),    dimension(:,:) :: mask
real, intent(inout), dimension(:,:) :: dat
integer :: i, j, nlon, nlat
real :: avg, dif

! restores zonal mean in latitude rows with missing values

  nlon = size(dat,1);  nlat = size(dat,2)

  do j = 1, nlat
     if (xknt(j) == nlon .or. xknt(j) == 0) cycle
     avg = 0.
     do i = 1, nlon
        if (mask(i,j) < .01) cycle
        avg = avg + dat(i,j)
     enddo
     dif = (xsum(j)-avg)/real(xknt(j))
     do i = 1, nlon
        if (mask(i,j) < .01) cycle
        dat(i,j) = dat(i,j) + dif
     enddo
  enddo
  deallocate ( xknt, xsum )

end subroutine fix_missing

!=======================================================================

subroutine intrp ( mask, dat )
real, intent(in),    dimension(:) :: mask
real, intent(inout), dimension(:) :: dat

integer, dimension(size(dat,1)) :: m1, m2, mbas
real,    dimension(size(dat,1)) :: base, slop
integer :: nlon, last, nseg, i, m, n

!  fill in missing values by linear interpolating

  m1(:)=99999
  m2(:)=99999
  mbas(:)=99999
  base(:)=99999.
  slop(:)=99999.

  nlon = size(dat,1)
  last = -1
  nseg =  1

  do i = 1, nlon
     if (mask(i) < .01) then
         if (last == 1) then 
           m1(nseg) = i
           mbas(nseg) = i-1
           base(nseg) = dat(i-1)
         endif
         last = 0
     else
         if (last == 0) then
           m2(nseg) = i-1
           slop(nseg) = (dat(i)-base(nseg))/real(i-mbas(nseg))
           nseg = nseg+1
         endif
         last = 1
     endif
  enddo

  if ( m1(nseg) == 99999 .and. m2(nseg) == 99999) nseg = nseg-1
  if ( m1(1)    == 99999 .and. m2(nseg) == 99999) then
       m1(1) = 1
       m2(nseg) = nlon
       mbas(1) = mbas(nseg)-nlon
       base(1) = base(nseg)
       slop(1) = (dat(m2(1)+1)-base(1))/real(m2(1)+1-mbas(1))
       slop(nseg) = slop(1)
  endif
  if (m1(1) == 99999) then
      m1(1) = 1
      mbas(1) = 0
      base(1) = dat(nlon)
      slop(1) = (dat(m2(1)+1)-base(1))/real(m2(1)+1-mbas(1))
  endif
  if (m2(nseg) == 99999) then
      m2(nseg) = nlon
      slop(nseg) = (dat(1)-base(nseg))/real(nlon+1-mbas(nseg))
  endif

  do n = 1, nseg
  do m = m1(n), m2(n)
       dat(m) = base(n) + real(m-mbas(n))*slop(n)
  enddo
  enddo

end subroutine intrp

!#######################################################################
!   initialization routines
!#######################################################################
 
 subroutine polar_filter_init ( Control, Hgrid, nlev, reflat, weight, sigma, verbose )

!-----------------------------------------------------------------------
!  Hgrid  = horizontal grid constants
!  nlev   = number of vertical model levels
!             all input data must have this number of levels
!  reflat = reference latitude in degrees, default=60.,
!             poleward of this latitude the filter is applied
!  weight = weight to strengthen filter (not recommended)
!  sigma  = flag to improve optimization for sigma coordinate models
!             default=FALSE, set to TRUE for sigma models
!  verbose = not used?
!-----------------------------------------------------------------------

 type(pfilt_control_type), intent(inout)        :: Control
 type(horiz_grid_type),    intent(in)           :: Hgrid
 integer,                  intent(in)           :: nlev
 real,                     intent(in), optional :: reflat
 integer,                  intent(in), optional :: weight, verbose 
 logical,                  intent(in), optional :: sigma

 integer :: nlpf, i
 real    :: filter_lats
!-----------------------------------------------------------------------

 if (do_log) then
   call write_version_number (version,tagname)
   do_log = .false. 
 endif

  ! reference latitude
    filter_lats = 30.
    if (present(reflat)) then
         if ( reflat-epsilon(reflat) <= 90. .and. reflat >= 0. ) then
           filter_lats = 90. - reflat
         else    
           call error_mesg  ('polar_filter_init',      &
              'reflat must lie between 0 and 90.', FATAL)
        endif
    endif
  ! number of latitude rows of filtering per hemisphere
    nlpf = int(float(Hgrid%Tmp%jeg-Hgrid%Tmp%jsg+1)*filter_lats/180.)
    Control%nlpf = nlpf
    Control%nlev = nlev
  ! grid size (for error checking)
    Control%isize = Hgrid%isize
    Control%jsize = Hgrid%jsize

  ! Temp/Vel grid constants
    call set_index_type ( nlpf, nlev, Hgrid%dlm, Hgrid%dph, Hgrid%Tmp, &
                          Control%Tmp, weight, sigma )
    call set_index_type ( nlpf, nlev, Hgrid%dlm, Hgrid%dph, Hgrid%Vel, &
                          Control%Vel, weight, sigma )

  ! fourier transform initialization
    call fft_init (Hgrid%nlon)

  ! trig constants for converting u,v to polar stereographic
    allocate ( Control%slm(Control%Vel%isd:Control%Vel%ied), &
               Control%clm(Control%Vel%isd:Control%Vel%ied)  )
    do i = Control%Vel%isd, Control%Vel%ied
       Control%slm(i) = sin( Hgrid%Vel%tlm(i) )
       Control%clm(i) = cos( Hgrid%Vel%tlm(i) )
    enddo

  ! initialize performance clock
    if (do_clock_init) then
       id_pfiltr = mpp_clock_id ('BGRID: polar_filter (TOTAL)',  &
                        flags=MPP_CLOCK_SYNC, grain=CLOCK_ROUTINE)
       do_clock_init = .false.
    endif

end subroutine polar_filter_init
 
!=======================================================================

 subroutine set_index_type ( nlpf, nlev, dlm, dph, Grid, Index, weight, sigma )
 integer,                intent(in)  :: nlpf, nlev
 real,                   intent(in)  :: dlm, dph
 type(bgrid_type),       intent(in)  :: Grid
 type(pfilt_index_type), intent(inout) :: Index
 integer, optional,      intent(in)  :: weight
 logical, optional,      intent(in)  :: sigma

 integer :: leng, lenc, j, k
 real    :: hpi, rlat

 ! optional arguments
   Index%weight = 1;       if (present(weight)) Index%weight = weight
   Index%sigma = .false.;  if (present(sigma))  Index%sigma  = sigma

 ! east-most box in zonal row needs one more point
 ! this is used to store the latitude value
   Index % rowend = Grid%ie == Grid%ieg
 ! for full zonal rows
   leng = Grid%ieg - Grid%isg + 1
   lenc = leng/2+1
 
   Index % leng = leng
   Index % lenc = lenc

 ! copy indexing for this domain
   Index % is = Grid%is
   Index % ie = Grid%ie
   Index % js = Grid%js
   Index % je = Grid%je
   Index % isd = Grid%isd
   Index % ied = Grid%ied
   Index % jsd = Grid%jsd
   Index % jed = Grid%jed

 ! reference latitude (first lat w/o filtering)
   hpi = acos(0.0)
   rlat = hpi - (float(nlpf)+0.50)*dph
   Index % cph0 = cos(rlat)

   allocate ( Index%jpf(Index%js:Index%je), &
              Index%cph(Index%js:Index%je)  )
   do j = Index%js, Index%je
      Index%jpf(j) = 0   ! jpf non-zero for filtered rows
      if (j <= Grid%jsg+nlpf-1) Index%jpf(j) = -1
      if (j >= Grid%jeg-nlpf+1) Index%jpf(j) = +1
      Index%cph(j) = cos( Grid%tph(j) )
   enddo

 ! set trig constants for response function
   allocate ( Index%sklm(lenc) )
   do k = 2, lenc
      Index%sklm(k) = sin( 0.5*float(k-1)*dlm )
   enddo

 ! set domain types
   call pf_domain_init ( nlpf, nlev,  &
            Grid%Domain, Index%Local%Domain, Index%Zonal%Domain )
   call mpp_get_data_domain ( Index%Local%Domain,                &
                              Index%Local%isd, Index%Local%ied, &
                              Index%Local%jsd, Index%Local%jed  )
   call mpp_get_data_domain ( Index%Zonal%Domain,                &
                              Index%Zonal%isd, Index%Zonal%ied, &
                              Index%Zonal%jsd, Index%Zonal%jed  )
   call pf_domain_init ( nlpf, nlev*2,  &
            Grid%Domain, Index%Local2%Domain, Index%Zonal2%Domain )
   call mpp_get_data_domain ( Index%Local2%Domain,                 &
                              Index%Local2%isd, Index%Local2%ied, &
                              Index%Local2%jsd, Index%Local2%jed  )
   call mpp_get_data_domain ( Index%Zonal2%Domain,                 &
                              Index%Zonal2%isd, Index%Zonal2%ied, &
                              Index%Zonal2%jsd, Index%Zonal2%jed  )

 end subroutine set_index_type

!=======================================================================

subroutine pf_domain_init ( nlpf, nlev, Dom, Dom1, Dom2 )
integer       , intent(in)  :: nlpf, nlev 
type(domain2d), intent(in)  :: Dom  
type(domain2d), intent(inout) :: Dom1, Dom2

type(domain1d) :: Dx, Dy
integer, allocatable, dimension(:) :: xext, yext, yext1, ybeg, yend
logical, allocatable, dimension(:,:) :: mask 
integer :: layout(2), nlon, nlat, npes, ndiv, lin, rem, n
integer :: isg, ieg, jsg, jeg, rows
logical :: mask_rows

  npes = mpp_npes()

! get x/y extents
  call mpp_get_domain_components ( Dom, Dx, Dy )
  call mpp_get_layout ( Dom, layout )
  allocate ( xext(layout(1)), yext(layout(2)) )
  allocate ( ybeg(layout(2)), yend(layout(2)) )
  call mpp_get_compute_domains   ( Dx, size=xext )
  call mpp_get_compute_domains   ( Dy, size=yext )
  call mpp_get_compute_domains   ( Dy, begin=ybeg, end=yend )
  call mpp_get_global_domain     ( Dom, isg, ieg, jsg, jeg, &
                                   xsize=nlon, ysize=nlat   )
  ! compute number of filtered rows per y-axis processor
  mask_rows = .false.
  do n = 1, layout(2)
     rows = max(0, min(jsg+nlpf-1,yend(n))-max(jsg,ybeg(n))+1) + &
            max(0, min(jeg,yend(n))-max(jeg-nlpf+1,ybeg(n))+1)
     if (rows == 0) then
        yext(n) = 1  ! one row minimum
        mask_rows = .true.
     else
        yext(n) = rows*nlev
     endif
  enddo
  xext(layout(1)) = xext(layout(1)) + 1   ! extra column will hold cos(lat)
 !if (mpp_pe()==mpp_root_pe()) then
 ! print *, 'layout = ', layout
 ! print *, 'xext = ', xext
 ! print *, 'yext = ', yext
 !endif

! define new domain for polar filter
  call mpp_define_domains ( (/1,nlon+1,1,sum(yext)/), layout, Dom1, &
                           xflags = CYCLIC_GLOBAL_DOMAIN,           &
                           xextent = xext, yextent = yext )
                          !xextent = xext, yextent = yext, name = 'global' )
!------------------------------------------------------
! define zonal domain
  ndiv = npes
  if (mask_rows) ndiv = npes+1  ! setup extra division for masking rows
  layout = (/ 1, ndiv /)        ! one-dimensional decomposition
  allocate ( yext1(ndiv), mask(layout(1),layout(2)) )
  yext1 = 0
  mask  = .false.

! determine number of rows per PE in each hemisphere
! this can only be done when running on multiple PEs
  if (npes > 1) then
     ! sh
       lin = max(nlpf*nlev/((npes+1)/2), 1)
       rem = max(nlpf*nlev - lin*((npes+1)/2), 0)
       do n = 1, (npes+1)/2
         yext1(n) = lin
         mask (1,n) = .true.
       enddo
       do n = 1, rem
         yext1(n+1) = yext1(n+1) + 1
       enddo
     ! nh
       lin = max(nlpf*nlev/(npes/2), 1)
       rem = max(nlpf*nlev - lin*(npes/2), 0)
       do n = ndiv, ndiv-(npes/2)+1, -1
         yext1(n) = lin
         mask (1,n) = .true.
       enddo
       do n = ndiv, ndiv-rem+1, -1
         yext1(n) = yext1(n) + 1
       enddo
     ! non filtered row
       if (mask_rows) then
          lin = sum(yext1)
          yext1((npes+1)/2+1) = sum(yext) - lin
       endif
  else
! special case for NPES=1
! just copy original decomposition
     yext1 = yext
     mask  = .true.
  endif

 !if (mpp_pe()==mpp_root_pe()) then
 ! do n = npes, 0, -1
 !  print *, 'yext1, mask = ', yext1(n+1), mask(1,n+1)
 ! enddo
 !  print *, 'sum(yext1) = ', sum(yext1)
 !endif

  call mpp_define_domains ( (/1,nlon+1,1,sum(yext1)/), layout, Dom2, &
                           xflags = CYCLIC_GLOBAL_DOMAIN,            &
                           yextent = yext1, maskmap = mask )
                          !yextent = yext1, maskmap = mask, name = 'filter' )

  deallocate ( xext, yext, yext1, mask )

 end subroutine pf_domain_init

!#######################################################################

end module bgrid_polar_filter_mod


module topo_drag_mod

!=======================================================================
! TOPOGRAPHIC DRAG CLOSURE -- Garner (2005)
!=======================================================================

!-----------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!-----------------------------------------------------------------------

use          mpp_mod, only: input_nml_file
use          fms_mod, only: file_exist, open_namelist_file,            &
                            close_file, error_mesg, FATAL, NOTE,       &
                            mpp_pe, mpp_root_pe, stdout, stdlog,       &
                            check_nml_error, write_version_number
use       fms_io_mod, only: read_data, field_size
use       fms_io_mod, only: register_restart_field, restart_file_type
use       fms_io_mod, only: save_restart, restore_state
use    constants_mod, only: Grav, Cp_Air, Rdgas, Pi, Radian
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp, horiz_interp_del

implicit none

private

character(len=128) :: version = '$Id: topo_drag.F90,v 19.0 2012/01/06 20:27:27 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

logical :: module_is_initialized = .false.

! horizontal array size

integer :: nlon, nlat

! arrays defined by topo_drag_init:

real, allocatable, dimension(:,:) :: t11, t21, t12, t22    ! drag tensor
real, allocatable, dimension(:,:) :: hmin, hmax

! parameters:

real, parameter :: u0=1.0       ! arbitrary velocity scale for diagnostics
real, parameter :: xl=80.0e3    ! arbitrary horiz length scale for diagnostics
real, parameter :: ro=1.2       ! arbitrary density scale for diagnostics
real, parameter :: lapse=Grav/Cp_Air ! adiabatic lapse rate
real, parameter :: tiny=1.0e-8

real, parameter :: resolution=30.0 ! # of points per degree in topo datasets

integer, parameter :: ipts=360*resolution
integer, parameter :: jpts=180*resolution

!--- for netcdf restart
type(restart_file_type), save :: Top_restart

! parameters in namelist (topo_drag_nml):

real :: &
   frcrit=0.7   &      ! critical value of Froude number for nonlinear flow
  ,anonlin=7.0  &      ! amplitude of nonpropagating drag
  ,gamma=0.4    &      ! exponent in aspect ratio power law
  ,epsi=0.0     &      ! exponent in distribution power law
  ,beta=0.5     &      ! bluntness of topographic features
  ,h_frac=0.0   &      ! ratio of min to max subgrid mountain height
  ,zref_fac=1.0 &      ! adjusts level separating breaking/laminar flow
  ,no_drag_frac=0.05 & ! fraction of lower atmosphere with no breaking
  ,tboost=1.0 &        ! surface T boost to improve PBL height estimate
  ,pcut=0.0            ! high-level cutoff pressure for momentum forcing
logical :: &
   do_conserve_energy=.true. &! conserve total energy?
  ,keep_residual_flux=.true.  ! redistribute residual pseudomomentum?

NAMELIST /topo_drag_nml/                                               &
  frcrit, anonlin, beta, gamma, epsi,                                  &
  h_frac, zref_fac, no_drag_frac, tboost, pcut,                        &
  do_conserve_energy, keep_residual_flux

public topo_drag, topo_drag_init, topo_drag_end
public topo_drag_restart

contains

!#######################################################################

subroutine topo_drag                                                   &
                                     ( is, js, delt, uwnd, vwnd, atmp, &
                                           pfull, phalf, zfull, zhalf, &
                          dtaux, dtauy, dtemp, taux, tauy, taus, kbot )

integer, intent(in) :: is, js
real,    intent(in) :: delt
integer, intent(in), optional, dimension(:,:) :: kbot

! INPUT
! -----

! UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
! VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
! ATMP     Temperature at full levels (IDIM x JDIM x KDIM)
! PFULL    Pressure at full levels (IDIM x JDIM x KDIM)
! PHALF    Pressure at half levels (IDIM x JDIM x KDIM+1)
! ZFULL    Height at full levels (IDIM x JDIM x KDIM)
! ZHALF    Height at half levels (IDIM x JDIM x KDIM+1)

real, intent(in), dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:,:) :: pfull, phalf, zfull, zhalf

! OUTPUT
! ------

! DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
! DTEMP        Tendency of the temperature in K/s (IDIM x JDIM x KDIM)
! TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
! TAUS         clipped saturation momentum flux (IDIM x JDIM x KDIM) for diagnostics

real, intent(out), dimension(:,:)   :: taux, tauy
real, intent(out), dimension(:,:,:) :: dtaux, dtauy, dtemp, taus

integer, dimension(size(zfull,1),size(zfull,2)) :: ktop, kcut
real,    dimension(size(zhalf,1),size(zhalf,2),size(zhalf,3)) :: tausat

! work arrays

real, dimension(size(zfull,1),size(zfull,2)) :: taub, taul, taup, taun
real, dimension(size(zfull,1),size(zfull,2)) :: frulo, fruhi, frunl, rnorm

integer :: idim
integer :: jdim
integer :: k, kdim

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! estimate height of pbl

  call get_pbl ( atmp, zfull, pfull, phalf, ktop, kcut )

! calculate base flux

  call base_flux (                                                     &
                                             is, js, uwnd, vwnd, atmp, &
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                            zfull, zhalf, pfull, ktop )

! calculate saturation flux profile

  call satur_flux (                                                    &
                                                     uwnd, vwnd, atmp, &
                                       taux, tauy, taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, ktop, kcut )

! calculate momentum tendency

  call topo_drag_tend (                                                &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
                dtaux, dtauy, dtemp, zfull, zhalf, pfull, phalf, ktop )

! put saturation flux profile into 'taus' for diagnostics

  do k=1,kdim
     taus(:,:,k) = 0.5*rnorm(:,:)*(tausat(:,:,k) + tausat(:,:,k+1))
  enddo

! put total drag into 'taux,tauy' for diagnostics

  taup = taup - tausat(:,:,1)
  taub = (taup + taun)/taul
  taux = taux*taub
  tauy = tauy*taub

end subroutine topo_drag

!=======================================================================
                                  
subroutine base_flux (                                                 &
                                             is, js, uwnd, vwnd, atmp, &
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                            zfull, zhalf, pfull, ktop )

integer, intent(in) :: is, js
real, intent(in),  dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension(:,:,:) :: zfull, zhalf, pfull
real, intent(out), dimension(:,:)   :: taux, tauy
real, intent(out), dimension(:,:,:) :: dtaux, dtauy
real, intent(out), dimension(:,:)   :: taub, taul, taup, taun
real, intent(out), dimension(:,:)   :: frulo, fruhi, frunl, rnorm
integer, intent(in), dimension(:,:) :: ktop

integer :: i, idim, id
integer :: j, jdim, jd
integer :: k, kdim

real :: usat, bfreq2, bfreq, dphdz, vtau
real :: dzhalf, density
real :: frmin, frmax, frumin, frumax, fruclp, fru0, frusat
real :: rnormal, gterm

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! compute base flux

  do j=1,jdim
     jd = js+j-1
     do i=1,idim
        id = is+i-1
        k = ktop(i,j)

        dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
        density = (pfull(i,j,k) - pfull(i,j,k-1))/(Grav*dzhalf)
        bfreq2 = Grav*((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/  & 
                  (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
        bfreq = sqrt(max(tiny, bfreq2))

        taux(i,j) = (uwnd(i,j,k)*t11(id,jd) + vwnd(i,j,k)*t21(id,jd))  &
                                                         *bfreq*density
        tauy(i,j) = (uwnd(i,j,k)*t12(id,jd) + vwnd(i,j,k)*t22(id,jd))  &
                                                         *bfreq*density
        taub(i,j) = max(tiny, sqrt(taux(i,j)**2 + tauy(i,j)**2))

!       min/max Froude numbers based on surface flow

        vtau = max(tiny, -(uwnd(i,j,k)*taux(i,j)                       &
                         + vwnd(i,j,k)*tauy(i,j))/taub(i,j))
        frmax = hmax(id,jd)*bfreq / vtau
        frmin = hmin(id,jd)*bfreq / vtau

!       linear momentum flux associated with min/max Froude numbers

        dphdz = bfreq / vtau
        usat = density/ro * vtau / sqrt(dphdz*xl)
        frusat = frcrit*usat

        frumin = frmin*usat         ! linear momentum flux
        frumax = frmax*usat
        frumax = max(frumax,frumin + tiny)
        fruclp = min(frumax,max(frumin,frusat))
        fru0 = (u0/vtau)*usat

!       total drag in linear limit

        rnormal = (frumax**(2.0*gamma - epsi)                          &
                 - frumin**(2.0*gamma - epsi))/(2.0*gamma - epsi) 
        rnormal = fru0**gamma * ro/rnormal  

        gterm = (frumax**(gamma - epsi - beta)                         &
               - fruclp**(gamma - epsi - beta))/(gamma - epsi - beta)
        gterm = gterm*frusat**beta

        taul(i,j) =                                                    &
                    (frumax**(2.0 + gamma - epsi)                      &
                   - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi)

!       separate propagating and nonpropagating parts of total drag

        taup(i,j) =                                                    &
                 ( (fruclp**(2.0 + gamma - epsi)                       &       
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                                                    + frusat**2*gterm )

        taun(i,j) = anonlin*usat/(1.0 + beta)*                         &
                 ( (frumax**(1.0 + gamma - epsi)                       &
                  - fruclp**(1.0 + gamma - epsi))/(1.0 + gamma - epsi) &
                                                       - frusat*gterm )

        fruhi(i,j) = frumax
        frulo(i,j) = frumin
        frunl(i,j) = frusat
        rnorm(i,j) = rnormal

     enddo
  enddo

! wind component opposite the drag at full levels (stored as 'dtaux')

  do k=1,kdim
     do j=1,jdim
        do i=1,idim
           dtaux(i,j,k) =                                              &
             -(uwnd(i,j,k)*taux(i,j) + vwnd(i,j,k)*tauy(i,j))/taub(i,j)
        enddo
     enddo
  enddo

end subroutine base_flux

!=======================================================================

subroutine satur_flux (                                                &
                                                     uwnd, vwnd, atmp, &
                                       taux, tauy, taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, ktop, kcut )

real, intent(in),  dimension (:,:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension (:,:,:) :: dtaux, dtauy
real, intent(in),  dimension (:,:,:) :: zfull, pfull, phalf
real, intent(in),  dimension (:,:)   :: taux, tauy, taup
real, intent(out), dimension (:,:)   :: taub
real, intent(out), dimension (:,:,:) :: tausat
real, intent(in),  dimension (:,:)   :: frulo, fruhi, frunl
integer, intent(in), dimension (:,:) :: ktop, kcut

real, dimension(size(zfull,1),size(zfull,2)) :: usat

real :: dzhalf, gterm, gterm0, density
real :: bfreq2, bfreq, vtau, dphdz
real :: frumin, frumax, fruclp, frusat, frusat0, fruclp0

integer :: i, idim
integer :: j, jdim
integer :: k, kdim, k1

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! get vertical profile of propagating part of momentum flux

  usat = frunl/frcrit

  do k=kdim,2,-1
     do j=1,jdim
        do i=1,idim

!          buoyancy frequency, velocity and density at half levels

           dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
           density = (pfull(i,j,k) - pfull(i,j,k-1))/(Grav*dzhalf)
           bfreq2 = Grav*                                              &
                       ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ & 
                   (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           
           vtau = max(tiny, 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k)))

!          min/max and critical momentum flux values at half levels

           dphdz = bfreq / vtau
           usat(i,j) = min(usat(i,j),                                  &
                      sqrt(density/ro * vtau/sqrt(dphdz*xl)))
           frusat = frcrit*usat(i,j)

           frumin = frulo(i,j)
           frumax = fruhi(i,j)
           fruclp = min(frumax,max(frumin,frusat))
           frusat0 = frunl(i,j)
           fruclp0 = min(frumax,max(frumin,frusat0))

!          propagating part of momentum flux (from WKB or EP)

           gterm0 = (frumax**(gamma - epsi - beta)                     &
                - fruclp0**(gamma - epsi - beta))/(gamma - epsi - beta)
           gterm = (fruclp0**(gamma - epsi)                            &
                               - fruclp**(gamma - epsi))/(gamma - epsi)

           tausat(i,j,k) =                                             &
                 ( (fruclp**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                         + frusat**2.0*(gterm0*frusat0**beta + gterm) )
        enddo
     enddo
  enddo

! make propagating flux constant with height in zero-drag top layer

  tausat(:,:,1) = tausat(:,:,2)

  k1 = maxval(kcut)
  do k=3,k1
     where (k <= kcut)
        tausat(:,:,k) = tausat(:,:,1)
     endwhere
  enddo

! make propagating flux constant with height in zero-drag surface layer

  k1 = minval(ktop)
  do k=kdim+1,k1+1,-1
     where (k > ktop)
        tausat(:,:,k) = taup
     endwhere
  enddo

! redistribute residual forcing

  if ( keep_residual_flux ) then
     taub(:,:) = tausat(:,:,1)/(phalf(:,:,kdim+1) - phalf(:,:,1))
     do k=1,kdim
        tausat(:,:,k) = tausat(:,:,k)                                  &
                        - taub(:,:)*(phalf(:,:,kdim+1) - phalf(:,:,k))
     enddo
  endif

endsubroutine satur_flux

!=======================================================================

subroutine topo_drag_tend (                                            &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
                dtaux, dtauy, dtemp, zfull, zhalf, pfull, phalf, ktop )

real, intent(in) :: delt
real, intent(in), dimension(:,:,:)   :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:,:)   :: zfull, zhalf, pfull, phalf
real, intent(in), dimension(:,:)     :: taux, tauy, taul, taun
real, intent(in), dimension(:,:,:)   :: tausat
real, intent(inout),dimension(:,:,:) :: dtaux, dtauy, dtemp
integer, intent(in), dimension (:,:) :: ktop

real, parameter :: bfmin=0.7e-2, bfmax=1.7e-2  ! min/max buoyancy freq [1/s]
real, parameter :: vvmin=1.0                   ! minimum surface wind [m/s]

integer,dimension(size(zfull,1),size(zfull,2)) :: kref
real :: dzhalf, zlast, rscale, phase, bfreq, bfreq2, vtau
real :: gfac, gfac1, dp, weight, wtsum

integer :: i, idim
integer :: j, jdim
integer :: k, kdim, kr

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! find reference level for non-propagating drag (z ~ pi U/N)

  do j=1,jdim
     do i=1,idim
        k = ktop(i,j)
        phase = 0.0
        zlast = zhalf(i,j,k)
        do while (phase <= Pi*zref_fac .and. k > 1)
           k = k-1
           vtau = 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k))
           dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
           bfreq2 = Grav*                                              &
                       ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ &
                   (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           rscale = max(bfmin, min(bfmax, bfreq))/max(vvmin, vtau)
           dzhalf = zfull(i,j,k-1) - zlast
           phase = phase + dzhalf*rscale
           zlast = zfull(i,j,k-1)
        enddo
        kref(i,j) = k
     enddo
  enddo

! CALCULATE DECELERATION DUE TO PROPAGATING DRAG (~-rho^-1 dtau/dz)

  do k=1,kdim
     do j=1,jdim
        do i=1,idim
          dp = phalf(i,j,k+1) - phalf(i,j,k)
          gfac = tausat(i,j,k+1) - tausat(i,j,k)
          gfac1 = gfac*Grav/(dp*taul(i,j))
          dtaux(i,j,k) = gfac1*taux(i,j)
          dtauy(i,j,k) = gfac1*tauy(i,j)
        enddo
     enddo
  enddo

! CALCULATE DECELERATION DUE TO NON-PROPAGATING DRAG

  do j=1,jdim
     do i=1,idim
        kr = kref(i,j)
        dp = phalf(i,j,kdim+1) - phalf(i,j,kr)
        gfac = taun(i,j)*Grav/(dp*taul(i,j))
        wtsum = 0.0
        do k=kr,kdim
           weight = pfull(i,j,k) - phalf(i,j,kr)
           wtsum = wtsum + weight
        enddo
        do k=kr,kdim
           weight = pfull(i,j,k) - phalf(i,j,kr)
           gfac1 = gfac*weight/wtsum
           dtaux(i,j,k) = dtaux(i,j,k) + gfac1*taux(i,j)
           dtauy(i,j,k) = dtauy(i,j,k) + gfac1*tauy(i,j)
        enddo
     enddo
  enddo

! CALCULATE HEATING TO CONSERVE TOTAL ENERGY

  if (do_conserve_energy) then
     dtemp = -((uwnd + 0.5*delt*dtaux)*dtaux                           &
             + (vwnd + 0.5*delt*dtauy)*dtauy)/Cp_Air
  else
     dtemp = 0.0
  endif

end subroutine topo_drag_tend

!=======================================================================

subroutine get_pbl ( atmp, zfull, pfull, phalf, ktop, kcut )

integer, intent(out), dimension(:,:) :: ktop, kcut
real, intent(in), dimension(:,:,:)   :: atmp
real, intent(in), dimension(:,:,:)   :: zfull, pfull, phalf

real, dimension(size(pfull,1),size(pfull,2)) :: pbot, tbot, zbot

integer :: i, idim
integer :: j, jdim
integer :: k, kdim

  idim = size(atmp,1)
  jdim = size(atmp,2)
  kdim = size(atmp,3)

  do j=1,jdim
     do i=1,idim
        pbot(i,j) = (1.0 - no_drag_frac)*phalf(i,j,kdim+1)
        tbot(i,j) = atmp(i,j,kdim) + tboost
        zbot(i,j) = zfull(i,j,kdim)
     enddo
  enddo

! find highest model level in no-drag surface layer

  ktop = kdim-1

  do k=kdim-2,2,-1
     where ( pfull(:,:,k) <= pbot(:,:) .and.                           &
           tbot(:,:) - atmp(:,:,k) > lapse*(zfull(:,:,k) - zbot(:,:)) )
        ktop = k
     endwhere
  enddo

! find lowest model level in no-drag top layer

  kcut = 1

  do k=2,kdim
     where ( pfull(:,:,k) <= pcut )
        kcut = k
     endwhere
  enddo

end subroutine get_pbl

!=======================================================================

subroutine topo_drag_init (lonb, latb)

real, intent(in), dimension(:,:) :: lonb, latb

character(len=128) :: msg
character(len=64)  :: restart_file='topo_drag.res.nc'
character(len=64)  :: topography_file='INPUT/postopog_2min_hp150km.nc'
character(len=64)  :: dragtensor_file='INPUT/dragelements_2min_hp150km.nc'
character(len=3)   :: tensornames(4) = (/ 't11', 't21', 't12', 't22' /)

real, parameter :: bfscale=1.0e-2      ! buoyancy frequency scale [1/s]

real, allocatable, dimension(:)   :: xdat, ydat
real, allocatable, dimension(:,:) :: zdat, zout
type (horiz_interp_type) :: Interp
real :: exponent

integer :: n
integer :: io, ierr, unit_nml, logunit
integer :: i, j
integer :: siz(4)
integer :: id_restart

  if (module_is_initialized) return

  nlon = size(lonb,1)-1
  nlat = size(latb,2)-1

! read namelist

#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=topo_drag_nml, iostat=io)
   ierr = check_nml_error(io,'topo_drag_nml')
#else   
  if( file_exist( 'input.nml' ) ) then
     unit_nml = open_namelist_file ( )
     ierr = 1
     do while ( ierr /= 0 )
        read( unit_nml, nml = topo_drag_nml, iostat = io, end = 10 )
        ierr = check_nml_error (io, 'topo_drag_nml')
     end do
 10  call close_file ( unit_nml )
  endif
#endif

! write version number and namelist to logfile

  call write_version_number(version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe())                                       &
                                    write (logunit, nml=topo_drag_nml)

  allocate (t11(nlon,nlat))
  allocate (t21(nlon,nlat))
  allocate (t12(nlon,nlat))
  allocate (t22(nlon,nlat))
  allocate (hmin(nlon,nlat))
  allocate (hmax(nlon,nlat))

  if (gamma == beta + epsi) gamma = gamma + tiny

! read restart file

  id_restart = register_restart_field(Top_restart, restart_file, 't11', t11)
  id_restart = register_restart_field(Top_restart, restart_file, 't21', t21)
  id_restart = register_restart_field(Top_restart, restart_file, 't12', t12)
  id_restart = register_restart_field(Top_restart, restart_file, 't22', t22)
  id_restart = register_restart_field(Top_restart, restart_file, 'hmin', hmin)
  id_restart = register_restart_field(Top_restart, restart_file, 'hmax', hmax)
  restart_file = 'INPUT/'//trim(restart_file)
  if ( file_exist(restart_file) ) then

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading restart file: ",a40)' ) restart_file
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif
     call restore_state(Top_restart)

  else if (file_exist(topography_file) .and.                           &
           file_exist(dragtensor_file)) then

!    read and interpolate topography datasets

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading topography file: ",a)')              &
                                                        trim(topography_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     ! check for correct field size in topography
     call field_size (topography_file, 'hpos', siz)
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Field \"hpos\" in file '// &
                   trim(topography_file)//' has the wrong size', FATAL)
     endif
     
     allocate (xdat(ipts+1))
     allocate (ydat(jpts+1))
     allocate (zdat(ipts,jpts))
     allocate (zout(nlon,nlat))

     do i=1,ipts+1
        xdat(i) = (i-1)/resolution / Radian
     enddo
     do j=1,jpts+1
        ydat(j) = (-90.0 + (j-1)/resolution) / Radian
     enddo

     ! initialize horizontal interpolation
     ! Note: interp_method will be conservative for lat/lon grid
     !       and bilinear for all other grids
     call horiz_interp_init
     call horiz_interp_new ( Interp, xdat, ydat, lonb, latb )

     call read_data (topography_file, 'hpos', zdat, no_domain=.true.)

     exponent = 2.0 - gamma
     zdat = max(0.0, zdat)**exponent
     
     call horiz_interp ( Interp, zdat, zout )

     hmax = abs(zout)**(1.0/exponent) * sqrt(2.0/exponent)
     hmin = hmax*h_frac

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading drag tensor file: ",a)')             &
                                                        trim(dragtensor_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     ! check for correct field size in tensor file
     call field_size (dragtensor_file, tensornames(1), siz)
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Field \"'//tensornames(1)// &
         '\" in file '//trim(dragtensor_file)//' has the wrong size', FATAL)
     endif

     do n=1,4
        call read_data (dragtensor_file, tensornames(n), zdat, no_domain=.true.)
        call horiz_interp ( Interp, zdat, zout )
       !call horiz_interp ( zdat, xdat, ydat, lonb, latb, zout,        &
       !                                 interp_method='conservative' )
        if ( tensornames(n) == 't11' ) then
           t11 = zout/bfscale
        else if ( tensornames(n) == 't21' ) then
           t21 = zout/bfscale
        else if ( tensornames(n) == 't12' ) then
           t12 = zout/bfscale
        else if ( tensornames(n) == 't22' ) then
           t22 = zout/bfscale
        endif
     enddo

     deallocate (zdat, zout)
     call horiz_interp_del ( Interp )

  else

     call ERROR_MESG ('topo_drag_init',                                &
                'No sub-grid orography available for topo_drag', FATAL)

  endif

  module_is_initialized = .true.

end subroutine topo_drag_init

!=======================================================================

subroutine topo_drag_end

! writes static arrays to restart file

  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('topo_drag_mod', 'Writing netCDF formatted restart file: RESTART/topo_drag.res.nc', NOTE)
  endif

  call topo_drag_restart

  module_is_initialized = .false.

end subroutine topo_drag_end

!#######################################################################
! <SUBROUTINE NAME="topo_drag_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine topo_drag_restart(timestamp)
   character(len=*), intent(in), optional :: timestamp

   call save_restart(Top_restart, timestamp)

end subroutine topo_drag_restart
! </SUBROUTINE>

!#######################################################################

endmodule topo_drag_mod

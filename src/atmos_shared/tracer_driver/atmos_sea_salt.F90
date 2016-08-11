module atmos_sea_salt_mod
! <DESCRIPTION>
!   This module evaluates the change of mass mixing ratio of sea salt
!   particles due to their emission at the ocean surface, and the removal by
!   gravitational settling. The sea salt particles are transported as dry
!   particles. Therefore, some conversion from wet to dry and vice-et-versa
!   are consdered.
!   The size distribution of sea salt ranges from 0.1 to 10 um (dry radius)
!   and is divided into 5 bins. For each bin, the volume size distribution
!   dV/dlnr is considered constant.
! </DESCRIPTION>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Paul Ginoux
! </CONTACT>
! 
!
!-----------------------------------------------------------------------
use mpp_mod, only: input_nml_file 
use              fms_mod, only : file_exist, &
                                 write_version_number, &
                                 close_file,              &
                                 mpp_pe, &
                                 mpp_root_pe, &
                                 open_namelist_file, file_exist,    &
                                 check_nml_error, error_mesg,  &
                                 stdlog
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_diag_field
use   tracer_manager_mod, only : get_tracer_index
use    field_manager_mod, only : MODEL_ATMOS
use        constants_mod, only : PI, GRAV,RDGAS, DENS_H2O, PSTD_MKS, WTMAIR
implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_sea_salt_sourcesink, atmos_sea_salt_init, atmos_sea_salt_end

!---- version number -----
character(len=128) :: version = '$Id: atmos_sea_salt.F90,v 20.0 2013/12/13 23:23:59 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


!-----------------------------------------------------------------------
!----------- namelist -------------------
character(len=80) :: scheme = " "
real, save :: coef1
real, save :: coef2
!---------------------------------------------------------------------
real :: coef_emis1=-999.
real :: coef_emis2=-999.
real :: critical_sea_fraction = 0.0   ! sea-salt aerosol production
                                      ! occurs in grid cells with
                                      ! ocn_flx_fraction .gt. this value


namelist /ssalt_nml/  scheme, coef_emis1, coef_emis2, &
                      critical_sea_fraction

!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=7), parameter :: mod_name = 'tracers'

integer :: nseasalt=0  ! tracer number for sea_salt
!--- identification numbers for  diagnostic fields and axes ----

integer :: id_SS_emis(5), id_SS_setl(5)
integer :: id_ocn_flx_fraction

logical :: module_is_initialized=.FALSE.
logical :: used
integer :: logunit

      real, parameter :: small_value = 1.e-20
      real, parameter :: mtv  = 1.    ! factor connversion for mixing ratio
      integer, parameter :: nrh= 65   ! number of RH in look-up table
      integer, parameter :: nr = 10   ! number of integration points 
                                      ! The difference with nr=100 & nr=5 < 1e-3
      real, dimension(nrh) :: rho_table, growth_table
!! Sea salt hygroscopic growth factor from 35 to 99% RH
!! We start at the deliquescence point of sea-salt for RH=37%
!! Any lower RH doesn't affect dry properties
!! Reference: Tang et al., JGR, v102(D19), 23,269-23,275, 1997. 
      data growth_table/1.000, 1.000, 1.396, &
       1.413, 1.428, 1.441, 1.454, 1.466, 1.478, 1.490, 1.501, 1.512, &
       1.523, 1.534, 1.545, 1.555, 1.566, 1.577, 1.588, 1.599, 1.610, &
       1.621, 1.632, 1.644, 1.655, 1.667, 1.679, 1.692, 1.704, 1.717, &
       1.730, 1.743, 1.757, 1.771, 1.786, 1.801, 1.816, 1.832, 1.849, &
       1.866, 1.884, 1.903, 1.923, 1.944, 1.966, 1.990, 2.014, 2.041, &
       2.069, 2.100, 2.134, 2.170, 2.210, 2.255, 2.306, 2.363, 2.430, &
       2.509, 2.605, 2.723, 2.880, 3.087, 3.402, 3.919, 5.048/
!! Seal salt density for 65 RH values from 35% to 99% [g/cm3]
      data rho_table/2.160, 2.160, 1.490, &
       1.475, 1.463, 1.452, 1.441, 1.432, 1.422, 1.414, 1.406, 1.398, &
       1.390, 1.382, 1.375, 1.368, 1.361, 1.354, 1.347, 1.341, 1.334, &
       1.328, 1.322, 1.315, 1.309, 1.303, 1.297, 1.291, 1.285, 1.279, &
       1.273, 1.267, 1.261, 1.255, 1.249, 1.243, 1.237, 1.231, 1.225, &
       1.219, 1.213, 1.207, 1.201, 1.195, 1.189, 1.183, 1.176, 1.170, &
       1.163, 1.156, 1.150, 1.142, 1.135, 1.128, 1.120, 1.112, 1.103, &
       1.094, 1.084, 1.074, 1.063, 1.051, 1.038, 1.025, 1.011/

real :: betha

!-----------------------------------------------------------------------

contains

!#######################################################################
!<SUBROUTINE NAME="atmos_sea_salt_sourcesink">
!<OVERVIEW>
! The routine that calculate the emission and settling of sea_salt.
!</OVERVIEW>
 subroutine atmos_sea_salt_sourcesink (i_SS,ra,rb,ssaltref,ssaltden, &
                                       lon, lat, ocn_flx_fraction, pwt, &
                                       zhalf, pfull, w10m, t, rh, &
                                       seasalt, seasalt_dt, dt, SS_setl, &
                                       SS_emis, Time, Time_next, &
                                       is,ie,js,je, kbot)

!-----------------------------------------------------------------------
   integer, intent(in)                 :: i_SS
   real, intent(in)                    :: ra
   real, intent(in)                    :: rb
   real, intent(in)                    :: dt
   real, intent(in)                    :: ssaltref
   real, intent(in)                    :: ssaltden
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: ocn_flx_fraction
   real, intent(in),  dimension(:,:)   :: w10m
   real, intent(out), dimension(:,:)   :: SS_setl, SS_emis
   real, intent(in),  dimension(:,:,:) :: pwt, seasalt
   real, intent(in),  dimension(:,:,:) :: zhalf, pfull, t, rh
   real, intent(out), dimension(:,:,:) :: seasalt_dt
   integer, intent(in)                 :: is, ie, js, je
   type(time_type), intent(in) :: Time, Time_next     
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------

real, dimension(size(seasalt,3)) :: SS_conc0, SS_conc1
integer  i, j, k, id, jd, kd, kb, ir, irh

!----------------------------------------------
!     Sea-Salt parameters
!----------------------------------------------
      real, dimension(size(seasalt,3)) :: setl

      integer :: istep, nstep
      real, parameter :: ptmb = 0.01     ! pascal to mb
      integer, parameter :: nstep_max = 5  !Maximum number of cyles for settling
      real :: step
      real :: rhb
      real :: rho_wet_salt, viscosity, free_path, C_c
      real :: rho_air, Bcoef
      real :: r, dr, rmid, rwet, seasalt_flux
      real :: a1, a2
      real, dimension(size(pfull,3))  :: vdep
!-----------------------------------------------------------------------

      id=size(seasalt,1); jd=size(seasalt,2); kd=size(seasalt,3)

      SS_emis(:,:)      = 0.0
      SS_setl(:,:)      = 0.0
      seasalt_dt(:,:,:) = 0.0

!-------------------------------------------------------------------------
! Smith et al. (1993) derived an expression for the sea-salt flux
! by assuming that particle size spectra measured at a height of 10 m
! on the ocast of the Outer Hevrides islands with prevailing winds from
! the ocean present a balance between production and loss. They approximate
! the flux by the sum of two lognormal distributions.
! 
! seasalt_flux is the flux of dry particles by bubble bursting .
! Monahan et al. (1986) established a formula for the flux of wet
! particles by bubble bursting. The radius needs to be converted
! into dry assuming an equilibrium relative humidity above water RH=80%
! Using Tang (1996) formula the sea-salt hygroscopic growth betha=2.009
! at RH=80 (cf. Fitzegrald and Hoppel, JGR, v103 D13, 16085-16102, 1998)
! Also, Monahan et al. (1986) formula is defined for radius in units of
! micrometers
! The flux of particles is converted into mass flux. The cube of the
! radius disappear by cancelation with r^3 in Monahan formula. However
! there is 1/betha3 remaining multiply by betha from the integrand, such
! that finally 1/betha^2 remained. 
!-------------------------------------------------------------------------

      if (scheme .ne. "Smith") then !  ie, Monahan et . (1986)
        r = ra* 1.e6
        dr= (rb - ra)/float(nr)* 1.e6
        seasalt_flux=0.
        do ir=1,nr
          rmid=r+dr*0.5   ! Dry radius
          r=r+dr
          Bcoef=(coef1-alog10(betha*rmid))/coef2
          seasalt_flux = seasalt_flux + &
             1.373*4./3.*pi*ssaltden/betha**2*1.e-18* &
             (1.+0.057*(betha*rmid)**1.05)*dr*      &
             10**(1.19*exp(-(Bcoef**2)))
        enddo
          do j=1,jd
            do i=1,id
              if (ocn_flx_fraction (i,j).gt.critical_sea_fraction) then
                SS_emis(i,j) = seasalt_flux*(ocn_flx_fraction (i,j))*   &
                                                           w10m(i,j)**3.41
              endif
            enddo
          enddo

      else  !  scheme .ne. smith
        if (mpp_pe() == mpp_root_pe())   &
          write (logunit,*) "Smith parameterization for sea-salt production"
          do j=1,jd
            do i=1,id
              if (ocn_flx_fraction (i,j).gt.critical_sea_fraction) then
!------------------------------------------------------------------
!    Surface emission of sea salt  (Smith et al. (1993))
!------------------------------------------------------------------
                seasalt_flux            = 0.0
                a1=exp(0.155*w10m(i,j)+5.595)
                a2=exp(2.2082*sqrt(w10m(i,j))-3.3986)
                r = ra* 1.e6
                dr= (rb - ra)/float(nr)* 1.e6
                do ir=1,nr
                  rmid=r+dr*0.5   ! Dry radius
                  r=r+dr
                  seasalt_flux            = seasalt_flux            + &
                     4.188e-18*rmid**3*ssaltden*betha*( &
                    + coef1*a1*exp(-3.1*(alog(betha*rmid/2.1))**2) &
                    + coef2*a2*exp(-3.3*(alog(betha*rmid/9.2))**2) )
                enddo
                  SS_emis(i,j) = seasalt_flux*ocn_flx_fraction (i,j)
              endif  
            enddo
          enddo
      endif

      if (present(kbot)) then
        do j=1,jd
          do i=1,id
            kb = kbot(i,j)
            seasalt_dt(i,j,kb)=amax1(0.,SS_emis(i,j)/pwt(i,j,kb)*mtv)
          enddo
        enddo
      else  
        do j=1,jd
          do i=1,id
            seasalt_dt(i,j,kd) =  amax1(0., SS_emis(i,j)/pwt(i,j,kd)*mtv)
          end do
        end do
      endif  

! Send the emission data to the diag_manager for output.
      if (id_SS_emis(i_SS) > 0 ) then
        used = send_data ( id_SS_emis(i_SS), SS_emis, Time_next, &
              is_in=is,js_in=js )
      endif
      if (id_ocn_flx_fraction > 0) then
        used = send_data ( id_ocn_flx_fraction, ocn_flx_fraction, Time_next, &
             is_in=is,js_in=js )
      endif

!------------------------------------------
!       Solve at the model TOP (layer plev-10)
!------------------------------------------
        do j=1,jd
          do i=1,id
            if (present(kbot)) then
              kb=kbot(i,j)
            else
              kb=kd
            endif
!
! Determine the maximum timestep to avoid particles settling more than 1 layer
!
            nstep=1
            do k=1,kb
              rhb=amin1(0.99,rh(i,j,k))
              rhb=amax1(0.001,rhb)
              irh=max0(1,int(rhb*100.-34.))
              rho_wet_salt=rho_table(irh)*1000. !Density of wet sea-salt [kg/m3]
              rwet=ssaltref*growth_table(irh) ! Radius of wet sea-salt [m]
              viscosity = 1.458E-6 * t(i,j,k)**1.5/(t(i,j,k)+110.4)
              free_path=6.6e-8*t(i,j,k)/293.15*(PSTD_MKS/pfull(i,j,k))
              C_c=1. + free_path/ssaltref* &            ! Slip correction [none]
                    (1.257+0.4*exp(-1.1*ssaltref/free_path))
              Vdep(k)=2./9.*C_c*GRAV*rho_wet_salt*rwet**2./viscosity
              step = (zhalf(i,j,k)-zhalf(i,j,k+1)) / vdep(k) / 2.
              nstep = max(nstep, int( dt/ step) )
!!! To avoid spending too much time on cycling the settling in case
!!! of very large particles falling through a tiny layer, impose
!!! maximum speed for the selected nstep_max. This is not physically
!!! correct, but as these particles are very large there will be removed
!!! fast enough to not change significantly their lifetime. The proper
!!! way would be to implement semi-lagrangian technique.
              if (nstep.gt.nstep_max) then
                nstep = nstep_max
                vdep(k)=(zhalf(i,j,k)-zhalf(i,j,k+1))*nstep / 2. /dt
              endif
            enddo
            step = dt / nstep

            SS_conc1(:) = seasalt(i,j,:) 
            do istep = 1, nstep
              SS_conc0(:) = SS_conc1(:) 
              do k=1,kb
                rho_air = pfull(i,j,k)/t(i,j,k)/RDGAS ! Air density [kg/m3]
                if (SS_conc0(k).gt.0.) then
!!!               settling flux [kg/m2/s]
                  setl(k)=SS_conc0(k)*rho_air/mtv*vdep(k) 
                else
                  setl(k)=0.
                endif
              enddo
              SS_setl(i,j)=SS_setl(i,j)+setl(kb)*step
              SS_conc1(1) = SS_conc0(1) - setl(1)/pwt(i,j,1)*mtv * step
              SS_conc1(2:kb)= SS_conc0(2:kb) &
              + ( setl(1:kb-1) - setl(2:kb) )/pwt(i,j,2:kb)*mtv * step
              where (SS_conc1 < 0 ) SS_conc1=0.0
            enddo
            seasalt_dt(i,j,:)=seasalt_dt(i,j,:)+ (SS_conc1(:)-seasalt(i,j,:))/dt
            SS_setl(i,j)=SS_setl(i,j)/dt
          enddo
        enddo 

! Send the settling data to the diag_manager for output.
      if (id_SS_setl(i_SS) > 0 ) then
        used = send_data ( id_SS_setl(i_SS), SS_setl, Time_next, &
              is_in=is,js_in=js )
      endif


!-----------------------------------------------------------------------

 end subroutine atmos_sea_salt_sourcesink
!</SUBROUTINE>


!#######################################################################
!<SUBROUTINE NAME="atmos_sea_salt_init">
!<OVERVIEW>
! The constructor routine for the sea_salt module.
!</OVERVIEW>
 subroutine atmos_sea_salt_init (lonb, latb, axes, Time, mask)
!-----------------------------------------------------------------------
real, intent(in),    dimension(:,:)                 :: lonb, latb
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
integer :: n, m
!
!-----------------------------------------------------------------------
!
      integer  unit,ierr,io
      character(len=1) :: numb(5)
      data numb/'1','2','3','4','5'/

      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=ssalt_nml, iostat=io)
        ierr = check_nml_error(io,'ssalt_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=ssalt_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'ssalt_nml')
        end do
10      call close_file (unit)
#endif
      endif
!--------- write version and namelist to standard log ------------
      call write_version_number(version, tagname)
      logunit=stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=ssalt_nml )
      if (scheme .eq. "Smith") then
        if (coef_emis1 .le. -990) then
          coef1 = 1.0
        else
          coef1 = coef_emis1
        endif
        if (coef_emis2 .le. -990) then
          coef2 = 1.0
        else
          coef2 = coef_emis2
        endif
      else
        if (coef_emis1 .le. -990) then
          coef1 = 0.38
        else
          coef1 = coef_emis1
        endif
        if (coef_emis2 .le. -990) then
          coef2 = 0.65
        else
          coef2 = coef_emis2
        endif
      endif

      betha = growth_table(46)  ! Growth factor at 80% RH

!----- set initial value of sea_salt ------------
       do m=1,5
       n = get_tracer_index(MODEL_ATMOS,'seasalt'//numb(m))
       if (n>0) then
         nseasalt=n
         if (nseasalt > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) 'Sea-salt',nseasalt
         if (nseasalt > 0 .and. mpp_pe() == mpp_root_pe()) write (logunit,30) 'Sea-salt',nseasalt
       endif
!

! Register a diagnostic field : emission of seasalt
       id_SS_emis(m) = register_diag_field ( mod_name,             &
                     'ssalt'//numb(m)//'_emis', axes(1:2),Time,          &
                     'ssalt'//numb(m)//'_emis', 'kg/m2/s',      &
                     missing_value=-999.  )

! Register a diagnostic field : settling of seasalt
       id_SS_setl(m) = register_diag_field ( mod_name,              &
                     'ssalt'//numb(m)//'_setl', axes(1:2),Time,        &
                     'ssalt'//numb(m)//'_setl', 'kg/m2/s',      &
                     missing_value=-999.  )
      enddo

      id_ocn_flx_fraction  = register_diag_field   &
                  ( mod_name,'salt_flux_area_frac', &
                   axes(1:2),Time,'fractional_area_with_salt_flux','-', &
                                                     missing_value=-999.)
      module_is_initialized = .TRUE.
 
  30        format (A,' was initialized as tracer number ',i2)

!-----------------------------------------------------------------------

 end subroutine atmos_sea_salt_init
!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_sea_salt_end">
!<OVERVIEW>
!  The destructor routine for the sea_salt module.
!</OVERVIEW>
 subroutine atmos_sea_salt_end
 
      module_is_initialized = .FALSE.

 end subroutine atmos_sea_salt_end
!</SUBROUTINE>


end module atmos_sea_salt_mod

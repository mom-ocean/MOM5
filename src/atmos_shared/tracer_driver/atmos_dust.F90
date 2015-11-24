module atmos_dust_mod
! <DESCRIPTION>
!   This module evaluates the change of mass mixing ratio for mineral dust
!   particles due to their emission from preferential sources, and the removal
!   by gravitational settling. The dust particles are transported as dry
!   particles. No hygroscopic growth is considered.
!   The size distribution of sea salt ranges from 0.1 to 10 um (dry radius)
!   and is divided into 5 bins. For each bin, the volume size distribution
!   dV/dlnr is considered constant.
! </DESCRIPTION>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
use              fms_mod, only : file_exist, &
                                 write_version_number, &
                                 mpp_pe, &
                                 mpp_root_pE, &
                                 close_file,           &
                                 open_namelist_file, file_exist,    &
                                 check_nml_error, error_mesg,  &
                                 stdlog
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_diag_field
use   tracer_manager_mod, only : get_tracer_index, &
                                 set_tracer_atts
use    field_manager_mod, only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                 dry_deposition
use interpolator_mod,    only:  interpolate_type, interpolator_init, &
                                obtain_interpolator_time_slices, &
                                unset_interpolator_time_flag, &
                                interpolator, interpolator_end, &
                                CONSTANT, INTERP_WEIGHTED_P
use     constants_mod, only : PI, GRAV, RDGAS, DENS_H2O, PSTD_MKS, WTMAIR


implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_dust_sourcesink, atmos_dust_init, atmos_dust_end, &
        atmos_dust_time_vary, atmos_dust_endts

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: ndust=0  ! tracer number for dust
!--- identification numbers for  diagnostic fields and axes ----

integer :: id_dust_emis(5), id_dust_setl(5)
integer :: id_dust_source

!--- Arrays to help calculate tracer sources/sinks ---
type(interpolate_type),save         ::  dust_source_interp


logical :: module_is_initialized=.FALSE.
logical :: used

real, save :: u_ts
real, save :: ch

!---------------------------------------------------------------------
!-------- namelist  ---------
character(len=32)  :: dust_source_filename = 'dust_source_1x1.nc'
character(len=32)  :: dust_source_name(1) = 'source'
real :: uthresh=-999.
real :: coef_emis =-999.

namelist /dust_nml/  dust_source_filename, dust_source_name, uthresh, coef_emis

!---- version number -----
character(len=128) :: version = '$Id: atmos_dust.F90,v 20.0 2013/12/13 23:23:54 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################
!<SUBROUTINE NAME="atmos_dust_sourcesink">
!<OVERVIEW>
! The routine that calculate the sources and sinks of dust.
!</OVERVIEW>
 subroutine atmos_dust_sourcesink (i_DU,ra,rb,dustref,dustden, &
       lon, lat, frac_land, pwt, &
       zhalf, pfull, w10m, t, rh, &
       dust, dust_dt, dust_emis, dust_setl, Time, Time_next, is,ie,js,je,kbot)

!-----------------------------------------------------------------------
   integer, intent(in)                 :: i_DU
   real, intent(in)                    :: ra
   real, intent(in)                    :: rb
   real, intent(in)                    :: dustref
   real, intent(in)                    :: dustden
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: frac_land
   real, intent(in),  dimension(:,:)   :: w10m
   real, intent(out),  dimension(:,:)  :: dust_setl, dust_emis
   real, intent(in),  dimension(:,:,:) :: pwt, dust
   real, intent(in),  dimension(:,:,:) :: zhalf, pfull, t, rh
   real, intent(out), dimension(:,:,:) :: dust_dt
   type(time_type), intent(in) :: Time, Time_next     
   integer, intent(in),  dimension(:,:), optional :: kbot
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
integer  i, j, k, id, jd, kd, kb
!----------------------------------------------
!     Dust parameters
!----------------------------------------------
      real, dimension(5) ::   frac_s

      real, dimension(size(dust,3)) :: setl
      real, dimension(size(dust,1),size(dust,2)) :: u_ts_2d, source

      real, parameter :: small_value = 1.e-20
      real, parameter :: mtcm = 100.            ! meter to cm
      real, parameter :: mtv  = 1. ! factor conversion for mixing ratio of dust
      real, parameter :: ptmb = 0.01     ! pascal to mb

      real :: rhb, rcm
      real :: ratio_r, rho_wet_dust, viscosity, free_path, C_c, vdep
      real :: rho_air
      real :: rwet
!-----------------------------------
!    SET-Up  DATA
!-----------------------------------

!yim: per pag 2/1/08
!     data frac_s/0.1,0.225,0.225,0.225,0.225/
      data frac_s/0.05,0.1125,0.225,0.225,0.225/

!-----------------------------------------------------------------------

      id=size(dust,1); jd=size(dust,2); kd=size(dust,3)

     u_ts_2d(:,:) = u_ts 
!----------- compute dust emission ------------
      dust_emis(:,:)   = 0.0
      dust_setl(:,:)   = 0.0
      dust_dt(:,:,:) = 0.0

!----------- dust sources on local grid
     source(:,:)=0.0
     call interpolator(dust_source_interp, Time, source, &
                       trim(dust_source_name(1)), is, js)
! Send the dust source data to the diag_manager for output.
     if (id_dust_source > 0 ) &
          used = send_data ( id_dust_source, source , Time_next )

      where ( frac_land.gt.0.1 .and. w10m .gt. u_ts_2d )
          dust_emis = CH * frac_s(i_DU)*source * frac_land &
             * w10m**2 * (w10m - u_ts_2d)
      endwhere
      dust_dt(:,:,kd)=dust_dt(:,:,kd)+dust_emis(:,:)/pwt(:,:,kd)*mtv

! Send the emission data to the diag_manager for output.
      if (id_dust_emis(i_DU) > 0 ) then
        used = send_data ( id_dust_emis(i_DU), dust_emis, Time_next, &
              is_in=is,js_in=js )
      endif

         rcm=dustref*mtcm            ! Particles radius in centimeters
!------------------------------------------
!       Solve at the model TOP (layer plev-10)
!------------------------------------------
      do j=1,jd
        do i=1,id
          setl(:)=0.
          if (present(kbot)) then
              kb=kbot(i,j)
          else
             kb=kd
          endif
          do k=1,kb
              rhb=amin1(0.99,rh(i,j,k))
              rhb=amax1(0.01,rhb)
!----------------------------------------------------------
!     Aerosol growth with relative humidity
!----------------------------------------------------------

            rwet=dustref  ! Add any particle growth here
            ratio_r=(dustref/rwet)**3.   ! Ratio dry over wet radius cubic power
            rho_wet_dust=ratio_r*dustden+(1.-ratio_r)*DENS_H2O     ! Density of wet aerosol [kg/m3]
            viscosity = 1.458E-6 * t(i,j,k)**1.5/(t(i,j,k)+110.4)     ! Dynamic viscosity
            free_path=6.6e-8*t(i,j,k)/293.15*(PSTD_MKS/pfull(i,j,k))
            C_c=1. + free_path/dustref* &          ! Slip correction [none]
                  (1.257+0.4*exp(-1.1*dustref/free_path))
            Vdep=2./9.*C_c*GRAV*rho_wet_dust*rwet**2./viscosity   ! Settling velocity [m/s]
            rho_air = pfull(i,j,k)/t(i,j,k)/RDGAS      ! Air density [kg/m3]
            if (dust(i,j,k).gt.0.) then
              setl(k)=dust(i,j,k)*rho_air/mtv*vdep    ! settling flux [kg/m2/s]
            endif
          enddo
          dust_dt(i,j,1)=dust_dt(i,j,1)-setl(1)/pwt(i,j,1)*mtv
          dust_dt(i,j,2:kb)=dust_dt(i,j,2:kb) &
             + ( setl(1:kb-1) - setl(2:kb) )/pwt(i,j,2:kb)*mtv
          dust_setl(i,j)=setl(kb)
        enddo
      enddo 

! Send the settling data to the diag_manager for output.
      if (id_dust_setl(i_DU) > 0 ) then
        used = send_data ( id_dust_setl(i_DU), dust_setl, Time_next, &
              is_in=is,js_in=js )
      endif


!-----------------------------------------------------------------------

 end subroutine atmos_dust_sourcesink
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_dust_init">
!<OVERVIEW>
! The constructor routine for the dust module.
!</OVERVIEW>
 subroutine atmos_dust_init (lonb, latb, axes, Time, mask)
!-----------------------------------------------------------------------
real, intent(in),    dimension(:,:)               :: lonb, latb
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
integer :: n, m, logunit
!
!-----------------------------------------------------------------------
!
      integer  unit,ierr, io
      character(len=1)  :: numb(5)
      data numb/'1','2','3','4','5'/


      if (module_is_initialized) return

      call write_version_number(version, tagname)
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=dust_nml, iostat=io)
        ierr = check_nml_error(io,'dust_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=dust_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'dust_nml')
        end do
10      call close_file (unit)
#endif
      endif

      if (uthresh .le. -990) then
        u_ts = 0.
      else
        u_ts=uthresh
      endif
      if (coef_emis .le. -990) then
        ch = 1.0e-10
      else
        ch = coef_emis
      endif
!----- set initial value of dust ------------
    logunit=stdlog()
    do m=1,5

       n = get_tracer_index(MODEL_ATMOS,'dust'//numb(m))
       if (n>0) then
         ndust=n
         call set_tracer_atts(MODEL_ATMOS,'dust'//numb(m),'dust'//numb(m),'mmr')
         if (ndust > 0 .and. mpp_pe() == mpp_root_pe()) then
                write (*,30) 'dust'//numb(m),ndust
                write (logunit,30) 'dust '//numb(m),ndust
         endif       
       endif


  30        format (A,' was initialized as tracer number ',i2)
! Register a diagnostic field : emission of dust
     id_dust_emis(m) = register_diag_field ( mod_name,            &
                     'dust'//numb(m)//'_emis', axes(1:2),Time,  &
                     'dust'//numb(m)//'_emis', 'kg/m2/s',       &
                     missing_value=-999.  )

! Register a diagnostic field : settling of dust
     id_dust_setl(m) = register_diag_field ( mod_name,            &
                     'dust'//numb(m)//'_setl', axes(1:2),Time,  &
                     'dust'//numb(m)//'_setl', 'kg/m2/s',       &
                     missing_value=-999.  )
enddo
!
     id_dust_source  = register_diag_field ( mod_name,             &
                      'DU_source',axes(1:2),Time,                    &
                      'DU_source', 'none')

     call interpolator_init (dust_source_interp, trim(dust_source_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = dust_source_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )


     call write_version_number(version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

 end subroutine atmos_dust_init
!</SUBROUTINE>

!######################################################################

subroutine atmos_dust_time_vary (Time)


type(time_type), intent(in) :: Time

      call obtain_interpolator_time_slices (dust_source_interp, Time)


end subroutine atmos_dust_time_vary 


!######################################################################

subroutine atmos_dust_endts              


      call unset_interpolator_time_flag (dust_source_interp)


end subroutine atmos_dust_endts 



!#######################################################################
!<SUBROUTINE NAME="atmos_dust_end">
!<OVERVIEW>
!  The destructor routine for the dust module.
!</OVERVIEW>
 subroutine atmos_dust_end

      call interpolator_end (dust_source_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_dust_end
!</SUBROUTINE>

end module atmos_dust_mod

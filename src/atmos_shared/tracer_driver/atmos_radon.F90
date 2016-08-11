module atmos_radon_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Larry Horowitz
! </REVIEWER>


! <OVERVIEW>
!     This code allows the implementation of an extremely simplified 
!     radon tracer in the FMS framework.
!
!    It should be taken as the implementation of a very simple tracer 
!   which bears some characteristics of radon.
! </OVERVIEW>

! <DESCRIPTION>
!   This module presents an implementation of a tracer.
!   It should be taken as representing radon only in a rudimentary manner.
! </DESCRIPTION>

!-----------------------------------------------------------------------

use              fms_mod, only : file_exist, &
                                 write_version_number, &
                                 mpp_pe, &
                                 mpp_root_pe, &
                                 error_mesg, &
                                 FATAL,WARNING, NOTE, &
                                 stdlog
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_static_field
use   tracer_manager_mod, only : get_tracer_index
use    field_manager_mod, only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                 dry_deposition


implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_radon_sourcesink, atmos_radon_init, atmos_radon_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
integer  :: ncopies_radon = 9
 
namelist /atmos_radon_nml/  &
                            ncopies_radon


!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

!--- identification numbers for  diagnostic fields and axes ----

integer :: id_emiss

logical :: module_is_initialized=.FALSE.


!---- version number -----
character(len=128) :: version = '$Id: atmos_radon.F90,v 19.0 2012/01/06 20:30:50 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################
!<SUBROUTINE NAME="atmos_radon_sourcesink">
!<OVERVIEW>
! The routine that calculate the sources and sinks of radon.
!</OVERVIEW>
!<DESCRIPTION>
! This is a very rudimentary implementation of radon.
!
! It is assumed that the Rn222 flux is 3.69e-21 kg/m*m/sec over land 
! for latitudes < 60N
!
!   Between 60N and 70N the source  = source * .5
!
!  Rn222 has a half-life time of 3.83 days, which corresponds to an 
!  e-folding time of 5.52 days.
!
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_radon_sourcesink (lon, lat, land, pwt, radon, radon_dt, 
!                              Time, kbot)
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     Land/sea mask.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="radon" TYPE="real" DIM="(:,:,:)">
!     The array of the radon mixing ratio.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

!   <OUT NAME="radon_dt" TYPE="real" DIM="(:,:,:)">
!     The array of the tendency of the radon mixing ratio.
!   </OUT>
 subroutine atmos_radon_sourcesink (lon, lat, land, pwt, radon, radon_dt,  &
                              Time, kbot)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:,:) :: pwt, radon
   real, intent(out), dimension(:,:,:) :: radon_dt
     type(time_type), intent(in) :: Time     
integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(radon,1),size(radon,2),size(radon,3)) ::  &
         source, sink
logical, dimension(size(radon,1),size(radon,2)) ::  maskeq,masknh
   real  radon_flux, dtr, deg60, deg70, deg300, deg336
integer  i,j,kb,id,jd,kd
!-----------------------------------------------------------------------

      id=size(radon,1); jd=size(radon,2); kd=size(radon,3)

      dtr=acos(0.0)/90.
      deg60=60.*dtr; deg70=70.*dtr; deg300=300.*dtr; deg336=336.*dtr

!----------- compute radon source ------------
!
!  rn222 flux is 3.69e-21 kg/m*m/sec over land for latitudes lt 60n
!   between 60n and 70n the source  = source * .5
!
!  molecular wt. of air is 28.9644 gm/mole
!  molecular wt. of radon is 222 gm/mole
!  scaling facter to get reasonable mixing ratio is 1.e+21
!
!  source = 3.69e-21 * g * 28.9644 * 1.e+21/(pwt * 222.) or
!
!  source = g * .4814353 / pwt
!
!  must initialize all rn to .001
!

      radon_flux = 3.69e-21 * 28.9644 * 1.e+21 / 222.
      source = 0.0
      maskeq = (land > 0.5) .and. lat > -deg60 .and. lat < deg60
      masknh = (land > 0.5) .and. lat >= deg60 .and. lat < deg70

      if (present(kbot)) then
          do j=1,jd
          do i=1,id
             kb=kbot(i,j)
             if (maskeq(i,j)) source(i,j,kb)=radon_flux/pwt(i,j,kb)
             if (masknh(i,j)) source(i,j,kb)=0.5*radon_flux/pwt(i,j,kb)
          enddo
          enddo
      else
          where (maskeq) source(:,:,kd)=radon_flux/pwt(:,:,kd)
          where (masknh) source(:,:,kd)=0.5*radon_flux/pwt(:,:,kd)
          where (masknh .and. lon > deg300 .and. lon < deg336)  &
               source(:,:,kd)=0.0
      endif


!------- compute radon sink --------------
!
!  rn222 has a half-life time of 3.83days 
!   (corresponds to an e-folding time of 5.52 days)
!
!  sink = 1./(86400.*5.52) = 2.09675e-6
!

    where (radon(:,:,:) >= 0.0)
       sink(:,:,:) = -2.09675e-6*radon(:,:,:)
    elsewhere
       sink(:,:,:) = 0.0
    endwhere

!------- tendency ------------------

      radon_dt=source+sink
      

!-----------------------------------------------------------------------

 end subroutine atmos_radon_sourcesink
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_radon_init">
!<OVERVIEW>
! The constructor routine for the radon module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the radon module.
!</DESCRIPTION>
!<TEMPLATE>
!call radon_init (r, mask, axes, Time)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
 subroutine atmos_radon_init (r, axes, Time, nradon, mask)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
real,             intent(inout), dimension(:,:,:,:) :: r
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
integer, dimension(:), pointer                         :: nradon
real, intent(in), dimension(:,:,:), optional        :: mask

integer :: n
!
!-----------------------------------------------------------------------
!
      integer  logunit
      character(len=64) ::  search_name
      character(len=4) ::  chname
      integer :: nn

      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number(version, tagname)
      logunit=stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=atmos_radon_nml )
 
      if (ncopies_radon > 9) then
        call error_mesg ('atmos_radonm_mod', &
          'currently no more than 9 copies of the radon tracer '//&
                                               'are allowed', FATAL)
      endif
      allocate (nradon(ncopies_radon))
      nradon = -1

      do nn=1,ncopies_radon
        write (chname,'(i1)') nn
        if (nn > 1) then
          search_name = 'radon_'// trim(chname)
        else
          search_name = 'radon'
        endif
!----- set initial value of radon ------------

       n = get_tracer_index(MODEL_ATMOS,search_name)
       if (n>0) then
         nradon(nn)=n
         if (nradon(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) trim(search_name), nradon(nn)
         if (nradon(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (logunit,30) trim(search_name), nradon(nn)
       endif

      end do

  30        format (A,' was initialized as tracer number ',i2)
!

! Register a static field for the emissions of your tracer
     id_emiss = register_static_field ( 'tracers',                    &
                     'rnemiss', axes(1:2),       &
                     'rnemiss', 'g/m2/s')

      module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

 end subroutine atmos_radon_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_radon_end">
!<OVERVIEW>
!  The destructor routine for the radon module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_radon_end
!</TEMPLATE>
 subroutine atmos_radon_end
 
      module_is_initialized = .FALSE.

 end subroutine atmos_radon_end
!</SUBROUTINE>


end module atmos_radon_mod




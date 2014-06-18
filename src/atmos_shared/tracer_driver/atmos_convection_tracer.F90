module atmos_convection_tracer_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Richard Hemler
! </CONTACT>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!                    
! </REVIEWER>


! <OVERVIEW>
!     This code allows the incorporation of an arbitrarily-specified   
!     tracer for testing within the donner_deep module.
!
!    This module is to serve as a testbed for assessing convective 
!    transport of tracers. 
! </OVERVIEW>

! <DESCRIPTION>
!   This module presents an implementation of an arbirary tracer, 
!   including its convective transport by the donner_deep module.
! </DESCRIPTION>

!-----------------------------------------------------------------------

use              fms_mod,       only : file_exist, &
                                       write_version_number, &
                                       error_mesg, &
                                       FATAL,WARNING,NOTE, &
                                       mpp_pe, mpp_root_pe, stdlog
use     time_manager_mod,       only : time_type
use     diag_manager_mod,       only : send_data,            &
                                       register_static_field
use   tracer_manager_mod,       only : get_tracer_index
use    field_manager_mod,       only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                       dry_deposition


implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_cnvct_tracer_sourcesink,  &
        atmos_convection_tracer_init,        &
        atmos_convection_tracer_end

!-----------------------------------------------------------------------
!----------- namelist -------------------

integer  :: ncopies_cnvct_trcr = 9

namelist /atmos_convection_tracer_nml/  &
                                        ncopies_cnvct_trcr

!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'


!--- identification numbers for  diagnostic fields and axes ----

integer :: id_emiss

logical :: module_is_initialized=.FALSE.


!---- version number -----
character(len=128) :: version = '$Id: atmos_convection_tracer.F90,v 19.0 2012/01/06 20:30:16 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################
!<SUBROUTINE NAME="atmos_cnvct_tracer_sourcesink">
!<OVERVIEW>
! The routine that calculate the sources and sinks of the 
! convection tracer.
!</OVERVIEW>
!<DESCRIPTION>
! This is an implementation of an arbitrarily-specified tracer.
! At this time it is assumed to have no source or sink.
!
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_cnvct_tracer_sourcesink (lon, lat, land, pwt, convtr, 
!                                         convtr_dt, Time, is, ie, 
!                                         js, je, kbot)
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
!   <IN NAME="convtr" TYPE="real" DIM="(:,:,:)">
!     The array of the convection tracer mixing ratio.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="is, ie, js, je" TYPE="integer">
!     Local domain boundaries.
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

!   <OUT NAME="convtr_dt" TYPE="real" DIM="(:,:,:)">
!     The array of the tendency of the convection tracer mixing ratio.
!   </OUT>
 subroutine atmos_cnvct_tracer_sourcesink (lon, lat, land, pwt,&
                                                convtr, convtr_dt,  &
                                                Time, is, ie, js, je, &
                                                kbot)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:,:) :: pwt, convtr
   real, intent(out), dimension(:,:,:) :: convtr_dt
     type(time_type), intent(in) :: Time     
   integer,           intent(in)       :: is, ie, js, je
integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(convtr,1),size(convtr,2),size(convtr,3)) ::  &
         source, sink
!-----------------------------------------------------------------------


!------  define source and sink of convection_tracer -------
!
!   it is currently assumed that the convection tracer has no source
!   or sink

      source = 0.
      sink   = 0.

!------- tendency ------------------

      convtr_dt = source + sink
      

!-----------------------------------------------------------------------

 end subroutine atmos_cnvct_tracer_sourcesink
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_convection_tracer_init">
!<OVERVIEW>
! The constructor routine for the convection tracer module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the convection tracer module.
!</DESCRIPTION>
!<TEMPLATE>
!call convection_tracer_init (r, phalf, mask, axes, Time)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!      pressure at model interface levels
!   </IN>
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
 subroutine atmos_convection_tracer_init (r, phalf, axes, Time, &
                                          nconvect, mask)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
real,             intent(inout), dimension(:,:,:,:) :: r
real,             intent(in),    dimension(:,:,:)   :: phalf
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
integer, dimension(:), pointer                         :: nconvect
real, intent(in), dimension(:,:,:), optional        :: mask

integer :: n
character(len=64) ::  search_name (10)
character(len=4) ::  chname
integer :: nn
!
!-----------------------------------------------------------------------
!
      real, dimension (size(r,1), size(r,2), size(r,3)) :: xgcm, pfull

      real :: xba = 1.0
      integer :: nlev, k, logunit
      character(len=64) :: filename

      nlev = size(r,3)

!---------------------------------------------------------------------
      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)
      logunit=stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=atmos_convection_tracer_nml )

!----- set initial value of convection tracer ------------

       if (ncopies_cnvct_trcr > 9) then
         call error_mesg ('atmos_convection_tracer_mod', &
         'currently no more than 9 copies of the convection tracer '//&
                                             'are allowed', FATAL)
       endif
       allocate (nconvect(ncopies_cnvct_trcr))
       nconvect = -1
 
       
        do nn=1,ncopies_cnvct_trcr
          write (chname,'(i1)') nn
          if (nn > 1) then
          search_name(nn) = 'cnvct_trcr_'// trim(chname)
          else
          search_name(nn) = 'cnvct_trcr'
          endif

       n = get_tracer_index(MODEL_ATMOS,search_name(nn) )
       if (n>0) then
         nconvect(nn)=n
         if (nconvect(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) trim(search_name(nn))  ,nconvect(nn)
         if (nconvect(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (logunit,30) trim(search_name(nn))  ,nconvect(nn)
       endif

      end do

  30        format (A,' was initialized as tracer number ',i2)
!

! Register a static field for the emissions of your tracer
     id_emiss = register_static_field ( 'tracers',                    &
                     'rnemiss', axes(1:2),       &
                     'rnemiss', 'g/m2/s')

!---------------------------------------------------------------------
!    if a convection_tracer.res file exists, it will have been prev-
!    iously processed. there is no need to do anything here.
!---------------------------------------------------------------------
      do nn = 1, ncopies_cnvct_trcr 
        if (nconvect(nn) > 0) then
          filename = 'INPUT/tracer_' //trim(search_name(nn)) // '.res'
          if (file_exist (filename)) then

!--------------------------------------------------------------------
!    if a .res file does not exist, initialize the convection_tracer.
!--------------------------------------------------------------------
          else   
            do k=1, nlev
              pfull(:,:,k) = 0.5*(phalf(:,:,k) + phalf(:,:,k+1))
            end do
            do k=1,nlev
              xgcm(:,:,k) = xba*  &
                           exp((pfull(:,:,k) - pfull(:,:,1))/       &
                               (pfull(:,:,1) - pfull(:,:,nlev)))
            end do
            do k=1,nlev
              r(:,:,nlev+1-k, nconvect(nn)) = xgcm(:,:,k)
            end do
          endif  ! (file_exist) 
        endif
      end do



      module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

 end subroutine atmos_convection_tracer_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_convection_tracer_end">
!<OVERVIEW>
!  The destructor routine for the convection tracer module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine marks the module as uninitialized and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_convection_tracer_end
!</TEMPLATE>
 subroutine atmos_convection_tracer_end
 
      module_is_initialized = .FALSE.

 end subroutine atmos_convection_tracer_end
!</SUBROUTINE>


end module atmos_convection_tracer_mod




!FDOC_TAG_GFDL
                module original_fms_rad_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!


!-----------------------------------------------------------------------
!                 radiation interface module 
!-----------------------------------------------------------------------

use           mcm_lw_mod, only: mcm_lw_init

use    mcm_sw_driver_mod, only: mcm_sw_driver_init

use            fsrad_mod, only: rdparm_init, fsrad, co2_data

use           clouds_mod, only: clouds, clouds_init, clouds_end

use     time_manager_mod, only: time_type

use              mpp_mod, only: input_nml_file
use              fms_mod, only: fms_init, FATAL, &
                                close_file, &
                                open_namelist_file,    &
                                check_nml_error, file_exist,       &
                                error_mesg, &
                                mpp_pe, mpp_root_pe, &
                                write_version_number

use    rad_utilities_mod, only: rad_utilities_init, &
                                radiative_gases_type, &
                                cldrad_properties_type, &
                                cld_specification_type, &
                                astronomy_type, &
                                atmos_input_type, &
                                surface_type, &
                                Sw_control, &
                                rad_output_type, &
                                fsrad_output_type

implicit none 
private 

!----------- public interfaces in this module -----------------------

public    original_fms_rad_init, original_fms_rad_end, original_fms_rad

!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------
character(len=128) :: version = '$Id: original_fms_rad.F90,v 19.0 2012/01/06 20:20:49 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!   ---- list of restart versions readable by this module ----
!   (sorry, but restart version 1 will not be readable by this module)
      integer, dimension(2) :: restart_versions = (/ 2, 3 /)
!-----------------------------------------------------------------------

     logical :: module_is_initialized = .false.

!-----------------------------------------------------------------------
!   the following code converts rvco2, the volume mixing ratio of co2
!   (used by lwrad) into rco2, the mass mixing ratio (used by swrad).

!   real, parameter :: ratco2mw=1.519449738  ! unused variable (pjp)

!   real, parameter :: rvco2 = 3.3e-4        ! unused variable (pjp)
!   real, parameter :: rco2 = rvco2*ratco2mw ! unused variable (pjp)


!--------- flags -------------------------------------------------------

    logical ::  do_rad
    logical ::  use_rad_date

!-----------------------------------------------------------------------
!  ------- allocatable global data saved by this module -------
!
!   tdt_rad      = radiative (sw + lw) heating rate
!   flux_sw_surf = net (down-up) sw flux at surface
!   flux_lw_surf = downward lw flux at surface
!   coszen_angle = cosine of the zenith angle 
!                  (used for the last radiation calculation)

    type(rad_output_type),save          ::  Rad_output
    real, allocatable, dimension(:,:)   ::  solar_save

!   real, allocatable, dimension(:,:,:) :: tdt_rad
!   real, allocatable, dimension(:,:)   :: flux_sw_surf, &
!                                          flux_lw_surf, &
!                                          coszen_angle, &
!                                          solar_save, &
!                                          fracday

!-----------------------------------------------------------------------
!   ------------------- time step constant --------------------------
!
!   rad_alarm     = time interval until the next radiation calculation
!                   (integer in seconds)
!   rad_time_step = radiation time step (integer in seconds)
!   total_pts = number of grid boxes to be processed every time step
!               (note: all grid boxes must be processed every time step)
!   num_pts   = counter for current number of grid boxes processed
!             (when num_pts=0 or num_pts=total_pts certain things happen)

           integer :: num_pts, total_pts
!          integer :: rad_alarm, rad_time_step, old_time_step
           integer :: rad_alarm, rad_time_step

!-----------------------------------------------------------------------
!------- private allocatable array for time averaged input data --------

      real, allocatable, dimension(:,:,:)   :: psum, tsum, qsum
      real, allocatable, dimension(:,:)     :: asum, csum, ssum,   &
                                               fsum, rrsum
   integer, allocatable, dimension(:,:)     :: nsum
      real, allocatable, dimension(:,:,:)   :: gas_component_sum
      real, allocatable, dimension(:,:,:,:) :: cld_component_sum

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------
!
!   dt_rad  =  radiative time step in seconds
!
!   offset  =  offset for radiative time step (in seconds) 
!              note if offset=1 radition will be done on the first
!              time step after t=0, dt_rad, 2*dt_rad, etc.
!              for now do not change this value

!! the following not needed as of 11-13-01 :
! NOTE:  when reading vers = 1 sea_esf_rad.resfile, offset MUST be < the
!        model timestep
!! end of not needed
!
!   solar_constant = solar constant in watts/m2
!
!   rad_date       = fixed date (year,month,day,hour,min,sec)
!                    for radiation (solar info, ozone, clouds)
!
!   co2std = standard co2 vol. mixing ratio (either 300 or 330 ppmv)
!
!   ratio  = co2 vol. mixing ratio in units of the standard vol.
!            mixing ratio (must lie between 0.5 and 4.0)
!   rsd (repeat same day) : call radiation for the specified rad_date,
!                           running through the diurnal cycle

!  integer :: dt_rad=43200, offset=1
!  logical ::                                         &
!             do_clear_sky_diag=.false., do_average=.false., &
!      do_average_gases = .false., do_average_clouds=.false.
!  logical :: rsd=.false.
!!!   real :: solar_constant = 1.96  !(1367.44w/m2)
      real :: solar_constant = 1365.
!  integer, dimension(6) :: rad_date = (/ 0, 0, 0, 0, 0, 0 /)
      logical :: do_mcm_radiation=.false.
      real :: co2std = 330., ratio = 1.0
!     integer :: jpt = -35, ipt = -35
       real :: lat_diag = -1000.
       real :: long_diag = -1000.
!     logical :: calc_hemi_integrals = .false.
!     logical :: renormalize_sw_fluxes=.false.
!     logical :: do_bounds_chk        =.false.

! character(len=24) :: zenith_spec = '      '
!  character(len=16)    :: rad_step_physics='default'

!logical         :: all_column_radiation = .true.
!logical         :: all_level_radiation = .true.
!integer         :: topmost_radiation_level=-99

namelist /original_fms_rad_nml/ &
!                               offset, do_average, &
!                               do_average_gases, do_average_clouds, &
!                               do_clear_sky_diag,          &
!                               zenith_spec,  &
!                               rad_step_physics, &
!                               calc_hemi_integrals, &
                                lat_diag, long_diag, &
!                               renormalize_sw_fluxes, &
!                               solar_constant, rad_date,   &
                                solar_constant,             &
!                               co2std, ratio, jpt, ipt, &
                                co2std, ratio, &
!                               do_bounds_chk, &
!  all_level_radiation, &
!  topmost_radiation_level, &
!  all_column_radiation, &
!                               rsd
                                do_mcm_radiation

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_alb_sfc, id_coszen

integer, dimension(2) :: id_tdt_sw,   id_tdt_lw,  &
                         id_swdn_toa, id_swup_toa, id_olr, &
                         id_swup_sfc, id_swdn_sfc,         &
                         id_lwup_sfc, id_lwdn_sfc

!character(len=9), parameter :: mod_name = 'radiation'
character(len=16), parameter :: mod_name = 'radiation'

real :: missing_value = -999.
logical  :: do_sea_esf_rad
logical  :: do_clear_sky_pass

!-------------------------------------------------------------------- 
!   list of restart versions readable by this module ----
!-------------------------------------------------------------------- 
integer, dimension(2) :: restart_versions_sea = (/ 2, 3 /)

!--------------------------------------------------------------------
!   tdt_rad = the current radiative (sw + lw) heating rate
!--------------------------------------------------------------------

!---------------------------------------------------------------------
! diagnostics field informations 
!---------------------------------------------------------------------
integer ::             id_cosz, id_fracday 

logical   :: do_diurnal, do_annual, do_daily_mean

integer            ::  id, jd, kmin, kmax
integer            ::  ipgl, jpgl
integer            ::  vers
integer            ::  ksrad, kerad ! vertical indices on radiation grid
                                    ! ksrad ==1

integer            ::  ks, ke ! vertical indices in model grid coords
                              !  ks is topmost_radiation_level
integer            ::  iomsgs
real               ::  rh2o_lower_limit_orig=3.0E-06
real               ::  rh2o_lower_limit_seaesf=2.0E-07
real               ::  rh2o_lower_limit

     integer :: israd, ierad, jsrad, jerad

!---------------------------------------------------------------------
!---------------------------------------------------------------------

                         contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!####################################################################

! <SUBROUTINE NAME="original_fms_rad_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call original_fms_rad_init
!  </TEMPLATE>
! </SUBROUTINE>
!
   subroutine original_fms_rad_init ( lonb, latb, pref, axes, Time , &
!      co2std, ratio, kmax,  &
!                     kmax,  &
                      kmax)
!      lat_diag, long_diag, ipgl, jpgl)
!                           ipgl, jpgl)

!-----------------------------------------------------------------------
           integer, intent(in)  :: kmax
!   real, intent(in) :: co2std, ratio
           real, intent(in), dimension(:,:) :: lonb, latb
           real, intent(in), dimension(:,:) :: pref
        integer, intent(in), dimension(4)   :: axes
type(time_type), intent(in)                 :: Time
!real, intent(inout) :: lat_diag, long_diag
!integer, intent(out) :: ipgl, jpgl
!-----------------------------------------------------------------------
      integer :: unit, io, ierr


!     character(len=4) :: chvers

      real, dimension(size (lonb,1)-1) :: idiff
      real, dimension(size (latb,2)-1) :: jdiff
      integer :: i, j, minindx
      real  :: mindiff
      integer :: id, jd

      if(module_is_initialized) return
      call fms_init
      call rad_utilities_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=original_fms_rad_nml, iostat=io)
      ierr = check_nml_error(io,'original_fms_rad_nml')
#else   
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
          read  (unit, nml=original_fms_rad_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'original_fms_rad_nml')
        enddo
  10    call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      if ( mpp_pe() == mpp_root_pe() ) then
        call write_version_number(version, tagname)
        write (unit, nml=original_fms_rad_nml)
      endif


    jd = size(latb,2)-1
    id = size(lonb,1)-1
    
!--------------------------------------------------------------------
!  define location of radiation diagnostics column, if present
!--------------------------------------------------------------------

     if (lat_diag /= -1000. .and. long_diag /= -1000.) then
       if ( (lat_diag < -90. .or. lat_diag > 90.) .or.  &
           (long_diag <    0. .or. long_diag > 360.) ) then 
          call error_mesg ('original_fms_rad_mod', & 
       ' bad values specified for lat_diag or long_diag', FATAL) 
       endif
       
!      jd = size(latb,2)-1
!      id = size(lonb,1) - 1
       lat_diag = lat_diag*4.0*atan(1.0)/180.0
       long_diag = long_diag*4.0*atan(1.0)/180.0
       if (  (lat_diag >= latb(1,1) .and. lat_diag <= latb(1,jd+1))   .and. &
          (long_diag >= lonb(1,1) .and. long_diag <= lonb(id+1,1))  ) then
         do j=1,jd
           jdiff(j) = abs ( 0.5*(latb(1,j)+latb(1,j+1)) - lat_diag)  
         end do
          mindiff = 10.0
          minindx = 1
         do j=1,jd
          if (jdiff(j) .lt. mindiff) then
            mindiff = jdiff(j)
            minindx = j
          endif 
          end do
          jpgl = minindx
         do i=1,id
           idiff(i) = abs ( 0.5*(lonb(i,1)+lonb(i+1,1)) - long_diag)  
         end do
          mindiff = 10.0
          minindx = 1
         do i=1,id
          if (idiff(i) .lt. mindiff) then
            mindiff = idiff(i)
            minindx = i
          endif 
          end do
          ipgl = minindx
         else
            jpgl = 0  
            ipgl = 0  
         endif
         else
            jpgl = 0  
            ipgl = 0  
         endif

!        print *,  get_my_pe() , ipgl, jpgl

!-----------------------------------------------------------------------

         call rdparm_init (kmax)

         if(do_mcm_radiation) then
           call mcm_lw_init (id,jd,kmax)
           call mcm_sw_driver_init(kmax)
         else
           call co2_data (co2std, ratio, pref)
         endif

         call clouds_init ( lonb, latb, axes, Time )

         module_is_initialized = .true.

   end subroutine original_fms_rad_init

!#######################################################################

! <SUBROUTINE NAME="original_fms_rad_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call original_fms_rad_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine original_fms_rad_end

!   --- terminate other modules ---

if(.not.module_is_initialized) return
 
     call clouds_end

!-----------------------------------------------------------------------
module_is_initialized = .false.

end subroutine original_fms_rad_end
      
!###################################################################

!subroutine original_fms_rad (is, ie, js, je, kerad, ipgl, jpgl,   &
!subroutine original_fms_rad (lat_in, lon_in, is, ie, js, je, kerad,   &
!subroutine original_fms_rad (is, ie, js, je, kerad, lat_in, lon_in, &
! <SUBROUTINE NAME="original_fms_rad">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call original_fms_rad (is, ie, js, je, phalf, lat_in, lon_in, &
!                do_clear_sky_pass, &
!                Rad_time, Time_diag, Atmos_input, &
!                Surface, &
!                Astro, Rad_gases, Cldrad_props, Cld_spec, &
!                Fsrad_output, mask, kbot) 
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <IN NAME="phalf" TYPE="real">
! 
!  </IN>
!  <IN NAME="lat_in" TYPE="real">
! 
!  </IN>
!  <IN NAME="lon_in" TYPE="real">
! 
!  </IN>
!  <IN NAME="do_clear_sky_pass" TYPE="logical">
! 
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
! 
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
! 
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
! 
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
! 
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
! 
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
! 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
! 
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
! 
!  </IN>
!  <INOUT NAME="Fsrad_output" TYPE="fsrad_output_type">
! 
!  </INOUT>
!  <IN NAME="mask" TYPE="real">
! 
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine original_fms_rad (is, ie, js, je, phalf, lat_in, lon_in, &
                             do_clear_sky_pass, &
                             Rad_time, Time_diag, Atmos_input, &
                             Surface, &
                             Astro, Rad_gases, Cldrad_props, Cld_spec, &
                             Fsrad_output, mask, kbot) 

!--------------------------------------------------------------------
!integer, intent(in) :: is, ie, js, je, kerad, ipgl, jpgl
!integer, intent(in) :: is, ie, js, je, kerad              
integer, intent(in) :: is, ie, js, je
real, intent(in), dimension(:,:,:) :: phalf
real, dimension(:,:), intent(in) :: lat_in, lon_in
type(time_type), intent(in) :: Rad_time, Time_diag
logical, intent(in) :: do_clear_sky_pass
type(atmos_input_type), intent(in) :: Atmos_input
type(surface_type), intent(in) :: Surface
type(astronomy_type), intent(in) :: Astro        
type(radiative_gases_type), intent(in) :: Rad_gases
type(cldrad_properties_type), intent(in) :: Cldrad_props
type(cld_specification_type), intent(in) :: Cld_spec
type(fsrad_output_type), intent(inout)  :: Fsrad_output
real, dimension(:,:,:), intent(in), optional :: mask
integer, dimension(:,:), intent(in), optional  :: kbot
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
! local variables

!    real,dimension(size( lat_in,1),size( lat_in,2), ke+1     )   ::   &
!                    press_in, temp_in, pflux_in, tflux_in, &
!      real, dimension(ie-is+1, je-js+1, kerad+1) :: &
      real, dimension(ie-is+1, je-js+1, size(Atmos_input%press,3)) :: &
                    cldamt,emcld, cuvrf,cirrf, cirab

!    real,dimension(size( lat_in,1),size( lat_in,2), ke               )   ::  &
!                   rh2o_in,  deltaz_in, cloud_water_in, cloud_ice_in,&
!                   q,qo3

!    real,dimension(size( lat_in,1),size( lat_in,2)     )   ::  &
!                      land_in, tsfc, asfc_in, psfc 
!     real, dimension(ie-is+1, je-js+1, kerad+1) :: temp, press, &
     real, dimension(ie-is+1, je-js+1, size(Atmos_input%press,3)) :: temp, press, &
                                                pflux, tflux

!    real, dimension(ie-is+1, je-js+1, kerad+1) :: temp, press, &
!                                               pflux, tflux
!    real, dimension(ie-is+1, je-js+1, kerad  ) :: rh2o        , &
     real, dimension(ie-is+1, je-js+1, size(Atmos_input%press,3)-1  ) :: rh2o        , &
                                                   q, qo3
!    integer,dimension(size(temp_in,1),size(temp_in,2), size(temp_in,3) ) ::  &
!  integer, dimension(ie-is+1, je-js+1, kerad+1  ) :: &
   integer, dimension(ie-is+1, je-js+1, size(Atmos_input%press,3)  ) :: &
                    ktop,kbtm, ktopsw, kbtmsw

!    integer,dimension(size(temp_in,1),size(temp_in,2))   :: nclds
      integer, dimension(ie-is+1, je-js+1 ) ::   nclds

     integer ::                   ip, jp
     logical :: no_clouds



     real :: rrsun

      real, dimension(ie-is+1, je-js+1) :: &
!               lat, lon, land, asfc, cosz, fracday, solar, &
                          land, asfc, cosz, fracday, solar, &
!23             lat ,lon ,land, asfc, cosz, fracday, solar, &
                 psfc, tsfc
!     real, dimension(ie-is+1, je-js+1, kerad  ) :: &
!     real, dimension(ie-is+1, je-js+1, size(Atmos_input%press,3)  ) :: &
!                             cloud_water, cloud_ice, deltaz, qo3_in, &
!                       qo3_out          

      integer :: ierad, jerad, kerad
      real :: rvco2

           ierad = ie - is + 1
           jerad = je -js + 1
           kerad = size(Atmos_input%press,3) - 1
!   print *, get_my_pe(), is, ie, js, je, kerad, do_clear_sky_pass
!   print *, 'lat', get_my_pe(), lat_in
!   print *, 'lon', get_my_pe(), lon_in
!   print *, 'ipgl,jpgl', get_my_pe(), ipgl, jpgl
!--------------------------------------------------------------------
!    define variables on the radiation grid which will be passed into
!    the radiation routines.
!---------------------------------------------------------------------
               press = Atmos_input%press
               temp  = Atmos_input%temp 
               rh2o  = Atmos_input%rh2o 
               asfc  = Surface%asfc
               psfc  = Atmos_input%psfc 
               tsfc  = Atmos_input%tsfc 
               pflux = Atmos_input%pflux
               tflux = Atmos_input%tflux
!              deltaz = Atmos_input%deltaz
               land  = Surface%land
!              cloud_water = Atmos_input%cloud_water
!              cloud_ice = Atmos_input%cloud_ice

!-------------------------------------------------------------------
!    make mods necessary for use with original fms radiation code:
!    1) the value expected for cosz in the diurnally-varying case is the
!    product of cosz and fracday; and 2) the solar constant is included 
!    in the solar array, in contrast to the sea-esf radiation, where it 
!    is resident in the shortwave module and is applied in that module.
!--------------------------------------------------------------------
!              cosz = Astro%cosz
!              if (Astro%do_diurnal) then
               if (Sw_control%do_diurnal) then
                 cosz = Astro%cosz*Astro%fracday
               else
                 cosz = Astro%cosz
               endif
               fracday = Astro%fracday
!       solar = Astro%solar
               solar = Astro%solar*solar_constant
               rrsun = Astro%rrsun

               qo3 = Rad_gases%qo3
               rvco2 = Rad_gases%rrvco2
     
        allocate (Fsrad_output%tdtsw    ( ierad, jerad, kerad) )
        allocate (Fsrad_output%tdtlw    ( ierad, jerad, kerad) )
        allocate (Fsrad_output%swdns    ( ierad, jerad       ) )
        allocate (Fsrad_output%swups    ( ierad, jerad       ) )
        allocate (Fsrad_output%lwdns    ( ierad, jerad       ) )
        allocate (Fsrad_output%lwups    ( ierad, jerad       ) )
        allocate (Fsrad_output%swin     ( ierad, jerad       ) )
        allocate (Fsrad_output%swout    ( ierad, jerad       ) )
        allocate (Fsrad_output%olr      ( ierad, jerad       ) )

        Fsrad_output%tdtsw   = 0.
        Fsrad_output%tdtlw  = 0.
        Fsrad_output%swdns   = 0.
        Fsrad_output%swups  = 0.
        Fsrad_output%lwdns   = 0.
        Fsrad_output%lwups   = 0.
        Fsrad_output%swin     = 0.
        Fsrad_output%swout  = 0.
        Fsrad_output%olr     = 0.
     if (do_clear_sky_pass) then
        allocate (Fsrad_output%tdtsw_clr( ierad, jerad, kerad) )
        allocate (Fsrad_output%tdtlw_clr( ierad, jerad, kerad) )
        allocate (Fsrad_output%swdns_clr( ierad, jerad       ) )
        allocate (Fsrad_output%swups_clr( ierad, jerad       ) )
        allocate (Fsrad_output%lwdns_clr( ierad, jerad       ) )
        allocate (Fsrad_output%lwups_clr( ierad, jerad       ) )
        allocate (Fsrad_output%swin_clr ( ierad, jerad       ) )
        allocate (Fsrad_output%swout_clr( ierad, jerad       ) )
        allocate (Fsrad_output%olr_clr  ( ierad, jerad       ) )

        Fsrad_output%tdtsw_clr = 0.
        Fsrad_output%tdtlw_clr =0.
        Fsrad_output%swdns_clr = 0.
        Fsrad_output%swups_clr = 0.
        Fsrad_output%lwdns_clr = 0.
        Fsrad_output%lwups_clr = 0.
        Fsrad_output%swin_clr  = 0.
        Fsrad_output%swout_clr = 0.
        Fsrad_output%olr_clr   = 0.
     endif

!----------------------------------------------------------------------
!   determine if radiation diagnostics column is present in current 
!   physics window. if so, determine its coordinates in the 
!   physics_window space. if not present, set ip and jp to 0.
!----------------------------------------------------------------------
!        call define_diag_column (is, ie, js, je, ipgl, jpgl, lat,   &
        call define_diag_column (is, ie, js, je,             lat_in,   &
                                 lon_in, ip, jp)
       
!-----------------------------------------------------------------------
!--------------- loop for clear sky diagnostics option -----------------

      if (do_clear_sky_pass) then
        no_clouds = .true.
!  redefine values to be input to clouds to be the time-averaged ones
!  now available, when that option is selected.
!  this change MAY affect isccp diagnostics from strat_cloud_mod.

     q = rh2o/(1.+rh2o)

!      call clouds (is,js, no_clouds, Rad_time, Time_diag, lat, land, tsfc, &
      call clouds (is,js, no_clouds, Rad_time, Time_diag, lat_in, land, tsfc, &
!                   press(:,:,1:kmax), pflux, temp(:,:,1:kmax), q,  &
                    press(:,:,1:kerad), pflux, temp(:,:,1:kerad), q,  &
                            cosz ,nclds, ktopsw, kbtmsw, ktop, kbtm,   &
                    cldamt, cuvrf, cirrf, cirab, emcld, mask, kbot)
!-----------------------------------------------------------------------
!----------------------------- radiation -------------------------------

     if (present(kbot)) then
!        call fsrad (ip,jp,press,temp,rh2o,Rad_gases%qo3,  &
         call fsrad (ip,jp,press,temp,rh2o,          qo3,  &
                     phalf,do_mcm_radiation,               &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                            cosz ,          solar,                            &
                     Fsrad_output%swin_clr,Fsrad_output%swout_clr, &
                     Fsrad_output%olr_clr,Fsrad_output%swups_clr, &
                     Fsrad_output%swdns_clr,Fsrad_output%lwups_clr, &
                     Fsrad_output%lwdns_clr,  &
                     Fsrad_output%tdtsw_clr,Fsrad_output%tdtlw_clr, &
                     kbot,psfc)
     else
!        call fsrad (ip,jp,press,temp,rh2o,Rad_gases%qo3,  &
         call fsrad (ip,jp,press,temp,rh2o,          qo3,  &
                     phalf,do_mcm_radiation,               &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                  cosz ,solar,                            &
                     Fsrad_output%swin_clr,Fsrad_output%swout_clr, &
                     Fsrad_output%olr_clr,Fsrad_output%swups_clr, &
                     Fsrad_output%swdns_clr,Fsrad_output%lwups_clr, &
                     Fsrad_output%lwdns_clr,  &
                     Fsrad_output%tdtsw_clr,Fsrad_output%tdtlw_clr )
     endif

  endif  ! (do_clear_sky_pass)

        no_clouds = .false.

!  redefine values to be input to clouds to be the time-averaged ones
!  now available, when that option is selected.
!  this change MAY affect isccp diagnostics from strat_cloud_mod.

     q = rh2o/(1.+rh2o)

!     call clouds (is,js, no_clouds, Rad_time, Time_diag, lat, land, tsfc, &
      call clouds (is,js, no_clouds, Rad_time, Time_diag, lat_in, land, tsfc, &
!                   press(:,:,1:kmax), pflux, temp(:,:,1:kmax), q,  &
                    press(:,:,1:kerad), pflux, temp(:,:,1:kerad), q,  &
                            cosz ,nclds, ktopsw, kbtmsw, ktop, kbtm,   &
                    cldamt, cuvrf, cirrf, cirab, emcld, mask, kbot)
 
!-----------------------------------------------------------------------
!----------------------------- radiation -------------------------------

     if (present(kbot)) then
!        call fsrad (ip,jp,press,temp,rh2o,Rad_gases%qo3,  &
         call fsrad (ip,jp,press,temp,rh2o,          qo3,  &
                     phalf,do_mcm_radiation,               &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                  cosz ,solar,                            &
                     Fsrad_output%swin,Fsrad_output%swout, &
                     Fsrad_output%olr,Fsrad_output%swups, &
                     Fsrad_output%swdns,Fsrad_output%lwups, &
                     Fsrad_output%lwdns,  &
                     Fsrad_output%tdtsw,Fsrad_output%tdtlw, kbot,psfc)
     else
!        call fsrad (ip,jp,press,temp,rh2o,Rad_gases%qo3,  &
         call fsrad (ip,jp,press,temp,rh2o,          qo3,  &
                     phalf,do_mcm_radiation,               &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                            cosz ,solar,                            &
                     Fsrad_output%swin,Fsrad_output%swout, &
                     Fsrad_output%olr,Fsrad_output%swups, &
                     Fsrad_output%swdns,Fsrad_output%lwups, &
                     Fsrad_output%lwdns,  &
                     Fsrad_output%tdtsw,Fsrad_output%tdtlw)
     endif

     if (do_clear_sky_pass) then
       Fsrad_output%npass = 2
      else

       Fsrad_output%npass = 1
      endif

!---------------------------------------------------------------------

end subroutine original_fms_rad

!#####################################################################

!subroutine define_diag_column (is, ie, js, je, ipgl, jpgl, lat,  &
! <SUBROUTINE NAME="define_diag_column">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_diag_column (is, ie, js, je,             lat,  &
!                lon, ip, jp) 
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <IN NAME="lat" TYPE="real">
! 
!  </IN>
!  <IN NAME="lon" TYPE="real">
! 
!  </IN>
!  <OUT NAME="ip" TYPE="integer">
! 
!  </OUT>
!  <OUT NAME="jp" TYPE="integer">
! 
!  </OUT>
! </SUBROUTINE>
!
subroutine define_diag_column (is, ie, js, je,             lat,  &
                               lon, ip, jp) 

!--------------------------------------------------------------------
!integer, intent(in) :: is, ie, js, je, ipgl, jpgl
integer, intent(in) :: is, ie, js, je
real, dimension(:,:), intent(in) :: lat, lon
integer, intent(out) :: ip, jp
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   determine if radiation diagnostics column is present in current 
!   physics window. if so, determine its coordinates in the 
!   physics_window space. if not present, set ip and jp to 0.
!----------------------------------------------------------------------
       
     if (jpgl == 0 .and. ipgl == 0) then
       ip = 0
       jp = 0
     else
       if  (   (is <= ipgl .and. ie >= ipgl) .and.   &
               (js <= jpgl .and. je >= jpgl) ) then
          ip = ipgl - is + 1
          jp = jpgl - js + 1
          print *, 'long and lat of point for radiation diagnostics', &
                    lon(ip,jp)*45./atan(1.0), lat(ip,jp)*45./atan(1.0)      
        else
          ip = 0
          jp = 0
        endif
      endif

end subroutine define_diag_column 

!######################################################################
      
                 end module original_fms_rad_mod


module rh_clouds_mod

!=======================================================================
!
!                          RH_CLOUDS MODULE
!
!=======================================================================

use mpp_mod,    only : input_nml_file
use fms_mod,    only : error_mesg, FATAL, file_exist,    &
                       check_nml_error, open_namelist_file,       &
                       close_file, mpp_pe, mpp_root_pe, &
                       write_version_number, stdlog
use fms_io_mod, only : restore_state, &
                       register_restart_field, restart_file_type, &
                       save_restart, get_mosaic_tile_file

!=======================================================================

implicit none
private

public  rh_clouds, rh_clouds_init, rh_clouds_end,  &
        rh_clouds_sum, rh_clouds_avg, do_rh_clouds

!=======================================================================

!
!  The public interface in this module are RH_CLOUDS_INIT 
!                                          RH_CLOUDS
!                                          
!  SUBROUTINE RH_CLOUDS_INIT
!  -- no input or output -- initializes module by reading namelist 
!  
!  SUBROUTINE RH_CLOUDS(RH, P_FULL, P_SURF, ZENITH, DEG_LAT,
!          N_CLOUD, TOP, BOT, CLDAMT,
!          ALB_UV,ALB_NIR,ABS_UV,ABS_NIR,EMISS)
!
!  input -- 
!
!     real, rh(:,:,:)     -- relative humidity(nlon,nlat,nlev)  
!                            third index runs from top of atmosphere to bottom 
!     real, p_full(:,:,:) -- pressure(nlon, nlat, nlev) at rh levels 
!     real, p_surf(:,:)   -- surface pressure(nlon,nlat)
!                            p_full and p_surf must be in same units
!     real, zenith(:,:)   -- cosine of zenith angle (nlon, nlat) 
!     real, deg_lat(:)    -- latitude in degrees (nlon,nlat) 
!
!  output --
!
!     integer, n_cloud(:,:)  -- number of distinct clouds (nlon, nlat)
!                  
!     integer, top(:,:,:) -- 
!     integer, bot(:,:,:) -- (nlon,nlat,max) --  max must be
!                            larger than any value of n_cloud 
!                            max = (size(rh,3)+1)/2 is safe
!                            the n'th cloud at horizontal location i,j
!                            fills the levels from top(i,j,n) to bot(i,j,n)
!                            inclusive (i.e., if the cloud is only one 
!                            level thick, then top = bot)
!                            cloud numbering starts from top of atmosphere
!     real, cldamt       --  horizontal area covered by nth cloud
!                            each cloud covers total area currently
!                            (i.e. cldamt = 1.)
!
!     real, alb_uv(:,:,:) -- (lon, lat, max)
!                            short wave albedo for each cloud
!     real, alb_nir(:,:,:) - near infrared albedo for each cloud
!     real, abs_uv(:,:,:) -- (lon, lat, max)
!                            short wave absorption coeff for each cloud = 0
!     real, abs_nir(:,:,:) - near infrared absorption coeff for each cloud
!     real, emiss(:.:.:) --(lon, lat, max)
!                        infrared emissivity for each cloud
!=======================================================================
!----------------- data for rh averaging code --------------------------

    real,    allocatable, dimension (:,:,:) :: rhsum
    integer, allocatable, dimension (:,:)   :: nsum

!-----------------------------------------------------------------------

interface rh_clouds
    module procedure  rh_clouds_3d, rh_clouds_2d, rh_clouds_1d
end interface

!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: rh_clouds.F90,v 19.0 2012/01/06 20:12:00 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

! sigma coordinate boundary between high and middle clouds varies linearly
!   in latitude between high_middle_pole at the pole and high_middle_eq
!   at the equator, and similarly for the boundary between middle and low
!   clouds

real :: high_middle_pole = 0.7
real :: high_middle_eq   = 0.4
real :: middle_low_pole  = 0.85
real :: middle_low_eq    = 0.7

! cloud is present when relative humidity >= rh_crit, which varies liearly
!   in sigma from rh_crit_top at sigma = 0 to rh_crit_bot at sigma = 1

real :: rh_crit_bot    = 1.00
real :: rh_crit_top    = 0.90


!  near infrared absorption coeffs 

real :: high_abs    = 0.04
real :: middle_abs  = 0.30
real :: low_abs     = 0.40

! infrared emissivities

real :: high_emiss    = 0.6
real :: middle_emiss  = 1.0
real :: low_emiss     = 1.0

real :: tuning_coeff_low_cld = 1.0

!  flag for time averaging rh

logical :: do_average = .false.
logical :: do_mcm_no_clouds_top = .false.
logical :: do_mcm_crit_rh = .false.

! albedos are computed from a table look-up as function of zenith angle

namelist /rh_clouds_nml/ high_middle_pole, high_middle_eq, &
                         middle_low_pole , middle_low_eq, &
                         rh_crit_bot, rh_crit_top, &
                         high_abs, middle_abs, low_abs, &
                         high_emiss, middle_emiss, low_emiss, &
                         do_average, tuning_coeff_low_cld, &
                         do_mcm_no_clouds_top, do_mcm_crit_rh

!=======================================================================

!  OTHER MODULE VARIABLES

!--- for netcdf restart
type(restart_file_type), pointer, save :: RH_restart => NULL()

logical :: module_is_initialized = .false.

contains

!#######################################################################

subroutine rh_clouds_init (nlon, nlat, nlev)

integer, intent(in) :: nlon, nlat, nlev

integer :: unit, ierr, io, logunit, id_restart

      if (module_is_initialized) return

!------------------- read namelist input -------------------------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=rh_clouds_nml, iostat=io)
      ierr = check_nml_error(io,'rh_clouds_nml')
#else   
      if (file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=rh_clouds_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'rh_clouds_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit, nml=rh_clouds_nml)
      endif

!---------- initialize for rh cloud averaging -------------------------

      allocate (rhsum(nlon,nlat,nlev), nsum(nlon,nlat))
      allocate(RH_restart)

      id_restart = register_restart_field(RH_restart, 'rh_clouds.res.nc', 'nsum',  nsum)
      id_restart = register_restart_field(RH_restart, 'rh_clouds.res.nc', 'rhsum', rhsum)
      if (file_exist('INPUT/rh_clouds.res.nc')) then
        call restore_state(RH_restart)
      else if (file_exist('INPUT/rh_clouds.res')) then
        call error_mesg ('rh_clouds_init', &
                         'Native restart files no longer supported.', FATAL)
!        unit = open_restart_file ('INPUT/rh_clouds.res', action='read')
!        call read_data (unit, nsum)
!        call read_data (unit, rhsum)
!        call close_file (unit)
      else
        rhsum = 0.0;  nsum = 0
      endif


      module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine rh_clouds_init

!#######################################################################

subroutine rh_clouds_end


    call save_restart(RH_restart)
    module_is_initialized = .false.

end subroutine rh_clouds_end

!#######################################################################

 function do_rh_clouds ( ) result (answer)
   logical :: answer

!  returns logical value for whether rh_clouds has been initialized
!  presumably if initialized then rh_cloud will be used

   answer = module_is_initialized

 end function do_rh_clouds

!#######################################################################

 subroutine rh_clouds_sum (is, js, rh)

!-----------------------------------------------------------------------
   integer, intent(in)                   :: is, js
      real, intent(in), dimension(:,:,:) :: rh
!-----------------------------------------------------------------------
   integer :: ie, je

   ie = is + size(rh,1) - 1
   je = js + size(rh,2) - 1

!--------- use time-averaged or instantaneous clouds -----------

   if (do_average) then
       nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
      rhsum(is:ie,js:je,:) = rhsum(is:ie,js:je,:) + rh(:,:,:)
   else
       nsum(is:ie,js:je)   =  1
      rhsum(is:ie,js:je,:) = rh(:,:,:)
   endif

!-----------------------------------------------------------------------

 end subroutine rh_clouds_sum

!#######################################################################

 subroutine rh_clouds_avg (is, js, rh, ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(out), dimension(:,:,:) :: rh
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
   integer ::ie, je, num, k
!-----------------------------------------------------------------------

   if (size(rh,3) .ne. size(rhsum,3)) call error_mesg ( &
                              'rh_clouds_avg in rh_clouds_mod',  &
                              'input argument has the wrong size',FATAL)

   ie = is + size(rh,1) - 1
   je = js + size(rh,2) - 1
   num = count(nsum(is:ie,js:je) == 0)

   if (num > 0) then

!     ----- no average, return error flag -----

!!!    call error_mesg ('rh_clouds_avg in rh_clouds_mod',  &
!!!                     'dividing by a zero counter', FATAL)
       ierr = 1

   else

!      ----- compute average -----

       do k = 1, size(rh,3)
          rh(:,:,k) = rhsum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
       enddo
       ierr = 0

   endif

    nsum(is:ie,js:je)   = 0
   rhsum(is:ie,js:je,:) = 0.0
     
!-----------------------------------------------------------------------

 end subroutine rh_clouds_avg

!#######################################################################

subroutine rh_clouds_3d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:,:,:)   :: rh, p_full
real   , intent(in) , dimension(:,:)     :: p_surf, zenith, deg_lat
integer, intent(out), dimension(:,:,:)   :: top, bot
integer, intent(out), dimension(:,:)     :: n_cloud
real   , intent(out), dimension(:,:,:)   :: cldamt,emiss
real   , intent(out), dimension(:,:,:)   :: alb_uv,alb_nir,abs_uv,abs_nir


integer :: i, j, k, n, nlev, max_n_cloud
real   , dimension(size(rh,1),size(rh,2))            :: rh_crit
logical, dimension(size(rh,1),size(rh,2),size(rh,3)) :: cloud
real, dimension(size(rh,1),size(rh,2),2) ::  high_alb, middle_alb, low_alb
real, dimension(size(rh,1),size(rh,2)) :: high_middle, middle_low
real :: sig_bot

! dummy checks

 if (.not.module_is_initialized) call error_mesg( 'RH_CLOUDS in RH_CLOUD        S_MOD', &
        'module not initialized', FATAL)
 if (size(zenith,1).ne.size(rh,1)) &
              call error_mesg( 'RH_CLOUDS in RH_CLOUDS_MOD', &
             'dimensions of zenith and top do not match', FATAL)
 if (size(zenith,2).ne.size(rh,2)) &
              call error_mesg( 'RH_CLOUDS in RH_CLOUDS_MOD', &
             'dimensions of zenith and top do not match', FATAL)
 if (size(deg_lat,1).ne.size(rh,1)) &
              call error_mesg( 'RH_CLOUDS in RH_CLOUDS_MOD', &
             'dimension of deg_lat and top do not match', FATAL)
 if (size(deg_lat,2).ne.size(rh,2)) &
              call error_mesg( 'RH_CLOUDS in RH_CLOUDS_MOD', &
             'dimension of deg_lat and top do not match', FATAL)

cloud = .false.
nlev  = size(rh,3)

do k = 1, nlev
  if ( do_mcm_crit_rh ) then
    rh_crit = rh_crit_top + (rh_crit_bot - rh_crit_top)*p_full(:,:,k)/p_full(:,:,nlev)
  else
    rh_crit = rh_crit_top + (rh_crit_bot - rh_crit_top)*p_full(:,:,k)/p_surf
  endif
  where(rh(:,:,k) >= rh_crit) cloud(:,:,k) = .true.
enddo

if ( do_mcm_no_clouds_top ) cloud(:,:,1) = .false.

n_cloud = 0

do j = 1, size(rh,2)
  do i = 1, size(rh,1)
    if(cloud(i,j,1)) then
      n_cloud(i,j) = 1
      top(i,j,n_cloud(i,j)) = 1
    end if
    do k = 2, nlev
      if(.not.cloud(i,j,k).and.cloud(i,j,k-1)) then
         bot(i,j,n_cloud(i,j)) = k-1
      else if(cloud(i,j,k).and..not.cloud(i,j,k-1)) then
        n_cloud(i,j) = n_cloud(i,j) + 1
        top(i,j,n_cloud(i,j)) = k
      end if
    end do
    if(cloud(i,j,nlev)) bot(i,j,n_cloud(i,j)) = nlev
  end do
end do

max_n_cloud = maxval(n_cloud)
if(size(top,3).lt.max_n_cloud) call error_mesg( 'RH_CLOUDS in RH_CLOUDS_MOD',&
             'third dimension of top not large enough', FATAL)
call cloud_bounds(high_middle, middle_low, deg_lat)
call cloud_albedo(zenith, high_alb, middle_alb, low_alb)


abs_uv(:,:,:) = 0.0

do j = 1, size(top,2)
  do i = 1, size(top,1)
    do n = 1, n_cloud(i,j)


      !set cloud amount
      cldamt(i,j,n) = 1.

      sig_bot = p_full(i,j,bot(i,j,n))/p_surf(i,j)
      
      !guarantee some transmission to the clouds
      !by reducing the actual cloud reflectance in uv and nir band
      ! this break is necessary to avoid the rest of the
      ! radiation code from breaking up.
      if (sig_bot.le.high_middle(i,j)) then      ! high cloud   !
        alb_uv(i,j,n) = high_alb(i,j,1)
        alb_nir(i,j,n) = high_alb(i,j,2)
        abs_nir(i,j,n) = MIN(0.99-alb_nir(i,j,n),high_abs)
        emiss(i,j,n) = high_emiss

      elseif (sig_bot.gt.high_middle(i,j) .and. &
                sig_bot.le.middle_low(i,j))  then  ! middle cloud !
        alb_uv(i,j,n) = middle_alb(i,j,1)
        alb_nir(i,j,n) = middle_alb(i,j,2)
        abs_nir(i,j,n) = MIN(0.99-alb_nir(i,j,n),middle_abs)
        emiss(i,j,n) = middle_emiss

      elseif (sig_bot.gt.middle_low(i,j))  then  ! low cloud    !
        alb_uv(i,j,n) = low_alb(i,j,1)
        alb_nir(i,j,n) = low_alb(i,j,2)
        abs_nir(i,j,n) = MIN(0.99-alb_nir(i,j,n),low_abs)
        emiss(i,j,n) = low_emiss

      endif
        
    end do
  end do
end do

end subroutine rh_clouds_3d

!#######################################################################

subroutine cloud_albedo(zenith, high, middle, low)

real, intent(in), dimension(:,:)    :: zenith
real, intent(out), dimension(:,:,:) :: high, middle, low

integer, parameter :: num_angles = 17
integer, parameter :: num_bands  = 2

real, dimension(num_angles, 2) :: high_cloud, middle_cloud, low_cloud
real, dimension(num_angles, 2) :: low_cloud_tun

real, dimension(size(zenith,1),size(zenith,2)) :: z

real    :: pi, del, x, r
integer :: n, i, j, ind
integer :: n1, n2

! high cloud albedos for zenith angles from 0-80 deg. every 5 degs.
!    first for band =1, then band = 2

data high_cloud  &
 /.04,.05,.05,.05,.06,.06,.07,.07,.08,.11,.13,.16,.21,.28,.39,.48,.61, &
  .04,.05,.05,.05,.06,.06,.07,.07,.08,.10,.11,.14,.19,.26,.35,.44,.55/

! middle cloud albedos 

data middle_cloud &
 /.18,.18,.19,.20,.21,.23,.24,.26,.29,.33,.37,.42,.47,.55,.64,.71,.79, &
  .14,.14,.15,.16,.17,.18,.18,.20,.23,.25,.29,.32,.37,.43,.50,.55,.61/

! low cloud albedos 

data low_cloud &
 /.50,.50,.51,.51,.52,.53,.54,.56,.58,.62,.65,.67,.69,.73,.78,.82,.86, &
  .42,.42,.43,.43,.44,.45,.46,.48,.50,.52,.55,.57,.59,.63,.66,.70,.74/

pi = 4.0*atan(1.0)
z  = acos(zenith)*180.0/pi
del = 90.0/float(num_angles+1)

  do n1 = 1,num_angles
    do n2 = 1,num_bands
      low_cloud_tun(n1,n2) = tuning_coeff_low_cld*low_cloud(n1,n2)
    end do
  end do

! if zenith angle >= 80 degrees, use albedos for zenith angle = 80
do j = 1, size(zenith,2)
  do i = 1, size(zenith,1)
    if (z(i,j) .ge. 80.0) then
      do n = 1,num_bands
          high(i,j,n) =   high_cloud(num_angles,n)
        middle(i,j,n) = middle_cloud(num_angles,n)
           low(i,j,n) = low_cloud_tun(num_angles,n)
      end do
    else
      x = z(i,j)/del
      ind = floor(x)
      r = x - ind
      ind = ind + 1
      do n = 1,num_bands
          high(i,j,n) =   high_cloud(ind,n) &
                   + r*(  high_cloud(ind+1,n) -   high_cloud(ind,n))
        middle(i,j,n) = middle_cloud(ind,n) &
                   + r*(middle_cloud(ind+1,n) - middle_cloud(ind,n))
           low(i,j,n) =    low_cloud_tun(ind,n) &
                   + r*(   low_cloud_tun(ind+1,n) - low_cloud_tun(ind,n))
      end do
    end if
  end do
end do


end subroutine cloud_albedo

!#######################################################################

subroutine cloud_bounds(high_middle, middle_low, deg_lat)

real, intent(in) , dimension(:,:) :: deg_lat
real, intent(out), dimension(:,:) :: high_middle, middle_low

real,dimension(size(deg_lat,1),size(deg_lat,2)) :: x

   x = (90.0 - abs(deg_lat))/90.

   high_middle = high_middle_pole + x*(high_middle_eq - high_middle_pole)
   middle_low  = middle_low_pole  + x*(middle_low_eq  - middle_low_pole )

return
end subroutine cloud_bounds


!#######################################################################
!  THE FOLLOWING CODE ALLOWS RH_CLOUDS TO BE USED IN 2D AND 1D MODELS
!#######################################################################

subroutine rh_clouds_2d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:,:)   :: rh, p_full
real   , intent(in) , dimension(:)     :: p_surf,zenith,deg_lat
integer, intent(out), dimension(:,:)   :: top, bot
integer, intent(out), dimension(:)     :: n_cloud
real   , intent(out), dimension(:,:)   :: cldamt,alb_uv,alb_nir,abs_uv
real   , intent(out), dimension(:,:)   :: abs_nir,emiss

real, dimension(1, size(rh,1),size(rh,2)) :: rh_3d, p_full_3d
real, dimension(1, size(rh,1)           ) :: p_surf_3d, zenith_3d,deg_lat_3d
real, dimension(1, size(rh,1),size(rh,2)) :: cldamt_3d,alb_uv_3d,alb_nir_3d
real, dimension(1, size(rh,1),size(rh,2)) :: abs_uv_3d,abs_nir_3d,emiss_3d
integer, dimension(1, size(rh,1), size(rh,2)) :: top_3d, bot_3d
integer, dimension(1, size(rh,1)        ) :: n_cloud_3d

!assign variables to 3d matrices
rh_3d(1,:,:) = rh
p_full_3d(1,:,:) = p_full
p_surf_3d(1,:) = p_surf
zenith_3d(1,:) = zenith
deg_lat_3d(1,:) = deg_lat

call rh_clouds_3d(rh_3d, p_full_3d, p_surf_3d, zenith_3d, deg_lat_3d,&
   n_cloud_3d,top_3d,bot_3d,cldamt_3d,alb_uv_3d,alb_nir_3d,abs_uv_3d,&
   abs_nir_3d,emiss_3d)

!patch back results to 2d matrices
n_cloud = n_cloud_3d(1,:)
top = top_3d(1,:,:)
bot = bot_3d(1,:,:)
cldamt = cldamt_3d(1,:,:)
alb_uv = alb_uv_3d(1,:,:)
alb_nir = alb_nir_3d(1,:,:)
abs_uv = abs_uv_3d(1,:,:)
abs_nir = abs_nir_3d(1,:,:)
emiss = emiss_3d(1,:,:)

end subroutine rh_clouds_2d

!#######################################################################

subroutine rh_clouds_1d(rh, p_full, p_surf, zenith, deg_lat,&
            n_cloud,top,bot,cldamt,alb_uv,alb_nir,abs_uv,abs_nir,emiss)

real   , intent(in) , dimension(:)   :: rh, p_full
real   , intent(in)                  :: p_surf,zenith,deg_lat
integer, intent(out), dimension(:)   :: top, bot
integer, intent(out)                 :: n_cloud
real   , intent(out), dimension(:)   :: cldamt,alb_uv,alb_nir,abs_uv
real   , intent(out), dimension(:)   :: abs_nir,emiss

real, dimension(1, 1,size(rh(:))) :: rh_3d, p_full_3d
real, dimension(1, 1            ) :: p_surf_3d, zenith_3d,deg_lat_3d
real, dimension(1, 1,size(rh(:))) :: cldamt_3d,alb_uv_3d,alb_nir_3d
real, dimension(1, 1,size(rh(:))) :: abs_uv_3d,abs_nir_3d,emiss_3d
integer, dimension(1, 1, size(rh(:))) :: top_3d, bot_3d
integer, dimension(1, 1             ) :: n_cloud_3d

!assign variables to 3d matrices
rh_3d(1,1,:) = rh
p_full_3d(1,1,:) = p_full
p_surf_3d(1,1) = p_surf
zenith_3d(1,1) = zenith
deg_lat_3d(1,1) = deg_lat

call rh_clouds_3d(rh_3d, p_full_3d, p_surf_3d, zenith_3d, deg_lat_3d,&
   n_cloud_3d,top_3d,bot_3d,cldamt_3d,alb_uv_3d,alb_nir_3d,abs_uv_3d,&
   abs_nir_3d,emiss_3d)

!patch back results to 2d matrices
n_cloud = n_cloud_3d(1,1)
top = top_3d(1,1,:)
bot = bot_3d(1,1,:)
cldamt = cldamt_3d(1,1,:)
alb_uv = alb_uv_3d(1,1,:)
alb_nir = alb_nir_3d(1,1,:)
abs_uv = abs_uv_3d(1,1,:)
abs_nir = abs_nir_3d(1,1,:)
emiss = emiss_3d(1,1,:)

end subroutine rh_clouds_1d


!#######################################################################



end module rh_clouds_mod


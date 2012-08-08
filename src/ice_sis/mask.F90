!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! mask_mod reads the land/sea mask from netCDF file - Mike Winton (Michael.Winton)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module mask_mod

use rot_mod, only: rot_intrp_sclr
use fms_mod,   only: error_mesg, FATAL

implicit none
include 'netcdf.inc'
private

public :: lonb, latb, mask, lonb_uv, latb_uv, mask_uv, make_mask

real,    allocatable, dimension(:  ) :: lonb, latb, lonb_uv, latb_uv
logical, allocatable, dimension(:,:) :: mask, mask_uv
real,    allocatable, dimension(:,:) :: rmask

contains

!
! make_v_grid - makes B-grid longitudes, latitudes and mask for velocity-cells
!
subroutine make_v_grid(mask_t, mask_v, lon_t, lat_t, lon_v, lat_v)
  logical, dimension (:,:),                           intent(in   ) :: mask_t
  logical, dimension (size(mask_t,1),size(mask_t,2))                :: mask_v
  real,    dimension (size(mask_t,1)+1),              intent(in   ) :: lon_t
  real,    dimension (size(mask_t,2)+1),              intent(in   ) :: lat_t
  real,    dimension (size(lon_t(:)))                                  :: lon_v
  real,    dimension (size(lat_t(:)))                                  :: lat_v

  integer, dimension(size(mask_t,1)) :: wes
  integer, dimension(size(mask_t,2)) :: sou
  integer                               i,j
  real,    dimension(size(mask_t,1)) :: lon_tmp
  real                                  pi

  sou = (/(j,j=1,size(sou(:)))/); wes = (/(i,i=1,size(wes(:)))/)
  sou = cshift(sou, -1)      ; wes = cshift(wes, -1)

  mask_v = mask_t(:,:)   .and. mask_t(wes,:) .and. &
           mask_t(:,sou) .and. mask_t(wes,sou)

  pi      = 4*atan(1.0)
  lon_tmp  = lon_t(1:size(lon_tmp(:)))
  lon_tmp  = lon_tmp-0.5*mod(lon_tmp+2*pi-cshift(lon_tmp,-1),2*pi)
  lon_v(1:size(lon_tmp(:))) = lon_tmp
  lon_v(size(lon_v(:)))     = lon_v(1)

  lat_v    = (lat_t+cshift(lat_t,-1))/2
  lat_v(1) = lat_t(1)              ! southernmost velocities not computed anyway
end subroutine make_v_grid

subroutine make_mask(lon_start, lon_end, lon_incr, lat_start, lat_end, lat_incr)
real :: lon_start, lon_end, lon_incr, lat_start, lat_end, lat_incr
integer      :: dims(4), start(4), count(4), rcode
integer      :: i, j, id, jd, im, jm 
integer      :: ncid, varid, xv_id, yv_id
real         :: pi
real, dimension(:  ), allocatable :: xd, yd, xt, yt
real, dimension(:,:), allocatable :: map

   pi = 4*atan(1.0)

   rcode = nf_open('INPUT/map.nc',0,ncid)
   if (rcode/=0) call error_mesg ('mask_mod', &
                                  'cannot open INPUT/map.nc', FATAL)
   rcode = nf_inq_varid(ncid, 'MAP', varid)
   rcode = nf_inq_vardimid(ncid, varid, dims)
   rcode = nf_inq_dimlen(ncid, dims(1), id)
   rcode = nf_inq_dimlen(ncid, dims(2), jd)

   allocate(xd(id), yd(jd), map(id,jd))

   start = 1; count = 1; count(1) = id;
   rcode = nf_get_vara_double(ncid, dims(1), start, count, xd)
   start = 1; count = 1; count(1) = jd;
   rcode = nf_get_vara_double(ncid, dims(2), start, count, yd)
   start = 1; count = 1; count(1) = id; count(2) = jd;
   rcode = nf_get_vara_double(ncid, varid, start, count, map)

   im = nint((lon_end-lon_start)/lon_incr)
   jm = nint((lat_end-lat_start)/lat_incr)
   allocate(rmask(im,jm), mask(im,jm), mask_uv(im,jm), xt(im), yt(jm) )
   allocate(lonb(im+1), lonb_uv(im+1), latb(jm+1), latb_uv(jm+1))
  
   xt = lon_start+(/(i-0.5,i=1,im)/)*lon_incr
   yt = lat_start+(/(j-0.5,j=1,jm)/)*lat_incr

   if (maxval(abs(yt))<45.0) then ! rotate poles onto equator
     call rot_intrp_sclr(map, xd, yd, id, jd, rmask, xt, yt, size(xt(:)),  &
                         size(yt(:)), 2*atan(1.0), -2*atan(1.0), atan(1.0) )
   else
     call rot_intrp_sclr(map, xd, yd, id, jd, rmask, xt, yt, size(xt(:)), &
                         size(yt(:)), 0.0, 0.0, 0.0                       )
   endif
   deallocate(map, xd, yd, xt, yt)

   lonb = (lon_start+(/(i,i=0,im)/)*lon_incr)*pi/180
   latb = (lat_start+(/(j,j=0,jm)/)*lat_incr)*pi/180
   mask = rmask > 0.5

   call make_v_grid(mask, mask_uv, lonb, latb, lonb_uv, latb_uv)

end subroutine make_mask

end module mask_mod

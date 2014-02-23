module land_io_mod

use constants_mod,     only : PI
use fms_mod,           only : file_exist, error_mesg, FATAL, stdlog, mpp_pe, &
     mpp_root_pe, write_version_number, string, check_nml_error, close_file

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif


use horiz_interp_mod,  only : horiz_interp_type, &
     horiz_interp_new, horiz_interp_del, &
     horiz_interp

use land_numerics_mod, only : nearest, bisect
use nf_utils_mod,      only : nfu_validtype, nfu_get_dim, nfu_get_dim_bounds, &
     nfu_get_valid_range, nfu_is_valid, nfu_inq_var, nfu_get_var

implicit none
private

! ==== public interface ======================================================
public :: init_cover_field
public :: read_field
public :: read_land_io_namelist

public :: print_netcdf_error

public :: input_buf_size
! ==== end of public interface ===============================================

interface read_field
   module procedure read_field_N_2D, read_field_N_3D
   module procedure read_field_I_2D, read_field_I_3D
end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_io_mod', &
     version     = '$Id: land_io.F90,v 20.0 2013/12/13 23:29:51 fms Exp $', &
     tagname     = '$Name: tikal $'

logical :: module_is_initialized = .false.

character(len=64)  :: interp_method = "conservative"
integer :: input_buf_size = 65536 ! input buffer size for tile and cohort reading
namelist /land_io_nml/ interp_method, input_buf_size

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_land_io_namelist()
  integer :: io, ierr, unit


  module_is_initialized = .TRUE.

  ! [1] print out version number
  call write_version_number (version, tagname)

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_io_nml, iostat=io)
     ierr = check_nml_error(io, 'land_io_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_io_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_io_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif   
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_io_nml)
     call close_file (unit)
  endif

  if(trim(interp_method) .NE. "conservative" .AND. trim(interp_method) .NE. "conserve_great_circle") then
     call error_mesg ( module_name,'interp_method should be "conservative" or "conserve_great_circle"', FATAL)
  endif
  
  if (input_buf_size <= 0) then
     call error_mesg ( module_name,'input_buf_size must be larger than zero', FATAL)
  endif

end subroutine read_land_io_namelist


! ============================================================================
! This procedure creates and initializes a field of fractional coverage.
subroutine init_cover_field( &
     cover_to_use, filename, cover_field_name, frac_field_name, &
     lonb, latb, uniform_cover, input_cover_types, frac)
  character(len=*), intent(in) :: cover_to_use
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: cover_field_name, frac_field_name
  real            , intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  integer         , intent(in) :: uniform_cover
  integer         , intent(in) :: input_cover_types(:)
  real            , intent(out):: frac(:,:,:) ! output-global map of soil fractional coverage

  ! ---- local vars ---------------------------------------------------------
  integer :: i,j,k     ! iterators
  integer :: cover_id
  real    :: maxfrac, total

  if( .not. module_is_initialized ) &
       call error_mesg(module_name,'land_io_init is not called', FATAL)
 
  frac = 0
  
  if (cover_to_use == 'multi-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
  else if (cover_to_use=='single-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do j = 1,size(frac,2)
     do i = 1,size(frac,1)
        total = sum(frac(i,j,:))
        if (total <= 0) cycle ! operate on valid input data points only
        maxfrac=0 ; cover_id=1
        do k = 1,size(frac,3)
           if(frac(i,j,k).gt.maxfrac) then
              maxfrac=frac(i,j,k)
              cover_id=k
           endif
        enddo
        ! set all fractions except dominant fraction to zero
        frac(i,j,:) = 0.0
        frac(i,j,cover_id) = total
     enddo
     enddo
  else if (cover_to_use == 'uniform') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do j = 1,size(frac,2)
     do i = 1,size(frac,1)
        total = sum(frac(i,j,:))
        if (total <= 0) cycle ! operate on valid input data points only
        ! set all fractions except dominant fraction to zero
        frac(i,j,:) = 0.0
        frac(i,j,uniform_cover) = total
     enddo
     enddo
  else
     call error_mesg ( module_name,'illegal value of cover_to_use '//cover_to_use, FATAL )
  endif

end subroutine init_cover_field


! ============================================================================
subroutine read_cover_field(file, cover_field_name, frac_field_name,&
     lonb, latb, input_cover_types, frac)
  character(len=*)  , intent(in)  :: file            ! file to read from
  character(len=*)  , intent(in)  :: cover_field_name, frac_field_name
  real              , intent(in)  :: lonb(:,:),latb(:,:) ! boundaries of the model grid
  real              , intent(out) :: frac(:,:,:)     ! resulting fractions
  integer, optional , intent(in)  :: input_cover_types(:)

  ! --- local vars
  integer :: ncid, varid

  if (.not.file_exist(file)) &
       call error_mesg(module_name,'input file "'//trim(file)//'" does not exist',FATAL)

  __NF_ASRT__( nf_open(file, NF_NOWRITE, ncid) )
  if(nf_inq_varid(ncid,cover_field_name,varid)==NF_NOERR) then
     call do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
  else if ( nf_inq_varid(ncid,frac_field_name,varid)==NF_NOERR) then
     call do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)
  else
     call error_mesg(module_name,&
          'neither "'//trim(cover_field_name)//'" nor "'//&
          frac_field_name//'" is present in input file "'//trim(file)//'"' ,&
          FATAL)
  endif
  __NF_ASRT__( nf_close(ncid) )

end subroutine read_cover_field

! ============================================================================
subroutine do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: ncid, varid
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat ! size of input map
  integer :: k
  integer, allocatable :: in_cover(:,:)
  real, allocatable    :: in_lonb(:), in_latb(:), x(:,:)
  type(horiz_interp_type) :: interp
  integer :: vardims(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v
  integer :: in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: iret ! result of netcdf calls

  ! find out dimensions, etc
  __NF_ASRT__( nf_inq_vardimid(ncid,varid,vardims) )
  ! get size of the longitude and latitude axes
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(1), nlon) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(2), nlat) )
  allocate ( in_lonb (nlon+1), in_latb (nlat+1) )
  __NF_ASRT__( nfu_get_dim_bounds(ncid, vardims(1), in_lonb) )
  __NF_ASRT__( nfu_get_dim_bounds(ncid, vardims(2), in_latb) )
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  ! to minimize the i/o and work done by horiz_interp, find the boundaries
  ! of latitude belt in input data that covers the entire latb array
  in_j_start=bisect(in_latb, minval(latb))
  in_j_count=bisect(in_latb, maxval(latb))-in_j_start+1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_cover_field','input latitude start index ('&
                     //string(in_j_start)//') is out of bounds', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_cover_field','input latitude count ('&
                     //string(in_j_count)//') is too large (start index='&
                     //string(in_j_start)//')', FATAL)

  ! allocate input data buffers
  allocate ( x(nlon,in_j_count), in_cover(nlon,in_j_count) )

  ! read input data
  iret = nf_get_vara_int(ncid,varid, &
                    (/1,in_j_start/), (/nlon,in_j_count/), in_cover)
  __NF_ASRT__( iret )
  __NF_ASRT__( nfu_get_valid_range(ncid,varid,v) )

  call horiz_interp_new(interp, in_lonb,in_latb(in_j_start:in_j_start+in_j_count), &
       lonb,latb, interp_method=trim(interp_method))
  frac=0
  do k = 1,size(input_cover_types(:))
     x=0
     where(nfu_is_valid(in_cover,v).and.in_cover==input_cover_types(k)) x = 1

     call horiz_interp(interp,x,frac(:,:,k))
  enddo

  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine do_read_cover_field


! ============================================================================
subroutine do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: ncid, varid
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:,:)

  ! ---- local vars
  integer :: nlon, nlat, ntypes, k, cover
  real, allocatable :: in_frac(:,:,:)
  real, allocatable :: in_lonb(:), in_latb(:)
  real, allocatable :: in_mask(:,:)
  type(horiz_interp_type) :: interp
  type(nfu_validtype) :: v
  integer :: vardims(NF_MAX_VAR_DIMS)
  integer :: in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: iret ! result of netcdf calls

  ! find out dimensions, etc
  __NF_ASRT__( nf_inq_vardimid(ncid,varid,vardims) )
  ! get size of the longitude and latitude axes
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(1), nlon) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(2), nlat) )
  __NF_ASRT__( nf_inq_dimlen(ncid, vardims(3), ntypes))
  allocate ( in_lonb(nlon+1), in_latb(nlat+1) )
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(1), in_lonb))
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(2), in_latb))
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  ! find the boundaries of latitude belt in input data that covers the 
  ! entire latb array
  in_j_start=bisect(in_latb, minval(latb))
  in_j_count=bisect(in_latb, maxval(latb))-in_j_start+1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_fraction_field','input latitude start index ('&
                     //string(in_j_start)//') is out of bounds', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_fraction_field','input latitude count ('&
                     //string(in_j_count)//') is too large (start index='&
                     //string(in_j_start)//')', FATAL)

  allocate( in_mask(nlon,in_j_count), in_frac(nlon,in_j_count,ntypes) )

  ! read input data
  iret = nf_get_vara_double(ncid, varid, &
          (/1,in_j_start,1/), (/nlon,in_j_count,ntypes/), in_frac)
  __NF_ASRT__( iret ) 
  __NF_ASRT__( nfu_get_valid_range(ncid,varid,v) )

  frac = 0
  do k = 1,size(input_cover_types)
     cover = input_cover_types(k)
     if (cover<1.or.cover>ntypes) then
        cycle ! skip all invalid indices in the array of input cover types
     endif
     in_mask = 0.0
     where(nfu_is_valid(in_frac(:,:,cover),v)) in_mask = 1.0
     call horiz_interp_new(interp, &
          in_lonb,in_latb(in_j_start:in_j_start+in_j_count), lonb,latb,&
          interp_method=trim(interp_method), mask_in=in_mask)
     call horiz_interp(interp,in_frac(:,:,cover),frac(:,:,k))
     call horiz_interp_del(interp)
  enddo

  ! clean up memory
  deallocate(in_lonb, in_latb, in_frac, in_mask)

end subroutine do_read_fraction_field


! ============================================================================
subroutine read_field_N_2D(filename, varname, lon, lat, data, interp)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional :: interp

  ! ---- local vars ----------------------------------------------------------
  real    :: data3(size(data,1),size(data,2),1)

  call read_field_N_3D(filename, varname, lon, lat, data3, interp)
  data = data3(:,:,1)

end subroutine read_field_N_2D

! ============================================================================
subroutine read_field_N_3D(filename, varname, lon, lat, data, interp)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  real, intent(in)  :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional :: interp

  ! ---- local vars ----------------------------------------------------------
  integer :: ncid
  integer :: iret

  iret = nf_open(filename,NF_NOWRITE,ncid)
  if(iret/=NF_NOERR) then
     call error_mesg('read_field','Can''t open netcdf file "'//trim(filename)//'"',FATAL)
  endif
  call read_field_I_3D(ncid, varname, lon, lat, data, interp)
  __NF_ASRT__( nf_close(ncid) )

end subroutine read_field_N_3D

! ============================================================================
subroutine read_field_I_2D(ncid, varname, lon, lat, data, interp)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:)
  character(len=*), intent(in), optional  :: interp
  ! ---- local vars
  real    :: data3(size(data,1),size(data,2),1)

  call read_field_I_3D(ncid, varname, lon, lat, data3, interp)
  data = data3(:,:,1)

end subroutine read_field_I_2D

! ============================================================================
subroutine read_field_I_3D(ncid, varname, lon, lat, data, interp)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real, intent(in) :: lon(:,:),lat(:,:)
  real, intent(out) :: data(:,:,:)
  character(len=*), intent(in), optional  :: interp

  ! ---- local vars ----------------------------------------------------------
  integer :: nlon, nlat, nlev ! size of input grid
  integer :: varndims ! number of variable dimension
  integer :: vardims(NF_MAX_VAR_DIMS) ! IDs of variable dimension
  integer :: dimlens(NF_MAX_VAR_DIMS) ! sizes of respective dimensions
  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:)
  real,    allocatable :: x(:,:,:) ! input buffer
  logical, allocatable :: mask(:,:,:) ! mask of valid values
  real,    allocatable :: rmask(:,:,:) ! real mask for interpolator
  character(len=20) :: interpolation 
  integer :: i,j,k,imap,jmap !
  type(nfu_validtype) :: v

  interpolation = "bilinear"
  if(present(interp)) interpolation = interp
  
  ! get the dimensions of our variable
  __NF_ASRT__( nfu_inq_var(ncid,varname,ndims=varndims,dimids=vardims,dimlens=dimlens) )
  if(varndims<2.or.varndims>3) then
     call error_mesg('read_field','variable "'//trim(varname)//'" is '//string(varndims)//&
          'D, but only reading 2D or 3D variables is supported', FATAL)
  endif
  nlon = dimlens(1) ; nlat = dimlens(2)
  nlev = 1; 
  if (varndims==3) nlev=dimlens(3)
  if(nlev/=size(data,3)) then
     call error_mesg('read_field','3rd dimension length of the variable "'&
          //trim(varname)//'" ('//trim(string(nlev))//') is different from the expected size of data ('// &
          trim(string(size(data,3)))//')', FATAL)
  endif

  allocate (                 &
       in_lon  (nlon),   in_lat  (nlat),   &
       in_lonb (nlon+1), in_latb (nlat+1), &
       x       (nlon, nlat, nlev) ,&
       mask    (nlon, nlat, nlev) , rmask(nlon, nlat, nlev) )

  ! read boundaries of the grid cells in longitudinal direction
  __NF_ASRT__(nfu_get_dim(ncid, vardims(1), in_lon))
  __NF_ASRT__(nfu_get_dim(ncid, vardims(2), in_lat))
  in_lon = in_lon*PI/180.0; in_lat = in_lat*PI/180.0
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(1), in_lonb))
  __NF_ASRT__(nfu_get_dim_bounds(ncid, vardims(2), in_latb))
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
  __NF_ASRT__(nfu_get_valid_range(ncid,varname,v))
  ! read input data
  __NF_ASRT__( nfu_get_var(ncid,varname,x) ) ! assuming real is real*8
  mask = nfu_is_valid(x,v)
  rmask = 1.0
  where(.not.mask) rmask = 0.0

  select case(trim(interpolation))
  case ("bilinear")
     do k = 1,size(data,3)
        call horiz_interp(x(:,:,k), in_lonb, in_latb, lon,lat, data(:,:,k), mask_in=rmask(:,:,k), &
             interp_method='bilinear')
     enddo
  case ("nearest")
     do k = 1,size(data,3)
     do j = 1,size(data,2)
     do i = 1,size(data,1)
        call nearest (mask(:,:,k), in_lon, in_lat, lon(i,j), lat(i,j), imap, jmap)
        data(i,j,k) = x(imap,jmap,k)
     enddo
     enddo
     enddo
  case default
     call error_mesg(module_name, interpolation//" is not a valid interpolation method",FATAL)
  end select

  deallocate(in_lonb, in_latb, in_lon, in_lat, x, mask, rmask)

end subroutine read_field_I_3D

! ============================================================================
subroutine print_netcdf_error(ierr, file, line)
  ! prints out NetCDF library error message, including file name and line number
  integer,          intent(in) :: ierr ! error code
  character(len=*), intent(in) :: file ! name of the file
  integer,          intent(in) :: line ! number of line in the file

  ! ---- local vars
  character(len=1024) :: mesg

  if (ierr.ne.NF_NOERR) then
     write(mesg, "('File ',a,' Line ',i4.4,' :: ',a)") &
          trim(file),line,trim(NF_STRERROR(ierr))
     call error_mesg('NetCDF', mesg, FATAL)
  endif
end subroutine print_netcdf_error

end module

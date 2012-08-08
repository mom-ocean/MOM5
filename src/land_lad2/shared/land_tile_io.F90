module land_tile_io_mod

use mpp_mod, only : mpp_send, mpp_recv, mpp_sync
use mpp_mod, only : COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
use mpp_mod, only : COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
use fms_mod, only : error_mesg, FATAL, mpp_pe, get_mosaic_tile_file
use fms_io_mod, only : get_instance_filename
use time_manager_mod, only : time_type
use data_override_mod, only : data_override

use nf_utils_mod, only : nfu_inq_dim, nfu_inq_var, nfu_def_dim, nfu_def_var, &
     nfu_get_var, nfu_put_att
use land_io_mod, only : print_netcdf_error, read_field
use land_tile_mod, only : land_tile_type, land_tile_list_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), &
     get_elmt_indices
use land_data_mod, only  : lnd
use land_utils_mod, only : put_to_tiles_r0d_fptr

implicit none
private

! ==== public interfaces =====================================================
! restart i/o subroutines: those use CF "compression by gathering" technique
! to pack tile data.
public :: create_tile_out_file
public :: read_tile_data_r0d_fptr,  read_tile_data_r1d_fptr
public :: read_tile_data_i0d_fptr
public :: write_tile_data_r0d_fptr, write_tile_data_r1d_fptr
public :: write_tile_data_i0d_fptr

! data override subroutines
public :: override_tile_data_r0d_fptr

! auxiliary subroutines
public :: get_tile_by_idx
public :: print_netcdf_error

public :: read_field

public :: get_input_restart_name

public :: sync_nc_files ! synchronizes writer and reader processors
! ==== end of public interfaces ==============================================
interface create_tile_out_file
   module procedure create_tile_out_file_idx
   module procedure create_tile_out_file_fptr
end interface

interface read_tile_data_r1d_fptr
   module procedure read_tile_data_r1d_fptr_all
   module procedure read_tile_data_r1d_fptr_idx
end interface

interface write_tile_data_r1d_fptr
   module procedure write_tile_data_r1d_fptr_all
   module procedure write_tile_data_r1d_fptr_idx
end interface

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_tile_io_mod', &
     version     = '$Id: land_tile_io.F90,v 19.0.6.2 2012/05/14 19:16:06 Zhi.Liang Exp $', &
     tagname     = '$Name: siena_201207 $'

! name of the "compressed" dimension (and dimension variable) in the output 
! netcdf files -- that is, the dimensions written out using compression by 
! gathering, as described in CF conventions. See subroutines write_tile_data,
! read_tile_data, read_unpack_tile_data, write_cohort_data
character(len=*),   parameter :: tile_index_name   = 'tile_index'
integer, parameter :: INPUT_BUF_SIZE=1024 ! size of the input buffer for tile input 


! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! =============================================================================
! given a generic name of the restart file, checks if a file with one of the 
! possible restarts file names exists, and if it does returns the tile-qualified 
! (or tile- and processor-qualified) name of the restart.
subroutine get_input_restart_name(name, restart_exists, actual_name)
  character(*), intent(in)  :: name        ! "generic" name of the restart
  logical     , intent(out) :: restart_exists ! TRUE if any file found
  character(*), intent(out) :: actual_name ! name of the found file, if any

  ! ---- local vars
  character(6) :: PE_suffix ! PE number
  
  ! Build the restart file name.
  call get_instance_filename(trim(name), actual_name)
  call get_mosaic_tile_file(trim(actual_name),actual_name,.false.,lnd%domain)
  ! we can't use fms file_exist function here, because it lies: it checks not
  ! just the original name, but the name with PE suffix, and returns true if
  ! either of those exist
  inquire (file=trim(actual_name), exist=restart_exists)
  if (.not.restart_exists) then
     ! try the name with current PE number attached
     write(PE_suffix,'(".",I4.4)') lnd%io_id
     actual_name = trim(actual_name)//trim(PE_suffix)
     inquire (file=trim(actual_name), exist=restart_exists)
  endif

end subroutine

! =============================================================================
! this subroutine creates netcdf file for output of tiled data using "compression
! by gathering," as described in CF conventions, and creates coordinate system
! necessary for write_tile_data subroutines.
! In particular:
!  "compressed" dimension and integer variable with appropriate attributes, with 
! variable filled with packing indices
!   horizontal dimensions "lat" and "lon" with associated variable, and boundaries, 
! describing global grid
!   dimension "tile," without any associated variable. length of this dimension is
! equal to current global max of number of tiles per grid cell
!
! The file is actually created only by root processor of our io_domain; the rest 
! of the processors just open the created file in NOWRITE mode. 
subroutine create_tile_out_file_idx(ncid, name, glon, glat, tidx, tile_dim_length, reserve)
  integer          , intent(out) :: ncid      ! resulting NetCDF id
  character(len=*) , intent(in)  :: name      ! name of the file to create
  real             , intent(in)  :: glon(:)   ! longitudes of the grid centers
  real             , intent(in)  :: glat(:)   ! latitudes of the grid centers
  integer          , intent(in)  :: tidx(:)   ! integer compressed index of tiles
  integer          , intent(in)  :: tile_dim_length ! length of tile axis
  integer, optional, intent(in)  :: reserve   ! amount of space to reserve for 
  ! header expansion. This subroutine and following calls to write_tile_data
  ! will work even if this is set to 0, but for efficiency it is useful to
  ! specify some non-zero value, so that netcdf library do not have to rewrite
  ! entire file each time a new variable is added. Default value (8K) should do
  ! fine in most cases.

  ! ---- local vars
  integer        :: reserve_  ! local value of space to reserve at the end of NetCDF header
  character(256) :: full_name ! full name of the file, including the processor number
  character(6)   :: PE_suffix ! PE number
  integer, allocatable :: ntiles(:) ! list of land tile numbers for each of PEs in io_domain
  integer, allocatable :: tidx2(:)  ! array of tile indices from all PEs in io_domain
  integer :: p ! io_domain PE iterator 
  integer :: k ! current index in tidx2 array for receive operation

  ! form the full name of the file
  call get_instance_filename(trim(name), full_name)
  call get_mosaic_tile_file(trim(full_name),full_name,.false.,lnd%domain)
  if (lnd%append_io_id) then
      write(PE_suffix,'(".",I4.4)') lnd%io_id
  else
      PE_suffix = ''
  endif

  full_name = trim(full_name)//trim(PE_suffix)

  if(tile_dim_length<=0) &
    call error_mesg('create_tile_out_file','tile axis length must be positive', FATAL)

  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if current PE doesn't do io, we just send the data to the processor that
     ! does
     call mpp_send(size(tidx), plen=1,          to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
     call mpp_send(tidx(1),    plen=size(tidx), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
  else
     ! gather an array of tile sizes from all processors in our io_domain
     allocate(ntiles(size(lnd%io_pelist)))
     ntiles(1) = size(tidx)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ntiles(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
     ! gather tile indices from all processors in our io_domain
     allocate(tidx2(sum(ntiles(:))))
     tidx2(1:ntiles(1))=tidx(:)
     k=ntiles(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(tidx2(k), from_pe=lnd%io_pelist(p), glen=ntiles(p), tag=COMM_TAG_2)
        k = k+ntiles(p)
     enddo

     ! create netcdf file
#ifdef use_netCDF3
     __NF_ASRT__(nf_create(full_name,NF_CLOBBER,ncid))
#elif use_LARGEFILE
     __NF_ASRT__(nf_create(full_name,ior(NF_64BIT_OFFSET,NF_CLOBBER),ncid))
#else
     __NF_ASRT__(nf_create(full_name,ior(NF_NETCDF4,NF_CLASSIC_MODEL),ncid))
#endif

     ! create lon, lat, dimensions and variables
     __NF_ASRT__(nfu_def_dim(ncid,'lon' ,glon(:) ,'longitude','degrees_east'))
     __NF_ASRT__(nfu_def_dim(ncid,'lat' ,glat(:) ,'latitude','degrees_north'))

     __NF_ASRT__(nfu_def_dim(ncid,'tile',tile_dim_length))
     ! the size of tile dimension really does not matter for the output, but it does
     ! matter for uncompressing utility, since it uses it as a size of the array to
     ! unpack to
     ! create tile index dimension and variable
     __NF_ASRT__(nfu_def_dim(ncid,tile_index_name,tidx2,'compressed land point index'))
     __NF_ASRT__(nfu_put_att(ncid,tile_index_name,'compress','tile lat lon'))
     __NF_ASRT__(nfu_put_att(ncid,tile_index_name,'valid_min',0))
     ! release the data we no longer need
     deallocate(ntiles,tidx2)

     ! determine the local value of space reserved in the header; by default 16K
     reserve_ = 1024*16
     if(present(reserve)) reserve_ = reserve

     ! end definition mode, reserving some space for future additions
     ! this call also commits the changes to the disk
     __NF_ASRT__(nf__enddef(ncid,reserve_,4,0,4))
     ! arguments are ncid,h_minfree,v_align,v_minfree,r_align; default is (ncid,0,4,0,4).
     ! The above call reserves some space at the end of the netcdf header for
     ! future expansion without library's having to rewrite the entire file. See 
     ! manual pages netcdf(3f) or netcdf(3) for more information.
  endif
  ! make sure send-receive operations and file creation have finished
  call mpp_sync()
  ! open file on non-writing processors to have access to the tile index
  if (mpp_pe()/=lnd%io_pelist(1)) then
     __NF_ASRT__(nf_open(full_name,NF_NOWRITE,ncid))
  endif
end subroutine create_tile_out_file_idx

! =============================================================================
subroutine create_tile_out_file_fptr(ncid, name, glon, glat, tile_exists, &
     tile_dim_length, reserve, created)
  integer          , intent(out) :: ncid      ! resulting NetCDF id
  character(len=*) , intent(in)  :: name      ! name of the file to create
  real             , intent(in)  :: glon(:)   ! longitudes of the grid centers
  real             , intent(in)  :: glat(:)   ! latitudes of the grid centers
  integer          , intent(in)  :: tile_dim_length ! length of tile axis
  integer, optional, intent(in)  :: reserve   ! amount of space to reserve for
  logical, optional, intent(out) :: created   ! indicates wether the file was 
      ! created; it is set to false if no restart needs to be written, in case 
      ! the total number of qualifying tiles in this domain is equal to zero
  ! the following interface describes the "detector function", which is passed 
  ! through the argument list and must return true for any tile to be written 
  ! to the specific restart, false otherwise
  interface
     logical function tile_exists(tile)
        use land_tile_mod, only : land_tile_type
        type(land_tile_type), pointer :: tile
     end function tile_exists
  end interface

  ! ---- local vars
  type(land_tile_enum_type) :: ce, te ! tile list elements
  integer, allocatable :: idx(:)   ! integer compressed index of tiles
  integer :: i,j,k,n

  ! count total number of tiles in this domain
  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
  te = tail_elmt (lnd%tile_map)
  n  = 0
  do while (ce/=te)
     if (tile_exists(current_tile(ce))) n = n+1
     ce=next_elmt(ce)
  end do
  
  ! calculate compressed tile index to be written to the restart file;
  allocate(idx(max(n,1))); idx(:) = -1 ! set init value to a known invalid index
  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
  n = 1
  do while (ce/=te)
     call get_elmt_indices(ce,i,j,k)
     
     if(tile_exists(current_tile(ce))) then
        idx (n) = (k-1)*lnd%nlon*lnd%nlat + (j-1)*lnd%nlon + (i-1)
        n = n+1
     endif
     ce=next_elmt(ce)
  end do
  ! create tile output file, defining horizontal coordinate and compressed
  ! dimension
  call create_tile_out_file_idx(ncid, name, glon, glat, idx, tile_dim_length, reserve)

  if (present(created)) created = .true.
  
end subroutine

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array returns a pointer to the tile
! corresponding to the compressed index, or NULL is the index is outside 
! current domain, or such tile does not exist.
subroutine get_tile_by_idx(idx,nlon,nlat,tiles,is,js,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: nlon, nlat
   integer, intent(in) :: is, js
   type(land_tile_list_type), intent(in) :: tiles(is:,js:)
   type(land_tile_type), pointer :: ptr

   ! ---- local vars
   integer :: i,j,k
   type(land_tile_enum_type) :: ce,te
   
   ptr=>null()

   if(idx<0) return ! negative indices do not correspond to any tile

   ! given tile idx, calculate global lon, lat, and tile indices
   k = idx
   i = modulo(k,nlon)+1 ; k = k/nlon
   j = modulo(k,nlat)+1 ; k = k/nlat
   ! do nothing if the indices is outside of our domain
   if (i<is.or.i>is+size(tiles,1)-1) return ! skip points outside of domain
   if (j<js.or.j>js+size(tiles,2)-1) return ! skip points outside of domain
   ! loop through the list of tiles at the given point to find k+1st tile
   ce = first_elmt (tiles(i,j))
   te = tail_elmt  (tiles(i,j))
   do while(ce/=te.and.k>0)
     ce=next_elmt(ce); k = k-1
   enddo
   ! set the returned pointer   
   ! NOTE that if (ce==te) at the end of the loop (that is, there are less
   ! tiles in the list then requested by the idx), current_tile(ce) returns
   ! NULL
   ptr=>current_tile(ce)
   
end subroutine


! ============================================================================
! given the netcdf file id, name of the variable, and accessor subroutine, this
! subroutine reads integer 2D data (a scalar value per grid cell, that's why 
! there is 0d in the name of this subroutine) and assigns the input values to 
! each tile in respective grid cell.  
subroutine read_tile_data_i0d_fptr(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   ! subroutine returning the pointer to the data to be written
   interface ; subroutine fptr(tile, ptr)
      use land_tile_mod, only : land_tile_type
      type(land_tile_type), pointer :: tile ! input
      integer             , pointer :: ptr  ! returned pointer to the data
    end subroutine fptr 
   end interface
   
   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_i0d_fptr'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index 
   integer, allocatable :: x1d(:) ! storage for the data
   integer :: i,j,bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile   
   integer, pointer :: ptr
   
   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize = min(INPUT_BUF_SIZE,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_int(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(INPUT_BUF_SIZE,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                              lnd%is,lnd%js, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine


! ============================================================================
subroutine read_tile_data_r0d_fptr(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   ! subroutine returning the pointer to the data to be written
   interface ; subroutine fptr(tile, ptr)
      use land_tile_mod, only : land_tile_type
      type(land_tile_type), pointer :: tile ! input
      real                , pointer :: ptr  ! returned pointer to the data
   end subroutine fptr
   end interface
   
   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r0d_fptr'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index 
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i, j, bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile   
   real, pointer :: ptr
   
   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(INPUT_BUF_SIZE,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                              lnd%is,lnd%js, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine

! ============================================================================
subroutine read_tile_data_r1d_fptr_all(ncid,name,fptr)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   ! subroutine returning the pointer to the data to be written
   interface ; subroutine fptr(tile, ptr)
      use land_tile_mod, only : land_tile_type
      type(land_tile_type), pointer :: tile ! input
      real                , pointer :: ptr(:) ! returned pointer to the data
   end subroutine fptr
   end interface
   
   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r1d_fptr'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(2) ! IDs of the variable dimensions
   integer :: dimlen(2) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:)   ! storage for compressed index 
   real   , allocatable :: x1d(:)   ! storage for the data
   integer :: i, j, bufsize
   integer :: varid,idxid
   integer :: start(2), count(2) ! input slab parameters
   type(land_tile_type), pointer :: tileptr ! pointer to tile   
   real, pointer :: ptr(:)
   
   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=2) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 2-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(INPUT_BUF_SIZE,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize*dimlen(2)))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! set up slab parameters
      start(1) = j ; count(1) = min(bufsize,dimlen(1)-j+1)
      start(2) = 1 ; count(2) = dimlen(2)
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,start(1),count(1),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,start,count,x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                              lnd%is,lnd%js, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr(:) = x1d(i:count(1)*count(2):count(1))
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine


! ============================================================================
subroutine read_tile_data_r1d_fptr_idx (ncid,name,fptr,index)
   integer     , intent(in) :: ncid ! netcdf file id
   character(*), intent(in) :: name ! name of the variable to read
   ! subroutine returning the pointer to the data to be written
   interface ; subroutine fptr(tile, ptr)
      use land_tile_mod, only : land_tile_type
      type(land_tile_type), pointer :: tile ! input
      real                , pointer :: ptr(:) ! returned pointer to the data
   end subroutine fptr 
   end interface
   integer    , intent(in) :: index ! index where to read the data
   
   ! ---- local constants
   character(*), parameter :: module_name='read_tile_data_r0d_fptr'
   ! ---- local vars
   integer :: ndims     ! number of the variable dimensions
   integer :: dimids(1) ! IDs of the variable dimensions
   integer :: dimlen(1) ! size of the variable dimensions
   character(NF_MAX_NAME) :: idxname ! name of the index variable
   integer, allocatable :: idx(:) ! storage for compressed index 
   real   , allocatable :: x1d(:) ! storage for the data
   integer :: i,j,bufsize
   integer :: varid, idxid
   type(land_tile_type), pointer :: tileptr ! pointer to tile   
   real, pointer :: ptr(:)
   
   ! get the number of variable dimensions, and their lengths
   __NF_ASRT__(nfu_inq_var(ncid,name,id=varid,ndims=ndims,dimids=dimids,dimlens=dimlen))
   if(ndims/=1) then
      call error_mesg(module_name,'variable "'//trim(name)//'" has incorrect number of dimensions -- must be 1-dimensional', FATAL)
   endif
   ! get the name of compressed dimension and ID of corresponding variable
   __NF_ASRT__(nfu_inq_dim(ncid,dimids(1),idxname))
   __NF_ASRT__(nfu_inq_var(ncid,idxname,id=idxid))
   ! allocate input buffers for compression index and the variable
   bufsize=min(INPUT_BUF_SIZE,dimlen(1))
   allocate(idx(bufsize),x1d(bufsize))
   ! read the input buffer-by-buffer
   do j = 1,dimlen(1),bufsize
      ! read the index variable
      __NF_ASRT__(nf_get_vara_int(ncid,idxid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),idx))
      ! read the data
      __NF_ASRT__(nf_get_vara_double(ncid,varid,(/j/),(/min(bufsize,dimlen(1)-j+1)/),x1d))
      ! distribute the data over the tiles
      do i = 1, min(bufsize,dimlen(1)-j+1)
         call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                              lnd%is,lnd%js, tileptr)
         call fptr(tileptr, ptr)
         if(associated(ptr)) ptr(index) = x1d(i)
      enddo
   enddo
   ! release allocated memory
   deallocate(idx,x1d)
end subroutine


! ============================================================================
! The subroutines write_tile_data_* below write tiled data (that is, data provided
! in arrays congruous with current tiling) to NetCDF files using "compression
! by gathering" (see CF conventions). They assume that the compressed dimension
! is already created, has certain name (see parameter "tile_index_name" at the beginning
! of this file), and has length equal to the number of actually used tiles in the 
! current domain.

! ============================================================================
! writes out 1-d integer tiled data using "compression by gathering"
subroutine write_tile_data_i1d(ncid,name,data,mask,long_name,units)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name
  integer         , intent(inout) :: data(:) ! data to write
  integer         , intent(inout) :: mask(:) ! mask of valid data
  character(len=*), intent(in), optional :: units, long_name
  ! data and mask are "inout" to save the memory on send-receive buffers. On the
  ! root io_domain PE mask is destroyed and data is filled with the information 
  ! from other PEs in our io_domain. On other PEs these arrays reman intact.

  ! local vars
  integer :: varid,iret,p
  integer, allocatable :: buffer(:) ! data buffers

  ! if our PE doesn't do io (that is, it isn't the root io_domain processor),  
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(data(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_3)
     call mpp_send(mask(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_4)
  else
     ! gather data and masks from all processors in io_domain
     allocate(buffer(size(data)))
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(1), glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_3)
        call mpp_recv(mask(1),   glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_4)
        where (mask>0) data = buffer
     enddo
     ! clean up allocated memory
     deallocate(buffer)
   
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        __NF_ASRT__(nfu_def_var(ncid,name,NF_INT,(/tile_index_name/),long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
     __NF_ASRT__(nf_put_var_int(ncid,varid,data))
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine


! ============================================================================
! writes out 1-d real tiled data using "compression by gathering"
subroutine write_tile_data_r1d(ncid,name,data,mask,long_name,units)
  integer         , intent(in) :: ncid    ! netcdf ID
  character(len=*), intent(in) :: name    ! name of the variable
  real            , intent(inout) :: data(:) ! data to write
  integer         , intent(inout) :: mask(:) ! mask of valid data
  character(len=*), intent(in), optional :: units, long_name ! attributes
  ! data and mask are "inout" to save the memory on send-receive buffers. On the
  ! root io_domain PE mask is destroyed and data is filled with the information 
  ! from other PEs in our io_domain. On other PEs these arrays reman intact.

  ! ---- local vars
  integer :: varid,iret,p
  real,    allocatable :: buffer(:) ! data buffer

  ! if our PE doesn't do io (that is, it isn't the root io_domain processor),  
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(data(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_5)
     call mpp_send(mask(1), plen=size(data), to_pe=lnd%io_pelist(1), tag=COMM_TAG_6)
  else
     ! gather data and masks from the processors in io_domain
     allocate(buffer(size(data)))
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(1), glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_5)
        call mpp_recv(mask(1),   glen=size(data), from_pe=lnd%io_pelist(p), tag=COMM_TAG_6)
        where(mask>0) data = buffer
     enddo
     ! clean up allocated memory
     deallocate(buffer)
     
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,(/tile_index_name/),long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
     __NF_ASRT__(nf_put_var_double(ncid,varid,data))
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine


! ============================================================================
! writes out 2-d real tiled data using "compression by gathering". The dimension
! of the data is (tile,z), and both tile and z dimensions are assumed to be 
! already created
subroutine write_tile_data_r2d(ncid,name,data,mask,zdim,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  real            , intent(inout) :: data(:,:) ! (tile,z)
  integer         , intent(inout) :: mask(:) ! mask of valid data
  character(len=*), intent(in), optional :: units, long_name
  ! data and mask are "inout" to save the memory on send-receive buffers. On the
  ! root io_domain PE mask is destroyed and data is filled with the information 
  ! from other PEs in our io_domain. On other PEs these arrays reman intact.

  ! local vars
  integer :: varid,iret,p,i
  character(NF_MAX_NAME)::dimnames(2)
  real, allocatable :: buffer(:,:) ! send/receive buffer

  ! if our PE doesn't do io (that is, it isn't the root io_domain processor),  
  ! simply send the data and mask of valid data to the root IO processor
  if (mpp_pe()/=lnd%io_pelist(1)) then
     call mpp_send(data(1,1), plen=size(data),   to_pe=lnd%io_pelist(1), tag=COMM_TAG_7)
     call mpp_send(mask(1),   plen=size(data,1), to_pe=lnd%io_pelist(1), tag=COMM_TAG_8)
  else
     allocate(buffer(size(data,1),size(data,2)))
     ! gather data and masks from the processors in our io_domain
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(buffer(1,1), glen=size(data),   from_pe=lnd%io_pelist(p), tag=COMM_TAG_7)
        call mpp_recv(mask(1),     glen=size(data,1), from_pe=lnd%io_pelist(p), tag=COMM_TAG_8)
        do i=1,size(data,1)
           if(mask(i)>0) data(i,:) = buffer(i,:)
        enddo
     enddo
     ! clean up allocated memory
     deallocate(buffer)
   
     ! create variable, if it does not exist
     if(nf_inq_varid(ncid,name,varid)/=NF_NOERR) then
        dimnames(1) = tile_index_name
        dimnames(2) = zdim
        __NF_ASRT__(nfu_def_var(ncid,name,NF_DOUBLE,dimnames,long_name,units,varid))
     endif
     ! write data
     iret = nf_enddef(ncid) ! ignore errors: its OK if file is in data mode already
     __NF_ASRT__(nf_put_var_double(ncid,varid,data))
  endif
  ! wait for all PEs to finish: necessary because mpp_send doesn't seem to 
  ! copy the data, and therefore on non-root io_domain PE there would be a chance
  ! that the data and mask are destroyed before they are actually sent.
  call mpp_sync()
end subroutine


! ============================================================================
subroutine write_tile_data_i0d_fptr(ncid,name,fptr,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(tile, ptr)
     use land_tile_mod, only : land_tile_type
     type(land_tile_type), pointer :: tile ! input
     integer             , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr
  end interface
  
  ! ---- local vars
  integer, allocatable :: idx(:)    ! index dimension
  integer, allocatable :: data(:)   ! data to be written
  integer, allocatable :: mask(:)   ! mask of valid data
  type(land_tile_type), pointer :: tileptr ! pointer to tile
  integer, pointer     :: ptr ! pointer to the tile data
  integer :: ntiles  ! total number of tiles (length of compressed dimension)
  integer :: i

  ! get the size of the output array. Note that at this point the variable
  ! might not yet exist, so we can't use nfu_inq_var
  __NF_ASRT__(nfu_inq_dim(ncid,tile_index_name,len=ntiles))

  ! allocate data
  allocate(data(ntiles),idx(ntiles), mask(ntiles))
  ! fill the data with initial values. This is for the case when some of the
  ! compressed tile indices are invalid, so that corresponding indices of the
  ! array are skipped in the loop below. The invalid indices occur when a restart
  ! is written for the domain where no tiles exist, e.g. the ocean-covered 
  ! region
  data = NF_FILL_INT
  mask = 0

  ! read tile index
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,tile_index_name,idx))

  ! gather data into an array along the tile dimension. It is assumed that 
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                          lnd%is,lnd%js, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) then
        data(i) = ptr
        mask(i) = 1
     endif
  enddo

  ! write data
  call write_tile_data_i1d(ncid,name,data,mask,long_name,units)
  
  ! release allocated memory
  deallocate(data,idx,mask)
end subroutine


! ============================================================================
subroutine write_tile_data_r0d_fptr(ncid,name,fptr,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(tile, ptr)
     use land_tile_mod, only : land_tile_type
     type(land_tile_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr
  end interface
  
  ! ---- local vars
  integer, allocatable :: mask(:)   ! mask of valid data
  integer, allocatable :: idx(:)    ! index dimension
  real   , allocatable :: data(:)   ! data to be written
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr ! pointer to the tile data
  integer :: ntiles  ! total number of tiles (length of compressed dimension)
  integer :: i
  
  ! get the size of the output array. Note that at this point the variable
  ! might not yet exist, so we can't use nfu_inq_var
  __NF_ASRT__(nfu_inq_dim(ncid,tile_index_name,len=ntiles))

  ! allocate data
  allocate(data(ntiles),idx(ntiles),mask(ntiles))
  data = NF_FILL_DOUBLE
  mask = 0

  ! read tile index
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,tile_index_name,idx))

  ! gather data into an array along the tile dimension. It is assumed that 
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                          lnd%is,lnd%js, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) then
        data(i) = ptr
        mask(i) = 1
     endif
  enddo

  ! write data
  call write_tile_data_r1d(ncid,name,data,mask,long_name,units)

  ! free allocated memory
  deallocate(data,idx,mask)
end subroutine


! ============================================================================
subroutine write_tile_data_r1d_fptr_all(ncid,name,fptr,zdim,long_name,units)
  integer         , intent(in) :: ncid ! netcdf id
  character(len=*), intent(in) :: name ! name of the variable to write
  character(len=*), intent(in) :: zdim ! name of the z-dimension
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(tile, ptr)
     use land_tile_mod, only : land_tile_type
     type(land_tile_type), pointer :: tile ! input
     real                , pointer :: ptr(:) ! returned pointer to the data
  end subroutine fptr
  end interface
  
  ! ---- local vars
  integer, allocatable :: idx(:)    ! index dimension
  real   , allocatable :: data(:,:) ! data to be written
  integer, allocatable :: mask(:)   ! mask of valid data
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr(:) ! pointer to the tile data
  integer :: ntiles  ! total number of tiles (length of compressed dimension)
  integer :: i
  integer :: nlev ! number of levels of the output variable
  
  ! get the size of the output array. Note that at this point the variable
  ! might not yet exist, so we can't use nfu_inq_var
  __NF_ASRT__(nfu_inq_dim(ncid,tile_index_name,len=ntiles))
  __NF_ASRT__(nfu_inq_dim(ncid,zdim,len=nlev))

  ! allocate data
  allocate(data(ntiles,nlev),idx(ntiles),mask(ntiles))
  data = NF_FILL_DOUBLE
  mask = 0

  ! read tile index
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,tile_index_name,idx))

  ! gather data into an array along the tile dimension. It is assumed that 
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                          lnd%is,lnd%js, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) then
        data(i,:) = ptr(:)
        mask(i) = 1
     endif
  enddo

  ! write data
  call write_tile_data_r2d(ncid,name,data,mask,zdim,long_name,units)
  
  ! free allocated memory
  deallocate(data,idx)
  
end subroutine


! ============================================================================
subroutine write_tile_data_r1d_fptr_idx(ncid,name,fptr,index,long_name,units)
  integer         , intent(in) :: ncid  ! netcdf id
  character(len=*), intent(in) :: name  ! name of the variable to write
  integer         , intent(in) :: index ! index of the fptr array element to 
                                        ! write out
  character(len=*), intent(in), optional :: units, long_name
  ! subroutine returning the pointer to the data to be written
  interface ; subroutine fptr(tile, ptr)
     use land_tile_mod, only : land_tile_type
     type(land_tile_type), pointer :: tile ! input
     real                , pointer :: ptr(:) ! returned pointer to the data
  end subroutine fptr 
  end interface
  
  ! ---- local vars
  integer, allocatable :: idx(:)    ! index dimension
  real   , allocatable :: data(:)   ! data to be written
  integer, allocatable :: mask(:)   ! mask of valid data
  type(land_tile_type), pointer :: tileptr ! pointer to tiles
  real   , pointer :: ptr(:) ! pointer to the tile data
  integer :: ntiles  ! total number of tiles (length of compressed dimension)
  integer :: i

  ! get the size of the output array. Note that at this point the variable
  ! might not yet exist, so we can't use nfu_inq_var
  __NF_ASRT__(nfu_inq_dim(ncid,tile_index_name,len=ntiles))

  ! allocate data
  allocate(data(ntiles),idx(ntiles),mask(ntiles))
  data = NF_FILL_DOUBLE
  mask = 0

  ! read tile index
  i = nf_enddef(ncid) ! ignore errors (file may be in data mode already)
  __NF_ASRT__(nfu_get_var(ncid,tile_index_name,idx))

  ! gather data into an array along the tile dimension. It is assumed that 
  ! the tile dimension spans all the tiles that need to be written.
  do i = 1, size(idx)
     call get_tile_by_idx(idx(i),lnd%nlon,lnd%nlat,lnd%tile_map,&
                          lnd%is,lnd%js, tileptr)
     call fptr(tileptr, ptr)
     if(associated(ptr)) then
        data(i) = ptr(index)
        mask(i) = 1
     endif
  enddo

  ! write data
  call write_tile_data_r1d(ncid,name,data,mask,long_name,units)

  ! free allocated memory
  deallocate(data,idx)
end subroutine


! ============================================================================
subroutine override_tile_data_r0d_fptr(fieldname,fptr,time,override)
  character(len=*), intent(in)   :: fieldname ! field to override
  type(time_type),  intent(in)   :: time      ! model time
  logical, optional, intent(out) :: override  ! true if the field has been 
                                              ! overridden successfully
  ! subroutine returning the pointer to the data to be overridden
  interface ; subroutine fptr(tile, ptr)
     use land_tile_mod, only : land_tile_type
     type(land_tile_type), pointer :: tile ! input
     real                , pointer :: ptr  ! returned pointer to the data
  end subroutine fptr
  end interface

  ! ---- local vars
  real    :: data2D(lnd%is:lnd%ie,lnd%js:lnd%je) ! storage for the input data
  logical :: override_
  
  call data_override('LND',fieldname,data2D, time, override_ )
  if(present(override)) override=override_
  if(.not.override_) return ! do nothing if the field was not overridden 

  ! distribute the data over the tiles
  call put_to_tiles_r0d_fptr(data2d,lnd%tile_map,fptr)
  
end subroutine

! =============================================================================
! given netcdf ID, synchronizes the definitions between writing and reading 
! processors
subroutine sync_nc_files(ncid)
  integer, intent(in) :: ncid

  integer :: iret

  if(mpp_pe()==lnd%io_pelist(1)) then
     iret = nf_enddef(ncid)
     ! commit possible definition changes and data to the disk
     __NF_ASRT__(nf_sync(ncid))
  endif
  call mpp_sync()
  if(mpp_pe()/=lnd%io_pelist(1)) then
     ! synchronize in-memory data structures with the changes on the disk
    __NF_ASRT__(nf_sync(ncid))
  endif
end subroutine sync_nc_files

end module

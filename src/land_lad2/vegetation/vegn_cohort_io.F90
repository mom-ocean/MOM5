module cohort_io_mod

use fms_mod,          only : error_mesg, FATAL, WARNING
use mpp_mod,          only : mpp_pe, mpp_max, mpp_send, mpp_recv, mpp_sync, &
                             COMM_TAG_1, COMM_TAG_2
use nf_utils_mod,     only : nfu_inq_dim, nfu_get_var, nfu_put_var, &
     nfu_get_rec, nfu_put_rec, nfu_def_dim, nfu_def_var, nfu_put_att, &
     nfu_inq_var
use land_io_mod,      only : print_netcdf_error
use land_tile_mod,    only : land_tile_type, land_tile_list_type, &
     land_tile_enum_type, first_elmt, tail_elmt, next_elmt, get_elmt_indices, &
     current_tile, operator(/=)
use land_tile_io_mod, only : get_tile_by_idx, sync_nc_files

use vegn_cohort_mod, only: vegn_cohort_type
use land_data_mod, only : lnd

implicit none
private

! ==== public interfaces =====================================================
! input
public :: read_create_cohorts
public :: read_cohort_data_r0d_fptr
public :: read_cohort_data_i0d_fptr
! output
public :: create_cohort_dimension
public :: write_cohort_data_r0d_fptr
public :: write_cohort_data_i0d_fptr
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'cohort_io_mod', &
     version     = '$Id: vegn_cohort_io.F90,v 19.0.4.2 2012/05/14 19:18:34 Zhi.Liang Exp $', &
     tagname     = '$Name: siena_201207 $'
! name of the "compressed" dimension (and dimension variable) in the output 
! netcdf files -- that is, the dimensions written out using compression by 
! gathering, as described in CF conventions.
character(len=*),   parameter :: cohort_index_name   = 'cohort_index'
integer, parameter :: INPUT_BUF_SIZE = 1024 ! max size of the input buffer for
                                     ! cohort input

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ============================================================================
! given compressed index, sizes of the global grid, 2D array of tile lists
! and the lower boundaries of this array, returns a pointer to the cohort
! corresponding to the compressed index, or NULL is the index is outside 
! current domain, or such tile does not exist, or such cohort does not exist.
subroutine get_cohort_by_idx(idx,nlon,nlat,ntiles,tiles,is,js,ptr)
   integer, intent(in) :: idx ! index
   integer, intent(in) :: nlon, nlat, ntiles
   integer, intent(in) :: is, js
   type(land_tile_list_type), intent(in) :: tiles(is:,js:)
   type(vegn_cohort_type), pointer :: ptr
   
   ! ---- local vars
   integer :: tile_idx, k
   type(land_tile_type), pointer :: tile
   
   ptr=>NULL()
   
   tile_idx = modulo(idx,nlon*nlat*ntiles)
   call get_tile_by_idx(tile_idx,nlon,nlat,tiles,is,js,tile)
   if(associated(tile)) then
      if (associated(tile%vegn)) then
         k = idx/(nlon*nlat*ntiles) ! calculate cohort index within a tile
         ptr=>tile%vegn%cohorts(k+1)
      endif
   endif

end subroutine

! ============================================================================
subroutine read_create_cohorts(ncid)
  integer, intent(in) :: ncid

  integer :: ncohorts ! total number of cohorts in restart file
  integer :: nlon, nlat, ntiles ! size of respective dimensions
 
  integer, allocatable :: idx(:)
  integer :: i,j,t,k,m, n, nn, idxid
  integer :: bufsize
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type), pointer :: tile
  character(len=64) :: info ! for error message

  ! get the size of dimensions
  nlon = lnd%nlon ; nlat = lnd%nlat
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))

  ! read the cohort index
  __NF_ASRT__(nfu_inq_dim(ncid,cohort_index_name,len=ncohorts))
  __NF_ASRT__(nfu_inq_var(ncid,cohort_index_name,id=idxid))
  bufsize = min(INPUT_BUF_SIZE,ncohorts)
  allocate(idx(bufsize))
  
  do nn = 1, ncohorts, bufsize
     __NF_ASRT__(nf_get_vara_int(ncid,idxid,nn,min(bufsize,ncohorts-nn+1),idx))
     
     do n = 1,min(bufsize,ncohorts-nn+1)
        if(idx(n)<0) cycle ! skip illegal indices
        k = idx(n)
        i = modulo(k,nlon)+1   ; k = k/nlon
        j = modulo(k,nlat)+1   ; k = k/nlat
        t = modulo(k,ntiles)+1 ; k = k/ntiles
        k = k+1
        
        if (i<lnd%is.or.i>lnd%ie) cycle ! skip points outside of domain
        if (j<lnd%js.or.j>lnd%je) cycle ! skip points outside of domain
        
        ce = first_elmt(lnd%tile_map(i,j))
        do m = 1,t-1
           ce=next_elmt(ce)
        enddo
        tile=>current_tile(ce)
        if(.not.associated(tile%vegn)) then
           info = ''
           write(info,'("(",3i3,")")')i,j,t
           call error_mesg('read_create_cohort',&
                'vegn tile'//trim(info)//' does not exist, but is necessary to create a cohort', &
                WARNING)
        else
           tile%vegn%n_cohorts = tile%vegn%n_cohorts + 1
        endif
     enddo
  enddo

  ! go through all tiles in the domain and allocate requested numner of cohorts
  ce = first_elmt(lnd%tile_map); te = tail_elmt(lnd%tile_map)
  do while (ce/=te)
     tile=>current_tile(ce); ce = next_elmt(ce)
     if(.not.associated(tile%vegn))cycle
     allocate(tile%vegn%cohorts(tile%vegn%n_cohorts))
  enddo

  ! clean up memory
  deallocate(idx)
end subroutine read_create_cohorts


! ============================================================================
! creates cohort dimension, if necessary, in the output restart file. NOTE 
! that this subroutine should be called even if restart has not been created
! (because, for example, there happen to be no vegetation in a certain domain),
! for the reason that it calls mpp_max, and that should be called for each
! processor to work.
subroutine create_cohort_dimension(ncid)
  integer, intent(in) :: ncid


  ! ---- local vars
  type(land_tile_enum_type) :: ce, te ! tile list enumerators
  type(land_tile_type), pointer :: tile
 
  integer, allocatable :: idx(:)   ! integer compressed index of tiles
  integer :: i,j,k,c,n,ntiles,max_cohorts,p
  integer :: iret
  integer, allocatable :: ncohorts(:) ! array of idx sizes from all PEs in io_domain
  integer, allocatable :: idx2(:) ! array of cohort indices from all PEs in io_domain

  ! count total number of cohorts in compute domain and max number of
  ! of cohorts per tile
  ce = first_elmt(lnd%tile_map)
  te = tail_elmt (lnd%tile_map)
  n  = 0
  max_cohorts = 0
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn))then
        n = n+tile%vegn%n_cohorts 
        max_cohorts = max(max_cohorts,tile%vegn%n_cohorts)
     endif
     ce=next_elmt(ce)
  enddo

  call mpp_max(max_cohorts)

  ! get the size of the tile dimension from the file
  __NF_ASRT__(nfu_inq_dim(ncid,'tile',len=ntiles))
  
  ! calculate compressed cohort index to be written to the restart file
  allocate(idx(max(n,1))) ; idx(:) = -1
  ce = first_elmt(lnd%tile_map, lnd%is, lnd%js)
  n = 1
  do while (ce/=te)
     tile=>current_tile(ce)
     if(associated(tile%vegn)) then
        call get_elmt_indices(ce,i,j,k)
        do c = 1,tile%vegn%n_cohorts
           idx (n) = &
                (c-1)*lnd%nlon*lnd%nlat*ntiles + &
                (k-1)*lnd%nlon*lnd%nlat + &
                (j-1)*lnd%nlon + &
                (i-1)        
           n = n+1
        enddo
     endif
     ce=next_elmt(ce)
  end do

  if (mpp_pe()/=lnd%io_pelist(1)) then
     ! if this processor is not doing io (that is, it's not root io_domain
     ! processor), simply send the data to the root io_domain PE
     call mpp_send(size(idx), plen=1,         to_pe=lnd%io_pelist(1), tag=COMM_TAG_1)
     call mpp_send(idx(1),    plen=size(idx), to_pe=lnd%io_pelist(1), tag=COMM_TAG_2)
  else
     ! gather the array of cohort index sizes
     allocate(ncohorts(size(lnd%io_pelist)))
     ncohorts(1) = size(idx)
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(ncohorts(p), from_pe=lnd%io_pelist(p), glen=1, tag=COMM_TAG_1)
     enddo
     ! gather cohort index from the processors in our io_domain
     allocate(idx2(sum(ncohorts(:))))
     idx2(1:ncohorts(1))=idx(:)
     k=ncohorts(1)+1
     do p = 2,size(lnd%io_pelist)
        call mpp_recv(idx2(k), from_pe=lnd%io_pelist(p), glen=ncohorts(p), tag=COMM_TAG_2)
        k = k+ncohorts(p)
     enddo
     ! create cohort dimension in the output file
     iret = nf_redef(ncid)
     __NF_ASRT__(nfu_def_dim(ncid,'cohort',max_cohorts))
     ! create cohort index
     __NF_ASRT__(nfu_def_dim(ncid,cohort_index_name,idx2,'compressed vegetation cohort index'))
     __NF_ASRT__(nfu_put_att(ncid,cohort_index_name,'compress','cohort tile lat lon'))
     __NF_ASRT__(nfu_put_att(ncid,cohort_index_name,'valid_min',0))
     ! deallocate the data we no longer need
     deallocate(ncohorts,idx2)
     ! leave the define mode to commit the new definitions to the disk
     iret = nf_enddef(ncid)
  endif
  call sync_nc_files(ncid)
end subroutine create_cohort_dimension


#define F90_TYPE real
#define NF_TYPE NF_DOUBLE
#define NF_FILL_VALUE NF_FILL_DOUBLE
#define READ_0D_FPTR read_cohort_data_r0d_fptr
#define WRITE_0D_FPTR write_cohort_data_r0d_fptr
#include "vegn_cohort_io.inc"

#define F90_TYPE integer
#define NF_TYPE NF_INT
#define NF_FILL_VALUE NF_FILL_INT
#define READ_0D_FPTR read_cohort_data_i0d_fptr
#define WRITE_0D_FPTR write_cohort_data_i0d_fptr
#include "vegn_cohort_io.inc"

end module cohort_io_mod

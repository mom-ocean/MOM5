! -*-f90-*-
! $Id: read_remap_cohort_data.inc,v 17.0 2009/07/21 03:03:09 fms Exp $

! some sanity checks
#ifndef F90_TYPE
#error F90_TYPE is not defined: must be one of FORTRAN 90 types
#endif

#ifndef READ_REMAP_SUB
#error name of READ_REMAP_SUB is not defined
#endif

! ============================================================================
subroutine READ_REMAP_SUB(ncid,name,fptr,map_i,map_j,rec)
  integer           , intent(in) :: ncid ! netcdf id
  character(len=*)  , intent(in) :: name ! name of the variable to read
  integer           , intent(in) :: map_i(lnd%is:,lnd%js:) ! re-mapping index
  integer           , intent(in) :: map_j(lnd%is:,lnd%js:) ! re-mapping index
  integer, optional , intent(in) :: rec  ! record number (in case there are 
                                         ! several in the file) 
  ! subroutine returning the pointer to the data to be written
  interface
     subroutine fptr(cohort, ptr)
       use vegn_cohort_mod, only : vegn_cohort_type
       type(vegn_cohort_type), pointer :: cohort ! input
       F90_TYPE              , pointer :: ptr    ! returned pointer to the data
     end subroutine fptr
  end interface

  ! ---- local constants
  character(*), parameter :: module_name = "read_remap_cohort_data"

  ! ---- local vars
  integer :: i,j,k,n,ii,jj,ndims, iret
  integer :: rec_     ! record number
  type(land_tile_enum_type) :: ce, te
  type(land_tile_type)   , pointer :: tile
  type(vegn_cohort_type) , pointer :: cohort
  F90_TYPE, pointer :: ptr ! pointer to the individual cohort data
  F90_TYPE, allocatable :: data(:,:,:,:) ! buffer for input data
  logical,  allocatable :: mask(:,:,:,:) ! validity mask for input data
  logical :: has_records, is_compressed
  integer :: dimlens(NF_MAX_VAR_DIMS)
  type(nfu_validtype) :: v

  ! assign the internal record number
  if(present(rec)) then
     rec_ = rec
  else
     rec_ = 1
  endif

  ! get the size of dimensions
  iret=nfu_inq_compressed_var(ncid, name, ndims=ndims, dimlens=dimlens,&
       has_records=has_records, is_compressed=is_compressed)
  __NF_ASRT__(iret)

  ! calculate the dimensions of input buffers, based on the dimensions of
  ! input variable
  if(has_records)ndims = ndims-1
  do i = ndims+1,4
     dimlens(i) = 1
  enddo

  ! allocate input buffers
  allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  allocate(mask(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
  !             lon        lat        tile       cohort

  mask = .FALSE.
  __NF_ASRT__(nfu_get_compressed_rec(ncid,name,rec_,data,mask))
  if (.not.is_compressed) then
     __NF_ASRT__( nfu_get_valid_range(ncid,name,v) )
     mask=nfu_is_valid(data,v)
  endif

  ! distribute data over cohorts. NOTE that this is slightly different from the restart
  ! reading procedure. On reading the restart, all the tiles are counted in sequence,
  ! while here only tne vegetation tiles.
  do j = lnd%js, lnd%je
  do i = lnd%is, lnd%ie
     ii = map_i(i,j); jj = map_j(i,j)
     if ((ii.le.0).or.(jj.le.0)) cycle ! skip un-mapped points
     if (.not.any(mask(ii,jj,:,:))) cycle ! skip points where there is no data 

     ce = first_elmt (lnd%tile_map(i,j))
     te = tail_elmt  (lnd%tile_map(i,j))
     k = 1
tile_loop:  do while(ce/=te.and.k<=dimlens(3))
        tile=>current_tile(ce); ce=next_elmt(ce);
        if (.not.associated(tile%vegn)) cycle
        ! find index of the next valid tile in the input data
        do while(.not.any(mask(ii,jj,k,:)))
           k=k+1 ! go to the next tile if there's no data (i.e. all mask 
                 ! values are false for this tile)
           if(k>dimlens(3)) exit tile_loop 
        enddo
        
        do n = 1,min(size(tile%vegn%cohorts(:)),dimlens(4))
           cohort=>tile%vegn%cohorts(n)
           call fptr(cohort,ptr)
           if(associated(ptr).and.mask(ii,jj,k,n)) ptr = data(ii,jj,k,n)
        enddo
        k = k+1 ! go to the next tile in input data
     enddo tile_loop
  enddo
  enddo
  
  ! free allocated memory
  deallocate(data,mask)

end subroutine 

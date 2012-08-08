! -*-f90-*-
! $Id: getput_compressed.inc,v 17.0 2009/07/21 03:02:47 fms Exp $

! some sanity checks
#ifndef F90_TYPE
#error F90_TYPE is not defined: must be one of FORTRAN 90 types
#endif

#ifndef NF_TYPE
#error NF_TYPE is not defined: must be netcdf type name corresponding to F90_TYPE
#endif

! #### Basic macro definition #################################################

! macro definition for concatenation -- for construction of names based on the
! names of the operations, types, and dimension numbers
#define CONCAT3(op,T,D) op##T##D
#define CONCAT2(op,T)   op##T

! names of the functions we define in this module
#define GET_VAR(T,D) CONCAT3(get_compressed_var_,T,D)
#define PUT_VAR(T,D) CONCAT3(put_compressed_var_,T,D)
#define GET_REC(T,D) CONCAT3(get_compressed_rec_,T,D)

! define names of the corresponding netcdf functions. The two-stage definition is
! necessary because of the preprocessor argument pre-scan rules. See, for example,
! http://gcc.gnu.org/onlinedocs/cpp/Argument-Prescan.html
#define NF_GET_VAR_(T) CONCAT2(nf_get_var_,T)
#define NF_GET_VAR_T NF_GET_VAR_(NF_TYPE) 

#define NF_PUT_VAR_(T) CONCAT2(nf_put_var_,T)
#define NF_PUT_VAR_T NF_PUT_VAR_(NF_TYPE) 

! #### Interface definition ###################################################
! define specific names of the subroutines
#define GET_VAR_D1N GET_VAR(NF_TYPE, D1N)
#define GET_VAR_D1I GET_VAR(NF_TYPE, D1I)

#define PUT_VAR_D1N PUT_VAR(NF_TYPE, D1N)
#define PUT_VAR_D1I PUT_VAR(NF_TYPE, D1I)

#define GET_REC_D4N GET_REC(NF_TYPE,D4N)
#define GET_REC_D3N GET_REC(NF_TYPE,D3N)
#define GET_REC_D2N GET_REC(NF_TYPE,D2N)
#define GET_REC_D1N GET_REC(NF_TYPE,D1N)
#define GET_REC_D1I GET_REC(NF_TYPE,D1I)


#ifdef __INTERFACE_SECTION__

interface nfu_get_compressed_var
   module procedure GET_VAR_D1I, GET_VAR_D1N
end interface
interface nfu_put_compressed_var
   module procedure PUT_VAR_D1I, PUT_VAR_D1N
end interface
interface nfu_get_compressed_rec
   module procedure GET_REC_D1N,GET_REC_D2N,GET_REC_D3N,GET_REC_D4N,GET_REC_D1I
end interface
#endif
! #### END of interface definition ############################################


! #### Implementation definition ##############################################
#ifdef __BODY_SECTION__ 

! ===========================================================================
function GET_VAR_D1N(ncid,name,data,mask) result (iret)
  integer      , intent(in)    :: ncid
  character(*) , intent(in)  :: name
  F90_TYPE     , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = GET_VAR(NF_TYPE,D1I)(ncid,varid,data,mask)
7 return
end function 

! ===========================================================================
function GET_VAR_D1I(ncid,varid,data,mask) result (iret)
  integer          , intent(in)    :: ncid,varid
  F90_TYPE         , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  F90_TYPE, allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=diminfo(i)%length,name=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=dimlen,name=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var(ncid,dimname,diminfo(i)%idx),iret,7)
        ! calculate corresponding stride in output (unpacked) array
        length = 1
        do n = 1,cndims
           length = length*cdimlens(n)
        enddo
     else
        length = diminfo(i)%length
     endif
     diminfo(i)%stride = stride
     stride = stride*length
  enddo
        
  ! get the entire variable
  allocate(buffer(varsize))
  __NF_TRY__(NF_GET_VAR_T(ncid,varid,buffer),iret,7)

  ! move the data to the output buffer
  idx(:) = 0
  do i = 1,size(buffer)
     ! calculate destination index
     ii = 1
     do n = 1,ndims
        if(associated(diminfo(n)%idx)) then
           if(diminfo(n)%idx(idx(n)) >= 0)then
              ii = ii+diminfo(n)%idx(idx(n))*diminfo(n)%stride
           else
              ii = -1 ! set a value flagging an invalid point
              exit    ! from index loop
           endif
        else
           ii = ii+idx(n)*diminfo(n)%stride
        endif
     enddo

     ! if index is negative, skip an invalid point
     if (ii > 0) then
        data(ii) = buffer(i)
        if(present(mask))mask(ii) = .true.
     endif

     ! increment indices
     do n = 1,ndims
        idx(n) = idx(n)+1
        if(idx(n)<diminfo(n)%length)exit
        idx(n) = 0
     enddo
  enddo

7 continue
  ! clean up memory
  do i = 1,size(diminfo)
     if(associated(diminfo(i)%idx)) &
          deallocate(diminfo(i)%idx)
  enddo
  if (allocated(buffer)) &
       deallocate(buffer)
end function


! ===========================================================================
function PUT_VAR_D1N (ncid,name,src) result (iret)
  integer     , intent(in)    :: ncid
  character(*), intent(in)    :: name
  F90_TYPE    , intent(inout) :: src(*)      ! data to write
  integer :: iret

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = PUT_VAR(NF_TYPE,D1I)(ncid,varid,src)
7 return
end function

! ===========================================================================
function PUT_VAR_D1I(ncid,varid,src) result (iret)
  integer  , intent(in)    :: ncid,varid
  F90_TYPE , intent(inout) :: src(*)      ! data to write
  integer :: iret

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  F90_TYPE, allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=diminfo(i)%length,name=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=dimlen,name=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var(ncid,dimname,diminfo(i)%idx),iret,7)
        ! calculate corresponding stride in output (unpacked) array
        length = 1
        do n = 1,cndims
           length = length*cdimlens(n)
        enddo
     else
        length = diminfo(i)%length
     endif
     diminfo(i)%stride = stride
     stride = stride*length
  enddo
        
  ! get the entire variable
  allocate(buffer(varsize))

  ! move the data to the output buffer
  idx(:) = 0
  do i = 1,size(buffer)
     ! calculate destination index
     ii = 1
     do n = 1,ndims
        if(associated(diminfo(n)%idx)) then
           ii = ii+diminfo(n)%idx(idx(n))*diminfo(n)%stride
        else
           ii = ii+idx(n)*diminfo(n)%stride
        endif
     enddo

     buffer(i) = src(ii)

     ! increment indices
     do n = 1,ndims
        idx(n) = idx(n)+1
        if(idx(n)<diminfo(n)%length)exit
        idx(n) = 0
     enddo
  enddo

  __NF_TRY__(NF_PUT_VAR_T(ncid,varid,buffer),iret,7)

7 continue
  ! clean up memory
  do i = 1,size(diminfo)
     if(associated(diminfo(i)%idx)) &
          deallocate(diminfo(i)%idx)
  enddo
  if (allocated(buffer)) &
       deallocate(buffer)
end function


! ===========================================================================
function GET_REC_D4N(ncid,name,rec,data,mask) result (iret)
  integer          , intent(in)    :: ncid   ! ID of the netcdf file
  character(*)     , intent(in)    :: name   ! name of the variable
  integer          , intent(in)    :: rec    ! number of the record to get
  F90_TYPE         , intent(inout) :: data(:,:,:,:)
  logical, optional, intent(inout) :: mask(:,:,:,:)
  integer :: iret ! return value

  if (present(mask)) then
     iret = GET_REC_D1N(ncid,name,rec,data,mask)
  else
     iret = GET_REC_D1N(ncid,name,rec,data)
  endif
end function 

! ===========================================================================
function GET_REC_D3N(ncid,name,rec,data,mask) result (iret)
  integer          , intent(in)    :: ncid   ! ID of the netcdf file
  character(*)     , intent(in)    :: name   ! name of the variable
  integer          , intent(in)    :: rec    ! number of the record to get
  F90_TYPE         , intent(inout) :: data(:,:,:)
  logical, optional, intent(inout) :: mask(:,:,:)
  integer :: iret ! return value

  if (present(mask)) then
     iret = GET_REC_D1N(ncid,name,rec,data,mask)
  else
     iret = GET_REC_D1N(ncid,name,rec,data)
  endif
end function 

! ===========================================================================
function GET_REC_D2N(ncid,name,rec,data,mask) result (iret)
  integer          , intent(in)    :: ncid   ! ID of the netcdf file
  character(*)     , intent(in)    :: name   ! name of the variable
  integer          , intent(in)    :: rec    ! number of the record to get
  F90_TYPE         , intent(inout) :: data(:,:)
  logical, optional, intent(inout) :: mask(:,:)
  integer :: iret ! return value

  if (present(mask)) then
     iret = GET_REC_D1N(ncid,name,rec,data,mask)
  else
     iret = GET_REC_D1N(ncid,name,rec,data)
  endif
end function 

! ===========================================================================
function GET_REC_D1N(ncid,name,rec,data,mask) result (iret)
  integer          , intent(in)    :: ncid   ! ID of the netcdf file
  character(*)     , intent(in)    :: name   ! name of the variable
  integer          , intent(in)    :: rec    ! number of the record to get
  F90_TYPE         , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = GET_REC_D1I(ncid,varid,rec,data,mask)
7 return

end function 

! ===========================================================================
function GET_REC_D1I(ncid,varid,rec,data,mask) result (iret)
  integer          , intent(in)    :: ncid,varid
  integer          , intent(in)    :: rec    ! number of the record to get
  F90_TYPE         , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: recsize
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  F90_TYPE, allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride
  logical :: has_records

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,recsize=recsize,has_records=has_records)
  __NF_TRY__(iret,iret,7)

  ! adjust the number of dimensions for one record: if there is a record dimension,
  ! it has to be the last
  if (has_records)ndims=ndims-1
  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=diminfo(i)%length,name=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),len=dimlen,name=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var(ncid,dimname,diminfo(i)%idx),iret,7)
        ! calculate corresponding stride in output (unpacked) array
        length = 1
        do n = 1,cndims
           length = length*cdimlens(n)
        enddo
     else
        length = diminfo(i)%length
     endif
     diminfo(i)%stride = stride
     stride = stride*length
  enddo
        
  ! get the entire variable
  allocate(buffer(recsize))
  __NF_TRY__(nfu_get_rec(ncid,varid,rec,buffer),iret,7)

  ! move the data to the output buffer
  idx(:) = 0
  do i = 1,size(buffer)
     ! calculate destination index
     ii = 1
     do n = 1,ndims
        if(associated(diminfo(n)%idx)) then
           if(diminfo(n)%idx(idx(n)) >= 0)then
              ii = ii+diminfo(n)%idx(idx(n))*diminfo(n)%stride
           else
              ii = -1 ! set a value flagging an invalid point
              exit    ! from index loop
           endif
        else
           ii = ii+idx(n)*diminfo(n)%stride
        endif
     enddo

     ! if index is negative, skip an invalid point
     if (ii > 0) then
        data(ii) = buffer(i)
        if(present(mask))mask(ii) = .true.
     endif

     ! increment indices
     do n = 1,ndims
        idx(n) = idx(n)+1
        if(idx(n)<diminfo(n)%length)exit
        idx(n) = 0
     enddo
  enddo

7 continue
  ! clean up memory
  do i = 1,size(diminfo)
     if(associated(diminfo(i)%idx)) &
          deallocate(diminfo(i)%idx)
  enddo
  if (allocated(buffer)) &
       deallocate(buffer)
end function


#endif

! #### Undefine all auxilary definitions #######################################
#undef GET_VAR
#undef GET_REC
#undef PUT_VAR
#undef GET_REC_D4N
#undef CONCAT2
#undef CONCAT3
#undef NF_GET_VAR_
#undef NF_GET_VAR_T
#undef NF_PUT_VAR_
#undef NF_PUT_VAR_T

#undef GET_VAR_D1N
#undef GET_VAR_D1I 
#undef PUT_VAR_D1N 
#undef PUT_VAR_D1I 
#undef GET_REC_D4N 
#undef GET_REC_D3N 
#undef GET_REC_D2N 
#undef GET_REC_D1N 
#undef GET_REC_D1I 

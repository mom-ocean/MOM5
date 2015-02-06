!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
!
!-----------------------------------------------------------------------
#define __NF_TRY__(err_code,iret,LABEL) iret=err_code;\
call cdfe(iret,"",__LINE__,verb);\
if(iret/=NF_NOERR)goto LABEL


module nfu_compress_mod

  use nfu_mod

implicit none
private

! ==== public interface ======================================================
public :: nfu_inq_compressed_dim, nfu_inq_compressed_var
public :: nfu_get_compressed_var_r8, nfu_get_compressed_var_int
public :: nfu_put_compressed_var_r8, nfu_put_compressed_var_int
public :: nfu_get_compressed_var_r8n
public :: nfu_put_compressed_var_r8n
! ==== end of public interface ===============================================

! ==== interfaces for overloaded functions ===================================
interface nfu_inq_compressed_dim
   module procedure inq_compressed_dim_n
   module procedure inq_compressed_dim_i
end interface

interface nfu_inq_compressed_var
   module procedure inq_compressed_var_n
   module procedure inq_compressed_var_i
end interface

#include <netcdf.inc>

integer :: verb = 0


! ---- private type - used to hold dimension/packing information during unpacking
! (see get_compressed_var_i_r8)
type diminfo_type
   integer, pointer :: idx(:)=>NULL() ! packing information
   integer :: length  ! size of the dimension in the input array
   integer :: stride  ! stide along the dimension in the output array
end type 


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ===========================================================================
function nfu_get_compressed_var_r8n(ncid,name,data,mask) result (iret)
  integer          , intent(in)    :: ncid
  character(*)     , intent(in)  :: name
  real(kind=8)     , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nfu_get_compressed_var_r8(ncid,varid,data,mask)
7 return

end function 

! ===========================================================================
function nfu_get_compressed_var_r8(ncid,varid,data,mask) result (iret)
  integer          , intent(in)    :: ncid,varid
  real(kind=8)     , intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)
  integer :: iret

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  real(kind=8),allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=diminfo(i)%length,dimname=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=dimlen,dimname=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var_int(ncid,dimname,diminfo(i)%idx),iret,7)
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
  __NF_TRY__(nf_get_var_double(ncid,varid,buffer),iret,7)

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
function nfu_get_compressed_var_int(ncid,varid,data,mask) result (iret)
  integer :: iret
  integer, intent(in)    :: ncid,varid
  integer, intent(inout) :: data(*)
  logical, optional, intent(inout) :: mask(*)

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  integer, allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=diminfo(i)%length,dimname=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=dimlen,dimname=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var_int(ncid,dimname,diminfo(i)%idx),iret,7)
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
  __NF_TRY__(nf_get_var_int(ncid,varid,buffer),iret,7)

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
        if(present(mask)) mask(ii) = .true.
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
function nfu_put_compressed_var_r8n(ncid,name,src) result (iret)
  integer     , intent(in)    :: ncid
  character(*), intent(in)    :: name
  real(kind=8), intent(inout) :: src(*)      ! data to write
  integer :: iret

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nfu_put_compressed_var_r8(ncid,varid,src)
7 return
end function

! ===========================================================================
function nfu_put_compressed_var_r8(ncid,varid,src) result (iret)
  integer     , intent(in)    :: ncid,varid
  real(kind=8), intent(inout) :: src(*)      ! data to write
  integer :: iret

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  real(kind=8),allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=diminfo(i)%length,dimname=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=dimlen,dimname=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var_int(ncid,dimname,diminfo(i)%idx),iret,7)
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

  __NF_TRY__(nf_put_var_double(ncid,varid,buffer),iret,7)

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
function nfu_put_compressed_var_int(ncid,varid,src) result (iret)
  integer :: iret
  integer, intent(in)    :: ncid,varid
  integer, intent(inout) :: src(*)

  integer :: ndims,dimids(NF_MAX_VAR_DIMS),dimlen
  integer :: varsize ! total size of the compressed variable
  integer :: cndims, cdimids(NF_MAX_VAR_DIMS),cdimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname

  integer, allocatable :: buffer(:)
  integer :: i, ii, n, length, idx(NF_MAX_VAR_DIMS)
  integer :: stride

  type(diminfo_type) :: diminfo(NF_MAX_VAR_DIMS)

  ! get the information for the compressed variable
  iret = nfu_inq_var(ncid,varid,ndims=ndims,dimids=dimids,varsize=varsize)
  __NF_TRY__(iret,iret,7)

  ! get the compressed dimensions
  stride = 1
  do i = 1,ndims
     __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=diminfo(i)%length,dimname=dimname),iret,7)
     if(nfu_inq_compressed_dim(ncid,dimids(i),&
          ndims=cndims,dimids=cdimids,dimlens=cdimlens)==NF_NOERR) then
        ! it is a compressed dimension; get dimension itself and calculate
        ! get the dimension (that is, compression information)
        __NF_TRY__(nfu_inq_dim(ncid,dimids(i),dimlen=dimlen,dimname=dimname),iret,7)
        allocate(diminfo(i)%idx(0:dimlen-1))
        __NF_TRY__(nfu_get_var_int(ncid,dimname,diminfo(i)%idx),iret,7)
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
     ! increment indices
     do n = 1,ndims
        idx(n) = idx(n)+1
        if(idx(n)<diminfo(n)%length)exit
        idx(n) = 0
     enddo
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
  enddo

  __NF_TRY__(nf_get_var_int(ncid,varid,buffer),iret,7)

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
function inq_compressed_dim_n(ncid,name,ndims,dimids,dimlens,dimid) result (iret)
  integer :: iret
  integer, intent(in)  :: ncid
  character(*), intent(in) :: name
  integer, intent(out), optional :: ndims
  integer, intent(out), optional :: dimids(:)
  integer, intent(out), optional :: dimlens(:)
  integer, intent(out), optional :: dimid

  integer :: dimid_

  __NF_TRY__(nf_inq_dimid(ncid,name,dimid_),iret,7)
  if(present(dimid)) dimid = dimid_
  __NF_TRY__(inq_compressed_dim_i(ncid,dimid_,ndims,dimids,dimlens),iret,7)
7 return
end function

! ===========================================================================
function inq_compressed_dim_i(ncid,dimid,ndims,dimids,dimlens,dimname) result (iret)
  integer :: iret
  integer, intent(in)  :: ncid,dimid
  integer, intent(out), optional :: ndims
  integer, intent(out), optional :: dimids(:)
  integer, intent(out), optional :: dimlens(:)
  character(*), intent(out), optional :: dimname
  
  character(NF_MAX_NAME) :: dimname_
  character(1024) :: compress ! should be more than enough to hold the compression info
  integer :: dimlen,dimid0,n,is,ie

  __NF_TRY__(nfu_inq_dim(ncid,dimid,dimname=dimname_),iret,7)
  if(present(dimname)) dimname = dimname_
  compress = ''
  __NF_TRY__(nfu_get_att(ncid,dimname_,'compress',compress),iret,7)

  ! parse the description of the compression
  ie = len_trim(compress)
  n = 0
  do while(ie>0)
     is = scan(compress(1:ie),' ',back=.true.)
     if(is==ie) then
        ! skip space runs
     else
        n = n+1
        iret = nfu_inq_dim(ncid,compress(is+1:ie),dimlen=dimlen,dimid=dimid0)
        __NF_TRY__(iret,iret,7)
        if(present(dimids)) dimids(n) = dimid0
        if(present(dimlens)) dimlens(n) = dimlen
     endif
     ie = is-1
  enddo
  if(present(ndims))ndims=n
7 return
end function

! ============================================================================
function inq_compressed_var_n(ncid, name, id, xtype, ndims, dimids, dimlens, natts, &
     is_dim, has_records, varsize, recsize, nrec, is_compressed) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(*),intent(in) :: name
  integer, intent(out), optional :: id
  integer, intent(out), optional :: xtype
  integer, intent(out), optional :: ndims
  integer, intent(out), optional :: dimids(:)
  integer, intent(out), optional :: dimlens(:)
  integer, intent(out), optional :: natts
  logical, intent(out), optional :: is_dim ! true if variable is a dimension variable
  logical, intent(out), optional :: has_records ! true if variable depends on record dimension
  integer, intent(out), optional :: varsize ! total size of the variable
  integer, intent(out), optional :: recsize ! size of a single record
  integer, intent(out), optional :: nrec    ! number of records
  logical, intent(out), optional :: is_compressed ! true if variable is actually compressed

  integer :: vid
  character(len=NF_MAX_NAME) :: vname

  __NF_TRY__(nf_inq_varid(ncid,name,vid),iret,7)
  if(present(id)) id = vid
  iret = inq_compressed_var_i(ncid,vid,vname,xtype,ndims,dimids,dimlens,natts,&
       is_dim,has_records,varsize,recsize,nrec,is_compressed)

7 return  
end function

! ============================================================================
function inq_compressed_var_i(ncid, vid, name, xtype, ndims, dimids, dimlens, &
     natts, is_dim, has_records, varsize, recsize, nrec, is_compressed) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  integer, intent(in) :: vid
  character(*),intent(out), optional :: name
  integer, intent(out), optional :: xtype
  integer, intent(out), optional :: ndims
  integer, intent(out), optional :: dimids(:)
  integer, intent(out), optional :: dimlens(:)
  integer, intent(out), optional :: natts
  logical, intent(out), optional :: is_dim ! true if variable is a dimension variable
  logical, intent(out), optional :: has_records ! true if variable depends on record dimension
  integer, intent(out), optional :: varsize ! total size of the variable
  integer, intent(out), optional :: recsize ! size of a single record
  integer, intent(out), optional :: nrec    ! number of records
  logical, intent(out), optional :: is_compressed ! true if variable is actually compressed

  
  integer :: nd0, dids0(NF_MAX_VAR_DIMS),dlens0(NF_MAX_VAR_DIMS)
  integer :: nd1, dids1(NF_MAX_VAR_DIMS),dlens1(NF_MAX_VAR_DIMS)
  integer :: i,n,unlimdim,vsize,rsize

  iret =  nfu_inq_var(ncid, vid, name, xtype, nd0, dids0, dlens0, natts, &
     is_dim, has_records, varsize, recsize, nrec)

  nd1=1
  if(present(is_compressed)) is_compressed=.false.
  do i = 1, nd0
     if(nfu_inq_compressed_dim(ncid,dids0(i),&
          ndims=n,dimids=dids1(nd1:),dimlens=dlens1(nd1:))==NF_NOERR) then
        nd1 = nd1+n
        if(present(is_compressed)) is_compressed=.true.
     else
        dlens1(nd1) = dlens0(i)
        dids1(nd1) = dids0(i)
        nd1 = nd1+1
     endif
  enddo
  nd1 = nd1-1

  if(present(ndims))   ndims   = nd1
  if(present(dimids))  dimids  = dids1
  if(present(dimlens)) dimlens = dlens1
  if(present(varsize).or.present(recsize)) then
     __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
     vsize = 1; rsize=1
     do i = 1,nd1
        vsize = vsize*dlens1(i)
        if(dids1(i)/=unlimdim)&
             rsize = rsize*dlens1(i)
     enddo
     if (present(varsize)) varsize=vsize
     if (present(recsize)) recsize=rsize
  end if
7 return

end function

end module nfu_compress_mod

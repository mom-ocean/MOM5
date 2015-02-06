module nfc_mod

  use nfu_mod

implicit none
private

! ==== public interface ======================================================
public :: nfu_inq_compressed_dim, nfu_inq_compressed_var
public :: nfu_get_compressed_var
public :: nfu_put_compressed_var
public :: nfu_get_compressed_rec
! ==== end of public interface ===============================================

! ==== interfaces for overloaded functions ===================================
#define __INTERFACE_SECTION__
interface nfu_inq_compressed_dim
   module procedure inq_compressed_dim_n, inq_compressed_dim_i
end interface

interface nfu_inq_compressed_var
   module procedure inq_compressed_var_n, inq_compressed_var_i
end interface

#define F90_TYPE real(8)
#define NF_TYPE  double
#include "getput_compressed.inc"

#define F90_TYPE integer
#define NF_TYPE  int
#include "getput_compressed.inc"

#undef __INTERFACE_SECTION__
! ---- module constants ------------------------------------------------------
character(len=*), parameter :: &
     version = '$Id: nfc.F90,v 20.0 2013/12/13 23:30:40 fms Exp $', &
     tagname = '$Name: tikal $'

! ---- private type - used to hold dimension/packing information during unpacking
! (see get_compressed_var_i_r8)
type diminfo_type
   integer, pointer :: idx(:)=>NULL() ! packing information
   integer :: length  ! size of the dimension in the input array
   integer :: stride  ! stide along the dimension in the output array
end type 

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_TRY__(err_code,iret,LABEL)iret=err_code;if(iret/=NF_NOERR)goto LABEL
contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


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
  integer :: dimlen,dimid0,varid,n,is,ie

  __NF_TRY__(nfu_inq_dim(ncid,dimid,name=dimname_),iret,7)
  if(present(dimname)) dimname = dimname_
  compress = ''
  __NF_TRY__(nf_inq_varid(ncid,dimname_,varid),iret,7)
  __NF_TRY__(nf_get_att_text(ncid,varid,'compress',compress),iret,7)

  ! parse the description of the compression
  ie = len_trim(compress)
  n = 0
  do while(ie>0)
     is = scan(compress(1:ie),' ',back=.true.)
     if(is==ie) then
        ! skip space runs
     else
        n = n+1
        iret = nfu_inq_dim(ncid,compress(is+1:ie),len=dimlen,dimid=dimid0)
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
  logical :: compressed

  iret =  nfu_inq_var(ncid, vid, name, xtype, nd0, dids0, dlens0, natts, &
     is_dim, has_records, varsize, recsize, nrec)

  nd1=1
  if(present(is_compressed)) is_compressed=.false.
  do i = 1, nd0
     __NF_TRY__(nfu_inq_dim(ncid,dids0(i),is_compressed=compressed),iret,7)
     if (compressed) then
        iret = nfu_inq_compressed_dim(ncid,dids0(i),&
                             ndims=n,dimids=dids1(nd1:),dimlens=dlens1(nd1:))
        if (iret/=NF_NOERR) goto 7
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

#define __BODY_SECTION__
#define F90_TYPE real(8)
#define NF_TYPE  double
#include "getput_compressed.inc"

#define F90_TYPE integer
#define NF_TYPE  int
#include "getput_compressed.inc"

end module nfc_mod

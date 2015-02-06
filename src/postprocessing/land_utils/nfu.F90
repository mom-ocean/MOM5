!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
!
!-----------------------------------------------------------------------
#define __NF_TRY__(err_code,iret,LABEL) iret=err_code;\
call cdfe(iret,"",__LINE__,verb);\
if(iret/=NF_NOERR)goto LABEL

module nfu_mod

implicit none
private

! ==== public interface ======================================================
public :: nfu_def_var, nfu_def_dim
public :: nfu_inq_var, nfu_inq_dim, nfu_inq_att
public :: nfu_get_var_r8, nfu_get_var_r4, nfu_get_var_int
public :: nfu_put_var_r8, nfu_put_var_r4, nfu_put_var_int
public :: nfu_get_rec_r8, nfu_get_rec_r4
public :: nfu_put_rec_r8, nfu_put_rec_r4
public :: nfu_get_att, nfu_put_att, nfu_append_att, nfu_copy_att
public :: nfu_clone_dim, nfu_clone_var
public :: nfu_copy_var_data
public :: nfu_check_err
public :: nfu_get_valid_range, is_valid, validtype

! ---- error codes
public :: nfu_ediffdimsize,nfu_elimunlim

public :: cdfe
! ==== end of public interface ===============================================

! ==== module constants ======================================================
integer, parameter :: nfu_ediffdimsize = -1001
integer, parameter :: nfu_elimunlim    = -1002

type validtype
   logical :: hasmax = .false.
   logical :: hasmin = .false.
!   real(kind=8) :: max=HUGE(max),min=-HUGE(min)
   real(kind=8) :: max=0,min=0
end type

! ==== interfaces for overloaded functions ===================================
interface nfu_def_var
   module procedure nfu_def_var_0
   module procedure nfu_def_var_1
end interface

interface nfu_inq_var
   module procedure inq_var_i
   module procedure inq_var_n
end interface

interface nfu_inq_dim
   module procedure inq_dim_i
   module procedure inq_dim_n
end interface

interface nfu_inq_att
   module procedure inq_att_i_n
   module procedure inq_att_n_n
   module procedure inq_att_i_i
   module procedure inq_att_n_i
end interface

interface nfu_get_rec_r8
   module procedure get_r8_rec
   module procedure get_r8_n_rec
end interface

interface nfu_get_rec_r4
   module procedure get_r4_rec
   module procedure get_r4_n_rec
end interface

interface nfu_put_rec_r8
   module procedure put_r8_rec
   module procedure put_r8_n_rec_1d
   module procedure put_r8_n_rec_2d
end interface

interface nfu_put_rec_r4
   module procedure put_r4_rec
   module procedure put_r4_n_rec
end interface

interface nfu_get_att
   module procedure get_att_text
   module procedure get_att_int
   module procedure get_att_r4
   module procedure get_att_r8
   module procedure get_att_int_1
   module procedure get_att_r4_1
   module procedure get_att_r8_1
end interface

interface nfu_put_att
   module procedure put_att_text
   module procedure put_att_int
   module procedure put_att_r8
   module procedure put_att_r4
   module procedure put_att_int_1
   module procedure put_att_r8_1
   module procedure put_att_r4_1
end interface

interface nfu_append_att
   module procedure append_att_text_i
   module procedure append_att_text_n
end interface

interface nfu_clone_dim
  module procedure clone_dim_n
  module procedure clone_dim_i
end interface

interface nfu_clone_var
   module procedure clone_var_n
   module procedure clone_var_i
end interface

interface nfu_copy_var_data
   module procedure copy_var_data_i
   module procedure copy_var_data_n
end interface

! ==== module data ==========================================================
integer :: verb = 0 ! verbosity level

#include <netcdf.inc>

contains

! ============================================================================
function inq_dim_n(ncid,name,dimlen,is_unlim,dimid) result (iret)
  integer, intent(in) :: ncid
  character(*), intent(in) :: name
  integer, optional, intent(out) :: dimid
  logical, optional, intent(out) :: is_unlim
  integer, optional, intent(out) :: dimlen
  integer :: iret

  integer :: dimid_
  __NF_TRY__(nf_inq_dimid(ncid,name,dimid_),iret,7)
  if(present(dimid)) dimid = dimid_
  __NF_TRY__(inq_dim_i(ncid,dimid_,dimlen,is_unlim),iret,7)
7 return
end function


! ============================================================================
function inq_dim_i(ncid,dimid,dimlen,is_unlim,dimname) result (iret)
  integer, intent(in) :: ncid
  integer, intent(in) :: dimid
  character(*), optional, intent(out) :: dimname
  logical     , optional, intent(out) :: is_unlim
  integer     , optional, intent(out) :: dimlen
  integer :: iret
  
  integer :: unlimdimid

  if(present(dimname)) then
     __NF_TRY__(nf_inq_dimname(ncid,dimid,dimname),iret,7)
  end if
  if(present(is_unlim)) then
     __NF_TRY__(nf_inq_unlimdim(ncid,unlimdimid),iret,7)
     is_unlim = (dimid==unlimdimid)
  end if
  if(present(dimlen))then
     __NF_TRY__(nf_inq_dimlen(ncid,dimid,dimlen),iret,7)
  end if

7 return
end function

! ============================================================================
function inq_var_n(ncid, name, id, xtype, ndims, dimids, dimlens, natts, &
     is_dim, has_records, varsize, recsize, nrec) result(iret)
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

  integer :: vid
  character(len=NF_MAX_NAME) :: vname

  __NF_TRY__(nf_inq_varid(ncid,name,vid),iret,7)
  if(present(id)) id = vid
  iret = inq_var_i(ncid,vid,vname,xtype,ndims,dimids,dimlens,natts,&
       is_dim,has_records,varsize,recsize,nrec)

7 return  
end function

! ============================================================================
function inq_var_i(ncid, vid, name, xtype, ndims, dimids, dimlens,natts, &
     is_dim, has_records, varsize, recsize, nrec) result(iret)
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

  integer :: vxtype, vndims, vdimids(NF_MAX_VAR_DIMS), vnatts
  integer :: vsize, vrecsize
  integer :: unlimdim, did, dlen, i
  character(len=NF_MAX_NAME) :: vname

  __NF_TRY__(nf_inq_var(ncid,vid,vname,vxtype,vndims,vdimids,vnatts),iret,7)
  if (present(name)) name = vname
  if (present(xtype)) xtype = vxtype
  if (present(ndims)) ndims = vndims
  if (present(dimids)) dimids(1:min(vndims,size(dimids))) = &
       vdimids(1:min(vndims,size(dimids)))
  if (present(natts)) natts = vnatts
  if (present(is_dim)) then
     is_dim = (nf_inq_dimid(ncid,vname,did)==NF_NOERR)
  endif
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  if (present(has_records)) then
     has_records = ANY(vdimids(1:vndims)==unlimdim)
  endif
  if (present(varsize).or.present(recsize).or.present(dimlens)) then
     vsize = 1; vrecsize=1
     do i = 1,vndims
        __NF_TRY__(nf_inq_dimlen(ncid,vdimids(i),dlen),iret,7)
        vsize = vsize*dlen
        if(vdimids(i)/=unlimdim) vrecsize=vrecsize*dlen
        if(present(dimlens)) dimlens(i)=dlen
     enddo
     if(present(varsize)) varsize=vsize
     if(present(recsize)) recsize=vrecsize 
  endif
  if(present(nrec)) then
     nrec=1
     if(unlimdim/=-1.and.ANY(vdimids(1:vndims)==unlimdim)) then
        __NF_TRY__(nf_inq_dimlen(ncid,unlimdim,nrec),iret,7)
     endif
  endif

7 return  
end function

! ============================================================================
function inq_att_i_n(ncid, varid, att, xtype, len, attid) result (iret)
  integer     , intent(in) :: ncid
  integer     , intent(in) :: varid
  character(*), intent(in) :: att
  integer, optional, intent(out) :: xtype
  integer, optional, intent(out) :: len
  integer, optional, intent(out) :: attid
  integer :: iret

  integer :: xtype_, len_

  __NF_TRY__(nf_inq_att(ncid,varid,att,xtype_,len_),iret,7)
  if(present(attid)) then
     __NF_TRY__(nf_inq_attid(ncid,varid,att,attid),iret,7)
  endif
  if(present(xtype)) xtype = xtype_
  if(present(len))   len   = len_
  
7 return
end function

! ============================================================================
function inq_att_n_n(ncid, var, att, xtype, len, attid) result (iret)
  integer     , intent(in) :: ncid
  character(*), intent(in) :: var
  character(*), intent(in) :: att
  integer, optional, intent(out) :: xtype
  integer, optional, intent(out) :: len
  integer, optional, intent(out) :: attid
  integer :: iret


  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,var,varid),iret,7)
  __NF_TRY__(inq_att_i_n(ncid,varid,att,xtype,len,attid),iret,7)
7 return
end function

! ============================================================================
function inq_att_i_i(ncid, varid, attid, xtype, len, name) result (iret)
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(in) :: attid
  integer, optional, intent(out) :: xtype
  integer, optional, intent(out) :: len
  character(*), optional, intent(out) :: name
  integer :: iret

  character(NF_MAX_NAME) :: name_

  __NF_TRY__(nf_inq_attname(ncid,varid,attid,name_),iret,7)
  __NF_TRY__(inq_att_i_n(ncid,varid,name_,xtype,len),iret,7)
  if(present(name)) name = name_
7 return
end function


! ============================================================================
function inq_att_n_i(ncid, var, attid, xtype, len, name) result (iret)
  integer, intent(in) :: ncid
  character(*) :: var
  integer, intent(in) :: attid
  integer, optional, intent(out) :: xtype
  integer, optional, intent(out) :: len
  character(*), optional, intent(out) :: name
  integer :: iret

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,var,varid),iret,7)
  __NF_TRY__(inq_att_i_i(ncid,varid,attid,xtype,len,name),iret,7)
7 return
end function

! ============================================================================
function nfu_def_dim(ncid,name,size,xtype,long_name,units,bounds,dimid,varid) &
     result (iret)
  integer         , intent(in) :: ncid  ! id of NetCDF file to create 
  character(len=*), intent(in) :: name  ! name of the dimension
  integer         , intent(in) :: size  ! size of the dimension
  integer,optional, intent(in) :: xtype ! external type of the associated variable
  character(len=*), intent(in), optional :: &
       long_name, &
       units,     &
       bounds
  integer,optional,intent(out) :: dimid,varid
  integer :: iret

  integer :: did,vid

  iret=nf_redef(ncid)

  did = -1; vid = -1
  __NF_TRY__(nf_def_dim(ncid,name,size,did),iret,7)
  if(present(xtype)) then
     __NF_TRY__(nf_def_var(ncid,name,xtype,1,(/did/),vid),iret,7)
     if (present(long_name)) then
        __NF_TRY__(nf_put_att_text(ncid,vid,'long_name',len(long_name),long_name),iret,7)
     endif
     if (present(units)) then
        __NF_TRY__(nf_put_att_text(ncid,vid,'units',len(units),units),iret,7)
     endif
     if (present(bounds)) then
        __NF_TRY__(nf_put_att_text(ncid,vid,'bounds',len(bounds),bounds),iret,7)
     endif
  endif
  if(present(dimid))dimid=did
  if(present(varid))varid=vid
7 return
end function

! ============================================================================
! define variable using its symbolic name and symbolic names of the dimensions
function nfu_def_var_0(ncid,name,xtype,dimc,dimv) result(iret)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name       ! name of the variable
  integer         , intent(in) :: xtype      ! external type of the var
  integer         , intent(in) :: dimc       ! number of dimensions
  character(len=*), intent(in) :: dimv(dimc) ! vector of dimension names 
  integer :: iret

  integer :: i,dimids(NF_MAX_VAR_DIMS),varid
  do i=1,dimc
     __NF_TRY__(nf_inq_dimid(ncid,dimv(i),dimids(i)),iret,7)
  enddo
  __NF_TRY__(nf_def_var(ncid,name,xtype,dimc,dimids,varid),iret,7)

7 return
end function


function nfu_def_var_1(ncid,name,xtype,dimv) result(iret)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name       ! name of the variable
  integer         , intent(in) :: xtype      ! external type of the var
  character(len=*), intent(in) :: dimv(:)    ! vector of dimension names 
  integer :: iret

  integer :: i,dimids(NF_MAX_VAR_DIMS),varid
  integer :: dimc

  dimc = size(dimv)
  do i=1,dimc
     __NF_TRY__(nf_inq_dimid(ncid,dimv(i),dimids(i)),iret,7)
  enddo
  __NF_TRY__(nf_def_var(ncid,name,xtype,dimc,dimids,varid),iret,7)

7 return
end function

! ===========================================================================
function clone_var_n(incid,vname,oncid,name) result(iret)
  integer :: iret ! return value
  integer      , intent(in) :: incid, oncid
  character*(*), intent(in) :: vname
  character*(*), intent(in), optional :: name

  integer :: varid
  __NF_TRY__(nf_inq_varid(incid,vname,varid),iret,7)
  __NF_TRY__(clone_var_i(incid,varid,oncid,name),iret,7)
7 return
end function

! ===========================================================================
function clone_var_i(incid,ivarid,oncid,name) result(iret)
  integer :: iret ! return value
  integer      , intent(in) :: incid, oncid
  integer      , intent(in) :: ivarid
  character*(*), intent(in), optional :: name
  
  character(len=NF_MAX_NAME) :: vname ! name of the var to be cloned
  character(len=NF_MAX_NAME) :: oname ! resulting name of the variable
  character(len=NF_MAX_NAME) :: str   ! name of the dimension or attribute 
  integer :: xtype,ndims,dimids(NF_MAX_VAR_DIMS),natts,i,ovarid

  ! get the info about variable to be cloned
  __NF_TRY__(nf_inq_var(incid,ivarid,vname,xtype,ndims,dimids,natts),iret,7)
  ! get respective dimension ids in output file
  do i = 1,ndims
     __NF_TRY__(nf_inq_dimname(incid,dimids(i),str),iret,7)
     __NF_TRY__(nf_inq_dimid(oncid,str,dimids(i)),iret,7)
  enddo
  ! set the name of the variable clone
  if (present(name)) then
     oname = name
  else
     oname = vname
  end if
  ! define clone variable
  __NF_TRY__(nf_def_var(oncid,oname,xtype,ndims,dimids,ovarid),iret,7)
  ! copy all attributes
  do i = 1,natts
     __NF_TRY__(nf_inq_attname(incid,ivarid,i,str),iret,7)
     __NF_TRY__(nf_copy_att(incid,ivarid,str,oncid,ovarid),iret,7)
  enddo
7 return
  ! note taht according to docs allocatable array is supposed to
  ! be deallocated on exit from procedure
end function

! ===========================================================================
function clone_dim_n(incid,iname,oncid,name) result(iret)
  integer :: iret
  integer, intent(in) :: incid,oncid
  character*(*),intent(in) :: iname
  character*(*),intent(in),optional :: name
  
  integer :: dimid
  
  __NF_TRY__(nf_inq_dimid(incid,iname,dimid),iret,7)
  iret = clone_dim_i(incid,dimid,oncid,name)
7 return  
end function

! ===========================================================================
function clone_dim_i(incid,dimid,oncid,name) result(iret)
  integer :: iret
  integer, intent(in) :: incid,dimid,oncid
  character*(*), intent(in), optional :: name ! name of the dimension copy
    
  character(len=NF_MAX_NAME) :: dname
  integer :: unlimid,len,newdimid,newlen,newunlimid

  ! get the name of the dimension
  __NF_TRY__(nf_inq_dim(incid,dimid,dname,len),iret,7)
  ! get the id of unlimited dimension in source file
  __NF_TRY__(nf_inq_unlimdim(incid,unlimid),iret,7)
  if(dimid==unlimid) len = NF_UNLIMITED
  ! if new output name is specified, replace the name of the input dim
  ! with it
  if(present(name)) dname = name
  ! check if the desired dimension already exists in the dest file
  if(nf_inq_dimid(oncid,dname,newdimid)==NF_NOERR) then
    ! dimension already exists in output file, make sure the size is the same,
    ! or both are unlimited dimensions
    __NF_TRY__(nf_inq_dimlen(oncid,newdimid,newlen),iret,7)
    __NF_TRY__(nf_inq_unlimdim(oncid,newunlimid),iret,7)
    
    if((dimid==unlimid).neqv.(newdimid==newunlimid)) then
       __NF_TRY__(nfu_elimunlim,iret,7)
    endif
    if(dimid/=unlimid.and.len/=newlen) then
       __NF_TRY__(nfu_ediffdimsize,iret,7)
    endif
  else
    __NF_TRY__(nf_def_dim(oncid,dname,len,newdimid),iret,7)
  endif
7 return  
end function

! ============================================================================
function copy_var_data_n(ncid,varname,ncid1,name) result (iret)
  integer :: iret
  integer     , intent(in) :: ncid
  character(*), intent(in) :: varname
  integer     , intent(in) :: ncid1
  character(*), intent(in), optional :: name ! name of output variable
  
  integer :: varid, varid1
  
  __NF_TRY__(nf_inq_varid(ncid,varname,varid),iret,7)
  if(present(name)) then
     __NF_TRY__(nf_inq_varid(ncid1,name,varid1),iret,7)
  else
     __NF_TRY__(nf_inq_varid(ncid1,varname,varid1),iret,7)
  endif
  
  __NF_TRY__(copy_var_data_i(ncid,varid,ncid1,varid1),iret,7)
7 continue
end function
  

function copy_var_data_i(ncid,varid,ncid1,varid1) result (iret)
  integer :: iret
  integer, intent(in) :: ncid,  varid
  integer, intent(in) :: ncid1, varid1
  
  integer :: recsize,nrec,rec
  real(kind=8), allocatable :: buffer(:)
  
  __NF_TRY__(nfu_inq_var(ncid,varid,recsize=recsize,nrec=nrec),iret,7)
  allocate(buffer(recsize))
  
  if(ncid/=ncid1.or.varid/=varid1) then
    do rec=1,nrec
       __NF_TRY__(nfu_get_rec_r8(ncid,varid,rec,buffer),iret,7)
       __NF_TRY__(nfu_put_rec_r8(ncid1,varid1,rec,buffer),iret,7)
    enddo
  endif

7 continue
  if (allocated(buffer)) deallocate(buffer)
end function


! ============================================================================
function nfu_get_var_r8(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer, intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in)  :: name
  real(kind=8), intent(out) :: var(*) ! storage for the variable     [out]

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_get_var_double(ncid,varid,var)
7 return
end function

function nfu_get_var_r4(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer, intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in)  :: name
  real(kind=4), intent(out) :: var(*) ! storage for the variable     [out]

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_get_var_real(ncid,varid,var)
7 return
end function

function nfu_get_var_int(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer, intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in)  :: name
  integer     , intent(out) :: var(*) ! storage for the variable

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_get_var_int(ncid,varid,var)
7 return
end function

! ============================================================================
function nfu_put_var_r8(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer     , intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in) :: name
  real(kind=8), intent(in) :: var(*)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_put_var_double(ncid,varid,var)
7 return
end function

function nfu_put_var_r4(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer     , intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in) :: name
  real(kind=4), intent(in) :: var(*)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_put_var_real(ncid,varid,var)
7 return
end function

function nfu_put_var_int(ncid,name,var) result(iret)
  integer :: iret ! return value
  integer     , intent(in) :: ncid   ! id of netcdf data set
  character(*), intent(in) :: name
  integer,      intent(in) :: var(*)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = nf_put_var_int(ncid,varid,var)
7 return
end function

! ============================================================================
function get_r8_n_rec(ncid, name, rec, data) result(iret)
  integer :: iret
  integer     , intent(in)  :: ncid
  character(*), intent(in)  :: name
  integer     , intent(in)  :: rec
  real(kind=8), intent(out) :: data(*)

  integer :: varid
  
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret=get_r8_rec(ncid, varid, rec, data)
7 return

end function

! ============================================================================
function get_r4_n_rec(ncid, name, rec, data) result(iret)
  integer :: iret
  integer     , intent(in)  :: ncid
  character(*), intent(in)  :: name
  integer     , intent(in)  :: rec
  real(kind=4), intent(out) :: data(*)

  integer :: varid
  
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret=get_r4_rec(ncid, varid, rec, data)
7 return

end function

! ============================================================================
! reads double precision record from the file
function get_r8_rec(ncid, varid, rec, var) result(iret)
  integer :: iret
  integer, intent(in) :: ncid   ! id of netcdf data set
  integer, intent(in) :: varid  ! id of the variable           [in]
  integer, intent(in) :: rec       ! number of the record to get  [in]
  real(kind=8), intent(out) :: var(*) ! storage for the variable     [out]
 
  integer :: dimids(NF_MAX_VAR_DIMS), ndims, unlimdim
  integer :: start(NF_MAX_VAR_DIMS)
  integer :: count(NF_MAX_VAR_DIMS)
  integer :: i
      
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  __NF_TRY__(nf_inq_varndims(ncid,varid,ndims),iret,7)
  __NF_TRY__(nf_inq_vardimid(ncid,varid,dimids),iret,7)

  do i = 1, ndims
     if (dimids(i).eq.unlimdim) then
        start(i) = rec
        count(i) = 1
     else
        start(i) = 1
        __NF_TRY__(nf_inq_dimlen(ncid,dimids(i),count(i)),iret,7)
     endif
     ! write(*,*) i, dimids(i), start(i), count(i)
  enddo
  iret = nf_get_vara_double(ncid,varid,start,count,var)

7 return
end function

! ============================================================================
! reads single precision record from the file
function get_r4_rec(ncid, varid, rec, var) result(iret)
  integer :: iret
  integer, intent(in) :: ncid   ! id of netcdf data set
  integer, intent(in) :: varid  ! id of the variable           [in]
  integer, intent(in) :: rec       ! number of the record to get  [in]
  real(kind=4), intent(out) :: var(*) ! storage for the variable     [out]
 
  integer :: dimids(NF_MAX_VAR_DIMS), ndims, unlimdim
  integer :: start(NF_MAX_VAR_DIMS)
  integer :: count(NF_MAX_VAR_DIMS)
  integer :: i
      
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  __NF_TRY__(nf_inq_varndims(ncid,varid,ndims),iret,7)
  __NF_TRY__(nf_inq_vardimid(ncid,varid,dimids),iret,7)

  do i = 1, ndims
     if (dimids(i).eq.unlimdim) then
        start(i) = rec
        count(i) = 1
     else
        start(i) = 1
        __NF_TRY__(nf_inq_dimlen(ncid,dimids(i),count(i)),iret,7)
     endif
     ! write(*,*) i, dimids(i), start(i), count(i)
  enddo
  iret = nf_get_vara_real(ncid,varid,start,count,var)

7 return
end function

! ============================================================================
function put_r4_n_rec(ncid, name, rec, data) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(*), intent(in) :: name
  integer, intent(in) :: rec
  real(4), intent(in) :: data(*)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = put_r4_rec(ncid, varid, rec, data)
7 return

end function

! ============================================================================
function put_r8_n_rec_2d(ncid, name, rec, data) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(*), intent(in) :: name
  integer, intent(in) :: rec
  real(8), intent(in) :: data(:,:)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = put_r8_rec(ncid, varid, rec, data)
7 return

end function

! ============================================================================
function put_r8_n_rec_1d(ncid, name, rec, data) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(*), intent(in) :: name
  integer, intent(in) :: rec
  real(8), intent(in) :: data(*)

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = put_r8_rec(ncid, varid, rec, data)
7 return

end function

! ============================================================================
function put_r4_rec(ncid, varid, rec, data) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(in) :: rec
  real(4), intent(in) :: data(*)
      
  integer :: dimids(NF_MAX_VAR_DIMS), ndims
  integer :: unlimdim, unlimlen
  integer :: start(NF_MAX_VAR_DIMS)
  integer :: count(NF_MAX_VAR_DIMS)
  integer :: i
      
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  __NF_TRY__(nf_inq_varndims(ncid,varid,ndims),iret,7)
  __NF_TRY__(nf_inq_vardimid(ncid,varid,dimids),iret,7)

  do i = 1, ndims
     if (dimids(i).eq.unlimdim) then
        start(i) = rec
        count(i) = 1
     else
        start(i) = 1
        __NF_TRY__(nf_inq_dimlen(ncid,dimids(i),count(i)),iret,7)
     endif
!          write(*,*) i, dimids(i), start(i), count(i)
  enddo
  iret = nf_put_vara_real(ncid,varid,start,count,data)
7 return
end function

! ============================================================================
function put_r8_rec(ncid, varid, rec, data) result(iret)
  integer :: iret
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(in) :: rec
  real(8), intent(in) :: data(*)
      
  integer :: dimids(NF_MAX_VAR_DIMS), ndims
  integer :: unlimdim, unlimlen
  integer :: start(NF_MAX_VAR_DIMS)
  integer :: count(NF_MAX_VAR_DIMS)
  integer :: i
      
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  __NF_TRY__(nf_inq_varndims(ncid,varid,ndims),iret,7)
  __NF_TRY__(nf_inq_vardimid(ncid,varid,dimids),iret,7)

  do i = 1, ndims
     if (dimids(i).eq.unlimdim) then
        start(i) = rec
        count(i) = 1
     else
        start(i) = 1
        __NF_TRY__(nf_inq_dimlen(ncid,dimids(i),count(i)),iret,7)
     endif
!          write(*,*) i, dimids(i), start(i), count(i)
  enddo
  iret = nf_put_vara_double(ncid,varid,start,count,data)
7 return
end function

! ===========================================================================
function get_att_text(ncid, name, att, text) result(iret)
  integer :: iret
  integer, intent(in) :: ncid         ! id of the NetCDF file [in]
  character(*),intent(in) :: name   ! name of the variable  [in]
  character(*),intent(in) :: att    ! attribute name [in]
  character(*),intent(out) :: text   ! attribute value [out]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_get_att_text(ncid, varid, att, text)
      
7 return
end function

! ===========================================================================
function get_att_int(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  integer     ,intent(inout) :: d(:)   ! attribute value [out]
  
  integer :: varid,len
  integer, allocatable :: data(:)

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(nf_inq_attlen(ncid,varid,att,len),iret,7)
  allocate(data(len))
  iret = nf_get_att_int(ncid, varid, att, data)
  d(1:min(size(d),size(data))) = data(1:min(size(d),size(data)))
7 return
end function

! ===========================================================================
function get_att_r4(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  real(4)     ,intent(inout) :: d(:)   ! attribute value [out]
  
  integer :: varid,len
  real(4), allocatable :: data(:)

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(nf_inq_attlen(ncid,varid,att,len),iret,7)
  allocate(data(len))
  iret = nf_get_att_real(ncid, varid, att, data)
  d(1:min(size(d),size(data))) = data(1:min(size(d),size(data)))
7 return

end function

! ===========================================================================
function get_att_r8(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  real(8)     ,intent(out) :: d(:)   ! attribute value [out]
  
  integer :: varid,len
  real(8), allocatable :: data(:)

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(nf_inq_attlen(ncid,varid,att,len),iret,7)
  allocate(data(len))
  iret = nf_get_att_double(ncid, varid, att, data)
  d(1:min(size(d),size(data))) = data(1:min(size(d),size(data)))
7 return
end function


! ===========================================================================
function get_att_int_1(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  integer     ,intent(out) :: d    ! attribute value [out]

  integer :: data(1)

  iret = get_att_int(ncid,name,att,data)
  d = data(1)
end function

! ===========================================================================
function get_att_r4_1(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  real(4)     ,intent(out) :: d    ! attribute value [out]

  real(4) :: data(1)

  iret = get_att_r4(ncid,name,att,data)
  d = data(1)
end function

! ===========================================================================
function get_att_r8_1(ncid, name, att, d) result(iret)
  integer :: iret
  integer, intent(in)      :: ncid   ! id of the NetCDF file [in]
  character(*),intent(in)  :: name   ! name of the variable  [in]
  character(*),intent(in)  :: att    ! attribute name [in]
  real(8)     ,intent(out) :: d    ! attribute value [out]

  real(8) :: data(1)

  iret = get_att_r8(ncid,name,att,data)
  d = data(1)
end function


! ===========================================================================
function put_att_text(ncid,name,att,text) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  character*(*),intent(in) :: text   ! text to put in attribute [in]
  
  integer :: varid
      
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_text(ncid,varid,att,len(text),text)
      
7 return
end function

! ===========================================================================
function put_att_int(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  integer      ,intent(in) :: data(:)! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_int(ncid,varid,att,NF_INT,size(data),data)
      
7 return
end function

! ===========================================================================
function put_att_r4(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  real(4)      ,intent(in) :: data(:)! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_real(ncid,varid,att,NF_REAL,size(data),data)
      
7 return
end function

! ===========================================================================
function put_att_r8(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  real(8)      ,intent(in) :: data(:)! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_double(ncid,varid,att,NF_DOUBLE,size(data),data)
      
7 return
end function

! ===========================================================================
function put_att_int_1(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  integer      ,intent(in) :: data   ! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_int(ncid,varid,att,NF_INT,1,data)
      
7 return
end function

! ===========================================================================
function put_att_r4_1(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  real(4)      ,intent(in) :: data   ! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_real(ncid,varid,att,NF_REAL,1,data)
      
7 return
end function

! ===========================================================================
function put_att_r8_1(ncid,name,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: name   ! name of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  real(8)      ,intent(in) :: data   ! text to put in attribute [in]
  
  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = nf_put_att_double(ncid,varid,att,NF_DOUBLE,1,data)
      
7 return
end function

! ===========================================================================
! appends specified text to the attribute
function append_att_text_n(ncid,var,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  character*(*),intent(in) :: var  ! id of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  character*(*),intent(in) :: data   ! text to put in attribute [in]

  integer :: varid

  __NF_TRY__(nf_inq_varid(ncid,var,varid),iret,7)
  __NF_TRY__(append_att_text_i(ncid,varid,att,data),iret,7)
7 return
end function

! ===========================================================================
! appends specified text to the attribute
function append_att_text_i(ncid,varid,att,data) result(iret)
  integer :: iret
  integer      ,intent(in) :: ncid   ! id of the NetCDF file [in]
  integer      ,intent(in) :: varid  ! id of the variable  [in]
  character*(*),intent(in) :: att    ! attribute name [in]
  character*(*),intent(in) :: data   ! text to put in attribute [in]

  integer :: i,n ! original length of the attribute
  character, allocatable :: text(:)

  if(nf_inq_attlen(ncid,varid,att,n)/=NF_NOERR)then
     iret = nf_put_att_text(ncid,varid,att,len(data),data)
  else
     allocate(text(n+len(data)))
     __NF_TRY__(nf_get_att_text(ncid,varid,att,text),iret,7)
     ! erase trailing zero byte (normally appended by C), if necessary
     if(text(n)==char(0))n = n-1
     ! is there a better way to copy string to array of chars?
     ! or, even better, to allocate a string of specified length?
     do i = 1,len(data)
        text(n+i) = data(i:i)
     enddo
     __NF_TRY__(nf_put_att_text(ncid,varid,att,n+len(data),text),iret,7)
  endif
7 return
end function

! ===========================================================================
function nfu_copy_att(incid,iname,oncid,oname,aname) result(iret)
  integer :: iret
  integer, intent(in) :: incid,oncid
  character*(*), intent(in) :: iname,oname
  character*(*), intent(in) :: aname

  integer :: ivarid,ovarid

  __NF_TRY__(nf_inq_varid(incid,iname,ivarid),iret,7)
  __NF_TRY__(nf_inq_varid(oncid,oname,ovarid),iret,7)
  __NF_TRY__(nf_copy_att(incid,ivarid,aname,oncid,ovarid),iret,7)
7 return
end function

! ========================================================================
! based on presence/absence of attributes, defines valid range or missing 
! value. For details, see section 8.1 of NetCDF User Guide
function nfu_get_valid_range(ncid, name, v) result (iret)
  integer     , intent(in) :: ncid
  character(*), intent(in) :: name
  type(validtype),intent(out) :: v ! validator

  integer :: iret
  
  integer :: var_T, valid_T, scale_T, T ! types variable and of attributes
  real(kind=8) :: scale, offset, fill, r(2)
  
  ! find the type of the variable
  __NF_TRY__(nfu_inq_var(ncid,name,xtype=var_T),iret,7)

  ! find the widest type of scale and offset; note that the code
  ! uses assumption that NetCDF types are arranged in th order of rank,
  ! that is NF_BYTE < NF_CHAR < NF_SHORT < NF_INT < NF_FLOAT < NF_DOUBLE
  scale = 1; offset = 0;
  scale_T = 0
  if(nfu_inq_att(ncid,name,'scale_factor',xtype=T)==NF_NOERR) then
     __NF_TRY__(nfu_get_att(ncid,name,'scale_factor',scale),iret,7)
     scale_T = T
  endif
  if(nfu_inq_att(ncid,name,'add_offset',xtype=T)==NF_NOERR) then
     __NF_TRY__(nfu_get_att(ncid,name,'add_offset',offset),iret,7)
     scale_T = max(scale_T,T)
  endif
     
  ! examine possible range attributes
  valid_T = 0; v%hasmax=.false. ; v%hasmin=.false.
  if (nfu_inq_att(ncid,name,'valid_range',xtype=T)==NF_NOERR) then
     __NF_TRY__(nfu_get_att(ncid,name,'valid_range',r),iret,7)
     v%min = r(1)      ; v%max = r(2)
     v%hasmax = .true. ; v%hasmin = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,name,'valid_max',xtype=T)==NF_NOERR) then
     __NF_TRY__(nfu_get_att(ncid,name,'valid_max',v%max),iret,7)
     v%hasmax = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,name,'valid_min',xtype=T)==NF_NOERR) then
     __NF_TRY__(nfu_get_att(ncid,name,'valid_min',v%min),iret,7)
     v%hasmin = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,name,'missing_value',xtype=T)==NF_NOERR) then
     ! here we always scale, since missing_value is supposed to be in 
     ! external representation
     __NF_TRY__(nfu_get_att(ncid,name,'missing_value',v%min),iret,7)
     v%min = v%min*scale + offset
  else
     ! as a last resort, define range based on _FillValue
     ! get fill value and its type: from var, from file, or default
     if(nfu_get_att(ncid,name,'_FillValue',fill)/=NF_NOERR) then
        if(nf_get_att_double(ncid,NF_GLOBAL,'_FillValue',fill)/=NF_NOERR) then
           select case(var_T)
           case(NF_CHAR)
              fill = NF_FILL_CHAR
           case(NF_BYTE)
              fill = NF_FILL_BYTE
           case(NF_SHORT)
              fill = NF_FILL_SHORT
           case(NF_INT)
              fill = NF_FILL_INT
           case(NF_REAL)
              fill = NF_FILL_REAL
           case(NF_DOUBLE)
              fill = NF_FILL_DOUBLE
           end select
        endif
     endif
     if(fill>0) then
        ! if _FillValue is positive, then it defines valid maximum
        v%hasmax = .true.
        v%max = fill
        select case(T)
        case (NF_BYTE,NF_CHAR,NF_SHORT,NF_INT)
           v%max = v%max-1
        case (NF_FLOAT)
           v%max = nearest(nearest(real(v%max,4),-1.0),-1.0)
        case (NF_DOUBLE)
           v%max = nearest(nearest(real(v%max,8),-1.0),-1.0)
        end select
     else
        ! if _FillValue is negative or zero, then it defines valid minimum
        v%hasmin = .true.
        v%min = fill
        select case(T)
        case (NF_BYTE,NF_CHAR,NF_SHORT,NF_INT)
           v%min = v%min+1
        case (NF_FLOAT)
           v%min = nearest(nearest(real(v%min,4),+1.0),+1.0)
        case (NF_DOUBLE)
           v%min = nearest(nearest(real(v%min,8),+1.0),+1.0)
        end select
     endif
     ! NOTE: if we go through _FillValue branch, valid_T is 0, so values
     ! are always scaled, as it should be because _FillValue is in external 
     ! representation
  endif
  ! If valid_range is the same type as scale_factor (actually the wider of
  ! scale_factor and add_offset) and this is wider than the external data, then it
  ! will be interpreted as being in the units of the internal (unpacked) data.
  ! Otherwise it is in the units of the external (packed) data.
  if(.not.((valid_T == scale_T).and.(scale_T>var_T))) then
     v%min = v%min*scale + offset
     v%max = v%max*scale + offset
  endif
7 return
end function
   
! ========================================================================
elemental function is_valid(x, v) result (lret)
  real           , intent(in) :: x ! real value to be examined
  type(validtype), intent(in) :: v ! validator
  logical :: lret

!  if (x is NaN) then
!     lret = .false.
!  else 
  if (v%hasmin.or.v%hasmax) then
     lret = ((.not.v%hasmin).or.v%min<=x).and.((.not.v%hasmax).or.x<=v%max)
  else
     lret = x/=v%min
  endif
end function

! ===========================================================================
subroutine nfu_check_err(code,file,line)
  integer      , intent(in) :: code  ! error code
  character*(*), intent(in) :: file  ! file name
  integer      , intent(in) :: line  ! line number

  call cdfe(code,file,line,1)
  if(code/=NF_NOERR) call exit(code)

end subroutine

! ===========================================================================
! prints error message
subroutine cdfe(code, file, line, verb)
  integer      , intent(in) :: code  ! error code
  character*(*), intent(in) :: file  ! file name
  integer      , intent(in) :: line  ! line number
  integer      , intent(in) :: verb  ! verbosity level
      
  if((code/=NF_NOERR).and.(verb>0)) then
     write(*,'("File ",a," ; Line ",i5," # ",a)') file,line,nfu_strerror(code)
!     write(*,*) file,line,nfu_strerror(code)
  endif
end subroutine

function nfu_strerror(code) result(string)
  character*(80) :: string
  integer, intent(in) :: code

  select case (code)
  case (nfu_ediffdimsize)
     string = 'dimension already exists and has different size'
  case (nfu_elimunlim)
     string = 'dimension already exists and has different UNLIMITED status'
  case default
     string = nf_strerror(code)
  end select
end function


end module nfu_mod

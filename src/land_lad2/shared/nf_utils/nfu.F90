module nfu_mod

implicit none
private

! ==== public interfaces =====================================================
public :: nfu_inq_dim, nfu_inq_var, nfu_inq_att
public :: nfu_def_dim, nfu_def_var
public :: nfu_put_att
public :: nfu_get_dim, nfu_get_dim_bounds
public :: nfu_put_var, nfu_put_rec
public :: nfu_get_var, nfu_get_rec

public :: nfu_get_valid_range, nfu_is_valid, nfu_validtype, nfu_validtype2ascii
! ==== end of public interfaces ==============================================

#define __INTERFACE_SECTION__

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
interface nfu_inq_var
   module procedure inq_var_i
   module procedure inq_var_n
end interface
interface nfu_def_dim
   module procedure def_dim_0
   module procedure def_dim_r
   module procedure def_dim_i
end interface
interface nfu_def_var
   module procedure def_var_i, def_var_n, def_var_scalar
end interface
interface nfu_put_att
   module procedure put_att_text_i
   module procedure put_att_text_n
   module procedure put_att_int_i
   module procedure put_att_int_n
end interface

#define F90_TYPE integer
#define NF_TYPE  int
#include "getput.inc"

#define F90_TYPE real(8)
#define NF_TYPE  double
#include "getput.inc"

interface nfu_get_valid_range
   module procedure get_valid_range_i
   module procedure get_valid_range_n
end interface
interface nfu_is_valid
   module procedure nfu_is_valid_i
   module procedure nfu_is_valid_r
end interface
#undef __INTERFACE_SECTION__
! ---- module constants ------------------------------------------------------
character(len=*), parameter :: &
     module_name = 'nf_utils_mod', &
     version     = '$Id: nfu.F90,v 17.0 2009/07/21 03:02:52 fms Exp $', &
     tagname     = '$Name: siena_201207 $'

! ---- module types ----------------------------------------------------------
type nfu_validtype
   private
   logical :: hasmax = .false.
   logical :: hasmin = .false.
!   real(kind=8) :: max=HUGE(max),min=-HUGE(min)
   real(kind=8) :: max=0,min=0
end type

! ---- module variables ------------------------------------------------------
logical :: module_is_initialized =.FALSE.

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_TRY__(err_code,iret,LABEL)iret=err_code;if(iret/=NF_NOERR)goto LABEL
contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#define __BODY_SECTION__
! ============================================================================
function inq_dim_n(ncid, name, len, dimid) result (iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(*),intent(in) :: name
  integer, intent(out), optional :: len
  integer, intent(out), optional :: dimid

  integer :: id

  __NF_TRY__(nf_inq_dimid(ncid,name, id),iret,7)
  if(present(dimid))dimid = id
  if(present(len)) &
       iret = nf_inq_dimlen(ncid,id,len)
7 return
end function

! ============================================================================
function inq_dim_i(ncid, id, name, len) result (iret)
  integer :: iret
  integer, intent(in) :: ncid
  integer, intent(in) :: id
  character(*), intent(out), optional :: name
  integer     , intent(out), optional :: len

  if(present(name)) then
     __NF_TRY__(nf_inq_dimname(ncid,id,name),iret,7)
  endif
  if(present(len)) then
     __NF_TRY__(nf_inq_dimlen(ncid,id,len),iret,7)
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
function def_dim_0(ncid,name,size,xtype,long_name,units,edges,dimid,varid) &
     result (iret)
  integer         , intent(in) :: ncid  ! id of NetCDF file to create 
  character(len=*), intent(in) :: name  ! name of the dimension
  integer         , intent(in) :: size  ! size of the dimension
  integer,optional, intent(in) :: xtype ! external type of the associated variable
  character(len=*), intent(in), optional :: &
       long_name, &
       units,     &
       edges
  integer,optional,intent(out) :: dimid,varid
  integer :: iret

  integer :: did,vid

  iret = nf_redef(ncid)

  did = -1; vid = -1
  __NF_TRY__(nf_def_dim(ncid,name,size,did),iret,7)
  if(present(xtype)) then
     __NF_TRY__(nf_def_var(ncid,name,xtype,1,(/did/),vid),iret,7)
     if (present(long_name)) then
        __NF_TRY__(nfu_put_att(ncid,vid,'long_name',long_name),iret,7)
     endif
     if (present(units)) then
        __NF_TRY__(nfu_put_att(ncid,vid,'units',units),iret,7)
     endif
     if (present(edges)) then
        __NF_TRY__(nfu_put_att(ncid,vid,'edges',edges),iret,7)
     endif
  endif
  if(present(dimid))dimid=did
  if(present(varid))varid=vid
7 return
end function

! ============================================================================
function def_dim_r(ncid,name,data,long_name,units,edges,dimid,varid) result (iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  real            , intent(in) :: data(:)
  character(len=*), intent(in), optional :: long_name, units, edges
  integer,optional,intent(out) :: dimid,varid
  
  integer :: vid
  iret = nf_redef(ncid)
  
  __NF_TRY__(def_dim_0(ncid,name,size(data),NF_DOUBLE,long_name,units,edges,dimid,varid=vid),iret,7)
  iret = nf_enddef(ncid)
  iret = nf_put_var_double(ncid,vid,data)
  if(present(varid)) varid = vid
7 return
end function

! ============================================================================
function def_dim_i(ncid,name,data,long_name,units,edges,dimid,varid) result (iret)
  integer :: iret
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  integer         , intent(in) :: data(:)
  character(len=*), intent(in), optional :: long_name, units, edges
  integer,optional,intent(out) :: dimid,varid
  
  integer :: vid
  iret = nf_redef(ncid)
  
  __NF_TRY__(def_dim_0(ncid,name,size(data),NF_INT,long_name,units,edges,dimid,varid=vid),iret,7)
  iret = nf_enddef(ncid)
  iret = nf_put_var_int(ncid,vid,data)
  if(present(varid)) varid = vid
7 return
end function

! ============================================================================
function def_var_n(ncid,name,xtype,dims,long_name,units,varid) result(iret)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name       ! name of the variable
  integer         , intent(in) :: xtype      ! external type of the var
  character(len=*), intent(in) :: &
       dims(:)       ! vector of dimension names 
  character(len=*), intent(in), optional :: &
       long_name, &  ! name of the variable
       units         ! name of the variable
  integer         , intent(out), optional :: &
       varid   ! ID of the defined variable
  integer :: iret

  ! ---- local vars
  integer :: dimc,dimids(NF_MAX_VAR_DIMS)
  integer :: i

  dimc = size(dims)
  do i = 1,dimc
     __NF_TRY__(nf_inq_dimid(ncid,dims(i),dimids(i)),iret,7)
  enddo
  iret=def_var_i(ncid,name,xtype,dimids(1:dimc),long_name,units,varid)

7 return
end function

! ============================================================================
function def_var_scalar(ncid,name,xtype,long_name,units,varid) result(iret)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name       ! name of the variable
  integer         , intent(in) :: xtype      ! external type of the var
  character(len=*), intent(in), optional :: &
       long_name, &  ! name of the variable
       units         ! name of the variable
  integer         , intent(out), optional :: &
       varid   ! ID of the defined variable
  integer :: iret

  ! ---- local vars
  integer :: varid_

  iret = nf_redef(ncid); ! ignore errors here since file can be in define mode already
  __NF_TRY__(nf_def_var(ncid,name,xtype,0,(/1/),varid_),iret,7)
  if(present(varid)) varid = varid_
  if(present(long_name)) then
     __NF_TRY__(nfu_put_att(ncid,varid_,'long_name',long_name),iret,7)
  endif
  if(present(units)) then
     __NF_TRY__(nfu_put_att(ncid,varid_,'units',units),iret,7)
  endif
  
7 return
end function

! ============================================================================
function def_var_i(ncid,name,xtype,dimids,long_name,units,varid) result(iret)
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: name       ! name of the variable
  integer         , intent(in) :: xtype      ! external type of the var
  integer         , intent(in) :: &
       dimids(:)     ! vector of dimension ids
  character(len=*), intent(in), optional :: &
       long_name, &  ! name of the variable
       units         ! name of the variable
  integer         , intent(out), optional :: &
       varid   ! ID of the defined variable
  integer :: iret

  ! ---- local vars
  integer :: dimc,varid_

  dimc = size(dimids)
  iret = nf_redef(ncid); ! ignore errors here since file can be in define mode already
  __NF_TRY__(nf_def_var(ncid,name,xtype,dimc,dimids,varid_),iret,7)
  if(present(varid)) varid = varid_
  if(present(long_name)) then
     __NF_TRY__(nfu_put_att(ncid,varid_,'long_name',long_name),iret,7)
  endif
  if(present(units)) then
     __NF_TRY__(nfu_put_att(ncid,varid_,'units',units),iret,7)
  endif
  
7 return
end function

! ============================================================================
function put_att_text_i(ncid,varid,name,text) result (iret)
  integer :: iret
  integer         , intent(in) :: ncid,varid
  character(len=*), intent(in) :: name,text
  
  iret = nf_redef(ncid)
  iret = nf_put_att_text(ncid,varid,name,len(text),text)
end function

! ============================================================================
function put_att_text_n(ncid,varname,name,text) result (iret)
  integer :: iret
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: varname,name,text
  
  integer :: varid
  
  __NF_TRY__(nf_inq_varid(ncid,varname,varid),iret,7)
  iret = nf_redef(ncid)
  iret = nf_put_att_text(ncid,varid,name,len(text),text)
7 return
end function

! ============================================================================
function put_att_int_i(ncid,varid,name,value) result (iret)
  integer :: iret
  integer         , intent(in) :: ncid,varid
  character(len=*), intent(in) :: name
  integer         , intent(in) :: value
  
  iret = nf_redef(ncid)
  iret = nf_put_att_int(ncid,varid,name,NF_INT,1,value)
end function

! ============================================================================
function put_att_int_n(ncid,varname,name,value) result (iret)
  integer :: iret
  integer         , intent(in) :: ncid
  character(len=*), intent(in) :: varname,name
  integer         , intent(in) :: value
  
  integer :: varid
  
  __NF_TRY__(nf_inq_varid(ncid,varname,varid),iret,7)
  iret = nf_redef(ncid)
  iret = nf_put_att_int(ncid,varid,name,NF_INT,1,value)
7 return
end function

! ============================================================================
function nfu_get_dim(ncid, dimid, x) result(iret)
  integer, intent(in) :: ncid,dimid
  real   , intent(out) :: x(:)
  integer :: iret
  
  integer :: varid
  character(len=NF_MAX_NAME) :: name
  
  __NF_TRY__(nf_inq_dimname(ncid,dimid,name),iret,7)
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(nf_get_var_double(ncid,varid,x),iret,7)

7 return
end function

! ============================================================================
function nfu_get_dim_bounds(ncid,dimid,edges)result(iret)
  integer, intent(in) :: ncid,dimid
  real   , intent(out) :: edges(:)
  integer :: iret
    
  ! ---- local vars
  character(len=NF_MAX_NAME) :: name, edges_name
  real    :: x(size(edges)-1)
  integer :: varid, len
  
  __NF_TRY__( nf_inq_dimname(ncid,dimid,name),iret,7 )
  __NF_TRY__( nf_inq_dimlen(ncid,dimid,len),iret,7 )
  __NF_TRY__( nf_inq_varid(ncid,name,varid),iret,7 )
  edges_name = " "
  if (nf_get_att_text(ncid,varid,'edges',edges_name)==NF_NOERR) then
     __NF_TRY__(nf_inq_varid(ncid,edges_name,varid),iret,7)
     __NF_TRY__(nf_get_var_double(ncid,varid,edges),iret,7)
  else
     __NF_TRY__( nf_get_var_double(ncid,varid,x),iret,7 )
     edges(2:len) = (x(1:len-1)+x(2:len))/2
     edges(1) = x(1)-(edges(2)-x(1))
     edges(len+1) = x(len)+(x(len)-edges(len))
  endif
7 return
end function





! ============================================================================
! nfu_get_var interface
! ============================================================================
#define F90_TYPE integer
#define NF_TYPE  int
#include "getput.inc"

#define F90_TYPE real(8)
#define NF_TYPE  double
#include "getput.inc"



function get_valid_range_n(ncid, varname, v) result (iret)
  integer           , intent(in)  :: ncid
  character(*)      , intent(in)  :: varname
  type(nfu_validtype), intent(out) :: v ! validator

  integer :: iret
  integer :: varid

  __NF_TRY__(nfu_inq_var(ncid,varname,id=varid),iret,7)
  iret = get_valid_range_i(ncid, varid, v)

7 return
end function

! ========================================================================
! based on presence/absence of attributes, defines valid range or missing 
! value. For details, see section 8.1 of NetCDF User Guide
function get_valid_range_i(ncid, varid, v) result (iret)
  integer           , intent(in)  :: ncid
  integer           , intent(in)  :: varid
  type(nfu_validtype), intent(out) :: v ! validator

  integer :: iret
  
  integer :: var_T, valid_T, scale_T, T ! types variable and of attributes
  real(kind=8) :: scale, offset, fill, r(2)
  
  ! find the type of the variable
  __NF_TRY__(nfu_inq_var(ncid,varid,xtype=var_T),iret,7)

  ! find the widest type of scale and offset; note that the code
  ! uses assumption that NetCDF types are arranged in th order of rank,
  ! that is NF_BYTE < NF_CHAR < NF_SHORT < NF_INT < NF_FLOAT < NF_DOUBLE
  scale = 1; offset = 0;
  scale_T = 0
  if(nfu_inq_att(ncid,varid,'scale_factor',xtype=T)==NF_NOERR) then
     __NF_TRY__(nf_get_att_double(ncid,varid,'scale_factor',scale),iret,7)
     scale_T = T
  endif
  if(nfu_inq_att(ncid,varid,'add_offset',xtype=T)==NF_NOERR) then
     __NF_TRY__(nf_get_att_double(ncid,varid,'add_offset',offset),iret,7)
     scale_T = max(scale_T,T)
  endif
     
  ! examine possible range attributes
  valid_T = 0; v%hasmax=.false. ; v%hasmin=.false.
  if (nfu_inq_att(ncid,varid,'valid_range',xtype=T)==NF_NOERR) then
     __NF_TRY__(nf_get_att_double(ncid,varid,'valid_range',r),iret,7)
     v%min = r(1)      ; v%max = r(2)
     v%hasmax = .true. ; v%hasmin = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,varid,'valid_max',xtype=T)==NF_NOERR) then
     __NF_TRY__(nf_get_att_double(ncid,varid,'valid_max',v%max),iret,7)
     v%hasmax = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,varid,'valid_min',xtype=T)==NF_NOERR) then
     __NF_TRY__(nf_get_att_double(ncid,varid,'valid_min',v%min),iret,7)
     v%hasmin = .true.
     valid_T = max(valid_T,T)
  else if(nfu_inq_att(ncid,varid,'missing_value',xtype=T)==NF_NOERR) then
     ! here we always scale, since missing_value is supposed to be in 
     ! external representation
     __NF_TRY__(nf_get_att_double(ncid,varid,'missing_value',v%min),iret,7)
     v%min = v%min*scale + offset
  else
     ! as a last resort, define range based on _FillValue
     ! get fill value and its type: from var, from file, or default
     if(nf_get_att_double(ncid,varid,'_FillValue',fill)/=NF_NOERR) then
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
elemental function nfu_is_valid_r(x, v) result (lret)
  real               , intent(in) :: x ! real value to be examined
  type(nfu_validtype), intent(in) :: v ! validator
  logical :: lret

!  if (x is NaN) then
!     lret = .false.
!  else 
  if (v%hasmin.or.v%hasmax) then
     lret = .not.(((v%hasmin).and.x<v%min).or.((v%hasmax).and.x>v%max))
  else
     lret = (x /= v%min)
  endif
end function

! ========================================================================
elemental function nfu_is_valid_i(x, v) result (lret)
  integer            , intent(in) :: x ! real value to be examined
  type(nfu_validtype), intent(in) :: v ! validator
  logical :: lret
  
  lret = nfu_is_valid_r(real(x),v)
end function

! ========================================================================
function nfu_validtype2ascii(v) result (string)
  character(len=64) :: string
  type(nfu_validtype), intent(in) :: v

  if(v%hasmin.and.v%hasmax) then
     write(string,'("[",g,",",g,"]")') v%min, v%max
  else if (v%hasmin) then
     write(string,'("[",g,")")') v%min
  else if (v%hasmax) then
     write(string,'("(",g,"]")') v%max
  else 
     write(string,'("/=",g)') v%min
  endif
end function

end module nfu_mod

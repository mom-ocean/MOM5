! -*-f90-*-
! $Id: getput.inc,v 19.0 2012/01/06 20:42:09 fms Exp $

! ============================================================================
! nfu_get_var, nfu_get_rec interface and implementation
! ============================================================================

! some sanity checks
#ifndef F90_TYPE
#error F90_TYPE is not defined: must be one of FORTRAN 90 types
#endif

#ifndef NF_TYPE
#error NF_TYPE is not defined: must be netcdf type name corresponding to F90_TYPE
#endif

! macro definition for concatenation -- for construction of names based on the
! names of the operations, types, and dimension numbers
#if defined(__GFORTRAN__)

#define CONCAT1(op) op
#define CONCAT2(op,T) CONCAT1(op)T
#define CONCAT3(op,T,D) CONCAT2(op,T)D

#else

#define CONCAT3(op,T,D) op##T##D
#define CONCAT2(op,T)   op##T

#endif

! names of the functions we define
#define GET_VAR(T,D) CONCAT3(get_var_,T,D)
#define PUT_VAR(T,D) CONCAT3(put_var_,T,D)
#define GET_REC(T,D) CONCAT3(get_rec_,T,D)
#define PUT_REC(T,D) CONCAT3(put_rec_,T,D)

! define names of the corresponding netcdf functions. The two-stage definition is
! necessary because of the preprocessor argument pre-scan rules. See, for example,
! http://gcc.gnu.org/onlinedocs/cpp/Argument-Prescan.html
#define NF_GET_VAR_(T) CONCAT2(nf_get_var_,T)
#define NF_GET_VAR_T NF_GET_VAR_(NF_TYPE) 

#define NF_PUT_VAR_(T) CONCAT2(nf_put_var_,T)
#define NF_PUT_VAR_T NF_PUT_VAR_(NF_TYPE) 

#define NF_GET_VARA_(T) CONCAT2(nf_get_vara_,T)
#define NF_GET_VARA_T NF_GET_VARA_(NF_TYPE) 

#define NF_PUT_VARA_(T) CONCAT2(nf_put_vara_,T)
#define NF_PUT_VARA_T NF_PUT_VARA_(NF_TYPE) 

! #### Interface definition ###################################################
! define specific names of the subroutines
#define GET_VAR_D0 GET_VAR(NF_TYPE,D0)
#define GET_VAR_D1 GET_VAR(NF_TYPE,D1)
#define GET_VAR_D2 GET_VAR(NF_TYPE,D2)
#define GET_VAR_D3 GET_VAR(NF_TYPE,D3)

#define PUT_VAR_D0 PUT_VAR(NF_TYPE,D0)
#define PUT_VAR_D1 PUT_VAR(NF_TYPE,D1)
#define PUT_VAR_D2 PUT_VAR(NF_TYPE,D2)
#define PUT_VAR_D3 PUT_VAR(NF_TYPE,D3)

#define GET_REC_D0N GET_REC(NF_TYPE,D0N)
#define GET_REC_D1N GET_REC(NF_TYPE,D1N)
#define GET_REC_D2N GET_REC(NF_TYPE,D2N)
#define GET_REC_D3N GET_REC(NF_TYPE,D3N)

#define GET_REC_D0I GET_REC(NF_TYPE,D0I)
#define GET_REC_D1I GET_REC(NF_TYPE,D1I)
#define GET_REC_D2I GET_REC(NF_TYPE,D2I)
#define GET_REC_D3I GET_REC(NF_TYPE,D3I)

#define PUT_REC_D0N PUT_REC(NF_TYPE,D0N)
#define PUT_REC_D1N PUT_REC(NF_TYPE,D1N)
#define PUT_REC_D2N PUT_REC(NF_TYPE,D2N)
#define PUT_REC_D3N PUT_REC(NF_TYPE,D3N)

#define PUT_REC_D0I PUT_REC(NF_TYPE,D0I)
#define PUT_REC_D1I PUT_REC(NF_TYPE,D1I)
#define PUT_REC_D2I PUT_REC(NF_TYPE,D2I)
#define PUT_REC_D3I PUT_REC(NF_TYPE,D3I)

#ifdef __INTERFACE_SECTION__
! nfu_get_var interface
interface nfu_get_var
   module procedure GET_VAR_D0, GET_VAR_D1, GET_VAR_D2, GET_VAR_D3
end interface
interface nfu_put_var
   module procedure PUT_VAR_D0, PUT_VAR_D1, PUT_VAR_D2, PUT_VAR_D3
end interface
interface nfu_get_rec
   module procedure GET_REC_D0N, GET_REC_D1N, GET_REC_D2N, GET_REC_D3N
   module procedure GET_REC_D0I, GET_REC_D1I, GET_REC_D2I, GET_REC_D3I
end interface
interface nfu_put_rec
   module procedure PUT_REC_D0N, PUT_REC_D1N, PUT_REC_D2N, PUT_REC_D3N
   module procedure PUT_REC_D0I, PUT_REC_D1I, PUT_REC_D2I, PUT_REC_D3I
end interface
#endif

! #### END of interface definition ############################################


! #### Implementation definition ##############################################

#ifdef __BODY_SECTION__ 
! ============================================================================
! nfu_get_var implemenatation
! ============================================================================
function GET_VAR_D0(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(inout) :: var   ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = NF_GET_VAR_T(ncid,varid,var)
7 return
end function
! ============================================================================
function GET_VAR_D1(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE   , intent(inout) :: var(*) ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  iret = NF_GET_VAR_T(ncid,varid,var)
7 return
end function
! ============================================================================
function GET_VAR_D2(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(inout) :: var(:,:) ! storage for the variable
  integer :: iret ! return value

  iret = GET_VAR_D1(ncid,name,var)
7 return
end function
! ============================================================================
function GET_VAR_D3(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(inout) :: var(:,:,:) ! storage for the variable
  integer :: iret ! return value

  iret = GET_VAR_D1(ncid,name,var)
7 return
end function

! ============================================================================
function PUT_VAR_D0(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE   , intent(in) :: var    ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = NF_PUT_VAR_T(ncid,varid,var)
7 return
end function
! ============================================================================
function PUT_VAR_D1(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(in) :: var(*) ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = NF_PUT_VAR_T(ncid,varid,var)
7 return
end function
! ============================================================================
function PUT_VAR_D2(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(in) :: var(:,:) ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = NF_PUT_VAR_T(ncid,varid,var)
7 return
end function
! ============================================================================
function PUT_VAR_D3(ncid,name,var) result(iret)
  integer     , intent(in) :: ncid   ! id of netcdf file
  character(*), intent(in) :: name   ! name of the variable
  F90_TYPE    , intent(in) :: var(:,:,:) ! storage for the variable
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid, name, varid), iret, 7)
  iret = NF_PUT_VAR_T(ncid,varid,var)
7 return
end function

! ============================================================================
! nfu_get_rec implementation
! ============================================================================
function GET_REC_D0N(ncid,name,rec,var) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  character(*) , intent(in)    :: name   ! name of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var    ! storage for the variable
  integer :: iret ! return value

  F90_TYPE :: var1(1)
  __NF_TRY__(GET_REC_D1N(ncid,name,rec,var1),iret,7)
  var=var1(1)
7 return
end function
! ============================================================================
function GET_REC_D1N(ncid,name,rec,var,start,count) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  character(*) , intent(in)    :: name   ! name of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var(*) ! storage for the variable
  integer, optional, intent(in) :: start(:), count(:) ! slab to read
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(GET_REC_D1I(ncid,varid,rec,var,start,count),iret,7)
7 return
end function
! ============================================================================
function GET_REC_D2N(ncid,name,rec,var,start,count) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  character(*) , intent(in)    :: name   ! name of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var(:,:) ! storage for the variable
  integer, optional, intent(in) :: start(2), count(2) ! slab to read
  integer :: iret ! return value

  iret=GET_REC_D1N(ncid,name,rec,var,start,count)
end function
! ============================================================================
function GET_REC_D3N(ncid,name,rec,var,start,count) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  character(*) , intent(in)    :: name   ! name of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var(:,:,:) ! storage for the variable
  integer, optional, intent(in) :: start(3), count(3) ! slab to read
  integer :: iret ! return value

  iret=GET_REC_D1N(ncid,name,rec,var)
end function
! ============================================================================
function GET_REC_D0I(ncid,varid,rec,var,start) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  integer      , intent(in)    :: varid  ! id of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var    ! storage for the variable
  integer, optional, intent(in) :: start(:) ! slab to read
  integer :: iret ! return value

  F90_TYPE :: var1(1)
  integer :: count_(NF_MAX_VAR_DIMS)

  count_(:) = 1
  __NF_TRY__(GET_REC_D1I(ncid,varid,rec,var1,start,count_),iret,7)
  var=var1(1)
7 return
end function
! ============================================================================
function GET_REC_D1I(ncid,varid,rec,var,start,count) result(iret)
  integer      , intent(in)     :: ncid   ! id of netcdf file
  integer      , intent(in)     :: varid  ! id of the variable
  integer      , intent(in)     :: rec    ! number of the record to get
  F90_TYPE     , intent(inout)  :: var(*) ! storage for the variable
  integer      , intent(in), optional :: start(:), count(:) ! definition of
                                          ! the slab to read
  integer :: iret ! return value

  integer :: dimids(NF_MAX_VAR_DIMS), ndims, unlimdim
  integer :: start_(NF_MAX_VAR_DIMS)
  integer :: count_(NF_MAX_VAR_DIMS)
  integer :: i
      
  __NF_TRY__(nf_inq_unlimdim(ncid,unlimdim),iret,7)
  __NF_TRY__(nf_inq_varndims(ncid,varid,ndims),iret,7)
  __NF_TRY__(nf_inq_vardimid(ncid,varid,dimids),iret,7)

  do i = 1, ndims
     if (dimids(i).eq.unlimdim) then
        start_(i) = rec
        count_(i) = 1
     else
        start_(i) = 1
        __NF_TRY__(nf_inq_dimlen(ncid,dimids(i),count_(i)),iret,7)
        if (present(start)) then
           start_(i) = start(i)
           count_(i) = count_(i)-start_(i)+1
        endif
        if (present(count)) then
           count_(i) = count(i)
        endif 
     endif
     ! write(*,*) i, dimids(i), start_(i), count_(i)
  enddo
  iret = NF_GET_VARA_T(ncid,varid,start_,count_,var)

7 return
end function
! ============================================================================
function GET_REC_D2I(ncid,varid,rec,var,start,count) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  integer      , intent(in)    :: varid  ! id of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var(:,:)    ! storage for the variable
  integer      , intent(in), optional :: start(2), count(2) ! definition of
                                         ! the slab to read
  integer :: iret ! return value

  iret=GET_REC_D1I(ncid,varid,rec,var,start,count)
end function
! ============================================================================
function GET_REC_D3I(ncid,varid,rec,var,start,count) result(iret)
  integer      , intent(in)    :: ncid   ! id of netcdf file
  integer      , intent(in)    :: varid  ! id of the variable
  integer      , intent(in)    :: rec    ! number of the record to get
  F90_TYPE     , intent(inout) :: var(:,:,:)    ! storage for the variable
  integer      , intent(in), optional :: start(3), count(3) ! definition of
                                          ! the slab to read
  integer :: iret ! return value

  iret=GET_REC_D1I(ncid,varid,rec,var,start,count)
end function

! ============================================================================
! nfu_put_rec implementation
! ============================================================================
function PUT_REC_D0N(ncid,name,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  character(*) , intent(in) :: name   ! name of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var    ! data to write
  integer :: iret ! return value

  F90_TYPE :: var1(1)
  var1(1)=var
  iret = PUT_REC_D1N(ncid,name,rec,var1)
7 return
end function
! ============================================================================
function PUT_REC_D1N(ncid,name,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  character(*) , intent(in) :: name   ! name of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(*) ! data to write
  integer :: iret ! return value

  integer :: varid
  __NF_TRY__(nf_inq_varid(ncid,name,varid),iret,7)
  __NF_TRY__(PUT_REC_D1I(ncid,varid,rec,var),iret,7)
7 return
end function
! ============================================================================
function PUT_REC_D2N(ncid,name,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  character(*) , intent(in) :: name   ! name of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(:,:) ! data to write
  integer :: iret ! return value

  iret=PUT_REC_D1N(ncid,name,rec,var)
end function
! ============================================================================
function PUT_REC_D3N(ncid,name,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  character(*) , intent(in) :: name   ! name of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(:,:,:) ! data to write
  integer :: iret ! return value

  iret=PUT_REC_D1N(ncid,name,rec,var)
end function
! ============================================================================
function PUT_REC_D0I(ncid,varid,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  integer      , intent(in) :: varid  ! id of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var    ! data to write
  integer :: iret ! return value

  F90_TYPE :: var1(1)
  var1(1)=var
  iret = PUT_REC_D1I(ncid,varid,rec,var1)
end function
! ============================================================================
function PUT_REC_D1I(ncid,varid,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  integer      , intent(in) :: varid  ! id of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(*) ! data to write
  integer :: iret ! return value

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
  i = nf_enddef(ncid) ! ignore errors here (the file may be in define mode already)
  iret = NF_PUT_VARA_T(ncid,varid,start,count,var)

7 return
end function
! ============================================================================
function PUT_REC_D2I(ncid,varid,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  integer      , intent(in) :: varid  ! id of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(:,:) ! data to write
  integer :: iret ! return value

  iret=PUT_REC_D1I(ncid,varid,rec,var)
end function
! ============================================================================
function PUT_REC_D3I(ncid,varid,rec,var) result(iret)
  integer      , intent(in) :: ncid   ! id of netcdf file
  integer      , intent(in) :: varid  ! id of the variable
  integer      , intent(in) :: rec    ! number of the record to get
  F90_TYPE     , intent(in) :: var(:,:,:) ! data to write
  integer :: iret ! return value

  iret=PUT_REC_D1I(ncid,varid,rec,var)
end function
#endif
! #### End of implementation definition ######################################
#undef CONCAT3
#undef CONCAT2

#undef GET_VAR
#undef PUT_VAR
#undef GET_REC
#undef PUT_REC

#undef NF_GET_VAR_
#undef NF_GET_VAR_T

#undef NF_PUT_VAR_
#undef NF_PUT_VAR_T

#undef NF_GET_VARA_
#undef NF_GET_VARA_T

#undef GET_VAR_D0
#undef GET_VAR_D1
#undef GET_VAR_D2
#undef GET_VAR_D3

#undef PUT_VAR_D0
#undef PUT_VAR_D1
#undef PUT_VAR_D2
#undef PUT_VAR_D3

#undef GET_REC_D0N
#undef GET_REC_D1N
#undef GET_REC_D2N
#undef GET_REC_D3N

#undef GET_REC_D0I
#undef GET_REC_D1I
#undef GET_REC_D2I
#undef GET_REC_D3I

#undef PUT_REC_D0N
#undef PUT_REC_D1N
#undef PUT_REC_D2N
#undef PUT_REC_D3N

#undef PUT_REC_D0I
#undef PUT_REC_D1I
#undef PUT_REC_D2I
#undef PUT_REC_D3I

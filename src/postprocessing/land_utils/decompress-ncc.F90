!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
!
! This program reads several input netcdf files, presumably containing "compressed
! by gathering" data, and combines them into a single output file
!
!-----------------------------------------------------------------------

#define __NF_ASRT__(ierr) call nfu_check_err(ierr,__FILE__,__LINE__)
program decompress

  use nfu_mod
  use nfu_compress_mod
  implicit none
  include 'netcdf.inc'

  integer, parameter :: PATH_MAX = 1024 ! max len of the file name; 
  integer, parameter :: HEADERPAD = 16384 ! Use mpp_io large headers;
  integer            :: blksz = 65536  ! blksz must be writable for nf__create

  character(PATH_MAX), allocatable :: files(:) ! names of all files on the command line
  character(PATH_MAX)              :: outfile  ! name of the output file
  integer :: nfiles    ! number of files on command line
  integer :: debug = 0 ! debug output verbosity level
  integer, allocatable :: input(:)             ! netcdf IDs of input files
  integer :: i,iret,ncid,dimid,varid,varid1,xtype,dimlen,vsize,ndims,nvars,natts
  integer :: in_format, cmode
  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname,varname,attname
  logical :: is_dim, is_compressed, has_records, add_missing_value = .FALSE.
  real   , allocatable :: buffer(:)
  logical, allocatable :: mask(:)
  real :: missing

  ! get command line options and list of files
  call parse_command_line() ! modigies global data!

  call assert(nfiles>0,'at least one input file must be specified')
  if(debug>0) then
     do i = 1,nfiles
        write(*,'("input file",i3,":",a)')i, '"'//trim(files(i))//'"'
     enddo
     write(*,'("output file:",a)')'"'//trim(outfile)//'"'
  endif

  ! open all input files
  allocate(input(nfiles))
  do i = 1,nfiles
     __NF_ASRT__(nf_open(files(i),NF_NOWRITE,input(i)))
     __NF_ASRT__(nf_inq_format(input(i),in_format))
  enddo

  if (in_format==NF_FORMAT_NETCDF4) then
     cmode = NF_NETCDF4
  elseif (in_format==NF_FORMAT_NETCDF4_CLASSIC) then
     cmode=IOR(NF_NETCDF4,NF_CLASSIC_MODEL)
  elseif (in_format==NF_FORMAT_64BIT) then
     cmode=IOR(NF_CLOBBER,NF_64BIT_OFFSET)
     if(debug>0)write(*,'("output file is 64-bit netcdf")')
  elseif (in_format==NF_FORMAT_CLASSIC) then
     cmode=IOR(NF_CLOBBER,NF_CLASSIC_MODEL)
     if(debug>0)write(*,'("output file is 32-bit netcdf")')
  else
     call assert(.false.,'Unknown netCDF format')
  endif

  ! create output file
  __NF_ASRT__(nf__create(outfile,cmode,0,blksz,ncid))

  ! Create netcdf structure in the output NetCDF file, using last input file
  ! as a template.

  ! clone all dimensions except compressed ones; compressed are just skipped
  __NF_ASRT__(nf_inq_ndims(input(nfiles),ndims))
  do dimid = 1,ndims
     __NF_ASRT__(nfu_inq_dim(input(nfiles),dimid,dimname=dimname,dimlen=dimlen,is_unlim=has_records))
     if(nfu_inq_att(input(nfiles),dimname,'compress')==NF_NOERR) cycle
     if(has_records)&
          dimlen=NF_UNLIMITED
     if(debug>0)&
          write(*,*)'defining dimension "'//trim(dimname)//'" with length',dimlen
     __NF_ASRT__(nf_def_dim(ncid,dimname,dimlen,i)) ! i is just a space for id
  enddo

  ! clone all variable definitions, replacing compressed dimensions with sets
  ! of uncompressed ones
  __NF_ASRT__(nf_inq_nvars(input(nfiles),nvars))
  do varid = 1,nvars
     iret = nfu_inq_compressed_var(input(nfiles),varid,varname,xtype,ndims,dimids,&
          natts=natts,is_dim=is_dim, is_compressed=is_compressed)
     if(debug>0)&
          write(*,*)'defining variable "'//trim(varname)//'"'
     __NF_ASRT__(iret) ! just because the line is going to be too long
     if(is_dim.and.is_compressed) cycle
     do i = 1,ndims
        __NF_ASRT__(nf_inq_dimname(input(nfiles),dimids(i),dimname))
        __NF_ASRT__(nf_inq_dimid(ncid,dimname,dimids(i)))
     enddo
     __NF_ASRT__(nf_def_var(ncid,varname,xtype,ndims,dimids,varid1))
     do i=1,natts
        __NF_ASRT__(nf_inq_attname(input(nfiles),varid,i,attname))
        __NF_ASRT__(nf_copy_att(input(nfiles),varid,attname,ncid,varid1))
     enddo
     if(add_missing_value.and.is_compressed) then
        ! check if missing_value or _FillValue attributes are present 
        if ( nf_inq_atttype(input(nfiles),varid,'missing_value',iret)/=NF_NOERR .and. &
             nf_inq_atttype(input(nfiles),varid,'_FillValue',iret)/=NF_NOERR ) then
           ! if not, define the missing value attribute
           select case(xtype)
           case(NF_DOUBLE) 
              missing = NF_FILL_DOUBLE
           case(NF_FLOAT)
              missing = NF_FILL_FLOAT
           case(NF_INT)
              missing = NF_FILL_INT
           end select
           ! and add it to the output variable
           __NF_ASRT__(nf_put_att_double(ncid,varid1,'missing_value',xtype,1,missing))
        endif
     endif
  enddo

  ! clone all global attributes
  __NF_ASRT__(nf_inq_natts(input(nfiles),natts))
  do i = 1,natts
     __NF_ASRT__(nf_inq_attname(input(nfiles),NF_GLOBAL,i,attname))
     __NF_ASRT__(nf_copy_att(input(nfiles),NF_GLOBAL,attname,ncid,NF_GLOBAL))
  enddo

  ! ---- end of definition stage
  __NF_ASRT__(nf__enddef(ncid,HEADERPAD,4,0,4))


  ! grow unlimited dimension, if necessary
  do varid=1,nvars
     __NF_ASRT__(nfu_inq_compressed_var(input(nfiles),varid,name=varname,dimlens=dimlens,has_records=has_records))
     if(has_records) then
        ! just write an integer at the very end of the variable -- that extends the 
        ! record dimensions as well
        __NF_ASRT__(nf_inq_varid(ncid,varname,varid))
        __NF_ASRT__(nf_put_var1_int(ncid,varid,dimlens,0))
        exit ! this loop
     endif
  enddo

  ! gather and copy data
  __NF_ASRT__(nf_inq_nvars(ncid,nvars))
  do varid = 1,nvars
!     __NF_ASRT__(nfu_inq_compressed_var(ncid,varid,name=varname,varsize=vsize))
     __NF_ASRT__(nfu_inq_var(ncid,varid,varname,xtype,varsize=vsize))
     if(debug>0) &
          write(*,*)'processing var "'//trim(varname)//'"'
     allocate(buffer(vsize),mask(vsize))
     mask(:) = .false.

     ! obtain the missing value 
     if(nfu_get_att(ncid,varname,'missing_value',missing)==NF_NOERR) then
        ! do nothing, the value is already in the "missing" variable
     else if(nfu_get_att(ncid,varname,'_FillValue',missing)==NF_NOERR) then
        ! do nothing, the value is already in the "missing" variable
     else
        ! get fill value for the type instead of the missing value
        select case(xtype)
        case(NF_DOUBLE) 
           missing = NF_FILL_DOUBLE
        case(NF_FLOAT)
           missing = NF_FILL_FLOAT
        case(NF_INT)
           missing = NF_FILL_INT
        end select
     endif
     ! fill the buffer with the missing value
     buffer=missing

     ! read the variable
     do i=1,nfiles
        __NF_ASRT__(nfu_get_compressed_var_r8n(input(i),varname,buffer,mask))
     enddo
     ! write the variable
     __NF_ASRT__(nfu_put_var_r8(ncid,varname,buffer))
     deallocate(buffer,mask)
  enddo
     
  __NF_ASRT__(nf_close(ncid))

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ---- parses command line arguments, getting options and gathering list of
! file names
! NOTE: updates global variables.
subroutine parse_command_line()
  character(PATH_MAX) :: arg, param

  integer :: nargs     ! number of command-line arguments
  logical :: do_interpret_arguments
  integer :: i, iostat

  integer, external :: iargc

  nargs = iargc()
  if(nargs==0) then
     call usage()
     call exit(1)
  endif
  
  allocate(files(nargs))  ! allocate storage for all file names
  do_interpret_arguments = .true.
  add_missing_value = .false.
  i=1        ! counter of all command-line arguments
  nfiles = 0 ! counter of input files
  do while (i<=nargs)
     call getarg(i,arg)
     if(debug>1) write(*,*)'argument ',i, trim(arg)
     if(arg(1:1)=='-'.and.do_interpret_arguments) then
        select case(trim(arg))
        case('--')
           do_interpret_arguments = .false.

        case('-D','--debug-level')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           call getarg(i+1,param)
           read(param,*,iostat=iostat) debug
           call assert(iostat==0,trim(arg)//' flag must be followed by integer verbosity level')
           i=i+1

        case('-m','--add-missing-value')
	   add_missing_value = .TRUE.

        case ('-h','-?','--help')
           call usage()
           call exit(1)

        case default
           call usage()
           call assert(.false.,'argument "'//trim(arg)//'" is illegal')
        end select
     else
        ! argument is either input or output file name, add it to the list
        nfiles = nfiles+1
        files(nfiles) = arg
     endif
     i = i+1
  enddo
  if (nfiles>0) then
     outfile = files(nfiles)
     nfiles  = nfiles-1
  endif
end subroutine


! ---- prints usage information
subroutine usage()
! this program reads several input netcdf files, presumably containing "compressed
! by gathering" data, and combines them into a single output file
  character(len=PATH_MAX) :: name
  call getarg(0,name)
  write(*,'(a)')'Converts one or several compressed-by-gathering netcdf file into'
  write(*,'(a)')'one regular netcdf. Normally used to convert lm3 restarts into a'
  write(*,'(a)')'form suitable for visualization applications.'
  write(*,'(a)')
  write(*,'(a)')'Usage:'
  write(*,'(a)')'  '//trim(name)//' [-D debug-level] [-m] in.nc [...] out.nc'
  write(*,'(a)')
  write(*,'(a)')'-D debug-level   Specifies level of debug output verbosity'
  write(*,'(a)')'-m               Forces adding a missing_value attribute to the variables'
  write(*,'(a)')'                 that do not have it'
  write(*,'(a)')'in.nc            Input file name(s)'
  write(*,'(a)')'out.nc           Output file name'
  write(*,'(a)')
  write(*,'(a)')'WARNING: output file is overwritten.'
end subroutine


! ===========================================================================
! ---- prints error message an exits if condition is not satisfied
subroutine assert(cond,message)
  logical     , intent(in) :: cond    ! condition to check
  character(*), intent(in) :: message ! error message to print if condition is not satisfied
  
  if(.not.cond) then
     write(*,*) 'ERROR :: ',trim(message)
     call exit(1)
  endif
end subroutine

end program decompress

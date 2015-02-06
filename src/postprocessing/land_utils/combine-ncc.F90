!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
!
! This program reads several input netcdf files, presumably containing 
! "compressed by gathering" data, and combines them into a single output file
!
!-----------------------------------------------------------------------

#define __NF_ASRT__(ierr) call nfu_check_err(ierr,__FILE__,__LINE__)
program combine_res

  use nfu_mod
  use nfu_compress_mod
  implicit none
  include 'netcdf.inc'

  integer, parameter :: PATH_MAX = 1024 ! max len of the file name; 
  integer, parameter :: HEADERPAD = 16384 ! Use mpp_io large headers; 
  integer            ::  blksz = 65536  ! blksz must be writable for nf__create

  character(PATH_MAX), allocatable :: files(:) ! names of all files on the command line
  character(PATH_MAX)              :: outfile  ! name of the output file
  integer :: nfiles    ! number of files on command line
  integer :: debug = 0 ! debug output verbosity level
  integer, allocatable :: input(:)             ! netcdf IDs of input files
  integer :: i,ncid,dimid,varid,dimlen,vsize,ndims,nvars,ngatts
  integer :: dimlens(NF_MAX_DIMS)
  logical :: has_records
  integer :: in_format ! format of input files
  integer :: cmode     ! mode for output file creation
  character(NF_MAX_NAME) :: dimname,varname,attname
  real   , allocatable :: buffer(:)
  logical, allocatable :: mask(:)

  ! get command line options and list of files
  call parse_command_line() ! modigies global data!

  call assert(nfiles>0,'at least one input file must be specified')
  if(debug>0) then
     do i = 1,nfiles
        write(*,'("input file",i3,":",a)')i, '"'//trim(files(i))//'"'
     enddo
     write(*,'("output file:",a)')'"'//trim(outfile)//'"'
  endif

  ! open all input files and determine the creation mode of output file:
  ! if any of the input files is 64-bit then the output is 64-bit as well,
  ! otherwise it's 32-bit
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

  ! mpp_io supports environment variables to set these. For Riga, we'll simply use the defaults`

  __NF_ASRT__(nf__create(outfile,cmode,0,blksz,ncid))

  ! Create netcdf structure in the output NetCDF file, using last input file
  ! as a template.

  ! clone all dimensions; for compressed dimensions calculate the length
  __NF_ASRT__(nf_inq_ndims(input(nfiles),ndims))
  do dimid = 1,ndims
     __NF_ASRT__(nfu_inq_dim(input(nfiles),dimid,dimname=dimname,dimlen=dimlen,is_unlim=has_records))
     if(nfu_inq_att(input(nfiles),dimname,'compress')==NF_NOERR) then
        __NF_ASRT__(nfu_inq_compressed_var(input(nfiles),dimname,varsize=vsize))
        allocate(buffer(vsize),mask(vsize))
        mask(:) = .false.
        do i=1,nfiles
           __NF_ASRT__(nfu_get_compressed_var_r8n(input(i),dimname,buffer,mask))
        enddo
        dimlen = max(count(mask),1)
        ! can't have 0-length dimension, since it is (mis-)understood by netcdf as 
        ! a record one.
        deallocate(buffer,mask)
     endif
     if(debug>0)&
          write(*,*)'defining dimension "'//trim(dimname)//'" with length',dimlen
     if(has_records)then
        dimlen = NF_UNLIMITED
     endif
     __NF_ASRT__(nf_def_dim(ncid,dimname,dimlen,i)) ! i is just a space for id

  enddo

  ! clone all variable definitions
  __NF_ASRT__(nf_inq_nvars(input(nfiles),nvars))
  do i = 1,nvars
     __NF_ASRT__(nfu_clone_var(input(nfiles),i,ncid))
     ! NOTE: since cloning of variable definition relies on dimension names,
     ! each variable tile and compressed dimensions automaticaly get the right
     ! size, as defined while creating dimensions in the output file
  enddo

  ! clone all global attributes
  __NF_ASRT__(nf_inq_natts(input(nfiles),ngatts))
  do i = 1,ngatts
     __NF_ASRT__(nf_inq_attname(input(nfiles),NF_GLOBAL,i,attname))
     __NF_ASRT__(nf_copy_att(input(nfiles),NF_GLOBAL,attname,ncid,NF_GLOBAL))
  enddo

  ! ---- end of definition stage
  __NF_ASRT__(nf__enddef(ncid,HEADERPAD,4,0,4))

  ! grow unlimited dimension, if necessary
  do varid = 1,nvars
     __NF_ASRT__(nfu_inq_var(input(nfiles),varid,dimlens=dimlens,has_records=has_records))
     if(has_records)then
        __NF_ASRT__(nf_put_var1_int(ncid,varid,dimlens,0))
        exit ! loop
     endif
  enddo


  ! gather and copy data
  do varid = 1,nvars
     __NF_ASRT__(nfu_inq_compressed_var(ncid,varid,name=varname,varsize=vsize))
     if(debug>0) &
          write(*,*)'processing var "'//trim(varname)//'"'
     allocate(buffer(vsize),mask(vsize))
     mask(:) = .false.
     do i=1,nfiles
        __NF_ASRT__(nfu_get_compressed_var_r8n(input(i),varname,buffer,mask))
     enddo
     if (count(mask)>0) then
        __NF_ASRT__(nfu_put_var_r8(ncid,varname,pack(buffer,mask)))
     endif
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
  if(debug>0) &
       write(*,*) nfiles, ' input files'
end subroutine


! ---- prints usage information
subroutine usage()
  character(len=PATH_MAX) :: name
  call getarg(0,name)
  write(*,'(a)')'Combines several compressed-by-gathering netcdf files into one.'
  write(*,'(a)')'Normally used to combine lm3 restarts generated by each processor.'
  write(*,'(a)')
  write(*,'(a)')'Usage:'
  write(*,'(a)')'  '//trim(name)//' [-D debug-level] in.nc [...] out.nc'
  write(*,'(a)')
  write(*,'(a)')'-D debug-level   Specifies level of debug output verbosity'
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

end program combine_res

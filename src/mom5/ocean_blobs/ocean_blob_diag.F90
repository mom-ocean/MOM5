module ocean_blob_diag_mod
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Controls the diagnostic output from individual blobs.
!</OVERVIEW>
!
!<DESCRIPTION>
! Controls the diagnostics from individual blobs.  Blob diagnostics
! are snapshots at full E system time steps of blob properties.  There is
! no averaging of individual blob properties. The properties that can be 
! output are set in a diagnostic file.  The name of the diagnostic file 
! is specified in the namelist.
!</DESCRIPTION>
!
!<INFO>
!<NAMELIST NAME="ocean_blob_diag_nml">
!
!  <DATA NAME ="blob_diagnostics" TYPE="LOGICAL">
!    Logical as to whether diagnostics should be saved or not.
!    Default is .false.
!  </DATA>
!
!  <DATA NAME ="diag_table" TYPE="character">
!    Name of file to look for blob diagnostic information.
!    Default is "blob_diag_table"
!  </DATA>
!
!  <DATA NAME ="dump_num" TYPE="integer">
!    The number of entried to keep in memory before writing
!    them to file.  The higher the number, the more memory that
!    the module will take up, but, it should lower the frequency of
!    IO operations.
!    Default is 2000000
!  </DATA>
!
!  <DATA NAME ="frequency" TYPE="integer">
!    The frequency (in number of E system time steps) that blob
!    diagnostics should be saved.
!    Default is 1
!  </DATA>
!
!</NAMELIST>
!</INFO>
!
use constants_mod, only: c2dbars
use fms_mod,       only: stderr, stdout, stdlog, error_mesg, FATAL, WARNING
use fms_mod,       only: close_file, check_nml_error, open_namelist_file
use fms_mod,       only: mpp_error
use mpp_mod,       only: mpp_pe, mpp_sum

use ocean_blob_util_mod,  only: put_att, def_var, inq_var, put_int, put_double, hashfun
use ocean_blob_util_mod,  only: free_blob_memory
use ocean_parameters_mod, only: DEPTH_BASED
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, blob_grid_type
use ocean_types_mod,      only: ocean_blob_type, ocean_time_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type

implicit none

include 'netcdf.inc'
include 'mpif.h'

private

! Logicals for optional diagnostics
logical :: usei
logical :: usej
logical :: usek
logical :: usedepth
logical :: usest
logical :: usemass
logical :: usedensity
logical :: useent
logical :: usedet
logical :: userichardson
logical :: usegprime
logical :: usevolume
logical :: useu
logical :: usev
logical :: usew
logical :: useage
logical, allocatable, dimension(:) :: usetracer(:)
logical, allocatable, dimension(:) :: usefield(:)

character(len=31) :: filename
character(len=31) :: filename_modeltime
integer :: num_prog_tracers
integer :: index_temp
integer :: length
integer :: nblobs
integer :: vert_coordinate_class
real    :: dtime

integer :: isc,iec,jsc,jec,nk

type, private :: diag_blob_type
   integer :: type
   integer :: hash
   integer :: number
   real    :: time
   real    :: lat
   real    :: lon
   real    :: geodepth
   integer :: i
   integer :: j
   integer :: k
   real    :: depth
   real    :: st
   real    :: mass
   real    :: density
   real    :: ent
   real    :: det
   real    :: richardson
   real    :: gprime
   real    :: volume
   real    :: u
   real    :: v
   real    :: w
   real    :: age
   real, allocatable :: tracer(:)
   real, allocatable :: field(:)
   type(diag_blob_type), pointer :: next
end type diag_blob_type

type(diag_blob_type), pointer :: diag_head

public blob_diag_init
public blob_diag
public open_netcdf_file
public blob_diag_end

! namelist variables and default values
logical           :: blob_diagnostics = .false.  
character(len=31) :: diag_table       = "blob_diag_table"
integer           :: dump_num         = 2000000
integer           :: frequency        = 1                 !How many steps from one snapshot to the next

namelist /ocean_blob_diag_nml/ blob_diagnostics, diag_table, dump_num, frequency

contains

!######################################################################
! <SUBROUTINE NAME="ocean_blob_diag_init">
!
! <DESCRIPTION>
! Initialises the blob diagnostic module.
! </DESCRIPTION>
!
  subroutine blob_diag_init(T_prog, Grid, Domain, num_tracers, dtimein, &
                            itemp, ver_coordinate_class, do_diag)
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)
  type(ocean_grid_type),        intent(in)  :: Grid
  type(ocean_domain_type),      intent(in)  :: Domain
  integer,                      intent(in)  :: num_tracers
  real,                         intent(in)  :: dtimein
  integer,                      intent(in)  :: itemp
  integer,                      intent(in)  :: ver_coordinate_class
  logical,                      intent(out) :: do_diag
  

  integer :: ioun, ierr, io_status
  integer :: stdoutunit,stdlogunit 
  integer :: rootid, blobdim
  integer :: typeid, entryid, hashid, numberid, timeid
  integer :: latid, lonid, geodepthid
  integer :: iid, jid, kid
  integer :: depthid, stid, massid, volumeid, densityid, gprimeid
  integer :: entid, detid, richardsonid
  integer :: uid, vid, wid, ageid
  integer :: retval
  integer :: i,j,k,n,m
  integer :: maxhash
  logical :: validvar
  character(len=128) :: varname
  integer,allocatable,dimension(:) :: tracerid, fieldid

  stdoutunit = stdout(); stdlogunit=stdlog()

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_blob_diag_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_blob_diag_nml)  
  write (stdlogunit,ocean_blob_diag_nml)
  ierr = check_nml_error(io_status,'ocean_blob_diag_nml')
  call close_file (ioun)

  do_diag = blob_diagnostics

  if (.not. blob_diagnostics) return

  num_prog_tracers      = num_tracers
  vert_coordinate_class = ver_coordinate_class
  dtime                 = dtimein
  index_temp            = itemp
  diag_head => NULL()

  allocate(usetracer(num_prog_tracers))
  allocate(usefield(num_prog_tracers))

  allocate(tracerid(num_prog_tracers))
  allocate(fieldid(num_prog_tracers))

  usei        = .false.
  usej        = .false.
  usek        = .false.
  usedepth    = .false.
  usest       = .false.
  usemass     = .false.
  usedensity  = .false.
  useent      = .false.
  usedet      = .false.
  userichardson=.false.
  usegprime   = .false.
  usevolume   = .false.
  useu        = .false.
  usev        = .false.
  usew        = .false.
  useage      = .false.
  usetracer(:)  = .false.
  usefield(:)   = .false.
  
  write(stdoutunit,'(/,a)') 'Opening blob diagnostic table '//trim(diag_table)
  open(unit=1, iostat=ioun, file=trim(diag_table))
  if (ioun>0) then
     call error_mesg('ocean_blob_diag_mod, ocean_blob_diag_init, ',&
          'error opening blob diagnostic table: '//trim(diag_table), FATAL)
  endif

  write(stdoutunit,'(a)') 'Will write the following optional blob diagnostics:'
  m=0
  fileread: do
     read(unit=1, fmt=*, iostat=ioun) varname
     if (ioun>0) then
        ! there was an error: raise a warning but continue reading the file
        call error_mesg('ocean_blob_diag_mod, ocean_blob_diag_init, ',&
             'error reading blob diagnostic table: returned error '&
             //char(ioun)//' for variable '//trim(varname), WARNING)
        cycle fileread
     elseif (ioun<0) then
        ! end of the file: exit the loop
        exit fileread
     else
        ! line read ok: keep going
        write(stdoutunit,'(a)') trim(varname)
     endif
     m=m+1
     validvar=.false.

     if(trim(varname)=='entry')    validvar=.true.
     if(trim(varname)=='type')     validvar=.true.
     if(trim(varname)=='hash')     validvar=.true.
     if(trim(varname)=='number')   validvar=.true.
     if(trim(varname)=='time')     validvar=.true.
     if(trim(varname)=='lat')      validvar=.true.
     if(trim(varname)=='lon')      validvar=.true.
     if(trim(varname)=='geodepth') validvar=.true.
     if(trim(varname)=='i') then
        usei=.true.
        validvar=.true.
     endif
     if(trim(varname)=='j') then
        usej=.true.
        validvar=.true.
     endif
     if(trim(varname)=='k') then
        usek=.true.
        validvar=.true.
     endif
     if(trim(varname)=='depth') then
        usedepth=.true.
        validvar=.true.
     endif
     if(trim(varname)=='st') then
        usest=.true.
        validvar=.true.
     endif
     if(trim(varname)=='mass') then
        usemass=.true.
        validvar=.true.
     endif
     if(trim(varname)=='density') then
        usedensity=.true.
        validvar=.true.
     endif
     if(trim(varname)=='ent') then
        useent=.true.
        validvar=.true.
     endif
     if(trim(varname)=='det') then
        usedet=.true.
        validvar=.true.
     endif
     if(trim(varname)=='richardson') then
        userichardson = .true.
        validvar=.true.
     endif
     if(trim(varname)=='gprime') then
        usegprime=.true.
        validvar=.true.
     endif
     if(trim(varname)=='volume') then
        usevolume=.true.
        validvar=.true.
     endif
     if(trim(varname)=='u') then
        useu=.true.
        validvar=.true.
     endif
     if(trim(varname)=='v') then
        usev=.true.
        validvar=.true.
     endif
     if(trim(varname)=='w') then
        usew=.true.
        validvar=.true.
     endif
     if(trim(varname)=='blob_age') then
        useage=.true.
        validvar=.true.
     endif
     do n=1,num_prog_tracers
        if (n==index_temp) then
           if(trim(varname)=='heat') then
              usetracer(n)=.true.
              validvar=.true.
           endif
        else
           if(trim(varname)==trim(T_prog(n)%name)//'_cont') then
              usetracer(n)=.true.
              validvar=.true.
           endif
        endif
        if(trim(varname)==trim(T_prog(n)%name)) then
           usefield(n)=.true.
           validvar=.true.
        endif
     enddo
     if (.not.validvar) call error_mesg('ocean_blob_diag_mod, ocean_blob_diag_init, ', & 
          trim(varname)//' is not a recognised blob diagnostic field', WARNING)
  enddo fileread

  close(unit=1) 
  if(m==0) call error_mesg('ocean_blob_diag_mod, ocean_blob_diag_init, ',&
       'blob diagnostics table '//trim(diag_table)//' is empty', WARNING)

  write(filename, '("ocean_blobs.nc.",i4.4)') mpp_pe()

  call create_netcdf_file(filename, rootid, .false.)

  ! Write the maximum number of hash numbers
  isc=Domain%isc; iec=Domain%iec; jsc=Domain%jsc; jec=Domain%jec
  nk = Grid%nk
  maxhash = 0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if (Grid%tmask(i,j,k) > 0.) maxhash = maxhash + 1
        enddo
     enddo
  enddo
  retval = nf_put_att_int(rootid, NF_GLOBAL, 'nhash', NF_INT, 1, maxhash)
  if (retval/=NF_NOERR) call handle_error(retval, 'blob_diag_init, nhash, nf_put_att_int')

  ! Define the dimensions
  retval = nf_def_dim(rootid, 'entry', NF_UNLIMITED, blobdim)
  if (retval/=NF_NOERR) call handle_error(retval, 'blob_diag_init, entry nf_def_dim')

  ! Define the variables
  entryid    = def_var(rootid, 'entry',    NF_INT,    blobdim)
  typeid     = def_var(rootid, 'type',     NF_INT,    blobdim)
  hashid     = def_var(rootid, 'hash',     NF_INT,    blobdim)
  numberid   = def_var(rootid, 'number',   NF_INT,    blobdim)
  timeid     = def_var(rootid, 'time',     NF_DOUBLE, blobdim)
  latid      = def_var(rootid, 'lat',      NF_DOUBLE, blobdim)
  lonid      = def_var(rootid, 'lon',      NF_DOUBLE, blobdim)
  geodepthid = def_var(rootid, 'geodepth', NF_DOUBLE, blobdim)

  call put_att(rootid, entryid,    'name',     'entry number')
  call put_att(rootid, typeid,     'name',     'blob type')
  call put_att(rootid, hashid,     'name',     'blob hash')
  call put_att(rootid, numberid,   'name',     'blob number')
  call put_att(rootid, timeid,     'name',     'time')
  call put_att(rootid, timeid,     'units',    'days')
  call put_att(rootid, latid,      'name',     'blob latitude')
  call put_att(rootid, latid,      'units',    'degrees N')
  call put_att(rootid, lonid,      'name',     'blob longitude')
  call put_att(rootid, lonid,      'units',    'degrees E')
  call put_att(rootid, geodepthid, 'name',     'blob depth relative to z=0')
  call put_att(rootid, geodepthid, 'units',    'm')
  call put_att(rootid, geodepthid, 'positive', 'up')
  
  if (usei) then
     iid = def_var(rootid, 'i', NF_INT, blobdim)
     call put_att(rootid, iid, 'name', 'zonal grid cell number')
  endif

  if (usej) then
     jid = def_var(rootid, 'j', NF_INT, blobdim)
     call put_att(rootid, jid, 'name', 'meridional grid cell number')
  endif

  if (usek) then
     kid = def_var(rootid, 'k', NF_INT, blobdim)
     call put_att(rootid, kid, 'name', 'vertical grid cell number')
  endif

  if (usedepth) then
     depthid = def_var(rootid, 'depth', NF_DOUBLE, blobdim)
     call put_att(rootid, depthid, 'name',     'blob depth relative to free surface')
     call put_att(rootid, depthid, 'units',    'm')
     call put_att(rootid, depthid, 'positive', 'up')
  endif

  if (usest) then
     stid = def_var(rootid, 'st', NF_DOUBLE, blobdim)
     call put_att(rootid, stid, 'name', 'vertical position in native coordiante')
  endif

  if (usemass) then
     massid = def_var(rootid, 'mass', NF_DOUBLE, blobdim)
     call put_att(rootid, massid, 'name',  'blob mass')
     call put_att(rootid, massid, 'units', 'kg')
  endif

  if (useent) then
     entid = def_var(rootid, 'ent', NF_DOUBLE, blobdim)
     call put_att(rootid, entid, 'name', 'blob entrainment velocity')
     call put_att(rootid, entid, 'units', 'm/s')
  endif

  if (usedet) then
     detid = def_var(rootid, 'det', NF_DOUBLE, blobdim)
     call put_att(rootid, detid, 'name', 'blob detrainment velocity')
     call put_att(rootid, detid, 'units', 'm/s')
  endif

  if (userichardson) then
     richardsonid = def_var(rootid, 'richardson', NF_DOUBLE, blobdim)
     call put_att(rootid, richardsonid, 'name', 'blob Richardson number')
     call put_att(rootid, richardsonid, 'units', 'dimensionless')
  endif

  if (usedensity) then
     densityid = def_var(rootid, 'density', NF_DOUBLE, blobdim)
     call put_att(rootid, densityid, 'name',  'blob density')
     call put_att(rootid, densityid, 'units', 'kg/m^3')
  endif

  if (usegprime) then
     gprimeid = def_var(rootid, 'gprime', NF_DOUBLE, blobdim)
     call put_att(rootid, gprimeid, 'name',  'blob reduced gravity')
     call put_att(rootid, gprimeid, 'units', 'm/s^2')
  endif

  if (usevolume) then
     volumeid = def_var(rootid, 'volume', NF_DOUBLE, blobdim)
     call put_att(rootid, volumeid, 'name',  'blob volume')
     call put_att(rootid, volumeid, 'units', 'm^3')
  endif

  if (useu) then
     uid = def_var(rootid, 'u', NF_DOUBLE, blobdim)
     call put_att(rootid, uid, 'name',  'blob zonal velocity')
     call put_att(rootid, uid, 'units', 'm/s')
  endif

  if (usev) then
     vid = def_var(rootid, 'v', NF_DOUBLE, blobdim)
     call put_att(rootid, vid, 'name',  'blob meridional velocity')
     call put_att(rootid, vid, 'units', 'm/s')
  endif

  if (usew) then
     wid = def_var(rootid, 'w', NF_DOUBLE, blobdim)
     call put_att(rootid, wid, 'name',  'blob vertical velocity')
     call put_att(rootid, wid, 'units', 'm/s')
  endif

  if(useage) then
     ageid = def_var(rootid, 'blob_age', NF_DOUBLE, blobdim)
     call put_att(rootid, ageid, 'name', 'blob age')
     call put_att(rootid, ageid, 'units', 'yr')
  endif

  do n=1,num_prog_tracers
     if (usetracer(n)) then
        if(n==index_temp) then
           tracerid(n) = def_var(rootid, 'heat', NF_DOUBLE, blobdim)
           call put_att(rootid, tracerid(n), 'name',  'blob heat content')
           call put_att(rootid, tracerid(n), 'units', 'J')
        else
           tracerid(n) = def_var(rootid, trim(T_prog(n)%name)//'_cont', NF_DOUBLE, blobdim)
           call put_att(rootid, tracerid(n), 'name',  'blob '//trim(T_prog(n)%name)//' content')
           call put_att(rootid, tracerid(n), 'units', 'kg')
        endif
     endif
     
     if (usefield(n)) then
        fieldid(n) = def_var(rootid, trim(T_prog(n)%name), NF_DOUBLE, blobdim)
        if(n==index_temp) then
           call put_att(rootid, fieldid(n), 'name', trim(T_prog(n)%name))
           call put_att(rootid, fieldid(n), 'units', 'C')
        else
           call put_att(rootid, fieldid(n), 'name',  'blob '//trim(T_prog(n)%name)//' concentration')
           call put_att(rootid, fieldid(n), 'units', 'kg/kg')
        endif
     endif
     
  enddo

  call close_netcdf_file(rootid)

  nblobs = 0
  length = 0

end subroutine blob_diag_init
! </SUBROUTINE> NAME="ocean_blob_diag_init"

!######################################################################
! <SUBROUTINE NAME="blob_diag">
!
! <DESCRIPTION>
! Accummulates the blob diagnostics by creating a linked list of
! diagnostic blobs.  The blobs are kept in the linked list until
! there are more than dump_num of them.  Then, they are written
! (using write_blobs), and erased from memory.
! </DESCRIPTION>
!
subroutine blob_diag(Time, head, T_prog, blob_type)
  type(ocean_time_type),        intent(in) :: Time
  type(ocean_blob_type),        pointer    :: head
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: blob_type

  type(ocean_blob_type), pointer :: this=>NULL()
  type(diag_blob_type), pointer  :: blob=>NULL()
  integer :: n, lenflag
  real :: day
  
  if (.not. blob_diagnostics) return

  if (mod(Time%itt,frequency) /= 0) return

  day = real(Time%itt0)*dtime/86400.
  lenflag = 0

  if(associated(head)) then
     this=>head
     blobcycle : do
        
        allocate(blob)
        allocate(blob%tracer(num_prog_tracers))
        allocate(blob%field(num_prog_tracers))

        blob%type     = blob_type
        blob%hash     = this%hash
        blob%number   = this%number
        blob%time     = day
        blob%lat      = this%lat
        blob%lon      = this%lon
        blob%geodepth = this%geodepth
        blob%i        = this%i
        blob%j        = this%j
        blob%k        = this%k
        blob%depth    = this%depth
        blob%st       = this%st
        blob%mass     = this%mass
        blob%density  = this%density
        blob%ent      = this%ent
        blob%det      = this%det
        blob%richardson=this%richardson
        blob%gprime   = this%gprime
        blob%volume   = this%volume
        blob%u        = this%v(1)
        blob%v        = this%v(2)
        blob%w        = this%v(3)
        blob%age      = this%age
        do n=1,num_prog_tracers
           blob%tracer(n)  = this%tracer(n)
           blob%field(n)   = this%field(n)
        enddo

        if(associated(diag_head)) then
           blob%next => diag_head
        else
           blob%next => NULL()
        endif
        diag_head => blob
        nullify(blob)

        length = length + 1

        this=>this%next
        if (.not.associated(this)) exit blobcycle
     enddo blobcycle
  endif
  
  ! So that we don't use too much memory, if any PE has more blobs
  ! than dump_num, we write all the diagnostic blobs (on all pes)
  ! and deallocate them from memory, and then start afresh with an 
  ! empty list.
  if(length>dump_num) lenflag = 1

  call mpp_sum(lenflag)
  if(lenflag>0) then
     call write_blobs(T_prog(:))
     length = 0
  endif

end subroutine blob_diag
! </SUBROUTINE> NAME="ocean_blob_diag_init"

!######################################################################
! <FUNCTION NAME="varid">
!
! <DESCRIPTION>
! Reads the variable id of a netcdf file.
! </DESCRIPTION>
!
integer function varid(rootid, name)
  integer,          intent(in)  :: rootid
  character(len=*), intent(in)  :: name

  integer :: retval

  retval = nf_inq_varid(rootid, trim(name), varid)
  if(retval/=NF_NOERR) then 
     call handle_error(retval, 'varid, '//trim(name)//' nf_inq_varid')
  endif

end function varid
! </FUNCTION> NAME="varid"


!######################################################################
! <SUBROUTINE NAME="blob_diag_end">
!
! <DESCRIPTION>
! </DESCRIPTION>
!
subroutine blob_diag_end(T_prog)
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  
  if (.not. blob_diagnostics) return

  call write_blobs(T_prog(:))

end subroutine blob_diag_end
! </SUBROUTINE> NAME="ocean_blob_diag_end"

!######################################################################
! <SUBROUTINE NAME="write_blobs">
!
! <DESCRIPTION>
! Write the diagnostics of individual blobs.
! </DESCRIPTION>
!
subroutine write_blobs(T_prog)
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  type(diag_blob_type), pointer :: prev=>NULL()
  type(diag_blob_type), pointer :: blob=>NULL()
  integer :: n
  integer :: rootid, retval
  integer :: entryid, typeid, hashid, numberid, timeid
  integer :: latid, lonid, geodepthid
  integer :: iid, jid, kid
  integer :: depthid, stid, massid, volumeid, densityid, gprimeid
  integer :: entid, detid, richardsonid
  integer :: uid, vid, wid, ageid
  integer :: tracerid(num_prog_tracers)
  integer :: fieldid(num_prog_tracers)

  call open_netcdf_file(filename, rootid)
  
  retval = nf_inq_dimid(rootid, 'entry', entryid)
  if(retval/=NF_NOERR) call handle_error(retval, 'blob_diag, entry nf_inq_dimid')
  
  typeid     = varid(rootid, 'type')
  hashid     = varid(rootid, 'hash')
  numberid   = varid(rootid, 'number')
  timeid     = varid(rootid, 'time')
  latid      = varid(rootid, 'lat')
  lonid      = varid(rootid, 'lon')
  geodepthid = varid(rootid, 'geodepth')
  if (usei)       iid       = varid(rootid, 'i')
  if (usej)       jid       = varid(rootid, 'j')
  if (usek)       kid       = varid(rootid, 'k')
  if (usedepth)   depthid   = varid(rootid, 'depth')
  if (usest)      stid      = varid(rootid, 'st')
  if (usemass)    massid    = varid(rootid, 'mass')
  if (usedensity) densityid = varid(rootid, 'density')
  if (useent)     entid     = varid(rootid, 'ent')
  if (usedet)     detid     = varid(rootid, 'det')
  if (userichardson) richardsonid = varid(rootid, 'richardson')
  if (usegprime)  gprimeid  = varid(rootid, 'gprime')
  if (usevolume)  volumeid  = varid(rootid, 'volume')
  if (useu)       uid       = varid(rootid, 'u')
  if (usev)       vid       = varid(rootid, 'v')
  if (usew)       wid       = varid(rootid, 'w')
  if (useage)     ageid     = varid(rootid, 'blob_age')
  do n=1,num_prog_tracers
     if (n==index_temp) then
        if (usetracer(n))  tracerid(n)  = varid(rootid, 'heat')
     else
        if (usetracer(n))  tracerid(n)  = varid(rootid, trim(T_prog(n)%name)//'_cont')
     endif
     if (usefield(n)) fieldid(n) = varid(rootid, trim(T_prog(n)%name))
  enddo
  
  if (associated(diag_head)) then
     blob=>diag_head
     blobcycle: do
        nblobs=nblobs+1
        call put_int(   rootid, entryid,    nblobs, nblobs)
        call put_int(   rootid, typeid,     nblobs, blob%type)
        call put_int(   rootid, hashid,     nblobs, blob%hash)
        call put_int(   rootid, numberid,   nblobs, blob%number)
        call put_double(rootid, timeid,     nblobs, blob%time)
        call put_double(rootid, latid,      nblobs, blob%lat)
        call put_double(rootid, lonid,      nblobs, blob%lon)
        call put_double(rootid, geodepthid, nblobs, blob%geodepth)
        if (usest) then
           if (vert_coordinate_class==DEPTH_BASED) then
              call put_double(rootid, stid, nblobs, blob%st)
           else
              call put_double(rootid, stid, nblobs, blob%st*c2dbars)
           endif
        endif
        if (usei)       call put_int(   rootid, iid,        nblobs, blob%i)
        if (usej)       call put_int(   rootid, jid,        nblobs, blob%j)
        if (usek)       call put_int(   rootid, kid,        nblobs, blob%k)
        if (usedepth)   call put_double(rootid, depthid,    nblobs, blob%depth)
        if (usemass)    call put_double(rootid, massid,     nblobs, blob%mass)
        if (usedensity) call put_double(rootid, densityid,  nblobs, blob%density)
        if (useent)     call put_double(rootid, entid,      nblobs, blob%ent)
        if (usedet)     call put_double(rootid, detid,      nblobs, blob%det)
        if (userichardson) call put_double(rootid, richardsonid, nblobs, blob%richardson)
        if (usegprime)  call put_double(rootid, gprimeid,   nblobs, blob%gprime)
        if (usevolume)  call put_double(rootid, volumeid,   nblobs, blob%volume)
        if (useu)       call put_double(rootid, uid,        nblobs, blob%u)
        if (usev)       call put_double(rootid, vid,        nblobs, blob%v)
        if (usew)       call put_double(rootid, wid,        nblobs, blob%w)
        if (useage)     call put_double(rootid, ageid,      nblobs, blob%age)
        do n=1,num_prog_tracers
           if (usetracer(n))  call put_double(rootid, tracerid(n),  nblobs, blob%tracer(n))
           if (usefield(n))   call put_double(rootid, fieldid(n),   nblobs, blob%field(n))
        enddo

        prev=>blob
        blob=>blob%next
        if(allocated(prev%tracer))  deallocate(prev%tracer)
        if(allocated(prev%field))   deallocate(prev%field)
        deallocate(prev)
        nullify(prev)
        if(.not.associated(blob)) exit blobcycle
     enddo blobcycle
     nullify(diag_head)
  endif

  call close_netcdf_file(rootid)

end subroutine write_blobs
! </SUBROUTINE> NAME="write_blobs"


!######################################################################
! <SUBROUTINE NAME="handle_error">
!
! <DESCRIPTION>
! Handles any errors from the reading/writing of netcdf files.  It
! should (hopefully) provide some sort of useful idea of what went 
! wrong.
! </DESCRIPTION>
!
subroutine handle_error(retval, operation)
  integer,          intent(in) :: retval
  character(len=*), intent(in) :: operation
  
  print '(a,i3)', 'ocean_blob_diag_mod, '&
       //trim(operation)//': failed with error number ', retval
  if(retval==nf_emaxname)   print *, '==> NF_MAX_NAME exceeded'
  if(retval==nf_enameinuse) print *, '==> Name already in use'
  if(retval==nf_ebadname)   print *, '==> Attribute or variable name contains illegal characters'
  if(retval==nf_ebadid)     print *, '==> Invalid ncid'
  if(retval==nf_einval)     print *, '==> Size is invalid'
  if(retval==nf_enomem)     print *, '==> Out of memory'
  if(retval==nf_ebadtype)   print *, '==> Cannot find the type id'
  if(retval==nf_ebaddim)    print *, '==> Bad dimension id'
  if(retval==nf_enotvar)    print *, '==> Variable not found'
  
  call error_mesg('ocean_blob_diag_mod, '//trim(operation), &
       ' netcdf function returned a failure!', FATAL)

end subroutine handle_error
! </SUBROUTINE> NAME="handle_error"

!######################################################################
! <SUBROUTINE NAME="create_netcdf_file">
!
! <DESCRIPTION>
! Creates a new netcdf file.
! </DESCRIPTION>
!
subroutine create_netcdf_file(filename, fileid, parallel)
  character(len=31), intent(in)  :: filename
  integer,           intent(out) :: fileid
  logical,           intent(in)  :: parallel
  integer :: mret
  integer :: stderrunit

  stderrunit = stderr()

  if (parallel) then
     call error_mesg('ocean_blob_diag_mod: ', 'Parallel access not supported!', FATAL)
  else
     mret = nf_create(filename, NF_CLOBBER, fileid)
  endif
  if (mret .ne. NF_NOERR) then
     write(stderrunit, '(/,a)') 'ocean_blob_diag_mod, create_netcdf_file: NF_CREATE failed'
     call error_mesg('ocean_blob_diag_mod,  create_netcdf_file', &
          'netcdf function NF_CREATE returned a failure!', FATAL)
  endif
end subroutine create_netcdf_file
! </SUBROUTINE> NAME="create_netcdf_file"


!######################################################################
! <SUBROUTINE NAME="open_netcdf_file">
!
! <DESCRIPTION>
! Opens an existing netcdf file.
! </DESCRIPTION>
!
subroutine open_netcdf_file(filename, fileid)
  character(len=31), intent(in)    :: filename
  integer,           intent(inout) :: fileid
  integer :: mret
  integer :: stderrunit

  stderrunit = stderr()

  mret = nf_open(filename, NF_WRITE, fileid)
  if (mret .ne. NF_NOERR) then
     write(stderrunit, '(/,a)') 'ocean_blob_diag_mod, open_netcdf_file: NF_OPEN failed'
     call error_mesg('ocean_blob_diag_mod, open_netcdf_file', &
          'netcdf function NF_OPEN returned a failure!', FATAL)
  endif
end subroutine open_netcdf_file
! </SUBROUTINE> NAME="open_netcdf_file"


!######################################################################
! <SUBROUTINE NAME="close_netcdf_file">
!
! <DESCRIPTION>
! Closes an existing netcdf file.
! </DESCRIPTION>
!
subroutine close_netcdf_file(fileid)
  integer, intent(inout) :: fileid
  integer :: mret
  integer :: stderrunit

  stderrunit = stderr()

  mret = nf_close(fileid)
  if (mret .ne. NF_NOERR) then
     write(stderrunit, '(/,a)') 'ocean_blob_diag_mod, close_netcdf_file: NF_CLOSE failed'
     call error_mesg('ocean_blob_diag_mod,  close_netcdf_file', &
          'netcdf function NF_CLOSE returned a failure!', FATAL)
  endif
end subroutine close_netcdf_file
! </SUBROUTINE> NAME="close_netcdf_file"


end module ocean_blob_diag_mod

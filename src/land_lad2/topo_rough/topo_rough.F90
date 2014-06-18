module topo_rough_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Sergey Malyshev
! </CONTACT>

  use time_manager_mod,   only : time_type
  use mpp_domains_mod,    only : domain2d

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

  use fms_mod,            only : write_version_number, error_mesg, FATAL, NOTE, &
       open_restart_file, set_domain, read_data, &
       write_data, close_file, file_exist, check_nml_error, mpp_pe, &
       mpp_root_pe, stdlog
  use diag_manager_mod,   only : register_static_field, send_data
  use topography_mod,     only : get_topog_stdev

implicit none
private
! ==== public interface ======================================================
public :: topo_rough_init
public :: topo_rough_end
public :: update_topo_rough
! ==== end of public interface ===============================================


! <NAMELIST NAME = "topo_rough_nml">
!   <DATA NAME="use_topo_rough" TYPE="logical" DEFAULT="false">
!     If true, the topographic momentum drag scaling scheme is used
!   </DATA>
!   <DATA NAME="max_topo_rough" TYPE="real" DEFAULT="100" UNITS="m">
!     Maximum of topographic "roughness length" used for momentum drag scaling
!   </DATA>
!   <DATA NAME="topo_rough_factor" TYPE="real" DEFAULT="1.0">
!     Scaling factor to convert topography variance to topographic 
!     "roughness length"
!   </DATA>
!   <DATA NAME="topo_rough_source" TYPE="caharacter(len=16)" DEFAULT="'computed'">
!     Source of the sub-grid topography variance data for topographic momentum drag scaling. 
!     'computed' means that the variance is calculated based on high-resolution 
!     topography data. 'input' means that the data will be provided in specified file
!     (NetCDF of IEEE binary)
!   </DATA>
!   <DATA NAME="topo_rough_file" TYPE="character(len=256)" DEFAULT="INPUT/mg_drag.data.nc">
!     Name of the file to be used as an input for sub-grid topography variance data. 
!     The file can be either NetCDF (in this case variable name can also be specified), or
!     IEEE.
!   </DATA>
!   <DATA NAME="topo_rough_var" TYPE="character(len=128)" DEFAULT="ghprime">
!     Name of the NetCDF variable to be used as a topography variance field. Ignored if
!     the file specified in topo_rough_file is not NetCDF file.
!   </DATA>
! </NAMELIST>

logical     :: use_topo_rough    = .false.
real        :: max_topo_rough    = 100 ! m
real        :: topo_rough_factor = 1.0
character(len=16) :: topo_rough_source = 'computed'
character(len=256):: topo_rough_file   = 'INPUT/mg_drag.data.nc'
character(len=128):: topo_rough_var    = 'ghprime'

namelist/topo_rough_nml/ use_topo_rough, topo_rough_factor, max_topo_rough, &
     topo_rough_source, topo_rough_file, topo_rough_var

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name   = 'she_topo_rough', &
     diag_mod_name = 'topo_rough', &
     version       = '$Id: topo_rough.F90,v 19.0 2012/01/06 20:43:20 fms Exp $', &
     tagname       = '$Name: tikal $'

! ==== module private data ===================================================
real, allocatable, save ::topo_stdev(:,:)
logical :: module_is_initialized = .FALSE.

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'

contains ! ###################################################################

subroutine topo_rough_init(time, lonb, latb, domain, id_lon,id_lat)
  type(time_type), intent(in) :: time            ! current time
  type(domain2d) , intent(in) :: domain          ! our domain
  real           , intent(in) :: latb(:,:),lonb(:,:) ! boundaries of the grid cells
  integer        , intent(in) :: id_lon,id_lat   ! IDs of diagnostic axes
!   <ERROR MSG="could not read topography data" STATUS="FATAL">
!     get_topog_stdev failed to provide topography variance data.
!   </ERROR>  
!   <ERROR MSG="input file for for topography standard deviation ... does not exist" STATUS="FATAL">
!     topo_rough_source is set to 'input', but input file name either
!     not specified or specified incorrectly, so the program cannot 
!     find it.
!   </ERROR>
!   <ERROR MSG="... is not a valid value for topo_rough_source" STATUS="FATAL">
!     specified value of namelist parameter topo_rough_source is invalid; 
!     valid values are 'computed' or 'input'.
!   </ERROR>
  ! --- local vars
  integer :: ierr,io,unit
  integer :: id
  logical :: used, got_stdev

  ! write the version and tagname to the logfile
  call write_version_number(version, tagname)

  ! read and write (to logfile) namelist variables
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=topo_rough_nml, iostat=io)
  ierr = check_nml_error(io, 'topo_rough_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=topo_rough_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'topo_rough_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=topo_rough_nml)
  endif

  ! allocate topo_stdev according to specified domain
  allocate(topo_stdev(size(lonb,1)-1, size(lonb,2)-1))

  if (use_topo_rough) then

     if(trim(topo_rough_source) == 'computed') then
        call error_mesg('topo_rough_init','computing topography standard deviation',NOTE)
        got_stdev = get_topog_stdev(lonb,latb,topo_stdev)
        if (.not.got_stdev) &
             call error_mesg ('topo_rough_init', &
             'could not read topography data', FATAL)
     else if (trim(topo_rough_source)=='input') then
        call error_mesg('topo_rough_init','reading topography standard deviation from "'&
             //trim(topo_rough_file)//'"',NOTE)
        if(.not.file_exist(topo_rough_file,domain))&
             call error_mesg('topo_rough_init',            &
             'input file for topography standard deviation "'// &
             trim(topo_rough_file)//'" does not exist', FATAL)
        
        call set_domain(domain)
        call read_data(topo_rough_file,topo_rough_var,topo_stdev)
     else
        call error_mesg('topo_rough_init','"'//trim(topo_rough_source)//&
             '" is not a valid value for topo_rough_source', FATAL)
     endif
     topo_stdev = min(topo_stdev*topo_rough_factor,max_topo_rough)
  else
     topo_stdev = 0.0
  endif

  ! diag output : send topo_stdev to diagnostics
  id = register_static_field(diag_mod_name,'topo_rough',(/id_lon,id_lat/), &
       'momentum drag coefficient scaling lenght','m',missing_value=-1.0 )
  if(id > 0) &
       used = send_data(id,topo_stdev,time)
  module_is_initialized = .TRUE.
end subroutine topo_rough_init

! ============================================================================
subroutine topo_rough_end()
  deallocate(topo_stdev)
  module_is_initialized = .FALSE.
end subroutine

! ============================================================================
subroutine update_topo_rough(topo_rough)
  real, intent(out) :: topo_rough(:,:,:)

  ! ---- local vars
  integer :: k

  ! just assign standard deviation (scaled and trimmed according to namelist 
  ! parameters) to the output field 
  do k = 1, size(topo_rough,3)
     topo_rough(:,:,k) = topo_stdev(:,:)
  enddo
end subroutine

end module topo_rough_mod

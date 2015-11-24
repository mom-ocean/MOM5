module ocean_grids_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang 
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> M.J. Harrison 
!</REVIEWER>
!
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
! S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! Set up the ocean model grid spacing 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module sets up the ocean model grid based on information read in 
! from the grid_spec.nc file. It translates the generic names from the 
! grid_spec.nc file to the names used by MOM. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, A. Rosati, and R.C. Pacanowski 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_grids_nml">
!
!  <DATA NAME="read_rho0_profile" TYPE="logical">
!  To read in an initial rho0(k) profile to assist in defining the 
!  initial settings for the pressure increments dst, for use in 
!  setting the pressure-based vertical coordinate grids.  Ideally,
!  this profile is determined by the level averaged density in 
!  the initial conditions. Note that it is essential to have 
!  rho0_profile have a sensible value at all depths even if there
!  is no water there, since there are places where we divide by 
!  rho0_profile in rock.  Also, be mindful that with denser water 
!  at depth, the pressure levels will be coarser at depth than if 
!  using the trivial density profile rho0(k)=rho0. 
!  This option is experimental, so it is recommended that user
!  maintain the default read_rho0_profile=.false.  
!  </DATA> 
!
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is do_bitwise_exact_sum=.false.
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging. Note that most of the debugging stuff 
!  has been removed, but keep flag around in case need in future.
!  </DATA> 
!
!  <DATA NAME="verbose_init" TYPE="logical">
!  Prints out lots of initial checksums.  Useful to have on, so 
!  defaulted to true. 
!   </DATA> 
!
! </NAMELIST>
!
use constants_mod,     only: pi, radian, radius, epsln, c2dbars, deg_to_rad
use diag_manager_mod,  only: diag_axis_init, set_diag_global_att
use diag_manager_mod,  only: register_static_field, send_data
use fms_mod,           only: write_version_number
use fms_mod,           only: read_data, write_data
use fms_mod,           only: open_namelist_file, close_file, check_nml_error
use fms_mod,           only: file_exist, field_exist, field_size, get_global_att_value
use mpp_domains_mod,   only: mpp_update_domains, domain2d
use mpp_domains_mod,   only: mpp_global_sum, CORNER, EAST, NORTH, XUPDATE
use mpp_domains_mod,   only: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM
use mpp_domains_mod,   only: mpp_define_domains, mpp_get_domain_extents, mpp_deallocate_domain
use mpp_domains_mod,   only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_mod,           only: input_nml_file, stdout, stdlog, mpp_error, mpp_sum, FATAL, NOTE, mpp_pe
use mosaic_mod,        only: get_mosaic_ntiles, get_mosaic_ncontacts, get_mosaic_contact

use ocean_domains_mod,    only: get_local_indices, get_global_indices, get_halo_sizes, get_domain_offsets
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, PSIGMA
use ocean_parameters_mod, only: PRESSURE_BASED, DEPTH_BASED
use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
use ocean_parameters_mod, only: missing_value, rho0, grav
use ocean_types_mod,      only: ocean_grid_type, ocean_time_type, ocean_domain_type
use ocean_workspace_mod,  only: wrk1_2d 
use ocean_util_mod,       only: diagnose_2d_u, diagnose_3d_u, diagnose_2d, diagnose_3d, write_chksum_2d

implicit none

private

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS
integer, parameter :: xhalo=1
integer, parameter :: yhalo=1
#else
integer, private   :: xhalo
integer, private   :: yhalo
#endif

! for diagnostics 
integer :: id_ht          =-1
integer :: id_hu          =-1
integer :: id_dht_dx      =-1
integer :: id_dht_dy      =-1
integer :: id_gradH       =-1
integer :: id_kmt         =-1
integer :: id_kmu         =-1
integer :: id_tmask       =-1
integer :: id_umask       =-1
integer :: id_tmasken(2)  =-1

integer :: id_dxt=-1
integer :: id_dxu=-1
integer :: id_dyt=-1
integer :: id_dyu=-1
integer :: id_dxtn=-1
integer :: id_dyte=-1

logical :: used

character(len=128) :: version=&
     '$Id: ocean_grids.F90,v 20.0 2013/12/14 00:10:44 fms Exp $'
character (len=128) :: tagname = &
     '$Name: tikal $'

public ocean_grids_init
public set_ocean_grid_size
public set_ocean_hgrid_arrays
public set_ocean_vgrid_arrays
public init_grids_diag
public update_boundaries

private axes_info


! for setting dst with pressure-based vertical coordinate models 
real, allocatable, dimension(:) :: rho0_profile 

integer :: vert_coordinate 
integer :: vert_coordinate_class
integer :: horz_grid 

logical :: module_is_initialized = .FALSE.

! for bitwise exact global sums independent of PE number
logical :: do_bitwise_exact_sum = .false.
integer :: global_sum_flag

!mosaic         
integer             :: grid_version 
integer, parameter  :: VERSION_0 = 0  ! grid file with field geolon_t
integer, parameter  :: VERSION_1 = 1  ! grid file with field x_T
integer, parameter  :: VERSION_2 = 2  ! mosaic file

! nml variables 
logical :: debug_this_module = .false.
logical :: verbose_init      = .true.
logical :: read_rho0_profile = .false. 
logical :: write_grid        = .false.

namelist /ocean_grids_nml/ debug_this_module, verbose_init, read_rho0_profile, &
                           do_bitwise_exact_sum, write_grid

! grid file name
character(len=256) :: ocean_hgrid        ! will be set in set_ocean_grid_size
character(len=256) :: ocean_vgrid = 'INPUT/ocean_vgrid.nc'

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_grids_init">
!
! <DESCRIPTION>
! Initialize the grids module. 
! </DESCRIPTION>
!
subroutine ocean_grids_init(ver_coordinate, ver_coordinate_class, hor_grid, debug)

  integer, intent(in)           :: ver_coordinate
  integer, intent(in)           :: ver_coordinate_class
  integer, intent(in)           :: hor_grid 
  logical, intent(in), optional :: debug

  integer :: ioun, io_status, ierr
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error from ocean_grids_mod (ocean_grids_init): module has been initialized')
  endif 

  module_is_initialized = .TRUE.

  vert_coordinate       = ver_coordinate 
  vert_coordinate_class = ver_coordinate_class
  horz_grid             = hor_grid  

  call write_version_number(version, tagname)

  if (PRESENT(debug)) debug_this_module = debug

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_grids_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_grids_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_grids_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_grids_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_grids_nml)  
  write (stdlogunit,ocean_grids_nml)

  if(do_bitwise_exact_sum) then
     global_sum_flag = BITWISE_EXACT_SUM
  else
     global_sum_flag = NON_BITWISE_EXACT_SUM
  endif


end subroutine ocean_grids_init
! </SUBROUTINE> NAME="ocean_grids_init"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_grid_size">
!
! <DESCRIPTION>
! Set the ocean grid size.  Model expects the grid specification file
! to be called grid_spec.nc.  
! </DESCRIPTION>
!
subroutine set_ocean_grid_size(Grid, grid_file, grid_name)

  type(ocean_grid_type), intent(inout)        :: Grid
  character(len=*),      intent(in), optional :: grid_file
  character(len=*),      intent(in), optional :: grid_name

  integer                                   :: siz(4)
  integer                                   :: m
  integer                                   :: ntiles, ncontacts
  integer, dimension(2)                     :: tile1, tile2
  integer, dimension(2)                     :: istart1, iend1, jstart1, jend1
  integer, dimension(2)                     :: istart2, iend2, jstart2, jend2
  character(len=256)                        :: grd_file, ocean_mosaic, attvalue
  real                                      :: f_plane_latitude
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  Grid%ni=0 ; Grid%nj=0 ; Grid%nk=0 
  Grid%mosaic = .false.
  grd_file = "INPUT/grid_spec.nc"
  if(present(grid_file)) grd_file = grid_file
  if (.not. file_exist(grd_file) ) then 
       call mpp_error(FATAL, '==> Error from ocean_grids_mod(set_ocean_grid_size): '// &
                             'grid specification file '//grd_file//' does not exist')
  endif 
  !
  !  Determine if the grid is mosaic file
  !
  if(field_exist(grd_file, 'ocn_mosaic_file') .or. field_exist(grd_file, 'gridfiles') ) then ! read from mosaic file
     write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from mosaic version grid'
     grid_version = VERSION_2
     Grid%mosaic = .true.
     if( field_exist(grd_file, 'ocn_mosaic_file') ) then ! coupler mosaic
        call read_data(grd_file, "ocn_mosaic_file", ocean_mosaic)
        ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
     else
        ocean_mosaic = trim(grd_file)
     end if
     ntiles = get_mosaic_ntiles(ocean_mosaic)
     if(ntiles .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                   'ntiles should be 1 for ocean mosaic, contact developer')
! Notify diag_manager of mosaic grid
     call set_diag_global_att ('ocean','mosaic','1')
     call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
     ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)
     call field_size(ocean_hgrid, 'x', siz)
     
     if(mod(siz(1),2) .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
          'x-size of x in file '//trim(ocean_hgrid)//' should be 2*ni+1')
     if(mod(siz(2),2) .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
          'y-size of x in file '//trim(ocean_hgrid)//' should be 2*nj+1')
     Grid%ni = siz(1)/2 
     Grid%nj = siz(2)/2 
     call field_size(ocean_vgrid , "zeta", siz) 
     if(mod(siz(1),2) .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
          'size of dimension zeta in file '//trim(ocean_vgrid)//' should be 2*nk+1')
     Grid%nk  = siz(1)/2            
  else  if(field_exist(grd_file, 'x_T')) then
     ocean_hgrid = grd_file
     write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from new version grid'
     grid_version = VERSION_1
     call field_size( ocean_hgrid, 'x_T', siz)
     Grid%ni = siz(1)
     Grid%nj = siz(2) 
     call field_size( ocean_hgrid, "zt", siz) 
     Grid%nk  = siz(1)     
  else if(field_exist(grd_file, 'geolon_t')) then
     ocean_hgrid = grd_file
     write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from old version grid'
     grid_version = VERSION_0 
     call field_size( ocean_hgrid, 'geolon_t', siz)  
     Grid%ni = siz(1)
     Grid%nj = siz(2)
     call field_size( ocean_hgrid, "zt", siz) 
     Grid%nk  = siz(1)       
     else
     call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                     'x_T, geolon_t, ocn_mosaic_file, gridfiles does not exist in file ' //trim(grd_file))
     endif

  if (Grid%ni == 0 .or. Grid%nj == 0 .or. Grid%nk == 0) then
     write(stdoutunit,*) '==>Error reading grid information from ',trim(grd_file),'. Make sure file exists'
     call mpp_error(FATAL,'==>Error reading grid information from grid file.  Are you sure file exists?')
  endif

  Grid%cyclic_x=.false.;Grid%cyclic_y=.false.;Grid%tripolar=.false.
  Grid%f_plane=.false.; Grid%beta_plane=.false.
  Grid%f_plane_latitude = -999.0

  if(grid_version == VERSION_2) then
     !z1l: f_plane, beta_plane area not supported in mosaic grid. Need to think about to implement this. 
     if(field_exist(ocean_mosaic, "contacts") ) then
        ncontacts = get_mosaic_ncontacts(ocean_mosaic)
        if(ncontacts < 1) call mpp_error(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                       'number of contacts should be larger than 0 when field contacts exist in file '//trim(ocean_mosaic) )
        if(ncontacts > 2) call mpp_error(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                       'number of contacts should be no larger than 2')
        call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
             istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
             istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
        do m = 1, ncontacts
           if(istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
              if(istart2(m) .NE. iend2(m) ) call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): only cyclic condition is allowed for x-boundary")
             Grid%cyclic_x = .true.
           else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
              if(jstart2(m) .NE. jend2(m) ) call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): "//&
                   "only cyclic/folded-north condition is allowed for y-boundary")
              if( jstart1(m) == jstart2(m) ) then ! folded north
             Grid%tripolar = .true.
         else
                 Grid%cyclic_y = .true.
         endif 
           else 
              call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): invalid boundary contact")
         end if
  end do
     end if
  else
     if( get_global_att_value(ocean_hgrid, "x_boundary_type", attvalue) ) then
        if(attvalue == 'cyclic') Grid%cyclic_x = .true.
     end if
     if( get_global_att_value(ocean_hgrid, "y_boundary_type", attvalue) ) then
        if(attvalue == 'cyclic') then
           Grid%cyclic_y = .true.
        else if(attvalue == 'fold_north_edge') then
           Grid%tripolar = .true.
        end if
     end if
     if(get_global_att_value(ocean_hgrid, "f_plane", attvalue) ) then
        if(attvalue == 'y') Grid%f_plane = .true.
     end if
     if(get_global_att_value(ocean_hgrid, "beta_plane", attvalue) ) then
        if(attvalue == 'y') Grid%beta_plane = .true.
     end if
     if( get_global_att_value(ocean_hgrid, "f_plane_latitude", f_plane_latitude) ) then
        Grid%f_plane_latitude = f_plane_latitude
     end if
  end if

  if(Grid%cyclic_x) then
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is cyclic')
  else
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is solid_walls')
  end if

  if(Grid%tripolar) then
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is fold_north_edge')
  else if(Grid%cyclic_y) then
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is cyclic')
  else
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is solid_walls')
  endif

  if(Grid%f_plane .or. Grid%beta_plane) then
     if(abs(Grid%f_plane_latitude) > 90.0 ) then
        write(stdoutunit,*)  "==>Error from ocean_grids_mod(set_ocean_grid_size): "// &
                 "f_plane_latitude must be between -90 and 90 degrees. Please check your grid file."
        call mpp_error(FATAL, "==>Error from ocean_grids_mod(set_ocean_grid_size): "// &
                 "f_plane_latitude must be between -90 and 90 degrees. Please check your grid file.")
     end if
     write(stdoutunit,'(a,f8.2)') '==>f-plane/beta-plane latitude = ',Grid%f_plane_latitude
  end if

  if(Grid%f_plane) then 
    write(stdoutunit,*) '==>NOTE: Setting geometry and Coriolis according to f-plane' 
  endif 
  if(Grid%beta_plane) then 
    write(stdoutunit,*) '==>NOTE: Setting geometry and Coriolis according to beta-plane' 
  endif 

  if(Grid%tripolar) then 
    write (stdoutunit,'(1x,/a)')  ' ==> Note: Energy conversion errors are nontrivial when using tripolar=.true.'
    write (stdoutunit,'(7x,a)')   'The cause is related to the need to update redundantly computed information' 
    write (stdoutunit,'(7x,a)')   'across the Arctic bipolar fold in a bit-wise exact manner for terms contributing'
    write (stdoutunit,'(7x,a)')   'to the energy conversion analysis.  The extra code and mpp calls have not been'
    write (stdoutunit,'(7x,a)')   'implemented.'
  endif 

  if (PRESENT(grid_name)) then
     Grid%name = grid_name
  else
     Grid%name = 'ocean'
  endif

end subroutine set_ocean_grid_size
! </SUBROUTINE> NAME="set_ocean_grid_size"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_hgrid_arrays">
!
! <DESCRIPTION>
! Define horizontal (and some vertical) grid arrays.
!
!---------------------------------------------------------------------------------------------------------------------
! Grid%       grid_spec      grid_spec     grid_spec                       Description
!  var        field          field         field
!             VERSION_0      VERSION_1     VERSION_2(mosaic) 
!--------------------------------------------------------------------------------------------------------------------- 
!                                                    ocean_vgrid.nc
!                                                    k=1,nk
! zt          zt             zt            zeta(2k-1)
! zw          zw             zb            zeta(2k)
!
!                                                    ocean_hgrid.nc
!                                                    i=1,ni
!                                                    j=1,nj
! grid_x_t    gridlon_t      grid_x_T      x(2i  ,2)
! grid_x_u    gridlon_vert_t grid_x_C      x(2i+1,1)
! grid_y_t    gridlat_t      grid_y_T      y(ni/4,2j)
! grid_y_u    gridlat_vert_t grid_y_C      y(ni/4,2j+1)
!
!T
! xt(i,j)     geolon_t(i,j)  x_T(i,j)      x(2i,2j)
! yt          geolat_t       y_T           y(2i,2j)
! dtw         dtw            ds_01_11_T    dx(2i-1,2j)                 distance to western face of t cell
! dte         dte            ds_11_21_T    dx(2i,2j)                   distance to eastern face of t cell
! dts         dts            ds_10_11_T    dy(2i,2j-1)                 distance to southern face of t cell 
! dtn         dtn            ds_11_12_T    dy(2i,2j)                   distance to northern face of t cell
! 
! dxt         dxt            ds_01_21_T    dx(2i,2j)    +dx(2i-1,2j)   width of t cell   
! dxtn        dxtn           ds_02_22_T    dx(2i-1,2j+1)+dx(2i,2j+1)   width of northern face of t cell
! dxte        dxte           ds_00_20_C    dx(2i,2j)    +dx(2i+1,2j)   distance to adjacent t cell to the east!
! dyt         dyt            ds_10_12_T    dy(2i,2j)    +dy(2i,2j-1)   height of t cell
! dytn        dytn           ds_00_02_C    dy(2i,2j)    +dy(2i,2j+1)   distance to adjacent t cell to the north!
! dyte        dyte           ds_20_22_T    dy(2i+1,2j-1)+dy(2i+1,2j)   height of eastern face of t cell 
!
!C 
! NOTE: The "first" (I,J) C-cell is the one shifted NE of the "first" (I,J) T-cell
!            
!
! xu          geolon_c       x_C           x(2i+1,2j+1)
! yu          geolat_c       y_c           y(2i+1,2j+1)
! dxu         dxu            ds_01_21_C    dx(2i+1,2j+1)+dx(2i,2j+1)   width of u cell
! dxun        dxun           ds_02_22_C    dx(2i,2j+2)+dx(2i+1,2j+2)   width of northern face of u cell
! dyu         dyu            ds_10_12_C    dy(2i+1,2j+1)+dy(2i+1,2j)   height of u cell
! dyue        dyue           ds_20_22_C    dy(2i+2,2j)+dy(2i+2,2j+1)   height of eastern face of u cell
!
! dyun        dyun           ds_11_12_C    dy(2i+1,2j+1)+dy(2i+1,2j+2) distance to adjacent u cell to the north 
!                            +ds_10_11_C(i,j+1)                         satisfies sum rule dyte(i,j)=dyun(i,j-1) 
! dxue        dxue           ds_11_21_C    dx(2i+1,2j+1)+dx(2i+2,2j+1) distance to adjacent u cell to the east! 
!                            +ds_01_11_C(i+1,j)   
!
! duw         duw            ds_01_11_C    dx(2i,2j+1)                 distance to western face of u cell
! due         due            ds_11_21_C    dx(2i+1,2j+1)               distance to eastern face of u cell
! dus         dus            ds_10_11_C    dy(2i+1,2j)                 distance to southern face of u cell
! dun         dun            ds_11_12_C    dy(2i+1,2j+1)               distance to northern face of u cell
!  
! sin_rot     sin_rot        angle_C       sin(angle_dx(2*i+1,2*j+1)   sin of rotation angle at corner cell centers
! cos_rot     cos_rot        angle_C       cos(angle_dx(2*i+1,2*j+1)   cos of rotation angle at corner cell centers
!
!Following are the available fields in mosaic files
!--------------------------------------------------------
!Mosaic file     fields
!--------------------------------------------------------
!ocean_hgrid.nc  x, y, dx, dy, angle_dx, area
!ocean_vgrid.nc  zeta
!topog.nc        depth
!
! </DESCRIPTION>
subroutine set_ocean_hgrid_arrays(Domain, Grid)

  type(ocean_domain_type), intent(inout) :: Domain
  type(ocean_grid_type),   intent(inout) :: Grid

  real, dimension(:),   allocatable          :: data
  real, dimension(:,:), allocatable          :: tmp
  real, dimension(:,:), allocatable          :: tmp_local
  real, dimension(:,:), allocatable          :: tmp1_local
  real, dimension(:,:), allocatable          :: tmp2_local
  real, dimension(:,:), allocatable          :: tmpx
  real, dimension(:,:), allocatable          :: tmpy
  integer, dimension(:),allocatable          :: xext, yext
  type(ocean_domain_type)                    :: Domain2  

  integer :: isd2, ied2, jsd2, jed2
  integer :: isc2, iec2, jsc2, jec2
  integer :: isg2, ieg2, jsg2, jeg2
  integer :: isc1, iec1, jsc1, jec1
  integer :: i, j, k, lon_Tripol
  integer :: ioff, joff
  integer :: start(4), nread(4)

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! set the grid points coordinates (degrees) and grid spacing (degrees)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Domain, isg, ieg, jsg, jeg)
  call get_halo_sizes(Domain, xhalo, yhalo)

  ni = Grid%ni;nj = Grid%nj; nk = Grid%nk
 
  allocate (Grid%xt(isd:ied,jsd:jed))
  allocate (Grid%yt(isd:ied,jsd:jed))
  allocate (Grid%xu(isd:ied,jsd:jed))
  allocate (Grid%yu(isd:ied,jsd:jed))
  allocate (Grid%grid_x_t(ni))
  allocate (Grid%grid_y_t(nj))
  allocate (Grid%grid_x_u(ni))
  allocate (Grid%grid_y_u(nj))

  allocate (Grid%zt(nk))
  allocate (Grid%zw(nk))

  allocate(Grid%phit(isd:ied,jsd:jed))
  allocate(Grid%phiu(isd:ied,jsd:jed))
  allocate(Grid%h1t(isd:ied,jsd:jed))
  allocate(Grid%h2t(isd:ied,jsd:jed))
  allocate(Grid%h1u(isd:ied,jsd:jed))
  allocate(Grid%h2u(isd:ied,jsd:jed))

  allocate (Grid%dxt(isd:ied,jsd:jed))
  allocate (Grid%dxu(isd:ied,jsd:jed))
  allocate (Grid%dyt(isd:ied,jsd:jed))
  allocate (Grid%dyu(isd:ied,jsd:jed))
  allocate (Grid%dat(isd:ied,jsd:jed))
  allocate (Grid%dat_frac(isd:ied,jsd:jed))
  allocate (Grid%dau(isd:ied,jsd:jed))

  allocate (Grid%dxtn(isd:ied,jsd:jed))
  allocate (Grid%dytn(isd:ied,jsd:jed))
  allocate (Grid%dxte(isd:ied,jsd:jed))
  allocate (Grid%dyte(isd:ied,jsd:jed))
  allocate (Grid%dxun(isd:ied,jsd:jed))
  allocate (Grid%dyun(isd:ied,jsd:jed))
  allocate (Grid%dxue(isd:ied,jsd:jed))
  allocate (Grid%dyue(isd:ied,jsd:jed))

  allocate (Grid%dxtr(isd:ied,jsd:jed))
  allocate (Grid%dxur(isd:ied,jsd:jed))
  allocate (Grid%dytr(isd:ied,jsd:jed))
  allocate (Grid%dyur(isd:ied,jsd:jed))
  allocate (Grid%dxuer(isd:ied,jsd:jed))
  allocate (Grid%dyuer(isd:ied,jsd:jed))
  allocate (Grid%dxunr(isd:ied,jsd:jed))
  allocate (Grid%dyunr(isd:ied,jsd:jed))
  allocate (Grid%dxter(isd:ied,jsd:jed))
  allocate (Grid%dyter(isd:ied,jsd:jed))
  allocate (Grid%dxtnr(isd:ied,jsd:jed))
  allocate (Grid%dytnr(isd:ied,jsd:jed))
  allocate (Grid%dyue_dxuer(isd:ied,jsd:jed))
  allocate (Grid%dxun_dyunr(isd:ied,jsd:jed))
  allocate (Grid%datr(isd:ied,jsd:jed))
  allocate (Grid%daur(isd:ied,jsd:jed))
  allocate (Grid%dater(isd:ied,jsd:jed))
  allocate (Grid%datnr(isd:ied,jsd:jed))

  allocate (Grid%dxt_dxter(isd:ied,jsd:jed))
  allocate (Grid%dxte_dxtr(isd:ied,jsd:jed))
  allocate (Grid%dxt_dxtnr(isd:ied,jsd:jed))
  allocate (Grid%dxtn_dxtr(isd:ied,jsd:jed))
  allocate (Grid%dyt_dyter(isd:ied,jsd:jed))
  allocate (Grid%dyte_dytr(isd:ied,jsd:jed))
  allocate (Grid%dyt_dytnr(isd:ied,jsd:jed))
  allocate (Grid%dytn_dytr(isd:ied,jsd:jed))

  allocate (Grid%dh1dy(isd:ied,jsd:jed))
  allocate (Grid%dh2dx(isd:ied,jsd:jed))

  allocate ( Grid%duw(isd:ied,jsd:jed) )
  allocate ( Grid%due(isd:ied,jsd:jed) )
  allocate ( Grid%dus(isd:ied,jsd:jed) )
  allocate ( Grid%dun(isd:ied,jsd:jed) )
  allocate ( Grid%dtw(isd:ied,jsd:jed) )
  allocate ( Grid%dte(isd:ied,jsd:jed) )
  allocate ( Grid%dts(isd:ied,jsd:jed) )
  allocate ( Grid%dtn(isd:ied,jsd:jed) )

  allocate(Grid%sin_rot(isd:ied,jsd:jed))
  allocate(Grid%cos_rot(isd:ied,jsd:jed))

  allocate (Grid%obc_tmask(isd:ied,jsd:jed))
  allocate (Grid%obc_umask(isd:ied,jsd:jed))

#endif

  !--- initialize grid data
  
  Grid%xt=0.0;    Grid%yt=0.0;     Grid%grid_x_t=0.0; Grid%grid_y_t=0.0 
  Grid%xu=0.0;    Grid%yu=0.0;     Grid%grid_x_u=0.0; Grid%grid_y_u=0.0
  Grid%dxtn=1.0;  Grid%dytn=1.0;   Grid%dxte=1.0;     Grid%dyte=1.0
  Grid%dxun=1.0;  Grid%dyun=1.0;   Grid%dxue=1.0;     Grid%dyue=1.0
  Grid%duw=0.0 ;  Grid%due=0.0 ;   Grid%dus=0.0 ;     Grid%dun=0.0 
  Grid%dtw=0.0 ;  Grid%dte=0.0 ;   Grid%dts=0.0;      Grid%dtn=0.0
  Grid%zt=0.0;    Grid%zw=0.0;     Grid%sin_rot=0.0;  Grid%cos_rot=0.0  
  Grid%dxt=epsln; Grid%dxu=epsln ; Grid%dyt=epsln ;   Grid%dyu=epsln
  Grid%dh1dy=0.0; Grid%dh2dx=0.0

  call get_domain_offsets(Domain,ioff, joff)  
  select case( grid_version )
  case( VERSION_0 )
     call read_data(ocean_hgrid, "zt",             Grid%zt,       no_domain = .true.)
     call read_data(ocean_hgrid, "zw",             Grid%zw,       no_domain = .true.)
     call read_data(ocean_hgrid, "gridlon_t",      Grid%grid_x_t, no_domain = .true.)
     call read_data(ocean_hgrid, "gridlat_t",      Grid%grid_y_t, no_domain = .true.)
     allocate(data(ni+1))
     call read_data(ocean_hgrid, "gridlon_vert_t", data,          no_domain = .true.)
     Grid%grid_x_u = data(2:ni+1)
     deallocate(data)
     allocate(data(nj+1))
     call read_data(ocean_hgrid, "gridlat_vert_t", data,          no_domain = .true.)
     Grid%grid_y_u = data(2:nj+1)
     deallocate(data)
  case( VERSION_1 )
     call read_data(ocean_hgrid, "zt",             Grid%zt,       no_domain = .true.)
     call read_data(ocean_hgrid, "zb",             Grid%zw,       no_domain = .true.)
     call read_data(ocean_hgrid, "grid_x_T",      Grid%grid_x_t, no_domain = .true.)
     call read_data(ocean_hgrid, "grid_y_T",      Grid%grid_y_t, no_domain = .true.)
     call read_data(ocean_hgrid, "grid_x_C",      Grid%grid_x_u, no_domain = .true.)
     call read_data(ocean_hgrid, "grid_y_C",      Grid%grid_y_u, no_domain = .true.)
  case( VERSION_2 )
     allocate(data(2*nk+1) )  
     call read_data(ocean_vgrid, "zeta", data, no_domain=.TRUE.)
     do k=1,nk
        Grid%zt(k) = data(2*k)  
        Grid%zw(k) = data(2*k+1)  
     enddo
     deallocate(data)
     !--- set up domain for supergrid.
     allocate( xext(Domain%layout(1)), yext(Domain%layout(2)) )
     call mpp_get_domain_extents(domain%domain2d, xext, yext)
     xext = 2*xext
     yext = 2*yext

     call mpp_define_domains((/1,2*ni,1,2*nj/),domain%layout, Domain2%domain2d, maskmap=Domain%maskmap&
           , xflags = Domain%xflags, yflags = Domain%yflags, xhalo=2, yhalo=2,name='ocean_super_grid'&
           , xextent = xext, yextent = yext, symmetry = .true. )

     call mpp_get_compute_domain(Domain%domain2d, isc1, iec1, jsc1, jec1)
     call mpp_get_compute_domain(Domain2%domain2d, isc2, iec2, jsc2, jec2)
     call mpp_get_data_domain(Domain2%domain2d, isd2, ied2, jsd2, jed2)
     call mpp_get_global_domain(Domain2%domain2d, isg2, ieg2, jsg2, jeg2)
     !--- The following adjust is for the consideration of static memory.
     ioff = isc1 - isc; joff = jsc1 - jsc
     isc2 = isc2 - 2*ioff; iec2 = iec2 - 2*ioff
     jsc2 = jsc2 - 2*joff; jec2 = jec2 - 2*joff
     isd2 = isd2 - 2*ioff; ied2 = ied2 - 2*ioff
     jsd2 = jsd2 - 2*joff; jed2 = jed2 - 2*joff

     if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec ) then
        call mpp_error(FATAL, 'ocean_grids_mod (set_ocean_hgrid_arrays): supergrid compute domain is not set properly')
     endif     
     if(isd2 .NE. 2*isd-1 .OR. ied2 .NE. 2*ied .OR. jsd2 .NE. 2*jsd-1 .OR. jed2 .NE. 2*jed ) then
        print*, isd, ied, jsd, jed, isd2, ied2, jsd2, jed2
        call mpp_error(FATAL, 'ocean_grids_mod (set_ocean_hgrid_arrays): supergrid data domain is not set properly')
     endif  

     Domain2%isc = isc2; Domain2%iec = iec2
     Domain2%jsc = jsc2; Domain2%jec = jec2
     Domain2%isd = isd2; Domain2%ied = ied2
     Domain2%jsd = jsd2; Domain2%jed = jed2
     Domain2%isg = isg2; Domain2%ieg = ieg2
     Domain2%jsg = jsg2; Domain2%jeg = jeg2
     Domain2%xhalo = 2;  Domain2%yhalo = 2
     Domain2%ioff = 2*ioff; Domain2%joff = 2*joff

     start = 1; nread = 1
     start(2) = 2; nread(1) = 2*ni+1; nread(2) = 2
     allocate(tmpx(2*ni+1,2), tmpy(2,2*nj+1) )
     call read_data(ocean_hgrid, "x", tmpx, start, nread, no_domain=.true.)

     do i = 1, ni
        Grid%grid_x_t(i) = tmpx(2*i,  1) 
        Grid%grid_x_u(i) = tmpx(2*i+1,2) 
     enddo     
     lon_Tripol = ni/4 ! 90 for 1 degree grid
     start = 1; nread = 1
     start(1) = 2*lon_Tripol+1; nread(1) = 2; nread(2) = 2*nj+1
     call read_data(ocean_hgrid, "y", tmpy, start, nread, no_domain=.true.)

     do j = 1, nj  
        Grid%grid_y_t(j) = tmpy(1 , 2*j)
        Grid%grid_y_u(j) = tmpy(1 , 2*j+1)
     enddo
     deallocate(tmpx, tmpy)
     allocate(tmpx(isc2:iec2+1, jsc2:jec2+1))
     allocate(tmpy(isc2:iec2+1, jsc2:jec2+1))
     call read_data(ocean_hgrid, "x", tmpx, Domain2%domain2d, position=CORNER)
     call read_data(ocean_hgrid, "y", tmpy, Domain2%domain2d, position=CORNER)
  end select

  !--- Grid%xt
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolon_t', Grid%xt, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'x_T', Grid%xt, Domain%domain2d)   
  case(VERSION_2)
     Grid%xt(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%xt,0,0,0,0)

  !--- Grid%yt
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolat_t', Grid%yt, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'y_T', Grid%yt, Domain%domain2d)   
  case(VERSION_2)
     Grid%yt(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%yt,0,0,0,0)

  !--- Grid%xu
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolon_c', Grid%xu, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'x_C', Grid%xu, Domain%domain2d)   
  case(VERSION_2)
     Grid%xu(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%xu,-1,-1,0,0)

  !--- Grid%yu
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolat_c', Grid%yu, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'y_C', Grid%yu, Domain%domain2d)   
  case(VERSION_2)
     Grid%yu(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%yu,-1,-1,0,0)

  if (Grid%beta_plane .or. Grid%f_plane) then
     Grid%phiu(:,:) = Grid%f_plane_latitude/radian
     Grid%phit(:,:) = Grid%phiu(:,:)
  else
     Grid%phit(:,:) = Grid%yt(:,:)/radian
     Grid%phiu(:,:) = Grid%yu(:,:)/radian
  endif

 Grid%h1t(:,:) = radius*cos(Grid%phit(:,:))
 where (cos(Grid%phit) == 0.0) Grid%h1t = radius*abs(epsln)
 Grid%h1u(:,:) = radius*cos(Grid%phiu(:,:))
 where (cos(Grid%phiu) == 0.0) Grid%h1u = radius*abs(epsln)
 Grid%h2t(:,:) = radius
 Grid%h2u(:,:) = radius

 ! set cell widths (meters) and area (meters^2)
 if(grid_version == VERSION_2) then
    deallocate(tmpx, tmpy)
    allocate(tmpx(isd2:ied2,jsd2:jed2+1), tmpy(isd2:ied2+1,jsd2:jed2))
    tmpx = 0.0
    tmpy = 0.0
    call read_data(ocean_hgrid, "dx", tmpx, Domain2%domain2d, position=NORTH)
    call read_data(ocean_hgrid, "dy", tmpy, Domain2%domain2d, position=EAST)   
    call update_boundaries(Domain2,Grid,tmpx,0,-1,0,1) 
    call update_boundaries(Domain2,Grid,tmpy,-1,0,1,0)
 end if

  ! --- Grid%dxt
  select case(grid_version) 
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxt', Grid%dxt, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_01_21_T', Grid%dxt, Domain%domain2d)   
  case(VERSION_2)
     Grid%dxt(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2,2*jsc:2*jec:2) + tmpx(2*isc:2*iec:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxt,0,0,0,0,default_data=epsln)

  !--- Grid%dyt
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dyt', Grid%dyt, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_10_12_T', Grid%dyt, Domain%domain2d)   
  case(VERSION_2)
     Grid%dyt(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc-1:2*jec-1:2) + tmpy(2*isc:2*iec:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dyt,0,0,0,0,default_data=epsln)

  ! --- Grid%dxu
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxu', Grid%dxu, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_01_21_C', Grid%dxu, Domain%domain2d)   
  case(VERSION_2)
     Grid%dxu(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxu,-1,-1,0,0,default_data=epsln)

  !--- Grid%dyu
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dyu', Grid%dyu, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_10_12_C', Grid%dyu, Domain%domain2d)   
  case(VERSION_2)
     Grid%dyu(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc:2*jec:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dyu,-1,-1,0,0,default_data=epsln)

  Grid%dat(:,:)  = Grid%dxt(:,:)*Grid%dyt(:,:)
  Grid%dau(:,:)  = Grid%dxu(:,:)*Grid%dyu(:,:)

  ! set lengths at edges of grid cells (meters)
  allocate(tmp1_local(isd:ied,jsd:jed), tmp2_local(isd:ied,jsd:jed))

  ! --- Grid%dxtn
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxtn', Grid%dxtn, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_02_22_T', Grid%dxtn, Domain%domain2d)
  case(VERSION_2)
     Grid%dxtn(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc:2*iec:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxtn,0,-1,0,0,default_data=1.0)

  ! --- Grid%dytn
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dytn', Grid%dytn, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_00_02_C', Grid%dytn, Domain%domain2d)
  case(VERSION_2)
     Grid%dytn(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc:2*jec:2) + tmpy(2*isc:2*iec:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dytn,0,-1,0,0,default_data=1.0)

  ! --- Grid%dxte
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxte', Grid%dxte, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_00_20_C', Grid%dxte, Domain%domain2d)
  case(VERSION_2)
     Grid%dxte(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc:2*jec:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxte,-1,0,0,0,default_data=1.0) 

  ! --- Grid%dyte
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dyte', Grid%dyte, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_20_22_T', Grid%dyte, Domain%domain2d)
  case(VERSION_2)
     Grid%dyte(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc-1:2*jec-1:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dyte,-1,0,0,0,default_data=1.0)

  ! --- Grid%dxun
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxun', Grid%dxun, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_02_22_C', Grid%dxun, Domain%domain2d)
  case(VERSION_2)
     Grid%dxun(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc+2:2*jec+2:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc+2:2*jec+2:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxun,-1,-2,0,0,default_data=1.0)

  ! --- Grid%dyun
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dyun', Grid%dyun, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_10_11_C', tmp2_local, Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_12_C', tmp1_local, Domain%domain2d)
     call update_boundaries(Domain,Grid,tmp2_local,-1,-1,0,0,field_ref=tmp1_local(isc:iec,jsc:jec) )
     Grid%dyun(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc:iec,jsc+1:jec+1) 
  case(VERSION_2)
     Grid%dyun(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc+2:2*jec+2:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dyun,-1,-2,0,0,default_data=1.0)

  ! --- Grid%dxue
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dxue', Grid%dxue, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_01_11_C', tmp2_local, Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_21_C', tmp1_local, Domain%domain2d)
     call update_boundaries(Domain,Grid,tmp2_local,-1,-1,0,0,field_ref=tmp1_local(isc:iec,jsc:jec),default_data=1.0)
     Grid%dxue(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc+1:iec+1,jsc:jec) 
  case(VERSION_2)
     Grid%dxue(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc+2:2*iec+2:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dxue,-2,-1,0,0,default_data=1.0) 

  ! --- Grid%dyue
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dyue', Grid%dyue, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_20_22_C', Grid%dyue, Domain%domain2d)
  case(VERSION_2)
     Grid%dyue(isc:iec,jsc:jec) = tmpy(2*isc+2:2*iec+2:2,2*jsc:2*jec:2) + tmpy(2*isc+2:2*iec+2:2,2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dyue,-2,-1,0,0,default_data=1.0) 

  ! set reciprocals and related quantities

  do j=jsd,jed
     do i=isd,ied
        Grid%dxtr(i,j)       = 1.0/(Grid%dxt(i,j)+epsln)
        Grid%dxur(i,j)       = 1.0/(Grid%dxu(i,j)+epsln)
        Grid%dytr(i,j)       = 1.0/(Grid%dyt(i,j)+epsln)
        Grid%dyur(i,j)       = 1.0/(Grid%dyu(i,j)+epsln)
        Grid%dxuer(i,j)      = 1.0/(Grid%dxue(i,j)+epsln)
        Grid%dyuer(i,j)      = 1.0/(Grid%dyue(i,j)+epsln)
        Grid%dxunr(i,j)      = 1.0/(Grid%dxun(i,j)+epsln)
        Grid%dyunr(i,j)      = 1.0/(Grid%dyun(i,j)+epsln)
        Grid%dxter(i,j)      = 1.0/(Grid%dxte(i,j)+epsln)
        Grid%dyter(i,j)      = 1.0/(Grid%dyte(i,j)+epsln)
        Grid%dxtnr(i,j)      = 1.0/(Grid%dxtn(i,j)+epsln)
        Grid%dytnr(i,j)      = 1.0/(Grid%dytn(i,j)+epsln)
        Grid%dyue_dxuer(i,j) = Grid%dyue(i,j)/(Grid%dxue(i,j)+epsln)
        Grid%dxun_dyunr(i,j) = Grid%dxun(i,j)/(Grid%dyun(i,j)+epsln)
        Grid%datr(i,j)       = 1.0/(Grid%dat(i,j)+epsln)
        Grid%daur(i,j)       = 1.0/(Grid%dau(i,j)+epsln)
        Grid%dater(i,j)      = 1.0/(Grid%dxte(i,j)*Grid%dyte(i,j)+epsln)
        Grid%datnr(i,j)      = 1.0/(Grid%dxtn(i,j)*Grid%dytn(i,j)+epsln)

        Grid%dxt_dxter(i,j)  = Grid%dxt(i,j)/(Grid%dxte(i,j)+epsln)
        Grid%dxte_dxtr(i,j)  = Grid%dxte(i,j)/(Grid%dxt(i,j)+epsln)
        Grid%dxt_dxtnr(i,j)  = Grid%dxt(i,j)/(Grid%dxtn(i,j)+epsln)
        Grid%dxtn_dxtr(i,j)  = Grid%dxtn(i,j)/(Grid%dxt(i,j)+epsln)
        Grid%dyt_dyter(i,j)  = Grid%dyt(i,j)/(Grid%dyte(i,j)+epsln)
        Grid%dyte_dytr(i,j)  = Grid%dyte(i,j)/(Grid%dyt(i,j)+epsln)
        Grid%dyt_dytnr(i,j)  = Grid%dyt(i,j)/(Grid%dytn(i,j)+epsln)
        Grid%dytn_dytr(i,j)  = Grid%dytn(i,j)/(Grid%dyt(i,j)+epsln)

     enddo
  enddo

  ! set quantities which account for rotation of basis vectors

 allocate(tmp_local(isd:ied,jsd:jed))
 tmp_local=1.0
 
 ! expanded version of dh2dx(:,:) = BDX_EU(tmp(:,:))
 do j=jsc,jec
    do i=isc,iec
       Grid%dh2dx(i,j) = (Grid%dyue(i,j)*tmp_local(i,j) - Grid%dyue(i-1,j)*tmp_local(i-1,j))*Grid%daur(i,j)
    enddo
 enddo
 if (.not. Grid%tripolar .and. jec+joff == nj) then
   do i=isc,iec
     Grid%dh2dx(i,jec) = 0.0
   enddo
 endif

 call update_boundaries(Domain,Grid,Grid%dh2dx,-1,-1,0,0)

 ! expanded version of dh1dy(:,:) = BDY_NU(tmp(:,:))
 do j=jsc,jec
    do i=isc,iec
       Grid%dh1dy(i,j) = (Grid%dxun(i,j)*tmp_local(i,j) - Grid%dxun(i,j-1)*tmp_local(i,j-1))*Grid%daur(i,j)
    enddo
 enddo
 if (.not. Grid%tripolar .and. jec+joff == nj) then
   do i=isc,iec
     Grid%dh1dy(i,jec) = 0.0
   enddo
 endif

 call update_boundaries(Domain,Grid,Grid%dh1dy,-1,-1,0,0)

 ! build distances from grid point to cell faces {meters}
  ! --- Grid%duw
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'duw', Grid%duw, Domain%domain2d )
     call read_data(ocean_hgrid, 'due', tmp_local, Domain%domain2d  )
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_01_11_C', Grid%duw, Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_21_C', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%duw(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc+1:2*jec+1:2)
     tmp_local(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%duw,-1,-1,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%due
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'due', Grid%due, Domain%domain2d )
     call read_data(ocean_hgrid, 'duw', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_11_21_C', Grid%due, Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_01_11_C', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%due(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
     tmp_local(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%due,-1,-1,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dus
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dus', Grid%dus(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dun', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_10_11_C', Grid%dus(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_12_C', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dus(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc:2*jec:2)
     tmp_local(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dus,-1,-1,0,0,field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dun
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dun', Grid%dun(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dus', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_11_12_C', Grid%dun(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_10_11_C', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dun(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
     tmp_local(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dun,-1,-1,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dtw
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dtw', Grid%dtw(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dte', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_01_11_T', Grid%dtw(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_21_T', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dtw(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2, 2*jsc:2*jec:2)
     tmp_local(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dtw,0,0,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dte
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dte', Grid%dte(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dtw', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_11_21_T', Grid%dte(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_01_11_T', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dte(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc:2*jec:2)
     tmp_local(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2, 2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dte,0,0,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dts
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dts', Grid%dts(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dtn', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_10_11_T', Grid%dts(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_11_12_T', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dts(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc-1:2*jec-1:2)
     tmp_local(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc:2*jec:2)
  end select
  call update_boundaries(Domain,Grid,Grid%dts,0,0,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! --- Grid%dtn
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'dtn', Grid%dtn(isc:iec,jsc:jec), Domain%domain2d )
     call read_data(ocean_hgrid, 'dts', tmp_local, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'ds_11_12_T', Grid%dtn(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(ocean_hgrid, 'ds_10_11_T', tmp_local, Domain%domain2d)
  case(VERSION_2)
     Grid%dtn(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc:2*jec:2)
     tmp_local(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc-1:2*jec-1:2)
     deallocate(tmpx, tmpy)
  end select
  call update_boundaries(Domain,Grid,Grid%dtn,0,0,0,0, field_ref=tmp_local(isc:iec,jsc:jec))

  ! calculate rotation angles on velocity points
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'sin_rot', Grid%sin_rot, Domain%domain2d)
     call read_data(ocean_hgrid, 'cos_rot', Grid%cos_rot, Domain%domain2d)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'angle_C', tmp_local, Domain%domain2d)
     tmp_local     = tmp_local*deg_to_rad
     Grid%sin_rot(isc:iec,jsc:jec) = sin(tmp_local(isc:iec,jsc:jec))
     Grid%cos_rot(isc:iec,jsc:jec) = cos(tmp_local(isc:iec,jsc:jec))
  case(VERSION_2)
     allocate(tmp(isc2:iec2+1,jsc2:jec2+1) )
     call read_data(ocean_hgrid, "angle_dx", tmp, Domain2%domain2d, position=CORNER)
     tmp = tmp*deg_to_rad
     do j = jsc, jec
        do i = isc, iec
           Grid%sin_rot(i,j) = sin(tmp(2*i+1,2*j+1) )
           Grid%cos_rot(i,j) = cos(tmp(2*i+1,2*j+1) )
        end do
     end do
     deallocate(tmp)
     call mpp_deallocate_domain(Domain2%domain2d)
  end select
  call update_boundaries(Domain,Grid,Grid%sin_rot,-1,-1,0,0)
  call update_boundaries(Domain,Grid,Grid%cos_rot,-1,-1,0,0)


 ! ensure that all grid factors are set properly at the southern boundary.
 ! many grid factors with index j=0 have no affect within the computational domain 
 ! Even so, for completeness they are explicity set below.
 tmp_local=1.0
 if (jsc+joff == 1) then
   
     call mpp_error(NOTE, '==>Note from ocean_grids_mod (set_ocean_hgrid_arrays): altering U-grid arrays at j=0')

     ! phiu(:,0), h1u(:,0), and dxu(:,0) are within the computational domain and are
     ! not arbitrary. their values follow directly from grid_spec definitions.
     ! dyu(:,0) is arbitrary and has been set equal to dyu(:,1) by update_boundaries.
     ! However, a symmetric ocean domain requires that dyu be symmetric so dyu(:,0) = dyu(:,nj) 

     if (Grid%beta_plane .or. Grid%f_plane) then
         Grid%phiu(:,0) = Grid%f_plane_latitude/radian
     else
         Grid%phiu(:,0) = Grid%yu(:,1)/radian - Grid%dyte(:,1)/radius
     endif
     Grid%phiu(:,0) = max(min(Grid%phiu(:,0),90.0/radian),-90.0/radian)
     where (cos(Grid%phiu(:,0)) == 0.0)
         Grid%h1u(:,0) = radius*abs(epsln)
     elsewhere
         Grid%h1u(:,0) = radius*cos(Grid%phiu(:,0))
     end where
     
     Grid%dxu(:,0)  = (Grid%dxu(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
     Grid%dau(:,0)  = Grid%dxu(:,0)*Grid%dyu(:,0)
     
     Grid%dxur(:,0) = 1.0/(Grid%dxu(:,0)+epsln)
     Grid%dyur(:,0) = 1.0/(Grid%dyu(:,0)+epsln)
     Grid%daur(:,0) = 1.0/(Grid%dau(:,0)+epsln)

     Grid%dun(:,0)  = Grid%dyte(:,1) - Grid%dus(:,1)
     Grid%dus(:,0)  = Grid%dau(:,0)/(Grid%dxu(:,0)+epsln) - Grid%dun(:,0)
     Grid%due(:,0)  = (Grid%due(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
     Grid%duw(:,0)  = Grid%dxu(:,0) - Grid%due(:,0)

     Grid%dxue(:,0)  = (Grid%dxue(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
     Grid%dxun(:,0)  = (Grid%dxun(:,1)/Grid%h1t(:,2))*Grid%h1t(:,1)
     Grid%dxuer(:,0) = 1.0/(Grid%dxue(:,0)+epsln)
     Grid%dxunr(:,0) = 1.0/(Grid%dxun(:,0)+epsln)
                             
     Grid%dyun(:,0)  = Grid%dyte(:,1)
     Grid%dyunr(:,0) = 1.0/(Grid%dyun(:,0)+epsln)
     Grid%dxun_dyunr(:,0) = Grid%dxun(:,0)/(Grid%dyun(:,0)+epsln)

     ! not correct, but any use should be zeroed out by the land mask at j=0
     Grid%dyue(:,0)  = Grid%dyu(:,0) 
     Grid%dyuer(:,0) = 1.0/(Grid%dyue(:,0)+epsln) 
     Grid%dyue_dxuer(:,0) = Grid%dyue(:,0)/(Grid%dxue(:,0)+epsln)

     Grid%dh2dx(:,0) = 0.0
     Grid%dh1dy(:,1) = (Grid%dxun(:,1)*tmp_local(:,1) - Grid%dxun(:,0)*tmp_local(:,0))*Grid%daur(:,1)
     Grid%dh1dy(:,0) = 0.0

 endif

  deallocate(tmp_local, tmp1_local, tmp2_local)
 
  if (debug_this_module .or. verbose_init) then
     call write_chksum_2d('xt', Grid%xt(COMP))
     call write_chksum_2d('xu', Grid%xu(COMP))
     call write_chksum_2d('yt', Grid%yt(COMP))
     call write_chksum_2d('yu', Grid%yu(COMP))
     call write_chksum_2d('dxt', Grid%dxt(COMP))
     call write_chksum_2d('dxu', Grid%dxu(COMP))
     call write_chksum_2d('dyt', Grid%dyt(COMP))
     call write_chksum_2d('dyu', Grid%dyu(COMP))
     call write_chksum_2d('dat', Grid%dat(COMP))
     call write_chksum_2d('dau', Grid%dau(COMP))
     call write_chksum_2d('dxtn', Grid%dxtn(COMP))
     call write_chksum_2d('dytn', Grid%dytn(COMP))
     call write_chksum_2d('dxte', Grid%dxte(COMP))
     call write_chksum_2d('dyte', Grid%dyte(COMP))
     call write_chksum_2d('dxun', Grid%dxun(COMP))
     call write_chksum_2d('dyun', Grid%dyun(COMP))
     call write_chksum_2d('dxue', Grid%dxue(COMP))
     call write_chksum_2d('dyue', Grid%dyue(COMP))
     call write_chksum_2d('dtn+dts', Grid%dtn(COMP) + Grid%dts(COMP))
     call write_chksum_2d('dun+dus', Grid%dun(COMP) + Grid%dus(COMP))
     call write_chksum_2d('dte+dtw', Grid%dte(COMP) + Grid%dtw(COMP))
     call write_chksum_2d('due+duw', Grid%due(COMP) + Grid%duw(COMP))
     call write_chksum_2d('dte', Grid%dte(COMP))
     call write_chksum_2d('dtw', Grid%dtw(COMP))
     call write_chksum_2d('due', Grid%due(COMP))
     call write_chksum_2d('duw', Grid%duw(COMP))
     call write_chksum_2d('sin_rot', Grid%sin_rot(COMP))
     call write_chksum_2d('cos_rot', Grid%cos_rot(COMP))
  endif

  if(write_grid) then
     call write_data("ocean_grid_check.nc", "xt", Grid%xt(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "yt", Grid%yt(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "xu", Grid%xu(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "yu", Grid%yu(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "zt", Grid%zt,                  no_domain=.true.)
     call write_data("ocean_grid_check.nc", "zw", Grid%zw,                  no_domain=.true.)
     call write_data("ocean_grid_check.nc", "dxt", Grid%dxt(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dyt", Grid%dyt(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dxu", Grid%dxu(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dyu", Grid%dyu(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dat", Grid%dat(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dau", Grid%dau(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dxtn", Grid%dxtn(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dytn", Grid%dytn(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dxte", Grid%dxte(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dyte", Grid%dyte(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dxun", Grid%dxun(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dyun", Grid%dyun(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dxue", Grid%dxue(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dyue", Grid%dyue(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dtn", Grid%dtn(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dts", Grid%dts(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dun", Grid%dun(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dus", Grid%dus(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dte", Grid%dte(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "dtw", Grid%dtw(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "due", Grid%due(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "duw", Grid%duw(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "sin_rot", Grid%sin_rot(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "cos_rot", Grid%cos_rot(isc:iec,jsc:jec), Domain%domain2d)
  end if

end subroutine set_ocean_hgrid_arrays
! </SUBROUTINE> NAME="set_ocean_hgrid_arrays"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_vgrid_arrays">
!
! <DESCRIPTION>
! Compute vertical (and some horizontal) grids for ocean model. 
! Also compute axes information for diagnostic manager.  
! </DESCRIPTION>
!
subroutine set_ocean_vgrid_arrays (Domain, Grid, obc)

  type(ocean_domain_type), intent(inout)        :: Domain
  type(ocean_grid_type),   intent(inout)        :: Grid
  logical,                 intent(in), optional :: obc

  real    :: ht_fax, ht_fay, ht_jp1_fax, ht_ip1_fay
  real    :: tcella_sphere, ucella_sphere
  real    :: wet_t_points, wet_u_points

  integer :: i, j, k, kzt, kzu, kmin

  logical :: have_obc=.false.

  real    :: sigma_sign, zwnk, zwnk_r

  integer :: stdoutunit 
  stdoutunit=stdout() 

  
  if (PRESENT(obc)) have_obc = obc
  
#ifndef MOM_STATIC_ARRAYS
  allocate (Grid%tcella(nk))
  allocate (Grid%ucella(nk))
  allocate (Grid%mask(isd:ied,jsd:jed,nk))
  allocate (Grid%tmask(isd:ied,jsd:jed,nk))
  allocate (Grid%umask(isd:ied,jsd:jed,nk))
  allocate (Grid%tmasken(isd:ied,jsd:jed,nk,2))
  allocate (Grid%tmask_depth(isd:ied,jsd:jed,nk))
  allocate (Grid%umask_depth(isd:ied,jsd:jed,nk))

  allocate (Grid%dht_dx(isd:ied,jsd:jed))
  allocate (Grid%dht_dy(isd:ied,jsd:jed))
  allocate (Grid%gradH(isd:ied,jsd:jed))

  allocate (Grid%dzt(nk))
  allocate (Grid%dztlo(nk))
  allocate (Grid%dztup(nk))
  allocate (Grid%dzw(0:nk))

  allocate (Grid%st(nk))
  allocate (Grid%sw(nk))
  allocate (Grid%dst(nk))
  allocate (Grid%dstlo(nk))
  allocate (Grid%dstup(nk))
  allocate (Grid%dsw(0:nk))

  allocate (Grid%fracdz(nk,0:1)) 
  allocate (Grid%dzwr(0:nk))

#endif

  ! set depth arrays
  Grid%dzw(0)         = Grid%zt(1)
  Grid%dzw(1)         = Grid%zt(2)-Grid%zt(1)
  Grid%dzwr(0)        = 1.0/Grid%dzw(0)
  Grid%dzwr(1)        = 1.0/Grid%dzw(1)
  Grid%dzt(1)         = Grid%zw(1)
  Grid%fracdz(1,0)    = Grid%zt(1)/Grid%dzt(1)  
  Grid%fracdz(1,1)    = (Grid%zw(1) - Grid%zt(1))/Grid%dzt(1)  
  do k=2,nk-1
     Grid%dzt(k)      = Grid%zw(k)-Grid%zw(k-1)
     Grid%dzw(k)      = Grid%zt(k+1)-Grid%zt(k)
     Grid%dzwr(k)     = 1.0/Grid%dzw(k)
     Grid%fracdz(k,0) = (Grid%zt(k) - Grid%zw(k-1))/Grid%dzt(k)
     Grid%fracdz(k,1) = (Grid%zw(k) - Grid%zt(k))/Grid%dzt(k)
  enddo
  Grid%dzt(nk)        = Grid%zw(nk)-Grid%zw(nk-1)
  Grid%dzw(nk)        = Grid%zw(nk) - Grid%zt(nk)  
  Grid%dzwr(nk)       = 1.0/Grid%dzw(nk)
  Grid%fracdz(nk,0)   = (Grid%zt(nk) - Grid%zw(nk-1))/Grid%dzt(nk)
  Grid%fracdz(nk,1)   = (Grid%zw(nk) - Grid%zt(nk))/Grid%dzt(nk)

  ! half-thicknesses 
  k=1
  Grid%dztlo(k) = Grid%zw(k)-Grid%zt(k)
  Grid%dztup(k) = Grid%dzw(k-1)
  do k=2,nk
     Grid%dztlo(k) = Grid%zw(k)-Grid%zt(k)
     Grid%dztup(k) = Grid%zt(k)-Grid%zw(k-1)
  enddo


 ! construct T cell and U cell land/sea masks
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          if (Grid%kmt(i,j) .ge. k) then
             Grid%tmask(i,j,k) = 1.0
          else
             Grid%tmask(i,j,k) = 0.0
          endif
          if (Grid%kmu(i,j) .ge. k) then
             Grid%umask(i,j,k) = 1.0
          else
             Grid%umask(i,j,k) = 0.0
          endif
       enddo
    enddo
 enddo

 if(horz_grid == MOM_CGRID) then 
    Grid%mask(:,:,:) = Grid%tmask(:,:,:)
    do k=1,nk
       do j=jsd,jec
          do i=isd,iec
             Grid%tmasken(i,j,k,1) = min(Grid%tmask(i,j,k),Grid%tmask(i+1,j,k))
             Grid%tmasken(i,j,k,2) = min(Grid%tmask(i,j,k),Grid%tmask(i,j+1,k))
          enddo
       enddo
    enddo
    call mpp_update_domains(Grid%tmasken(:,:,:,1),Domain%domain2d)
    call mpp_update_domains(Grid%tmasken(:,:,:,2),Domain%domain2d)
 else 
    Grid%mask(:,:,:) = Grid%umask(:,:,:)
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             Grid%tmasken(i,j,k,1) = Grid%umask(i,j,k)
             Grid%tmasken(i,j,k,2) = Grid%umask(i,j,k)
          enddo
       enddo
    enddo
 endif 


 ! Compute masks as if depth or pressure were the vertical coordinates.
 ! This new mask is used for case when have vert_coordinate=ZSIGMA or PSIGMA,
 ! in which case tmask and umask are 1 for all k=1,nk except where ht=0.0.
 ! We need tmask_depth and umask_depth for vertical remapping in case when 
 ! remap fields into depth or pressure space from terrain following s-space. 
 Grid%tmask_depth(:,:,:) = 0.0
 Grid%umask_depth(:,:,:) = 0.0

 k=1
 do j=jsd,jed
    do i=isd,ied
       if (Grid%ht(i,j) > 0.0) then
           Grid%tmask_depth(i,j,k) = 1.0
       endif
       if (Grid%hu(i,j) > 0.0) then
           Grid%umask_depth(i,j,k) = 1.0
       endif
    enddo
 enddo
 do k=1,nk-1
    do j=jsd,jed
       do i=isd,ied
          if (Grid%ht(i,j) >= Grid%zw(k)) then
              Grid%tmask_depth(i,j,k+1) = 1.0
          endif
          if (Grid%hu(i,j) >= Grid%zw(k)) then
              Grid%umask_depth(i,j,k+1) = 1.0
          endif
       enddo
    enddo
 enddo

  ! compute surface area and volume of ocean (T cells and U cells)
  Grid%tcellv    = 0.0 !total ocean volume on T cells (assuming eta=0)
  Grid%ucellv    = 0.0 !total ocean volume on U cells (assuming eta=0)
  Grid%tcella(:) = 0.0 !total ocean surface area on T cells in level k
  Grid%ucella(:) = 0.0 !total ocean surface area on U cells in level k

  ! exclude points at open boundaries from diagnostics
  do j=jsc,jec
    do i=isc,iec
       kzt = Grid%kmt(i,j)
       if(have_obc) kzt = kzt * int(Grid%obc_tmask(i,j))
       if (kzt .gt. 0) then
           do k=1,kzt
              Grid%tcella(k) = Grid%tcella(k) + Grid%dat(i,j)
           enddo
           Grid%tcellv = Grid%tcellv + Grid%dat(i,j)*Grid%ht(i,j)
       endif
    enddo
 enddo

 do j=jsc,jec
    do i=isc,iec
       kzu = Grid%kmu(i,j)
       if(have_obc) kzu = kzu * int(Grid%obc_umask(i,j))
       if (kzu .gt. 0) then
           do k=1,kzu
              Grid%ucella(k) = Grid%ucella(k) + Grid%dau(i,j)
           enddo
           Grid%ucellv = Grid%ucellv + Grid%dau(i,j)*Grid%hu(i,j)
       endif
    enddo
 enddo

 ! compute area of full domain, including land regions 
 tcella_sphere=0.0
 ucella_sphere=0.0
 do j=jsc,jec
    do i=isc,iec
       tcella_sphere = tcella_sphere + Grid%dat(i,j)
       ucella_sphere = ucella_sphere + Grid%dau(i,j)
    enddo
 enddo

 call mpp_sum (Grid%tcella, nk)
 call mpp_sum (Grid%ucella, nk)
 call mpp_sum (Grid%tcellv)
 call mpp_sum (Grid%ucellv)
 call mpp_sum (tcella_sphere)
 call mpp_sum (ucella_sphere)

 ! information about model grid size  
 Grid%total_t_points = Grid%ni*Grid%nj*Grid%nk
 if(global_sum_flag==BITWISE_EXACT_SUM) then 
     Grid%wet_t_points = nint(mpp_global_sum(Domain%domain2d,Grid%tmask(:,:,:), BITWISE_EXACT_SUM)) 
     Grid%wet_u_points = nint(mpp_global_sum(Domain%domain2d,Grid%umask(:,:,:), BITWISE_EXACT_SUM)) 
 else 
     wet_t_points = real(mpp_global_sum(Domain%domain2d,Grid%tmask(:,:,:), NON_BITWISE_EXACT_SUM), KIND=4)
     wet_u_points = real(mpp_global_sum(Domain%domain2d,Grid%umask(:,:,:), NON_BITWISE_EXACT_SUM), KIND=4)
     Grid%wet_t_points = nint(wet_t_points)
     Grid%wet_u_points = nint(wet_u_points)
 endif

 if (debug_this_module .or. verbose_init) then
     write (stdoutunit,'(/a,i12)')   '  Number of wet ocean tracer points          =', Grid%wet_t_points
     write (stdoutunit,'(a,i12)')    '  Number of wet ocean velocity points        =', Grid%wet_u_points
     write (stdoutunit,'(a,i12)')    '  Number of computed ocean tracer points     =', Grid%total_t_points

     write (stdoutunit,'(/a,es24.17,a)') &
     '  Wet ocean volume with eta_t=0.0   (T-cells)  =', Grid%tcellv,    ' m^3 (not bit reproducible)'
     write (stdoutunit,'(a,es24.17,a)') &
      '  Ocean surface area               (T-cells)  =', Grid%tcella(1), ' m^2 (not bit reproducible)'

     write (stdoutunit,'(/a,es24.17,a)') &
     '  Wet ocean volume with eta_u=0.0   (U-cells)  =', Grid%ucellv,    ' m^3 (not bit reproducible)'
     write (stdoutunit,'(a,es24.17,a)') &
     '  Wet ocean surface area            (U-cells)  =', Grid%ucella(1), ' m^2 (not bit reproducible)'

     write (stdoutunit,'(/a,es24.17,a)') &
     '  Wet ocean + masked-out (land) surface area (T-cells) =', tcella_sphere,  ' m^2 (not bit reproducible)'
     write (stdoutunit,'(a,es24.17,a/)') &
     '  Wet ocean + masked-out (land) surface area (U-cells) =', ucella_sphere,  ' m^2 (not bit reproducible)'
 endif
 
 ! compute bitwise reproducible surface areas
 if(global_sum_flag==BITWISE_EXACT_SUM) then 
     if(have_obc) then 
         Grid%tcellsurf = &
         mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1)*Grid%obc_tmask(:,:), BITWISE_EXACT_SUM)
         Grid%ucellsurf = &
         mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1)*Grid%obc_umask(:,:), BITWISE_EXACT_SUM)
     else 
         Grid%tcellsurf = &
          mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1), BITWISE_EXACT_SUM)
         Grid%ucellsurf = &
          mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1), BITWISE_EXACT_SUM)
     endif
 else 
     if(have_obc) then 
         Grid%tcellsurf = &
         real(mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1)*Grid%obc_tmask(:,:), NON_BITWISE_EXACT_SUM), KIND=4) 
         Grid%ucellsurf = &
         real(mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1)*Grid%obc_umask(:,:), NON_BITWISE_EXACT_SUM), KIND=4) 
     else 
         Grid%tcellsurf = &
         real(mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1), NON_BITWISE_EXACT_SUM), KIND=4) 
         Grid%ucellsurf = &
         real(mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1), NON_BITWISE_EXACT_SUM), KIND=4) 
     endif
 endif


 ! compute fraction of total wet area occupied by a single tracer cell
 do j=jsc,jec
    do i=isc,iec
       Grid%dat_frac(i,j) = Grid%dat(i,j)/Grid%tcellsurf
    enddo
 enddo

  ! compute topographic slopes at U-cell. Note: cannot use 
  ! operators as ocean_operators_init has not yet been called.
  do j=jsc,jec
    do i=isc,iec
      ht_fax           = 0.5*( Grid%ht(i,j) + Grid%ht(i+1,j))
      ht_fay           = 0.5*( Grid%ht(i,j) + Grid%ht(i,j+1))
      ht_jp1_fax       = 0.5*( Grid%ht(i,j+1) + Grid%ht(i+1,j+1))
      ht_ip1_fay       = 0.5*( Grid%ht(i+1,j) + Grid%ht(i+1,j+1))
      Grid%dht_dx(i,j) = (ht_ip1_fay-ht_fay)*Grid%dxur(i,j)*Grid%umask(i,j,1)
      Grid%dht_dy(i,j) = (ht_jp1_fax-ht_fax)*Grid%dyur(i,j)*Grid%umask(i,j,1)
      Grid%gradH(i,j)  = sqrt(Grid%dht_dx(i,j)**2 + Grid%dht_dy(i,j)**2)
    enddo
  enddo
  call update_boundaries(Domain,Grid,Grid%dht_dx,0,0,0,0)
  call update_boundaries(Domain,Grid,Grid%dht_dy,0,0,0,0)
  call update_boundaries(Domain,Grid,Grid%gradH,0,0,0,0)

  allocate(rho0_profile(nk))
  rho0_profile(:) = rho0
  if(read_rho0_profile) then 
      write(stdoutunit,'(a)') '==>Warning: ocean_grids_mod: read rho0 profile to set dst grid spacing.'
      write(stdoutunit,'(a)') '            This option is experimental, and NOT generally supported.'
      if(vert_coordinate_class==DEPTH_BASED) then 
          write(stdoutunit,'(a)') '   Since using DEPTH_BASED vertical coordinates, rho0_profile = rho0.'
      endif
      if(vert_coordinate_class==PRESSURE_BASED) then 
         call read_data('INPUT/rho0_profile.nc','rho0_profile', rho0_profile)
      endif 
  endif 

  ! compute static s-grid information 
  if(vert_coordinate == GEOPOTENTIAL .or. vert_coordinate == ZSTAR) then 
    do k=1,nk
      Grid%st(k)    = Grid%zt(k)
      Grid%sw(k)    = Grid%zw(k)
      Grid%dst(k)   = Grid%dzt(k)
      Grid%dstlo(k) = Grid%dztlo(k)
      Grid%dstup(k) = Grid%dztup(k)
    enddo
    do k=0,nk
      Grid%dsw(k) = Grid%dzw(k)
    enddo

  elseif(vert_coordinate == PRESSURE .or. vert_coordinate == PSTAR) then 
    do k=1,nk
      Grid%st(k)    =  rho0_profile(k)*grav*Grid%zt(k)
      Grid%sw(k)    =  rho0_profile(k)*grav*Grid%zw(k)
      Grid%dst(k)   = -rho0_profile(k)*grav*Grid%dzt(k)
      Grid%dstlo(k) = -rho0_profile(k)*grav*Grid%dztlo(k)
      Grid%dstup(k) = -rho0_profile(k)*grav*Grid%dztup(k)
    enddo
    do k=0,nk
      kmin = max(1,k)
      Grid%dsw(k)   = -rho0_profile(kmin)*grav*Grid%dzw(k)
    enddo

  elseif(vert_coordinate == ZSIGMA .or. vert_coordinate == PSIGMA) then 

     ! assume initial bottom pressure = rho0*g*ht, which means that  
     ! PSIGMA and ZSIGMA initialization are same, to within sign. 
     ! -1 <= ZSIGMA <= 0, so dst and dsw are positive 
     !  0 <= PSIGMA <= 1, so dst and dsw are negative 
     ! Recall s-coordinate is dimensionless when use ZSIGMA or PSIGMA.

     if(vert_coordinate==ZSIGMA) then  
         sigma_sign = -1.0
     elseif(vert_coordinate==PSIGMA) then
         sigma_sign =  1.0
     endif

     zwnk   = Grid%zw(nk)
     zwnk_r = 1.0/(epsln + zwnk)

     do k=1,nk
        Grid%st(k) = sigma_sign*Grid%zt(k)*zwnk_r
        Grid%sw(k) = sigma_sign*Grid%zw(k)*zwnk_r
     enddo

     ! set s-grid increments from st and sw
     Grid%dsw(0) = - Grid%st(1)
     Grid%dsw(1) = -(Grid%st(2)-Grid%st(1))
     Grid%dst(1) = - Grid%sw(1)
     do k=2,nk-1
        Grid%dst(k) = -(Grid%sw(k)  - Grid%sw(k-1))
        Grid%dsw(k) = -(Grid%st(k+1)- Grid%st(k))
     enddo
     Grid%dst(nk) = -(Grid%sw(nk)-Grid%sw(nk-1))
     Grid%dsw(nk) = -(Grid%sw(nk)-Grid%st(nk))

     k=1
     Grid%dstlo(k) = Grid%sw(k)-Grid%st(k)
     Grid%dstup(k) = Grid%dsw(k-1)
     do k=2,nk
        Grid%dstlo(k) = Grid%sw(k)-Grid%st(k)
        Grid%dstup(k) = Grid%st(k)-Grid%sw(k-1)
     enddo

  endif 

  ! set up the axis definition for variables
  call axes_info(Domain, Grid)


end subroutine set_ocean_vgrid_arrays
! </SUBROUTINE> NAME="set_ocean_vgrid_arrays"

!#######################################################################
! <SUBROUTINE NAME="init_grids_diag">
!
! <DESCRIPTION>
! Initialize some grid diagnostics. 
! </DESCRIPTION>
subroutine init_grids_diag(Grid, Time)

  type(ocean_grid_type), intent(inout) :: grid
  type(ocean_time_type), intent(in) :: Time

  integer :: id_geolon_t
  integer :: id_geolat_t
  integer :: id_geolon_uv
  integer :: id_geolat_uv
  integer :: id_area_t
  integer :: id_area_u

  ! output grid information to the diagnostic manager

  id_geolon_t = register_static_field('ocean_model','geolon_t', Grid%tracer_axes(1:2), &
                                      'tracer longitude','degrees_E', range=(/-281.,361./))
  if (id_geolon_t > 0)  used = send_data(id_geolon_t,Grid%xt(isc:iec,jsc:jec),  Time%model_time)

  id_geolat_t = register_static_field('ocean_model','geolat_t', Grid%tracer_axes(1:2), &
                                       'tracer latitude','degrees_N', range=(/-91.,91./))  
  if (id_geolat_t > 0)  used = send_data(id_geolat_t,Grid%yt(isc:iec,jsc:jec),  Time%model_time)

  id_geolon_uv = register_static_field('ocean_model','geolon_c', Grid%vel_axes_uv(1:2), &
                                       'uv longitude','degrees_E', range=(/-281.,361./))
  if (id_geolon_uv > 0) used = send_data(id_geolon_uv,Grid%xu(isc:iec,jsc:jec), Time%model_time)

  id_geolat_uv = register_static_field('ocean_model','geolat_c',Grid%vel_axes_uv(1:2), &
                                       'uv latitude','degrees_N', range=(/-91.,91./))
  if (id_geolat_uv > 0) used = send_data(id_geolat_uv,Grid%yu(isc:iec,jsc:jec), Time%model_time)

  id_area_t = register_static_field('ocean_model','area_t',Grid%tracer_axes(1:2), 'tracer cell area',&
                                    'm^2', range=(/0.,1.e15/))
  if (id_area_t > 0)    used = send_data(id_area_t,Grid%dat(isc:iec,jsc:jec),   Time%model_time)

  id_area_u = register_static_field('ocean_model','area_u',Grid%vel_axes_uv(1:2), 'velocity cell area',&
                                    'm^2', range=(/0.,1.e15/))  
  if (id_area_u > 0)    used = send_data(id_area_u,Grid%dau(isc:iec,jsc:jec),   Time%model_time)    

  id_kmt = register_static_field ('ocean_model', 'kmt', Grid%tracer_axes(1:2), &
           'number of depth levels on t-grid', 'dimensionless',                &
            missing_value=missing_value, range=(/-1e1,1e9/))
  if (id_kmt > 0) then 
      wrk1_2d(:,:) = Grid%kmt(:,:)
      call diagnose_2d(Time, Grid, id_kmt, wrk1_2d(:,:))
  endif 
  id_kmu = register_static_field ('ocean_model', 'kmu', Grid%vel_axes_uv(1:2), &
           'number of depth levels on u-grid', 'dimensionless',                &
            missing_value=missing_value, range=(/-1e1,1e9/))
  if (id_kmu > 0) then 
      wrk1_2d(:,:) = Grid%kmu(:,:)
      call diagnose_2d(Time, Grid, id_kmu, wrk1_2d(:,:))
  endif 

  id_ht  = register_static_field ('ocean_model', 'ht', Grid%tracer_axes(1:2),      &
                                  'ocean depth on t-cells', 'm',                   &
                                  missing_value=missing_value, range=(/-1e9,1e9/), &
                                  standard_name='sea_floor_depth_below_geoid')
  call diagnose_2d(Time, Grid, id_ht, Grid%ht(:,:))

  id_hu  = register_static_field ('ocean_model', 'hu', Grid%vel_axes_uv(1:2), 'ocean depth on u-cells', &
                                  'm',missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_hu, Grid%hu(:,:))

  id_dht_dx  = register_static_field ('ocean_model', 'dht_dx', Grid%vel_axes_uv(1:2), 'd(ht)/dx on u-cells', &
                                      'm/m', missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_dht_dx, Grid%dht_dx(:,:))

  id_dht_dy  = register_static_field ('ocean_model', 'dht_dy', Grid%vel_axes_uv(1:2), 'd(ht)/dy on u-cells', &
                                      'm/m', missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_dht_dy, Grid%dht_dy(:,:))

  id_gradH = register_static_field ('ocean_model', 'gradH', Grid%vel_axes_uv(1:2), '|grad bottom on U-cell|', &
                                      'm/m', missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_gradH, Grid%gradH(:,:))

  id_dxt = register_static_field ('ocean_model', 'dxt', Grid%tracer_axes(1:2), 'ocean dxt on t-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d(Time, Grid, id_dxt, Grid%dxt(:,:))

  id_dxu = register_static_field ('ocean_model', 'dxu', Grid%vel_axes_uv(1:2), 'ocean dxu on u-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_dxu, Grid%dxu(:,:))

  id_dyt = register_static_field ('ocean_model', 'dyt', Grid%tracer_axes(1:2), 'ocean dyt on t-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d(Time, Grid, id_dyt, Grid%dyt(:,:))

  id_dyu = register_static_field ('ocean_model', 'dyu', Grid%vel_axes_uv(1:2), 'ocean dyu on u-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d_u(Time, Grid, id_dyu, Grid%dyu(:,:))

  id_dxtn = register_static_field ('ocean_model', 'dxtn', Grid%tracer_axes(1:2), 'ocean dxtn on t-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d(Time, Grid, id_dxtn, Grid%dxtn(:,:))

  id_dyte = register_static_field ('ocean_model', 'dyte', Grid%tracer_axes(1:2), 'ocean dyte on t-cells', 'm',&
                                  missing_value=missing_value, range=(/-1e9,1e9/))
  call diagnose_2d(Time, Grid, id_dyte, Grid%dyte(:,:))

  id_tmask = register_static_field ('ocean_model', 'tmask', Grid%tracer_axes(1:3), 'tracer mask', 'dimensionless',&
                                  missing_value=missing_value, range=(/-1e1,1e1/))
  call diagnose_3d(Time, Grid, id_tmask, Grid%tmask(:,:,:))

  id_umask = register_static_field ('ocean_model', 'umask', Grid%vel_axes_uv(1:3), 'velocity mask', 'dimensionless',&
                                  missing_value=missing_value, range=(/-1e1,1e1/))
  call diagnose_3d_u(Time, Grid, id_umask, Grid%umask(:,:,:))

  id_tmasken(1) = register_static_field ('ocean_model', 'tmasken1', Grid%tracer_axes(1:3), &
                 'min(tmask(i),tmask(i+1))', 'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))
  call diagnose_3d(Time, Grid, id_tmasken(1), Grid%tmasken(:,:,:,1))

  id_tmasken(2) = register_static_field ('ocean_model', 'tmasken2', Grid%tracer_axes(1:3), &
                 'min(tmask(j),tmask(j+1))', 'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))
  call diagnose_3d(Time, Grid, id_tmasken(2), Grid%tmasken(:,:,:,2))


end subroutine init_grids_diag
! </SUBROUTINE> NAME="init_grids_diag"



!#######################################################################
! <SUBROUTINE NAME="axes_info">
!
! <DESCRIPTION>
! Set up axes definitions. 
! </DESCRIPTION>
!
subroutine axes_info(Domain, Grid)

  type(ocean_domain_type), intent(in)    :: Domain
  type(ocean_grid_type),   intent(inout) :: Grid

  integer :: id_xt, id_yt, id_xu, id_yu
  integer :: id_zt_bounds, id_zt, id_zw_bounds, id_zw
  integer :: id_st_bounds, id_st, id_sw_bounds, id_sw
  integer :: k

  real, allocatable, dimension(:) :: zt_bounds
  real, allocatable, dimension(:) :: zw_bounds
  real, allocatable, dimension(:) :: st_bounds
  real, allocatable, dimension(:) :: sw_bounds

  real :: meter2dbar 
  meter2dbar = c2dbars*rho0*grav 


  ! horizontal axes
  ! Note: 1d grid arrays are not in general representative of the grid.
  ! Need to generalize diagnostics to employ CF conventions.

  id_xt = diag_axis_init ('xt_'//trim(Grid%name),Grid%grid_x_t,'degrees_E','x','tcell longitude',&
       set_name='ocean', Domain2=Domain%domain2d, aux='geolon_t')
  id_yt = diag_axis_init ('yt_'//trim(Grid%name),Grid%grid_y_t,'degrees_N','y','tcell latitude',&
       set_name='ocean', Domain2=Domain%domain2d, aux='geolat_t')
  id_xu = diag_axis_init ('xu_'//trim(Grid%name),Grid%grid_x_u,'degrees_E','x','ucell longitude',&
       set_name='ocean', Domain2=Domain%domain2d, aux='geolon_c')
  id_yu = diag_axis_init ('yu_'//trim(Grid%name),Grid%grid_y_u,'degrees_N','y','ucell latitude',&
       set_name='ocean', Domain2=Domain%domain2d, aux='geolat_c')


  ! vertical axes (needs to be re-thought for partial cells and general vertical coordinates)

  allocate(zt_bounds(Grid%nk+1))
  allocate(zw_bounds(Grid%nk+1))
  zt_bounds(1) = Grid%zt(1)-Grid%dzt(1)/2.0
  zw_bounds(1) = Grid%zw(1)-Grid%dzw(1)/2.0
  do k=2,Grid%nk+1
     zt_bounds(k)=zt_bounds(k-1)+Grid%dzt(k-1)
     zw_bounds(k)=zw_bounds(k-1)+Grid%dzw(k-1)
  enddo

  allocate(st_bounds(Grid%nk+1))
  allocate(sw_bounds(Grid%nk+1))

  if(vert_coordinate_class==DEPTH_BASED) then 
      ! for depth-based, st and sw are > 0 and dst and dsw are > 0
      st_bounds(1) = Grid%st(1)-Grid%dst(1)/2.0
      sw_bounds(1) = Grid%sw(1)-Grid%dsw(1)/2.0
      do k=2,Grid%nk+1
         st_bounds(k)=st_bounds(k-1)+Grid%dst(k-1)
         sw_bounds(k)=sw_bounds(k-1)+Grid%dsw(k-1)
      enddo
  else
      ! for pressure-based, st and sw are > 0 but dst and dsw are < 0
      st_bounds(1) = Grid%st(1)+Grid%dst(1)/2.0
      sw_bounds(1) = Grid%sw(1)+Grid%dsw(1)/2.0
      do k=2,Grid%nk+1
         st_bounds(k)=st_bounds(k-1)-Grid%dst(k-1)
         sw_bounds(k)=sw_bounds(k-1)-Grid%dsw(k-1)
      enddo
  endif


  if(vert_coordinate==GEOPOTENTIAL) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), zt_bounds, 'meters', 'z',&
           'tcell depth edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), Grid%zt, 'meters', 'z','tcell depth',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), zw_bounds, 'meters', 'z',&
           'ucell depth edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), Grid%zw, 'meters', 'z','ucell depth',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), st_bounds, 'meters', 'z',&
           'tcell depth edges',direction=-1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), Grid%st, 'meters', 'z','tcell depth',&
           edges = id_st_bounds,direction=-1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), sw_bounds, 'meters', 'z',&
           'ucell depth edges',direction=-1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), Grid%sw, 'meters', 'z','ucell depth',&
           edges = id_sw_bounds,direction=-1, set_name='ocean')

  elseif(vert_coordinate==ZSTAR) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), zt_bounds, 'meters', 'z',&
           'tcell zstar depth edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), Grid%zt, 'meters', 'z','tcell zstar depth',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), zw_bounds, 'meters', 'z',&
           'ucell zstar depth edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), Grid%zw, 'meters', 'z','ucell zstar depth',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), st_bounds, 'meters', 'z',&
           'tcell zstar depth edges',direction=-1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), Grid%st, 'meters', 'z','tcell zstar depth',&
           edges = id_st_bounds,direction=-1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), sw_bounds, 'meters', 'z',&
           'ucell zstar depth edges',direction=-1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), Grid%sw, 'meters', 'z','ucell zstar depth',&
           edges = id_sw_bounds,direction=-1, set_name='ocean')

  elseif(vert_coordinate==ZSIGMA) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), zt_bounds, 'meters', 'z',&
           'tcell depth edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), Grid%zt, 'meters', 'z','tcell depth',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), zw_bounds, 'meters', 'z',&
           'ucell depth edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), Grid%zw, 'meters', 'z','ucell depth',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), st_bounds, 'dimensioness', 'z',&
           'tcell sigma-depth edges',direction=1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), Grid%st, 'dimensionless', 'z','tcell sigma-depth',&
           edges = id_st_bounds,direction=1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), sw_bounds, 'dimensionless', 'z',&
           'ucell sigma-depth edges',direction=1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), Grid%sw, 'dimensionless', 'z','ucell sigma-depth',&
           edges = id_sw_bounds,direction=1, set_name='ocean')


  elseif(vert_coordinate==PRESSURE) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), meter2dbar*zt_bounds, 'dbars', 'z',&
           'tcell pressure edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), meter2dbar*Grid%zt, 'dbars', 'z','tcell pressure',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), meter2dbar*zw_bounds, 'dbars', 'z',&
           'ucell pressure edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), meter2dbar*Grid%zw, 'dbars', 'z','ucell pressure',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), c2dbars*st_bounds, 'dbars', 'z',&
           'tcell pressure edges',direction=-1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), c2dbars*Grid%st, 'dbars', 'z','tcell pressure',&
           edges = id_st_bounds,direction=-1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), c2dbars*sw_bounds, 'dbars', 'z',&
           'ucell pressure edges',direction=-1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), c2dbars*Grid%sw, 'dbars', 'z','ucell pressure',&
           edges = id_sw_bounds,direction=-1, set_name='ocean')

  elseif(vert_coordinate==PSTAR) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), meter2dbar*zt_bounds, 'dbars', 'z',&
           'tcell pstar edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), meter2dbar*Grid%zt, 'dbars', 'z','tcell pstar',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), meter2dbar*zw_bounds, 'dbars', 'z',&
           'ucell pstar edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), meter2dbar*Grid%zw, 'dbars', 'z','ucell pstar',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), c2dbars*st_bounds, 'dbars', 'z',&
           'tcell pstar edges',direction=-1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), c2dbars*Grid%st, 'dbars', 'z','tcell pstar',&
           edges = id_st_bounds,direction=-1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), c2dbars*sw_bounds, 'dbars', 'z',&
           'ucell pstar edges',direction=-1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), c2dbars*Grid%sw, 'dbars', 'z','ucell pstar',&
           edges = id_sw_bounds,direction=-1, set_name='ocean')


  elseif(vert_coordinate==PSIGMA) then 

      id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), meter2dbar*zt_bounds, 'dbars', 'z',&
           'tcell pressure edges',direction=-1, set_name='ocean')

      id_zt = diag_axis_init ('zt_'//trim(Grid%name), meter2dbar*Grid%zt, 'dbars', 'z','tcell pressure',&
           edges = id_zt_bounds,direction=-1, set_name='ocean')

      id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), meter2dbar*zw_bounds, 'dbars', 'z',&
           'ucell pressure edges',direction=-1, set_name='ocean')

      id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), meter2dbar*Grid%zw, 'dbars', 'z','ucell pressure',&
           edges = id_zw_bounds,direction=-1, set_name='ocean')


      id_st_bounds = diag_axis_init ('st_edges_'//trim(Grid%name), st_bounds, 'dimensionless', 'z',&
           'tcell pressure-sigma edges',direction=-1, set_name='ocean')

      id_st = diag_axis_init ('st_'//trim(Grid%name), Grid%st, 'dimensionless', 'z',&
           'tcell pressure-sigma', edges = id_st_bounds,direction=-1, set_name='ocean')

      id_sw_bounds = diag_axis_init ( 'sw_edges_'//trim(Grid%name), sw_bounds, 'dimensionless', 'z', &
            'ucell pressure-sigma edges',direction=-1, set_name='ocean')

      id_sw = diag_axis_init ( 'sw_'//trim(Grid%name), Grid%sw, 'dimensionless', 'z',&
           'ucell pressure-sigma', edges = id_sw_bounds, direction=-1, set_name='ocean')

  endif


  ! attributes for variables

  Grid%tracer_axes        = (/ id_xt, id_yt, id_st /)
  Grid%tracer_axes_flux_x = (/ id_xu, id_yt, id_st /)
  Grid%tracer_axes_flux_y = (/ id_xt, id_yu, id_st /)
  Grid%tracer_axes_wt     = (/ id_xt, id_yt, id_sw /)

  Grid%tracer_axes_depth        = (/ id_xt, id_yt, id_zt /)
  Grid%tracer_axes_flux_x_depth = (/ id_xu, id_yt, id_zt /)
  Grid%tracer_axes_flux_y_depth = (/ id_xt, id_yu, id_zt /)
  Grid%tracer_axes_wt_depth     = (/ id_xt, id_yt, id_zw /)

  Grid%vel_axes_uv              = (/ id_xu, id_yu, id_st /)
  Grid%vel_axes_wu_depth        = (/ id_xu, id_yu, id_zw /)
  Grid%vel_axes_wt_depth        = (/ id_xt, id_yt, id_zw /)
  Grid%vel_axes_uv_depth        = (/ id_xu, id_yu, id_zt /)
  Grid%vel_axes_flux_x_depth    = (/ id_xt, id_yu, id_zt /)
  Grid%vel_axes_flux_y_depth    = (/ id_xu, id_yt, id_zt /)

  Grid%vel_axes_wu              = (/ id_xu, id_yu, id_sw /)
  Grid%vel_axes_wt              = (/ id_xt, id_yt, id_sw /)
  Grid%vel_axes_flux_x          = (/ id_xt, id_yu, id_st /)
  Grid%vel_axes_flux_y          = (/ id_xu, id_yt, id_st /)

  if(horz_grid == MOM_BGRID) then 
     Grid%vel_axes_u = (/ id_xu, id_yu, id_st /)
     Grid%vel_axes_v = (/ id_xu, id_yu, id_st /)
  else 
     Grid%vel_axes_u = (/ id_xu, id_yt, id_st /)
     Grid%vel_axes_v = (/ id_xt, id_yu, id_st /)
  endif 

  deallocate(zt_bounds, zw_bounds)
  deallocate(st_bounds, sw_bounds)

end subroutine axes_info
! </SUBROUTINE> NAME="axes_info"


!#######################################################################
! <SUBROUTINE NAME="update_boundaries">
!
! <DESCRIPTION>
! Set halo points at model boundaries equal to values at boundaries 
! if no grid connectivity exists. If model is connected along
! boundary (e.g., tripolar) then use mpp_update_domains.
! </DESCRIPTION>
!
subroutine update_boundaries(Domain, Grid, field, i_offset, j_offset, ishift, jshift, field_ref, default_data)  
type(ocean_domain_type),                            intent(inout) :: Domain
type(ocean_grid_type),                              intent(in)    :: Grid
real, dimension(Domain%isd:,Domain%jsd:),           intent(inout) :: field
integer,                                            intent(in)    :: i_offset, j_offset
integer,                                            intent(in)    :: ishift, jshift
real, dimension(Domain%isc:,Domain%jsc:), optional, intent(inout) :: field_ref
real,                                     optional, intent(in)    :: default_data

integer                                                      :: i,j
real, dimension(Domain%isd:Domain%ied+ishift,Domain%jsd:Domain%jed+jshift) :: tmp
real, dimension(:,:), allocatable                                          :: tmp2
integer                                                      :: isc, iec, jsc, jec
integer                                                      :: isd, ied, jsd, jed
integer                                                      :: isg, ieg, jsg, jeg
integer                                                      :: xhalo, yhalo
real                                                         :: default_value

 isc = Domain%isc; iec = Domain%iec
 jsc = Domain%jsc; jec = Domain%jec
 isd = Domain%isd; ied = Domain%ied
 jsd = Domain%jsd; jed = Domain%jed
 isg = Domain%isg; ieg = Domain%ieg
 jsg = Domain%jsg; jeg = Domain%jeg
 xhalo = Domain%xhalo; yhalo = Domain%yhalo
 default_value = 0
 if(present(default_data)) default_value = default_data

 if(ishift == 0 .AND. jshift == 0) then
    call mpp_update_domains(field,Domain%domain2d)
 else
    if(ishift == 1 .AND. jshift == 0) then
       allocate(tmp2(isd:ied+ishift,jsd:jed+jshift) )
       tmp2 = field
       call mpp_update_domains(tmp2, Domain%domain2d, position = EAST)
    else if(ishift == 0 .AND. jshift == 1) then
       allocate(tmp2(isd+ishift:ied+ishift,jsd+jshift:jed+jshift) )
       tmp2(isd+ishift:ied+ishift,jsd+jshift:jed+jshift) = field(isd+ishift:ied+ishift,jsd+jshift:jed+jshift)
       call mpp_update_domains(tmp2, Domain%domain2d)
    else
       call mpp_error(FATAL, "ocean_grids_mod(update_boundaries): invalid value of ishift and/or jshift")
    endif
    field(isd+ishift:ied+ishift,jsd+jshift:jed+jshift) = tmp2(isd+ishift:ied+ishift,jsd+jshift:jed+jshift)
    deallocate(tmp2)
 endif

!    call mpp_update_domains(field(isd+ishift:ied+ishift,jsd+jshift:jed+jshift), Domain%domain2d)
 iec = iec+ishift; jec = jec+jshift
 ied = ied+ishift; jed = jed+jshift
 ieg = ieg+ishift; jeg = jeg+jshift


 if (Grid%tripolar) then
    tmp = default_value
    if(i_offset == 0 .AND. j_offset == 0) then
       if(present(field_ref)) then
          tmp(isc:iec,jsc:jec) = field_ref(isc:iec,jsc:jec)
       else
          tmp(isc:iec,jsc:jec) = field(isc:iec,jsc:jec)
       endif 
       call mpp_update_domains(tmp,Domain%domain2d)    
    else if(i_offset == -1 .AND. j_offset == -1) then
       if(present(field_ref)) then
          tmp(isc:iec,jsc:jec) = field_ref(isc:iec,jsc:jec)
       else
          tmp(isc:iec,jsc:jec) = field(isc:iec,jsc:jec)
       endif
       call mpp_update_domains(tmp,Domain%domain2d,position=CORNER) 
    else if(i_offset == -1 .AND. j_offset == 0) then
       tmp(isc:iec,jsc:jec) = field(isc:iec,jsc:jec)
       call mpp_update_domains(tmp,Domain%domain2d,position=EAST)
    else if(i_offset == 0 .AND. j_offset == -1) then
       tmp(isc:iec,jsc:jec) = field(isc:iec,jsc:jec)
       call mpp_update_domains(tmp,Domain%domain2d,position=NORTH)
    else if(i_offset == -1 .AND. j_offset == -2) then
       tmp(isc:iec,jsc+1:jec) = field(isc:iec,jsc:jec-1)
       call mpp_update_domains(tmp, Domain%domain2d,position=CORNER)
    else if(i_offset == -2 .AND. j_offset == -1) then 
       tmp(isc:iec,jsc:jec) = field(isc-1:iec-1,jsc:jec)
       call mpp_update_domains(tmp, Domain%domain2d,position=CORNER)
    else 
       call mpp_error(FATAL, "ocean_grids_mod(update_boundaries): invalid value of i_offset and/or j_offset")
    endif
     if(jec + Domain%joff == jeg) field(isd:ied,jec+1:jed) = tmp(isd:ied,jec+1:jed)
 endif

! fill e/w halo points along computational rows
if (iec + Domain%ioff == ieg .and. .not. Grid%cyclic_x) then
   do i=iec+1,iec+xhalo
      field(i,jsc:jec) = field(iec,jsc:jec)
   enddo
endif
if (isc + Domain%ioff == isg .and. .not. Grid%cyclic_x) then
   do i=isc-xhalo,isc-1
      field(i,jsc:jec) = field(isc,jsc:jec)
   enddo
endif

! fill southernmost halo rows
if (jsc + Domain%joff == jsg .and. .not. Grid%cyclic_y) then
   do j=jsc-yhalo,jsc-1
      field(:,j) = field(:,jsc)
   enddo
endif

! fill northernmost halo rows
if (jec + Domain%joff == jeg) then
   if (.not. Grid%tripolar .and. .not. Grid%cyclic_y) then
      do j=jec+1,jec+yhalo
         field(:,j) = field(:,jec)
      enddo
   endif
endif

! fill halo points at west boundary along south halo rows 
if (.not. (Grid%cyclic_x .or. Grid%cyclic_y) .and. isc + Domain%ioff == isg .and. jsc + Domain%joff /= jsg ) then 
   do j=jsc-yhalo,jsc-1
      do i=isc-xhalo,isc-1
         field(i,j) = field(isc,j)
      enddo
   enddo
endif

! fill halo points at west boundary along north halo rows
if (.not. (Grid%cyclic_x .or. Grid%cyclic_y) .and. isc + Domain%ioff == isg .and. jec + Domain%joff /= jeg ) then 
   do j=jec+1,jec+yhalo
      do i=isc-xhalo,isc-1
         field(i,j) = field(isc,j)
      enddo
   enddo
endif

! fill halo points at east boundary along south halo rows
if (.not. (Grid%cyclic_x .or. Grid%cyclic_y) .and. iec + Domain%ioff == ieg .and. jsc + Domain%joff /= jsg) then 
   do j=jsc-yhalo,jsc-1
      do i=iec+1,iec+xhalo
         field(i,j) = field(iec,j)
      enddo
   enddo
endif

! fill halo points east boundary along north halo rows (I know, it is confusing !!)
if (.not. (Grid%cyclic_x .or. Grid%cyclic_y) .and. iec + Domain%ioff == ieg .and. jec + Domain%joff /= jeg) then 
   do j=jec+1,jec+yhalo
      do i=iec+1,iec+xhalo
         field(i,j) = field(iec,j)
      enddo
   enddo
endif

return

end subroutine update_boundaries
! </SUBROUTINE> NAME="update_boundaries"

end module ocean_grids_mod

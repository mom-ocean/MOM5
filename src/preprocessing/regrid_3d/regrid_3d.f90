program regrid_3d
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  !
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Bonnie Samuels </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Zhi Liang</CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">M.J. Harrison</REVIEWER>
  !
  !<DESCRIPTION>
  ! regrid 3-d lat-lon gridded data to logically rectangular grid
  ! described by grid descriptor file. Applies only to scalar fields
  ! No missing points allowed on input grid.
  !
  !</DESCRIPTION>

  use mpp_mod,          only : mpp_error, mpp_pe,  mpp_npes, mpp_root_pe
  use mpp_mod,          only : FATAL, WARNING, stdout, stdlog, mpp_chksum
  use mpp_io_mod,       only : mpp_open, mpp_close, mpp_read, mpp_write, mpp_write_meta
  use mpp_io_mod,       only : mpp_copy_meta, axistype, fieldtype, atttype
  use mpp_io_mod,       only : mpp_get_atts, mpp_get_info, mpp_get_fields, mpp_get_times
  use mpp_io_mod,       only : mpp_get_axes, mpp_get_axis_data
  use mpp_io_mod,       only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
  use mpp_domains_mod,  only : mpp_update_domains, mpp_define_domains, mpp_global_field
  use mpp_domains_mod,  only : domain2d, mpp_define_layout, mpp_get_compute_domain
  use mpp_domains_mod,  only : mpp_domains_set_stack_size
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type
  use axis_utils_mod,   only : get_axis_cart, interp_1d
  use fms_mod,          only : fms_init, fms_end, open_namelist_file, close_file, file_exist
  use fms_mod,          only : check_nml_error, write_version_number, lowercase
  use constants_mod,    only : constants_init, PI

  implicit none

  integer, parameter :: max_fields       = 10
  integer, parameter :: max_ntimes_saved = 12

  !--- namelist interface
  !<NAMELIST NAME="regrid_3d_nml">
  ! <DATA NAME="src_file" TYPE="character(len=128)" DEFAULT="src_file.nc">
  !  Name of input file containing grid and data to be regridded.
  ! </DATA>
  ! <DATA NAME="numfields"  TYPE="integer" DEFAULT="2">
  !  Number of fields.
  ! </DATA>
  ! <DATA NAME="src_field_name" TYPE="character(len=128), dimension(max_fields)" >
  !  Name of input field(s). default is (/'temp', 'salt'/)
  ! </DATA>
  ! <DATA NAME="dest_field_name" TYPE="character(len=128), dimension(max_fields)" >
  !  Name of output field(s). If it is not specified in the namelist, it will 
  !  get the value from src_field_name
  ! </DATA>
  ! <DATA NAME="dest_grid" TYPE="character(len=128)" DEFAULT="dest_grid.nc">
  !  Name of grid descriptor file containing target grid information.
  ! </DATA>
  ! <DATA NAME="dest_file"  TYPE="character(len=128)" DEFAULT="dest_file.nc">
  !  Name of output file.
  ! </DATA>
  ! <DATA NAME="num_nbrs"  TYPE="integer" DEFAULT="10">
  ! Number of nearest neighbors for regridding  
  ! </DATA>
  ! <DATA NAME="max_dist"  TYPE="integer" DEFAULT="0.17" UNITS="radians">
  !  Maximum radial influence for regridding.
  ! </DATA>
  ! <DATA NAME="scale_factor" TYPE="real ">
  ! scaling factor for data (e.g. -1 to flip sign or 0.01 to convert from centimeters)
  ! </DATA>
  ! <DATA NAME="stop_crit"  TYPE="character(len=1),dimension(2)" DEFAULT="0.001">
  !  The stopping criteria when extrapping data onto missing points.
  ! </DATA>
  ! <DATA NAME="use_source_vertical_grid" TYPE="logical" DEFAULT=".false.">
  !  when use_source_vertical_grid is set to true, the destination data will 
  !  have the same vertical level as the source data. When use_source_vertical_grid 
  !  is false, the vertical grid of destination data will come from dest_grid. 
  !  A linear vertical interpolation will be done when the source vertical is different
  !  from destination vertical grid.
  ! </DATA>
  ! <DATA NAME="apply_mask"  TYPE="logical" DEFAULT="true">
  !  flag to indicate if the land/sea mask of source/destination grid will be applied 
  !  on the output dest_file. When apply_mask is false, the destination data will be 
  !  global data, i.e. no missing value in the destination data file. When apply_mask 
  !  is true, mask will be applied to the destination data. The mask can be either 
  !  source grid or destination grid determined by nml use_source_vertical_grid. 
  !  When use_source_vertical_grid is true, source grid mask will be applied, otherwise
  !  destination grid mask will be applied.
  ! </DATA>
  ! <DATA NAME="interp_method"  TYPE= "character(len=20)" >
  !  specifying the remapping method when remampping data onto current grid.
  !  Its value can be "spherical" or " bilinear". "spherical" interpolation is a 
  !  inverse distance weighted interpolation algorithm. Default value is "bilinear". 
  !  "bilinear" interpolation is recommanded, since bilinear interpolation will provide 
  !  more smooth results than "spherical" interpolation (especially when interpolating 
  !  from coarse grid to fine grid). Plus bilinear interpolation is much more efficiency 
  !  than "spherical interpolation". 
  ! </DATA>
  ! <DATA NAME="ntimes_saved" TYPE= "integer" >
  !  number of time levels to be saved. Its value has to be less than or equal to the number 
  !  of time levels in the source data file.
  ! </DATA>
  ! <DATA NAME="timelevel_saved" TYPE= "integer(max_ntimes_saved)" >
  !  specify the selection of time levels to be saved. The number of elements to be specified 
  !  should be equal to ntimes_saved. 
  ! </DATA>
  ! <DATA NAME="debug" TYPE="logical">
  ! For Debugging. Set true to print out chksum information for debugging reproducing ability 
  ! accross processors. default is false.
  ! </DATA>
  !</NAMELIST>
  character(len=128) :: src_file                          = 'src_file.nc'
  character(len=128) :: dest_grid                         = 'dest_grid.nc'
  character(len=128) :: dest_file                         = 'dest_file.nc'
  integer            :: numfields                         = 2
  character(len=128) :: src_field_name(max_fields)        = ''
  character(len=128) :: dest_field_name(max_fields)       = ''
  real               :: stop_crit(max_fields)             = 0.005 
  integer            :: num_nbrs                          = 5
  real               :: max_dist                          = 0.1
  real               :: scale_factor(max_fields)          = 1.0
  logical            :: use_source_vertical_grid          = .FALSE.
  logical            :: apply_mask                        = .TRUE.
  character(len=32)  :: interp_method                     = "bilinear"
  logical            :: debug                             = .FALSE.
  integer            :: ntimes_saved                      = 0
  integer            :: timelevel_saved(max_ntimes_saved) = 0

  namelist /regrid_3d_nml/ src_file, src_field_name, dest_field_name, numfields, dest_file, &
                           dest_grid, scale_factor, num_nbrs, max_dist, stop_crit,          &
                           interp_method, apply_mask, use_source_vertical_grid, debug,      &
                           ntimes_saved, timelevel_saved

  !---------------------------------------------------------------------
  integer            :: ni_src, nj_src, nk_src, ni_dst, nj_dst, nk_dst, ntime_src
  type(axistype)     :: depth_axis, time_axis, axes_dst(2)
  type(fieldtype)    :: field_lon_dst, field_lat_dst
  type(fieldtype)    :: dest_field(max_fields)
  integer            :: dest_unit, dst_grid_unit         
  logical            :: time_axis_exists = .false.
  real, parameter    :: tol = 1.e-10   ! tolerance for detecting missing values
  real, parameter    :: max_val=1.e20
  real, parameter    :: rel_coef = 0.9
  integer, parameter :: max_iter = 2000
  logical            :: is_cyclic = .true.  ! we suppose the source data is always global data.
  real               :: missing(max_fields)
  real               :: D2R
  type(fieldtype), dimension(:), allocatable :: src_field
  integer                                    :: src_unit
  real, dimension(:),            allocatable :: time_in
  real, dimension(:),            allocatable :: depth_src, depth_dst
  real, dimension(:),            allocatable :: lon_src, lat_src
  real, dimension(:,:),          allocatable :: lon_dst, lat_dst
  real, dimension(:,:,:),        allocatable :: mask_dst
  !--- version information variables
  character(len=128) :: version='CVS $Id: regrid_3d.f90,v 20.0 2013/12/14 00:31:13 fms Exp $'
  character(len=128) :: tagname='Tag $Name: tikal $'

  ! --- Begin of the program

  ! --- call fms_init, which will call mpp_init, mpp_io_init, mpp_domains_init
  call fms_init
  call constants_init

  !--- call regrid_3d initialization routine
  call regrid_3d_init

  !--- read the dest_grid file
  call read_dst_grid

  !--- read src_file
  call read_src_file

  !--- set up metadata of output file
  call setup_meta()

  !--- remap data from src grid to dest grid and write out data to output file
  call process_data()

  call regrid_3d_end

  call fms_end

contains

  !#####################################################################
  ! --- read the namelist and write the version and namelist to logfile. Also
  ! --- write the namelist to standard output
  subroutine regrid_3d_init

    integer :: io_status, unit, ierr, n

    D2R = PI/180.0
    ! --- the default src_field_name is 'temp and 'salt'
    src_field_name(1) = 'temp' 
    src_field_name(2) = 'salt'

    ! --- read namelist ------------------------------------------------
    if(file_exist('input.nml')) then
       unit = open_namelist_file()
       read (unit,regrid_3d_nml,IOSTAT=io_status)
       write (stdout(),'(/)')
       write (stdout(),regrid_3d_nml)  
       ierr = check_nml_error(io_status, 'regrid_3d_nml')
       call close_file(unit)
    else
       call mpp_error(FATAL, 'regrid_3d: file input.nml does not exist' )
    endif

    if(ntimes_saved > 0) then
       if(ntimes_saved > max_ntimes_saved) call mpp_error(FATAL, 'regrid_3d: nml ntimes_saved is greater than ' &
               //'max_ntimes_saved, increase max_ntimes_saved or decrease ntimes_saved')
       do n = 1, ntimes_saved
          if(timelevel_saved(n) .le. 0) call mpp_error(FATAL, 'regrid_3d: nml timelevel_saved is not ' &
              //'specified properly, modify timelevel_saved')
       enddo
    endif

    if(numfields .gt. max_fields) call mpp_error(FATAL, &
             'regrid_3d: numfields should be less than max_fields')
    if (numfields .le. 0) call mpp_error(FATAL,'regrid_3d: No field specified')
    !--- if dest_field_name is not defined in the namelist, get it from src_field_name
    do n = 1, numfields
       if(trim(dest_field_name(n)) == '') dest_field_name(n) = trim(src_field_name(n))
    enddo 
 
    !--- write version information
    call write_version_number(version, tagname)

  end subroutine regrid_3d_init

  !#####################################################################
  !--- open grid file and store grid info
  subroutine read_dst_grid

    integer                                    :: ndim, nvar, natt, ntime, i, j, k
    integer                                    :: len1, siz_in(3)
    logical                                    :: found_xt, found_yt, found_kmt
    character(len=32)                          :: name
    real, allocatable, dimension(:,:)          :: kmt
    type(axistype), allocatable, dimension(:)  :: axes
    type(fieldtype), allocatable, dimension(:) :: fields

    if(.not. file_exist(trim(dest_grid)) ) &
         call mpp_error(FATAL, 'regrid_3d: file '//trim(dest_grid)//' does not exist')

    call mpp_open(dst_grid_unit, trim(dest_grid),&
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    call mpp_get_info(dst_grid_unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar), axes(ndim) )
    call mpp_get_axes(dst_grid_unit, axes)
    call mpp_get_fields(dst_grid_unit,fields)

    !--------------------------------------------------------------------
    ! get output grid information
    !--------------------------------------------------------------------
    ni_dst=0; nj_dst=0; nk_dst=0
    do i=1,ndim
       call mpp_get_atts(axes(i),name=name,len=len1)
       select case (trim(name))
       case ('grid_x_T')
          ni_dst = len1
       case ('grid_y_T')
          nj_dst = len1
       case ('zt')
          nk_dst = len1
          allocate(depth_dst(nk_dst))
          call mpp_get_axis_data(axes(i),depth_dst)
          if(.not.use_source_vertical_grid) depth_axis = axes(i)         
       end select
    enddo
    if(ni_dst==0) call mpp_error(FATAL,'regrid_3d: file '//trim(dest_grid)//' does not contain axis grid_x_T')
    if(nj_dst==0) call mpp_error(FATAL,'regrid_3d: file '//trim(dest_grid)//' does not contain axis grid_y_T')
    if(nk_dst==0) call mpp_error(FATAL,'regrid_3d: file '//trim(dest_grid)//' does not contain axis zt')

    allocate(lon_dst(ni_dst,nj_dst), lat_dst(ni_dst,nj_dst), kmt(ni_dst,nj_dst) )
    found_xt = .FALSE.;  found_yt = .FALSE.; found_kmt = .false.
    do i=1,nvar
       call mpp_get_atts(fields(i),name=name,ndim=ndim)
       select case (trim(name))
       case ('x_T')
          found_xt = .true.
          call mpp_read(dst_grid_unit,fields(i),lon_dst)
          field_lon_dst = fields(i)
          call mpp_get_atts(fields(i),axes=axes_dst)
       case ('y_T')
          found_yt = .true.
          call mpp_read(dst_grid_unit,fields(i),lat_dst)
          field_lat_dst = fields(i)           
       case ('num_levels')
          found_kmt = .true.
          call mpp_read(dst_grid_unit,fields(i),kmt )
       end select
    enddo
    if(.not.found_kmt) call mpp_error(FATAL,'regrid_3d: field num_levels is not in the file '//trim(dest_grid) )
    if(.not.found_xt) call mpp_error(FATAL,'regrid_3d: field x_T is not in the file '//trim(dest_grid) )
    if(.not.found_yt) call mpp_error(FATAL,'regrid_3d: field y_T is not in the file '//trim(dest_grid) )

    allocate(mask_dst(ni_dst,nj_dst,nk_dst) )

    do k = 1, nk_dst
       do j = 1, nj_dst
          do i = 1, ni_dst
             if(kmt(i,j) .ge. k) then
                mask_dst(i,j,k) = 1.0
             else
                mask_dst(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo
    
    deallocate(fields, axes, kmt)    

  end subroutine read_dst_grid

  !#####################################################################
  !--- read source grid informationfrom src_file
  subroutine read_src_file

    integer                                    :: unit, ndim, nvar, natt, n
    integer                                    :: nt, i, j, k, jj, len1
    logical                                    :: found_src_field(numfields)
    character(len=1)                           :: cart
    character(len=32)                          :: name, units
    type(axistype), allocatable, dimension(:)  :: axes
    type(fieldtype), allocatable, dimension(:) :: fields
    real                                       :: scale, add

    if(.not. file_exist(trim(src_file)) ) &
         call mpp_error(FATAL, 'regrid_3d: file '//trim(src_file)//' does not exist')

    call mpp_open(src_unit, trim(src_file),&
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(src_unit, ndim, nvar, natt, ntime_src)

    allocate(fields(nvar))
    call mpp_get_fields(src_unit, fields)
    allocate(src_field(numfields))
    found_src_field = .FALSE.
    do n=1,numfields
       do i=1,nvar
          call mpp_get_atts(fields(i),name=name)
          if (lowercase(trim(src_field_name(n))) == lowercase(trim(name))) then
             src_field(n) = fields(i)
             found_src_field(n) = .TRUE.
             write(stdout(),*) 'Interpolating src field : ',trim(name), ' grid ', trim(dest_grid)
          endif
       end do
    enddo
    do n=1,numfields 
       if(.not. found_src_field(n)) call mpp_error(FATAL, 'regrid_3d: field '&
            //trim(src_field_name(n))//' is not in the file '//trim(src_file) )
    enddo
    !--- get the src grid
    call mpp_get_atts(src_field(1),ndim=ndim)
    allocate(axes(ndim))
    call mpp_get_atts(src_field(1),axes=axes)
    ni_src=0; nj_src=0; nk_src=0; ntime_src = 1
    do j=1,ndim
       call mpp_get_atts(axes(j),len=len1,units=units)
       call get_axis_cart(axes(j),cart)
       select case (cart)
       case ('X')
          ni_src = len1
          allocate(lon_src(ni_src))
          call mpp_get_axis_data(axes(j),lon_src)
       case('Y')
          nj_src = len1
          allocate(lat_src(nj_src))
          call mpp_get_axis_data(axes(j),lat_src)
       case('Z')
          nk_src = len1
          allocate(depth_src(nk_src))
          if(use_source_vertical_grid) depth_axis = axes(j)
          call mpp_get_axis_data(axes(j),depth_src)
          if (trim(units) == 'cm') then
             depth_src(:) = depth_src(:)*.01
          endif
       case ('T')
          ntime_src = len1
          time_axis_exists = .true.
          allocate(time_in(ntime_src))
          call mpp_get_times(src_unit, time_in)
          time_axis = axes(j)
       end select
    enddo
    if(ni_src==0) call mpp_error(FATAL,'regrid_3d: file ' &
         //trim(src_file)//' does not contain axis with cartesian attributes = "X" ')
    if(nj_src==0) call mpp_error(FATAL,'regrid_3d: file '&
         //trim(src_file)//' does not contain axis with cartesian attributes = "Y" ')
    if(nk_src==0) call mpp_error(FATAL,'regrid_3d: file '&
         //trim(src_file)//' does not contain axis with cartesian attributes = "Z" ')
    !--- if use selective time level, make sure the time level selection is suitable.
    if(ntimes_saved > ntime_src) call mpp_error(FATAL, 'regrid_3d: nml ntime_src is greater than ' &
             //'number of time levels in the file '// trim(src_file) )
    do n = 1, ntimes_saved
       if(timelevel_saved(n) > ntime_src) call mpp_error(FATAL, 'regrid_3d: some entry in nml timelevel_saved ' &
             // 'is greater than number of time levels in the file '// trim(src_file) )
    enddo

    !--- get the missing value
    do n=1, numfields
       call mpp_get_atts(src_field(n),missing=missing(n),scale=scale, add=add )
       if(scale .NE. 1 .OR. add .NE. 0) missing(n) = missing(n)*scale + add
    enddo

    deallocate(fields, axes)    

  end subroutine read_src_file

  !#####################################################################
  !--- setup metadata of output file
  subroutine setup_meta

    integer            :: i, j, n, nt
    character(len=32)  :: units
    character(len=128) :: longname
    !--------------------------------------------------------------------
    ! write output file metadata
    !--------------------------------------------------------------------

    call mpp_open(dest_unit, trim(dest_file),MPP_OVERWR,MPP_NETCDF,threading=MPP_SINGLE,&
         fileset=MPP_SINGLE)
    !
    ! write axis metadata
    !
    do i=1,size(axes_dst(:))
       call mpp_copy_meta(dest_unit,axes_dst(i))
    end do

    call mpp_copy_meta(dest_unit,depth_axis)

    if (time_axis_exists) then
       call mpp_copy_meta(dest_unit, time_axis)
    endif

    !
    ! write variable metadata
    !
    call mpp_copy_meta(dest_unit,field_lon_dst,axes=axes_dst)
    call mpp_copy_meta(dest_unit,field_lat_dst,axes=axes_dst)

    do n = 1, numfields
       call mpp_get_atts(src_field(n),units=units,longname=longname)
       if (time_axis_exists) then 
          call mpp_write_meta(dest_unit, dest_field(n), (/axes_dst(1),axes_dst(2), depth_axis,time_axis/), &
               trim(dest_field_name(n)), units,longname, missing=missing(n), pack=1)
       else
          call mpp_write_meta(dest_unit, dest_field(n), (/axes_dst(1),axes_dst(2), depth_axis/), &
               trim(dest_field_name(n)), units,longname, missing=missing(n), pack=1)
       endif
    enddo
    ! write axis data

    do i=1,size(axes_dst(:))
       call mpp_write(dest_unit,axes_dst(i))
    end do

    call mpp_write(dest_unit,depth_axis)

    ! write variable data

    call mpp_write(dest_unit, field_lon_dst,lon_dst)
    call mpp_write(dest_unit, field_lat_dst,lat_dst)

  end subroutine setup_meta


  !#####################################################################
  !--- remap data from src grid to destination grid.
  subroutine process_data

    integer                             :: ndivs, i, j, k, n, nt, stat, ntimelevel, l
    integer                             :: isc, iec, jsc, jec, layout(2) = (/1,0/)
    integer                             :: kstart, kend
    real                                :: tmp_x, tmp_y
    type(domain2d)                      :: Domain
    type(horiz_interp_type)             :: Interp
    real, dimension(:,:,:), allocatable :: tmp1, tmp2, tmp, mask_src, data_dst, data_src
    real, dimension(:,:,:), allocatable :: depth_src_3d, depth_dst_3d

    !--- decompose model grid points
    !--- mapping can get expensive so we distribute the task at this level
    ndivs = mpp_npes()
    call mpp_define_layout ((/1,ni_dst,1,nj_dst/), ndivs, layout)
    call mpp_define_domains((/1,ni_dst,1,nj_dst/),layout, Domain,xhalo=0,yhalo=0)  
    call mpp_get_compute_domain (Domain, isc, iec, jsc, jec)

    call horiz_interp_new(Interp, lon_src*D2R, lat_src*D2R, lon_dst(isc:iec,jsc:jec)*D2R, &
         lat_dst(isc:iec,jsc:jec)*D2R, interp_method = trim(interp_method),                &
         num_nbrs = num_nbrs, max_dist=max_dist, grid_at_center = .true. )

    allocate(tmp1(ni_src,nj_src, nk_src))
    allocate(data_src(ni_src,nj_src, nk_src))
    if( use_source_vertical_grid) then 
       allocate(data_dst(ni_dst,nj_dst,nk_src), tmp(isc:iec,jsc:jec, nk_src) )
       call mpp_domains_set_stack_size(2*ni_dst*nj_dst*nk_src)
    else
       call mpp_domains_set_stack_size(2*ni_dst*nj_dst*nk_dst)
       allocate(tmp2(isc:iec,jsc:jec, nk_src) )
       allocate(data_dst(ni_dst,nj_dst,nk_dst), tmp(isc:iec,jsc:jec, nk_dst) )
       allocate(depth_src_3d(isc:iec,jsc:jec, nk_src))
       allocate(depth_dst_3d(isc:iec,jsc:jec,nk_dst))
       do k=1,nk_src
          depth_src_3d(:,:,k) = depth_src(k)
       enddo
       do k=1,nk_dst
          depth_dst_3d(:,:,k) = depth_dst(k)
       enddo
    endif


    if(ntimes_saved > 0) then
       ntimelevel = ntimes_saved       
    else
       ntimelevel = ntime_src
    endif

    write(stdout(),*)' There are ', ntimelevel, ' time steps to be saved to output file.'

    ! for vertical interpolation, Set value of levels shallower than the shallowest source level to 
    ! be the source value at shallowest level. 
    ! Set value of levels deeper than the deepest source level to 
    ! be the source value at deepest level. 

    if( .not. use_source_vertical_grid) then 
       do kstart = 1, nk_dst
          if( depth_dst(kstart) .GE. depth_src(1) ) exit 
       enddo
       do kend = nk_dst, 1, -1
          if( depth_dst(kend) .LE. depth_src(nk_src) ) exit 
       enddo
       if( kstart > 1 ) then
          write(stdout(),*)"NOTE from regrid_3d: the value from level 1 to level ", kstart-1, &
                 " will be set to the value at the shallowest source levle."
       endif
       if( kend < nk_dst ) then
          write(stdout(),*)"NOTE from regrid_3d: the value from level ", kend+1, " to level ", nk_dst, &
                 " will be set to the value at the deepest source levle."
       endif
    endif

    do l = 1, ntimelevel 
       if(ntimes_saved > 0) then
          nt = timelevel_saved(l)
       else
          nt = l
       endif
       write(stdout(),*)'**************At time step ', l       
       do n = 1, numfields
          !--- read source data
          call mpp_read(src_unit,src_field(n),data_src,nt)
          if( n==1 .and. nt == 1) then
             allocate(mask_src(ni_src,nj_src,nk_src))
             mask_src = 1.0
             do k = 1, nk_src
                do j = 1, nj_src
                   do i = 1, ni_src
                      if (abs(data_src(i,j,k) - missing(1)) <= tol) mask_src(i,j,k) = 0.
                      if (abs(data_src(i,j,k)) > max_val) mask_src(i,j,k) = 0.0
                   enddo
                enddo
             enddo
          endif
          !--- if scale_factor is not 1, multiple scale_factor to each data
          if(scale_factor(n) .ne. 1) then
             do k =1, nk_src
                do j =1, nj_src
                   do i =1, ni_src
                      if(mask_src(i,j,k) >0.5) data_src(i,j,k) = data_src(i,j,k)*scale_factor(n)
                   enddo
                enddo
             enddo
          endif

          !--- begin to regrid data
          if( use_source_vertical_grid) then 
             if(apply_mask) then
                do k = 1, nk_src
                   call horiz_interp(Interp, data_src(:,:,k), tmp(:,:,k),  &
                        mask_in = mask_src(:,:,k), missing_value=missing(n) )
                enddo
             else
                if(any(mask_src == 0.0) ) then ! do laplace extrap if needed.
                   call extrap(data_src(:,:,:), tmp1(:,:,:), stop_crit(n), &
                      missing(n)*scale_factor(n), is_cyclic ) 
                   do k = 1, nk_src 
                      call horiz_interp(Interp, tmp1(:,:,k), tmp(:,:,k))
                   enddo
                else
                   do k = 1, nk_src
                      call horiz_interp(Interp, data_src(:,:,k), tmp(:,:,k) )
                   enddo
                endif
             endif
          else
             if(any(mask_src == 0.0) ) then ! do laplace extrap if needed.
                call extrap(data_src(:,:,:), tmp1(:,:,:), stop_crit(n), &
                      missing(n)*scale_factor(n), is_cyclic ) 
                do k = 1, nk_src
                   call horiz_interp(Interp, tmp1(:,:,k), tmp2(:,:,k) )
                enddo
             else
                do k = 1, nk_src
                   call horiz_interp(Interp, data_src(:,:,k), tmp2(:,:,k) )
                enddo
             endif
             
             do k = 1, kstart-1
                tmp(:,:,k) = tmp2(:,:,1)
             enddo
                
             do k = kend+1, nk_dst
                tmp(:,:,k) = tmp2(:,:,nk_dst)
             enddo

             call interp_1d(depth_src_3d,depth_dst_3d(:,:,kstart:kend),tmp2,tmp(:,:,kstart:kend))

             do k = 1, nk_dst
                if(apply_mask) then
                   do j = jsc, jec
                      do i = isc, iec
                         if(mask_dst(i,j,k) < 0.5) tmp(i,j,k) = missing(n)
                      enddo
                   enddo
                endif
             enddo
          endif

          !--- get global data
          call mpp_global_field(Domain,tmp, data_dst(:,:,:) )
          if(mpp_pe()==mpp_root_pe()) then  ! --- write out data from root pe
             if (time_axis_exists) then 
                call mpp_write(dest_unit, dest_field(n), data_dst,time_in(nt))
             else
                call mpp_write(dest_unit, dest_field(n), data_dst)  
             endif

             if(debug) then            !--- the chksum is for debugging purpose

                write(stdout(),*)'NOTE: At time level ',nt,', Chksum for ',trim(src_field_name(n)), &
                     ' after-regrid data ', mpp_chksum(data_dst, (/mpp_root_pe()/) )
             endif
          endif
       enddo
    enddo

    if(allocated(tmp)) deallocate(tmp)
    if(allocated(tmp1)) deallocate(tmp1)
    if(allocated(data_src)) deallocate(data_src)
    if(allocated(data_dst)) deallocate(data_dst)
    if(allocated(mask_src)) deallocate(mask_src)
    if(allocated(tmp2)) deallocate(tmp2)
    if(allocated(depth_src_3d)) deallocate(depth_src_3d)
    if(allocated(depth_dst_3d)) deallocate(depth_dst_3d)

    !--- write out chksum for parallel checking

    call mpp_close(dst_grid_unit)
    call mpp_close(dest_unit)
    call mpp_close(src_unit)
    call horiz_interp_end

  end subroutine process_data

  !#####################################################################
  !--- release the memory
  subroutine regrid_3d_end

  deallocate(depth_src, depth_dst, lon_src, lat_src, lon_dst, lat_dst, mask_dst)
  if(time_axis_exists) deallocate(time_in)

  end subroutine regrid_3d_end


  !#####################################################################
  subroutine extrap(data_in, data_out, crit, missing_value, is_cyclic )
    real, dimension(:,:,:),  intent(in) :: data_in
    real, dimension(:,:,:), intent(out) :: data_out
    real,                    intent(in) :: crit, missing_value 
    logical,                 intent(in) :: is_cyclic 
    real                                :: resmax, initial_guess = 0.0
    integer                             :: ni, nj, nk, i, j, k, n

    real, dimension(0:size(data_in,1)+1, 0:size(data_in,2)+1) :: tmp
    real, dimension(size(data_in,1), size(data_in,2) )        :: sor, res

    real, dimension(size(data_in,1), size(data_in,2) )        :: cfn, cfs, cfe, cfw, cfc
    real                                                      :: latp, latm, cstr, csm, csj
    real, dimension(size(data_in,1))                          :: dxu, dxt
    real, dimension(size(data_in,2))                          :: dyu, dyt

    ni = size(data_in,1)
    nj = size(data_in,2)
    nk = size(data_in,3)

    ! construct grid factors for a sphere
    do j=1,nj-1
      dyu(j) = lat_src(j+1)-lat_src(j)
    enddo
    dyu(nj) = dyu(nj-1)
    do j=2,nj
      dyt(j) = 0.5*(dyu(j)+dyu(j-1))
    enddo
    dyt(1) = dyt(2)
    do i=1,ni-1
      dxu(i) = lon_src(i+1)-lon_src(i)
    enddo
    dxu(ni) = dxu(ni-1)
    do i=2,ni
      dxt(i) = 0.5*(dxu(i)+dxu(i-1))
    enddo
    dxt(1) = dxt(2)
    do j= 1, nj
      if (j == nj) then
        latp = lat_src(j) + 0.5*(lat_src(j) - lat_src(j-1))
      else
        latp = 0.5*(lat_src(j) + lat_src(j+1))
      endif
      if (j == 1) then
        latm = lat_src(j) - 0.5*(lat_src(j+1) - lat_src(j))
      else
        latm = 0.5*(lat_src(j) + lat_src(j-1))
      endif
      csj  = cos(latp*pi/180.0)
      csm  = cos(latm*pi/180.0)
      cstr = 1.0/cos(lat_src(j)*pi/180.)       
      do i= 1,ni
        cfn(i,j) = csj*cstr/(dyt(j)*dyu(j))
        cfs(i,j) = csm*cstr/(dyt(j)*dyu(max(j-1,1)))
        cfe(i,j) = cstr**2/(dxu(i)*dxt(i))
        cfw(i,j) = cstr**2/(dxu(max(i-1,1))*dxt(i))
        cfc(i,j) = 1.0/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
        cfn(i,j) = cfn(i,j)*cfc(i,j)
        cfs(i,j) = cfs(i,j)*cfc(i,j)
        cfe(i,j) = cfe(i,j)*cfc(i,j)
        cfw(i,j) = cfw(i,j)*cfc(i,j)         
      enddo
    enddo

    tmp = 0.0
    do j= 1, nj
       do i= 1, ni
          if(abs(data_in(i,j,1) - missing_value) <= tol) then
             tmp(i,j) = initial_guess
          endif
       enddo
    enddo

    do k = 1, nk
       do j= 1, nj
          do i= 1, ni
             if(abs(data_in(i,j,k) - missing_value) <= tol ) then
                sor(i,j) = rel_coef
             else
                tmp(i,j) = data_in(i,j,k)
                sor(i,j) = 0.0
             endif
          enddo
       enddo

       call fill_boundaries(tmp, is_cyclic)

       ! iterate
       n=1
       do            
          resmax=0.0
          do j= 1, nj
             do i= 1, ni
                res(i,j) = cfw(i,j)*tmp(i-1,j) + cfe(i,j)*tmp(i+1,j) + cfs(i,j)*tmp(i,j-1) + cfn(i,j)*tmp(i,j+1) - tmp(i,j)
             enddo
          enddo

          do j= 1, nj
             do i= 1, ni
                res(i,j) = res(i,j)*sor(i,j)
                tmp(i,j)=tmp(i,j)+res(i,j)
                resmax = max(abs(res(i,j)),resmax)
             enddo
          enddo

          if(resmax .le. crit .or. n > max_iter) then
             data_out(:,:,k) = tmp(1:ni,1:nj)
             exit
          endif

          !--- update boundaries

          call fill_boundaries(tmp, is_cyclic)

          n=n+1

       enddo

       if(mpp_pe() == mpp_root_pe() ) write(stdout(),'(a,i6,a)') 'Stopped after ',n,' iterations'
       if(mpp_pe() == mpp_root_pe() ) write(stdout(),'(a,f10.4)') 'maxres= ',resmax
    enddo

  end subroutine extrap

  !#####################################################################


  subroutine fill_boundaries(data, is_cyclic)

    real, dimension(0:,0:), intent(inout) :: data
    logical,                   intent(in) :: is_cyclic
    integer :: i,j, ni, nj

    ni = size(data,1) - 2
    nj = size(data,2) - 2

    if(is_cyclic) then
       data(0,1:nj) = data(ni,1:nj)
       data(ni+1,1:nj) = data(1,1:nj)
    endif

    return

  end subroutine fill_boundaries

  !#####################################################################

end program regrid_3d








program regrid_2d
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">M.J. Harrison </CONTACT>
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Zhi Liang</CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">Bonnie Samuels</REVIEWER>
  !
  !<DESCRIPTION>
  ! regrid 2-d lat-lon gridded data to logically rectangular grid
  ! described by grid descriptor file.  If two fields are specified 
  ! for regridding, it is assumed that they represent vector components.
  ! Rotation to the local grid direction on the target grid will be performed.
  !
  !</DESCRIPTION>
  !

  use mpp_mod,          only : mpp_error, mpp_pe,  mpp_npes, mpp_root_pe, mpp_chksum
  use mpp_mod,          only : FATAL, WARNING, stdout, stdlog
  use mpp_io_mod,       only : mpp_open, mpp_close, mpp_read, mpp_write, mpp_write_meta
  use mpp_io_mod,       only : mpp_copy_meta, axistype, fieldtype, atttype
  use mpp_io_mod,       only : mpp_get_atts, mpp_get_info, mpp_get_fields, mpp_get_times
  use mpp_io_mod,       only : mpp_get_axes, mpp_get_axis_data
  use mpp_io_mod,       only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
  use mpp_io_mod,       only : mpp_get_att_name, mpp_get_att_char, mpp_get_att_type, mpp_get_att_real
  use mpp_domains_mod,  only : mpp_update_domains, mpp_define_domains, mpp_global_field
  use mpp_domains_mod,  only : domain2d, mpp_define_layout, mpp_get_compute_domain
  use mpp_domains_mod,  only : mpp_domains_set_stack_size
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_end, horiz_interp_type
  use axis_utils_mod,   only : get_axis_cart
  use fms_mod,          only : fms_init, fms_end, open_namelist_file, close_file, file_exist
  use fms_mod,          only : check_nml_error, write_version_number, lowercase
  use constants_mod,    only : constants_init, PI

  implicit none

  !--- namelist interface
  !<NAMELIST NAME="regrid_2d_nml">
  ! <DATA NAME="src_file" TYPE="character(len=128)" DEFAULT="src_file.nc">
  !  Name of input file containing grid and data to be regridded.
  ! </DATA>
  ! <DATA NAME="src_field_name" TYPE="character(len=128),dimension(2)" >
  !  Name of input field(s).
  ! </DATA>
  ! <DATA NAME="dest_field_name" TYPE="character(len=128),dimension(2)" >
  !  Name of output field(s). If it is not specified in the namelist, it will 
  !  get the value from src_field_name
  ! </DATA>
  ! <DATA NAME="dest_grid" TYPE="character(len=128)" DEFAULT="dest_grid.nc">
  !  Name of grid descriptor file containing target grid information.
  ! </DATA>
  ! <DATA NAME="dest_file"  TYPE="character(len=128)" DEFAULT="dest_file.nc">
  !  Name of output file.
  ! </DATA>
  ! <DATA NAME="numfields"  TYPE="integer" DEFAULT="1">
  !  Number of fields (1 or 2). If numfields=2, then the fields are assumed
  !  to represent vector components.
  ! </DATA>
  ! <DATA NAME="dest_grid_type"  TYPE="character(len=1)" DEFAULT="T">
  !  Name of target grid type.  Valid choices are (T)racer and (C)orner
  ! </DATA>
  ! <DATA NAME="stop_crit"  TYPE="character(len=1),dimension(2)" DEFAULT="0.001">
  !  The stopping criteria when extrapping data onto missing points.
  ! </DATA>
  ! <DATA NAME="apply_dest_mask"  TYPE="logical" DEFAULT="true">
  !  flag to indicate if the land/sea mask of destination grid will be applied 
  !  on the output dest_file. when true, land/sea mask of destination grid will 
  !  be applied on the output dest_file. When false, the output data will be 
  !  global data. Default is true. 
  ! </DATA>
  ! <DATA NAME="vector_field"  TYPE="logical" DEFAULT="False">
  !  True if fields are vector components.
  ! </DATA>
  ! <DATA NAME="level"  TYPE="integer" DEFAULT="1">
  !  Vertical level from the source grid if one exists.  
  ! </DATA>
  ! <DATA NAME="num_nbrs"  TYPE="integer" DEFAULT="10">
  ! Number of nearest neighbors for regridding  
  ! </DATA>
  ! <DATA NAME="max_dist"  TYPE="integer" DEFAULT="0.17" UNITS="radians">
  !  Maximum radial influence for regridding.
  ! </DATA>
  ! <DATA NAME="scale_factor" TYPE="real,dimension(2) ">
  ! scaling factor for data (e.g. -1 to flip sign or 0.01 to convert from centimeters)
  ! </DATA>
  !<DATA NAME="interp_method"  TYPE= "character(len=20)" >
  !  specifying the remapping method when remampping data onto current grid.
  !  Its value can be "spherical" or " bilinear". "spherical" interpolation is a 
  !  inverse distance weighted interpolation algorithm. Default value is "bilinear". 
  !  "bilinear" interpolation is recommanded, since bilinear interpolation will provide 
  !  more smooth results than "spherical" interpolation (especially when interpolating 
  !  from coarse grid to fine grid). Plus bilinear interpolation is much more efficiency 
  !  than "spherical interpolation". 
  ! </DATA>
  !</NAMELIST>
  character(len=128) :: src_file = 'src_file.nc'
  character(len=128) :: dest_grid = 'dest_grid.nc'
  character(len=128) :: dest_file = 'dest_file.nc'
  character(len=128) :: src_field_name(2) 
  character(len=128) :: dest_field_name(2) = ''
  real               :: stop_crit(2)       = 0.005 
  character(len=1)   :: dest_grid_type = 'T'
  integer            :: numfields = 1
  integer            :: level = 1
  integer            :: num_nbrs = 5
  real               :: max_dist = 0.1
  real               :: scale_factor(2) = 1.0
  logical            :: vector_field=.false.
  logical            :: apply_dest_mask = .true.
  character(len=20)  :: interp_method = "bilinear"

  namelist /regrid_2d_nml/ src_file, src_field_name, dest_field_name, dest_file, dest_grid, &
       dest_grid_type, numfields, scale_factor, vector_field,level, num_nbrs, max_dist,     &
       apply_dest_mask, stop_crit, interp_method

  !---------------------------------------------------------------------
  integer            :: src_unit, dst_grid_unit
  integer            :: ni_src, nj_src, ni_dst, nj_dst, ntime_src
  type(axistype)     :: time_axis, axes_dst(2)
  type(fieldtype)    :: field_lon_dst, field_lat_dst, src_field(2)
  character(len=16)  :: y_boundary_type
  character(len=128) :: title
  logical            :: time_axis_exists = .false.
  real, parameter    :: tol = 1.e-10   ! tolerance for detecting missing values
  real, parameter    :: max_val=1.e20
  real, parameter    :: rel_coef = 0.9
  integer, parameter :: max_iter = 2000
  logical            :: is_cyclic = .true.  ! we suppose the source data is always global data.
  real               :: missing(2) = -1e20
  real               :: D2R
  logical            :: do_extrap
  type(atttype), dimension(:), allocatable :: global_atts
  real, dimension(:),          allocatable :: time_in
  real, dimension(:),          allocatable :: lon_src, lat_src
  real, dimension(:,:),        allocatable :: lon_dst, lat_dst
  real, dimension(:,:),        allocatable :: sin_rot, cos_rot
  real, dimension(:,:),        allocatable :: mask_dst
  real, dimension(:,:,:,:),    allocatable :: data_dst, data_src

  !--- version information variables
  character(len=128) :: version='CVS $Id: regrid_2d.f90,v 20.0 2013/12/14 00:31:09 fms Exp $'
  character(len=128) :: tagname='Tag $Name: tikal $'

  ! --- Begin of the program

  ! --- call fms_init, which will call mpp_init, mpp_io_init, mpp_domains_init
  call fms_init
  call constants_init

  !--- call regrid_2d initialization routine
  call regrid_2d_init

  !--- read the dest_grid file
  call read_dst_grid

  !--- read src_file
  call read_src_file

  !--- remap data from src grid to dest grid
  call remap_data

  !--- write out  metadata and data to output file dest_file
  if(mpp_pe() == mpp_root_pe() ) call write_dst_file

  call regrid_2d_end

  call fms_end

contains

  !#####################################################################
  ! --- read the namelist and write the version and namelist to logfile. Also
  ! --- write the namelist to standard output
  subroutine regrid_2d_init

    integer :: io_status, unit, ierr, n

    D2R = PI/180.0
    ! --- read namelist ------------------------------------------------

    if(file_exist('input.nml')) then
       unit = open_namelist_file()
       read (unit,regrid_2d_nml,IOSTAT=io_status)
       write (stdout(),'(/)')
       write (stdout(),regrid_2d_nml)  
       ierr = check_nml_error(io_status, 'regrid_2d_nml')
       call close_file(unit)
    else
       call mpp_error(FATAL, 'regrid_2d: file input.nml does not exist' )
    endif

    if (numfields.ne.2.and.vector_field) &
         call mpp_error(FATAL,'regrid_2d: 2 components of vector field not specified')
    if (numfields == 2.and..not.vector_field) &
         call mpp_error(FATAL,'regrid_2d: only 1 scalar at a time')
    if (numfields .gt. 2) call mpp_error(FATAL,'regrid_2d: more than 2 fields specified')
    if (numfields .le. 0) call mpp_error(FATAL,'regrid_2d: No field specified')
    if (dest_grid_type .ne. 'T' .and. dest_grid_type .ne. 'C') &
         call mpp_error(FATAL,'regrid_2d: nml dest_grid_type = '//dest_grid_type//', should be either T or C')

    !--- if dest_field_name is not defined in the namelist, get it from src_field_name
    do n = 1, numfields
       if(trim(dest_field_name(n)) == '') dest_field_name(n) = trim(src_field_name(n))
    enddo 
 
    !--- write version information
    call write_version_number(version, tagname)

  end subroutine regrid_2d_init

  !#####################################################################
  !--- open grid file and store grid info
  subroutine read_dst_grid

    integer                                    :: ndim, nvar, natt, ntime, i
    integer                                    :: len1, siz_in(3)
    logical                                    :: found_xt, found_yt, found_wet
    logical                                    :: found_xc, found_yc, found_angle
    character(len=32)                          :: name
    real, dimension(:,:),          allocatable :: angle
    type(axistype), allocatable, dimension(:)  :: axes
    type(fieldtype), allocatable, dimension(:) :: fields

    if(.not. file_exist(trim(dest_grid)) ) &
         call mpp_error(FATAL, 'regrid_2d: file '//trim(dest_grid)//' does not exist')

    call mpp_open(dst_grid_unit, trim(dest_grid),&
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    call mpp_get_info(dst_grid_unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar), global_atts(natt), axes(ndim) )
    call mpp_get_axes(dst_grid_unit, axes)
    call mpp_get_atts(dst_grid_unit, global_atts)
    call mpp_get_fields(dst_grid_unit,fields)

    do i=1,natt
       if (trim(mpp_get_att_name(global_atts(i))) == 'y_boundary_type') &
               y_boundary_type = trim(mpp_get_att_char(global_atts(i)))
       if (trim(mpp_get_att_name(global_atts(i))) == 'filename') &
               title = trim(mpp_get_att_char(global_atts(i)))
    enddo

    write(stdout(),*) 'Retrieving input grid information ...'
    write(stdout(),*) 'Input grid name: ',trim(title)
    write(stdout(),*) 'y_boundary_type: ',trim(y_boundary_type)

    !--------------------------------------------------------------------
    ! get output grid information
    !--------------------------------------------------------------------
    ni_dst=0; nj_dst=0
    do i=1,ndim
       call mpp_get_atts(axes(i),name=name,len=len1)
       select case (trim(name))
       case ('grid_x_T')
          ni_dst = len1
       case ('grid_y_T')
          nj_dst = len1
       end select
    enddo
    if(ni_dst==0) call mpp_error(FATAL,'regrid_2d: file '//trim(dest_grid)//' does not contain axis grid_x_T')
    if(nj_dst==0) call mpp_error(FATAL,'regrid_2d: file '//trim(dest_grid)//' does not contain axis grid_y_T')

    allocate(lon_dst(ni_dst,nj_dst), lat_dst(ni_dst,nj_dst) )
    allocate(mask_dst(ni_dst,nj_dst) )

    found_xt = .FALSE.;  found_yt = .FALSE.; found_wet   = .FALSE.
    found_xc = .FALSE.;  found_yc = .FALSE.; found_angle = .FALSE.
    do i=1,nvar
       call mpp_get_atts(fields(i),name=name,ndim=ndim)
       !--- get the land/sea mask
       if(trim(name) == 'wet') then
          found_wet = .true.
          call mpp_read(dst_grid_unit,fields(i),mask_dst)
       endif
       if(dest_grid_type == 'T') then
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
       end select
       else if ( dest_grid_type == 'C' ) then
          select case (trim(name))
          case ('x_C')
             found_xc = .true.
             call mpp_read(dst_grid_unit,fields(i),lon_dst)
             field_lon_dst = fields(i)
             call mpp_get_atts(fields(i),axes=axes_dst)
          case('y_C')
             found_yc = .true.
             call mpp_read(dst_grid_unit,fields(i),lat_dst)
             field_lat_dst = fields(i)
          case('angle_C')
             found_angle = .true.
             allocate(angle(ni_dst,nj_dst), sin_rot(ni_dst,nj_dst), cos_rot(ni_dst,nj_dst) )
             call mpp_read(dst_grid_unit,fields(i),angle)
             sin_rot = sin(angle*D2R)
             cos_rot = cos(angle*D2R)
             deallocate(angle)
          end select
       endif
    enddo
    if(dest_grid_type == 'T') then
       if(.not.found_xt) call mpp_error(FATAL,'regrid_2d: field x_T is not in the file '//trim(dest_grid) )
       if(.not.found_yt) call mpp_error(FATAL,'regrid_2d: field y_T is not in the file '//trim(dest_grid) )    
             else
       if(.not.found_xc) call mpp_error(FATAL,'regrid_2d: field x_C is not in the file '//trim(dest_grid) )
       if(.not.found_yc) call mpp_error(FATAL,'regrid_2d: field y_C is not in the file '//trim(dest_grid) )    
       if(.not.found_angle) call mpp_error(FATAL,'regrid_2d: field angle_C is not in the file '//trim(dest_grid) )
             endif
    if(.not.found_wet) call mpp_error(FATAL,'regrid_2d: field wet is not in the file '//trim(dest_grid) )
    
    if(.not. apply_dest_mask) mask_dst = 1.0  ! will get global data

    deallocate(fields, axes)    

  end subroutine read_dst_grid

  !#####################################################################
  !--- read source grid and source data from src_file
  subroutine read_src_file

    integer                                    :: ndim, nvar, natt, n
    integer                                    :: nt, i, j, k, jj, len1, nk_src
    logical                                    :: flip_y, found_src_field(2)
    character(len=1)                           :: cart
    character(len=32)                          :: name
    real, dimension(:),           allocatable  :: tmp1d 
    real, dimension(:,:),         allocatable  :: tmp 
    real, dimension(:,:,:),       allocatable  :: tmp3d 
    type(axistype), allocatable, dimension(:)  :: axes
    type(fieldtype), allocatable, dimension(:) :: fields


    if(.not. file_exist(trim(src_file)) ) &
         call mpp_error(FATAL, 'regrid_2d: file '//trim(src_file)//' does not exist')

    call mpp_open(src_unit, trim(src_file),&
         action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(src_unit, ndim, nvar, natt, ntime_src)

    allocate(fields(nvar))
    call mpp_get_fields(src_unit, fields)

    if (numfields > nvar) call mpp_error(FATAL,'not enough fields in file')
    found_src_field = .FALSE.
    do n=1,numfields
       do i=1,nvar
          call mpp_get_atts(fields(i),name=name)
          if (lowercase(trim(src_field_name(n))) == lowercase(trim(name))) then
             found_src_field(n) = .true.
             src_field(n) = fields(i)
             write(stdout(),*) 'Interpolating src field : ',trim(name), ' to ocean model grid'
          endif
       end do
    enddo

    do n=1,numfields 
       if(.not. found_src_field(n)) call mpp_error(FATAL, 'regrid_2d: field '&
            //trim(src_field_name(n))//' is not in the file '//trim(src_file) )
    enddo

    !--- get the src grid
    call mpp_get_atts(src_field(1),ndim=ndim)
    allocate(axes(ndim))
    call mpp_get_atts(src_field(1),axes=axes)
    flip_y = .false.
    ni_src=0; nj_src=0; nk_src=1; ntime_src = 1
    do j=1,ndim
       call mpp_get_atts(axes(j),len=len1)
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
          if (lat_src(2) < lat_src(1)) then
             flip_y = .true.
             allocate(tmp1d(nj_src))
             tmp1d=lat_src
             do jj=1,nj_src
                lat_src(jj)=tmp1d(nj_src-jj+1)
             enddo
             deallocate(tmp1d)
          endif
       case('Z')
          nk_src = len1
       case ('T')
          ntime_src = len1
          time_axis_exists = .true.
          allocate(time_in(ntime_src))
          call mpp_get_times(src_unit, time_in)
          time_axis = axes(j)
       end select
    enddo

    if(ni_src==0) call mpp_error(FATAL,'regrid_2d: file ' &
         //trim(src_file)//' does not contain axis with cartesian attributes = "X" ')
    if(nj_src==0) call mpp_error(FATAL,'regrid_2d: file '&
         //trim(src_file)//' does not contain axis with cartesian attributes = "Y" ')

    allocate(data_src(ni_src,nj_src,2,ntime_src) )
    allocate(tmp3d(ni_src,nj_src,nk_src), tmp(ni_src,nj_src))

    !--------------------------------------------------------------------
    ! loop through source fields 
    !--------------------------------------------------------------------
    do n=1, numfields
       !--- get the missing value
       call mpp_get_atts(src_field(n),missing=missing(n) )
       if (level > nk_src) call mpp_error(FATAL,'selected level exceeds size of input array')
       if (nk_src > 1) write(*,*) 'warning: selecting level ',level,' from 3d array'
       do nt=1,ntime_src
          call mpp_read(src_unit,src_field(n), tmp3d,nt)

          tmp=tmp3d(:,:,level)
          !--- set up the mask of source data
          if(n==1 .and. nt == 1) then
             do_extrap = .false.
             do j = 1, nj_src
                do i = 1, ni_src
                   if( abs(tmp(i,j) - missing(n)) <= tol ) then
                      do_extrap = .true.
                   endif
                enddo
             enddo
          endif
           
          if (flip_y) then
             do j =1, nj_src
                data_src(:,j,n,nt)=tmp(:,nj_src-j+1)* scale_factor(n)
             enddo
          else
             data_src(:,:,n,nt) = tmp(:,:)* scale_factor(n)
          endif
       enddo

    enddo

    deallocate(fields, axes, tmp, tmp3d)   


  end subroutine read_src_file

  !#####################################################################
  !--- remap data from src grid to destination grid.
  subroutine remap_data

    integer                           :: ndivs, i, j, n, nt, l
    integer                             :: isc, iec, jsc, jec, layout(2) = (/1,0/)
    real                                :: tmp_x, tmp_y
    type(domain2d)                      :: Domain
    type(horiz_interp_type)             :: Interp
    real, dimension(:,:), allocatable :: tmp_src, tmp_dst

    !--- decompose model grid points
    !--- mapping can get expensive so we distribute the task at this level
    ndivs = mpp_npes()
    call mpp_define_layout ((/1,ni_dst,1,nj_dst/), ndivs, layout)
    call mpp_define_domains((/1,ni_dst,1,nj_dst/),layout, Domain,xhalo=0,yhalo=0)  
    call mpp_get_compute_domain (Domain, isc, iec, jsc, jec)

    call horiz_interp_new(Interp, lon_src*D2R, lat_src*D2R, lon_dst(isc:iec,jsc:jec)*D2R, &
         lat_dst(isc:iec,jsc:jec)*D2R, interp_method = trim(interp_method),                &
         num_nbrs = num_nbrs, max_dist=max_dist, grid_at_center = .true. )

    allocate(data_dst(ni_dst,nj_dst,numfields,ntime_src))
    allocate(tmp_src(ni_src, nj_src), tmp_dst(isc:iec,jsc:jec) )

    call mpp_domains_set_stack_size(2*ni_dst*nj_dst)
    do n=1, numfields
       do l=1,ntime_src
          if(do_extrap) then
             call extrap(data_src(:,:,n,l), tmp_src(:,:), stop_crit(n), missing(n)*scale_factor(n), is_cyclic ) 
             call horiz_interp(Interp, tmp_src(:,:), tmp_dst )
          else
             call horiz_interp(Interp, data_src(:,:,n,l), tmp_dst)
          endif
          !--- apply the mask to destination data
          do j = jsc, jec
             do i = isc, iec
                if(mask_dst(i,j) < 0.5) tmp_dst(i,j) = missing(n)
             enddo
          enddo
          call mpp_global_field(Domain,tmp_dst, data_dst(:,:,n,l) )
       enddo
    enddo

    deallocate(tmp_src, tmp_dst)

    call horiz_interp_end

    if (vector_field) then
          do nt=1,ntime_src
          do j=1,nj_dst
             do i=1,ni_dst
                if(mask_dst(i,j) >= 0.5) then 
                   tmp_x = cos_rot(i,j)*data_dst(i,j,1,nt)+sin_rot(i,j)*data_dst(i,j,2,nt)
                   tmp_y = -1.*sin_rot(i,j)*data_dst(i,j,1,nt)+cos_rot(i,j)*data_dst(i,j,2,nt)
                   data_dst(i,j,1,nt) = tmp_x
                   data_dst(i,j,2,nt) = tmp_y
                endif
             enddo
          enddo
          if (y_boundary_type == 'fold_north_edge' .and. dest_grid == 'c') then
             do i=1,ni_dst/2
                if(mask_dst(i,nj_dst) >= 0.5) then 
                   data_dst(i,nj_dst,1,nt) = -1.*data_dst(ni_dst-i,nj_dst,1,nt)
                   data_dst(i,nj_dst,2,nt) = -1.*data_dst(ni_dst-i,nj_dst,2,nt)
                endif
             enddo
          endif
       enddo
    endif

    !--- apply the destination mask to the 

    !--- write out chksum for parallel checking
    if(mpp_pe()==mpp_root_pe()) then
       write(stdout(),*)'NOTE: Chksum for after-regrid data ', mpp_chksum(data_dst, (/mpp_root_pe()/) )
    endif

  end subroutine remap_data

  !#####################################################################
  !--- write remapped data to file dest_file
  subroutine write_dst_file

    type(fieldtype)    :: dest_field(2)
    integer            :: unit, i, j, n
    character(len=32)  :: units
    character(len=128) :: longname
    !--------------------------------------------------------------------
    ! write output file metadata
    !--------------------------------------------------------------------

    call mpp_open(unit, trim(dest_file),MPP_OVERWR,MPP_NETCDF,threading=MPP_SINGLE,&
         fileset=MPP_SINGLE)
    !
    ! write global atts, will be changed in the future. The reason is we 
    ! need to rename regrid_2d.f90 to regrid_2d.F90
    !
!!$    do i=1,size(global_atts(:))
!!$       if (global_atts(i)%name == 'title') then
!!$          global_atts(i)%catt = 'SBC field interpolated to ocean grid'
!!$          global_atts(i)%len = 36
!!$       endif
!!$    end do

    !
    ! write axis metadata
    !
    do i=1,size(axes_dst(:))
       call mpp_copy_meta(unit,axes_dst(i))
    end do

    if (time_axis_exists) then
       call mpp_copy_meta(unit, time_axis)
    endif

    !
    ! write variable metadata
    !
    call mpp_copy_meta(unit,field_lon_dst,axes=axes_dst)
    call mpp_copy_meta(unit,field_lat_dst,axes=axes_dst)

    do n = 1, numfields
       call mpp_get_atts(src_field(n),units=units,longname=longname)
       if (time_axis_exists) then 
          call mpp_write_meta(unit, dest_field(n), (/axes_dst(1),axes_dst(2), time_axis/), &
               trim(dest_field_name(n)), units,longname, missing=missing(n), pack=1)
       else
          call mpp_write_meta(unit, dest_field(n), (/axes_dst(1),axes_dst(2)/), &
               trim(dest_field_name(n)), units,longname, missing=missing(n), pack=1)
       endif
    enddo

    ! write axis data

    do i=1,size(axes_dst(:))
       call mpp_write(unit,axes_dst(i))
    end do

    ! write variable data

    call mpp_write(unit, field_lon_dst,lon_dst)
    call mpp_write(unit, field_lat_dst,lat_dst)

    do j=1, ntime_src
       do n=1,numfields
    if (time_axis_exists) then 
             call mpp_write(unit, dest_field(n), data_dst(:,:,n,j),time_in(j))
          else
             call mpp_write(unit, dest_field(n), data_dst(:,:,n,1)) 
          endif
          enddo
       enddo

    call mpp_close(dst_grid_unit)
    call mpp_close(src_unit)
    call mpp_close(unit)

  end subroutine write_dst_file

  !#####################################################################
  !--- release the memory
  subroutine regrid_2d_end

    deallocate( lon_src, lat_src, lon_dst, lat_dst)
  deallocate(data_dst, data_src, mask_dst)
  if(time_axis_exists) deallocate(time_in)
    if(dest_grid_type == 'C') deallocate(sin_rot, cos_rot)
  end subroutine regrid_2d_end

  !#####################################################################
  subroutine extrap(data_in, data_out, crit, missing_value, is_cyclic )
    real, dimension(:,:),    intent(in) :: data_in
    real, dimension(:,:),   intent(out) :: data_out
    real,                    intent(in) :: crit, missing_value 
    logical,                 intent(in) :: is_cyclic 
    real                                :: resmax, initial_guess = 0.0
    integer                             :: ni, nj, i, j, n

    real, dimension(0:size(data_in,1)+1, 0:size(data_in,2)+1) :: tmp
    real, dimension(size(data_in,1), size(data_in,2) )        :: sor, res

    real, dimension(size(data_in,1), size(data_in,2) )        :: cfn, cfs, cfe, cfw, cfc
    real                                                      :: latp, latm, cstr, csm, csj
    real, dimension(size(data_in,1))                          :: dxu, dxt
    real, dimension(size(data_in,2))                          :: dyu, dyt

    ni = size(data_in,1)
    nj = size(data_in,2)

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
          if(abs(data_in(i,j) - missing_value) <= tol ) then
             tmp(i,j) = initial_guess
             sor(i,j) = rel_coef
          else
             tmp(i,j)=data_in(i,j)
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
          data_out(:,:) = tmp(1:ni,1:nj)
          exit
       endif

       !--- update boundaries

       call fill_boundaries(tmp, is_cyclic)

       n=n+1

    enddo

    if(mpp_pe() == mpp_root_pe() ) write(stdout(),'(a,i6,a)') 'Stopped after ',n,' iterations'
    if(mpp_pe() == mpp_root_pe() ) write(stdout(),'(a,f10.4)') 'maxres= ',resmax

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

end program regrid_2d








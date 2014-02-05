program edit_grid
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>

  !<OVERVIEW>
  ! edit grid topography.
  !</OVERVIEW>

  !<DESCRIPTION>
  ! This program can edit the topography of input grid_spec file "orig_grid" 
  ! according to the ascii input file "grid_edits". Then it will output the 
  ! new grid_spec file "mod_grid". The program read file "grid_edits" line
  ! by line. Each line contains grid points position and new topography value 
  ! of those grid points. The grid points position is specified by the grid 
  ! index. You can specify a point or a region at each line. For example,
  ! <PRE>
  ! 100, 60, 0 
  ! will set the depth at point (100,60) to 0. 
  ! 40:45, 30:34, 1000
  ! will set the depth at region ( index i from 54 to 50 and j from 30 to 34 ) to 1000.
  ! </PRE>
  ! 
  !</DESCRIPTION>

  use mpp_mod,       only : mpp_error, FATAL, NOTE, mpp_pe, mpp_npes
  use mpp_io_mod,    only : mpp_open, mpp_close, mpp_read, mpp_write, mpp_write_meta
  use mpp_io_mod,    only : MPP_RDONLY, MPP_NETCDF, MPP_SINGLE, MPP_ASCII, MPP_OVERWR
  use mpp_io_mod,    only : axistype, fieldtype, atttype, mpp_get_atts, mpp_get_info
  use mpp_io_mod,    only : mpp_get_axes, mpp_get_fields, mpp_get_axis_data, mpp_copy_meta
  use mpp_io_mod,    only : mpp_get_att_char, mpp_get_att_name, mpp_get_att_real_scalar
  use mpp_io_mod,    only : mpp_get_att_type, mpp_get_att_real
  use fms_mod,       only : fms_init, fms_end, file_exist, close_file, stdout
  use fms_mod,       only : open_namelist_file, check_nml_error, write_version_number
  use constants_mod, only : constants_init
  use topog_mod,     only : process_topo, show_deepest, set_topog_nml

  implicit none
#include <netcdf.inc>

  !--- namelist interface
 !
  !<NAMELIST NAME="edit_grid_nml">
  !<DATA NAME= "mod_grid" TYPE="character(len=128)" >
  ! original grid file
  !</DATA>
  !<DATA NAME= "orig_grid" TYPE="character(len=128)" >
  ! output grid file after modification.
  !</DATA>
  !<DATA NAME="grid_edits" TYPE="character(len=128)">
  !  input text file. Each line is in the format as 
  !  "is:ie, js:je, depth", which means set the depth at region
  !  (is:ie, js:je) to value "depth". is and ie can be equal or 
  !  different. js and je can be same or different.
  !</DATA>
  !<DATA NAME="debug" TYPE="logical" >
  ! Control standard output. Default value is false.
  !</DATA>
  !</NAMELIST>

  character(len=128) :: orig_grid       = 'orig_grid'       ! original grid file
  character(len=128) :: mod_grid        = 'mod_grid     '   ! modified grid file
  character(len=128) :: grid_edits      = 'grid_edits.txt'  ! test file used to edit grid
  logical            :: debug           = .FALSE.           ! will fail at compilation if included 
                                                            ! in the namelist

  namelist /edit_grid_nml/ orig_grid, mod_grid, grid_edits, debug

  !--- version information variables -----------------------------------
  character(len=128) :: version = '$Id: edit_grid.F90,v 19.0 2012/01/06 22:07:48 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'
  !---------------------------------------------------------------------
  logical :: tripolar_grid        =.false. ! indicate the grid is tripolar grid or not.
  logical :: cyclic_x             =.false. ! true indicate cyclic in x-direction
  logical :: cyclic_y             =.false. ! true indicate cyclic in y-direction
  logical :: full_cell            =.false. ! do not generate partial bottom cells 
  logical :: fill_isolated_cells  =.false. ! Do not allow non-advective tracer cells
  logical :: dont_change_landmask =.false. ! Do not change land/sea mask when filling isolated cells
  logical :: fill_shallow         =.false. ! Make cells less than minimum depth land
  logical :: deepen_shallow       =.false. ! Make cells less than minimum depth equal to minimum depth
  logical :: round_shallow        =.false. ! Make cells land if depth is less than 1/2 mimumim depth, 
                                           ! otherwise make ocean
  logical :: adjust_topo          =.false. ! adjust topography 
  logical :: fill_first_row       =.false. ! make first row of ocean model all land points for ice model
  integer :: kmt_min              = 2      ! minimum number of vertical levels

  real, dimension(:),      allocatable :: zw         ! vertical grid at T-cell bound
  real, dimension(:,:),    allocatable :: ht         ! depth at T-cell center
  real, dimension(:,:),    allocatable :: hu         ! depth at UV-cell center  
  real, dimension(:,:),    allocatable :: geolon_t   ! geographical longitude of T-cell center
  real, dimension(:,:),    allocatable :: geolat_t   ! geographical latitude of T-cell center
  real, dimension(:,:),    allocatable :: dxte, dytn 
  real, dimension(:,:),    allocatable :: wet        ! land/sea mask
  real, dimension(:,:),    allocatable :: wet_c        ! land/sea mask  
  integer, dimension(:,:), allocatable :: kmt        ! number of vertical levels
  integer, dimension(:,:), allocatable :: kmu        ! number of vertical levels  
  character(len=64)                    :: topography ! type of topography

  type(axistype), dimension(:), allocatable  :: axes_in, axes_out
  type(fieldtype), dimension(:), allocatable :: fields_in, fields_out
  type(atttype), dimension(:), allocatable   :: global_atts

  integer :: ni, nj, nk  ! grid size
  integer :: i,j
  integer :: unit_input  ! corresponding to the file orig_grid.    
  logical :: hu_exist


  !--- begin of the program 
  call fms_init
  call constants_init

  call edit_grid_init

  !--- first read from orig_grid
  call read_orig_grid

  if(topography == 'idealized')  then
      call mpp_error(FATAL,'edit_grid: can not edit grid when topography == idealized') 
  endif

  !--- read the grid edit file and edit ht
  call reset_depth

  !--- set topog_mod namelist option
  call set_topog_nml(full_cell, fill_isolated_cells, dont_change_landmask, fill_shallow, &
       deepen_shallow, round_shallow, adjust_topo, fill_first_row, kmt_min, debug )

  ! Compare "ht" to bottom of deepest model level. 
  if(debug) call show_deepest(zw, ht)

  call process_topo(ht, kmt, zw, dxte, dytn, geolon_t, geolat_t, tripolar_grid, cyclic_x, cyclic_y)

  if( hu_exist ) then
! MJH
     do j=1,nj-1
        do i=1,ni-1
           hu(i,j) = min(ht(i,j),ht(i+1,j), ht(i,j+1),ht(i+1,j+1))
           kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))        
        enddo
        hu(ni,j) = hu(ni-1,j)
        kmu(ni,j) = kmu(ni-1,j)
     enddo
     hu(:,nj) = hu(:,nj-1)
     kmu(:,nj) = kmu(:,nj-1)

     wet_c=0.0
     where (hu .gt. 0.0) wet_c = 1.0
  endif

  !--- define land/sea mask according to depth
  wet = 0.0
  where(ht(1:ni,1:nj) .gt. 0.0 ) wet = 1.0

  call write_mod_grid ()

  call fms_end

contains

  !#####################################################################
  ! --- read the namelist and write the version information to logfile. Also
  ! --- write the namelist to standard output
  subroutine edit_grid_init

    integer :: io_status, unit, ierr

    !--- This program do not support parallel, so always run on 1 pe
    if(mpp_npes() .gt. 1) call mpp_error(FATAL,'program edit_grid:  set npes = 1 in the runscripts')

    ! provide for namelist over-ride of defaults 
    if(file_exist('input.nml')) then
       unit = open_namelist_file()
       read (unit,edit_grid_nml,IOSTAT=io_status)
       write (stdout(),'(/)')
       write (stdout(),edit_grid_nml)  
       ierr = check_nml_error(io_status, 'edit_grid_nml')
       call close_file(unit)
    else
       call mpp_error(NOTE, 'edit_grid: file input.nml does not exist' )
    endif

    !--- write out version information ---------------------------------
    call write_version_number(version,tagname)

  end subroutine edit_grid_init

  !#####################################################################
  !--- read the grid information from the file orig_grid 
  subroutine read_orig_grid
    integer            :: ndim, nvar, natt, ntime, i, len
    character(len=128) :: name
    logical            :: zw_found, gridlon_found, gridlat_found
    logical            :: ht_found, hu_found, dxte_found, kmu_found
    logical            :: dytn_found, xt_found, yt_found
    real, dimension(:,:), allocatable :: tmp2d
    
    if(.not. file_exist(trim(orig_grid))) &
         call mpp_error(FATAL, 'file '//trim(orig_grid)//' does not exist')

    call mpp_open(unit_input,trim(orig_grid),MPP_RDONLY,MPP_NETCDF,threading=MPP_SINGLE,&
         fileset=MPP_SINGLE)
    call mpp_get_info(unit_input,ndim,nvar,natt,ntime)
    allocate (global_atts(natt))
    call mpp_get_atts(unit_input,global_atts)
    allocate(axes_in(ndim), axes_out(ndim) )
    call mpp_get_axes(unit_input,axes_in)
    axes_out = axes_in
    allocate(fields_in(nvar), fields_out(nvar))
    call mpp_get_fields(unit_input,fields_in)
    fields_out = fields_in

    do i=1,natt
       select case (trim(mpp_get_att_name(global_atts(i))))
       case ('x_boundary_type')
          if (trim(mpp_get_att_char(global_atts(i))) == 'cyclic') then
             cyclic_x = .true.
          else
             cyclic_x = .false.
          endif
       case ('y_boundary_type')
          if (trim(mpp_get_att_char(global_atts(i))) == 'fold_north_edge') then
             tripolar_grid = .true.
          else if(trim(mpp_get_att_char(global_atts(i))) == 'cyclic') then
             cyclic_y = .true.
          endif
       case ('topography')
          topography = mpp_get_att_char(global_atts(i))
       case ('full_cell')
          full_cell = .true.
       case ('fill_isolated_cells')
          fill_isolated_cells = .true.
       case ('dont_change_landmask')
          dont_change_landmask = .true.
       case ('fill_shallow')
          fill_shallow = .true.
       case ('deepen_shallow')
          deepen_shallow = .true.
       case ('round_shallow')
          round_shallow = .true.
       case ('kmt_min')
          kmt_min = mpp_get_att_real_scalar(global_atts(i))
       case ('adjust_topo')
          adjust_topo = .true.
       case ('fill_first_row')
          fill_first_row = .true.
       end select
    enddo

    do i=1,ndim
       call mpp_get_atts(axes_in(i),name=name,len=len)
       select case (trim(name))
       case ('zb') 
          allocate(zw(len))
          nk = len
          call mpp_get_axis_data(axes_in(i),zw)
          zw_found = .true.
       case ('grid_x_T')
          ni = len
          gridlon_found = .true.
       case ('grid_y_T')
          nj = len
          gridlat_found=.true.
       end select
    enddo

    if (.not.zw_found )      call mpp_error(FATAL,'edit_grid: axis zb not found in the file '//trim(orig_grid) )
    if (.not.gridlon_found ) call mpp_error(FATAL,'edit_grid: axis grid_x_T not found in the file '//trim(orig_grid) )
    if (.not.gridlat_found ) call mpp_error(FATAL,'edit_grid: axis grid_y_T not found in the file '//trim(orig_grid) )

    allocate( ht(0:ni+1,0:nj+1), hu(ni,nj),kmt(0:ni+1,0:nj+1), kmu(ni,nj), wet(ni,nj), wet_c(ni,nj)  )
    allocate( geolon_t(ni,nj), geolat_t(ni,nj), dxte(ni,nj), dytn(ni,nj) ) 
    ht = 0.0
    hu = 0.0
    kmt = 0
    kmu = 0

    ht_found = .false.
    hu_found = .false.
    dxte_found = .false.
    dytn_found = .false.
    xt_found = .false.
    yt_found = .false.
    kmu_found = .false.
    do i=1,nvar
       call mpp_get_atts(fields_in(i),name=name)
       select case (trim(name))
       case('depth_t')
          call mpp_read(unit_input,fields_in(i),ht(1:ni,1:nj) )
          ht_found = .true.
       case ('ds_01_21_E')
          call mpp_read(unit_input,fields_in(i),dxte)
          dxte_found = .true.
       case ('ds_10_12_N')
          call mpp_read(unit_input,fields_in(i),dytn)
          dytn_found = .true.
       case ('x_T')
          call mpp_read(unit_input,fields_in(i),geolon_t(1:ni,1:nj))
          xt_found = .true.
       case ('y_T')
          call mpp_read(unit_input,fields_in(i),geolat_t(1:ni,1:nj))
          yt_found = .true.
      case ('depth_c')
          call mpp_read(unit_input,fields_in(i),hu(1:ni,1:nj))          
          hu_found = .true.
      case ('num_levels_c')
          allocate(tmp2d(ni,nj))
          call mpp_read(unit_input,fields_in(i),tmp2d(1:ni,1:nj))          
          kmu(1:ni,1:nj) = tmp2d(1:ni,1:nj)
          deallocate(tmp2d)
          kmu_found = .true.
       end select
    enddo

    if (.not.ht_found)   call mpp_error(FATAL,'edit_grid: depth_t not found in the file '//trim(orig_grid) )
    if (.not.dxte_found) call mpp_error(FATAL,'edit_grid: ds_01_21_E not found in the file '//trim(orig_grid) )
    if (.not.dytn_found) call mpp_error(FATAL,'edit_grid: ds_10_12_N not found in the file '//trim(orig_grid) )
    if (.not.xt_found)   call mpp_error(FATAL,'edit_grid: x_T not found in the file '//trim(orig_grid) )
    if (.not.yt_found)   call mpp_error(FATAL,'edit_grid: y_T not found in the file '//trim(orig_grid) )
    if ( hu_found .neqv. kmu_found ) &
       call mpp_error(FATAL,'edit_grid: depth_c and num_levels_c should co-exist in the file '//trim(orig_grid) )

    hu_exist = hu_found
    if( hu_exist ) then
       write(stdout(),*)"NOTE from edit_grid: field depth_c and num_levels exist in the file "//trim(orig_grid)//".", &
            "Please be aware that the output topography should only be used in radiating open boundary condition model"
    endif

  end subroutine read_orig_grid

  !#####################################################################
  !--- write out the grid information after modification according to 
  !    grid_edits to file mod_grid.
  subroutine write_mod_grid
    integer                              :: unit, i, n
    real, dimension(:,:) ,   allocatable :: tmp2d
    real, dimension(:,:,:) , allocatable :: tmp3d
    integer, dimension(4)                :: siz_in = 1
    character(len=128)                   :: name

    call mpp_open(unit,trim(mod_grid),MPP_OVERWR,MPP_NETCDF,threading=MPP_SINGLE,&
         fileset=MPP_SINGLE)

    ! ----------------------------------------------------------------------
    ! write global atts
    ! ----------------------------------------------------------------------
    do i=1,size(global_atts(:))
       if (mpp_get_att_type(global_atts(i)) == NF_FLOAT) &
            call mpp_write_meta(unit,mpp_get_att_name(global_atts(i)), &
                 rval=mpp_get_att_real(global_atts(i)))
       if (mpp_get_att_type(global_atts(i)) == NF_CHAR) then
          if (trim(mpp_get_att_name(global_atts(i))) == 'filename') then
             call mpp_write_meta(unit,mpp_get_att_name(global_atts(i)),cval=trim(mod_grid))
          else
             call mpp_write_meta(unit,mpp_get_att_name(global_atts(i)), cval=trim(mpp_get_att_char(global_atts(i))))
          end if

       endif
    enddo

    !-----------------------------------------------------------------------
    ! write axis metadata
    ! ----------------------------------------------------------------------
    do i=1,size(axes_out(:))
       call mpp_copy_meta(unit,axes_out(i))
    enddo

    ! ----------------------------------------------------------------------
    ! write variable metadata
    ! ----------------------------------------------------------------------
    do i=1,size(fields_out(:))
       call mpp_copy_meta(unit,fields_out(i))
    enddo

    ! ----------------------------------------------------------------------
    ! write axis data
    ! ----------------------------------------------------------------------
    do i=1,size(axes_out(:))
       call mpp_write(unit,axes_out(i))
    enddo

    do i=1,size(fields_out(:))
       call mpp_get_atts(fields_out(i),name=name,siz=siz_in)
       select case (trim(name))
       case ('depth_t')
           call mpp_write(unit,fields_out(i),ht(1:ni,1:nj))
       case ('depth_c')
          call mpp_write(unit,fields_out(i),hu(1:ni,1:nj))           
       case ('num_levels')
          allocate(tmp2d(ni,nj))
          tmp2d = kmt(1:ni,1:nj)
          call mpp_write(unit,fields_out(i),tmp2d )
          deallocate(tmp2d)
       case ('num_levels_c')
           allocate(tmp2d(ni,nj))
          tmp2d = kmu(1:ni,1:nj)
          call mpp_write(unit,fields_out(i),tmp2d )
          deallocate(tmp2d)
       case ('wet')
           call mpp_write(unit,fields_out(i),wet)
       case ('wet_c')
           call mpp_write(unit,fields_out(i),wet_c)
       case default
          if(siz_in(3) > 1 ) then
             allocate(tmp3d(siz_in(1),siz_in(2), siz_in(3)) )
             call mpp_read(unit_input,fields_in(i),tmp3d)
             call mpp_write(unit,fields_out(i),tmp3d)
             deallocate(tmp3d)
          else 
             allocate(tmp2d(siz_in(1),siz_in(2)) )    
             call mpp_read(unit_input,fields_in(i),tmp2d)
             call mpp_write(unit,fields_out(i),tmp2d)
             deallocate(tmp2d)
          endif
       end select
    enddo

    call mpp_close(unit)
    call mpp_close(unit_input)

  end subroutine write_mod_grid

  !#####################################################################
  !--- reset depth according to file grid_edits
  subroutine reset_depth

    integer            :: unit, i, j, is, ie, js, je
    character(len=128) :: txt
    logical            :: flag
    real               :: ht_new

    if(.not. file_exist(trim(grid_edits))) &
         call mpp_error(FATAL, 'file '//trim(grid_edits)//' does not exist')

    call mpp_open(unit,trim(grid_edits),MPP_RDONLY,MPP_ASCII,threading=MPP_SINGLE,fileset=MPP_SINGLE)

    do 
       read(unit,'(a)',end=99,err=99) txt
       call parse_edits(txt,is,ie,js,je,ht_new,flag)
       if (flag) then
          do j=js,je
             do i=is,ie
                if (i < 0 .or. i > ni .or. j < 0 .or. j > nj) call mpp_error(FATAL,'indices exceed grid bounds')
                if (ht_new < 0.0 .or. ht_new > zw(nk)) call mpp_error(FATAL,'depth exceeds allowable depth')
                write(stdout(),*) 'Resetting depth at location (i,j) : ',i,j,' = ',ht_new
                ht(i,j) = ht_new
             enddo
          enddo
       endif
    enddo

99  call mpp_close(unit)

  end subroutine reset_depth

  !#####################################################################
  !--- read each line of the file grid_edits
  subroutine parse_edits(txt,is,ie,js,je,depth,flag)

    character(len=*), intent(in) :: txt
    integer, intent(inout)       :: is,ie,js,je
    real, intent(inout)          :: depth
    logical, intent(inout)       :: flag

    integer :: i1,i2, i3
    character(len=128) :: txt2

    flag = .true.
    i1 = scan(txt,',')
    txt2 = txt(1:i1-1)

    if (i1 <= 0) goto 90

    i2 = scan(txt2,':')
    if (i2 <= 0) then
       read(txt2,*,err=90) is
       ie = is
    else
       read(txt2(1:i2-1),*,err=90) is
       read(txt2(i2+1:),*,err=90) ie
    endif

    i2 = scan(txt(i1+1:),',')
    txt2 = txt(i1+1:i1+i2-1)
   if (i2 <= 0) goto 90

    i3 = scan(txt2,':')
     if (i3 <= 0) then
       read(txt2,*,err=90) js
       je = js
    else
       read(txt2(1:i3-1),*,err=90) js
       read(txt2(i3+1:),*,err=90) je
    endif

    read(txt(i1+i2+1:),*,err=90) depth

    return

90  continue

    flag = .false.

    return
  end subroutine parse_edits

  !#####################################################################

end program edit_grid


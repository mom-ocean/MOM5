program compare_grid
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
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>

  !<OVERVIEW>
  ! compare depth and land/sea mask of two grid_spec file. 
  !</OVERVIEW>

  !<DESCRIPTION>
  ! This program reads in two grid descriptor files (generated via ocean_grid_generator) 
  ! and creates a text file output listing line-by-line differences between the
  ! two files. Output file format is the same as the grid_edits file used by
  ! edit_grid.F90. These two files should have same grid size. 
  !
  ! Originally developed by Jeffery B. Greenblatt on 12/11/2001 at Princeton University     
  !</DESCRIPTION>

  use mpp_mod,    only : mpp_error, FATAL, NOTE, mpp_pe, mpp_init, mpp_exit
  use mpp_io_mod, only : mpp_open, mpp_read, mpp_close, fieldtype
  use mpp_io_mod, only : mpp_get_info, mpp_get_atts, mpp_get_fields
  use mpp_io_mod, only : MPP_RDONLY, MPP_NETCDF, MPP_SINGLE, MPP_OVERWR, MPP_ASCII 
  use fms_mod,    only : fms_init, fms_end, file_exist, stdout, close_file
  use fms_mod,    only : open_namelist_file, check_nml_error, write_version_number

  implicit none

  !--- namelist information --------------------------------------------
  !
  !<NAMELIST NAME="compare_grid_nml">
  !<DATA NAME="grid_file_1" TYPE="character(len=128)">
  !  First grid files to be compared with grid_file_2.
  !</DATA>
  !<DATA NAME="grid_file_2" TYPE="character(len=128)">
  !  Second grid files to be compared with grid_file_1.
  !</DATA>
  !<DATA NAME="grid_edits" TYPE="character(len=128)">
  !  output text file. Each line is in the format as 
  !  "i, j, depth_new, #was depth_old ". depth_new is 
  !  the depth at point (i,j) of grid_file_2 and depth_old 
  !  is the depth at point (i,j) of grid_file_1. 
  !</DATA>
  !<DATA NAME="mask_diff" TYPE="character(len=128)">
  !  output text file. Each line is in the format as 
  !  "i, j, wet_new, #was wet_old ". wet_new is 
  !  the land/sea mask at point (i,j) of grid_file_2 and wet_old 
  !  is the land/sea mask at point (i,j) of grid_file_1. 
  !</DATA>
  !</NAMELIST>
  character(len=128) :: grid_file_1 = 'grid_file_1'
  character(len=128) :: grid_file_2 = 'grid_file_2'
  character(len=128) :: grid_edits  = 'grid_edits.txt'
  character(len=128) :: mask_diff   = 'mask_diff.txt'

  namelist /compare_grid_nml/ grid_file_1, grid_file_2, grid_edits, mask_diff

  !--- version information variables -----------------------------------
  character(len=128) :: version = '$Id: compare_grid.f90,v 11.0 2004/09/28 20:07:16 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  !--- compare_grid_type
  type compare_grid_type
     real,    dimension(:,:), pointer :: ht => NULL()
     integer, dimension(:,:), pointer :: wet => NULL()
     real,    dimension(:,:), pointer :: geolon_t => NULL()
     real,    dimension(:,:), pointer :: geolat_t => NULL()
     integer                          :: ni, nj
  end type compare_grid_type

  !--- other variables
  type(compare_grid_type) :: grid_1
  type(compare_grid_type) :: grid_2

  call fms_init ()

  call compare_grid_init ()

  ! --- Read grid files ------------------------------------------------

  call read_grid(grid_file_1, grid_1)
  call read_grid(grid_file_2, grid_2)

  ! Ensure we are comparing same-size arrays
  call grid_check()

  ! Now compare selected arrays and write out differences
  call grid_compare()

  call fms_end()

contains 
  !#####################################################################
  ! --- read the namelist and write the version information to logfile. Also
  ! --- write the namelist to standard output
  subroutine compare_grid_init

    integer :: io_status, unit, ierr

 ! --- provide for namelist over-ride of defaults ---------------------
  if(file_exist('input.nml')) then
     unit = open_namelist_file()
     read (unit,compare_grid_nml,IOSTAT=io_status)
     write (stdout(),'(/)')
     write (stdout(),compare_grid_nml)  
     ierr = check_nml_error(io_status, 'compare_grid_nml')
     call close_file(unit)
  else
     call mpp_error(NOTE, 'file input.nml does not exist' )
  endif

   !--- write out version information ---------------------------------
    call write_version_number(version, tagname)

  end subroutine compare_grid_init

  !#####################################################################

  !--- read the grid information from file.
  subroutine read_grid(file, grid)
    character(len=*),           intent(in) :: file   
    type(compare_grid_type), intent(inout) :: grid

    !--- Local variables
    integer            :: i, j, unit, len, ndim, nvar, natt, ntime
    integer            :: siz_in(3)
    logical            :: ht_found, wet_found, xt_found, yt_found
    character(len=128) :: name
    type(fieldtype), dimension(:), allocatable :: fields
    real, dimension(:,:), allocatable :: tmp 

    if(.not. file_exist(trim(file)) ) &
        call mpp_error(FATAL, 'compare_grid: file '//trim(file)//' does not exist')

    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_SINGLE,&
         fileset=MPP_SINGLE)
    call mpp_get_info(unit,ndim,nvar,natt,ntime)

    allocate(fields(nvar))
    call mpp_get_fields(unit,fields)

    ht_found=.false.
    wet_found=.false.  

    do i=1,nvar
       call mpp_get_atts(fields(i),name=name,siz=siz_in)
       select case (trim(name))
       case('depth_t')
          grid%ni = siz_in(1); grid%nj = siz_in(2)
          allocate(grid%ht(grid%ni,grid%nj))
          call mpp_read(unit,fields(i),grid%ht)
          ht_found = .true.
       case ('wet')
          allocate(grid%wet(siz_in(1), siz_in(2)),tmp(siz_in(1), siz_in(2)))
          call mpp_read(unit,fields(i),tmp)
          grid%wet = tmp
          deallocate(tmp)
          wet_found = .true.    
       case('x_T')
          allocate(grid%geolon_t(siz_in(1), siz_in(2)))
          call mpp_read(unit,fields(i),grid%geolon_t)
          xt_found = .true.
       case('y_T')
          allocate(grid%geolat_t(siz_in(1), siz_in(2)))
          call mpp_read(unit,fields(i),grid%geolat_t)
          yt_found = .true.
       end select
    enddo

    if (.not. ht_found) call mpp_error(FATAL,'compare_grid: depth_t not found in the file '//trim(file) )
    if (.not. wet_found) call mpp_error(FATAL,'compare_grid: wet not found in the file '//trim(file) )
    if (.not. xt_found) call mpp_error(FATAL,'compare_grid: x_T not found in the file '//trim(file) )
    if (.not. yt_found) call mpp_error(FATAL,'compare_grid: x_T not found in the file '//trim(file) )

    call mpp_close(unit)

  end subroutine read_grid
 
  !#####################################################################
  ! compare grid size and geographic grid.
  subroutine grid_check
    integer :: i, j

    if(grid_1%ni .ne. grid_2%ni .or. grid_1%nj .ne. grid_2%nj ) then
       call mpp_error(FATAL,'compare_grid: grid sizes of '//trim(grid_file_1)//  &
            ' and'//trim(grid_file_2)//' do not match')
    endif

    !--- make sure it is the same grid
    do j = 1, grid_1%nj
       do i = 1, grid_1%ni
          if(grid_1%geolon_t(i,j) .ne. grid_2%geolon_t(i,j) .or.  grid_1%geolat_t(i,j) .ne. grid_2%geolat_t(i,j)) then
               write(stdout(),*) 'grid_1%geolon_t(',i,',',j,') = ', grid_1%geolon_t(i,j),  &
                                'grid_2%geolon_t(',i,',',j,') = ', grid_2%geolon_t(i,j)
               write(stdout(),*) 'grid_1%geolat_t(',i,',',j,') = ', grid_1%geolat_t(i,j),  &
                                'grid_2%geolat_t(',i,',',j,') = ', grid_2%geolat_t(i,j)
               call mpp_error(FATAL, 'file '//trim(grid_file_1)//' and file '//trim(grid_file_2)// &
               ' do not have the same geographical grid' )
          endif
       enddo
    enddo

  end subroutine grid_check

  !#####################################################################
  !--- compare the two grid file and write out the difference to file grid_edits
  subroutine grid_compare

  logical :: flag
  integer :: i, j, unit

  !--- check the depth_t
  call mpp_open(unit,trim(grid_edits),MPP_OVERWR,MPP_ASCII,threading=&
       MPP_SINGLE,fileset=MPP_SINGLE, nohdrs = .TRUE.)

  flag = .false.
  write(unit, *) '#   i     j           depth_new'
  write(unit, *) '#'
  do i = 1, grid_1%ni
     do j = 1, grid_2%nj
        if (grid_1%ht(i,j) .ne. grid_2%ht(i,j)) then
           write(unit, '(i6,a,i6,a,f24.15,a,f24.15)') i,',',j,',',grid_2%ht(i,j),'   # was ',grid_1%ht(i,j)
           flag = .true.
        endif
     enddo
  enddo
  if (.not. flag) write(unit,*) '# NO CHANGES in ht'
  call mpp_close(unit)

  !--- check land_sea mask
  call mpp_open(unit,trim(mask_diff),MPP_OVERWR,MPP_ASCII,threading=&
       MPP_SINGLE,fileset=MPP_SINGLE)

  write(unit, *) '#   i     j    mask_new'
  write(unit, *) '#'
  flag = .false.
  do i = 1, grid_1%ni
     do j = 1, grid_1%nj
        if (grid_1%wet(i,j) .ne. grid_2%wet(i,j)) then
           write(unit, '(i6,a,i6,a,i6,a,i6)') i,',',j,',',grid_2%wet(i,j),'        # was ', grid_1%wet(i,j)
           flag = .true.
        endif
     enddo
  enddo
  if (.not. flag) write(unit,*) '# NO CHANGES in wet'
  call mpp_close(unit)

  end subroutine grid_compare

  !#####################################################################

end program compare_grid








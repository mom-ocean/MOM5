module grids_util_mod
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>


  !<OVERVIEW>
  ! <TT>grids_util_mod</TT> contains some public interface used by several modules in
  !  generate_ocean_grid package. 
  !</OVERVIEW>

  use fms_mod,       only : stdout, string
  use mpp_mod,       only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_chksum, FATAL, NOTE 
  use mpp_io_mod,    only : mpp_open, MPP_OVERWR, MPP_NETCDF, MPP_SINGLE, axistype, fieldtype
  use mpp_io_mod,    only : mpp_write_meta, mpp_write, mpp_close
  use constants_mod, only : PI

  implicit none
  private

  public :: make_axis, gcell, get_file_unit, write_field_meta, write_field_data, set_grid

  !--- public interface ------------------------------------------------
  ! <INTERFACE NAME="write_field_data">
  !   <OVERVIEW>
  !     Write data to corresponding grid file. 
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     For the purpose of generating higher resolution grid. There is a
  !     2 GB limit of each grid file. The interface will allow the grid information
  !     be stored at multiple files.
  !   </DESCRIPTION>
  !   <IN NAME="filename" TYPE= "character(len=*)">
  !    The name of the grid file to be generated. If the grid file is 
  !    over 2 GB limit, it will break into several files with file name filename,
  !    filename2, filename3 ....
  !   </IN>
  !   <IN NAME="fieldname" TYPE= "character(len=*)">
  !    name of the field to be written into the file filename
  !   </IN>
  !   <IN NAME="fielddata" TYPE= "real, dimension(:,:) or dimension(:,:,:)">
  !    data of fieldname to be written to the file filename.
  !   </IN>
  ! </INTERFACE>
  !
  interface write_field_data
     module procedure  write_field_data_2d
     module procedure  write_field_data_3d
  end interface write_field_data

  real, dimension(:), allocatable :: xt, yt, xc, yc
  integer                         :: ni, nj               ! grid size
  logical                         :: is_root_pe = .false. ! to indicate if current pe is root pe.

  real, parameter :: max_file_size = 4294967295.  !The maximum size in bytes for any one file. 
                                                  ! With NetCDF3, this  should be 2 Gb or less. 
                                                  ! the reason using real number is the maximum integer.
  real, parameter    :: tolr = 1.0e-4

  integer, parameter :: max_fields = 300    ! can be increased if needed
  integer, parameter :: max_files = 50      ! can be increased if needed
  integer            :: num_files = 0       ! if >= 0, means a series file will be created.
  integer            :: num_fields = 0      ! number of fields written to the file.


  character(len=128)   :: opened_files(0:max_files) = ''
  integer              :: files_unit(0:max_files) = 0
  type(axistype),save  :: axis_xt(0:max_files), axis_yt(0:max_files)
  type(axistype),save  :: axis_xc(0:max_files), axis_yc(0:max_files)
  type(axistype),save  :: axis_v(0:max_files)
  type(fieldtype),save :: flds(max_fields)
  character(len=128)   :: flds_name(max_fields)
  integer              :: file_num(max_fields) = -1

contains

   !#######################################################################
  ! <SUBROUTINE NAME="gcell">
  !   <OVERVIEW>
  !     grid cell construction.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     A domain is composed of one or more regions: Build "num" T cells with resolution 
  !     "deltat(n) n=1,num" within the domain composed of regions bounded by "bounds".
  !     Also construct "num" C-cells of resolution "deltau(n) n=1,num" with the relation 
  !     between T and U cells given by: deltat(n) = 0.5*(deltau(n-1) + deltau(n)).
  !     Resolution may be constant or smoothly varying within each region AND there must 
  !     be an integral number of grid cells within each region. The domain is the sum of all regions.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call gcell (maxlen, n_bounds, bounds, d_bounds, nbpts, num, deltat, deltau, stretch)
  !   </TEMPLATE>
  !   <IN NAME="maxlen" TYPE="integer" >
  !     maximum length of "deltat" and "deltau"
  !   </IN>
  !   <IN NAME="n_bounds" TYPE="integer" >
  !     number of bounds needed to define the regions
  !   </IN>
  !   <IN NAME="bounds" TYPE="real, dimension(n_bounds)" >
  !     latitude, longitude, or depth at each bound
  !   </IN>
  !   <IN NAME="d_bounds" TYPE="real, dimension(n_bounds)" >
  !     delta (resolution) at each of the "bounds"
  !   </IN>
  !   <IN NAME="nbpts" TYPE="integer" >
  !     number of extra boundary cells to add to the domain. (usually one at the beginning and end)
  !   </IN>
  !   <IN NAME="stretch" TYPE="real" >
  !     stretching factor for last region (should only be used in the vertical) to provide
  !     increased stretching of grid points. "stretch" = 1.0 gives no increased stretching.
  !     "stretch" = 1.2 gives increased stretching...etc
  !   </IN>
  !   <IN NAME="debug" TYPE="logical, optional" >
  !      flag that controls standard output.
  !   </IN>
  !   <OUT NAME="num" TYPE="integer" >
  !     total number of grid cells within the domain
  !   </OUT>
  !   <OUT NAME="deltat" TYPE="real, dimension(1-nbpts:maxlen)" >
  !     resolution of T grid cells: n=1,num
  !   </OUT>
  !   <OUT NAME="deltau" TYPE="real, dimension(1-nbpts:maxlen)" >
  !     resolution of C grid cells: n=1,num
  !   </OUT>
  ! </SUBROUTINE>
  subroutine gcell (region, bnd1, bnd2, dbnd1, dbnd2, num, ntotal, maxlen, deltau)
    integer, intent(in)             :: region
    real,    intent(in)             :: bnd1, bnd2, dbnd1, dbnd2
    integer, intent(out)            :: num
    integer, intent(in)             :: ntotal, maxlen
    real, dimension(1:maxlen), intent(out) :: deltau

    real    :: avg_res, chg_res, wid, an, err, d_new, del
    integer :: m, i

    avg_res = 0.5*(dbnd1 + dbnd2)
    chg_res = dbnd2 - dbnd1
    wid = abs(bnd2 - bnd1)
    an  = wid/avg_res
    m   = nint(an)
    err = wid - avg_res * m

    if( abs(err) > tolr ) then
       write(stdout(),*) "==>Error: non integral number of cells in region #", region, &
            ", average resolution within region =",avg_res,                            &
            ", this implies ",an," grid cells, Change grid specifications ",           &
            "via table, Here is some help..."
       if(m>1) then
          d_new = (2.0*wid)/(m-1) - dbnd1;
          write(stdout(),*) "Note: to get ",m-1," grid cells within region ",region, &
               ", change resolution from ",dbnd2," to ",d_new
       endif

       d_new = (2.0*wid)/m - dbnd1
       write(stdout(),*)"Note: to get ",m," grid cells within region ",region, &
            ", change resolution from ",dbnd2," to ",d_new      
       d_new = (2.0*wid)/(m+1) - dbnd1
       write(stdout(),*)"Note: to get ",m+1," grid cells within region ",region, &
            ", change resolution from ",dbnd2," to ",d_new

       call mpp_error(FATAL);
    endif

    if(ntotal+m > maxlen) call mpp_error(FATAL, "grids_util_mod: maxlen exceeded in gcell, "// &
                                      "increase size of maxlen in hgrid_mod or vgrid_mod")

    ! Calculate resolution of corner cells: "deltau"
    ! T grid points will be centered in these cells

    do i = 1, m
       del = avg_res - 0.5*chg_res*cos((PI/m)*(i-0.5))
       deltau(i) = del
    enddo

    num = m

    return

  end subroutine gcell

!///////////////////////////////////////////////////////////
!/
!/ Make_axis
!/ define horizontal grid resolution and axis data
!/
!///////////////////////////////////////////////////////////

subroutine make_axis(cart, maxlen, num_regions,bounds,delta,corners,centers, &
                     num, square_grid, extend_square_grid )
    character(len=1),       intent(in) :: cart
    integer,                intent(in) :: maxlen
    integer,             intent(inout) :: num_regions
    real,  dimension(:), intent(inout) :: bounds, delta
    real, dimension(0:), intent(inout) :: corners, centers
    integer,             intent(inout) :: num
    logical,                intent(in) :: square_grid, extend_square_grid

  integer:: i, n, m, ii
  real   :: delta_corners(maxlen)

  if(cart == 'Y') then  
     if(square_grid) then
        call iso_grid (maxlen, num_regions, bounds, delta, num, corners, centers, extend_square_grid)
        return
     endif
   endif

   !--- bounds should increase monotonically
   do n = 2, num_regions
      if( bounds(n-1) .gt. bounds(n)) then
         do m=1,num_regions
            write (stdout(),'(i3,f10.5)') m, bounds(m)
         end do
         call mpp_error(FATAL, &
              ' grids_util_mod: longitude boundaries with cart = '// cart //&
              'do not increase monotonically ')
      endif
   enddo

  num = 0;
  do n = 2, num_regions
    write(stdout(),*) " region # ",n," going from ",bounds(n-1)," (res=",delta(n-1), &
           ") to ",bounds(n)," (res=",delta(n),")"
    
    call gcell (n, bounds(n-1), bounds(n), delta(n-1), delta(n), m, num, maxlen, delta_corners)

    !*********************************************************************
    !   Build the grid points on a "B" grid. The T and C
    !   cells are staggered in the horizontal but at the same level in
    !   the vertical. However, the W cells (for vertical
    !   advection velocities at the bottoms of the C and T cells) are
    !   staggered in the vertical with respect to C and T cells.
    !   (no extra boundary point at the start)
    !*********************************************************************/
    corners(num) = bounds(n-1)

    do i = 1, m
      ii = num+i
      corners(ii) = corners(ii-1) + delta_corners(i)
    enddo
    num = num+m
  enddo

  corners(num) = bounds(num_regions)
  corners(num+1) = 2*corners(num) - corners(num-1)
  centers(0) = bounds(1) - 0.5*delta(1)

  do i = 1, num+1
     centers(i) = 2.0*corners(i-1) - centers(i-1)
  enddo 

  return   

end subroutine make_axis

  !#######################################################################
  subroutine iso_grid (maxlen, nylats, y_lat, dy_lat, nj, corners, centers, extend_square_grid)
    ! compute latitude resolution of grid cells to match the convergence
    ! of meridians.
    ! 
    ! inputs:
    !
    ! maxlen   = maximum length of "dytdeg0" and "dyudeg0"
    ! nylats   = should equal 2 (defines one region)
    ! dy_lat   = latitudinal resolution of grid cell on equator
    ! y_lat(1) = southern boundary of the domain (it will be adjusted
    !            to fit an integral number of cells)
    ! y_lat(2) = northern boundary of the domain (it will be adjusted
    !            to fit an integral number of cells)
    !
    ! outputs:
    !
    ! nj      = number of grid cells
    ! yt0     = latitude of point within T cell (degrees)
    ! yu0     = latitude of point within U cell (degrees)

    integer,                   intent(in)    :: maxlen, nylats
    integer,                   intent(out)   :: nj
    real, dimension(nylats),   intent(inout) :: y_lat
    real, dimension(nylats),   intent(in)    :: dy_lat
    real, dimension(0:maxlen), intent(out)   :: corners, centers
    logical,                   intent(in)    :: extend_square_grid
    !--- local variables -------------------------------------------------
    integer,                    parameter :: jmaxlat = 1000, lenjmax = 2*jmaxlat+1
    real, dimension(-jmaxlat-1:jmaxlat+1) :: dusq, dtsq, usq, tsq 
    real,                    dimension(2) :: y_bound, dy_bound
    integer                               :: n, n1, n2, nps, npn, j1, jmts, jmtn, j
    real                                  :: pole1, pole2, wid, stretch, D2R
    !---------------------------------------------------------------------

    D2R = PI/180.

    if (nylats .gt. 2) then
       call mpp_error(FATAL,'grids_util_mod: when hgrid_nml "square_grid"=true, hgrid_nml "nylats"= ' &
                 //trim(string(nylats))//' should eqaul 2')
    endif

    if (dy_lat(1) .ne. dy_lat(2)) then
       call mpp_error(FATAL,'grids_util_mod: nml dy_lat(1)= '//trim(string(dy_lat(1)))// &
             ' must equal dy_lat(2)= '//trim(string(dy_lat(2)))//' when hgrid_nml "square_grid"=true')
    endif
    if ((abs(y_lat(1)) > 89.9999) .or. (abs(y_lat(2)) > 89.9999)) then
       call mpp_error(FATAL,'grids_util_mod: Cannot specify hgrid_nml "ylat" = 90 deg ' &
                     //'when hgrid_nml "square_grid"=true ')
    endif

    ! build a square grid

    usq(0)  = 0.0
    dusq(0) = dy_lat(1)
    tsq(1)  = 0.5*dusq(0)
    tsq(0)  = -tsq(1)
    dtsq(1) = dusq(0)*cos(tsq(1)*D2R)
    dtsq(0) = dtsq(1)
    do n=1,jmaxlat
       dusq(n)   = 2.0*dtsq(n) - dusq(n-1)
       usq(n)    = tsq(n) + 0.5*dusq(n)
       if (tsq(n) .lt. 90.0) then
          tsq(n+1)  = tsq(n) + dusq(n)
          dtsq(n+1) = dusq(0)*cos(tsq(n+1)*D2R)
       else
          tsq(n+1)  = tsq(n)
          dtsq(n+1) = 0.0
       endif
    enddo
    do n=-1,-jmaxlat,-1
       dusq(n) = dusq(-n)
       usq(n)  = -usq(-n)
       dtsq(n) = dtsq(-(n-1)) 
       tsq(n)  = -tsq(-(n-1))
    enddo

    ! pick out cells between bounding latitudes

    n1 = -jmaxlat
    n2 = jmaxlat
    do n=0,jmaxlat
       if (usq(n) > y_lat(2)) then
          n2 = n
          exit
       endif
    enddo
    do n=0,-jmaxlat,-1
       if (usq(n) < y_lat(1)) then
          n1 = n+1
          exit
       endif
    enddo

    ! re-define bounding latitudes to match square boundaries
    y_lat(1) = usq(n1-1)
    y_lat(2) = usq(n2)  

    do n=n1,n2
       write(stdout(),*) 'n=',n,' tsq=',tsq(n),', dyt=',dtsq(n)
    enddo
    if (n1 .eq. -jmaxlat .or. n2 .eq. jmaxlat) then
       call mpp_error(FATAL,'hgrid_mod: Need to increase jmaxlat to reach the max latitude')
    endif

    if (extend_square_grid) then

       ! extend square grid in southern hemisphere to South pole

       ! set south pole (pole1) to skip over land in antarctica
       pole1 = -80.0
       stretch = 1.0
       wid = (y_lat(1) - pole1)
       nps = nint(wid/dusq(n1)) - 1
       y_bound(1) = pole1
       y_bound(2) = y_lat(1)
       dy_bound(1) = 2.0*wid/nps - dusq(n1)
       dy_bound(2) = dusq(n1)
       j1 = n1-nps
       call gcell(1, y_bound(1), y_bound(2), dy_bound(1), dy_bound(2), jmts, 0, maxlen, dtsq(j1) ) 
       dusq(j1) = 2*dtsq(j1)-dy_bound(1)
       do n = 1, jmts-1
          dusq(j1+n) = 2*dtsq(j1+n)-dusq(j1+n-1)
       enddo

       ! construct latitude of T and U grid points
       do n=n1-1, n1-jmts+1, -1
          tsq(n) = tsq(n+1) - dusq(n)
          usq(n) = tsq(n+1) - 0.5*dusq(n)
          write(stdout(),*)'extended n=',n,' usq=',usq(n),'tsq=',tsq(n), 'dyt=',dusq(n)
       enddo

       n1 = n1 - jmts + 1

       ! extend square grid to North pole

       pole2 = 90.0
       npn = nint((pole2 - y_lat(2))/dusq(n2)) - 1
       wid = (pole2 - y_lat(2))
       y_bound(1) = y_lat(2)
       y_bound(2) = pole2
       dy_bound(1) = dusq(n2)
       dy_bound(2) = 2.0*wid/npn - dusq(n2)

       call gcell(1, y_bound(1), y_bound(2), dy_bound(1), dy_bound(2), jmtn, 0, maxlen, dtsq(n2+1) ) 
       dusq(n2+1) = 2*dtsq(n2+1)-dy_bound(1)

       do n = 1, jmtn-1
          dusq(n2+n) = 2*dtsq(n2+n)-dusq(n2+n-1)
       enddo

       ! construct latitude of T and U grid points
       do n = n2+1, n2+jmtn-1
          tsq(n) = tsq(n-1) + dusq(n)
          usq(n) = tsq(n) + 0.5*dusq(n)
          write(stdout(),*)'extended n=',n,' usq=',usq(n),'tsq=',tsq(n), 'dyt=',dusq(n)
       enddo
       n2 = n2 + jmtn - 1
       ! re-define bounding latitudes to match extended boundaries

       y_lat(1) = pole1
       y_lat(2) = pole2 
    endif
    nj = n2 - n1 + 1

    do j=1,nj
       corners(j)     = usq(j+n1-1)
       centers(j)     = tsq(j+n1-1)
    enddo
    
    corners(0) = 2*centers(1) - corners(1)
    centers(0) = 2*corners(0) - centers(1)

    centers(nj+1) = 2*corners(nj)-centers(nj)
    corners(nj+1) = 2*centers(nj+1)- corners(nj)

  end subroutine iso_grid


  !#######################################################################
  ! This function is only for global meta and vgrid. So get_file_unit should 
  ! always linked to file opened_file(0)
  ! <FUNCTION NAME="get_file_unit">
  !   <OVERVIEW>
  !    returns the io unit corresponding to filename.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    If the file filename is already open, return the io unit of this 
  !    opened file. Otherwise will open the file and return the io unit.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     get_file_unit(filename)
  !   </TEMPLATE>
  !   <IN NAME="filename" TYPE= "character(len=*)">
  !    The name of the grid file to be generated.
  !   </IN>

  function get_file_unit(filename)   
    character(len=*), intent(in) :: filename
    integer                      :: get_file_unit

    if(.not. is_root_pe) return

    if(trim(opened_files(0)) == trim(filename) ) then 
          get_file_unit = files_unit(0)
          return
    endif

    !--- if file is not opened, open the file

    call mpp_open(get_file_unit, trim(filename), MPP_OVERWR,MPP_NETCDF, &
          threading=MPP_SINGLE,  fileset=MPP_SINGLE )

   files_unit(num_files) = get_file_unit
   opened_files(num_files) = trim(filename)

  end function get_file_unit
  ! </FUNCTION>
  !#####################################################################
  ! <SUBROUTINE NAME="write_field_meta">
  !   <OVERVIEW>
  !    Write meta data of a field to a netcdf file.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    It will check if the grid file will over the 2 GB limit. If do, will open 
  !    a new file with name filename? (? is 1, 2, 3 ....) and write axis metadata
  !    to the new file.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !    call write_field_meta(filename, fieldname, units, field_longname, fielddim, x_pos, y_pos)
  !   </TEMPLATE>
  !   <IN NAME="filename" TYPE= "character(len=*)">
  !    The name of the grid file to be generated. If the grid file is 
  !    over 2 GB limit, it will break into several files with file name filename,
  !    filename1, filename2, filename3 ....
  !   </IN>
  !   <IN NAME="fieldname" TYPE= "character(len=*)">
  !    name of the field to be written into the file filename
  !   </IN>
  !   <IN NAME="units" TYPE= "character(len=*)">
  !    units of field fieldname.
  !   </IN>
  !   <IN NAME="field_longname" TYPE= "character(len=*)">
  !    longname of fielname.
  !   </IN>
  !   <IN NAME="fielddim" TYPE= "integer">
  !    Indicate the dimension of fieldname. fielddim should be either 2 or 3. 
  !   </IN>
  !   <IN NAME="x_pos, y_pos" TYPE= "character(len=1)">
  !    To indicate the cell position. its value can be "T" or "C".
  !   </IN>
  subroutine write_field_meta(filename, fieldname, units, field_longname, fielddim, x_pos, y_pos)
    character(len=*),           intent(in) :: filename, fieldname
    character(len=*),           intent(in) :: units, field_longname
    integer,                    intent(in) :: fielddim
    character(len=1), intent(in), optional :: x_pos, y_pos

    character(len=128) :: curr_file, basename
    logical            :: is_first_field 
    integer            :: unit, fieldsize, length
    type(axistype)     :: axis_x, axis_y
    logical, save      :: do_write_axis = .true.
    real, save         :: size_in_file = 20000.0

    if(.not. is_root_pe) return

    curr_file = filename
    num_fields = num_fields + 1
    flds_name(num_fields) = trim(fieldname)

    if(num_fields ==1) is_first_field = .true.

    if(size_in_file .gt. max_file_size) then   ! 
       num_files = num_files + 1
       if(num_files .gt. max_files) call mpp_error(FATAL, &
           'grids_util_mod: number of files is over max_files, increase max_files ')

    endif

    if(num_files .gt. 0) then
       length = len_trim(filename)
       if(filename(length-2:length) == '.nc') then
          basename = filename(1:length-3)
       else
          basename = filename(1:length)
       endif
       if(num_files < 10)  then   ! num_files should be less than 100
          write(curr_file, '(a,i1,a)')trim(basename), num_files, '.nc'
       else
          write(curr_file, '(a,i2,a)')trim(basename), num_files, '.nc'
       endif
    endif

    if(trim(opened_files(num_files)) == curr_file) then  ! the file is already opened
       unit = files_unit(num_files)
    else   ! need to open the file and write out axis meta.
       call mpp_open(unit, trim(curr_file), MPP_OVERWR,MPP_NETCDF, &
            threading=MPP_SINGLE,  fileset=MPP_SINGLE )

       opened_files(num_files) = trim(curr_file)
       files_unit(num_files) = unit
       do_write_axis = .true.
    endif
    file_num(num_fields) = num_files
    if(do_write_axis) then
       !--- write axis meta -------------------------------------------
       call mpp_write_meta(unit, axis_xt(num_files),'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
            cartesian ='X', data = xt )
       call mpp_write_meta(unit, axis_yt(num_files),'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
            cartesian ='Y', data = yt )
       call mpp_write_meta(unit, axis_xc(num_files),'grid_x_C','degree_east','Nominal Longitude of C-cell center', &
            cartesian ='X', data = xc )
       call mpp_write_meta(unit, axis_yc(num_files),'grid_y_C','degree_north','Nominal Latitude of C-cell center', &
            cartesian ='Y', data = yc )
       call mpp_write_meta(unit, axis_v(num_files), 'vertex', 'none ',        &
            'Vertex position from southwest couterclockwise', data = (/1.,2.,3.,4./) )   
       do_write_axis = .false.
       !--- open a new file, set size_in_file to 20000
       size_in_file = 20000.

    endif

    if(present(x_pos)) then
      if (x_pos .ne. 'T' .and. x_pos .ne. 'C')  &
        call mpp_error(FATAL,'grids_util_mod: x_pos should be either "T" or "C" ')
    endif

    if(present(y_pos)) then
        if (y_pos .ne. 'T' .and. y_pos .ne. 'C') &
        call mpp_error(FATAL,'grids_util_mod: x_pos should be either "T" or "C" ')
    endif

    axis_x = axis_xt(num_files)
    if(present(x_pos))then
       if( x_pos == 'C') axis_x = axis_xc(num_files)
    endif
    axis_y = axis_yt(num_files)
    if(present(y_pos))then
       if(y_pos == 'C') axis_y = axis_yc(num_files)
    endif
    !--- write out field meta
    if(fielddim==2) then
       fieldsize = 8*ni*nj
       call mpp_write_meta(unit, flds(num_fields), (/axis_x, axis_y/), fieldname, units, field_longname, pack=1) 
    else if(fielddim==3) then
       fieldsize = 32*ni*nj
       call mpp_write_meta(unit, flds(num_fields), (/axis_x, axis_y,axis_v(num_files)/), &
            fieldname, units, field_longname, pack=1) 
    else
       call mpp_error(FATAL, 'grids_util_mod: dimension of field '//trim(fieldname)//' should be either 2 or 3')
    endif

    size_in_file = size_in_file + fieldsize


  end subroutine write_field_meta
  !  </SUBROUTINE>
  !#####################################################################
  !  <SUBROUTINE NAME="write_field_data_2d" INTERFACE="write_field_data">
  !  <IN NAME="filename" TYPE= "character(len=*)"> </IN>
  !   <IN NAME="fieldname" TYPE= "character(len=*)"> </IN>
  !   <IN NAME="fielddata" TYPE= "real, dimension(:,:)"> </IN>
  !<PUBLICROUTINE INTERFACE="write_field_data">
  subroutine write_field_data_2d(filename, fieldname, fielddata)
  !</PUBLICROUTINE>
    character(len=*),     intent(in) :: filename, fieldname
    real, dimension(:,:), intent(in) :: fielddata             
    integer                          :: n, k
    logical                          :: is_first_field, is_last_field

    if(.not. is_root_pe) return

    is_first_field = .false.
    is_last_field = .false.
    do n = 1, num_fields
       if(trim(flds_name(n)) == trim(fieldname)) then
          k = file_num(n)
          if(n ==1) then
             is_first_field = .true.
          else
             if(file_num(n) .ne. file_num(n-1)) then
                is_first_field = .true.
             endif
          endif
          if(n == max_fields) then
             is_last_field = .true.
          else
             if(file_num(n) .ne. file_num(n+1)) is_last_field = .true.
          endif

          ! if current field is the first field in a file, need to writ out axis data
          if(is_first_field ) then 
             call mpp_write(files_unit(k), axis_xt(k) )
             call mpp_write(files_unit(k), axis_yt(k) )
             call mpp_write(files_unit(k), axis_xc(k) )
             call mpp_write(files_unit(k), axis_yc(k) )
             call mpp_write(files_unit(k), axis_v(k) )
          endif

          !--- data should be already on global domain
          if(size(fielddata,1) .ne. ni .or. size(fielddata,2) .ne. nj ) then
             call mpp_error(FATAL,'grids_util_mod: data of field '//trim(fieldname)//' is not on global domain')
          endif

          call mpp_write(files_unit(k), flds(n), fielddata)
          if(is_last_field) call mpp_close(files_unit(k))

          return
       endif
    enddo

    !--- if fieldname is not found in the array flds_name, abort the program
    call mpp_error(FATAL,'grids_util_mod: file '//trim(filename)// &
         ' does not contain meta data of field '//trim(fieldname))

  end subroutine write_field_data_2d
  !  </SUBROUTINE>

  !#####################################################################
  !  <SUBROUTINE NAME="write_field_data_3d" INTERFACE="write_field_data">
  !   <IN NAME="fielddata" TYPE= "real, dimension(:,:,:)"> </IN>
  subroutine write_field_data_3d(filename, fieldname, fielddata)
    character(len=*),       intent(in) :: filename, fieldname
    real, dimension(:,:,:), intent(in) :: fielddata             
    integer                            :: n, k
    logical                            :: is_first_field, is_last_field

    if(.not. is_root_pe) return

    is_first_field = .false.
    is_last_field = .false.
    do n = 1, num_fields
       if(trim(flds_name(n)) == trim(fieldname)) then
             k = file_num(n)
          if(n ==1) then
             is_first_field = .true.
          else
             if(file_num(n) .ne. file_num(n-1)) then
                is_first_field = .true.
             endif
          endif 
        if(n == max_fields) then
           is_last_field = .true.
        else
           if(file_num(n) .ne. file_num(n+1)) is_last_field = .true.
        endif
           
          ! if current field is the first field in a file, need to writ out axis data
          if(is_first_field) then 
             call mpp_write(files_unit(k), axis_xt(k) )
             call mpp_write(files_unit(k), axis_xt(k) )
             call mpp_write(files_unit(k), axis_xc(k) )
             call mpp_write(files_unit(k), axis_yc(k) )
             call mpp_write(files_unit(k), axis_v(k) )
          endif
          call mpp_write(files_unit(k), flds(n), fielddata)
          if(is_last_field) call mpp_close(files_unit(k))

          return
        endif
     enddo

     !--- if fieldname is not found in the array flds_name, abort the program
     call mpp_error(FATAL,'grids_util_mod: file '//trim(filename)// &
          ' does not contain meta data of field '//trim(fieldname))

  end subroutine write_field_data_3d
  !  </SUBROUTINE>

  !#####################################################################
  ! <SUBROUTINE NAME="set_grid">

  !   <OVERVIEW>
  !    set the axis grid information.
  !   </OVERVIEW>     
  !   <TEMPLATE>
  !    call set_grid(grid_xt, grid_yt, grid_xc, grid_yc)
  !   </TEMPLATE>
  !   <IN NAME="grid_xt, grid_yt" TYPE="real, dimension(:)">
  !    longitude and latitude of the T-cell grid.
  !   </IN>
  !   <IN NAME="grid_xc, grid_yc" TYPE="real, dimension(:)">
  !    longitude and latitude of the C-cell grid.
  !   </IN>
  subroutine set_grid(grid_xt, grid_yt, grid_xc, grid_yc)
    real, dimension(:), intent(in) :: grid_xt, grid_yt, grid_xc, grid_yc


    ni = size(grid_xt(:))
    nj = size(grid_yt(:))

    allocate(xt(ni), yt(nj), xc(ni), yc(nj))
    xt = grid_xt
    yt = grid_yt
    xc = grid_xc
    yc = grid_yc 

    if(mpp_pe() == mpp_root_pe()) is_root_pe = .true.

  end subroutine set_grid
  !  </SUBROUTINE>
  !#####################################################################

end module grids_util_mod

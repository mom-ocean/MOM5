program regrid
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
  !
  ! <OVERVIEW>
  !   Since most analyses are placed on a spherical latitude-longitude grid and
  !   most global ocean models configured from mom4 are run with tripolar grids.
  !   A regridding tool is needed to regrid data from tripolar grid onto latitude-longitude 
  !   grid. This program map data from any logical rectangular grid 
  !   (tripolar or latitude-longitude) to a spherical latitude-longitude grid.
  !   When the data is on Tracer cell (T-cell), the interpolation will be a conservative 
  !   scheme. When the data is on other position (C-cell, E-cell or N-cell), regridding 
  !   is accomplished non-conservatively using a nearest neighbor distance weighting algorithm.
  ! </OVERVIEW>

  !<DESCRIPTION>
  ! Before using this regridding tool, you need to use preprocessing tool 
  ! make_xgrids to generate the grid file which contains source grid, destination 
  ! grid and exchange grid information between source grid and destination grid. 
  ! These exchange grid information is needed when doing conservative remapping.
  ! Suppose the file name of your source grid is src_grid.nc and the file name of
  ! your destination grid is dst_grid.nc (dst grid should be spherical 
  ! latitude-lontitude grid), the command will be 
  ! "make_xgrids -o src_grid.nc -a dst_grid.nc".
  ! Before using make_xgrids, you need to make sure src_grid.nc does not contain
  ! any exchange grid information. Otherwise you are not going to get the 
  ! desired exchange grid information. If your src_grid.nc do contains exchange
  ! grid information, you can remove all those fields by using ncks. Those fields are 
  ! AREA_ATMxOCN, DI_ATMxOCN, DJ_ATMxOCN, I_ATM_ATMxOCN, J_ATM_ATMxOCN, I_OCN_ATMxOCN, 
  ! J_OCN_ATMxOCN, AREA_ATMxLND, DI_ATMxLND, DJ_ATMxLND, I_ATM_ATMxLND, J_ATM_ATMxLND,
  ! I_LND_ATMxLND, J_LND_ATMxLND, AREA_LNDxOCN, DI_LNDxOCN, DJ_LNDxOCN, 
  ! I_LND_LNDxOCN, J_LND_LNDxOCN, I_OCN_LNDxOCN, J_OCN_LNDxOCN, xba, yba, xta, yta, 
  ! AREA_ATM, xbl, ybl, xtl, ytl, AREA_LND, AREA_LND_CELL, xto, yto and AREA_OCN. 
  ! After the grid file is generated, the grid file should be passed into the program 
  ! through nml option "grid_spec_file". 
  !
  ! This program expects to read data from a netcdf file, which is specfied
  ! by the namelist variable "src_data". The number of data to be remapped is 
  ! specified by num_flds. The name of field to be remapped is specified 
  ! by the namelist variable "fld_name". The output file is a netcdf file specified
  ! by the namelist variable "dst_data". Each field can be a scalar variable or
  ! a vector field, which is specified by vector_fld. The vector field should
  ! always be paired together. The data will be always on the source vertical grid.
  ! Previous experiences show that linear vertical interpolation will create lots of noises.
  ! If we find better vertical interpolation algorithm, we may implement it in the
  ! future. 
  !</DESCRIPTION>

  use mpp_mod,          only : mpp_npes, mpp_pe, mpp_root_pe, mpp_error, FATAL, NOTE, WARNING
  use mpp_mod,          only : stdout, stdlog, mpp_chksum
  use mpp_domains_mod,  only : mpp_define_layout, mpp_define_domains, domain2d
  use mpp_domains_mod,  only : mpp_global_field, mpp_get_compute_domain 
  use mpp_domains_mod,  only : mpp_domains_set_stack_size
  use fms_mod,          only : fms_init, open_namelist_file, close_file, fms_end
  use fms_mod,          only : check_nml_error, write_version_number, lowercase, uppercase
  use fms_io_mod,       only : fms_io_exit
  use constants_mod,    only : constants_init, PI, epsln
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type

  implicit none

#include "netcdf.inc"

  integer, parameter :: max_flds     = 20
  integer, parameter :: max_dims      = 20
  integer, parameter :: max_atts     = 20
  real,    parameter :: small        = 1.0e-6

  !--- namelist interface ----------------------------------------------
  !<NAMELIST NAME="regrid_nml">
  ! <DATA NAME="num_flds"  TYPE="integer" DEFAULT="0">
  !  Number of fields. 
  ! </DATA>
  ! <DATA NAME="src_data" TYPE="character(len=128)" >
  !  Name of input file containing to-be-remapped data.
  ! </DATA>
  ! <DATA NAME="dst_data" TYPE="character(len=128)" >
  !  Name of output file containing after-remapped data.
  ! </DATA>
  ! <DATA NAME="grid_spec_file" TYPE="character(len=128)" >
  !  Name of grid descriptor file containing source and target grid information, 
  !  exchange grid information between source and target grid. This grid file can 
  !  be created using preprocessing tool make_xgrids.
  ! </DATA>
  ! <DATA NAME="fld_name" TYPE="character(len=128), dimension(max_flds)" >
  !  Name of field to be regridded in input file.
  ! </DATA>
  ! <DATA NAME="fld_pos"  TYPE="character(len=1),dimension(max_flds)" DEFAULT="T">
  !  Name of grid where the field located.  Valid choices are (T)racer, (C)orner, (E)ast and (N)orth.
  ! </DATA>
  ! <DATA NAME="vector_field"  TYPE="logical,dimension(max_flds)" DEFAULT="False">
  !  True if fields are vector components. All the vector field should be paired together.
  !  That is, if vector_field(n) is .true., then vector_field(n+1) should be true.
  ! </DATA>
  ! <DATA NAME="num_nbrs" TYPE="integer">
  !  Number of nearest neighbors for regridding.
  ! </DATA>
  ! <DATA NAME="max_dist" TYPE="real" UNITS="radians">
  !  Maximum region of influence around destination grid points.
  ! </DATA>
  ! <DATA NAME="debug" TYPE="logical">
  ! For Debugging. Set true to print out chksum information for debugging reproducing ability 
  ! accross processors. default is false.
  ! </DATA>
  !</NAMELIST>

  integer            :: num_flds                  = 0
  character(len=128) :: src_data                  = 'INPUT/src_data.nc' 
  character(len=128) :: grid_spec_file            = 'INPUT/grid_spec.nc'
  character(len=128) :: dst_data                  = 'OUTPUT/dst_data.nc'
  character(len=16)  :: fld_name(max_flds)        = ''
  character(len=1)   :: fld_pos(max_flds)         = ''     ! with value 'T','C','E','N'
  logical            :: vector_fld(max_flds)      = .false. ! True if fields are vector components.
  integer            :: num_nbrs                  = 5
  real               :: max_dist                  = 0.05
  logical            :: debug                     = .false.

  namelist /regrid_nml/ num_flds, src_data, dst_data, grid_spec_file, fld_name, &
       fld_pos, vector_fld, num_nbrs, max_dist, debug

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: regrid.F90,v 14.0 2007/03/15 22:45:28 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  !--- variables about source data or on source grid
  integer                             :: ncid_src                     ! ncid corresponding to source file.
  integer                             :: ni_src, nj_src, nk_src       ! source grid size
  integer                             :: ntime_src                    ! number of time levels of src data
  logical                             :: src_is_cyclic                ! indicate if the source is cyclic
  logical                             :: src_is_tripolar              ! indicate if the source grid is tripolar
  real, dimension(max_flds)           :: missing_value =  -1e20       ! missing value for source field.
  logical                             :: has_taxis,tavg_exist         ! indicate if all the field has time axis
  logical, dimension(max_flds)        :: has_zaxis     = .false.      ! indicate if the field has z axis
  real, dimension(:),     allocatable :: time_value                   ! time axis value
  real, dimension(:),     allocatable :: t1_avg, t2_avg, dt_avg       ! time average data
  real, dimension(:,:),   allocatable :: area_ocn                     ! ocean area
  real, dimension(:),     allocatable :: zt_src                       ! source vertical grid
  real, dimension(:,:),   allocatable :: geolon_t, geolat_t           ! geographic grid location on T-cell
  real, dimension(:,:),   allocatable :: geolon_e, geolat_e           ! geographic grid location on E-cell
  real, dimension(:,:),   allocatable :: geolon_c, geolat_c           ! geographic grid location on C-cell
  real, dimension(:,:),   allocatable :: geolon_n, geolat_n           ! geographic grid location on N-cell
  real, dimension(:,:),   allocatable :: sinrot, cosrot               ! rotation

  !--- variables about destination data
  integer                             :: ni_dst, nj_dst                   ! destination grid size
  integer                             :: ncid_dst                         ! ncid corresponding to destination file.
  integer                             :: id_t1_dst, id_t2_dst, id_dt_dst
  integer                             :: id_time_dst, id_dst_fld(max_flds)
  character(len=128)                  :: dims_name(max_dims)              ! output file dimension name
  integer                             :: dims_id(max_dims)  
  character(len=1)                    :: dims_cart(max_dims)              ! 'X', 'Y', 'Z' or 'T'
  character(len=1)                    :: dims_pos(max_dims)               ! 'T', 'C', 'E' or 'N'   
  real, dimension(:), allocatable     :: xt_dst, yt_dst, xb_dst, yb_dst   !lat-lon destination grid
  type(domain2D), save                :: Domain                           ! destination grid domain.
  integer                             :: isc, iec, jsc, jec               ! compute domain decompsition of dst grid.
  integer                             :: ndims_dst                        ! number of dimensions in the output file

  !--- exchange grid information
  integer, dimension(:),    allocatable :: i_atm_atmxocn
  integer, dimension(:),    allocatable :: j_atm_atmxocn
  integer, dimension(:),    allocatable :: i_ocn_atmxocn
  integer, dimension(:),    allocatable :: j_ocn_atmxocn
  real, dimension(:),       allocatable :: area_atmxocn

  !--- other variables
  logical                             :: is_root_pe                    ! if current pe is root pe.
  type(horiz_interp_type),     target :: Interp_c, Interp_e, Interp_n  ! horiz_interp data at C,E,N-cell
  integer                             :: jstart                        ! starting index to do interpolation



  !---  begin of the program -------------------------------------------

  call regrid_init()

  call read_src_grid()

  call read_dst_grid()

  call setup_interp()

  call read_xgrid()

  call get_src_info()

  call setup_meta()

  call process_data ()

  call regrid_end ()

contains

  !#####################################################################
  subroutine regrid_init
    integer          :: unit, io_status, ierr, n
    character(len=1) :: direction

    call fms_init()
    call constants_init()

    !--- reading namelist ------------------------------------------------
    unit = open_namelist_file()
    read  (unit, regrid_nml,iostat=io_status)
    write (stdout(),'(/)')
    write (stdout(), regrid_nml)
    ierr = check_nml_error(io_status,'regrid_nml')
    call close_file (unit)

    !--- write out version information ---------------------------------
    call write_version_number(version, tagname)

    if(num_flds == 0) call error_handler('regrid: nml num_fiels = 0 should be a positive number')
    if(num_flds .gt. max_flds) call error_handler('regrid: nml num_fiels is greater than maximum'// &
         ' allowed fields max_flds. Modify num_flds or increase "max_flds" at top of the program' )

    do n = 1, num_flds
       if(fld_pos(n) .ne. 'T' .and. fld_pos(n) .ne. 'C' .and. fld_pos(n) .ne. 'E' .and. fld_pos(n) .ne. 'N' ) &
            call error_handler('regrid: fld_pos should be "T", "C", "E" or "N" ')
    enddo

    !--- when all the data are on T-cell, should always use 1 processor. Since no parallization is 
    !--- implemeneted.
    if(all(fld_pos(1:num_flds) == 'T') .and. mpp_npes() > 1) then
       write(stdout(),*)'WARNING: All the fields are on T-cell but npes is greater than 1, '//&
                        'should change processor count to 1'
    endif

    !--- number of vector_fld is true should be even and should be paired adjacent.
    !--- The first one of the pair should be u-component and the second one 
    !--- should be v-component.
    n = 1
    do while ( n .le. num_flds )
       if(vector_fld(n)) then
          if(fld_pos(n) .ne. 'C') call error_handler('regrid:When the field is a vector field, we assume '// &
             'this field is located at C-cell center. ')
          n = n + 1
          if( n .gt. num_flds .or. ( .NOT. vector_fld(n) ) ) then
             call error_handler('regrid: vector field should be paired together')
          endif
       endif
       n = n+1
    enddo

    !--- write out field information to be remapped.
    Write(stdout(),*) '*****************************************************************'
    write(stdout(),*) 'number of fields to be remapped are: ', num_flds
    write(stdout(),*) 'The source data file is: ', trim(src_data)
    write(stdout(),*) 'The grid spec file that contains exchange grid is: ', trim(grid_spec_file)
    Write(stdout(),*) 'The output data will be located on source vertical grid, '

    direction = 'x'
    do n = 1, num_flds
       if(vector_fld(n)) then
          write(stdout(),*) '** field ',n,': ',trim(fld_name(n)), ' at '//  &
               trim(fld_pos(n))//'-cell center is the '//direction//'-component of a vector field'         
          if(direction == 'x') then
             direction = 'y'
          else
             direction = 'x'
          endif
       else
          write(stdout(),*) '**field ',n,': ',trim(fld_name(n)), ' at '//  &
               trim(fld_pos(n))//'-cell center is a scalar variable'            
       endif
    enddo

    is_root_pe = mpp_pe() == mpp_root_pe()

  end subroutine regrid_init

  !#################################################################
  !--- read the src grid information from file grid_spec_file. The grid information
  !--- will include boundary condition, horizontal and vertical grid information. 
  !--- No need to read source grid land/sea mask, since we can always obtain source grid land/sea
  !--- mask from source data (the value of the point is equal to missing value will be 
  !--- land point). 

  subroutine read_src_grid
    integer              :: i
    integer              :: rcode, ncid, natt
    character(len=64)    :: name, catt
    real,    allocatable :: tmp(:,:)
    logical              :: is_new_grid
    real                 :: D2R

    D2R = PI/180.

    rcode = nf_open(trim(grid_spec_file), NF_NOWRITE, ncid)    
    call error_handler('error in open file '//trim(grid_spec_file), rcode)

    !--- get boundary condition ( cyclic or tripolar )
    !--- the boundary condition are specified by global attribute 
    !--- 'x_boundary_type' and 'y_boundary_type'
    rcode =  nf_inq_natts(ncid, natt)
    call error_handler('error in inquring ngatts of file '//trim(grid_spec_file), rcode )
    do i = 1, natt
       rcode = nf_inq_attname(ncid,NF_GLOBAL,i,name )
       call error_handler('error in inquring att name of file '//trim(grid_spec_file), rcode )
       catt = ''
       select case( trim(name) )
       case('x_boundary_type')
          rcode = nf_get_att_text(ncid,NF_GLOBAL,name,catt)
          call error_handler('error in inquring x_boundary_type value of file ' &
                              //trim(grid_spec_file), rcode )
          if(trim(catt) == 'cyclic') then
             src_is_cyclic = .true.
             Write(stdout(),*)' NOTE: x_boundary_type of grid '//trim(grid_spec_file)//' is cyclic'
             src_is_cyclic = .false.
          else
             Write(stdout(),*)' NOTE: x_boundary_type of '//trim(grid_spec_file)//' is solid_walls'
          endif
       case ('y_boundary_type')
          rcode = nf_get_att_text(ncid,NF_GLOBAL,name,catt)
          call error_handler('error in inquring y_boundary_type value of file ' &
                   //trim(grid_spec_file), rcode )
          if (trim(catt) == 'fold_north_edge') then
             src_is_tripolar = .true.
             Write(stdout(),*)' NOTE: y_boundary_type of '//trim(grid_spec_file)//' is tripolar'
          else
             src_is_tripolar = .false.
             Write(stdout(),*)' NOTE: y_boundary_type of '//trim(grid_spec_file)//' is solid_walls' 
          endif
       end select
    end do

    !--- To determine the grid is in new grid format or old grid format
    is_new_grid = check_is_new_grid(ncid)

    !--- get the vertical grid information
    nk_src = get_dimlen(ncid, 'zt')
    allocate(zt_src(nk_src) ) 
    call get_var_real_1d(ncid, 'zt', zt_src )

    if(is_new_grid) then
       ni_src = get_dimlen(ncid, 'grid_x_T')
       nj_src = get_dimlen(ncid, 'grid_y_T')
    else
       ni_src = get_dimlen(ncid, 'gridlon_t')
       nj_src = get_dimlen(ncid, 'gridlat_t')
    endif

    allocate(geolon_t(ni_src,nj_src), geolat_t(ni_src,nj_src) )
    allocate(geolon_c(ni_src,nj_src), geolat_c(ni_src,nj_src) )
    allocate(geolon_e(ni_src,nj_src), geolat_e(ni_src,nj_src) )
    allocate(geolon_n(ni_src,nj_src), geolat_n(ni_src,nj_src) )
    allocate(sinrot  (ni_src,nj_src), cosrot  (ni_src,nj_src) )  
    allocate(tmp     (ni_src,nj_src) )

    if(is_new_grid) then
       call get_var_real_2d(ncid, 'x_T', geolon_t)
       call get_var_real_2d(ncid, 'y_T', geolat_t)
       call get_var_real_2d(ncid, 'x_C', geolon_c)
       call get_var_real_2d(ncid, 'y_C', geolat_c)
       call get_var_real_2d(ncid, 'x_T', geolon_e)
       call get_var_real_2d(ncid, 'y_T', geolat_e)
       call get_var_real_2d(ncid, 'x_C', geolon_n)
       call get_var_real_2d(ncid, 'y_C', geolat_n)
       call get_var_real_2d(ncid, 'angle_C', tmp)
       sinrot = sin(tmp*D2R)
       cosrot = cos(tmp*D2R)
    else
       call get_var_real_2d(ncid, 'geolon_t', geolon_t)
       call get_var_real_2d(ncid, 'geolat_t', geolat_t)
       call get_var_real_2d(ncid, 'geolon_c', geolon_c)
       call get_var_real_2d(ncid, 'geolat_c', geolat_c)
       call get_var_real_2d(ncid, 'geolon_e', geolon_e)
       call get_var_real_2d(ncid, 'geolat_e', geolat_e)
       call get_var_real_2d(ncid, 'geolon_n', geolon_n)
       call get_var_real_2d(ncid, 'geolat_n', geolat_n)       
       call get_var_real_2d(ncid, 'sin_rot',  sinrot  )
       call get_var_real_2d(ncid, 'cos_rot',  cosrot  )
    endif       

    rcode = nf_close(ncid)

    return

  end subroutine read_src_grid


  !#################################################################
  !--- read the destination horizontal grid information from file "grid_spec_file".

  subroutine read_dst_grid
    integer              :: ncid, rcode

    !--- first reading grid information from file "grid_spec_file"
    rcode = nf_open(trim(grid_spec_file), NF_NOWRITE, ncid)    
    call error_handler('error in open file '//trim(grid_spec_file), rcode)

    !--- get the destination grid size
    ni_dst = get_dimlen(ncid, 'xta')
    nj_dst = get_dimlen(ncid, 'yta')   

    allocate(xt_dst(ni_dst), yt_dst(nj_dst), xb_dst(ni_dst+1), yb_dst(nj_dst+1) )
    call get_var_real_1d(ncid, 'xta', xt_dst)
    call get_var_real_1d(ncid, 'yta', yt_dst)  
    call get_var_real_1d(ncid, 'xba', xb_dst)
    call get_var_real_1d(ncid, 'yba', yb_dst)

    rcode = nf_close(ncid)

    !--- For the conservative interpolation, it is very complicated to implmeent parallization
    !--- also the conservative interpolation is efficiency enough for postprocessing purpose.
    !--- That's the reason we didn't implement parallization for the conservative scheme.


    return

  end subroutine read_dst_grid

  !#####################################################################
  !--- In many situation, the source grid match the destination grid in the
  !--- region southern of some latitude, indexed with jstart. 
  !--- Interpolation is necessary  only for the region northen of jstart.
  !--- this will save lots of interpolation time. When some field is located
  !--- on C, E, or N-cell, call horiz_interp_new to set up the interpolation.
  !--- in this case, parallization is implemented.
  subroutine setup_interp

    integer :: i, j, layout(2)
    real    :: D2R

    D2R = PI/180.

    !--- find the starting point to do interpolatioon jstart
    if(ni_src == ni_dst) then
       J_LOOP: do j = 1, nj_dst
          do i = 1, ni_dst
             if( abs( geolon_t(i,j) - xt_dst(i) )   .gt. small ) exit J_LOOP
             if( abs( geolat_t(i,j) - yt_dst(j) )   .gt. small ) exit J_LOOP
             if( abs( geolon_c(i,j) - xb_dst(i+1) ) .gt. small ) exit J_LOOP
             if( abs( geolat_c(i,j) - yb_dst(j+1) ) .gt. small ) exit J_LOOP
          enddo
       enddo J_LOOP
       jstart = j
    else
       jstart = 1
    endif

    if(jstart .ne. 1) write(stdout(),*)'NOTE: horizontal interpolation will be done ', &
                                       'starting at lat(',jstart,')= ',yt_dst(jstart)

    !--- set up domain for parallization if needed
    layout = (/1,0/)
    call mpp_define_layout((/1,ni_dst,jstart,nj_dst/), mpp_npes(), layout)
    call mpp_define_domains((/1,ni_dst,jstart,nj_dst/),layout, Domain )
    call mpp_get_compute_domain(Domain,isc,iec,jsc,jec) 

    !--- call horiz_interp_new to calculate the reampping weight in horizontal direction

    if( any(fld_pos == 'C') ) then
       call horiz_interp_new(Interp_c, geolon_c*D2R, geolat_c*D2R, xb_dst(isc+1:iec+1)*D2R, &
            yb_dst(jsc+1:jec+1)*D2R, interp_method="spherical",   &
            num_nbrs=num_nbrs, max_dist=max_dist, src_modulo = src_is_cyclic)
    endif

    if( any(fld_pos == 'E') ) then
       call horiz_interp_new(Interp_e, geolon_e*D2R, geolat_e*D2R, xb_dst(isc+1:iec+1)*D2R, &
            yt_dst(jsc:jec)*D2R, interp_method="spherical",   &
            num_nbrs=num_nbrs, max_dist=max_dist, src_modulo = src_is_cyclic)
    endif

    if( any(fld_pos == 'N') ) then
       call horiz_interp_new(Interp_n, geolon_n*D2R, geolat_n*D2R, xt_dst(isc:iec)*D2R, &
            yb_dst(jsc+1:jec+1)*D2R, interp_method="spherical",   &
            num_nbrs=num_nbrs, max_dist=max_dist, src_modulo = src_is_cyclic)
    endif

    return

  end subroutine setup_interp

  !#################################################################
  ! read needed exchange grid information from file grid_spec_file
  ! for the purpose of first-order conservative interpolation for tracer data
  ! If needed, we may adding the option to do second-order conservative
  ! interpolation in the future. The algorithm will be very similiar to
  ! the algorithm used in coupler flux exchange, which is implemented in
  ! shared code xgrid.f90.

  subroutine read_xgrid()
    integer              :: rcode, ncid, ni_lon, nj_lat, ncells, i, j

    !--- It is needed only when doing conservative interpolation, which means
    !--- there is a field located on T-cell.
    if(.not. any(fld_pos == 'T')) return


    !--- read xgrid from the grid_spec.nc file --------------------------
    rcode = nf_open(grid_spec_file, NF_NOWRITE, ncid)
    call error_handler('error in open file '//grid_spec_file, rcode)

    !--- get exchange grid information
    ncells = get_dimlen(ncid, 'I_ATM_ATMxOCN')
    allocate(i_atm_atmxocn(ncells), j_atm_atmxocn(ncells), &
             i_ocn_atmxocn(ncells), j_ocn_atmxocn(ncells), &
             area_atmxocn(ncells) )

    call get_var_int_1d (ncid, 'I_ATM_ATMxOCN', i_atm_atmxocn)
    call get_var_int_1d (ncid, 'J_ATM_ATMxOCN', j_atm_atmxocn )
    call get_var_int_1d (ncid, 'I_OCN_ATMxOCN', i_ocn_atmxocn )
    call get_var_int_1d (ncid, 'J_OCN_ATMxOCN', j_ocn_atmxocn )
    call get_var_real_1d(ncid, 'AREA_ATMxOCN',  area_atmxocn )

    rcode = nf_close(ncid)
    !--- this variable will be used for the conservative interpolation.
    allocate(area_ocn(ni_dst, jstart:nj_dst))

  end subroutine read_xgrid

  !#####################################################################
  !--- This routine will read source data information from file src_data.
  !--- This information includes missing_value for each field, if any field
  !--- contains attribute "time_avg_info". Get the time information. 
  !--- Since the metadata of destination data will be determined by the
  !--- source data, get the metadata information for the dst_data. 
  subroutine get_src_info
    integer            :: rcode, ndims, natt, natt_dim, id_fld, dimids(4), id_dim
    integer            :: num_has_taxis, i, j, m, n
    character(len=128) :: name, cart, att_name
    logical            :: found_dim    

    rcode = nf_open(trim(src_data), NF_NOWRITE, ncid_src) 
    call error_handler('error in opening file '//trim(src_data), rcode) 

    num_has_taxis  = 0
    ndims_dst      = 0
    tavg_exist     = .false.
    ntime_src      = 1

    do n = 1, num_flds
       rcode = nf_inq_varid(ncid_src, fld_name(n), id_fld)
       call error_handler('error in inquring id of field '//trim(fld_name(n))// &
                 ' of file '//trim(src_data), rcode)
       rcode = nf_inq_varndims( ncid_src, id_fld, ndims )
       call error_handler('error in inquring ndims of field '//trim(fld_name(n))// &
                 ' of file '//trim(src_data), rcode)
       rcode = nf_inq_vardimid( ncid_src, id_fld, dimids )
       call error_handler('error in inquring dimids of field '//trim(fld_name(n))// &
                 ' of file '//trim(src_data), rcode)
       rcode = nf_inq_varnatts(ncid_src, id_fld, natt )
       call error_handler('error in inquring natt of field '//trim(fld_name(n))// &
                 ' of file '//trim(src_data), rcode)

       !--- get the missing value and check if the field has tag_exist attribute.

       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_fld,i,name )          
          call error_handler('error in inquring attribute name', rcode)
          if(trim(name) == 'time_avg_info') then
             tavg_exist = .true.
          else if(trim(name) == 'missing_value') then
             rcode = nf_get_att_double(ncid_src, id_fld, name, missing_value(n) )
          endif
       enddo

       !--- to check if each field has z or time axis
       do m = 1, ndims
          rcode = nf_inq_dimname ( ncid_src, dimids(m), name )
          call error_handler('error in inquring dimension name of dimid', rcode)
          rcode = nf_inq_varid(ncid_src, name, id_dim)                    
          call error_handler('error in inquring variable id of '//trim(name), rcode)
          if(lowercase(trim(name)) == 'time') then
             cart = 'T'
             if(num_has_taxis == 0) then
                ntime_src = get_dimlen(ncid_src, trim(name) )
                allocate(time_value(ntime_src) )
                call get_var_real_1d(ncid_src, trim(name), time_value)
             endif
             num_has_taxis = num_has_taxis + 1   
          else
             rcode = nf_inq_varnatts(ncid_src, id_dim, natt_dim)
             call error_handler('error in inquring number of attributes of variable '//trim(name), rcode)
             do i = 1, natt_dim
                rcode = nf_inq_attname(ncid_src, id_dim,i,att_name )
                call error_handler('error in inquring attribute name ', rcode)
                if(trim(att_name) == 'cartesian_axis' ) then
                   rcode = nf_get_att_text(ncid_src,id_dim,att_name,cart)
                   call error_handler('error in get attribute value of attribute '//trim(att_name), rcode)
                   select case(cart(1:1))
                   case('X','Y')
                      ! do nothing 
                   case('Z') 
                      has_zaxis(n) = .true.
                   case default
                      call error_handler('The cartesian_axis of each axis should be X, Y or Z ' )
                   end select
                   exit
                endif
             enddo
          endif

          found_dim = .false.
          do j = 1, ndims_dst
             if(trim(dims_name(j)) == trim(name)) then
                found_dim = .true.
                exit
             endif
          enddo
          if(.not. found_dim) then
             ndims_dst = ndims_dst + 1
             dims_name(ndims_dst) = trim(name)
             dims_pos(ndims_dst) = fld_pos(n)
             dims_cart(ndims_dst) = cart(1:1)
             dims_id(ndims_dst)   = id_dim
          endif
       enddo
    enddo

    if(num_has_taxis == 0) then
       has_taxis = .false.
    else if(num_has_taxis == num_flds) then
       has_taxis = .true.
    else
       call error_handler('Either all or none of the fields has time axis ' )
    endif

    !--- suppose the time_avg_info is specified by average_T1, average_T2, average_DT.
    if(tavg_exist) then
       !--- some files don't have 'average_T1', 'average_T2' and 'average_DT' even though 
       !--- some variable do have time_avg_info attribute.  We suppose 'average_T1', 
       !--- 'average_T2' and 'average_DT' always co-exist.
       rcode = nf_inq_varid(ncid_src, 'average_T1', id_fld)
       if(rcode == 0) then
          allocate( t1_avg(ntime_src), t2_avg(ntime_src), dt_avg(ntime_src) )
          call get_var_real_1d(ncid_src, 'average_T1', t1_avg )
          call get_var_real_1d(ncid_src, 'average_T2', t2_avg )
          call get_var_real_1d(ncid_src, 'average_DT', dt_avg )
       else
          tavg_exist = .false.
       endif
    endif

  end subroutine get_src_info

  !#####################################################################
  !--- This routine will read source information.
  !--- set up axis and field meta data for destination data file
  !--- This routine is kind of long. We may change the design in the future
  !--- to have shorter subroutine for easy reading.

  subroutine setup_meta
       integer            :: rcode, natt, type, ndims, natt_dim, id_fld, dimids(4)
       integer            :: dim_time_dst, dims_dst(max_dims), id_axes(max_dims)
       integer            :: i, j, n, len, fld_dims(4), start(4), nwrite(4)
       character(len=128) :: name
       character(len=512) :: catt

    !--- simply return if current pe is not root pe.
    if(.not. is_root_pe) return

    rcode = nf_create(trim(dst_data),NF_WRITE, ncid_dst)
    call error_handler('error in creating file '//trim(dst_data), rcode) 

    !--- read the global attributes from source data and write it to output file
    rcode =  nf_inq_natts(ncid_src, natt )
    call error_handler('error in inquring natts of file '//trim(src_data), rcode)
    do i = 1, natt
       rcode = nf_inq_attname(ncid_src,NF_GLOBAL,i,name )
       call error_handler('error in inquring attribute name of file '//trim(src_data), rcode)
       rcode = nf_inq_att(ncid_src,NF_GLOBAL,name,type,len)
       call error_handler('error in inquring attribute of file '//trim(src_data), rcode)
       if(type == NF_CHAR ) then      ! only keep character global attributes
          if(len .le. 512) then
             if(trim(name) == 'filename') then
                catt = trim(dst_data)
             else
                rcode = nf_get_att_text(ncid_src,NF_GLOBAL,name,catt)
                call error_handler('error in getting attribute text of file '//trim(src_data), rcode)
             endif
             if(is_root_pe) then
                rcode = nf_put_att_text(ncid_dst, NF_GLOBAL, name, len, catt )
                call error_handler('error in putting  attribute to file '//trim(dst_data), rcode)
             endif
          else
             write(stdout(),*)'GLOBAL ATT '//trim(name)//' too long - not reading this metadata'
          endif
       endif
    enddo

    rcode = nf_inq_ndims(ncid_src, ndims)
    call error_handler('error in inquring ndims of file '//trim(src_data), rcode)
    do j = 1, ndims_dst
       if(lowercase(trim(dims_name(j)))=='time' ) then
          rcode = nf_def_dim(ncid_dst, dims_name(j), NF_UNLIMITED, dims_dst(j))
          call error_handler('error in defining time dimension of file '//trim(dst_data), rcode)
          rcode = nf_def_var(ncid_dst, dims_name(j), NF_DOUBLE, 1, dims_dst(j), id_axes(j))
          call error_handler('error in defining time variable of file '//trim(dst_data), rcode)
          id_time_dst = id_axes(j)
          dim_time_dst = dims_dst(j)
       else
          select case(dims_cart(j))
          case('X')
             len = ni_dst
          case('Y')
             len = nj_dst
          case('Z')
             len = nk_src             
          end select
          rcode = nf_def_dim(ncid_dst, dims_name(j), len, dims_dst(j))
          call error_handler('error in defining dimension '//trim(dims_name(j))//' of file '//trim(dst_data), rcode)
          rcode = nf_def_var(ncid_dst, dims_name(j), NF_DOUBLE, 1, dims_dst(j), id_axes(j) )
          call error_handler('error in defining axis var '//trim(dims_name(j))//' of file '//trim(dst_data), rcode)
       endif

       !--- write out axis attribute
       rcode = nf_inq_varnatts(ncid_src, dims_id(j), natt_dim)
       call error_handler('error in inquring number of attributes', rcode)
       do i = 1, natt_dim 
          rcode = nf_inq_attname(ncid_src, dims_id(j), i, name)             
          rcode = nf_copy_att(ncid_src,dims_id(j), name, ncid_dst, id_axes(j))
          call error_handler('error in copy attribute '//trim(name), rcode)            
       enddo
    enddo

    !--- read attribue of each field and write them to dst_data file
    do n = 1, num_flds
       rcode = nf_inq_varid(ncid_src,fld_name(n), id_fld)
       call error_handler('error in inquring id of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       rcode = nf_inq_varndims( ncid_src, id_fld, ndims )
       call error_handler('error in inquring ndims of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       rcode = nf_inq_vardimid( ncid_src, id_fld, dimids )
       call error_handler('error in inquring dimids of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       do i = 1, ndims
          rcode = nf_inq_dimname ( ncid_src, dimids(i), name )
          call error_handler('error in inquring dimension name of dimid', rcode)
          do j = 1, ndims_dst
             if(trim(name) == trim(dims_name(j))) then
                fld_dims(i) = dims_dst(j)
                exit
             endif
          enddo
       enddo
       call define_field_metadata(ncid_src, ncid_dst, fld_name(n), ndims, fld_dims, id_dst_fld(n))       
    enddo

    !--- define time avg info, suppose the time_avg_info is specified by average_T1, average_T2, average_DT.
    if(tavg_exist) then
       fld_dims(1) = dim_time_dst
       call define_field_metadata(ncid_src, ncid_dst, 'average_T1', 1, fld_dims, id_t1_dst)
       call define_field_metadata(ncid_src, ncid_dst, 'average_T2', 1, fld_dims, id_t2_dst)
       call define_field_metadata(ncid_src, ncid_dst, 'average_DT', 1, fld_dims, id_dt_dst)
    endif

    rcode = nf_enddef(ncid_dst)

    !--- write axis data to dst_data
    start = 1; nwrite = 1
    do j = 1, ndims_dst
       select case(dims_cart(j))
       case('X')
          nwrite(1) = ni_dst
          if(dims_pos(j) == 'T' .or. dims_pos(j) == 'N') then
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, xt_dst )
          else
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, xb_dst(2:) )             
          endif
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       case('Y')
          nwrite(1) = nj_dst
          if(dims_pos(j) == 'T' .or. dims_pos(j) == 'E') then
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, yt_dst )
          else
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, yb_dst(2:) )             
          endif
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       case('Z')
          nwrite(1) = nk_src
          rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, zt_src )
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       end select
    enddo

  end subroutine setup_meta

  !#####################################################################
    !-------------------------------------------------------------------
    !--- read src data, remap onto current grid and write output -------
  subroutine process_data
    integer :: n, m, l, rcode, nk, nf 
    integer :: start(4), nread(4), nwrite(4), id_fld(max_flds)
    real, dimension(:,:,:,:), allocatable :: data_src, data_dst
    real, dimension(:,:,:),   allocatable :: wrk1, wrk2

    do n = 1, num_flds
       rcode = nf_inq_varid(ncid_src, fld_name(n), id_fld(n) )
       call error_handler('error in inquring varid of '//trim(fld_name(n))// &
                          ' of file '//trim(src_data), rcode)
    enddo

    write(stdout(),*)'***************************************'
    write(stdout(),*)' ntime_src is ', ntime_src

    do m = 1, ntime_src
       write(stdout(),*)'************************************************'
       write(stdout(),*)'*** At time level ', m
       n = 1
       !--- write out time information 
       if(is_root_pe) then
          if(has_taxis) then       
             rcode = nf_put_var1_double(ncid_dst, id_time_dst, m, time_value(m))
             call error_handler('error in put time data to file '//trim(src_data), rcode )
             if(tavg_exist) then
                rcode = nf_put_var1_double(ncid_dst, id_t1_dst, m, t1_avg(m))
                call error_handler('error in put average_T1 data to file '//trim(src_data), rcode )
                rcode = nf_put_var1_double(ncid_dst, id_t2_dst, m, t2_avg(m))
                call error_handler('error in put average_T2 data to file '//trim(src_data), rcode )
                rcode = nf_put_var1_double(ncid_dst, id_dt_dst, m, dt_avg(m))
                call error_handler('error in put average_DT data to file '//trim(src_data), rcode )
             endif
          endif
       endif

       do while ( n .le. num_flds )
          nf = 1
          if(vector_fld(n) .and. src_is_tripolar) nf = 2
          nk = 1
          if(has_zaxis(n)) nk = nk_src

          allocate(data_src(ni_src, nj_src, nk,nf), data_dst(ni_dst, nj_dst, nk,nf) )
          allocate(wrk1(isc:iec,jsc:jec,nk), wrk2(1:ni_dst,jstart:nj_dst,nk) )

          start = 1; nread = 1
          nread(1) = ni_src; nread(2) = nj_src
          if(has_zaxis(n)) then
             nread(3) = nk; start(4) = m
          else
             start(3) = m
          endif
          rcode = nf_get_vara_double(ncid_src,id_fld(n), start, nread, data_src(:,:,:,1) ) 
          call error_handler('error in get data 1 from file '//trim(src_data), rcode )

          !--- if the field is a vector field and src_grid is tripolar, read the next field and 
          !--- rotate data onto spherical grid.
          if( vector_fld(n) .and. src_is_tripolar ) then
             rcode = nf_get_vara_double(ncid_src,id_fld(n+1), start, nread, data_src(:,:,:,2) ) 
             call error_handler('error in get data 2 from file '//trim(src_data), rcode )

             call rotate_data(data_src(:,jstart:,:,:), cosrot(:,jstart:), sinrot(:,jstart:), missing_value(n) )
          endif

          do l = 1, nf
             !---- doing horizontal interpolaiton
             if(fld_pos(n) == 'T' ) then !--- this case data is already global data
                call conserve_interp(data_src(:,:,:,1), wrk2(:,:,:), missing_value(1))
             else
                call mpp_domains_set_stack_size(2*ni_dst*nj_dst*nk)

                call nonconserve_interp(data_src(:,:,:,l), wrk1(:,:,:), fld_pos(n), missing_value(1) )
                !---- get the global data
                call mpp_global_field(Domain, wrk1, wrk2 )
             endif

             !--- assign data to global region
             if(jstart == 1) then
                data_dst(:,:,:,l) = wrk2
             else
                data_dst(:,1:jstart-1,:,l) = data_src(:,1:jstart-1,:,l)
                data_dst(:,jstart:,:,l)    = wrk2(:,:,:)
             endif
          enddo

          !--- write out data from root_pe
          if(is_root_pe) then
             do l = 1, nf
                start = 1; nwrite = 1
                nwrite(1) = ni_dst; nwrite(2) = nj_dst
                if(has_zaxis(n)) then
                   nwrite(3) = nk; start(4) = m
                   rcode = nf_put_vara_double(ncid_dst, id_dst_fld(n+l-1), start, nwrite, data_dst(:,:,:,l) ) 
                   call error_handler('error in putting '//trim(fld_name(n+l-1))//' in file '//trim(dst_data), rcode )
                else
                   start(3) = m
                   rcode = nf_put_vara_double(ncid_dst, id_dst_fld(n+l-1), start, nwrite, data_dst(:,:,1,l) ) 
                   call error_handler('error in putting '//trim(fld_name(n+l-1))//' in file '//trim(dst_data), rcode )
                endif
                !--- for debugging reproducing ability accross processor count.
                if(debug) write(stdout(),*)'NOTE:after regrid, chksum for '//trim(fld_name(n+l-1))// &
                              ' at time level ',m,' ====> ',mpp_chksum(data_dst(:,:,:,l), (/mpp_root_pe()/) )
             enddo
          endif
          n = n+nf
          deallocate(data_src, data_dst, wrk1, wrk2)
       enddo
    enddo

    rcode = nf_close(ncid_src)
    if(is_root_pe) rcode = nf_close(ncid_dst)

  end subroutine process_data

  !#####################################################################

  subroutine regrid_end


    call fms_io_exit
    call fms_end()

  end subroutine regrid_end

  !#####################################################################
  ! error handling routine.
  subroutine error_handler(mesg, status)
    character(len=*),  intent(in) :: mesg
    integer, optional, intent(in) :: status
    character(len=256) :: msg

    if(present(status)) then
       if(status == 0) return
       msg = nf_strerror(status)
       msg = trim(mesg)//': '// trim(msg)
    else
       msg = trim(mesg)
    endif

    call mpp_error(FATAL,trim(msg) )

  end subroutine error_handler


  !#####################################################################

  ! get the dimension length of any one dimensional variable
  function get_dimlen(ncid, name)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer                      :: get_dimlen
    integer                      :: varid, rcode, dims(1)

    get_dimlen = 0

    rcode = nf_inq_varid(ncid, trim(name), varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_inq_vardimid(ncid, varid, dims)
    call error_handler('error in inquiring dimension id of '//trim(name), rcode)

    rcode = nf_inq_dimlen(ncid, dims(1), get_dimlen)
    call error_handler('error in inquiring dimension length of '//trim(name), rcode)

  end function get_dimlen

  !#####################################################################
  ! read the 1d integer data from netcdf file.
  subroutine get_var_int_1d(ncid, name, data)
    integer,                intent(in) :: ncid
    character(len=*),       intent(in) :: name
    integer, dimension(:), intent(out) :: data
    integer                            :: rcode, varid

    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)
    rcode = nf_get_var_int(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)

  end subroutine get_var_int_1d

  !#####################################################################
  ! read the 1d real data from netcdf file.
  subroutine get_var_real_1d(ncid, name, data)
    integer,             intent(in) :: ncid
    character(len=*),    intent(in) :: name
    real, dimension(:), intent(out) :: data
    integer                         :: rcode, varid

    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_get_var_double(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)


  end subroutine get_var_real_1d

  !#####################################################################
  ! read the 2d real data from netcdf file.
  subroutine get_var_real_2d(ncid, name, data)
    integer,               intent(in) :: ncid
    character(len=*),      intent(in) :: name
    real, dimension(:,:), intent(out) :: data
    integer                           :: rcode, varid


    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_get_var_double(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)

  end subroutine get_var_real_2d


  !#####################################################################

  !--- rotate the data when the field is a vector field.
  subroutine rotate_data(data, rotcos, rotsin, missing)
    real, dimension(:,:,:,:), intent(inout) :: data
    real, dimension(:,:),     intent(in)    :: rotsin, rotcos
    real,                        intent(in) :: missing

    integer                       :: i, j, k
    real                          :: temp1, temp2

    if(size(data,4).ne. 2 ) call error_handler('number of data should be paired when do rotation. ')

    do k = 1, size(data,3)
       do j = 1, size(data,2)
          do i = 1, size(data,1)
             if(data(i,j,k,1) /= missing ) then ! assume data(:,:,:,1) is not missing, then
                ! data(:,:,:,2) is not missing
                temp1           = data(i,j,k,1)*rotcos(i,j)-data(i,j,k,2)*rotsin(i,j)
                temp2           = data(i,j,k,1)*rotsin(i,j)+data(i,j,k,2)*rotcos(i,j)
                data(i,j,k,1)   = temp1
                data(i,j,k,2)   = temp2
             endif
          enddo
       enddo
    enddo

  end subroutine rotate_data

  !#####################################################################
  subroutine conserve_interp(input, output, missing)
    real, dimension(:,:,:),        intent(in) :: input
    real, dimension(:,jstart:,:), intent(out) :: output
    real,                          intent(in) :: missing

    integer                                   :: i, j, k

    !--- put data on exchange grid
    output = 0.0
    area_ocn  = 0.0

    do k = 1, size(input,3)
       area_ocn = 0.0
       do i = 1, size(i_atm_atmxocn)
          if(j_atm_atmxocn(i) >= jstart ) then
             if(input(i_ocn_atmxocn(i),j_ocn_atmxocn(i),k) /= missing) then
                output(i_atm_atmxocn(i),j_atm_atmxocn(i),k) = output(i_atm_atmxocn(i),j_atm_atmxocn(i),k) &
                     + area_atmxocn(i)*input(i_ocn_atmxocn(i),j_ocn_atmxocn(i),k)
                area_ocn(i_atm_atmxocn(i),j_atm_atmxocn(i)) = area_ocn(i_atm_atmxocn(i),j_atm_atmxocn(i)) &
                     + area_atmxocn(i)
             endif
          endif
       enddo

       do j = jstart, nj_dst
          do i = 1, ni_dst
             if(area_ocn(i,j) > epsln ) then
                output(i,j,k) = output(i,j,k)/area_ocn(i,j)
             else
                output(i,j,k) = missing
             endif
          enddo
       enddo
    enddo

  end subroutine conserve_interp

  !#####################################################################
    subroutine nonconserve_interp(input, output, position, missing)
    real, dimension(:,:,:),  intent(in) :: input
    real, dimension(:,:,:), intent(out) :: output
    character(len=1),        intent(in) :: position
    real,                    intent(in) :: missing

    integer                                       :: i, j, k
    real, dimension(size(input,1), size(input,2)) :: mask
    type(horiz_interp_type), pointer              :: Interp=>NULL()

    select case(position)
    case('C')
       Interp => Interp_c
    case('E')
       Interp => Interp_e
    case('N')
       Interp => Interp_n
    end select

    do k = 1, size(input,3)
       !--- define land/sea mask
       do j = 1, size(input,2)
          do i = 1, size(input,1)
             if(input(i,j,k) == missing) then
                mask(i,j) = 0.0
             else
                mask(i,j) = 1.0
             endif
          enddo
       enddo
       !--- remap data onto destination horizontal grid through horiz_interp
       call horiz_interp(Interp, input(:,:,k), output(:,:,k), &
            mask_in=mask, missing_value=missing )
    enddo

    return

  end subroutine nonconserve_interp

  !#####################################################################
  !--- To determine the grid is in new grid format or old grid format
  function check_is_new_grid(ncid)
    integer, intent(in) :: ncid
    logical             :: check_is_new_grid
    integer             :: rcode, id_fld

    rcode = nf_inq_varid(ncid, 'x_T', id_fld)
    if(rcode == 0) then
       check_is_new_grid = .true.
    else
       rcode = nf_inq_varid(ncid, 'geolon_t', id_fld)
       call error_handler('check_is_new_grid: error in inquring id of geolon_t ', rcode)       
       check_is_new_grid = .false.
    endif

  end function check_is_new_grid

  !#####################################################################
  !--- define a field metadata and copy the attribute from source data to destination data
  subroutine define_field_metadata(ncid_in, ncid_out, name, ndim, dims, id_fld_out)
    integer, intent(in)          :: ncid_in, ncid_out, ndim, dims(:)
    character(len=*), intent(in) :: name
    integer, intent(out)         :: id_fld_out

    integer            :: rcode, id_fld_in, natt, i
    character(len=128) :: attname

    rcode = nf_def_var(ncid_out, trim(name), NF_DOUBLE, ndim, dims(1:ndim), id_fld_out)
    call error_handler('define_field_metadata: error in defining field '//trim(name), rcode )
    rcode = nf_inq_varid(ncid_in, trim(name), id_fld_in)
    call error_handler('define_field_metadata: error in inquring field of '//trim(name), rcode) 
    rcode = nf_inq_varnatts(ncid_in, id_fld_in, natt)
    call error_handler('define_field_metadata: error in inquring natt of '//trim(name), rcode)       

    do i = 1, natt
       rcode = nf_inq_attname(ncid_in, id_fld_in, i, attname)
       call error_handler('error in inquring '//trim(name)//' attribute '//trim(attname), rcode) 
       if(trim(attname) /= '_FillValue') then 
          rcode = nf_copy_att(ncid_in,id_fld_in, attname, ncid_out, id_fld_out)
          call error_handler('error in copy '//trim(name)//' attribute '//trim(attname), rcode) 
       endif
    enddo

  end subroutine define_field_metadata

  !#####################################################################
  
end program regrid

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
  !   This program remap data from logically rectangular grid to logically rectangular grid.
  ! </OVERVIEW>

  !<DESCRIPTION>
  ! This program expects to read data from a netcdf file, which is specfied
  ! by the namelist variable "src_data". The number of data to be remapped is 
  ! specified by num_flds. The name of field to be remapped is specified 
  ! by the namelist variable "fld_name". The source data should be on the source
  ! grid which is specified by namelist variable src_grid. The destination grid is
  ! specified by nml dst_grd. The output file is a netcdf file specified
  ! by the namelist variable "dst_data". Each field can be a scalar variable or
  ! a vector field, which is specified by vector_fld. The vector field should
  ! always be paired together. A laplacian extrapolation will be performed when
  ! there is any missing value in the source data to interpolate data onto missing points.
  !</DESCRIPTION>

  use mpp_mod,          only : mpp_npes, mpp_pe, mpp_root_pe, mpp_error
  use mpp_mod,          only : FATAL, WARNING, stdout, mpp_chksum
  use mpp_domains_mod,  only : mpp_define_layout, mpp_define_domains, domain2d
  use mpp_domains_mod,  only : mpp_global_field, mpp_domains_set_stack_size, mpp_get_compute_domain 
  use mpp_domains_mod,  only : mpp_update_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
  use fms_mod,          only : fms_init, open_namelist_file, close_file, fms_end
  use fms_mod,          only : check_nml_error, write_version_number, lowercase, uppercase
  use fms_io_mod,       only : fms_io_exit
  use constants_mod,    only : constants_init, PI
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type
  use axis_utils_mod,   only : interp_1d

  implicit none

#include "netcdf.inc"

  integer, parameter :: max_flds              = 20
  integer, parameter :: max_iter              = 2000
  integer, parameter :: max_atts              = 20
  real,    parameter :: rel_coef              = 0.6
  real,    parameter :: epsln                 = 1.0e-10

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
  ! <DATA NAME="dst_grid" TYPE="character(len=128)" >
  !  Name of grid descriptor file containing target grid information.
  ! </DATA>
  ! <DATA NAME="src_grid" TYPE="character(len=128)" >
  !  Name of grid descriptor file containing source grid information.
  ! </DATA>
  ! <DATA NAME="fld_name" TYPE="character(len=128), dimension(max_flds)" >
  !  Name of runoff field in input file.
  ! </DATA>
  ! <DATA NAME="fld_pos"  TYPE="character(len=1),dimension(max_flds)" DEFAULT="T">
  !  Name of grid where the field located.  Valid choices are (T)racer, (C)orner, (E)ast and (N)orth.
  ! </DATA>
  ! <DATA NAME="vector_field"  TYPE="logical,dimension(max_flds)" DEFAULT="False">
  !  True if fields are vector components. All the vector field should be paired together.
  !  That is, if vector_field(n) is .true., then vector_field(n+1) should be true.
  ! </DATA>
  ! <DATA NAME="stop_crit"  TYPE="character(len=1),dimension(max_flds)" DEFAULT="0.001">
  !  The stopping criteria when extrapping data onto missing points.
  ! </DATA>
  ! <DATA NAME="num_nbrs" TYPE="integer">
  !  Number of nearest neighbors for regridding.
  ! </DATA>
  ! <DATA NAME="max_dist" TYPE="real" UNITS="radians">
  !  Maximum region of influence around destination grid points.
  ! </DATA>
  ! <DATA NAME="use_source_vertical_grid" TYPE="logical" DEFAULT=".false.">
  !  when use_source_vertical_grid is set to true, the destination data will 
  !  have the same vertical level as the source data. When use_source_vertical_grid 
  !  is false, the vertical grid of destination data will come from dest_grid. 
  !  a linear vertical interpolation will be done when the source vertical is different
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
  !<DATA NAME="interp_method"  TYPE= "character(len=20)" >
  ! specifying the remapping method when remampping data onto current grid.
  ! Its value can be "spherical" or " bilinear". "spherical" interpolation is a 
  ! inverse distance weighted interpolation algorithm. Default value is "bilinear". 
  ! "bilinear" interpolation is recommanded, since bilinear interpolation will provide 
  ! more smooth results than "spherical" interpolation (especially when interpolating 
  ! from coarse grid to fine grid). Plus bilinear interpolation is much more efficiency 
  ! than "spherical interpolation". Since bilinear interpolation suppose the source grid 
  ! is a lat-lon grid, "spherical" need to be used if the source grid is not a lat-lon grid.
  ! </DATA>
  ! <DATA NAME="debug" TYPE="logical">
  ! For Debugging. Set true to print out chksum information for debugging reproducing ability 
  ! accross processors. default is false.
  ! </DATA>
  !</NAMELIST>

  integer            :: num_flds                  = 0
  character(len=128) :: src_data                  = '' 
  character(len=128) :: src_grid                  = ''
  character(len=128) :: dst_data                  = ''
  character(len=128) :: dst_grid                  = ''
  character(len=16)  :: fld_name(max_flds)        = ''
  character(len=1)   :: fld_pos(max_flds)         = ''      ! with value 'T','C','E','N'
  logical            :: vector_fld(max_flds)      = .false. ! True if fields are vector components.
  real               :: stop_crit(max_flds)       = 0.001 
  integer            :: num_nbrs                  = 5
  real               :: max_dist                  = 0.1
  logical            :: use_source_vertical_grid  = .FALSE.
  logical            :: apply_mask                = .TRUE.
  character(len=32)  :: interp_method             = "bilinear"
  logical            :: debug                     = .false.

  namelist /regrid_nml/ num_flds, src_data, src_grid, dst_data, dst_grid, fld_name, &
       fld_pos, vector_fld, num_nbrs, max_dist, stop_crit, interp_method,           &
       apply_mask, use_source_vertical_grid, debug


  !--- data type for easy data management ------------------------------
  type cell_type
     real, dimension(:,:),   pointer :: geolon=>NULL(), geolat=>NULL() ! geographical grid
     real, dimension(:,:,:), pointer :: mask  =>NULL()                 ! mask at each vertical level
     real, dimension(:,:),   pointer :: sinrot=>NULL(), cosrot=>NULL() ! rotation
  end type cell_type

  type regrid_type
     type(cell_type)          :: T,C,E,N            ! T,C,E,N-cell
     integer                  :: ni, nj,nk          ! grid size
     integer                  :: isc, iec, jsc, jec ! compute domain decompsition src grid
     integer                  :: isd, ied, jsd, jed ! data domain decompsition on dst grid
     type(domain2d)           :: domain             ! domain of the grid ( only destination grid will use it).
     logical                  :: is_cyclic          ! indicate if the grid is cyclic
     logical                  :: is_tripolar        ! indicate if the grid is tripolar
     real, pointer            :: zt(:) =>NULL()     ! vertical grid 
     real, pointer            :: grid_xt(:)=>NULL() ! T-cell longitude grid
     real, pointer            :: grid_yt(:)=>NULL() ! T-cell latitude grid
     real, pointer            :: grid_xc(:)=>NULL() ! C-cell longitude grid
     real, pointer            :: grid_yc(:)=>NULL() ! C-cell latitude grid
  end type regrid_type

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: regrid.F90,v 20.0 2013/12/14 00:31:05 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  !--- other variables
  integer                             :: ntime_src                     ! number of time levels of src data
  logical                             :: do_vertical_interp            ! indicate if vertical interpolation is needed.
  logical                             :: is_root_pe                    ! if current pe is root pe.
  type(regrid_type)                   :: Src                           ! regrid_type data of source
  type(regrid_type)                   :: Dst                           ! regrid_type data of destination
  type(horiz_interp_type),     target :: Interp_T, Interp_C            ! horiz_interp data at T,C-cell
  type(horiz_interp_type),     target :: Interp_E, Interp_N            ! horiz_interp data at E,N-cell
  real, dimension(:),     allocatable :: missing_value                 ! missing value for source field.
  real, dimension(:),     allocatable :: time_value                    ! time axis value
  real, dimension(:),     allocatable :: t1_avg, t2_avg, dt_avg        ! time average data
  real, dimension(:,:),   allocatable :: ht, hc, he, hn                ! topography
  logical, dimension(:),  allocatable :: has_zaxis                     ! indicate if the field has z axis
  logical                             :: has_taxis,tavg_exist          ! indicate if all the field has time axis
  integer                             :: ncid_dst                      ! ncid corresponding to 
  integer                             :: id_t1_dst, id_t2_dst, id_dt_dst
  integer                             :: id_time_dst, id_dst_fld(max_flds)

  !---  begin of the program -------------------------------------------

  call regrid_init()

  call read_grid()

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
    call write_version_number( version, tagname )

    if(num_flds == 0) call error_handler('regrid: nml num_fiels = 0 should be a positive number')
    if(num_flds .gt. max_flds) call error_handler('regrid: nml num_fiels is greater than maximum'// &
         ' allowed fields max_flds. Modify num_flds or increase "max_flds" at top of the program' )

    do n = 1, num_flds
       if(fld_pos(n) .ne. 'T' .and. fld_pos(n) .ne. 'C' .and. fld_pos(n) .ne. 'E' .and. fld_pos(n) .ne. 'N' ) &
            call error_handler('regrid: fld_pos should be "T", "C", "E" or "N" ')
    enddo

    !--- number of vector_fld is true should be even and should be paired adjacent.
    !--- The first one of the pair should be u-component and the second one 
    !--- should be v-component.
    n = 1
    do while ( n .le. num_flds )
       if(vector_fld(n)) then
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
    write(stdout(),*) 'The source grid file is: ', trim(src_grid)
    write(stdout(),*) 'The destination grid file is: ', trim(dst_grid)
    if(use_source_vertical_grid) then
       write(stdout(),*) 'use_source_vertical_grid is set to true. The destination data will have '// &
                         'same vertical level as the source data '
    else
       write(stdout(),*) 'use_source_vertical_grid is set to false. vertical interpolation will ', &
                         'be performed if the source and destination has different vertical grid. '
    endif

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

    !--- memory allocation
    allocate(has_zaxis(num_flds), missing_value(num_flds) )
    has_zaxis = .false.
    missing_value = -1e20

    is_root_pe = mpp_pe() == mpp_root_pe()

  end subroutine regrid_init

  !#################################################################
  !--- read the src grid and dst grid and land/sea mask of dst grid.

  subroutine read_grid
    integer  :: isc, iec, jsc, jec
    integer  :: k, npes, layout(2)
    real     :: D2R
    logical  :: use_src_vgrid

    D2R = PI/180.0

    !--- read some general grid informaiton for source grid and destination grid
    call get_grid_info(src_grid, Src, .TRUE.)
    call get_grid_info(dst_grid, Dst, .FALSE.)

    do_vertical_interp = .false.
    if(use_source_vertical_grid) then !--- modify the destination vertical information for the future use.
       !--- In this case the destination vertical grid should be the same as source.
       Dst%nk = Src%nk
       deallocate(Dst%zt)
       allocate(Dst%zt(Dst%nk))
       Dst%zt = Src%zt
    else     !--- when use_source_vertical_grid is false, vertical interpolation will be done if the source vertical
       !--- grid didn't match destination vertical grid.
       if(Src%nk .ne. Dst%nk) then
          do_vertical_interp = .true.
       else
          do k = 1, Src%nk
             if( abs(Src%zt(k) - Dst%zt(k)) .gt. epsln ) then
                do_vertical_interp = .true.
                exit
             endif
          enddo
       endif
       if(do_vertical_interp) then
          write(stdout(),*)'NOTE: Since source vertical grid does not match ' //&
               'destination vertical grid and nml use_source_vertical_grid is set ' //&
               'to false, vertical interpolation will be performed. '
       else
          write(stdout(),*)'NOTE: Even though nml use_source_vertical_grid is set to false, ', &
               'but since source vertical grid match destination vertical grid, ',             &
               'no vertical interpolation will be performed. '
       endif
    endif

    isc= Dst%isc; iec= Dst%iec; jsc= Dst%jsc; jec= Dst%jec;

    !--- when source is tripolar grid, 'spherical' interpolation will be used 
    !--- no matter your choice of "interp_method".
    if(Src%is_tripolar .and. trim(interp_method) == "bilinear") then
       write(stdout(),*) "WARNING: Since src grid is tripolar grid, even though you chose ", &
             "bilinear option, spherical option will be used."
    endif

    !--- call horizinterp_new to calculate the reampping weight in horizontal direction
    if( any(fld_pos == 'T') ) then
       if(Src%is_tripolar) then
          call horiz_interp_new(Interp_T, Src%T%geolon*D2R, Src%T%geolat*D2R,           &
               Dst%T%geolon(isc:iec,jsc:jec)*D2R, Dst%T%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method='spherical', num_nbrs=num_nbrs, max_dist=max_dist,          &
               src_modulo = Src%is_cyclic)
       else
          call horiz_interp_new(Interp_T, Src%T%geolon(:,1)*D2R, Src%T%geolat(1,:)*D2R, &
               Dst%T%geolon(isc:iec,jsc:jec)*D2R, Dst%T%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method=trim(interp_method), num_nbrs=num_nbrs, max_dist=max_dist,  &
               src_modulo = Src%is_cyclic, grid_at_center = .true.)
       endif
    endif

    if( any(fld_pos == 'C') ) then
       if(Src%is_tripolar) then
          call horiz_interp_new(Interp_C, Src%C%geolon*D2R, Src%C%geolat*D2R,           &
               Dst%C%geolon(isc:iec,jsc:jec)*D2R, Dst%C%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method='spherical', num_nbrs=num_nbrs, max_dist=max_dist,          &
               src_modulo = Src%is_cyclic )
       else
          call horiz_interp_new(Interp_C, Src%C%geolon(:,1)*D2R, Src%C%geolat(1,:)*D2R, &
               Dst%C%geolon(isc:iec,jsc:jec)*D2R, Dst%C%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method=trim(interp_method), num_nbrs=num_nbrs, max_dist=max_dist,  &
               src_modulo = Src%is_cyclic, grid_at_center = .true.)
       endif
    endif

    if( any(fld_pos == 'E') ) then
       if(Src%is_tripolar) then
          call horiz_interp_new(Interp_E, Src%E%geolon*D2R, Src%E%geolat*D2R,           &
               Dst%E%geolon(isc:iec,jsc:jec)*D2R, Dst%E%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method='spherical', num_nbrs=num_nbrs, max_dist=max_dist,          &
               src_modulo = Src%is_cyclic)
       else
          call horiz_interp_new(Interp_E, Src%E%geolon(:,1)*D2R, Src%E%geolat(1,:)*D2R, &
               Dst%E%geolon(isc:iec,jsc:jec)*D2R, Dst%E%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method=trim(interp_method), num_nbrs=num_nbrs, max_dist=max_dist,  &
               src_modulo = Src%is_cyclic, grid_at_center = .true.)
       endif
    endif

    if( any(fld_pos == 'N') ) then
       if(Src%is_tripolar) then
          call horiz_interp_new(Interp_N, Src%N%geolon*D2R, Src%N%geolat*D2R,           &
               Dst%N%geolon(isc:iec,jsc:jec)*D2R, Dst%N%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method='spherical', num_nbrs=num_nbrs, max_dist=max_dist,          &
               src_modulo = Src%is_cyclic)
       else
          call horiz_interp_new(Interp_N, Src%N%geolon(:,1)*D2R, Src%N%geolat(1,:)*D2R, &
               Dst%N%geolon(isc:iec,jsc:jec)*D2R, Dst%N%geolat(isc:iec,jsc:jec)*D2R,     &
               interp_method=trim(interp_method), num_nbrs=num_nbrs, max_dist=max_dist,  &
               src_modulo = Src%is_cyclic, grid_at_center = .true.)
       endif
    endif

  end subroutine read_grid

  !#####################################################################

  subroutine get_grid_info(gridfile, Regrid, is_src_grid)
    character(len=*),  intent(in)      :: gridfile
    type(regrid_type), intent(inout)   :: Regrid
    logical,              intent(in)   :: is_src_grid
    integer                            :: rcode, ncid, natt, id_zt, id_kmt, dims(2)
    integer                            :: id_grid_xt, id_grid_yt, id_grid_xc, id_grid_yc
    integer                            :: npes, i, j, k, layout(2)
    integer                            :: isc, iec, jsc, jec, isd, ied, jsd, jed
    logical                            :: is_new_grid
    character(len=128)                 :: name, catt
    integer,allocatable,dimension(:,:) :: kmt, kmc, kme, kmn
    real,   allocatable,dimension(:,:) :: tmp

    rcode = nf_open(trim(gridfile), NF_NOWRITE, ncid)    
    call error_handler('error in open file '//trim(gridfile), rcode)

    Regrid%ni = 0; Regrid%nj = 0; Regrid%nk = 0    
    !--- get vertical grid information
    rcode = nf_inq_varid(ncid, 'zt', id_zt)
    call error_handler('error in inquring id of field zt from file '//trim(gridfile), rcode)
    rcode = nf_inq_vardimid(ncid, id_zt, dims)
    call error_handler('error in inquring dims of field zt from file '//trim(gridfile), rcode)
    rcode = nf_inq_dimlen(ncid, dims(1), Regrid%nk)
    call error_handler('error in inquring dimlen of field zt from file '//trim(gridfile), rcode)
    allocate(Regrid%zt(Regrid%nk) ) 
    rcode = nf_get_var_double(ncid, id_zt, Regrid%zt) 
    call error_handler('error in getting data of zt from file '//trim(gridfile), rcode)

    !--- get horizontal grid size 
    rcode = nf_inq_varid(ncid, 'num_levels', id_kmt)
    if(rcode == 0) then
       is_new_grid = .true.
    else
       is_new_grid = .false.
       rcode = nf_inq_varid(ncid, 'kmt', id_kmt)
    endif

    call error_handler('Can not find kmt/num_levels in file '//trim(gridfile), rcode)

    rcode = nf_inq_vardimid(ncid, id_kmt, dims)
    call error_handler('error in inquring dims of field num_levels/kmt from file '//trim(gridfile), rcode)
    rcode = nf_inq_dimlen(ncid, dims(1), Regrid%ni )
    call error_handler('error in inquring dimlen ni of field num_levels/kmt from file '//trim(gridfile), rcode)
    rcode = nf_inq_dimlen(ncid, dims(2), Regrid%nj )
    call error_handler('error in inquring dimlen nj of field num_levels/kmt from file '//trim(gridfile), rcode)

    !--- get global attributes ( cyclic or tripolar )
    rcode =  nf_inq_natts(ncid, natt)
    call error_handler('error in inquring natts of file '//trim(gridfile), rcode )
    do i = 1, natt
       rcode = nf_inq_attname(ncid,NF_GLOBAL,i,name )
       call error_handler('error in inquring att name of file '//trim(gridfile), rcode )
       catt = ''
       select case( trim(name) )
       case('x_boundary_type')
          rcode = nf_get_att_text(ncid,NF_GLOBAL,name,catt)
          call error_handler('error in inquring x_boundary_type value of file '//trim(gridfile), rcode )
          if(trim(catt) == 'cyclic') then
             Regrid%is_cyclic = .true.
             Write(stdout(),*)' NOTE: x_boundary_type of grid '//trim(gridfile)//' is cyclic'
          else
             Write(stdout(),*)' NOTE: x_boundary_type of '//trim(gridfile)//' is solid_walls'
          endif
       case ('y_boundary_type')
          rcode = nf_get_att_text(ncid,NF_GLOBAL,name,catt)
          call error_handler('error in inquring y_boundary_type value of file '//trim(gridfile), rcode )
          if (trim(catt) == 'fold_north_edge') then
             Regrid%is_tripolar = .true.
             Write(stdout(),*)' NOTE: y_boundary_type of '//trim(gridfile)//' is tripolar'
          else
             Write(stdout(),*)' NOTE: y_boundary_type of '//trim(gridfile)//' is solid_walls' 
          endif
       end select
    end do

    !--- define domain decompsition for destination grid ---------------
    if(.not. is_src_grid) then
       npes = mpp_npes()
       layout = (/1,0/)
       call mpp_define_layout((/1,Regrid%ni,1,Regrid%nj/),npes,layout)

       if(Regrid%is_tripolar) then
          call mpp_define_domains((/1,Regrid%ni,1,Regrid%nj/),layout, Regrid%domain, xflags = CYCLIC_GLOBAL_DOMAIN, &
               yflags = FOLD_NORTH_EDGE, xhalo=1, yhalo=1)
       else if (Regrid%is_cyclic) then
          call mpp_define_domains((/1,Regrid%ni,1,Regrid%nj/),layout, Regrid%domain, xflags = cyclic_global_domain, &
               xhalo=1, yhalo=1)
       else
          call mpp_define_domains((/1,Regrid%ni,1,Regrid%nj/),layout, Regrid%domain,  xhalo=1, yhalo=1)
       endif
       call mpp_get_compute_domain(Regrid%domain,Regrid%isc,Regrid%iec,Regrid%jsc,Regrid%jec)
    endif

    !--- when the grid is source grid, need to get global grid information.
    !--- when the grid is destination grid, only need to get local grid information.
    if(is_src_grid) then
       isc = 1; iec = Regrid%ni
       jsc = 1; jec = Regrid%nj
    else
       isc = Regrid%isc  ; iec = Regrid%iec
       jsc = Regrid%jsc  ; jec = Regrid%jec
    endif
    isd = isc-1       ; ied = iec+1
    jsd = jsc-1       ; jed = jec+1

    ! get number vertical level
    allocate(kmt(isd:ied,jsd:jed), kmc(isc:iec,jsc:jec), kme(isc:iec,jsc:jec), kmn(isc:iec,jsc:jec)  )
    allocate(tmp(Regrid%ni,Regrid%nj) )    

    rcode = nf_get_var_double(ncid, id_kmt, tmp) 
    call error_handler('error in getting value of kmt from file '//trim(gridfile), rcode) 
    kmt = 0
    kmt(isc:iec,jsc:jec) = tmp(isc:iec,jsc:jec)
    if(is_src_grid) then
       kmt(isd,jsc:jec) = kmt(iec,jsc:jec)
       kmt(ied,jsc:jec) = kmt(isc,jsc:jec)
    else
       call mpp_update_domains(kmt,Dst%domain)    
    endif

    !--- calculate the mask at T,C,E,N-cell if needed
    if( any(fld_pos == 'T') ) then
       allocate(Regrid%T%mask(isc:iec,jsc:jec,Regrid%nk))
       do k=1,Regrid%nk
          do j=jsc,jec
             do i=isc,iec
                if (kmt(i,j) .ge. k) then
                   Regrid%T%mask(i,j,k) = 1.0
                else
                   Regrid%T%mask(i,j,k) = 0.0
                endif
             enddo
          enddo
       enddo
    endif

    if( any(fld_pos == 'C') ) then
       do j = jsc, jec
          do i = isc, iec
             kmc(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
          enddo
       enddo

       allocate(Regrid%C%mask(isc:iec,jsc:jec,Regrid%nk))
       do k=1,Regrid%nk
          do j=jsc,jec
             do i=isc,iec
                if (kmc(i,j) .ge. k) then
                   Regrid%C%mask(i,j,k) = 1.0
                else
                   Regrid%C%mask(i,j,k) = 0.0
                endif
             enddo
          enddo
       enddo
    endif

    if( any(fld_pos == 'E') ) then
       do j = jsc, jec
          do i = isc, iec
             kme(i,j) = min(kmt(i,j), kmt(i+1,j))
          enddo
       enddo

       allocate(Regrid%E%mask(isc:iec,jsc:jec,Regrid%nk))
       do k=1,Regrid%nk
          do j=jsc,jec
             do i=isc,iec
                if (kme(i,j) .ge. k) then
                   Regrid%E%mask(i,j,k) = 1.0
                else
                   Regrid%E%mask(i,j,k) = 0.0
                endif
             enddo
          enddo
       enddo
    endif

    if( any(fld_pos == 'N') ) then
       do j = jsc, jec
          do i = isc, iec
             kmn(i,j) = min(kmt(i,j), kmt(i,j+1))
          enddo
       enddo

       allocate(Regrid%N%mask(isc:iec,jsc:jec,Regrid%nk))
       do k=1,Regrid%nk
          do j=jsc,jec
             do i=isc,iec
                if (kmn(i,j) .ge. k) then
                   Regrid%N%mask(i,j,k) = 1.0
                else
                   Regrid%N%mask(i,j,k) = 0.0
                endif
             enddo
          enddo
       enddo
    endif

    !--- get the destination axis grid information ---------------------------------
    if(.not. is_src_grid) then
       allocate(Regrid%grid_xt(Regrid%ni),Regrid%grid_yt(Regrid%nj) )
       allocate(Regrid%grid_xc(Regrid%ni),Regrid%grid_yc(Regrid%nj) )  
       if(is_new_grid) then
          rcode = nf_inq_varid(ncid, 'grid_x_T', id_grid_xt)
          call error_handler('error in inquring variable id of field grid_x_T', rcode )
          rcode = nf_inq_varid(ncid, 'grid_y_T', id_grid_yt)
          call error_handler('error in inquring variable id of field grid_y_T', rcode )
          rcode = nf_inq_varid(ncid, 'grid_x_C', id_grid_xc)
          call error_handler('error in inquring variable id of field grid_x_C', rcode )
          rcode = nf_inq_varid(ncid, 'grid_y_C', id_grid_yc)
          call error_handler('error in inquring variable id of field grid_y_C', rcode )
       else
          rcode = nf_inq_varid(ncid, 'gridlon_t', id_grid_xt)
          call error_handler('error in inquring variable id of field gridlon_t', rcode )
          rcode = nf_inq_varid(ncid, 'gridlat_t', id_grid_yt)
          call error_handler('error in inquring variable id of field gridlat_t', rcode )
          rcode = nf_inq_varid(ncid, 'gridlon_c', id_grid_xc)
          call error_handler('error in inquring variable id of field gridlon_c', rcode )
          rcode = nf_inq_varid(ncid, 'gridlat_c', id_grid_yc)
          call error_handler('error in inquring variable id of field gridlat_c', rcode )
       endif
       rcode =  nf_get_var_double(ncid, id_grid_xt, Regrid%grid_xt)
       call error_handler('error in getting data of variable grid_x_T/gridlon_t', rcode )       
       rcode =  nf_get_var_double(ncid, id_grid_yt, Regrid%grid_yt)
       call error_handler('error in getting data of variable grid_y_T/gridlat_t', rcode )    
       rcode =  nf_get_var_double(ncid, id_grid_xc, Regrid%grid_xc)
       call error_handler('error in getting data of variable grid_x_C/gridlon_c', rcode )       
       rcode =  nf_get_var_double(ncid, id_grid_yc, Regrid%grid_yc)
       call error_handler('error in getting data of variable grid_y_C/gridlat_c', rcode ) 
    endif

    if(any(fld_pos == 'T') ) then
       call get_cell_info(Regrid%T,'T',ncid,Regrid%ni,Regrid%nj,isc,iec,jsc,jec,is_new_grid)
    endif

    if(any(fld_pos == 'C') ) then
       call get_cell_info(Regrid%C,'C',ncid,Regrid%ni,Regrid%nj,isc,iec,jsc,jec,is_new_grid)
    endif

    if(any(fld_pos == 'E') ) then
       call get_cell_info(Regrid%E,'E',ncid,Regrid%ni,Regrid%nj,isc,iec,jsc,jec,is_new_grid)
    endif

    if(any(fld_pos == 'N') ) then
       call get_cell_info(Regrid%N,'N',ncid,Regrid%ni,Regrid%nj,isc,iec,jsc,jec,is_new_grid)
    endif


    deallocate(kmt, kmc, kme, kmn,tmp)

    rcode = nf_close(ncid)

  end subroutine get_grid_info

  !#####################################################################
  subroutine get_cell_info(Cell, type, ncid, ni, nj, isc, iec,jsc,jec,is_new_grid )     
    type(cell_type), intent(inout) :: Cell
    character(len=1),   intent(in) :: type
    integer,            intent(in) :: ncid, ni, nj
    integer,            intent(in) :: isc, iec, jsc, jec
    logical,            intent(in) :: is_new_grid
    character(len=1)               :: lc_type, uc_type
    integer                        :: rcode, id_x, id_y, id_angle, id_sinrot, id_cosrot
    real, allocatable              :: tmp(:,:)

    allocate(Cell%geolon(isc:iec, jsc:jec), Cell%geolat(isc:iec, jsc:jec) )
    allocate(Cell%sinrot(isc:iec, jsc:jec), Cell%cosrot(isc:iec, jsc:jec) )
    allocate(tmp(ni,nj) )

    lc_type = lowercase(type)
    uc_type = uppercase(type)

    if(is_new_grid) then
       rcode = nf_inq_varid(ncid, 'x_'//uc_type, id_x)  
       call error_handler('error in inquring variable id of field x_'//uc_type, rcode )
       rcode = nf_inq_varid(ncid, 'y_'//uc_type, id_y)    
       call error_handler('error in inquring variable id of field y_'//uc_type, rcode )
       rcode = nf_inq_varid(ncid, 'angle_'//uc_type, id_angle)
       call error_handler('error in inquring variable id of field angle_'//uc_type, rcode )
    else
       rcode = nf_inq_varid(ncid, 'geolon_'//lc_type, id_x)
       call error_handler('error in inquring variable id of field geolon_'//lc_type, rcode )
       rcode = nf_inq_varid(ncid, 'geolat_'//lc_type, id_y) 
       call error_handler('error in inquring variable id of field geolat_'//lc_type, rcode )
       if(lc_type == 'c') then
          rcode = nf_inq_varid(ncid, 'sin_rot', id_sinrot)
          call error_handler('error in inquring variable id of field sin_rot', rcode )
          rcode = nf_inq_varid(ncid, 'cos_rot', id_cosrot)  
          call error_handler('error in inquring variable id of field cos_rot', rcode )
       else
          write(stdout(),*)'==> NOTE: No rotation information avaible for '// &
               uc_type//'-cell, will set sinrot = 0 and cosrot = 1' 
       endif
    endif

    rcode =  nf_get_var_double(ncid, id_x, tmp)
    call error_handler('error in getting data of variable x_'//uc_type//'/geolon_'//lc_type, rcode )
    Cell%geolon(isc:iec, jsc:jec) = tmp(isc:iec,jsc:jec)
    rcode =  nf_get_var_double(ncid, id_y, tmp)
    call error_handler('error in getting data of variable y_'//uc_type//'/geolat_'//lc_type, rcode )
    Cell%geolat(isc:iec, jsc:jec) = tmp(isc:iec,jsc:jec)

    if(is_new_grid) then
       rcode =  nf_get_var_double(ncid, id_angle, tmp) 
       call error_handler('error in getting data of variable angle_'//uc_type, rcode )
       Cell%sinrot(isc:iec,jsc:jec) = sin(tmp(isc:iec,jsc:jec))
       Cell%cosrot(isc:iec,jsc:jec) = cos(tmp(isc:iec,jsc:jec))
    else
       if(lc_type == 'c') then
          rcode =  nf_get_var_double(ncid, id_sinrot, tmp) 
          call error_handler('error in getting data of variable sin_rot', rcode )
          Cell%sinrot(isc:iec,jsc:jec) = tmp(isc:iec,jsc:jec)
          rcode =  nf_get_var_double(ncid, id_cosrot, tmp) 
          call error_handler('error in getting data of variable cos_rot', rcode )
          Cell%cosrot(isc:iec,jsc:jec) = tmp(isc:iec,jsc:jec)
       else
          Cell%sinrot = 0
          Cell%cosrot = 1
       endif
    endif

    deallocate(tmp)

  end subroutine get_cell_info

  !#####################################################################
  !--- set up axis and field meta data for destination data file

  subroutine setup_meta
    integer            :: ncid_src, rcode, id_fld, id_dim, id_time, dims(4)
    integer            :: type, len, ndims, natt_dim, natt, dimids(4), ndims_dst
    integer            :: n, m, i, j, num_has_taxis, start(4), nwrite(4)
    integer            :: id_t1_src, id_t2_src, id_dt_src, dim_time_dst
    integer, parameter :: max_dim = 10
    integer            :: dims_dst(max_dim), id_axes(max_dim), fld_dims(4), dims_id(max_dim)
    character(len=1)   :: dims_pos(max_dim), dims_cart(max_dim)
    character(len=128) :: name, att_name, cart, dims_name(max_dim)
    character(len=512) :: catt
    logical            :: found_dim = .false.
    logical            :: has_missing(max_flds)
    integer            ::fsize = 65536, inital = 0

    !--- open the dst_data file from root_pe.
    if(is_root_pe) then
#ifdef use_netCDF3
        rcode = NF__CREATE(trim(dst_data), NF_CLOBBER, inital, fsize, ncid_dst)
#else  ! netcdf 4
        rcode = NF__CREATE(trim(dst_data), IOR(NF_NETCDF4,NF_CLASSIC_MODEL), inital, fsize, ncid_dst )
#endif
       call error_handler('error in creating file '//trim(dst_data), rcode) 
    endif

    !--- get src_data information ( ntime, global attributes,  axis and field meta information )
    rcode = nf_open(trim(src_data), NF_NOWRITE, ncid_src) 
    call error_handler('error in opening file '//trim(src_data), rcode) 

    !--- get missing value, tavg_information, has_zaxis or has_taxis
    num_has_taxis = 0
    ndims_dst     = 0
    tavg_exist = .false.
    ntime_src = 1
    do n = 1, num_flds
       rcode = nf_inq_varid(ncid_src,fld_name(n), id_fld)
       call error_handler('error in inquring id of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       rcode = nf_inq_varndims( ncid_src, id_fld, ndims )
       call error_handler('error in inquring ndims of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       rcode = nf_inq_vardimid( ncid_src, id_fld, dimids )
       call error_handler('error in inquring dimids of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       rcode = nf_inq_varnatts(ncid_src, id_fld, natt )
       call error_handler('error in inquring natt of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       !--- to check if each field has z or time axis

       do m = 1, ndims
          rcode = nf_inq_dimname ( ncid_src, dimids(m), name )
          call error_handler('error in inquring dimension name of dimid', rcode)
          rcode = nf_inq_varid(ncid_src, name, id_dim)                    
          call error_handler('error in inquring variable id of '//trim(name), rcode)
          rcode = nf_inq_varnatts(ncid_src, id_dim, natt_dim)
          call error_handler('error in inquring number of attributes of variable '//trim(name), rcode)
          if(lowercase(trim(name)) == 'time') then
             cart = 'T'
             if(num_has_taxis == 0) then
                rcode = nf_inq_dimlen(ncid_src, dimids(m), ntime_src)
                call error_handler('error in inquring time dimension from file '//trim(src_data), rcode)
                allocate(time_value(ntime_src) )
                rcode = nf_get_var_double(ncid_src, id_dim, time_value) 
                call error_handler('error in inquring time value from file '//trim(src_data), rcode)
             endif
             num_has_taxis = num_has_taxis + 1   
          else             
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

       has_missing(n) = .false.
       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_fld,i,name )          
          call error_handler('error in inquring attribute name', rcode)
          if(trim(name) == 'time_avg_info') then
             tavg_exist = .true.
          else if(trim(name) == 'missing_value') then
             has_missing(n) = .true.
             rcode = nf_get_att_double(ncid_src, id_fld, name, missing_value(n) )
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
       allocate( t1_avg(ntime_src), t2_avg(ntime_src), dt_avg(ntime_src) )
       rcode = nf_inq_varid(ncid_src, 'average_T1', id_t1_src)
       call error_handler('error in inquring variable id of average_T1', rcode)
       rcode = nf_get_var_double(ncid_src,id_t1_src, t1_avg )
       call error_handler('error in getting average_T1 data', rcode)
       rcode = nf_inq_varid(ncid_src, 'average_T2', id_t2_src)
       call error_handler('error in inquring variable id of average_T2', rcode)
       rcode = nf_get_var_double(ncid_src,id_t2_src, t2_avg )
       call error_handler('error in getting average_T2 data', rcode)
       rcode = nf_inq_varid(ncid_src, 'average_DT', id_dt_src)
       call error_handler('error in inquring variable id of average_Dt', rcode)
       rcode = nf_get_var_double(ncid_src,id_dt_src, dt_avg )
       call error_handler('error in getting average_DT data', rcode)
    endif

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

    !--- read dimension information from src_data and write them to dst_data file.
    if( .not. is_root_pe) then
       rcode = nf_close(ncid_src)
       return
    endif

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
             len = Dst%ni
          case('Y')
             len = Dst%nj
          case('Z')
             len = Dst%nk             
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
       rcode = nf_def_var(ncid_dst, fld_name(n), NF_DOUBLE, ndims, fld_dims(1:ndims), id_dst_fld(n) )
       !--- copy the source field attribute to destination field.
       rcode = nf_inq_varnatts(ncid_src, id_fld, natt )
       call error_handler('error in inquring natt of field '//trim(fld_name(n))// ' of file '//trim(src_data), rcode)
       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_fld, i, name)
          call error_handler('error in inquring attribute '//trim(name), rcode)  
          if(trim(name) .ne. '_FillValue') then
             rcode = nf_copy_att(ncid_src,id_fld, trim(name), ncid_dst, id_dst_fld(n))
             call error_handler('error in copy attribute '//trim(name)//' of variable '//trim(fld_name(n)), rcode)  
          endif
       enddo
       if(.not. has_missing(n)) then
          rcode = nf_put_att_double(ncid_dst, id_dst_fld(n), 'missing_value',  NF_DOUBLE, 1, missing_value(n:n) )
          call error_handler('error in put missing_value attribute of var '//trim(fld_name(n))// &
                        ' of file '//trim(dst_data), rcode)   
       endif
       call error_handler('error in defining var '//trim(fld_name(n))//' of file '//trim(dst_data), rcode)       
    enddo

    !--- define time avg info, suppose the time_avg_info is specified by average_T1, average_T2, average_DT.
    if(tavg_exist) then
       rcode = nf_def_var(ncid_dst, 'average_T1', NF_DOUBLE, 1, dim_time_dst, id_t1_dst)
       rcode = nf_inq_varnatts(ncid_src, id_t1_src, natt )
       call error_handler('error in inquring natt of average_T1 of file '//trim(src_data), rcode)       
       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_t1_src, i, name)
          call error_handler('error in inquring average_T1 attribute '//trim(name), rcode)  
          rcode = nf_copy_att(ncid_src,id_t1_src, name, ncid_dst, id_t1_dst)
          call error_handler('error in copy average_T1 attribute '//trim(name), rcode) 
       enddo
       rcode = nf_def_var(ncid_dst, 'average_T2', NF_DOUBLE, 1, dim_time_dst, id_t2_dst)
       rcode = nf_inq_varnatts(ncid_src, id_t2_src, natt )
       call error_handler('error in inquring natt of average_T2 of file '//trim(src_data), rcode)       
       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_t2_src, i, name)
          call error_handler('error in inquring average_T2 attribute '//trim(name), rcode)  
          rcode = nf_copy_att(ncid_src,id_t2_src, name, ncid_dst, id_t2_dst)
          call error_handler('error in copy average_T2 attribute '//trim(name), rcode) 
       enddo
       rcode = nf_def_var(ncid_dst, 'average_DT', NF_DOUBLE, 1, dim_time_dst, id_dt_dst)
       rcode = nf_inq_varnatts(ncid_src, id_dt_src, natt )
       call error_handler('error in inquring natt of average_DT of file '//trim(src_data), rcode)       
       do i = 1, natt
          rcode = nf_inq_attname(ncid_src, id_dt_src, i, name)
          call error_handler('error in inquring average_DT attribute '//trim(name), rcode)  
          rcode = nf_copy_att(ncid_src,id_dt_src, name, ncid_dst, id_dt_dst)
          call error_handler('error in copy average_DT attribute '//trim(name), rcode) 
       enddo

    endif

    rcode = nf_enddef(ncid_dst)
    rcode = nf_close(ncid_src)

    !--- write axis data to dst_data
    start = 1; nwrite = 1
    do j = 1, ndims_dst
       select case(dims_cart(j))
       case('X')
          nwrite(1) = Dst%ni
          if(dims_pos(j) == 'T' .or. dims_pos(j) == 'N') then
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, Dst%grid_xt )
          else
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, Dst%grid_xc )             
          endif
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       case('Y')
          nwrite(1) = Dst%nj
          if(dims_pos(j) == 'T' .or. dims_pos(j) == 'E') then
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, Dst%grid_yt )
          else
             rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, Dst%grid_yc )             
          endif
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       case('Z')
          nwrite(1) = Dst%nk
          rcode = nf_put_vara_double(ncid_dst,id_axes(j), start, nwrite, Dst%zt )
          call error_handler('error in putting data of '//trim(dims_name(j))//' to file '//trim(dst_data), rcode)
       end select
    enddo

  end subroutine setup_meta

  !#####################################################################

  subroutine process_data
    integer                               :: ncid_src, rcode, i, j, n, m, k, l, nf, nk_src, nk_dst
    integer                               :: start(4), nread(4), nwrite(4), id_fld(num_flds)
    integer                               :: isc, iec, jsc, jec
    integer                               :: kstart, kend
    real                                  :: temp1, temp2
    real, dimension(:,:,:,:), allocatable :: data_src, data_dst
    real, dimension(:,:,:),   allocatable :: tmp1, tmp2, tmp_src
    real, dimension(:,:,:),   allocatable :: depth_src_3d, depth_dst_3d
    real, dimension(:,:),     allocatable :: mask_in
    type(horiz_interp_type),      pointer :: Interp   => NULL()
    real, dimension(:,:),         pointer :: sinrot   => NULL(), cosrot => NULL()
    real, dimension(:,:,:),       pointer :: mask     => NULL()
    logical                               :: is_first = .true.

    isc = Dst%isc; iec = Dst%iec; jsc = Dst%jsc; jec = Dst%jec;    

    !-------------------------------------------------------------------
    !--- read src data, remap onto current grid and write output -------
    rcode = nf_open(trim(src_data), NF_NOWRITE, ncid_src)
    do n = 1, num_flds
       rcode = nf_inq_varid(ncid_src, fld_name(n), id_fld(n) )
       call error_handler('error in inquring varid of '//trim(fld_name(n))//' of file '//trim(src_data), rcode)
    enddo

    write(stdout(),*)'***************************************'
    write(stdout(),*)' ntime_src is ', ntime_src

    if(do_vertical_interp) then
       allocate(depth_src_3d(isc:iec,jsc:jec,Src%nk), depth_dst_3d(isc:iec,jsc:jec,Dst%nk))
       do k=1, Src%nk
          depth_src_3d(:,:,k) = Src%zt(k)
       enddo
       do k=1, Dst%nk
          depth_dst_3d(:,:,k) = Dst%zt(k)
       enddo
       ! for vertical interpolation, Set value of levels shallower than the shallowest source level to 
       ! be the source value at shallowest level. 
       ! Set value of levels deeper than the deepest source level to 
       ! be the source value at deepest level. 
       
       do kstart = 1, Dst%nk
          if( Dst%zt(kstart) .GE. Src%zt(1) ) exit 
       enddo
       do kend = Dst%nk, 1, -1
          if( Dst%zt(kend) .LE. Src%zt(Src%nk) ) exit 
       enddo
       if( kstart > 1 ) then
          write(stdout(),*)"NOTE from regrid_3d: the value from level 1 to level ", kstart-1, &
               " will be set to the value at the shallowest source levle."
       endif
       if( kend < Dst%nk ) then
          write(stdout(),*)"NOTE from regrid_3d: the value from level ", kend+1, " to level ", Dst%nk, &
               " will be set to the value at the deepest source levle."
       endif

    endif

    allocate(mask_in(Src%ni, Src%nj) )

    do m = 1, ntime_src
       write(stdout(),*)'******************* at time level', m
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
          if(vector_fld(n) .and. Src%is_tripolar) nf = 2
          nk_src = 1
          if(has_zaxis(n)) nk_src = Src%nk
          nk_dst = 1
          if(has_zaxis(n)) nk_dst = Dst%nk

          write(stdout(),*)'************************************************************'
          write(stdout(),*)'The following are for the # ',n,' field '//trim(fld_name(n))
          write(stdout(),*)'has_zaxis is ', has_zaxis(n)
          write(stdout(),*)'nk_dst is  ', nk_dst
          allocate(data_src(Src%ni, Src%nj, nk_src,nf), data_dst(Dst%ni, Dst%nj, nk_dst,nf) )
          allocate(tmp_src(Src%ni, Src%nj, nk_src) )
          allocate(tmp1(isc:iec,jsc:jec,nk_src), tmp2(isc:iec,jsc:jec,nk_dst) )

          select case(fld_pos(n))
          case('T')
             Interp => Interp_T
          case('C')
             Interp => Interp_C
          case('E')
             Interp => Interp_E
          case('N')
             Interp => Interp_N
          end select

          start = 1; nread = 1
          nread(1) = Src%ni; nread(2) = Src%nj;
          if(has_zaxis(n)) then
             nread(3) = nk_src; start(4) = m
          else 
             start(3) = m
          endif
          rcode = nf_get_vara_double(ncid_src,id_fld(n), start, nread, data_src(:,:,:,1) ) 
          call error_handler('error in getting '//trim(fld_name(n))//' from file '//trim(src_data), rcode )

          !--- if the field is a vector field and src_grid is tripolar, read the next field and 
          !--- rotate data onto spherical grid.
          if( vector_fld(n) .and. Src%is_tripolar ) then
             select case(fld_pos(n))
             case('T')
                sinrot => Src%T%sinrot; cosrot => Src%T%cosrot
             case('C')
                sinrot => Src%C%sinrot; cosrot => Src%C%cosrot
             case('E')
                sinrot => Src%E%sinrot; cosrot => Src%E%cosrot
             case('N')
                sinrot => Src%N%sinrot; cosrot => Src%N%cosrot
             end select

             rcode = nf_get_vara_double(ncid_src,id_fld(n+1), start, nread, data_src(:,:,:,2) ) 

             do k = 1, nk_src
                do j = 1, Src%nj 
                   do i = 1, Src%ni
                      if(data_src(i,j,k,1) /= missing_value(n)) then
                         temp1             = data_src(i,j,k,1)*cosrot(i,j)-data_src(i,j,k,2)*sinrot(i,j)
                         temp2             = data_src(i,j,k,1)*sinrot(i,j)+data_src(i,j,k,2)*cosrot(i,j)
                         data_src(i,j,k,1) = temp1
                         data_src(i,j,k,2) = temp2
                      endif
                   enddo
                enddo
             enddo
          endif

          !--- extrap data onto land grid 
          do l = 1, nf
             if(use_source_vertical_grid) then !--- the source grid mask will decide the destination data mask.
                if(apply_mask) then
                   do k = 1, size(data_src,3)
                      where(data_src(:,:,k,l) == missing_value(n) )
                         mask_in = 0.0
                      elsewhere
                         mask_in = 1.0
                      end where
                      call horiz_interp(Interp, data_src(:,:,k,l), tmp2(:,:,k), mask_in=mask_in, &
                             missing_value=missing_value(n) )
                   enddo
                else
                   if(any(data_src(:,:,:,l) == missing_value(n)))  then
                      if(Src%is_tripolar .and. is_first) then
                         is_first = .false.
                         write(stdout(),*)'WARNING: when the source is tripolar, ',      &
                             'the laplacian extrapolation will not very accurate at the tripolar region. ', &
                             'Even the source grid is lat-lon grid, the extrapolation will not also not  ', &
                             'very accurate if the lat-lon grid is not uniform. '
                      endif
                      call extrap(data_src(:,:,:,l), tmp_src, stop_crit(n), missing_value(n), fld_pos(n) ) 
                      !--- remap data onto destination horizontal grid
                      do k = 1, size(tmp_src,3)
                         call horiz_interp(Interp, tmp_src(:,:,k), tmp2(:,:,k))
                      enddo
                   else
                      !--- remap data onto destination horizontal grid
                      do k = 1, size(data_src,3)
                         call horiz_interp(Interp, data_src(:,:,k,l), tmp2(:,:,k))
                      enddo
                   endif
                endif
             else
                if(any(data_src(:,:,:,l) == missing_value(n)))  then
                   if(Src%is_tripolar .and. is_first) then
                      is_first = .false.
                      write(stdout(),*)'WARNING: when the source is tripolar, ',      &
                           'the laplacian extrapolation will not very accurate at the tripolar region. ', &
                           'Even the source grid is lat-lon grid, the extrapolation will not also not  ', &
                           'very accurate if the lat-lon grid is not uniform. '
                   endif
                   call extrap(data_src(:,:,:,l), tmp_src, stop_crit(n), missing_value(n), fld_pos(n) ) 
                   !--- remap data onto destination horizontal grid
                   do k = 1, size(tmp_src,3)
                      call horiz_interp(Interp, tmp_src(:,:,k), tmp1(:,:,k))
                   enddo
                else
                   !--- remap data onto destination horizontal grid
                   do k = 1, size(data_src,3)
                      call horiz_interp(Interp, data_src(:,:,k,l), tmp1(:,:,k))
                   enddo
                endif

                if(has_zaxis(n) .and. do_vertical_interp ) then
                   do k = 1, kstart-1
                      tmp2(:,:,k) = tmp1(:,:,1)
                   enddo

                   do k = kend+1, nk_dst
                      tmp2(:,:,k) = tmp1(:,:,nk_dst)
                   enddo
                   call interp_1d(depth_src_3d,depth_dst_3d(:,:,kstart:kend),tmp1,tmp2(:,:,kstart:kend))
                else
                   tmp2 = tmp1
                endif

                !--- apply mask if needed,
                if(apply_mask) then
                   select case(fld_pos(n))
                   case('T')
                      mask => Dst%T%mask
                   case('C')
                      mask => Dst%C%mask
                   case('E')
                      mask => Dst%E%mask
                   case('N')
                      mask => Dst%N%mask
                   end select
                   do k = 1, nk_dst
                      do j = jsc, jec
                         do i = isc, iec
                            if(mask(i,j,k) < 0.5) tmp2(i,j,k) = missing_value(n)
                         enddo
                      enddo
                   enddo
                endif
             endif
             !--- get global data ------------------------------
             call mpp_domains_set_stack_size(2*Dst%ni*Dst%nj*nk_dst)
             call mpp_global_field(Dst%domain,tmp2, data_dst(:,:,:,l) )                
          enddo

      !--- rotate the data if the field is a vector field and dst_grid is tripolar.
          if( vector_fld(n) .and. Dst%is_tripolar ) then
             select case(fld_pos(n))
             case('T')
                sinrot => Dst%T%sinrot; cosrot => Dst%T%cosrot
             case('C')
                sinrot => Dst%C%sinrot; cosrot => Dst%C%cosrot
             case('E')
                sinrot => Dst%E%sinrot; cosrot => Dst%E%cosrot
             case('N')
                sinrot => Dst%C%sinrot; cosrot => Dst%C%cosrot
             end select
             do k = 1, nk_dst
                do j = 1, Dst%nj
                   do i = 1, Dst%ni
                      if(data_dst(i,j,k,1) /= missing_value(n)) then
                         temp1             =  data_dst(i,j,k,1)*cosrot(i,j)+data_dst(i,j,k,2)*sinrot(i,j)
                         temp2             = -data_dst(i,j,k,1)*sinrot(i,j)+data_dst(i,j,k,2)*cosrot(i,j) 
                         data_dst(i,j,k,1) =  temp1
                         data_dst(i,j,k,2) =  temp2
                      endif
                   enddo
                enddo
             enddo

             select case(fld_pos(n))
             case('C')
                do k = 1, nk_dst
                   do i=1,Dst%ni/2                   
                      if(data_dst(Dst%ni-i,Dst%nj,k,1) /= missing_value(n)) then
                         data_dst(i,Dst%nj,k,1) = -1.*data_dst(Dst%ni-i,Dst%nj,k,1)
                         data_dst(i,Dst%nj,k,2) = -1.*data_dst(Dst%ni-i,Dst%nj,k,2)
                      endif
                   enddo
                enddo
             case('N')
                do k = 1, nk_dst
                   do i=1,Dst%ni/2  
                      if(data_dst(Dst%ni-i+1,Dst%nj,k,1) /= missing_value(n)) then
                         data_dst(i,Dst%nj,k,1) = -1.*data_dst(Dst%ni-i+1,Dst%nj,k,1)
                         data_dst(i,Dst%nj,k,2) = -1.*data_dst(Dst%ni-i+1,Dst%nj,k,2)
                      endif
                   enddo
                enddo
             end select
          endif
          !--- write out data from root_pe
          if(is_root_pe) then
             if(debug) write(stdout(),*)'NOTE: after regrid data chksum :', mpp_chksum(data_dst, (/mpp_root_pe()/) )
             do l = 1, nf
                start = 1; nwrite = 1
                nwrite(1) = Dst%ni; nwrite(2) = Dst%nj;
                if(has_zaxis(n)) then
                   nwrite(3) = nk_dst
                   start(4) = m
                rcode = nf_put_vara_double(ncid_dst, id_dst_fld(n+l-1), start, nwrite, data_dst(:,:,:,l) ) 
                call error_handler('error in putting '//trim(fld_name(n+l-1))//' in file '//trim(dst_data), rcode )
                else
                   start(3) = m
                rcode = nf_put_vara_double(ncid_dst, id_dst_fld(n+l-1), start, nwrite, data_dst(:,:,1,l) ) 
                call error_handler('error in putting '//trim(fld_name(n+l-1))//' in file '//trim(dst_data), rcode )
                endif
!                rcode = nf_put_vara_double(ncid_dst, id_dst_fld(n+l-1), start, nwrite, data_dst(:,:,:,l) ) 
!                call error_handler('error in putting '//trim(fld_name(n+l-1))//' in file '//trim(dst_data), rcode )
             enddo
          endif
          n = n+nf
          deallocate(data_src, data_dst, tmp_src, tmp1, tmp2)
       enddo
    enddo

    rcode = nf_close(ncid_src)
    if(is_root_pe) rcode = nf_close(ncid_dst)
    if(do_vertical_interp) deallocate(depth_src_3d, depth_dst_3d)
    deallocate(mask_in)

  end subroutine process_data

  !#####################################################################

  subroutine regrid_end


    call fms_io_exit
    call fms_end()

  end subroutine regrid_end

  !#####################################################################
  subroutine extrap(data_in, data_out, crit, missing, pos)
    real, dimension(:,:,:),  intent(in) :: data_in
    real, dimension(:,:,:), intent(out) :: data_out
    real,                    intent(in) :: crit, missing 
    character(len=1),        intent(in) :: pos     
    real                                :: resmax, initial_guess = 0.0
    integer                             :: ni,nj,nk, i, j, k, n
    real, dimension(0:size(data_in,1)+1, 0:size(data_in,2)+1) :: tmp
    real, dimension(size(data_in,1), size(data_in,2) )        :: sor, res

    ni = size(data_in,1)
    nj = size(data_in,2)
    nk = size(data_in,3)

    tmp = 0.0

    do j= 1, nj
       do i= 1, ni
          if(data_in(i,j,1) == missing) then
             tmp(i,j) = initial_guess
          endif
       enddo
    enddo


    do k=1,nk
       do j= 1, nj
          do i= 1, ni
             if(data_in(i,j,k) == missing) then
                sor(i,j) = rel_coef
             else
                tmp(i,j)=data_in(i,j,k)
                sor(i,j) = 0.0
             endif
          enddo
       enddo

       call fill_boundaries(tmp,pos)

       ! iterate
       n=1
       do            
          resmax=0.0
          do j= 1, nj
             do i= 1, ni
                res(i,j) = 0.25*(tmp(i-1,j)+tmp(i+1,j)+tmp(i,j-1)+tmp(i,j+1)) -tmp(i,j)
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

          call fill_boundaries(tmp,pos)

          n=n+1

       enddo

       write(stdout(),'(a,i4,a,i4)') 'At k-level= ',k
       write(stdout(),'(a,i6,a)') 'Stopped after ',n,' iterations'
       write(stdout(),'(a,f10.4)') 'maxres= ',resmax

    enddo

  end subroutine extrap

  !##########################################################

  subroutine fill_boundaries(data, pos)

    real, dimension(0:,0:), intent(inout) :: data
    character(len=1),          intent(in) :: pos
    integer :: i,j, ni, nj

    ni = size(data,1) - 2
    nj = size(data,2) - 2

    if(Dst%is_cyclic) then
       data(0,1:nj) = data(ni,1:nj)
       data(ni+1,1:nj) = data(1,1:nj)
    endif

    if(Dst%is_tripolar) then
       do i = 1, ni
          select case(pos)
          case('T')
             data(i,nj+1) = data(ni-i+1,nj)
          case('C')
             data(i,nj+1) = data(ni-i,nj-1)
          case('E')
             data(i,nj+1) = data(ni-i+1,nj-1)
          case('N')
             data(i,nj+1) = data(ni-i,nj)
          end select
       enddo
    endif

    return

  end subroutine fill_boundaries

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

end program regrid

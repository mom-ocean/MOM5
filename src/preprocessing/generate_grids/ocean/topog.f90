module topog_mod
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
  !
  !<OVERVIEW>
  ! <TT>topog_mod</TT> Generate topography for ocean model. 
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  !    The topography can be idealized or remapped from some source topography data.
  !    The type of topography is specified by the namelist variable "topography" and
  !    "topog_depend_on_vgrid". See the documentation of namelist variable "topography"
  !    and "topog_depend_on_vgrid" for details.
  !</DESCRIPTION>
  !
  !<NAMELIST NAME="topog_nml">
  !<DATA NAME="topography" TYPE="character(len=24)">
  !<PRE>
  ! rectangular_basin : Constructing a rectangular basin with a flat bottom
  ! bowl              : From "Simulation of density-driven frictional downslope 
  !                     flow in  z-coordinate mocean models"  Winton et al. 
  !                     JPO, Vol 28, No 11, 2163-2174,  November 1998
  ! gaussian          : sets "kmt" to a gaussian bump on a sloping bottom.
  ! dome              : similar (not identical) to DOME configuration of Legg etal Ocean Modelling (2005) 
  ! idealized         : generates an "idealized" not very realistic topography.
  ! all_land          : constructing a all land topography.
  ! from_file         : Remap the topography onto the current grid from some source data file.
  !</PRE>
  !</DATA>
  !<DATA NAME="topog_depend_on_vgrid" TYPE = "logical" >
  ! when topography /= "from_file", topog_depend_on_vgrid must be true (default value). 
  ! When it is false, topography is obtained by a simple remapping onto current grid.
  !</DATA>
  !<DATA NAME="topog_file" TYPE="character(len=128)">
  ! name of topograhy file (e.g. scripps, navy_topo, ...)
  !</DATA>
  !<DATA NAME="topog_field" TYPE="character(len=24)">
  ! name of topography field in file
  !</DATA>
  !<DATA NAME="flat_bottom" TYPE="logical">
  ! generate flat bottom over ocean points. Default value is false.
  !</DATA>
  !<DATA NAME="full_cell" TYPE="logical">
  ! do not generate partial bottom cells. Default value is false.
  !</DATA>
  !<DATA NAME="fill_isolated_cells" TYPE="logical">
  ! Do not allow non-advective tracer cells (strongly recommended). Default value is true.
  !</DATA>
  !<DATA NAME="dont_change_landmask" TYPE="logical">
  ! Do not change land/sea mask when filling isolated cells. Default value is false.
  !</DATA>
  !<DATA NAME="fill_shallow" TYPE="logical">
  ! Make cells less than minimum depth land. Default value is false. 
  !</DATA>
  !<DATA NAME="fill_first_row" TYPE="logical">
  ! if true make first row of ocean model all land points for ice model when
  ! topography is "from_file". It will do nothing when topography is not "from_file".
  ! Default value is true.
  !</DATA>
  !<DATA NAME="deepen_shallow" TYPE="logical">
  ! Make cells less than minimum depth equal to minimum depth. Default value is false.
  !</DATA>
  !<DATA NAME="round_shallow" TYPE="logical">
  ! Make cells land if depth is less than 1/2 mimumim depth, otherwise make ocean. Default value is false.
  !</DATA>
  !<DATA NAME="gauss_amp" TYPE="real" >
  ! height of gaussian bump as percentage of ocean depth
  !</DATA>
  !<DATA NAME="gauss_scale" TYPE="real">
  ! width of gaussian bump as percentag e of basin width
  !</DATA>
  !<DATA NAME="slope_x" TYPE="real" UNITS="(m/deg)">
  ! rise of the ocean floor to the east for the gaussian bump
  !</DATA>
  !<DATA NAME="slope_y" TYPE="real" UNITS="(m/deg)">
  ! rise of the ocean floor to the north for the gaussian bump
  !</DATA>
  !<DATA NAME="bowl_south" TYPE="real" UNITS="degrees">
  ! southern boundary of Winton bowl
  !</DATA>
  !<DATA NAME="bowl_north" TYPE="real" UNITS="degrees">
  ! northern boundary of Winton bowl
  !</DATA>
  !<DATA NAME="bowl_west" TYPE="real" UNITS="degrees">
  ! western boundary of Winton bowl
  !</DATA>
  !<DATA NAME="bowl_east" TYPE="real" UNITS="degrees">
  ! eastern boundary of Winton bowl
  !</DATA>
  !<DATA NAME="bowl_min_depth" TYPE="real" UNITS="meters">
  ! minimum depth of Winton bowl 
  !</DATA>
  !<DATA NAME="bowl_max_depth" TYPE="real" UNITS="meters">
  ! maximum depth of Winton bowl 
  !</DATA>
  !
  !<DATA NAME="dome_slope" TYPE="real" UNITS="radians">
  ! Slope for the dome configuration.  Default = 0.01
  !</DATA>
  !<DATA NAME="dome_bottom" TYPE="real" UNITS="m">
  ! Bottom of the dome configuration.  Default=3600.0
  !</DATA>
  !<DATA NAME="dome_embayment_west" TYPE="real" UNITS="dimensionless">
  ! western edge of dome embayment. Default=19.0
  !</DATA>
  !<DATA NAME="dome_embayment_east" TYPE="real" UNITS="dimensionless">
  ! eastern edge of dome embayment. Default=21.0
  !</DATA>
  !<DATA NAME="dome_embayment_south" TYPE="real" UNITS="dimensionless">
  ! southern edge of dome embayment. Default=69.0
  !</DATA>
  !<DATA NAME="dome_embayment_depth" TYPE="real" UNITS="m">
  ! Depth of the embayment. Default=600.0
  !</DATA>
  !<DATA NAME="kmt_min" TYPE="integer">
  ! minimum number of vertical levels
  !</DATA>
  !<DATA NAME="filter_topog" TYPE="logical">
  ! apply filter to topography. Default value is false.
  !</DATA>
  !<DATA NAME="num_filter_pass" TYPE="integer">
  ! number of passes of spatial filter
  !</DATA>
  !<DATA NAME="adjust_topo" TYPE="logical">
  ! adjust topography (enforce_min_depth;remove_isolated_cells;restrict_partial_cells)
  ! Strongly recommended. Default value is true.
  !</DATA>
  !<DATA NAME="fraction_full_cell" TYPE="real" >
  ! Fraction of the associated full cell that a corresponding partial cell thickness 
  ! is no smaller than.  That is, we maintain 
  ! partial_cell_min_dht(i,j,k) = fraction_full_cell*full_cell_dzt(k)
  ! If fraction_full_cell=0.0, then partial_cell_min_dht = min(zw(1), 50.0)
  !</DATA>
  !<DATA NAME="scale_factor" TYPE="real" >
  ! scaling factor for topography data (e.g. -1 to flip sign or 0.01 to convert from centimeters)
  !</DATA>
  !<DATA NAME="smooth_topo_allow_deepening" TYPE="logical">
  ! allow filter to deepen cells. Default value is false.
  !</DATA>
  !<DATA NAME="interp_method"  TYPE= "character(len=64)" >
  ! specifying the remapping method when remampping topography from source data to current grid.
  ! Its value can be "spherical", " bilinear" or "conservative". Default value is "bilinear". when the source
  ! topography is on the regular grid (nml src_is_spherical is true), "bilinear" interpolation 
  ! is recommanded, since bilinear interpolation will provide more smooth results than 
  ! "spherical" interpolation (especially when interpolating from coarse grid to fine grid). 
  ! Plus bilinear interpolation is much more efficiency than "spherical interpolation". 
  ! When the source data is on non-regular grid (nml src_is_spherical is false), "bilinear" 
  ! interpolation may not work well because the destination is not inside the source grid, 
  ! in this case, you need to set interp_method to "spherical".
  !</DATA>
  !<DATA NAME="num_nbrs" TYPE="integer">
  ! Number of nearest neighbors for regridding.
  !</DATA>
  !<DATA NAME="max_dist" TYPE="real" UNITS="radians">
  ! Maximum region of influence around destination grid points.
  !</DATA>
  !<DATA NAME="src_is_spherical" TYPE="logical" >
  ! Determine if the source grid is spherical grid or not. If true, source grid is spherical grid,
  ! otherwise not. Default value is .true. When src_is_spherical is .true., lon_field and lat_field
  ! need to be set.
  !</DATA>
  !<DATA NAME="lon_field" TYPE="character(len=24)">
  ! name of geographic longitude field in source file
  !</DATA>
  !<DATA NAME="lat_field" TYPE="character(len=24)">
  ! name of geographic latitude field in source file
  !</DATA>
  !<DATA NAME="output_corner_topog" TYPE="logical" >
  ! When true, will write topography information on the C-cell cneter to the grid file.
  ! Default value is false.
  !</DATA>
  !<DATA NAME="check_mask" TYPE="logical" >
  ! When true, check the possible masking ( all-land region) for certain layout. The print out message 
  ! will provide coupler_nml ( or ocean_solo_nml) entry : nmask, layout_mask and mask_list.
  ! Default value is false.
  !</DATA>
  !<DATA NAME="open_very_this_cell" TYPE="logical">
  !  When set to false, check which change is larger, opening or closing the cell, and 
  !  to do that with smaller  effect in depth_t. Default is true, will always 
  !  opening the cell.
  !</DATA>
  !<DATA NAME="min_thickness" TYPE="real" UNITS="METERS">
  !  minimum vertical thickness allowed. with default value 0.1. Increase or decrease this number as needed.
  !<DATA NAME="debug" TYPE="logical" >
  ! Control standard output. Default value is true so to show lots of information.
  !</DATA>
  !</NAMELIST>
  use mpp_mod,          only : mpp_error,FATAL, NOTE, mpp_npes, mpp_pe, mpp_root_pe, mpp_chksum
  use mpp_domains_mod,  only : domain2d, mpp_define_domains, mpp_define_layout, mpp_global_field
  use mpp_domains_mod,  only : mpp_get_compute_domain
  use mpp_io_mod,       only : axistype, fieldtype, mpp_write_meta, mpp_write, mpp_open
  use mpp_io_mod,       only : MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
  use mpp_io_mod,       only : mpp_get_info, mpp_get_fields, mpp_get_atts, mpp_get_axis_data
  use mpp_io_mod,       only : mpp_read, mpp_close
  use fms_mod,          only : lowercase, open_namelist_file, close_file, check_nml_error
  use fms_mod,          only : stdout, write_version_number, file_exist, string
  use axis_utils_mod,   only : nearest_index, get_axis_cart, get_axis_bounds
  use grids_type_mod,   only : hgrid_data_type, vgrid_data_type, topog_data_type
  use grids_util_mod,   only : write_field_meta, write_field_data, get_file_unit
  use horiz_interp_mod, only : horiz_interp
  use constants_mod,    only : pi
  use check_mask_mod,   only : checking_mask

  implicit none
  private

  !--- namelist variables ----------------------------------------------
  character(len=128) :: topog_file = 'scripps'
  character(len=24)  :: topog_field = 'topog'
  character(len=24)  :: topography = 'from_file'
  logical            :: topog_depend_on_vgrid = .TRUE.
  logical :: flat_bottom=.false.          ! sets ocean bottom to the deepest level when true
  logical :: full_cell=.false.            ! sets bottom cell thickness equal to thickness of vertical level when true 
  logical :: fill_isolated_cells=.true. 
  logical :: fill_first_row = .true.
  logical :: dont_change_landmask=.false. ! when filling isolated cells
  real    :: gauss_amp=0.50               ! height of gaussian bump as a percent of max ocean depth 
  real    :: gauss_scale=0.25             ! half width of gaussian bump as a percent of basin width
  real    :: slope_x=0.0                  ! (m/deg) rise of the ocean floor to the east for the gaussian bump
  real    :: slope_y=0.0                  ! (m/deg) rise of the ocean floor to the north for the gaussian bump 
  real    :: bowl_south=60.0              ! southern boundary of Winton bowl (deg)
  real    :: bowl_north=70.0              ! northern boundary of Winton bowl (deg)
  real    :: bowl_west=0.0                ! western boundary of Winton bowl (deg)
  real    :: bowl_east=20.0               ! eastern boundary of Winton bowl (deg)
  real    :: bowl_min_depth=500.0         ! unit is m
  real    :: bowl_max_depth=3000.0        ! unit is m
  real    :: read_epsln=.1                

  real    :: dome_slope=0.01              ! slope for the dome configuration (radians)
  real    :: dome_bottom=3600.0           ! bottom depth (m) of the dome configuration 
  real    :: dome_embayment_west=19.0     ! western edge of dome embayment (degrees)
  real    :: dome_embayment_east=21.0     ! eastern edge of dome embayment (degrees)
  real    :: dome_embayment_south=69.0    ! southern edge of the dome embayment (degrees)
  real    :: dome_embayment_depth=600.0   ! depth of the dome embayment (m)

  real    :: deg2metre=111.324e3          ! convert degrees latitude to metres  

  logical :: fill_shallow = .false.
  logical :: deepen_shallow = .false.
  logical :: round_shallow=.false.
  integer :: kmt_min = 2, num_filter_pass = 1, num_nbrs = 10
  logical :: filter_topog=.false.
  logical :: adjust_topo = .true.
  real    :: fraction_full_cell=0.20      ! partial cell thickness >= fraction_full_cell*dzt 
  real    :: scale_factor=1               ! scaling factor for topography (e.g. -1 to flip sign)
  logical :: smooth_topo_allow_deepening = .false.
  real    :: max_dist=0.1
  character(len=64) :: interp_method = 'bilinear'
  logical           :: src_is_spherical = .true.
  character(len=24) :: lon_field        = 'x_T'
  character(len=24) :: lat_field        = 'y_T'
  logical :: output_corner_topog        = .false.
  logical :: check_mask                 = .false.
  logical :: open_very_this_cell        = .true.
  real    :: min_thickness              = 0.1
  logical :: debug = .true.

  namelist /topog_nml/ topography, topog_depend_on_vgrid, topog_file, topog_field, &
       flat_bottom, full_cell, fill_isolated_cells, dont_change_landmask,          &
       fill_shallow, deepen_shallow, round_shallow, gauss_amp, gauss_scale,        &
       slope_x, slope_y, bowl_south, bowl_north, bowl_east, bowl_west,             &
       bowl_max_depth, bowl_min_depth, kmt_min, filter_topog, num_filter_pass,     &
       scale_factor, adjust_topo, smooth_topo_allow_deepening,num_nbrs,max_dist,   &
       debug, fill_first_row, fraction_full_cell, interp_method, src_is_spherical, &
       lon_field, lat_field, output_corner_topog, check_mask, open_very_this_cell, &
       min_thickness, dome_slope, dome_bottom, dome_embayment_west,                &
       dome_embayment_east, dome_embayment_south, dome_embayment_depth

!--- namelist variables for OBC ----------------------------------------------

  integer, parameter :: max_obc = 4      ! maximum number of open boundaries (increase if want more)
  integer                               :: nobc= 0
  character(len=10), dimension(max_obc) :: direction='none' ! open directions; to be consistent with is, ie, js, je
  integer, dimension(max_obc)           :: is=-999, ie=-999, js=-999, je=-999     ! boundary position
  integer, dimension(max_obc)           :: nsmooth=2        ! number of points to be smoothed
  character(len=32), dimension(max_obc) :: name='none' 
  
  namelist /obc_nml/ nobc, direction, is, ie, js, je, nsmooth, name

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: topog.f90,v 20.0 2013/12/14 00:30:50 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'  
  !--- other variables
  logical              :: module_is_initialized = .false.
  real                 :: grid_tol = 1.e-2
  type(fieldtype),save :: fld_depth, fld_num_levels, fld_wet
  real                 :: small = 0.
  logical              :: have_obc = .false.

  !---public interface -------------------------------------------------

  public :: generate_topog, topog_init, topog_end, process_topo, write_topog_global_meta 
  public :: write_topog_field_meta, write_topog_data, show_deepest,  set_topog_nml

contains

  !#######################################################################
  !<SUBROUTINE NAME="topog_init">
  !   <OVERVIEW>
  !    Initialization routine.
  !   </OVERVIEW>
  !
  !<DESCRIPTION>
  !   Read topography namelist.
  !</DESCRIPTION>
  !   <TEMPLATE>
  !     call topog_init(Topog, Hgrid)
  !   </TEMPLATE>
  !<IN NAME="Hgrid" TYPE="type(hgrid_data_type)">
  !   A derived-type variable that contains horizontal grid information.
  !</IN>
  !<INOUT NAME="Topog" TYPE="type(topog_data_type)">
  !   A derived-type variable that contains topography.
  !</INOUT>
  !</SUBROUTINE>

  subroutine topog_init(Topog, Hgrid)
    type(hgrid_data_type),   intent(in)  :: Hgrid
    type(topog_data_type), intent(inout) :: Topog

    integer :: unit, ierr, io, ni, nj

    !---- read namelist --------------------------------------------------
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
       read  (unit, nml=topog_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'topog_nml')  ! also initializes nml error codes
    enddo
10  call close_file (unit)
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
       read  (unit, nml=obc_nml, iostat=io, end=20)
       ierr = check_nml_error(io,'obc_nml')  ! also initializes nml error codes
    enddo
20  call close_file (unit)

    !--- write version info and namelist to logfile ----------------------

    call write_version_number(version,tagname)
    write (stdout(), nml=topog_nml)
    write (stdout(), nml=obc_nml)

    if(trim(topography) == 'from_file') then
       if(trim(interp_method) .ne. 'spherical' .and. trim(interp_method) .ne. 'bilinear' &
                                          .and. trim(interp_method) .ne. 'conservative') &
          call mpp_error(FATAL, 'topog_mod: interp_method = '// trim(interp_method)//' is not a valid option')
       if(src_is_spherical) then
          if( trim(interp_method) == 'spherical' ) write(stdout(),*)" NOTE from topog_mod: ", &
              "when src_is_spherical is true, interp_methos set to 'bilinear' is recommanded"
       endif
    endif

    if(.NOT. topog_depend_on_vgrid .AND. (trim(topography) .NE. 'from_file' .and. trim(topography) .NE. 'all_land') ) &
         call mpp_error(FATAL, &
        'topog_mod: when topography /= "from_file" or "all_land", topog_depend_on_vgrid should be true' )

    ni = Hgrid%ni
    nj = Hgrid%nj

    allocate(Topog%depth_t(ni,nj), Topog%num_levels(ni,nj))
    if(output_corner_topog) allocate(Topog%depth_c(ni, nj), Topog%num_levels_c(ni, nj))

    if (nobc > 0) have_obc = .true.

    module_is_initialized = .true.

    return

  end subroutine topog_init

  !#######################################################################
  !<SUBROUTINE NAME="generate_topog">
  !<OVERVIEW>
  !  generate topography data.
  !</OVERVIEW>
  !<DESCRIPTION>
  !Call horiz_interp to calculate regridded topography.
  !Perform topography checks
  !</DESCRIPTION>
  !   <TEMPLATE>
  !     call generate_topog(Hgrid, Topog, Vgrid)
  !   </TEMPLATE>
  !<IN NAME="Hgrid" TYPE="type(hgrid_data_type)">
  !  A derived-type variable that contains horizontal grid information.
  !</IN>
  !<IN NAME="Vgrid" TYPE="type(vgrid_data_type), optional">
  !  A derived-type variable that contains vertical grid information.
  !</IN>
  !<INOUT NAME="Topog" TYPE="type(topog_data_type)">
  !  A derived-type variable that contains topography data.
  !</INOUT>   
  !</SUBROUTINE>

  subroutine generate_topog( Topog, Hgrid, Vgrid)
    type(topog_data_type), intent(inout) :: Topog
    type(hgrid_data_type), intent(in)    :: Hgrid
    type(vgrid_data_type), intent(in), optional  :: Vgrid

    real, dimension(:,:), allocatable    :: ht
    real, dimension(:,:), allocatable    :: hu
    real, dimension(:,:), allocatable    :: dxte, dytn
    real, dimension(:,:), allocatable    :: lonv, latv
    real, dimension(:,:,:), allocatable  :: xv, yv
    integer, dimension(:,:), allocatable :: kmt
    integer, dimension(:,:), allocatable :: kmu
    logical                              :: tripolar_grid, cyclic_x, cyclic_y
    integer                              :: pe, ni ,nj, i, j

    if(topog_depend_on_vgrid) then
       if(.not. present(Vgrid)) call mpp_error(FATAL, 'topog_mod: when topog_depend_on_vgrid is true, Vgrid is needed')
    else
       if(present(Vgrid)) call mpp_error(FATAL, 'topog_mod: when topog_depend_on_vgrid is false, Vgrid is not needed')
    endif

    tripolar_grid = Hgrid%tripolar_grid
    cyclic_x      = Hgrid%cyclic_x
    cyclic_y      = Hgrid%cyclic_y

    ni = Hgrid%ni; nj = Hgrid%nj
    allocate(ht(0:ni+1, 0:nj+1), kmt(0:ni+1,0:nj+1), Topog%wet(ni,nj))
    kmt=0
    ht = 0.0
    if(output_corner_topog) then
    allocate(hu(0:ni+1, 0:nj+1), kmu(0:ni+1,0:nj+1), Topog%wet_c(ni,nj))
       kmu=0
       hu = 0.0
    end if

    select case (trim(topography))
    case('rectangular_basin')
       call rectangular_basin(Vgrid%zb, ht(1:ni,1:nj) )
       kmt(:,:) = size(Vgrid%zb)
    case('bowl')             
       call bowl(Hgrid%T%x(1:ni,1:nj), Hgrid%T%y(1:ni,1:nj), Vgrid%zb, ht(1:ni,1:nj) )  
    case('gaussian')         
       call gaussian(Hgrid%T%x(1:ni,1:nj), Hgrid%T%y(1:ni,1:nj), Vgrid%zb, ht(1:ni,1:nj) ) 
    case('dome')         
       call dome(Hgrid%T%x(1:ni,1:nj), Hgrid%T%y(1:ni,1:nj), Vgrid%zb, ht(1:ni,1:nj) ) 
    case('idealized')
       call idealized (Hgrid%T%x(1:ni,1), Hgrid%T%y(1,1:nj), Vgrid%zb, ht(1:ni,1:nj), kmt(1:ni,1:nj)  )
    case('all_land') 
       ht(:,:) = 0.0
       kmt(:,:) = 0
    case ( 'from_file' )
       if(interp_method == "conservative") then
          allocate(lonv(ni+1, nj+1), latv(ni+1, nj+1) )
          allocate(xv(ni,nj,4), yv(ni,nj,4) )
          call mpp_global_field(Hgrid%Domain, Hgrid%T%x_vert, xv)
          call mpp_global_field(Hgrid%Domain, Hgrid%T%y_vert, yv)
          lonv(1:ni,1:nj) = xv(1:ni,1:nj,1)
          lonv(ni+1,1:nj) = xv(ni,  1:nj,2)
          lonv(1:ni,nj+1) = xv(1:ni,  nj,3)
          lonv(ni+1,nj+1) = xv(ni,    nj,4)
          latv(1:ni,1:nj) = yv(1:ni,1:nj,1)
          latv(ni+1,1:nj) = yv(ni,  1:nj,2)
          latv(1:ni,nj+1) = yv(1:ni,  nj,3)
          latv(ni+1,nj+1) = yv(ni,    nj,4)
          call get_topog_from_file(lonv, latv, ht(1:ni,1:nj) )
          deallocate(xv, yv, lonv, latv)
       else
          call get_topog_from_file(Hgrid%T%x(1:ni,1:nj), Hgrid%T%y(1:ni,1:nj), ht(1:ni,1:nj) )
       end if
         ht = scale_factor*ht         ! flip sign or change units (user specified) 
         where (ht < 0) ht = 0        ! set elevation of land points to zero
    case default 
       call mpp_error(FATAL,'topog_mod: '//trim(topography)//' is not a valid option of nml "topography" ')
    end select

    if (filter_topog) call filter_topo(ht(1:ni,1:nj), num_filter_pass)

    !--- process topography if needed
    if(topog_depend_on_vgrid) then

       ! Compare "ht" to bottom of deepest model level. 
       if(debug) call show_deepest(Vgrid%zb, ht(1:ni,1:nj) )

       ! make first row of ocean model all land points for ice model
       if(fill_first_row .and. trim(topography) == "from_file" )  ht(:,1) = 0.0 
    
       !--- for idealized topog, do not need to process topography
       if(trim(topography) .ne. "all_land" .and. trim(topography) .ne. "idealized" ) then
          ! Discretize "kmt" based on "ht" and zw(1..nk). 
          allocate(dxte(ni,nj), dytn(ni,nj))
          call mpp_global_field(Hgrid%Domain, Hgrid%E%ds_01_21, dxte)
          call mpp_global_field(Hgrid%Domain, Hgrid%N%ds_10_12, dytn)
          call process_topo(ht, kmt, Vgrid%zb, dxte, dytn,    &
               Hgrid%T%x(1:ni,1:nj), Hgrid%T%y(1:ni,1:nj) , tripolar_grid, cyclic_x, cyclic_y)
          deallocate(dxte, dytn)
       endif
    endif
    if (have_obc) call prep_obc_topog(kmt, ht, Hgrid%ni, Hgrid%nj)

    if(output_corner_topog) then
       do j=1,nj
          do i=1,ni
             kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
             hu (i,j) = min(ht (i,j), ht (i+1,j), ht (i,j+1), ht (i+1,j+1))
          enddo
       enddo
    end if

    !--- define topography data -----------------------------------------
    Topog%depth_t(1:ni,1:nj)      = ht(1:ni,1:nj)
    if(topog_depend_on_vgrid) Topog%num_levels(1:ni,1:nj) = kmt(1:ni,1:nj)
    if(output_corner_topog) then
    Topog%depth_c(1:ni,1:nj)      = hu(1:ni,1:nj)
       if(topog_depend_on_vgrid) Topog%num_levels_c(1:ni,1:nj) = kmu(1:ni,1:nj)
    end if

    !--- define land/sea mask --------------------------------------
    Topog%wet = 0.0
    where(Topog%depth_t .gt. 0.0) Topog%wet = 1.0
    if(output_corner_topog) then
    Topog%wet_c = 0.0
    where(Topog%depth_c .gt. 0.0) Topog%wet_c = 1.0
    end if

    if(check_mask) call checking_mask(Topog%wet, cyclic_x, cyclic_y, tripolar_grid, have_obc)

    !--- chcksum ---------------------------------------------------
    if(debug) then
       pe = mpp_pe()
       write (stdout(),'(/(a,I20/))')'topog chksum: depth_t = ',   mpp_chksum(Topog%depth_t, (/pe/))
       if(topog_depend_on_vgrid) write (stdout(),'(/(a,I20/))')'topog chksum: num_levels = ', &
              mpp_chksum(Topog%num_levels, (/pe/))
       write (stdout(),'(/(a,I20/))')'topog chksum: wet = ',       mpp_chksum(Topog%wet, (/pe/))
       if(output_corner_topog) then
       write (stdout(),'(/(a,I20/))')'topog chksum: depth_c = ',   mpp_chksum(Topog%depth_c, (/pe/))
          if(topog_depend_on_vgrid) write (stdout(),'(/(a,I20/))')'topog chksum: num_levels_c = ', &
               mpp_chksum(Topog%num_levels_c, (/pe/))
       write (stdout(),'(/(a,I20/))')'topog chksum: wet_c = ',       mpp_chksum(Topog%wet_c, (/pe/))
    endif
        endif

    deallocate(ht, kmt)
    if(output_corner_topog) deallocate(hu, kmu)

    return

  end subroutine generate_topog

  !#######################################################################
  ! <SUBROUTINE NAME="write_topog_meta">

  !   <OVERVIEW>
  !     Write out topography meta data.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_topog_meta(unit, axis_x, axis_y)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="axis_x, axis_y" TYPE="type(axistype)">
  !     axis of T-cell center
  !   </IN>
  ! </SUBROUTINE>


  subroutine write_topog_global_meta(file)
    character(len=*), intent(in) :: file    
    real                         :: ttmp
    integer                      :: unit
        
    if(mpp_pe() .ne. mpp_root_pe() ) return

    !--- get the io unit associated with file
    unit = get_file_unit(file)

    call mpp_write_meta(unit,'topography', cval=trim(topography) )
    if(trim(topography) == 'from_file') then
       call mpp_write_meta(unit,'input_file',cval=trim(topog_file))
       call mpp_write_meta(unit,'input_field',cval=trim(topog_field))
    end if

    if (full_cell) call mpp_write_meta(unit,'full_cell',cval='y')
    if (fill_isolated_cells) call mpp_write_meta(unit,'fill_isolated_cells',cval='y')
    if (fill_first_row .and. trim(topography) == "from_file" )  &
             call mpp_write_meta(unit,'fill_first_row',cval='y')
    if (dont_change_landmask) call mpp_write_meta(unit,'dont_change_landmask',cval='y')
    if (fill_shallow) call mpp_write_meta(unit,'fill_shallow',cval='y')
    if (deepen_shallow) call mpp_write_meta(unit,'deepen_shallow',cval='y')
    if (round_shallow) call mpp_write_meta(unit,'round_shallow',cval='y')
    if (adjust_topo) call mpp_write_meta(unit,'adjust_topo',cval='y')
    if (filter_topog) then
       call mpp_write_meta(unit,'filter_topog',cval='y')
       ttmp=num_filter_pass
       call mpp_write_meta(unit,'num_filter_pass',rval=ttmp)
    endif

  end subroutine write_topog_global_meta

  !#####################################################################

  subroutine write_topog_field_meta(file)
    character(len=*),      intent(in) :: file

    if(mpp_pe() .ne. mpp_root_pe() ) return

    call write_field_meta(file,'depth_t','meters', 'topographic depth of T-cell', 2)
    if(output_corner_topog) call write_field_meta(file,'depth_c','meters', 'topographic depth of C-cell', 2, 'C', 'C')
    if(topog_depend_on_vgrid) &
    call write_field_meta(file,'num_levels', 'none', 'number of vertical T-cells', 2)
    if(output_corner_topog .and. topog_depend_on_vgrid) &
    call write_field_meta(file,'num_levels_c', 'none', 'number of vertical C-cells', 2, 'C', 'C')
    call write_field_meta(file, 'wet', 'none','land/sea flag (0=land) for T-cell', 2)
    if(output_corner_topog) call write_field_meta(file, 'wet_c', 'none','land/sea flag (0=land) for C-cell', 2, 'C', 'C')

    return
  end subroutine write_topog_field_meta

  !#######################################################################
  ! <SUBROUTINE NAME="write_topog_data">
  !   <OVERVIEW>
  !     write the topography data to netcdf file
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_topog_data (unit,Topog)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Topog" TYPE="type(topog_data_type)">
  !     A derived-type variable that contains topography data.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_topog_data(file, Topog)
    character(len=*),      intent(in) :: file
    type(topog_data_type), intent(in) :: Topog

    call write_field_data(file, 'depth_t', Topog%depth_t)
    if(output_corner_topog) call write_field_data(file, 'depth_c', Topog%depth_c)
    if(topog_depend_on_vgrid) call write_field_data(file, 'num_levels', Topog%num_levels)
    if(output_corner_topog .and. topog_depend_on_vgrid) call write_field_data(file, 'num_levels_c', Topog%num_levels_c)
    call write_field_data(file, 'wet', Topog%wet)
    if(output_corner_topog) call write_field_data(file, 'wet_c', Topog%wet_c)

    return
  end subroutine write_topog_data

  !#######################################################################
  ! <SUBROUTINE NAME="topog_end">

  !   <OVERVIEW>
  !     Destruction routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "topog_data_type" variables.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call topog_end ( Topog )
  !   </TEMPLATE>
  !   <INOUT NAME="Topog" TYPE="type(topog_data_type)">
  !     A derived-type variable that contains topography data.
  !   </INOUT>
  ! </SUBROUTINE>

  subroutine topog_end(Topog)
    type(topog_data_type), intent(inout) :: Topog

    deallocate( Topog%depth_t, Topog%num_levels, Topog%wet)
    if(output_corner_topog) deallocate( Topog%depth_c, Topog%num_levels_c, Topog%wet_c )
    module_is_initialized = .false.

    return

  end subroutine topog_end

  !#######################################################################
  !--- get the minimum depth of surrounding cells
  function hmin(ht, i,j)
    real, dimension(0:,0:), intent(in) :: ht
    integer,                intent(in) :: i, j
    real :: hmin

    hmin = min(ht(i,j), ht(i+1,j), ht(i,j+1), ht(i+1,j+1))
    return

  end function hmin

  !#######################################################################
  !--- get the minimum number of vertical levels of surrounding cells
  function kmin(kmt, i,j)
    integer, dimension(0:,0:), intent(in) :: kmt
    integer,                   intent(in) :: i, j
    integer :: kmin

    kmin = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
    return

  end function kmin

  !#######################################################################
  !--- analyze topographic slopes and write to standard output.
  subroutine analyze_topographic_slopes(dxte, dytn, geolon_t, geolat_t, ht, kmt, nk)
    real, dimension(:,:), intent(in)    :: dxte, dytn
    real, dimension(:,:), intent(in)    :: geolon_t, geolat_t
    real, dimension(:,:), intent(in)    :: ht
    integer, dimension(:,:), intent(in) :: kmt
    integer,              intent(in)    :: nk

    integer,dimension(:), allocatable :: slope
    integer :: i, j, k, nbin, it, iu, jt, ju, ni, nj
    real    :: sum_t, sum_u, slope_t, slope_u

    ni = size(ht,1); nj = size(ht,2)
    allocate (slope(nk))
    slope(:) = 0
    do j=1,nj-1
       do i=1,ni-1
          if (kmt(i+1,j) > 0 .and. kmt(i,j) > 0) then
             nbin = abs(kmt(i+1,j) - kmt(i,j)) + 1
             slope(nbin) = slope(nbin) + 1
          endif
          if (kmt(i,j+1) > 0 .and. kmt(i,j) > 0) then
             nbin = abs(kmt(i,j+1) - kmt(i,j)) + 1
             slope(nbin) = slope(nbin) + 1
          endif
       enddo
    enddo

    ! not checking the slopes at the boundaries here since it is unclear how
    ! to handle tripolar grid.  need to resolve !!

    write (stdout(),'(//,10x,a/)') 'How well does the specified resolution resolve Topography?'
    sum_t = 0.0; sum_u = 0.0
    do k=1,nk
       sum_t = sum_t + slope(k)
       sum_u = sum_u + slope(k)
    enddo

    if (full_cell) then
       write (stdout(),'(/a)') ' Note: Since full_cell = true, topography cannot be well &
            & resolved and this analysis is inappropriate.'
       write (stdout(),'(a/)') '       Set full_cell = false to see the analysis.'
       write (stdout(),'(/a)') ' Warning: full_cell = true is only meant for testing the partial cell algorithm.'
       write (stdout(),'(a/)') '          Its use is not recommended for physically realistic experiments.'
    else
       if (slope(1) /= 0) write (stdout(),'(1x,i8,a,f7.3,a)') slope(1), ' T-cell slopes (',100.0*slope(1)/sum_t,&
            '%) are resolved without any change in vertical levels.'
       if(nk .gt. 2) then
          if (slope(2) /= 0) write (stdout(),'(1x,i8,a,f7.3,a)') slope(2), ' T-cell slopes (',100.0*slope(2)/sum_t,&
               '%) are resolved by a change of 1 vertical level.'
       endif
       do k=3,nk
          if (slope(k) /= 0) write (stdout(),'(1x,i8,a,f7.3,a,i3,a)') slope(k), ' T-cell slopes (',100.0*slope(k)/sum_t,&
               '%) require a change of ',k-1,' vertical levels and are thus unresolved.'
       enddo
    endif
    deallocate (slope)

    ! compute maximum slopes

    slope_t = 0.0
    slope_u = 0.0
    iu=1; ju=1; it=1; jt=1
    ! these are approximate slopes where grids are unevenly spaced.
    ! fix this
    do j=1,nj-1
       do i=1,ni-1
          if (ht(i,j) /= 0.0 .and. ht(i+1,j) /= 0.0) then
             if (slope_t < abs((ht(i+1,j)-ht(i,j))/dxte(i,j))) then
                slope_t = abs((ht(i+1,j)-ht(i,j))/dxte(i,j))
                it = i; jt = j
             endif
          endif
          if (ht(i,j) /= 0.0 .and. ht(i,j+1) /= 0.0) then
             if (slope_t < abs((ht(i,j+1)-ht(i,j))/dytn(i,j))) then
                slope_t = abs((ht(i,j+1)-ht(i,j))/dytn(i,j))
                it = i; jt = j
             endif
          endif
       enddo
    enddo

    write (stdout(),'(/,a,f8.6,a,i4,a1,i4,a,f7.2,a,f7.2,a)') ' The maximum T-cell slope =',slope_t,&    
         ' at (i,j) = (',it,',',jt,'),  (lon,lat) = (',geolon_t(it,jt),',',geolat_t(it,jt),')'
  end subroutine analyze_topographic_slopes

  !#######################################################################
  !--- Constructing a rectangular basin with a flat bottom
  subroutine rectangular_basin(zw, ht)
    real, dimension(:),      intent(in) :: zw
    real, dimension(:,:), intent(inout) :: ht

    integer :: i, j, ni, nj, nk

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)

    write (stdout(),'(/a/)')' Constructing a rectangular basin with a flat bottom.'

    do j=1,nj
       do i=1,ni
          ht(i,j) = zw(nk)
       enddo
    enddo
  end subroutine rectangular_basin

  !#######################################################################
  !--- Constructing a gaussian bump
  subroutine gaussian(geolon_t, geolat_t, zw, ht)
    real, dimension(:,:),   intent(in)  :: geolon_t, geolat_t
    real, dimension(:),     intent(in)  :: zw
    real, dimension(:,:), intent(inout) :: ht

    real    :: bump_height, bump_scale, xcent, ycent, arg, bottom
    real    :: xe, xw, ys, yn
    integer :: i, j, ni, nj, nk

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)

    xw = geolon_t(1,1); ys = geolat_t(1,1); xe = xw; yn = ys
    do j=1,nj
       do i=1,ni
          xw = min(geolon_t(i,j),xw); xe = max(geolon_t(i,j),xe)
          ys = min(geolat_t(i,j),ys); yn = max(geolat_t(i,j),yn)
       enddo
    enddo

    bump_height = gauss_amp*zw(nk)
    bump_scale = gauss_scale*min(xe-xw, yn-ys)
    xcent = 0.5*(xe+xw)
    ycent = 0.5*(yn+ys)
    write (stdout(),'(/a,f8.2,a,f8.2,a)')&
         ' Constructing a gaussian bump of height =', bump_height,' meters with a scale width of ', bump_scale,' degrees.'
    write (stdout(),'(a,f6.2,a,f6.2,a)') ' The bump is centered at (lon,lat) = (',xcent,',',ycent,') deg.'
    write (stdout(),'(a,f8.2,a,f8.2,a/)') ' The ocean floor rises with a slope of ',slope_x, &
         ' meters/deg towards the east and ', slope_y,' meters/deg to the north.'
    do j=1,nj
       do i=1,ni
          arg = ((geolon_t(i,j)-xcent)**2 + (geolat_t(i,j) - ycent)**2)
          bottom = zw(nk) - bump_height*exp(-arg/bump_scale**2)
          bottom = bottom - slope_x*(geolon_t(i,j)-xw)- slope_y*(geolat_t(i,j)-ys)
          ht(i,j) = max(bottom,zw(2))  
       enddo
    enddo

    return

  end subroutine gaussian



  !#######################################################################
  !--- Construct a domain similar to the DOME configuration used in 
  !    Legg, Hallberg, and Girton (2005) Ocean Modelling.  
  !
  !
  !                       |     |
  !                       |     |
  !                      /|     |
  !                     / |     |
  !     up             /  |     |
  !     ^             /   |     |
  !     |            /    |     |
  !     |           /     |     |
  !                /      |     |
  !---------------/-------|-----|
  !     --->north 
  !
  !
  subroutine dome (geolon_t, geolat_t, zw, ht)
    real, dimension(:,:),   intent(in)  :: geolon_t, geolat_t
    real, dimension(:),     intent(in)  :: zw
    real, dimension(:,:), intent(inout) :: ht

    real    :: xe, xw, ys, yn
    real    :: xe_embay, xw_embay, yn_embay, ys_embay
    real    :: yn_slope, ys_slope
    real    :: height 
    integer :: inear, jnear
    integer :: i, j, ni, nj, nk

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)

    ! determine full domain extents 
    xw = geolon_t(1,1); ys = geolat_t(1,1); xe = xw; yn = ys
    do j=1,nj
       do i=1,ni
          xw = min(geolon_t(i,j),xw); xe = max(geolon_t(i,j),xe)
          ys = min(geolat_t(i,j),ys); yn = max(geolat_t(i,j),yn)
       enddo
    enddo

    ! define embayment location
    xw_embay = dome_embayment_west 
    xe_embay = dome_embayment_east 
    yn_embay = yn
    ys_embay = dome_embayment_south

    ! find grid (lat,lon) corresponding to embayment specifications 
    do j=1,1 
      do i=1,ni-1
         if(geolon_t(i,j) <= xw_embay .and. geolon_t(i+1,j) >= xw_embay) then 
           xw_embay = geolon_t(i,j)
         endif 
         if(geolon_t(i,j) <= xe_embay .and. geolon_t(i+1,j) >= xe_embay) then 
           xe_embay = geolon_t(i,j)
         endif 
      enddo
    enddo
    do j=1,nj-1
      do i=1,1
         if(geolat_t(i,j) <= ys_embay .and. geolat_t(i,j+1) >= ys_embay) then 
           ys_embay = geolat_t(i,j)
         endif 
      enddo
    enddo

    ! define latitudes sloping surface starts and finishes 
    yn_slope = ys_embay
    ys_slope = yn_slope - (dome_bottom-dome_embayment_depth)/(dome_slope*deg2metre) 
    do j=1,nj-1
      do i=1,1
         if(geolat_t(i,j) <= ys_slope .and. geolat_t(i,j+1) >= ys_slope) then 
           ys_slope = geolat_t(i,j)
         endif 
      enddo
    enddo

    write (stdout(),'(/a)')     ' Constructing DOME configuration topography'
    write(stdout(),'(/a,f10.4)')' western boundary (longitude) of full domain = ',xw
    write(stdout(),'(a,f10.4)') ' eastern boundary (longitude) of full domain = ',xe
    write(stdout(),'(a,f10.4)') ' southern boundary (latitude) of full domain = ',ys
    write(stdout(),'(a,f10.4)') ' northern boundary (latitude) of full domain = ',yn

    write(stdout(),'(a,f10.4)') ' depth (m) of the embayment                  = ',dome_embayment_depth
    write(stdout(),'(a,f10.4)') ' depth (m) of the domain bottom              = ',dome_bottom
    write(stdout(),'(a,f10.4)') ' western boundary (longitude) of embayment   = ',xw_embay 
    write(stdout(),'(a,f10.4)') ' eastern boundary (longitude) of embayment   = ',xe_embay 
    write(stdout(),'(a,f10.4)') ' southern boundary (latitude) of embayment   = ',ys_embay
    write(stdout(),'(a,f10.4)') ' northern boundary (latitude) of embayment   = ',yn_embay
 
    write(stdout(),'(/a,f10.4)')' slope of the dome configuration             = ',dome_slope
    write(stdout(),'(a,f10.4)') ' northern boundary (latitude) of slope       = ',yn_slope 
    write(stdout(),'(a,f10.4)') ' southern boundary (latitude) of slope       = ',ys_slope 


    ! define flat bottom
    do j=1,nj
       do i=1,ni
          ht(i,j) = dome_bottom 
       enddo
    enddo

    ! define sloping bottom 
    do j=1,nj
       do i=1,ni
          if(geolat_t(i,j) >= ys_slope .and. geolat_t(i,j) < yn_slope) then 
              height  = dome_slope*deg2metre*(geolat_t(i,j)-ys_slope)
              ht(i,j) = dome_bottom - height
          endif
       enddo
    enddo
    
    ! define shelf
    do j=1,nj
       do i=1,ni
          if(geolat_t(i,j) >= yn_slope) then 
              ht(i,j) = 0.0
          endif
       enddo
    enddo

    ! define embayment 
    do j=1,nj
       do i=1,ni
          if(geolon_t(i,j) >= xw_embay .and. geolon_t(i,j) <= xe_embay .and. &
             geolat_t(i,j) >= ys_embay .and. geolat_t(i,j) <= yn_embay) then
             ht(i,j) = dome_embayment_depth
          endif
       enddo
    enddo
  

    return

  end subroutine dome



  !#######################################################################
  ! From "Simulation of density-driven frictional downslope flow in z-coordinate mocean models"
  ! Winton et al. JPO, Vol 28, No 11, 2163-2174, November 1998
  subroutine bowl(geolon_t, geolat_t, zw, ht)
    real, dimension(:,:),   intent(in)  :: geolon_t, geolat_t
    real, dimension(:),     intent(in)  :: zw
    real, dimension(:,:), intent(inout) :: ht

    real    :: bottom
    integer :: i, j, ni, nj

    ni = size(ht,1); nj = size(ht,2)

    write (stdout(),'(/a)') ' Constructing a Winton bowl.'
    do j=1,nj
       do i=1,ni
          if (geolon_t(i,j) <= bowl_west .or. geolon_t(i,j) >= bowl_east &
               .or. geolat_t(i,j) <= bowl_south .or. geolat_t(i,j) >= bowl_north) then           
             bottom = bowl_min_depth
          else
             bottom = bowl_min_depth + bowl_max_depth*(1.0-exp(-((geolat_t(i,j)-bowl_south)/2.0)**2))&
                  *(1.0-exp(-((geolat_t(i,j)-bowl_north)/2.0)**2))&
                  *(1.0-exp(-((geolon_t(i,j)-bowl_west)/4.0)**2))&
                  *(1.0-exp(-((geolon_t(i,j)-bowl_east)/4.0)**2))
          endif
          ht(i,j) = max(bottom,zw(2))  
       enddo
    enddo
  end subroutine bowl

  !#######################################################################
  ! limit the minimum number of levels. kmt_min should be >= 2
  subroutine enforce_min_depth(zw, ht, kmt)
    real, dimension(:),         intent(in) :: zw 
    real, dimension(:,:),    intent(inout) :: ht
    integer, dimension(:,:), intent(inout) :: kmt

    integer :: n, i, j, kmt_shallow, ni, nj

    ni = size(ht,1); nj = size(ht,2)

    if(debug) write (stdout(),'(/a,i2)') ' Enforcing the minimum number of ocean cells in the vertical to be ', kmt_min
    if(size(zw) .lt. kmt_min) call mpp_error(FATAL,'topog_mod: number of vertical cells= '// &
         trim(string(size(zw)))//' is less than nml "kmt_min"= '//trim(string(kmt_min)) ) 
    n = 0
    do i=1,ni
       do j=1,nj
          if (kmt(i,j) /= 0 .and. kmt(i,j) < kmt_min) then
             n = n + 1
             kmt_shallow = kmt(i,j)
             if (round_shallow .or. (.not. deepen_shallow .and. .not. fill_shallow)) then
                if (zw(kmt(i,j)) < 0.5*zw(kmt_min)) then
                   if (debug) write(stdout(),'(a,3i6,a)') 'Making location i,j,kmt= ',i,j,kmt(i,j),' land'
                   kmt(i,j) = 0
                   ht(i,j)  = 0.0
                else
                   if (debug) write(stdout(),'(a,3i6,a)') 'Setting location i,j,kmt= ', &
                        i,j,kmt(i,j),'  to minimum ocean depth'
                   kmt(i,j) = kmt_min
                   ht(i,j)  = zw(kmt_min)
                end if
             endif
             if (fill_shallow) then
                if (debug) write(stdout(),'(a,3i6,a)') 'Making location i,j,kmt= ',i,j,kmt(i,j),' land'
                kmt(i,j) = 0
                ht(i,j)  = 0.0
             endif
             if (deepen_shallow) then
                if (debug) write(stdout(),'(a,3i6,a)') 'Setting location i,j,kmt= ',i,j,kmt(i,j),'  to minimum ocean depth'
                kmt(i,j) = kmt_min
                ht(i,j)  = zw(kmt_min)
             endif

          endif
       enddo
    enddo

    if(debug) then
    if (n > 0) then
          write (stdout(),'(a,i5,a/)')  '  -> Modified', n,' shallow cells'
    else
          write (stdout(),'(a/)') '  -> No modifications needed'
       endif
    endif

  end subroutine enforce_min_depth

  !#######################################################################
  !--- remove isolated cells
  subroutine remove_isolated_cells(ht, kmt, tripolar_grid)
    real, dimension(0:,0:),    intent(inout) :: ht
    integer, dimension(0:,0:), intent(inout) :: kmt
    logical,                   intent(in   ) :: tripolar_grid

    integer :: i, j, k, n, ni, nj, nj_max
    real    :: tmp

    ni = size(ht,1) - 2; nj = size(ht,2) - 2
    ! when tripolar_grid is true, do not check at j = nj (folded region).
    ! Will upgrade it if needed.
    if (tripolar_grid) then
      nj_max = nj-1
    else
      nj_max = nj
    endif

    if (fill_isolated_cells) then
       ! fill isolated cells (potholes and trenches) at all levels in kmt.
       ! filled kmt array is the maximum of the surrounding kmt levels.
    n=0
       if(debug) write (stdout(),'(a)') ' Searching for isolated T cells...'
        do j=1,nj_max
          do i=1,ni
            k = max(kmin(kmt,i,j), kmin(kmt,i-1,j), kmin(kmt,i,j-1), kmin(kmt,i-1,j-1))
             if ((k > 0) .or. .not. dont_change_landmask) then
                if (kmt(i,j) /= k) then
                   n = n + 1
              tmp = max(hmin(ht,i,j), hmin(ht,i-1,j), hmin(ht,i,j-1),hmin(ht,i-1,j-1))
                   if (debug) write(stdout(),'(a,3i6,f10.3,a,i4,f10.3)') 'Resetting location i,j,kmt,ht= ', &
                        i,j,kmt(i,j),ht(i,j),'  to ',  k,tmp
                   ht(i,j) = tmp
                kmt(i,j) = k
              endif
            endif
          enddo
        enddo
       if(debug) then
          if (n > 0) then
             write (stdout(),*) ' -> Found ',n,' and filled them in.'
        else
          write (stdout(),*) ' -> None Found.'
        endif
       endif
      else
       if(debug) write (stdout(),'(/a)') 'Note: Not filling isolated T cells...'
    endif

  end subroutine remove_isolated_cells

  !#######################################################################  
  !--- Restricting partial cells
  subroutine restrict_partial_cells(zw, ht, kmt)
    real, dimension(:), intent(in)      :: zw  
    real, dimension(:,:), intent(inout) :: ht
    integer, dimension(:,:), intent(inout) :: kmt

    real, dimension(size(zw)) :: p_cell_min
    real    :: tmp, tmp1, tmp2
    integer :: i, j, n, k, ni, nj, nk, itmp

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)
    n = 0

    if (fraction_full_cell==0.0) then 
        do k=1,nk
           p_cell_min(k) = min(50.0,zw(1))
      enddo
    else
      p_cell_min(1) = fraction_full_cell*zw(1)
      if (nk >= 2) then
        do k=2,nk
           p_cell_min(k) = max(fraction_full_cell*(zw(k)-zw(k-1)),min_thickness) 
        enddo
      endif
    endif

    if(debug) then 
      write (stdout(),'(/a)') 'Partial cell restriction' 
    do k=1,nk
        write (stdout(),'(a,i6,a,e12.6,a)') ' Restricting partial cell # ',k,' to a min thickness of ', p_cell_min(k),' m.'
    enddo
    endif 

    if( open_very_this_cell  ) then
       do j=1,nj       
          do i=1,ni
             k = kmt(i,j)
             if (k > 1) then
                ! Use a "fuzzy" small epsilon value in comparison, to prevent
                ! single-precision comparison from re-editing cells
                ! that are "exactly" at the minimum partial cell thickness.
                if ((ht(i,j)-zw(k-1))*(1.+small) < p_cell_min(k)) then
                   tmp = zw(k-1) + p_cell_min(k)
                   if (debug) write(stdout(),'(a,2i4,f12.5,a,f10.5)') 'Resetting depth i,j,ht= ',i,j,ht(i,j),'  to ',&
                        tmp
                   ht(i,j) = tmp
                   n = n + 1
                endif
             endif
          enddo
       enddo
    else
    do j=1,nj       
      do i=1,ni
        k = kmt(i,j)
        if (k > 1) then
          ! Use a "fuzzy" small epsilon value in comparison, to prevent
          ! single-precision comparison from re-editing cells
          ! that are "exactly" at the minimum partial cell thickness.
          if ((ht(i,j)-zw(k-1))*(1.+small) < p_cell_min(k)) then
                   tmp1 = ht(i,j)-zw(k-1)
                   tmp2 = zw(k-1) + p_cell_min(k) - ht(i,j)
                   if (tmp2 .gt. tmp1) then
                      itmp = kmt(i,j) -1 
                      tmp = zw(k-1)
                      if (debug) write(stdout(),'(a,2i4,i4,a,i4)') 'Resetting kmt i,j,kmt= ',i,j,kmt(i,j),'  to ',&
                           itmp
                      kmt(i,j) = itmp
                   else  
             tmp = zw(k-1) + p_cell_min(k)
                   endif
                   if (debug) write(stdout(),'(a,2i4,f12.5,a,f10.5)') 'Resetting depth i,j,ht= ',i,j,ht(i,j),'  to ',&
                        tmp
             ht(i,j) = tmp
             n = n + 1
          endif
        endif
      enddo
    enddo
    end if
    if(debug) then
    if (n > 0) then
          write (stdout(),*) ' -> Found ',n, ' cells with too thin partical cell thickness, and so reset depth ht for these cells.'
    else
          write (stdout(),*) ' -> No cells were found whose partial cells were too thin.'
       endif
    endif
  end subroutine restrict_partial_cells

  !#######################################################################
    !     smooth topographic depth "d" with "num_pass" applications of a 2D
    !     version of a shapiro filter (weights = 1/4, 1/2, 1/4) . 
    !     allow filtering to decrease the bottom depth but not increase it.
    !     do not allow original geometry to change.
    !     note: depth "d" should be on a grid of uniformly constant spacing
  subroutine filter_topo (d, num_pass)
    real, dimension(:,:),  intent(inout) :: d
    integer,                  intent(in) :: num_pass
    real, dimension(size(d,1),size(d,2)) :: rmask, s
    real,           dimension(-1:1,-1:1) :: f
    integer                              :: n, i, j, ip, jp, ni, nj
    real                                 :: d_old

    ni = size(d,1); nj = size(d,2)
    ! 2D symmetric filter weights

    f(-1,-1) = 1.0/16.0
    f( 0,-1) = 1.0/8.0
    f( 1,-1) = 1.0/16.0
    f(-1, 0) = 1.0/8.0
    f( 0, 0) = 1.0/4.0
    f( 1, 0) = 1.0/8.0
    f(-1, 1) = 1.0/16.0
    f( 0, 1) = 1.0/8.0
    f( 1, 1) = 1.0/16.0

    ! geometry mask

    where (d(:,:) == 0.0)
       rmask(:,:) = 0.0
    elsewhere
       rmask(:,:) = 1.0
    endwhere

    s=d

    do n=1,num_pass
       do j=2,nj-1
          !      if (j >= js .and. j <= je) then
          do i=2,ni-1
             s(i,j) = 0.0
             d_old  = d(i,j)
               do ip=-1,1
                  do jp=-1,1
                   if (rmask(i+ip,j+jp) .eq. 0.0) then
                        s(i,j) = s(i,j) + d(i,j)*f(ip,jp)
                     else
                        s(i,j) = s(i,j) + d(i+ip,j+jp)*f(ip,jp)
                     endif
                  enddo
               enddo
             if (.not. smooth_topo_allow_deepening) then
                if (s(i,j) .gt. d_old) then
                   s(i,j) = d_old
                endif
             endif
          enddo
          do i=2,ni-1
             s(i,j) = s(i,j)*rmask(i,j)
          enddo
          !      endif
       enddo

            d(:,:) = s(:,:)
    enddo
  end subroutine filter_topo

  !#######################################################################
  !--- print out deepest topography
  subroutine show_deepest(zw, ht )
    real, dimension(:),   intent(in) :: zw
    real, dimension(:,:), intent(in) :: ht

    integer :: i, j, ni, nj, nk
    real    :: deepest

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)
    deepest = 0.0
    do j=1,nj
       do i=1,ni
          if (ht(i,j) /= 0.0) then
             deepest  = max(ht(i,j),deepest)
          endif
       enddo
    enddo

    if ((deepest - zw(nk)) > 1.0) then
       write (stdout(),'(/a,f8.1,a,f8.1,a/a/)')'  Warning: Topography reaches to a depth of ',deepest,    &
            ' m. The deepest model level only reaches ', zw(nk), ' m.     Re-think the vertical grid  '// &
            'specification if the idea is to accurately capture the specified topography.'
    elseif (nk .gt. 1 .and. deepest <= zw(nk-1)) then
       write (stdout(),'(/a,f8.1,a,f8.1,a,i3,a/a)')&
            '   Warning: Topography reaches to a depth of ',deepest,' m. The deepest model level reaches ', &
            zw(nk),' m.            Fewer than ',nk,' vertical levels are required. '//                      &
            'The specified number of vertical levels is wasteful.'
    else
       write (stdout(),'(/a/)')&
            ' Note: The vertical grid specification contains the correct number of &
            &levels to efficiently contain the specified topography.'
    endif

  end subroutine show_deepest

  !#######################################################################
  !--- processing topography
  subroutine process_topo(ht, kmt, zw, dxte, dytn, geolon_t, geolat_t, tripolar_grid, cyclic_x, cyclic_y)

    real, dimension(0:,0:),    intent(inout) :: ht
    integer, dimension(0:,0:), intent(inout) :: kmt
    logical,                      intent(in) :: tripolar_grid, cyclic_x, cyclic_y
    real, dimension(:),           intent(in) :: zw
    real, dimension(:,:),         intent(in) :: dxte, dytn
    real, dimension(:,:),         intent(in) :: geolon_t, geolat_t

    integer :: i,j, ni, nj, nk
    real    :: ht_prev

    ni = size(ht,1)-2; nj = size(ht,2)-2; nk = size(zw)
    do j=1,nj
      do i=1,ni
        if (ht(i,j) >  0.0) then
          kmt(i,j) = nearest_index (ht(i,j), zw(1:nk))
          if (zw(kmt(i,j)) < ht(i,j)) then
                if ((ht(i,j)-zw(kmt(i,j))) < min_thickness) then
                   ht_prev = ht(i,j)
                   ht(i,j) = zw(kmt(i,j))
                   write (stdout(),*) ' Warning: very thin partial cell at ',i,j, &
                       ' (possibly because of netcdf-accuracy) is changed from ',ht_prev,' to ',ht(i,j)
                   write (stdout(),*) ' If this is not wanted, set "min_thickness" to zero in topog_nml'
                else
            if (kmt(i,j) == nk) then
              ht(i,j) = zw(kmt(i,j))
            else
              kmt(i,j) = kmt(i,j) + 1
            endif
          endif
             endif
        endif
      enddo
    enddo

    if (full_cell) then
       write (stdout(),'(/a/a/)')' Warning: Replacing partial bottom cells with full cell thicknesses.'
       do j=1,nj
          do i=1,ni
             if (kmt(i,j) /= 0) then
                ht(i,j) = zw(kmt(i,j))
             endif
          enddo
       enddo
    endif

    if (flat_bottom) then 
       write (stdout(),'(/a/a/)')' Warning: Replacing the ocean topography with a flat bottom where kmt(i,j)=nk.'
       do j=1,nj
          do i=1,ni
             if (kmt(i,j) /= 0) then
                kmt(i,j) = nk
                ht(i,j)  = zw(nk)
             endif
          enddo
       enddo
    endif

   if (tripolar_grid) then
       do i=1,ni
          kmt(i,nj+1) = kmt(ni-i+1,nj)
          ht(i,nj+1)  = ht(ni-i+1,nj)
       enddo
    else if(cyclic_y) then
       kmt(:,0)    = kmt(:,nj)
       kmt(:,nj+1) = kmt(:,1)
       ht(:,0)     = ht(:,nj)
       ht(:,nj+1)  = ht(:,1)       
    endif

   if (cyclic_x) then
       kmt(0,:)    = kmt(ni,:)
       kmt(ni+1,:) = kmt(1,:)
       ht(0,:)     = ht(ni,:)
       ht(ni+1,:)  = ht(1,:)
    endif

    if(debug)then
       write (stdout(),*) ' Checksum: Original ht =', sum( ht(1:ni, 1:nj) )
       write (stdout(),*) ' Checksum: Original kmt=', sum (kmt(1:ni, 1:nj) )
    endif

    if (adjust_topo) then
       call remove_isolated_cells(ht, kmt, tripolar_grid)
       call restrict_partial_cells(zw, ht(1:ni,1:nj), kmt(1:ni,1:nj) )
       call enforce_min_depth(zw, ht(1:ni,1:nj), kmt(1:ni,1:nj))
       if(debug) then
          write (stdout(),*) ' Checksum: Final kmt=', sum(kmt(1:ni, 1:nj))
          write (stdout(),*) ' Checksum: Final  ht=', sum( ht(1:ni, 1:nj))
       endif
    endif

    ! do some analysis.
    if(debug) call analyze_topographic_slopes(dxte, dytn, geolon_t, geolat_t, ht(1:ni,1:nj), kmt(1:ni,1:nj), nk) 

  end subroutine process_topo

  !#######################################################################
   ! construct a highly "idealized" world ... piece by piece
    !
    ! note: the purpose of this geometry/topography is to automatically
    !       map into arbitrary resolution as grid dimensions are
    !       changed, thereby facilitating the implementation
    !       and verification of the model on various computer platforms
    !       without referencing the topographic database.  Although it
    !       somewhat resembles the real world, it is NOT realistic.
    ! Note: this routine needs to be re-thought for generalized curvilinear coordinates
  subroutine idealized (xt,yt, zw, ht, kmt)
    real, dimension(:),    intent(in)      :: xt
    real, dimension(:),    intent(in)      :: yt
    real, dimension(:), intent(in)         :: zw
    real, dimension(:,:), intent(inout)    :: ht
    integer, dimension(:,:), intent(inout) :: kmt

    integer :: nj2, ni2, i, j, level, ni, nj, nk
    real    :: bot, arg

    ni = size(ht,1); nj = size(ht,2); nk = size(zw)
    kmt(:,:) = nk

    ! antarctica

    call setkmt (kmt(1:ni,1:nj), xt, yt, -90.0, 0.0, 360.0, -80.0, 0.0, 360.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -80.0, 360.0-25.0, 360.0, -70.0, 360.0, 360.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -80.0, 0.0, 360.0, -70.0, 0.0, 170.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -80.0, 360.0-135.0, 360.0-60.0, -68.0, 360.0-75.0, 360.0-60.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -70.0, 0.0, 155.0, -67.0, 50.0, 145.0, 0)

    ! australia

    call setkmt (kmt(1:ni,1:nj), xt, yt, -35.0, 116.0, 120.0, -31.0, 114.0, 130.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -38.0, 140.0, 151.0, -31.0, 130.0, 151.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -31.0, 115.0, 153.0, -20.0, 113.0, 149.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -20.0, 113.0, 149.0, -11.0, 131.0, 143.0, 0)

    ! south america

    call setkmt (kmt(1:ni,1:nj), xt, yt, -50.0, 360.0-74.0, 360.0-68.0, -40.0, 360.0-73.0, 360.0-62.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -40.0, 360.0-73.0, 360.0-62.0, -20.0, 360.0-70.0, 360.0-40.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -20.0, 360.0-70.0, 360.0-40.0, -16.0, 360.0-81.0, 360.0-35.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -16.0, 360.0-81.0, 360.0-35.0, 0.0, 360.0-80.0, 360.0-50.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 0.0, 360.0-80.0, 360.0-50.0, 11.0, 360.0-75.0, 360.0-60.0, 0)

    ! central america

    call setkmt (kmt(1:ni,1:nj), xt, yt, 6.0, 360.0-78.0, 360.0-75.0, 20.0, 360.0-105.0, 360.0-97.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 20.0, 360.0-105.0, 360.0-97.0, 30.0, 360.0-115.0, 360.0-94.0, 0)

    ! north america

    call setkmt (kmt(1:ni,1:nj), xt, yt, 25.0, 360.0-82.0, 360.0-80.0, 30.0, 360.0-85.0, 360.0-81.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 30.0, 360.0-115.0, 360.0-80.0, 40.0, 360.0-124.0, 360.0-74.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 40.0, 360.0-124.0, 360.0-74.0, 50.0, 360.0-124.0, 360.0-57.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 50.0, 360.0-124.0, 360.0-57.0, 60.0, 360.0-140.0, 360.0-64.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 60.0, 360.0-165.0, 360.0-64.0, 65.0, 360.0-140.0, 360.0-64.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 65.0, 360.0-140.0, 360.0-64.0, 70.0, 360.0-162.0, 360.0-72.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 70.0, 360.0-162.0, 360.0-140.0, 72.0, 360.0-157.0, 360.0-157.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 70.0, 360.0-130.0, 360.0-70.0, 75.0, 360.0-120.0, 360.0-80.0, 0)

    ! greenland
          
    call setkmt (kmt(1:ni,1:nj), xt, yt, 60.0, 360.0-45.0, 360.0-45.0, 75.0, 360.0-58.0, 360.0-19.0, 0)

    ! africa

    call setkmt (kmt(1:ni,1:nj), xt, yt, -35.0, 19.0, 28.0, 6.0, 8.0, 50.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 6.0, 0.0, 50.0, 18.0, 0.0, 56.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 18.0, 0.0, 56.0, 26.0, 0.0, 59.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 6.0, 360.0-10.0, 360.0, 18.0, 360.0-18.0, 360.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 18.0, 360.0-18.0, 360.0, 26.0,  360.0-15.0, 360.0, 0)

    ! northern africa and europe and asia
    
    call setkmt (kmt(1:ni,1:nj), xt, yt, 26.0, 360.0-15.0, 360.0, 40.0, 360.0-7.0, 360.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 40.0, 360.0-7.0, 360.0, 50.0, 360.0, 360.0, 0)
    
    call setkmt (kmt(1:ni,1:nj), xt, yt, 8.0, 77.0, 78.0, 26.0, 65.0, 90.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 4.0, 99.0, 100.0, 26.0, 90.0, 115.0, 0)

    call setkmt (kmt(1:ni,1:nj), xt, yt, 26.0, 0.0, 126.0, 40.0, 0.0, 122.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 40.0, 0.0, 130.0, 50.0, 0.0, 140.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 50.0, 0.0, 140.0, 60.0, 8.0, 140.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 60.0, 8.0, 163.0, 65.0, 13.0, 180.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 65.0, 13.0, 188.0, 70.0, 20.0, 180.0, 0)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 70.0, 70.0, 180.0, 75.0, 90.0, 100.0, 0)

    ! add an "idealized" undulating topography

    bot = zw(nk)
    nj2 = nj+2
    ni2 = ni+2
      do j=1,nj
        do i=1,ni
          if (kmt(i,j) .ne. 0) then
             arg = bot*(1-0.4*abs(cos(((j+1)*pi)/nj2)*sin(((i+1)*2*pi)/ni2)))
             kmt(i,j) = nearest_index (arg, zw)
          endif
        enddo
      enddo
      
    ! add "idealized" ridges

    level = nearest_index (0.666*zw(nk), zw)

    call setkmt (kmt(1:ni,1:nj), xt, yt, -20.0, 360.0-20.0, 360.0-10.0, 30.0, 360.0-45.0, 360.0-35.0, level)
    call setkmt (kmt(1:ni,1:nj), xt, yt, 30.0, 360.0-45.0, 360.0-35.0, 60.0, 360.0-20.0,  360.0-30.0, level)
    call setkmt (kmt(1:ni,1:nj), xt, yt, -60.0,360.0-100.0, 360.0-130.0, 40.0, 360.0-160.0, 180.0, level)

    level = nearest_index (0.5*zw(nk), zw)

    call setkmt (kmt(1:ni,1:nj), xt, yt, -50.0, 360.0-120.0, 360.0-120.0, 30.0, 190.0, 190.0, level)

    ! set ht to full cell depth.
    
    do j=1,nj
        do i=1,ni
          if (kmt(i,j) .ne. 0) then
             ht(i,j) = zw(kmt(i,j))
          else
            ht(i,j)  = 0.0
          endif
        enddo
      enddo

  end subroutine idealized

  !#####################################################################
    !     set the topography mask "kmt(i,j)" = "num" within the area of
    !     the trapezoid bounded by vertices:
    !     (alat1,slon1), (alat1,elon1), (alat2,slon2), and (alat2,elon2)
    !
    !     inputs:
    !
    !     xt = longitudes of T points in degrees
    !     yt = latitudes of T points in degrees
    !     alat1 = southern latitude of trapezoid (degrees)
    !     slon1 = starting longitude of southern edge of trapezoid (deg)
    !     elon1 = ending longitude of southern edge of trapezoid (deg)
    !     alat2 = northern latitude of trapezoid (degrees)
    !     slon2 = starting longitude of northern edge of trapezoid (deg)
    !     elon2 = ending longitude of northern edge of trapezoid (deg)
    !     num   = number of vertical levels
    !
    !     outputs:
    !
    !     kmt   = mask of vertical levels on model T points
  subroutine setkmt (kmt, xt, yt, alat1, slon1, elon1, alat2, slon2, elon2, num)
    real, intent(in),        dimension(:)  :: xt
    real, intent(in),        dimension(:)  :: yt
    integer, intent(inout), dimension(:,:) :: kmt 
    real,                       intent(in) :: alat1, slon1, elon1, alat2, slon2, elon2
    integer,                    intent(in) :: num
    integer                                :: j1, j2, js, je, i1, i2, is1, ie1
    integer                                :: is2, ie2, is, ie, i, j
    real                                   :: rdj

    j1 = nearest_index (alat1, yt )
    j2 = nearest_index (alat2, yt )
    js = min (j1,j2)
    je = max (j1,j2)

    i1  = nearest_index (slon1, xt)
    i2  = nearest_index (elon1, xt)
    is1 = min (i1,i2)
    ie1 = max (i1,i2)

    i1  = nearest_index (slon2, xt )
    i2  = nearest_index (elon2, xt )
    is2 = min (i1,i2)
    ie2 = max (i1,i2)

    is = is1
    ie = ie1

    ! fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
    ! the nudging of 1.e-5 is to insure that the test case resolution
    ! generates the same topography and geometry on various computers.

    if (js .eq. je) then
       rdj = 1.0
      else
       rdj = 1.0/(je-js)
      endif 
    do j=js,je
       do i=is,ie
          kmt(i,j) = num
    enddo
       is = nint(rdj*((j-js)*is2 + (je-j)*is1) + 1.0e-5)
       ie = nint(rdj*((j-js)*ie2 + (je-j)*ie1) + 1.0e-5)
    enddo
  end subroutine setkmt

  !#######################################################################
  !--- reading data from source data file topog_file and remap it onto current grid
  subroutine get_topog_from_file(lon_dst, lat_dst, data_dst)
    real,dimension(:,:), intent(in)  :: lon_dst, lat_dst
    real,dimension(:,:), intent(out) :: data_dst

    !--- local variables -------------------------------------------------
    integer           :: isc, iec, jsc, jec, layout(2), ndivs, unit, ndim, nvar, natt, ntime, n, len
    integer           :: i, j, nlon_src, nlat_src, nlon_dst, nlat_dst
    real              :: missing, D2R
    logical           :: do_remap = .true.
    logical           :: found_lon, found_lat
    character(len=1)  :: cart
    character(len=64) :: name
    type(domain2D)    :: domain  
    type(fieldtype)   :: data_field
    real, dimension(:,:), allocatable          :: tmp, global_tmp
    real, dimension(:),   allocatable          :: lon, lat, lonb, latb
    real, dimension(:,:), allocatable          :: data_src, mask_src, lon_src, lat_src
    type(fieldtype), allocatable, dimension(:) :: fields
    type(axistype), allocatable, dimension(:)  :: axes
    real  :: min_lat, max_lat
    integer :: jstart, jend

    !-------------------------------------------------------------------
    D2R = PI/180.
    layout = (/ 0, 0 /)

    !--- First read data from topog_file -------------------------------
    if(.not. file_exist(trim(topog_file))) &
         call mpp_error(FATAL,'topog_mod: file '//trim(topog_file)//' does not exist')
    call mpp_open(unit, trim(topog_file), action=MPP_RDONLY, form=MPP_NETCDF, &
         threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit,fields)

    write (stdout(),'(//,10x,a/)') 'Reading topography from file'

    nlon_src=0; nlat_src=0
    do n=1,nvar
       call mpp_get_atts(fields(n), name=name)
       if( trim(lowercase(name)) .eq. trim(lowercase(topog_field))) then
          call mpp_get_atts(fields(n), ndim=ndim)
          allocate(axes(ndim))
          call mpp_get_atts(fields(n), axes=axes)
          do i=1,ndim
             call get_axis_cart(axes(i),cart)
             call mpp_get_atts(axes(i),len=len)
             select case(cart)
             case('X')
                nlon_src = len
                allocate(lon(nlon_src) )
                call mpp_get_axis_data(axes(i),lon)
             case('Y')
                nlat_src = len
                allocate(lat(nlat_src))
                call mpp_get_axis_data(axes(i),lat)
             case('N')
                call mpp_error(FATAL,'topog_mod: axis cart of field '//trim(topog_field) &
                         //' is not correct, check file '//trim(topog_file)//' to make sure field ' &
                         //trim(topog_field)// ' has suitable attribute units or cartesian_axis')
             end select
          enddo
          data_field = fields(n)
       endif
    enddo

    if (nlon_src == 0 .or. nlat_src == 0 ) call mpp_error(FATAL,'topog_mod: field '//trim(topog_field)// &
         ' not found in file '//trim(topog_file) )

    !--- we assume the source grid is uniform in the zonal direction
    nlon_dst = size(data_dst,1)
    nlat_dst = size(data_dst,2)
    if( trim(interp_method) == 'conservative' ) then
       allocate(lonb(nlon_src+1), latb(nlat_src+1) )
       do i = 2, nlon_src
          lonb(i) = (lon(i-1) + lon(i))*0.5;
       end do
       lonb(1) = 2*lon(1) - lonb(2)
       lonb(nlon_src+1) = 2*lon(nlon_src) - lonb(nlon_src)
       do j = 2, nlat_src
          latb(j) = (lat(j-1) + lat(j))*0.5;
       end do
       latb(1) = 2*lat(1) - latb(2)
       latb(nlat_src+1) = 2*lat(nlat_src) - latb(nlat_src)
       ! lon_dst and lat_dst is at C-cell
       nlon_dst = nlon_dst - 1
       nlat_dst = nlat_dst - 1
    end if

    !--- if the src grid is not spherical grid, we need to get the geographical grid
    allocate(lon_src(nlon_src,nlat_src), lat_src(nlon_src,nlat_src) )
    if(src_is_spherical) then
       do i = 1, nlon_src
          lon_src(i,:) = lon(i)
       enddo
       do j = 1, nlat_src
          lat_src(:,j) = lat(j)
       enddo
    else !--- if the src grid is not spherical grid, we need to get the geographical grid
       found_lon = .false.; found_lat = .false.
       do n=1,nvar
          call mpp_get_atts(fields(n), name=name)
          if( trim(lowercase(name)) .eq. trim(lowercase(lon_field))) then
             call mpp_read(unit, fields(n), lon_src)
             found_lon = .true.
          else  if( trim(lowercase(name)) .eq. trim(lowercase(lat_field))) then            
             call mpp_read(unit, fields(n), lat_src)
             found_lat = .true.
          endif
       enddo
       if(.not.found_lon)call mpp_error(FATAL,"topog_mod: src_is_spherical is false, but field " //&
                           trim(lon_field)//" does not exist in the file "//trim(topog_file) )
       if(.not.found_lat)call mpp_error(FATAL,"topog_mod: src_is_spherical is false, but field " //&
                           trim(lat_field)//" does not exist in the file "//trim(topog_file) )
    endif

    !--- read the source data --------------------------------------------
    allocate(data_src(nlon_src,nlat_src), mask_src(nlon_src,nlat_src) )
    call mpp_read(unit, data_field, data_src)
    mask_src=1.0
    call mpp_get_atts(data_field,missing=missing)
    where(data_src == missing)          mask_src = 0.0
    where(data_src*scale_factor <= 0.0) mask_src = 0.0 

    call mpp_close(unit)
    deallocate(fields, axes)

    ! decompose model grid points
    ndivs = mpp_npes()
    call mpp_define_layout((/1,nlon_dst,1,nlat_dst/),ndivs,layout)
    call mpp_define_domains((/1,nlon_dst,1,nlat_dst/),layout,domain,xhalo=0,yhalo=0)  
    call mpp_get_compute_domain(domain,isc, iec, jsc, jec)

    ! --- check if a remap is needed or not
    if ( nlon_src == nlon_dst .and. nlat_src == nlat_dst .and. interp_method .NE. "conservative" ) then
       do_remap = .false.
       do j=1,nlat_dst
          do i=1,nlon_dst
             if (abs(lon_dst(i,j) -lon_src(i,j)) > grid_tol) then
                do_remap = .true.
                exit
             endif
          enddo
       enddo

       if (.not. do_remap) then
          do j=1,nlat_dst
             do i=1,nlon_dst
                if (abs(lat_dst(i,j)-lat_src(i,j)) > grid_tol) then
                   do_remap = .true.
                   exit
                endif
             enddo
          enddo
       endif
    endif

    if (do_remap) then
       ! use a temporary array for the regridding. use a global communication call to
       ! send the result on the global array back to all PEs

       allocate(tmp(isc:iec,jsc:jec),global_tmp(nlon_dst, nlat_dst))
       tmp = 0.0; data_dst=0.0
       select case(trim(interp_method))
       case('spherical')
             write (stdout(),'(//,10x,a/)') 'Interpolating topography spherical'
             call horiz_interp(data_src(:,:), lon_src*D2R, lat_src*D2R, lon_dst(isc:iec,jsc:jec)*D2R, &
                  lat_dst(isc:iec,jsc:jec)*D2R, tmp, mask_in=mask_src, interp_method = "spherical",   &
                  num_nbrs=num_nbrs, max_dist=max_dist)
       case('bilinear')       
          if(src_is_spherical) then
             write (stdout(),'(//,10x,a/)') 'Interpolating topography from spherical grid with method bilinear'
             call horiz_interp(data_src, lon*D2R, lat*D2R, lon_dst(isc:iec,jsc:jec)*D2R,   &
               lat_dst(isc:iec,jsc:jec)*D2R, tmp, mask_in=mask_src,                        &
               interp_method = "bilinear", grid_at_center = .true. )
          else
             write (stdout(),'(//,10x,a/)') 'Interpolating topography from non-spherical grid with method bilinear'
             call horiz_interp(data_src, lon_src*D2R, lat_src*D2R, lon_dst(isc:iec,jsc:jec)*D2R,   &
               lat_dst(isc:iec,jsc:jec)*D2R, tmp, mask_in=mask_src, interp_method = "bilinear")
          endif
       case('conservative')
          !--- find the starting and ending index of source grid to improve efficiency.
          min_lat = minval(lat_dst(isc:iec+1,jsc:jec+1))
          max_lat = maxval(lat_dst(isc:iec+1,jsc:jec+1))
          do j = 2, nlat_src+1
             jstart = j-1
             if(latb(j) > min_lat) exit
          enddo
          do j = nlat_src, 1, -1
             jend = j+1
             if(latb(j) < max_lat) exit
          enddo
          !--- To avoid truncation error, extend one more point.
          jstart = max(1, jstart-1)
          jend   = min(nlat_src+1, jend+1)

          if(src_is_spherical) then
             write (stdout(),'(//,10x,a/)') 'Interpolating topography from spherical grid with method conservative'
             call horiz_interp(data_src(:, jstart:jend-1), lonb*D2R, latb(jstart:jend)*D2R, lon_dst(isc:iec+1,jsc:jec+1)*D2R,   &
               lat_dst(isc:iec+1,jsc:jec+1)*D2R, tmp, mask_in=mask_src(:, jstart:jend-1), interp_method = "conservative" )
             deallocate(lonb, latb)
          else
             call mpp_error(FATAL,'topog_mod: for conservative interpolation, the source should be regular lon-lat grid.');
          end if
       case default
          call mpp_error(FATAL,'topog_mod: nml interp_method should be either "spherical", "bilinear" or "conservative" ')
       end select
       call mpp_global_field(domain, tmp,global_tmp)
       data_dst = global_tmp
       deallocate(tmp, global_tmp)
    else
       write (stdout(),'(//,10x,a/)') 'Topography is on-grid, no remapping'
       data_dst(:,:)=data_src(:,:)*mask_src(:,:)
    endif

    return

  end subroutine get_topog_from_file

  !#####################################################################
  subroutine set_topog_nml(is_full_cell, is_fill_isolated_cells, is_dont_change_landmask, &
       is_fill_shallow, is_deepen_shallow,is_round_shallow, is_adjust_topo,               &
       is_fill_first_row, min_kmt, verbose )    
    logical, intent(in)  :: is_full_cell, is_fill_isolated_cells, is_dont_change_landmask
    logical, intent(in)  :: is_fill_shallow, is_deepen_shallow, is_round_shallow
    logical, intent(in)  :: is_adjust_topo, is_fill_first_row, verbose
    integer, intent(in)  :: min_kmt

    full_cell            = is_full_cell
    fill_isolated_cells  = is_fill_isolated_cells
    dont_change_landmask = is_dont_change_landmask
    fill_shallow         = is_fill_shallow
    deepen_shallow       = is_deepen_shallow
    round_shallow        = is_round_shallow
    adjust_topo          = is_adjust_topo
    fill_first_row       = is_fill_first_row
    kmt_min              = min_kmt
    debug                = verbose
    small                = 1.e-6

    return

  end subroutine set_topog_nml

  !#####################################################################
  subroutine prep_obc_topog(kmt, ht, ni, nj)  
  ! Topography should be smooth neer open boundaries. Kmt and ht are set
  ! to the same value at the boundary point and two additional ocean points 
  ! near the boundary. To avoid steps, kmt and ht from the second inner point 
  ! near the boundary is used. The first oint outside the boundary is set to the
  ! depth at the boundary. In this case also velocity points at the B-grid are 
  ! well defined. If the boundary is at index i=1 or j=1, the buffer rows and 
  ! columnes at i=0 or j=0 are used. These points will not be stored in the 
  ! gridfile and the corresponding settings are evetually overwritten in other
  ! subroutines. So it is recommended, to use i>1 and j>1 as boundaries. 
    integer, dimension(0:,0:), intent(inout) :: kmt
    real,    dimension(0:,0:), intent(inout) :: ht
    integer,                   intent(in)    :: ni, nj
    integer :: i, j, n, ib, jb, ichange, kmt_temp
    real    :: ht_temp

    if (nobc == 0) return
      
    write(stdout(),'(/a/)') 'Checking bathymetry at open boundaries'
    do n=1, nobc
       write(stdout(),*) 'boundary', n, ' direction ',trim(direction(n))
       select case( trim(direction(n)) )
       case('west')
         ichange=0 
         if (is(n) .ne. ie(n)) call mpp_error(FATAL,'is must be equal to ie at a western boundary')
         ib = is(n)
!        If there is land near the boundary close the row.
!        If there is no land near the boundary set kmt and ht in
!        this row to ht(ib+nsmooth,jrow).
!        If this is land, the row is closed automatically.
         do j=js(n), je(n)
           kmt_temp = min(kmt(ib,j),kmt(ib+1,j),kmt(ib+2,j))
           if(kmt_temp.eq.0) then
             ht_temp = 0.
           else  
             ht_temp = ht (ib+nsmooth(n),j)
             kmt_temp= kmt(ib+nsmooth(n),j)
           endif 
           do i=is(n)-1,is(n)+2
! show all changes made to kmt and/or ht, use appropriate output format
             if (kmt(i,j).ne.kmt_temp .or. ht(i,j).ne.ht_temp) then
               ichange=ichange+1 
               write(stdout(),'(6x,2(a,i4),2(a,f10.2),a,i4)') &
                   'ht(',i,',',j,') =',ht_temp,' ! changed from',ht(i,j), ', kmt_new =',kmt_temp
               kmt(i,j) = kmt_temp
               ht(i,j)  = ht_temp
             endif
           enddo
         enddo
         if (ichange == 0) then
           write(stdout(),'(a)') 'No bathymetry changes at western obc'
         endif
       case('east')
         ichange=0 
         if (is(n) .ne. ie(n)) call mpp_error(FATAL,'is must be equal to ie at a eastern boundary')
         ib = is(n)
!        If there is land near the boundary close the row.
!        If there is no land near the boundary set kmt and ht in
!        this row to ht(ib-nsmooth,jrow).
!        If this is land, the row is closed automatically.
         do j=js(n), je(n)
           kmt_temp = min(kmt(ib,j),kmt(ib-1,j),kmt(ib-2,j))
           if(kmt_temp.eq.0) then
             ht_temp = 0.
           else  
             ht_temp = ht (ib-nsmooth(n),j)
             kmt_temp= kmt(ib-nsmooth(n),j)
           endif  
           do i=is(n)-2, min(is(n)+2,ni+1)
! show all changes made to kmt and/or ht, use appropriate output format
             if (kmt(i,j).ne.kmt_temp .or. ht(i,j).ne.ht_temp) then
               ichange=ichange+1 
               write(stdout(),'(6x,2(a,i4),2(a,f10.2),a,i4)')  &
                   'ht(',i,',',j,') =',ht_temp,' ! changed from',ht(i,j), ', kmt_new =',kmt_temp
               kmt(i,j) = kmt_temp
               ht(i,j)  = ht_temp
             endif
           enddo
         enddo
         if (ichange == 0) then
           write(stdout(),'(a)') 'No bathymetry changes at eastern obc'
         endif
       case('south')
         ichange=0 
         if (js(n) .ne. je(n)) call mpp_error(FATAL,'js must be equal to je at a southern boundary')
         jb = js(n)
!        If there is land near the boundary close the column.
!        If there is no land near the boundary set kmt and ht in
!        this row to ht(i,jb+2).
!        If this is land, the row is closed automatically.
         do i=is(n), ie(n)
           kmt_temp = min(kmt(i,jb),kmt(i,jb+1),kmt(i,jb+2))
           if(kmt_temp.eq.0) then
             ht_temp = 0.
           else  
             ht_temp = ht (i,jb+nsmooth(n))
             kmt_temp= kmt(i,jb+nsmooth(n))
           endif  
           do j=js(n)-1,js(n)+2
! show all changes made to kmt and/or ht, use appropriate output format
           if (kmt(i,j).ne.kmt_temp .or. ht(i,j).ne.ht_temp) then
               ichange=ichange+1 
               write(stdout(),'(6x,2(a,i4),2(a,f10.2),a,i4)')  &
                   'ht(',i,',',j,') =',ht_temp,' ! changed from',ht(i,j), ', kmt_new =',kmt_temp
               kmt(i,j) = kmt_temp
               ht(i,j)  = ht_temp
             endif
           enddo
         enddo
         if (ichange == 0) then
           write(stdout(),'(a)') 'No bathymetry changes at southern obc'
         endif
       case('north')
         ichange=0 
         if (js(n) .ne. je(n)) call mpp_error(FATAL,'js must be equal to je at a northern boundary')
         jb = js(n)
!        If there is land near the boundary close the column.
!        If there is no land near the boundary set kmt and ht in
!        this row to ht(i,jb-nsmooth).
!        If this is land, the row is closed automatically.
         do i=is(n), ie(n)
           kmt_temp = min(kmt(i,jb),kmt(i,jb-1),kmt(i,jb-2))
           if(kmt_temp.eq.0) then
             ht_temp = 0.
           else  
             ht_temp = ht (i,jb-nsmooth(n))
             kmt_temp= kmt(i,jb-nsmooth(n))
           endif  
           do j=js(n)-2,min(js(n)+2, nj+1)
! show all changes made to kmt and/or ht, use appropriate output format
             if (kmt(i,j).ne.kmt_temp .or. ht(i,j).ne.ht_temp) then
               ichange=ichange+1 
               write(stdout(),'(6x,2(a,i4),2(a,f10.2),a,i4)')  &
                   'ht(',i,',',j,') =',ht_temp,' ! changed from',ht(i,j), ', kmt_new =',kmt_temp
               kmt(i,j) = kmt_temp
               ht(i,j)  = ht_temp
             endif
           enddo
         enddo
         if (ichange == 0) then
           write(stdout(),'(a)') 'No bathymetry changes at northern obc'
         endif
       case default
          call mpp_error(FATAL,'each element of nml direction should be west, east, south or north')
       end select    
    enddo ! loop over all boundaries
    
    return
  end subroutine prep_obc_topog
  !#####################################################################

  end module topog_mod

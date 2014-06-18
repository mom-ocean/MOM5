module ocean_topog_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Set up ocean bottom topography. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Set up ocean bottom topography. Reads information from grid specification file. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_topog_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.
!  </DATA> 
!  <DATA NAME="flat_bottom" TYPE="logical">
!  For debugging, it is often useful to over-ride the grid spec file
!  and simply make the domain flat bottom. 
!  </DATA> 
!  <DATA NAME="flat_bottom_kmt" TYPE="integer">
!  Number of depth levels to use for the flat_bottom option.
!  </DATA> 
!  <DATA NAME="flat_bottom_ht" TYPE="real">
!  Depth to make the flat_bottom. 
!  </DATA>
!  <DATA NAME="min_thickness" TYPE="real">
!  min_thickness is only used for Mosaic grid. Since there is no kmt available
!  in mosaic grid, need to set min_thickness to configure kmt based on ht and zw.
!  Default min_thickness=1.0 metre.
!  </DATA> 
!  <DATA NAME="kmt_recompute" TYPE="logical">
!  To recompute the kmt array based on min_thickness.  This step is not recommended
!  in general, since it can modify the kmt array which may be in the grid spec file.  
!  But it may be of use for specialized situations, such as when you wish to use 
!  the same topography file with a refined vertical resolution.  
!  </DATA> 
!  <DATA NAME="kmt_recompute_offset" TYPE="integer">
!  To recompute the kmt array based on min_thickness, with an offset
!  determined by kmt_recompute_offset.  Default kmt_recompute_offset=0.
!  </DATA> 
! </NAMELIST>
!
use fms_mod,         only: open_namelist_file, close_file, check_nml_error
use fms_mod,         only: field_exist, read_data, file_exist, write_data
use mpp_domains_mod, only: mpp_update_domains
use mpp_mod,         only: input_nml_file, mpp_error, mpp_min, mpp_pe
use mpp_mod,         only: FATAL, WARNING, NOTE, stdout, stdlog
use axis_utils_mod,  only: nearest_index

use ocean_domains_mod,    only: get_local_indices, get_domain_offsets
use ocean_parameters_mod, only: TERRAIN_FOLLOWING, grav
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_util_mod,       only: write_chksum_2d, write_chksum_2d_int

implicit none
private

character(len=256) :: version='CVS $Id: ocean_topog.F90,v 20.0 2013/12/14 00:12:35 fms Exp $'
character(len=256) :: tagname='Tag $Name: tikal $'

#include <ocean_memory.h>

! for output
integer :: writeunit=6

logical :: flat_bottom          = .false.
integer :: flat_bottom_kmt      = 50
real    :: flat_bottom_ht       = 5500.0
real    :: min_thickness        = 1.0
integer :: kmt_recompute_offset = 0
logical :: kmt_recompute        = .false.
logical :: write_topog          = .false.
logical :: debug_this_module    = .true.

namelist /ocean_topog_nml/ flat_bottom, flat_bottom_kmt, flat_bottom_ht, write_topog, &
                           min_thickness, kmt_recompute, kmt_recompute_offset,        &
                           debug_this_module

public ocean_topog_init

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_topog_init">
!
! <DESCRIPTION>
! Initialize the ocean bottom topography.
!
! There are two reasons to prefer land be at j=1. 
!
! A/ MOM employs a northest B-grid.  To construct 
! horizontal remapping operators, we need information 
! about grid factors one row outside of the global 
! domain boundaries.  In particular, we need j=0
! grid information. However, when constructing the grid spec
! file, we assume nothing about the region outside 
! the global domain.  So MOM's requirement of j=0
! grid information necessitates extrapolation.  
! This extrapolation is done inside ocean_grids.F90, and
! it can lead to non-symmetric values of grid spacing
! for the region j=0 and j=nj, even if the domain global 
! limits are symmetric across the equator.  
!
! B/ The FMS Sea Ice Simulator (SIS) requires land to be present
! for all points at jsc+joff=1.  If this is not the case, then 
! the ocean model cannot be coupled to SIS.  The SIS requirement 
! of all land at jsc+joff=1 is related to the use of a northeast 
! B-grid convention.  To couple the models, the ocean grid
! must be generated with fill_first_row=.true. 
!
! </DESCRIPTION>
!
subroutine ocean_topog_init (Domain, Grid, grid_file, vert_coordinate_type)

  type(ocean_domain_type), intent(inout)        :: Domain
  type(ocean_grid_type),   intent(inout)        :: Grid
  character(len=*),        intent(in), optional :: grid_file
  integer,                 intent(in)           :: vert_coordinate_type

  character(len=128)  :: ocean_topog = "INPUT/topog.nc"
  character(len=128)  :: grd_file
  logical             :: land_jeq1
  logical             :: error_flag=.false.
  real                :: min_depth, min_depth0, fudge 
  integer             :: imin, jmin, kmin
  integer             :: ioun, io_status, ierr
  integer             :: i, j, ioff, joff
  integer             :: k, kb                        

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  grd_file = "INPUT/grid_spec.nc"
  if(present(grid_file)) grd_file = grid_file

  write( stdlogunit,'(/a/)') trim(version)

  joff=0
  ioff=0
#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (Grid%kmt(isd:ied,jsd:jed))
  allocate (Grid%ht(isd:ied,jsd:jed))
  allocate (Grid%htr(isd:ied,jsd:jed))
  allocate (Grid%kmu(isd:ied,jsd:jed))
  allocate (Grid%hu(isd:ied,jsd:jed))
#else
  call get_domain_offsets(Domain,ioff,joff)
#endif
  

  ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_topog_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_topog_nml')
#else
  ioun = open_namelist_file()
  read (ioun,ocean_topog_nml,IOSTAT=io_status)
  ierr = check_nml_error(io_status, 'ocean_topog_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_topog_nml)  
  write (stdlogunit,ocean_topog_nml)

  Grid%kmt=0;Grid%kmu=0;Grid%ht=0;Grid%hu=0
  if(field_exist(grd_file, 'depth_t') ) then ! new grid file
     call read_data(grd_file, "depth_t", Grid%ht(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(grd_file, "num_levels", Grid%kmt(isc:iec,jsc:jec), Domain%domain2d)
     if(field_exist(grd_file, 'depth_c')) then
        if(.not. field_exist(grd_file, 'num_levels_c')) call mpp_error(FATAL, &
           'ocean_topog_mod: depth_c exist but num_levels_c does not exist in file '//trim(grd_file) )
        call read_data(grd_file, "depth_c", Grid%hu(isc:iec,jsc:jec), Domain%domain2d)
        call read_data(grd_file, "num_levels_c", Grid%kmu(isc:iec,jsc:jec), Domain%domain2d)
     elseif(field_exist(grd_file, 'num_levels_c')) then
        call mpp_error(FATAL, &
           'ocean_topog_mod: num_levels_c exist but depth_c does not exist in file '//trim(grd_file) )
     endif

  elseif(field_exist(grd_file, 'ht') ) then ! old grid file
     call read_data(grd_file, "ht", Grid%ht(isc:iec,jsc:jec), Domain%domain2d)
     call read_data(grd_file, "kmt", Grid%kmt(isc:iec,jsc:jec), Domain%domain2d)

  elseif(file_exist(ocean_topog)) then
     call read_data(ocean_topog, 'depth', Grid%ht(isc:iec,jsc:jec), Domain%domain2d)

     !--- calculate kmt based on ht and zw.
     do j=jsc,jec
        do i=isc,iec
           if(Grid%ht(i,j) <= 0.0) then
              Grid%kmt(i,j) = 0
           else
              Grid%kmt(i,j) = nearest_index(Grid%ht(i,j), Grid%zw)
              if( Grid%ht(i,j) >= Grid%zw(Grid%kmt(i,j)) + min_thickness) then 
                 if(Grid%kmt(i,j) < nk) then
                    Grid%kmt(i,j) = Grid%kmt(i,j) + 1
                 endif
              endif
           endif
        enddo
     enddo
  else
     call mpp_error(FATAL, 'ocean_topog_mod: depth_t and ht do not exist in file '//trim(grd_file)// &
                    ', also file '//trim(ocean_topog)//' does not exist')
  endif

  if(kmt_recompute) then 
      call mpp_error(NOTE, 'ocean_topog_mod: recomputing kmt array given ht, zw, and min_thickness.')
      do j=jsc,jec
         do i=isc,iec
            if(Grid%kmt(i,j) > 1) then
                kb = Grid%kmt(i,j)              
                kbloop:          do k=2,kb
                   if( Grid%zw(k) + min_thickness > Grid%ht(i,j) ) then 
                       Grid%kmt(i,j) = k-kmt_recompute_offset
                       if(debug_this_module) then 
                          write(*,'(a,3i6,f12.3,a,3i6,f12.3)') 'kmt_recompute: resetting i,j,kmt_orig,ht= ',&
                          i,j,kb,Grid%ht(i,j),'  to ',  i,j,Grid%kmt(i,j),Grid%ht(i,j)
                       endif 
                       exit kbloop
                   endif
                enddo kbloop
            endif
         enddo
      enddo
  endif

  call mpp_update_domains(Grid%ht, Domain%domain2d)
  call mpp_update_domains(Grid%kmt,Domain%domain2d)

  if(write_topog) then
     call write_data("ocean_grid_check.nc", "ht", Grid%ht(isc:iec,jsc:jec), Domain%domain2d)
     call write_data("ocean_grid_check.nc", "kmt", Grid%kmt(isc:iec,jsc:jec), Domain%domain2d)
  end if

  ! ensure there are enough kmt levels for nontrivial vertical processes to be computed 
  error_flag=.false. 
  do j=jsc,jec
     do i=isc,iec
        if(Grid%kmt(i,j)==1) then 
            write(writeunit,'(a,i4,a,i4,a)') &
            '==>ocean_topog_init: WARNING: kmt(',i+Domain%ioff,',',j+Domain%joff, &
            ') = 1 is not permitted. kmt>2 recommended for wet ocean domain.'
            error_flag=.true.
        endif
        if(Grid%kmt(i,j)==2) then 
            write(writeunit,'(a,i4,a,i4,a)') &
            '==>ocean_topog_init: WARNING: kmt(',i+Domain%ioff,',',j+Domain%joff, &
            ') = 2 is generally not permitted. kmt>2 recommended for wet ocean domain.'
            error_flag=.true.
        endif
     enddo
  enddo
  if(error_flag) then 
      call mpp_error(WARNING,&
     '==>ocean_topog_init: Regions with 0 < kmt < 2 have been found. Recommend kmt > 2 for wet ocean regions.')
  endif


  ! for terrain following vertical coordinates, all levels are squeazed into wet regions 
  if(vert_coordinate_type == TERRAIN_FOLLOWING) then
      call mpp_error(NOTE,&
      '==>ocean_topog_init: terrain following coordinate places all nk levels into wet regions.')
      do j=jsd,jed
         do i=isd,ied
            if(Grid%kmt(i,j) > 0) then  
                Grid%kmt(i,j) = nk
            endif
         enddo
      enddo
  endif

  ! ensure that land points have ht = 0.0.
  ! this step may be needed if read ht from land/sea dataset, 
  ! or if ht contains "missing values" from a visualization package. 
  do j=jsd,jed
     do i=isd,ied
        if(Grid%ht(i,j) <= 0.0) Grid%ht(i,j) = 0.0
     enddo
  enddo

  if(flat_bottom) then 
      call mpp_error(NOTE,&
      '==>ocean_topog_init: flat_bottom=.true. so will set all cells that are not land equal to their max depth.')
      do j=jsd,jed
         do i=isd,ied
            if(Grid%kmt(i,j) > 0) then  
               Grid%kmt(i,j) = flat_bottom_kmt
               Grid%ht(i,j)  = flat_bottom_ht
            endif  
         enddo
      enddo
  endif


  ! if num_levels_c exist in the grid file, depths at u-cell points will be get from grid file,
  ! otherwise construct depths at u-cell points as minimum of surrounding t cells
  if(.NOT. field_exist(grd_file, 'num_levels_c')) then
     do j=jsc,jec
        do i=isc,iec
           Grid%kmu(i,j) = min(Grid%kmt(i,j), Grid%kmt(i+1,j), Grid%kmt(i,j+1), Grid%kmt(i+1,j+1))
           Grid%hu(i,j)  = min(Grid%ht(i,j), Grid%ht(i+1,j), Grid%ht(i,j+1), Grid%ht(i+1,j+1))
        enddo
     enddo
  end if
  call mpp_update_domains(Grid%kmu,Domain%domain2d)
  call mpp_update_domains(Grid%hu,Domain%domain2d)

  ! check whether any ocean occupies j=1. 
  if (jsc+joff==1) then
    land_jeq1 = .false. 
    do i=isc,iec
      if(Grid%kmt(i,1) > 0) land_jeq1 = .true. 
    enddo
    if(land_jeq1) then 
      call mpp_error(NOTE,&
      '==>ocean_topog_init: Ocean grid has water at j=1.  Recommend placing land at this latitude row.')
    endif 
  endif 

  ! set inverse depth on t-cells 
  Grid%htr = 0.0 
  do j=jsd,jed
     do i=isd,ied
        if(Grid%kmt(i,j) > 0) then    
            Grid%htr(i,j) = 1.0/Grid%ht(i,j)
        endif
     enddo
  enddo

  ! find shallowest water and provide caveat for overly shallow regions
  imin= 0 ; jmin=0 ; kmin=0 ; min_depth = 1.e10
  do j=jsc,jec
    do i=isc,iec
      if(Grid%kmt(i,j) > 0 .and. Grid%ht(i,j) < min_depth) then 
        imin      = i
        jmin      = j
        kmin      = Grid%kmt(i,j)
        min_depth = Grid%ht(i,j)
      endif 
    enddo
  enddo 
  fudge      = 1 + 1.e-12*mpp_pe() ! to distinguish processors when min_depth is independent of processor
  min_depth  = min_depth*fudge 
  min_depth0 = min_depth
  call mpp_min(min_depth)  

  if(min_depth0 == min_depth) then 
      write(writeunit,'(/a,f12.5)')' The shallowest wet ocean model grid cell has depth (meters) ', min_depth  
      write(writeunit,'(a,i3,a,i3,a,i3,a)')' and this occurs at (i,j,k) = (',&
                      imin+Domain%ioff,',',jmin+Domain%joff ,',',kmin,')'  
      write(writeunit,'(a,f10.4,a,f10.4,a,f12.5,a/)')' which has (long,lat,depth) =  (' &
           ,Grid%xt(imin,jmin),',',Grid%yt(imin,jmin), ',', Grid%ht(imin,jmin),')'
      if(min_depth < 50.0) then 
          write(writeunit,'(a)')' Beware that shallow regions (e.g., those shallower than 50m) may be subject' 
          write(writeunit,'(a)')' to numerical problems if strong surface forcing is not mixed vertically.' 
          write(writeunit,'(a)')' Such problems may occur especially in shallow regions with kmt==2.'
          write(writeunit,'(a)')' Current speeds and/or tracer deviations may become large due to the deposition'
          write(writeunit,'(a)')' of wind and/or buoyancy over just a small upper ocean region. Such problems'
          write(writeunit,'(a)')' can be resolved by adding sufficient vertical mixing in these regions.'
          write(writeunit,'(a)')' Such happens in Nature due to tides and breaking surface waves.'
      endif
  endif

  write(stdoutunit,*) 'Topography checksums' 
  call write_chksum_2d('ht', Grid%ht(COMP))
  call write_chksum_2d('hu', Grid%hu(COMP))
  call write_chksum_2d('htr', Grid%htr(COMP))
  call write_chksum_2d_int('kmu', Grid%kmu(COMP))
  call write_chksum_2d_int('kmt', Grid%kmt(COMP))
  write(stdoutunit,*) ' ' 
  

end subroutine ocean_topog_init
! </SUBROUTINE> NAME="ocean_topog_init"


end module ocean_topog_mod


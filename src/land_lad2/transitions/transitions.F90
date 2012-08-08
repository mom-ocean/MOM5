#include <fms_platform.h>

module land_transitions_mod

#include "../shared/debug.inc"

use constants_mod, only : PI

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod, only : write_version_number, string, error_mesg, FATAL, WARNING, NOTE, &
     mpp_pe, write_version_number, file_exist, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII
use time_manager_mod, only : time_type, set_date, get_date, set_time, &
     operator(+), operator(-), operator(>), operator(<), operator(<=), operator(/), &
     operator(//), days_in_year, print_date, increment_date, get_time, &
     valid_calendar_types, get_calendar_type
use get_cal_time_mod, only : get_cal_time
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del, horiz_interp
use time_interp_mod, only : time_interp
use diag_manager_mod, only : register_diag_field, send_data

use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_valid_range, nfu_is_valid

use vegn_data_mod, only : &
     N_LU_TYPES, LU_NTRL, LU_SCND, landuse_name, landuse_longname

use cana_tile_mod, only : cana_tile_heat
use snow_tile_mod, only : snow_tile_heat
use vegn_tile_mod, only : vegn_tile_heat
use soil_tile_mod, only : soil_tile_heat

use land_tile_mod, only : &
     land_tile_type, land_tile_list_type, land_tile_enum_type, new_land_tile, delete_land_tile, &
     first_elmt, tail_elmt, next_elmt, operator(/=), operator(==), current_tile, &
     land_tile_list_init, land_tile_list_end, &
     empty, erase, remove, insert, land_tiles_can_be_merged, merge_land_tiles, &
     get_tile_water, land_tile_carbon, land_tile_heat
use land_tile_io_mod, only : print_netcdf_error

use land_data_mod, only : &
     land_data_type, lnd
use vegn_tile_mod, only : &
     vegn_tile_type, vegn_tran_priority
use vegn_harvesting_mod, only : &
     vegn_cut_forest

use land_debug_mod, only : set_current_point, is_watch_point, get_current_point
     
implicit none
private

! ==== public interface =====================================================
public :: land_transitions_init
public :: land_transitions_end
public :: save_land_transitions_restart

public :: land_transitions
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter   :: &
     version = '$Id: transitions.F90,v 19.0 2012/01/06 20:43:22 fms Exp $', &
     tagname = '$Name: siena_201207 $', &
     module_name = 'land_transitions_mod', &
     diag_mod_name = 'landuse'
! selectors for overshoot handling options, for efficiency
integer, parameter :: &
     OPT_IGNORE = 0, &
     OPT_STOP   = 1, &
     OPT_REPORT = 2

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== data types ===========================================================
type :: tran_type
   integer :: donor    = 0  ! kind of donor tile
   integer :: acceptor = 0  ! kind of acceptor tile
   real    :: frac     = 0  ! area of transition
end type tran_type

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid ! netcd id of the input file
integer :: input_ids (N_LU_TYPES,N_LU_TYPES) ! id's of input transition rate fields
integer :: diag_ids  (N_LU_TYPES,N_LU_TYPES)
real, allocatable :: buffer_in(:,:) ! buffer for input data
type(horiz_interp_type), save :: interp
type(time_type), allocatable :: time_in(:) ! time axis in input data
type(time_type) :: time0 ! time of previous transition calculations
integer :: overshoot_opt ! selector for overshoot handling options, for efficiency
integer :: conservation_opt ! selector for non-conservation handling options, for efficiency

! ---- namelist variables ---------------------------------------------------
logical :: do_landuse_change = .FALSE. ! if true, then the landuse changes with time
character(len=512) :: input_file = ''
! sets how to handle transition overshoot: that is, the situation when transition 
! is larger than available area of the given land use type.
character(len=16) :: overshoot_handling = 'report' ! or 'stop', or 'ignore'
real :: overshoot_tolerance = 1e-4 ! tolerance interval for overshoots 
! specifies how to handle non-conservation
character(len=16) :: conservation_handling = 'stop' ! or 'report', or 'ignore'

namelist/landuse_nml/input_file, do_landuse_change, &
     overshoot_handling, overshoot_tolerance, &
     conservation_handling
     

contains ! ###################################################################

! ============================================================================
subroutine land_transitions_init(id_lon, id_lat)
  integer, intent(in) :: id_lon, id_lat ! the IDs of land diagnostic axes

  ! ---- local vars
  logical        :: grid_initialized = .false.
  integer        :: len, unit, ierr, io
  integer        :: year,month,day,hour,min,sec
  integer        :: k1,k2,i
  real,allocatable :: lon_in(:,:),lat_in(:,:)
  character(len=12) :: fieldname
  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  integer :: nrec ! number of records in the input file
  real, allocatable :: time(:)  ! real values of time coordinate
  real, allocatable :: mask_in(:,:) ! valid data mask on the input data grid
  type(nfu_validtype) :: v ! valid values range
  integer :: timedim ! id of the record (time) dimension
  integer :: timevar ! id of the time variable
  character(len=NF_MAX_NAME) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units ot time in the file
  character(len=24) :: calendar ! model calendar

  if(module_is_initialized) return

  call horiz_interp_init
  call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=landuse_nml, iostat=io)
  ierr = check_nml_error(io, 'landuse_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=landuse_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'landuse_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=landuse_nml)
  endif

  ! read restart file, if any
  if (file_exist('INPUT/landuse.res')) then
     call error_mesg('land_transitions_init',&
          'reading restart "INPUT/landuse.res"',&
          NOTE)
     call mpp_open(unit,'INPUT/landuse.res', action=MPP_RDONLY, form=MPP_ASCII)
     read(unit,*) year,month,day,hour,min,sec
     time0 = set_date(year,month,day,hour,min,sec)
     call mpp_close(unit)
  else
     call error_mesg('land_transitions_init',&
          'cold-starting land transitions',&
          NOTE)
     time0 = set_date(0001,01,01);
  endif
  
  ! parse the overshoot handling option
  if (trim(overshoot_handling)=='stop') then
     overshoot_opt = OPT_STOP
  else if (trim(overshoot_handling)=='ignore') then
     overshoot_opt = OPT_IGNORE
  else if (trim(overshoot_handling)=='report') then
     overshoot_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','overshoot_handling value "'//&
          trim(overshoot_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  ! parse the non-conservation handling option
  if (trim(conservation_handling)=='stop') then
     conservation_opt = OPT_STOP
  else if (trim(conservation_handling)=='ignore') then
     conservation_opt = OPT_IGNORE
  else if (trim(conservation_handling)=='report') then
     conservation_opt = OPT_REPORT
  else
     call error_mesg('land_transitions_init','conservation_handling value "'//&
          trim(conservation_handling)//'" is illegal, use "stop", "report", or "ignore"',&
          FATAL)
  endif

  if (do_landuse_change) then
     if (trim(input_file)=='') call error_mesg('landuse_init', &
          'do_landuse_change is requested, but landuse transition file is not specified', &
          FATAL)

     ierr=nf_open(input_file,NF_NOWRITE,ncid)

     if(ierr/=NF_NOERR) call error_mesg('landuse_init', &
          'do_landuse_change is requested, but landuse transition file "'// &
          trim(input_file)//'" could not be opened because '//nf_strerror(ierr), FATAL)
  
     ! initialize array of input field ids
     input_ids(:,:) = 0
     do k1 = 1,size(input_ids,1)
     do k2 = 1,size(input_ids,2)
        ! construct a name of input field and register the field
        fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
        if(trim(fieldname)=='2') cycle ! skip unspecified tiles
        
        ierr = nfu_inq_var(ncid,fieldname, id=input_ids(k1,k2))
        if (ierr/=NF_NOERR) then
           if (ierr==NF_ENOTVAR) then
              input_ids(k1,k2)=-1 ! to indicate that input field is not present in the input data
           else
              call error_mesg('landuse_init',&
                   'error initializing field "'//trim(fieldname)//&
                   '" from file "'//trim(input_file)//'" : '//&
                   nf_strerror(ierr), &
                   FATAL)
           endif
        endif
        ! initialize the input data grid and horizontal interpolator
        if ((.not.grid_initialized).and.(input_ids(k1,k2)>0)) then
           ! we assume that all transition rate fields are specified on the same grid, 
           ! in both horizontal and time "directions". Therefore there is a single grid
           ! for all fields, initialized only once.
           
           __NF_ASRT__(nfu_inq_var(ncid,fieldname,dimids=dimids,dimlens=dimlens,nrec=nrec))
           ! allocate temporary variables
           allocate(lon_in(dimlens(1)+1,1), &
                    lat_in(1,dimlens(2)+1), &
                    time(nrec), mask_in(dimlens(1),dimlens(2)) )
           ! allocate module data
           allocate(buffer_in(dimlens(1),dimlens(2)),time_in(nrec))

           ! get the boundaries of the horizontal axes and initialize horizontal
           ! interpolator
           __NF_ASRT__(nfu_get_dim_bounds(ncid, dimids(1), lon_in(:,1)))
           __NF_ASRT__(nfu_get_dim_bounds(ncid, dimids(2), lat_in(1,:)))

           ! get the first record from variable and obtain the mask of valid data
           ! assume that valid mask does not change with time
           __NF_ASRT__(nfu_get_rec(ncid,fieldname,1,buffer_in))
           ! get the valid range for the variable
           __NF_ASRT__(nfu_get_valid_range(ncid,fieldname,v))
           ! get the mask
           where (nfu_is_valid(buffer_in,v))
              mask_in = 1
           elsewhere
              mask_in = 0
           end where

           ! add mask_in and mask_out to this call
           call horiz_interp_new(interp, lon_in*PI/180,lat_in*PI/180, &
                lnd%lonb, lnd%latb, &
                interp_method='conservative',&
                mask_in=mask_in, is_latlon_in=.TRUE. )
           
           ! get the time axis 
           __NF_ASRT__(nf_inq_unlimdim(ncid, timedim))
           __NF_ASRT__(nfu_get_dim(ncid, timedim, time))
           ! get units of time
           __NF_ASRT__(nf_inq_dimname(ncid, timedim, timename))
           __NF_ASRT__(nf_inq_varid(ncid,timename, timevar))
           timeunits = ' '
           __NF_ASRT__(nf_get_att_text(ncid,timevar,'units',timeunits))
           ! get model calendar
           calendar=valid_calendar_types(get_calendar_type())

           ! loop through the time axis and get time_type values in time_in
           do i = 1,size(time)
              time_in(i) = get_cal_time(time(i),timeunits,calendar)
           end do
           
           grid_initialized = .true.
           ! get rid of allocated data
           deallocate(lon_in,lat_in,time,mask_in)
        endif
     enddo
     enddo
  endif ! do_landuse_changes

  ! initialize diagnostics
  diag_ids(:,:) = 0

  do k1 = 1,size(diag_ids,1)
  do k2 = 1,size(diag_ids,2)
     ! skip unnamed tiles
     if(landuse_name(k1)=='')cycle
     if(landuse_name(k2)=='')cycle
     ! construct a name of input field and register the field
     fieldname = trim(landuse_name(k1))//'2'//trim(landuse_name(k2))
     diag_ids(k1,k2) = register_diag_field(diag_mod_name,fieldname,(/id_lon,id_lat/), lnd%time, &
          'rate of transition from '//trim(landuse_longname(k1))//' to '//trim(landuse_longname(k2)),& 
          units='1/year', missing_value=-1.0)
  enddo
  enddo

  module_is_initialized=.TRUE.

end subroutine


! ============================================================================
subroutine land_transitions_end()
  
  if (do_landuse_change) &
       call horiz_interp_del(interp)
  if(allocated(time_in)) &
       deallocate(time_in,buffer_in)
  module_is_initialized=.FALSE.

end subroutine


! ============================================================================
subroutine save_land_transitions_restart(timestamp)
  character(*), intent(in) :: timestamp ! timestamp to add to the file name
  
  integer :: unit,year,month,day,hour,min,sec

  call mpp_open( unit, 'RESTART/'//trim(timestamp)//'landuse.res', nohdrs=.TRUE. )
  if (mpp_pe() == mpp_root_pe()) then
     call get_date(time0, year,month,day,hour,min,sec)
     write(unit,'(6i6,8x,a)') year,month,day,hour,min,sec, &
          'Time of previous landuse transition calculation'
  endif
  call mpp_close(unit)

end subroutine save_land_transitions_restart


! ============================================================================
! performs transitions between tiles, e.g. conversion of forests to crop, etc.
subroutine land_transitions (time)
  type(time_type), intent(in) :: time 

  ! ---- local vars.
  integer :: i,j,k1,k2
  type(tran_type), pointer :: transitions(:,:,:) => NULL()
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1

  if (.not.do_landuse_change) &
       return ! do nothing if landuse change not requested
  ! NB: in this case file/interp/data are not initialized, so it is
  ! not even possible to use the code below

  call get_date(time,             year0,month0,day0,hour,minute,second)
  call get_date(time-lnd%dt_slow, year1,month1,day1,hour,minute,second)
  if(year0 == year1) &
!!$  if(day0 == day1) &
       return ! do nothing during a year 

  ! get transition rates for current time: read map of transitions, and accumulate
  ! as many layers in array of transitions as necessary. Note that "transitions"
  ! array gets reallocated inside get_transitions as necessary, it has only as many
  ! layers as the max number of transitions occuring at a point at the time.
  do k1 = 1,N_LU_TYPES
  do k2 = 1,N_LU_TYPES
     call get_transitions(time0,time,k1,k2,transitions)
  enddo
  enddo

  ! perform the transitions
  do j = lnd%js,lnd%je
  do i = lnd%is,lnd%ie
     if(empty(lnd%tile_map(i,j))) cycle ! skip cells where there is no land
     ! set current point for debugging
     call set_current_point(i,j,1)
     ! transiton land area between different tile types
     call land_transitions_0d(lnd%tile_map(i,j), &
          transitions(i,j,:)%donor, &
          transitions(i,j,:)%acceptor,&
          transitions(i,j,:)%frac )
  enddo
  enddo
  
  ! deallocate array of transitions
  if (associated(transitions)) deallocate(transitions)
  
  ! store current time for future reference
  time0=time

end subroutine land_transitions


! =============================================================================
! performs tile transitions in a given grid cell
subroutine land_transitions_0d(d_list,d_kinds,a_kinds,area)
  type(land_tile_list_type), intent(inout) :: d_list ! list of tiles
  integer, intent(in) :: d_kinds(:) ! array of donor tile kinds
  integer, intent(in) :: a_kinds(:) ! array of acceptor tile kinds
  real   , intent(in) :: area(:)    ! array of areas changing from donor tiles to acceptor tiles
  
  ! ---- local vars
  integer :: i
  type(land_tile_type), pointer :: ptr
  type(land_tile_list_type) :: a_list
  type(land_tile_enum_type) :: ts, te
  real :: atot ! total fraction of tiles that can be involved in transitions
  ! variable used for conservation check:
  real :: lmass0, fmass0, cmass0, heat0, &
       soil_heat0, vegn_heat0, cana_heat0, snow_heat0 ! pre-transition values 
  real :: lmass1, fmass1, cmass1, heat1, &
       soil_heat1, vegn_heat1, cana_heat1, snow_heat1 ! post-transition values
  real :: lm, fm ! buffers for transition calulations

  ! conservation check code, part 1: calculate the pre-transition grid
  ! cell totals
  lmass0 = 0 ; fmass0 = 0 ; cmass0 = 0 ; heat0 = 0
  soil_heat0 = 0 ;  vegn_heat0 = 0 ; cana_heat0 = 0 ; snow_heat0 = 0
  ts = first_elmt(d_list) ; te=tail_elmt(d_list)
  do while (ts /= te)
     ptr=>current_tile(ts); ts=next_elmt(ts)
     call get_tile_water(ptr,lm,fm)
     lmass0 = lmass0 + lm*ptr%frac ; fmass0 = fmass0 + fm*ptr%frac

     heat0  = heat0  + land_tile_heat  (ptr)*ptr%frac
     cmass0 = cmass0 + land_tile_carbon(ptr)*ptr%frac
     
     if(associated(ptr%soil)) soil_heat0 = soil_heat0 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat0 = vegn_heat0 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat0 = cana_heat0 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat0 = snow_heat0 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo

  ! calculate the area that can participate in land transitions
  atot = 0 ; ts = first_elmt(d_list) ; te=tail_elmt(d_list)
  do while (ts /= te)
     ptr=>current_tile(ts); ts=next_elmt(ts)
     if (.not.associated(ptr%vegn)) cycle
     atot = atot + ptr%frac
  enddo

  if (is_watch_point()) then
     write(*,*)'### land_transitions_0d: input parameters ###'
     do i = 1, size(d_kinds)
        __DEBUG4__(i,d_kinds(i),a_kinds(i),area(i))
        ! write(*,'(a,i2.2,100(2x,a,g))')'i='i,d_kinds(i)d_kinds(i),a_kinds(i),)
     enddo

     write(*,*)'### land_transitions_0d: land fractions before transitions (initial state) ###'
     ts = first_elmt(d_list) ; te=tail_elmt(d_list)
     do while (ts /= te)
        ptr=>current_tile(ts); ts=next_elmt(ts)
        if (.not.associated(ptr%vegn)) cycle
        write(*,*)'landuse=',ptr%vegn%landuse,' area=',ptr%frac
     enddo
     write(*,'(a,g)')'total area=',atot
  endif

  ! split each donor tile and gather the parts that undergo a 
  ! transition into a separate list. Note that the kind of the landuse is
  ! changed during this transition, including forest harvesting if necessary.
  ! This has to occur at some time before the tiles are merged, and it seems
  ! to be the most convenient place as both original and final landuse kind
  ! is known for each part.
  call land_tile_list_init(a_list)
  do i = 1,size(d_kinds)
     call split_changing_tile_parts(d_list,d_kinds(i),a_kinds(i),area(i)*atot,a_list)
     ! the factor atot normalizes the transitions to the total area in the grid cell
     ! available for the land use, that is, the area of land excluding lakes and glaciers
  enddo
  if (is_watch_point()) then
     write(*,*)'### land_transitions_0d: land fractions after splitting changing parts ###'
     atot = 0 ; ts = first_elmt(d_list) ; te=tail_elmt(d_list)
     do while (ts /= te)
        ptr=>current_tile(ts); ts=next_elmt(ts)
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g,2x))')'   donor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     ts = first_elmt(a_list); te=tail_elmt(a_list)
     do while (ts /= te)
        ptr=>current_tile(ts); ts=next_elmt(ts)
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g,2x))')'acceptor: landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     write(*,'(a,g)')'total area=',atot
  endif

  ! move all tiles from the donor list to the acceptor list -- this will ensure
  ! that all the tiles that can be merged at this time will be
  te = tail_elmt(d_list) 
  do 
     ts=first_elmt(d_list)
     if(ts==te) exit ! reached the end of the list
     ptr=>current_tile(ts)
     if(ptr%frac <= 0.0) then
        call erase(ts) ! if area of the tile is zero, free it
     else
        ! othervise, move it to a_list
        call remove(ts)
        call insert(ptr,a_list)
     endif
  enddo
  ! d_list is empty at this point

  ! merge all generated tiles into the source (donor) list
  te = tail_elmt(a_list) 
  do
     ts=first_elmt(a_list)
     if(ts==te) exit ! break out of loop
     ptr=>current_tile(ts)
     call remove(ts)
     call land_tile_merge(ptr,d_list)
  enddo
  ! a_list is empty at this point
  call land_tile_list_end(a_list)

  if (is_watch_point()) then
     write(*,*)'### land_transitions_0d: land fractions final state ###'
     atot = 0
     ts = first_elmt(d_list); te=tail_elmt(d_list)
     do while (ts /= te)
        ptr=>current_tile(ts); ts=next_elmt(ts)
        if (.not.associated(ptr%vegn)) cycle
        write(*,'(2(a,g,2x))')'landuse=',ptr%vegn%landuse,' area=',ptr%frac
        atot = atot + ptr%frac
     enddo
     write(*,'(a,g)')'total area=',atot
  endif

  ! conservation check part 2: calculate grid cell totals in final state, and 
  ! compare them with pre-transition totals
  lmass1 = 0 ; fmass1 = 0 ; cmass1 = 0 ; heat1 = 0
  soil_heat1 = 0 ;  vegn_heat1 = 0 ; cana_heat1 = 0 ; snow_heat1 = 0
  ts = first_elmt(d_list) ; te=tail_elmt(d_list)
  do while (ts /= te)
     ptr=>current_tile(ts); ts=next_elmt(ts)
     call get_tile_water(ptr,lm,fm)
     lmass1 = lmass1 + lm*ptr%frac ; fmass1 = fmass1 + fm*ptr%frac

     heat1  = heat1  + land_tile_heat  (ptr)*ptr%frac
     cmass1 = cmass1 + land_tile_carbon(ptr)*ptr%frac

     if(associated(ptr%soil)) soil_heat1 = soil_heat1 + soil_tile_heat(ptr%soil)*ptr%frac
     if(associated(ptr%vegn)) vegn_heat1 = vegn_heat1 + vegn_tile_heat(ptr%vegn)*ptr%frac
     if(associated(ptr%cana)) cana_heat1 = cana_heat1 + cana_tile_heat(ptr%cana)*ptr%frac
     if(associated(ptr%snow)) snow_heat1 = snow_heat1 + snow_tile_heat(ptr%snow)*ptr%frac
  enddo
  call check_conservation ('liquid water', lmass0, lmass1, 1e-6)
  call check_conservation ('frozen water', fmass0, fmass1, 1e-6)
  call check_conservation ('carbon'      , cmass0, cmass1, 1e-6)
  call check_conservation ('canopy air heat content', cana_heat0 , cana_heat1 , 1e-6)
  call check_conservation ('vegetation heat content', vegn_heat0 , vegn_heat1 , 1e-6)
  call check_conservation ('snow heat content',       snow_heat0 , snow_heat1 , 1e-6)
  call check_conservation ('soil heat content',       soil_heat0 , soil_heat1 , 1e-6)
  call check_conservation ('heat content', heat0 , heat1 , 1e-6)

end subroutine 


! ==============================================================================
! given a pointer to a tile and a tile list, insert the tile into the list so that
! if tile can be merged with any one already present, it is merged; otherwise 
! the tile is inserted into the list
subroutine land_tile_merge(tile, list)
  type(land_tile_type), pointer :: tile
  type(land_tile_list_type), intent(inout) :: list

  ! ---- local vars
  type(land_tile_type), pointer :: ptr
  type(land_tile_enum_type) :: ct,et
  
  ! try to find a tile that we can merge to
  ct = first_elmt(list) ; et = tail_elmt(list)
  do while(ct/=et)
     ptr=>current_tile(ct) ; ct = next_elmt(ct)
     if (land_tiles_can_be_merged(tile,ptr)) then
        call merge_land_tiles(tile,ptr)
        call delete_land_tile(tile)
        return ! break out of the subroutine
     endif
  enddo
  ! we reach here only if no suitable files was found in the list
  ! if no suitable tile was found, just insert given tile into the list
  call insert(tile,list)
end subroutine land_tile_merge


! =============================================================================
! splits changing parts of donor tiles into a separate tile list, performing
! land use changes in the process
subroutine split_changing_tile_parts(d_list,d_kind,a_kind,dfrac,a_list)
  type(land_tile_list_type), intent(in) :: d_list ! list of donor tiles
  integer, intent(in) :: d_kind ! kind of donor tiles
  integer, intent(in) :: a_kind ! kind of acceptor tiles
  real,    intent(in) :: dfrac  ! fraction of land area that changes kind
  type(land_tile_list_type), intent(inout) :: a_list ! list of acceptors
  
  ! ---- local vars
  type(land_tile_enum_type) :: ct, et
  type(land_tile_type), pointer :: temp
  type(land_tile_type), pointer :: tile
  real :: area, darea, area0, area1
  real :: x0,x1,x2 ! values of transition intensity
  real, parameter :: eps = 1e-6 ! area calculation precision
  integer :: iter
  real :: factor = 1.6
  integer :: severity ! severity of overshoot errors
  integer :: i,j,k,face ! coordinates of current point, for overshoot diagnostics

  ! calculate total area of the tiles that should be transitioned to another kind
  area = 0
  ct = first_elmt(d_list); et=tail_elmt(d_list)
  do while (ct /= et)
     tile=>current_tile(ct); ct=next_elmt(ct)
     if (.not.associated(tile%vegn)) cycle
     if (tile%vegn%landuse == d_kind)  &
          area = area + tile%frac
  enddo

  ! check for overshoot situtaion: that is, a case where the transition area is
  ! larger than the available area
  if(overshoot_opt /= OPT_IGNORE.and.dfrac>area+overshoot_tolerance) then
     severity = WARNING
     if (overshoot_opt==OPT_STOP) severity = FATAL
     call get_current_point(i,j,k,face)
     call error_mesg('she_landuse',&
          'transition at ('//trim(string(i))//','//trim(string(j))//&
          ',face='//trim(string(face))//&
          ') from "'//trim(landuse_name(d_kind))// &
          '" to "'  //trim(landuse_name(a_kind))//&
          '" ('//trim(string(dfrac))//') is larger than area of "'&
          //trim(landuse_name(d_kind))//'" ('//trim(string(area))//')', &
          severity)
  endif

  ! if area of the tiles of requested kind is zero we cannot transition
  ! anything, so just return
  if (area==0) return
       
  ! transition cannot be more than current total area of specified kind
  darea = min(dfrac, area)

  ! solve equation to get transition intensity
  ! (1) bracket transition intensity interval so that requested area is within it
  x0=0.0; area0 = total_transition_area(d_list, d_kind, a_kind, x0)
  x1=1.0; area1 = total_transition_area(d_list, d_kind, a_kind, x1)
  iter = 0
  do
     if ((area0<=darea).and.(area1>=darea)) exit
     if (area0>darea) then
        x0 = x0-(x1-x0)*factor
        area0 = total_transition_area(d_list, d_kind, a_kind, x0)
     else
        x1 = x1+(x1-x0)*factor
        area1 = total_transition_area(d_list, d_kind, a_kind, x1)
     endif
     iter = iter+1
     if (iter>50) then
        call error_mesg('veg_tile_transitions',&
             'cannot braket transition intensity interval after 50 iterations',&
             FATAL) 
     endif
  enddo

  ! find solution for transition intensity by binary search
  do iter = 1,50
     x2 = (x0+x1)/2
     area = total_transition_area(d_list, d_kind, a_kind, x2)
     if (abs(x1-x2)<eps) exit
     if (area>darea) then
        x1=x2
     else
        x0=x2
     endif
  enddo

  ! do tile transitions to destination list
  ct = first_elmt(d_list); et=tail_elmt(d_list)
  do while (ct /= et)
     tile=>current_tile(ct) ; ct=next_elmt(ct)
     if(.not.associated(tile%vegn)) cycle
     if(tile%vegn%landuse == d_kind) then
        darea = vegn_tran_priority(tile%vegn, a_kind, x2)
        if(darea > 0) then
           ! make a copy 
           temp => new_land_tile(tile)
           temp%frac = tile%frac*darea
           tile%frac = tile%frac*(1.0-darea)
           ! convert land use type of the tile:
           ! cut the forest, if necessary
           if(temp%vegn%landuse==LU_NTRL.or.temp%vegn%landuse==LU_SCND) &
                call vegn_cut_forest(temp%vegn, a_kind)
           ! change landuse type of the tile
           temp%vegn%landuse = a_kind
           ! add the new tile to the resulting list 
           call insert(temp, a_list) ! insert tile into output list
        endif
     endif
  enddo

end subroutine split_changing_tile_parts


! ============================================================================
! calculates total area (fraction of grid cell area) participating in 
! vegetation transition from src_kind to dst_kind for given transition
! intensity tau
function total_transition_area(list,src_kind,dst_kind,tau) result (total_area)
  real :: total_area
  type(land_tile_list_type), intent(in) :: list ! list of tiles
  integer , intent(in) :: src_kind, dst_kind ! source and destination kinds
  real    , intent(in) :: tau                ! transition intensity

  ! ---- local vars
  type(land_tile_enum_type) :: ct, et
  type(land_tile_type), pointer :: tile

  total_area = 0
  ct = first_elmt(list) ; et = tail_elmt(list) 
  do while (ct/=et)
     tile=>current_tile(ct);  ct = next_elmt(ct)
     if (.not.associated(tile%vegn)) cycle ! skip non-vegetated tiles
     if(tile%vegn%landuse == src_kind) &
          total_area = total_area + tile%frac*vegn_tran_priority(tile%vegn,dst_kind,tau)
  enddo
  
end function


! ============================================================================

subroutine get_transitions(time0,time1,k1,k2,tran)
  type(time_type), intent(in) :: time0       ! time of previous calculation of 
    ! transitions (the integral transitinos will be calculated between time0 
    ! and time)
  type(time_type), intent(in) :: time1       ! current time
  integer, intent(in) :: k1,k2               ! kinds of tiles
  type(tran_type), pointer :: tran(:,:,:)    ! transition info

  ! ---- local vars
  integer :: i,j,k,sec,days
  type(tran_type), pointer :: ptr(:,:,:) => NULL()
  real    :: frac(lnd%is:lnd%ie,lnd%js:lnd%je)
  real    :: part_of_year
  logical :: used

  ! allocate array of transitions, if necessary
  if (.not.associated(tran)) then
     allocate(tran(lnd%is:lnd%ie,lnd%js:lnd%je,1) )
  end if

  ! get transition rate for this specific transition
  frac(:,:) = 0.0
  if(input_ids(k1,k2)>0) then
     call integral_transition(time0,time1,input_ids(k1,k2),frac)
     
     do j = lnd%js,lnd%je
     do i = lnd%is,lnd%ie
        if(frac(i,j) == 0) cycle ! skip points where transition rate is zero
        ! find the first empty transition element for the current indices
        k = 1
        do while ( k <= size(tran,3) )
           if(tran(i,j,k)%donor == 0) exit
           k = k+1
        enddo

        if (k>size(tran,3)) then
           ! if there is no room, make the array of transitions larger
           allocate(ptr(lnd%is:lnd%ie,lnd%js:lnd%je,size(tran,3)*2))
           ptr(:,:,1:size(tran,3)) = tran
           deallocate(tran)
           tran => ptr
           nullify(ptr)
        end if

        ! store the transition element
        tran(i,j,k) = tran_type(k1,k2,frac(i,j))
     enddo
     enddo
  endif

  ! send transition data to diagnostics
  if(diag_ids(k1,k2)>0) then
     call get_time(time1-time0, sec,days)
     part_of_year = (days+sec/86400.0)/days_in_year(time0)
     used = send_data(diag_ids(k1,k2),frac/part_of_year,time1)
  endif

end subroutine get_transitions


! ===========================================================================
! given beginnig and end of the time period, and the id of the external
! tile transition field, calculates total transition during the specified period.
! The transition rate data are assumed to be in fraction of land area per year,
! timestamped at the beginning of the year
subroutine integral_transition(t1, t2, id, frac)
  type(time_type), intent(in)  :: t1,t2 ! time boundaries
  integer        , intent(in)  :: id    ! id of the field
  real           , intent(out) :: frac(:,:)

  ! ---- local vars
  integer :: n ! size of time axis
  type(time_type) :: ts,te
  integer         :: i1,i2
  real :: w  ! time interpolation weight
  real :: dt ! current time interval, in years
  real :: sum(size(frac,1),size(frac,2))
  integer :: i,j

  ! adjust the integration limits, in case they are out of range   
  n = size(time_in)
  ts = t1; 
  if (ts<time_in(1)) ts = time_in(1)
  if (ts>time_in(n)) ts = time_in(n) 
  te = t2
  if (te<time_in(1)) te = time_in(1)
  if (te>time_in(n)) te = time_in(n) 

  call time_interp(ts, time_in, w, i1,i2)
  __NF_ASRT__(nfu_get_rec(ncid,id,i1,buffer_in))

  frac = 0;
  call horiz_interp(interp,buffer_in,frac)
  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  sum = -frac*w*dt
  do while(time_in(i2)<=te)
     __NF_ASRT__(nfu_get_rec(ncid,id,i1,buffer_in))
     call horiz_interp(interp,buffer_in,frac)
     dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
     sum = sum+frac*dt
     i2 = i2+1
     i1 = i2-1
     if(i2>size(time_in)) exit ! from loop
  enddo

  call time_interp(te,time_in,w,i1,i2)
  __NF_ASRT__(nfu_get_rec(ncid,id,i1,buffer_in))
  call horiz_interp(interp,buffer_in,frac)
  dt = (time_in(i2)-time_in(i1))//set_time(0,days_in_year((time_in(i2)+time_in(i1))/2))
  frac = sum+frac*w*dt
  do i = 1,size(frac,1)
  do j = 1,size(frac,2)
     if(frac(i,j)<0) then
        call error_mesg('get','transition rate is below zero',FATAL)
     endif
  enddo
  enddo
end subroutine


! ==============================================================================
! checks conservation and aborts with fatal error if tolerance is exceeded
subroutine check_conservation(name, d1, d2, tolerance)
  character(*), intent(in) :: name ! name of the component
  real, intent(in) :: d1,d2 ! values to check
  real, intent(in) :: tolerance ! tolerance of the test

  integer :: curr_i, curr_j, face
  integer :: severity ! severity of the generated message
  character(256) :: message

  if (conservation_opt == OPT_IGNORE) return ! do nothing

  severity = WARNING
  if (overshoot_opt==OPT_STOP) severity = FATAL
  
  if (abs(d1-d2)>tolerance) then
     call get_current_point(i=curr_i,j=curr_j,face=face)
     write(message,'(a,3(x,a,i4), 2(x,a,g))')&
          'conservation of '//trim(name)//' is violated', &
          'at i=',curr_i,'j=',curr_j,'face=',face, &
          'value before=', d1, 'after=', d2
     call error_mesg('land_transitions',message,severity)
  endif
end subroutine 

end module

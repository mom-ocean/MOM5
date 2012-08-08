module static_vegn_mod

use constants_mod,      only : pi
use mpp_mod,            only : mpp_max, mpp_sum
use time_manager_mod,   only : time_type, set_date, time_type_to_real, &
     get_calendar_type, valid_calendar_types, operator(-), get_date
use get_cal_time_mod,   only : get_cal_time

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use fms_mod,            only : write_version_number, error_mesg, FATAL, NOTE, &
     mpp_pe, file_exist, close_file, check_nml_error, stdlog, &
     mpp_root_pe, get_mosaic_tile_file
use time_interp_mod,    only : time_interp
use diag_manager_mod,   only : get_base_date

use nf_utils_mod,       only : nfu_inq_dim, nfu_get_dim, nfu_def_dim, &
     nfu_inq_compressed_var, nfu_get_compressed_rec, nfu_validtype, &
     nfu_get_valid_range, nfu_is_valid, nfu_put_rec, nfu_put_att
use land_data_mod,      only : lnd
use land_io_mod,        only : print_netcdf_error
use land_numerics_mod,  only : nearest
use land_tile_io_mod,   only : create_tile_out_file,sync_nc_files
use land_tile_mod,      only : land_tile_type, land_tile_enum_type, first_elmt, &
     tail_elmt, next_elmt, current_tile, operator(/=), nitems
use vegn_cohort_mod,    only : vegn_cohort_type
use cohort_io_mod,      only : create_cohort_dimension, &
     write_cohort_data_i0d_fptr, write_cohort_data_r0d_fptr


implicit none
private

! ==== public interface =====================================================
public :: read_static_vegn_namelist
public :: static_vegn_init
public :: static_vegn_end

public :: read_static_vegn
public :: write_static_vegn
! ==== end of public interface ==============================================

! ==== module constants =====================================================
character(len=*), parameter :: &
     module_name = 'static_vegn_mod', &
     version     = '$Id: vegn_static_override.F90,v 19.0 2012/01/06 20:44:38 fms Exp $', &
     tagname     = '$Name: siena_201207 $'

! ==== module data ==========================================================
logical :: module_is_initialized = .FALSE.
integer :: ncid  ! netcdf id of the input file
integer :: ncid2 ! netcdf id of the output file
type(time_type),allocatable :: time_line(:) ! time line of input data
type(time_type)             :: ts,te        ! beginning and end of time interval
integer, allocatable :: map_i(:,:), map_j(:,:)! remapping arrays: for each of the
     ! land grid cells in current domain they hold indices of corresponding points 
     ! in the input grid.
type(time_type) :: base_time ! model base time for static vegetation output

! ---- namelist variables ---------------------------------------------------
logical :: use_static_veg = .FALSE.
character(len=512) :: input_file = & ! name of input file for static vegetation
     "INPUT/static_veg_data.nc"
character(len=10)  :: timeline   = 'normal' ! type of timeline ('normal' or 'loop')
integer, dimension(6) :: &
     start_loop = (/1,1,1,0,0,0/), & ! beginning of the time loop
     end_loop   = (/1,1,1,0,0,0/)    ! end of the time loop
logical :: fill_land_mask = .FALSE. ! if true, all the vegetation points on the
     ! map are filled with the information from static vegetation data, using
     ! nearest point remap; otherwise only the points that overlap with valid
     ! static vegetation data are overriden.
logical :: write_static_veg = .FALSE. ! if true, the state of vegetation is saved 
     ! periodically for future use as static vegetation input
character(16) :: static_veg_freq = 'daily' ! or 'monthly', or 'annual'
     ! specifies the frequency for writing the static vegetation data file
namelist/static_veg_nml/use_static_veg,input_file,timeline,start_loop,end_loop,&
     fill_land_mask, write_static_veg, static_veg_freq

! ==== NetCDF declarations ==================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

contains

! ===========================================================================
subroutine read_static_vegn_namelist(static_veg_used)
  logical, intent(out) :: static_veg_used

  ! ---- local vars
  integer :: unit, ierr, io

  call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=static_veg_nml, iostat=io)
  ierr = check_nml_error(io, 'static_veg_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=static_veg_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'static_veg_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=static_veg_nml)
  endif

  if (    (trim(static_veg_freq)=='daily') &
      .or.(trim(static_veg_freq)=='monthly') &
      .or.(trim(static_veg_freq)=='annual') ) then
     ! static_veg_freq is OK -- do nothing
  else
     call error_mesg('static_vegn_init','option static_veg_freq="'&
          //trim(static_veg_freq)&
          //'" is invalid, use "daily", "monthly", or "annual"', FATAL)
  endif

  static_veg_used = use_static_veg
end subroutine read_static_vegn_namelist


! ===========================================================================
subroutine static_vegn_init()

  ! ---- local vars
  integer :: unlimdim, timelen, timeid
  integer :: i,j,k,iret
  character(len=NF_MAX_NAME) :: dimname  ! name of the dimension variable : time, lon, and lat
  integer                    :: ndims    ! rank of input vars
  integer                    :: dimids (NF_MAX_VAR_DIMS) ! netcdf IDs of input var dimensions
  integer                    :: dimlens(NF_MAX_VAR_DIMS) ! sizes of respective dimensions
  real, allocatable          :: t(:)     ! temporary real timeline
  character(len=256)         :: units    ! units of time in the file
  character(len=256)         :: calendar ! calendar of the data
  real, allocatable          :: in_lon(:)! longitude coordinates in input file
  real, allocatable          :: in_lat(:)! latitude coordinates in input file
  logical, allocatable       :: mask(:,:)! mask of valid points in input data 
  integer, allocatable       :: data(:,:,:,:) ! temprary array used to calculate the mask of
                                         ! valid input data
  logical                    :: has_records ! true if input variable has records
  integer :: tile_dim_length ! length of tile dimension in output files 
                             ! global max of number of tiles per gridcell
  integer :: year, month, day, hour, minute, sec ! components of base date
  character(len=1024) :: actual_input_file
  logical :: input_is_multiface ! TRUE if the input files are face-specific

  if(module_is_initialized) return


  if(use_static_veg) then

     ! SET UP LOOP BOUNDARIES
     ts = set_date(start_loop(1),start_loop(2),start_loop(3), start_loop(4),start_loop(5),start_loop(6))
     te = set_date(end_loop(1)  ,end_loop(2)  ,end_loop(3)  , end_loop(4)  ,end_loop(5)  ,end_loop(6)  )
     
     ! OPEN INPUT FILE
     if (nf_open(input_file,NF_NOWRITE,ncid)/=NF_NOERR) then
        if(lnd%nfaces==1) then
           ! for 1-face grid we can't use multi-face input, even if it exists
           call error_mesg('static_vegn_init','input file "'//trim(input_file)&
                   //'" does not exist', FATAL)
        else
           ! if there's more then one face, try opening face-specific input
           call get_mosaic_tile_file(trim(input_file),actual_input_file,.FALSE.,lnd%domain)
           if (nf_open(actual_input_file,NF_NOWRITE,ncid)/=NF_NOERR) then
              call error_mesg('static_vegn_init','Neither "'//trim(input_file)&
                   //'" nor "'//trim(actual_input_file)//'" files exist', FATAL)
           else
              call error_mesg('static_vegn_init','Reading face-specific vegetation file "'&
                   //trim(actual_input_file)//'"', NOTE)
              input_is_multiface = .TRUE.
           endif
        endif
     else
        call error_mesg('static_vegn_init','Reading global static vegetation file "'&
             //trim(input_file)//'"', NOTE)
        input_is_multiface = .FALSE.
     endif
     
     ! READ TIME AXIS DATA
     __NF_ASRT__(nf_inq_unlimdim( ncid, unlimdim ))
     __NF_ASRT__(nf_inq_dimname ( ncid, unlimdim, dimname ))
     __NF_ASRT__(nf_inq_varid   ( ncid, dimname, timeid ))
     __NF_ASRT__(nf_inq_dimlen( ncid, unlimdim, timelen ))
     allocate (time_line(timelen), t(timelen))
     __NF_ASRT__(nf_get_var_double (ncid, timeid, t ))
     
     ! GET UNITS OF THE TIME
     units = ' '
     __NF_ASRT__(nf_get_att_text(ncid, timeid,'units',units))
     
     ! GET CALENDAR OF THE DATA
     calendar = ' '
     iret = nf_get_att_text(ncid, timeid, 'calendar',calendar)
     if(iret/=NF_NOERR) &
          iret = nf_get_att_text(ncid, timeid,'calendar_type',calendar)
     if(iret/=NF_NOERR) &
          calendar='JULIAN' ! use model calendar? how to get the name of the model calendar?
       
     ! CONVERT TIME TO THE FMS TIME_TYPE AND STORE IT IN THE TIMELINE FOR THE
     ! DATA SET
     do i = 1, size(t)
        ! set the respective value in the timeline
        time_line(i) = get_cal_time(t(i),units,calendar)
     enddo

     ! READ HORIZONTAL COORDINATES
     iret = nfu_inq_compressed_var(ncid,'species',ndims=ndims,dimids=dimids,dimlens=dimlens,&
          has_records=has_records)
     __NF_ASRT__(iret)
     allocate(in_lon(dimlens(1)),in_lat(dimlens(2)))
     __NF_ASRT__(nfu_get_dim(ncid,dimids(1),in_lon)) ! get longitude
     __NF_ASRT__(nfu_get_dim(ncid,dimids(2),in_lat)) ! get latitude
     in_lon = in_lon*PI/180.0 ; in_lat = in_lat*PI/180.0

     ! COMPUTE INDEX REMAPPING ARRAY
     allocate(map_i(lnd%is:lnd%ie,lnd%js:lnd%je))
     allocate(map_j(lnd%is:lnd%ie,lnd%js:lnd%je))
     allocate(mask(size(in_lon),size(in_lat)))

     if(fill_land_mask) then
        ! CALCULATE THE DIMENSIONS OF THE BUFFER FOR THE INPUT DATA
        if (has_records) ndims=ndims-1
        do i = ndims+1,4
           dimlens(i) = 1
        enddo
        ! READ THE FIRST RECORD AND CALCULTE THE MASK OF THE VALID INPUT DATA
        allocate(data(dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
        !             lon        lat        tile       cohort
        data(:,:,:,:) = -1
        __NF_ASRT__(nfu_get_compressed_rec(ncid,'species',1,data))
        do j = 1,size(data,2)
        do i = 1,size(data,1)
           mask(i,j) = any(data(i,j,:,:)>=0)
        enddo
        enddo
        deallocate(data)
     else
        mask(:,:) = .TRUE.
     endif

     if(input_is_multiface) then
        ! check that the sizes of input data and the model data are the same
        if(dimlens(1)/=lnd%nlon.or.dimlens(2)/=lnd%nlat) then
           call error_mesg('static_vegn_init','size of face-specific static vegetation '&
                //'data isn''t the same as the size of the mosaic face', FATAL)
        endif
        ! in case of multi-face input, we don't do any remapping
        do j = lnd%js,lnd%je
        do i = lnd%is,lnd%ie
           map_i(i,j) = i; map_j(i,j) = j
        enddo
        enddo
     else
        ! do the nearest-point remapping
        do j = lnd%js,lnd%je
        do i = lnd%is,lnd%ie
           call nearest(mask,in_lon,in_lat,lnd%lon(i,j),lnd%lat(i,j),map_i(i,j),map_j(i,j))
        enddo
        enddo
     endif

     deallocate (in_lon,in_lat,mask)
     deallocate(t)
  endif

  if(write_static_veg) then
     ! create output file for static vegetation

     ! count all land tiles and determine the lenght of tile dimension
     ! sufficient for the current domain
     tile_dim_length = 0
     do j = lnd%js, lnd%je
     do i = lnd%is, lnd%ie
        k = nitems(lnd%tile_map(i,j))
        tile_dim_length = max(tile_dim_length,k)
     enddo
     enddo
   
     ! [1.1] calculate the tile dimension length by taking the max across all domains
     call mpp_max(tile_dim_length)

     call create_tile_out_file(ncid2,'static_veg_out.nc', &
          lnd%coord_glon, lnd%coord_glat, vegn_tile_exists, tile_dim_length)
     ! create compressed dimension for vegetation cohorts
     call create_cohort_dimension(ncid2)
     ! get the base date of the simulation
     call get_base_date(year,month,day,hour,minute,sec)
     base_time = set_date(year, month, day, hour, minute, sec)
     if(mpp_pe()==lnd%io_pelist(1)) then
        ! create time axis, on root IO processors only
        units = ' '
        write(units, 11) year, month, day, hour, minute, sec
11      format('days since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
        __NF_ASRT__(nfu_def_dim(ncid2,'time',NF_UNLIMITED,NF_DOUBLE,units=trim(units)))
        ! add calendar attribute to the time axis
        iret=nfu_put_att(ncid2,'time','calendar',&
             trim(valid_calendar_types(get_calendar_type())))
        __NF_ASRT__(iret)
     endif
     call sync_nc_files(ncid2)
  endif
  module_is_initialized = .true.

end subroutine static_vegn_init

! ===========================================================================
subroutine static_vegn_end()
  if(use_static_veg) then
     __NF_ASRT__(nf_close(ncid))
     deallocate(time_line,map_i,map_j)
  endif
  if(write_static_veg) then
     __NF_ASRT__(nf_close(ncid2))
  endif
  module_is_initialized = .false.
end subroutine static_vegn_end

! ===========================================================================
subroutine read_static_vegn (time)
  type(time_type), intent(in)    :: time

  ! ---- local vars 
  integer :: index1, index2 ! result of time interpolation (only index1 is used)
  real    :: weight         ! another result of time interp, not used
  character(len=256) :: err_msg

  if(.not.use_static_veg)return;

  !   time_interp to find out the index of the current time interval
  if (timeline == 'loop') then
     err_msg = ''
     call time_interp(time, ts, te, time_line, weight, index1, index2, &
                      correct_leap_year_inconsistency=.true.)
     if(err_msg /= '') then
       call error_mesg('subroutine read_static_vegn',trim(err_msg), FATAL)
     endif
  else if (timeline == 'normal') then
     call time_interp(time, time_line, weight, index1, index2)
  else
     call error_mesg(module_name,'timeline option "'//trim(timeline)// &
          '" is incorrect, use "normal" or "loop"', FATAL)
  endif

  ! read the data into cohort variables
  call read_remap_cohort_data_i0d_fptr(ncid, 'species' , cohort_species_ptr , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bl'      , cohort_bl_ptr      , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'blv'     , cohort_blv_ptr     , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'br'      , cohort_br_ptr      , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bsw'     , cohort_bsw_ptr     , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bwood'   , cohort_bwood_ptr   , map_i, map_j, index1)
  call read_remap_cohort_data_r0d_fptr(ncid, 'bliving' , cohort_bliving_ptr , map_i, map_j, index1)
  call read_remap_cohort_data_i0d_fptr(ncid, 'status'  , cohort_status_ptr  , map_i, map_j, index1)

  ! derived variables will be updated in update_land_bc_fast
end subroutine read_static_vegn


! ===========================================================================
subroutine write_static_vegn()

  real :: t ! time in output units
  integer :: rec ! number of record to write
  ! components of the date
  integer :: second, minute, hour, day0, day1, month0, month1, year0, year1

  if(.not.write_static_veg) return;

  ! get components of calendar dates for this and previous time step
  call get_date(lnd%time,             year0,month0,day0,hour,minute,second)
  call get_date(lnd%time-lnd%dt_fast, year1,month1,day1,hour,minute,second)

  if (     (trim(static_veg_freq)=='daily'  .and.  day1/=day0)   &
       .or.(trim(static_veg_freq)=='monthly'.and.month1/=month0) &
       .or.(trim(static_veg_freq)=='annual' .and. year1/=year0) )&
       then
     ! sync output files with the disk so that every processor sees the same 
     ! information, number of records being critical here
     call sync_nc_files(ncid2)
     ! get the current number of records in the output file
     __NF_ASRT__(nfu_inq_dim(ncid2,'time',rec))
     rec = rec+1
     ! create new record in the output file and store current value of time
     if(mpp_pe()==lnd%io_pelist(1)) then
        t = (time_type_to_real(lnd%time)-time_type_to_real(base_time))/86400
        __NF_ASRT__(nfu_put_rec(ncid2,'time',rec,t))
     endif
     ! write static vegetation data
     call write_cohort_data_i0d_fptr(ncid2,'species', cohort_species_ptr, &
          'vegetation species',record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'bl',      cohort_bl_ptr, &
          'biomass of leaves per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'blv',     cohort_blv_ptr, &
          'biomass of virtual leaves (labile store) per individual','kg C/m2',record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'br',      cohort_br_ptr, &
          'biomass of fine roots per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'bsw',     cohort_bsw_ptr, &
          'biomass of sapwood per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'bwood',   cohort_bwood_ptr, &
          'biomass of heartwood per individual','kg C/m2', record=rec)
     call write_cohort_data_r0d_fptr(ncid2,'bliving', cohort_bliving_ptr, &
          'total living biomass per individual','kg C/m2', record=rec)
     call write_cohort_data_i0d_fptr(ncid2,'status',  cohort_status_ptr, &
          'leaf status', record=rec)
  endif
end subroutine write_static_vegn


! ============================================================================
#define F90_TYPE       integer
#define READ_REMAP_SUB read_remap_cohort_data_i0d_fptr
#include "read_remap_cohort_data.inc"

#define F90_TYPE       real
#define READ_REMAP_SUB read_remap_cohort_data_r0d_fptr
#include "read_remap_cohort_data.inc"

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;\
end subroutine

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(integer,status)

end module static_vegn_mod

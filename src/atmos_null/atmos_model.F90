module atmos_model_mod

!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang
!</CONTACT>
!
!<OVERVIEW>
! Null atmosphere model. 
!</OVERVIEW>
!<DESCRIPTION>
! Null atmosphere model. 
!</DESCRIPTION>
!
!<NAMELIST NAME="atmos_model_nml">
!  <DATA NAME="layout" TYPE="integer">
!  Processor domain layout for atmos model. 
!  </DATA> 
!  <DATA NAME="mask_table" TYPE="character">
!   A text file to specify n_mask, layout and mask_list to reduce number of processor
!   usage by masking out some domain regions which contain all land points. 
!   The default file name of mask_table is "INPUT/atmos_mask_table". Please note that 
!   the file name must begin with "INPUT/". The first 
!   line of mask_table will be number of region to be masked out. The second line 
!   of the mask_table will be the layout of the model. User need to set atmos_model_nml
!   variable layout to be the same as the second line of the mask table.
!   The following n_mask line will be the position of the processor to be masked out.
!   The mask_table could be created by tools check_mask. 
!   For example the mask_table will be as following if n_mask=2, layout=4,6 and 
!   the processor (1,2) and (3,6) will be masked out. 
!     2
!     4,6
!     1,2
!     3,6
!  </DATA>
!

!</NAMELIST>




use mpp_mod,           only : mpp_npes, mpp_pe, mpp_error, FATAL, mpp_chksum
use mpp_mod,           only : input_nml_file

use mpp_domains_mod,   only : domain2d
use mpp_domains_mod,   only : mpp_define_layout, mpp_define_domains
use mpp_domains_mod,   only : CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain
use mpp_domains_mod,   only : mpp_get_compute_domain, mpp_get_tile_id
use mpp_domains_mod,   only : mpp_get_current_ntile
use fms_mod,           only : field_exist, read_data, field_size, stdout
use fms_mod,           only : open_namelist_file, check_nml_error, file_exist, close_file
use fms_io_mod,        only : parse_mask_table
use time_manager_mod,  only : time_type
use coupler_types_mod, only : coupler_2d_bc_type
use diag_manager_mod,  only : diag_axis_init
use constants_mod,     only : cp_air, hlv
use mosaic_mod,        only : get_mosaic_ntiles
use xgrid_mod,         only : grid_box_type
use grid_mod,          only : get_grid_ntiles, define_cube_mosaic
use grid_mod,          only : get_grid_size, get_grid_cell_vertices
use grid_mod,          only : get_grid_cell_centers
use diag_integral_mod, only : diag_integral_init
use tracer_manager_mod,only : register_tracers
use field_manager_mod, only : MODEL_LAND

implicit none
private

public atmos_data_type
public atmos_model_end
public atmos_model_init
public ice_atmos_boundary_type
public land_ice_atmos_boundary_type
public land_atmos_boundary_type
public surf_diff_type
public update_atmos_model_down
public update_atmos_model_up
public atm_stock_pe
public atmos_model_restart
public atmos_data_type_chksum
public lnd_ice_atm_bnd_type_chksum, lnd_atm_bnd_type_chksum
public ice_atm_bnd_type_chksum

!<PUBLICTYPE >
! This type should be defined in one spot and "used" from there
type surf_diff_type
  real, pointer, dimension(:,:)   :: dtmass  => NULL()
  real, pointer, dimension(:,:)   :: dflux_t => NULL()
  real, pointer, dimension(:,:)   :: delta_t => NULL()
  real, pointer, dimension(:,:)   :: delta_u => NULL()
  real, pointer, dimension(:,:)   :: delta_v => NULL()
  real, pointer, dimension(:,:,:) :: dflux_tr => NULL()   ! tracer flux tendency
  real, pointer, dimension(:,:,:) :: delta_tr => NULL()   ! tracer tendency
  real, pointer, dimension(:,:)     :: sst_miz => NULL()
end type surf_diff_type
!</PUBLICTYPE >

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real, pointer, dimension(:,:) :: lon_bnd  => NULL() ! local longitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: lat_bnd  => NULL() ! local latitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: t_bot    => NULL() ! temperature at lowest model level
     real, pointer, dimension(:,:,:) :: tr_bot => NULL() ! tracers at lowest model level, including specific humidity
     real, pointer, dimension(:,:) :: z_bot    => NULL() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => NULL() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => NULL() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => NULL() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => NULL() ! surface pressure 
     real, pointer, dimension(:,:) :: slp      => NULL() ! sea level pressure 
     real, pointer, dimension(:,:) :: gust     => NULL() ! gustiness factor
     real, pointer, dimension(:,:) :: coszen   => NULL() ! cosine of the zenith angle
     real, pointer, dimension(:,:) :: flux_sw  => NULL() ! net shortwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: flux_sw_dir            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_dif            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dir   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dif   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dir =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dif =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dir        =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dif        =>NULL()
     real, pointer, dimension(:,:) :: flux_lw  => NULL() ! net longwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: lprec    => NULL() ! mass of liquid precipitation since last time step (Kg/m2)
     real, pointer, dimension(:,:) :: fprec    => NULL() ! ass of frozen precipitation since last time step (Kg/m2)
     logical,pointer,dimension(:,:):: maskmap  => NULL() ! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used.indicate if a domain region will be loaded.
     type (surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer, pointer              :: pelist(:) =>NULL() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     type(coupler_2d_bc_type)      :: fields   ! array of fields used for additional tracers
     type(grid_box_type)           :: grid
 end type
!</PUBLICTYPE >

!<PUBLICTYPE >
type land_ice_atmos_boundary_type
   ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
   real, dimension(:,:),   pointer :: t              =>NULL() ! surface temperature for radiation calculations
   real, dimension(:,:),   pointer :: albedo         =>NULL() ! surface albedo for radiation calculations
   real, dimension(:,:),   pointer :: albedo_vis_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_vis_dif =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dif =>NULL()
   real, dimension(:,:),   pointer :: land_frac      =>NULL() ! fraction amount of land in a grid box 
   real, dimension(:,:),   pointer :: dt_t           =>NULL() ! temperature tendency at the lowest level
   real, dimension(:,:,:), pointer :: dt_tr          =>NULL() ! tracer tendency at the lowest level, including specific humidity
   real, dimension(:,:),   pointer :: u_flux         =>NULL() ! zonal wind stress
   real, dimension(:,:),   pointer :: v_flux         =>NULL() ! meridional wind stress
   real, dimension(:,:),   pointer :: dtaudu         =>NULL() ! derivative of wind stress w.r.t. the lowest level wind speed
   real, dimension(:,:),   pointer :: dtaudv         =>NULL() ! derivative of wind stress w.r.t. the lowest level wind speed
   real, dimension(:,:),   pointer :: u_star         =>NULL() ! friction velocity
   real, dimension(:,:),   pointer :: b_star         =>NULL() ! bouyancy scale
   real, dimension(:,:),   pointer :: q_star         =>NULL() ! moisture scale
   real, dimension(:,:),   pointer :: rough_mom      =>NULL() ! surface roughness (used for momentum)
   real, dimension(:,:,:), pointer :: data =>NULL() !collective field for "named" fields above
   real, dimension(:,:),   pointer :: frac_open_sea  =>null() ! non-seaice fraction (%)
   integer :: xtype             !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type :: land_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from land alone to atmos (none at present)
end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
type :: ice_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from ice alone to atmos (none at present)
end type ice_atmos_boundary_type
!</PUBLICTYPE >
  
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmos_model.F90,v 20.0 2013/12/13 23:08:53 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!---- atmos_model_nml
integer :: layout(2)
! mask_table contains information for masking domain ( n_mask, layout and mask_list).
character(len=128) :: mask_table = "INPUT/atmos_mask_table"

namelist /atmos_model_nml/ layout, mask_table


contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation, 
!   vertical diffusion of momentum, tracers, and heat/moisture.
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </IN>

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
  type (atmos_data_type), intent(in) :: Atmos

  return

end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero. 
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </IN>
subroutine update_atmos_model_up( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!-----------------------------------------------------------------------

   type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
   type (atmos_data_type), intent(in) :: Atmos
   
   return

end subroutine update_atmos_model_up
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file. 
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in)          :: Time_init, Time, Time_step

real, dimension(:,:), allocatable     :: glon, glat, glon_bnd, glat_bnd
integer, dimension(:), allocatable    :: tile_ids
real, dimension(:,:),  allocatable    :: area
integer                               :: is, ie, js, je, i
integer                               :: nlon, nlat, ntile, tile
integer                               :: ntracers, ntprog, ndiag
integer                               :: namelist_unit, io, ierr, stdoutunit

!--- read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, atmos_model_nml, iostat=io)
#else
   namelist_unit = open_namelist_file()
   ierr=1
   do while (ierr /= 0)
     read(namelist_unit, nml=atmos_model_nml, iostat=io, end=20)
     ierr = check_nml_error (io, 'atmos_model_nml')
   enddo
   20 call close_file (namelist_unit)
#endif

stdoutunit = stdout()
!---- set the atmospheric model time ------

Atmos % Time_init = Time_init
Atmos % Time      = Time
Atmos % Time_step = Time_step
 
call get_grid_ntiles('ATM',ntile)
call get_grid_size('ATM',1,nlon,nlat)


if(file_exist(mask_table)) then
   if(ntile > 1) then
      call mpp_error(FATAL, &
        'atmos_model_init: file '//trim(mask_table)//' should not exist when ntile is not 1')
   endif
   if(layout(1) == 0 .OR. layout(2) == 0 ) call mpp_error(FATAL, &
        'atmos_model_init: atmos_model_nml layout should be set when file '//trim(mask_table)//' exists')

   write(stdoutunit, '(a)') '==> NOTE from atmos_model_init:  reading maskmap information from '//trim(mask_table)
   allocate(Atmos%maskmap(layout(1), layout(2)))
   call parse_mask_table(mask_table, Atmos%maskmap, "Atmos model")
else
   if(layout(1)*layout(2) .NE. mpp_npes()/ntile) call mpp_define_layout((/1,nlon,1,nlat/), mpp_npes()/ntile, layout)
endif


if(ntile ==1) then
   if( ASSOCIATED(Atmos%maskmap) ) then
      call mpp_define_domains((/1,nlon,1,nlat/), layout, Atmos%domain, &
           xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, maskmap = Atmos%maskmap , name='atmos model')
   else
      call mpp_define_domains((/1,nlon,1,nlat/), layout, Atmos%domain, &
           xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, name='atmos model')
   end if
else
   call define_cube_mosaic('ATM', Atmos%domain, layout, halo=1)
endif

call mpp_get_compute_domain(Atmos%domain,is,ie,js,je)

allocate ( glon_bnd(nlon+1,nlat+1))
allocate ( glat_bnd(nlon+1,nlat+1))
allocate ( glon(nlon, nlat))
allocate ( glat(nlon, nlat))
allocate ( Atmos%lon_bnd(ie-is+2,je-js+2) )
allocate ( Atmos%lat_bnd(ie-is+2,je-js+2) )
allocate ( area(ie-is+1,je-js+1) )

allocate(tile_ids(mpp_get_current_ntile(Atmos%domain)))
tile_ids = mpp_get_tile_id(Atmos%domain)
tile = tile_ids(1)
deallocate(tile_ids)

call get_grid_cell_vertices('ATM',tile,glon_bnd,glat_bnd)
call get_grid_cell_centers ('ATM',tile,glon, glat)

Atmos % lon_bnd(:,:) = glon_bnd(is:ie+1, js:je+1)*atan(1.0)/45.0
Atmos % lat_bnd(:,:) = glat_bnd(is:ie+1, js:je+1)*atan(1.0)/45.0

if(ntile==1) then
   Atmos%axes(1) = diag_axis_init('lon',glon(:,1),'degrees_E','X','longitude',&
        set_name='atmos',domain2 = Atmos%domain)

   Atmos%axes(2) = diag_axis_init('lat',glat(1,:),'degrees_N','Y','latitude',&
        set_name='atmos',domain2 = Atmos%domain)  
else
   Atmos%axes(1) = diag_axis_init('lon',(/(real(i),i=1,nlon)/),'degrees_E','X','longitude',&
        set_name='atmos',domain2 = Atmos%domain)

   Atmos%axes(2) = diag_axis_init('lat',(/(real(i),i=1,nlat)/),'degrees_N','Y','latitude',&
        set_name='atmos',domain2 = Atmos%domain)  
endif

call register_tracers(MODEL_LAND, ntracers, ntprog, ndiag)
allocate ( Atmos % t_bot(is:ie,js:je) )
allocate ( Atmos % tr_bot(is:ie,js:je,ntprog) )
allocate ( Atmos % z_bot(is:ie,js:je) )
allocate ( Atmos % p_bot(is:ie,js:je) )
allocate ( Atmos % u_bot(is:ie,js:je) )
allocate ( Atmos % v_bot(is:ie,js:je) )
allocate ( Atmos % p_surf(is:ie,js:je) )
allocate ( Atmos % slp(is:ie,js:je) )
allocate ( Atmos % gust(is:ie,js:je) )
allocate ( Atmos % coszen(is:ie,js:je) )
allocate ( Atmos % flux_sw(is:ie,js:je) )
allocate ( Atmos % flux_sw_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_vis_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_vis_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_total_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_total_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis_dif(is:ie,js:je) )
allocate ( Atmos % flux_lw(is:ie,js:je) )
allocate ( Atmos % lprec(is:ie,js:je) )
allocate ( Atmos % fprec(is:ie,js:je) )

Atmos % t_bot=273.0
Atmos % tr_bot = 0.0
Atmos % z_bot = 10.0
Atmos % p_bot = .99e5
Atmos % u_bot = 0.0
Atmos % v_bot = 0.0
Atmos % p_surf = 1.e5
Atmos % slp = 1.e5
Atmos % gust = 0.0
Atmos % coszen = 0.0
Atmos % flux_sw = 0.0
Atmos % flux_lw = 0.0
Atmos % flux_sw_dir = 0.0
Atmos % flux_sw_dif = 0.0 
Atmos % flux_sw_down_vis_dir = 0.0 
Atmos % flux_sw_down_vis_dif = 0.0 
Atmos % flux_sw_down_total_dir = 0.0
Atmos % flux_sw_down_total_dif = 0.0
Atmos % flux_sw_vis = 0.0 
Atmos % flux_sw_vis_dir = 0.0 
Atmos % flux_sw_vis_dif = 0.0
Atmos % lprec = 0.0
Atmos % fprec = 0.0

allocate ( Atmos % Surf_diff % dtmass(is:ie, js:je) )
allocate ( Atmos % Surf_diff % dflux_t(is:ie, js:je) )
allocate ( Atmos % Surf_diff % delta_t(is:ie, js:je) )
allocate ( Atmos % Surf_diff % delta_u(is:ie, js:je) )
allocate ( Atmos % Surf_diff % delta_v(is:ie, js:je) )
allocate ( Atmos % Surf_diff % dflux_tr(is:ie, js:je, ntprog) )
allocate ( Atmos % Surf_diff % delta_tr(is:ie, js:je, ntprog) )

Atmos % Surf_diff % dtmass   = 0.0
Atmos % Surf_diff % dflux_t  = 0.0
Atmos % Surf_diff % delta_t  = 0.0
Atmos % Surf_diff % delta_u  = 0.0
Atmos % Surf_diff % delta_v  = 0.0
Atmos % Surf_diff % dflux_tr = 0.0
Atmos % Surf_diff % delta_tr = 0.0

area = 0.0

allocate ( Atmos % grid % dx    (   is:ie  , js:je+1))
allocate ( Atmos % grid % dy    (   is:ie+1, js:je  ))
allocate ( Atmos % grid % area  (   is:ie  , js:je  ))
allocate ( Atmos % grid % edge_w(            js:je+1))
allocate ( Atmos % grid % edge_e(            js:je+1))
allocate ( Atmos % grid % edge_s(   is:ie+1         ))
allocate ( Atmos % grid % edge_n(   is:ie+1         ))
allocate ( Atmos % grid % en1   (3, is:ie  , js:je+1))
allocate ( Atmos % grid % en2   (3, is:ie+1, js:je  ))
allocate ( Atmos % grid % vlon  (3, is:ie  , js:je  ))
allocate ( Atmos % grid % vlat  (3, is:ie  , js:je  ))
     
Atmos % grid % dx    = 1.0
Atmos % grid % dy    = 1.0
Atmos % grid % area  = 1.0
Atmos % grid % edge_w= 0.0
Atmos % grid % edge_e= 1.0
Atmos % grid % edge_s= 0.0
Atmos % grid % edge_n= 1.0
Atmos % grid % en1   = 0.0
Atmos % grid % en2   = 0.0
Atmos % grid % vlon  = 0.0
Atmos % grid % vlat  = 0.0

call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                         Atmos % lon_bnd(:,:),  &
                         Atmos % lat_bnd(:,:), area)

end subroutine atmos_model_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </IN>

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(in) :: Atmos

return

end subroutine atmos_model_end
! </SUBROUTINE>

  !#######################################################################
  ! <SUBROUTINE NAME="atmos_model_restart">
  ! <DESCRIPTION>
  ! dummy routines.
  ! </DESCRIPTION>
  subroutine atmos_model_restart(Atmos, timestamp)
    type (atmos_data_type),   intent(inout) :: Atmos
    character(len=*),  intent(in)           :: timestamp


  end subroutine atmos_model_restart
  ! </SUBROUTINE>

subroutine atm_stock_pe (Atm, index, value)

type (atmos_data_type), intent(inout) :: Atm
integer,                intent(in)    :: index
real,                   intent(out)   :: value

value = 0.0

end subroutine atm_stock_pe

!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm 
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    integer :: n, outunit

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd               )
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd               )
  write(outunit,100) ' atm%t_bot                  ', mpp_chksum(atm%t_bot                 )
  do n = 1, size(atm%tr_bot,3)
  write(outunit,100) ' atm%tr_bot(:,:,n)          ', mpp_chksum(atm%tr_bot(:,:,n)         )
  enddo
  write(outunit,100) ' atm%z_bot                  ', mpp_chksum(atm%z_bot                 )
  write(outunit,100) ' atm%p_bot                  ', mpp_chksum(atm%p_bot                 )
  write(outunit,100) ' atm%u_bot                  ', mpp_chksum(atm%u_bot                 )
  write(outunit,100) ' atm%v_bot                  ', mpp_chksum(atm%v_bot                 )
  write(outunit,100) ' atm%p_surf                 ', mpp_chksum(atm%p_surf                )
  write(outunit,100) ' atm%slp                    ', mpp_chksum(atm%slp                   )
  write(outunit,100) ' atm%gust                   ', mpp_chksum(atm%gust                  )
  write(outunit,100) ' atm%coszen                 ', mpp_chksum(atm%coszen                )
  write(outunit,100) ' atm%flux_sw                ', mpp_chksum(atm%flux_sw               )
  write(outunit,100) ' atm%flux_sw_dir            ', mpp_chksum(atm%flux_sw_dir           )
  write(outunit,100) ' atm%flux_sw_dif            ', mpp_chksum(atm%flux_sw_dif           )
  write(outunit,100) ' atm%flux_sw_down_vis_dir   ', mpp_chksum(atm%flux_sw_down_vis_dir  )
  write(outunit,100) ' atm%flux_sw_down_vis_dif   ', mpp_chksum(atm%flux_sw_down_vis_dif  )
  write(outunit,100) ' atm%flux_sw_down_total_dir ', mpp_chksum(atm%flux_sw_down_total_dir)
  write(outunit,100) ' atm%flux_sw_down_total_dif ', mpp_chksum(atm%flux_sw_down_total_dif)
  write(outunit,100) ' atm%flux_sw_vis            ', mpp_chksum(atm%flux_sw_vis           )
  write(outunit,100) ' atm%flux_sw_vis_dir        ', mpp_chksum(atm%flux_sw_vis_dir       )
  write(outunit,100) ' atm%flux_sw_vis_dif        ', mpp_chksum(atm%flux_sw_vis_dif       )
  write(outunit,100) ' atm%flux_lw                ', mpp_chksum(atm%flux_lw               )
  write(outunit,100) ' atm%lprec                  ', mpp_chksum(atm%lprec                 )
  write(outunit,100) ' atm%fprec                  ', mpp_chksum(atm%fprec                 )
!  call surf_diff_type_chksum(id, timestep, atm%surf_diff)

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine lnd_ice_atm_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_ice_Atm_bnd_type):: ', id, timestep
100 FORMAT("CHECKSUM::",A32," = ",Z20)
    write(outunit,100) 'lnd_ice_atm_bnd_type%t             ',mpp_chksum(bnd_type%t              )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo        ',mpp_chksum(bnd_type%albedo         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dir',mpp_chksum(bnd_type%albedo_vis_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dir',mpp_chksum(bnd_type%albedo_nir_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dif',mpp_chksum(bnd_type%albedo_vis_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dif',mpp_chksum(bnd_type%albedo_nir_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%land_frac     ',mpp_chksum(bnd_type%land_frac      )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_t          ',mpp_chksum(bnd_type%dt_t           )
    do n = 1, size(bnd_type%dt_tr,3)
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_tr(:,:,n)  ',mpp_chksum(bnd_type%dt_tr(:,:,n)   )
    enddo
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_flux        ',mpp_chksum(bnd_type%u_flux         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%v_flux        ',mpp_chksum(bnd_type%v_flux         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudu        ',mpp_chksum(bnd_type%dtaudu         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudv        ',mpp_chksum(bnd_type%dtaudv         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_star        ',mpp_chksum(bnd_type%u_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%b_star        ',mpp_chksum(bnd_type%b_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%q_star        ',mpp_chksum(bnd_type%q_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%rough_mom     ',mpp_chksum(bnd_type%rough_mom      )
!    write(outunit,100) 'lnd_ice_atm_bnd_type%data          ',mpp_chksum(bnd_type%data           )

end subroutine lnd_ice_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'lnd_atm_bnd_type%data',mpp_chksum(bnd_type%data)

100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine lnd_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call ice_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!


subroutine ice_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(ice_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'ice_atm_bnd_type%data',mpp_chksum(data_type%data)

100 FORMAT("CHECKSUM::",A32," = ",Z20)


end subroutine ice_atm_bnd_type_chksum
! </SUBROUTINE>

end module atmos_model_mod

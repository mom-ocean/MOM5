module land_model_mod ! This is the null version

!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang
!</CONTACT>
!
!<OVERVIEW>
! Null land model. 
!</OVERVIEW>
!<DESCRIPTION>
! Null land model. 
!</DESCRIPTION>
!
!<NAMELIST NAME="land_model_nml">
!  <DATA NAME="layout" TYPE="integer">
!  Processor domain layout for land model. 
!  </DATA> 
!  <DATA NAME="mask_table" TYPE="character">
!   A text file to specify n_mask, layout and mask_list to reduce number of processor
!   usage by masking out some domain regions which contain all land points. 
!   The default file name of mask_table is "INPUT/land_mask_table". Please note that 
!   the file name must begin with "INPUT/". The first 
!   line of mask_table will be number of region to be masked out. The second line 
!   of the mask_table will be the layout of the model. User need to set land_model_nml
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


use  mpp_mod,           only : mpp_pe, mpp_chksum, mpp_root_pe
use  mpp_mod,           only : input_nml_file

use  mpp_domains_mod,   only : CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domain
use  mpp_domains_mod,   only : domain2d, mpp_define_layout, mpp_define_domains
use  mpp_domains_mod,   only : mpp_get_ntile_count, mpp_get_tile_id, mpp_get_current_ntile

use time_manager_mod,   only : time_type

use diag_manager_mod,   only : diag_axis_init

use tracer_manager_mod, only : register_tracers

use field_manager_mod,  only : MODEL_LAND

use          fms_mod,  only : write_version_number, error_mesg, FATAL, mpp_npes, stdout
use          fms_mod,  only : open_namelist_file, check_nml_error, file_exist, close_file
use       fms_io_mod,  only : parse_mask_table

use         grid_mod,  only : get_grid_ntiles, get_grid_size, define_cube_mosaic
use         grid_mod,  only : get_grid_cell_vertices, get_grid_cell_centers
use         grid_mod,  only : get_grid_cell_area, get_grid_comp_area

implicit none
private

! ==== public interfaces =====================================================

public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! saves the land model restart(s)
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler
public :: Lnd_stock_pe          ! return stocks of conservative quantities
public land_data_type_chksum, atm_lnd_bnd_type_chksum

! ==== end of public interfaces ==============================================

character(len=*), parameter :: &
     version = '$Id: land_model.F90,v 20.0 2013/12/13 23:31:21 fms Exp $', &
     tagname = '$Name: tikal $'

type :: atmos_land_boundary_type
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperture of precipitation, degK
        sw_flux_down_vis_dir   => NULL(), & ! visible direct 
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature 
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value, 
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   type(domain2d) :: domain  ! our computation domain
   logical, pointer :: maskmap(:,:) 
   integer, pointer, dimension(:) :: pelist
end type land_data_type

!---- land_model_nml
integer :: layout(2)
! mask_table contains information for masking domain ( n_mask, layout and mask_list).
character(len=128) :: mask_table = "INPUT/land_mask_table"

namelist /land_model_nml/ layout, mask_table


contains

! ============================================================================
subroutine land_model_init (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow, &
                            glon_bnd, glat_bnd, Domain_in)

  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step
  real, dimension(:,:), optional :: glon_bnd, glat_bnd
  type(domain2d), optional :: Domain_in

 ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid
  integer :: is,ie,js,je,id_lon,id_lat,i
  type(domain2d), save :: domain
  real, allocatable, dimension(:,:)  :: garea, gcellarea, gfrac
  real, allocatable, dimension(:,:)  :: glon, glat
  integer, allocatable, dimension(:) :: tile_ids
  integer :: ntracers, ntprog, ndiag, face, npes_per_tile
  integer :: namelist_unit, io, ierr, stdoutunit

!--- read namelist
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, land_model_nml, iostat=io)
#else
   namelist_unit = open_namelist_file()
   ierr=1
   do while (ierr /= 0)
     read(namelist_unit, nml=land_model_nml, iostat=io, end=20)
     ierr = check_nml_error (io, 'land_model_nml')
   enddo
   20 call close_file (namelist_unit)
#endif

  stdoutunit = stdout()

  call write_version_number (version, tagname)

  ! define the processor layout information according to the global grid size 
  call get_grid_ntiles('LND',ntiles)
  call get_grid_size('LND',1,nlon,nlat)
 
  if(file_exist(mask_table)) then
     if(ntiles > 1) then
        call error_mesg('land_model_init', &
          'file '//trim(mask_table)//' should not exist when ntiles is not 1', FATAL)
     endif
     if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('land_model_init', &
          'land_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

     write(stdoutunit, '(a)') '==> NOTE from land_model_init:  reading maskmap information from '//trim(mask_table)
     allocate(land2cplr%maskmap(layout(1), layout(2)))
     call parse_mask_table(mask_table, land2cplr%maskmap, "Land model")
  else
     if(layout(1)*layout(2) .NE. mpp_npes()/ntiles) call mpp_define_layout((/1,nlon,1,nlat/), mpp_npes()/ntiles, layout)
  endif

  if(ntiles ==1) then
     if( ASSOCIATED(land2cplr%maskmap) ) then
        call mpp_define_domains((/1,nlon,1,nlat/), layout, domain, &
             xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, maskmap = land2cplr%maskmap , name='land model')
     else
        call mpp_define_domains((/1,nlon,1,nlat/), layout, domain, &
             xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, name='land model')
     end if
  else
     call define_cube_mosaic('LND', domain, layout, halo=1)
  endif

  land2cplr%domain = domain

  npes_per_tile = mpp_npes()/ntiles
  face = (mpp_pe()-mpp_root_pe())/npes_per_tile + 1
  allocate(garea(nlon,nlat), gcellarea(nlon,nlat), gfrac(nlon,nlat))
  call get_grid_cell_area    ('LND',face,gcellarea)
  call get_grid_comp_area    ('LND',face,garea)
  gfrac = garea/gcellarea
  

  call mpp_get_compute_domain(domain, is,ie,js,je)
  allocate(land2cplr%tile_size(is:ie,js:je,1))
  land2cplr%tile_size(is:ie,js:je,1) = gfrac(is:ie,js:je)
  deallocate(gfrac, garea, gcellarea)

  allocate(tile_ids(mpp_get_current_ntile(domain)))
  tile_ids = mpp_get_tile_id(domain)
  allocate(glon(nlon,nlat), glat(nlon,nlat))
  call get_grid_cell_centers ('LND', tile_ids(1), glon,  glat)

  if(mpp_get_ntile_count(domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define longitude axes and its edges
     id_lon  = diag_axis_init('lon', glon(:,1), 'degrees_E', 'X', &
          long_name='longitude', set_name='land', domain2=domain )

     ! define latitude axes and its edges
     id_lat = diag_axis_init ('lat', glat(1,:), 'degrees_N', 'Y', &
          long_name='latitude',  set_name='land', domain2=domain )
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
          long_name='T-cell longitude', set_name='land', domain2=domain, aux='geolon_t' )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
          long_name='T-cell latitude',  set_name='land', domain2=domain, aux='geolat_t' )
  endif
  land2cplr%axes = (/id_lon,id_lat/)

  call register_tracers(MODEL_LAND, ntracers, ntprog, ndiag)

  allocate(land2cplr%mask          (is:ie,js:je,1))
  allocate(land2cplr%t_surf        (is:ie,js:je,1))
  allocate(land2cplr%t_ca          (is:ie,js:je,1))
  allocate(land2cplr%tr            (is:ie,js:je,1,ntprog))
  allocate(land2cplr%albedo        (is:ie,js:je,1))
  allocate(land2cplr%albedo_vis_dir(is:ie,js:je,1))
  allocate(land2cplr%albedo_nir_dir(is:ie,js:je,1))
  allocate(land2cplr%albedo_vis_dif(is:ie,js:je,1))
  allocate(land2cplr%albedo_nir_dif(is:ie,js:je,1))
  allocate(land2cplr%rough_mom     (is:ie,js:je,1))
  allocate(land2cplr%rough_heat    (is:ie,js:je,1))
  allocate(land2cplr%rough_scale   (is:ie,js:je,1))
  allocate(land2cplr%discharge     (is:ie,js:je))
  allocate(land2cplr%discharge_heat(is:ie,js:je))
  allocate(land2cplr%discharge_snow(is:ie,js:je))
  allocate(land2cplr%discharge_snow_heat(is:ie,js:je))

  land2cplr%mask              = .FALSE.
  land2cplr%t_surf            = 280.0
  land2cplr%t_ca              = 280.0
  land2cplr%tr                = 0.0
  land2cplr%albedo            = 0.0
  land2cplr%albedo_vis_dir    = 0.0
  land2cplr%albedo_nir_dir    = 0.0
  land2cplr%albedo_vis_dif    = 0.0
  land2cplr%albedo_nir_dif    = 0.0
  land2cplr%rough_mom         = 0.0
  land2cplr%rough_heat        = 0.0
  land2cplr%rough_scale       = 1.0
  land2cplr%discharge         = 0.0
  land2cplr%discharge_heat    = 0.0
  land2cplr%discharge_snow    = 0.0
  land2cplr%discharge_snow_heat = 0.0

  allocate(cplr2land%t_flux(is:ie,js:je,1) )
  allocate(cplr2land%lw_flux(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux(is:ie,js:je,1) )
  allocate(cplr2land%lprec(is:ie,js:je,1) )
  allocate(cplr2land%fprec(is:ie,js:je,1) )
  allocate(cplr2land%tprec(is:ie,js:je,1) )
  allocate(cplr2land%dhdt(is:ie,js:je,1) )
  allocate(cplr2land%dhdq(is:ie,js:je,1) )
  allocate(cplr2land%drdt(is:ie,js:je,1) )
  allocate(cplr2land%p_surf(is:ie,js:je,1) )
  allocate(cplr2land%tr_flux(is:ie,js:je,1,ntprog) )
  allocate(cplr2land%dfdtr(is:ie,js:je,1,ntprog) )
  allocate(cplr2land%lwdn_flux(is:ie,js:je,1) )
  allocate(cplr2land%swdn_flux(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_vis_dir(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_total_dir(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_vis_dif(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_total_dif(is:ie,js:je,1) )
  allocate(cplr2land%cd_t(is:ie,js:je,1) )
  allocate(cplr2land%cd_m(is:ie,js:je,1) )
  allocate(cplr2land%bstar(is:ie,js:je,1) )
  allocate(cplr2land%ustar(is:ie,js:je,1) )
  allocate(cplr2land%wind(is:ie,js:je,1) )
  allocate(cplr2land%z_bot(is:ie,js:je,1) )
  allocate(cplr2land%drag_q(is:ie,js:je,1) )

  cplr2land%t_flux                 = 0.0
  cplr2land%lw_flux                = 0.0
  cplr2land%sw_flux                = 0.0
  cplr2land%lprec                  = 0.0
  cplr2land%fprec                  = 0.0
  cplr2land%tprec                  = 0.0
  cplr2land%dhdt                   = 0.0
  cplr2land%dhdq                   = 0.0
  cplr2land%drdt                   = 0.0
  cplr2land%p_surf                 = 1.0e5
  cplr2land%tr_flux                = 0.0
  cplr2land%dfdtr                  = 0.0
  cplr2land%lwdn_flux              = 0.0
  cplr2land%swdn_flux              = 0.0
  cplr2land%sw_flux_down_vis_dir   = 0.0
  cplr2land%sw_flux_down_total_dir = 0.0
  cplr2land%sw_flux_down_vis_dif   = 0.0
  cplr2land%sw_flux_down_total_dif = 0.0
  cplr2land%cd_t                   = 0.0
  cplr2land%cd_m                   = 0.0
  cplr2land%bstar                  = 0.0
  cplr2land%ustar                  = 0.0
  cplr2land%wind                   = 0.0
  cplr2land%z_bot                  = 0.0
  cplr2land%drag_q                 = 0.0

  deallocate(glon, glat, tile_ids)

end subroutine land_model_init

! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr
  
end subroutine land_model_end

! ============================================================================
subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp
  
end subroutine land_model_restart

! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  call error_mesg('update_land_model_fast','Should not be calling null version of update_land_model_fast',FATAL)

end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  call error_mesg('update_land_model_slow','Should not be calling null version of update_land_model_slow',FATAL)

end subroutine update_land_model_slow

! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd 
integer             , intent(in)  :: index
real                , intent(out) :: value ! Domain water (Kg) or heat (Joules)

value = 0.0
 
end subroutine Lnd_stock_pe
! ============================================================================

!#######################################################################
! <SUBROUTINE NAME="atm_lnd_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_land_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_land_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atm_lnd_bnd_type_chksum(id, timestep, albt)
! </TEMPLATE>

! <IN NAME="albt" TYPE="type(atmos_land_boundary_type)">
!   Derived-type variable that contains fields in the atmos_land_boundary_type.
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
subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_land_boundary_type), intent(in) :: albt
 integer ::   n, outunit
 
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(atmos_land_boundary_type):: ', id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%swdn_flux               ', mpp_chksum( albt%swdn_flux)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%dhdq                  ', mpp_chksum( albt%dhdq)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)
    do n = 1,size(albt%tr_flux,4)
    write(outunit,100) 'albt%tr_flux               ', mpp_chksum( albt%tr_flux(:,:,:,n))
    enddo
    do n = 1,size(albt%dfdtr,4)
    write(outunit,100) 'albt%dfdtr                 ', mpp_chksum( albt%dfdtr(:,:,:,n))
    enddo

100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum

! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="land_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call land_data_type_chksum(id, timestep, land)
! </TEMPLATE>

! <IN NAME="land" TYPE="type(land_data_type)">
!   Derived-type variable that contains fields in the land_data_type.
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

subroutine land_data_type_chksum(id, timestep, land)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_data_type), intent(in) :: land
    integer ::   n, outunit
 
    outunit = stdout()

    write(outunit,*) 'BEGIN CHECKSUM(land_data_type):: ', id, timestep
    write(outunit,100) 'land%tile_size         ',mpp_chksum(land%tile_size)
    write(outunit,100) 'land%t_surf            ',mpp_chksum(land%t_surf)
    write(outunit,100) 'land%t_ca              ',mpp_chksum(land%t_ca)
    write(outunit,100) 'land%albedo            ',mpp_chksum(land%albedo)
    write(outunit,100) 'land%albedo_vis_dir    ',mpp_chksum(land%albedo_vis_dir)
    write(outunit,100) 'land%albedo_nir_dir    ',mpp_chksum(land%albedo_nir_dir)
    write(outunit,100) 'land%albedo_vis_dif    ',mpp_chksum(land%albedo_vis_dif)
    write(outunit,100) 'land%albedo_nir_dif    ',mpp_chksum(land%albedo_nir_dif)
    write(outunit,100) 'land%rough_mom         ',mpp_chksum(land%rough_mom)
    write(outunit,100) 'land%rough_heat        ',mpp_chksum(land%rough_heat)
    write(outunit,100) 'land%rough_scale       ',mpp_chksum(land%rough_scale)

    do n = 1, size(land%tr,4)
    write(outunit,100) 'land%tr                ',mpp_chksum(land%tr(:,:,:,n))
    enddo
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)
    write(outunit,100) 'land%discharge_heat    ',mpp_chksum(land%discharge_heat)


100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine land_data_type_chksum

! </SUBROUTINE>
end module land_model_mod

! ============================================================================
! 
! ============================================================================
module land_types_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Christopher Milly
! </CONTACT> 

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Elena Shevliakova
! </REVIEWER>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Sergey Malyshev
! </REVIEWER>

! <OVERVIEW>
!   Initializes land data types and allocates boundary data.
! </OVERVIEW>

! <DESCRIPTION>
!   This module defines two derived types. atmos_land_boundary_type represents
!   data passed from the atmosphere to the surface and the derivatives of the
!   fluxes. Data describing the land is land_data_type. This module contains
!   routines for initializing land type data and allocating and deallocating
!   data for the specified domain and number of tiles.
! </DESCRIPTION>

  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain
  use fms_mod, only : write_version_number, stdout
  use mpp_mod,                 only: mpp_chksum

implicit none
private

! ==== public interfaces =====================================================
public :: land_types_init
public :: atmos_land_boundary_type
public :: land_data_type

public :: allocate_boundary_data
public :: deallocate_boundary_data
public :: land_data_type_chksum, atm_lnd_bnd_type_chksum
! ==== end of public interfaces ==============================================

! <TYPE NAME="atmos_land_boundary_type">
type :: atmos_land_boundary_type
!   <DESCRIPTION>
!     Data passed from the atmosphere to the surface and derivatives
!     of the fluxes. This data is provided by the flux_exchange.
!   </DESCRIPTION>

!   <DATA NAME="t_flux" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Sensible heat flux
!   </DATA>

!   <DATA NAME="q_flux" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!     Water vapor flux
!   </DATA>

!   <DATA NAME="lw_flux" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Net longwave flux
!   </DATA>

!   <DATA NAME="sw_flux" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Net shortwave flux
!   </DATA>

!   <DATA NAME="sw_flux_down_vis_dir" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Downward visible band direct shortwave flux
!   </DATA>

!   <DATA NAME="sw_flux_down_total_dir" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Total direct downward shortwave flux
!   </DATA>

!   <DATA NAME="sw_flux_down_vis_dif" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Downward visible diffuse shortwave flux
!   </DATA>

!   <DATA NAME="sw_flux_down_total_diff" UNITS="W/m2" TYPE="real, pointer" DIM="3">
!     Total downward diffuse shortwave flux
!   </DATA>
 
!   <DATA NAME="lprec" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!     Liquid precipitation
!   </DATA>

!   <DATA NAME="fprec" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!     Solid precipitation (snowfall)
!   </DATA>

!   <DATA NAME="tprec" UNITS="deg K" TYPE="real, pointer" DIM="3">
!     Temperature of precipitation
!   </DATA>

   ! data passed from the atmosphere to the surface
   real, dimension(:,:,:), pointer :: &
        t_flux =>NULL(),  &
        lw_flux =>NULL(), &
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux =>NULL(), &
        sw_flux_down_vis_dir =>NULL(), &
        sw_flux_down_total_dir =>NULL(), &
        sw_flux_down_vis_dif =>NULL(), &
        sw_flux_down_total_dif =>NULL(), &
        lprec =>NULL(),   &
        fprec =>NULL(),   &
        tprec =>NULL()

!   <DATA NAME="dhdt" UNITS="W/m2/K" TYPE="real, pointer" DIM="3">
!     Derivative of sensible heat over land surface temperature
!   </DATA>

!   <DATA NAME="dedt" UNITS="kg/m2/s/K" TYPE="real, pointer" DIM="3">
!     Derivative of evaporation over land surface temperature 
!   </DATA>

!   <DATA NAME="dedq" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!      Derivative of evaporation over near-surface specific humidity
!   </DATA>

!   <DATA NAME="drdt" UNITS="W/m2/K" TYPE="real, pointer" DIM="3">
!     Derivative of LW radiation over land surface temperature
!   </DATA>

   ! derivatives of the fluxes
   real, dimension(:,:,:), pointer :: &
        dhdt =>NULL(),    &
        drdt =>NULL()

!   <DATA NAME="drag_q" UNITS="m/s" TYPE="real, pointer" DIM="3">
!     Product of drag coefficient for water vapor by wind
!   </DATA>

!   <DATA NAME="p_surf" UNITS="N/m2" TYPE="real, pointer" DIM="3">
!     Surface pressure
!   </DATA>

   real, dimension(:,:,:), pointer :: &
        cd_m => NULL(),      &   ! drag coefficient for momentum, dimensionless
        cd_t => NULL(),      &   ! drag coefficient for tracers, dimensionless
        ustar => NULL(),     &   ! turbulent wind scale, m/s
        bstar => NULL(),     &   ! turbulent bouyancy scale, m/s
        wind => NULL(),      &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot => NULL(),     &   ! height of the bottom atmos. layer above surface, m
        drag_q =>NULL(),     &   ! product of cd_q by wind
        p_surf =>NULL()          ! surface pressure

#ifdef LAND_BND_TRACERS
   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile,tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! tendency, including evap over surface specific humidity
#else
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        q_flux => NULL(),    &   ! water vapor flux
        dedt => NULL(),      &   ! evap over surface temperature (assuming saturation @ surf.)
        dedq => NULL()           ! evap over surface specufuc humidity
#endif 


!   <DATA NAME="xtype" TYPE="integer">
!     REGRID, REDIST or DIRECT
!   </DATA>

   integer :: xtype             !REGRID, REDIST or DIRECT

end type atmos_land_boundary_type
! </TYPE>

! <TYPE NAME="land_data_type">
!   <DESCRIPTION>
!     Data describing the land.
!   </DESCRIPTION>
type :: land_data_type

!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     The computational domain
!   </DATA>

   type(domain2d) :: domain  ! our computation domain

!   <DATA NAME="tile_size" TYPE="real, pointer" DIM="3">
!     Fractional coverage of cell by tile, dimensionless
!   </DATA>

!   <DATA NAME="t_surf" UNITS="K" TYPE="real, pointer" DIM="3">
!     Ground surface temperature
!   </DATA>

!   <DATA NAME="t_ca" UNITS="K" TYPE="real, pointer" DIM="3">
!     Canopy air temperature
!   </DATA>

!   <DATA NAME="q_ca" UNITS="kg/kg" TYPE="real, pointer" DIM="3">
!     Canopy air specific humidity
!   </DATA>

!   <DATA NAME="albedo" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land albedo
!   </DATA>

!   <DATA NAME="albedo_vis_dir" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land albedo for direct visible
!   </DATA>
 
!   <DATA NAME="albedo_nir_dir" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land albedo for direct nir
!   </DATA>

!   <DATA NAME="albedo_vis_dif" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land albedo for diffuse visible
!   </DATA>

!   <DATA NAME="albedo_nir_dif" TYPE="real, pointer" DIM="3">
!     Snow-adjusted land albedo for diffuse nir
!   </DATA>

!   <DATA NAME="rough_mom" UNITS="m" TYPE="real, pointer" DIM="3">
!     Momentum roughness length
!   </DATA>

!   <DATA NAME="rough_heat" UNITS="m" TYPE="real, pointer" DIM="3">
!     Roughness length for tracers (heat and water)
!   </DATA>

!   <DATA NAME="rough_scale" UNITS="m" TYPE="real, pointer" DIM="3">
!     Roughness length for topographic momentum drag coefficient scaling
!   </DATA>
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size =>NULL(),       & ! fractional coverage of cell by tile, dimensionless
        t_surf =>NULL(),          & ! ground surface temperature, degK
        t_ca =>NULL(),            & ! canopy air temperature, degK
        albedo =>NULL(),          & ! snow-adjusted land albedo
        albedo_vis_dir =>NULL(),  & ! albedo for direct visible-band radiation
        albedo_nir_dir =>NULL(),  & ! albedo for direct nir-band radiation 
        albedo_vis_dif =>NULL(),  & ! albedo for diffuse visible-band radiation 
        albedo_nir_dif =>NULL(),  & ! albedo for diffuse nir-band radiation
        rough_mom =>NULL(),       & ! momentum roughness length, m
        rough_heat =>NULL(),      & ! roughness length for tracers (heat and water), m
        rough_scale =>NULL()        ! roughness length for drag scaling

#ifdef LAND_BND_TRACERS
   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr   => NULL()               ! tracers, including canopy air specific humidity, kg/kg
#else
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        q_ca => NULL()               ! canopy air specific humidity, kg/kg
#endif


!   <DATA NAME="discharge" UNITS="kg/m2/s" TYPE="real, pointer" DIM="2">
!     Outflow of fresh water from river mouths into the ocean (per unit area
!     of the ocean part of the grid cell)
!   </DATA>

!   <DATA NAME="discharge_heat" UNITS="W/m2" TYPE="real, pointer" DIM="2">
!     Sensible heat of discharge  (0 C datum)
!   </DATA>

!   <DATA NAME="discharge_snow" UNITS="kg/m2/s" TYPE="real, pointer" DIM="2">
!     Snow analogue of discharge
!   </DATA>

!   <DATA NAME="discharge_snow_heat" UNITS="W/m2" TYPE="real, pointer" DIM="2">
!     Sensible heat of discharge_snow  (0 C datum)
!   </DATA>

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
        discharge           => NULL(),  & ! flux from surface drainage network out of land model
        discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
        discharge_snow      => NULL(),  & ! snow analogue of discharge
        discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

!   <DATA NAME="mask" TYPE="logical, pointer" DIM="3">
!     Land mask; true if land
!   </DATA>

   logical, pointer, dimension(:,:,:):: &
        mask =>NULL()        ! true if land

!   <DATA NAME="maskmap" TYPE="logical, pointer" DIM="2">
!     A pointer to an array indicating which logical processors are actually used for
!     the ocean code. The other logical processors would be all land points and
!     are not assigned to actual processors. This need not be assigned if all logical
!     processors are used. This variable is dummy and need not to be set, 
!     but it is needed to pass compilation.
!   </DATA>

   logical, pointer, dimension(:,:) :: maskmap =>NULL()  ! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. 

!   <DATA NAME="axes(2)" TYPE="integer" DIM="2">
!     Axes IDs for diagnostics
!   </DATA>

   integer :: axes(2)      ! axes IDs for diagnostics  
   logical :: pe

   integer, pointer, dimension(:) :: pelist

end type land_data_type
! </TYPE>

logical :: module_is_initialized =.FALSE.
character(len=128) :: version = '$Id: land_types.F90,v 20.0 2013/12/13 23:28:40 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================

! <SUBROUTINE NAME="land_types_init">

!   <OVERVIEW>
!     Initializes the land types.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Initializes the land types module and writes the CVS version
!     and tagname
!   </DESCRIPTION>

!   <TEMPLATE>
!     call land_types_init()
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="">
subroutine land_types_init()
!   </PUBLICROUTINE>

  module_is_initialized =.TRUE. 
  call write_version_number(version, tagname)

end subroutine land_types_init
! </SUBROUTINE>


! <SUBROUTINE NAME="allocate_boundary_data">

!   <OVERVIEW>
!     Allocates the boundary data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Allocates data for the specified domain and number of tiles.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call allocate_boundary_data (bnd, domain, n_tiles)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine allocate_boundary_data (a2l, bnd, domain, n_tiles, n_tracers)
  type(atmos_land_boundary_type), intent(inout) :: a2l
  type(land_data_type), intent(inout) :: bnd     ! data to allocate
  type(domain2d),       intent(in)  :: domain  ! domain to allocate for
  integer,              intent(in)  :: n_tiles ! number of tiles
  integer,              intent(in)  :: n_tracers ! number of tracers
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: is,ie,js,je ! boundaries of the domain

  ! get the size of our computation domain
  call mpp_get_compute_domain ( domain, is,ie,js,je )

  ! allocate data according to the domain boundaries
  allocate ( &
       bnd % tile_size  (is:ie,js:je,n_tiles), & 
       bnd % t_surf     (is:ie,js:je,n_tiles), &
       bnd % t_ca       (is:ie,js:je,n_tiles), &
       bnd % albedo     (is:ie,js:je,n_tiles), & 
       bnd % albedo_vis_dir (is:ie,js:je,n_tiles), &
       bnd % albedo_nir_dir (is:ie,js:je,n_tiles), &
       bnd % albedo_vis_dif (is:ie,js:je,n_tiles), &
       bnd % albedo_nir_dif (is:ie,js:je,n_tiles), &
       bnd % rough_mom  (is:ie,js:je,n_tiles), & 
       bnd % rough_heat (is:ie,js:je,n_tiles), & 
       bnd % rough_scale(is:ie,js:je,n_tiles), & 

       bnd % discharge             (is:ie,js:je),   & 
       bnd % discharge_heat        (is:ie,js:je),   & 
       bnd % discharge_snow        (is:ie,js:je),   &
       bnd % discharge_snow_heat   (is:ie,js:je),   &

       bnd % mask       (is:ie,js:je,n_tiles) &
       )
#ifdef LAND_BND_TRACERS
  allocate( bnd % tr  (is:ie,js:je,n_tiles,n_tracers) )
#else
  allocate( bnd % q_ca(is:ie,js:je,n_tiles) )
#endif  

  allocate( a2l % t_flux  (is:ie,js:je,n_tiles) )
  allocate( a2l % lw_flux (is:ie,js:je,n_tiles) )
  allocate( a2l % sw_flux (is:ie,js:je,n_tiles) )
  allocate( a2l % lprec   (is:ie,js:je,n_tiles) )
  allocate( a2l % fprec   (is:ie,js:je,n_tiles) )
  allocate( a2l % tprec   (is:ie,js:je,n_tiles) )
  allocate( a2l % dhdt    (is:ie,js:je,n_tiles) )
  allocate( a2l % drdt    (is:ie,js:je,n_tiles) )
  allocate( a2l % drag_q  (is:ie,js:je,n_tiles) )
  allocate( a2l % p_surf  (is:ie,js:je,n_tiles) )
  allocate( a2l % sw_flux_down_vis_dir   (is:ie,js:je,n_tiles) )
  allocate( a2l % sw_flux_down_total_dir (is:ie,js:je,n_tiles) )
  allocate( a2l % sw_flux_down_vis_dif   (is:ie,js:je,n_tiles) )
  allocate( a2l % sw_flux_down_total_dif (is:ie,js:je,n_tiles) )

#ifdef LAND_BND_TRACERS
  allocate( a2l % tr_flux(is:ie,js:je,n_tiles,n_tracers) )
  allocate( a2l % dfdtr(is:ie,js:je,n_tiles,n_tracers) )
#else
  allocate( a2l % q_flux(is:ie,js:je,n_tiles) )
  allocate( a2l % dedt(is:ie,js:je,n_tiles) )
  allocate( a2l % dedq(is:ie,js:je,n_tiles) )
#endif

! set up initial values (discharge_heat and discharge_snow_heat are set to zero and never change)
  bnd % discharge_heat      = 0.0
  bnd % discharge_snow_heat = 0.0

end subroutine allocate_boundary_data
! </SUBROUTINE>


! <SUBROUTINE NAME="deallocate_boundary_data">

!   <OVERVIEW>
!     Deallocates the boundary data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Deallocates the boundary data.
!     
!   </DESCRIPTION>

!   <TEMPLATE>
!     call deallocate_boundary_data ( bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine deallocate_boundary_data ( a2l, bnd )
  type(atmos_land_boundary_type), intent(inout) :: a2l
  type(land_data_type), intent(inout) :: bnd  ! data to deallocate
!   </PUBLICROUTINE>

  deallocate ( &
       bnd % tile_size  , & 
       bnd % t_surf     , &
       bnd % t_ca       , &
       bnd % albedo     , & 
       bnd % albedo_vis_dir , &
       bnd % albedo_nir_dir , &
       bnd % albedo_vis_dif , &
       bnd % albedo_nir_dif , &
       bnd % rough_mom  , & 
       bnd % rough_heat , & 
       bnd % rough_scale, &

       bnd % discharge      ,      & 
       bnd % discharge_heat ,      & 
       bnd % discharge_snow ,      &
       bnd % discharge_snow_heat , &

       bnd % mask        &
       )
#ifdef LAND_BND_TRACERS
  if(associated( bnd%tr ))             deallocate( bnd%tr )
#else
  if(associated( bnd%q_ca ))           deallocate( bnd%q_ca )
#endif

  deallocate( a2l % t_flux   )
  deallocate( a2l % lw_flux  )
  deallocate( a2l % sw_flux  )
  deallocate( a2l % lprec    )
  deallocate( a2l % fprec    )
  deallocate( a2l % tprec    )
  deallocate( a2l % dhdt     )
  deallocate( a2l % drdt     )
  deallocate( a2l % drag_q   )
  deallocate( a2l % p_surf   )
  deallocate( a2l % sw_flux_down_vis_dir    )
  deallocate( a2l % sw_flux_down_total_dir  )
  deallocate( a2l % sw_flux_down_vis_dif    )
  deallocate( a2l % sw_flux_down_total_dif  )

#ifdef LAND_BND_TRACERS
  deallocate( a2l % tr_flux )
  deallocate( a2l % dfdtr )
#else
  deallocate( a2l % q_flux )
  deallocate( a2l % dedt )
  deallocate( a2l % dedq )
#endif

end subroutine deallocate_boundary_data
! </SUBROUTINE>

subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_land_boundary_type), intent(in) :: albt
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) "BEGIN CHECKSUM(atmos_land_boundary_type):: ", id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    if ( associated(albt%lwdn_flux) ) write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    if ( associated(albt%cd_m) )  write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    if ( associated(albt%cd_t) )  write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    if ( associated(albt%ustar) ) write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    if ( associated(albt%bstar) ) write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    if ( associated(albt%wind) )  write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    if ( associated(albt%z_bot) ) write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)
#ifdef LAND_BND_TRACERS
    do n = 1,size(albt%tr_flux,4)
    write(outunit,100) 'albt%tr_flux               ', mpp_chksum( albt%tr_flux(:,:,:,n))
    enddo
    do n = 1,size(albt%dfdtr,4)
    write(outunit,100) 'albt%dfdtr                 ', mpp_chksum( albt%dfdtr(:,:,:,n))
    enddo
#else
    write(outunit,100) 'albt%q_flux                ', mpp_chksum( albt%q_flux)
    write(outunit,100) 'albt%dedt                  ', mpp_chksum( albt%dedt)
    write(outunit,100) 'albt%dedq                  ', mpp_chksum( albt%dedq)
#endif 

100 FORMAT("   CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum

subroutine land_data_type_chksum(id, timestep, land)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_data_type), intent(in) :: land
    integer ::   n, outunit
    
    outunit = stdout()

    write(outunit,*) "BEGIN CHECKSUM(land_type):: ", id, timestep
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

#ifdef LAND_BND_TRACERS
    do n = 1, size(land%tr,4)
    write(outunit,100) 'land%tr                ',mpp_chksum(land%tr(:,:,:,n))
    enddo
#else
    write(outunit,100) 'land%q_ca              ',mpp_chksum(land%q_ca)
#endif
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)


100 FORMAT("   CHECKSUM::",A32," = ",Z20)
end subroutine land_data_type_chksum


end module land_types_mod

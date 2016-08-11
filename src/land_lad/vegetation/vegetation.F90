! ============================================================================
! vegetation-related processes
! ============================================================================
module vegetation_mod

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
!   Module containing processes relating to vegetation.
! </OVERVIEW>

! <DESCRIPTION>
!   Vegetation type describing vegetation characteristics is defined.
!   Vegetation data is allocated and the initial value for canopy air
!   humidity is set up. Updates vegetation and boundary data on the slow
!   and fast time-scales. Deallocates vegetation data, empties memory and
!   cleans up, if necessary.
! </DESCRIPTION>

  ! things from the other modules we use in the interface part of this module
  use time_manager_mod, only: time_type, get_time
  use mpp_domains_mod,  only: domain2d, mpp_get_compute_domain
  use land_types_mod,   only: land_data_type
  use soil_mod,         only: soil_type
  use constants_mod,    only: rdgas, rvgas
  use fms_mod,          only: write_version_number, error_mesg, FATAL,      &
                              file_exist, open_restart_file, close_file,    &
                              read_data, write_data, set_domain, mpp_pe,    &
                              open_namelist_file, check_nml_error, mpp_pe,  &
                              NOTE, mpp_root_pe, mpp_error, stdlog
  use fms_io_mod,       only: get_restart_io_mode
  use sat_vapor_pres_mod, only: escomp

  use field_manager_mod,  only: MODEL_LAND
  use tracer_manager_mod, only: get_tracer_index, NO_TRACER

implicit none
private

! ==== public interfaces =====================================================
public :: vegetation_type

public :: vegetation_init              ! initialize vegetation data

public :: vegetation_end               ! finish using vegetation data

public :: update_vegetation_slow       ! slow time-scale vegetation update
public :: update_vegetation_fast_up    ! fast time-scale update of veg. state
public :: update_vegetation_fast_down  ! fast time-scale update of veg. state
public :: update_vegetation_bnd_fast
public :: update_vegetation_bnd_slow

public :: vegetation_stock_pe          ! calculate and return total amount of
                                       ! requested quantitiy per PE
! ==== end of public interfaces ==============================================

! <TYPE NAME="vegetation_type">
type vegetation_type
   private

!   <DESCRIPTION>
!      Private data type describing vegetation characteristics.
!   </DESCRIPTION>

   type(domain2d) :: domain   ! computational domain
   real           :: dt       ! fast time step, s

!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     Computational domain
!   </DATA>
!   <DATA NAME="dt" UNITS="s" TYPE="real">
!     Fast time step
!   </DATA>

   integer :: is,ie,js,je  ! computational domain bounds
   integer :: n_tiles      ! number of tiles
!   <DATA NAME="is" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="ie" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="js" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="je" TYPE="integer">
!     Computational domain bounds
!   </DATA>
!   <DATA NAME="n_tiles" TYPE="integer">
!     Number of tiles
!   </DATA>

   logical, pointer, dimension(:,:,:) :: &
        mask =>NULL(),    & ! land mask
        bonedry =>NULL()    ! true if the "bone dry" conditions occur (evaporation
                          ! during time step larger than available water)
!   <DATA NAME="mask" TYPE="logical, pointer" DIM="3">
!     Land mask
!   </DATA>
!   <DATA NAME="bonedry" TYPE="logical, pointer" DIM="3">
!     True if the "bone dry" conditions occur (evaporation during time step
!     is larger than the available water)
!   </DATA>

   real, pointer, dimension(:,:,:) :: &
        q_ca =>NULL()       ! specific humidity of canopy air
!   <DATA NAME="q_ca" UNITS="kg/kg" TYPE="real, pointer" DIM="3">
!     Specific humidity of canopy air
!   </DATA>

   real, pointer, dimension(:,:,:) :: &
        t_surf =>NULL(),  & ! soil surface temperature, degK
        evap =>NULL(),    & ! explicit estimate of the water vapor flux
        c_surf =>NULL(),  & ! conductance between surface and canopy air
        dqsatdt =>NULL(), & ! derivative of sat. humidity over T surface
        e_q =>NULL(),     & ! implicit scheme coefficient
        f_q =>NULL(),     & ! implicit scheme coefficient
        beta =>NULL()       ! water availability for evaporation
!   <DATA NAME="t_surf" UNITS="K" TYPE="real, pointer" DIM="3">
!     Soil surface temperature
!   </DATA>
!   <DATA NAME="evap" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!     Explicit estimate of the water vapor flux
!   </DATA>
!   <DATA NAME="c_surf" UNITS="kg/m2/s" TYPE="real, pointer" DIM="3">
!     Conductance between surface and canopy air
!   </DATA>
!   <DATA NAME="dqsatdt" UNITS="kg/kg/degK" TYPE="real, pointer" DIM="3">
!     Derivative of sat. humidity over T surface
!   </DATA>
!   <DATA NAME="e_q" TYPE="real, pointer" DIM="3">
!     Implicit scheme coefficient
!   </DATA>
!   <DATA NAME="f_q" TYPE="real, pointer" DIM="3">
!     Implicit scheme coefficient
!   </DATA>
!   <DATA NAME="beta" TYPE="real, pointer" DIM="3">
!     Water availability for evaporation
!   </DATA>

end type vegetation_type
! </TYPE>

! some names, for information only
logical :: module_is_initialized =.FALSE.
character(len=*), private, parameter :: module_name = 'vegetation_mod'
character(len=128), private, parameter :: version     = '$Id: vegetation.F90,v 15.0 2007/08/14 04:00:20 fms Exp $'
character(len=128), private, parameter :: tagname        = '$Name: tikal $'

! module constants
real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622
real, parameter :: d608 = d378/d622

! <NAMELIST NAME="vegetation_nml">
!   <DATA NAME="klev" TYPE="integer" DEFAULT="0">
!     Soil level at which to specify frozen-soil factor for modifying beta.
!   </DATA>
!   <DATA NAME="do_netcdf_restart" TYPE="logical" DEFAULT=".true.">
!     Do netcdf restart.
!   </DATA>
! ---- namelist variables and their default values ---------------------------
integer :: klev    = 0 
logical :: do_netcdf_restart = .true.

namelist /vegetation_nml/ klev, do_netcdf_restart
! </NAMELIST>

! for diagnostics only
integer :: i 
integer, parameter :: iwatch = 38, jwatch=60


! module variables
integer :: isphum ! index of specific humidity in the tracer table

contains


! <SUBROUTINE NAME="vegetation_init">

!   <OVERVIEW>
!     Initializes vegetation data.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Allocates vegetation data and sets up initial value for canopy air
!     humidity.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call vegetation_init &
!     ( veg, gblon, gblat, garea, gfrac, time, dt_fast, dt_slow, domain, &
!     frac, mask, id_lon, id_lat )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine vegetation_init &
     ( veg, gblon, gblat, garea, gfrac, time, dt_fast, dt_slow, domain, &
     frac, mask, id_lon, id_lat )

  type(vegetation_type),intent(inout) :: veg        ! state of a particular
                                                    ! realization to initialize
  real,                 intent(in)    :: gblon(:,:) ! longitude corners of the
                                                    ! grid cells
  real,                 intent(in)    :: gblat(:,:) ! latitude corners of the
                                                    ! grid cells
  real,                 intent(in)    :: garea(:,:) ! grid cell area 
  real,                 intent(in)    :: gfrac(:,:) ! fraction of grid cell
                                                    ! covered by land 
  type(time_type),      intent(in)    :: time       ! current time
  type(time_type),      intent(in)    :: dt_fast    ! fast time step
  type(time_type),      intent(in)    :: dt_slow    ! slow time step
  type(domain2d),       intent(in)    :: domain     ! our domain
  real,                 intent(in)    :: frac(:,:,:)! fractional area of tiles
  logical,              intent(in)    :: mask(:,:,:)! land mask
  integer,              intent(in)    :: id_lon     ! ID of X (longitude) diag
                                                    ! axis
  integer,              intent(in)    :: id_lat     ! ID of Y (latitude) diag
                                                    ! axis
!   </PUBLICROUTINE>

  ! ---- local vars ---------------------------------------------------------
  integer :: unit, ierr, io     ! restart file unit
  integer :: sec, day ! components of time

  if ( file_exist( 'input.nml' ) ) then
     unit = open_namelist_file ( )
     ierr = 1
     do while ( ierr /= 0 )
        read ( unit,  nml = vegetation_nml, iostat = io, end = 10 )
        ierr = check_nml_error ( io, 'vegetation_nml' )
     enddo
10   continue
     call close_file (unit)
  endif
  call get_restart_io_mode(do_netcdf_restart)

! write version and tag information to logfile
  call write_version_number(version, tagname) 
  !  write the namelist to a log file
  if( mpp_pe()==0 ) then
     unit = stdlog( )
     write (unit, nml=vegetation_nml)
     call close_file (unit)
  endif
  ! copy specified domain to our data
  veg % domain = domain

  ! get the size of our domain
  call mpp_get_compute_domain ( veg%domain, veg%is, veg%ie, veg%js, veg%je )
  veg%n_tiles = size(frac, 3)

  ! setup time-related data
  call get_time(dt_fast, sec, day); veg%dt = day*86400.0+sec

  ! allocate data
  allocate ( &
       veg%mask    (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%q_ca    (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%evap    (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%t_surf  (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%c_surf  (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%dqsatdt (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%e_q     (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%f_q     (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%beta    (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles), &
       veg%bonedry (veg%is:veg%ie, veg%js:veg%je, veg%n_tiles)  )

  veg%mask = mask

  ! set up initial value for canopy air humidity
  veg%q_ca = 0
  call set_domain ( veg%domain )

  if (file_exist('INPUT/vegetation.res.nc') .or. file_exist('INPUT/vegetation.res.tile1.nc') ) then
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('vegetation_mod', &
          'Reading NetCDF formatted restart file: INPUT/vegetation.res.nc', NOTE)
     call read_data( 'INPUT/vegetation.res.nc', 'q_ca', veg%q_ca)
  else
     if (file_exist('INPUT/vegetation.res')) then
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('vegetation_mod', &
             'Reading native formatted restart file.', NOTE)
        unit = open_restart_file ( 'INPUT/vegetation.res', 'read')
        call read_data  ( unit, veg%q_ca )
        call close_file ( unit )
     endif
  endif 

! initialize tracers
#ifdef LAND_BND_TRACERS
  isphum = get_tracer_index ( MODEL_LAND, 'sphum' )
  if (isphum==NO_TRACER) then
     call error_mesg('vegetation_init','no required "sphum" tracer',FATAL)
  endif
#else
  isphum = NO_TRACER
#endif

  module_is_initialized =.TRUE.

end subroutine vegetation_init
! </SUBROUTINE>


! <SUBROUTINE NAME="vegetation_end">

!   <OVERVIEW>
!      Deallocates vegetation data; empty memory and do clean-up, if
!      necessary.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Deallocates vegetation data; empty memory and do clean-up, if
!      necessary.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call vegetation_end ( veg )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine vegetation_end ( veg )

  type(vegetation_type), intent(inout) :: veg ! data to finish using
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: unit ! restart file unit 

  ! save restart file
  call set_domain ( veg%domain )
  if( do_netcdf_restart) then
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('vegetation_mod', &
          'Writing NetCDF formatted restart file: RESTART/vegetation.res.nc', NOTE)
     call write_data('RESTART/vegetation.res.nc', 'q_ca', veg%q_ca)
  else
     if (mpp_pe() == mpp_root_pe()) call mpp_error ('vegetation_mod', &
          'Writing native formatted restart file.', NOTE)
     unit = open_restart_file ( 'RESTART/vegetation.res', 'write' )
     call write_data ( unit, veg%q_ca )
     call close_file ( unit )
  endif

  ! deallocate data
  deallocate (        &
       veg%mask,      &
       veg%q_ca,      &
       veg%evap,      &
       veg%t_surf,    &
       veg%c_surf,    &
       veg%dqsatdt,   &
       veg%e_q,       &
       veg%f_q,       &
       veg%beta,      &
       veg%bonedry    )

    module_is_initialized =.FALSE.

end subroutine vegetation_end
! </SUBROUTINE>


! <SUBROUTINE NAME="update_vegetation_slow">

!   <OVERVIEW>
!     Slow time-scale vegetation update.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Slow time-scale vegetation update.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_vegetation_slow ( veg )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_vegetation_slow ( veg )

  type(vegetation_type), intent(inout) :: veg ! data to update
!   </PUBLICROUTINE>
  
  ! call diagnostics_slow ( veg )

end subroutine update_vegetation_slow
! </SUBROUTINE>


! <SUBROUTINE NAME="update_vegetation_fast_down">

!   <OVERVIEW>
!     Fast time-scale vegetation update given soil data inputs.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Fast time-scale vegetation update given soil data inputs. Calculates
!     water availability for evapotranspiration. Calculates saturated specific
!     humidity at the surface and its derivative over the surface temperature.
!     Air density is calculated here using surface q and T; in principle we
!     should use canopy q and T, but surface values are the only ones we have
!     available in this particular implementation. 
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_vegetation_fast_down( veg, soil, evap, dedq, drag_q, psurf, evap1, dedt)
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_vegetation_fast_down( veg, soil, evap, dedq, drag_q, psurf, evap1, dedt)

  type(vegetation_type), intent(inout) :: veg            ! data to update
  type(soil_type),       intent(in)    :: soil           ! soil data inputs
  real, intent(in)  :: &
       evap   (veg%is:veg%ie,veg%js:veg%je,veg%n_tiles),&! evaporation from the
                                                         ! surface into the atm
       drag_q (:,:,:),&                                  ! drag coefficient
       dedq   (:,:,:),&                                  ! derivative of evap
                                                         ! over q_ca
       psurf  (veg%is:veg%ie,veg%js:veg%je,veg%n_tiles)  ! surface pressure

  real, intent(out) :: &
       dedt   (:,:,:), &                                 ! derivative of evap
                                                         ! over T
       evap1  (veg%is:veg%ie,veg%js:veg%je,veg%n_tiles)  ! evaporation from
                                                         ! stomatal into sfc
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  real, dimension(veg%is:veg%ie,veg%js:veg%je,veg%n_tiles) :: &
       qsat,    &
       qsat1,   &
       rho,     &     ! canopy air density 
       fusion_factor  ! frozen-soil factor for modifying beta


  real, parameter :: del_temp = 0.1 ! temperature increment for q_sat derivative calc.

  ! save soil surface temperature for future use in update_vegetation_fast_up
  where (veg%mask)
     veg%t_surf = soil%temp(:,:,:,1)
  elsewhere
     veg%t_surf = 200.0 ! safeguards to avoid errors in escomp
  endwhere

  if (klev == 0) then
      fusion_factor= 1.
  else
      where (soil%mask)
        fusion_factor= 1. - (soil%fusion(:,:,:,klev)/soil%max_fusion(klev))
      endwhere
  endif

  ! calculate water availability for evapotranspiration
  if (soil%conserve_glacier_mass) then
      where(soil%mask .and. (.not.soil%snow>0))
        where (soil%glacier)
          veg%beta = 0.0
        elsewhere
          veg%beta = max(0.0, min(1.0, soil%water/(0.75*soil%max_water)))*  &
             fusion_factor 
        endwhere
      elsewhere
        veg%beta = 1
      endwhere
    else
      where(soil%mask .and. (.not.soil%glacier) .and. (.not.soil%snow>0))
         veg%beta = max(0.0, min(1.0, soil%water/(0.75*soil%max_water)))*  &
            fusion_factor
      elsewhere
         veg%beta = 1
      endwhere
    endif


  ! calculate saturated specific humidity at the surface and its derivative
  ! over T sfc.
  call escomp(veg%t_surf,         qsat ) ! calculate sat water vapor pressure
  call escomp(veg%t_surf+del_temp,qsat1) ! calculate sat water vapor pressure

  evap1 = 0  ! + slm, Mar 29 2002, for some reason if I do not do it, model crashes in bone_dry 
  where (soil%mask)
     ! conversion is placed here because psurf does not have to be defined 
     ! outside of land mask
     qsat   = d622*qsat /(psurf-d378*qsat ) ! convert pres. to spec. humidity
     qsat1  = d622*qsat1/(psurf-d378*qsat1) ! convert pres. to spec. humidity

     ! air density is calculated here using surface q and T; in principle
     ! we should use canopy q and T, but surface values are the only ones we
     ! have available in this particular implementation, and that should do
     rho    = psurf / (rdgas* veg%t_surf  * (1.0 + d608*veg%q_ca))

     veg%dqsatdt = (qsat1 - qsat)/del_temp
     where (soil%stomatal > 0.0)
        veg%c_surf = rho * veg%beta/(soil%stomatal+(1-veg%beta)/drag_q)
        veg%evap      = veg%c_surf * (qsat - veg%q_ca)
        evap1         = veg%evap + veg%c_surf*(evap - veg%evap)/(veg%c_surf + dedq)
        dedt          = veg%dqsatdt*veg%c_surf*(1.0-veg%c_surf/(veg%c_surf + dedq))
     elsewhere
        veg%e_q       = 1 - dedq / ( rho*drag_q )
        veg%f_q       = evap * veg%e_q / ( rho*drag_q*( 1 - veg%e_q ) ) 
        evap1         = evap + dedq * (1-veg%beta) * veg%f_q / (1-(1-veg%beta)*veg%e_q)
        dedt          =        dedq * veg%beta * veg%dqsatdt / (1-(1-veg%beta)*veg%e_q)
     endwhere

  endwhere
! pcm fixed for glacier-mass conservation. however, this will not catch cases
!     where explicit evap1 is ok but implicit evap_new is excessive. i think it
!     might be cleaner to delete this section and instead allow negative stores of
!     snow and/or soil water, as slm once suggested
  veg%bonedry = .false.
  if (soil%conserve_glacier_mass) then
     where (soil%mask)
        where ( evap1*veg%dt > soil%water+soil%snow )
           veg%bonedry = .true.
           evap1      = (soil%water+soil%snow)/veg%dt
           veg%evap   = evap1
           dedt       = 0
        endwhere
     endwhere
  else
     where (soil%mask)
        where ( (.not.soil%glacier) .and. (evap1*veg%dt > soil%water+soil%snow ) )
           veg%bonedry = .true.
           evap1      = (soil%water+soil%snow)/veg%dt
           veg%evap   = evap1
           dedt       = 0
        endwhere
     endwhere
  endif

 
  ! diagnostic output
  if (veg%is<=iwatch.and.iwatch<=veg%ie.and.veg%js<=jwatch.and.jwatch<=veg%je) then
!!$     do i = 1, size(veg%mask,3)
!!$        write(*,'(i2,2L2,100g14.4)')  i, &
!!$             veg%mask     (iwatch,jwatch,i), &
!!$             soil%glacier (iwatch,jwatch,i), &
!!$             soil%stomatal(iwatch,jwatch,i), &
!!$             soil%water   (iwatch,jwatch,i), &
!!$             veg%beta     (iwatch,jwatch,i), &
!!$             veg%t_surf   (iwatch,jwatch,i), &
!!$             veg%q_ca     (iwatch,jwatch,i), &
!!$             qsat         (iwatch,jwatch,i), &
!!$             evap         (iwatch,jwatch,i), &
!!$             evap1        (iwatch,jwatch,i)
!!$             dedt         (iwatch,jwatch,i), &
!!$             psurf        (iwatch,jwatch,i), &
!!$             veg%dqsatdt  (iwatch,jwatch,i)
!!$     enddo
  endif

end subroutine update_vegetation_fast_down
! </SUBROUTINE>


! <SUBROUTINE NAME="update_vegetation_fast_up">

!   <OVERVIEW>
!     Fast time-scale vegetation update given soil data inputs.
!   </OVERVIEW>

!   <DESCRIPTION>
!     Fast time-scale vegetation update given soil data inputs.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_vegetation_fast_up( veg, soil, drag_q, evap, dedq )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_vegetation_fast_up( veg, soil, drag_q, evap, dedq )

  type(vegetation_type), intent(inout) :: veg  ! data to update
  type(soil_type),       intent(in)    :: soil ! soil data inputs
  real, intent(in)  :: &
       drag_q (:,:,:),  &  ! drag coefficient for atmosphere (above vegetation)
       evap   (:,:,:),  &  ! evaporation from surface into the atmosphere
       dedq   (:,:,:)      ! derivative of evap over q_ca
!   </PUBLICROUTINE>

  ! -----local variables-------------------------------------------------------
  real, dimension(veg%is:veg%ie,veg%js:veg%je,veg%n_tiles) :: &
       delta_q_ca, &    ! change in the surface humidity
       delta_t_surf     ! change in surface temperature

  where (veg%mask)
     delta_t_surf  = soil%temp(:,:,:,1) - veg%t_surf
     where ( veg%bonedry )
        delta_q_ca = (veg%evap-evap)/dedq
     elsewhere
        where ( soil%stomatal > 0.0 )
           delta_q_ca = (veg%evap - evap + veg%c_surf*veg%dqsatdt*delta_t_surf) &
                /(veg%c_surf+dedq)
        elsewhere     ! bare soil or glacier (beta=1)
           delta_q_ca = (veg%beta*veg%dqsatdt*delta_t_surf+(1-veg%beta)*veg%f_q) &
                /(1-(1-veg%beta)*veg%e_q)
        endwhere
     endwhere
     veg%q_ca = veg%q_ca + delta_q_ca
 endwhere
 
end subroutine update_vegetation_fast_up
! </SUBROUTINE>


! <SUBROUTINE NAME="update_vegetation_bnd_fast">

!   <OVERVIEW>
!      Updates vegetation boundary data on the fast time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Updates vegetation boundary data on the fast time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_vegetation_bnd_fast ( veg, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE INTERFACE="">
subroutine update_vegetation_bnd_fast ( veg, bnd )

  type(vegetation_type), intent(in)    :: veg ! vegetation data
  type(land_data_type),  intent(inout) :: bnd ! boundary data
!   </PUBLICROUTINE>

#ifdef LAND_BND_TRACERS
  bnd%tr(:,:,:,isphum) = veg%q_ca
#else
  bnd%q_ca = veg%q_ca
#endif

end subroutine update_vegetation_bnd_fast
! </SUBROUTINE>


! <SUBROUTINE NAME="update_vegetation_bnd_slow">

!   <OVERVIEW>
!      Updates vegetation boundary data on the slow time-scale.
!   </OVERVIEW>

!   <DESCRIPTION>
!      Updates vegetation boundary data on the slow time-scale.
!   </DESCRIPTION>

!   <TEMPLATE>
!     call update_vegetation_bnd_slow ( veg, bnd )
!   </TEMPLATE>

!   <PUBLICROUTINE>
subroutine update_vegetation_bnd_slow ( veg, bnd )

  type(vegetation_type), intent(in)    :: veg ! vegetation data
  type(land_data_type),  intent(inout) :: bnd ! boundary data
!   </PUBLICROUTINE>

end subroutine update_vegetation_bnd_slow
! </SUBROUTINE>

subroutine vegetation_stock_pe ( veg, index, value )
  type(vegetation_type), intent(in) :: veg ! vegetation state
  integer , intent(in)  :: index ! ID of the stock to calculate
  real    , intent(out) :: value ! calculated value of the stock
  
  value = 0 ! vegetation doesn't have any capacity
end subroutine vegetation_stock_pe

end module vegetation_mod

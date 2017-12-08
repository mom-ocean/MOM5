module ocean_shortwave_mod
#include <fms_platform.h>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module sets up the shortwave routines. 
!</OVERVIEW>
!
!<DESCRIPTION>
! There are two shortwave routines available.  The more complete one is 
! from GFDL, and the streamlined and simpler one is from CSIRO.  
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_shortwave_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!  Must be .true. to run with module. Default is false.
!  </DATA> 
!
!  <DATA NAME="use_shortwave_gfdl=" TYPE="logical">
!  Must be .true. to run with the GFDL shortwave module. 
!  Default is true.  
!  </DATA> 
!  <DATA NAME="use_shortwave_csiro=" TYPE="logical">
!  Must be .true. to run with the CSIRO shortwave module. 
!  Default is false.  
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln 
use diag_manager_mod, only: register_diag_field
use fms_mod,          only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,          only: stdout, stdlog, FATAL, WARNING, NOTE
use mpp_mod,          only: input_nml_file, mpp_error, mpp_max, mpp_min

use ocean_domains_mod,          only: get_local_indices
use ocean_parameters_mod,       only: cp_ocean, missing_value, rho0r
use ocean_shortwave_csiro_mod,  only: ocean_shortwave_csiro_init, sw_source_csiro
use ocean_shortwave_jerlov_mod, only: ocean_shortwave_jerlov_init, sw_source_jerlov
use ocean_shortwave_gfdl_mod,   only: ocean_shortwave_gfdl_init, sw_source_gfdl
use ocean_types_mod,            only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,            only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,            only: ocean_thickness_type, ocean_options_type, ocean_density_type
use ocean_tpm_util_mod,         only: otpm_set_diag_tracer
use ocean_workspace_mod,        only: wrk1, wrk2, wrk3, wrk4, wrk5
use ocean_util_mod,             only: diagnose_2d, diagnose_3d, diagnose_sum
use ocean_tracer_util_mod,      only: diagnose_3d_rho

implicit none

private

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$Id: ocean_shortwave.F90,v 20.0 2013/12/14 00:16:14 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

public  ocean_shortwave_init
public  ocean_irradiance_init
public  sw_source
private sw_source_ext
private watermass_diag_init
private watermass_diag

! internally set for computing watermass diagnostics
logical :: compute_watermass_diag = .false. 

! for global area normalization
real    :: cellarea_r

! for diagnostics 
integer :: id_sw_frac       =-1
integer :: id_sw_heat       =-1
integer :: id_sw_heat_on_nrho=-1
integer :: id_irradiance    =-1

integer :: id_neut_rho_sw          =-1
integer :: id_wdian_rho_sw         =-1
integer :: id_tform_rho_sw         =-1
integer :: id_neut_rho_sw_on_nrho  =-1
integer :: id_wdian_rho_sw_on_nrho =-1
integer :: id_tform_rho_sw_on_nrho =-1
integer :: id_eta_tend_sw_pen      =-1
integer :: id_eta_tend_sw_pen_glob =-1


logical :: used

! inverse of specific heat
real :: cp_r

! for irradiance index
integer :: index_irr

! for FAFMIP redistributed heat tracer
integer :: index_redist_heat = -1 

logical :: module_is_initialized  = .false.
logical :: use_this_module        = .false.
logical :: use_shortwave_gfdl     = .true.
logical :: use_shortwave_csiro    = .false.
logical :: use_shortwave_jerlov   = .false.
logical :: use_shortwave_ext      = .false.

namelist /ocean_shortwave_nml/ use_this_module, use_shortwave_gfdl,       &
                               use_shortwave_csiro, use_shortwave_jerlov, &
                               use_shortwave_ext

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_shortwave_init">
!
! <DESCRIPTION>
! Initialization for the shortwave module.
! </DESCRIPTION>
  subroutine ocean_shortwave_init(Grid, Domain, Time, Dens, T_prog, vert_coordinate, Ocean_options)

    type(ocean_grid_type),        intent(in), target :: Grid
    type(ocean_domain_type),      intent(in), target :: Domain
    type(ocean_time_type),        intent(in)         :: Time 
    type(ocean_density_type),     intent(in)         :: Dens
    type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
    integer,                      intent(in)         :: vert_coordinate
    type(ocean_options_type),     intent(inout)      :: Ocean_options

    integer :: unit, io_status, ierr
    integer :: num_schemes=0
    integer :: n, num_prog_tracers
    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) return
    
    module_is_initialized = .TRUE.

    call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_shortwave_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_shortwave_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_shortwave_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_shortwave_nml')
    call close_file(unit)
#endif
    write(stdoutunit,'(/)')
    write(stdoutunit,ocean_shortwave_nml)    
    write(stdlogunit,ocean_shortwave_nml)

    Dom => Domain
    Grd => Grid

    cp_r = 1.0/cp_ocean
    cellarea_r = 1.0/(epsln + Grd%tcellsurf)

#ifndef MOM_STATIC_ARRAYS    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif 

    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING shortwave_mod.')
    else 
      call mpp_error(NOTE, '==>Note: NOT using shortwave_mod.')
      Ocean_options%shortwave = 'Did NOT use the shortwave penetration option.'
      return 
    endif 

    if(use_shortwave_gfdl) then 
       call ocean_shortwave_gfdl_init (Grid, Domain, Time, vert_coordinate, Ocean_options)
       num_schemes = num_schemes+1
    elseif(use_shortwave_csiro) then 
       call ocean_shortwave_csiro_init(Grid, Domain, Time, Ocean_options)       
       num_schemes = num_schemes+1
    elseif(use_shortwave_jerlov) then 
       call ocean_shortwave_jerlov_init(Grid, Domain, vert_coordinate, Ocean_options)       
       num_schemes = num_schemes+1
    elseif(use_shortwave_ext) then 
      call mpp_error(NOTE, &
      '==>Note: using shortwave_ext. Obtaining radiation from another model.')
       num_schemes = num_schemes+1
    endif 
    if(num_schemes==0) then 
      call mpp_error(WARNING,&
      '==>shortwave_mod: no shortwave scheme selected, yet using shortwave_mod. Is this correct?.')
    endif 
    if(num_schemes > 1) then 
      call mpp_error(FATAL,&
      '==>shortwave_mod: choose only ONE of the shortwave schemes: GFDL, CSIRO, JERLOV, or External.')
    endif

    ! to determine index for FAFMIP temperature tracer     
    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'redist_heat') index_redist_heat = n
    enddo

    ! for diagnostics      
    id_sw_frac = register_diag_field ('ocean_model', 'sw_frac',                       &
         Grid%tracer_axes(1:3), Time%model_time, 'fraction of shortwave penetrating', &
         'dimensionless', missing_value=-1e10, range=(/-1.e10,1.e10/))

    id_sw_heat = register_diag_field ('ocean_model', 'sw_heat',                   & 
         Grid%tracer_axes(1:3), Time%model_time, 'penetrative shortwave heating', &
         'W/m^2', missing_value=-1e10, range=(/-1.e10,1.e10/),                    &
         standard_name='downwelling_shortwave_flux_in_sea_water')

    id_irradiance = register_diag_field ('ocean_model', 'irradiance', &
         Grid%tracer_axes(1:3), Time%model_time, 'irradiance', &
         'W/m^2',missing_value=missing_value, range=(/-1e6,1e6/))

    call watermass_diag_init(Time, Dens)


end subroutine ocean_shortwave_init
! </SUBROUTINE> NAME="ocean_shortwave_init"



!#######################################################################
! <SUBROUTINE NAME="ocean_irradiance_init">
!
! <DESCRIPTION>
!
! Initialize the irradiance diagnostic tracer.
!
! </DESCRIPTION>

subroutine ocean_irradiance_init

  character(len=48), parameter :: sub_name = 'ocean_irradiance_init'
  character(len=48), parameter :: mod_name = 'ocean_shortwave_mod'

  ! set the irradiance diagnostic tracer
  index_irr = otpm_set_diag_tracer('irr',                             &
       longname = 'Irradiance', units = 'Watts/m^2',                  &
       const_init_tracer = .true., const_init_value = 0.0,            &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')

end subroutine ocean_irradiance_init
! </SUBROUTINE> NAME="ocean_irradiance_init">



!#######################################################################
! <SUBROUTINE NAME="sw_source">
!
! <DESCRIPTION>
!
! Choose either of the GFDL, CSIRO, JERLOV or External sw_source methods.
!
! </DESCRIPTION>
subroutine sw_source (Time, Thickness, Dens, T_prog, T_diag, swflx, swflx_vis, Temp, sw_frac_zt, opacity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:),     intent(in)    :: swflx
  real, dimension(isd:,jsd:),     intent(in)    :: swflx_vis
  type(ocean_prog_tracer_type),   intent(inout) :: Temp
  real, dimension(isd:,jsd:,:),   intent(inout) :: sw_frac_zt
  real, dimension(isd:,jsd:,:),   intent(inout) :: opacity 

  integer :: i,j,k,tau

  if (.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_shortwave_mod (sw_source): module must be initialized ')
  endif 

  tau = Time%tau

  if (Temp%name /= 'temp') then 
    call mpp_error(FATAL,&
    '==>Error in ocean_shortwave_mod (sw_source): invalid tracer for sw_source')
  endif 

  ! initialize wrk1 array used to diagnose heating from shortwave 
  Temp%wrk1(:,:,:) = 0.0

  if(use_shortwave_gfdl) then 
     call sw_source_gfdl (Time, Thickness, T_diag(:), swflx, swflx_vis, index_irr, Temp, sw_frac_zt, opacity)
  elseif(use_shortwave_csiro) then 
     call sw_source_csiro (Time, Thickness, T_diag(:), swflx, index_irr, Temp, sw_frac_zt)
  elseif(use_shortwave_jerlov) then 
     call sw_source_jerlov (Thickness, T_diag(:), swflx, swflx_vis, index_irr, Temp, sw_frac_zt, opacity)
  elseif(use_shortwave_ext) then
     call sw_source_ext(Time, Thickness, T_diag(:),  swflx, Temp, sw_frac_zt)
  endif 

  ! add heating rate to thickness*density weighted temperature
  ! tendency. cp_r factor needed to convert from W/m^2 to (kg/m^3)*(m/s)*degC.
  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        Temp%th_tendency(i,j,k) = Temp%th_tendency(i,j,k) + Temp%wrk1(i,j,k)*cp_r
      enddo
    enddo
  enddo


  ! add heating rate to thickness*density weighted redistributed temperature
  ! for the case when including the FAFMIP redistributed heat tracer.
  if(index_redist_heat > 0) then 
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
          T_prog(index_redist_heat)%th_tendency(i,j,k) = T_prog(index_redist_heat)%th_tendency(i,j,k) + Temp%wrk1(i,j,k)*cp_r
        enddo
      enddo
    enddo
  endif 


#ifndef MOM_STATIC_ARRAYS
  if(_ALLOCATED(Temp%radiation)) then
#endif
    do k=1,nk-1
      do j=jsc,jec
        do i=isc,iec
          Temp%radiation(i,j,k) = Temp%wrk1(i,j,k)
        enddo
      enddo
    enddo
#ifndef MOM_STATIC_ARRAYS
  endif
#endif


  ! diagnostics 

  call diagnose_3d(Time, Grd, id_sw_frac, sw_frac_zt(:,:,:))
  call diagnose_3d(Time, Grd, id_sw_heat, Temp%wrk1(:,:,:))
  call diagnose_3d(Time, Grd, id_irradiance, T_diag(index_irr)%field(:,:,:))

  call watermass_diag(Time, Temp, Dens)


end subroutine sw_source
! </SUBROUTINE> NAME="sw_source"


!#######################################################################
! <SUBROUTINE NAME="sw_source_ext">
!
! <DESCRIPTION>
!
! Example of a routine that applies an externally supplied shortwave 
! heating rate (i.e. top minus bottom radiation flux in W/m^2). Users
! should modify this routine for their own purposes.
!
! </DESCRIPTION>

subroutine sw_source_ext (Time, Thickness, T_diag, swflx, Temp, sw_frac_zt)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:) ,   intent(in)    :: swflx
  type(ocean_prog_tracer_type),  intent(inout) :: Temp
  real, dimension(isd:,jsd:,:),  intent(inout) :: sw_frac_zt

  integer :: k, j, i
  real :: sw_heat_rate(isd:ied,jsd:jed,1:nk) 

  ! load sw_heat_rate here
  ! sw_heat_rate = ...
  ! 
  ! for our testing purposes we'll just invoke one of our methods
  ! Note here that index_irr is "global"
  call sw_source_csiro (Time, Thickness, T_diag, swflx, index_irr, Temp, sw_frac_zt)
  ! the heating rate is stored in 
  sw_heat_rate = Temp%wrk1

  ! load sw_heat_rate
  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        Temp%wrk1(i,j,k) = sw_heat_rate(i,j,k)
      enddo
    enddo
  enddo

end subroutine sw_source_ext
! </SUBROUTINE> NAME="sw_source_ext">


!#######################################################################
! <SUBROUTINE NAME="watermass_diag_init">
!
! <DESCRIPTION>
! Initialization of watermass diagnostic output files. 
! </DESCRIPTION>
!
subroutine watermass_diag_init(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens

  integer :: stdoutunit
  stdoutunit=stdout()
  compute_watermass_diag = .false. 

  id_sw_heat_on_nrho = register_diag_field ('ocean_model', 'sw_heat_on_nrho',   &
       Dens%neutralrho_axes(1:3), Time%model_time, 'penetrative shortwave heating binned to neutral density', &
       'W/m^2', missing_value=-1e10, range=(/-1.e20,1.e20/))
  if(id_sw_heat_on_nrho > 0) compute_watermass_diag = .true.

  id_neut_rho_sw = register_diag_field ('ocean_model', 'neut_rho_sw',&
    Grd%tracer_axes(1:3), Time%model_time,                           &
    'update of locally ref potrho from shortwave penetration',       &
    '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sw > 0) compute_watermass_diag = .true. 

  id_neut_rho_sw_on_nrho = register_diag_field ('ocean_model',                            &
   'neut_rho_sw_on_nrho', Dens%neutralrho_axes(1:3), Time%model_time,                     &
   'update of locally ref potrho from shortwave penetration binned to neutral rho layers',&
   '(kg/m^3)/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_neut_rho_sw_on_nrho > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sw = register_diag_field ('ocean_model', 'wdian_rho_sw',&
    Grd%tracer_axes(1:3), Time%model_time,                             &
    'dianeutral mass transport due to penetrative shortwave heating',  &
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sw > 0) compute_watermass_diag = .true. 

  id_wdian_rho_sw_on_nrho = register_diag_field ('ocean_model',                                   &
    'wdian_rho_sw_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                            &
    'dianeutral mass transport due to penetrative shortwave heating binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_wdian_rho_sw_on_nrho > 0) compute_watermass_diag = .true. 

  id_tform_rho_sw = register_diag_field ('ocean_model',        &
    'tform_rho_sw',Grd%tracer_axes(1:3), Time%model_time,      &
    'watermass transform due to penetrative shortwave heating',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sw > 0) compute_watermass_diag = .true. 

  id_tform_rho_sw_on_nrho = register_diag_field ('ocean_model',                             &
    'tform_rho_sw_on_nrho',Dens%neutralrho_axes(1:3), Time%model_time,                      &
    'watermass transform due to penetrative shortwave heating binned to neutral rho layers',&
    'kg/sec', missing_value=missing_value, range=(/-1.e20,1.e20/))
  if(id_tform_rho_sw_on_nrho > 0) compute_watermass_diag = .true. 

  id_eta_tend_sw_pen= -1          
  id_eta_tend_sw_pen= register_diag_field ('ocean_model','eta_tend_sw_pen', &
       Grd%tracer_axes(1:2), Time%model_time,                               &
       'non-Bouss steric sea level tendency from swflx penetration', 'm/s', &
       missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_sw_pen > 0) compute_watermass_diag = .true. 

  id_eta_tend_sw_pen_glob= -1          
  id_eta_tend_sw_pen_glob= register_diag_field ('ocean_model', 'eta_tend_sw_pen_glob',&
       Time%model_time,                                                               &
       'global mean non-bouss steric sea level tendency from swflx penetration',      &
       'm/s', missing_value=missing_value, range=(/-1e10,1.e10/))
  if(id_eta_tend_sw_pen_glob > 0) compute_watermass_diag = .true. 


  if(compute_watermass_diag) then 
      write(stdoutunit,'(/a/)') &
           '==>Note: running ocean_shortwave_mod w/ compute_watermass_diag=.true.'  
  endif

end subroutine watermass_diag_init
! </SUBROUTINE> NAME="watermass_diag_init"


!#######################################################################
! <SUBROUTINE NAME="watermass_diag">
!
! <DESCRIPTION>
! Diagnose effects from shortwave heating on watermass transformation.  
! </DESCRIPTION>
!
subroutine watermass_diag(Time, Temp, Dens)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_prog_tracer_type),   intent(in) :: Temp
  type(ocean_density_type),       intent(in) :: Dens

  integer :: i,j,k,tau
  real, dimension(isd:ied,jsd:jed) :: eta_tend

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_shortwave (watermass_diag): module needs initialization ')
  endif 

  if (.not. compute_watermass_diag) return

  tau = Time%tau 

  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  wrk5(:,:,:) = 0.0

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           wrk1(i,j,k) = Grd%tmask(i,j,k)*Dens%drhodT(i,j,k)*cp_r*Temp%wrk1(i,j,k)
           wrk2(i,j,k) = wrk1(i,j,k)*Dens%rho_dztr_tau(i,j,k)
           wrk3(i,j,k) = wrk2(i,j,k)*Dens%stratification_factor(i,j,k)
           wrk4(i,j,k) = wrk1(i,j,k)*Grd%dat(i,j)*Dens%watermass_factor(i,j,k)
           wrk5(i,j,k) =-wrk1(i,j,k)/(epsln+Dens%rho(i,j,k,tau)**2) ! for eta_tend
        enddo
     enddo
  enddo

  call diagnose_3d(Time, Grd, id_neut_rho_sw, wrk2(:,:,:))
  call diagnose_3d(Time, Grd, id_wdian_rho_sw, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_tform_rho_sw, wrk4(:,:,:))
  call diagnose_3d_rho(Time, Dens, id_neut_rho_sw_on_nrho, wrk2)
  call diagnose_3d_rho(Time, Dens, id_wdian_rho_sw_on_nrho, wrk3)
  call diagnose_3d_rho(Time, Dens, id_tform_rho_sw_on_nrho, wrk4)
  call diagnose_3d_rho(Time, Dens, id_sw_heat_on_nrho, Temp%wrk1)

  if(id_eta_tend_sw_pen > 0 .or. id_eta_tend_sw_pen_glob > 0) then
      eta_tend(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eta_tend(i,j) = eta_tend(i,j) + wrk5(i,j,k)
            enddo
         enddo
      enddo
      call diagnose_2d(Time, Grd, id_eta_tend_sw_pen, eta_tend(:,:))
      call diagnose_sum(Time, Grd, Dom, id_eta_tend_sw_pen_glob, eta_tend, cellarea_r)
  endif

end subroutine watermass_diag
! </SUBROUTINE>  NAME="watermass_diag"


end module ocean_shortwave_mod

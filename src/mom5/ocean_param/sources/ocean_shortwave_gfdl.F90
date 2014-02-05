module ocean_shortwave_gfdl_mod
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John P. Dunne
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<REVIEWER EMAIL="russell.fiedler@csiro.au"> Russell Fiedler 
!</REVIEWER>
!
!<OVERVIEW>
! This module returns thickness weighted and density weighted 
! temperature tendency [deg C *m/sec *kg/m^3] from penetrative 
! shortwave heating.
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness and density weighted tendency [deg C *m/sec *kg/m^3]
! of temperature associated with penetrative shortwave heating in the upper
! ocean. Generally penetration is taken as a function of monthly optical 
! properties of the upper ocean, where optical properties are read 
! in from a file of climatological data or from an ecosystem model. 
!
! Presently there is account taken only of chlorophyll-a on the optical
! properties of ocean water.  Other particulates can be added so to 
! have a more complete picture of the ocean optical properties.
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Jerlov (1968): Optical Oceanography, Elsevier Press
! </REFERENCE>
!
! <REFERENCE>
! Morel and Antoine (1994), Heating rate in the upper ocean 
! in relation to its bio-optical state.
! Journal of Physical Oceanography vol 24 pages 1652-1664
! </REFERENCE>
!
! <REFERENCE>
! Manizza, M., C Le Quere, A. J. Watson, and E. T. Buitenhuis (2005)
! Bio-optical feedbacks among phytoplankton, upper ocean physics and 
! sea-ice in a global model. Geophys. Res. Let. 32, L05603, 
! doi:10.1029/2004GL020778
! </REFERENCE>
!
! <REFERENCE>
! Paulson and Simpson (1977)
! Irradiance measurements in the upper ocean
! Journal of Physical Oceanography vol 7 pages 952-956
! </REFERENCE>
!
! <REFERENCE>
! Rosati and Miyakoda (1988)
! A General Circulation Model for Upper Ocean Simulation
! Journal of Physical Oceanography vol 18 pages 1601-1626.
! </REFERENCE>
!
! <NOTE>
! Optimized for vector peformance by R. Fiedler (russell.fiedler@csiro.au)
! June 2003 on the Australian NEC computer. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_shortwave_gfdl_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be .true. to run with module. Default is false.
!  </DATA> 
!
!  <DATA NAME="use_sw_morel_mom4p0" TYPE="logical">
!  For backward compatibility with older simulations using 
!  MOM4.0. The new subroutine removes some confusing and unnecessary
!  logic to recompute a vertical k-index.  The differences 
!  between the old and new approach are nonzero and so will
!  result in bitwise changes to the simulation, but these changes
!  are deemed to be trivial.  Default use_sw_morel_mom4p0=.false.
!  </DATA> 
!
!  <DATA NAME="read_chl" TYPE="logical">
!  If .true. then read in climatological data of chlorophyll-a.
!  </DATA> 
!
!  <DATA NAME="optics_morel_antoine" TYPE="logical">
!  For using the Morel and Antoine optics.  This was the default in 
!  MOM4.0 for use with chlorophyll data. This scheme is NOT available
!  in MOM4p1 for use with the prognostic biology models, since it has 
!  been improved by the Manizza scheme. 
!  Default optics_morel_antoine=.false.
!  </DATA> 
!
!  <DATA NAME="optics_manizza" TYPE="logical">
!  For using the Manizza optics with chlorophyll data. Note that 
!  when running with a prognostic biology model, GFDL scientists use the 
!  Manizza optics.  
!  Default optics_manizza=.false.
!  </DATA> 
!
!  <DATA NAME="sw_frac_top" TYPE="real">
!  The fraction of shortwave radiation that should be incorporated into 
!  the sw_source array at k=1.  The generic treatment in MOM is to assume
!  that shortwave radiation is already contained inside the 
!  T_prog(index_temp)%stf field. Hence, to avoid   
!  double counting, sw_frac(k=0)=sw_frac_top should=0.0.
!  If one removes shortwave from stf, then set sw_frac_top=1.0.
!  </DATA> 
!  <DATA NAME="zmax_pen" UNITS="meter" TYPE="real">
!   Maximum depth of penetration of shortwave radiation. 
!   Below this depth, shortwave penetration is exponentially 
!   small and so is ignored.  This option formerly was useful,
!   since computation of exponentials expensive.  But with more 
!   modern computers, exponentials are cheap, so the default 
!   has been changed from 200 to 1e6, making this option irrelevant.  
!   But the option remains both for legacy purposes, and for those 
!   computers where exponentials are not cheap.  
!   Default zmax_pen=1e6.  
!  </DATA>
!  <DATA NAME="chl_default" UNITS="mg/m^3" TYPE="real">
!   Default concentration chl_default=0.08 roughly yields Jerlov Type 1A optics.
!  </DATA>
!  <DATA NAME="enforce_sw_frac" TYPE="logical">
!  To ensure the shortwave fraction is monotonically decreasing with depth. 
!  Applied only if optics_morel=.true.
!  Default enforce_sw_frac=.true. 
!  </DATA> 
!  <DATA NAME="sw_morel_fixed_depths" TYPE="logical">
!  To compute penetration assuming fixed depths via Grd%zw(k) depths.
!  This is strictly incorrect when have undulating free surface and/or 
!  generatlized vertical coordinates.  This option is here for purposes
!  of legacy, as this was done in MOM4.0 versions. The default is 
!  sw_morel_fixed_depths=.false.
!  </DATA> 
!
!  <DATA NAME="override_f_vis" TYPE="logical">
!  To fix the fraction of incoming shortwave assigned to the visible at 0.57.
!  </DATA> 
!  <DATA NAME="optics_for_uniform_chl" TYPE="logical">
!  To set the coefficients for optical model assuming the chlorophyll 
!  has a uniform distribution.  
!  Default optics_for_uniform_chl=.false.
!  </DATA> 

!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>

use axis_utils_mod,           only: frac_index
use constants_mod,            only: epsln, c2dbars
use diag_manager_mod,         only: register_diag_field
use field_manager_mod,        only: fm_get_index
use fms_mod,                  only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,                  only: stdout, stdlog, FATAL, NOTE
use mpp_mod,                  only: input_nml_file, mpp_error, mpp_max, mpp_min
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_ROUTINE
use time_interp_external_mod, only: time_interp_external, init_external_field

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: cp_ocean, rho_cp, missing_value
use ocean_parameters_mod,     only: PRESSURE, PSTAR, GEOPOTENTIAL, ZSTAR 
use ocean_types_mod,          only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,          only: ocean_thickness_type, ocean_options_type
use ocean_workspace_mod,      only: wrk1, wrk2, wrk3, wrk4 
use ocean_util_mod,           only: diagnose_2d

implicit none

private

! for diagnostics 
integer :: id_sat_chl  =-1
integer :: id_f_vis    =-1
logical :: used

! for vertical coordinate 
integer :: vert_coordinate 

! clock ids
integer :: id_sw_morel
integer :: id_sw_morel_mom4p0

integer :: index_chl
integer :: sbc_chl
logical :: verbose_flag=.false.

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS
  real, dimension(isd:ied,jsd:jed)      :: sw_fk_zt    ! sw (radiation) fractional decay on t grid
  real, dimension(isd:ied,jsd:jed)      :: sw_fk_zw    ! sw (radiation) fractional decay on w grid
  real, dimension(isd:ied,jsd:jed)      :: sat_chl     ! chlorophyll concentration (mg/m^3) (a proxy for color)   
  real, dimension(isd:ied,jsd:jed,0:nk) :: sw_frac_zw  ! fractional short wave radiation on w-points     
#else  
  real, allocatable, dimension(:,:)   :: sw_fk_zt    ! sw (radiation) fractional decay on t grid
  real, allocatable, dimension(:,:)   :: sw_fk_zw    ! sw (radiation) fractional decay on w grid
  real, allocatable, dimension(:,:)   :: sat_chl     ! chlorophyll concentration (mg/m^3) (a proxy for color)   
  real, allocatable, dimension(:,:,:) :: sw_frac_zw  ! fractional short wave radiation on w-points  
#endif

! coefficients for the optical model
real, dimension(6) :: V1_coef
real, dimension(6) :: V2_coef
real, dimension(6) :: Z1_coef
real, dimension(6) :: Z2_coef

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$Id: ocean_shortwave_gfdl.F90,v 20.0 2013/12/14 00:16:18 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'
character(len=48), parameter          :: mod_name = 'ocean_shortwave_gfdl_mod'

public ocean_shortwave_gfdl_init
public sw_source_gfdl
private sw_morel
private sw_morel_mom4p0

logical :: use_this_module        = .false.
logical :: use_sw_morel_mom4p0    = .false.
logical :: read_chl               = .false.
logical :: optics_morel_antoine   = .false.
logical :: optics_manizza         = .false.
logical :: module_is_initialized  = .FALSE.
logical :: debug_this_module      = .false. 
logical :: enforce_sw_frac        = .true. 
logical :: override_f_vis         = .true. 
logical :: sw_morel_fixed_depths  = .false. 
logical :: optics_for_uniform_chl = .false. 


! (mg/m^3) default concentration 0.08 roughly yields Jerlov Type 1A optics
real :: chl_default = 0.08  

! maximum depth (m) of solar penetration. 
! below, penetration is exponentially small and so is ignored
real :: zmax_pen    = 1e6

! set to 1.0 if do not have shortwave radiation inside of T_prog(index_temp)%stf.
real :: sw_frac_top = 0.0   

namelist /ocean_shortwave_gfdl_nml/ use_this_module, read_chl, chl_default,       &
                                    zmax_pen, sw_frac_top, debug_this_module,     &
                                    enforce_sw_frac, override_f_vis,              &
                                    sw_morel_fixed_depths, optics_for_uniform_chl,&
                                    optics_morel_antoine, optics_manizza  

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_shortwave_gfdl_init">
!
! <DESCRIPTION>
! Initialization for the shorwave module
! </DESCRIPTION>
  subroutine ocean_shortwave_gfdl_init(Grid, Domain, Time, ver_coordinate, Ocean_options)

    type(ocean_grid_type),     intent(in), target :: Grid
    type(ocean_domain_type),   intent(in), target :: Domain
    type(ocean_time_type),     intent(in)         :: Time
    integer,                   intent(in)         :: ver_coordinate
    type(ocean_options_type),  intent(inout)      :: Ocean_options

    integer :: unit, io_status, ierr, i, j
#ifdef MOM_STATIC_ARRAYS    
    real, dimension(isd:ied,jsd:jed)   :: chl_data
#else
    real, allocatable, dimension(:,:)  :: chl_data
#endif
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


    if ( module_is_initialized ) return
    
    module_is_initialized = .TRUE.
    vert_coordinate = ver_coordinate

    call write_version_number( version, tagname )

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_shortwave_gfdl_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_shortwave_gfdl_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_shortwave_gfdl_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_shortwave_gfdl_nml')
    call close_file(unit)
#endif
    write(stdoutunit,'(/)')
    write(stdoutunit,ocean_shortwave_gfdl_nml)    
    write(stdlogunit,ocean_shortwave_gfdl_nml)

    Dom => Domain
    Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif 
    
    if(use_this_module) then 
        call mpp_error(NOTE, '==>Note: USING shortwave_gfdl_mod.')

        if(optics_morel_antoine) then 
            Ocean_options%shortwave = 'Used GFDL shortwave penetration module & Morel-Antoine optics.'
            write(stdoutunit,'(a)') &
            '=>Note: Using shortwave penetration with GFDL formulaton & Morel-Antoine optics.'
        endif
        if(optics_manizza) then 
            Ocean_options%shortwave = 'Used GFDL shortwave penetration module with Manizza etal optics.'
            write(stdoutunit,'(a)') &
            '=>Note: Using shortwave penetration with GFDL formulaton & Manizza etal optics.'
        endif

        if(.not. optics_morel_antoine .and. .not. optics_manizza) then 
           call mpp_error(FATAL, '==>Error: Must choose an optics scheme for shortwave penetration.')
        endif 
        if(optics_morel_antoine .and. optics_manizza) then 
           call mpp_error(FATAL, '==>Error: Must choose just one optics scheme for shortwave penetration.')
        endif 

    else 
        call mpp_error(NOTE, '==>Note: NOT using shortwave_gfdl_mod.')
        Ocean_options%shortwave = 'Did NOT use any shortwave penetration option.'
        return 
    endif

    if(use_sw_morel_mom4p0) then 
        call mpp_error(NOTE, '==>Note: Using use_sw_morel_mom4p0=.true. Recommend use_sw_morel_mom4p0=.false.')
    endif

    if(optics_morel_antoine) then 

        if(enforce_sw_frac) then  
            write(stdoutunit,'(a)') &
                 '==>Note: enforce_sw_frac=.true. enforcing monotonic decrease of sw_frac with depth.'
        else 
            write(stdoutunit,'(a)') &
                 '==>Note: enforce_sw_frac=.false. non-monotonic sw_frac w/ some penetration profiles.'
        endif

        if(sw_morel_fixed_depths) then
            write(stdoutunit,'(a)') &
                 ' ==>Warning: sw_morel_fixed_depths=.true. is unsuitable for time varying thicknesses.'
            write(stdoutunit,'(a)') &
                 '             Time varying thicknesses are the norm in MOM, so recommend'
            write(stdoutunit,'(a)') &
                 '             setting sw_morel_fixed_depths=.false.  However, to reproduce MOM4.0'
            write(stdoutunit,'(a)') &
                 '             algorithm, then set sw_morel_fixed_depths=.true.' 
        endif

    endif


#ifndef MOM_STATIC_ARRAYS    
    allocate( chl_data(isd:ied,jsd:jed))
    allocate( sat_chl(isd:ied,jsd:jed))
    allocate( sw_fk_zt(isd:ied,jsd:jed))
    allocate( sw_fk_zw(isd:ied,jsd:jed))
    allocate( sw_frac_zw(isd:ied,jsd:jed,0:nk))
#endif
    chl_data(:,:)     = 0.0 
    sat_chl(:,:)      = 0.0
    sw_fk_zt(:,:)     = 0.0 
    sw_fk_zw(:,:)     = 0.0 
    sw_frac_zw(:,:,:) = 0.0

    ! set clock ids     
    id_sw_morel        = mpp_clock_id('(Ocean shortwave morel pen) '       ,grain=CLOCK_ROUTINE)
    id_sw_morel_mom4p0 = mpp_clock_id('(Ocean shortwave morel pen-mom4p0)' ,grain=CLOCK_ROUTINE)

    sat_chl(:,:)   = chl_default*Grd%tmask(:,:,1)

    index_chl = fm_get_index('/ocean_mod/diag_tracers/chl')
    
    ! for reading chlorophyll climatology data  
    if(read_chl) then 

        call mpp_error(NOTE, &
             '==>Note: Reading in chlorophyll-a from data file for shortwave penetration.')

        ! get the unit number "sbc_chl" for reading chl data 
        sbc_chl = init_external_field('INPUT/chl','chl',domain=Domain%domain2d)
        if (sbc_chl == -1) then 
            call mpp_error(FATAL, &
                 '==>Error in ocean_shortwave_gfdl_mod: failure to find sbc_chl data file')
        endif

        ! update chl in case of restart
        chl_data = 0.0
        call time_interp_external(sbc_chl, Time%model_time, chl_data, verbose=debug_this_module)
        do j=jsc,jec
           do i=isc,iec
              sat_chl(i,j) = Grd%tmask(i,j,1)*max(0.0,chl_data(i,j))
           enddo
        enddo

        id_sat_chl = register_diag_field ('ocean_model', 'sat_chl', &
             Grid%tracer_axes(1:2), Time%model_time, 'Chlorophyll', &
             'mg/m^3',missing_value=missing_value, range=(/-10.0,10.0/))

    elseif (index_chl > 0) then
        if(optics_morel_antoine) then 
            call mpp_error(FATAL, &
            '==>Error: optics_morel_antoine in ocean_shortwave_gfdl is not coded for prognostic chlorophyll.')
        endif
        call mpp_error(NOTE, &
             '==>Note: Using prognostic chlorophyll from biology module for shortwave penetration.')
    else
        call mpp_error(NOTE, &
             '==>Note: Setting chl=chl_default in shortwave_gfdl_init.')

    endif  ! endif for read_chl 


    if(sw_frac_top==0.0) then 
        write(stdoutunit,'(a)') &
             '=>Note: computing solar shortwave penetration. Assume stf has sw-radiation field'
        write(stdoutunit,'(a)') &
             '  included.  Hence, solar shortwave penetration effects placed in sw_source will '
        write(stdoutunit,'(a)') &
             '  subtract out the effects of shortwave at k=1 to avoid double-counting.'
    elseif(sw_frac_top==1.0) then 
        write(stdoutunit,'(a)') &
             '=>Note: computing solar shortwave penetration. Assume stf does not have sw-radiation'
        write(stdoutunit,'(a)') &
             ' field included.  Shortwave penetration effects are placed completely in sw_source.'
        write(stdoutunit,'(a)') &
             ' This is not the standard approach used in MOM.'
    elseif(sw_frac_top/=1.0 .and. sw_frac_top/=0.0) then 
        write(stdoutunit,'(a)') &
             '=>Note: Computing solar shortwave penetration. Assume a portion of sw-effects are'
        write(stdoutunit,'(a)') &
             '  included in stf and a portion in sw_source.  Are you sure you wish to do this?'
    endif


    if(optics_for_uniform_chl) then 

        ! Defining coefficients for optical model taken from Morel and Antoine (1994).

        if(optics_manizza) then 
            call mpp_error(FATAL, &
            '==>Error: optics_manizza in ocean_shortwave_gfdl is not coded for optics_for_uniform_chl.')
        endif

        ! These coefficients represent a uniform distribution of chlorophyll-a 
        ! through the water column. These may be more appropriate when using
        ! chl-a from surface-viewing satellite data. 
        write(stdoutunit,'(a)') &
        ' ==>Note: Setting optical model coefficients assuming uniform chl distribution.'
        V1_coef = (/0.353,  -0.047,  0.083,  0.047, -0.011, -0.009/)
        V2_coef = (/0.647,   0.047, -0.083, -0.047,  0.011,  0.009/)
        Z1_coef = (/1.662,  -0.605,  0.128, -0.033, -0.051, -0.004/)
        Z2_coef = (/8.541,  -8.924,  4.020, -0.077, -0.536,  0.055/)

    else 

        ! These coefficients represent a non-uniform distribution of chlorophyll-a 
        ! through the depth of the water column. 
        write(stdoutunit,'(a)') &
        ' ==>Note: Setting optical model coefficients assuming nonuniform chl distribution.'
        V1_coef = (/0.321,   0.008,   0.132,   0.038, -0.017,  -0.007/)
        V2_coef = (/0.679,  -0.008,  -0.132,  -0.038,  0.017,   0.007/)
        Z1_coef = (/1.540,  -0.197,   0.166,  -0.252, -0.055,   0.042/)
        Z2_coef = (/7.925,  -6.644,   3.662,  -1.815, -0.218,   0.502/)

    endif 

    id_f_vis = register_diag_field ('ocean_model', 'f_vis',           &
               Grid%tracer_axes(1:2), Time%model_time,                &
               'Fraction of incoming shortwave in the visible range', &
               'dimensionless',missing_value=missing_value, range=(/-10.0,10.0/))

end subroutine ocean_shortwave_gfdl_init
! </SUBROUTINE> NAME="ocean_shortwave_gfdl_init"


!#######################################################################
! <SUBROUTINE NAME="sw_source_gfdl">
!
! <DESCRIPTION>
! Add short wave penetrative heating to T_prog(index_temp)%th_tendency.
!
! Note that the divergence of shortwave for the first
! level "div_sw(0)" is compensating for the effect of having
! the shortwave component already included in the total
! surface tracer flux "stf(i,j,temp)"
!
! If the shortwave penetration routine is activated but Chlorophyll
! is not being read from data, then that implies that an ecological 
! model is being used to determine chlorophyll concentration.  
! In this case, the shortwave penetration is calcualted using the
! algorithm of
!
! Manizza, M., C Le Quere, A. J. Watson, and E. T. Buitenhuis (2005)
! Bio-optical feedbacks among phytoplankton, upper ocean physics and 
! sea-ice in a global model. Geophys. Res. Let. 32, L05603, 
! doi:10.1029/2004GL020778.
!
! This algorithm assumes that all infrared light is absorbed in the 
! top level.  It separates visible light into equal portions of 
! red and blue bands, treating separately absorption by water and
! chlorophyll.
!
! If the Chlorophyll is read from data, then we generally use the 
! Morel and Antoine optics scheme. Here, we take their approach 
! for computing a vertical profile based on the surface Chlorophyll. 
! However, one may also wish to use the Manizza scheme with 
! surface Chlorophyll data. In this case, we assume the surface
! Chlorophyll concentration is the same throughout the depth.
! This assumption is not generally good, but it does provide
! for a simple means of using Manizza etal scheme with Chlorophyll
! data. Note that GFDL scientists prefer Manizza etal for use with 
! prognostic 3d models.
!
! NOTE: Determine depths to T-points and W-points.
! This code is needed in particular for GEOPOTENTIAL, since
! depth_zwt and depth_zt for this coordinate do not include
! the surface height undulations. For the shortwave calculation,
! we wish to include the depth level undulations, unless enable
! sw_morel_fixed_depths=.true. 
!
! </DESCRIPTION>
subroutine sw_source_gfdl(Time, Thickness, T_diag, swflx, swflx_vis, index_irr, Temp, sw_frac_zt, opacity)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:),     intent(in)    :: swflx
  real, dimension(isd:,jsd:),     intent(in)    :: swflx_vis
  integer,                        intent(in)    :: index_irr
  type(ocean_prog_tracer_type),   intent(inout) :: Temp
  real, dimension(isd:,jsd:,:),   intent(inout) :: sw_frac_zt
  real, dimension(isd:,jsd:,:),   intent(inout) :: opacity

  real, dimension(isd:ied,jsd:jed)  :: f_vis
  real, dimension(isd:ied,jsd:jed)  :: zt_sw
  real, dimension(isd:ied,jsd:jed)  :: zw_sw
  real, dimension(isd:ied,jsd:jed)  :: chl_data

  real    :: div_sw 
  integer :: i, j, k
  integer :: tau

  ! zero out the wrk1 array used to diagnose heating from shortwave 
  Temp%wrk1(:,:,:) = 0.0

  if (.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_shortwave_gfdl_mod (sw_source_gfdl): module must be initialized ')
  endif 

  if (Temp%name /= 'temp') then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_pen_mod (sw_source): invalid tracer for sw_source')
  endif 

  tau = Time%tau 

  f_vis(:,:)         = 0.0
  zt_sw(:,:)         = 0.0
  zw_sw(:,:)         = 0.0
  chl_data(:,:)      = 0.0
  sw_frac_zw(:,:,:)  = 0.0  


  ! fill sat_chl with time interpolated chlorophyll data 
  if (read_chl) then
    chl_data=0.0
    call time_interp_external(sbc_chl, Time%model_time, chl_data, verbose=debug_this_module)
    do j=jsc,jec
      do i=isc,iec
        sat_chl(i,j) = Grd%tmask(i,j,1)*max(0.0,chl_data(i,j))
      enddo
    enddo
    call diagnose_2d(Time, Grd, id_sat_chl, sat_chl(:,:))
  endif

  ! F_vis is the amount of light in the shortwave verses the long wave. 
  ! F_vis=0.54 on sunny days and F_vis=0.60 on cloudy days. 
  ! In the GFDL atmospheric model, F_vis varies much more than 
  ! 0.54 <= F_vis <= 0.60.  So it is useful to include it explicitly
  ! in this module for those cases when using a realistic atmospheric
  ! radiation module.  Otherwise, we set f_vis=0.57
  if(override_f_vis) then
     do j=jsc,jec
       do i=isc,iec
       f_vis(i,j) = 0.57
      enddo
    enddo
  else

    !Use max() to ensure f_vis > 0
    do j=jsc,jec
      do i=isc,iec
        f_vis(i,j) = max(swflx_vis(i,j)/(epsln + swflx(i,j)),epsln) 
      enddo
    enddo

  endif
  call diagnose_2d(Time, Grd, id_f_vis, f_vis(:,:))

  ! zero out the fractional decay 
  sw_fk_zt(:,:) = 0.0
  sw_fk_zw(:,:) = 0.0

  ! only compute 3-D sw_fract for ocean regions shallower than zmax_pen 
  zw_sw=0.0
  zt_sw=0.0
  sw_frac_zw(:,:,0) = sw_frac_top


  ! for reading chlorophyll from data and using optics_morel_antoine
  if (read_chl .and. optics_morel_antoine) then

      ! determine depths to T-points and W-points
      wrk3(:,:,:) = 0.0
      wrk4(:,:,:) = 0.0
      if(vert_coordinate==GEOPOTENTIAL) then 
          do j=jsd,jed
             do i=isd,ied
                wrk3(i,j,1) = Thickness%dzwt(i,j,0)
                wrk4(i,j,1) = Thickness%dzt(i,j,1)
             enddo
          enddo
          do k=2,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk3(i,j,k) = wrk3(i,j,k-1) + Thickness%dzwt(i,j,k-1)
                   wrk4(i,j,k) = wrk4(i,j,k-1) + Thickness%dzt(i,j,k)
                enddo
             enddo
          enddo
      else 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk3(i,j,k) = Thickness%depth_zt(i,j,k)
                   wrk4(i,j,k) = Thickness%depth_zwt(i,j,k)
                enddo
             enddo
          enddo
      endif

      ! for using the older mom4p0 sw_morel (not recommended)
      if(use_sw_morel_mom4p0) then 

          do k=1,nk-1
             if(sw_morel_fixed_depths) then 
                 if(Grd%zw(k) <= zmax_pen) then 
                     zw_sw(isc:iec,jsc:jec) = Grd%zw(k)
                     call sw_morel_mom4p0(T_diag, zw_sw, f_vis, sw_fk_zw)
                     zt_sw(isc:iec,jsc:jec) = Grd%zt(k)
                     call sw_morel_mom4p0(T_diag, zt_sw, f_vis, sw_fk_zt)
                 else
                     sw_fk_zt(:,:) = 0.0
                     sw_fk_zw(:,:) = 0.0
                 endif
             else
                 do j=jsc,jec
                    do i=isc,iec
                       zw_sw(i,j) = wrk4(i,j,k)
                       zt_sw(i,j) = wrk3(i,j,k)
                    enddo
                 enddo
                 call sw_morel_mom4p0(T_diag, zw_sw, f_vis, sw_fk_zw)
                 call sw_morel_mom4p0(T_diag, zt_sw, f_vis, sw_fk_zt)  
             endif
             sw_frac_zt(:,:,k) = sw_fk_zt(:,:)
             sw_frac_zw(:,:,k) = sw_fk_zw(:,:)
          enddo

      else 

          do k=1,nk-1
             if(sw_morel_fixed_depths) then 
                 if(Grd%zw(k) <= zmax_pen) then 
                     zw_sw(isc:iec,jsc:jec) = Grd%zw(k)
                     call sw_morel(T_diag, zw_sw, f_vis, sw_fk_zw, k)
                     zt_sw(isc:iec,jsc:jec) = Grd%zt(k)
                     call sw_morel(T_diag, zt_sw, f_vis, sw_fk_zt, k)
                 else
                     sw_fk_zt(:,:) = 0.0
                     sw_fk_zw(:,:) = 0.0
                 endif
             else
                 do j=jsc,jec
                    do i=isc,iec
                       zw_sw(i,j) = wrk4(i,j,k)
                       zt_sw(i,j) = wrk3(i,j,k)
                    enddo
                 enddo
                 call sw_morel(T_diag, zw_sw, f_vis, sw_fk_zw, k)
                 call sw_morel(T_diag, zt_sw, f_vis, sw_fk_zt, k)  
             endif
             sw_frac_zt(:,:,k) = sw_fk_zt(:,:)
             sw_frac_zw(:,:,k) = sw_fk_zw(:,:)
          enddo

      endif  ! endif for use_sw_morel_mom4p0


      if(enforce_sw_frac) then   
          do k=2,nk-1
             do j=jsc,jec
                do i=isc,iec
                   sw_frac_zt(i,j,k) = min(sw_frac_zt(i,j,k),sw_frac_zt(i,j,k-1))
                enddo
             enddo
          enddo
      endif

  endif  ! endif for (read_chl .and. optics_morel_antoine)


  ! choose the Manizza etal optics for computing shortwave penetration 
  if (optics_manizza) then


      ! for restart reproducibility, need to use time
      ! independent depth_swt array.  more accurate is to use 
      ! depth_zwt, but use of depth_zwt breaks restart 
      ! reproducibility.   
      if(vert_coordinate==PRESSURE .or. vert_coordinate==PSTAR) then 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk4(i,j,k) = Thickness%depth_swt(i,j,k)*c2dbars
                enddo
             enddo
          enddo
      elseif(vert_coordinate==GEOPOTENTIAL .or. vert_coordinate==ZSTAR) then 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk4(i,j,k) = Thickness%depth_swt(i,j,k)
                enddo
             enddo
          enddo
      else 
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk4(i,j,k) = Thickness%depth_zwt(i,j,k)
                enddo
             enddo
          enddo
      endif


      ! wrk3 will contain the chl concentration array 
      wrk3(:,:,:) = 0.0

      ! assume chl data to be constant with depth (not always a good assumption)
      if(read_chl) then 
          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   wrk3(i,j,k) = sat_chl(i,j)
                enddo
             enddo
          enddo

      ! take chl from prognostic 3d ecosystem model 
      else 
          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   wrk3(i,j,k) = T_diag(index_chl)%field(i,j,k)
                enddo
             enddo
          enddo
      endif

      wrk1(:,:,:) = 0.0  ! temporary for red_frac_zw
      wrk2(:,:,:) = 0.0  ! temporary for blu_frac_zw

      ! compute shortwave penetration using Manizza etal optics 
      do j=jsc,jec
         do i=isc,iec
            if (Grd%tmask(i,j,1) > 0.0) then

                k=1
                wrk1(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.225+0.037*wrk3(i,j,k)**0.629))
                wrk2(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.0232+0.074*wrk3(i,j,k)**0.674))
                sw_frac_zw(i,j,k) = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
                sw_frac_zt(i,j,k) = 0.5*(f_vis(i,j) + sw_frac_zw(i,j,1))

                do k=2,nk-1
                   if (wrk4(i,j,k) < zmax_pen) then
                       wrk1(i,j,k) = wrk1(i,j,k-1) &
                            *exp(-Thickness%dzt(i,j,k)*(0.2250+0.037*wrk3(i,j,k)**0.629))
                       wrk2(i,j,k) = wrk2(i,j,k-1) &
                            *exp(-Thickness%dzt(i,j,k)*(0.0232+0.074*wrk3(i,j,k)**0.674))
                       sw_frac_zw(i,j,k)  = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
                       sw_frac_zt(i,j,k)  = &
                            0.5*(sw_frac_zw(i,j,k-1) + sw_frac_zw(i,j,k))
                   endif
                enddo

            endif
         enddo
      enddo

  endif  ! endif for optics_manizza

  
  ! when chlorophyll is being read through T_diag (not from climatology).
  ! note that this irradiance tracer is only used in MOM for purpose of 
  ! diagnostics, such as in in the residency time module.  It is NOT used
  ! by the ocean biogeochemistry, which actually computes its own irradiance
  ! as a function of the opacity.  
  if (index_irr > 0) then 
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               T_diag(index_irr)%field(i,j,k) = swflx(i,j) * sw_frac_zt(i,j,k) 
            enddo
         enddo
      enddo
  endif

  ! compute and load heating rate, to be added to   
  ! Temp%th_tendency inside ocean_shortwave.F90. 
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec           
           div_sw           = sw_frac_zw(i,j,k-1) - sw_frac_zw(i,j,k)*Grd%tmask(i,j,k+1)
           Temp%wrk1(i,j,k) = Grd%tmask(i,j,k)*swflx(i,j)*div_sw
        enddo
     enddo
  enddo

  ! Compute opacity, which is used by biogeochemistry. 
  ! We split off the k=1 level, since sw_frac_zw(k=0)=sw_frac_top,
  ! which is typically set to sw_frac_top=0.0 for purposes of 
  ! accounting (as swflx is also in stf(index_temp). A value 
  ! of sw_frac_zw(k=0)=0.0 to compute opacity would result in a
  ! negative opacity at k=1, which is not physical. Instead, for 
  ! purposes of opacity calculation, we need sw_frac_zw(k=0)=1.0. 
  k=1
  do j=jsc,jec
     do i=isc,iec           
        opacity(i,j,k) = -log( sw_frac_zw(i,j,k)/(f_vis(i,j)+epsln) + epsln) &
                          /(Thickness%dzt(i,j,k) + epsln)
     enddo
  enddo
  do k=2,nk-1
     do j=jsc,jec
        do i=isc,iec           
           opacity(i,j,k) = -log( sw_frac_zw(i,j,k)/(sw_frac_zw(i,j,k-1)+epsln) + epsln) &
                             /(Thickness%dzt(i,j,k) + epsln)
        enddo
     enddo
  enddo


end subroutine sw_source_gfdl
! </SUBROUTINE> NAME="sw_source_gfdl"


!#######################################################################
! <SUBROUTINE NAME="sw_morel">
!
! <DESCRIPTION>
!  Solar shortwave energy penetrates below the ocean surface and is aborbed 
!  by water and organic matter (both particulate and dissolved). This
!  routine estimates fraction of shortwave penetration using chlorophyll-a.
!  Absorbtion of shortwave radiation in the water assumes energy partitions
!  between three exponentials:
!
!  The first exponential is for wavelength > 0.75 um (microns) and assumes a
!  single attenuation of 0.267 m if the "zenith_angle" is 0.  Presently the 
!  code assumes a zero zenith angle, but this could be modified easily. 
!
!  The second and third exponentials represent a parameterization of the
!  attenuation coeficient for light between 300 um and 750 um in the following
!  form:
!
!	E(z) = E(0) * [V1 *  exp(z/efold1) + V2 * exp(z/efold2)]
!       with z < 0 the ocean depth 
!
!  Here, V1+V2=1 represent the partitioning between long (V1) and short (V2)
!  wavelengths between 300 um and 750 um. Thoughout most of the ocean V1<0.5
!  and V2>0.5. The "efold1" and "efold2" are the efolding depth of the long and short
!  visable and ultra violet light. Throughout most of the ocean efold1 should not exceed 3 m
!  while the efold2 will vary between 30 m in oligotrophic waters and 4 m in coastal
!  regions. All of these constants are based on satellite estimates of chlorophyll a and
!  taken from Morel and Antoine (JPO 1994, (24) 1652-1664).
!
!  If the thickness of the first ocean level "dzt(1)" is 50 meters,
!  then shortwave penetration does not do much. However, for higher
!  vertical resolution, such as dzt(1) = 10 meters commonly used
!  in ocean climate models, the effect of shortwave heating can
!  be significant. This can be particularly noticable in the summer
!  hemisphere.
!
! </DESCRIPTION>
!
! <INFO>
!
! <NOTE>
!  The terms contributing to sw_fk(i,j) are depth independent
!  when chl is depth independent.  However, we anticipate implementing 
!  a biological model, whereby chl will be depth dependent.  
! </NOTE> 
!
! <NOTE>
!  Simpson and Dickey (1981) and others have argued between one and 
!  two exponentials for light between 300 um and 750 um.  
!  With vertical grid resolution of 5 meters or finer
!  for the upper 20 meters, the second exponential will make a difference.
!  We anticipate using such resolutions, and so have implemented both 
!  exponentials. 
! </NOTE> 
!
! </INFO>
!
subroutine sw_morel (T_diag, z_sw, f_vis, sw_fk, ksw)

  type(ocean_diag_tracer_type), dimension(:), intent(in)    :: T_diag
  real, dimension(isd:,jsd:),                 intent(in)    :: z_sw     ! vertical depth
  real, dimension(isd:,jsd:),                 intent(in)    :: f_vis    ! fraction of incoming sw that is visible
  real, dimension(isd:,jsd:),                 intent(inout) :: sw_fk    ! sw fractional decay
  integer,                                    intent(in)    :: ksw      ! index of the k-level 

  ! parameters for optical model using C=log10(chl)
  real    :: V1, V2, efold1, efold2 
  real    :: C, C2, C3, C4, C5
  real    :: zenith_angle
  real    :: swmax, swmin, chmax, chmin

  integer :: i, j, kswp1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_sw_morel)  

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_gfdl_mod (sw_morel): module must be initialized')
  endif 

  zenith_angle = 0.0
  kswp1 = min(ksw+1,nk)

  ! compute shortwave fraction based on triple exponential 
  ! note that 0.02 and 60.0 provide a floor and ceiling which 
  ! keep sw_fk between 0.0 and 1.0.  These values are dependent
  ! on details of the coefficients used in the exponential. 
  ! If the coefficients change, then the floor/ceiling needs 
  ! to be reevaluated. shortwave fraction set to zero for 
  ! depths greater than zmax_pen.
  do j=jsc,jec
     do i=isc,iec

        if ( z_sw(i,j) > zmax_pen .or. Grd%tmask(i,j,kswp1) == 0.0) then
            sw_fk(i,j) = 0.0

        else

            if (.not. read_chl .and. index_chl .gt. 0) then      
              ! chlorophyll-a from biology.
              C = log10(min(max(T_diag(index_chl)%field(i,j,ksw),0.02),60.0))
            else
              ! chlorophyll-a read from climatology
              C = log10(min(max(sat_chl(i,j),0.02),60.0))
            endif
            C2 = C  * C
            C3 = C2 * C
            C4 = C3 * C
            C5 = C4 * C
            V1 = V1_coef(1) + V1_coef(2) * C + V1_coef(3) * C2  &
                 + V1_coef(4) * C3 + V1_coef(5) * C4 + V1_coef(6) * C5
            V2 = V2_coef(1) + V2_coef(2) * C + V2_coef(3) * C2  &
                 + V2_coef(4) * C3 + V2_coef(5) * C4 + V2_coef(6) * C5
            efold1 = Z1_coef(1) + Z1_coef(2) * C + Z1_coef(3) * C**2  &
                 + Z1_coef(4) * C3 + Z1_coef(5) * C4 + Z1_coef(6) * C5
            efold2 = Z2_coef(1) + Z2_coef(2) * C + Z2_coef(3) * C**2  &
                 + Z2_coef(4) * C3 + Z2_coef(5) * C4 + Z2_coef(6) * C5

            sw_fk(i,j) = (1-f_vis(i,j)) * exp( -z_sw(i,j)/(0.267 * cos(zenith_angle)) )   &
                 + f_vis(i,j)  * ( V1 * exp( -z_sw(i,j)/efold1 )                 &
                 + V2 * exp( -z_sw(i,j)/efold2 ) )

        endif

     enddo  ! i-loop finish 
  enddo  ! j-loop finish

  if(debug_this_module) then 
      swmax=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_max(swmax);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): max sw_fk=',swmax
      swmin=minval(sw_fk(isc:iec,jsc:jec))
      call mpp_min(swmin);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): min sw_fk=',swmin
      if (index_chl .le. 0) then
        chmax=maxval(sat_chl(isc:iec,jsc:jec))
        call mpp_max(chmax);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): max chl=',chmax
        chmin=minval(sat_chl(isc:iec,jsc:jec)) 
        call mpp_min(chmin);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): min chl=',chmin
      endif
  endif
  call mpp_clock_end(id_sw_morel)

end subroutine sw_morel
! </SUBROUTINE> NAME="sw_morel"



!#######################################################################
! <SUBROUTINE NAME="sw_morel_mom4p0">
!
! <DESCRIPTION>
! As in sw_morel, but uses the MOM4.0 algorithm to re-compute a k-level. 
! The recomputation is not needed, it can be costly, and produces
! no physically significant differences.  This routine is 
! retained for legacy only and it is not otherwise recommended.  
! </DESCRIPTION>
!
subroutine sw_morel_mom4p0 (T_diag, z_sw, f_vis, sw_fk)

  type(ocean_diag_tracer_type), dimension(:), intent(in)    :: T_diag
  real, dimension(isd:,jsd:),                 intent(in)    :: z_sw     ! vertical depth
  real, dimension(isd:,jsd:),                 intent(in)    :: f_vis    ! fraction of incoming sw that is visible
  real, dimension(isd:,jsd:),                 intent(inout) :: sw_fk    ! sw fractional decay

  ! flag for setting k of bottom of boundary layer           
  logical :: keep_going(isc:iec)                            
  ! introducing kb as 1d vector improves vector performance  
  integer, dimension(isc:iec) :: kb   
  integer, dimension(isc:iec) :: kb_old   

  ! parameters for optical model using C=log10(chl)
  real    :: V1, V2, efold1, efold2 
  real    :: C, C2, C3, C4, C5
  real    :: zenith_angle
  real    :: swmax, swmin, chmax, chmin

  integer :: i, j, k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_sw_morel_mom4p0)  

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_gfdl_mod (sw_morel_mom4p0): module must be initialized')
  endif 

  zenith_angle = 0.0


  ! compute shortwave fraction based on triple exponential 
  ! note that 0.02 and 60.0 provide a floor and ceiling which 
  ! keep sw_fk between 0.0 and 1.0.  These values are dependent
  ! on details of the coefficients used in the exponential. 
  ! If the coefficients change, then the floor/ceiling needs 
  ! to be reevaluated. shortwave fraction set to zero for 
  ! depths greater than zmax_pen.

  do j=jsc,jec

     ! compute kb 
     ! (this is the irrelevant part of this subroutine)
     k = 0
     kb(:)=nk                                   
     keep_going(:) = .true.
     do while (k < nk .and. any(keep_going(:))) 
        k = k+1
        do  i=isc,iec
           if (z_sw(i,j) <= Grd%zw(k) .and.  keep_going(i)) then   
               kb(i) = k
               keep_going(i) = .false.
           endif
        enddo
     enddo
     do i=isc,iec
       kb(i)=min(kb(i)+1,nk)
     enddo

     ! check with older (non-vectorized) method for kb calculation
     if(debug_this_module) then 
         do i=isc,iec 
            kb_old(i) = ceiling(frac_index(z_sw(i,j), (/0.,Grd%zw(1:nk)/)))
            kb_old(i) = min(kb_old(i),nk)
         enddo
         do i=isc,iec 
            if(kb(i) - kb_old(i) /= 0) then 
                write(stdoutunit,*) &
                'In sw_morel, kb computed two ways: kbnew(',i+Dom%ioff,',',j+Dom%joff,')= ', kb(i), &
                                  ' kbold(',i+Dom%ioff,',',j+Dom%joff,')= ',kb_old(i)
            endif
         enddo
     endif

     do i=isc,iec

        if ( z_sw(i,j) > zmax_pen .or. Grd%tmask(i,j,kb(i)) == 0) then
            sw_fk(i,j) = 0.0
        else

            if (.not. read_chl .and. index_chl .gt. 0) then      ! chlorophyll-a from biology.
              C = log10(min(max(T_diag(index_chl)%field(i,j,kb(i)),0.02),60.0))
            else                                                 ! chlorophyll-a read from climatology
              C = log10(min(max(sat_chl(i,j),0.02),60.0))
            endif
            C2 = C  * C
            C3 = C2 * C
            C4 = C3 * C
            C5 = C4 * C
            V1 = V1_coef(1) + V1_coef(2) * C + V1_coef(3) * C2  &
                 + V1_coef(4) * C3 + V1_coef(5) * C4 + V1_coef(6) * C5
            V2 = V2_coef(1) + V2_coef(2) * C + V2_coef(3) * C2  &
                 + V2_coef(4) * C3 + V2_coef(5) * C4 + V2_coef(6) * C5
            efold1 = Z1_coef(1) + Z1_coef(2) * C + Z1_coef(3) * C**2  &
                 + Z1_coef(4) * C3 + Z1_coef(5) * C4 + Z1_coef(6) * C5
            efold2 = Z2_coef(1) + Z2_coef(2) * C + Z2_coef(3) * C**2  &
                 + Z2_coef(4) * C3 + Z2_coef(5) * C4 + Z2_coef(6) * C5

            sw_fk(i,j) = (1-f_vis(i,j)) * exp( -z_sw(i,j)/(0.267 * cos(zenith_angle)) )   &
                 + f_vis(i,j)  * ( V1 * exp( -z_sw(i,j)/efold1 )                 &
                 + V2 * exp( -z_sw(i,j)/efold2 ) )
        endif

     enddo  ! i-loop finish 

  enddo  ! j-loop finish

  if(debug_this_module) then 
      swmax=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_max(swmax);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): max sw_fk=',swmax
      swmin=minval(sw_fk(isc:iec,jsc:jec))
      call mpp_min(swmin);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): min sw_fk=',swmin
      if (index_chl .le. 0) then
        chmax=maxval(sat_chl(isc:iec,jsc:jec))
        call mpp_max(chmax);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): max chl=',chmax
        chmin=minval(sat_chl(isc:iec,jsc:jec)) 
        call mpp_min(chmin);write(stdoutunit,*)'In ocean_shortwave_gfdl (sw_morel): min chl=',chmin
      endif
  endif
  call mpp_clock_end(id_sw_morel_mom4p0)

end subroutine sw_morel_mom4p0
! </SUBROUTINE> NAME="sw_morel_mom4p0"


end module ocean_shortwave_gfdl_mod

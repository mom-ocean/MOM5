module ocean_shortwave_csiro_mod
!
!<CONTACT EMAIL="Russell.Fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! This module returns thickness and density weighted temperature 
! tendency [kg/m^3 * deg C *m/sec] from penetrative shortwave heating.
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness and density weighted tendency [deg C *m/sec *kg/m^3]
! of temperature associated with penetrative shortwave heating in the upper
! ocean. Generally penetration is taken as a function of monthly optical 
! properties of the upper ocean, where optical properties are read 
! in from a file of climatological data.
!
! This module ussumes a simple single exponential decay law. The e-folding 
! depth may vary spatially and temporaly.  This routine is commonly 
! used by researchers at CSIRO Marine and Atmospheric Research in 
! Australia.  It has been optimized for vector peformance in 
! June 2003 on the Australian NEC computer. 
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Jerlov (1968)
! Optical Oceanography
! Elsevier Press
! </REFERENCE>
!
! <REFERENCE>
! Morel and Antoine (1994)
! Heating rate in the upper ocean in relation to its bio-optical state 
! Journal of Physical Oceanography vol 24 pages 1652-1664
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
! </INFO>
!
!<NAMELIST NAME="ocean_shortwave_csiro_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!  Must be .true. to run with module. Default is false.
!  </DATA> 
!
!  <DATA NAME="read_depth" TYPE="logical">
!  If .true. then read in e folding depth for radiation attenuation. 
!  </DATA> 
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
!   small and so is ignored.
!  </DATA>
!  <DATA NAME="depth_default" UNITS="mg/m^3" TYPE="real">
!   Default efolding depth = 20m.
!  </DATA>
!  <DATA NAME="enforce_sw_frac" TYPE="logical">
!  To ensure the shortwave fraction is monotonically decreasing with depth. 
!  </DATA> 
!  <DATA NAME="sw_pen_fixed_depths" TYPE="logical">
!  To compute penetration assuming fixed depths via Grd%zw(k) depths.
!  This is strictly incorrect when have undulating free surface or 
!  generatlized vertical coordinates.  This option is here for purposes
!  of legacy, as this was done in MOM4.0 versions. The default is .false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>

use axis_utils_mod,           only: frac_index
use field_manager_mod,        only: fm_get_index
use fms_mod,                  only: write_version_number, open_namelist_file
use fms_mod,                  only: close_file, check_nml_error
use fms_mod,                  only: stdout, stdlog, FATAL, NOTE
use mpp_mod,                  only: input_nml_file, mpp_error, mpp_max, mpp_min
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_ROUTINE
use time_interp_external_mod, only: time_interp_external, init_external_field

use ocean_domains_mod,        only: get_local_indices
use ocean_types_mod,          only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_thickness_type, ocean_options_type
use ocean_types_mod,          only: ocean_diag_tracer_type

implicit none

private

! clock ids
integer :: id_sw_pen

integer :: sbc_opt
logical :: verbose_flag=.false.

#include <ocean_memory.h>
  

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$Id: ocean_shortwave_csiro.F90,v 20.0 2013/12/14 00:16:16 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'

! F_vis is the amount of light in the shortwave verses the long wave.
! F_vis=0.54 on sunny days and F_vis=0.60 on cloudy days. 
real, public :: F_vis=0.57

real, public, allocatable, dimension(:,:) :: ssw_atten_depth  !  Attenuation depth for ssw (m)

public  ocean_shortwave_csiro_init
public  sw_source_csiro
private sw_pen

logical :: use_this_module        = .false.
logical :: read_depth             = .false.
logical :: module_is_initialized  = .FALSE.
logical :: debug_this_module      = .false. 
logical :: enforce_sw_frac        = .true. 
logical :: sw_pen_fixed_depths    = .false. 

real :: depth_default = 20.0  ! (m) default  efolding length
real :: zmax_pen      = 120.0 ! maximum depth (m) of solar penetration. 
                              ! below, penetration is exponentially small and so is ignored
real :: sw_frac_top   = 0.0   ! set to 1.0 if do not have shortwave radiation inside of T_prog(index_temp)%stf.

namelist /ocean_shortwave_csiro_nml/ use_this_module, read_depth, depth_default, &
                               zmax_pen, sw_frac_top, debug_this_module,         &
                               enforce_sw_frac, sw_pen_fixed_depths 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_shortwave_csiro_init">
!
! <DESCRIPTION>
! Initialization for the shortwave module
! </DESCRIPTION>
  subroutine ocean_shortwave_csiro_init(Grid, Domain, Time, Ocean_options)

    type(ocean_grid_type),    intent(in), target :: Grid
    type(ocean_domain_type),  intent(in), target :: Domain
    type(ocean_time_type),    intent(in)         :: Time
    type(ocean_options_type), intent(inout)      :: Ocean_options

    integer :: unit, io_status, ierr, i, j
#ifdef MOM_STATIC_ARRAYS    
    real, dimension(isd:ied,jsd:jed)   :: depth_data
#else
    real, allocatable, dimension(:,:)  :: depth_data
#endif
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


    if ( module_is_initialized ) return
    
    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_shortwave_csiro_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_shortwave_csiro_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_shortwave_csiro_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_shortwave_csiro_nml')
    call close_file(unit)
#endif
    write (stdoutunit,'(/)')
    write(stdoutunit,ocean_shortwave_csiro_nml)    
    write(stdlogunit,ocean_shortwave_csiro_nml)

    Dom => Domain
    Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif 

#ifndef MOM_STATIC_ARRAYS    
    allocate( depth_data(isd:ied,jsd:jed))
#endif
    allocate(ssw_atten_depth(isd:ied,jsd:jed))
    ssw_atten_depth(:,:) = depth_default

    ! even for use_this_module=.false., we need to have ssw_atten_depth
    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING shortwave_csiro_mod.')
      Ocean_options%shortwave = 'Used the shortwave penetration from the CSIRO formulaton.'
    else 
      call mpp_error(NOTE, '==>Note: NOT using shortwave_csiro_mod.')
      Ocean_options%shortwave = 'Did NOT use any shortwave penetration option.'
      return 
    endif 

    ! set clock ids     
    id_sw_pen = mpp_clock_id('(Ocean shortwave penetrate) ' ,grain=CLOCK_ROUTINE)

    if(read_depth) then 

        call mpp_error(NOTE,&
         '==>Note: Reading in attenuation length scale from data file for shortwave penetration.')

        sbc_opt = init_external_field('INPUT/ssw_atten_depth','ssw_atten_depth',domain=Domain%domain2d)
        if (sbc_opt == -1) call mpp_error(FATAL,&
        '==>Error in ocean_shortwave_csiro_mod: failure to find ssw_atten_depth data file')

        ! update depth in case of restart
        depth_data = 0.0
        call time_interp_external(sbc_opt, Time%model_time, depth_data, verbose=debug_this_module)
        do j=jsc,jec
           do i=isc,iec
              ssw_atten_depth(i,j) = depth_data(i,j)
           enddo
        enddo

    else
        call mpp_error(NOTE, '==>Note: Setting depth=depth_default in shortwave_csiro_init.')
    endif

    if(sw_frac_top==0.0) then 
        write(stdoutunit,*) &
        '=>Note: computing solar shortwave penetration. Assume stf has sw-radiation field'
        write(stdoutunit,*) &
        '  included.  Hence, solar shortwave penetration effects placed in sw_source will '
        write(stdoutunit,*) &
        '  subtract out the effects of shortwave at k=1 to avoid double-counting.'
    elseif(sw_frac_top==1.0) then 
        write(stdoutunit,*) &
        '=>Note: computing solar shortwave penetration. Assume stf does not have sw-radiation'
        write(stdoutunit,*) &
        ' field included.  Shortwave penetration effects are placed completely in sw_source.'
        write(stdoutunit,*) &
        ' This is not the standard approach used in MOM.'
    elseif(sw_frac_top/=1.0 .and. sw_frac_top/=0.0) then 
        write(stdoutunit,*) &
        '=>Note: Computing solar shortwave penetration. Assume a portion of sw-effects are'
        write(stdoutunit,*) &
        '  included in stf and a portion in sw_source.  Are you sure you wish to do this?'
    endif

    if(enforce_sw_frac) then  
        write(stdoutunit,*) &
        '=>Note: enforce_sw_frac=.true. enforcing monotonic decrease of sw_frac with depth.'
    else 
        write(stdoutunit,*) &
        '=>Note: enforce_sw_frac=.false. non-monotonic sw_frac w/ some penetration profiles.'
    endif

    if(sw_pen_fixed_depths) then
        write(stdoutunit,*) &
        ' ==>Warning: sw_pen_fixed_depths=.true. is unsuitable for time varying thicknesses.'
        write(stdoutunit,*)&
        '             Time varying thicknesses are the norm in MOM, so recommend'
        write(stdoutunit,*) &
        '             setting sw_pen_fixed_depths=.false.  However, to reproduce MOM4.0'
        write(stdoutunit,*) &
        '             algorithm, then set sw_pen_fixed_depths=.true.' 
    endif

#ifndef MOM_STATIC_ARRAYS
    deallocate(depth_data)
#endif


end subroutine ocean_shortwave_csiro_init
! </SUBROUTINE> NAME="ocean_shortwave_csiro_init"



!#######################################################################
! <SUBROUTINE NAME="sw_source_csiro">
!
! <DESCRIPTION>
! Add short wave penetrative heating to T_prog(index_temp)%th_tendency.
!
! Note that the divergence of shortwave for the first
! level "div_sw(0)" is compensating for the effect of having
! the shortwave component already included in the total
! surface tracer flux "stf(i,j,temp)"
!
! </DESCRIPTION>
subroutine sw_source_csiro (Time, Thickness, T_diag, swflx, index_irr, Temp, sw_frac_zt)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:) ,   intent(in)    :: swflx
  integer,                        intent(in)    :: index_irr
  type(ocean_prog_tracer_type),  intent(inout) :: Temp
  real, dimension(isd:,jsd:,:),  intent(inout) :: sw_frac_zt

  real, dimension(isd:ied,jsd:jed,0:nk)  :: sw_frac_zw
  real, dimension(isd:ied,jsd:jed)       :: zt_sw
  real, dimension(isd:ied,jsd:jed)       :: zw_sw
  real, dimension(isd:ied,jsd:jed)       :: depth_data
  real, dimension(isd:ied,jsd:jed)       :: sw_fk_zt  ! sw radiation fractional deacy on t-grid 
  real, dimension(isd:ied,jsd:jed)       :: sw_fk_zw  ! sw radiation fractional deacy on w-grid 

  real    :: div_sw 
  integer :: i, j, k

  ! zero out the wrk1 array used to diagnose heating from shortwave 
  Temp%wrk1(:,:,:) = 0.0

  if (.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_csiro_mod (sw_source_csiro): module must be initialized ')
  endif 

  if (read_depth) then
    depth_data=0.0
    call time_interp_external(sbc_opt, Time%model_time, depth_data, verbose=debug_this_module)
    do j=jsc,jec
      do i=isc,iec
        ssw_atten_depth(i,j) = depth_data(i,j)
      enddo
    enddo
  endif

  ! zero out the fractional decay 
  sw_fk_zt(:,:) = 0.0
  sw_fk_zw(:,:) = 0.0

  ! only compute 3-D sw_fract for ocean regions shallower than zmax_pen 
  zw_sw=0.0
  zt_sw=0.0
  sw_frac_zw(:,:,0) = sw_frac_top

  do k=1,nk-1

     if(sw_pen_fixed_depths) then 

         if(Grd%zw(k) <= zmax_pen) then 
             zw_sw(isc:iec,jsc:jec) = Grd%zw(k)
             call sw_pen(zw_sw, sw_fk_zw)
             zt_sw(isc:iec,jsc:jec) = Grd%zt(k)
             call sw_pen(zt_sw, sw_fk_zt)
         else
             sw_fk_zt(:,:) = 0.0
             sw_fk_zw(:,:) = 0.0
         endif

     else 

         do j=jsc,jec
            do i=isc,iec
               zw_sw(i,j) = Thickness%depth_zwt(i,j,k)
               zt_sw(i,j) = Thickness%depth_zt(i,j,k)
            enddo
         enddo
         call sw_pen(zw_sw, sw_fk_zw)
         call sw_pen(zt_sw, sw_fk_zt)  

     endif

     sw_frac_zt(:,:,k) = sw_fk_zt(:,:)
     sw_frac_zw(:,:,k) = sw_fk_zw(:,:)

  enddo

  if(enforce_sw_frac) then   
      do k=2,nk-1
         do j=jsc,jec
            do i=isc,iec
               sw_frac_zt(i,j,k) = min(sw_frac_zt(i,j,k),sw_frac_zt(i,j,k-1))
            enddo
         enddo
      enddo
  endif

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

  ! compute and load heating rate.
  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        div_sw           = sw_frac_zw(i,j,k-1) - sw_frac_zw(i,j,k)*Grd%tmask(i,j,k+1)
        Temp%wrk1(i,j,k) = Grd%tmask(i,j,k)*swflx(i,j)*div_sw
      enddo
    enddo
  enddo

end subroutine sw_source_csiro
! </SUBROUTINE> NAME="sw_source_csiro"


!#######################################################################
! <SUBROUTINE NAME="sw_pen">
!
! <DESCRIPTION>
!  Absorbtion of shortwave radiation in the water assumes energy partitions
!  represented by a single exponential:
!
!  The exponentialsrepresents a parameterization of the
!  attenuation coefficient for light between 300 um and 750 um in the following
!  form:
!
!	E(z) = E(0) * exp(z/efold))
!       with z < 0 the ocean depth 
!
!  The "efold" s the efolding depth of the long and short
!  visable and ultra violet light.
!  efold will vary between 30 m in oligotrophic waters and 4 m in coastal
!  regions. 
!
!  If the thickness of the first ocean level "dzt(1)" is 50 meters,
!  then shortwave penetration does not do much. However, for finer 
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
! </NOTE> 
!
! <NOTE>
!  Simpson and Dickey (1981) and others have argued between one and 
!  two exponentials for light between 300 um and 750 um.  
!  With vertical grid resolution of 5 meters or finer
!  for the upper 20 meters, a second exponential will make a difference.
! </NOTE> 
!
! </INFO>
!
subroutine sw_pen (z_sw, sw_fk)

  real, intent(in),    dimension(isd:,jsd:) :: z_sw     ! vertical depth
  real, intent(inout), dimension(isd:,jsd:) :: sw_fk    ! sw fractional decay

! flag for setting k of bottom of boundary layer           
  logical :: keep_going(isd:iec)                            
! introducing kb as 1d vector improves vector performance  
  integer, dimension(isc:iec) :: kb   
  integer, dimension(isc:iec) :: kb_old   

  real    :: swmax, swmin
  integer :: i, j, k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_sw_pen)  

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_csiro_mod (sw_pen): module must be initialized')
  endif 

  ! split the i,j loops to allow vectorization in j-loop 
  do j=jsc,jec

     ! compute kb
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
                'In sw_pen, kb computed two ways: kbnew(',i+Dom%ioff,',',j+Dom%joff,')= ', kb(i), &
                ' kbold(',i+Dom%ioff,',',j+Dom%joff,')= ',kb_old(i)
            endif
         enddo
     endif

     ! compute shortwave fraction based on sinle exponential 
     ! note that 0.02 and 60.0 provide a floor and ceiling which 
     ! keep sw_fk between 0.0 and 1.0.  These values are dependent
     ! on details of the coefficients used in the exponential. 
     ! If the coefficients change, then the floor/ceiling needs 
     ! to be reevaluated. shortwave fraction set to zero for 
     ! depths greater than zmax_pen.
     do i=isc,iec

        if ( z_sw(i,j) > zmax_pen .or. Grd%tmask(i,j,kb(i)) == 0) then
            sw_fk(i,j) = 0.0
        else

            sw_fk(i,j) =   F_vis  * exp( -z_sw(i,j)/ssw_atten_depth(i,j) ) 
        endif

     enddo  ! i-loop finish 

  enddo  ! j-loop finish

  if(debug_this_module) then 
      swmax=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_max(swmax);write(stdoutunit,*)'In ocean_shortwave_csiro (sw_pen): max sw_fk=',swmax
      swmin=maxval(sw_fk(isc:iec,jsc:jec))
      call mpp_min(swmin);write(stdoutunit,*)'In ocean_shortwave_csiro (sw_pen): min sw_fk=',swmin
  endif
  call mpp_clock_end(id_sw_pen)

end subroutine sw_pen
! </SUBROUTINE> NAME="sw_pen"


end module ocean_shortwave_csiro_mod

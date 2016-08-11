module ocean_shortwave_jerlov_mod
!
!<CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt
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
! ocean. 
!
! This module ussumes a simple double exponential decay law. The e-folding 
! depth may vary spatially and temporaly.  
!
! The exponentials represents a parameterization of the
! attenuation coefficient for light between 300 um and 750 um in the following
! form:
!
!      E(z) = E(0) * (r1*exp(-z/efold1) + (1-r1)*exp(-z/efold2))
!      with z > 0 the ocean depth 
!
! The "efold" is the efolding depth of the long and short
! visable and ultra violet light.
! efold will vary between 30 m in oligotrophic waters and 4 m in coastal
! regions. 
!
! If the thickness of the first ocean level "dzt(1)" is 50 meters,
! then shortwave penetration does not do much. However, for finer 
! vertical resolution, such as dzt(1) = 10 meters commonly used
! in ocean climate models, the effect of shortwave heating can
! be significant. This can be particularly noticable in the summer
! hemisphere.
!
! Radiation at the bottom is set to zero, hence the remaining radiation
! at the bottom of the deepest ocean cells is totally absorbed 
! by these cells. This implies, that partial cells need not to
! be considered explicitly. Radiation at tracer depth within these 
! cells is not set to zero. This differs from 
! ocean_shortwave_gfdl and ocean_shortwave_csiro and reduces slightly the
! bouyancy forcing to the mixing layer in the kpp-scheme, if surface mixing 
! goes down to the bottom. 
! 
! </DESCRIPTION>
!
! <INFO>
!
! <NOTE>
!  The efolding depth is depth independent.
! </NOTE> 
!
! <NOTE>
!  Simpson and Dickey (1981) and others have argued between one and 
!  two exponentials for light between 300 um and 750 um.  
!  With vertical grid resolution of 5 meters or finer
!  for the upper 20 meters, a second exponential will make a difference.
! </NOTE> 
!
! <REFERENCE>
! Jerlov (1968)
! Optical Oceanography
! Elsevier Press
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
!<NAMELIST NAME="ocean_shortwave_jerlov_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!   Must be .true. to run with module. Default is false.
!  </DATA> 
!  <DATA NAME="sw_frac_top" TYPE="real">
!   The fraction of shortwave radiation that should be incorporated into 
!   the sw_source array at k=1.  The generic treatment in MOM is to assume
!   that shortwave radiation is already contained inside the 
!   T_prog(index_temp)%stf field. Hence, to avoid   
!   double counting, sw_frac(k=0)=sw_frac_top should=0.0.
!   If one removes shortwave from stf, then set sw_frac_top=1.0.
!  </DATA> 
!  <DATA NAME="f_vis_in" TYPE="real">
!   F_vis is the amount of light in the shortwave versus the long wave.
!   F_vis=0.54 on sunny days and F_vis=0.60 on cloudy days. With 
!   override_f_vis  = .true. F_vis is defined from f_vis_in.
!   We believe, that this effect is in the first exponential in 
!   Paulson and Simpson (1977). The default is f_vis_in=1., instead of .57
!   but it is still possible to define this quantity. 
!  </DATA> 
!  <DATA NAME="rpart_in" TYPE="real">
!   rpart_in = (0..1)   
!  </DATA> 
!  <DATA NAME="coef1_in" UNITS="meter" TYPE="real">
!  </DATA> 
!  <DATA NAME="coef2_in" UNITS="meter"  TYPE="real">
!  </DATA> 
!  <DATA NAME="override_coeff" TYPE="logical">
!   With override_coeff  = .true. rpart_in, coef1_in, coef2_in specify
!   the parameters for the double exponential. The default is .false..
!  </DATA> 
!  <DATA NAME="override_f_vis" TYPE="logical">
!   With override_f_vis  = .true. F_vis is defined from f_vis_in,
!   otherwise it is the shortwave versus the long wave amount of light.
!   The default is .true.
!  </DATA> 
!  <DATA NAME="zmax_pen" UNITS="meter" TYPE="real">
!   Maximum depth of penetration of shortwave radiation. 
!   Below this depth, shortwave penetration is exponentially 
!   small and so is ignored.
!  </DATA>
!  <DATA NAME="baltic_optics, jerlov_1, jerlov_2, jerlov_3, jerlov_1a, jerlov_1b" TYPE="logical">
!   Logical switch to select a watertype. Default=.false.. The model stops, if none is selected
!   and override_coeff=.false..  
!  </DATA>
!  <DATA NAME="enforce_sw_frac" TYPE="logical">
!  To ensure the shortwave fraction is monotonically decreasing with depth. 
!  </DATA> 
!  <DATA NAME="sw_pen_fixed_depths" TYPE="logical">
!  To compute penetration assuming fixed depths via Grd%zw(k) depths.
!  This is strictly incorrect when have undulating free surface or 
!  generatlized vertical coordinates.  This option is here for purposes
!  of legacy, as this was done in MOM4.0 versions. It saves some compute time
!  if the surface elevation is small compared with the upper cells' thickness.
!  The default is .false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>

use constants_mod,            only: epsln
use fms_mod,                  only: write_version_number, open_namelist_file
use fms_mod,                  only: close_file, check_nml_error
use fms_mod,                  only: stdout, stdlog, FATAL, NOTE
use mpp_mod,                  only: input_nml_file, mpp_error, mpp_max, mpp_min
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_ROUTINE
use time_interp_external_mod, only: time_interp_external, init_external_field

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: GEOPOTENTIAL 
use ocean_types_mod,          only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,          only:  ocean_options_type, ocean_diag_tracer_type
use ocean_workspace_mod,      only: wrk1, wrk3, wrk4 

implicit none

private

! for vertical coordinate 
integer :: vert_coordinate 

! clock ids
integer :: id_sw_pen

logical :: verbose_flag=.false.

#include <ocean_memory.h>
  
#ifdef MOM_STATIC_ARRAYS
  real, dimension(isd:ied,jsd:jed,0:nk) :: sw_frac_zw      ! fractional short wave radiation on w-points     
  real, dimension(0:nk)                 :: sw_frac_zw_fix  ! fractional short wave radiation on w-points
                                                           ! with sw_pen_fixed_depths   
  real, dimension(1:nk)                 :: sw_frac_zt_fix  ! fractional short wave radiation on t-points
                                                           ! with sw_pen_fixed_depths   
#else  
  real, allocatable, dimension(:,:,:) :: sw_frac_zw        ! fractional short wave radiation on w-points  
  real, allocatable, dimension(:)     :: sw_frac_zw_fix    ! fractional short wave radiation on w-points  
                                                           ! with sw_pen_fixed_depths   
  real, allocatable, dimension(:)     :: sw_frac_zt_fix    ! fractional short wave radiation on t-points  
                                                           ! with sw_pen_fixed_depths   
#endif

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

character(len=128)  :: version='$Id: ocean_shortwave_jerlov.F90,v 20.0 2013/12/14 00:16:20 fms Exp $'
character (len=128) :: tagname = '$Name: tikal $'
character(len=48), parameter          :: mod_name = 'ocean_shortwave_jerlov_mod'

real         :: rpart, rscl1, rscl2, rscl_ir

public  ocean_shortwave_jerlov_init
public  sw_source_jerlov

real    :: f_vis_in               = 1.       
logical :: use_this_module        = .false.
logical :: read_depth             = .false.
logical :: module_is_initialized  = .FALSE.
logical :: debug_this_module      = .false. 
logical :: enforce_sw_frac        = .false. 
logical :: override_f_vis         = .true. 
logical :: override_coeff         = .false. 
logical :: sw_pen_fixed_depths    = .false. 
logical :: jerlov_1               = .false. 
logical :: jerlov_2               = .false. 
logical :: jerlov_3               = .false. 
logical :: jerlov_1a              = .false. 
logical :: jerlov_1b              = .false. 
logical :: baltic_optics          = .false.
real    :: rpart_in               = 0.58  !Jerlov I
real    :: coef1_in               = 0.35  !Jerlov I
real    :: coef2_in               = 23.   !Jerlov I

real :: zmax_pen      = 120.0 ! maximum depth (m) of solar penetration. 
                              ! below, penetration is exponentially small and so is ignored
real :: sw_frac_top   = 0.0   ! set to 1.0 if do not have shortwave radiation inside of T_prog(index_temp)%stf.

namelist /ocean_shortwave_jerlov_nml/ use_this_module,                   &
                               zmax_pen, sw_frac_top, debug_this_module, &  
                               enforce_sw_frac, override_f_vis,          &
                               override_coeff,                           & 
                               sw_pen_fixed_depths , baltic_optics,      &
                               jerlov_1, jerlov_2, jerlov_3, jerlov_1a,  &
                               jerlov_1b, f_vis_in,                      &
                               rpart_in, coef1_in, coef2_in 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_shortwave_jerlov_init">
!
! <DESCRIPTION>
! Initialization for the shortwave module
! </DESCRIPTION>
  subroutine ocean_shortwave_jerlov_init(Grid, Domain, ver_coordinate, Ocean_options)

    type(ocean_grid_type),    intent(in), target :: Grid
    type(ocean_domain_type),  intent(in), target :: Domain
    integer,                  intent(in)         :: ver_coordinate
    type(ocean_options_type), intent(inout)      :: Ocean_options

    real    :: coef1, coef2, efold1, efold2, depth
    integer :: unit, io_status, ierr, k, nopt

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) return
    
    module_is_initialized = .TRUE.
    vert_coordinate = ver_coordinate

    call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_shortwave_jerlov_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_shortwave_jerlov_nml')
#else
    unit = open_namelist_file()
    read(unit, ocean_shortwave_jerlov_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_shortwave_jerlov_nml')
    call close_file(unit)
#endif
    write(stdoutunit,'(/)')
    write(stdoutunit,ocean_shortwave_jerlov_nml)    
    write(stdlogunit,ocean_shortwave_jerlov_nml)

    Dom => Domain
    Grd => Grid

#ifndef MOM_STATIC_ARRAYS    
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif 
    
    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING shortwave_jerlov_mod.')
      Ocean_options%shortwave = 'Used the shortwave penetration from the Jerlov formulaton.'
    else 
      call mpp_error(NOTE, '==>Note: NOT using shortwave_jerlov_mod.')
      Ocean_options%shortwave = 'Did NOT use any shortwave penetration option.'
      return 
    endif 

#ifndef MOM_STATIC_ARRAYS    
    allocate( sw_frac_zw(isd:ied,jsd:jed,0:nk))
#endif

    ! set clock ids     
    id_sw_pen = mpp_clock_id('(Ocean shortwave penetrate) ' ,grain=CLOCK_ROUTINE)
    
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
    
    nopt = 0
    if(jerlov_1) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Jerlov I water.'
        rpart = 0.58
        coef1 = 0.35
        coef2 = 23.
        nopt = nopt + 1
    endif
    if(jerlov_1a) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Jerlov IA water.'
        rpart = 0.62
        coef1 = 0.60
        coef2 = 20.
        nopt = nopt + 1
    endif
    if(jerlov_1b) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Jerlov IB water.'
        rpart = 0.67
        coef1 = 1.0
        coef2 = 17.
        nopt = nopt + 1
    endif
    if(jerlov_2) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Jerlov II water.'
        rpart = 0.77
        coef1 = 1.5
        coef2 = 14.
        nopt = nopt + 1
    endif
    if(jerlov_3) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Jerlov III water.'
        rpart = 0.78
        coef1 = 1.4
        coef2 = 7.9
        nopt = nopt + 1
    endif
    if(baltic_optics) then 
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients assuming Baltic water.'
        rpart = 0.521
        coef1 = 0.15
        coef2 = 3.3
        nopt = nopt + 1
    endif
    if(override_coeff) then
        write(stdoutunit,*)' ==>Note: Setting optical model coefficients from namelist.'
        rpart = rpart_in
        coef1 = coef1_in
        coef2 = coef2_in
        nopt = nopt + 1
    endif 
    if (nopt == 0 ) &
        call mpp_error(FATAL, &
        'FATAL ==>ocean_shortwave_jerlov_init: No water type specified for short wave radiation.')
    if (nopt >= 2 ) &
        call mpp_error(FATAL, &
        'FATAL ==>ocean_shortwave_jerlov_init: More then one water type specified for short wave radiation.')
    
    rscl1 = -1./coef1
    rscl2 = -1./coef2
    
    write(stdoutunit,*)'=>Note: computing solar shortwave penetration with Jerlov exponentials '
    write(stdoutunit,*)'        R  = ',rpart
    write(stdoutunit,*)'        z1 = ',coef1
    write(stdoutunit,*)'        z2 = ',coef2
    if(override_f_vis) then
      write(stdoutunit,*)'=>Note: using vis/(IR+vis) fraction f_vis = ',f_vis_in 
      write(stdoutunit,*)'  A value different from 1. is not consistent with the '
      write(stdoutunit,*)'  original work of Paulson and Simpson (1977).'
    endif
    
    if(sw_pen_fixed_depths) then 
#ifndef MOM_STATIC_ARRAYS    
      allocate (sw_frac_zw_fix(0:nk))
      allocate (sw_frac_zt_fix(1:nk))
#endif
      sw_frac_zw_fix = 0.
      sw_frac_zt_fix = 0.
      do k=1,nk
         depth = Grd%zt(k)
         if ( depth <= zmax_pen ) then
            efold1   = depth * rscl1; efold2   = depth * rscl2
            sw_frac_zt_fix(k) = rpart*exp(efold1) + (1.-rpart)*exp(efold2) 
         endif
      enddo
      do k=1,nk-1
         depth = Grd%zw(k)
         if ( depth <= zmax_pen ) then
            efold1   = depth * rscl1; efold2   = depth * rscl2
            sw_frac_zw_fix(k) = rpart*exp(efold1) + (1.-rpart)*exp(efold2) 
         endif
      enddo
    endif  

end subroutine ocean_shortwave_jerlov_init
! </SUBROUTINE> NAME="ocean_shortwave_jerlov_init"



!#######################################################################
! <SUBROUTINE NAME="sw_source_jerlov">
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
subroutine sw_source_jerlov (Thickness, T_diag, swflx, swflx_vis, index_irr, Temp, sw_frac_zt, opacity)

  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_diag_tracer_type),  intent(inout) :: T_diag(:)
  real, dimension(isd:,jsd:) ,   intent(in)    :: swflx
  real, dimension(isd:,jsd:),    intent(in)    :: swflx_vis
  integer,                       intent(in)    :: index_irr
  type(ocean_prog_tracer_type),  intent(inout) :: Temp
  real, dimension(isd:,jsd:,:),  intent(inout) :: sw_frac_zt
  real, dimension(isd:,jsd:,:),  intent(inout) :: opacity

  real, dimension(isd:ied,jsd:jed)       :: f_vis
    
  real    :: efold1, efold2, efold_ir, zenith_angle, depth
  real    :: swmax, swmin

  real    :: div_sw 
  integer :: i, j, k, kp1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_sw_pen)  

  ! zero out the wrk1 array used to diagnose heating from shortwave 
  Temp%wrk1(:,:,:) = 0.0

  if (.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_jerlov_mod (sw_source_jerlov): module must be initialized ')
  endif 
  if (Temp%name /= 'temp') then 
    call mpp_error(FATAL, &
    '==>Error in ocean_shortwave_pen_mod (sw_source): invalid tracer for sw_source')
  endif 

  ! only compute 3-D sw_fract for ocean regions shallower than zmax_pen 
  f_vis(:,:)        = 0.0
  sw_frac_zw(:,:,:) = 0.
  sw_frac_zt(:,:,:) = 0.
  sw_frac_zw(:,:,0) = sw_frac_top

  if(override_f_vis) then
     do j=jsc,jec
       do i=isc,iec
         f_vis(i,j) = f_vis_in
       enddo
    enddo
  else
    do j=jsc,jec
      do i=isc,iec
        f_vis(i,j) = max(swflx_vis(i,j)/(epsln + swflx(i,j)),epsln)
      enddo
    enddo
  endif
  
! If the zenith angle comes from astronomy, move these lines into the i,j-loops
  zenith_angle = 0.0
!  rscl_ir = -1./(0.267 * cos(zenith_angle))
  rscl_ir = -1./0.267 
  
  if(sw_pen_fixed_depths) then 

     do k=1,nk
        depth = Grd%zt(k)
        efold_ir = exp( depth * rscl_ir )
        do j=jsc,jec
           do i=isc,iec
              if ( depth <= zmax_pen .and. Grd%tmask(i,j,k) /= 0) then
                 sw_frac_zt(i,j,k) = (1-f_vis(i,j)) *efold_ir    &
                          + f_vis(i,j)  * sw_frac_zt_fix(k)
              endif
           enddo  ! i-loop finish 
        enddo  ! j-loop finish 
     enddo
     do k=1,nk-1
        kp1 = k+1
        depth = Grd%zw(k)
        efold_ir = exp( depth * rscl_ir )
        do j=jsc,jec
           do i=isc,iec
              if ( depth <= zmax_pen .and. Grd%tmask(i,j,kp1) /= 0) then
                 sw_frac_zw(i,j,k) = (1-f_vis(i,j)) * efold_ir   &
                          + f_vis(i,j)  * sw_frac_zw_fix(k)
              endif
           enddo  ! i-loop finish 
        enddo  ! j-loop finish 
     enddo

  else

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
     
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              depth = wrk3(i,j,k)
              if ( depth <= zmax_pen .and. Grd%tmask(i,j,k) /= 0) then
                 efold1   = depth * rscl1
                 efold2   = depth * rscl2
                 efold_ir = depth * rscl_ir
                 sw_frac_zt(i,j,k) = (1-f_vis(i,j)) * exp( efold_ir )   &
                          + f_vis(i,j)  * ( rpart*exp(efold1) + (1.-rpart)*exp(efold2) )
              endif
           enddo  ! i-loop finish 
        enddo  ! j-loop finish 
     enddo

!The deepest interface is never to another ocean cell
!Hence, sw_frac_zw(i,j,nk) = 0 in any case
     do k=1,nk-1 
        kp1 = k+1
        do j=jsc,jec
           do i=isc,iec
              depth = wrk4(i,j,k)
              if ( depth <= zmax_pen .and. Grd%tmask(i,j,kp1) /= 0) then
                 efold1   = depth * rscl1
                 efold2   = depth * rscl2
                 efold_ir = depth * rscl_ir
                 sw_frac_zw(i,j,k) = (1-f_vis(i,j)) * exp( efold_ir )   &
                          + f_vis(i,j)  * ( rpart*exp(efold1) + (1.-rpart)*exp(efold2) )
              endif
           enddo  ! i-loop finish 
        enddo  ! j-loop finish 
     enddo
  endif

  if(enforce_sw_frac) then   
      do k=2,nk
         do j=jsc,jec
            do i=isc,iec
               sw_frac_zt(i,j,k) = min(sw_frac_zt(i,j,k),sw_frac_zt(i,j,k-1))
            enddo
         enddo
      enddo
  endif

  if (index_irr > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               T_diag(index_irr)%field(i,j,k) = swflx(i,j) * sw_frac_zt(i,j,k) 
            enddo
         enddo
      enddo
  endif

  ! compute and load heating rate.
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        div_sw           = sw_frac_zw(i,j,k-1) - sw_frac_zw(i,j,k)
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
  if(debug_this_module) then 
     do k=2,nk
        swmax=maxval(sw_frac_zw(isc:iec,jsc:jec,k))
        call mpp_max(swmax);write(stdoutunit,*)'In ocean_shortwave_jerlov : max sw_fk_zw=',swmax
        swmin=maxval(sw_frac_zw(isc:iec,jsc:jec,k))
        call mpp_min(swmin);write(stdoutunit,*)'In ocean_shortwave_jerlov : min sw_fk_zw=',swmin
     enddo
  endif
  
  call mpp_clock_end(id_sw_pen)

end subroutine sw_source_jerlov
! </SUBROUTINE> NAME="sw_source_jerlov"

end module ocean_shortwave_jerlov_mod

module idealized_ic_mod 
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Z. Liang  </CONTACT>
  !<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies </REVIEWER>

  !<OVERVIEW>
  ! For preparing idealized mom4 initial conditions
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This program prepares idealized initial conditions (active and passive tracer). 
  ! Results are netcdf files that are then read into mom4.
  ! Various idealized options are available as selected by a namelist.
  !</DESCRIPTION>
  !

  use fms_mod,         only : write_version_number, stdout, open_namelist_file
  use fms_mod,         only : file_exist, check_nml_error, read_data
  use mpp_mod,         only : mpp_chksum, mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, NOTE 
  use mpp_io_mod,      only : mpp_open, mpp_close, mpp_write_meta, mpp_write, mpp_get_axis_data
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_get_atts, mpp_get_axis_data, axistype
  use mpp_io_mod,      only : MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_OVERWR, MPP_NETCDF, fieldtype
  use mpp_domains_mod, only : domain2d, mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_global_field, mpp_domains_set_stack_size
  use constants_mod,   only : RADIAN, grav, rho0, rho0r

  implicit none
  private

  integer, parameter :: max_tracer = 10

  !--- namelist ----------------------------------------------------------
  !<NAMELIST NAME="idealized_ic_nml">
  !  <DATA NAME="temp_type" TYPE="character(len=128)" >
  !  Control the iealized initial temperature condition options. 
  !  There are various options available and the default value is 
  !  "constant_temp_ic". When temp_type is
  !  1. = "constant_temp_ic", use spatially constant initial potential temperature.
  !  2. = "exponential_temp_ic", use initial potential temperature that is 
  !        exponential in the vertical.
  !  3. = "equatorial_temp_ic", use initial temp condition for idealized equatorial studies.
  !  4. = "shelfbowl_temp_ic", Initial conditions for Winton etal shelf-bowl test case.  
  !        Use tanh transition between cold shelf and warm bowl waters, instead of Heaviside. 
  !  5. = "zonal_levitus_temp_ic", use zonal average Levitus temp as initial conditions. 
  !  6. = "zonal_levitus_temp_dome_ic" zonal Levitus, except constant_temp_value in embayment
  !  7. = "temp_for_dome_ic" linear stratification with special treatment within embayment
  !  8. = "ideal_thermocline", idealized thermocline with profile from Smith and Vallis (2001). 
  !  </DATA>
  !
  !  <DATA NAME="salt_type" TYPE="character(len=128)" >
  !  Control the idealized initial salinity condition options. 
  !  The following options for salt_type are available, with default 
  !  "constant_salt_ic".
  !  "constant_salt_ic", use spatially constant initial salinity.
  !  "exponential_salt_ic", use initial salinity that is exponential in the vertical.
  !  "salinity_profile_ic", use initial salinity condition as set by profile in function salt0
  !  "salinity_for_dome_ic", salinity=1.0 for inflow in dome embayment and 0.0 elsewhere.
  !  "salinity_layer_ic"  , salinity=1.0 in a layer and 0.0 elsewhere. For passive salinity experiments. 
  !  "salinity_patch_ic"  , salinity=1.0 in a patch and 0.0 elsewhere. For passive salinity experiments. 
  !  "shelfbowl_salinity_ic", salinity=1.0 on shelf, and 0.0 elsewhere. For passive salinity experiments. 
  !  </DATA>
  !
  !  <DATA NAME="age_type" TYPE="character(len=128)" >
  !  Control the idealized initial age condition options. 
  !  There is only one option now available, and the default value is 
  !  "constant_age_ic". 
  !  </DATA>
  !
  !  <DATA NAME="constant_temp_value" TYPE="real">
  !  Value for uniform initial temp.
  !  </DATA>
  !  <DATA NAME="constant_salt_value" TYPE="real">
  !  Value for uniform initial salinity.
  !  </DATA>
  !
  !  <DATA NAME="efold_depth" TYPE="real" UNITS="metre">
  !  The efolding depth used for exponential temp or salinity profile. 
  !  Default=1000.0
  !  </DATA>
  !
  !  <DATA NAME="ideal_thermocline_deltaT" TYPE="real" UNITS="degC">
  !  Difference in temperature between bottom and surface for 
  !  the ideal thermocline initial conditions. Default
  !  ideal_thermocline_deltaT=20.0
  !  </DATA>
  !  <DATA NAME="ideal_thermocline_scale_thick" TYPE="real" UNITS="dimensionless">
  !  Dimensionless scale thickness for the ideal thermocline profile.  The dimensionful
  !  scale thickness is Hscale = H*ideal_thermocline_scale_thick, with H the ocean bottom.
  !  Default ideal_thermocline_scale_thick=0.10  (range of 0.05 to 0.15 recommended by 
  !  Smith and Vallis (2001).
  !  </DATA>
  !  <DATA NAME="ideal_thermocline_rho0" TYPE="real" UNITS="kg/m3">
  !  Density scale for use with the ideal thermocline configuration.
  !  Default ideal_thermocline_rho0=1035. 
  !  </DATA>
  !  <DATA NAME="ideal_thermocline_alpha" TYPE="real" UNITS="dimensionless">
  !  Linear scaling for thermocline. Default is ideal_thermocline_alpha=0.0005.
  !  </DATA>
  !  <DATA NAME="ideal_thermocline_alpha_eos" TYPE="real" UNITS="kg/(m3*C)">
  !  Linear equation of state parameter for ideal thermocline, assuming 
  !  rho = rho0 - ideal_thermocline_alpha_eos*theta. Default 
  !  ideal_thermocline_alpha_eos=0.255.
  !  </DATA>
  !  <DATA NAME="ideal_thermocline_offset" TYPE="real" UNITS="C">
  !  Offset to make the temperature profile realistic. 
  !  Default ideal_thermocline_offset=22.0.
  !  </DATA>
  !
  !  <DATA NAME="linear_theta_strat" TYPE="real" UNITS="degreesC/metre">
  !  The linear vertical stratification for theta for reference in DOME configuration. 
  !  </DATA>
  !
  !  <DATA NAME="dome_debug" TYPE="logical">
  ! For debugging DOME conditions. 
  !  </DATA>
  !  <DATA NAME="dome_embayment_west" TYPE="real" UNITS="dimensionless">
  ! western edge of dome embayment. Default=18.75
  !  </DATA>
  !  <DATA NAME="dome_embayment_east" TYPE="real" UNITS="dimensionless">
  ! eastern edge of dome embayment. Default=20.75
  !  </DATA>
  !  <DATA NAME="dome_embayment_north" TYPE="real" UNITS="dimensionless">
  ! northern edge of dome embayment. Default=69.75.  This is a resolution dependent
  ! value, and should be determined for each grid used. 
  !  </DATA>
  !  <DATA NAME="dome_embayment_south" TYPE="real" UNITS="dimensionless">
  ! southern edge of dome embayment. Default=69.0. This is a resolution dependent
  ! value, and should be determined for each grid used. 
  !  </DATA>
  !  <DATA NAME="dome_embayment_depth" TYPE="real" UNITS="m">
  ! Depth of the embayment. Default=600.0
  !  </DATA>
  !  <DATA NAME="dome_bottom" TYPE="real" UNITS="m">
  ! Depth of the deepest water in DOME configuration. Default=3600.0
  !  </DATA>
  !  <DATA NAME="dome_embayment_interface" TYPE="real" UNITS="m">
  ! Depth (m) determining the inflow mass flux.  Default=300.0
  !  </DATA>
  !  <DATA NAME="dome_embayment_coriolis" TYPE="real" UNITS="1/s">
  ! reference Coriolis parameter for determining Rossby radius in 
  ! the DOME inflow.  Default taken at 70N.
  !  </DATA>
  !  <DATA NAME="dome_crit_richardson" TYPE="real" UNITS="1/s">
  ! Critical Richardson number used for dome. When running with 
  ! coarse models, dome_crit_richardson=1/1000 is what Legg recommends,
  ! whereas in refined models (finer than 50km), dome_crit_richardson=1/3.
  ! Default is dome_crit_richardson=.001.
  !  </DATA>
  !  <DATA NAME="dome_drho_inflow" TYPE="real" UNITS="kg/m3">
  ! Density difference between inflow and reference for DOME 
  ! configuration. Default dome_drho_ref=2.0
  !  </DATA>
  !
  !  <DATA NAME="salinity_layer_ztop" TYPE="real">
  !  Depth at the top of the salinity layer. 
  !  </DATA>
  !  <DATA NAME="salinity_layer_zbot" TYPE="real">
  !  Depth at the bottom of the salinity layer. 
  !  </DATA>
  !  <DATA NAME="salinity_patch_ztop" TYPE="real">
  !  Depth at the top of the salinity patch. 
  !  </DATA>
  !  <DATA NAME="salinity_patch_zbot" TYPE="real">
  !  Depth at the bottom of the salinity patch. 
  !  </DATA>
  !  <DATA NAME="salinity_patch_ratio1" TYPE="real">
  !  For setting position of salinity patch.
  !  </DATA>
  !  <DATA NAME="salinity_patch_ratio2" TYPE="real">
  !  For setting position of salinity patch.
  !  </DATA>
  !
  !  <DATA NAME="constant_age_value" TYPE="real">
  !  Value for uniform initial age.
  !  </DATA>
  !
  !  <DATA NAME="generate_zonal_velocity_ic" TYPE="logical" >
  !  For generating the zonal velocity at the t-cell. 
  !  </DATA>
  !  <DATA NAME="generate_merid_velocity_ic" TYPE="logical" >
  !  For generating the meridional velocity at the t-cell. 
  !  </DATA>
  !  <DATA NAME="zonal_velocity_name" TYPE="character(len=128)" >
  !  name array of the zonal velocity component
  !  to be generated. Its default value is 'ut_inflow'.
  !  </DATA>
  !  <DATA NAME="merid_velocity_name" TYPE="character(len=128)" >
  !  name array of the meridional velocity component
  !  to be generated. Its default value is 'vt'.
  !  </DATA>
  !  <DATA NAME="zonal_velocity_file" TYPE="character(len=128)" >
  !  zonal velocity output file. Default is 'ut.res.nc'
  !  </DATA>
  !  <DATA NAME="merid_velocity_file" TYPE="character(len=128)" >
  !  meridional velocity output file.  Default is 'vt.res.nc'
  !  </DATA>
  !
  !  <DATA NAME="num_active_tracer" TYPE="integer" >
  !  Number of active tracers will be generated. Its value should be 0, 1, 2.
  !  Its default value is 2. ( temp and salt )
  !  </DATA>
  !  <DATA NAME="active_tracer" TYPE="character(len=128),dimension(2)" >
  !  name array of the active tracers will be generated. 
  !  Its element value should be 'temp' or 'salt'
  !  </DATA>
  !  <DATA NAME="active_tracer_file" TYPE="character(len=128)" >
  !  active tracer output file.
  !  </DATA>
  !
  !  <DATA NAME="num_passive_tracer" TYPE="integer" >
  !  Number of passive tracers will be generated. Its value should be no more than max_tracer.
  !  Its default value is 1 (age). 
  !  </DATA>
  !  <DATA NAME="" TYPE="character(len=128),dimension(max_tracer)" >
  !  name array of the passive tracers will be generated. 
  !  Its element value should be 'age'.
  !  </DATA>
  !  <DATA NAME="passive_tracer_file" TYPE="character(len=128)" >
  !  passive tracer output file.
  !  </DATA>
  !
  !  <DATA NAME="grid_file" TYPE="character(len=128)" >
  !  grid descriptor file.
  !  </DATA>
  !  <DATA NAME="t1" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="t0" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="z0" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !  <DATA NAME="thk" TYPE="real">
  !  For setting idealized vertical profile 
  !  </DATA>
  !
  !  <DATA NAME="write_time_axis" TYPE="logical">
  !  For writing a time axis for the IC. This is useful if 
  !  wish to use the IC as a static dataset for sponges. 
  !  Default write_time_axis=.true. 
  !  </DATA>
  ! </NAMELIST>

  character(len=128) :: temp_type = 'constant_temp_ic'
  character(len=128) :: salt_type = 'constant_salt_ic'
  character(len=128) :: age_type  = 'constant_age_ic'

  logical :: dome_debug=.false.           ! for debugging 
  real :: dome_embayment_west=18.75       ! western edge of dome embayment (degrees)
  real :: dome_embayment_east=20.75       ! eastern edge of dome embayment (degrees)
  real :: dome_embayment_south=69.0       ! southern edge of the dome embayment (degrees)
  real :: dome_embayment_north=69.75      ! northern edge of the dome embayment (degrees)
  real :: dome_embayment_depth=600.0      ! depth of the dome embayment (m)
  real :: dome_bottom=3600.0              ! depth of deepest water in dome configuration (m)
  real :: dome_embayment_interface=300.0  ! depth (m) determining the inflow mass flux 
  real :: dome_embayment_coriolis=1.37e-4 ! reference Coriolis parameter 
  real :: dome_drho_inflow=2.0            ! increase in density going from reference to inflow 
  real :: dome_crit_richardson=0.001      ! critical Richardson number for DOME inflow 

  real :: deg2metre=111.324e3             ! convert degrees latitude to metres  

  real :: ideal_thermocline_deltaT      = 20.0
  real :: ideal_thermocline_scale_thick = 0.10
  real :: ideal_thermocline_rho0        = 1035.0
  real :: ideal_thermocline_alpha       = 5e-4
  real :: ideal_thermocline_alpha_eos   = 0.255
  real :: ideal_thermocline_offset      = 20.0

  real               :: constant_temp_value      = 12.0
  real               :: constant_salt_value      = 34.72

  real               :: linear_theta_strat       = 2.19e-3

  real               :: salinity_layer_ztop      = 160.0
  real               :: salinity_layer_zbot      = 165.0
  real               :: salinity_patch_ztop      = 160.0
  real               :: salinity_patch_zbot      = 165.0
  real               :: salinity_patch_ratio1    = .333333333
  real               :: salinity_patch_ratio2    = .666666666

  real               :: efold_depth=1000.0

  real               :: constant_age_value       = 0.0
  real               :: t0=7.5, t1=10.0, z0=30.0, thk=80.0
  integer            :: num_active_tracer = 2   
  character(len=128) :: grid_file = 'grid_spec.nc'

  logical            :: generate_zonal_velocity_ic=.false.
  logical            :: generate_merid_velocity_ic=.false.
  character(len=128) :: zonal_velocity_name = 'ut'
  character(len=128) :: merid_velocity_name = 'vt'
  character(len=128) :: zonal_velocity_file = 'ut.res.nc'
  character(len=128) :: merid_velocity_file = 'vt.res.nc'

  character(len=128) :: active_tracer(2)
  character(len=128) :: active_tracer_file = 'ocean_temp_salt.res.nc'

  integer            :: num_passive_tracer = 1
  character(len=128) :: passive_tracer(max_tracer)
  character(len=128) :: passive_tracer_file = 'ocean_age.res.nc'

  logical            :: write_time_axis=.false. 

  namelist /idealized_ic_nml/ temp_type, salt_type, constant_temp_value, constant_salt_value, &
       salinity_layer_ztop, salinity_layer_zbot, salinity_patch_ztop, salinity_patch_zbot,    &
       salinity_patch_ratio1, salinity_patch_ratio2,                                          &
       t1, t0, z0, thk,  num_active_tracer, active_tracer, active_tracer_file,                &
       num_passive_tracer, passive_tracer, passive_tracer_file, grid_file,                    &
       generate_zonal_velocity_ic, generate_merid_velocity_ic,                                &
       zonal_velocity_name, merid_velocity_name, zonal_velocity_file, merid_velocity_file,    &
       efold_depth, write_time_axis, linear_theta_strat,                                      &
       dome_embayment_west, dome_embayment_east, dome_embayment_south, dome_embayment_north,  &
       dome_embayment_depth, dome_bottom, dome_embayment_interface, dome_embayment_coriolis,  &
       dome_drho_inflow, dome_crit_richardson, dome_debug,                                    &
       ideal_thermocline_deltaT, ideal_thermocline_scale_thick,                               &
       ideal_thermocline_rho0, ideal_thermocline_alpha, ideal_thermocline_offset

  !--- version information 
  character(len=128) :: version = '$Id: idealized_ic.f90,v 14.0 2007/03/15 22:47:01 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 

  !--- output data -------------------------------------------------------
  real, dimension(:,:,:), allocatable :: temp, salt, age
  real, dimension(:,:,:), allocatable :: global_temp, global_salt, global_age
  real, dimension(:,:,:), allocatable :: ut, vt
  real, dimension(:,:,:), allocatable :: global_ut, global_vt
  !--- other variables
  logical                           :: generate_temp_ic     = .FALSE.
  logical                           :: generate_salt_ic     = .FALSE.
  logical                           :: generate_age_ic      = .FALSE.
  real, dimension(:,:), allocatable :: xt, yt, kmt 
  real, dimension(:),   allocatable :: zt, grid_x_t, grid_y_t
  logical                           :: module_is_initialized = .FALSE.
  integer                           :: ni, nj, nk, isc, iec, jsc, jec
  type(domain2d),save               :: domain
  !--- public interface --------------------------------------------------
  public  idealized_ic_init, write_idealized_ic_data, idealized_ic_end

contains

  !#######################################################################
  ! <SUBROUTINE NAME="idealized_ic_init">
  !
  ! <DESCRIPTION>
  ! Initialize the module generating ideal initial conditions.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine idealized_ic_init

    !--- local variables -------------------------------------------------
    integer :: unit, io_status, ierr, n
    integer :: chksum_t, chksum_s, chksum_a
    !---------------------------------------------------------------------

    if ( module_is_initialized ) call mpp_error(FATAL, 'idealized_ic_mod: attempted reinitialization')

    module_is_initialized = .TRUE.

    !--- default value of tracer_name
    active_tracer(1)  = 'temp'
    active_tracer(2)  = 'salt'
    passive_tracer(1) = 'age'

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, idealized_ic_nml,iostat=io_status, end = 10)
          ierr = check_nml_error(io_status,'idealized_ic_nml')
       enddo
10     call mpp_close (unit)
    else
       call mpp_error(NOTE,'idealized_ic_mod: file input.nml does not exist, use default namelist option')  
    endif

    if(num_active_tracer .gt. 2) call mpp_error(FATAL, &
            'idealized_ic_mod: num_active_tracer should be the value of 0, 1 or 2')

    if(num_passive_tracer .gt. max_tracer) call mpp_error(FATAL, &
            'idealized_ic_mod: num_passive_tracer should be no greater than max_tracer')   
 
    do n = 1, num_active_tracer
       select case(active_tracer(n))
       case('temp')
          generate_temp_ic = .true.
       case('salt')
          generate_salt_ic = .true.
       case default
          call mpp_error(FATAL, trim(active_tracer(n))//' is not a valid option of nml active_tracer')
       end select
    enddo

    do n = 1, num_passive_tracer
       select case(passive_tracer(n))
       case('age')
          generate_age_ic = .true.
       case default
          call mpp_error(FATAL, trim(passive_tracer(n))//' is not a valid option of nml passive_tracer')
       end select
    enddo

    !--- at least one of generate_temp_ic and generate_salt_ic should be true
    if(.not.generate_temp_ic .and. .not.generate_salt_ic .and. .not.generate_age_ic)  &
      call mpp_error(FATAL, &
         'idealized_ic_mod: nml "generate_temp_ic"="generate_salt_ic"="generate_salt_ic"=false.  At least one should be true')

    !--- write out version information and namelist option ---------------
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
       write (stdout(), nml= idealized_ic_nml)
    endif

    !--- get grid information
    call get_grid

    !--- allocate memory to module variables
    allocate(temp(isc:iec,jsc:jec,nk), salt(isc:iec,jsc:jec,nk), age(isc:iec,jsc:jec,nk) )
    allocate(global_temp(ni,nj,nk), global_salt(ni,nj,nk), global_age(ni,nj,nk) )

    ! idealized temperature initial conditions 
    if(generate_temp_ic) then
       call idealized_temp_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,temp,global_temp)
       chksum_t = mpp_chksum(temp)
       write (stdout(),'(/(a,I20/))')' Idealized_ic Checksum: T=', chksum_t 
    endif

    ! idealized salinity initial conditions 
    if(generate_salt_ic) then 
       call idealized_salt_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,salt,global_salt)
       chksum_s = mpp_chksum(salt)
       write (stdout(),'(/(a,I20/))') ' Idealized_ic Checksum: S=', chksum_s
    endif

    ! idealized age initial conditions 
    if(generate_age_ic) then 
       call idealized_age_ic()
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,age,global_age)
       chksum_a = mpp_chksum(age)
       write (stdout(),'(/(a,I20/))') ' Idealized_ic Checksum: age=', chksum_a
    endif

    if(generate_zonal_velocity_ic) then 
       allocate(ut(isc:iec,jsc:jec,nk))
       allocate(global_ut(ni,nj,nk))
       call idealized_zonal_velocity_ic() 
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,ut,global_ut)
       chksum_t = mpp_chksum(ut)
       write (stdout(),'(/(a,I20/))')' Idealized_ic Checksum: ut=', chksum_t 
    endif 
    if(generate_merid_velocity_ic) then 
       allocate(vt(isc:iec,jsc:jec,nk))
       allocate(global_vt(ni,nj,nk))
       call idealized_merid_velocity_ic() 
       call mpp_domains_set_stack_size(ni*nj*nk*2)
       call mpp_global_field(domain,vt,global_vt)
       chksum_t = mpp_chksum(vt)
       write (stdout(),'(/(a,I20/))')' Idealized_ic Checksum: vt=', chksum_t 
    endif 


    return

  end subroutine idealized_ic_init

  !#######################################################################
  !--- idealized temperature initial condition
  subroutine idealized_temp_ic

    real    :: t_value, fact
    integer :: i, j, k, kb

    select case (trim(temp_type))
    case('constant_temp_ic')
       t_value = constant_temp_value
       temp =  t_value
       write (stdout(),'(/a,(f6.2,a)/)') ' Note: constant temp IC is being set. T=',t_value,' deg C.'
    case('equatorial_temp_ic')
       write (stdout(),'(/a/)') ' Note: idealized equatorial temp IC is being set.'
       do k=1,nk
          if(z0 == 0.0) call mpp_error(FATAL,'idealized_ic_mod: nml "z0"= 0.0, but should be nonzero')
          t_value = t0*(1.0-tanh((zt(k)-thk)/z0)) + t1*(1.0-zt(k)/zt(nk))
          temp(:,:,k) = t_value
          write (stdout(),*) ' k=',k,' T=',t_value,' deg C.'
       enddo
    case('exponential_temp_ic')
       write (stdout(),'(/a/)') ' Note: exponential temperature IC is being set.'
       do k=1,nk
          t_value = 22.0*exp(-zt(k)/efold_depth)
          temp(:,:,k) = t_value
          write (stdout(),*) ' k=',k,' T=',t_value,' deg C.'
       enddo
    case('shelfbowl_temp_ic')
       write (stdout(),'(/a/)') ' Note: shelfbowl temperature ic being set, assuming south-edge of shelf is at 70N.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                t_value = 5.0*(1.0 + 0.5*(1.0-tanh(0.5*(yt(i,j)-70.0))) )
                temp(i,j,k)  = t_value
             enddo
          enddo
       enddo
    case('zonal_levitus_temp_dome_ic')
       write (stdout(),'(/a/)') ' Note: zonal averaged Levitus temp plus cold embayment IC set for depth < 2000m.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                temp(i,j,k)  = theta0 (yt(i,j), zt(k))
             enddo
          enddo
       enddo
       do j=jsc,jec
         do i=isc,iec
           if(yt(i,j) >= dome_embayment_south) then 
              temp(i,j,:) = constant_temp_value
           endif 
         enddo
       enddo
    case('temp_for_dome_ic')
       write (stdout(),'(/a/)') ' Note: linear stratification for DOME configuration with inflow condition.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                temp(i,j,k) = dome_conditions(xt(i,j), yt(i,j), zt(k), 'temp')
             enddo
          enddo
       enddo

    case('ideal_thermocline')
       write (stdout(),'(/a/)') ' Note: ideal_thermocline stratification specified for temp IC.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                temp(i,j,k) = ideal_thermocline(zt(k),zt(nk))
             enddo
          enddo
       enddo

    case('zonal_levitus_temp_ic')
       write (stdout(),'(/a/)') ' Note: zonal averaged Levitus temp IC set for depth < 2000m.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                temp(i,j,k)  = theta0 (yt(i,j), zt(k))
             enddo
          enddo
       enddo

    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(temp_type)//' is not a valid option of nml "temp_type" ')
    end select

    ! zero out temperature in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                temp(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

  end subroutine idealized_temp_ic

!#######################################################################
  ! idealized salinity initial conditions 
  subroutine idealized_salt_ic

    real    :: s_value, fact
    real    :: xwidth, ywidth, xwest, xeast, ysouth, ynorth
    integer :: i,j,k,kb

    select case (trim(salt_type))
    case('constant_salt_ic')
       s_value = constant_salt_value
       salt = s_value
       write (stdout(),'(/a,2(f6.2,a)/)') ' Note: constant salinity IC being set. s=',s_value,' psu.'
    case('exponential_salt_ic')
       write (stdout(),'(/a/)') ' Note: exponential salinity IC is being set.'
       do k=1,nk
          s_value = 22.0*exp(-zt(k)/efold_depth)
          salt(:,:,k) = s_value
          write (stdout(),*) ' k=',k,' S=',s_value, ' psu.'
       enddo
    case('salinity_profile_ic')
       write (stdout(),'(/a/)') ' Note: general salinity IC profile from salt0 is being set.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                salt(i,j,k)  = salt0 (zt(k))
             enddo
          enddo
       enddo
    case('salinity_for_dome_ic')
       write (stdout(),'(/a/)') ' Note: setting salinity according to DOME specification.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                salt(i,j,k) = dome_conditions(xt(i,j), yt(i,j), zt(k), 'salinity')
             enddo
          enddo
       enddo
    case('salinity_layer_ic')
        write (stdout(),'(/a/)') ' Note: setting salinity to 1.0 in a layer, 0.0 elsewhere.'
        salt(:,:,:) = 0.0
        do k=1,nk
           if(zt(k) >=salinity_layer_ztop .and. zt(k) <= salinity_layer_zbot) then 
               salt(:,:,k) = 1.0
           endif
        enddo
    case('salinity_patch_ic')
        write (stdout(),'(/a/)') ' Note: setting salinity to 1.0 in a patch, 0.0 elsewhere.'
        salt(:,:,:) = 0.0
        xwidth = abs(xt(ni,1)-xt(1,1))
        ywidth = abs(yt(1,nj)-yt(1,1)) 
        xwest  = xt(1,1) + salinity_patch_ratio1*xwidth
        xeast  = xt(1,1) + salinity_patch_ratio2*xwidth
        ysouth = yt(1,1) + salinity_patch_ratio1*ywidth
        ynorth = yt(1,1) + salinity_patch_ratio2*ywidth
        write(stdout(),*)' for salinity_patch_ic, xwest  = ', xwest
        write(stdout(),*)' for salinity_patch_ic, xeast  = ', xeast
        write(stdout(),*)' for salinity_patch_ic, ynorth = ', ynorth
        write(stdout(),*)' for salinity_patch_ic, ysouth = ', ysouth
        do k=1,nk
           if(zt(k) >=salinity_patch_ztop .and. zt(k) <= salinity_patch_zbot) then 
               do j=jsc,jec
                  do i=isc,iec
                      if(xt(i,j) >= xwest  .and. xt(i,j) <= xeast   .and. &
                         yt(i,j) >= ysouth .and. yt(i,j) <= ynorth ) then    
                         salt(i,j,k) = 1.0
                      endif 
                  enddo
               enddo
           endif
        enddo
    case('shelfbowl_salinity_ic')
       write (stdout(),'(/a/)') ' Note: shelfbowl salinity ic being set, assuming south-edge of shelf is at 70N.'
       do k=1,nk
          do j=jsc,jec
             do i=isc,iec
                if(yt(i,j)>=69.5) then 
                  s_value=constant_salt_value
                else
                  s_value=0.0
                endif
                salt(i,j,k) = s_value
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(salt_type)//' is not a valid option of nml "salt_type" ')
    end select

    ! zero out salinity in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                salt(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

end  subroutine idealized_salt_ic


!#######################################################################
  ! idealized age initial conditions 
  subroutine idealized_age_ic

    real    :: age_value, fact
    integer :: i,j,k,kb

    select case (trim(age_type))
    case('constant_age_ic')
       age_value = constant_age_value
       age = age_value
       write (stdout(),'(/a,2(f6.2,a)/)') ' Note: constant age IC being set. age=',age_value,' seconds.'
    case default
       call mpp_error(FATAL,'idealized_ic_mod: '//trim(age_type)//' is not a valid option of nml "age_type" ')
    end select

    ! zero out age in land points
    do j=jsc,jec
       do i=isc,iec
          kb = kmt(i,j)
          do k=1,nk
             if (k .gt. kb) then
                age(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

    return

end  subroutine idealized_age_ic

!#######################################################################
! idealized zonal velocity initial conditions (nothing yet here)
subroutine idealized_zonal_velocity_ic

  integer :: i,j,k

  write (stdout(),'(/a/)') ' Note: setting zonal velocity to zero.'
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           ut(i,j,k) = 0.0
        enddo
     enddo
  enddo

  return

end  subroutine idealized_zonal_velocity_ic



!#######################################################################
  ! idealized meridional velocity initial conditions 
subroutine idealized_merid_velocity_ic

  integer :: i,j,k

  write (stdout(),'(/a/)') ' Note: setting merid velocity according to DOME inflow conditions.'
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           vt(i,j,k) = dome_conditions(xt(i,j), yt(i,j), zt(k), 'vt')
        enddo
     enddo
  enddo

  return

end  subroutine idealized_merid_velocity_ic



!#######################################################################
  !--- read grid from grid_spec file.
  subroutine get_grid

    integer, dimension(2) :: layout  = (/ 1, 0 /)
    integer :: unit, npes, ndim, nvar, natt, ntime, len, i, j, k, kb, mk
    type(axistype), dimension(:), allocatable :: axes
    real, dimension(:,:), allocatable :: tmp, ht
    character(len=128) :: name

    !--- get grids information -------------------------------------------
    npes = mpp_npes()
    if(file_exist(trim(grid_file))) then
       call mpp_open(unit,trim(grid_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
       call mpp_get_info(unit, ndim, nvar, natt, ntime)
       allocate(axes(ndim))
       call mpp_get_axes(unit,axes)
       do i=1, ndim
          call mpp_get_atts(axes(i),name=name,len=len)
          select case(trim(name))
          case ('grid_x_T')
             ni = len
             allocate(grid_x_t(ni))
             call mpp_get_axis_data(axes(i),data=grid_x_t)
          case ('grid_y_T')
             nj = len
             allocate(grid_y_t(nj))
             call mpp_get_axis_data(axes(i),data=grid_y_t)        
          case ( 'zt' )
             nk = len
             allocate(zt(nk))
             call mpp_get_axis_data(axes(i),data=zt ) 
          end select
       enddo
       call mpp_close(unit)
       if(ni == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis grid_x_T')
       if(nj == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis grid_y_T')
       if(nk == 0) call mpp_error(FATAL,'idealized_ic_mod: '//trim(grid_file)//' does not contain axis zt')

       !--- define domain ------------------------------------------------
       call mpp_define_layout((/1,ni,1,nj/),npes,layout)
       call mpp_define_domains((/1,ni,1,nj/),layout, domain)
       call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
       allocate(xt(ni,nj), yt(ni,nj), kmt(isc:iec,jsc:jec),  &
            tmp(isc:iec,jsc:jec), ht(isc:iec,jsc:jec) )
       call read_data(trim(grid_file), 'x_T', xt)
       call read_data(trim(grid_file), 'y_T', yt)
       call read_data(trim(grid_file), 'depth_t', ht, domain)
       call read_data(trim(grid_file), 'num_levels', tmp,domain)
       kmt = tmp
       deallocate(tmp, axes)
    else
       call mpp_error(FATAL, 'idealized_ic_mod: file '//trim(grid_file)//' does not exist in your work directory')
    endif

  end subroutine get_grid

  !#######################################################################
  ! <SUBROUTINE NAME="write_idealized_ic_data">
  !
  ! <DESCRIPTION>
  !    Write out tracer data to netcdf file 
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine write_idealized_ic_data

    real, dimension(:,:,:,:),         allocatable :: tracer
    character(len=128), dimension(:), allocatable :: fld_name, fld_unit, fld_longname
    integer                                       :: n

    if(num_active_tracer .gt. 0) then
       allocate(tracer(ni,nj,nk,num_active_tracer), fld_name(num_active_tracer) )
       allocate( fld_unit(num_active_tracer), fld_longname(num_active_tracer) )
       n = 0
       if(generate_temp_ic) then
          n = n + 1
          tracer(:,:,:,n) = global_temp
          fld_name(n)     = 'temp'
          fld_unit(n)     = 'deg_c'
          fld_longname(n) = 'initial potential temp'
       endif
       if(generate_salt_ic) then
          n = n + 1
          tracer(:,:,:,n) = global_salt
          fld_name(n)     = 'salt'
          fld_unit(n)     = 'psu'
          fld_longname(n) = 'initial salinity'
       endif
       call write_tracer_data(active_tracer_file, tracer, fld_name, fld_unit, fld_longname )
       deallocate(tracer, fld_name, fld_unit, fld_longname)
    endif

    if(num_passive_tracer .gt. 0) then
       allocate(tracer(ni,nj,nk,num_passive_tracer), fld_name(num_passive_tracer) )
       allocate(fld_unit(num_passive_tracer), fld_longname(num_passive_tracer) )
       n = 0
       if(generate_age_ic) then
          n = n + 1
         tracer(:,:,:,n) = global_age
          fld_name(n)     = 'age_global'
          fld_unit(n)     = 'sec'
          fld_longname(n) = 'initial age'
       endif
       call write_tracer_data(passive_tracer_file, tracer, fld_name, fld_unit, fld_longname )
       deallocate(tracer, fld_name, fld_unit, fld_longname)
    endif

    if(generate_zonal_velocity_ic) then
       call write_velocity_data(zonal_velocity_file, global_ut, 'ut', &
                                'm/s','zonal velocity component at t-cell')
    endif
    if(generate_merid_velocity_ic) then
       call write_velocity_data(merid_velocity_file, global_vt, 'vt', &
                                'm/s','merid velocity component at t-cell')
    endif


    return

  end subroutine write_idealized_ic_data


  !#######################################################################
  !--- write velocity data to netcdf file, including writing out meta.

  subroutine write_velocity_data(file, velocity, fld_name, fld_unit, fld_longname)

    character(len=*),        intent(in)    :: file
    real, dimension(:,:,:),  intent(inout) :: velocity
    character(len=*),        intent(in)    :: fld_name, fld_unit, fld_longname

    integer          :: unit
    type(axistype)   :: axis_xt, axis_yt, axis_zt, axis_t
    type(fieldtype)  :: fld_x_T, fld_y_T
    type(fieldtype)  :: fld_velocity

    call mpp_open(unit,trim(file), MPP_OVERWR, MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)

    !--- write out meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_xt,'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = grid_x_t )
    call mpp_write_meta(unit, axis_yt,'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = grid_y_t )
    call mpp_write_meta(unit,axis_zt,'zt','meters','zt',&
         data=zt,cartesian='Z',sense=-1)
    call mpp_write_meta(unit, fld_x_T, (/axis_xt, axis_yt/), 'x_T', 'degree_east',  &
         'Geographic longitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_T, (/axis_xt, axis_yt/), 'y_T', 'degree_north',  &
         'Geographic latitude of T_cell centers', pack=1)

    call mpp_write_meta(unit,axis_t,'TIME','days since 0001-01-01 00:00:00','TIME',cartesian='T' )

    if(write_time_axis) then 
           call mpp_write_meta(unit, fld_velocity, (/axis_xt, axis_yt, axis_zt, axis_t/), &
                               trim(fld_name), trim(fld_unit), trim(fld_longname), pack=1)
    else 
           call mpp_write_meta(unit, fld_velocity, (/axis_xt, axis_yt, axis_zt/), &
                               trim(fld_name), trim(fld_unit), trim(fld_longname), pack=1)
    endif

    !--- write out data --------------------------------------------------
    call mpp_write(unit, axis_xt)
    call mpp_write(unit, axis_yt)
    call mpp_write(unit, axis_zt)
    call mpp_write(unit, fld_x_T, xt)
    call mpp_write(unit, fld_y_T, yt)
    call mpp_write(unit, fld_velocity, velocity(:,:,:) )
    call mpp_close(unit)

  end subroutine write_velocity_data



  !#######################################################################
  !--- write tracer data to netcdf file, including writing out meta.

  subroutine write_tracer_data(file, tracer, fld_name, fld_unit, fld_longname)

    character(len=*),               intent(in)   :: file
    real, dimension(:,:,:,:),    intent(inout)   :: tracer
    character(len=*), dimension(:), intent(in)   :: fld_name, fld_unit, fld_longname

    integer                                      :: unit, num_tracers, n
    type(axistype)                               :: axis_xt, axis_yt, axis_zt, axis_t
    type(fieldtype)                              :: fld_x_T, fld_y_T
    type(fieldtype),dimension(size(fld_name(:))) :: fld_tracer

    num_tracers = size(fld_name(:))

    call mpp_open(unit,trim(file), MPP_OVERWR, MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)

    !--- write out meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_xt,'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = grid_x_t )
    call mpp_write_meta(unit, axis_yt,'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = grid_y_t )
    call mpp_write_meta(unit,axis_zt,'zt','meters','zt',&
         data=zt,cartesian='Z',sense=-1)
    call mpp_write_meta(unit, fld_x_T, (/axis_xt, axis_yt/), 'x_T', 'degree_east',  &
         'Geographic longitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_T, (/axis_xt, axis_yt/), 'y_T', 'degree_north',  &
         'Geographic latitude of T_cell centers', pack=1)

    call mpp_write_meta(unit,axis_t,'TIME','days since 0001-01-01 00:00:00','TIME',cartesian='T' )

    if(write_time_axis) then 
        do n=1,num_tracers
           call mpp_write_meta(unit, fld_tracer(n), (/axis_xt, axis_yt, axis_zt, axis_t/), &
                               trim(fld_name(n)), trim(fld_unit(n)), trim(fld_longname(n)), pack=1)
        enddo
    else 
        do n=1,num_tracers
           call mpp_write_meta(unit, fld_tracer(n), (/axis_xt, axis_yt, axis_zt/), &
                               trim(fld_name(n)), trim(fld_unit(n)), trim(fld_longname(n)), pack=1)
        enddo
    endif

    !--- write out data --------------------------------------------------
    call mpp_write(unit, axis_xt)
    call mpp_write(unit, axis_yt)
    call mpp_write(unit, axis_zt)
    call mpp_write(unit, fld_x_T, xt)
    call mpp_write(unit, fld_y_T, yt)
    do n = 1, num_tracers    
       call mpp_write(unit, fld_tracer(n), tracer(:,:,:,n) )
    enddo
    call mpp_close(unit)

  end subroutine write_tracer_data

  !#######################################################################
 ! <SUBROUTINE NAME="idealized_ic_end">
  !
  ! <DESCRIPTION>
  !  Release memory.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine idealized_ic_end

    deallocate(temp, salt, age, xt, yt, kmt, zt, grid_x_t, grid_y_t)

    module_is_initialized = .false.
    return

  end subroutine idealized_ic_end


  !#######################################################################
  ! This subroutine returns estimates of global mean potential
  ! temperature for model initialization as a function of depth.
  ! it is used to produce a reference thermal stratification for the
  ! upper 2000m of the OCEAN`s test case.  below 2000m, the
  ! potential temperature returned is 2.0 degrees C.  surface
  ! values are set slightly above 18.4 degrees C at the reference
  ! latitude "reflat".
  ! the estimates are produced from a 7th order ploynomial fit to
  ! the annual mean world ocean potential temperature observations
  ! of Levitus (1982).
  !
  ! input [units]:
  !
  !   a latitdue (ydeg): [degrees] <BR/>
  !   a zt value (depth): [meters] <BR/>
  !
  ! output [units]:
  !
  !   potential temperature estimate (est): [degrees centigrade]
  !
  ! variables:
  !
  !   coeft     = coefficients for the polynomial fit of potential
  !               temperature vs. depth
  !
  !   reflat    = reference latitude at which observed surface
  !               temperatures approximately equal coeft(1)
  !
  !   factor    = the ratio of the cosine of the latitude requested
  !               ("ydeg") to the reference latitude ("reflat")
  !               used to scale the upper 2000 meters of the vertical
  !               temperature profile
  !
  !   tmin,tmax = the minumum and maximum potential temperatures
  !               allowed at the time of model initialization
  !
  ! reference:
  !   Levitus, S., Climatological atlas of the world ocean, NOAA
  ! Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.
  !
  ! Originally coded by Keith Dixon (Keith.Dixon)
  function theta0 (ydeg, depth)

    real, intent(in)   :: ydeg, depth
    integer, parameter :: ndeg=7
    real, save, dimension(ndeg+1) :: coeft = (/ 0.184231944E+02,-0.430306621E-01, 0.607121504E-04,-0.523806281E-07,&
                                                0.272989082E-10,-0.833224666E-14,0.136974583E-17,-0.935923382E-22 /)
    real, save :: tmin=2.0, tmax=25.0, reflat=34.0
    real       :: refcos, coslat, factor, z, est, bb, theta0, theta, t0, t1, h, z0
    integer    :: nn

    refcos = abs(cos(reflat/RADIAN))

    coslat = abs(cos(ydeg/RADIAN))
    factor = coslat/refcos

    if (depth > 2000.) then
       est = 2.0
    else
       est = 0.0
       bb = 1.0
       do nn=1,ndeg+1
          est = est + coeft(nn)*bb
          bb = bb*depth
       enddo
       est = est * factor
    endif

    if (est > tmax) est = tmax
    if (est < tmin) est = tmin

    theta0 = est

  end function theta0



  !#######################################################################
  ! This function returns ideal thermocline profile taken from 
  ! the appendix to Smith and Vallis (2001) JPO pages 554-571.  

  function ideal_thermocline(zdepth,zbottom)

    real, intent(in) :: zdepth
    real, intent(in) :: zbottom

    real :: ideal_thermocline
    real :: rho, delta_rho

    delta_rho = ideal_thermocline_alpha_eos*ideal_thermocline_deltaT

    rho = ideal_thermocline_rho0 &
    *(1.0 + (delta_rho/ideal_thermocline_rho0)*tanh(zdepth/zbottom/ideal_thermocline_scale_thick)**2) &
    *(1.0+ideal_thermocline_alpha*zdepth/zbottom) 
   
    ideal_thermocline = ideal_thermocline_offset + (ideal_thermocline_rho0-rho)/ideal_thermocline_alpha_eos

  
  end function ideal_thermocline
  !#######################################################################


  !#######################################################################
  ! This function returns initial and boundary conditions for DOME 
  ! as motivated by the work of Legg, Hallberg, and Girton (2006).
  ! It assumes a very coarse grid relative to the embayment size...
  ! e.g., four or less grid points within the embayment.  For this 
  ! case, we simply provide an exponential interface, with cold water
  ! beneath the interface with a southward velocity, and tagged with 
  ! salinity of one.  
  !
  ! For temperature, at the northern boundary of the embayment, we replace 
  ! a background temperature field with an inflow condition.
  ! 
  ! For salinity, it is given a value between zero and unity only on the
  ! northern boundary of the embayment, where the inflow water is 
  ! prescribed. 
  !
  ! For meridional velocity, there is a nonzero value only at the northern 
  ! row of the embayment with a specified profile. 
  ! 
  ! Picture at j=nj within the embayment:
  !
  !     ^            --------------------
  !     |            |                  |  ^
  !     |            |                  |  |
  !     |            |                  |  |dome_embayment_interface 
  !     |            |                  |  |
  !     |            |                  |  |
  !     |            |\    warm         |  v
  !     |            | \   reference    |  
  !     |            |  \               |
  !     |            |   \              |
  !dome_embayment_depth   \             |
  !     |            |     \<----exponential profile
  !     |            | cold \           |
  !     |            | inflow\          |
  !     |            |        \         |
  !     |            |         \        |
  !     |            |          \       |
  !     |            |           \      |
  !     |            |    X       \     |
  !     |            |             \    |
  !     v            |------------------|
  ! 
  !  A point in the inflow..X..has coordinates (inflow_x,inflow_z)
  !  with inflow_z the depth of the point, and inflow_x the distance
  !  from the western wall of the embayment. 
  ! 
  ! 


  function dome_conditions(xlon, ylat, zdepth, field_name)

    real,             intent(in) :: xlon
    real,             intent(in) :: ylat
    real,             intent(in) :: zdepth
    character(LEN=*), intent(in) :: field_name

    real :: dome_conditions
   
    real :: gprime   
    real :: deform 
    real :: vt_max 
    real :: interface_x 
    real :: interface_z 
    real :: inflow_x 
    real :: embayment_center 
    real :: tref
    real :: talpha
    real :: dtemp
    real :: temp_dome 
    real :: salinity_dome
    real :: vt_dome 

    ! values for reference background state 
    tref          = constant_temp_value - linear_theta_strat*zdepth
    temp_dome     = tref 
    salinity_dome = 0.0
    vt_dome       = 0.0

    ! compute properties at northern end of embayment 
    if(ylat  >=dome_embayment_north .and. &
       xlon  >=dome_embayment_west  .and. &
       xlon  <=dome_embayment_east) then 
!!$       xlon  <=dome_embayment_east  .and. &
!!$       zdepth<=dome_embayment_depth) then 

        ! reduced gravity between reference fluid and inflow fluid.
        ! assume reference fluid has density = rho0 
        gprime = dome_drho_inflow*rho0r*grav 

        ! maximum of inflow speed as determined by gprime and thickness of inflow
        vt_max = sqrt(gprime*dome_embayment_interface)

        ! deformation radius as determined by gprime, thickness of inflow, and coriolis 
        deform = vt_max/dome_embayment_coriolis

        ! distance in meters from the western wall of the embayment
        inflow_x = deg2metre*(xlon-dome_embayment_west)*cos(ylat)

        ! center of embayment 
        embayment_center =  dome_embayment_west + 0.5*(dome_embayment_east-dome_embayment_west) 
  
        ! depth of interface between reference and inflow water 
        interface_z = dome_embayment_interface + &
                      (dome_embayment_depth-dome_embayment_interface)*(1.0-exp(-inflow_x/deform))

        ! temperature difference between inflow (colder) and reference (warmer) 
        talpha = 2.0e-4
        dtemp  = dome_drho_inflow*rho0r/talpha

        ! values of the inflow along northern wall  
        if(zdepth >= interface_z .and. xlon <= embayment_center) then 
            temp_dome     = tref-dtemp
            salinity_dome = 1.0
            vt_dome       = -vt_max*exp(-inflow_x/deform)
        else 
            temp_dome     = tref  
            salinity_dome = 0.0
            vt_dome       = 0.0 
        endif

        if(dome_debug) then 
            write(*,*)' ' 
            write(*,'(a,f12.3,a,f12.3,a,f12.3,a)') &
                 ' in embayment with (x,y,z)     = (',xlon, ',', ylat, ',',zdepth,')' 
            write(*,'(a,f12.3,a,f12.3,a)') &
                 ' (inflow_x(km),interface_z(m)) = (',inflow_x/1000.0,',',interface_z,')'
            write(*,'(a,f12.3,a,f12.3,a,f12.3,a)') &
                 ' (temp,salinity,vt)            = (',temp_dome,',',salinity_dome,',',vt_dome,')'
        endif

    endif

    if(field_name==    'temp')     then 
       dome_conditions=temp_dome
    elseif(field_name=='salinity') then
       dome_conditions=salinity_dome
    elseif(field_name=='vt' )      then 
       dome_conditions=vt_dome
    endif     

  

  end function dome_conditions
  !#######################################################################





  !#######################################################################
  ! This function returns a salinity in psu
  function salt0 (depth)

    real, intent(in) :: depth
    real :: salt0

    salt0 = 32.0 + .001*depth 

  end function salt0
  !#######################################################################

end module idealized_ic_mod

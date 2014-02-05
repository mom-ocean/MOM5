module idealized_bc_mod 
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  ! This file is a part of MOM.                                                                 
  !                                                                      
  ! MOM is free software; you can redistribute it and/or modify it and  
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
  !!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Z. Liang  </CONTACT>
  !<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies </REVIEWER>

  !<OVERVIEW>
  ! For preparing idealized mom4 surface boundary conditions
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This program prepares idealized surface boundary conditions. The generated data
  ! contains some or all of the following, temperature, salinity, water flux,
  ! zonal and meridinal wind stress. The choice of output data is controlled by 
  ! namelist option.  Results are netcdf files that are then read into mom4.
  ! Various idealized options are available as selected by a namelist.
  !</DESCRIPTION>

  use fms_mod,         only : write_version_number, stdout, open_namelist_file
  use fms_mod,         only : file_exist, check_nml_error,  read_data
  use mpp_domains_mod, only : mpp_global_sum, domain2d, mpp_define_layout, mpp_global_field
  use mpp_domains_mod, only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_update_domains, BITWISE_EXACT_SUM, CYCLIC_GLOBAL_DOMAIN
  use mpp_io_mod,      only : mpp_open, mpp_close, axistype, fieldtype, MPP_RDONLY,MPP_NETCDF
  use mpp_io_mod,      only : MPP_MULTI, MPP_SINGLE, MPP_OVERWR, mpp_get_atts,  mpp_get_axis_data
  use mpp_io_mod,      only : mpp_get_info, mpp_get_axes, mpp_write_meta, mpp_write
  use mpp_mod,         only : mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, NOTE
  use mpp_mod,         only : mpp_min, mpp_max, mpp_chksum
  use constants_mod,   only : PI, RADIUS, EPSLN, RADIAN

  implicit none
  private

  !--- namelist ----------------------------------------------------------
  !<NAMELIST NAME="idealized_bc_nml">
  !  <DATA NAME="wind_type" TYPE="character(len=128)" >
  !  Control the iealized boundary wind stress condition options. 
  !  There are four options available and the default value is 
  !  "constant_tau". When temp_type is
  !  1. = "constant_tau", use space-time constant wind stress.
  !  2. = "cosine_zonal_winds", compute idealized winds using a cosine in latitude profile.  
  !        Makes sense only when running spherical coordinate model.
  !  3. = "frank_bryan_winds", compute idealized surface wind stress using zonal wind 
  !       profile originally  used by Frank Bryan. Makes sense only when running spherical coordinate model.
  !  4. = "frank_bryan_winds_compress", use the full Frank Bryan wind profiled over a 
  !        latitudinally compressed domain. Makes sense only when running spherical coordinate model. 
  !  </DATA>
  !  <DATA NAME="taux0"  TYPE="real">
  !  Constant zonal wind stress 
  !  </DATA>
  !  <DATA NAME="tauy0"  TYPE="real">
  !  Constant meridional wind stress 
  !  </DATA>
  !  <DATA NAME="qw0" UNITS="meter/sec" TYPE="real">
  !  Fresh water flux scaling parameter for idealized surface water flux 
  !  </DATA>
  !  <DATA NAME="generate_wind_bc" TYPE="logical">
  !   Control if wind stress data will be generated. When true (default value), idealized surface boundary
  !   wind stress data will be generated. If false, do not generate wind stress data. 
  !  </DATA>
  !  <DATA NAME="generate_temp_bc" TYPE="logical">
  !   Control if temperature data will be generated. When true (default value), idealized surface boundary
  !   temperature data will be generated. If false, do not generate temperature data. 
  !  </DATA>
  !  <DATA NAME="generate_salt_bc" TYPE="logical">
  !   Control if salinity data will be generated. When true (default value), idealized surface boundary
  !   salinity data will be generated. If false, do not generate salinity data. 
  !  </DATA>
  !  <DATA NAME="generate_water_bc" TYPE="logical">
  !   Control if water flux data will be generated. When true (default value), idealized surface boundary
  !   water flux data will be generated. If false, do not generate water flux data. 
  !  </DATA>
  ! </NAMELIST>
  character(len=128) :: wind_type = 'constant_tau'
  character(len=128) :: temp_type = 'idealized_temp_restore'
  character(len=128) :: salt_type = 'idealized_salt_restore'
  character(len=128) :: water_type = 'idealized_water_flux'
  logical            :: cyclic = .true.
  logical            :: generate_wind_bc = .true.
  logical            :: generate_temp_bc = .true.
  logical            :: generate_salt_bc = .true.
  logical            :: generate_water_bc = .true.
  real               :: qw0=-3.17e-8  ! 1m/year fresh water velocity (converted to m/s) for idealized
  real               :: taux0=0.0, tauy0=0.0  ! constant zonal/meridional windstress (newton/m^2) for idealized
  namelist /idealized_bc_nml/ wind_type, taux0, tauy0, qw0,  generate_wind_bc, &
       generate_temp_bc, generate_salt_bc, generate_water_bc, cyclic

  !--- output data variables -------------------------------------------
  real, dimension(:,:), allocatable :: taux, tauy, temp, salt, water
  real, dimension(:,:), allocatable :: global_taux, global_tauy, global_temp, global_salt, global_water
  !--- output file name
  character(len=128) :: temp_file  = 'temp_sfc_restore.nc'
  character(len=128) :: salt_file  = 'salt_sfc_restore.nc'
  character(len=128) :: wind_file  = 'tau.nc'
  character(len=128) :: water_file = 'water_flux.nc'

  !--- version information 
  character(len=128) :: version = '$Id: idealized_bc.f90,v 17.0 2009/07/21 03:22:49 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 

  !--- other module variables
  real, dimension(:,:),    allocatable :: xu, yu, xt, yt, h1t, area_t
  real, dimension(:,:),    allocatable :: umask, tmask
  integer, dimension(:,:), allocatable :: kmt, kmu
  real, dimension(:),      allocatable :: grid_x_t, grid_y_t, grid_x_c, grid_y_c
  integer                              :: isc, iec, jsc, jec, ni, nj
  real,                      parameter :: c2mks = 0.1   ! converts dynes/cm^2 to newton/m^2 
  type(domain2d),                 save :: domain
  character(len=128)                   :: grid_file = 'INPUT/grid_spec.nc'
  logical                              :: module_is_initialized = .FALSE.

  !--- public interface
  public :: idealized_bc_init, write_idealized_bc_data, idealized_bc_end

contains

  !#######################################################################
  ! <SUBROUTINE NAME="idealized_bc_init">
  !
  ! <DESCRIPTION>
  ! Initialize the module generating ideal surface boundary conditions.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine idealized_bc_init

    !--- local variables -----------------------------------------------
    integer :: unit, ierr, io_status, chksum_temp, chksum_salt
    integer :: chksum_water, chksum_taux, chksum_tauy
    !-------------------------------------------------------------------
    if ( module_is_initialized ) call mpp_error(FATAL, 'idealized_bc_mod: attempted reinitialization')
    module_is_initialized = .TRUE.

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1
       do while (ierr /= 0)
          read  (unit, idealized_bc_nml,iostat=io_status, end = 10)
          ierr = check_nml_error(io_status,'idealized_bc_nml')
       enddo
10     call mpp_close (unit)
    else
       call mpp_error(NOTE,'idealized_bc_mod: file input.nml does not exist, use default namelist option')  
    endif

    !--- at least one of generate_temp_bc, generate_salt_bc, generate_water_bc and
    !--- generate_wind_bc should be true
    if(.not.generate_temp_bc .and. .not.generate_salt_bc .and. .not.generate_wind_bc .and. .not.generate_water_bc)  &
         call mpp_error(FATAL, 'idealized_bc_mod: at least one of nml generate_temp_bc, generate_salt_bc,'//  &
                         ' generate_wind_bc, generate_water_bc should be true')

    !--- write out version information and namelist option ---------------
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
       write (stdout(), nml= idealized_bc_nml)
    endif

    !--- get grid infomration
    call get_grid

    !--- allocate memory to module variables
    allocate(temp(isc:iec,jsc:jec), salt(isc:iec,jsc:jec),  &
         taux(isc:iec,jsc:jec), tauy(isc:iec,jsc:jec), water( isc:iec,jsc:jec) )
    allocate(global_temp(ni,nj), global_salt(ni,nj), global_taux(ni,nj), &
         global_tauy(ni,nj), global_water(ni,nj) )

    ! idealized wind boundary conditions 
    if(generate_wind_bc) then
       call idealized_wind_bc()
       call mpp_global_field(domain, taux, global_taux)
       call mpp_global_field(domain, tauy, global_tauy)
       chksum_taux = mpp_chksum(taux)
       chksum_tauy = mpp_chksum(tauy)
       write (stdout(),'(/(a,I20/))')' Idealized_bc condition Checksum: taux=', chksum_taux
       write (stdout(),'(/(a,I20/))')' Idealized_bc condition Checksum: tauy=', chksum_tauy
    endif

    ! idealized temperature boundary conditions
    if(generate_temp_bc) then
       call idealized_temp_bc()
       call mpp_global_field(domain, temp, global_temp)
       chksum_temp = mpp_chksum(temp) 
       write (stdout(),'(/(a,I20/))')' Boundary condition Checksum: temp=', chksum_temp
    endif
    ! idealized salinity boundary conditions
    if(generate_salt_bc) then
       call idealized_salt_bc()
       call mpp_global_field(domain, salt, global_salt)
       chksum_salt = mpp_chksum(salt)
       write (stdout(),'(/(a,I20/))')' Boundary condition Checksum: salt=', chksum_salt
    endif
    ! idealized water flux boundary conditions
    if(generate_water_bc) then
       call idealized_water_bc()
       call mpp_global_field(domain, water, global_water)
       chksum_water = mpp_chksum(water)
       write (stdout(),'(/(a,I20/))')' Boundary condition Checksum: water=', chksum_water
    endif

    return

  end subroutine idealized_bc_init

  !#######################################################################
  !--- idealized wind stress boundary condition
  subroutine idealized_wind_bc( )

    real    :: yunj, yu1, y_lat, arg_lat, frank_scale
    integer :: i, j

    select case(trim(wind_type))
    case ('constant_tau')
       write(stdout(),*)'idealized_wind_stress being used with constant_tau.'
       taux(:,:) = taux0*umask(:,:)
       tauy(:,:) = tauy0*umask(:,:)
    case ('cosine_zonal_winds')
       write(stdout(),*)'idealized_wind_stress being used with cosine_zonal_winds.'
       yunj = yu(isc,jec)
       yu1  = yu(isc,jsc)
       do j=jsc,jec
          do i=isc,iec
             yu1  = min(yu(i,j),yu1)
             yunj = max(yu(i,j),yunj)
          enddo
       enddo
       call mpp_min (yu1)
       call mpp_max (yunj)
       do j=jsc,jec
          do i=isc,iec
             arg_lat   = (yu(i,j)-yu1)/(yunj-yu1)*2.0*PI
             taux(i,j) = -c2mks*cos(arg_lat)*umask(i,j)
             tauy(i,j) = 0.0 
          enddo
       enddo
    case ('frank_bryan_winds')
       write(stdout(),*)'idealized_wind_stress being used with frank_bryan_winds.'
       do j=jsc,jec
          do i=isc,iec
             y_lat     = abs(yu(i,j))/RADIAN
             taux(i,j) = c2mks*(0.8*(-sin(6.0*y_lat)-1.0)+0.5*(tanh(5.0*PI-10.0*y_lat)+tanh(10.0*y_lat)))
             taux(i,j) = taux(i,j)*umask(i,j)
             tauy(i,j) = 0.0 
          enddo
       enddo
    case ('frank_bryan_winds_compress')
       write(stdout(),*)'idealized_wind_stress being used with frank_bryan_winds_compress.'
       do j=jsc,jec
          do i=isc,iec
             frank_scale = 90.0/yu(i,nj)
             y_lat      = frank_scale*abs(yu(i,j))/RADIAN
             taux(i,j) = c2mks*(0.8*(-sin(6.0*y_lat)-1.0)+0.5*(tanh(5.0*PI-10.0*y_lat)+tanh(10.0*y_lat)))
             taux(i,j) = taux(i,j)*umask(i,j)
             tauy(i,j) = 0.0 
          enddo
       enddo
    case default
       call mpp_error(FATAL,'idealized_bc_mod: '//trim(wind_type)//' is not a valid option of nml " wind_type" ')
    end select

  end subroutine idealized_wind_bc

  !#######################################################################
    !--- idealized temperature boundary condition
  subroutine idealized_temp_bc

    select case(trim(temp_type))
    case('idealized_temp_restore')
       write(stdout(),*)'idealized_temp_restore gives an idealized profile for restoring surface temperature.'
       temp(isc:iec,jsc:jec) = 26.93 - 0.3618*abs(yt(isc:iec,jsc:jec)) 
    case default
       call mpp_error(FATAL,'idealized_bc_mod: '//trim(temp_type)//' is not a valid option of nml "temp_type"')
    end select

  end subroutine idealized_temp_bc

  !#######################################################################
  !--- idealized profile for restoring surface salinity
  subroutine idealized_salt_bc

    !--- idealized salinity
    select case(trim(salt_type))
    case('idealized_salt_restore')
       write(stdout(),*)'idealized_salt_restore gives an idealized profile for restoring surface salinity.'
       salt(isc:iec,jsc:jec) = 24.+12.*sin(2.*abs(yt(isc:iec,jsc:jec))/RADIAN)
    case default
       call mpp_error(FATAL,'idealized_bc_mod: '//trim(salt_type)//' is not a valid option of nml salt_type')
    end select

  end subroutine idealized_salt_bc

  !#######################################################################
  !--- idealized water flux
  subroutine idealized_water_bc

    real    :: yt1, ytnj, surf_area, total_sfwf
    integer :: i, j     

    select case(trim(water_type))
    case('idealized_water_flux')
       write(stdout(),*)'idealized_water_flux gives an idealized profile for surface fresh water flux.'
       ! bounding latitudes for ideal water flux
       ytnj = yt(isc,jsc); yt1 = yt(isc,jsc)
       do j=jsc,jec
          do i=isc,iec
             yt1  = min(yt(i,j),yt1)
             ytnj = max(yt(i,j),ytnj)
          enddo
       enddo
       call mpp_min (yt1)
       call mpp_max (ytnj)

       surf_area=0.0 ; total_sfwf=0.0
       water(:,:) = tmask(:,:)*(qw0*RADIUS/h1t(:,:))*(1.0-2.0*(yt(isc:iec,jsc:jec)-yt1)/(ytnj-yt1))
       total_sfwf = mpp_global_sum (domain, area_t(:,:)*water(:,:), BITWISE_EXACT_SUM)        
       surf_area  = mpp_global_sum (domain, area_t(:,:)*tmask(:,:),  BITWISE_EXACT_SUM) 
       total_sfwf = total_sfwf/(epsln+surf_area)
       water(:,:) = tmask(:,:)*(water(:,:)-total_sfwf)
    case default
       call mpp_error(FATAL,'idealized_bc_mod: '//trim(water_type)//' is not a valid option of nml water_type')
    end select

  end subroutine idealized_water_bc

  !#######################################################################
  ! <SUBROUTINE NAME="idealized_bc_end">
  !
  ! <DESCRIPTION>
  !  Release memory.
  ! </DESCRIPTION>
  ! </SUBROUTINE>

  subroutine idealized_bc_end

    deallocate(taux, tauy, temp, salt, water, kmt, kmu, xu, yu, xt, yt, h1t)
    deallocate(area_t, umask, tmask, grid_x_t, grid_y_t, grid_x_c, grid_y_c)

    module_is_initialized = .false.
    return

  end subroutine idealized_bc_end

  !#######################################################################
  !--- read grid from grid_spec file.
  subroutine get_grid

    integer,                     dimension(2) :: layout  = (/ 1, 0 /)
    integer                                   :: unit, npes, ndim, nvar, natt, ntime, len
    integer                                   :: mi, mj, i, j, k, kb, isd, ied, jsd, jed
    type(axistype), dimension(:), allocatable :: axes
    real, dimension(:,:),         allocatable :: tmp, ht
    real, dimension(:),           allocatable :: dzt, fracdz 
    real, dimension(:,:,:),       allocatable :: dht, ztp
    character(len=128)                        :: name

    !--- get grids information -------------------------------------------
    npes = mpp_npes()
    ni = 0; nj = 0
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
          case ('grid_x_C')
             mi = len
             allocate(grid_x_c(mi))
             call mpp_get_axis_data(axes(i),data=grid_x_c)
          case ('grid_y_C')
             mj = len
             allocate(grid_y_c(mj))
             call mpp_get_axis_data(axes(i),data=grid_y_c)        
          end select
       enddo
       call mpp_close(unit)
       if(ni == 0) call mpp_error(FATAL,'idealized_bc_mod: '//trim(grid_file)//' does not contain axis grid_x_T')
       if(nj == 0) call mpp_error(FATAL,'idealized_bc_mod: '//trim(grid_file)//' does not contain axis grid_y_T')
       if(mi .ne. ni .or. mj .ne. mj) call mpp_error(FATAL,'idealized_bc_mod: in '//trim(grid_file) &
            //' T-grid and C-grid should have the same size')

       !--- define domain ------------------------------------------------
       call mpp_define_layout((/1,ni,1,nj/),npes,layout)
       if(cyclic) then
          call mpp_define_domains((/1,ni,1,nj/),layout, domain, xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1)
       else
          call mpp_define_domains((/1,ni,1,nj/),layout, domain, xhalo=1, yhalo=1)
       endif
       call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
       call mpp_get_data_domain(domain, isd, ied, jsd, jed)
       allocate(xt(ni,nj), yt(ni,nj), xu(ni,nj),  &
            yu(ni,nj), kmt(isc:iec,jsc:jec), kmu(isc:iec,jsc:jec),   &
            h1t(isc:iec,jsc:jec), area_t(isc:iec,jsc:jec), tmask(isc:iec,jsc:jec), &
            umask(isc:iec,jsc:jec), tmp(isd:ied,jsd:jed) )
       tmp = 10000.0  ! dummy large number
       call read_data(trim(grid_file), 'x_T', xt)
       call read_data(trim(grid_file), 'y_T', yt)
       call read_data(trim(grid_file), 'x_C', xu)
       call read_data(trim(grid_file), 'y_C', yu )
       call read_data(trim(grid_file), 'num_levels', tmp(isc:iec,jsc:jec), domain)
       call read_data(trim(grid_file), 'area_T',area_t, domain)
       call mpp_update_domains(tmp,domain)
       kmt(isc:iec,jsc:jec) = tmp(isc:iec,jsc:jec)
       do j=jsc,jec
          do i=isc,iec
             kmu(i,j) = min(tmp(i,j), tmp(i+1,j), tmp(i,j+1), tmp(i+1,j+1))
          enddo
       enddo
       tmask = 0.0
       where(kmt .ge. 1) tmask = 1.0
       umask = 0.0
       where(kmu .ge. 1) umask = 1.0

       do j=jsc,jec
          do i=isc,iec
             if (cos(yt(i,j)/RADIAN) == 0.0) then
                h1t(i,j) = radius*abs(EPSLN)
             else
                h1t(i,j) = radius*cos(yt(i,j)/RADIAN)
             endif
          enddo
       enddo
    else
       call mpp_error(FATAL, 'idealized_bc_mod: file '//trim(grid_file)//' does not exist in your work directory')
    endif

    ! release memory
    deallocate(axes, tmp)

  end subroutine get_grid

  !#######################################################################
  ! <SUBROUTINE NAME="write_idealized_bc_data">
  !
  ! <DESCRIPTION>
  !    Write the idealized boundary condition to netcdf file.
  ! </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine write_idealized_bc_data

    if(generate_temp_bc) call write_flux_data(temp_file, global_temp, 'temp', 'deg_c', 'surface restoring temp')
    if(generate_salt_bc) call write_flux_data(salt_file, global_salt, 'salt', 'g/kg', 'surface restoring salinity')
    if(generate_water_bc) call write_flux_data(water_file, global_water, 'lprec', 'm/s', 'surface liquid water flux')
    if(generate_wind_bc) call write_wind_data()

    return

  end subroutine write_idealized_bc_data

  !#######################################################################
  !--- write tracer data to output file
  subroutine write_flux_data(file, flux, fld_name, fld_unit, fld_longname)

    character(len=*), intent(in)       :: file, fld_name, fld_unit, fld_longname
    real, dimension(:,:),intent(inout) :: flux

    integer         :: unit
    type(axistype)  :: axis_x, axis_y, axis_t
    type(fieldtype) :: fld_x_t, fld_y_t, fld_flux

    call mpp_open(unit,trim(file), MPP_OVERWR, MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)

    !--- write out meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_x,'grid_x_T','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = grid_x_t )
    call mpp_write_meta(unit, axis_y,'grid_y_T','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = grid_y_t )
    call mpp_write_meta(unit,axis_t,'TIME','days since 0001-01-01 00:00:00','TIME',cartesian='T' )
    call mpp_write_meta(unit, fld_x_t, (/axis_x, axis_y/), 'x_T', 'degree_east',  &
         'Geographic longitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_t, (/axis_x, axis_y/), 'y_T', 'degree_north',  &
         'Geographic latitude of T_cell centers', pack=1)
    call mpp_write_meta(unit, fld_flux, (/axis_x, axis_y, axis_t/), trim(fld_name), &
         trim(fld_unit), trim(fld_longname), pack=1)

    !--- write out data --------------------------------------------------
    call mpp_write(unit, axis_x)
    call mpp_write(unit, axis_y)
    call mpp_write(unit, fld_x_t, xt)
    call mpp_write(unit, fld_y_t, yt)
    call mpp_write(unit, fld_flux, flux, tstamp=1.)
    call mpp_close(unit)

  end subroutine write_flux_data
  !#######################################################################
  !--- write stress wind to the output file wind_file
  subroutine write_wind_data( )

    integer         :: unit
    type(axistype)  :: axis_x, axis_y, axis_t
    type(fieldtype) :: fld_x_c, fld_y_c, fld_taux, fld_tauy

    call mpp_open(unit,trim(wind_file), MPP_OVERWR, MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)

    !--- write out meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_x,'grid_x_C','degree_east','Nominal Longitude of C-cell center', &
         cartesian ='X', data = grid_x_c )
    call mpp_write_meta(unit, axis_y,'grid_y_C','degree_north','Nominal Latitude of C-cell center', &
         cartesian ='Y',  data = grid_y_c )
    call mpp_write_meta(unit,axis_t,'TIME','days since 0001-01-01 00:00:00','TIME', cartesian='T' )
    call mpp_write_meta(unit, fld_x_c, (/axis_x, axis_y/), 'x_C', 'degree_east',  &
         'Geographic longitude of C_cell centers', pack=1)
    call mpp_write_meta(unit, fld_y_c, (/axis_x, axis_y/), 'y_C', 'degree_north',  &
         'Geographic latitude of C_cell centers', pack=1)
    call mpp_write_meta(unit, fld_taux, (/axis_x, axis_y, axis_t/), 'taux', &
         'N/m^2', 'ZONAL WIND STRESS', pack=1)
    call mpp_write_meta(unit, fld_tauy, (/axis_x, axis_y, axis_t/), 'tauy', &
         'N/m^2', 'MERIDIONAL WIND STRESS', pack=1)
    !--- write out data --------------------------------------------------
    call mpp_write(unit, axis_x)
    call mpp_write(unit, axis_y)
    call mpp_write(unit, fld_x_c, xu)
    call mpp_write(unit, fld_y_c, yu)
    call mpp_write(unit, fld_taux,  global_taux, tstamp=1.)
    call mpp_write(unit, fld_tauy,  global_tauy, tstamp=1.)
    call mpp_close(unit)

  end subroutine write_wind_data

  !#######################################################################

end module idealized_bc_mod

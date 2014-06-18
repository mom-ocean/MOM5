program regrid_hlfs

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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT>
  !
  ! <OVERVIEW>
  !   This program is to calculate surface upward laten heat flux on the atmosphere grid.
  ! </OVERVIEW>

  !<DESCRIPTION>
  ! to compile, in the command line, type 
  !   "
  !</DESCRIPTION>

  implicit none
#include "netcdf.inc"

  real, parameter   :: epsln = 1.e-20
  integer           :: stdout = 6

  character(len=24) :: grid_file = "input/grid_spec.nc"
  character(len=24) :: fld_in_from_lnd = 'latent'
  character(len=24) :: fld_in_from_ice =  'LH'
  character(len=24) :: src_from_lnd,    src_from_ice
  character(len=24) :: fld_out_from_ice = 'LH_from_ice'
  character(len=24) :: fld_out_from_tot = 'sfc_latent_tot'
  character(len=24) :: output_from_ice, output_from_tot
  character(len=32) :: tunits = " "

  logical           :: time_bounds_exist
  integer           :: ni_ice, nj_ice, ntime
  integer           :: ni_atm, nj_atm
  integer           :: ncid_ice, ncid_tot
  integer           :: id_time_ice, id_timeb_ice
  integer           :: id_time_tot, id_timeb_tot
  integer           :: id_lh_ice, id_lh_tot
  integer           :: id_t1_ice, id_t2_ice, id_dt_ice
  integer           :: id_t1_tot, id_t2_tot, id_dt_tot
  real              :: missing_lnd, missing_ice


  real, dimension(:),    allocatable :: lon, lat, lonb, latb
  integer, dimension(:), allocatable :: i_atm_atmxocn, j_atm_atmxocn 
  integer, dimension(:), allocatable :: i_ocn_atmxocn, j_ocn_atmxocn
  real, dimension(:),    allocatable :: area_atmxocn
  real, dimension(:,:),  allocatable :: area_atm, area_ocn 
  real, dimension(:),    allocatable :: time_value, t1, t2, dt
  real, dimension(:,:),  allocatable :: time_bounds

  call regrid_hfls_init

  !---- read grid information
  call read_grid

  !--- get source data information
  call get_src_info

  !--- setup meta data for the output file
  call setup_meta

  !--- processing data and write out the data sets
  call process_data()


contains

  !#####################################################################
  subroutine regrid_hfls_init
    integer :: unit_begin, unit_end, unit
    logical :: opened, file_exist

    src_from_lnd = 'input/'//trim(fld_in_from_lnd)//'.nc'
    src_from_ice = 'input/'//trim(fld_in_from_ice)//'.nc'
    output_from_ice = trim(fld_out_from_ice)//'.nc'
    output_from_tot = trim(fld_out_from_tot)//'.nc'

  end subroutine regrid_hfls_init

  !#####################################################################
  subroutine read_grid

    integer :: rcode, ncid, ncells

    rcode = nf_open(grid_file, NF_NOWRITE, ncid)
    call error_handler('error in open file '//grid_file, rcode)


    !--- read the ocean grid information ( source grid )
    ni_ice = get_dimlen(ncid, 'gridlon_t')
    nj_ice = get_dimlen(ncid, 'gridlat_t')

    !--- read the atmos grid information ( destination grid )
    ni_atm = get_dimlen(ncid, 'xta')
    nj_atm = get_dimlen(ncid, 'yta')

    allocate(lon(ni_atm), lat(nj_atm), lonb(ni_atm+1), latb(nj_atm+1) )

    call get_var_real_1d(ncid, 'xta', lon)
    call get_var_real_1d(ncid, 'yta', lat)  
    call get_var_real_1d(ncid, 'xba', lonb)
    call get_var_real_1d(ncid, 'yba', latb)  

    !--- get exchange grid information
    ncells = get_dimlen(ncid, 'I_ATM_ATMxOCN')
    allocate(i_atm_atmxocn(ncells), j_atm_atmxocn(ncells)    )
    allocate(i_ocn_atmxocn(ncells), j_ocn_atmxocn(ncells)    )
    allocate(area_atmxocn(ncells),  area_atm(ni_atm, nj_atm) )

    call get_var_int_1d(ncid, 'I_ATM_ATMxOCN', i_atm_atmxocn)
    call get_var_int_1d(ncid, 'J_ATM_ATMxOCN',j_atm_atmxocn )
    call get_var_int_1d(ncid, 'I_OCN_ATMxOCN',i_ocn_atmxocn )
    call get_var_int_1d(ncid, 'J_OCN_ATMxOCN',j_ocn_atmxocn )
    call get_var_real_1d(ncid, 'AREA_ATMxOCN',area_atmxocn )
    call get_var_real_2d(ncid, 'AREA_ATM',area_atm )

    allocate(area_ocn(ni_atm,nj_atm) )

    rcode = nf_close(ncid)


  end subroutine read_grid

  !#####################################################################

  subroutine get_src_info

    integer :: rcode, ncid, varid

    rcode = nf_open(trim(src_from_lnd), NF_NOWRITE, ncid) 
    call error_handler('error in opening file '//trim(src_from_lnd), rcode) 

    ! check if time_bounds exist or not
    rcode = nf_inq_varid(ncid, 'time_bounds', varid)
    if(rcode == 0) then
       time_bounds_exist = .true.;
    else
       time_bounds_exist = .false.;
       write(stdout,*)" time average information does not exist, ", &
                   "There will be not time_bounds, average_T1, average_T2 ", &
                   "and average_DT in the output file"
    endif

    ntime = get_dimlen(ncid, 'time')
    allocate(time_value(ntime), time_bounds(2, ntime) )
    allocate(t1(ntime), t2(ntime), dt(ntime) )
    call get_var_real_1d(ncid, 'time', time_value)
    if(time_bounds_exist) then
       call get_var_real_1d(ncid, 'average_T1', t1)
       call get_var_real_1d(ncid, 'average_T2', t2)
       call get_var_real_1d(ncid, 'average_DT', dt)
       call get_var_real_2d(ncid, 'time_bounds', time_bounds)
    endif

    !--- get the time units
    call get_time_units(ncid)

    !--- get the missing value
    call get_missing_value(ncid, fld_in_from_lnd, missing_lnd )

    rcode = nf_close(ncid) 

    rcode = nf_open(trim(src_from_ice), NF_NOWRITE, ncid) 
    call error_handler('error in opening file '//trim(src_from_ice), rcode) 

    call get_missing_value(ncid, fld_in_from_ice, missing_ice )

  end subroutine get_src_info

  !#####################################################################
  subroutine setup_meta

    integer          :: rcode, dims(6), id_xt, id_yt, id_xb, id_yb, id_nv
    character(len=4) :: units

    if(tunits(1:4) == 'days') then
       units = 'days'
    else
       call error_handler('time units should be days, need to update the source code')
    endif

    !--- setup for the part from ice model
    rcode = nf_create(trim(output_from_ice),NF_WRITE, ncid_ice)
    call error_handler('error in creating file '//trim(output_from_ice), rcode) 

    !--- define dimension
    dims(1) = define_dim(ncid_ice, 'xt',ni_atm)
    dims(2) = define_dim(ncid_ice, 'yt',nj_atm)
    dims(3) = define_dim(ncid_ice, 'time', NF_UNLIMITED)
    dims(4) = define_dim(ncid_ice, 'xb', ni_atm+1)
    dims(5) = define_dim(ncid_ice, 'yb', nj_atm+1)
    dims(6) = define_dim(ncid_ice, 'nv', 2)

    !--- define variable
    id_xt        = define_var_1d(ncid_ice, 'xt', NF_FLOAT, dims(1) )
    call put_att_text(ncid_ice, id_xt, 'long_name', 'longitude')
    call put_att_text(ncid_ice, id_xt, 'units', 'degrees_E')
    call put_att_text(ncid_ice, id_xt, 'cartesian_axis', 'X')
    call put_att_text(ncid_ice, id_xt, 'edges', 'xb')

    id_yt        = define_var_1d(ncid_ice, 'yt', NF_FLOAT, dims(2) )
    call put_att_text(ncid_ice, id_yt, 'long_name', 'latitude')
    call put_att_text(ncid_ice, id_yt, 'units', 'degrees_N')
    call put_att_text(ncid_ice, id_yt, 'cartesian_axis', 'Y')
    call put_att_text(ncid_ice, id_yt, 'edges', 'yb')

    id_xb        = define_var_1d(ncid_ice, 'xb', NF_FLOAT, dims(4) )
    call put_att_text(ncid_ice, id_xb, 'long_name', 'longitude edges')
    call put_att_text(ncid_ice, id_xb, 'units', 'degrees_E')
    call put_att_text(ncid_ice, id_xb, 'cartesian_axis', 'X')

    id_yb        = define_var_1d(ncid_ice, 'yb', NF_FLOAT, dims(5) )
    call put_att_text(ncid_ice, id_yb, 'long_name', 'latitude edges')
    call put_att_text(ncid_ice, id_yb, 'units', 'degrees_N')
    call put_att_text(ncid_ice, id_yb, 'cartesian_axis', 'Y')

    id_time_ice  = define_var_1d(ncid_ice, 'time', NF_DOUBLE, dims(3) )
    call put_att_text(ncid_ice, id_time_ice, 'long_name', 'time')
    call put_att_text(ncid_ice, id_time_ice, 'units', trim(tunits) )
    call put_att_text(ncid_ice, id_time_ice, 'cartesian_axis', 'T')
    call put_att_text(ncid_ice, id_time_ice, 'calendar_type', 'NOLEAP')
    call put_att_text(ncid_ice, id_time_ice, 'calendar', 'NOLEAP')
    if(time_bounds_exist) then
       call put_att_text(ncid_ice, id_time_ice, 'bounds', 'time_bounds')
    endif

    id_nv        = define_var_1d(ncid_ice, 'nv', NF_FLOAT, dims(6) )
    call put_att_text(ncid_ice, id_nv, 'long_name', 'vertex number')
    call put_att_text(ncid_ice, id_nv, 'units', 'none')
    call put_att_text(ncid_ice, id_nv, 'cartesian_axis', 'N')

    if(time_bounds_exist) then
       id_t1_ice    = define_var_1d(ncid_ice, 'average_T1', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_ice, id_t1_ice, 'long_name', 'Start time for average period')
       call put_att_text(ncid_ice, id_t1_ice, 'units', trim(tunits) )

       id_t2_ice    = define_var_1d(ncid_ice, 'average_T2', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_ice, id_t2_ice, 'long_name', 'End time for average period')
       call put_att_text(ncid_ice, id_t2_ice, 'units',  trim(tunits) )

       id_dt_ice    = define_var_1d(ncid_ice, 'average_DT', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_ice, id_dt_ice, 'long_name', 'Length of average period')
       call put_att_text(ncid_ice, id_dt_ice, 'units', trim(units))

       id_timeb_ice = define_var_2d(ncid_ice, 'time_bounds', NF_DOUBLE, (/dims(6), dims(3)/) )
       call put_att_text(ncid_ice, id_timeb_ice, 'long_name', 'time axis boundaries')
       call put_att_text(ncid_ice, id_timeb_ice, 'units', trim(units) )
    endif

    id_lh_ice    = define_var_3d(ncid_ice, trim(fld_out_from_ice), NF_FLOAT, dims(1:3) )
    call put_att_text(ncid_ice, id_lh_ice, 'long_name', 'latent heat flux')
    call put_att_text(ncid_ice, id_lh_ice, 'units', 'W/m^2')
    call put_att_real(ncid_ice, id_lh_ice, 'missing_value', missing_ice)
    call put_att_text(ncid_ice, id_lh_ice, 'cell_methods', 'time: mean')

    if(time_bounds_exist) then
       call put_att_text(ncid_ice, id_lh_ice, 'time_avg_info', 'average_T1,average_T2,average_DT')
    endif

    rcode = nf_enddef(ncid_ice)

    !--- write out axis data
    call put_var_real_1d(ncid_ice, id_xt, lon)
    call put_var_real_1d(ncid_ice, id_yt, lat)
    call put_var_real_1d(ncid_ice, id_xb, lonb)
    call put_var_real_1d(ncid_ice, id_yb, latb)
    call put_var_real_1d(ncid_ice, id_nv, (/1.,2./) )

    !--- setup for the part from addition of ice and land model
    rcode = nf_create(trim(output_from_tot),NF_WRITE, ncid_tot)
    call error_handler('error in creating file '//trim(output_from_tot), rcode) 

    !--- define dimension
    dims(1) = define_dim(ncid_tot, 'lon',ni_atm)
    dims(2) = define_dim(ncid_tot, 'lat',nj_atm)
    dims(3) = define_dim(ncid_tot, 'time', NF_UNLIMITED)
    dims(4) = define_dim(ncid_tot, 'lonb', ni_atm+1)
    dims(5) = define_dim(ncid_tot, 'latb', nj_atm+1)
    dims(6) = define_dim(ncid_tot, 'nv', 2)

    !--- define variable
    id_xt        = define_var_1d(ncid_tot, 'lon', NF_FLOAT, dims(1) )
    call put_att_text(ncid_tot, id_xt, 'long_name', 'longitude')
    call put_att_text(ncid_tot, id_xt, 'units', 'degrees_E')
    call put_att_text(ncid_tot, id_xt, 'cartesian_axis', 'X')
    call put_att_text(ncid_tot, id_xt, 'edges', 'lonb')

    id_yt        = define_var_1d(ncid_tot, 'lat', NF_FLOAT, dims(2) )
    call put_att_text(ncid_tot, id_yt, 'long_name', 'latitude')
    call put_att_text(ncid_tot, id_yt, 'units', 'degrees_N')
    call put_att_text(ncid_tot, id_yt, 'cartesian_axis', 'Y')
    call put_att_text(ncid_tot, id_yt, 'edges', 'latb')

    id_xb        = define_var_1d(ncid_tot, 'lonb', NF_FLOAT, dims(4) )
    call put_att_text(ncid_tot, id_xb, 'long_name', 'longitude edges')
    call put_att_text(ncid_tot, id_xb, 'units', 'degrees_E')
    call put_att_text(ncid_tot, id_xb, 'cartesian_axis', 'X')

    id_yb        = define_var_1d(ncid_tot, 'latb', NF_FLOAT, dims(5) )
    call put_att_text(ncid_tot, id_yb, 'long_name', 'latitude edges')
    call put_att_text(ncid_tot, id_yb, 'units', 'degrees_N')
    call put_att_text(ncid_tot, id_yb, 'cartesian_axis', 'Y')

    id_time_tot  = define_var_1d(ncid_tot, 'time', NF_DOUBLE, dims(3) )
    call put_att_text(ncid_tot, id_time_tot, 'long_name', 'time')
    call put_att_text(ncid_tot, id_time_tot, 'units', trim(tunits) )
    call put_att_text(ncid_tot, id_time_tot, 'cartesian_axis', 'T')
    call put_att_text(ncid_tot, id_time_tot, 'calendar_type', 'NOLEAP')
    call put_att_text(ncid_tot, id_time_tot, 'calendar', 'NOLEAP')
    if(time_bounds_exist) then
       call put_att_text(ncid_tot, id_time_tot, 'bounds', 'time_bounds')
    endif

    id_nv        = define_var_1d(ncid_tot, 'nv', NF_FLOAT, dims(6) )
    call put_att_text(ncid_tot, id_nv, 'long_name', 'vertex number')
    call put_att_text(ncid_tot, id_nv, 'units', 'none')
    call put_att_text(ncid_tot, id_nv, 'cartesian_axis', 'N')

    if(time_bounds_exist) then
       id_t1_tot    = define_var_1d(ncid_tot, 'average_T1', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_tot, id_t1_tot, 'long_name', 'Start time for average period')
       call put_att_text(ncid_tot, id_t1_tot, 'units', trim(tunits) )

       id_t2_tot    = define_var_1d(ncid_tot, 'average_T2', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_tot, id_t2_tot, 'long_name', 'End time for average period')
       call put_att_text(ncid_tot, id_t2_tot, 'units', trim(tunits) )

       id_dt_tot    = define_var_1d(ncid_tot, 'average_DT', NF_DOUBLE, dims(3) )
       call put_att_text(ncid_tot, id_dt_tot, 'long_name', 'Length of average period')
       call put_att_text(ncid_tot, id_dt_tot, 'units', trim(units))

       id_timeb_tot = define_var_2d(ncid_tot, 'time_bounds', NF_DOUBLE, (/dims(6), dims(3)/) )
       call put_att_text(ncid_tot, id_timeb_tot, 'long_name', 'time axis boundaries')
       call put_att_text(ncid_tot, id_timeb_tot, 'units', trim(units))
    endif

    id_lh_tot    = define_var_3d(ncid_tot, trim(fld_out_from_tot), NF_FLOAT, dims(1:3) )    
    call put_att_text(ncid_tot, id_lh_tot, 'long_name', 'latent heat flux')
    call put_att_text(ncid_tot, id_lh_tot, 'units', 'W/m^2')
    call put_att_real(ncid_tot, id_lh_tot, 'missing_value', missing_lnd)
    call put_att_text(ncid_tot, id_lh_tot, 'cell_methods', 'time: mean')
    if(time_bounds_exist) then
       call put_att_text(ncid_tot, id_lh_tot, 'time_avg_info', 'average_T1,average_T2,average_DT')
    endif

    rcode = nf_enddef(ncid_tot)

    !--- write out axis data
    call put_var_real_1d(ncid_tot, id_xt, lon)
    call put_var_real_1d(ncid_tot, id_yt, lat)
    call put_var_real_1d(ncid_tot, id_xb, lonb)
    call put_var_real_1d(ncid_tot, id_yb, latb)
    call put_var_real_1d(ncid_tot, id_nv, (/1.,2./) )

    return

  end subroutine setup_meta

  !#####################################################################

  subroutine process_data

    integer                         :: rcode, m, i, j
    real, dimension(ni_ice, nj_ice) :: data_src_ice
    real, dimension(ni_atm, nj_atm) :: data_src_lnd, data_dst_ice, data_dst_tot

    do m = 1, ntime

       write(stdout,*)'******* processing at time step: ', m

       !--- read input data
       call get_var_level_2d(src_from_ice, fld_in_from_ice, m, data_src_ice )
       call get_var_level_2d(src_from_lnd, fld_in_from_lnd, m, data_src_lnd )    

       !--- data conversion
       data_dst_ice = 0.0
       area_ocn     = 0.0
       do i = 1, size(i_atm_atmxocn)
          if(data_src_ice(i_ocn_atmxocn(i),j_ocn_atmxocn(i)) .ne. missing_ice ) then
             data_dst_ice(i_atm_atmxocn(i),j_atm_atmxocn(i)) = data_dst_ice(i_atm_atmxocn(i),j_atm_atmxocn(i)) &
                  + data_src_ice(i_ocn_atmxocn(i),j_ocn_atmxocn(i)) * area_atmxocn(i)
             area_ocn(i_atm_atmxocn(i),j_atm_atmxocn(i)) = area_ocn(i_atm_atmxocn(i),j_atm_atmxocn(i)) + area_atmxocn(i)
          endif
       enddo

       do j = 1, nj_atm
          do i = 1, ni_atm
             if(area_ocn(i,j) > epsln ) then
                data_dst_ice(i,j) = data_dst_ice(i,j)/area_ocn(i,j)
             else
                data_dst_ice(i,j) = missing_ice
             endif
          enddo
       enddo

       do j = 1, nj_atm
          do i = 1, ni_atm
             if(data_dst_ice(i,j) .ne. missing_ice ) then
                data_dst_tot(i,j) = data_dst_ice(i,j) * area_ocn(i,j)
             else
                data_dst_tot(i,j) = 0.0
             endif

             if(data_src_lnd(i,j) .ne. missing_lnd ) then
                data_dst_tot(i,j) = data_dst_tot(i,j) + data_src_lnd(i,j) * (area_atm(i,j)-area_ocn(i,j))
             endif
             data_dst_tot(i,j) = data_dst_tot(i,j)/area_atm(i,j)
          enddo
       enddo

       !--- write out data from ice model
       call put_var_level_2d(ncid_ice, id_lh_ice, m, data_dst_ice )
       call put_var_level_0d(ncid_ice, id_time_ice, m, time_value(m) ) 
       if(time_bounds_exist) then
          call put_var_level_0d(ncid_ice, id_t1_ice, m, t1(m) ) 
          call put_var_level_0d(ncid_ice, id_t2_ice, m, t2(m) ) 
          call put_var_level_0d(ncid_ice, id_dt_ice, m, dt(m) ) 
          call put_var_level_1d(ncid_ice, id_timeb_ice, m, time_bounds(:,m) )      
       endif

       !--- write out the total hfls
       call put_var_level_2d(ncid_tot, id_lh_tot, m, data_dst_tot )
       call put_var_level_0d(ncid_tot, id_time_tot, m, time_value(m) ) 
       if(time_bounds_exist) then
          call put_var_level_0d(ncid_tot, id_t1_tot, m, t1(m) ) 
          call put_var_level_0d(ncid_tot, id_t2_tot, m, t2(m) ) 
          call put_var_level_0d(ncid_tot, id_dt_tot, m, dt(m) ) 
          call put_var_level_1d(ncid_tot, id_timeb_tot, m, time_bounds(:,m) )   
       endif

    enddo

    rcode = nf_close(ncid_ice)
    rcode = nf_close(ncid_tot)

  end subroutine process_data

  !#####################################################################
  ! get the dimension length of any one dimensional variable
  function get_dimlen(ncid, name)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer                      :: get_dimlen
    integer                      :: varid, rcode, dims(1)

    rcode = nf_inq_varid(ncid, trim(name), varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_inq_vardimid(ncid, varid, dims)
    call error_handler('error in inquiring dimension id of '//trim(name), rcode)

    rcode = nf_inq_dimlen(ncid, dims(1), get_dimlen)
    call error_handler('error in inquiring dimension length of '//trim(name), rcode)

  end function get_dimlen

  !#####################################################################
  ! read the 1d integer data from netcdf file.
  subroutine get_var_int_1d(ncid, name, data)
    integer,                intent(in) :: ncid
    character(len=*),       intent(in) :: name
    integer, dimension(:), intent(out) :: data
    integer                            :: rcode, varid

    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)
    rcode = nf_get_var_int(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)

  end subroutine get_var_int_1d

  !#####################################################################
  ! read the 1d real data from netcdf file.
  subroutine get_var_real_1d(ncid, name, data)
    integer,             intent(in) :: ncid
    character(len=*),    intent(in) :: name
    real, dimension(:), intent(out) :: data
    integer                         :: rcode, varid

    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_get_var_double(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)


  end subroutine get_var_real_1d

  !#####################################################################
  ! read the 2d real data from netcdf file.
  subroutine get_var_real_2d(ncid, name, data)
    integer,               intent(in) :: ncid
    character(len=*),      intent(in) :: name
    real, dimension(:,:), intent(out) :: data
    integer                           :: rcode, varid


    rcode = nf_inq_varid(ncid, name, varid)
    call error_handler('error in inquiring variable id of '//trim(name), rcode)

    rcode = nf_get_var_double(ncid, varid, data)
    call error_handler('error in reading data of '//trim(name), rcode)

  end subroutine get_var_real_2d

  !#####################################################################
  subroutine get_missing_value(ncid, name, missing)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: name
    real,             intent(in) :: missing

    integer                      :: rcode, id_fld

    rcode = nf_inq_varid(ncid,trim(name), id_fld)
    call error_handler('error in inquring id of field '//trim(name), rcode)

    rcode = nf_get_att_double(ncid, id_fld, 'missing_value', missing )
    call error_handler('error in get field '//trim(name)//' missing value' , rcode)    
    rcode = nf_close(ncid)

  end subroutine get_missing_value

  !#####################################################################
  subroutine get_time_units(ncid)
     integer,           intent(in) :: ncid

    integer                      :: rcode, id_fld

    rcode = nf_inq_varid(ncid,'time', id_fld)
    call error_handler('error in inquring id of time', rcode)

    rcode = nf_get_att_text(ncid, id_fld, 'units', tunits )
    call error_handler('error in get time units' , rcode)    

  end subroutine get_time_units

  !#####################################################################
  function define_dim(ncid, name, size )
    integer,          intent(in) :: ncid, size
    character(len=*), intent(in) :: name
    integer                      :: define_dim

    integer                      :: rcode

    rcode = nf_def_dim(ncid, trim(name), size, define_dim)
    call error_handler('error in defining dimension '//trim(name) , rcode)   

  end function define_dim

  !#####################################################################
  function define_var_1d(ncid, name, type, dim )
    integer,          intent(in) :: ncid, type, dim
    character(len=*), intent(in) :: name
    integer                      :: define_var_1d

    integer                      :: rcode

     rcode = nf_def_var(ncid, trim(name), type, 1, dim, define_var_1d)
    call error_handler('error in defining variable '//trim(name) , rcode)  


  end function define_var_1d

  !#####################################################################
  function define_var_2d(ncid, name, type, dim )
    integer,               intent(in) :: ncid, type
    character(len=*),      intent(in) :: name
    integer, dimension(:), intent(in) :: dim
    integer                           :: define_var_2d

    integer                           :: rcode

     rcode = nf_def_var(ncid, trim(name), type, 2, dim, define_var_2d)
    call error_handler('error in defining variable '//trim(name) , rcode)  


  end function define_var_2d

  !#####################################################################

  function define_var_3d(ncid, name, type, dim )
    integer,               intent(in) :: ncid, type
    character(len=*),      intent(in) :: name
    integer, dimension(:), intent(in) :: dim
    integer                           :: define_var_3d

    integer                           :: rcode

     rcode = nf_def_var(ncid, trim(name), type, 3, dim, define_var_3d)
    call error_handler('error in defining variable '//trim(name) , rcode)  


  end function define_var_3d

  !#####################################################################

  subroutine get_var_level_2d(file, field, level, data )
     character(len=*), intent(in) :: file, field
     integer,          intent(in) :: level
     real,            intent(out) :: data(:,:)
    
     integer :: start(4), nread(4), ncid, rcode, id_fld, nlon, nlat

    nlon = size(data,1)
    nlat = size(data,2)

    rcode = nf_open(trim(file), NF_NOWRITE, ncid) 
    call error_handler('error in opening file '//trim(file), rcode) 

    rcode = nf_inq_varid(ncid, trim(field), id_fld )
    call error_handler('error in inquiring variable id of '//trim(field), rcode) 

    start = 1; nread = 1
    nread(1) = nlon; nread(2) = nlat
    start(3) = level
    rcode = nf_get_vara_double(ncid, id_fld, start, nread, data)
    call error_handler('Error in reading the data ', rcode)

    rcode = nf_close(ncid)

  end subroutine get_var_level_2d

  !#####################################################################
  subroutine put_att_text(ncid, id_fld, name, cval)
    integer,          intent(in) :: ncid, id_fld
    character(len=*), intent(in) :: name, cval

    integer                      :: rcode 

    rcode = nf_put_att_text(ncid,id_fld, trim(name),len_trim(cval), trim(cval) )
    call error_handler('error in putting attribute '//trim(name) , rcode)  

  end subroutine put_att_text
 
  !#####################################################################

  subroutine put_att_real(ncid, id_fld, name, rval)
    integer,          intent(in) :: ncid, id_fld
    character(len=*), intent(in) :: name
    real,             intent(in) :: rval

    integer                      :: rcode 

    rcode = nf_put_att_double(ncid,id_fld, trim(name), NF_DOUBLE, 1, rval )
    call error_handler('error in putting attribute '//trim(name) , rcode)  

  end subroutine put_att_real
 
  !#####################################################################

  subroutine put_var_real_1d(ncid, varid, data)
     integer,            intent(in) :: ncid, varid
     real, dimension(:), intent(in) :: data
     integer :: rcode, start(4), nwrite(4)

     start = 1; nwrite = 1;
     nwrite(1) = size(data)
     rcode = nf_put_vara_double(ncid, varid, start, nwrite, data)
     call error_handler('Error in put real 1d data', rcode )
     

  end subroutine  put_var_real_1d

  !#####################################################################
  subroutine put_var_level_0d(ncid, id_fld, level, data)
    integer, intent(in) :: ncid, id_fld, level
    real,    intent(in) :: data

    integer :: rcode

    rcode = nf_put_var1_double(ncid, id_fld, level, data)
    call error_handler('error in putting scalar data', rcode)

  end subroutine put_var_level_0d

  !#####################################################################

  subroutine put_var_level_1d(ncid, id_fld, level, data)
    integer, intent(in) :: ncid, id_fld, level
    real,    intent(in) :: data(:)

    integer :: nwrite(4), start(4), rcode

    start  = 1;  nwrite = 1 
    start(2)  = level; nwrite(1) = size(data)

    rcode = nf_put_vara_double(ncid, id_fld, start, nwrite, data)
    call error_handler('error in putting 1d data', rcode)

  end subroutine put_var_level_1d

  !#####################################################################
  subroutine put_var_level_2d(ncid, id_fld, level, data)
    integer, intent(in) :: ncid, id_fld, level
    real,    intent(in) :: data(:,:)

    integer             :: rcode, start(4), nwrite(4)

    start = 1; nwrite = 1;
    start(3) = level
    nwrite(1) = size(data,1); nwrite(2) = size(data,2)

    rcode = nf_put_vara_double(ncid, id_fld, start, nwrite, data)
    call error_handler('error in putting 2d data ', rcode)

  end subroutine put_var_level_2d

  !#####################################################################
  ! error handling routine.
  subroutine error_handler(mesg, status)
    character(len=*),  intent(in) :: mesg
    integer, optional, intent(in) :: status
    character(len=256) :: msg


    if(present(status)) then
       if(status == 0) return
       msg = nf_strerror(status)
       msg = trim(mesg)//': '// trim(msg)
    else
       msg = trim(mesg)
    endif

    write(stdout,*)'Error: '//trim(msg) 
    call abort()

  end subroutine error_handler

  !#####################################################################

end program regrid_hlfs

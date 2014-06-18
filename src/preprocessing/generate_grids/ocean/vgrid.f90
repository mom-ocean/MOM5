!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vgrid_mod
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
  ! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">Z. Liang </CONTACT>
  ! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">S. M. Griffies</REVIEWER>
  !

  !<OVERVIEW>
  ! <TT>vgrid_mod</TT> Generate vertical grid. 
  !</OVERVIEW>

  !<DESCRIPTION>
  !   The grid file contains the following information,

  ! <PRE>
  !
  !       zt = depth of tracer points
  !       zb = depth of tracer_boundaries
  !
  !              +---------------+
  !              |               |
  !              |               |
  !              |               |
  !              |               |
  !              |      +zt_k    |
  !              |               |
  !              |               |
  !              |               |
  !              +------+zb_k----+
  !
  ! </PRE>
  !</DESCRIPTION>

  use mpp_mod,        only : mpp_pe, mpp_root_pe, mpp_error, FATAL, NOTE,  mpp_chksum
  use mpp_mod,        only : lowercase
  use mpp_io_mod,     only : MPP_NETCDF, MPP_RDONLY, MPP_ASCII, MPP_MULTI, MPP_SINGLE
  use mpp_io_mod,     only : mpp_open, mpp_write_meta, mpp_write, axistype, mpp_close
  use mpp_io_mod,     only : mpp_get_info, mpp_get_atts, mpp_get_axes, mpp_get_axis_data 
  use fms_mod,        only : write_version_number, open_namelist_file, string
  use fms_mod,        only : file_exist, close_file, check_nml_error, stdlog, stdout
  use grids_type_mod, only : vgrid_data_type
  use grids_util_mod, only : make_axis, get_file_unit
  use constants_mod,  only : PI

  implicit none
  private

  integer, parameter :: maxlen=10000,maxbounds=11
  !------ namelist interface ---------------------------------------------
  !------ specify a spherical grid resolution in depth
  integer                    :: nzdepths = 0 
  real, dimension(maxbounds) :: z_depth, dz_depth
  logical                    :: read_my_grid = .false.
  character(len=128)         :: my_grid_file = 'my_vgrid'
  logical                    :: debug = .false.
  character(len=24)          :: z_axis_t   = 'zt_k'
  character(len=24)          :: z_axis_b   = 'zw_k'
  integer                    :: z_axis_b_offset = 1
  !
  !<NAMELIST NAME="vgrid_nml">
  !<DATA NAME="nzdepths" TYPE="integer">
  ! number of depth regions for varying resolution
  ! </DATA>
  ! <DATA NAME="z_depth" TYPE="real" DIM="(nxlons)" UNITS="meters">
  ! boundaries for defining depth regions of varying resolution
  ! </DATA>
  ! <DATA NAME="dz_depth" TYPE="real" DIM="(nxlons)" UNITS="meters">
  ! nominal resolution of depth regions
  ! </DATA>
  ! <DATA NAME="stretch_z" TYPE="real">
  !  stretch factor of vertical grids.
  ! </DATA>
  !<DATA NAME="read_my_grid" TYPE="logical">
  ! read ASCII grid information for supplying user-defined grids. 
  !</DATA>
  !<DATA NAME="my_grid_file" TYPE="character(len=128)">
  ! Name of ASCII or netcdf user grid file
  !</DATA>
  !<DATA NAME="z_axis_t" TYPE="character(len=24)">
  ! Name of z_t axis, if the file is netcdf
  !</DATA>
  !<DATA NAME="z_axis_b" TYPE="character(len=24)">
  ! Name of z_b axis, if the file is netcdf
  !</DATA>
  !<DATA NAME="z_axis_b_offset" TYPE="integer">
  ! offset of z_b axis, if the file is netcdf
  ! 1 corresponds to an axis starting at k=0, z_b=0 (mom3)
  !</DATA>
  !<DATA NAME="debug" TYPE="logical">
  ! control standard output.
  !</DATA>
  !</NAMELIST>
  namelist /vgrid_nml/ nzdepths, z_depth, dz_depth, read_my_grid, my_grid_file, debug,&
                       z_axis_b_offset, z_axis_t, z_axis_b
  !-----------------------------------------------------------------------
  !--------private data---------------------------------------------------
  integer                         :: nk   !number of grid cells in depth
  real, dimension(:), allocatable :: zt0, zb0
  type(axistype),save             :: axis_zt, axis_zb
  logical                         :: module_is_initialized = .false.
  !---------version information-------------------------------------------
  character(len=128) :: version = '$Id: vgrid.f90,v 13.0 2006/03/28 21:45:06 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $' 
  !---------public interface----------------------------------------------
  public :: generate_vgrid, vgrid_init, vgrid_end, write_vgrid_meta, write_vgrid_data


contains

  !#######################################################################
  ! <SUBROUTINE NAME="vgrid_init" >
  !   <OVERVIEW>
  !    Initialization routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    Read namelist, write out version and namelist informaiton and generate depth resolution.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call vgrid_init ( )
  !   </TEMPLATE>
  ! </SUBROUTINE>

  subroutine vgrid_init

    integer :: unit, ierr, io, k, pe
    integer :: ndim, nvar, natt, ntime, len, i, nkb
    character(len=128)          :: txt
    character(len=128)          :: name
    type(axistype), allocatable :: axes(:)
    logical                     :: found_z_t, found_z_b

    z_depth(:)  = 0.0
    dz_depth(:) = 0.0

    !---- read namelist --------------------------------------------------
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
       read  (unit, nml=vgrid_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'vgrid_nml') 
    enddo
10  call close_file (unit)

    !--- write version info and namelist to logfile ----------------------

    call write_version_number(version,tagname)
    if (mpp_pe() == mpp_root_pe()) then
       write (stdout(), nml=vgrid_nml)
    endif

    allocate (zt0(0:maxlen), zb0(0:maxlen))

    if (read_my_grid) then
       nk=0
       nkb=0
       if(.not. file_exist(trim(my_grid_file))) &
            call mpp_error(FATAL,'vgrid_mod: file '//trim(my_grid_file)//' does not exist')
       len = len_trim(my_grid_file)
       if(my_grid_file(len-2:len) == '.nc') then
          call mpp_open(unit,trim(my_grid_file),action=MPP_RDONLY,form=MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
          call mpp_get_info(unit, ndim, nvar, natt, ntime)
          allocate(axes(ndim))
          call mpp_get_axes(unit, axes)
          found_z_t = .false.; found_z_b = .false.
          do i=1, ndim
             call mpp_get_atts(axes(i),name=name,len=len)  
             if (trim(lowercase(name)) .eq. trim(lowercase(z_axis_t))) then
               found_z_t = .true.
               nk = len
               call mpp_get_axis_data(axes(i), zt0(1:nk) )
             endif  
             if (trim(lowercase(name)) .eq. trim(lowercase(z_axis_b))) then
               found_z_b = .true.
               nkb = len
               call mpp_get_axis_data(axes(i), zb0(1-z_axis_b_offset:nkb-z_axis_b_offset) )
             endif  
          enddo
          if(.not. found_z_t) call mpp_error(FATAL,'axis '//trim(z_axis_t)//' does not exist in file '//trim(my_grid_file) )
          if(.not. found_z_b) call mpp_error(FATAL,'axis '//trim(z_axis_b)//' does not exist in file '//trim(my_grid_file) )
          deallocate(axes)
          call mpp_close(unit)
       else
          call mpp_open(unit,trim(my_grid_file),action=MPP_RDONLY,form=MPP_ASCII,  &
               threading=MPP_MULTI,fileset=MPP_SINGLE)
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) nk
          read(unit,*) zt0(1:nk)
          read(unit,*) txt ! header line (ignored) 
          read(unit,*) 
          read(unit,*) zb0(1:nk)
          call mpp_close(unit)
       endif
       if (nk == 0) call mpp_error(FATAL,'vgrid_mod: error reading file '//trim(my_grid_file) )

    else
       write (stdout(),'(//,36x,a,/)') 'V G R I D   G E N E R A T I O N'
       !--- check namelist
       if (nzdepths .lt. 2 ) then 
          call mpp_error(FATAL, 'vgrid_mod: nml "nzdepths" = '//trim(string(nzdepths))//' should be no less than 2')
       endif
       if (nzdepths .gt. maxbounds) then
          call mpp_error(FATAL, 'vgrid_mod: maxbouds= '//trim(string(maxbounds))// &
               ' is less than nml nzdepths= '//trim(string(nzdepths)) )
       endif

       write (stdout(),'(/a)') 'Generating the vertical resolution for the computational domain:'
       call make_axis('Z',maxlen, nzdepths,z_depth,dz_depth,zb0,zt0, nk,.false., .false.)

    endif

    ! print out vgrid coordinates

    write (stdout(),9103) nk
    write (stdout(),9002) (zt0(k),k=1,nk)
    write (stdout(),9104) nk
    write (stdout(),9002) (zb0(k),k=1,nk)

    ! Compute a grid checksum
    if(debug) then
       write (stdout(),'(/)')
       pe = mpp_pe()
       write (stdout(),*) 'Vertical Grid checksum = ',  mpp_chksum(zt0(1:nk),(/pe/))+&
            mpp_chksum(zb0(1:nk),(/pe/))
       write (stdout(),'(/)')
    endif
    module_is_initialized = .true.

    return

9002 format (1x,10f10.2)
9101 format (/,  a,i4,' in units of ',a,' as follows:')
9103 format (/,' Depth to T and U grid points (m): zt(k) k=1,',i3)
9104 format (/,' Depth to W grid points (m): zb(k) k=1,',i3)

  end subroutine vgrid_init

  !#######################################################################
  ! <SUBROUTINE NAME="generate_vgrid" >
  !   <OVERVIEW>
  !    Generate vertical grid.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call generate_vgrid (Vgrid)
  !   </TEMPLATE>
  !   <INOUT NAME="Vgrid" TYPE="vgrid_data_type">
  !     A derived-type variable that contains vertical grid information.
  !   </INOUT>
  ! </SUBROUTINE>

  subroutine generate_vgrid (Vgrid)

    type(vgrid_data_type), intent(inout) :: Vgrid

    !--- allocate memory to data type ------------------------------------

    allocate(Vgrid%zt(nk), Vgrid%zb(nk))

    Vgrid%zt(:) = zt0(1:nk)   
    Vgrid%zb(:) = zb0(1:nk)

    deallocate(zt0, zb0)

  end subroutine generate_vgrid

  !#######################################################################
  ! <SUBROUTINE NAME="write_vgrid_data">
  !   <OVERVIEW>
  !     write the vertical grid data to netcdf file
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_vgrid_data (unit)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_vgrid_data (file)
    character(len=*), intent(in) :: file
    integer                      :: unit

    if(mpp_pe() .ne. mpp_root_pe() ) return

    unit = get_file_unit(file)

    call mpp_write(unit,axis_zt)
    call mpp_write(unit,axis_zb)

  end subroutine write_vgrid_data

  !#######################################################################
  ! <SUBROUTINE NAME="write_vgrid_meta">

  !   <OVERVIEW>
  !     Write out vertical grid meta data.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     call write_vgrid_meta(unit, Vgrid)
  !   </TEMPLATE>
  !   <IN NAME="unit" TYPE="integer" >
  !     The unit corresponding the output netcdf file. Always is returned by mpp_open.
  !   </IN>
  !   <IN NAME="Vgrid" TYPE="vgrid_data_type">
  !     A derived-type variable that contains vertical grid information.
  !   </IN>
  ! </SUBROUTINE>

  subroutine write_vgrid_meta(file, Vgrid)
    character(len=*),      intent(in) :: file
    type(vgrid_data_type), intent(in) :: Vgrid
    integer                           :: unit

    if(mpp_pe() .ne. mpp_root_pe() ) return

    unit = get_file_unit(file)

    call mpp_write_meta(unit,axis_zt,'zt','meters','zt',&
         data=Vgrid%zt(1:nk),cartesian='z',sense=-1)
    call mpp_write_meta(unit,axis_zb,'zb','meters','zb',&
         data=Vgrid%zb,cartesian='z',sense=-1)

    return

  end subroutine write_vgrid_meta

  !#######################################################################
  ! <SUBROUTINE NAME="vgrid_end">

  !   <OVERVIEW>
  !     Destruction routine.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Deallocates memory used by "vgrid_data_type" variables.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call vgrid_end ( Vgrid )
  !   </TEMPLATE>
  !   <INOUT NAME="Vgrid" TYPE="vgrid_data_type">
  !     A derived-type variable that contains vertical grid information.
  !   </INOUT>
  ! </SUBROUTINE>
  subroutine vgrid_end(Vgrid)
    type(vgrid_data_type), intent(inout) :: Vgrid

    deallocate(Vgrid%zt, Vgrid%zb )
    module_is_initialized = .false.

  end subroutine vgrid_end

  !#######################################################################


end module vgrid_mod

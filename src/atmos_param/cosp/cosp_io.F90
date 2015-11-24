#include "cosp_defs.H"
#ifdef COSP_GFDL
 
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------

! $Id: cosp_io.F90,v 20.0 2013/12/13 23:15:43 fms Exp $
! $Name: tikal $
! cosp_version = 1.3.2

#endif

! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2008 - A. Bodas-Salcedo - Initial version
! Oct 2008 - S. Bony - In nc_write_cosp_1d and nc_write_cosp_2d :
!                      the label of layered cloud fractions was wrong -> corrected
!                      (before: low was actually mid, mid was high, high was total,
!                      total was low)
! Sep 2009 - A. Bodas-Salcedo - CMIP5 variable names implemented
!

!#include "cosp_defs.h" 
#ifndef COSP_GFDL
#include "cosp_defs.h" 
#endif
MODULE MOD_COSP_IO
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
#ifndef COSP_GFDL
  USE cmor_users_functions
#endif
  USE netcdf
  use MOD_COSP_Modis_Simulator
#ifdef COSP_GFDL
  use mpp_mod, only: input_nml_file
  use fms_mod, only: open_namelist_file, open_file, close_file,   &
                     file_exist, mpp_pe, mpp_root_pe,   &
                     error_mesg, FATAL, &
                     check_nml_error, write_version_number, stdlog
#endif
  
  IMPLICIT NONE
  
#ifdef COSP_GFDL
!---------------------------------------------------------------------
!----------- version number for this module --------------------------
character(len=128)  :: version =  '$Id: cosp_io.F90,v 20.0 2013/12/13 23:15:43 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'
#endif

  ! Types to be used as arrays of pointers
  TYPE var1d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(3)
     integer :: dimssz(2)
     real,pointer,dimension(:) :: pntr
  END TYPE
  TYPE var2d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(4)
     integer :: dimssz(3)
     real,pointer,dimension(:,:) :: pntr
  END TYPE
  TYPE var3d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(5)
     integer :: dimssz(4)
     real,pointer,dimension(:,:,:) :: pntr
  END TYPE
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR1D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR1D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(3)
     integer,intent(in) :: dimssz(2)
     real,dimension(:),target,intent(in) :: pntr
     type(var1d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR2D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR2D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(4)
     integer,intent(in) :: dimssz(3)
     real,dimension(:,:),target,intent(in) :: pntr
     type(var2d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR3D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR3D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(5)
     integer,intent(in) :: dimssz(4)
     real,dimension(:,:,:),target,intent(in) :: pntr
     type(var3d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR3D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE MAP_POINT_TO_LL---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_POINT_TO_LL(Nx,Ny,geomode,x1,x2,x3,x4,y2,y3,y4,y5)
     ! Input arguments
     integer,intent(in) :: Nx,Ny,geomode
     real,intent(in),optional :: x1(:),x2(:,:),x3(:,:,:), &
                                 x4(:,:,:,:)
     real,intent(out),optional :: y2(:,:),y3(:,:,:), &
                                  y4(:,:,:,:),y5(:,:,:,:,:)
     ! Local variables
     integer :: Npoints
     integer :: px(Nx*Ny),py(Nx*Ny)
     integer :: i,j,k,l,m
     integer :: Ni,Nj,Nk,Nl
     integer :: Mi,Mj,Mk,Ml,Mm
     character(len=128) :: proname='MAP_POINT_TO_LL'

     Npoints = Nx*Ny
     
     px=0
     py=0
     ! Obtain pointers to do the mapping
     if (geomode == 2) then ! (lon,lat) mode
      do j=1,Ny
        do i=1,Nx
            k = (j-1)*Nx+i
            px(k) = i  
            py(k) = j  
        enddo
      enddo
     else if (geomode == 3) then ! (lon,lat) mode
      do j=1,Nx
        do i=1,Ny
            k = (j-1)*Ny+i
            px(k) = j
            py(k) = i  
        enddo
      enddo
     else
#ifdef COSP_GFDL
       call error_mesg ('cosp_io:map_point_to_ll',  &
                    ' -- '//trim(proname)//': geomode not supported, ', &
                                                             FATAL)
#else
       print *, ' -- '//trim(proname)//': geomode not supported, ',geomode
       stop
#endif
     endif

     if (present(x1).and.present(y2)) then
        Ni = size(x1,1)
        Mi = size(y2,1)
        Mj = size(y2,2)
        if (Mi*Mj /= Ni) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
           ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 1)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 1)'
          stop
#endif
        endif
        do i=1,Npoints
          y2(px(i),py(i)) = x1(i)
        enddo
     else if (present(x2).and.present(y3)) then
        Ni = size(x2,1)
        Nj = size(x2,2)
        Mi = size(y3,1)
        Mj = size(y3,2)
        Mk = size(y3,3)
        if (Mi*Mj /= Ni) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
           ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 2)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 2)'
          stop
#endif
        endif
        if (Nj /= Mk) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
           ' -- '//trim(proname)//': Nj /= Mk (opt 2)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 2)'
          stop
#endif
        endif
        do k=1,Mk
         do i=1,Npoints
            y3(px(i),py(i),k) = x2(i,k)
         enddo
        enddo
     else if (present(x3).and.present(y4)) then
        Ni = size(x3,1)
        Nj = size(x3,2)
        Nk = size(x3,3)
        Mi = size(y4,1)
        Mj = size(y4,2)
        Mk = size(y4,3)
        Ml = size(y4,4)
        if (Mi*Mj /= Ni) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
           ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 3)'
          stop
#endif 
        endif
        if (Nj /= Mk) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
           ' -- '//trim(proname)//': Nj /= Mk (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 3)'
          stop
#endif
        endif
        if (Nk /= Ml) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
                   ' -- '//trim(proname)//': Nk /= Ml (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nk /= Ml (opt 3)'
          stop
#endif
        endif
        do l=1,Ml
         do k=1,Mk
          do i=1,Npoints
            y4(px(i),py(i),k,l) = x3(i,k,l)
          enddo
         enddo
        enddo
     else if (present(x4).and.present(y5)) then
        Ni = size(x4,1)
        Nj = size(x4,2)
        Nk = size(x4,3)
        Nl = size(x4,4)
        Mi = size(y5,1)
        Mj = size(y5,2)
        Mk = size(y5,3)
        Ml = size(y5,4)
        Mm = size(y5,5)
        if (Mi*Mj /= Ni) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
                 ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 4)'
          stop
#endif
        endif
        if (Nj /= Mk) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
                   ' -- '//trim(proname)//': Nj /= Mk (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 4)'
          stop
#endif
        endif
        if (Nk /= Ml) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
                   ' -- '//trim(proname)//': Nk /= Ml (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nk /= Ml (opt 4)'
          stop
#endif
        endif
        if (Nl /= Mm) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_point_to_ll',  &
                   ' -- '//trim(proname)//': Nl /= Mm (opt 4)', &
                                                             FATAL)
#else
          print *, ' -- '//trim(proname)//': Nl /= Mm (opt 4)'
          stop
#endif
        endif
        do m=1,Mm
         do l=1,Ml
          do k=1,Mk
            do i=1,Npoints
                y5(px(i),py(i),k,l,m) = x4(i,k,l,m)
            enddo
          enddo
         enddo
        enddo
     else
#ifdef COSP_GFDL
        call error_mesg ('cosp_io:map_point_to_ll',  &
                 ' -- '//trim(proname)//': wrong option', &
                                                            FATAL)
#else
        print *, ' -- '//trim(proname)//': wrong option'
        stop
#endif
     endif

     
  END SUBROUTINE MAP_POINT_TO_LL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE MAP_LL_TO_POINT---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_LL_TO_POINT(Nx,Ny,Np,x2,x3,x4,x5,y1,y2,y3,y4)
     ! Input arguments
     integer,intent(in) :: Nx,Ny,Np
     real,intent(in),optional :: x2(:,:),x3(:,:,:), &
                                 x4(:,:,:,:),x5(:,:,:,:,:)
     real,intent(out),optional :: y1(:),y2(:,:),y3(:,:,:), &
                                 y4(:,:,:,:)
     ! Local variables
     integer :: px(Nx*Ny),py(Nx*Ny)
     integer :: i,j,k,l,m
     integer :: Ni,Nj,Nk,Nl,Nm
     integer :: Mi,Mj,Mk,Ml
     character(len=128) :: proname='MAP_LL_TO_POINT'
     
     px=0
     py=0
     if (Nx*Ny < Np) then
#ifdef COSP_GFDL
       call error_mesg ('cosp_io:map_ll_to_point',  &
                ' -- '//trim(proname)//': Nx*Ny < Np', &
                                                            FATAL)
#else
       print *, ' -- '//trim(proname)//': Nx*Ny < Np'
       stop
#endif
     endif
     do j=1,Ny
       do i=1,Nx
          k = (j-1)*Nx+i
          px(k) = i  
          py(k) = j  
       enddo
     enddo
     
     if (present(x2).and.present(y1)) then
        Ni = size(x2,1)
        Nj = size(x2,2)
        Mi = size(y1,1)
        if (Ni*Nj < Mi) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 1)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 1)'
          stop
#endif
        endif
        do j=1,Np
          y1(j) = x2(px(j),py(j))
        enddo
     else if (present(x3).and.present(y2)) then
        Ni = size(x3,1)
        Nj = size(x3,2)
        Nk = size(x3,3)
        Mi = size(y2,1)
        Mj = size(y2,2)
        if (Ni*Nj < Mi) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
             ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 2)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 2)'
          stop
#endif
        endif
        if (Nk /= Mj) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nk /= Mj (opt 2)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 2)'
          stop
#endif
        endif
        do k=1,Nk
          do j=1,Np
            y2(j,k) = x3(px(j),py(j),k)
          enddo
        enddo
     else if (present(x4).and.present(y3)) then
        Ni = size(x4,1)
        Nj = size(x4,2)
        Nk = size(x4,3)
        Nl = size(x4,4)
        Mi = size(y3,1)
        Mj = size(y3,2)
        Mk = size(y3,3)
        if (Ni*Nj < Mi) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
              ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 3)'
          stop
#endif
        endif
        if (Nk /= Mj) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nk /= Mj (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 3)'
          stop
#endif
        endif
        if (Nl /= Mk) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nl /= Mk (opt 3)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 3)'
          stop
#endif
        endif
        do l=1,Nl
         do k=1,Nk
          do j=1,Np
            y3(j,k,l) = x4(px(j),py(j),k,l)
          enddo
         enddo
        enddo
     else if (present(x5).and.present(y4)) then
        Ni = size(x5,1)
        Nj = size(x5,2)
        Nk = size(x5,3)
        Nl = size(x5,4)
        Nm = size(x5,5)
        Mi = size(y4,1)
        Mj = size(y4,2)
        Mk = size(y4,3)
        Ml = size(y4,4)
        if (Ni*Nj < Mi) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
               ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 4)', &
                                                             FATAL)
#else
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 4)'
          stop
#endif
        endif
        if (Nk /= Mj) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nk /= Mj (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 4)'
          stop
#endif
        endif
        if (Nl /= Mk) then
#ifdef COSP_GFDL
          call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nl /= Mk (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 4)'
          stop
#endif
        endif
        if (Nm /= Ml) then
#ifdef COSP_GFDL
           call error_mesg ('cosp_io:map_ll_to_point',  &
                   ' -- '//trim(proname)//': Nm /= Ml (opt 4)', &
                                                            FATAL)
#else
          print *, ' -- '//trim(proname)//': Nm /= Ml (opt 4)'
          stop
#endif
        endif
        do m=1,Nm
         do l=1,Nl
          do k=1,Nk
           do j=1,Np
            y4(j,k,l,m) = x5(px(j),py(j),k,l,m)
           enddo
          enddo
         enddo
        enddo
     else
#ifdef COSP_GFDL
        call error_mesg ('cosp_io:map_ll_to_point',  &
                 ' -- '//trim(proname)//': wrong option', &
                                                            FATAL)
#else
        print *, ' -- '//trim(proname)//': wrong option'
        stop
#endif
     endif
  
  END SUBROUTINE MAP_LL_TO_POINT
  
#ifndef COSP_GFDL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE NC_READ_INPUT_FILE -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_READ_INPUT_FILE(fname,Npnts,Nl,Nhydro,lon,lat,p,ph,z,zh,T,qv,rh,tca,cca, &
            mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
            fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,dem_c,skt,landmask,sfc_height, &
            mr_ozone,u_wind,v_wind,sunlit,emsfc_lw,mode,Nlon,Nlat,verbosity)
    
    !Arguments
    character(len=512),intent(in) :: fname ! File name
    integer,intent(in) :: Npnts,Nl,Nhydro
    real,dimension(Npnts),intent(out) :: lon,lat
    real,dimension(Npnts,Nl),target,intent(out) :: p,ph,z,zh,T,qv,rh,tca,cca, &
                  mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
                  fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
    real,dimension(Npnts,Nl,Nhydro),intent(out) :: Reff
    real,dimension(Npnts),intent(out) :: skt,landmask,sfc_height,u_wind,v_wind,sunlit
    real,intent(out) :: emsfc_lw
    integer,intent(out) :: mode,Nlon,Nlat
    integer,optional :: verbosity
        
    !Local variables
    integer :: Npoints,Nlevels,i,j,k
    character(len=128) :: vname
    integer,parameter :: NMAX_DIM=5
    integer :: vrank,vdimid(NMAX_DIM)
    character(len=256) :: dimname(NMAX_DIM) ! 256 hardcoded, instead of MAXNCNAM. This works for NetCDF 3 and 4.
    integer :: ncid,vid,ndims,nvars,ngatts,recdim,dimsize(NMAX_DIM)
    integer :: errst
    logical :: Llat,Llon,Lpoint
    integer :: Na,Nb,Nc,Nd,Ne
    real,dimension(Npnts) :: ll
    integer,dimension(:),allocatable :: plon,plat
    real,allocatable :: x1(:),x2(:,:),x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:) ! Temporary arrays
    character(len=64) :: routine_name='NC_READ_INPUT_FILE'
    character(len=128) :: errmsg,straux
    
    mode = 0
    Nlon = 0
    Nlat = 0
    
    Npoints = Npnts
    Nlevels = Nl
    
    ! Open file
    errst = nf90_open(fname, nf90_nowrite, ncid)
    if (errst /= 0) then
      errmsg="Couldn't open "//trim(fname)
      call cosp_error(routine_name,errmsg)
    endif
    
    ! Get information about dimensions. Curtain mode or lat/lon mode?
    Llat  =.false.
    Llon  =.false.
    Lpoint=.false.
    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    if (errst /= 0) then
      errmsg="Error in  nf90_inquire"
      call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    do i = 1,ndims
       errst = nf90_Inquire_Dimension(ncid,i,name=dimname(i),len=dimsize(i))
       if (errst /= 0) then
         write(straux, *)  i
         errmsg="Error in nf90_Inquire_Dimension, i: "//trim(straux)
         call cosp_error(routine_name,errmsg)
       endif
       if ((trim(dimname(i)).eq.'level').and.(Nlevels > dimsize(i))) then
         errmsg='Number of levels selected is greater than in input file '//trim(fname)
         call cosp_error(routine_name,errmsg)
       endif
       if (trim(dimname(i)).eq.'point') then
         Lpoint = .true.
         if (Npnts > dimsize(i)) then
           errmsg='Number of points selected is greater than in input file '//trim(fname)
           call cosp_error(routine_name,errmsg)
         endif
       endif
       if (trim(dimname(i)).eq.'lon') then
         Llon = .true.
         Nlon = dimsize(i)
       endif
       if (trim(dimname(i)).eq.'lat') then
         Llat = .true.
         Nlat = dimsize(i)
       endif
    enddo

    ! Get lon and lat
    if (Llon.and.Llat) then ! 2D mode
        if ((Npnts) > Nlon*Nlat) Npoints=Nlon*Nlat
        lon = R_UNDEF
        lat = R_UNDEF
        mode = 2 ! Don't know yet if (lon,lat) or (lat,lon) at this point
    else if (Lpoint) then ! 1D mode
        Nlon = Npoints
        Nlat = Npoints
        mode = 1
    else
        errmsg= trim(fname)//' file contains wrong dimensions'
        call cosp_error(routine_name,errmsg)
    endif
    errst = nf90_inq_varid(ncid, 'lon', vid)
    if (errst /= 0) then
        errmsg="Error in nf90_inq_varid, var: lon"
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/Nlon/))
    if (errst /= 0) then
        errmsg="Error in nf90_get_var, var: lon"
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_inq_varid(ncid, 'lat', vid)
    if (errst /= 0) then
        errmsg="Error in nf90_inq_varid, var: lat"
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/Nlat/))
    if (errst /= 0) then
        errmsg="Error in nf90_get_var, var: lat"
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    
    ! Get all variables
    do vid = 1,nvars
       vdimid=0
       errst = nf90_Inquire_Variable(ncid, vid, name=vname, ndims=vrank, dimids=vdimid)
       if (errst /= 0) then
         write(straux, *)  vid
         errmsg='Error in nf90_Inquire_Variable, vid '//trim(straux)
         call cosp_error(routine_name,errmsg,errcode=errst)
       endif
       ! Read in into temporary array of correct shape
       if (present(verbosity).and.(verbosity == 1)) print *, 'Reading '//trim(vname)//' ...'
       if (vrank == 1) then
          Na = dimsize(vdimid(1))
          allocate(x1(Na))
          errst = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/Na/))
       endif
       if (vrank == 2) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          allocate(x2(Na,Nb))
          errst = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/Na,Nb/))
       endif
       if (vrank == 3) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          allocate(x3(Na,Nb,Nc))
          errst = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/Na,Nb,Nc/))
          if ((mode == 2).or.(mode == 3)) then
            if ((Na == Nlon).and.(Nb == Nlat)) then
              mode = 2
            else if ((Na == Nlat).and.(Nb == Nlon)) then
              mode = 3
            else
              errmsg='Wrong mode for variable '//trim(vname)
              call cosp_error(routine_name,errmsg)
            endif
          endif
       endif
       if (vrank == 4) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          allocate(x4(Na,Nb,Nc,Nd))
          errst = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/Na,Nb,Nc,Nd/))
       endif
       if (vrank == 5) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          Ne = dimsize(vdimid(5))
          allocate(x5(Na,Nb,Nc,Nd,Ne))
          errst = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/Na,Nb,Nc,Nd,Ne/))
       endif
       if (errst /= 0) then
          write(straux, *)  vid
          errmsg='Error in nf90_get_var, vid '//trim(straux)
          call cosp_error(routine_name,errmsg,errcode=errst)
       endif
       ! Map to the right input argument
       select case (trim(vname))
       case ('pfull')
         if (Lpoint) then
           p(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=p)
         endif
       case ('phalf')
         if (Lpoint) then
           ph(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=ph)
         endif
       case ('height')
         if (Lpoint) then
           z(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=z)
         endif
       case ('height_half')
         if (Lpoint) then
           zh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=zh)
         endif
       case ('T_abs')
         if (Lpoint) then
           T(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=T)
         endif
       case ('qv')
         if (Lpoint) then
           qv(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=qv)
         endif
       case ('rh')
         if (Lpoint) then
           rh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=rh)
         endif
       case ('tca')
         if (Lpoint) then
           tca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=tca)
         endif
         tca = tca
       case ('cca')
         if (Lpoint) then
           cca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=cca)
         endif
         cca = cca
       case ('mr_lsliq')
         if (Lpoint) then
           mr_lsliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsliq)
         endif
       case ('mr_lsice')
         if (Lpoint) then
           mr_lsice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsice)
         endif
       case ('mr_ccliq')
         if (Lpoint) then
           mr_ccliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccliq)
         endif
       case ('mr_ccice')
         if (Lpoint) then
           mr_ccice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccice)
         endif
       case ('fl_lsrain')
         if (Lpoint) then
           fl_lsrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsrain)
         endif
       case ('fl_lssnow')
         if (Lpoint) then
           fl_lssnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lssnow)
         endif
       case ('fl_lsgrpl')
         if (Lpoint) then
           fl_lsgrpl(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsgrpl)
         endif
       case ('fl_ccrain')
         if (Lpoint) then
           fl_ccrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccrain)
         endif
       case ('fl_ccsnow')
         if (Lpoint) then
           fl_ccsnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccsnow)
         endif
       case ('dtau_s')
         if (Lpoint) then
           dtau_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_s)
         endif
       case ('dtau_c')
         if (Lpoint) then
           dtau_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_c)
         endif
       case ('dem_s')
         if (Lpoint) then
           dem_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_s)
         endif
       case ('dem_c')
         if (Lpoint) then
           dem_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_c)
         endif
       case ('Reff')
         if (Lpoint) then
           Reff(1:Npoints,:,:) = x3(1:Npoints,1:Nlevels,:)
         else
           call map_ll_to_point(Na,Nb,Npoints,x4=x4,y3=Reff)
         endif
       case ('skt')
         if (Lpoint) then
           skt(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=skt)
         endif
       case ('landmask')
         if (Lpoint) then
           landmask(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=landmask)
         endif
       case ('orography')
         if (Lpoint) then
           sfc_height(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sfc_height)
         endif
       case ('mr_ozone')
         if (Lpoint) then
           mr_ozone(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ozone)
         endif
       case ('u_wind')
         if (Lpoint) then
           u_wind(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=u_wind)
         endif
       case ('v_wind')
         if (Lpoint) then
           v_wind(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=v_wind)
         endif
       case ('sunlit')
         if (Lpoint) then
           sunlit(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sunlit)
         endif
       end select
!        select case (trim(vname))
!        case ('pfull')
!          if (Lpoint) then
!            p(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=p)
!          endif
!        case ('phalf')
!          if (Lpoint) then
!            ph(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=ph)
!          endif
!        case ('zfull')
!          if (Lpoint) then
!            z(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=z)
!          endif
!        case ('zhalf')
!          if (Lpoint) then
!            zh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=zh)
!          endif
!        case ('ta')
!          if (Lpoint) then
!            T(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=T)
!          endif
!        case ('hus')
!          if (Lpoint) then
!            qv(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=qv)
!          endif
!        case ('hur')
!          if (Lpoint) then
!            rh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=rh)
!          endif
!        case ('cl')
!          if (Lpoint) then
!            tca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=tca)
!          endif
!          tca = tca/100.0
!        case ('clc')
!          if (Lpoint) then
!            cca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=cca)
!          endif
!          cca = cca/100.0
!        case ('clws')
!          if (Lpoint) then
!            mr_lsliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsliq)
!          endif
!        case ('clis')
!          if (Lpoint) then
!            mr_lsice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsice)
!          endif
!        case ('clwc')
!          if (Lpoint) then
!            mr_ccliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccliq)
!          endif
!        case ('clic')
!          if (Lpoint) then
!            mr_ccice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccice)
!          endif
!        case ('prsprof')
!          if (Lpoint) then
!            fl_lsrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsrain)
!          endif
!        case ('prsns')
!          if (Lpoint) then
!            fl_lssnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lssnow)
!          endif
!        case ('grplprof')
!          if (Lpoint) then
!            fl_lsgrpl(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsgrpl)
!          endif
!        case ('prcprof')
!          if (Lpoint) then
!            fl_ccrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccrain)
!          endif
!        case ('prsnc')
!          if (Lpoint) then
!            fl_ccsnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccsnow)
!          endif
!        case ('dtaus')
!          if (Lpoint) then
!            dtau_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_s)
!          endif
!        case ('dtauc')
!          if (Lpoint) then
!            dtau_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_c)
!          endif
!        case ('dems')
!          if (Lpoint) then
!            dem_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_s)
!          endif
!        case ('demc')
!          if (Lpoint) then
!            dem_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_c)
!          endif
!        case ('reff')
!          if (Lpoint) then
!            Reff(1:Npoints,:,:) = x3(1:Npoints,1:Nlevels,:)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x4=x4,y3=Reff)
!          endif
!        case ('ts')
!          if (Lpoint) then
!            skt(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=skt)
!          endif
!        case ('landmask')
!          if (Lpoint) then
!            landmask(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=landmask)
!          endif
!        case ('orog')
!          if (Lpoint) then
!            sfc_height(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sfc_height)
!          endif
!        case ('mrozone')
!          if (Lpoint) then
!            mr_ozone(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ozone)
!          endif
!        case ('uas')
!          if (Lpoint) then
!            u_wind(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=u_wind)
!          endif
!        case ('vas')
!          if (Lpoint) then
!            v_wind(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=v_wind)
!          endif
!        case ('sunlit')
!          if (Lpoint) then
!            sunlit(1:Npoints) = x1(1:Npoints)
!          else
!            call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sunlit)
!          endif
!        end select
       ! Free memory
       if (vrank == 1) deallocate(x1)
       if (vrank == 2) deallocate(x2)
       if (vrank == 3) deallocate(x3)
       if (vrank == 4) deallocate(x4)
       if (vrank == 5) deallocate(x5)
    enddo
       
    ! SFC emissivity
    errst = nf90_inq_varid(ncid, 'emsfc_lw', vid)
    if (errst /= 0) then
        errmsg='Error in nf90_inq_varid, var: emsfc_lw'
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_get_var(ncid, vid, emsfc_lw)
    if (errst /= 0) then
        errmsg='Error in nf90_get_var, var: emsfc_lw'
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif

    
    ! Fill in the lat/lon vectors with the right values for 2D modes
    ! This might be helpful if the inputs are 2D (gridded) and 
    ! you want outputs in 1D mode
    allocate(plon(Npoints),plat(Npoints))
    if (mode == 2) then !(lon,lat)
      ll = lat
      do j=1,Nb
        do i=1,Na
          k = (j-1)*Na + i
          plon(k) = i  
          plat(k) = j
        enddo
      enddo
      lon(1:Npoints) = lon(plon(1:Npoints))
      lat(1:Npoints) = ll(plat(1:Npoints))
    else if (mode == 3) then !(lat,lon)
      ll = lon
      do j=1,Nb
        do i=1,Na
          k = (j-1)*Na + i
          lon(k) = ll(j)
          lat(k) = lat(i)
        enddo
      enddo
      lon(1:Npoints) = ll(plon(1:Npoints))
      lat(1:Npoints) = lat(plat(1:Npoints))
    endif
    deallocate(plon,plat)
    
    ! Close file
    errst = nf90_close(ncid)
    if (errst /= 0) then
        errmsg='Error in nf90_close'
        call cosp_error(routine_name,errmsg,errcode=errst)
    endif

  END SUBROUTINE NC_READ_INPUT_FILE
#endif

#ifndef COSP_GFDL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE NC_WRITE_COSP_1D ---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef RTTOV
  SUBROUTINE NC_WRITE_COSP_1D(cmor_nl,wmode,cfg,vgrid,gb,sg,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
  SUBROUTINE NC_WRITE_COSP_1D(cmor_nl,wmode,cfg,vgrid,gb,sg,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif
   ! Input arguments
   character(len=*),intent(in) :: cmor_nl
   character(len=*) :: wmode ! Writing mode 'replace' or 'append'
   type(cosp_config),intent(in) :: cfg
   type(cosp_vgrid),intent(in) :: vgrid
   type(cosp_gridbox),intent(in) :: gb
   type(cosp_subgrid),intent(in) :: sg
   type(cosp_sgradar),intent(in) :: sgradar ! Subgrid radar
   type(cosp_sglidar),intent(in) :: sglidar ! Subgrid lidar
   type(cosp_isccp),intent(in)   :: isccp   ! ISCCP outputs
   type(cosp_misr),intent(in)    :: misr    ! MISR outputs
   type(cosp_modis),intent(in)   :: modis   ! MODIS outputs
#ifdef RTTOV 
   type(cosp_rttov),intent(in)   :: rttov   ! RTTOV outputs 
#endif
   type(cosp_radarstats),intent(in) :: stradar ! Summary statistics from radar simulator
   type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator
   !--- Local variables ---
   integer,parameter :: n3d = 8
#ifdef RTTOV
   integer,parameter :: n2d = 7
#else
   integer,parameter :: n2d = 6
#endif
   integer,parameter :: n1d = 29
   integer :: error_flag,i,j, nc_action
   logical :: lfound
   character(len=512) :: inpath,outpath,start_date,model_id,experiment_id,institution,institute_id,source,calendar, &
                 contact,history,comment,table,parent_experiment_id,parent_experiment_rip,forcing
   character(len=2056) :: references
   double precision :: branch_time
   integer :: initialization_method,physics_version
   integer :: realization,Npoints,Ncolumns,Nlevels,Nlvgrid,maxtsteps,Nchannels
   namelist/CMOR/inpath,outpath,start_date,model_id,experiment_id,branch_time,parent_experiment_id,parent_experiment_rip, &
              forcing,institution,institute_id,source,calendar,realization,initialization_method,physics_version, &
              contact,history,comment,references,table,maxtsteps
   real,dimension(:),allocatable :: profile_ax,column_ax,dbze_ax,channel_ax
   real,dimension(:,:),allocatable :: dbze_bounds,vgrid_bounds,mgrid_bounds,sratio_bounds
   integer :: profile_axid,column_axid,height_axid,dbze_axid,height_mlev_axid,sratio_axid, &
              tau_axid,pressure2_axid,time_axid,sza_axid,MISR_CTH_axid,channel_axid
   integer :: grid_id, latvar_id, lonvar_id, blatvar_id, blonvar_id
   integer :: var3d_id(n3d),var2d_id(n2d),var1d_id(n1d)
   type(var1d) :: v1d(n1d)
   type(var2d) :: v2d(n2d)
   type(var3d) :: v3d(n3d)
   integer :: d2(2),d3(3),d4(4),d5(5)
   logical,dimension(n1d) :: lout1d = .false.
   logical,dimension(n2d) :: lout2d = .false.
   logical,dimension(n3d) :: lout3d = .false.
   double precision :: tbnds(2,1)
   character(len=64) :: pro_name = 'NC_WRITE_COSP_1D'

   ! Initialise values of axis ids that may not be used
   height_mlev_axid = 0
   column_axid      = 0
   channel_axid     = 0
   tau_axid         = 0
   pressure2_axid   = 0
   misr_cth_axid    = 0

   Npoints   = gb%Npoints
   Ncolumns  = gb%Ncolumns
   Nlevels   = gb%Nlevels
   Nchannels = gb%Nchan
   Nlvgrid   = vgrid%Nlvgrid
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Allocate memory and initialise
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   tbnds(:,1) = gb%time_bnds
   allocate(profile_ax(Npoints),column_ax(Ncolumns),dbze_ax(DBZE_BINS), &
            dbze_bounds(2,DBZE_BINS),vgrid_bounds(2,Nlvgrid),mgrid_bounds(2,Nlevels), &
            sratio_bounds(2,SR_BINS), channel_ax(Nchannels))
   ! Profile
   do i=1,Npoints
     profile_ax(i) = i
   enddo
   ! Column
   do i=1,Ncolumns
     column_ax(i) = i
   enddo
   ! Channels
   channel_ax = float(gb%ichan)
   ! Radar Ze
   do i=1,DBZE_BINS
     dbze_ax(i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 0.5)
   enddo
   do i=1,DBZE_BINS
     dbze_bounds(1,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 1)
     dbze_bounds(2,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*i
   enddo
   ! Height of model levels
   do i=1,Nlevels
     mgrid_bounds(1,i) = vgrid%mzl(i)
     mgrid_bounds(2,i) = vgrid%mzu(i)
   enddo
   ! Height of std grid
   do i=1,Nlvgrid
     vgrid_bounds(1,i) = vgrid%zl(i)
     vgrid_bounds(2,i) = vgrid%zu(i)
   enddo
   ! Lidar scattering ratio bounds (They are output by cosp_cfad_sr->diag_lidar in lmd_ipsl_stats.f90)
   sratio_bounds(2,:)         = stlidar%srbval(:) ! srbval contains the upper limits from lmd_ipsl_stats.f90
   sratio_bounds(1,2:SR_BINS) = stlidar%srbval(1:SR_BINS-1)
   sratio_bounds(1,1)         = 0.0
   sratio_bounds(2,SR_BINS)   = 1.e5 ! This matches with Chepfer et al., JGR, 2009. However, it is not consistent 
                                     ! with the upper limit in lmd_ipsl_stats.f90, which is LIDAR_UNDEF-1=998.999
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Read namelist with information for CMOR output file
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   open(10,file=cmor_nl,status='old')
   read(10,nml=cmor)
   close(10)

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Specify path for tables and set up other CMOR options
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef USE_CMOR1
   error_flag = cmor_setup(inpath=trim(inpath),netcdf_file_action=trim(wmode))
#else
   nc_action = CMOR_APPEND_3
   if (trim(wmode) == 'replace') nc_action = CMOR_REPLACE_3
   error_flag = cmor_setup(inpath=trim(inpath),netcdf_file_action=nc_action,create_subdirectories=0)
#endif

   print *, '---------------Define dataset'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define dataset as output from COSP, and other model details
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef USE_CMOR1
   error_flag = cmor_dataset(outpath=trim(outpath),experiment_id=trim(experiment_id),institution=trim(institution), &
                    source=trim(source),calendar=trim(calendar),realization=realization,contact=trim(contact), &
                    history=trim(history),comment=trim(comment),references=trim(references))
#else
   error_flag = cmor_dataset(outpath=trim(outpath),experiment_id=trim(experiment_id),institution=trim(institution), &
                    source=trim(source),calendar=trim(calendar),realization=realization,contact=trim(contact), &
                    history=trim(history),comment=trim(comment),references=trim(references),model_id=trim(model_id), &
                    branch_time=branch_time,parent_experiment_id=trim(parent_experiment_id),forcing=trim(forcing), &
                    institute_id=trim(institute_id),parent_experiment_rip=trim(parent_experiment_rip), &
                    initialization_method=initialization_method,physics_version=physics_version)
   error_flag = cmor_set_cur_dataset_attribute('cosp_version',trim(COSP_VERSION))
#endif
   
   print *, '---------------Define axis'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define axis
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   profile_axid = cmor_axis(table=table, table_entry='location', units='1', length=Npoints, coord_vals=profile_ax)
   height_axid  = cmor_axis(table=table, table_entry='alt40', units='m', length=Nlvgrid, &
                                                       coord_vals=vgrid%z,cell_bounds=vgrid_bounds)
   dbze_axid    = cmor_axis(table=table, table_entry='dbze', units='dBZ', length=DBZE_BINS, &  ! '10 lg(re 1.e-18 m3)'
                            coord_vals=dbze_ax,cell_bounds=dbze_bounds)
   sratio_axid  = cmor_axis(table=table, table_entry='scatratio', units='1', length=SR_BINS, &
                            coord_vals=(sratio_bounds(1,:)+sratio_bounds(2,:))/2.0,cell_bounds=sratio_bounds)
   sza_axid  = cmor_axis(table=table, table_entry='sza5', units='degree', length=PARASOL_NREFL, coord_vals=PARASOL_SZA)
   time_axid  = cmor_axis(table=table, table_entry='time1', units='days since '//trim(start_date), length=maxtsteps)
   if ((trim(table) /= 'CMIP5_cf3hr').and.(trim(table) /= 'TAMIP_3hr_curt')) then
        column_axid  = cmor_axis(table=table, table_entry='column', units='1', length=Ncolumns, coord_vals=column_ax)
        channel_axid = cmor_axis(table=table, table_entry='channel', units='1', length=Nchannels, coord_vals=channel_ax)
        height_mlev_axid  = cmor_axis(table=table, table_entry='alevel', &
                                      units='1', length=Nlevels, coord_vals=vgrid%mz,cell_bounds=mgrid_bounds)
        tau_axid     = cmor_axis(table=table, table_entry='tau', units='1', length=7, &
                                    coord_vals=ISCCP_TAU,cell_bounds=ISCCP_TAU_BNDS)
        pressure2_axid = cmor_axis(table=table, table_entry='plev7', units='Pa', length=7, &
                                    coord_vals=ISCCP_PC,cell_bounds=ISCCP_PC_BNDS)
        MISR_CTH_axid = cmor_axis(table=table, table_entry='cth16', units='m', length=MISR_N_CTH, &
                                    coord_vals=MISR_CTH,cell_bounds=MISR_CTH_BNDS)
   endif
                            
   print *, '---------------Define grid'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define grid
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   grid_id = cmor_grid((/profile_axid, time_axid/))
   latvar_id = cmor_time_varying_grid_coordinate(grid_id,'latitude','degrees_north',R_UNDEF)
   lonvar_id = cmor_time_varying_grid_coordinate(grid_id,'longitude','degrees_east' ,R_UNDEF)
   if (grid_id > 0) then
        print *,  '*********** Error, grid_id: ', grid_id
        stop
   endif
   
   print *, '---------------Fill in var info'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Fill in variable info and associate pointers
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables
   d3 = (/grid_id,0,0/)
   d2 = (/Npoints,0/)
   call construct_var1d('cllcalipso', d3, d2, stlidar%cldlayer(:,1),v1d(1),units='%') ! correction sb
   call construct_var1d('clmcalipso', d3, d2, stlidar%cldlayer(:,2),v1d(2),units='%') ! correction sb
   call construct_var1d('clhcalipso', d3, d2, stlidar%cldlayer(:,3),v1d(3),units='%') ! correction sb
   call construct_var1d('cltcalipso', d3, d2, stlidar%cldlayer(:,4),v1d(4),units='%') ! correction sb
   call construct_var1d('cltlidarradar', d3, d2, stradar%radar_lidar_tcc,v1d(5),units='%')
   call construct_var1d('cltisccp', d3, d2, isccp%totalcldarea,v1d(6),units='%')
   call construct_var1d('pctisccp', d3, d2, isccp%meanptop,v1d(7),units='Pa')
   call construct_var1d('tauisccp', d3, d2, isccp%meantaucld,v1d(8),units='1')
   call construct_var1d('albisccp', d3, d2, isccp%meanalbedocld,v1d(9),units='1')
   call construct_var1d('meantbisccp', d3, d2, isccp%meantb,v1d(10),units='K')
   call construct_var1d('meantbclrisccp', d3, d2, isccp%meantbclr,v1d(11),units='K')
   call construct_var1d('cltmodis', d3, d2, modis%Cloud_Fraction_Total_Mean,v1d(12),units='%')
   call construct_var1d('clwmodis', d3, d2, modis%Cloud_Fraction_Water_Mean,v1d(13),units='%')
   call construct_var1d('climodis', d3, d2, modis%Cloud_Fraction_Ice_Mean,  v1d(14),units='%')
   call construct_var1d('clhmodis', d3, d2, modis%Cloud_Fraction_High_Mean,v1d(15),units='%')
   call construct_var1d('clmmodis', d3, d2, modis%Cloud_Fraction_Mid_Mean,v1d(16),units='%')
   call construct_var1d('cllmodis', d3, d2, modis%Cloud_Fraction_Low_Mean,  v1d(17),units='%')
   call construct_var1d('tautmodis', d3, d2, modis%Optical_Thickness_Total_Mean,v1d(18),units='1')
   call construct_var1d('tauwmodis', d3, d2, modis%Optical_Thickness_Water_Mean,v1d(19),units='1')
   call construct_var1d('tauimodis', d3, d2, modis%Optical_Thickness_Ice_Mean,v1d(20),units='1')
   call construct_var1d('tautlogmodis', d3, d2, modis%Optical_Thickness_Total_LogMean,v1d(21),units='1')
   call construct_var1d('tauwlogmodis', d3, d2, modis%Optical_Thickness_Water_LogMean,v1d(22),units='1')
   call construct_var1d('tauilogmodis', d3, d2, modis%Optical_Thickness_Ice_LogMean,v1d(23),units='1')
   call construct_var1d('reffclwmodis', d3, d2, modis%Cloud_Particle_Size_Water_Mean,v1d(24),units='m')
   call construct_var1d('reffclimodis', d3, d2, modis%Cloud_Particle_Size_Ice_Mean,  v1d(25),units='m')
   call construct_var1d('pctmodis', d3, d2, modis%Cloud_Top_Pressure_Total_Mean, v1d(26),units='Pa')
   call construct_var1d('lwpmodis', d3, d2, modis%Liquid_Water_Path_Mean, v1d(27),units='kg m-2')
   call construct_var1d('iwpmodis', d3, d2, modis%Ice_Water_Path_Mean,    v1d(28),units='kg m-2')
   call construct_var1d('toffset', d3, d2, gb%toffset,    v1d(29),units='day')
   ! 2D variables
   d4 = (/grid_id,height_axid,0,0/)
   d3 = (/Npoints,Nlvgrid,0/)
   call construct_var2d('clcalipso', d4, d3, stlidar%lidarcld,v2d(1),units='%')
   call construct_var2d('clcalipso2',  d4, d3, stradar%lidar_only_freq_cloud,v2d(2),units='%')
   d4 = (/grid_id,height_mlev_axid,0,0/)
   d3 = (/Npoints,Nlevels,0/)
   ! reshape   d4 = (/profile_axid,height_mlev_axid,time_axid,0/)
   call construct_var2d('lidarBetaMol532', d4, d3, sglidar%beta_mol,v2d(3),units='m-1 sr-1')
   d4 = (/grid_id,column_axid,0,0/)
   ! reshape d4 = (/profile_axid,column_axid,time_axid,0/)
   d3 = (/Npoints,Ncolumns,0/)
   call construct_var2d('boxtauisccp', d4, d3, isccp%boxtau,v2d(4),units='1')
   call construct_var2d('boxptopisccp', d4, d3, isccp%boxptop,v2d(5),units='Pa')
   d4 = (/grid_id,sza_axid,0,0/)
   d3 = (/Npoints,PARASOL_NREFL,0/)
   call construct_var2d('parasolRefl', d4, d3, stlidar%parasolrefl,v2d(6),units='1')
#ifdef RTTOV 
   !reshape d4 = (/profile_axid,channel_axid,time_axid,0/) 
   d4 = (/grid_id,channel_axid,0,0/) 
   d3 = (/Npoints,Nchannels,0/) 
   call construct_var2d('tbrttov', d4, d3, rttov%tbs,v2d(7),units='K') 
#endif

   ! 3D variables
   ! reshape d5 = (/profile_axid,column_axid,height_mlev_axid,time_axid,0/)
   d5 = (/grid_id,column_axid,height_mlev_axid,0,0/)
   d4 = (/Npoints,Ncolumns,Nlevels,0/)
   call construct_var3d('dbze94', d5, d4, sgradar%Ze_tot,v3d(1),units='1')
   call construct_var3d('atb532', d5, d4, sglidar%beta_tot,v3d(2),units='m-1 sr-1')
   call construct_var3d('fracout', d5, d4, sg%frac_out,v3d(3),units='1')
   ! reshape d5 = (/profile_axid,dbze_axid,height_axid,time_axid,0/)
   d5 = (/grid_id,dbze_axid,height_axid,0,0/)
   d4 = (/Npoints,DBZE_BINS,Nlvgrid,0/)
   call construct_var3d('cfadDbze94', d5, d4, stradar%cfad_ze,v3d(4),units='1')
   ! reshape d5 = (/profile_axid,sratio_axid,height_axid,time_axid,0/)
   d5 = (/grid_id,sratio_axid,height_axid,0,0/)
   d4 = (/Npoints,SR_BINS,Nlvgrid,0/)
   call construct_var3d('cfadLidarsr532', d5, d4, stlidar%cfad_sr,v3d(5),units='1')
   ! reshape d5 = (/profile_axid,tau_axid,pressure2_axid,time_axid,0/)
   d5 = (/grid_id,tau_axid,pressure2_axid,0,0/)
   d4 = (/Npoints,7,7,0/)
   call construct_var3d('clisccp', d5, d4, isccp%fq_isccp,v3d(6),units='%')
   call construct_var3d('clmodis', d5, d4, modis%Optical_Thickness_vs_Cloud_Top_Pressure, v3d(7), units='%')
   ! reshape d5 = (/profile_axid,tau_axid,MISR_CTH_axid,time_axid,0/)
   d5 = (/grid_id,tau_axid,MISR_CTH_axid,0,0/)
   d4 = (/Npoints,7,MISR_N_CTH,0/)
   call construct_var3d('clMISR', d5, d4, misr%fq_MISR,v3d(8),units='%')
  
   print *, '---------------Find list'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Find list of outputs to be written
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do i=1,N_OUT_LIST
     lfound = .false.
     if (trim(cfg%out_list(i)) /= '') then
        do j=1,n1d
          if (trim(v1d(j)%name) == trim(cfg%out_list(i))) then
            lout1d(j) = .true.
            lfound = .true.
            exit
          endif
        enddo
        if (.not.lfound) then
          do j=1,n2d
            if (trim(v2d(j)%name) == trim(cfg%out_list(i))) then
                lout2d(j) = .true.
                lfound = .true.
                exit
            endif
          enddo
        endif
        if (.not.lfound) then
          do j=1,n3d
            if (trim(v3d(j)%name) == trim(cfg%out_list(i))) then
                lout3d(j) = .true.
                lfound = .true.
                exit
            endif
          enddo
        endif
     endif
   enddo
   
   print *, '---------------Define variables'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define variables. Fill in dimensions table first if needed
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables
   do i=1,n1d
     if (lout1d(i)) var1d_id(i) = cmor_variable(table=table, table_entry=v1d(i)%name, units=v1d(i)%units, &
                                 axis_ids=v1d(i)%dimsid(1:1), missing_value=R_UNDEF)
   enddo
   ! 2D variables
   do i=1,n2d
     if (lout2d(i)) var2d_id(i) = cmor_variable(table=table, table_entry=v2d(i)%name, units=v2d(i)%units, &
                                 axis_ids=v2d(i)%dimsid(1:2), missing_value=R_UNDEF)
   enddo
   ! 3D variables
   do i=1,n3d
     if (lout3d(i)) then 
          var3d_id(i) = cmor_variable(table=trim(table), table_entry=v3d(i)%name, units=v3d(i)%units, &
                                  axis_ids=v3d(i)%dimsid(1:3), missing_value=R_UNDEF)
     endif
   enddo
   
   print *, '---------------Write variables'
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Write variables to file
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables
   do i=1,n1d
      if (lout1d(i)) then
        error_flag = cmor_write(var_id=var1d_id(i), data=v1d(i)%pntr, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v1d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=lonvar_id, data=gb%longitude,store_with=var1d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lon for '//trim(v1d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=latvar_id, data=gb%latitude,store_with=var1d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lat for '//trim(v1d(i)%name)
           stop
        endif
      endif
   enddo
   ! 2D variables
   do i=1,n2d
      if (lout2d(i)) then
        error_flag = cmor_write(var_id=var2d_id(i), data=v2d(i)%pntr, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v2d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=lonvar_id, data=gb%longitude,store_with=var2d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lon for '//trim(v2d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=latvar_id, data=gb%latitude,store_with=var2d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lat for '//trim(v2d(i)%name)
           stop
        endif
      endif
   enddo
   ! 3D variables
   do i=1,n3d
      if (lout3d(i)) then
        error_flag = cmor_write(var_id=var3d_id(i), data=v3d(i)%pntr, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v3d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=lonvar_id, data=gb%longitude,store_with=var3d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lon for '//trim(v3d(i)%name)
           stop
        endif
        error_flag = cmor_write(var_id=latvar_id, data=gb%latitude,store_with=var3d_id(i), &
                           ntimes_passed=1,time_vals=(/gb%time/))
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing lat for '//trim(v3d(i)%name)
           stop
        endif
      endif
   enddo
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Close files
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   error_flag = cmor_close()
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Deallocate memory
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   deallocate(profile_ax,column_ax,dbze_ax,dbze_bounds,sratio_bounds, &
             vgrid_bounds,mgrid_bounds,channel_ax)
   
  END SUBROUTINE NC_WRITE_COSP_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE NC_WRITE_COSP_2D ---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef RTTOV
  SUBROUTINE NC_WRITE_COSP_2D(cmor_nl,wmode,cfg,vgrid,gb,sg,sgradar,sglidar,isccp,misr,modis,rttov, &
                              stradar,stlidar,geomode,Nlon,Nlat)
#else
  SUBROUTINE NC_WRITE_COSP_2D(cmor_nl,wmode,cfg,vgrid,gb,sg,sgradar,sglidar,isccp,misr,modis, &
                              stradar,stlidar,geomode,Nlon,Nlat)
#endif
   ! Input arguments
   character(len=*),intent(in) :: cmor_nl
   character(len=*)  :: wmode ! Writing mode 'replace' or 'append'
   type(cosp_config),intent(in) :: cfg
   type(cosp_vgrid),intent(in) :: vgrid
   type(cosp_gridbox),intent(in) :: gb
   type(cosp_subgrid),intent(in) :: sg
   type(cosp_sgradar),intent(in) :: sgradar ! Subgrid radar
   type(cosp_sglidar),intent(in) :: sglidar ! Subgrid lidar
   type(cosp_isccp),intent(in)   :: isccp   ! ISCCP outputs
   type(cosp_misr),intent(in)    :: misr    ! MISR outputs
   type(cosp_modis),intent(in)   :: modis   ! MODIS outputs
#ifdef RTTOV
   type(cosp_rttov),intent(in)   :: rttov   ! RTTOV outputs
#endif
   type(cosp_radarstats),intent(in) :: stradar ! Summary statistics from radar simulator
   type(cosp_lidarstats),intent(in) :: stlidar ! Summary statistics from lidar simulator
   integer,intent(in) :: geomode,Nlon,Nlat
   !--- Local variables ---
   integer,parameter :: n3d = 8
#ifdef RTTOV 
   integer,parameter :: n2d = 7 
#else 
   integer,parameter :: n2d = 6 
#endif
   integer,parameter :: n1d = 28
   integer :: error_flag,i,j,Npoints,Ncolumns,Nlevels,Nlvgrid,maxtsteps,Nchannels
   logical :: lfound
   real :: lon_ax(Nlon),lat_ax(Nlat)
   character(len=512) :: inpath,outpath,start_date,model_id,experiment_id,institution,institute_id,source,calendar, &
                 contact,history,comment,table,parent_experiment_id,parent_experiment_rip,forcing
   character(len=2056) :: references
   integer :: realization, nc_action,initialization_method,physics_version
   double precision :: branch_time
   namelist/CMOR/inpath,outpath,start_date,model_id,experiment_id,branch_time,parent_experiment_id,parent_experiment_rip, &
              forcing,institution,institute_id,source,calendar,realization,initialization_method,physics_version, &
              contact,history,comment,references,table,maxtsteps
   real,dimension(:),allocatable :: column_ax,dbze_ax,channel_ax
   real,dimension(:,:),allocatable :: dbze_bounds,vgrid_bounds,sratio_bounds, &
                                      lon_bounds,lat_bounds,mgrid_bounds
   integer :: column_axid,height_axid,dbze_axid,height_mlev_axid,sratio_axid, &
              tau_axid,pressure2_axid,lon_axid,lat_axid,time_axid,sza_axid,MISR_CTH_axid, &
              channel_axid
   integer :: var3d_id(n3d),var2d_id(n2d),var1d_id(n1d)   
   type(var1d) :: v1d(n1d)
   type(var2d) :: v2d(n2d)
   type(var3d) :: v3d(n3d)
   integer :: d2(2),d3(3),d4(4),d5(5)
   logical,dimension(n1d) :: lout1d = .false.
   logical,dimension(n2d) :: lout2d = .false.
   logical,dimension(n3d) :: lout3d = .false.
   real,allocatable :: y2(:,:),y3(:,:,:),y4(:,:,:,:)
   double precision :: tbnds(2,1)
   character(len=64) :: pro_name = 'NC_WRITE_COSP_2D'
   
   Npoints   = gb%Npoints
   Ncolumns  = gb%Ncolumns
   Nlevels   = gb%Nlevels
   Nchannels = gb%Nchan
   Nlvgrid   = vgrid%Nlvgrid

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Safety checks
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (Npoints > Nlon*Nlat) then
      Npoints = Nlon*Nlat
      print *, ' -- '//trim(pro_name)//' Warning: Npoints > Nlon*Nlat'
   endif
   if (Npoints < Nlon*Nlat) then
      print *, ' -- '//trim(pro_name)//': only Npoints >= Nlon*Nlat is supported'
      stop
   endif
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Allocate memory and compute axes and bounds
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   tbnds(:,1) = gb%time_bnds
   allocate(column_ax(Ncolumns),dbze_ax(DBZE_BINS),channel_ax(Nchannels), &
            dbze_bounds(2,DBZE_BINS),vgrid_bounds(2,Nlvgrid),mgrid_bounds(2,Nlevels),sratio_bounds(2,SR_BINS), &
            lon_bounds(2,Nlon),lat_bounds(2,Nlat))
   ! Column
   do i=1,gb%Ncolumns
     column_ax(i) = i
   enddo
   ! Channels
   channel_ax = float(gb%ichan)
   ! Radar Ze
   do i=1,DBZE_BINS
     dbze_ax(i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 0.5)
   enddo
   do i=1,DBZE_BINS
     dbze_bounds(1,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 1)
     dbze_bounds(2,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*i
   enddo
   ! Height of model levels
   do i=1,Nlevels
     mgrid_bounds(1,i) = vgrid%mzl(i)
     mgrid_bounds(2,i) = vgrid%mzu(i)
   enddo
   ! Height of std grid
   do i=1,Nlvgrid
     vgrid_bounds(1,i) = vgrid%zl(i)
     vgrid_bounds(2,i) = vgrid%zu(i)
   enddo
   ! Lidar scattering ratio bounds (They are output by cosp_cfad_sr->diag_lidar in lmd_ipsl_stats.f90)
   sratio_bounds(2,:)         = stlidar%srbval(:) ! srbval contains the upper limits from lmd_ipsl_stats.f90
   sratio_bounds(1,2:SR_BINS) = stlidar%srbval(1:SR_BINS-1)
   sratio_bounds(1,1)         = 0.0
   sratio_bounds(2,SR_BINS)   = 1.e5 ! This matches with Chepfer et al., JGR, 2009. However, it is not consistent 
                                     ! with the upper limit in lmd_ipsl_stats.f90, which is LIDAR_UNDEF-1=998.999
   ! Lat lon axes
   if (geomode == 2) then
      lon_ax = gb%longitude(1:Nlon)
      lat_ax = gb%latitude(1:Npoints:Nlon)
   else if (geomode == 3) then
      lon_ax = gb%longitude(1:Npoints:Nlat)
      lat_ax = gb%latitude(1:Nlat)
   else if (geomode == 4) then
      lon_ax = gb%longitude(1:Nlon)
      lat_ax = gb%latitude(1:Nlat)
   endif
   lon_bounds(1,2:Nlon) = (lon_ax(1:Nlon-1) + lon_ax(2:Nlon))/2.0
   lon_bounds(1,1) = lon_ax(1) - (lon_bounds(1,2) - lon_ax(1))
   lon_bounds(2,1:Nlon-1) = lon_bounds(1,2:Nlon)
   lon_bounds(2,Nlon) = lon_ax(Nlon) + (lon_ax(Nlon) - lon_bounds(2,Nlon-1))
   lat_bounds(1,2:Nlat) = (lat_ax(1:Nlat-1) + lat_ax(2:Nlat))/2.0
   lat_bounds(1,1) = lat_ax(1) - (lat_bounds(1,2) - lat_ax(1))
   lat_bounds(2,1:Nlat-1) = lat_bounds(1,2:Nlat)
   lat_bounds(2,Nlat) = lat_ax(Nlat) + (lat_ax(Nlat) - lat_bounds(2,Nlat-1))

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Read namelist with information for CMOR output file
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   open(10,file=cmor_nl,status='old')
   read(10,nml=cmor)
   close(10)
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Specify path for tables and set up other CMOR options
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef USE_CMOR1
   error_flag = cmor_setup(inpath=trim(inpath),netcdf_file_action=trim(wmode))
#else
   nc_action = CMOR_APPEND_3
   if (trim(wmode) == 'replace') nc_action = CMOR_REPLACE_3
   error_flag = cmor_setup(inpath=trim(inpath),netcdf_file_action=nc_action,create_subdirectories=0)
#endif
      
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define dataset as output from COSP, and other model details
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef USE_CMOR1
   error_flag = cmor_dataset(outpath=trim(outpath),experiment_id=trim(experiment_id),institution=trim(institution), &
                    source=trim(source),calendar=trim(calendar),realization=realization,contact=trim(contact), &
                    history=trim(history),comment=trim(comment),references=trim(references))
#else
   error_flag = cmor_dataset(outpath=trim(outpath),experiment_id=trim(experiment_id),institution=trim(institution), &
                    source=trim(source),calendar=trim(calendar),realization=realization,contact=trim(contact), &
                    history=trim(history),comment=trim(comment),references=trim(references),model_id=trim(model_id), &
                    branch_time=branch_time,parent_experiment_id=trim(parent_experiment_id),forcing=trim(forcing), &
                    institute_id=trim(institute_id),parent_experiment_rip=trim(parent_experiment_rip), &
                    initialization_method=initialization_method,physics_version=physics_version)
   error_flag = cmor_set_cur_dataset_attribute('cosp_version',trim(COSP_VERSION))
#endif
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define axis
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   lon_axid     = cmor_axis(table=table, table_entry='longitude', units='degrees_east', length=Nlon, coord_vals=lon_ax, &
                        cell_bounds = lon_bounds)
   lat_axid     = cmor_axis(table=table, table_entry='latitude', units='degrees_north', length=Nlat, coord_vals=lat_ax, &
                        cell_bounds = lat_bounds)
   column_axid  = cmor_axis(table=table, table_entry='column', units='1', length=Ncolumns, coord_vals=column_ax)
   channel_axid = cmor_axis(table=table, table_entry='channel', units='1', length=Nchannels, coord_vals=channel_ax)
   height_axid  = cmor_axis(table=table, table_entry='alt40', units='m', length=Nlvgrid, &
                                                       coord_vals=vgrid%z,cell_bounds=vgrid_bounds)
   dbze_axid    = cmor_axis(table=table, table_entry='dbze', units='dBZ', length=DBZE_BINS, &
                            coord_vals=dbze_ax,cell_bounds=dbze_bounds)
   height_mlev_axid  = cmor_axis(table=table, table_entry='alevel', units='1', length=Nlevels, &
                                                       coord_vals=vgrid%mz,cell_bounds=mgrid_bounds)
   sratio_axid  = cmor_axis(table=table, table_entry='scatratio', units='1', length=SR_BINS, &
                            coord_vals=(sratio_bounds(1,:)+sratio_bounds(2,:))/2.0,cell_bounds=sratio_bounds)
   tau_axid     = cmor_axis(table=table, table_entry='tau', units='1', length=7, &
                            coord_vals=ISCCP_TAU,cell_bounds=ISCCP_TAU_BNDS)
   pressure2_axid = cmor_axis(table=table, table_entry='plev7', units='Pa', length=7, &
                            coord_vals=ISCCP_PC,cell_bounds=ISCCP_PC_BNDS)
   sza_axid   = cmor_axis(table=table, table_entry='sza5', units='degree', length=PARASOL_NREFL, coord_vals=PARASOL_SZA)
   MISR_CTH_axid = cmor_axis(table=table, table_entry='cth16', units='m', length=MISR_N_CTH, &
                            coord_vals=MISR_CTH,cell_bounds=MISR_CTH_BNDS)
   time_axid  = cmor_axis(table=table, table_entry='time1', units='days since '//trim(start_date), length=maxtsteps)
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Fill in variable info and associate pointers
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables
   d3 = (/lon_axid,lat_axid,time_axid/)
   d2 = (/Nlon,Nlat/)
   call construct_var1d('cllcalipso',     d3, d2, stlidar%cldlayer(:,1),v1d(1),units='%') ! correction sb
   call construct_var1d('clmcalipso',     d3, d2, stlidar%cldlayer(:,2),v1d(2),units='%') ! correction sb
   call construct_var1d('clhcalipso',     d3, d2, stlidar%cldlayer(:,3),v1d(3),units='%') ! correction sb
   call construct_var1d('cltcalipso',     d3, d2, stlidar%cldlayer(:,4),v1d(4),units='%') ! correction sb
   call construct_var1d('cltlidarradar',  d3, d2, stradar%radar_lidar_tcc,v1d(5),units='%')
   call construct_var1d('cltisccp',       d3, d2, isccp%totalcldarea,v1d(6),units='%')
   call construct_var1d('pctisccp',       d3, d2, isccp%meanptop,v1d(7),units='Pa')
   call construct_var1d('tauisccp',       d3, d2, isccp%meantaucld,v1d(8),units='1')
   call construct_var1d('albisccp',       d3, d2, isccp%meanalbedocld,v1d(9),units='1')
   call construct_var1d('meantbisccp',    d3, d2, isccp%meantb,v1d(10),units='K')
   call construct_var1d('meantbclrisccp', d3, d2, isccp%meantbclr,v1d(11),units='K')
   call construct_var1d('cltmodis', d3, d2, modis%Cloud_Fraction_Total_Mean,v1d(12),units='%')
   call construct_var1d('clwmodis', d3, d2, modis%Cloud_Fraction_Water_Mean,v1d(13),units='%')
   call construct_var1d('climodis', d3, d2, modis%Cloud_Fraction_Ice_Mean,  v1d(14),units='%')
   call construct_var1d('clhmodis', d3, d2, modis%Cloud_Fraction_High_Mean,v1d(15),units='%')
   call construct_var1d('clmmodis', d3, d2, modis%Cloud_Fraction_Mid_Mean,v1d(16),units='%')
   call construct_var1d('cllmodis', d3, d2, modis%Cloud_Fraction_Low_Mean,  v1d(17),units='%')
   call construct_var1d('tautmodis', d3, d2, modis%Optical_Thickness_Total_Mean,v1d(18),units='1')
   call construct_var1d('tauwmodis', d3, d2, modis%Optical_Thickness_Water_Mean,v1d(19),units='1')
   call construct_var1d('tauimodis', d3, d2, modis%Optical_Thickness_Ice_Mean,v1d(20),units='1')
   call construct_var1d('tautlogmodis', d3, d2, modis%Optical_Thickness_Total_LogMean,v1d(21),units='1')
   call construct_var1d('tauwlogmodis', d3, d2, modis%Optical_Thickness_Water_LogMean,v1d(22),units='1')
   call construct_var1d('tauilogmodis', d3, d2, modis%Optical_Thickness_Ice_LogMean,v1d(23),units='1')
   call construct_var1d('reffclwmodis', d3, d2, modis%Cloud_Particle_Size_Water_Mean,v1d(24),units='m')
   call construct_var1d('reffclimodis', d3, d2, modis%Cloud_Particle_Size_Ice_Mean,  v1d(25),units='m')
   call construct_var1d('pctmodis', d3, d2, modis%Cloud_Top_Pressure_Total_Mean, v1d(26),units='Pa')
   call construct_var1d('lwpmodis', d3, d2, modis%Liquid_Water_Path_Mean, v1d(27),units='kg m-2')
   call construct_var1d('iwpmodis', d3, d2, modis%Ice_Water_Path_Mean,    v1d(28),units='kg m-2')
   ! 2D variables
   d4 = (/lon_axid,lat_axid,height_axid,time_axid/)
   d3 = (/Nlon,Nlat,Nlvgrid/)
   call construct_var2d('clcalipso',  d4, d3, stlidar%lidarcld,v2d(1),units='%')
   call construct_var2d('clcalipso2', d4, d3, stradar%lidar_only_freq_cloud,v2d(2),units='%')
   d4 = (/lon_axid,lat_axid,height_mlev_axid,time_axid/)
   d3 = (/Nlon,Nlat,Nlevels/)
   call construct_var2d('lidarBetaMol532', d4, d3, sglidar%beta_mol,v2d(3),units='m-1 sr-1')
   d4 = (/lon_axid,lat_axid,column_axid,time_axid/)
   d3 = (/Nlon,Nlat,Ncolumns/)
   call construct_var2d('boxtauisccp',  d4, d3, isccp%boxtau,v2d(4),units='1')
   call construct_var2d('boxptopisccp', d4, d3, isccp%boxptop,v2d(5),units='Pa')
   d4 = (/lon_axid,lat_axid,sza_axid,time_axid/)
   d3 = (/Nlon,Nlat,PARASOL_NREFL/)
   call construct_var2d('parasolRefl', d4, d3, stlidar%parasolrefl,v2d(6),units='1')
#ifdef RTTOV 
   d4 = (/lon_axid,lat_axid,channel_axid,time_axid/) 
   d3 = (/Nlon,Nlat,Nchannels/) 
   call construct_var2d('tbrttov', d4, d3, rttov%tbs,v2d(7),units='K') 
#endif

   ! 3D variables
   d5 = (/lon_axid,lat_axid,column_axid,height_mlev_axid,time_axid/)
   d4 = (/Nlon,Nlat,Ncolumns,Nlevels/)
   call construct_var3d('dbze94', d5, d4, sgradar%Ze_tot,v3d(1),units='1')
   call construct_var3d('atb532', d5, d4, sglidar%beta_tot,v3d(2),units='m-1 sr-1')
   call construct_var3d('fracout', d5, d4, sg%frac_out,v3d(3),units='1')
   d5 = (/lon_axid,lat_axid,dbze_axid,height_axid,time_axid/)
   d4 = (/Nlon,Nlat,DBZE_BINS,Nlvgrid/)
   call construct_var3d('cfadDbze94', d5, d4, stradar%cfad_ze,v3d(4),units='1')
   d5 = (/lon_axid,lat_axid,sratio_axid,height_axid,time_axid/)
   d4 = (/Nlon,Nlat,SR_BINS,Nlvgrid/)
   call construct_var3d('cfadLidarsr532', d5, d4, stlidar%cfad_sr,v3d(5),units='1')
   d5 = (/lon_axid,lat_axid,tau_axid,pressure2_axid,time_axid/)
   d4 = (/Nlon,Nlat,7,7/)
   call construct_var3d('clisccp', d5, d4, isccp%fq_isccp,v3d(6),units='%')
   call construct_var3d('clmodis', d5, d4, modis%Optical_Thickness_vs_Cloud_Top_Pressure, v3d(7), units='%')
!    call construct_var3d('clmodis',                                                   &
!                         (/lon_axid, lat_axid, modis_tau_axid, modis_p_axid, time_axid/), &
!                         (/Nlon, Nlat, numModisTauBins, numModisPressureBins/), & 
!                         modis%Optical_Thickness_vs_Cloud_Top_Pressure, v3d(7), units='%')
   d5 = (/lon_axid,lat_axid,tau_axid,MISR_CTH_axid,time_axid/)
   d4 = (/Nlon,Nlat,7,MISR_N_CTH/)
   call construct_var3d('clMISR', d5, d4, misr%fq_MISR,v3d(8),units='%')
      
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Find list of outputs to be written
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do i=1,N_OUT_LIST
     lfound = .false.
     if (trim(cfg%out_list(i)) /= '') then
        do j=1,n1d
          if (trim(v1d(j)%name) == trim(cfg%out_list(i))) then
            lout1d(j) = .true.
            lfound = .true.
            exit
          endif
        enddo
        if (.not.lfound) then
          do j=1,n2d
            if (trim(v2d(j)%name) == trim(cfg%out_list(i))) then
                lout2d(j) = .true.
                lfound = .true.
                exit
            endif
          enddo
        endif
        if (.not.lfound) then
          do j=1,n3d
            if (trim(v3d(j)%name) == trim(cfg%out_list(i))) then
                lout3d(j) = .true.
                lfound = .true.
                exit
            endif
          enddo
        endif
     endif
   enddo
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Define variables. Fill in dimensions table first if needed
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables
   do i=1,n1d
     if (lout1d(i)) var1d_id(i) = cmor_variable(table=table, table_entry=v1d(i)%name, units=v1d(i)%units, &
                                 axis_ids=v1d(i)%dimsid, missing_value=R_UNDEF)
   enddo
   ! 2D variables
   do i=1,n2d
     if (lout2d(i)) var2d_id(i) = cmor_variable(table=table, table_entry=v2d(i)%name, units=v2d(i)%units, &
                                 axis_ids=v2d(i)%dimsid, missing_value=R_UNDEF)
   enddo
   ! 3D variables
   do i=1,n3d
     if (lout3d(i)) var3d_id(i) = cmor_variable(table=table, table_entry=v3d(i)%name, units=v3d(i)%units, &
                                  axis_ids=v3d(i)%dimsid, missing_value=R_UNDEF)
   enddo
   
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Write variables to file
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 1D variables (2D output)
   do i=1,n1d
      if (lout1d(i)) then
        allocate(y2(v1d(i)%dimssz(1),v1d(i)%dimssz(2)))
        call map_point_to_ll(Nlon,Nlat,geomode,x1=v1d(i)%pntr,y2=y2) ! Regridding
        error_flag = cmor_write(var_id=var1d_id(i), data=y2, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v1d(i)%name)
           stop
        endif
        deallocate(y2)
      endif
   enddo
   ! 2D variables (3D output)
   do i=1,n2d
      if (lout2d(i)) then
        allocate(y3(v2d(i)%dimssz(1),v2d(i)%dimssz(2),v2d(i)%dimssz(3)))
        call map_point_to_ll(Nlon,Nlat,geomode,x2=v2d(i)%pntr,y3=y3) ! Regridding
        error_flag = cmor_write(var_id=var2d_id(i), data=y3, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v2d(i)%name)
           stop
        endif
        deallocate(y3)
      endif
   enddo
   ! 3D variables (4D output)
   do i=1,n3d
      if (lout3d(i)) then
        allocate(y4(v3d(i)%dimssz(1),v3d(i)%dimssz(2),v3d(i)%dimssz(3),v3d(i)%dimssz(4)))
        call map_point_to_ll(Nlon,Nlat,geomode,x3=v3d(i)%pntr,y4=y4) ! Regridding
        error_flag = cmor_write(var_id=var3d_id(i), data=y4, &
                           ntimes_passed=1,time_vals=(/gb%time/),time_bnds=tbnds)
        if (error_flag < 0) then
           print *,  trim(pro_name)//': Error writing '//trim(v3d(i)%name)
           stop
        endif
        deallocate(y4)
      endif
   enddo
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Close files
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   error_flag = cmor_close()
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Deallocate memory
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   deallocate(column_ax,dbze_ax,channel_ax,dbze_bounds,sratio_bounds, &
             vgrid_bounds,mgrid_bounds,lon_bounds,lat_bounds)
   
  END SUBROUTINE NC_WRITE_COSP_2D
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE READ_COSP_OUTPUT_NL(cosp_nl,cfg)
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(out) :: cfg
  ! Local variables
  integer :: i
  logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, &
             Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis

  namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp, &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp, &
             Lcltisccp,Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis
#ifdef COSP_GFDL
  integer :: unit, io, ierr, logunit
#endif

  do i=1,N_OUT_LIST
    cfg%out_list(i)=''
  enddo
#ifdef COSP_GFDL
!---------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=cosp_output, iostat=io)
   ierr = check_nml_error(io,"cosp_output")
#else
!---------------------------------------------------------------------
    if ( file_exist('input.nml')) then
     unit =  open_namelist_file ()
     ierr=1; do while (ierr /= 0)
     read  (unit, nml=cosp_output, iostat=io, end=10)
     ierr = check_nml_error(io,'cosp_output')
     enddo
10      call close_file (unit)
   endif
#endif
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
     call write_version_number(version, tagname)
     logunit = stdlog()
     if (mpp_pe() == mpp_root_pe() )    &
                        write (logunit, nml=cosp_output)


#else
  open(10,file=cosp_nl,status='old')
  read(10,nml=cosp_output)
  close(10)
#endif
  
  ! Deal with dependencies
  if (.not.Lradar_sim) then
    LcfadDbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
    Ldbze94        = .false.
    Lclcalipso2    = .false. ! Needs radar & lidar
  endif
  if (.not.Llidar_sim) then
    Latb532 = .false.
    LcfadLidarsr532 = .false.
    Lclcalipso2      = .false. ! Needs radar & lidar
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    Lcltlidarradar   = .false.
    LparasolRefl    = .false.
    LlidarBetaMol532 = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
  endif
  if (.not.Lisccp_sim) then
    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp        = .false.
    Lpctisccp       = .false.
    Ltauisccp       = .false.
    Lcltisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.
  endif
  if (.not.Lmisr_sim) then
    LclMISR = .false.
  endif
  if (.not.Lrttov_sim) then
    Ltbrttov = .false.
  endif
  if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
      (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
    Lfracout = .false.
  endif
  if (.not.Lmodis_sim) then
    Lcltmodis=.false.
    Lclwmodis=.false.
    Lclimodis=.false.
    Lclhmodis=.false.
    Lclmmodis=.false.
    Lcllmodis=.false.
    Ltautmodis=.false.
    Ltauwmodis=.false.
    Ltauimodis=.false.
    Ltautlogmodis=.false.
    Ltauwlogmodis=.false.
    Ltauilogmodis=.false.
    Lreffclwmodis=.false.
    Lreffclimodis=.false.
    Lpctmodis=.false.
    Llwpmodis=.false.
    Liwpmodis=.false.
    Lclmodis=.false.
  endif
  if (Lmodis_sim) Lisccp_sim = .true.
  
  cfg%Lstats = .false.
  if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.
  
  ! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmodis_sim = Lmodis_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim
  
  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
    cfg%Lwrite_output = .true.
  endif
  
  ! Output diagnostics
  i = 1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Latb532)          cfg%out_list(i) = 'atb532'
  i = i+1
  if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (LcfadDbze94)      cfg%out_list(i) = 'cfadDbze94'
  i = i+1
  if (LcfadLidarsr532)  cfg%out_list(i) = 'cfadLidarsr532'
  i = i+1
  if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lclisccp)         cfg%out_list(i) = 'clisccp'
  i = i+1
  if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
  i = i+1
  if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lpctisccp)        cfg%out_list(i) = 'pctisccp'
  i = i+1
  if (Ldbze94)          cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Lcltisccp)        cfg%out_list(i) = 'cltisccp'
  i = i+1
  if (Ltoffset)         cfg%out_list(i) = 'toffset'
  i = i+1
  if (LparasolRefl)     cfg%out_list(i) = 'parasolRefl'
  i = i+1
  if (LclMISR)          cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (Lfracout)         cfg%out_list(i) = 'fracout'
  i = i+1
  if (LlidarBetaMol532) cfg%out_list(i) = 'lidarBetaMol532'
  i = i+1
  if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'
  i = i+1
  if (Lcltmodis)        cfg%out_list(i) = 'cltmodis'
  i = i+1
  if (Lclwmodis)        cfg%out_list(i) = 'clwmodis'
  i = i+1
  if (Lclimodis)        cfg%out_list(i) = 'climodis'
  i = i+1
  if (Lclhmodis)        cfg%out_list(i) = 'clhmodis'
  i = i+1
  if (Lclmmodis)        cfg%out_list(i) = 'clmmodis'
  i = i+1
  if (Lcllmodis)        cfg%out_list(i) = 'cllmodis'
  i = i+1
  if (Ltautmodis)       cfg%out_list(i) = 'tautmodis'
  i = i+1
  if (Ltauwmodis)       cfg%out_list(i) = 'tauwmodis'
  i = i+1
  if (Ltauimodis)       cfg%out_list(i) = 'tauimodis'
  i = i+1
  if (Ltautlogmodis)    cfg%out_list(i) = 'tautlogmodis'
  i = i+1
  if (Ltauwlogmodis)    cfg%out_list(i) = 'tauwlogmodis'
  i = i+1
  if (Ltauilogmodis)    cfg%out_list(i) = 'tauilogmodis'
  i = i+1
  if (Lreffclwmodis)    cfg%out_list(i) = 'reffclwmodis'
  i = i+1
  if (Lreffclimodis)    cfg%out_list(i) = 'reffclimodis'
  i = i+1
  if (Lpctmodis)        cfg%out_list(i) = 'pctmodis'
  i = i+1
  if (Llwpmodis)        cfg%out_list(i) = 'lwpmodis'
  i = i+1
  if (Liwpmodis)        cfg%out_list(i) = 'iwpmodis'
  i = i+1
  if (Lclmodis)         cfg%out_list(i) = 'clmodis'
    
  if (i /= N_OUT_LIST) then
#ifdef COSP_GFDL
       call error_mesg ('cosp_io:read_cosp_output_nl',  &
              'COSP_IO: wrong number of output diagnostics', &
                                                            FATAL)
#else
     print *, 'COSP_IO: wrong number of output diagnostics'
     stop
#endif
  endif

  ! Copy diagnostic flags to cfg structure
  ! ISCCP simulator  
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lclisccp = Lclisccp
  cfg%Lpctisccp = Lpctisccp
  cfg%Ltauisccp = Ltauisccp
  cfg%Lcltisccp = Lcltisccp
  ! CloudSat simulator  
  cfg%Ldbze94 = Ldbze94
  cfg%LcfadDbze94 = LcfadDbze94
  ! CALIPSO/PARASOL simulator  
  cfg%LcfadLidarsr532 = LcfadLidarsr532
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lcltlidarradar = Lcltlidarradar
  cfg%LparasolRefl = LparasolRefl
  ! MISR simulator  
  cfg%LclMISR = LclMISR
  ! Other
  cfg%Ltoffset = Ltoffset
  cfg%Lfracout = Lfracout
  cfg%LlidarBetaMol532 = LlidarBetaMol532
  ! RTTOV
  cfg%Ltbrttov = Ltbrttov
  ! MODIS simulator  
  cfg%Lcltmodis=Lcltmodis
  cfg%Lclwmodis=Lclwmodis
  cfg%Lclimodis=Lclimodis
  cfg%Lclhmodis=Lclhmodis
  cfg%Lclmmodis=Lclmmodis
  cfg%Lcllmodis=Lcllmodis
  cfg%Ltautmodis=Ltautmodis
  cfg%Ltauwmodis=Ltauwmodis
  cfg%Ltauimodis=Ltauimodis
  cfg%Ltautlogmodis=Ltautlogmodis
  cfg%Ltauwlogmodis=Ltauwlogmodis
  cfg%Ltauilogmodis=Ltauilogmodis
  cfg%Lreffclwmodis=Lreffclwmodis
  cfg%Lreffclimodis=Lreffclimodis
  cfg%Lpctmodis=Lpctmodis
  cfg%Llwpmodis=Llwpmodis
  cfg%Liwpmodis=Liwpmodis
  cfg%Lclmodis=Lclmodis
 END SUBROUTINE READ_COSP_OUTPUT_NL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE ERROR_CONTROL ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ERROR(routine_name,message,errcode) 
    character(len = *), intent(in) :: routine_name
    character(len = *), intent(in) :: message
    integer,optional :: errcode
    
    write(6, *) " ********** Failure in ", trim(routine_name)
    write(6, *) " ********** ", trim(message)
    if (present(errcode)) write(6, *) " ********** errcode: ", errcode
    flush(6)
    stop
  END SUBROUTINE COSP_ERROR
    
END MODULE MOD_COSP_IO

#include "cosp_defs.H"
#ifdef COSP_GFDL

!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------

! $Id: load_mie_table.F90,v 20.0 2013/12/13 23:16:47 fms Exp $
! $Name: tikal $
! cosp_version = 1.3.2

#endif

  subroutine load_mie_table(mie_table_name,mt)
  use radar_simulator_types
#ifdef COSP_GFDL
  use mpp_mod,only: get_unit
#endif
  implicit none
  
! Purpose:
!   Loads the Mie table data
!   Part of Quickbeam v1.03
!   http://reef.atmos.colostate.edu/haynes/radarsim
!
! Inputs:  
!   [mie_table_name]   Mie table input file
!
! Outputs:
!   [mt]            structure of Mie table data
!
! Created from Quickbeam v1.02 08/24/2006 by Roger Marchand  

! ----- INPUT -----
  character*200, intent(in) :: mie_table_name

! ----- OUTPUT -----
  type(mie), intent(out) :: mt

! ----- INTERNAL -----  
#ifdef COSP_GFDL
  integer :: i, funit
#else
  integer :: i
#endif

  integer*4 :: dummy_in(4)
     
#ifdef COSP_GFDL
    funit = get_unit()
    open(funit,file=mie_table_name,action='read')
    read(funit,*) dummy_in
#else
    open(51,file=mie_table_name,action='read')
    read(51,*) dummy_in
#endif
 
     if(dummy_in(1).ne. mt_nfreq .or. &
        dummy_in(2).ne. mt_ntt .or. &
        dummy_in(3).ne. mt_nf .or. &
        dummy_in(4).ne. mt_nd) then

          print *,'Mie file is of size :',dummy_in(:)
          print *,'  expected a size of:',mt_nfreq, mt_ntt,mt_nf,mt_nf
          print *,'  change paramters in radar_simulator_types.f90 ?? '
          stop
     endif

#ifdef COSP_GFDL
    read(funit,*) mt%freq
    read(funit,*) mt%tt
    read(funit,*) mt%f
    read(funit,*) mt%phase
    read(funit,*) mt%D
    read(funit,*) mt%qext
    read(funit,*) mt%qbsca
#else
    read(51,*) mt%freq
    read(51,*) mt%tt
    read(51,*) mt%f
    read(51,*) mt%phase
    read(51,*) mt%D
    read(51,*) mt%qext
    read(51,*) mt%qbsca
#endif
    
#ifdef COSP_GFDL
    close(funit)
#else
    close(51)
#endif

! // create arrays of liquid/ice temperature
  cnt_liq = 0
  cnt_ice = 0
  do i=1,mt_ntt
    if (mt%phase(i) == 0) cnt_liq = cnt_liq + 1
    if (mt%phase(i) == 1) cnt_ice = cnt_ice + 1
  enddo
#ifndef COSP_GFDL
  allocate(mt_ttl(cnt_liq),mt_tti(cnt_ice))
#endif
  do i=1,cnt_liq
    mt_ttl(i) = mt%tt(i)
  enddo
  do i=1,cnt_ice
    mt_tti(i) = mt%tt(cnt_liq+i)
  enddo

  end subroutine load_mie_table

#include "cosp_defs.H"
#ifdef COSP_GFDL

!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------

! $Id: cosp_rttov_simulator.F90,v 20.0 2013/12/13 23:15:49 fms Exp $
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
! Aug 2008 - V. John - Initial version
! Feb 2009 - V. John - Trace gases and max number of profiles
!


!#include "cosp_defs.h"
#ifndef COSP_GFDL
#include "cosp_defs.h"
#endif
MODULE MOD_COSP_RTTOV_SIMULATOR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
#ifdef RTTOV
  USE MOD_COSP_RTTOV
#endif
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE COSP_RTTOV_SIMULATOR ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RTTOV_SIMULATOR(gbx,y)
  
  ! Arguments
  type(cosp_gridbox),intent(in)  :: gbx ! Gridbox info
  type(cosp_rttov),intent(inout) :: y   ! RTTOV output
  
  ! some local variables for profile conversions etc.
  real, parameter :: eps    =  0.622
  real, parameter :: Mdry   =  28.966
  real, parameter :: Mo3    =  47.9983
  real, parameter :: Mco2   =  44.0096
  real, parameter :: Mch4   =  16.0426
  real, parameter :: Mn2o   =  44.0129
  real, parameter :: Mco    =  28.0102
  integer, parameter :: MaxLim  =  100
  
  ! Local variables 
  integer :: Npoints
  real :: sh(gbx%Npoints, gbx%Nlevels)
  real :: pp(gbx%Npoints, gbx%Nlevels)
  real :: tt(gbx%Npoints, gbx%Nlevels)
  real :: o3(gbx%Npoints, gbx%Nlevels)

  real :: co2,ch4,n2o,co
  real :: tt_surf(gbx%Npoints) ! 1.5 m T
  real :: sh_surf(gbx%Npoints) ! 1.5 m q 
  integer :: nloop,rmod,il
  integer :: istart,istop
  integer :: nprof,nlevels
    
  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Reverting Levels from TOA to surface
  sh  = gbx%sh(:,Nlevels:1:-1) 
  pp  = gbx%p(:,Nlevels:1:-1) / 100.
  tt  = gbx%t(:,Nlevels:1:-1) 
  o3  = gbx%mr_ozone(:,Nlevels:1:-1)
  
  ! FIXME: 1.5 m T and q should be added to input
  tt_surf  =  tt(:, Nlevels)
  sh_surf  =  sh(:, Nlevels)
  
  !Converting Specific Humidity to PPMV
  sh  =  ( sh / ( sh + eps * ( 1. - sh ) ) ) * 1e6

  !Converting Mass mixing ratio of other trace gases to ppmv
  o3   =  ( Mdry / Mo3  ) *     o3  * 1e6
  co2  =  ( Mdry / Mco2 ) * gbx%co2 * 1e6
  ch4  =  ( Mdry / Mch4 ) * gbx%ch4 * 1e6  
  n2o  =  ( Mdry / Mn2o ) * gbx%n2o * 1e6
  co   =  ( Mdry / Mco  ) * gbx%co  * 1e6
  
  !! RTTOV can handle only about 100 profiles at a time (FIXME: Check this with Roger) 
  !! So we are putting a loop of 100 
  
  nloop  =  Npoints / MaxLim
  rmod   =  MOD( Npoints, MaxLim )
  
  if( rmod .ne. 0 ) then  
     nloop = nloop + 1
  endif
  
  !! looping over MaxLim number of profiles
  do il = 1, nloop
     istart  =  (il - 1) * MaxLim + 1
     istop   =  min(il * MaxLim, Npoints) 
     
     if( ( il .eq. nloop ) .and. ( rmod .ne. 0 ) ) then
        nprof   =  rmod
     else
        nprof   =  MaxLim
     endif
          
#ifdef RTTOV
     call  rttov_multprof(              &
          gbx%Nchan,                    &
          gbx%ichan,                    &
          gbx%surfem,                   &
          nprof,                        &
          Nlevels,                      &
          gbx%Plat,                     &
          gbx%Sat,                      &
          gbx%Inst,                     &
          gbx%ZenAng,                   &
          pp(istart:istop, :),          &
          tt(istart:istop, :),          &
          sh(istart:istop, :),          &
          o3(istart:istop, :),          &
          co2,                          &
          ch4,                          &
          n2o,                          &
          co,                           &
          gbx%sfc_height(istart:istop), &
          gbx%u_wind(istart:istop),     &
          gbx%v_wind(istart:istop),     &
          gbx%skt(istart:istop),        &
          gbx%psfc(istart:istop)/100.,  &
          tt_surf(istart:istop),        &
          sh_surf(istart:istop),        &
          gbx%land(istart:istop),       &
          gbx%latitude(istart:istop),   &
          y%tbs(istart:istop, :) )
#endif
  enddo
  
END SUBROUTINE COSP_RTTOV_SIMULATOR

END MODULE MOD_COSP_RTTOV_SIMULATOR

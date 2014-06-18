#include "cosp_defs.H"
#ifdef COSP_GFDL

!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------

! $Id: icarus.F90,v 20.0 2013/12/13 23:16:02 fms Exp $
! $Name: tikal $
! cosp_version = 1.3.2

#endif

#ifdef COSP_GFDL
      SUBROUTINE ICARUS(          &
     &     debug,                 &
     &     debugcol,              &
     &     npoints,               &
     &     sunlit,                &
     &     nlev,                  &
     &     ncol,                  &
     &     pfull,                 &
     &     phalf,                 &
     &     qv,                    &
     &     cc,                    &
     &     conv,                  &
     &     dtau_s,                &
     &     dtau_c,                &
     &     top_height,            &
     &     top_height_direction,  &
     &     overlap,               &
     &     frac_out,              &
     &     skt,                   &
     &     emsfc_lw,              &
     &     at,                    &
     &     dem_s,                 &
     &     dem_c,                 &
     &     fq_isccp,              &
     &     totalcldarea,          &
     &     meanptop,              &
     &     meantaucld,            &
     &     meanalbedocld,         &
     &     meantb,                &
     &     meantbclr,             &
     &     boxtau,                &
     &     boxptop,               &
     &     dtau_col,              &
     &     dem_col                &
     &)
#else
      SUBROUTINE ICARUS(
     &     debug,
     &     debugcol,
     &     npoints,
     &     sunlit,
     &     nlev,
     &     ncol,
     &     pfull,
     &     phalf,
     &     qv,
     &     cc,
     &     conv,
     &     dtau_s,
     &     dtau_c,
     &     top_height,
     &     top_height_direction,
     &     overlap,
     &     frac_out,
     &     skt,
     &     emsfc_lw,
     &     at,
     &     dem_s,
     &     dem_c,
     &     fq_isccp,
     &     totalcldarea,
     &     meanptop,
     &     meantaucld,
     &     meanalbedocld,
     &     meantb,
     &     meantbclr,
     &     boxtau,
     &     boxptop
     &)
#endif

#ifndef COSP_GFDL
!$Id: icarus.F90,v 20.0 2013/12/13 23:16:02 fms Exp $
#endif

! *****************************COPYRIGHT****************************
! (c) 2009, Lawrence Livermore National Security Limited Liability 
! Corporation.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security 
!       Limited Liability Corporation nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************
! *****************************COPYRIGHT*******************************

#ifdef COSP_GFDL
use mpp_mod,only: get_unit
use fms_mod,only: stdlog, error_mesg, FATAL
#endif
      implicit none

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolprint
      
!     -----
!     Input 
!     -----

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER sunlit(npoints) !  1 for day points, 0 for night time

      REAL pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL cc(npoints,nlev)   
                  !  input cloud cover in each model level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds

      REAL conv(npoints,nlev) 
                  !  input convective cloud cover in each model
                  !   level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by convective clouds

      REAL dtau_s(npoints,nlev) 
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box

      REAL dtau_c(npoints,nlev) 
                  !  mean 0.67 micron optical depth of convective
                !  clouds in each
                  !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                              !  optical depth to adjust cloud top pressure. Note
                              !  that this calculation is most appropriate to compare
                              !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                              !  3 = adjust top height using only the computed
                              !  infrared brightness temperature. Note that this
                              !  calculation is most appropriate to compare to ISCCP
                              !  IR only algortihm (i.e. you can compare to nighttime
                              !  ISCCP data with this option)

      INTEGER top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 1 = find the *lowest* altitude (highest pressure) level
				 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 2 = find the *highest* altitude (lowest pressure) level
				 ! with interpolated temperature equal to the radiance 
				 ! determined cloud-top temperature
				 ! 
				 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
				 !				 !
				 ! 1 = old setting: matches all versions of 
				 ! ISCCP simulator with versions numbers 3.5.1 and lower
				 !
				 ! 2 = default setting: for version numbers 4.0 and higher
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL skt(npoints)                 !  skin Temperature (K)
      REAL emsfc_lw                     !  10.5 micron emissivity of surface (fraction)                                            
      REAL at(npoints,nlev)                   !  temperature in each model level (K)
      REAL dem_s(npoints,nlev)                !  10.5 micron longwave emissivity of stratiform
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
      REAL dem_c(npoints,nlev)                  !  10.5 micron longwave emissivity of convective
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column

#ifdef COSP_GFDL
       REAL, optional :: dtau_col(npoints,ncol,nlev)
                              ! tau values obtained from model
                              ! stochastic columns

       REAL, optional :: dem_col(npoints,ncol,nlev)
                              ! lw emissivity values obtained
                              ! from model stochastic columns
 

#endif


!     ------
!     Output
!     ------

      REAL fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  NOTE: This diagnostic
					! does not count model clouds with tau < isccp_taumin
                              ! Thus this diagnostic does not equal the sum over all entries of fq_isccp.
			      ! However, this diagnostic does equal the sum over entries of fq_isccp with
			      ! itau = 2:7 (omitting itau = 1)
      
      
      ! The following three means are averages only over the cloudy areas with tau > isccp_taumin.  
      ! If no clouds with tau > isccp_taumin are in grid box all three quantities should equal zero.      
                              
      REAL meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                              
      REAL meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
      real meanalbedocld(npoints)        ! mean cloud albedo
                                        ! linear averaging in albedo performed
					
      real meantb(npoints)              ! mean all-sky 10.5 micron brightness temperature
      
      real meantbclr(npoints)           ! mean clear-sky 10.5 micron brightness temperature
      
      REAL boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column
                              
                                                                                          
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
                              !  emissivity in part of
                              !  gridbox under consideration

      REAL ptrop(npoints)
      REAL attrop(npoints)
      REAL attropmin (npoints)
      REAL atmax(npoints)
      REAL btcmin(npoints)
      REAL transmax(npoints)

      INTEGER i,j,ilev,ibox,itrop(npoints)
      INTEGER ipres(npoints)
      INTEGER itau(npoints),ilev2
      INTEGER acc(nlev,ncol)
      INTEGER match(npoints,nlev-1)
      INTEGER nmatch(npoints)
      INTEGER levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
      real fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      real taumin(npoints)
      real dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      real press(npoints), dpress(npoints), atmden(npoints)
      real rvh20(npoints), wk(npoints), rhoave(npoints)
      real rh20s(npoints), rfrgn(npoints)
      real tmpexp(npoints),tauwv(npoints)
      
#ifdef COSP_GFDL
      character(len=1) :: cchar(6),cchar_realtops(6)
#else
      character*1 cchar(6),cchar_realtops(6)
#endif
      integer icycle
      REAL tau(npoints,ncol)
      LOGICAL box_cloudy(npoints,ncol)
      REAL tb(npoints,ncol)
      REAL ptop(npoints,ncol)
      REAL emcld(npoints,ncol)
      REAL fluxtop(npoints,ncol)
      REAL trans_layers_above(npoints,ncol)
      real isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      REAL albedocld(npoints,ncol)
      real boxarea
      integer debug       ! set to non-zero value to print out inputs
                    ! with step debug
      integer debugcol    ! set to non-zero value to print out column
                    ! decomposition with step debugcol
      integer rangevec(npoints),rangeerror

      integer index1(npoints),num1,jj,k1,k2
#ifdef COSP_GFDL
      integer   funit, logunit
#endif
      real rec2p13,tauchk,logp,logp1,logp2,atd
      real output_missing_value

#ifdef COSP_GFDL
      character(len=10) :: ftn09
#else
      character*10 ftn09
#endif
      
      DATA isccp_taumin / 0.3 /
      DATA output_missing_value / -1.E+30 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!     ------ End duplicate definitions common to wrapper routine

      tauchk = -1.*log(0.9999999)
      rec2p13=1./2.13

      ncolprint=0

#ifdef COSP_GFDL
      logunit = stdlog()
#endif
      if ( debug.ne.0 ) then
          j=1
#ifdef COSP_GFDL
          write(logunit,'(a10)') 'j='
          write(logunit,'(8I10)') j
          write(logunit,'(a10)') 'debug='
          write(logunit,'(8I10)') debug
          write(logunit,'(a10)') 'debugcol='
          write(logunit,'(8I10)') debugcol
          write(logunit,'(a10)') 'npoints='
          write(logunit,'(8I10)') npoints
          write(logunit,'(a10)') 'nlev='
          write(logunit,'(8I10)') nlev
          write(logunit,'(a10)') 'ncol='
          write(logunit,'(8I10)') ncol
          write(logunit,'(a11)') 'top_height='
          write(logunit,'(8I10)') top_height
          write(logunit,'(a21)') 'top_height_direction='
          write(logunit,'(8I10)') top_height_direction
          write(logunit,'(a10)') 'overlap='
          write(logunit,'(8I10)') overlap
          write(logunit,'(a10)') 'emsfc_lw='
          write(logunit,'(8f10.2)') emsfc_lw
        do j=1,npoints,debug
          write(logunit,'(a10)') 'j='
          write(logunit,'(8I10)') j
          write(logunit,'(a10)') 'sunlit='
          write(logunit,'(8I10)') sunlit(j)
          write(logunit,'(a10)') 'pfull='
          write(logunit,'(8f10.2)') (pfull(j,i),i=1,nlev)
          write(logunit,'(a10)') 'phalf='
          write(logunit,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
          write(logunit,'(a10)') 'qv='
          write(logunit,'(8f10.3)') (qv(j,i),i=1,nlev)
          write(logunit,'(a10)') 'cc='
          write(logunit,'(8f10.3)') (cc(j,i),i=1,nlev)
          write(logunit,'(a10)') 'conv='
          write(logunit,'(8f10.2)') (conv(j,i),i=1,nlev)
          write(logunit,'(a10)') 'dtau_s='
          write(logunit,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
          write(logunit,'(a10)') 'dtau_c='
          write(logunit,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
          write(logunit,'(a10)') 'skt='
          write(logunit,'(8f10.2)') skt(j)
          write(logunit,'(a10)') 'at='
          write(logunit,'(8f10.2)') (at(j,i),i=1,nlev)
          write(logunit,'(a10)') 'dem_s='
          write(logunit,'(8f10.3)') (dem_s(j,i),i=1,nlev)
          write(logunit,'(a10)') 'dem_c='
          write(logunit,'(8f10.3)') (dem_c(j,i),i=1,nlev)
        enddo
#else
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'debug='
          write(6,'(8I10)') debug
          write(6,'(a10)') 'debugcol='
          write(6,'(8I10)') debugcol
          write(6,'(a10)') 'npoints='
          write(6,'(8I10)') npoints
          write(6,'(a10)') 'nlev='
          write(6,'(8I10)') nlev
          write(6,'(a10)') 'ncol='
          write(6,'(8I10)') ncol
          write(6,'(a11)') 'top_height='
          write(6,'(8I10)') top_height
	  write(6,'(a21)') 'top_height_direction='
          write(6,'(8I10)') top_height_direction
          write(6,'(a10)') 'overlap='
          write(6,'(8I10)') overlap
          write(6,'(a10)') 'emsfc_lw='
          write(6,'(8f10.2)') emsfc_lw
        do j=1,npoints,debug
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'sunlit='
          write(6,'(8I10)') sunlit(j)
          write(6,'(a10)') 'pfull='
          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
          write(6,'(a10)') 'phalf='
          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
          write(6,'(a10)') 'qv='
          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
          write(6,'(a10)') 'cc='
          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
          write(6,'(a10)') 'conv='
          write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_s='
          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_c='
          write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
          write(6,'(a10)') 'skt='
          write(6,'(8f10.2)') skt(j)
          write(6,'(a10)') 'at='
          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_s='
          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_c='
          write(6,'(8f10.3)') (dem_c(j,i),i=1,nlev)
        enddo
#endif
      endif

!     ---------------------------------------------------!

      if (ncolprint.ne.0) then
      do j=1,npoints,1000
#ifdef COSP_GFDL
        write(logunit,'(a10)') 'j='
        write(logunit,'(8I10)') j
#else
        write(6,'(a10)') 'j='
        write(6,'(8I10)') j
#endif
      enddo
      endif

      if (top_height .eq. 1 .or. top_height .eq. 3) then 

      do j=1,npoints 
          ptrop(j)=5000.
          attropmin(j) = 400.
          atmax(j) = 0.
          attrop(j) = 120.
          itrop(j) = 1
      enddo 

      do 12 ilev=1,nlev
        do j=1,npoints 
#ifdef COSP_GFDL
          if (pfull(j,ilev) .lt. 40000. .and. &
     &          pfull(j,ilev) .gt.  5000. .and. &
     &          at(j,ilev) .lt. attropmin(j)) then
#else
         if (pfull(j,ilev) .lt. 40000. .and.
     &          pfull(j,ilev) .gt.  5000. .and.
     &          at(j,ilev) .lt. attropmin(j)) then
#endif
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           end if
        enddo
12    continue

      do 13 ilev=1,nlev
        do j=1,npoints 
#ifdef COSP_GFDL
           if (at(j,ilev) .gt. atmax(j) .and.    &
     &              ilev  .ge. itrop(j)) atmax(j)=at(j,ilev)
#else
           if (at(j,ilev) .gt. atmax(j) .and.
     &              ilev  .ge. itrop(j)) atmax(j)=at(j,ilev)
#endif
        enddo
13    continue

      end if


      if (top_height .eq. 1 .or. top_height .eq. 3) then
          do j=1,npoints
              meantb(j) = 0.
	      meantbclr(j) = 0. 
          end do
      else
          do j=1,npoints
              meantb(j) = output_missing_value
       	      meantbclr(j) = output_missing_value
          end do
      end if
      
!     -----------------------------------------------------!

!     ---------------------------------------------------!

      do ilev=1,nlev
        do j=1,npoints

          rangevec(j)=0

          if (cc(j,ilev) .lt. 0. .or. cc(j,ilev) .gt. 1.) then
!           error = cloud fraction less than zero
!           error = cloud fraction greater than 1
            rangevec(j)=rangevec(j)+1
          endif

          if (conv(j,ilev) .lt. 0. .or. conv(j,ilev) .gt. 1.) then
!           ' error = convective cloud fraction less than zero'
!           ' error = convective cloud fraction greater than 1'
            rangevec(j)=rangevec(j)+2
          endif

          if (dtau_s(j,ilev) .lt. 0.) then
!           ' error = stratiform cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+4
          endif

          if (dtau_c(j,ilev) .lt. 0.) then
!           ' error = convective cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+8
          endif

          if (dem_s(j,ilev) .lt. 0. .or. dem_s(j,ilev) .gt. 1.) then
!             ' error = stratiform cloud emissivity less than zero'
!             ' error = stratiform cloud emissivity greater than 1'
            rangevec(j)=rangevec(j)+16
          endif

          if (dem_c(j,ilev) .lt. 0. .or. dem_c(j,ilev) .gt. 1.) then
!             ' error = convective cloud emissivity less than zero'
!             ' error = convective cloud emissivity greater than 1'
              rangevec(j)=rangevec(j)+32
          endif
        enddo

        rangeerror=0
        do j=1,npoints
            rangeerror=rangeerror+rangevec(j)
        enddo

#ifdef COSP_GFDL
        if (rangeerror.ne.0) then
              write (logunit,*) 'Input variable out of range'
              write (logunit,*) 'rangevec:'
              write (logunit,*) rangevec
              call error_mesg('ICARUS','Input variable out of range',FATAL)

        endif
#else
        if (rangeerror.ne.0) then 
              write (6,*) 'Input variable out of range'
              write (6,*) 'rangevec:'
              write (6,*) rangevec
              call flush(6)
              STOP
        endif
#endif
      enddo

!
!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
        do j=1,npoints 
            tau(j,ibox)=0.
          albedocld(j,ibox)=0.
          boxtau(j,ibox)=output_missing_value
          boxptop(j,ibox)=output_missing_value
          box_cloudy(j,ibox)=.false.
        enddo
15    continue

      !compute total cloud optical depth for each column     
#ifdef COSP_GFDL
    if (present(dtau_col)     ) then
        do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints
                tau(j,ibox)=tau(j,ibox) &
      &                     + dtau_col(j,ibox,ilev)
              enddo
            enddo ! ibox
        enddo ! ilev

     else
      do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints 
                 if (frac_out(j,ibox,ilev).eq.1) then
                        tau(j,ibox)=tau(j,ibox)   &
     &                     + dtau_s(j,ilev)
                 endif
                 if (frac_out(j,ibox,ilev).eq.2) then
                        tau(j,ibox)=tau(j,ibox)   &
     &                     + dtau_c(j,ilev)
                 end if
              enddo
            enddo ! ibox
      enddo ! ilev
     endif
#else
      do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints 
                 if (frac_out(j,ibox,ilev).eq.1) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_s(j,ilev)
                 endif
                 if (frac_out(j,ibox,ilev).eq.2) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_c(j,ilev)
                 end if
              enddo
            enddo ! ibox
      enddo ! ilev
#endif
          if (ncolprint.ne.0) then

#ifdef COSP_GFDL
              do j=1,npoints ,1000
                write(logunit,'(a10)') 'j='
                write(logunit,'(8I10)') j
                write(logunit,'(i2,1X,8(f7.2,1X))')  &
     &          ilev,                                &
     &          (tau(j,ibox),ibox=1,ncolprint)
              enddo
#else
              do j=1,npoints ,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write(6,'(i2,1X,8(f7.2,1X))') 
     &          ilev,
     &          (tau(j,ibox),ibox=1,ncolprint)
              enddo
#endif
          endif 
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      if (top_height .eq. 1 .or. top_height .eq. 3) then


        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644
        wtmh20 = 18.01534
        Navo = 6.023E+23
        grav = 9.806650E+02
        pstd = 1.013250E+06
        t0 = 296.
#ifdef COSP_GFDL
        if (ncolprint .ne. 0) &
     &         write(logunit,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
#else
        if (ncolprint .ne. 0) 
     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
#endif
        do 125 ilev=1,nlev
          do j=1,npoints 
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10.
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02*(at(j,ilev)-t0))
#ifdef COSP_GFDL
               tauwv(j) = wk(j)*1.e-20*(          &
     &           (0.0224697*rh20s(j)*tmpexp(j)) + &
     &                (3.41817e-7*rfrgn(j)) )*0.98
#else
               tauwv(j) = wk(j)*1.e-20*( 
     &           (0.0224697*rh20s(j)*tmpexp(j)) + 
     &                (3.41817e-7*rfrgn(j)) )*0.98
#endif
               dem_wv(j,ilev) = 1. - exp( -1. * tauwv(j))
          enddo
               if (ncolprint .ne. 0) then
               do j=1,npoints ,1000
#ifdef COSP_GFDL
               write(logunit,'(a10)') 'j='
               write(logunit,'(8I10)') j
               write(logunit,'(i2,1X,3(f8.3,3X))') ilev,                 &
     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100.), &
     &           tauwv(j),dem_wv(j,ilev)
#else
               write(6,'(a10)') 'j='
               write(6,'(8I10)') j
               write(6,'(i2,1X,3(f8.3,3X))') ilev,
     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100.),
     &           tauwv(j),dem_wv(j,ilev)
#endif
               enddo
             endif
125     continue

        !initialize variables
        do j=1,npoints 
          fluxtop_clrsky(j) = 0.
          trans_layers_above_clrsky(j)=1.
        enddo

        do ilev=1,nlev
          do j=1,npoints 
 
            ! Black body emission at temperature of the layer

              bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
              !bb(j)= 5.67e-8*at(j,ilev)**4

              ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

#ifdef COSP_GFDL
                fluxtop_clrsky(j) = fluxtop_clrsky(j) &
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j)
#else
                fluxtop_clrsky(j) = fluxtop_clrsky(j) 
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
#endif
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

#ifdef COSP_GFDL
                trans_layers_above_clrsky(j)= &
     &            trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))
#else
                trans_layers_above_clrsky(j)=
     &            trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))
#endif
                   

          enddo   
            if (ncolprint.ne.0) then
             do j=1,npoints ,1000
#ifdef COSP_GFDL
              write(logunit,'(a10)') 'j='
              write(logunit,'(8I10)') j
              write (logunit,'(a)') 'ilev:'
              write (logunit,'(I2)') ilev

              write (logunit,'(a)') &
     &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
              write (logunit,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j), &
     &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
#else
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write (6,'(a)') 
     &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j),
     &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
#endif
             enddo   
            endif

        enddo   !loop over level
        
        do j=1,npoints 
          !add in surface emission
          bb(j)=1/( exp(1307.27/skt(j)) - 1. )
          !bb(j)=5.67e-8*skt(j)**4

#ifdef COSP_GFDL
          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) &
     &     * trans_layers_above_clrsky(j)
#else
          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) 
     &     * trans_layers_above_clrsky(j)
#endif
     
          !clear sky brightness temperature
          meantbclr(j) = 1307.27/(log(1.+(1./fluxtop_clrsky(j))))
	  
        enddo

        if (ncolprint.ne.0) then
        do j=1,npoints ,1000
#ifdef COSP_GFDL
          write(logunit,'(a10)') 'j='
          write(logunit,'(8I10)') j
          write (logunit,'(a)') 'id:'
          write (logunit,'(a)') 'surface'
 
          write (logunit,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
          write (logunit,'(5(f7.2,1X))') emsfc_lw,100.*bb(j), &
     &      100.*fluxtop_clrsky(j), &
     &       trans_layers_above_clrsky(j), meantbclr(j)
#else
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
          write (6,'(5(f7.2,1X))') emsfc_lw,100.*bb(j),
     &      100.*fluxtop_clrsky(j),
     &       trans_layers_above_clrsky(j), meantbclr(j)
#endif
        enddo
      endif
    

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



        if (ncolprint.ne.0) then

        do j=1,npoints ,1000
#ifdef COSP_GFDL
            write(logunit,'(a10)') 'j='
            write(logunit,'(8I10)') j
            write (logunit,'(a)') 'ts:'
            write (logunit,'(8f7.2)') (skt(j),ibox=1,ncolprint)

            write (logunit,'(a)') 'ta_rev:'
            write (logunit,'(8f7.2)') &
     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
#else
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
    
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)') 
     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
#endif

        enddo
        endif 
        !loop over columns 
        do ibox=1,ncol
          do j=1,npoints
            fluxtop(j,ibox)=0.
            trans_layers_above(j,ibox)=1.
          enddo
        enddo

        do ilev=1,nlev
              do j=1,npoints 
                ! Black body emission at temperature of the layer

              bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
              !bb(j)= 5.67e-8*at(j,ilev)**4
              enddo

            do ibox=1,ncol
              do j=1,npoints 

#ifdef COSP_GFDL
         if (present(dem_col)     ) then
               ! emissivity for point in this layer
                if (frac_out(j,ibox,ilev).eq.1 .or. &
     &              frac_out(j,ibox,ilev).eq.2) then
                dem(j,ibox)= 1. - &
     &          ( (1. - dem_wv(j,ilev)) * (1. -  dem_col(j,ibox,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if

         else
              ! emissivity for point in this layer
                if (frac_out(j,ibox,ilev).eq.1) then
                dem(j,ibox)= 1. - &
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_s(j,ilev)) )
                else if (frac_out(j,ibox,ilev).eq.2) then
                dem(j,ibox)= 1. - &
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_c(j,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if

         endif
#else
              ! emissivity for point in this layer
                if (frac_out(j,ibox,ilev).eq.1) then
                dem(j,ibox)= 1. - 
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_s(j,ilev)) )
                else if (frac_out(j,ibox,ilev).eq.2) then
                dem(j,ibox)= 1. - 
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_c(j,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if
#endif
                

                ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

#ifdef COSP_GFDL
                fluxtop(j,ibox) = fluxtop(j,ibox) &
     &            + dem(j,ibox) * bb(j)           &
     &            * trans_layers_above(j,ibox)
#else
                fluxtop(j,ibox) = fluxtop(j,ibox) 
     &            + dem(j,ibox) * bb(j)
     &            * trans_layers_above(j,ibox) 
#endif
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

#ifdef COSP_GFDL
                trans_layers_above(j,ibox)= &
     &            trans_layers_above(j,ibox)*(1.-dem(j,ibox))
#else
                trans_layers_above(j,ibox)=
     &            trans_layers_above(j,ibox)*(1.-dem(j,ibox))
#endif

              enddo ! j
            enddo ! ibox

            if (ncolprint.ne.0) then
              do j=1,npoints,1000
#ifdef COSP_GFDL
              write (logunit,'(a)') 'ilev:'
              write (logunit,'(I2)') ilev

              write(logunit,'(a10)') 'j='
              write(logunit,'(8I10)') j
              write (logunit,'(a)') 'emiss_layer:'
              write (logunit,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
 
              write (logunit,'(a)') '100.*bb(j):'
              write (logunit,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
 
              write (logunit,'(a)') '100.*f:'
              write (logunit,'(8f7.2)') &
     &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)

              write (logunit,'(a)') 'total_trans:'
              write (logunit,'(8f7.2)') &
     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
#else
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*bb(j):'
              write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)') 
     &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)
        
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)') 
     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
#endif
            enddo
          endif

        enddo ! ilev


          do j=1,npoints 
            !add in surface emission
            bb(j)=1/( exp(1307.27/skt(j)) - 1. )
            !bb(j)=5.67e-8*skt(j)**4
          end do

        do ibox=1,ncol
          do j=1,npoints 

            !add in surface emission

#ifdef COSP_GFDL
            fluxtop(j,ibox) = fluxtop(j,ibox) &
     &         + emsfc_lw * bb(j)             &
     &         * trans_layers_above(j,ibox)
#else
            fluxtop(j,ibox) = fluxtop(j,ibox) 
     &         + emsfc_lw * bb(j) 
     &         * trans_layers_above(j,ibox) 
#endif
            
          end do
        end do

        !calculate mean infrared brightness temperature
        do ibox=1,ncol
          do j=1,npoints 
            meantb(j) = meantb(j)+1307.27/(log(1.+(1./fluxtop(j,ibox))))
	  end do
        end do
	  do j=1, npoints
	    meantb(j) = meantb(j) / real(ncol)
	  end do        

        if (ncolprint.ne.0) then

          do j=1,npoints ,1000
#ifdef COSP_GFDL
          write(logunit,'(a10)') 'j='
          write(logunit,'(8I10)') j
          write (logunit,'(a)') 'id:'
          write (logunit,'(a)') 'surface'

          write (logunit,'(a)') 'emiss_layer:'
          write (logunit,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
 
          write (logunit,'(a)') '100.*bb(j):'
          write (logunit,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
 
          write (logunit,'(a)') '100.*f:'
          write (logunit,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
 
          write (logunit,'(a)') 'meantb(j):'
          write (logunit,'(8f7.2)') (meantb(j),ibox=1,ncolprint)
#else
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*bb(j):'
          write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
          
	  write (6,'(a)') 'meantb(j):'
          write (6,'(8f7.2)') (meantb(j),ibox=1,ncolprint)
#endif
      
          end do
      endif
    
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
      !If it turns out that the cloud brightness temperature
      !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
      !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
          do j=1,npoints  
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1. /  ( exp(1307.27/(attrop(j)-5.)) - 1. ) 
          enddo 
        do ibox=1,ncol
          do j=1,npoints  
#ifdef COSP_GFDL
            transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
     &                /(fluxtop_clrsky(j)-btcmin(j))
#else
            transmax(j) = (fluxtop(j,ibox)-btcmin(j))
     &                /(fluxtop_clrsky(j)-btcmin(j))
#endif
          !note that the initial setting of tauir(j) is needed so that
          !tauir(j) has a realistic value should the next if block be
          !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1. * log(max(min(transmax(j),0.9999999),0.001))

          enddo 

          if (top_height .eq. 1) then
            do j=1,npoints  
#ifdef COSP_GFDL
              if (transmax(j) .gt. 0.001 .and. &
     &          transmax(j) .le. 0.9999999) then
#else
              if (transmax(j) .gt. 0.001 .and. 
     &          transmax(j) .le. 0.9999999) then
#endif
                fluxtopinit(j) = fluxtop(j,ibox)
              tauir(j) = tau(j,ibox) *rec2p13
              endif
            enddo
            do icycle=1,2
              do j=1,npoints  
                if (tau(j,ibox) .gt. (tauchk            )) then 
#ifdef COSP_GFDL
                if (transmax(j) .gt. 0.001 .and. &
     &            transmax(j) .le. 0.9999999) then
                  emcld(j,ibox) = 1. - exp(-1. * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -   &
     &              ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06, &
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27 &
     &              / (log(1. + (1./fluxtop(j,ibox))))              
#else
                if (transmax(j) .gt. 0.001 .and. 
     &            transmax(j) .le. 0.9999999) then
                  emcld(j,ibox) = 1. - exp(-1. * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -   
     &              ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06,
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27
     &              / (log(1. + (1./fluxtop(j,ibox))))
#endif
                  if (tb(j,ibox) .gt. 260.) then
                  tauir(j) = tau(j,ibox) / 2.56
                  end if                   
                end if
                end if
              enddo
            enddo
                
          endif
        
          do j=1,npoints
            if (tau(j,ibox) .gt. (tauchk            )) then 
                !cloudy box 
		!NOTE: tb is the cloud-top temperature not infrared brightness temperature 
		!at this point in the code
                tb(j,ibox)= 1307.27/ (log(1. + (1./fluxtop(j,ibox))))
                if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
                         tb(j,ibox) = attrop(j) - 5. 
                   tau(j,ibox) = 2.13*taumin(j)
                end if
            else
                !clear sky brightness temperature
                tb(j,ibox) = meantbclr(j)
            end if
          enddo ! j
        enddo ! ibox

        if (ncolprint.ne.0) then

          do j=1,npoints,1000
#ifdef COSP_GFDL
          write(logunit,'(a10)') 'j='
          write(logunit,'(8I10)') j
 
          write (logunit,'(a)') 'attrop:'
          write (logunit,'(8f7.2)') (attrop(j))
 
          write (logunit,'(a)') 'btcmin:'
          write (logunit,'(8f7.2)') (btcmin(j))

          write (logunit,'(a)') 'fluxtop_clrsky*100:'
          write (logunit,'(8f7.2)') &
     &      (100.*fluxtop_clrsky(j))

          write (logunit,'(a)') '100.*f_adj:'
          write (logunit,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'transmax:'
          write (logunit,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'tau:'
          write (logunit,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'emcld:'
          write (logunit,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'total_trans:'
          write (logunit,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'total_emiss:'
          write (logunit,'(8f7.2)') &
     &        (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'total_trans:'
          write (logunit,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)

          write (logunit,'(a)') 'ppout:'
          write (logunit,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
#else
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j

          write (6,'(a)') 'attrop:'
          write (6,'(8f7.2)') (attrop(j))
    
          write (6,'(a)') 'btcmin:'
          write (6,'(8f7.2)') (btcmin(j))
    
          write (6,'(a)') 'fluxtop_clrsky*100:'
          write (6,'(8f7.2)') 
     &      (100.*fluxtop_clrsky(j))

          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'transmax:'
          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') 
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)') 
     &        (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') 
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
#endif
          enddo ! j
      endif

      end if

!     ---------------------------------------------------!

!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      do 30 ibox=1,ncol
        !segregate according to optical thickness
        if (top_height .eq. 1 .or. top_height .eq. 3) then  
          !find level whose temperature
          !most closely matches brightness temperature
          do j=1,npoints 
            nmatch(j)=0
          enddo
          do 29 k1=1,nlev-1
	    if (top_height_direction .eq. 2) then
	      ilev = nlev - k1 
	    else
	      ilev = k1
	    end if
            !cdir nodep
            do j=1,npoints 
	     if (ilev .ge. itrop(j)) then
#ifdef COSP_GFDL
              if ((at(j,ilev)   .ge. tb(j,ibox) .and. &
     &          at(j,ilev+1) .le. tb(j,ibox)) .or. &
     &          (at(j,ilev) .le. tb(j,ibox) .and. &
     &          at(j,ilev+1) .ge. tb(j,ibox))) then
#else
              if ((at(j,ilev)   .ge. tb(j,ibox) .and. 
     &          at(j,ilev+1) .le. tb(j,ibox)) .or.
     &          (at(j,ilev) .le. tb(j,ibox) .and. 
     &          at(j,ilev+1) .ge. tb(j,ibox))) then 
#endif
                nmatch(j)=nmatch(j)+1
		match(j,nmatch(j))=ilev
              end if  
	     end if                         
            enddo
29        continue

          do j=1,npoints 
            if (nmatch(j) .ge. 1) then
              k1 = match(j,nmatch(j))
	      k2 = k1 + 1
              logp1 = log(pfull(j,k1))
              logp2 = log(pfull(j,k2))
	      atd = max(tauchk,abs(at(j,k2) - at(j,k1)))
              logp=logp1+(logp2-logp1)*abs(tb(j,ibox)-at(j,k1))/atd
              ptop(j,ibox) = exp(logp)
#ifdef COSP_GFDL
              if(abs(pfull(j,k1)-ptop(j,ibox)) .lt. &
     &            abs(pfull(j,k2)-ptop(j,ibox))) then
#else
	      if(abs(pfull(j,k1)-ptop(j,ibox)) .lt.
     &            abs(pfull(j,k2)-ptop(j,ibox))) then
#endif
                 levmatch(j,ibox)=k1
              else
                 levmatch(j,ibox)=k2
              end if   
            else
              if (tb(j,ibox) .le. attrop(j)) then
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              end if
              if (tb(j,ibox) .ge. atmax(j)) then
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              end if                                
            end if
          enddo ! j

        else ! if (top_height .eq. 1 .or. top_height .eq. 3) 
 
          do j=1,npoints     
            ptop(j,ibox)=0.
          enddo
          do ilev=1,nlev
            do j=1,npoints     
#ifdef COSP_GFDL
              if ((ptop(j,ibox) .eq. 0. ) &
     &           .and.(frac_out(j,ibox,ilev) .ne. 0)) then
#else
              if ((ptop(j,ibox) .eq. 0. )
     &           .and.(frac_out(j,ibox,ilev) .ne. 0)) then
#endif
                ptop(j,ibox)=phalf(j,ilev)
              levmatch(j,ibox)=ilev
              end if
            end do
          end do
        end if                            
          
        do j=1,npoints
          if (tau(j,ibox) .le. (tauchk            )) then
            ptop(j,ibox)=0.
            levmatch(j,ibox)=0      
          endif 
        enddo

30    continue
              
!
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!

      !compute isccp frequencies

      !reset frequencies
      do 38 ilev=1,7
      do 38 ilev2=1,7
        do j=1,npoints ! 
             if (sunlit(j).eq.1 .or. top_height .eq. 3) then 
                fq_isccp(j,ilev,ilev2)= 0.
	     else
	        fq_isccp(j,ilev,ilev2)= output_missing_value
	     end if
        enddo
38    continue

      !reset variables need for averaging cloud properties
      do j=1,npoints 
        if (sunlit(j).eq.1 .or. top_height .eq. 3) then 
             totalcldarea(j) = 0.
             meanalbedocld(j) = 0.
             meanptop(j) = 0.
             meantaucld(j) = 0.
	else
             totalcldarea(j) = output_missing_value
             meanalbedocld(j) = output_missing_value
             meanptop(j) = output_missing_value
             meantaucld(j) = output_missing_value
	end if
      enddo ! j

      boxarea = 1./real(ncol)
     
      do 39 ibox=1,ncol
        do j=1,npoints 

#ifdef COSP_GFDL
          if (tau(j,ibox) .gt. (tauchk            ) &
     &      .and. ptop(j,ibox) .gt. 0.) then
#else
          if (tau(j,ibox) .gt. (tauchk            )
     &      .and. ptop(j,ibox) .gt. 0.) then
#endif
              box_cloudy(j,ibox)=.true.
          endif

          if (box_cloudy(j,ibox)) then

              if (sunlit(j).eq.1 .or. top_height .eq. 3) then

                boxtau(j,ibox) = tau(j,ibox)

		if (tau(j,ibox) .ge. isccp_taumin) then
		   totalcldarea(j) = totalcldarea(j) + boxarea
		
                   !convert optical thickness to albedo
#ifdef COSP_GFDL
                    albedocld(j,ibox) &
     &             = (tau(j,ibox)**0.895)/((tau(j,ibox)**0.895)+6.82)

                   !contribute to averaging
                   meanalbedocld(j) = meanalbedocld(j)  &
     &                                +albedocld(j,ibox)*boxarea
#else
                   albedocld(j,ibox)
     &		   = (tau(j,ibox)**0.895)/((tau(j,ibox)**0.895)+6.82)
         
                   !contribute to averaging
                   meanalbedocld(j) = meanalbedocld(j) 
     &                                +albedocld(j,ibox)*boxarea
#endif

                end if

            endif

          endif

          if (sunlit(j).eq.1 .or. top_height .eq. 3) then 

           if (box_cloudy(j,ibox)) then
          
              !convert ptop to millibars
              ptop(j,ibox)=ptop(j,ibox) / 100.
            
              !save for output cloud top pressure and optical thickness
              boxptop(j,ibox) = ptop(j,ibox)
    
              if (tau(j,ibox) .ge. isccp_taumin) then
	      	meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea
              end if		

              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              !determine optical depth category
#ifdef COSP_GFDL
              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau(j)=1
              else if (tau(j,ibox) .ge. isccp_taumin &
     &          .and. tau(j,ibox) .lt. 1.3) then
                itau(j)=2
              else if (tau(j,ibox) .ge. 1.3 &
     &          .and. tau(j,ibox) .lt. 3.6) then
                itau(j)=3
              else if (tau(j,ibox) .ge. 3.6 &
     &          .and. tau(j,ibox) .lt. 9.4) then
                  itau(j)=4
              else if (tau(j,ibox) .ge. 9.4 &
     &          .and. tau(j,ibox) .lt. 23.) then
                  itau(j)=5
              else if (tau(j,ibox) .ge. 23. &
     &          .and. tau(j,ibox) .lt. 60.) then
                  itau(j)=6
              else if (tau(j,ibox) .ge. 60.) then
                  itau(j)=7
              end if

              !determine cloud top pressure category
              if (    ptop(j,ibox) .gt. 0.   &
     &          .and.ptop(j,ibox) .lt. 180.) then
                 ipres(j)=1
              else if(ptop(j,ibox) .ge. 180. &
     &          .and.ptop(j,ibox) .lt. 310.) then
                  ipres(j)=2
              else if(ptop(j,ibox) .ge. 310. &
     &          .and.ptop(j,ibox) .lt. 440.) then
                  ipres(j)=3
              else if(ptop(j,ibox) .ge. 440. &
     &          .and.ptop(j,ibox) .lt. 560.) then
                  ipres(j)=4
              else if(ptop(j,ibox) .ge. 560. &
     &          .and.ptop(j,ibox) .lt. 680.) then
                  ipres(j)=5
              else if(ptop(j,ibox) .ge. 680. &
     &          .and.ptop(j,ibox) .lt. 800.) then
                  ipres(j)=6
              else if(ptop(j,ibox) .ge. 800.) then
                  ipres(j)=7
              end if
#else
              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau(j)=1
              else if (tau(j,ibox) .ge. isccp_taumin
     &                                    
     &          .and. tau(j,ibox) .lt. 1.3) then
                itau(j)=2
              else if (tau(j,ibox) .ge. 1.3 
     &          .and. tau(j,ibox) .lt. 3.6) then
                itau(j)=3
              else if (tau(j,ibox) .ge. 3.6 
     &          .and. tau(j,ibox) .lt. 9.4) then
                  itau(j)=4
              else if (tau(j,ibox) .ge. 9.4 
     &          .and. tau(j,ibox) .lt. 23.) then
                  itau(j)=5
              else if (tau(j,ibox) .ge. 23. 
     &          .and. tau(j,ibox) .lt. 60.) then
                  itau(j)=6
              else if (tau(j,ibox) .ge. 60.) then
                  itau(j)=7
              end if

              !determine cloud top pressure category
              if (    ptop(j,ibox) .gt. 0.  
     &          .and.ptop(j,ibox) .lt. 180.) then
                  ipres(j)=1
              else if(ptop(j,ibox) .ge. 180.
     &          .and.ptop(j,ibox) .lt. 310.) then
                  ipres(j)=2
              else if(ptop(j,ibox) .ge. 310.
     &          .and.ptop(j,ibox) .lt. 440.) then
                  ipres(j)=3
              else if(ptop(j,ibox) .ge. 440.
     &          .and.ptop(j,ibox) .lt. 560.) then
                  ipres(j)=4
              else if(ptop(j,ibox) .ge. 560.
     &          .and.ptop(j,ibox) .lt. 680.) then
                  ipres(j)=5
              else if(ptop(j,ibox) .ge. 680.
     &          .and.ptop(j,ibox) .lt. 800.) then
                  ipres(j)=6
              else if(ptop(j,ibox) .ge. 800.) then
                  ipres(j)=7
              end if 
#endif

              !update frequencies
              if(ipres(j) .gt. 0.and.itau(j) .gt. 0) then
#ifdef COSP_GFDL
              fq_isccp(j,itau(j),ipres(j))= &
     &          fq_isccp(j,itau(j),ipres(j))+ boxarea
#else
              fq_isccp(j,itau(j),ipres(j))=
     &          fq_isccp(j,itau(j),ipres(j))+ boxarea
#endif
              end if

            end if

          end if
                       
        enddo ! j
39    continue
      
      !compute mean cloud properties
      do j=1,npoints 
        if (totalcldarea(j) .gt. 0.) then
	  ! code above guarantees that totalcldarea > 0 
	  ! only if sunlit .eq. 1 .or. top_height = 3 
	  ! and applies only to clouds with tau > isccp_taumin
          meanptop(j) = meanptop(j) / totalcldarea(j)
          meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
          meantaucld(j) = (6.82/((1./meanalbedocld(j))-1.))**(1./0.895)
	else
	  ! this code is necessary so that in the case that totalcldarea = 0.,
	  ! that these variables, which are in-cloud averages, are set to missing
	  ! note that totalcldarea will be 0. if all the clouds in the grid box have
	  ! tau < isccp_taumin 
	  meanptop(j) = output_missing_value
          meanalbedocld(j) = output_missing_value
          meantaucld(j) = output_missing_value
        end if
      enddo ! j
!
!     ---------------------------------------------------!

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      if (debugcol.ne.0) then
!     
         do j=1,npoints,debugcol

            !produce character output
            do ilev=1,nlev
              do ibox=1,ncol
                   acc(ilev,ibox)=0
              enddo
            enddo

            do ilev=1,nlev
              do ibox=1,ncol
                   acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
#ifdef COSP_GFDL
                   if (levmatch(j,ibox) .eq. ilev) &
     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
#else
                   if (levmatch(j,ibox) .eq. ilev) 
     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
#endif
              enddo
            enddo

             !print test

          write(ftn09,11) j
11        format('ftn09.',i4.4)
#ifdef COSP_GFDL 
          funit = get_unit()
          open(funit, FILE=ftn09, FORM='FORMATTED')

             write(funit,'(a1)') ' '
             write(funit,'(10i5)')  &
     &                  (ilev,ilev=5,nlev,5)
             write(funit,'(a1)') ' '

             do ibox=1,ncol
               write(funit,'(40(a1),1x,40(a1))')                &
     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) &
     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev)
             end do
             close(funit)
#else
          open(9, FILE=ftn09, FORM='FORMATTED')

             write(9,'(a1)') ' '
             write(9,'(10i5)')
     &                  (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '
             
             do ibox=1,ncol
               write(9,'(40(a1),1x,40(a1))')
     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) 
     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
             end do
             close(9)
#endif

             if (ncolprint.ne.0) then
#ifdef COSP_GFDL
               write(logunit,'(a1)') ' '
                    write(logunit,'(a2,1X,5(a7,1X),a50)') &
     &                  'ilev',                           &
     &                  'pfull','at',                     &
     &                  'cc*100','dem_s','dtau_s',        &
     &                  'cchar'

!               do 4012 ilev=1,nlev
!                    write(logunit,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(logunit,'(i2,1X,5(f7.2,1X),50(a1))')
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
               write (logunit,'(a)') 'skt(j):'
               write (logunit,'(8f7.2)') skt(j)

               write (logunit,'(8I7)') (ibox,ibox=1,ncolprint)
               
               write (logunit,'(a)') 'tau:'
               write (logunit,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)

               write (logunit,'(a)') 'tb:'
               write (logunit,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)

               write (logunit,'(a)') 'ptop:'
               write (logunit,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
#else
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)') 
     &                  'ilev',
     &                  'pfull','at',
     &                  'cc*100','dem_s','dtau_s',
     &                  'cchar'

!               do 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(6,'(i2,1X,5(f7.2,1X),50(a1))') 
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
               write (6,'(a)') 'skt(j):'
               write (6,'(8f7.2)') skt(j)
                                      
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
            
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
#endif
             endif 
    
        enddo
       
      end if 

      return
      end 



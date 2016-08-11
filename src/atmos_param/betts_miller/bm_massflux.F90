
module bm_massflux_mod

!----------------------------------------------------------------------
use            mpp_mod, only:  input_nml_file
use            fms_mod, only:  file_exist, error_mesg, open_namelist_file,  &
                               check_nml_error, mpp_pe, FATAL,  &
                               close_file, mpp_root_pe, write_version_number, stdlog
use sat_vapor_pres_mod, only:  escomp, descomp
use      constants_mod, only:  HLv,HLs,Cp_air,Grav,rdgas,rvgas, cp_vapor, kappa

implicit none
private
!----------------------------------------------------------------------
!  ---- public interfaces ----

   public  bm_massflux, bm_massflux_init, bm_massflux_end

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id: bm_massflux.F90,v 19.0 2012/01/06 20:01:33 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'

!-----------------------------------------------------------------------
!   ---- local/private data ----

    real, parameter :: d622 = rdgas/rvgas
    real, parameter :: d378 = 1.-d622

    logical :: module_is_initialized=.false.

!-----------------------------------------------------------------------
!   --- namelist ----

real    :: tau_bm=7200.
real    :: rhbm = .8
real    :: revap = 1.0
real    :: minstatstab = .1
logical :: do_shallower = .false.
logical :: do_changeqref = .false.

namelist /bm_massflux_nml/  tau_bm, rhbm, do_shallower, do_changeqref, revap, &
   minstatstab

!-----------------------------------------------------------------------
!           description of namelist variables
!
!  tau_bm    =  betts-miller relaxation timescale (seconds)
!  rhbm      = relative humidity that you're relaxing towards
!             In contrast to standard BM, this can be put as high as 100%
!             as long as the convection scheme is not overwhelmd by 
!             the large-scale condensation.
!  revap     = the precipiation efficiency reevaporation parameter, i.e., 
!              the moistening of the atmosphere per heating.  
!
!  not used:
!
!
!  do_shallower = do the shallow convection scheme where it chooses a smaller
!                 depth such that precipitation is zero
! 
!  do_changeqref = do the shallow convection scheme where if changes the 
!                  profile of both q and T in order make precip zero
!
!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine bm_massflux (dt, tin, qin, pfull, phalf, coldT, &
                           rain, snow, tdel, qdel, q_ref, bmflag, &
                           klzbs, t_ref, massflux,&
                           mask, conv)

!-----------------------------------------------------------------------
!
!                     Betts-Miller Convection Scheme
!
!-----------------------------------------------------------------------
!
!   input:  dt       time step in seconds
!           tin      temperature at full model levels
!           qin      specific humidity of water vapor at full
!                      model levels
!           pfull    pressure at full model levels
!           phalf    pressure at half (interface) model levels
!           coldT    should precipitation be snow at this point?
!   optional:
!           mask     optional mask (0 or 1.) 
!           conv     logical flag; if true then no betts-miller
!                       adjustment is performed at that grid-point or
!                       model level
!
!  output:  rain     liquid precipitation (kg/m2)
!           snow     frozen precipitation (kg/m2)
!           tdel     temperature tendency at full model levels
!           qdel     specific humidity tendency (of water vapor) at
!                      full model levels
!           bmflag   flag for which routines you're calling
!           klzbs    stored klzb values
!           massflux the massflux used to calculate the humidity adjustment
!
!-----------------------------------------------------------------------
!--------------------- interface arguments -----------------------------

   real   , intent(in) , dimension(:,:,:) :: tin, qin, pfull, phalf
   real   , intent(in)                    :: dt
   logical, intent(in) , dimension(:,:)   :: coldT
   real   , intent(out), dimension(:,:)   :: rain,snow, bmflag, klzbs
   real   , intent(out), dimension(:,:,:) :: tdel, qdel, q_ref, t_ref, massflux
   real   , intent(in) , dimension(:,:,:), optional :: mask
   logical, intent(in) , dimension(:,:,:), optional :: conv
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------

   logical :: avgbl
   real,dimension(size(tin,1),size(tin,2),size(tin,3)) ::  rin
   real,dimension(size(tin,1),size(tin,2))             ::  &
                     hlcp, precip, cape, cin
   real,dimension(size(tin,3))                         :: rpc, tpc, &
             statstab, theta
   real                                                ::  & 
       cape1, cin1, tot, deltak, deltaq, qrefint, deltaqfrac, deltaqfrac2, &
       ptopfrac
integer  i, j, k, ix, jx, kx, klzb, ktop, klcl

!modif omp:
real :: prec_rev, ratio, q_src, ratio_q, ratio_T, en_acc, ratio_ml
!-----------------------------------------------------------------------
!     computation of precipitation by betts-miller scheme
!-----------------------------------------------------------------------

      if (.not. module_is_initialized) call error_mesg ('bm_massflux',  &
                         'bm_massflux_init has not been called.', FATAL)

      ix=size(tin,1)
      jx=size(tin,2)
      kx=size(tin,3)
      avgbl = .false.

!----- compute proper latent heat --------------------------------------
!
!  omp: at this stage,uses only the latent heat of vaporization. 
!       --- this is used only by the energy budget. The lapse rate 
!       has its own value. (This should be made consistent.)  
!      
!
!      WHERE (coldT)
!           hlcp = HLs/Cp_air
!      ELSEWHERE
           hlcp = HLv/Cp_air
!      END WHERE

!------------------- calculate CAPE and CIN ----------------------------

! calculate r
       rin = qin/(1.0 - qin)
       do i=1,ix
          do j=1,jx

             cape1 = 0.
             cin1 = 0.
             tot = 0.
             klzb=0
             klcl = 0
! the bmflag is written out to show what aspects of the bm scheme is called
! bmflag = 0 is no cape, no convection
! bmflag = 1 cape but no precip: nothing happens
! bmflag = 2 is deep convection
! it's either 0 or 2 in this scheme


             bmflag(i,j) = 0.
             tpc = tin(i,j,:)
             rpc = rin(i,j,:)
! calculate cape, cin, level of zero buoyancy, and parcel properties w/ leo's
! code
!             call capecalc( Cp_air, cp_vapor, d622, kx, pfull(i,j,:),&
!                            rin(i,j,:), rdgas, hlv, rvgas, tin(i,j,:),&
!                            cape1, cin1, tot, tpc, rpc, klcl, klzb)
! new code (second order in delta ln p and exact LCL calculation)
             call capecalcnew( kx,  pfull(i,j,:),  phalf(i,j,:),&
                            Cp_air, rdgas, rvgas, hlv, kappa, tin(i,j,:), &
                            rin(i,j,:), avgbl, cape1, cin1, tpc, &
                            rpc, klzb,klcl)
! set values for storage
             cape(i,j) = cape1
             cin(i,j) = cin1
             klzbs(i,j) = klzb
             massflux(i,j,1:kx) = 0.
             if(cape1.gt.0.) then
!             if( (tot .gt. 0.) .and. (cape1.gt.0.) ) then 
                bmflag(i,j) = 1.

! reference temperature is just that of the parcel all the way up

                t_ref(i,j,:) = tpc

                do k=klzb,kx
!                   q_ref(i,j,k) = rpc(k)/(1. + rpc(k))

! use this for relaxation toward virtual adiabat
!                t_ref(i,j,k) = tpc(k) * ( 1+0.608 * q_ref(i,j,k))/ (1+ 0.608 * qin(i,j,k))


! modif omp: free troposphere is relaxed  toward saturated profile at CURRENT temperature
!                   call escomp(tin(i,j,k),es)
!                   rs=d622*es/(pfull(i,j,k)+(d622-1.)*es)
!                   q_sat(k) = rs/(1. + rs)


! This shouldn't happen, but set reference humidity to be positive just in 
! case
!                   if (q_ref(i,j,k).lt.0.) then
!                      q_ref(i,j,k) = 0.
!                      write (*,*) 'doh! q neg'
!                   end if

                end do

! set the tendencies to zero where you don't adjust
! set the reference profiles to be the original profiles (for diagnostic 
! purposes only --  you can think of this as what you're relaxing to in
! areas above the actual convection

                

                do k=1,max(klzb-1,1)
                   qdel(i,j,k) = 0.0
                   tdel(i,j,k) = 0.0
                   q_ref(i,j,k) = qin(i,j,k)
                   t_ref(i,j,k) = tin(i,j,k)
                end do




! initialize p to zero for the loop

!modif omp: my way...

! check on the LCL. if LCL is less than 50mb above the surface, 
! the 'lcl' is put at the first 
! level above 950mb (required for adjustment near surface) 
! (adjustment dries out ML if lcl too low).

                if (phalf(i,j,kx+1) -phalf(i,j,klcl+1) .LT. 2500) then
                   do k = kx-1,1,-1
                      if (phalf(i,j,kx+1) -phalf(i,j,k+1) .GT. 2500) then
                         klcl = k
                         go to 11
                      end if
                   end do
                end if

! alternatively, just ensure that there is one layer where convection occurs.
!if (klcl .eq. kx ) then
!   klcl = kx -1
!end if


11  continue


                precip(i,j) = 0.0
                prec_rev = 0.0
                en_acc = 0.0
                q_src = 0.0

!adjustment above PBL. 1st compute difference between the 
! profile and computed tendecies

! We need the potential temperature for the static stability.  
                theta(1:kx) = tin(i,j,1:kx)* &
                   (phalf(i,j,kx+1)/pfull(i,j,1:kx))**kappa
! This is actually the static stability * dz (all you need).  
                statstab(1:kx-1) = theta(1:kx-1) - theta(2:kx)
! Since we're going to divide by the static stability, make sure it's not
! too close to zero (minstatstab is a nml parameter).
                do k=1,kx-1
                statstab(k) = max(statstab(k),minstatstab)
!               if (statstab(k).eq.minstatstab) then
!                write (*,*) 'ouch', k
!               end if 
                enddo
! This should never happen, but just in case, make sure the level of zero 
! buoyancy isn't at the very top.  
                if (klzb.eq.1) then
                   klzb = 2
!                   write (*,*) 'doh, klzb = 1'
                end if 
! T is adjusted toward pseudo-adiabat -- this and the static stab determines 
! a mass flux.  
! q is adjusted based on mass flux and moisture stability
                do k=klzb, klcl
! tdel is the temperature adjustment when you multiply it by dt/tau
                   tdel(i,j,k) = - (tin(i,j,k) - t_ref(i,j,k))
! This is actually the downward mass flux/dz*tau.
                   massflux(i,j,k) = tdel(i,j,k)/statstab(k-1)
! Cap the massflux at a maximum value.  1.0 is picked to be 1/2 the 
! maximum allowed by the Courant criterion for tau=7200 s, dt=1800 s
                   if (abs(massflux(i,j,k)).gt.1.0) then
!                   write (*,*) massflux(i,j,k)
                   massflux(i,j,k) = 1.0*massflux(i,j,k)/abs(massflux(i,j,k))
!                   write (*,*) massflux(i,j,k), k, tdel(i,j,k), statstab(k-1)
                   endif
! qdel is the humidity adjustment when you multiply it by dt/tau
! This is just upwind differencing.
! Also, there is a Courant criterion for this advection equation.  The 
! advection velocity is massflux(i,j,k)*dz/tau.  The condition is:
! dt < dz/M = tau/massflux(i,j,k), i.e., massflux(i,j,k) < tau/dt
! The massflux can achieve values around 10 (minstatstab = .5, and 
! the temperature deficits of 5 deg).  
                   if (massflux(i,j,k).gt.0.) then 
                      qdel(i,j,k) = massflux(i,j,k)*(qin(i,j,k-1) - qin(i,j,k))
                   else
                      qdel(i,j,k) = massflux(i,j,k)*(qin(i,j,k) - qin(i,j,k+1))
                   end if 
! Additional moistening from precip eff if heating > 0
                   if (tdel(i,j,k).gt.0.0) then
                      qdel(i,j,k) = qdel(i,j,k) + revap*Cp_air/hlv*tdel(i,j,k)
                   endif 
!!                   qdel(i,j,k) = - (qin(i,j,k) - rhbm * q_ref(i,j,k))
!                   qdel(i,j,k) = - (qin(i,j,k) - rhbm * q_sat(k))
                   precip(i,j) = precip(i,j) + tdel(i,j,k)*(phalf(i,j,k+1) - &
                        phalf(i,j,k))/grav/hlcp(i,j)
                   prec_rev = prec_rev + qdel(i,j,k)*(phalf(i,j,k+1) - &
                        phalf(i,j,k))/grav
                end do


                do k = klcl + 1, kx
!no temperature change in the ML
                   tdel(i,j,k) = 0.0
! relax the ML temp. toward the pseudo-adiabat
!                tdel(i,j,k) = - (tin(i,j,k) - t_ref(i,j,k))
                   precip(i,j) = precip(i,j) + tdel(i,j,k)*(phalf(i,j,k+1) - &
                        phalf(i,j,k))/grav/hlcp(i,j)
                end do

! at this point: precip = total amount of precip required to brin the temperature to the m. adiabat.
!                prec_rev = total water deficit (respectively to saturation at CURRENT temperature)

! humidity adjustment in the ML

! modif df: the following calculations are no longer used, due to the modification
!           described below
!                en_acc = precip(i,j) + min(prec_rev, revap*precip(i,j))
!                en_acc = en_acc * grav /   &
!                   (phalf(i,j,klcl+1) - phalf(i,j,klzb))


! en_acc is the average energy difference (measured in qv) between the atmosphere and the state at which it would
! equilibrate for the current adiabat. 

! adjustment in the pbl
                do k = klcl + 1, kx

! modif df: changed this line (it was caused energy not to be conserved in the
!           case where the atmosphere is too moist, and en_acc is negative --
!           this allowed q_src to be negative, and the expression in that 
!           case was incorrect.  this could have been fixed by putting in a 
!           different expression, but it seems like the boundary layer should
!           be able to be a source of moisture even if the free troposphere 
!           is being dried; further the boundary layer should never be 
!           moistened (which was happening before).  therefore, this line is 
!           removed and replaced just by uniform drying of the boundary layer)
!                qdel(i,j,k) = max( - qin(i,j,k) ,-en_acc)
                   qdel(i,j,k) = -qin(i,j,k)
                   q_src = q_src - qdel(i,j,k)*(phalf(i,j,k+1) - &
                        phalf(i,j,k))/grav
                end do

! Determine the adjustment time-scale


                if (precip(i,j) .gt. 0.0) then 

! If precip > 0, then correct energy. 

                   bmflag(i,j) = 2.

                   ratio = min(dt / tau_bm, 1.0)
! modif df: due to the modification above, it's impossible for this to be negative.
!           therefore everything is commented out.  in case someone wants to use
!           this however, i put in the proper correction into the ratio_q line
!           (the added minus sign, which keeps the timescale positive in this case),
!           so this will conserve energy if you uncomment all the rele,vant lines.  
!                if (q_src .lt. 0.0) then
!                   ratio_ML = 0.0
!                   if (prec_rev .gt. 0.0) then
!should not be happening
!                      write (*,*) 'BUGGGGGGG', precip(i,j), prec_rev, en_acc, q_src
!                      write (*,*) 'QREF:',q_ref(i,j,klcl+1:kx),'QIN',qin(i,j,klcl+1:kx)
!                   end if
!                   ratio_T = ratio
!                   ratio_q = -precip(i,j) * ratio_T / prec_rev
!                else

!try convection limited by low level moisture supply
! modif df: a lot of this doesn't make sense to me, so this is commented out and 
!           changed to what's described in olivier's writeup
!                   if (prec_rev .lt. 0.0 ) then
!                      ratio_q = ratio
!                      ratio_ML = ratio
!                      ratio_T = (ratio_ML * q_src - ratio_q * prec_rev)  &
!                           / precip(i,j)    
!                      if (ratio_T .GT. ratio) then
!                         ratio_T = ratio
!                         ratio_ML = ( ratio_T * precip(i,j) + ratio_q * prec_rev) &
!                           /q_src
!                         if (ratio_ML .GT. ratio) then
!                            write(*,*)'danger: ', ratio, ratio_ML, precip(i,j), q_src, prec_rev
!                         end if
!                      end if
!                   else 
!                      ratio_ML = ratio
!                      ratio_T = ratio * q_src / ((1.0 + revap) * precip(i,j))
!!               ratio_q = revap * ratio_T
!                      ratio_q = precip(i,j) * revap * ratio_T / prec_rev
!                      if (ratio_q .gt. ratio) then
!                         ratio_q = ratio
!                         ratio_T = (ratio_ML * q_src - ratio_q * prec_rev)  &
!                             / precip(i,j)  
!                      end if
!
!
!                      if ( ratio_T  .gt. ratio) then
!                  
!                         ratio_T = ratio
!                         ratio_q = max(ratio, precip(i,j)*revap*ratio_T/prec_rev)
!                         ratio_ML = (ratio_T*precip(i,j) + ratio_q*prec_rev) &
!                           /q_src
!                      end if
!                   end if
! Adjust T and q, and conserve energy by draining the boundary layer by the 
! required amount.  
! Don't allow the ML to be moistened!  
                   if (precip(i,j)+prec_rev.gt.0.0) then
                      ratio_q = ratio
                      ratio_T = ratio
                      ratio_ML = ( ratio_T * precip(i,j) + ratio_q * prec_rev) &
                        /q_src
! If you're draining the boundary layer too fast (in a way that would allow 
! humidities to become negative), then limit the convection by reducing 
! the adjustment of temperature and humidity in concert.  
                      if (ratio_ML.gt.ratio) then 
                         ratio_ML = ratio
                         ratio_T = (ratio_ML * q_src )  &
                              / (precip(i,j) + prec_rev)
                         ratio_q = ratio_T
                      end if
                   else
! If the scheme tries to moisten the ML, don't allow this.  
                      ratio_ML = 0.0
                      ratio_T = ratio
                      ratio_q = -ratio_T*precip(i,j)/prec_rev
                   end if

!adjust the tendencies
                
                   do k = klzb,klcl
                      qdel(i,j,k)= qdel(i,j,k) * ratio_q  
                      tdel(i,j,k)= tdel(i,j,k) * ratio_T 
                      massflux(i,j,k) = massflux(i,j,k)*ratio_q/dt* &
                          (phalf(i,j,k+1) - phalf(i,j,k))/grav
                   end do
                   do k = klcl+1,kx
                      qdel(i,j,k)= qdel(i,j,k) * ratio_ML  
                      tdel(i,j,k)= tdel(i,j,k) * ratio_T
                   end do
                   precip(i,j) = precip(i,j) * ratio_T
                       





                else

!omp: This is Dargan's and not used in my version (put do_shallower and  
!do change_qref as .false.

! Shallow / non-precipitating adjustment from dargan. It has not been tested 
! with this version of the deep convection scehme. Use do_shallower = .false.
! and do_changeqref = .false. to turn off. In this case, no convection occurs 
! when the computed precipitation is 0. 


! If precip < 0, then do the shallow conv routine.
! First option: do_shallower = true
! This chooses the depth of convection based on choosing the height that 
! it can make precip zero, i.e., subtract off heights until that precip 
! becomes positive.  

                    if (do_shallower) then
! ktop is the new top of convection.  set this initially to klzb.
                       ktop = klzb
! Work your way down until precip is positive again.
                       do while ( (precip(i,j).lt.0) .and. (ktop.le.kx) )
                          precip(i,j) = precip(i,j) - qdel(i,j,ktop)* &
                                   (phalf(i,j,ktop) - phalf(i,j,ktop+1))/grav
                          ktop = ktop + 1
                       end do
! since there will be an overshoot (precip is going to be greater than zero 
! once we finish this), the actual new top of convection is somewhere between
! the current ktop, and one level above this.  set ktop to the level above.
                       ktop = ktop - 1
! Adjust the tendencies in the places above back to zero, and the reference 
! profiles back to the original t,q.
                       if (ktop.gt.klzb) then
                          qdel(i,j,klzb:ktop-1) = 0.
                          q_ref(i,j,klzb:ktop-1) = qin(i,j,klzb:ktop-1)
                          tdel(i,j,klzb:ktop-1) = 0.
                          t_ref(i,j,klzb:ktop-1) = tin(i,j,klzb:ktop-1)
                       end if
! Then make the change only a fraction of the new top layer so the precip is 
! identically zero.
! Calculate the fractional penetration of convection through that top layer.  
! This is the amount necessary to make precip identically zero.  
                       ptopfrac = precip(i,j)/(qdel(i,j,ktop)* &
                          (phalf(i,j,ktop+1) - phalf(i,j,ktop)))
! Reduce qdel in the top layer by this fraction. 
                       qdel(i,j,ktop) = ptopfrac*qdel(i,j,ktop)
! Set precip to zero
                       precip(i,j) = 0.
! A diagnostic which allows calculating precip to make sure it's zero.
                       do k=ktop,kx
                          precip(i,j)=precip(i,j)+qdel(i,j,k)* &
                                 (phalf(i,j,k) - phalf(i,j,k+1))/grav
                       end do
                       if (abs(precip(i,j)).gt.1.e-5) &
                           write(6,*) 'doh! precip.ne.0'
! Now change the reference temperature in such a way to make the net 
! heating zero.
                       deltak = 0.
                       if (ktop.lt.kx) then
! Integrate temperature tendency up to 1 level below top.
                          do k=ktop+1,kx
                             deltak = deltak + tdel(i,j,k)* &
                                 (phalf(i,j,k) - phalf(i,j,k+1))
                          end do
! Then for the top level, use only a fraction.
                          deltak = deltak + ptopfrac*tdel(i,j,ktop)* &
                               (phalf(i,j,ktop) - phalf(i,j,ktop+1))
! Normalize by the pressure difference.
                          deltak = deltak/(phalf(i,j,kx+1) - & 
                           phalf(i,j,ktop+1) + ptopfrac*(phalf(i,j,ktop+1) - &
                           phalf(i,j,ktop)))
! Subtract this value uniformly from tdel, and make the according change to 
! t_ref.
                          do k=ktop,kx
                             tdel(i,j,k) = tdel(i,j,k) + deltak
                             t_ref(i,j,k) = t_ref(i,j,k) + deltak*tau_bm/dt
                          end do
                       end if
                    else if(do_changeqref) then
! Change the reference profile of q by a certain fraction so that precip is 
! zero.  This involves calculating the total integrated q_ref dp (this is the
! quantity intqref), as well as the necessary change in q_ref (this is the 
! quantity deltaq).  Then the fractional change in q_ref at each level (the 
! quantity deltaqfrac) is 1-deltaq/intqref.  (have to multiply q_ref by 
! 1-deltaq/intqref at every level)  Then the change in qdel is 
! -deltaq/intqref*q_ref*dt/tau_bm.
! Change the reference profile of T by a uniform amount so that precip is zero.
                       deltak = 0.
                       deltaq = 0.
                       qrefint = 0.
                       do k=klzb,kx
! deltaq = a positive quantity (since int qdel is positive).  It's how 
! much q_ref must be changed by, in an integrated sense.  The requisite 
! change in qdel is this without the factors of tau_bm and dt.
                          deltaq = deltaq - qdel(i,j,k)*tau_bm/dt* &
                                    (phalf(i,j,k) - phalf(i,j,k+1))
! deltak = the amount tdel needs to be changed
                          deltak  = deltak  + tdel(i,j,k)* &
                                    (phalf(i,j,k) - phalf(i,j,k+1))
! qrefint = integrated value of qref
                          qrefint = qrefint - q_ref(i,j,k)* &
                                    (phalf(i,j,k) - phalf(i,j,k+1))
                       end do
! Normalize deltak by total pressure.
                       deltak  = deltak /(phalf(i,j,kx+1) - phalf(i,j,klzb))
! multiplying factor for q_ref is 1 + the ratio
                       deltaqfrac = 1. - deltaq/qrefint
! multiplying factor for qdel adds dt/tau_bm
                       deltaqfrac2 = - deltaq/qrefint*dt/tau_bm
! let's check that the precip really is zero as in the shallower scheme
                       precip(i,j) = 0.0
                       do k=klzb,kx
                          qdel(i,j,k) = qdel(i,j,k) + deltaqfrac2*q_ref(i,j,k)
                          q_ref(i,j,k) = deltaqfrac*q_ref(i,j,k)
                          tdel(i,j,k) = tdel(i,j,k) + deltak
                          t_ref(i,j,k) = t_ref(i,j,k) + deltak*tau_bm/dt
                          precip(i,j) = precip(i,j) + qdel(i,j,k)* &
                                 (phalf(i,j,k) - phalf(i,j,k+1))/grav
                       end do
                       if (abs(precip(i,j)).gt.1.e-5) &
                         write(6,*) 'doh! precip.ne.0)'
                    else
                       precip(i,j) = 0.
                       tdel(i,j,:) = 0.
                       qdel(i,j,:) = 0.
                    end if
                end if
             else 
                tdel(i,j,:) = 0.0
                qdel(i,j,:) = 0.0
                precip(i,j) = 0.0
                q_ref(i,j,:) = qin(i,j,:)
                t_ref(i,j,:) = tin(i,j,:)
             end if
          end do
       end do

       rain = precip
       snow = 0.
   

!-----------------------------------------------------------------------

   end subroutine bm_massflux

!#######################################################################

    subroutine capecalc(cpd, cpv, epsilo, nlev, pback, rback, rd, rl, rv, &
                        tback, xcape, cin, tot, tpcback, rpcback,         &
                        klclback, klzbback)
!
! modif omp: klcl is an additional output

!      Calculates convective available potential energy for a cloud whose
!      temperature follows a saturated adiabat.
!
!      On Input
!
!      cpd     specific heat of dry air at constant pressure (J/(kg K))
!      cpv     specific heat of water vapor
!      epsilo  ratio of molecular weights of water vapor to dry air
!      nlev    number of levels 
!      p       pressure (Pa)
!              Index 1 refers to level nearest earth's surface.
!      r       mixing ratio (kg(H2O)/kg)
!              Index 1 refers to level nearest earth's surface.
!      rd      gas constant for dry air (J/(kg K))
!      rl      latent heat of vaporization (J/kg)
!      rv      gas constant for water vapor (J/(kg K))
!      t       temperature (K)
!              Index 1 refers to level nearest earth's surface.
!
!     Output:
!   
!     tpc      parcel temperature (K)
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     rpc      parcel mixing ratio (kg(H2O)/kg)
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     cin      convective inhibition (J/kg)
!              energy required to lift parcel from level istart to
!              level of free convection
!     xcape    convective available potential energy (J/kg)
!              energy released as parcel moves from level of free
!              convection to level of zero buoyancy
!     tot      xcape+cin (J/kg)
!     klcl     first level above the LCL
!     klzb     the level where you hit LZB
!
!     For definitions of cin and xcape, see notes (4 Apr 95) (LAN Notes).
!
        implicit none
        integer, intent(in) :: nlev
        REAL, INTENT (IN),    DIMENSION(:) :: pback, tback, rback
        real, intent (in) :: cpd, cpv, epsilo, rd, rl, rv
        integer, intent(out) :: klzbback, klclback
        real, intent (out) :: xcape, cin, tot
        real, intent (out), dimension(nlev) :: tpcback, rpcback

      integer :: istart, ieq, klcl, k, klzb, klfc, ieqa, nlevm
      logical :: capepos
      real, dimension(nlev) :: p, r, t, tpc, rpc
      real :: ro, tc, tp, plcl, es, rs, rlcl, tlcl, pb, tb, rb, q, cp, dp, &
              dt1, plzb, qe, tve, pc, qs, tv, rc, fact1, fact2, fact3, &
              dtdp, rbc, rbe, qc, tvc, delt

!modif omp
      real :: rsb, tplcl

!      parameter(nlev=25,nlevm=nlev-1)
!      dimension t(nlev),r(nlev),p(nlev)
!      dimension tpc(nlev),rpc(nlev)
!      logical capepos
!
      capepos=.false.
      nlevm = nlev-1

     do k=1,nlev
        t(k) = tback(nlev+1-k)
        r(k) = rback(nlev+1-k)
        p(k) = pback(nlev+1-k)
        tpc(k) = t(k)
        rpc(k) = r(k)
     end do

!
!     Calculate LCL
!     istart-index of level whose mixing ratio is conserved as a parcel 
!            leaves it undergoing dry adiabatic ascent
!
      istart=1
      ro=r(istart)
      tc=t(istart)
      tp=tc
      plcl=0.
      do k=1,istart
        tpc(k)=t(k)
        rpc(k)=r(k)
      end do
      do k=istart,nlev
        call establ(es,tp)
        rs=epsilo*es/(p(k)+(epsilo-1.)*es)
!       write(6,*) 'k,tp,rs= ',k,tp,rs
        ieq=iequ(rs,ro)
        if (ieq .eq. 0) then
          plcl=p(k)
          rlcl=r(k)
          tlcl=t(k)
          klcl=k
          go to 11
        end if
        if (k .eq. istart) then
           if (ieq .lt. 0) then
              plcl=p(istart)
              tlcl=t(istart)
              rlcl=r(istart)
! omp: try this before, did not change much. 
! (The lowest level should not be supersaturated.)
!             rlcl=rs
              klcl=istart
              go to 11
           else
              go to 13
           end if
        end if
        if (k .gt. 1) then
           pb=(p(k)+p(k-1))/2.
           tb=(t(k)+t(k-1))/2.
           rb=(r(k)+r(k-1))/2.
           if (rs .lt. ro) then

              fact1 = (ro - rs)/(rsb-rs)
              fact2 = (rsb - ro)/(rsb-rs)
             plcl=fact1* p(k-1) + fact2 * p(k) 
             tlcl=fact1* t(k-1) + fact2 * t(k) 
             rlcl=rb
             tplcl = fact1* tpc(k-1) + fact2 * tpc(k) 
             klcl=k
             go to 11
           end if
        end if
        if (k .eq. nlev) go to 11

!     Convert mixing ratio to specific humidity.
!
 13   continue
         q=ro/(1.+ro)
         cp=cpd*(1.+((cpv/cpd)-1.)*q)
         dp=p(k+1)-p(k)
         dtdp=rd*tp/cp
!        write(6,*) 'dp,dtdp,pb= ',dp,dtdp,pb
         dt1=dtdp*alog((p(k)+dp)/p(k))
         tp=tp+dt1
         tpc(k+1)=tp
         rpc(k+1)=ro
         rsb = rs
      end do
 11   continue
      ieq=iequ(plcl,0.)
      if (ieq .eq. 0) then
         xcape = 0.
         cin = 0.
         tot = 0.
         tpcback = tback
         rpcback = rback
!        write(6,*) 'plcl=0'
!         stop
         call error_mesg ('bm_massflux:capecalc', 'ieq = 0', FATAL)
      end if
!
!     Calculate temperature along saturated adiabat, starting at p(kLCL).
!
!       write(6,*) 'plcl,klcl,tlcl,rlcl= ',plcl,klcl,tlcl,rlcl
!       write(6,*) 'p(klcl)= ',p(klcl)

!modif omp: first find saturated temp at level klcl
!           In the previous version, the parcel temperature
!           was obtained from a dry adiabat ascent all the way
!           to the level k. 




      tc = tplcl
      call establ(es,tc)
      pc = plcl
      rs=epsilo*es/(pc+(epsilo-1.)*es)
      qs=rs/(1.+rs)
      tv=tc*(1.+.61*qs)
      dp=p(klcl)-plcl
      rc=(1.-qs)*rd+qs*rv
!        write(6,*) 'tv= ',tv
      pb=(p(klcl)+plcl)/2.
      fact1=rd/cpd
      fact2=tv+(rl*qs/rc)
!         write(6,*) 'fact1,fact2,rc= ',fact1,fact2,rc
      fact1=fact1*fact2
      fact3=epsilo*(rl**2)*es/(cpd*pb*rv*(tv**2))
!        write(6,*) 'fact1,fact3= ',fact1,fact3
      fact3=1.+fact3
      dtdp=fact1/fact3
!         write(6,*) 'dtdp= ',dtdp
      tc=tc+dtdp*alog((pc+dp)/pc)
!         write(6,*) 'tc,t= ',tc,t(k+1)
      tpc(klcl)=tc
      rpc(klcl)=rs
!         write(6,*) 'p,r,rs= ',p(k+1),r(k+1),rs

!       tc=tpc(klcl)
!end modif omp


! omp note: the adiabat is computed by a forward integeration of the
! lapse rate. Thiscould be improved at coarse resolution by implementing 
! a 2nd or 3rd order Runge-Kunta scheme.

       plzb=0.
       do k=klcl,nlevm
          qe=r(k)/(1.+r(k))
          tve=t(k)*(1.+.61*qe)
          call establ(es,tc)
          pc=p(k)
          rs=epsilo*es/(pc+(epsilo-1.)*es)
          qs=rs/(1.+rs)
          tv=tc*(1.+.61*qs)
!          write(6,*) 'k,tv,tve= ',k,tv,tve
          ieq=iequ(tv,tve)
          if ((ieq .gt. 0) .and. (.not. capepos)) then
             capepos=.true.
          end if
          if ((ieq .lt. 0) .and. (capepos)) then
             klzb=k
             plzb=(p(k)+p(k-1))/2.
!             write(6,*) 'klzb,plzb,p(klzb)= ',klzb,plzb,p(klzb)
             go to 12
          end if
          dp=p(k+1)-p(k)
          rc=(1.-qs)*rd+qs*rv
!           write(6,*) 'tv= ',tv
          pb=(p(k)+p(k+1))/2.
          fact1=rd/cpd
          fact2=tv+(rl*qs/rc)
!          write(6,*) 'fact1,fact2,rc= ',fact1,fact2,rc
          fact1=fact1*fact2
          fact3=epsilo*(rl**2)*es/(cpd*pb*rv*(tv**2))
!          write(6,*) 'fact1,fact3= ',fact1,fact3
          fact3=1.+fact3
          dtdp=fact1/fact3
!          write(6,*) 'dtdp= ',dtdp
          tc=tc+dtdp*alog((pc+dp)/pc)
!          write(6,*) 'tc,t= ',tc,t(k+1)
          tpc(k+1)=tc
          rpc(k+1)=rs
!          write(6,*) 'p,r,rs= ',p(k+1),r(k+1),rs
       end do
 12    continue
      ieq=iequ(plzb,0.)
      if (ieq .eq. 0) then
         xcape = 0.
         cin = 0.
         tot = 0.
         tpcback = tback
         rpcback = rback
!         write(6,*) 'plzb=0'
         return
      end if
      cin=0.
      xcape=0.
      tot=0.
!
!     Calculate convective inhibition.
!
       klfc=0
       do k=istart,nlevm
          ieq=iequ(p(k),plzb)
          if (ieq .le. 0) then
!             write(6,*) 'cin= ',cin
!            write(6,*) 'cape = 0 NO LFC'
             return 
          end if
          rbc=(rpc(k)+rpc(k+1))/2.
          rbe=(r(k)+r(k+1))/2.
          qc=rbc/(1.+rbc)
          qe=rbe/(1.+rbe)
          tvc=tpc(k)*(1.+.61*qc)
          tve=t(k)*(1.+.61*qe)
!          write(6,*) 'k,tvc,tve= ',k,tvc,tve
          ieq=iequ(tvc,tve)
          ieqa=iequ(p(k),plcl)
          if ((ieq .le. 0) .or. (ieqa .ge. 0)) then
             delt=rd*(tvc-tve)*alog(p(k)/p(k+1))
             cin=cin-delt
          else
             klfc=k
             go to 14
          end if
       end do
 14    continue
!
!      Calculate convective available potential energy.
!
!       write(6,*) 'klfc,p(klfc)= ',klfc,p(klfc)
!
! omp note: CAPE is calculated using full levels. This can create a
! significant amount of flickering, espcially when the LCL and LZB
! switch from one model level to another. 
!
       if (klfc .eq. 0) then
          xcape=0.
!          write(6,*) 'klfc=0'
          return
       end if
       do k=klfc,klzb
          ieq=iequ(p(k+1),plzb)
          if (ieq .ge. 0) then
             rbc=(rpc(k)+rpc(k+1))/2.
             rbe=(r(k)+r(k+1))/2.
             qc=rbc/(1.+rbc)
             qe=rbe/(1.+rbe)
             tvc=tpc(k)*(1.+.61*qc)
             tve=t(k)*(1.+.61*qe)
             ieq=iequ(tvc,tve)
             if (ieq .gt. 0) then
                delt=rd*(tvc-tve)*alog(p(k)/p(k+1))
!                write(6,*) 'cape k,delt,xcape= ',k,delt,xcape
                xcape=xcape+delt
                if (xcape .lt. 0.) then
!                   write(6,*) 'xcape error'
                    call error_mesg ('bm_massflux:capecalc', &
                                     'xcape error', FATAL)
!                   stop
                end if
              end if
          end if
       end do
       tot=xcape-cin
!       write(6,*) 'cin= ',cin,' J/kg'
!       write(6,*) 'xcape= ',xcape,' J/kg'
!       write(6,*) 'tot= ',tot,' J/kg'
       do k=1,nlev
          tpcback(k) = tpc(nlev+1-k)
          rpcback(k) = rpc(nlev+1-k)
       end do
       klzbback = nlev + 1 - klzb
      klclback = nlev + 1 - klcl

       return
       end subroutine capecalc

!###############################################################
!all new cape calculation.

      subroutine capecalcnew(kx,p,phalf,cp_air,rdgas,rvgas,hlv,kappa,tin,rin,&
                             avgbl,cape,cin,tp,rp,klzb,klcl)

!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp_air      specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron, 
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!    avgbl       if true, the parcel is averaged in theta and r up to its LCL
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set 
!                to zero)
!    tp          Parcel temperature (set to the environmental temperature 
!                where no adjustment)
!    rp          Parcel specific humidity (set to the environmental humidity 
!                where no adjustment, and set to the saturation humidity at 
!                the parcel temperature below the LCL)
!    klzb        Level of zero buoyancy
!    klcl        Lifting condensation level
!
!    Algorithm: 
!    Start with surface parcel. 
!    Calculate the lifting condensation level (uses an analytic formula and a 
!       lookup table).  
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.  
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.  
      implicit none
      integer, intent(in)                    :: kx
      logical, intent(in)                    :: avgbl
      real, intent(in), dimension(:)         :: p, phalf, tin, rin
      real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp_air
      integer, intent(out)                   :: klzb, klcl
      real, intent(out), dimension(:)        :: tp, rp
      real, intent(out)                      :: cape, cin

      integer            :: k, klfc, klcl2
      logical            :: nocape
      real, dimension(kx)   :: theta
      real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
                               a, b, dtdlnp, thetam, rm, tlcl2, &
                               plcl2, plcl, plzb

      pstar = 1.e5

      nocape = .true.
      cape = 0.
      cin = 0.
      plcl = 0.
      plzb = 0.
      klfc = 0
      klcl = 0
      klzb = 0
      tp(1:kx) = tin(1:kx)
      rp(1:kx) = rin(1:kx)

! start with surface parcel
      t0 = tin(kx)
      r0 = rin(kx)
! calculate the lifting condensation level by the following:
! are you saturated to begin with?  
!      call establ(es,t0)
      call escomp(t0,es)
      rs = rdgas/rvgas*es/p(kx)
      if (r0.ge.rs) then
! if you¹re already saturated, set lcl to be the surface value.
         plcl = p(kx)
! the first level where you¹re completely saturated.
         klcl = kx
! saturate out to get the parcel temp and humidity at this level
! first order (in delta T) accurate expression for change in temp
         tp(kx) = t0 + (r0 - rs)/(cp_air/hlv + hlv*rs/rvgas/t0**2.)
!         call establ(es,tp(kx))
         call escomp(tp(kx),es)
         rp(kx) = rdgas/rvgas*es/p(kx)
      else
! if not saturated to begin with, use the analytic expression to calculate the 
! exact pressure and temperature where you¹re saturated.  
         theta0 = tin(kx)*(pstar/p(kx))**kappa
! the expression that we utilize is log(r/theta**(1/kappa)*pstar*rvgas/rdgas) =
! log(es/T**(1/kappa))
! The right hand side of this is only a function of temperature, therefore 
! this is put into a lookup table to solve for temperature.  
         if (r0.gt.0.) then
            value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas)
            call lcltabl(value,tlcl)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
! just in case plcl is very high up
            if (plcl.lt.p(1)) then
               plcl = p(1)
               tlcl = theta0*(plcl/pstar)**kappa
               write (*,*) 'hi lcl'
            end if
            k = kx
         else
! if the parcel sp hum is zero or negative, set lcl to 2nd to top level
            plcl = p(2)
            tlcl = theta0*(plcl/pstar)**kappa
!            write (*,*) 'zero r0', r0
            do k=2,kx
               tp(k) = theta0*(p(k)/pstar)**kappa
               rp(k) = 0.
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            end do
            go to 11
         end if
! calculate the parcel temperature (adiabatic ascent) below the LCL.  
! the mixing ratio stays the same
!!! the right command??
         do while (p(k).gt.plcl)
            tp(k) = theta0*(p(k)/pstar)**kappa
!            call establ(es,tp(k))
            call escomp(tp(k),es)
            rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            k = k-1
         end do
! first level where you're saturated at the level
         klcl = k
         if (klcl.eq.1) klcl = 2
! do a saturated ascent to get the parcel temp at the LCL.  
! use your 2nd order equation up to the pressure above.  
! moist adaibat derivatives: (use the lcl values for temp, humid, and 
! pressure)
         a = kappa*tlcl + hlv/cp_air*r0
         b = hlv**2.*r0/cp_air/rvgas/tlcl**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
! second order in p (RK2)
! first get temp halfway up 
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call escomp(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
         a = kappa*tp(klcl) + hlv/cp_air*rp(klcl)
         b = hlv**2./cp_air/rvgas*rp(klcl)/tp(klcl)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
!         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
!                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp_air/ &
!                   (1. + b)
! second order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
!             p(klcl)/plcl))**2.
!         call establ(es,tp(klcl))
         call escomp(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/p(klcl)
!         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
! CAPE/CIN stuff
         if ((tp(klcl).lt.tin(klcl)).and.nocape) then
! if you¹re not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(klcl) - &
                 tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
         else
! if you¹re buoyant, then add to cape
            cape = cape + rdgas*(tp(klcl) - &
                  tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it¹s the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = klcl
            endif
         end if
      end if
! then average the properties over the boundary layer if so desired.  to give 
! a new "parcel".  this may not be saturated at the LCL, so make sure you get 
! to a level where it is before moist adiabatic ascent!
!!!! take out all the below (between the exclamation points) if no avgbl !!!!
      if (avgbl) then
         theta(klcl:kx) = tin(klcl:kx)*(pstar/p(klcl:kx))**kappa
         thetam = 0.
         rm = 0.
         do k=klcl,kx
            thetam = thetam + theta(k)*(phalf(k+1) - phalf(k))
            rm = rm + rin(k)*(phalf(k+1) - phalf(k))
         end do
         thetam = thetam/(phalf(kx+1) - phalf(klcl))
         rm = rm/(phalf(kx+1) - phalf(klcl))
! check if you¹re saturated at the top level.  if not, then get a new LCL
         tp(klcl) = thetam*(p(klcl)/pstar)**kappa
!         call establ(es,tp(klcl))
         call escomp(tp(klcl),es)
         rs = rdgas/rvgas*es/p(klcl)
! if you¹re not saturated, get a new LCL
         if (rm.lt.rs) then
! reset CIN to zero.  
            cin = 0.
! again, use the analytic expression to calculate the exact pressure and 
! temperature where you¹re saturated.  
! the expression that we utilize is log(r/theta**(1/kappa)*pstar*rvgas/rdgas)=
! log(es/T**(1/kappa))
! The right hand side of this is only a function of temperature, therefore 
! this is put into a lookup table to solve for temperature.  
            value = log(thetam**(-1/kappa)*rm*pstar*rvgas/rdgas)
            call lcltabl(value,tlcl2)
            plcl2 = pstar*(tlcl2/thetam)**(1/kappa)
! just in case plcl is very high up
            if (plcl2.lt.p(1)) then
               plcl2 = p(1)
            end if
            k = kx
! calculate the parcel temperature (adiabatic ascent) below the LCL.  
! the mixing ratio stays the same
!!! the right command??
            do while (p(k).gt.plcl2) 
               tp(k) = thetam*(p(k)/pstar)**kappa
!               call establ(es,tp(k))
               call escomp(tp(k),es)
               rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
               k = k-1
            end do
! first level where you're saturated at the level
            klcl2 = k
            if(klcl2.eq.1) klcl2 = 2
! do a saturated ascent to get the parcel temp at the LCL.  
! use your 2nd order equation up to the pressure above.  
! moist adaibat derivatives: (use the lcl values for temp, humid, and 
! pressure)
            a = kappa*tlcl2 + hlv/cp_air*rm
            b = hlv**2.*rm/cp_air/rvgas/tlcl2**2.
            dtdlnp = a/(1. + b)
! first order in p
!            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
! second order in p (RK2)
! first get temp halfway up 
         tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)/2.
         if ((tp(klcl2).lt.173.16).and.nocape) go to 11
         call escomp(tp(klcl2),es)
         rp(klcl2) = rdgas/rvgas*es/(p(klcl2) + plcl2)*2.
         a = kappa*tp(klcl2) + hlv/cp_air*rp(klcl2)
         b = hlv**2./cp_air/rvgas*rp(klcl2)/tp(klcl2)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
!            d2tdlnp2 = (kappa + b - 1. - b/tlcl2*(hlv/rvgas/tlcl2 - &
!                          2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*rm/cp_air/ &
!                          (1. + b)
! second order in p
!            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2) + &
!               .5*d2tdlnp2*(log(p(klcl2)/plcl2))**2.
!            call establ(es,tp(klcl2))
            call escomp(tp(klcl2),es)
            rp(klcl2) = rdgas/rvgas*es/p(klcl2)
! CAPE/CIN stuff
            if ((tp(klcl2).lt.tin(klcl2)).and.nocape) then
! if you¹re not yet buoyant, then add to the CIN and continue
               cin = cin + rdgas*(tin(klcl2) - &
                    tp(klcl2))*log(phalf(klcl2+1)/phalf(klcl2))
            else
! if you¹re buoyant, then add to cape
               cape = cape + rdgas*(tp(klcl) - &
                     tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it¹s the first time buoyant, then set the level of free convection to k
               if (nocape) then
                  nocape = .false.
                  klfc = klcl2
               endif
            end if
         end if
      end if
!!!! take out all of the above (within the exclamations) if no avgbl !!!!
! then, start at the LCL, and do moist adiabatic ascent by the first order 
! scheme -- 2nd order as well
      do k=klcl-1,1,-1
         a = kappa*tp(k+1) + hlv/cp_air*rp(k+1)
         b = hlv**2./cp_air/rvgas*rp(k+1)/tp(k+1)**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
! second order in p (RK2)
! first get temp halfway up 
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call escomp(tp(k),es)
         rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
         a = kappa*tp(k) + hlv/cp_air*rp(k)
         b = hlv**2./cp_air/rvgas*rp(k)/tp(k)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
!         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - & 
!               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp_air*rp(k+1)/(1. + b)
! second order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
!             p(k)/p(k+1)))**2.
! if you're below the lookup table value, just presume that there's no way 
! you could have cape and call it quits
         if ((tp(k).lt.173.16).and.nocape) go to 11
!         call establ(es,tp(k))
         call escomp(tp(k),es)
         rp(k) = rdgas/rvgas*es/p(k)
         if ((tp(k).lt.tin(k)).and.nocape) then
! if you¹re not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
         elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
! if you have CAPE, and it¹s your first time being negatively buoyant, 
! then set the level of zero buoyancy to k+1, and stop the moist ascent
            klzb = k+1
            go to 11
         else
! if you¹re buoyant, then add to cape
            cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
! if it¹s the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = k
            endif
         end if
      end do
 11   if(nocape) then 
! this is if you made it through without having a LZB
! set LZB to be the top level.
         plzb = p(1)
         klzb = 0
         klfc = 0
         cin = 0.
         tp(1:kx) = tin(1:kx)
         rp(1:kx) = rin(1:kx)
      end if
!      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
!      write (*,*) 'tp, rp new', tp, rp
!       write (*,*) 'tp, new', tp
!       write (*,*) 'tin new', tin 
!       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
      end subroutine capecalcnew

! lookup table with e_s(T) using the new analytic expression
      subroutine establ2(es,t)
      

! Table of es values as a function of temperature.  
! Uses the analytic expression for e_s which assumes fixed latent heat 
!    coefficient.  
! Gives the values from -100 to 60 Celsius in 1 degree increments.  

      implicit none 
      real, intent(in)     :: t
      real, intent(out)    :: es

      integer              :: it
      real, dimension(161) :: table
      real                 :: t1, t2

      data table / 6.4876769e-03,   7.7642650e-03,   9.2730105e-03, & ! XXX Too many continuation lines 
       1.1052629e-02,   1.3147696e-02,   1.5609446e-02,   1.8496657e-02,   2.1876647e-02,   2.5826384e-02,   3.0433719e-02, &
       3.5798760e-02,   4.2035399e-02,   4.9272997e-02,   5.7658264e-02,   6.7357319e-02,   7.8557979e-02,   9.1472273e-02, &
       1.0633921e-01,   1.2342784e-01,   1.4304057e-01,   1.6551683e-01,   1.9123713e-01,   2.2062738e-01,   2.5416374e-01, &
       2.9237778e-01,   3.3586224e-01,   3.8527718e-01,   4.4135673e-01,   5.0491638e-01,   5.7686092e-01,   6.5819298e-01, &
       7.5002239e-01,   8.5357615e-01,   9.7020925e-01,   1.1014164e+00,   1.2488446e+00,   1.4143067e+00,   1.5997959e+00, &
       1.8075013e+00,   2.0398249e+00,   2.2993996e+00,   2.5891082e+00,   2.9121041e+00,   3.2718336e+00,   3.6720588e+00, &
       4.1168837e+00,   4.6107798e+00,   5.1586157e+00,   5.7656865e+00,   6.4377468e+00,   7.1810447e+00,   8.0023579e+00, &
       8.9090329e+00,   9.9090255e+00,   1.1010944e+01,   1.2224096e+01,   1.3558536e+01,   1.5025116e+01,   1.6635542e+01, &
       1.8402429e+01,   2.0339361e+01,   2.2460955e+01,   2.4782931e+01,   2.7322176e+01,   3.0096824e+01,   3.3126327e+01, &
       3.6431545e+01,   4.0034823e+01,   4.3960087e+01,   4.8232935e+01,   5.2880735e+01,   5.7932732e+01,   6.3420149e+01, &
       6.9376307e+01,   7.5836738e+01,   8.2839310e+01,   9.0424352e+01,   9.8634795e+01,   1.0751630e+02,   1.1711742e+02, &
       1.2748974e+02,   1.3868802e+02,   1.5077039e+02,   1.6379851e+02,   1.7783773e+02,   1.9295728e+02,   2.0923048e+02, &
       2.2673493e+02,   2.4555268e+02,   2.6577049e+02,   2.8748004e+02,   3.1077813e+02,   3.3576694e+02,   3.6255429e+02, &
       3.9125382e+02,   4.2198536e+02,   4.5487511e+02,   4.9005594e+02,   5.2766770e+02,   5.6785749e+02,   6.1078000e+02, &
       6.5659776e+02,   7.0548154e+02,   7.5761062e+02,   8.1317317e+02,   8.7236659e+02,   9.3539788e+02,   1.0024840e+03, &
       1.0738523e+03,   1.1497408e+03,   1.2303987e+03,   1.3160868e+03,   1.4070779e+03,   1.5036572e+03,   1.6061228e+03, &
       1.7147860e+03,   1.8299721e+03,   1.9520206e+03,   2.0812857e+03,   2.2181372e+03,   2.3629602e+03,   2.5161565e+03, &
       2.6781448e+03,   2.8493609e+03,   3.0302589e+03,   3.2213112e+03,   3.4230097e+03,   3.6358656e+03,   3.8604109e+03, &
       4.0971982e+03,   4.3468019e+03,   4.6098188e+03,   4.8868684e+03,   5.1785938e+03,   5.4856626e+03,   5.8087673e+03, &
       6.1486259e+03,   6.5059830e+03,   6.8816104e+03,   7.2763077e+03,   7.6909031e+03,   8.1262545e+03,   8.5832496e+03, &
       9.0628075e+03,   9.5658788e+03,   1.0093447e+04,   1.0646529e+04,   1.1226176e+04,   1.1833474e+04,   1.2469546e+04, &
       1.3135552e+04,   1.3832687e+04,   1.4562188e+04,   1.5325331e+04,   1.6123432e+04,   1.6957848e+04,   1.7829980e+04, &
       1.8741270e+04,   1.9693207e+04,   2.0687323e+04,   2.1725199e+04 /

      t1 = t
      if (t.lt.173.16) t1 = 173.16
      if (t.gt.333.16) t1 = 333.16
      it = floor(t1 - 173.16)
      t2 = 173.16 + it
      es = (t2 + 1.0 - t1)*table(it+1) + (t1 - t2)*table(it+2)
      end subroutine establ2

! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.  
! 
! the expression that we utilize is log(r/theta**(1/kappa)*pstar*rvgas/rdgas) = 
! log(es/T**(1/kappa))
! 
! Gives the values of the temperature for the following range: 
!   starts with -23, is uniformly distributed up to -10.4.  There are a 
! total of 127 values, and the increment is .1.  
!
      implicit none 
      real, intent(in)     :: value
      real, intent(out)    :: tlcl

      integer              :: ival
      real, dimension(127) :: lcltable
      real                 :: v1, v2

      data lcltable/  1.7364512e+02,   1.7427449e+02,   1.7490874e+02, &
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02, &
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02, &
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02, &
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02, &
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02, &
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02, &
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02, &
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02, &
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02, &
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02, &
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02, &
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02, &
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02, &
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02, &
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02, &
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02, &
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02, &
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02, &
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02, &
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02, &
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02, &
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02, &
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02, &
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02, &
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02, &
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02, &
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02, &
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02, &
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02, &
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02, &
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02 /

      v1 = value
      if (value.lt.-23.0) v1 = -23.0
      if (value.gt.-10.4) v1 = -10.4
      ival = floor(10.*(v1 + 23.0))
      v2 = -230. + ival
      v1 = 10.*v1
      tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(ival+2)


      end subroutine lcltabl



!#######################################################################


      SUBROUTINE ESTABL(ES,TP)
!
!   TABLE OF ES FROM -100 TO +60 C IN ONE-DEGREE INCREMENTS(ICE).
!
!   RAT GIVES THE RATIO OF ES(ICE)/ES(LIQUID)
!
!  es refers to liquid above 273, ice below 273
!
!

     implicit none
     real, intent(in)  :: TP
     real, intent(out) :: ES

     integer :: it
     real, dimension(161) :: table
     real :: ft, t2, tp1

!      DIMENSION TABLE(161)

      DATA TABLE/.01403,.01719,.02101,.02561,.03117,.03784, &
      .04584,.05542,.06685,.08049,.09672,.1160,.1388,.1658, &
      .1977,.2353,.2796,.3316,.3925,.4638,.5472,.6444,.7577, &
      .8894,1.042,1.22,1.425,1.622,1.936,2.252,2.615,3.032, &
      3.511,4.06,4.688,5.406,6.225,7.159,8.223,9.432,10.80, &
      12.36,14.13,16.12,18.38,20.92,23.80,27.03,30.67,34.76, &
      39.35,44.49,50.26,56.71,63.93,71.98,80.97,90.98,102.1, &
      114.5,128.3,143.6,160.6,179.4,200.2,223.3,248.8,276.9, &
      307.9,342.1,379.8,421.3,466.9,517.0,572.0,632.3,698.5, &
      770.9,850.2,937.0,1032.0,1146.6,1272.0,1408.1,1556.7, &
      1716.9,1890.3,2077.6,2279.6,2496.7,2729.8,2980.,3247.8, &
      3534.1,3839.8,4164.8,4510.5,4867.9,5265.1,5675.2,6107.8, &
      6566.2,7054.7,7575.3,8129.4,8719.2,9346.5,10013.,10722., &
      11474.,12272.,13119.,14017.,14969.,15977.,17044.,18173., &
      19367.,20630.,21964.,23373.,24861.,26430.,28086.,29831., &
      31671.,33608.,35649.,37796.,40055.,42430.,44927.,47551., &
      50307.,53200.,56236.,59422.,62762.,66264.,69934.,73777., &
      77802.,82015.,86423.,91034.,95855.,100890.,106160., &
      111660.,117400.,123400.,129650.,136170.,142980.,150070., &
      157460.,165160.,173180.,181530.,190220.,199260./ 

!      DATA TABLE/.01403,.01719,.02101,.02561,.03117,.03784,
!     A .04584,.05542,.06685,.08049,.09672,.1160,.1388,.1658,
!     B .1977,.2353,.2796,.3316,.3925,.4638,.5472,.6444,.7577,
!     C .8894,1.042,1.22,1.425,1.622,1.936,2.252,2.615,3.032,
!     D 3.511,4.06,4.688,5.406,6.225,7.159,8.223,9.432,10.80,
!     E 12.36,14.13,16.12,18.38,20.92,23.80,27.03,30.67,34.76,
!     F 39.35,44.49,50.26,56.71,63.93,71.98,80.97,90.98,102.1,
!     G 114.5,128.3,143.6,160.6,179.4,200.2,223.3,248.8,276.9,
!     H 307.9,342.1,379.8,421.3,466.9,517.0,572.0,632.3,698.5,
!     I 770.9,850.2,937.0,1032.0,1146.6,1272.0,1408.1,1556.7,
!     J 1716.9,1890.3,2077.6,2279.6,2496.7,2729.8,2980.,3247.8,
!     K 3534.1,3839.8,4164.8,4510.5,4867.9,5265.1,5675.2,6107.8,
!     L 6566.2,7054.7,7575.3,8129.4,8719.2,9346.5,10013.,10722.,
!     M 11474.,12272.,13119.,14017.,14969.,15977.,17044.,18173.,
!     N 19367.,20630.,21964.,23373.,24861.,26430.,28086.,29831.,
!     O 31671.,33608.,35649.,37796.,40055.,42430.,44927.,47551.,
!     P 50307.,53200.,56236.,59422.,62762.,66264.,69934.,73777.,
!     Q 77802.,82015.,86423.,91034.,95855.,100890.,106160.,
!     R 111660.,117400.,123400.,129650.,136170.,142980.,150070.,
!     S 157460.,165160.,173180.,181530.,190220.,199260./
!
      tp1 = tp
      IF (TP1 .LT. 173.16) GO TO 1
      IF (TP1 .LE. 333.16) GO TO 2
      TP1=333.16
      GO TO 2
 1    TP1=173.16
 2    IT= floor(TP1-173.16)
      FT= IT
      T2=173.16+FT
      ES=(T2+1.0-TP1)*TABLE(IT+1)+(TP1-T2)*TABLE(IT+2)
      ES=ES*.1
!
!     CONVERT FROM ES(LIQUID) TO ES(ICE)
!
!      R1=EXP(28.92-(6142./TP1))
!      R2=EXP(26.27-(5421./TP1))
!      RAT=R1/R2
!                        RAT=1.
!     ES=ES*R1/R2
      RETURN
      END subroutine establ

!#######################################################################


      integer function iequ(x,y)
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!

      real, intent(in) :: x, y
      real :: eps, epsm, d

      iequ=0
      eps=1.e-10
      epsm=-eps
      d=x-y
      if (d .gt. eps) iequ=10
      if (d .lt. epsm) iequ=-10
      return
      end function iequ


!#######################################################################

   subroutine bm_massflux_init ()

!-----------------------------------------------------------------------
!
!        initialization for bm_massflux
!
!-----------------------------------------------------------------------

  integer  unit,io,ierr, logunit

!----------- read namelist ---------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=bm_massflux_nml, iostat=io)
      ierr = check_nml_error(io,"bm_massflux_nml")
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=bm_massflux_nml, iostat=io, end=10)
            ierr = check_nml_error (io,'bm_massflux_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!---------- output namelist --------------------------------------------

      call write_version_number(version, tagname)
      if ( mpp_pe() == mpp_root_pe() ) then
           logunit = stdlog()
           write (logunit,nml=bm_massflux_nml)
      endif
      call close_file (unit)

      module_is_initialized =.true.

   end subroutine bm_massflux_init

!#######################################################################

   subroutine bm_massflux_end()
   
      module_is_initialized =.false.
   
   end subroutine bm_massflux_end

!#######################################################################

end module bm_massflux_mod



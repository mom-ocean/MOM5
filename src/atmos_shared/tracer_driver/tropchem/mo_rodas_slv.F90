      module mo_rodas_sol_mod

      use chem_mods_mod, only : clscnt5

      implicit none

      private
      public :: rodas_slv_init, rodas_sol

!     save

      integer :: ox_ndx
      integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                 c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, &
                 ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                 c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx, n_ndx
      real :: epsilon(max(1,clscnt5))
      real :: err_wghts(max(1,clscnt5))

character(len=128), parameter :: version     = '$Id: mo_rodas_slv.F90,v 19.0 2012/01/06 20:34:04 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine rodas_slv_init (retain_cm3_bugs)
!-----------------------------------------------------------------------      
!        ... initialize the implict solver
!-----------------------------------------------------------------------      

      use chem_mods_mod,  only : rodas
      use mo_grid_mod,    only : pcnstm1
      use mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx

      implicit none
      
      logical :: retain_cm3_bugs

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      real, parameter :: rel_err      = 1.e-2
      real, parameter :: high_rel_err = 1.e-3
      integer :: m
      real    :: eps(pcnstm1)
      real    :: wghts(pcnstm1)

      eps(:) = rel_err
!      ox_ndx = get_spc_ndx( 'OX' )
! for ox budget (jmao,1/7/2011)
   if (retain_cm3_bugs) then
      ox_ndx = get_spc_ndx( 'OX' )
   else
      ox_ndx = get_spc_ndx( 'O3' )
   endif
      if( ox_ndx > 0 ) then
         eps(ox_ndx) = high_rel_err
      else
         m = get_spc_ndx( 'O3' )
         if( m > 0 ) then
            eps(m) = high_rel_err
         end if
      end if
      m = get_spc_ndx( 'NO' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HNO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'N2O5' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'OH' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if

      wghts(:) = 1.
      n_ndx = get_spc_ndx( 'N' )
      if( n_ndx > 0 ) then
         wghts(n_ndx) = 0.
      end if
      do m = 1,rodas%clscnt
         epsilon(m)   = eps(rodas%clsmap(m))
         err_wghts(m) = wghts(rodas%clsmap(m))
      end do

      end subroutine rodas_slv_init

      subroutine rodas_sol( base_sol, reaction_rates, &
                            het_rates, extfrc, &
                            nstep, delt, &
                            plonl, plnplv )
!-----------------------------------------------------------------------
!              ... rodas_sol advances the volumetric mixing ratio
!           forward one time step via the implicit runge-kutta rosenbrock scheme
!-----------------------------------------------------------------------

      use chem_mods_mod,           only : rod_nzcnt, clscnt5, clsze, &
                                      rxntot, hetcnt, extcnt, rodas
      use mo_grid_mod,             only : pcnstm1
      use mo_indprd_mod,           only : indprd
      use mo_rodas_prod_loss_mod,  only : rodas_prod_loss
      use mo_rod_lin_matrix_mod,   only : rod_linmat
      use mo_rod_nln_matrix_mod,   only : rod_nlnmat
      use mo_rod_factor_mod,       only : rod_lu_fac
      use mo_rod_solve_mod,        only : rod_lu_slv

      implicit none

!-----------------------------------------------------------------------
!             ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::   nstep                     ! time step index (zero based)
      integer, intent(in) ::   plonl                     ! longitude tile dimension
      integer, intent(in) ::   plnplv                    ! plonl*plev
      real, intent(in)    ::   delt                      ! time step (s)
      real, intent(in)    ::   reaction_rates(plnplv,rxntot)
      real, intent(in)    ::   het_rates(plnplv,max(1,hetcnt)), &
                               extfrc(plnplv,max(1,extcnt))
      real, intent(inout) ::   base_sol(plnplv,pcnstm1)

!-----------------------------------------------------------------------
!             ... local variables
!-----------------------------------------------------------------------
      integer, parameter :: att_limit = 5
      real, parameter    :: hmin      = 1.
      real, parameter    :: min_val   = 1.e-30
      real, parameter    :: con3      = 8./3.

      integer ::   isec, j, k, m
      integer ::   lev, indx
      integer ::   attempts, failures, step_fail_cnt
      real    ::   con1, con2
      real, dimension(max(1,rod_nzcnt)) :: sys_jac, lin_jac
      real, dimension(max(1,clscnt5))   :: yn, prod, loss, &
                                                 u1, u2, u3, u4, &
                                                 ind_prd
      real, dimension(plnplv,max(1,clscnt5))  :: gl_ind_prd
      real, dimension(max(1,rxntot))    :: lrxt
      real, dimension(max(1,hetcnt))    :: lhet
      real, dimension(max(1,pcnstm1))   :: lsol, y_temp
      real, dimension(max(1,clscnt5))         :: spc_err
      real    ::   err
      real    ::   hfull, hinv, interval
      real    ::   h
      logical ::   interval_done

!-----------------------------------------------------------------------      
!        ... if there is "independent" production put it in the forcing
!        ... set the iteration invariant part of the function f(y)
!-----------------------------------------------------------------------      
      if( rodas%indprd_cnt /= 0 .or. extcnt > 0 ) then
         call indprd( 5, gl_ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,max(1,clscnt5)
            gl_ind_prd(:,m) = 0.
         end do
      end if
level_loop : &
!++lwh
!     do lev = 1,plev
      do lev = 1,plnplv/plonl
!--lwh
lon_tile_loop : &
         do isec = 1,plonl
            indx = (lev - 1)*plonl + isec
            h    = delt
            hinv = 1./h
            interval      = 0.
            step_fail_cnt = 0
            do m = 1,rxntot
               lrxt(m) = reaction_rates(indx,m) 
            end do
            if( hetcnt > 0 ) then
               do m = 1,hetcnt
                  lhet(m) = het_rates(indx,m) 
                end do
            end if
            if( rodas%indprd_cnt /= 0 .or. extcnt > 0 ) then
               do m = 1,max(1,clscnt5)
                  ind_prd(m) = gl_ind_prd(indx,m) 
               end do
            end if
!-----------------------------------------------------------------------      
!        ... full timestep loop
!-----------------------------------------------------------------------      
full_time_step_loop : &
            do
               interval_done = .false.
               failures      = 0
!-----------------------------------------------------------------------      
!        ... transfer from base to local work arrays
!-----------------------------------------------------------------------      
               do m = 1,pcnstm1
                  lsol(m)   = base_sol(indx,m) 
                  y_temp(m) = lsol(m)
               end do
!----------------------------------------------------------------------      
!        ... store values at t(n)
!-----------------------------------------------------------------------      
               do k = 1,clscnt5
                  j       = rodas%clsmap(k)
                  m       = rodas%permute(k)
                  yn(m) = lsol(j)
               end do
!-----------------------------------------------------------------------      
!        ... attemp step size
!-----------------------------------------------------------------------      
sub_step_loop : &
               do attempts = 1,att_limit
                  con1 = 2.*hinv
                  con2 = 4.*hinv
!-----------------------------------------------------------------------      
!        ... the linear component
!-----------------------------------------------------------------------      
                  if( rodas%lin_rxt_cnt > 0 ) then
                     call rod_linmat( lin_jac, lsol, lrxt, lhet )
                  end if
!-----------------------------------------------------------------------      
!        ... the non-linear component
!-----------------------------------------------------------------------      
                     call rod_nlnmat( sys_jac, lsol, lrxt, lin_jac, con1 )
!-----------------------------------------------------------------------      
!         ... factor the "system" matrix
!-----------------------------------------------------------------------      
                  call rod_lu_fac( sys_jac )
!-----------------------------------------------------------------------      
!           ... form dy/dt = prod - loss
!-----------------------------------------------------------------------      
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u1(m) = loss(m) - (prod(m) + ind_prd(m))
                     end do
                  else
                     do m = 1,clscnt5
                        u1(m) = loss(m) - prod(m)
                     end do
                  end if
                  do m = 1,clscnt5
                     u2(m) = u1(m)
                  end do
!-----------------------------------------------------------------------      
!           ... solve for the first intermediate
!-----------------------------------------------------------------------      
                  call rod_lu_slv( sys_jac, u1 )
!-----------------------------------------------------------------------      
!           ... solve for the second intermediate
!-----------------------------------------------------------------------      
                  do m = 1,clscnt5
                     u2(m) = u2(m) - con2*u1(m)
                  end do
                  call rod_lu_slv( sys_jac, u2 )
!-----------------------------------------------------------------------      
!           ... solve for the third intermediate
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(j) = yn(m) + 2.*u1(m)
                  end do
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u3(m) = loss(m) - (prod(m) + ind_prd(m) + hinv*(u1(m) - u2(m)))
                     end do
                  else
                     do m = 1,clscnt5
                        u3(m) = loss(m) - (prod(m) + hinv*(u1(m) - u2(m)))
                     end do
                  end if
                  call rod_lu_slv( sys_jac, u3 )
!-----------------------------------------------------------------------      
!           ... solve for the fourth intermediate
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(j) = yn(m) + 2.*u1(m) + u3(m)
                  end do
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u4(m) = loss(m) - (prod(m) + ind_prd(m) + hinv*(u1(m) - u2(m) - con3*u3(m)))
                     end do
                  else
                     do m = 1,clscnt5
                        u4(m) = loss(m) - (prod(m) + hinv*(u1(m) - u2(m) - con3*u3(m)))
                     end do
                  end if
                  call rod_lu_slv( sys_jac, u4 )
!-----------------------------------------------------------------------      
!           ... form y(n+1) from intermediates
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(j) = yn(m) + 2.*u1(m) + u3(m) + u4(m)
                  end do
!-----------------------------------------------------------------------      
!           ... form estimated trunc error
!-----------------------------------------------------------------------      
                  err = 0.
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     if( lsol(j) > min_val ) then
                        spc_err(k) = err_wghts(k) * u4(m) / (epsilon(k)*lsol(j))
                        err        = err + spc_err(k)*spc_err(k)
                        end if
                     end do
                  err = sqrt( err/real(clscnt5) )
                  if( err < 1. ) then
                     if( h == delt ) then
                        interval_done = .true.
                        hfull         = h
                        exit
                     end if
                     interval  = interval + h
                     h        = h * min( 10.,max( .1,1./(err**.33) ) )
                     h         = max( hmin,h )
                     hfull     = h
                     if( abs( interval - delt ) > 1.e-6*delt ) then
                        h    = min( delt-interval,h )
                        hinv = 1. / h
                     end if
                     exit
                  else
                     if( h == hmin ) then
                        interval      = interval + h
                        hfull         = h
                        step_fail_cnt = step_fail_cnt + 1
                        exit
                     end if
                     failures = failures + 1
                     if( attempts == att_limit ) then
                        interval = interval + h
                     end if
                     if( failures >= 2 ) then
                        h = .1 * h
                     else
                        h = h * min( 10.,max( .1,.5/(err**.33) ) )
                     end if
                     h = max( hmin,h )
                     h = min( delt-interval,h )
                     hinv = 1. / h
                     if( attempts == att_limit ) then
                        hfull         = h
                        step_fail_cnt = step_fail_cnt + 1
                        exit
                     end if
                     do m = 1,pcnstm1
                        lsol(m) = y_temp(m)
                     end do
                  end if
               end do sub_step_loop
               do m = 1,pcnstm1
                  base_sol(indx,m) = lsol(m)
               end do
               if( interval_done .or. abs( interval - delt ) <= 1.e-6*delt ) then
                  h = min( hfull,delt )
                  exit
               end if
            end do full_time_step_loop
         end do lon_tile_loop
      end do level_loop

      end subroutine rodas_sol

      end module mo_rodas_sol_mod

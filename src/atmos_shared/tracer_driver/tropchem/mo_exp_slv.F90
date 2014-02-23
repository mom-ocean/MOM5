      module MO_EXP_SOL_MOD

      implicit none

!     save

      private
      public :: exp_slv_init, exp_sol

      integer, parameter ::  inst = 1, avrg = 2

      integer ::  o3s_ndx, o3inert_ndx
      integer ::  oh_ndx, ho2_ndx, c2h4_ndx, c3h6_ndx, isop_ndx, &
                  mvk_ndx, macr_ndx, c10h16_ndx, no2_ndx, n2o5_ndx, &
                  no3_ndx, ox_ndx
      integer ::  jo1d_ndx, ox_l1_ndx, o1d_n2_ndx, o1d_o2_ndx, ox_l2_ndx, &
                  ox_l3_ndx, ox_l4_ndx, ox_l5_ndx, ox_l6_ndx, ox_l7_ndx, &
                  ox_l8_ndx, ox_l9_ndx, usr4_ndx, usr16_ndx, usr17_ndx
      logical ::  o3s_loss
      logical ::  class_hist_prod = .false.
      logical ::  class_hist_loss = .false.

character(len=128), parameter :: version     = '$Id: mo_exp_slv.F90,v 19.0 2012/01/06 20:33:52 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine exp_slv_init
!-----------------------------------------------------------------------      
!       ... Initialize the explicit solver
!-----------------------------------------------------------------------      

      use CHEM_MODS_MOD,   only : clscnt1, explicit
      use mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx

      implicit none

!-----------------------------------------------------------------------      
!       ... Local variables
!-----------------------------------------------------------------------      

      o3s_ndx     = get_spc_ndx( 'O3S' )
      o3inert_ndx = get_spc_ndx( 'O3INERT' )
      ox_ndx = get_spc_ndx( 'OX' )
      oh_ndx = get_spc_ndx( 'OH' )
      ho2_ndx = get_spc_ndx( 'HO2' )
      c2h4_ndx = get_spc_ndx( 'C2H4' )
      c3h6_ndx = get_spc_ndx( 'C3H6' )
      isop_ndx = get_spc_ndx( 'ISOP' )
      mvk_ndx = get_spc_ndx( 'MVK' )
      macr_ndx = get_spc_ndx( 'MACR' )
      c10h16_ndx = get_spc_ndx( 'C10H16' )
      no2_ndx = get_spc_ndx( 'NO2' )
      n2o5_ndx = get_spc_ndx( 'N2O5' )
      no3_ndx = get_spc_ndx( 'NO3' )

      jo1d_ndx = get_rxt_ndx( 'jo1d' )
      ox_l1_ndx = get_rxt_ndx( 'ox_l1' )
      ox_l2_ndx = get_rxt_ndx( 'ox_l2' )
      ox_l3_ndx = get_rxt_ndx( 'ox_l3' )
      ox_l4_ndx = get_rxt_ndx( 'ox_l4' )
      ox_l5_ndx = get_rxt_ndx( 'ox_l5' )
      ox_l6_ndx = get_rxt_ndx( 'ox_l6' )
      ox_l7_ndx = get_rxt_ndx( 'ox_l7' )
      ox_l8_ndx = get_rxt_ndx( 'ox_l8' )
      ox_l9_ndx = get_rxt_ndx( 'ox_l9' )
      o1d_n2_ndx = get_rxt_ndx( 'o1d_n2' )
      o1d_o2_ndx = get_rxt_ndx( 'o1d_o2' )
      usr4_ndx = get_rxt_ndx( 'usr4' )
      usr16_ndx = get_rxt_ndx( 'usr16' )
      usr17_ndx = get_rxt_ndx( 'usr17' )

!-----------------------------------------------------------------------      
!       ... Scan for class production to history file(s)
!-----------------------------------------------------------------------      
!     do file = 1,moz_file_cnt
!        do timetype = inst,avrg
!           if( hfile(file)%histout_cnt(14,timetype) > 0 ) then
!              il = hfile(file)%histout_ind(14,timetype)
!              iu = il + hfile(file)%histout_cnt(14,timetype) - 1
!              if( timetype == inst ) then
!                 if( ANY( hfile(file)%inst_map(il:iu)/1000 == 1 ) ) then
!                    class_hist_prod = .true.
!                    exit
!                 end if
!              else if( timetype == avrg ) then
!                 if( ANY( hfile(file)%timav_map(il:iu)/1000 == 1 ) ) then
!                    class_hist_prod = .true.
!                    exit
!                 end if
!              end if
!           end if
!        end do
!        if( class_hist_prod ) then
!           exit
!        end if
!     end do
!-----------------------------------------------------------------------      
!       ... Scan for class loss to history file(s)
!-----------------------------------------------------------------------      
!     do file = 1,moz_file_cnt
!        do timetype = inst,avrg
!           if( hfile(file)%histout_cnt(15,timetype) > 0 ) then
!              il = hfile(file)%histout_ind(15,timetype)
!              iu = il + hfile(file)%histout_cnt(15,timetype) - 1
!              if( timetype == inst ) then
!                 if( ANY( hfile(file)%inst_map(il:iu)/1000 == 1 ) ) then
!                    class_hist_loss = .true.
!                    exit
!                 end if
!              else if( timetype == avrg ) then
!                 if( ANY( hfile(file)%timav_map(il:iu)/1000 == 1 ) ) then
!                    class_hist_loss = .true.
!                    exit
!                 end if
!              end if
!           end if
!        end do
!        if( class_hist_loss ) then
!           exit
!        end if
!     end do

      end subroutine EXP_SLV_INIT

      subroutine EXP_SOL( base_sol, reaction_rates, &
                          het_rates, extfrc, &
                          nstep, delt, &
                          prod_out, loss_out,&
                          plonl, plnplv )
!-----------------------------------------------------------------------
!       ... Exp_sol advances the volumetric mixing ratio
!           forward one time step via the fully explicit
!           Euler scheme
!           Note : This code has o3inert and o3s as the last
!                  two class members;  neither has production
!                  or loss - some dimensionality below has been
!                  altered to acount for this
!-----------------------------------------------------------------------

      use chem_mods_mod,        only : clscnt1, explicit, extcnt, hetcnt, rxntot
      use MO_INDPRD_MOD,        only : INDPRD
      use MO_EXP_PROD_LOSS_MOD, only : EXP_PROD_LOSS
      use mo_grid_mod,          only : pcnstm1

      implicit none
!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  nstep            ! time step index
      integer, intent(in) ::  plonl            ! lon tile dim
      integer, intent(in) ::  plnplv           ! plonl*plev
      real, intent(in)    ::  delt             ! time step in seconds
      real, intent(in)    ::  reaction_rates(plnplv,max(1,rxntot))
      real, intent(in)    ::  het_rates(plnplv,max(1,hetcnt)), &
                              extfrc(plnplv,max(1,extcnt))
      real, intent(inout) ::  base_sol(plnplv,pcnstm1)
      real, intent(out), optional :: prod_out(plnplv,pcnstm1),loss_out(plnplv,pcnstm1)

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer  ::  k, l, m
      real, dimension(plnplv,max(1,clscnt1)) :: &
                   prod, &
                   loss, &
                   ind_prd

      if( explicit%indprd_cnt /= 0 ) then
!-----------------------------------------------------------------------      
!        ... Put "independent" production in the forcing
!-----------------------------------------------------------------------      
         call indprd( 1, ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,max(1,clscnt1)
            ind_prd(:,m) = 0.
         end do
      end if
!-----------------------------------------------------------------------      
!       ... Form F(y)
!-----------------------------------------------------------------------      
      call exp_prod_loss( prod, loss, base_sol, reaction_rates, het_rates )

!-----------------------------------------------------------------------      
!       ... Solve for the mixing ratio at t(n+1)
!-----------------------------------------------------------------------      
      do m = 1,clscnt1
         l             = explicit%clsmap(m)
         if( l /= o3s_ndx .and. l /= o3inert_ndx ) then
         base_sol(:,l) = base_sol(:,l) + delt * (prod(:,m) + ind_prd(:,m) - loss(:,m))
!++van : o3s is assigned in mo_chemdr.F90 so commenting these lines here
!        else if( l == o3s_ndx ) then
!-----------------------------------------------------------------------      
!       ... special code for o3s
! NB: The coefficients for O3S loss from rxn with ISOP, MVK, MACR, and C10H16
!     are unity. For the OX loss rate (in IMP_SOL) they are adjusted (downward)
!     to account for the regeneration of OX by these rxns. But here, we
!     consider this regenerated OX to be "tropospheric."  -- lwh 2/01
!     Also include O3S loss from NO2+OH, N2O5+aerosol, NO3+aerosol
!-----------------------------------------------------------------------      
!           do k = 1,plnplv
!              loss(k,m) = &
!                 reaction_rates(k,jo1d_ndx)*reaction_rates(k,ox_l1_ndx) &
!                 /(reaction_rates(k,o1d_n2_ndx) + reaction_rates(k,o1d_o2_ndx) &
!                   + reaction_rates(k,ox_l1_ndx)) &
!               + reaction_rates(k,ox_l2_ndx)*base_sol(k,oh_ndx) &
!               + reaction_rates(k,ox_l3_ndx)*base_sol(k,ho2_ndx) &
!               + reaction_rates(k,ox_l6_ndx)*base_sol(k,c2h4_ndx) &
!               + reaction_rates(k,ox_l4_ndx)*base_sol(k,c3h6_ndx) &
!               + reaction_rates(k,ox_l5_ndx)*base_sol(k,isop_ndx) &
!               + reaction_rates(k,ox_l7_ndx)*base_sol(k,mvk_ndx) &
!               + reaction_rates(k,ox_l8_ndx)*base_sol(k,macr_ndx) &
!               + reaction_rates(k,ox_l9_ndx)*base_sol(k,c10h16_ndx) &
!               + ((reaction_rates(k,usr4_ndx)*base_sol(k,no2_ndx)*base_sol(k,oh_ndx) &
!                  + 3.*reaction_rates(k,usr16_ndx)*base_sol(k,n2o5_ndx) &
!                  + 2.*reaction_rates(k,usr17_ndx)*base_sol(k,no3_ndx)) &
!                  / max( base_sol(k,ox_ndx), 1.e-20 ))
!              base_sol(k,l) = base_sol(k,l)*exp( -delt*loss(k,m) )
!              loss(k,m) = loss(k,m) * base_sol(k,l)
!           end do
         end if
         if( PRESENT( prod_out ) ) then
            prod_out(:,l) = prod(:,m) + ind_prd(:,m)
         end if
         if( PRESENT( loss_out ) ) then
            loss_out(:,l) = loss(:,m)
         end if
      end do

!-----------------------------------------------------------------------      
!       ... Check for explicit species production and loss output
!           First check instantaneous then time averaged
!-----------------------------------------------------------------------      
!     if( class_hist_prod ) then
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(14,1) > 0 ) then
!              do n = 1,hfile(file)%histout_cnt(14,1)
!                 class = hfile(file)%inst_map(hfile(file)%histout_ind(14,1)+n-1)/1000
!                 if( class == 1 ) then
!                    cls_ndx = mod( hfile(file)%inst_map(hfile(file)%histout_ind(14,1)+n-1),1000 )
!                    fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(14,1)+n-1)
!                    wrk(:)  = (prod(:,cls_ndx) + ind_prd(:,cls_ndx)) * hnm(:)
!                    call outfld( fldname, wrk, plonl, ip, lat, file )
!                 end if
!              end do
!           end if
!           if( hfile(file)%histout_cnt(14,2) > 0 ) then
!              do n = 1,hfile(file)%histout_cnt(14,2)
!                 class = hfile(file)%timav_map(hfile(file)%histout_ind(14,2)+n-1)/1000
!                 if( class == 1 ) then
!                    cls_ndx = mod( hfile(file)%timav_map(hfile(file)%histout_ind(14,2)+n-1),1000 )
!                    fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(14,2)+n-1)
!                    wrk(:) = (prod(:,cls_ndx) + ind_prd(:,cls_ndx)) * hnm(:)
!                    call outfld( fldname, wrk, plonl, ip, lat, file )
!                 end if
!              end do
!           end if
!        end do
!     end if
!     if( class_hist_loss ) then
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(15,1) > 0 ) then
!              do n = 1,hfile(file)%histout_cnt(15,1)
!                 class = hfile(file)%inst_map(hfile(file)%histout_ind(15,1)+n-1)/1000
!                 if( class == 1 ) then
!                    cls_ndx = mod( hfile(file)%inst_map(hfile(file)%histout_ind(15,1)+n-1),1000 )
!                    fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(15,1)+n-1)
!                    l       = explicit%clsmap(cls_ndx)
!                    wrk(:)  = loss(:,cls_ndx) * hnm(:)
!                    call outfld( fldname, wrk, plonl, ip, lat, file )
!                 end if
!              end do
!           end if
!           if( hfile(file)%histout_cnt(15,2) > 0 ) then
!              do n = 1,hfile(file)%histout_cnt(15,2)
!                 class = hfile(file)%timav_map(hfile(file)%histout_ind(15,2)+n-1)/1000
!                 if( class == 1 ) then
!                    cls_ndx = mod( hfile(file)%timav_map(hfile(file)%histout_ind(15,2)+n-1),1000 )
!                    fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(15,2)+n-1)
!                    l       = explicit%clsmap(cls_ndx)
!                    wrk(:)  = loss(:,cls_ndx) * hnm(:)
!                    call outfld( fldname, wrk, plonl, ip, lat, file )
!                 end if
!              end do
!           end if
!        end do
!     end if

      end subroutine EXP_SOL

      end module MO_EXP_SOL_MOD

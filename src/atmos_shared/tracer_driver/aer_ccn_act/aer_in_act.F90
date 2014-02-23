        module aer_in_act_mod

implicit none
      private 
      public Jhomo_wat, Jhete_dep, Jhomo_aer, &
                       aer_in_act_init, aer_in_act_end

!Parameters for look-up tables

integer, parameter :: tpDIM = 3 ! Dimension of temperature (K)
integer, parameter :: tp2DIM = 3 ! Dimension of temperature (K)
integer, parameter :: msDIM = 5 ! Dimension of scaled sub-micron dust mass
integer, parameter :: upDIM = 5 ! Dimension of updraft velocity (m/s)
real, parameter :: unitmass = 0.11 ! (ug/m3), one unit mass

real, dimension(tpDIM) :: tp = (/243.15, 253.15, 263.15/)
real, dimension(tp2DIM) :: tp2 = (/233.15, 223.15, 213.15/)   
real, dimension(tpDIM) :: a1 = (/1.416, 1.5809, 1.5171/)
real, dimension(tpDIM) :: a2 = (/3.7909, 6.4407, 8.184/)   
real, dimension(tp2DIM) :: b1 = (/0.1068, 0., -0.1016/)
real, dimension(tp2DIM) :: b2 = (/2.2905, 1.3846, 0.5668/)
real, dimension(tp2DIM) :: b3 = (/6.1191, 4.9161, 4.1684/)
real, dimension(tp2DIM) :: c1 = (/7.2759, 10.532, 17.628/)
real, dimension(tp2DIM) :: c2 = (/1.4301, 1.346, 1.3038/)
real, dimension(tp2DIM) :: d1 = (/0.2423, 0.4259, 0.6952/)
real, dimension(tp2DIM) :: d2 = (/-1.1814, -1.3778, -1.6693/)
real, dimension(tp2DIM) :: d3 = (/1.0016, 0.9998, 1.0114/)
real, dimension(tp2DIM) :: e1 = (/0.1412, 0.1541, 0.163/)
real, dimension(tp2DIM) :: e2 = (/0.3923, 0.3397, 0.3124/)
real, dimension(tp2DIM) :: e3 = (/-0.0819, -0.0269, 0.0009/)
real, dimension(tp2DIM) :: f = (/16., 19.6, 20.2/)
real, dimension(msDIM) :: ms = (/100.,10.,1.,0.1,0.01/) ! The last dimension is for mass threshold.
real, dimension(upDIM) :: up = (/0.1,0.05,0.01,0.005,0.001/)
real, dimension(tpDIM,msDIM,upDIM) :: crystal2

character(len=128) :: version = '$Id: aer_in_act.F90,v 19.0 2012/01/06 20:31:38 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
logical :: module_is_initialized  = .false.

contains

subroutine Jhete_dep (temp,Si,concen_dust_sub,crystal)
   real, intent(in) :: temp,Si,concen_dust_sub
   real, intent(out) :: crystal

   real dust

   if(.not. module_is_initialized) call aer_in_act_init()
   
   crystal=0.
   dust=concen_dust_sub/unitmass
   if (temp>243.15) then
      crystal=dust*min(exp(12.96*(Si-1)-0.639)*5.73e-4, 0.0091)
   else if (temp>223.15) then
      crystal=dust*min(exp(12.96*(Si-1.1))**0.3*9.09e-3, 0.0091)
   else if (Si<(6.5217e-3*(temp-273.15)+1.6276)) then
      crystal=min(max(exp(1.5*(Si-1.11))-1.,0.),0.05)*dust*18.12
   endif
end subroutine Jhete_dep

subroutine Jhomo_aer (dtcloud, temp, rh, updraft, aer, crystal)
   real, dimension(3), intent(inout) :: aer
   real, intent(in) :: dtcloud, temp, rh, updraft
   real, intent(out) :: crystal

   real totalnumber, ai, aw, J, Nhomo
   integer no_tp2
   
   crystal=0.
   
   if (temp<=237.15) then
      if (temp>228.15) then
         no_tp2 = 1
      else if (temp>218.15) then
         no_tp2 = 2
      else
         no_tp2 = 3
      endif
   Nhomo=c1(no_tp2)*updraft**c2(no_tp2)
!density sulfate 1.7418e3, sea-salt 2.17e3, OC 1.362e3      
!here only particles between 0.1 and 1 micron are considered. Smaller particles do not
!have a chance to freeze.
!Size dist. is the same as used in the prognostic paper.
   totalnumber=19.1*(aer(1)/1.48e-13 + aer(2)/1.85e-13 + aer(3)/1.16e-13)
   
   ai=exp((210368.+131.438*temp-3.32373e6/temp-41729.1*log(temp))/(8.314*temp))
   aw=min(rh-ai,0.34)
!homo. nucl. rate (cm-3 sec-1); from Koop et al. (2000)
   if(aw>0.26) then
      J=10.**(-906.7+8502*aw-26924*aw**2.+29180*aw**3.);
   else
      J=0;
   endif
!the volume-average diamter is 0.204 micron (4.46e-15 cm3)
   crystal=min(Nhomo,totalnumber*(1.-exp(-1.*J*dtcloud*4.46e-15)))
   endif
end subroutine Jhomo_aer

subroutine Jhomo_wat (T, J)
!return homogeneous nucleation rate constant J (cm-3 s-1) as a
!function of temperature T (K)
   real, intent(in) :: T
   real, intent(out) :: J
   real :: TT
   
   TT = T - 273.15

   if (TT>=-50 .and. TT<-30) then
      J=-0.0001536*TT**4-0.0265*TT**3
      J=J-1.7439*TT**2-52.6611*TT-606.3952
      J=10**J
   else if (TT<-50) then
      TT=-50;
      J=-0.0001536*TT**4-0.0265*TT**3
      J=J-1.7439*TT**2-52.6611*TT-606.3952
      J=10**J
   else
      J=0.
   endif
end subroutine Jhomo_wat

subroutine aer_in_act_init ()
   module_is_initialized  = .true.
end subroutine aer_in_act_init

subroutine aer_in_act_end ()
   module_is_initialized  = .false.
end subroutine aer_in_act_end

end module aer_in_act_mod

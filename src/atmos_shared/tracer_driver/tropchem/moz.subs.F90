      module mo_setrxt_mod

implicit none
      private
      public :: setrxt

character(len=128), parameter :: version     = '$Id: moz.subs.F90,v 19.0 2012/01/06 20:34:16 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      contains

!++lwh
      subroutine setrxt( rate, temp, m, plonl, plev, plnplv )
!--lwh

      use chem_mods_mod, only : rxntot
      use mo_jpl_mod,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... Dummy arguments
!-------------------------------------------------------
!++lwh
      integer, intent(in) :: plonl, plev, plnplv
!--lwh
      real, intent(in)    :: temp(plonl,plev), m(plonl,plev)
      real, intent(inout) :: rate(plonl,plev,rxntot)

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      real  ::  itemp(plonl,plev), exp_fac(plonl,plev)
      real, dimension(plonl,plev) :: ko, kinf

      rate(:,:,53) = 3.5e-12
      rate(:,:,56) = 0.
      rate(:,:,64) = 1.5e-10
      rate(:,:,74) = 1.1e-10
      rate(:,:,80) = 1.8e-12
      rate(:,:,82) = 1.8e-12
      rate(:,:,96) = 1e-12
      rate(:,:,103) = 2.e-13
      rate(:,:,104) = 6.8e-14
      rate(:,:,108) = 1.e-14
      rate(:,:,114) = 2.4e-12
      rate(:,:,117) = 1.4e-11
      rate(:,:,124) = 2.4e-12
      rate(:,:,127) = 1.4e-11
      rate(:,:,130) = 5.e-12
      rate(:,:,153) = 6.8e-13
      rate(:,:,156) = 2.4e-12
      rate(:,:,160) = 4.5e-11
      rate(:,:,161) = 1.3e-16
      rate(:,:,165) = 2.4e-12
      rate(:,:,175) = 4.e-14
      rate(:,:,176) = 3.e-12
      rate(:,:,177) = 1.e-11
      rate(:,:,224) = 3.17e-8
      itemp(:,:) = 1. / temp(:,:)
      rate(:,:,43) = 8e-12 * exp( -2060. * itemp(:,:) )
      rate(:,:,44) = 2.15e-11 * exp( 110. * itemp(:,:) )
      rate(:,:,45) = 3.3e-11 * exp( 55. * itemp(:,:) )
      rate(:,:,46) = 1.63e-10 * exp( 60. * itemp(:,:) )
      exp_fac(:,:) = exp( 20. * itemp(:,:) )
      rate(:,:,47) = 6.3e-11 * exp_fac(:,:)
      rate(:,:,48) = 4.7e-11 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250. * itemp(:,:) )
      rate(:,:,49) = 3.5e-12 * exp_fac(:,:)
      rate(:,:,81) = 4.8e-11 * exp_fac(:,:)
      rate(:,:,50) = 3e-12 * exp( -1500. * itemp(:,:) )
      rate(:,:,51) = 5.1e-12 * exp( 210. * itemp(:,:) )
      exp_fac(:,:) = exp( -2450. * itemp(:,:) )
      rate(:,:,52) = 1.2e-13 * exp_fac(:,:)
      rate(:,:,212) = 8.5e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170. * itemp(:,:) )
      rate(:,:,59) = 1.5e-11 * exp_fac(:,:)
      rate(:,:,193) = 1.8e-11 * exp_fac(:,:)
      rate(:,:,61) = 1.3e-12 * exp( 380. * itemp(:,:) )
      rate(:,:,63) = 2.45e-12 * exp( -1775. * itemp(:,:) )
      exp_fac(:,:) = exp( 300. * itemp(:,:) )
      rate(:,:,65) = 2.8e-12 * exp_fac(:,:)
      rate(:,:,150) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,66) = 6.03e-13 * exp( -453. * itemp(:,:) )
      rate(:,:,67) = 2.30e-14 * exp( 677. * itemp(:,:) )
      rate(:,:,68) = 4.1e-13 * exp( 750. * itemp(:,:) )
      exp_fac(:,:) = exp( 200. * itemp(:,:) )
      rate(:,:,69) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,76) = 3e-11 * exp_fac(:,:)
      rate(:,:,89) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,105) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,128) = 2.3e-11 * exp_fac(:,:)
      rate(:,:,148) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,152) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,171) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,208) = 5.5e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900. * itemp(:,:) )
      rate(:,:,70) = 3.4e-13 * exp_fac(:,:)
      rate(:,:,85) = 6.5e-15 * exp_fac(:,:)
      rate(:,:,91) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,71) = 5.5e-12 * exp( 125. * itemp(:,:) )
      rate(:,:,75) = 2.2e-11 * exp( 120. * itemp(:,:) )
      rate(:,:,77) = 1.7e-12 * exp( -940. * itemp(:,:) )
      rate(:,:,78) = 1e-14 * exp( -490. * itemp(:,:) )
      rate(:,:,83) = 2.8e-12 * exp( -1800. * itemp(:,:) )
      rate(:,:,86) = 4.6e-13 * exp( -1156. * itemp(:,:) )
      exp_fac(:,:) = exp( 180. * itemp(:,:) )
      rate(:,:,87) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,107) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,113) = 2.2e-12 * exp_fac(:,:)
      rate(:,:,145) = 4.2e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 700. * itemp(:,:) )
      rate(:,:,88) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,102) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,115) = 8.e-13 * exp_fac(:,:)
      rate(:,:,125) = 8.e-13 * exp_fac(:,:)
      rate(:,:,146) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,151) = 8.6e-13 * exp_fac(:,:)
      rate(:,:,157) = 8.e-13 * exp_fac(:,:)
      rate(:,:,166) = 8.e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270. * itemp(:,:) )
      rate(:,:,90) = 5.6e-12 * exp_fac(:,:)
      rate(:,:,92) = 8.1e-12 * exp_fac(:,:)
      rate(:,:,195) = 7.4e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040. * itemp(:,:) )
      rate(:,:,94) = 4.3e-13 * exp_fac(:,:)
      rate(:,:,131) = 4.30e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500. * itemp(:,:) )
      rate(:,:,95) = 2.0e-12 * exp_fac(:,:)
      rate(:,:,98) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,181) = 1.87e-13 * exp_fac(:,:)
      rate(:,:,99) = 1.05e-14 * exp( -2000. * itemp(:,:) )
      rate(:,:,100) = 8.7e-12 * exp( -1070. * itemp(:,:) )
      rate(:,:,101) = 2.6e-12 * exp( 365. * itemp(:,:) )
      rate(:,:,109) = 1.6e11 * exp( -4150. * itemp(:,:) )
      rate(:,:,110) = 1.2e-14 * exp( -2630. * itemp(:,:) )
      rate(:,:,111) = 2.54e-11 * exp( 410. * itemp(:,:) )
      rate(:,:,112) = 1.55e-11 * exp( -540. * itemp(:,:) )
      exp_fac(:,:) = exp( 400. * itemp(:,:) )
      rate(:,:,116) = 5.e-13 * exp_fac(:,:)
      rate(:,:,126) = 5.e-13 * exp_fac(:,:)
      rate(:,:,167) = 5.e-13 * exp_fac(:,:)
      rate(:,:,118) = 4.13e-12 * exp( 452. * itemp(:,:) )
      rate(:,:,119) = 7.52e-16 * exp( -1521. * itemp(:,:) )
      exp_fac(:,:) = exp( 175. * itemp(:,:) )
      rate(:,:,120) = 1.86e-11 * exp_fac(:,:)
      rate(:,:,163) = 1.86e-11 * exp_fac(:,:)
      rate(:,:,121) = 4.4e-15 * exp( -2500. * itemp(:,:) )
      exp_fac(:,:) = exp( 360. * itemp(:,:) )
      rate(:,:,122) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,123) = 1.3e-13 * exp_fac(:,:)
      rate(:,:,129) = 5.3e-12 * exp_fac(:,:)
      rate(:,:,155) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,164) = 2.7e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 640. * itemp(:,:) )
      rate(:,:,132) = 1.3e-12 * exp_fac(:,:)
      rate(:,:,168) = 1.3e-12 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530. * itemp(:,:) )
      rate(:,:,133) = 4.6e-12 * exp_fac(:,:)
      rate(:,:,134) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,137) = 1.2e-11 * exp( 444. * itemp(:,:) )
      rate(:,:,138) = 9.9e-15 * exp( -730. * itemp(:,:) )
      rate(:,:,139) = 5.6e-11 * exp( -650. * itemp(:,:) )
      rate(:,:,142) = 1.5e-11 * exp( -3600. * itemp(:,:) )
      rate(:,:,143) = 2.1e-11 * exp( 100. * itemp(:,:) )
      rate(:,:,144) = 8.7e-12 * exp( -615. * itemp(:,:) )
      rate(:,:,147) = 3.75e-13 * exp( -40. * itemp(:,:) )
      rate(:,:,154) = 3.03e-12 * exp( -446. * itemp(:,:) )
      rate(:,:,158) = 8.4e-13 * exp( 830. * itemp(:,:) )
      exp_fac(:,:) = exp( -1860. * itemp(:,:) )
      rate(:,:,159) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,162) = 1.4e-12 * exp_fac(:,:)
      rate(:,:,169) = 1.90e-12 * exp( 190. * itemp(:,:) )
      rate(:,:,172) = 2.9e-12 * exp( -345. * itemp(:,:) )
      rate(:,:,173) = 6.9e-12 * exp( -230. * itemp(:,:) )
      rate(:,:,179) = 1.1e-11 * exp( -240. * itemp(:,:) )
      rate(:,:,183) = 1.7e-12 * exp( -710. * itemp(:,:) )
      rate(:,:,184) = 1.4e-10 * exp( -470. * itemp(:,:) )
      rate(:,:,186) = 2.3e-11 * exp( -200. * itemp(:,:) )
      rate(:,:,187) = 2.8e-11 * exp( 85. * itemp(:,:) )
      exp_fac(:,:) = exp( 290. * itemp(:,:) )
      rate(:,:,188) = 6.4e-12 * exp_fac(:,:)
      rate(:,:,209) = 4.1e-13 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800. * itemp(:,:) )
      rate(:,:,190) = 2.9e-12 * exp_fac(:,:)
      rate(:,:,200) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,207) = 1.7e-11 * exp_fac(:,:)
      rate(:,:,191) = 7.3e-12 * exp( -1280. * itemp(:,:) )
      rate(:,:,192) = 2.6e-12 * exp( -350. * itemp(:,:) )
      exp_fac(:,:) = exp( 220. * itemp(:,:) )
      rate(:,:,194) = 2.7e-12 * exp_fac(:,:)
      rate(:,:,214) = 5.8e-12 * exp_fac(:,:)
      rate(:,:,196) = 8.1e-11 * exp( -30. * itemp(:,:) )
      exp_fac(:,:) = exp( 260. * itemp(:,:) )
      rate(:,:,202) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,204) = 8.8e-12 * exp_fac(:,:)
      rate(:,:,203) = 4.5e-12 * exp( 460. * itemp(:,:) )
      rate(:,:,205) = 1.2e-10 * exp( -430. * itemp(:,:) )
      rate(:,:,206) = 4.8e-12 * exp( -310. * itemp(:,:) )
      rate(:,:,210) = 6.0e-13 * exp( 230. * itemp(:,:) )
      rate(:,:,211) = 4.5e-14 * exp( -1260. * itemp(:,:) )

      itemp(:,:) = 300. * itemp(:,:)

      ko(:,:) = 2.e-30 * itemp(:,:)**4.4
      kinf(:,:) = 1.4e-12 * itemp(:,:)**.7
      call jpl( rate(1,1,54), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.8e-30 * itemp(:,:)**3.0
      kinf(:,:) = 2.8e-11
      call jpl( rate(1,1,57), m, .6, ko, kinf, plnplv )

      ko(:,:) = 2.0e-31 * itemp(:,:)**3.4
      kinf(:,:) = 2.9e-12 * itemp(:,:)**1.1
      call jpl( rate(1,1,60), m, .6, ko, kinf, plnplv )

      ko(:,:) = 5.9e-33 * itemp(:,:)**1.4
      kinf(:,:) = 1.1e-12 * itemp(:,:)**(-1.3)
      call jpl( rate(1,1,72), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.5e-13 * itemp(:,:)**(-0.6)
      kinf(:,:) = 2.1e9 * itemp(:,:)**(-6.1)
      call jpl( rate(1,1,73), m, .6, ko, kinf, plnplv )

      ko(:,:) = 8.e-27 * itemp(:,:)**3.5
      kinf(:,:) = 3.e-11
      call jpl( rate(1,1,84), m, .5, ko, kinf, plnplv )

      ko(:,:) = 9.7e-29 * itemp(:,:)**5.6
      kinf(:,:) = 9.3e-12 * itemp(:,:)**1.5
      call jpl( rate(1,1,93), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.e-28 * itemp(:,:)**4.5
      kinf(:,:) = 8.8e-12 * itemp(:,:)**0.85
      call jpl( rate(1,1,106), m, .6, ko, kinf, plnplv )

      ko(:,:) = 8.e-27 * itemp(:,:)**3.5
      kinf(:,:) = 3.e-11
      call jpl( rate(1,1,174), m, .5, ko, kinf, plnplv )

      ko(:,:) = 3.3e-31 * itemp(:,:)**4.3
      kinf(:,:) = 1.6e-12
      call jpl( rate(1,1,178), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 4.4e-32 * itemp(:,:)**1.3
      kinf(:,:) = 4.7e-11 * itemp(:,:)**0.2
      call jpl( rate(1,1,185), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 1.8e-31 * itemp(:,:)**3.4
      kinf(:,:) = 1.5e-11 * itemp(:,:)**1.9
      call jpl( rate(1,1,189), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 6.9e-31 * itemp(:,:)**1.0
      kinf(:,:) = 2.6e-11
      call jpl( rate(1,1,197), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 1.6e-32 * itemp(:,:)**4.5
      kinf(:,:) = 2.0e-12 * itemp(:,:)**2.4
      call jpl( rate(1,1,198), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 5.2e-31 * itemp(:,:)**3.2
      kinf(:,:) = 6.9e-12 * itemp(:,:)**2.9
      call jpl( rate(1,1,201), m, 0.6, ko, kinf, plnplv )

      ko(:,:) = 9.0e-32 * itemp(:,:)**1.5
      kinf(:,:) = 3.0e-11
      call jpl( rate(1,1,213), m, 0.6, ko, kinf, plnplv )

      end subroutine setrxt

      end module mo_setrxt_mod

      module mo_adjrxt_mod

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, plnplv )

      use chem_mods_mod, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in)    :: inv(plnplv,nfs)
      real, intent(in)    :: m(plnplv)
      real, intent(inout) :: rate(plnplv,rxntot)

!--------------------------------------------------------------------
!       ... Local variables
!--------------------------------------------------------------------

      rate(:, 44) = rate(:, 44) * inv(:, 2)
      rate(:, 45) = rate(:, 45) * inv(:, 3)
      rate(:, 54) = rate(:, 54) * inv(:, 1)
      rate(:, 55) = rate(:, 55) * inv(:, 1)
      rate(:, 57) = rate(:, 57) * inv(:, 1)
      rate(:, 60) = rate(:, 60) * inv(:, 1)
      rate(:, 62) = rate(:, 62) * inv(:, 1)
      rate(:, 72) = rate(:, 72) * inv(:, 1)
      rate(:, 84) = rate(:, 84) * inv(:, 1)
      rate(:, 93) = rate(:, 93) * inv(:, 1)
      rate(:, 97) = rate(:, 97) * inv(:, 1)
      rate(:,106) = rate(:,106) * inv(:, 1)
      rate(:,108) = rate(:,108) * inv(:, 3)
      rate(:,135) = rate(:,135) * inv(:, 1)
      rate(:,136) = rate(:,136) * inv(:, 1)
      rate(:,142) = rate(:,142) * inv(:, 3)
      rate(:,178) = rate(:,178) * inv(:, 1)
      rate(:,189) = rate(:,189) * inv(:, 1)
      rate(:,197) = rate(:,197) * inv(:, 1)
      rate(:,198) = rate(:,198) * inv(:, 1)
      rate(:,199) = rate(:,199) * inv(:, 1)
      rate(:,201) = rate(:,201) * inv(:, 1)
      rate(:,213) = rate(:,213) * inv(:, 1)
      rate(:, 42) = rate(:, 42) * inv(:, 3) * inv(:, 1)
      rate(:,185) = rate(:,185) * inv(:, 3) * inv(:, 1)
      rate(:, 43) = rate(:, 43) * m(:)
      rate(:, 46) = rate(:, 46) * m(:)
      rate(:, 47) = rate(:, 47) * m(:)
      rate(:, 48) = rate(:, 48) * m(:)
      rate(:, 49) = rate(:, 49) * m(:)
      rate(:, 50) = rate(:, 50) * m(:)
      rate(:, 51) = rate(:, 51) * m(:)
      rate(:, 52) = rate(:, 52) * m(:)
      rate(:, 53) = rate(:, 53) * m(:)
      rate(:, 54) = rate(:, 54) * m(:)
      rate(:, 56) = rate(:, 56) * m(:)
      rate(:, 57) = rate(:, 57) * m(:)
      rate(:, 58) = rate(:, 58) * m(:)
      rate(:, 59) = rate(:, 59) * m(:)
      rate(:, 60) = rate(:, 60) * m(:)
      rate(:, 61) = rate(:, 61) * m(:)
      rate(:, 63) = rate(:, 63) * m(:)
      rate(:, 64) = rate(:, 64) * m(:)
      rate(:, 65) = rate(:, 65) * m(:)
      rate(:, 66) = rate(:, 66) * m(:)
      rate(:, 67) = rate(:, 67) * m(:)
      rate(:, 68) = rate(:, 68) * m(:)
      rate(:, 69) = rate(:, 69) * m(:)
      rate(:, 70) = rate(:, 70) * m(:)
      rate(:, 71) = rate(:, 71) * m(:)
      rate(:, 72) = rate(:, 72) * m(:)
      rate(:, 73) = rate(:, 73) * m(:)
      rate(:, 74) = rate(:, 74) * m(:)
      rate(:, 75) = rate(:, 75) * m(:)
      rate(:, 76) = rate(:, 76) * m(:)
      rate(:, 77) = rate(:, 77) * m(:)
      rate(:, 78) = rate(:, 78) * m(:)
      rate(:, 79) = rate(:, 79) * m(:)
      rate(:, 80) = rate(:, 80) * m(:)
      rate(:, 81) = rate(:, 81) * m(:)
      rate(:, 82) = rate(:, 82) * m(:)
      rate(:, 83) = rate(:, 83) * m(:)
      rate(:, 84) = rate(:, 84) * m(:)
      rate(:, 85) = rate(:, 85) * m(:)
      rate(:, 86) = rate(:, 86) * m(:)
      rate(:, 87) = rate(:, 87) * m(:)
      rate(:, 88) = rate(:, 88) * m(:)
      rate(:, 89) = rate(:, 89) * m(:)
      rate(:, 90) = rate(:, 90) * m(:)
      rate(:, 91) = rate(:, 91) * m(:)
      rate(:, 92) = rate(:, 92) * m(:)
      rate(:, 93) = rate(:, 93) * m(:)
      rate(:, 94) = rate(:, 94) * m(:)
      rate(:, 95) = rate(:, 95) * m(:)
      rate(:, 96) = rate(:, 96) * m(:)
      rate(:, 98) = rate(:, 98) * m(:)
      rate(:, 99) = rate(:, 99) * m(:)
      rate(:,100) = rate(:,100) * m(:)
      rate(:,101) = rate(:,101) * m(:)
      rate(:,102) = rate(:,102) * m(:)
      rate(:,103) = rate(:,103) * m(:)
      rate(:,104) = rate(:,104) * m(:)
      rate(:,105) = rate(:,105) * m(:)
      rate(:,106) = rate(:,106) * m(:)
      rate(:,107) = rate(:,107) * m(:)
      rate(:,110) = rate(:,110) * m(:)
      rate(:,111) = rate(:,111) * m(:)
      rate(:,112) = rate(:,112) * m(:)
      rate(:,113) = rate(:,113) * m(:)
      rate(:,114) = rate(:,114) * m(:)
      rate(:,115) = rate(:,115) * m(:)
      rate(:,116) = rate(:,116) * m(:)
      rate(:,117) = rate(:,117) * m(:)
      rate(:,118) = rate(:,118) * m(:)
      rate(:,119) = rate(:,119) * m(:)
      rate(:,120) = rate(:,120) * m(:)
      rate(:,121) = rate(:,121) * m(:)
      rate(:,122) = rate(:,122) * m(:)
      rate(:,123) = rate(:,123) * m(:)
      rate(:,124) = rate(:,124) * m(:)
      rate(:,125) = rate(:,125) * m(:)
      rate(:,126) = rate(:,126) * m(:)
      rate(:,127) = rate(:,127) * m(:)
      rate(:,128) = rate(:,128) * m(:)
      rate(:,129) = rate(:,129) * m(:)
      rate(:,130) = rate(:,130) * m(:)
      rate(:,131) = rate(:,131) * m(:)
      rate(:,132) = rate(:,132) * m(:)
      rate(:,133) = rate(:,133) * m(:)
      rate(:,134) = rate(:,134) * m(:)
      rate(:,135) = rate(:,135) * m(:)
      rate(:,137) = rate(:,137) * m(:)
      rate(:,138) = rate(:,138) * m(:)
      rate(:,139) = rate(:,139) * m(:)
      rate(:,143) = rate(:,143) * m(:)
      rate(:,144) = rate(:,144) * m(:)
      rate(:,145) = rate(:,145) * m(:)
      rate(:,146) = rate(:,146) * m(:)
      rate(:,147) = rate(:,147) * m(:)
      rate(:,148) = rate(:,148) * m(:)
      rate(:,149) = rate(:,149) * m(:)
      rate(:,150) = rate(:,150) * m(:)
      rate(:,151) = rate(:,151) * m(:)
      rate(:,152) = rate(:,152) * m(:)
      rate(:,153) = rate(:,153) * m(:)
      rate(:,154) = rate(:,154) * m(:)
      rate(:,155) = rate(:,155) * m(:)
      rate(:,156) = rate(:,156) * m(:)
      rate(:,157) = rate(:,157) * m(:)
      rate(:,158) = rate(:,158) * m(:)
      rate(:,159) = rate(:,159) * m(:)
      rate(:,160) = rate(:,160) * m(:)
      rate(:,161) = rate(:,161) * m(:)
      rate(:,162) = rate(:,162) * m(:)
      rate(:,163) = rate(:,163) * m(:)
      rate(:,164) = rate(:,164) * m(:)
      rate(:,165) = rate(:,165) * m(:)
      rate(:,166) = rate(:,166) * m(:)
      rate(:,167) = rate(:,167) * m(:)
      rate(:,168) = rate(:,168) * m(:)
      rate(:,169) = rate(:,169) * m(:)
      rate(:,170) = rate(:,170) * m(:)
      rate(:,171) = rate(:,171) * m(:)
      rate(:,172) = rate(:,172) * m(:)
      rate(:,173) = rate(:,173) * m(:)
      rate(:,174) = rate(:,174) * m(:)
      rate(:,175) = rate(:,175) * m(:)
      rate(:,176) = rate(:,176) * m(:)
      rate(:,177) = rate(:,177) * m(:)
      rate(:,178) = rate(:,178) * m(:)
      rate(:,179) = rate(:,179) * m(:)
      rate(:,180) = rate(:,180) * m(:)
      rate(:,181) = rate(:,181) * m(:)
      rate(:,183) = rate(:,183) * m(:)
      rate(:,184) = rate(:,184) * m(:)
      rate(:,186) = rate(:,186) * m(:)
      rate(:,187) = rate(:,187) * m(:)
      rate(:,188) = rate(:,188) * m(:)
      rate(:,189) = rate(:,189) * m(:)
      rate(:,190) = rate(:,190) * m(:)
      rate(:,191) = rate(:,191) * m(:)
      rate(:,192) = rate(:,192) * m(:)
      rate(:,193) = rate(:,193) * m(:)
      rate(:,194) = rate(:,194) * m(:)
      rate(:,195) = rate(:,195) * m(:)
      rate(:,196) = rate(:,196) * m(:)
      rate(:,197) = rate(:,197) * m(:)
      rate(:,198) = rate(:,198) * m(:)
      rate(:,200) = rate(:,200) * m(:)
      rate(:,201) = rate(:,201) * m(:)
      rate(:,202) = rate(:,202) * m(:)
      rate(:,203) = rate(:,203) * m(:)
      rate(:,204) = rate(:,204) * m(:)
      rate(:,205) = rate(:,205) * m(:)
      rate(:,206) = rate(:,206) * m(:)
      rate(:,207) = rate(:,207) * m(:)
      rate(:,208) = rate(:,208) * m(:)
      rate(:,209) = rate(:,209) * m(:)
      rate(:,210) = rate(:,210) * m(:)
      rate(:,211) = rate(:,211) * m(:)
      rate(:,212) = rate(:,212) * m(:)
      rate(:,213) = rate(:,213) * m(:)
      rate(:,214) = rate(:,214) * m(:)
      rate(:,215) = rate(:,215) * m(:)
      rate(:,216) = rate(:,216) * m(:)
      rate(:,217) = rate(:,217) * m(:)
      rate(:,218) = rate(:,218) * m(:)
      rate(:,219) = rate(:,219) * m(:)
      rate(:,220) = rate(:,220) * m(:)
      rate(:,221) = rate(:,221) * m(:)
      rate(:,222) = rate(:,222) * m(:)
      rate(:,223) = rate(:,223) * m(:)

      end subroutine adjrxt

      end module mo_adjrxt_mod

      module mo_phtadj_mod

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, plnplv )

      use chem_mods_mod, only : nfs, phtcnt

      implicit none

!--------------------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in)    :: inv(plnplv,nfs)
      real, intent(in)    :: m(plnplv)
      real, intent(inout) :: p_rate(plnplv,phtcnt)

!--------------------------------------------------------------------
!       ... Local variables
!--------------------------------------------------------------------
      real    ::  im(plnplv)

      im(:) = 1. / m(:)
      p_rate(:,  1) = p_rate(:,  1)  * inv(:, 3) * im(:)

      end subroutine phtadj

      end module mo_phtadj_mod

      module mo_rxt_mod

      private
      public :: rxt_mod

      contains

      subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )

      use chem_mods_mod, only : rxntot, hetcnt, grpcnt

      implicit none

!---------------------------------------------------------------------------
!       ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) ::  plnplv
      real, intent(inout) ::  rate(plnplv,rxntot)
      real, intent(inout) ::  het_rates(plnplv,hetcnt)
      real, intent(in)    ::  grp_ratios(plnplv,grpcnt)


      end subroutine rxt_mod

      end module mo_rxt_mod

      module mo_make_grp_vmr_mod

      private
      public :: mak_grp_vmr

      contains

      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )

      use mo_grid_mod,   only : plev, pcnstm1
      use chem_mods_mod, only : grpcnt

      implicit none

!----------------------------------------------------------------------------
!        ... Dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)    :: vmr(plonl,plev,pcnstm1)
      real, intent(in)    :: group_ratios(plonl,plev,grpcnt)
      real, intent(out)   :: group_vmrs(plonl,plev,grpcnt)

!----------------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------------

      end subroutine mak_grp_vmr

      end module mo_make_grp_vmr_mod

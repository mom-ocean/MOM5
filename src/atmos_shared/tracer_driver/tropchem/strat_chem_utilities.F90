module strat_chem_utilities_mod

use       mpp_io_mod, only : mpp_open, mpp_close, MPP_RDONLY
use          mpp_mod, only : mpp_pe, mpp_root_pe, stdout
use          fms_mod, only : file_exist, open_namelist_file, close_file, &
                             error_mesg, FATAL
use    constants_mod, only : PI, DEG_TO_RAD, AVOGNO, PSTD_MKS, SECONDS_PER_DAY
use time_manager_mod, only : time_type, get_date, days_in_month, days_in_year, &
                             set_date, increment_time, set_time, &
                             operator(>), operator(<), operator(-), operator(+)
!++lwh
use interpolator_mod, only : interpolate_type, interpolator_init, &
                             obtain_interpolator_time_slices, &
                             unset_interpolator_time_flag, &
                             interpolator, interpolator_end, &
                             CONSTANT, INTERP_WEIGHTED_P
use time_interp_mod, only  : time_interp_init, time_interp
use diag_manager_mod, only : get_base_time

!--lwh

implicit none
private

integer, parameter :: nlon_input=144, nlat_input=90, nlev_input=48, &
                      nspecies_age=8, nspecies_lbc=15, &
                      ntime_tropc=151, year_start_tropc=1950, nspecies_tropc=9
! real, parameter :: agefact1 = 1.5, &
!                    agefact2 = 1.25
real, parameter :: clweight(7) = (/ 3., 2., 3., 4., 1., 3., 1. /)

real :: tropc(ntime_tropc,nspecies_tropc)
!++lwh
type(time_type) :: tropc_Time(ntime_tropc), cfc_entry, cfc_offset
logical :: time_varying_cfc_lbc, negative_offset_cfc

! real :: dfdage(nlat_input,nlev_input,nspecies_age)
! real :: lat_input(nlat_input)
real, parameter :: days_per_year = 365.25, &
                   tfact = 1./(days_per_year*SECONDS_PER_DAY)
! integer :: jstart
real :: age_factor, dclydt_factor
!--lwh

type psc_type
   private
   real, dimension(:,:), pointer :: &
      tice=>NULL(), wh2so4=>NULL(), am=>NULL(), aw=>NULL(), &
      aliq=>NULL(), rmean=>NULL(), asat=>NULL(), rnat=>NULL(), rice=>NULL()
   real, dimension(:,:,:), pointer :: &
      cond=>NULL()
   real :: adrop, anat, aice
end type psc_type

!++lwh
type(interpolate_type), save  :: dfdage_interp
character(len=32), save  :: dfdage_filename = "dfdage3.dat.nc"
character(len=32), dimension(nspecies_age), save :: dfdage_name = &      ! kerr
      (/ "dfdage_cfc11  ", "dfdage_cfc12  ", "dfdage_cfc113 ", "dfdage_ccl4   ", &
         "dfdage_ch3cl  ", "dfdage_ch3ccl3", "dfdage_hcfc22 ", "dfdage_bry    " /)
!--lwh

! For extra H2O calculation
real, dimension(:), allocatable :: ch4_value
type(time_type), dimension(:), allocatable :: ch4_time
logical :: fixed_ch4_lbc_time = .false.
type(time_type) :: ch4_entry

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public strat_chem_utilities_init, strat_chem_dcly_dt, strat_chem_get_aerosol, &
       strat_chem_dcly_dt_time_vary, strat_chem_dcly_dt_endts, &
       strat_chem_get_h2so4, strat_chem_get_psc, strat_chem_destroy_psc, &
       strat_chem_get_gamma, strat_chem_get_hetrates, strat_chem_psc_sediment, &
       strat_chem_get_extra_h2o, &
       psc_type


!---- version number -----
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''

logical :: module_is_initialized=.false.

CONTAINS

subroutine strat_chem_utilities_init( lonb, latb, age_factor_in, dclydt_factor_in, &
                                      set_min_h2o_strat, ch4_filename, ch4_scale_factor, &
                                      fixed_ch4_lbc_time_in, ch4_entry_in, &
                                      cfc_lbc_filename, time_varying_cfc_lbc_in, cfc_lbc_dataset_entry )

   implicit none
! dummy arguments
   real,             intent(in), dimension(:,:) :: lonb,latb
   real,             intent(in)                 :: age_factor_in, dclydt_factor_in
   logical,          intent(in)                 :: set_min_h2o_strat
   character(len=*), intent(in)                 :: ch4_filename
   real,             intent(in)                 :: ch4_scale_factor
   logical,          intent(in)                 :: fixed_ch4_lbc_time_in
   type(time_type),  intent(in)                 :: ch4_entry_in
   character(len=*), intent(in)                 :: cfc_lbc_filename
   logical,          intent(in)                 :: time_varying_cfc_lbc_in
   integer,          intent(in), dimension(:)   :: cfc_lbc_dataset_entry
   
! local variables
   real :: chlb_dummy(nlat_input,nspecies_lbc), &
           ozb_dummy(nlon_input, nlat_input, 12)
   integer :: unit, nc, n, year, outunit
   type(time_type) :: Model_init_time
   
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... initialize time_interp
!-----------------------------------------------------------------------
   call time_interp_init

!-----------------------------------------------------------------------
!     ... set scale factors for age of air and dcly/dt
!-----------------------------------------------------------------------
   age_factor = age_factor_in
   dclydt_factor = dclydt_factor_in

!-----------------------------------------------------------------------
!     ... read in chemical lower boundary 
!-----------------------------------------------------------------------
   call mpp_open( unit, 'INPUT/' // TRIM(cfc_lbc_filename),action=MPP_RDONLY )
   outunit= stdout()
   if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'reading INPUT/' // TRIM(cfc_lbc_filename)
   do nc = 1,15                                           
     read(unit,'(6E13.6)') chlb_dummy(:,nc)
   end do
   read(unit,'(6E13.6)') ozb_dummy
   read(unit,'(6e13.6)') tropc
   call mpp_close(unit)

!++lwh
!---------------------------------------------------------------------
!    convert the time stamps of the tropc series to time_type variables.     
!---------------------------------------------------------------------
   time_varying_cfc_lbc = time_varying_cfc_lbc_in
   Model_init_time = get_base_time()
   if ( cfc_lbc_dataset_entry(1) == 1 .and. &
        cfc_lbc_dataset_entry(2) == 1 .and. &
        cfc_lbc_dataset_entry(3) == 1 .and. &
        cfc_lbc_dataset_entry(4) == 0 .and. &
        cfc_lbc_dataset_entry(5) == 0 .and. &
        cfc_lbc_dataset_entry(6) == 0 ) then
      cfc_entry = Model_init_time
   else
      cfc_entry = set_date( cfc_lbc_dataset_entry(1), &
                            cfc_lbc_dataset_entry(2), &
                            cfc_lbc_dataset_entry(3), &
                            cfc_lbc_dataset_entry(4), &
                            cfc_lbc_dataset_entry(5), &
                            cfc_lbc_dataset_entry(6) )
   end if         
   if (time_varying_cfc_lbc) then
      cfc_offset = cfc_entry - Model_init_time
      if (Model_init_time > cfc_entry) then
         negative_offset_cfc = .true.
      else
         negative_offset_cfc = .false.
      end if
   end if

   do n = 1,ntime_tropc
      year = year_start_tropc + (n-1)
      tropc_Time(n) = set_date(year,1,1,0,0,0)
   end do
!--lwh


!++lwh -- replace dfdage ASCII file with NetCDF via interpolator

!  read in data for Cly and Bry computation
!  call mpp_open( unit, 'INPUT/ageair_fms_90.dat', action=MPP_RDONLY )
!  if (mpp_pe() == mpp_root_pe()) &
!     write(stdout(),*) 'strat_chem_utilities_init: Reading from INPUT/ageair_fms_90.dat'
!  read(unit,'(6e13.6)') age_dummy
!  read(unit,'(6e13.6)') dfdage
!  call mpp_close(unit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  THIS CODE WILL NOT WORK CORRECTLY FOR CUBED SPHERE  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  jstart = 0   
!  do j = 1,nlat_input
!     lat_input(j) = ( -90. + (180./(nlat_input-1))*(j-1) ) * DEG_TO_RAD
!     if (lat_input(j) >= latb(1,1) .and. lat_input(j) <= latb(1,2)) jstart = j
!  end do
!  if (mpp_pe() == mpp_root_pe()) &
!     write(stdout(),*) 'strat_chem_utilities_init: jstart=',jstart,' on PE ',mpp_pe()
!

   call interpolator_init (dfdage_interp,  &
                           dfdage_filename, lonb, latb,  &
                           data_names=dfdage_name(:),   &
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/) )  

   if (set_min_h2o_strat) then
      call strat_chem_extra_h2o_init( ch4_filename, ch4_scale_factor, &
                                      fixed_ch4_lbc_time_in, ch4_entry_in )
   end if

   module_is_initialized = .true.


end subroutine strat_chem_utilities_init



!#####################################################################

subroutine strat_chem_dcly_dt_time_vary(Time)


type(time_type), intent(in) :: Time

      call obtain_interpolator_time_slices (dfdage_interp, Time)

end subroutine strat_chem_dcly_dt_time_vary



!#####################################################################

subroutine strat_chem_dcly_dt_endts              

      call unset_interpolator_time_flag (dfdage_interp)

end subroutine strat_chem_dcly_dt_endts



!#####################################################################

!++lwh
subroutine strat_chem_dcly_dt(Time, phalf, is, js, age, cly, bry, dclydt, dbrydt)
!--lwh

implicit none

! dummy arguments

type(time_type),        intent(in)  :: Time
!++lwh
integer,                intent(in)  :: is, js
real, dimension(:,:,:), intent(in)  :: phalf, age, cly, bry
!--lwh
real, dimension(:,:,:), intent(out) :: dclydt, dbrydt

! local variables

integer :: iyear, imon, iday, ihour, imin, isec
real :: dt1, factor, extra_seconds
integer :: it1,it2
integer :: ic, i, j, k, il, jl, kl
real :: clytot, brytot
!++lwh
! real, dimension(size(age,1),size(age,2),nspecies_age) :: dfdtau
real, dimension(size(age,1),size(age,2),size(age,3),nspecies_age) :: dfdage
type(time_type) :: cfc_Time, cfc_base_Time
!--lwh
real, dimension(size(age,1),size(age,2),nspecies_tropc) :: cfc
character(len=256) :: err_msg

call get_date( Time, iyear, imon, iday, ihour, imin, isec )

il = size(age,1)
jl = size(age,2)
kl = size(age,3)

!-----------------------------------------------------------------------
!     ... Compute multiplying factor for missing CFCs, and include factor
!         for conversion of rates to a per second rate. 
!-----------------------------------------------------------------------

! time0 = iyear + REAL(imon-1)/12.
! it1 = INT(time0-year_start_tropc) + 1
! it1 = min(max(it1,1),ntime_tropc-1)
! it2 = it1+1
! dt1 = time0 - (year_start_tropc-1) - it1

! sum1 = 0.
! sum2 = 0.
! do ic = 1,7
!    sum1 = sum1 + clweight(ic) * tropc(it1,ic)
!    sum2 = sum2 + clweight(ic) * tropc(it2,ic)
! end do
! factor = ((1-dt1)*tropc(it1,8) + dt1*tropc(it2,8))*tfact / (sum1*(1-dt1) + sum2*dt1)

! monthfrac = REAL(iday)/REAL(days_in_month(Time))
! imon2 = MOD(imon,12) + 1


!++lwh -- Read in dfdage from NetCDF file via interpolator
call interpolator (dfdage_interp, Time, phalf, dfdage, dfdage_name(1), is, js)
!--lwh

if (time_varying_cfc_lbc) then
   if (negative_offset_cfc) then
      cfc_base_Time = Time - cfc_offset
   else
      cfc_base_Time = Time + cfc_offset
   end  if
else
   cfc_base_Time = cfc_entry
end if

level_loop: &
do k = 1,kl


!++lwh -- Switch to reading from NetCDF file via interpolator (above)
! Copy dfdage to locally indexed variable
!  do j = 1,jl
!  do ic = 1,8
!     dfdtau(:,j,ic) = dfdage(j+jstart-1+js-1,k,ic)
!  end do
!  end do
!--lwh

!-----------------------------------------------------------------------
!     ... Compute CFCs at time t - age
!-----------------------------------------------------------------------

   do j = 1,jl
   do i = 1,il
!++lwh
! Time-interpolate tropospheric CFC concentrations
!     time0 = iyear + (imon+monthfrac-1)/12. - age(i,j,k)*age_factor
!     it1 = INT(time0-year_start_tropc) + 1
!     it1 = min(max(it1,1),ntime_tropc-1)
!     it2 = it1 + 1
!     dt1 = time0 - (year_start_tropc-1) - it1
!     cfc(i,j,:) = tropc(it1,:)*(1-dt1) + tropc(it2,:)*dt1

      extra_seconds = age(i,j,k)*age_factor / tfact
      
      cfc_Time = increment_time( cfc_base_Time, -NINT(extra_seconds), 0)
      if (cfc_Time < tropc_Time(1)) then
         cfc_Time = tropc_Time(1)
      else if (cfc_Time > tropc_Time(ntime_tropc)) then
         cfc_Time = tropc_Time(ntime_tropc)
      end if
      call time_interp( cfc_Time, tropc_Time(:), dt1, it1, it2, err_msg=err_msg )
      if(err_msg /= '') then
         call error_mesg('strat_chem_dcly_dt', trim(err_msg) , FATAL)
      endif

      cfc(i,j,:) = tropc(it1,:)*(1-dt1) + tropc(it2,:)*dt1
!--lwh

      factor = cfc(i,j,8) / SUM(cfc(i,j,1:7)*clweight(1:7))
      dclydt(i,j,k) = 0.
      do ic = 1,7
         dclydt(i,j,k) = dclydt(i,j,k) &
                       + factor * 1.e-12 * tfact * dclydt_factor &
!++lwh
!                      * dfdtau(i,j,ic) * clweight(ic) * cfc(i,j,ic)
                       * dfdage(i,j,k,ic) * clweight(ic) * cfc(i,j,ic)
!--lwh
      end do
      clytot = 1.e-12*cfc(i,j,8)
      if (cly(i,j,k) >= clytot) dclydt(i,j,k) = 0.
!++lwh
!     dbrydt(i,j,k) = 1.e-12 * tfact * dfdtau(i,j,8)*cfc(i,j,9)
      dbrydt(i,j,k) = 1.e-12 * tfact * dfdage(i,j,k,8)*cfc(i,j,9)
!--lwh
      brytot = 1.e-12*cfc(i,j,9)
      if (bry(i,j,k) >= brytot) dbrydt(i,j,k) = 0.
   end do
   end do
   
end do level_loop


end subroutine strat_chem_dcly_dt


! <SUBROUTINE NAME="strat_chem_get_aerosol">
!   <OVERVIEW>
!     Estimate stratospheric aerosol surface area based on aerosol extinction.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Estimate stratospheric aerosol surface area based on aerosol extinction.
!     This subroutine is called from tropchem_driver. Aerosol extinction in
!     band 4 (centered at 1um) is saved as a diagnostic tracer in swresf.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strat_chem_get_aerosol (extinct, aerosol)
!   </TEMPLATE>
!   <IN NAME="extinct" TYPE="real" DIM="(:,:,:)">
!     Aerosol extinction for band 4 centered at 1 um (in units of 1/m)
!   </IN>
!   <OUT NAME="aerosol" TYPE="real" DIM="(:,:,:)">
!     Volcanic aerosol surface area (cm2/cm3)
!   </OUT>


subroutine strat_chem_get_aerosol( extinct, aerosol )

implicit none

! dummy arguments

real, dimension(:,:,:), intent(in)  :: extinct
real, dimension(:,:,:), intent(out) :: aerosol

! local variables

real, dimension(size(extinct,1),size(extinct,2)) :: extinct_local
integer :: k

do k = 1,size(extinct,3)

   extinct_local(:,:) = extinct(:,:,k) * 1.e3 ! convert to 1/km
   where( extinct_local(:,:) <= 0. )
      aerosol(:,:,k) = 0.
   elsewhere( extinct_local(:,:) <= 4.e-3 )
      aerosol(:,:,k) = 4.25e-6 * extinct_local(:,:)**0.68
   elsewhere( extinct_local(:,:) <= 2.e-2 )
      aerosol(:,:,k) = 1.223e-5 * extinct_local(:,:)**0.875
   elsewhere
      aerosol(:,:,k) = 2.e-5 * extinct_local(:,:)
   end where

end do

end subroutine strat_chem_get_aerosol


! <SUBROUTINE NAME="strat_chem_get_h2so4">
!   <OVERVIEW>
!     Set stratospheric H2SO4
!   </OVERVIEW>
!   <DESCRIPTION>
!     Set stratospheric H2SO4 as an analytical function of age, peaking at 0.5 ppbv
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strat_chem_get_h2so4( press, age, h2so4 )
!   </TEMPLATE>
!   <IN NAME="press" TYPE="real" DIM="(:,:)">
!     Pressure on model full levels (Pa)
!   </IN>
!   <IN NAME="age" TYPE="real" DIM="(:,:)">
!     Age-of-air tracer (yrs)
!   </IN>
!   <OUT NAME="h2so4" TYPE="real" DIM="(:,:)">
!     Sulfuric acid (H2SO4) VMR (mol/mol)
!   </OUT>
subroutine strat_chem_get_h2so4(press, age, h2so4)

implicit none

! Dummy arguments

real, dimension(:,:), intent(in)  :: press, age
real, dimension(:,:), intent(out) :: h2so4

! Local variables

integer :: i, k, il, kl
real :: x1, x2


il = size(press,1)
kl = size(press,2)

do k = 1,kl
   do i=1,il

        x1 = 3.*(log10(press(i,k)) - 3.5)**2
        x2 = 3./(10.*age(i,k) + 1.) * (age(i,k) - 1.)**2
        h2so4(i,k) = (0.01 + 0.49*exp(-x1-x2))*1.e-9

   end do
end do

end subroutine strat_chem_get_h2so4



! <SUBROUTINE NAME="strat_chem_get_psc">
!   <OVERVIEW>
!     Set up stratospheric PSCs
!   </OVERVIEW>
!   <DESCRIPTION>
!     Set up stratospheric PSCs
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strat_chem_get_psc( temp, press, hno3, h2o, h2so4, strat_aerosol, psc, psv_vmr_out )
!   </TEMPLATE>
!   <IN NAME="temp" TYPE="real" DIM="(:,:)">
!     Temperature (K)
!   </IN>
!   <IN NAME="press" TYPE="real" DIM="(:,:)">
!     Pressure on model full levels (Pa)
!   </IN>
!   <INOUT NAME="hno3" TYPE="real" DIM="(:,:)">
!     Nitric acid (HNO3) VMR (mol/mol)
!   </INOUT>
!   <INOUT NAME="h2o" TYPE="real" DIM="(:,:)">
!     Water (H2O) VMR (mol/mol)
!   </INOUT>
!   <IN NAME="h2so4" TYPE="real" DIM="(:,:)">
!     Sulfuric acid (H2SO4) VMR (mol/mol)
!   </IN>
!   <IN NAME="strat_aerosol" TYPE="real" DIM="(:,:)">
!     Stratospheric aerosol (liquid)
!   </IN>
!   <OUT NAME="psc" TYPE="psc_type">
!     PSC properties
!   </OUT>
!   <OUT NAME="psc_vmr_out" TYPE="real" DIM="(:,:)" OPTIONAL>
!     PSC VMR (mol/mol)
!   </OUT>

subroutine strat_chem_get_psc( temp, press, hno3, h2o, h2so4, strat_aerosol, psc, psc_vmr_out )

implicit none

! Dummy arguments

real, dimension(:,:), intent(in)    :: temp, press, h2so4, strat_aerosol
real, dimension(:,:), intent(inout) :: hno3, h2o
type(psc_type),       intent(out)   :: psc
real, dimension(:,:,:), intent(out), optional :: psc_vmr_out

! Local variables
integer :: il, kl, i, k
real :: rho(size(temp,1)), &
        DENS(size(temp,1)), &
        VHET(size(temp,1),4), &
        WF(size(temp,1),2)
real :: PH2O, &                    ! H2O partial pressure (atm)
        TSAT,ABT,AMT,BMT,PX,C2, &
        A1,A2,A3,A4,C3,P0H2O,Y1,Y2,T1,RMODE,RMODESAT,SANAT
!----------------------------------------------------------------------------
!
! AVGDR IS THE RECIPROCAL OF THE AVOGADRO CONSTANT
! RR IS THE GAS CONSTANT
!
!----------------------------------------------------------------------------
real, parameter :: AVGDR = 1/AVOGNO, &
                   RR = 8.3144
!----------------------------------------------------------------------------
!
!   ADROP, ANAT AND AICE ARE THE (FIXED) NUMBER OF DROPS OR PARTICLES PER CC
!   IN THE POLAR STRATOSPHERC CLOUDS. A LOG NORMAL SIZE DISTRIBUTION IS ASSUMED
!   WITH SIGMA = 1.8 GIVING LOG(SIGMA)**2 = 0.3455
!
!----------------------------------------------------------------------------
real, parameter :: SIGMAL=0.34549316


real, parameter ::  boltz = 1.38044e-16      ! erg/K


il = size(temp,1)
kl = size(temp,2)

allocate( psc%cond(il,kl,3), &
          psc%tice(il,kl), psc%wh2so4(il,kl), psc%am(il,kl), &
          psc%aw(il,kl), psc%aliq(il,kl), psc%rmean(il,kl), &
          psc%asat(il,kl), psc%rnat(il,kl), psc%rice(il,kl) )

!
!----------------------------------------------------------------------------
!
! COMPUTE SOLID CONDENSED MATTER: COND(IL,1:3) -- SAT, NAT, ICE
!
!----------------------------------------------------------------------------
psc%cond(:,:,1:3) = 0.
psc%adrop = 10.
psc%anat  = 10.
psc%aice  = 10.

do k = 1,kl
   do i=1,il
      rho(i) = 10. * press(i,k) / (boltz * temp(i,k))  ! molec/cm3
      PH2O = h2o(i,k)*press(i,k)/PSTD_MKS
      TSAT = 3236.0 / (11.502 - LOG10(PH2O*760.0))
      IF(temp(i,k) < TSAT) psc%cond(i,k,1) = h2so4(i,k)

      ABT = 39.1104 - 11397.0/temp(i,k) + 9.179E-3*temp(i,k)
      AMT = -2.7836 - 8.8E-4*temp(i,k)
      BMT = -2.1249 + LOG10(PH2O*PSTD_MKS)
      PX = AMT*BMT + ABT
      C2 = hno3(i,k) - 100./press(i,k) * 10.**PX
      IF(C2 > 0.) psc%cond(i,k,2) = C2

      psc%tice(i,k) = 2668.70/(10.4310 - LOG10(760.0*PH2O))
      A1 = 7.5502 - 2668.7/temp(i,k)
      C3 = h2o(i,k) - PSTD_MKS/press(i,k) * 10.**A1
!     C4 = PSTD_MKS/press(i,k) * 10.**A1
      IF(C3 > 0.) psc%cond(i,k,3) = C3

      hno3(i,k) = hno3(i,k) - psc%cond(i,k,2)
!     DYY(i,k,15) = DYY(i,k,15) - psc%cond(i,k,2)
!     CLIMIT = rho(i)*1. E-15
!     DYY(i,k,15) = MAX(DYY(i,k,15),CLIMIT)
      h2o(i,k) = h2o(i,k) - psc%cond(i,k,3)
   end do
!-------------------------------------------------------------------------
!
!  COMPUTE WEIGHT % H2SO4. FROM TABAZADEH ET AL., GRL, 24, 1931-1934, 1997.
!  TABLE A1 OF SHIA ET AL. JGR, 106, 24,529-24,274, 2001
!
!-------------------------------------------------------------------------
   do i = 1,IL
! Saturation water vapor pressure (in mbar)
      P0H2O = EXP(18.452406985 - 3505.1578807/temp(i,k) - &
         330918.55082/temp(i,k)**2 + 12725068.262/temp(i,k)**3)
! Water activity
      psc%aw(i,k) = h2o(i,k)*press(i,k)*0.01/P0H2O
!++lwh
      psc%aw(i,k) = MAX( MIN(psc%aw(i,k), 1.), 0.01 )
!--lwh
      IF (psc%aw(i,k) <= 0.05) THEN
         Y1 = 12.37208932*psc%aw(i,k)**(-0.16125516114) - 30.490657554*psc%aw(i,k) - 2.1133114241
         Y2 = 13.455394705*psc%aw(i,k)**(-0.1921312255) - 34.285174607*psc%aw(i,k) - 1.7620073078
      ELSEIF(psc%aw(i,k) < 0.85) THEN
         Y1 = 11.820654354*psc%aw(i,k)**(-0.20786404244) - 4.807306373*psc%aw(i,k) -  5.1727540348
         Y2 = 12.891938068*psc%aw(i,k)**(-0.23233847708) - 6.4261237757*psc%aw(i,k) - 4.9005471319
      ELSE
         Y1 = -180.06541028*psc%aw(i,k)**(-0.38601102592) - 93.317846778*psc%aw(i,k) + 273.88132245
         Y2 = -176.95814097*psc%aw(i,k)**(-0.36257048154) - 90.469744201*psc%aw(i,k) + 267.45509988
      ENDIF
! Sulfuric acid molality
!++lwh
!     psc%am(i,k) = Y1 + (temp(i,k) - 190.)*(Y2 - Y1)/70.
      psc%am(i,k) = (temp(i,k) - 190.)/70.
      psc%am(i,k) = MAX( MIN(psc%am(i,k), 1.), 0. )
      psc%am(i,k) = Y1 + psc%am(i,k)*(Y2 - Y1)
!--lwh
! Sulfuric acid weight percent
      psc%wh2so4(i,k) = 9800.*psc%am(i,k)/(98.*psc%am(i,k) + 1000.)
      WF(i,1) = 0.01*psc%wh2so4(i,k)
      WF(i,2) = 0.0
   end do
!---------------------------------------------------------------------------
!
!  COMPUTE DENSITY OF BINARY AEROSOL
!
!---------------------------------------------------------------------------
   CALL DENSITY(WF,temp(:,k),DENS)
!---------------------------------------------------------------------------
!
!  COMPUTE VOLUME OF BINARY AEROSOL/SAT/NAT/ICE
!  1.6, 1.35 and 0.928  are the densities of SAT, NAT and ICE
!
!---------------------------------------------------------------------------
   do i=1,IL
     T1 = h2so4(i,k)*press(i,k)/(rho(i)*temp(i,k)*RR)
     VHET(i,1) = T1*98.076E-6/(WF(i,1)*DENS(I))
     VHET(i,2) = psc%cond(i,k,1)*rho(i)*170.1*AVGDR/1.6
     VHET(i,3) = psc%cond(i,k,2)*rho(i)*117.1*AVGDR/1.35
     VHET(i,4) = psc%cond(i,k,3)*rho(i)*18.02*AVGDR/0.928
   end do
!---------------------------------------------------------------------------
!
!  COMPUTE PARTICLE PARAMETERS FROM WHICH THE HETEROGENEOUS REACTION RATES
!  ARE DETERMINED; ASSUME SURFACE AREA FROM SAT IS LIMITED BY NAT AMOUNT
!
!---------------------------------------------------------------------------
   A1 = EXP(-4.5*SIGMAL)
   A2 = EXP(0.5*SIGMAL)
   A3 = EXP(2.0*SIGMAL)
   A4 = 1.33333333*PI*psc%adrop
   do i=1,IL
     RMODE = (VHET(i,1)*A1/A4)**0.33333333
     psc%rmean(i,k) = MAX(RMODE*A2, 1.e-12)
!    psc%aliq(i,k) = 3.0*A4*(RMODE**2)*A3
     psc%aliq(i,k) = strat_aerosol(i,k)
     RMODESAT = (VHET(i,2)*A1/A4)**0.33333333
     psc%asat(i,k) = 3.0*A4*(RMODESAT**2)*A3
     psc%rnat(i,k) = (VHET(i,3)/(1.33333333*PI*psc%anat))**0.33333333
     psc%rice(i,k) = (VHET(i,4)/(1.33333333*PI*psc%aice))**0.33333333
     SANAT = 4.*PI * psc%anat * psc%rnat(i,k)**2
     psc%asat(i,k) = MAX(psc%asat(i,k) - SANAT, 0.)
!    psc%aliq(i,k) = MAX(psc%aliq(i,k) - SANAT, 0.)
   end do
end do

if (present(psc_vmr_out)) then
   psc_vmr_out(:,:,:) = psc%cond(:,:,:)
end if

end subroutine strat_chem_get_psc


subroutine strat_chem_destroy_psc( psc )

implicit none

! Dummy arguments

type(psc_type),   intent(inout)   :: psc


deallocate( psc%cond, psc%tice, psc%wh2so4, psc%am, &
            psc%aw, psc%aliq, psc%rmean, &
            psc%asat, psc%rnat, psc%rice )

end subroutine strat_chem_destroy_psc


subroutine strat_chem_get_gamma(temp, press, rho, hcl, clono2, psc, k, gamma)
!-------------------------------------------------------------------------
!
! SUBROUTINE TO CALCULATE REACTION PROBABILITIES ON SULPHATE AEROSOL
! BASED ON JPL'03 RECOMMENDATION
!
!-------------------------------------------------------------------------
IMPLICIT NONE


! dummy arguments
REAL, dimension(:),   intent(in)  :: temp, press, rho, hcl, clono2
type(psc_type),       intent(in)  :: psc
integer,              intent(in)  :: k
REAL, dimension(:,:), intent(out) :: gamma

! local variables

integer :: i,il
REAL, dimension(size(temp,1)) :: AMH2SO4, XH2SO4, VISC, AACID, TEMP2, RTT
REAL :: T2,Z1,Z2,Z3,RHOX,AA,X,T1,T3,AKH,AKH2O,AKHYDR,DIFF,SCLONO2, &
        CCLONO2,GAMMAB1,HHCL,AKHCL,Q1,RQ,A1,FCLONO2,GAMMARXN,      &
        GAMMABHCL,GAMMAS,FHCL,GAMMASP,GAMMABHCLP,GAMMAB,GCLONO2,   &
        SHOCL,HHOCL,FHOCL,WT,AK0,AK1,AK2,T0,HCLONO2,   &
        AMHCL,AKHOCL,CHOCL

!  The parameterisations used here break down below about 185K, so the
!  temperature is here limited to 185K and above (TEMP2).

   il = size(temp,1)

!-------------------------------------------------------------------------
!
! CALCULATE H2SO4 MOLARITY (AMH2SO4), MOLE FRACTION (XH2SO4),
!  VISCOSITY (VISC) AND ACID ACTIVITY (AACID)
!  TABLE A2 OF SHIA ET AL. JGR, 106, 24,529-24,274, 2001.
!
!-------------------------------------------------------------------------
   TEMP2(:) = MAX(temp(:),185.)
   RTT(:) = SQRT(temp(:))

long_loop: &
   do i = 1,il
      T2 = TEMP2(I)**2
      Z1 = 0.12364 - 5.6E-7*T2
      Z2 = -0.02954 + 1.814E-7*T2
      Z3 = 2.343E-3 - 1.487E-6*TEMP2(I) - 1.324E-8*T2
      RHOX = 1.0 + Z1*psc%am(i,k) + Z2*psc%am(i,k)**1.5 + Z3*psc%am(i,k)**2
      AMH2SO4(I) = RHOX*psc%wh2so4(i,k)/9.8
      XH2SO4(I) = psc%wh2so4(i,k)/(psc%wh2so4(i,k) + (100.0 - psc%wh2so4(i,k))*98.0/18.0)
      AA = 169.5 + 5.18*psc%wh2so4(i,k) - 0.0825*psc%wh2so4(i,k)**2 + 3.27E-3*psc%wh2so4(i,k)**3
      T0 = 144.11 + 0.166*psc%wh2so4(i,k) - 0.015*psc%wh2so4(i,k)**2 + 2.18E-4*psc%wh2so4(i,k)**3
      X = TEMP2(I)**(-1.43)
      VISC(I) = AA*X*EXP(448.0/(TEMP2(I) - T0))
      T1 = 60.51 - 0.095*psc%wh2so4(i,k) + 0.0077*psc%wh2so4(i,k)**2 - 1.61E-5*psc%wh2so4(i,k)**3
      T2 =  (-805.89 + 253.05*psc%wh2so4(i,k)**0.076)/RTT(I)
      T3 =   (1.76 + 2.52E-4*psc%wh2so4(i,k)**2)*RTT(I)
      AACID(I) = EXP(T1 + T2 - T3)
!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR CLONO2 + H2O AND CLONO2 + HCL AND
!  HENRY'S LAW COEFFICIENTS.
!  TABLE A3 OF SHIA ET AL. JGR, 106, 24,529-24,274, 2001.
!
!  The following formulation for the water activity is from Shi et al.,
!  but the differences between their parameterisation and that calculated
!  using the actual model H2O is not large.
!
!      awx = exp((-69.775*xh2so4(il) - 18253.7*xh2so4(il)**2 +   &
!           31072.2*xh2so4(il)**3 - 25668.8*xh2so4(il)**4)*      &
!           (1.0/temp(il) - 26.9033/(temp(il)**2)))
!      AKHYDR = AWx*(AKH2O + AKH*AACID(IL))
!
!-------------------------------------------------------------------------
      AKH = 1.22E12*EXP(-6200.0/TEMP2(I))
      AKH2O = 1.95E10*EXP(-2800.0/TEMP2(I))
      AKHYDR = psc%aw(i,k)*(AKH2O + AKH*AACID(I))
      DIFF = 5.0E-8*TEMP2(I)/VISC(I)
      SCLONO2 = 0.306 + 24.0/TEMP2(I)
      HCLONO2 = 1.6E-6*EXP(4710.0/TEMP2(I) - SCLONO2*AMH2SO4(I))
      CCLONO2 = 1474.*TEMP2(I)**0.5
      GAMMAB1 = (4.0*HCLONO2*0.082*TEMP2(I)/CCLONO2) * (DIFF*AKHYDR)**0.5

      HHCL = (0.094 - 0.61*XH2SO4(I) + 1.2*XH2SO4(I)**2)* &
         EXP(-8.68 + (8515. - 10718.*XH2SO4(I)**0.7)/TEMP2(I))
      AMHCL = HHCL*hcl(i)*press(I)/PSTD_MKS
      AKHCL = 7.9E11*AACID(I)*DIFF*AMHCL
      Q1 = (DIFF/(AKHYDR + AKHCL))**0.5
      RQ = psc%rmean(i,k)/Q1
      A1 = RQ + 0.312*RQ**2
      FCLONO2 = A1/(3.0 + A1)
      GAMMARXN = FCLONO2*GAMMAB1*(1.0 + AKHCL/AKHYDR)**0.5
      GAMMABHCL = GAMMARXN*AKHCL/(AKHCL + AKHYDR)
      GAMMAS = 66.12* EXP(-1374./TEMP2(I))*HCLONO2*AMHCL
      IF( hcl(i) /= 0. ) THEN
       FHCL = 1.0/(1.0 + 0.612*(GAMMAS + GAMMABHCL)*clono2(i)/hcl(i))
      ELSE
       FHCL = 0.0
      ENDIF
      GAMMASP = FHCL*GAMMAS
      GAMMABHCLP = FHCL*GAMMABHCL
      GAMMAB = GAMMABHCLP + GAMMARXN*AKHYDR/(AKHCL + AKHYDR)
      GCLONO2 = 1.0/(1.0 + 1.0/(GAMMASP + GAMMAB))
      gamma(i,5) = GCLONO2*(GAMMASP + GAMMABHCLP) / (GAMMASP + GAMMAB)
      gamma(i,4) = GCLONO2 - gamma(i,5)

!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR HOCL + HCL AND HENRY'S LAW COEFFICIENTS.
!  TABLE A4 OF SHIA ET AL. JGR, 106, 24,529-24,274, 2001.
!
!-------------------------------------------------------------------------
      SHOCL = 0.0776 + 59.18/TEMP2(I)
      HHOCL = 1.91E-6*EXP(5862.4/TEMP2(I) - SHOCL*AMH2SO4(I))
      DIFF = 6.4E-8 *TEMP2(I)/VISC(I)
      AKHOCL = 1.25E9*AACID(I)*DIFF*AMHCL
      CHOCL = 2009.*RTT(I)
      Q1 = (DIFF/AKHOCL)**0.5
      RQ = psc%rmean(i,k)/Q1
      A1 = RQ + 0.312*RQ**2
      FHOCL = A1/(3.0 + A1)
      GAMMARXN = (FHOCL*4.0*HHOCL*0.082*TEMP2(I)/CHOCL) * (DIFF*AKHOCL)**0.5
      gamma(i,1) = 1.0/(1.0 + 1.0/(GAMMARXN*FHCL))

!-------------------------------------------------------------------------
!
!  CALCULATE REACTION PROBABILITES FOR N2O5 + H2O
!  ROBINSON ET AL. JGR, 102, 3583-3601, 1997.
!
!-------------------------------------------------------------------------
      WT = MIN( psc%wh2so4(i,k), 80. )
      AK0 = -25.5265 - 0.133188*WT + 0.00930846*WT**2 - 9.0194E-5*WT**3
      AK1 = 9283.76 + 115.345*WT - 5.19258*WT**2 + 0.0483464*WT**3
      AK2 = -851801. - 22191.2*WT + 766.916*WT**2 - 6.85427*WT**3
      gamma(i,3) = EXP(AK0 + AK1/TEMP2(I) + AK2/(TEMP2(I)**2))
   end do long_loop

!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR
!       N2O5 + HCL
!       HOBR + HCL
!       HOCL + HBR
!  NO RECOMMENDATION IN JPL '03, ASSUMED ZERO
!
!-------------------------------------------------------------------------
   gamma(:,2) = 0.
   gamma(:,6) = 0.
   gamma(:,7) = 0.

!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR HOBR + HBR
!  ABBATT, JGR, 100, 14009-14017, 1995.
!
!-------------------------------------------------------------------------
   gamma(:,8) = 0.25

!-------------------------------------------------------------------------
!
!  REACTION PROBABILITES FOR BRONO2 + H2O
!  USE JPL '03 RECOMMENDATION (HANSON PERS. COMM.)
!
!-------------------------------------------------------------------------
   gamma(:,9) = 1.0/(1.2422 + 1.0/(0.114 + EXP(29.24 - 0.396*psc%wh2so4(:,k))))

end subroutine strat_chem_get_gamma

subroutine strat_chem_get_hetrates( temp, hcl, hbr, h2o, rho, psc, gamma, k, tstep, rates )

implicit none
!------------------------------------------------------------------------
!
!  This subroutine computes the equivalent 2nd order reaction rates for
!  the heterogeneous reactions on aerosol, nat and ice.
!
!------------------------------------------------------------------------


! dummy arguments

INTEGER,              intent(in)  :: k
REAL, dimension(:),   intent(in)  :: temp, &   ! temperature (K)
                                     hcl, &    ! HCl volume mixing ratio (VMR)
                                     hbr, &    ! HBr VMR
                                     h2o, &    ! water vapor VMR
                                     rho       ! atmospheric density (molec/cm3)
REAL, dimension(:,:), intent(in)  :: gamma     ! reaction probabilities
type(psc_type),       intent(in)  :: psc       ! polar stratospheric clouds (PSCs)
real,                 intent(in) ::  tstep     ! timestep (sec)
REAL, dimension(:,:), intent(out) :: rates     ! heterogeneous reaction rate constants (molec^-1 cm^3 s^-1)

! local variables

integer, parameter :: nhet_data = 9
INTEGER :: i, ic, IL, NHET
REAL, dimension(size(temp)) :: CHEMC, &        ! chemical species number density (molec/cm3)
                               G2NAT, &        ! reaction probabilities on NAT
                               DELT, &         ! temperature excess over ice saturation (K)
                               SICE, &         ! 
                               RTT, &
                               VEL1, &
                               VEL2
REAL :: VCONST, ANUM, ADEN, MW2, area, afac1, adsorb_frac

real, dimension(nhet_data), parameter :: &
         AMW = (/ 52.45, 108.00, 108.00, 97.45, 97.45, 96.91, 52.45, 96.91, 141.91 /), &
!------------------------------------------------------------------------
!  AMW = MOLECULAR WEIGHT OF GAS PHASE SPECIES
!------------------------------------------------------------------------
         GNAT = (/0.1, 3.0E-3, 4.0E-4, 4.0E-3, 0.2, 0.0, 0.0, 0.0, 0.0 /), &
         GICE = (/0.2, 0.03, 0.02, 0.3, 0.3, 0.3, 0.03, 0.1, 0.3 /)
integer, parameter :: &
         HET_H2O=0, HET_HCL=1, HET_HBR=2
integer, dimension(nhet_data), parameter :: &
         INN = (/ HET_HCL,HET_HCL,HET_H2O,HET_H2O,HET_HCL,HET_HCL,HET_HBR,HET_HBR,HET_H2O /)
real, parameter :: small_conc = 1.e-20, &
                   mw_hcl = 36.46, &
                   mw_hbr = 80.91, &
                   mw_h2o = 18.01, &
                   adsorb_sites = 1.e15 ! adsorption sites per cm^2
!------------------------------------------------------------------------
!  INN = INDEX NUMBER OF LIQUID/SOLID PHASE SPECIES (0= H2O)
!-------------------------------------------------------------------------
!     REACTION 70 HOCL + HCL --> H2O + CL2 (HETEROGENEOUS)
!     REACTION 71 N2O5 + HCL --> HNO3 + CLNO2 (HETEROGENEOUS)
!     REACTION 72 N2O5+H2O --> 2HNO3 (HETEROGENEOUS)
!     REACTION 73 CLONO2+H2O --> HOCL+HNO3 (HETEROGENEOUS)
!     REACTION 74 CLONO2+HCL --> CL2+HNO3 (HETEROGENEOUS)
!     REACTION 75 HOBR + HCL --> BRCL + H2O (HETEROGENEOUS)
!     REACTION 76 HOCL + HBR --> BRCL + H2O (HETEROGENEOUS)
!     REACTION 77 HOBR + HBR --> 2BR + H2O (HETEROGENEOUS)
!     REACTION 78 BRONO2 + H2O --> HOBR + HNO3 (HETEROGENEOUS)
!     (THE FIRST MOLECULE ON THE LEFT HAND SIDE IS IN THE GAS PHASE,
!     THE SECOND MOLECULE IS IN THE LIQUID/SOLID PHASE)
!-------------------------------------------------------------------------
   IL = size(temp)
   NHET = size(gamma,2)

   RTT(:) = sqrt(temp(:))
   DELT(:) = temp(:) - psc%tice(:,k)
   SICE(:) = 10**(2668.70*(1.0/temp(:) - 1.0/psc%tice(:,k)))
   SICE(:) = MIN( SICE(:), 3. )

   do ic = 1,NHET
      VCONST = sqrt(8.0*8.3144E7/(PI*AMW(ic)))
      select case (ic)
         case(4)
            G2NAT(:) = EXP(-9.03 + 2.81*SICE(:))
         case(5)
            G2NAT(:) = 1.0/(4.3478 + 1.4241*EXP(0.518*DELT(:)))
         case default
            G2NAT(:) = GNAT(ic)
      end select
      select case (INN(ic))
         case(HET_H2O)
            CHEMC(:) = h2o(:)*rho(:)
            mw2 = mw_h2o
         case(HET_HCL)
            CHEMC(:) = hcl(:)*rho(:)
            mw2 = mw_hcl
         case(HET_HBR)
            CHEMC(:) = hbr(:)*rho(:)
            mw2 = mw_hbr
      end select
      VEL1(:) = VCONST * RTT(:)
      VEL2(:) = sqrt(8.*8.3144E7/(pi*mw2)) * RTT(:)
      rates(:,ic) =  0.

      select case (INN(ic))
         case(HET_H2O)

            do i=1,IL
               ADEN = CHEMC(i)
               if (ADEN > small_conc) then
!-------------------------------------------------------------------------
!    Reactions on NAT
!-------------------------------------------------------------------------
                  area = 4.*PI * psc%anat * psc%rnat(i,k)**2
                  ANUM = 0.25*VEL1(i)*G2NAT(i) * area
                  IF(ANUM > 0.) rates(i,ic) =  ANUM/ADEN
!-------------------------------------------------------------------------
!    Reactions on ICE
!------------------------------------------------------------------------
                  area = 4.*PI * psc%aice * psc%rice(i,k)**2
                  ANUM = 0.25*VEL1(i)*GICE(ic) * area
                  IF(ANUM > 0.) rates(i,ic) = rates(i,ic) + ANUM/ADEN
!-------------------------------------------------------------------------
!    Reactions on LIQUID AEROSOL
!    aliq is the liquid surface area. const is sqrt(8R/(pi*mw)) for the
!    gaseous phase species, with the mean molecular speed equal to
!        const*sqrt(temp)
!-------------------------------------------------------------------------
                  ANUM = 0.25*VEL1(i)*gamma(i,ic) * psc%aliq(i,k)
                  if (ANUM > 0.) rates(i,ic) = rates(i,ic) + ANUM / ADEN
               end if
            end do

         case(HET_HCL,HET_HBR)

            do i=1,IL
!-------------------------------------------------------------------------
!    Reactions on NAT
!-------------------------------------------------------------------------
               area = 4.*PI * psc%anat * psc%rnat(i,k)**2
               ANUM = 0.25*VEL1(i)*G2NAT(i)
               if (ANUM>0. .and. area>0.) then
                  afac1 = 0.25*VEL2(i)*tstep*area
!                 adsorb_frac = MIN(0.5*afac1, 1.)
                  adsorb_frac = 1. - (1./afac1)*(1.-exp(-afac1))
                  rates(i,ic) =  ANUM*adsorb_frac/adsorb_sites
               end if
!-------------------------------------------------------------------------
!    Reactions on ICE
!------------------------------------------------------------------------
               area = 4.*PI * psc%aice * psc%rice(i,k)**2
               ANUM = 0.25*VEL1(i)*GICE(ic)
               if(ANUM>0. .and. area>0.) then
                  afac1 = 0.25*VEL2(i)*tstep*area
!                 adsorb_frac = MIN(0.5*afac1, 1.)
                  adsorb_frac = 1. - (1./afac1)*(1.-exp(-afac1))
                  rates(i,ic) = rates(i,ic) + ANUM*adsorb_frac/adsorb_sites
               end if
!-------------------------------------------------------------------------
!    Reactions on LIQUID AEROSOL
!    aliq is the liquid surface area. const is sqrt(8R/(pi*mw)) for the
!    gaseous phase species, with the mean molecular speed equal to
!        const*sqrt(temp)
!-------------------------------------------------------------------------
               area = psc%aliq(i,k)
               ANUM = 0.25*VEL1(i)*gamma(i,ic)
               if(ANUM>0. .and. area>0.) then
                  afac1 = 0.25*VEL2(i)*tstep*area
!                 adsorb_frac = MIN(0.5*afac1, 1.)
                  adsorb_frac = 1. - (1./afac1)*(1.-exp(-afac1))
                  rates(i,ic) = rates(i,ic) + ANUM*adsorb_frac/adsorb_sites
               end if
            end do

      end select

   end do

end subroutine strat_chem_get_hetrates


subroutine DENSITY(WF,T,DENS)
implicit none
!
!    Density of ternary solution in g cm-3
!

! dummy arguments
REAL, intent(in)  :: WF(:,:),T(:)
REAL, intent(out) :: DENS(:)

! local variables
INTEGER :: I
REAL :: W,WH,T2,V1,A1,A2,VS,VN,VMCAL

real, parameter :: &
       X(22) = (/ 2.393284E-02,-4.359335E-05,7.961181E-08,0.0,-0.198716351, &
                  1.39564574E-03,-2.020633E-06,0.51684706,-3.0539E-03,      &
                  4.505475E-06,-0.30119511,1.840408E-03,-2.7221253742E-06,  &
                 -0.11331674116,8.47763E-04,-1.22336185E-06,0.3455282,      &
                 -2.2111E-03,3.503768245E-06,-0.2315332,1.60074E-03,        &
                 -2.5827835E-06 /), &
       AMR(3) = (/ 0.05550622,0.01019576,0.01586899 /)

   DO I = 1,SIZE(T)
      W = WF(I,1) + WF(I,2)
      WH = 1.0 - W
      T2 = T(I)**2
      V1 = X(1) + X(2)*T(I) + X(3)*T2 + X(4)*T2*T(I)
      A1 = X(8) + X(9)*T(I) + X(10)*T2
      A2 = X(11) + X(12)*T(I) + X(13)*T2
      VS = X(5) + X(6)*T(I) + X(7)*T2 + A1*W + A2*W**2
      A1 = X(17) + X(18)*T(I) + X(19)*T2
      A2 = X(20) + X(21)*T(I) + X(22)*T2
      VN = X(14) + X(15)*T(I) + X(16)*T2 + A1*W + A2*W**2
      VMCAL = WH*V1*AMR(1) + VS*WF(I,1)*AMR(2) + VN*WF(I,2)*AMR(3)
      DENS(I) = 1.0E-3/VMCAL
   END DO

end subroutine DENSITY


subroutine strat_chem_psc_sediment( psc, pfull, dt, dpsc )      

implicit none
!------------------------------------------------------------------------
!
!  This subroutine calculates sedimentation rates of Type I and Type II
!  particles and vertically advects model NAT and ice
!
!------------------------------------------------------------------------


! dummy arguments
!
!  CALCULATES SEDIMENTATION RATES OF TYPE I AND TYPE 2 PARTICLES               
!  AND VERTICALLY ADVECTS MODEL NAT AND ICE                                    
!
REAL, dimension(:,:,:),   intent(in)  :: pfull
REAL, dimension(:,:,:,:), intent(in)  :: psc
REAL,                     intent(in)  :: dt
REAL, dimension(:,:,:,:), intent(out) :: dpsc

! local variables

real, dimension(SIZE(pfull,3)) :: ANAT, AICE, SNATS, SICES, F1, F2, ANAT2, AICE2
real :: PNAT,PICE,PNAT2,PICE2
real :: ANATMAX,AICEMAX
integer :: i, j, k, il, jl, kl
real :: temp, pfrac, dz, const, d1, d2, FIXNAT, FIXICE
!                                                                              
!  V1 = SEDIMENTATION VELOCITY (M/S) OF ICE PARTICLES
!  V2 = SEDIMENTATION VELOCITY OF NAT PARTICLES
!  R1, R2 = ASSUMED RADII
!  AM1, AM2 = MOLECULAR WEIGHTS
!  RHO1, RHO2 = DENSITIES OF THE PSCs (G/CM3)
!                                                                              
real, parameter :: V1=1.27E-2, V2=1.39E-4
real, parameter :: R1=7.0E-6, R2=0.5E-6
real, parameter :: AM1=18.0, AM2=117.0
real, parameter :: RHO1=0.928, RHO2=1.35
real, parameter :: RATIO = AM1*RHO2/(AM2*RHO1)*(R2/R1)**3

il = SIZE(pfull,1)
jl = SIZE(pfull,2)
kl = SIZE(pfull,3)

!                                                                              
!  CALCULATE FRACTION OF NAT PARTICLES USED AS TYPE 2 CORES (F1)               
!  AND FRACTION OF NAT PARTICLES THAT REMAIN AS TYPE 1 CORES (F2)              
!  DETERMINE MAXIMUM NAT AND ICE TO APPLY LIMITERS TO ADVECTED AMOUNTS
!                                                                              
Lat_loop : &
   DO j = 1,jl 
Lon_loop : &
   DO i = 1,il 
      ANAT(:) = psc(i,j,:,2)
      AICE(:) = psc(i,j,:,3)
      where(AICE(:) < 1.0E-18)
         AICE(:) = 0.
      end where
      where(ANAT(:) < 1.0E-18)
         ANAT(:) = 0.
         F1(:) = 0.
         F2(:) = 0.
      elsewhere
         F1(:) = AICE(:)*RATIO/ANAT(:)
         F1(:) = MIN(1.,F1(:))
         F2(:) = 1.0 - F1(:)
      end where
      ANATMAX = maxval(ANAT(:))
      AICEMAX = maxval(AICE(:))
!
! VERTICALLY ADVECT NAT AND ICE. NOTE THAT PART OF NAT IS ADVECTED
! AT TYPE 2 RATE AND THE REMAINDER AT TYPE 1 RATE. CALCULATE DESCENT IN
!  1 TIMESTEP; USE APPROXIMATE VERTICAL DISPLACEMENT BETWEEN LAYERS
!
      TEMP = 195.
      PNAT = 0.
      PICE = 0.
      DO k = 2,kl
         PFRAC = pfull(i,j,k)/pfull(i,j,k-1)
         DZ = 29.26*TEMP*LOG(PFRAC)
         CONST = dt/DZ                                              
         D1 = ANAT(k) - ANAT(k-1)
         D2 = AICE(k) - AICE(k-1)         
         SNATS(k) = -CONST*D1*(V1*F1(k) + V2*F2(k))
         SICES(k) = -CONST*D2*V1
         PNAT = PNAT + pfull(i,j,k)*ANAT(k)
         PICE = PICE + pfull(i,j,k)*AICE(k)
      END DO
!
!  set sedimented nat and ice to zero at top and bottom
!
      SNATS(1) = 0.0
      SICES(1) = 0.0 
      SNATS(kl) = 0.0
      SICES(kl) = 0.0
      ANAT2(:) = ANAT(:) + SNATS(:)
      AICE2(:) = AICE(:) + SICES(:)
!
!  APPLY LIMITERS TO NEW NAT AND ICE
!
      ANAT2(:) = MAX( MIN(ANAT2(:),ANATMAX), 0. )
      AICE2(:) = MAX( MIN(AICE2(:),AICEMAX), 0. )
!
! APPLY MASS FIXER
!
      PNAT2 = 0.0
      PICE2 = 0.0
      DO k = 1,kl
         PNAT2 = PNAT2 + pfull(i,j,k)*ANAT2(k)
         PICE2 = PICE2 + pfull(i,j,k)*AICE2(k)
      END DO
      IF(PNAT2 == 0.) THEN 
         FIXNAT = 1.0
      ELSE 
         FIXNAT = PNAT/PNAT2
      ENDIF
      IF(PICE2 == 0.) THEN
         FIXICE = 1.0
      ELSE
         FIXICE = PICE/PICE2
      ENDIF
      ANAT2(:) = ANAT2(:)*FIXNAT 
      AICE2(:) = AICE2(:)*FIXICE 
!
!  ADJUST NOY AND H2O TENDENCY FIELDS
!
!     ANOY(j,:) = ANOY(j,:) + (ANAT2(:) - ANAT(:))/dt 
!     AHNO3(j,:) = AHNO3(j,:) + (ANAT2(:) - ANAT(:))/dt 
!     AH2O(j,:) = AH2O(j,:) + (AICE2(:) - AICE(:))/dt 

!  ADJUST PSC FIELDS
      dpsc(i,j,:,1) = 0.
      dpsc(i,j,:,2) = ANAT2(:) - ANAT(:)
      dpsc(i,j,:,3) = AICE2(:) - AICE(:)


   end do Lon_loop
   end do Lat_loop

end subroutine strat_chem_psc_sediment


! <SUBROUTINE NAME="strat_chem_get_extra_h2o">
!   <OVERVIEW>
!     Set minimum allowed stratospheric water
!   </OVERVIEW>
!   <DESCRIPTION>
!     Constrain stratospheric H2O to be greater than or equal to 2*CH4
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strat_chem_get_extra_h2o( h2o, age, ch4, Time, extra_h2o )
!   </TEMPLATE>
!   <IN NAME="h2o" TYPE="real" DIM="(:,:)">
!     Total H2O volume mixing ratio (mol/mol)
!   </IN>
!   <IN NAME="age" TYPE="real" DIM="(:,:)">
!     Age-of-air tracer (yrs)
!   </IN>
!   <IN NAME="ch4" TYPE="real" DIM="(:,:)">
!     Methane volume mixing ratio (mol/mol)
!   </IN>
!   <IN NAME="age" TYPE="time_type">
!     Current model time
!   </IN>
!   <OUT NAME="extra_h2o" TYPE="real" DIM="(:,:)">
!     Additional stratospheric H2O VMR (mol/mol)
!   </OUT>
subroutine strat_chem_get_extra_h2o( h2o, age, ch4, Time, extra_h2o )

implicit none

! Dummy arguments

real, dimension(:,:), intent(in)  :: h2o, age, ch4
type(time_type),      intent(in)  :: Time
real, dimension(:,:), intent(out) :: extra_h2o

! Local variables

integer :: i, k, il, kl, index1, index2
real :: frac, ch4_trop, min_h2o
type(time_type) :: time_trop
character(len=256) :: err_msg

il = size(h2o,1)
kl = size(h2o,2)

do k = 1,kl
do i=1,il

   if (fixed_ch4_lbc_time) then
      time_trop = ch4_entry
   else
      time_trop = increment_time( Time, -NINT(age(i,k)/tfact), 0)
   end if
   call time_interp( time_trop, ch4_time(:), frac, index1, index2, err_msg=err_msg )
   if(err_msg /= '') then
      call error_mesg('strat_chem_get_extra_h2o', trim(err_msg) , FATAL)
   endif 
   ch4_trop = ch4_value(index1) + frac*(ch4_value(index2)-ch4_value(index1))
   min_h2o = 2. * MAX( 0., ch4_trop - ch4(i,k) )
   if (age(i,k) > 0.1) then
      extra_h2o(i,k) = MAX( 0., min_h2o - h2o(i,k) )
   else
      extra_h2o(i,k) = 0.
   end if

end do
end do

end subroutine strat_chem_get_extra_h2o



! <SUBROUTINE NAME="strat_chem_extra_h2o_init">
!   <OVERVIEW>
!     Initialize minimum stratospheric water calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialize constraint of stratospheric H2O to be greater than or equal to 2*CH4
!   </DESCRIPTION>
!   <TEMPLATE>
!     call strat_chem_extra_h2o_init()
!   </TEMPLATE>
!   <IN NAME="ch4_filename" TYPE="character">
!     Methane timeseries filename
!   </IN>
!   <IN NAME="ch4_scale_factor" TYPE="real">
!     Methane timeseries scale factor to convert to VMR (mol/mol)
!   </IN>
subroutine strat_chem_extra_h2o_init( ch4_filename, ch4_scale_factor, &
                                      fixed_ch4_lbc_time_in, ch4_entry_in )

implicit none

! Dummy arguments

character(len=*), intent(in) :: ch4_filename
real, intent(in)             :: ch4_scale_factor
logical, intent(in)          :: fixed_ch4_lbc_time_in
type(time_type), intent(in)  :: ch4_entry_in

! Local variables
character(len=64) :: filename
integer :: flb, series_length, n, year, diy
real :: extra_seconds
real, dimension(:), allocatable :: input_time
type(time_type) :: Year_t

fixed_ch4_lbc_time = fixed_ch4_lbc_time_in
ch4_entry = ch4_entry_in


filename = 'INPUT/' // trim(ch4_filename)
if( file_exist(filename) ) then
   flb = open_namelist_file( filename )
   read(flb, FMT='(i12)') series_length
   allocate( ch4_value(series_length), &
             input_time(series_length), &
             ch4_time(series_length) )
   do n = 1,series_length
      read (flb, FMT = '(2f12.4)') input_time(n), ch4_value(n)
   end do
   ch4_value(:) = ch4_value(:) * ch4_scale_factor
   call close_file( flb )
!---------------------------------------------------------------------
!    convert the time stamps of the series to time_type variables.     
!---------------------------------------------------------------------
   do n=1,series_length
      year = INT(input_time(n))
      Year_t = set_date(year,1,1,0,0,0)
      diy = days_in_year(Year_t)
      extra_seconds = (input_time(n) - year)*diy*SECONDS_PER_DAY 
      ch4_time(n) = Year_t + set_time(NINT(extra_seconds), 0)
   end do
   deallocate(input_time)
else
   call error_mesg ('strat_chem_extra_h2o_init', &
                    'Failed to find input file '//trim(filename), FATAL)
end if

end subroutine strat_chem_extra_h2o_init


end module strat_chem_utilities_mod

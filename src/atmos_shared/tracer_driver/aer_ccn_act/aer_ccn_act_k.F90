        module aer_ccn_act_k_mod

implicit none
private
    private    CalcG,  erff, CalcAlphaGamma, CalcBeta
      
    public aer_ccn_act_k, aer_ccn_act2_k, aer_ccn_act_wpdf_k,  &
           aer_ccn_act_k_init, aer_ccn_act_k_end, &
           dlocate, ghquad, aer_ccn_act_wpdf_m_k

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: aer_ccn_act_k.F90,v 19.0 2012/01/06 20:31:36 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!---------------- private data -------------------

  integer, parameter :: TY = 3  !  Number of aerosol types
  integer, parameter :: MD = 2  !  Number of lognormal modes

  real :: T = 283.15 !  Temperature (K)
  real :: P = 0.800e5 !  Pressure (Pa)
  real, parameter :: R = 8.314  !  Gas constant (J mol-1 K-1)
  real, parameter :: ZERO = 273.15 !  Zero degree C (K)
  real, parameter :: ATM = 1.01325e5 !  Standard atmosphere pressure (Pa)
  real, parameter :: PI = 3.1415926
  real, parameter :: eps = 1.e-5  ! epsilon used to prevent index
                                     ! calculation errors

  
!NO 1 Ammonium Sulfate  NO 2 Sea Salt NO 3 Organics
  
  real, dimension(TY) :: B_term = (/0.7822,0.6342,1.3764/) ! 2 * 0.3492/(Bprim)**(1/3)
  real, dimension(TY) :: Mass_scal = (/0.15896,0.198,0.1241/) ! scaling mass (ug m-3)

  real, dimension(MD) :: N = (/340., 60./) ! Total Number Concen (cm-3)
  real, dimension(MD) :: Dm = (/0.01, 0.07/) ! Geometric mean diameter (micron)
  real, dimension(MD) :: LNSIGMA = (/0.47, 0.6931/) ! ln( Sigma (St. Dev.) )

  logical :: nooc

!Parameters for look-up tables

  integer :: res, res2
  real    :: sul_concen, low_concen, high_concen
  real    :: lowup, highup, lowup2, highup2, lowmass2, &
             highmass2, lowmass3, highmass3,  &
             lowmass4, highmass4, lowmass5, highmass5, &
             lowT2, highT2

! real ::  lowup=0.3 !m/s
! real ::  highup=10.

! real ::  lowup2=0.0001 !m/s
! real ::  lowup2=0.05   !m/s
! real ::  highup2=0.3
! real ::  lowmass2=0.01 !ug m-3
! real ::  highmass2=1000.
! real ::  highmass2=100.
! real ::  lowmass3=0.01 !ug m-3
! real ::  highmass3=1000.
! real ::  highmass3=100.
! real :: lowT2=243.15 !K
! real :: highT2=308.15

real, dimension(:,:,:,:,:), allocatable  :: droplets, droplets2

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

logical :: module_is_initialized  = .false.
 
contains

subroutine aer_ccn_act_k_init      &
             (droplets_in, droplets2_in, res_in, res2_in, nooc_in,  &
              sul_concen_in, low_concen_in, high_concen_in, &
              lowup_in, highup_in, lowup2_in, highup2_in, lowmass2_in, &
              highmass2_in, lowmass3_in, highmass3_in,  &
              lowmass4_in, highmass4_in, lowmass5_in, highmass5_in, &
              lowT2_in, highT2_in)     

real, dimension(:,:,:,:,:), intent(in) :: droplets_in
real, dimension(:,:,:,:,:), intent(in) :: droplets2_in
integer, intent(in) :: res_in, res2_in
logical, intent(in) :: nooc_in
real, intent(in)    :: sul_concen_in, low_concen_in, high_concen_in
real, intent(in)    :: lowup_in, highup_in, lowup2_in, highup2_in,  &
                       lowmass2_in, highmass2_in, lowmass3_in,  &
                       highmass3_in, lowmass4_in, highmass4_in,  &
                       lowmass5_in, highmass5_in, lowT2_in, highT2_in

    if (module_is_initialized) return

    res = res_in
    res2 = res2_in

    allocate (droplets (res, res,res,res,res))
    allocate (droplets2 (res2, res2,res2,res2,res2))

    droplets = droplets_in
    droplets2 = droplets2_in
    nooc = nooc_in
    sul_concen = sul_concen_in
    low_concen = low_concen_in
    high_concen = high_concen_in
    lowup = lowup_in
    highup = highup_in
    lowup2 = lowup2_in
    highup2 = highup2_in
    lowmass3 = lowmass3_in
    highmass3 = highmass3_in
    lowmass4 = lowmass4_in
    highmass4 = highmass4_in
    lowmass5 = lowmass5_in
    highmass5 = highmass5_in
    lowmass2 = lowmass2_in
    highmass2 = highmass2_in
    lowT2 = lowT2_in
    highT2 = highT2_in

    module_is_initialized = .true.


end subroutine aer_ccn_act_k_init



subroutine aer_ccn_act_k (T1, P1, Updraft1, TotalMass, tym,      &      
                          Drop, ier, ermesg)
integer, intent(in) :: tym
real, dimension(tym), intent(inout) :: TotalMass
real, intent(in) :: T1, P1, Updraft1
real, intent(inout) :: Drop
integer, intent(out) :: ier
character(len=*), intent(out) :: ermesg
    
real tmass, tmass2, tmass3, tmass4, updr
integer nomass, nomass2, nomass3, nomass4, noup
        
  ier = 0

  if ( .not. module_is_initialized) then
    ier = 1
    ermesg = 'aer_ccn_act_k: module has not been initialized before &
                                &first call'
    return
  endif
  if (tym /=  4) then
    ier = 2
    ermesg = 'aer_ccn_act_k:dimension of TotalMass is incorrect'
    return
  endif
  
  if (nooc) TotalMass(4) = 0.

  tmass=(TotalMass(1)+TotalMass(3)+TotalMass(4))*1.e12
    
  if (Updraft1>lowup2 .and. tmass>lowmass2) then

    tmass3=TotalMass(2)*1.e12
    if(TotalMass(1)*1.e12<sul_concen) then
      if (tmass3<low_concen) then
        tmass=TotalMass(1)*1.e12
        tmass2=0.
      else if (tmass3>high_concen) then
        tmass2=TotalMass(1)*1.e12
        tmass=0.
      else
        tmass=TotalMass(1)*1.e12* &
                      (high_concen-tmass3)/(high_concen-low_concen)
        tmass2=TotalMass(1)*1.e12* &
                      (tmass3-low_concen)/(high_concen-low_concen)
      end if
    else
      tmass2=TotalMass(1)*1.e12
      tmass=0.
    end if
   
    tmass3=TotalMass(3)*1.e12
    tmass4=TotalMass(4)*1.e12
    
    if (Updraft1>highup2) then
    
      updr=max(min(Updraft1,highup-eps),lowup)
    
      noup= log(updr/lowup)/log(highup/lowup)*(res2-1.)

      tmass=max(min(tmass,highmass2-eps),lowmass2)
      nomass= log(tmass/lowmass2)/log(highmass2/lowmass2)*(res2-1.)

      tmass2=max(min(tmass2,highmass3-eps),lowmass3)
      nomass2= log(tmass2/lowmass3)/log(highmass3/lowmass3)*(res2-1.)
    
      tmass3=max(min(tmass3,highmass4-eps),lowmass4)
      nomass3= log(tmass3/lowmass4)/log(highmass4/lowmass4)*(res2-1.)
 
      tmass4=max(min(tmass4,highmass5-eps),lowmass5)
      nomass4= log(tmass4/lowmass5)/log(highmass5/lowmass5)*(res2-1.)

      Drop = 0.166667*(droplets2(noup+1,nomass4+1,nomass3+1,nomass2+1,nomass+1)+&
                       droplets2(noup+1,nomass4+1,nomass3+1,nomass2+1,nomass+2)+ &
                       droplets2(noup+1,nomass4+1,nomass3+1,nomass2+2,nomass+1)+ &
                       droplets2(noup+1,nomass4+1,nomass3+2,nomass2+1,nomass+1)+ &
                       droplets2(noup+1,nomass4+2,nomass3+1,nomass2+1,nomass+1)+ &
                       droplets2(noup+2,nomass4+1,nomass3+1,nomass2+1,nomass+1))
  
    else

      updr=max(min(Updraft1,highup2-eps),lowup2)
    
      noup= log(updr/lowup2)/log(highup2/lowup2)*(res-1.)

      tmass=max(min(tmass,highmass2-eps),lowmass2)
      nomass= log(tmass/lowmass2)/log(highmass2/lowmass2)*(res-1.)

      tmass2=max(min(tmass2,highmass3-eps),lowmass3)
      nomass2= log(tmass2/lowmass3)/log(highmass3/lowmass3)*(res-1.)
    
      tmass3=max(min(tmass3,highmass4-eps),lowmass4)
      nomass3= log(tmass3/lowmass4)/log(highmass4/lowmass4)*(res-1.)

      tmass4=max(min(tmass4,highmass5-eps),lowmass5)
      nomass4= log(tmass4/lowmass5)/log(highmass5/lowmass5)*(res-1.)

      Drop = 0.166667*(droplets(noup+1,nomass4+1,nomass3+1,nomass2+1,nomass+1)+&
                       droplets(noup+1,nomass4+1,nomass3+1,nomass2+1,nomass+2)+ &
                       droplets(noup+1,nomass4+1,nomass3+1,nomass2+2,nomass+1)+ &
                       droplets(noup+1,nomass4+1,nomass3+2,nomass2+1,nomass+1)+ &
                       droplets(noup+1,nomass4+2,nomass3+1,nomass2+1,nomass+1)+ &
                       droplets(noup+2,nomass4+1,nomass3+1,nomass2+1,nomass+1))
  
    endif
  
  
  else
  
    Drop=0.
  
  endif
  
end subroutine aer_ccn_act_k

subroutine aer_ccn_act2_k (T1, P1, Updraft1, TotalMass, tym, mu, &
                           airdens,Nc,qc,qt,qe,tc,te,Drop, ier, ermesg)

!T1 temperature (K)
!P1 pressure (Pa)
!Updraft1 updraft velocity (m/s)
!TotalMass aerosol mass ()
!mu entrainment coef. (/s)
!airdens air density (kg/m3 air)
!Nc droplet mixing ratio (#/kg air)
!qc in-cloud vapor mixing ratio (kg water/kg air)
!qt in-cloud total water mixing ratio qc + ql (kg water/kg air)
!qe environment vapor mixing ratio (kg water/kg air)
!tc in-cloud temperature (K)
!te environment temperature (K)
!Drop droplet number concentration (#/cc)

integer, intent(in) :: tym
real, dimension(tym), intent(in) :: TotalMass
real, intent(in) :: T1, P1, Updraft1, mu,airdens, Nc, qc, qt, qe, tc, te
real, intent(inout) :: Drop
integer, intent(out) :: ier
character(len=*), intent(out) :: ermesg

real :: G, alpha, gamma, Smax
real :: Diam, beta, Le_cpa, Dcut
integer :: i, j
        
  ier = 0
  if ( .not. module_is_initialized) then
    ier = 1
    ermesg = 'aer_ccn_act2_k: module has not been initialized before &
                                &first call'
    return
  endif
  if (tym /= 4) then
    ier = 2
    ermesg = ' aer_ccn_act2_k:dimension of TotalMass is incorrect'
    return
  endif

  Drop=0.
  
  if (Nc > 0.) then    
    T = T1
    P = P1

    call CalcAlphaGamma(alpha, gamma)
    call CalcBeta(beta, Le_cpa)
              
! Diam  average diameter of droplets (micron)
    if (qt > qc) then
      Diam= ((qt-qc)/Nc*1.91e15)**(1./3.)    
    else
      Diam= 20.
    endif

!set the upper and lower limits of Diam
    if (Diam < 10.) Diam=10.
  
    call calcG(Diam, G)
    
    Smax=(alpha-gamma*mu*(qt-qe)*airdens+beta*mu*(Le_cpa*(qc-qe)+(tc-te)))*Updraft1/ &
         (gamma)/(0.5*3.1415*1.e3*G*Diam*1.e-6*Nc*airdens)
  
    if (Smax>0.) then
      do i=1,TY
        Dcut=B_term(i)/T/(Smax**(2./3.))
        do j=1, MD
          Drop=Drop+TotalMass(i)/Mass_scal(i)*N(j)*0.5* &
               (1.-erff(log(Dcut/Dm(j))/LNSIGMA(j)*0.707107))      
        end do
      end do
    endif
  endif

end subroutine aer_ccn_act2_k

!-->cjg: addition
!
! Additional subroutines to compute CCN activation by integrating
! over an assumed subgrid-scale PDF of w

subroutine aer_ccn_act_wpdf_k (T, p, wm, wp2, totalmass, tym, drop,  &
                               ier, ermesg)

! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

integer,          intent(in)    :: tym
real,             intent(in)    :: T, p, wm, wp2
real,             intent(inout) :: totalmass(tym)
real,             intent(out)   :: drop
integer,          intent(out)   :: ier
character(len=*), intent(out)   :: ermesg

!  Parameters

real, parameter    :: wp2_eps = 0.0001 ! w variance threshold
integer, parameter :: npoints = 64     ! # for Gauss-Hermite quadrature
real, parameter    :: wmin =  0.0      ! min w for ccn_act
real, parameter    :: wmax = 10.0      ! max w for ccn_act

!  Internal

real(kind=8), dimension(npoints) :: x, w
integer init

logical lintegrate
integer i, ia, ib
real wtmp
real(kind=8) :: tmp, a, b, sum1, sum2

save init, x, w

data init/0/


  ier = 0
  if ( .not. module_is_initialized) then
    ier = 1
    ermesg = 'aer_ccn_act_wpdf _k: module has not been initialized &
                                             &before first call'
    return
  endif
  if (tym /=  4) then
    ier = 2
    ermesg = 'aer_ccn_act_wpdf_k:dimension of TotalMass is incorrect'
    return
  endif

! On first call, initialize arrays with abscissas and weights for
! integration

if ( init .eq. 0 ) then
  call ghquad( npoints, x, w )
  init = 1
endif

! Determine whether integration is needed to compute number
! of activated drops. lintegrate = .true. indicates that
! numerical integration is to be performed.

lintegrate = .false.
if ( wp2 .gt. wp2_eps ) then

  ! Integration bounds: from wmin to wmax (0 to 10 m/s)

  tmp = 1.0d0 / sqrt(2.0 * wp2 ) 
  a = (wmin - wm ) * tmp
  b = (wmax - wm ) * tmp

  ! Locate indices within integration bounds

  call dlocate( x, npoints, a, ia )
  call dlocate( x, npoints, b, ib )

  ! ia (ib) is zero if a (b) is smaller than the lowest abscissa.
  ! In that case, start the integration with the first abscissa.
!  ia = max(ia,1)
!  ib = max(ib,1)
  ia = min(max(ia,1),size(x))
  ib = min(max(ib,1),size(x))

  if ( ib .gt. ia ) lintegrate = .true.

endif

! Compute number of activated drops.

if (lintegrate) then

  ! Perform integration

  sum1 = 0.0d0
  sum2 = 0.0d0
  tmp = sqrt(2.0 * wp2 )
  do i=ia,ib
    wtmp = tmp * x(i) + wm
    call aer_ccn_act_k( T, p, wtmp, totalmass, TYm, drop, ier, ermesg )
    if (ier /= 0) return
    sum1 = sum1 + w(i)*drop
    sum2 = sum2 + w(i)
  enddo

  ! Normalize

  drop = sum1 / sum2

else

  ! No integration, use single point evaluation

  call aer_ccn_act_k( T, p, wm, totalmass, TYm, drop, ier, ermesg )
  if (ier /= 0) return

endif


end subroutine aer_ccn_act_wpdf_k

!----------------------------------------------------------------------
subroutine aer_ccn_act_wpdf_m_k (T, p, wm, wp2, offs, totalmass, tym, drop,  &
!cms                               ier, ermesg)
                                       ier, ermesg)

! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

integer,          intent(in)    :: tym
real,             intent(in)    :: T, p, wm, wp2
integer,          intent(in)    :: offs
real,             intent(inout) :: totalmass(tym)
real,             intent(out)   :: drop
integer,          intent(out)   :: ier
character(len=*), intent(out)   :: ermesg

!  Parameters

real, parameter    :: wp2_eps = 0.0001 ! w variance threshold
integer, parameter :: npoints = 64     ! # for Gauss-Hermite quadrature
real, parameter    :: wmin =  0.0      ! min w for ccn_act
real, parameter    :: wmax = 10.0      ! max w for ccn_act

!  Internal

real(kind=8), dimension(npoints) :: x, w
integer init

logical lintegrate
integer i, ia, ib

real wtmp
real(kind=8) :: tmp, a, b, sum1, sum2

save init, x, w

data init/0/


  ier = 0
  if ( .not. module_is_initialized) then
    ier = 1
    ermesg = 'aer_ccn_act_wpdf _k: module has not been initialized &
                                             &before first call'
    return
  endif
  if (tym /=  4) then
    ier = 2
    ermesg = 'aer_ccn_act_wpdf_k:dimension of TotalMass is incorrect'
    return
  endif

! On first call, initialize arrays with abscissas and weights for
! integration

if ( init .eq. 0 ) then
  call ghquad( npoints, x, w )
  init = 1
endif

! Determine whether integration is needed to compute number
! of activated drops. lintegrate = .true. indicates that
! numerical integration is to be performed.

lintegrate = .false.
if ( wp2 .gt. wp2_eps ) then

  ! Integration bounds: from wmin to wmax (0 to 10 m/s)

  tmp = 1.0d0 / sqrt(2.0 * wp2 ) 
  a = (wmin - wm ) * tmp
  b = (wmax - wm ) * tmp

  ! Locate indices within integration bounds

  call dlocate( x, npoints, a, ia )
  call dlocate( x, npoints, b, ib )

  ! ia (ib) is zero if a (b) is smaller than the lowest abscissa.
  ! In that case, start the integration with the first abscissa.
!  ia = max(ia,1)
!  ib = max(ib,1)

!cms++ 
  ia = ia + offs
!cms--

  ia = min(max(ia,1),size(x))
  ib = min(max(ib,1),size(x))

  if ( ib .gt. ia ) lintegrate = .true.

endif

! Compute number of activated drops.

if (lintegrate) then

  ! Perform integration

  sum1 = 0.0d0
  sum2 = 0.0d0
  tmp = sqrt(2.0 * wp2 )
  do i=ia,ib
    wtmp = tmp * x(i) + wm
    call aer_ccn_act_k( T, p, wtmp, totalmass, TYm, drop, ier, ermesg )
    if (ier /= 0) return
    sum1 = sum1 + w(i)*drop
    sum2 = sum2 + w(i)
  enddo

  ! Normalize

  drop = sum1 / sum2

else

  ! No integration, use single point evaluation

  call aer_ccn_act_k( T, p, wm, totalmass, TYm, drop, ier, ermesg )
  if (ier /= 0) return

endif


end subroutine aer_ccn_act_wpdf_m_k


subroutine dlocate(xx,n,x,j)

! Subroutine to locate the position of element in an ordered array

integer,      intent(in)  :: n
real(kind=8), intent(in)  :: x, xx(n)
integer,      intent(out) :: j
integer jl,jm,ju

jl=0
ju=n+1
do while (ju-jl.gt.1)
  jm=(ju+jl)/2
  if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
    jl=jm
  else
    ju=jm
  endif
enddo
j=jl

return
end subroutine dlocate

subroutine ghquad( n, x, w )

! This subroutine returns double precision abscissas [x(1:n)] and
! weights [w(1:n)] of a n-point Gauss-Hermite quadrature formula,
! with n = {8, 12, 16, 24, 32, 48, 64, 96}.
! 
! The values of the absicssas and weights were computed offline 
! using subroutine 'gaussq' from netlib.org

integer, intent(in)       :: n
real(kind=8), intent(out) :: x(n), w(n)

select case (n)

  case(  8)
    x(1:n/2) =                                      &
    (/                                              &
      -0.2930637420257D+01, -0.1981656756696D+01,   &
      -0.1157193712447D+01, -0.3811869902073D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.1996040722114D-03,  0.1707798300741D-01,   &
       0.2078023258149D+00,  0.6611470125582D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 12)
    x(1:n/2) =                                      &
    (/                                              &
      -0.3889724897870D+01, -0.3020637025121D+01,   &
      -0.2279507080501D+01, -0.1597682635153D+01,   &
      -0.9477883912402D+00, -0.3142403762544D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.2658551684356D-06,  0.8573687043588D-04,   &
       0.3905390584629D-02,  0.5160798561588D-01,   &
       0.2604923102642D+00,  0.5701352362625D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 16)
    x(1:n/2) =                                      &
    (/                                              &
      -0.4688738939306D+01, -0.3869447904860D+01,   &
      -0.3176999161980D+01, -0.2546202157847D+01,   &
      -0.1951787990916D+01, -0.1380258539199D+01,   &
      -0.8229514491447D+00, -0.2734810461382D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.2654807474011D-09,  0.2320980844865D-06,   &
       0.2711860092538D-04,  0.9322840086242D-03,   &
       0.1288031153551D-01,  0.8381004139899D-01,   &
       0.2806474585285D+00,  0.5079294790166D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 24)
    x(1:n/2) =                                      &
    (/                                              &
      -0.6015925561426D+01, -0.5259382927668D+01,   &
      -0.4625662756424D+01, -0.4053664402448D+01,   &
      -0.3520006813035D+01, -0.3012546137566D+01,   &
      -0.2523881017011D+01, -0.2049003573662D+01,   &
      -0.1584250010962D+01, -0.1126760817611D+01,   &
      -0.6741711070372D+00, -0.2244145474725D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.1664368496489D-15,  0.6584620243078D-12,   &
       0.3046254269988D-09,  0.4018971174941D-07,   &
       0.2158245704902D-05,  0.5688691636404D-04,   &
       0.8236924826884D-03,  0.7048355810073D-02,   &
       0.3744547050323D-01,  0.1277396217846D+00,   &
       0.2861795353464D+00,  0.4269311638687D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 32)
    x(1:n/2) =                                      &
    (/                                              &
      -0.7125813909831D+01, -0.6409498149270D+01,   &
      -0.5812225949516D+01, -0.5275550986516D+01,   &
      -0.4777164503503D+01, -0.4305547953351D+01,   &
      -0.3853755485471D+01, -0.3417167492819D+01,   &
      -0.2992490825002D+01, -0.2577249537732D+01,   &
      -0.2169499183606D+01, -0.1767654109463D+01,   &
      -0.1370376410953D+01, -0.9765004635897D+00,   &
      -0.5849787654359D+00, -0.1948407415694D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.7310676427384D-22,  0.9231736536518D-18,   &
       0.1197344017093D-14,  0.4215010211326D-12,   &
       0.5933291463397D-10,  0.4098832164771D-08,   &
       0.1574167792546D-06,  0.3650585129562D-05,   &
       0.5416584061820D-04,  0.5362683655280D-03,   &
       0.3654890326654D-02,  0.1755342883157D-01,   &
       0.6045813095591D-01,  0.1512697340766D+00,   &
       0.2774581423025D+00,  0.3752383525928D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 48)
    x(1:n/2) =                                      &
    (/                                              &
      -0.8975315081932D+01, -0.8310752190705D+01,   &
      -0.7759295519766D+01, -0.7266046554164D+01,   &
      -0.6810064578074D+01, -0.6380564096186D+01,   &
      -0.5971072225014D+01, -0.5577316981224D+01,   &
      -0.5196287718792D+01, -0.4825757228133D+01,   &
      -0.4464014546934D+01, -0.4109704603561D+01,   &
      -0.3761726490228D+01, -0.3419165969364D+01,   &
      -0.3081248988645D+01, -0.2747308624822D+01,   &
      -0.2416760904873D+01, -0.2089086660944D+01,   &
      -0.1763817579895D+01, -0.1440525220138D+01,   &
      -0.1118812152402D+01, -0.7983046277786D+00,   &
      -0.4786463375945D+00, -0.1594929358489D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.7935551460774D-35,  0.5984612693314D-30,   &
       0.3685036080151D-26,  0.5564577468902D-23,   &
       0.3188387323505D-20,  0.8730159601187D-18,   &
       0.1315159622658D-15,  0.1197589865479D-13,   &
       0.7046932581546D-12,  0.2815296537838D-10,   &
       0.7930467495165D-09,  0.1622514135896D-07,   &
       0.2468658993670D-06,  0.2847258691735D-05,   &
       0.2528599027748D-04,  0.1751504318012D-03,   &
       0.9563923198194D-03,  0.4153004911978D-02,   &
       0.1444496157498D-01,  0.4047967698460D-01,   &
       0.9182229707929D-01,  0.1692044719456D+00,   &
       0.2539615426648D+00,  0.3110010303780D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 64)
    x(1:n/2) =                                      &
    (/                                              &
      -0.1052612316796D+02, -0.9895287586830D+01,   &
      -0.9373159549647D+01, -0.8907249099965D+01,   &
      -0.8477529083380D+01, -0.8073687285010D+01,   &
      -0.7689540164040D+01, -0.7321013032781D+01,   &
      -0.6965241120551D+01, -0.6620112262636D+01,   &
      -0.6284011228775D+01, -0.5955666326799D+01,   &
      -0.5634052164350D+01, -0.5318325224633D+01,   &
      -0.5007779602199D+01, -0.4701815647408D+01,   &
      -0.4399917168228D+01, -0.4101634474567D+01,   &
      -0.3806571513945D+01, -0.3514375935741D+01,   &
      -0.3224731291992D+01, -0.2937350823005D+01,   &
      -0.2651972435431D+01, -0.2368354588632D+01,   &
      -0.2086272879882D+01, -0.1805517171466D+01,   &
      -0.1525889140210D+01, -0.1247200156943D+01,   &
      -0.9692694230712D+00, -0.6919223058100D+00,   &
      -0.4149888241211D+00, -0.1383022449870D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.5535706535857D-48,  0.1679747990108D-42,   &
       0.3421138011256D-38,  0.1557390624630D-34,   &
       0.2549660899113D-31,  0.1929103595465D-28,   &
       0.7861797788926D-26,  0.1911706883301D-23,   &
       0.2982862784280D-21,  0.3152254566504D-19,   &
       0.2351884710676D-17,  0.1280093391322D-15,   &
       0.5218623726591D-14,  0.1628340730710D-12,   &
       0.3959177766948D-11,  0.7615217250145D-10,   &
       0.1173616742322D-08,  0.1465125316476D-07,   &
       0.1495532936727D-06,  0.1258340251031D-05,   &
       0.8788499230850D-05,  0.5125929135786D-04,   &
       0.2509836985131D-03,  0.1036329099508D-02,   &
       0.3622586978534D-02,  0.1075604050988D-01,   &
       0.2720312895369D-01,  0.5873998196410D-01,   &
       0.1084983493062D+00,  0.1716858423491D+00,   &
       0.2329947860627D+00,  0.2713774249413D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case( 96)
    x(1:n/2) =                                      &
    (/                                              &
      -0.1311613002166D+02, -0.1252923074668D+02,   &
      -0.1204480674113D+02, -0.1161362082884D+02,   &
      -0.1121687519420D+02, -0.1084488807385D+02,   &
      -0.1049185453928D+02, -0.1015395127471D+02,   &
      -0.9828492746192D+01, -0.9513501769412D+01,   &
      -0.9207469353372D+01, -0.8909210581460D+01,   &
      -0.8617773244215D+01, -0.8332377307170D+01,   &
      -0.8052373330720D+01, -0.7777213031061D+01,   &
      -0.7506427894428D+01, -0.7239613294008D+01,   &
      -0.6976416464285D+01, -0.6716527240497D+01,   &
      -0.6459670819555D+01, -0.6205602024718D+01,   &
      -0.5954100706412D+01, -0.5704968013526D+01,   &
      -0.5458023340071D+01, -0.5213101801819D+01,   &
      -0.4970052133167D+01, -0.4728734920353D+01,   &
      -0.4489021106217D+01, -0.4250790715903D+01,   &
      -0.4013931763628D+01, -0.3778339308797D+01,   &
      -0.3543914636027D+01, -0.3310564538524D+01,   &
      -0.3078200688047D+01, -0.2846739077725D+01,   &
      -0.2616099526362D+01, -0.2386205234770D+01,   &
      -0.2156982386193D+01, -0.1928359784146D+01,   &
      -0.1700268521938D+01, -0.1472641679023D+01,   &
      -0.1245414039907D+01, -0.1018521831948D+01,   &
      -0.7919024787540D+00, -0.5654943662667D+00,   &
      -0.3392366188789D+00, -0.1130688831515D+00    &
    /)
    x(n/2+1:n) = -x(n/2:1:-1)
    w(1:n/2) =                                      &
    (/                                              &
       0.1315337147701D-74,  0.3480841387719D-68,   &
       0.4468702421681D-63,  0.1093382473752D-58,   &
       0.8733732835919D-55,  0.3022293931502D-51,   &
       0.5383222506970D-48,  0.5537064819783D-45,   &
       0.3568938561781D-42,  0.1531798964711D-39,   &
       0.4587401665859D-37,  0.9946800059220D-35,   &
       0.1608873207938D-32,  0.1989546894270D-30,   &
       0.1919970706795D-28,  0.1471241554974D-26,   &
       0.9085965647049D-25,  0.4580601098550D-23,   &
       0.1906260188688D-21,  0.6612951176481D-20,   &
       0.1928879326837D-18,  0.4766831171186D-17,   &
       0.1004905149157D-15,  0.1818169130493D-14,   &
       0.2838767395169D-13,  0.3843681724877D-12,   &
       0.4533303374584D-11,  0.4676010385985D-10,   &
       0.4233593783511D-09,  0.3375584731128D-08,   &
       0.2377366214757D-07,  0.1482969405005D-06,   &
       0.8213544547833D-06,  0.4048233178777D-05,   &
       0.1779180591822D-04,  0.6985432444541D-04,   &
       0.2454180927305D-03,  0.7726958216817D-03,   &
       0.2183131976097D-02,  0.5541631289439D-02,   &
       0.1265137511280D-01,  0.2600034027124D-01,   &
       0.4813990673109D-01,  0.8035395008064D-01,   &
       0.1209831168625D+00,  0.1643796305985D+00,   &
       0.2016130134792D+00,  0.2232700234975D+00    &
    /)
    w(n/2+1:n) = w(n/2:1:-1)

  case default
     stop 'ghquad: invalid value of n'

end select

return
end subroutine ghquad

!
!<--cjg: end of addition


subroutine CalcG(Dp, G)

real, intent(inout) :: Dp
real, intent(inout) :: G
real :: rhow = 1.0e3  ! density of water (Kg m-3)
real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
real :: alpc = 1.0  ! mass accomodation coef. 
real :: alpt = 0.97  ! thermal accomodation coef.
!real :: alpc = 0.042  ! mass accomodation coef. 
!real :: alpt = 1.  ! thermal accomodation coef.
!real :: alpc = 0.2  ! mass accomodation coef. 
!real :: alpt = 1.  ! thermal accomodation coef.
real :: delt = 2.16e-1 !thermal jump (micron)
real :: delv = 1.096e-1 !vapor jump (micron)
real vpres, Dv, ka, Le, mass, heat, TC
      
      
      Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
      Dv = Dv/(Dp/(Dp+delv*2.)+2*Dv/(alpc*(Dp*1e-6))*(2.*PI*Mw/R/T)**0.5)
      ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
      ka = ka/(Dp/(Dp+delt*2.)+2*ka/(alpt*(Dp*1e-6)*1.007e3*P)*(2.*PI*R*T/0.028965)**0.5)
      TC = T-ZERO
      vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
              +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure(Pa)
      Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1) 
      
      mass = rhow*R*T/(vpres*Dv*Mw)
      heat = Le*rhow/(ka*T)*(Le*Mw/T/R-1)
!      print *, Dv, vpres, Mw, rhow, mass, heat
      G = 4./(mass+heat) ! (m2 s-1)
end subroutine CalcG

recursive function erff(x) RESULT(y)

! Error function from Numerical Recipes.
! erf(x) = 1 - erfc(x)

real dumerfc, x
real t, z, y


  z = abs(x)
  t = 1.0 / ( 1.0 + 0.5 * z )

  dumerfc =     t * exp(-z * z - 1.26551223 + t *      &
            ( 1.00002368 + t * ( 0.37409196 + t *    &
            ( 0.09678418 + t * (-0.18628806 + t *    &
            ( 0.27886807 + t * (-1.13520398 + t *    &
            ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

  if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
 
  y = 1.0 - dumerfc

end function erff


subroutine CalcAlphaGamma(alpha, gamma)

  real, intent(inout) :: alpha, gamma
  real rhoa ! density of air (Kg m-3)
  real :: Cpa = 1.007e3 ! specific heat of air (J Kg-1 K-1)
  real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
  real :: Ma = 0.028965  ! molecular weight of air (Kg mol-1)
  real :: g = 9.815 ! gravitational acceleration (m s-2) 
  real vpres, Dv, ka, Le, TC
  
  rhoa = P*Ma/R/T  ! (Kg m-3)
  Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!  Dv = Dv/(1+2*Dv/(alpc*Dp)*(2.*PI*Mw/R/T)**0.5)
  ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
  TC = T-ZERO
  vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
          +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure (Pa)
  Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)
  alpha = g*Mw*Le/(Cpa*R*T**2.)-g*Ma/(R*T) ! (m-1)
  gamma = R*T/(vpres*Mw)+Mw*Le**2./(Cpa*P*Ma*T) ! (m3 Kg-1)
end subroutine CalcAlphaGamma

subroutine CalcBeta(beta, Le_cpa)

  real, intent(inout) :: beta, Le_cpa 
  real rhoa ! density of air (Kg m-3)
  real :: Cpa = 1.007e3 ! specific heat of air (J Kg-1 K-1)
  real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
  real :: Ma = 0.028965  ! molecular weight of air (Kg mol-1)
  real vpres, Dv, ka, Le, TC
  
  rhoa = P*Ma/R/T  ! (Kg m-3)
  Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!      Dv = Dv/(1+2*Dv/(alpc*Dp)*(2.*PI*Mw/R/T)**0.5)
  ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
  TC = T-ZERO
  vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
          +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure (Pa)
  Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)
  Le_cpa = Le/Cpa
  beta = Mw*Le/(R*T**2.) ! (K-1)
end subroutine CalcBeta

subroutine aer_ccn_act_k_end()

  module_is_initialized = .false.

end subroutine aer_ccn_act_k_end

end module aer_ccn_act_k_mod

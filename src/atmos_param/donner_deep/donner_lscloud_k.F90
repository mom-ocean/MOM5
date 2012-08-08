
!VERSION NUMBER:
!  $Id: donner_lscloud_k.F90,v 19.0 2012/01/06 20:07:28 fms Exp $

!module donner_lscloud_inter_mod

!#include "donner_types.h"

!end module donner_lscloud_inter_mod


!#####################################################################

subroutine don_l_lscloud_driver_k   &
         (isize, jsize, nlev_lsm, cloud_tracers_present, Param,  &
          Col_diag, pfull, temp, exit_flag,   &
          mixing_ratio, qlin, qiin, qain, phalf, Don_conv, &
          donner_humidity_factor, donner_humidity_area, dql, dqi, dqa, &
          mhalf_3d, &
          ermesg, error) 

!---------------------------------------------------------------------
!    subroutine don_l_lscloud_driver obtains variables needed by 
!    strat_cloud_mod that are dependent on the donner_deep parameter-
!    ization. specifically, the convective cell plus mesoscale anvil
!    cloud fraction (donner_humidity_area), the ratio of the large-scale 
!    specific humidity to the specific humidity in the environment out-
!    side of the convective system (donner_humidity_ratio), and the 
!    changes in cloud liquid, cloud ice and cloud area due to the con-
!    vective-system vertical mass flux and detrainment from the mesoscale
!    anvil to the large scale (dql, dqi, dqa) are passed out for use in 
!    strat_cloud_mod.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_param_type, &
                             donner_column_diag_type

implicit none

!---------------------------------------------------------------------
integer,                    intent(in)    :: isize, jsize, nlev_lsm
logical,                    intent(in)    :: cloud_tracers_present
type(donner_param_type),    intent(in)    :: Param
type(donner_column_diag_type),                    &
                            intent(in)    :: Col_diag
logical, dimension(isize,jsize), intent(in) :: exit_flag
real, dimension(isize,jsize,nlev_lsm),         &
                            intent(in)    :: pfull, temp, mixing_ratio, &
                                             qlin, qiin, qain
real, dimension(isize,jsize,nlev_lsm+1),        &
                            intent(in)    :: phalf 
type(donner_conv_type),     intent(inout) :: Don_conv
real, dimension(isize,jsize,nlev_lsm),           &
                            intent(out)   :: donner_humidity_factor,  &
                                             donner_humidity_area, dql, &
                                             dqi, dqa
real, dimension(isize, jsize, nlev_lsm+1),   &
                            intent(out)   :: mhalf_3d
character(len=*),           intent(out)   :: ermesg
integer,                    intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field on model half levels [ Pa ]
!     temp           temperature field at model full levels [ deg K ]
!     mixing_ratio   water vapor specific humidity at model full 
!                    levels [ kg(h2o) / kg(air) ]
!     qlin           large-scale cloud liquid specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qiin           large-scale cloud ice specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qain           large-scale cloud fraction [ fraction ]
!
!   intent(inout) variables:
!
!     Don_conv
!
!
!   intent(out) variables:
!
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to the
!                    specific humidity in the environment outside
!                    of the convective system [ dimensionless ]
!     donner_humidity_area
!                    fractional area of cell plus meso circulation
!                    associated with donner_deep_mod [ fraction ]
!     dql            increment to large-scale cloud liquid field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqi            increment to large-scale cloud ice field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqa            increment to large-scale cloud area field from
!                    donner_deep_mod [ fraction ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:


      real, dimension   &
            (isize, jsize, nlev_lsm  ) :: dmeso_3d

!---------------------------------------------------------------------
!   local variables:
!
!     dmeso_3d       detrainment rate from convective system 
!                    [ sec**(-1) ]
!     mhalf_3d       mass flux at model half-levels 
!                    [ kg / (m**2 sec) ]
!---------------------------------------------------------------------

      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    call define_donner_mass_flux to define the convective system 
!    detrainment rate (dmeso_3d) and the mass flux at model interface 
!    levels (mhalf_3d) that is associated with deep convection.
!---------------------------------------------------------------------
      call don_l_define_mass_flux_k    &
           (isize, jsize, nlev_lsm, pfull, phalf, Don_conv,   &
            dmeso_3d, mhalf_3d, exit_flag, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

 
!---------------------------------------------------------------------
!    call adjust_tiedtke_inputs to obtain the convective cloud area 
!    (donner_humidity_area) and the ratio of large-scale specific humid-
!    ity to the humidity in the environment of the convective system 
!    (donner_humidity_ratio).
!---------------------------------------------------------------------
      call don_l_adjust_tiedtke_inputs_k    &
           (isize, jsize, nlev_lsm, Param, Col_diag, pfull,temp,   &
            mixing_ratio, phalf, Don_conv, donner_humidity_factor, &
            donner_humidity_area,exit_flag,  ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
 
!---------------------------------------------------------------------
!    when strat_cloud is active, call strat_cloud_donner_tend to
!    define increments to cloudice, cloudwater and cloud area associated
!    with deep convective vertical mass flux and detrainment from the
!    mesoscale to the large-scale. 
!---------------------------------------------------------------------
      if (cloud_tracers_present) then
        call don_l_strat_cloud_donner_tend_k   &
             (isize, jsize, nlev_lsm, Param, Col_diag, dmeso_3d,   &
              Don_conv%xliq, Don_conv%xice, qlin, qiin, qain, mhalf_3d, &
              phalf, dql, dqi, dqa, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif

!---------------------------------------------------------------------


end subroutine don_l_lscloud_driver_k



!#####################################################################

subroutine don_l_define_mass_flux_k   &
         (isize, jsize, nlev_lsm, pfull, phalf, Don_conv, dmeso_3d,   &
          mhalf_3d, exit_flag, ermesg, error)

!---------------------------------------------------------------------
!    subroutine define_donner_mass_flux calculates the detrainment rate
!    of the mesoscale to the large-scale and the mass flux associated 
!    with the donner deep convection parameterization.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type

implicit none

!---------------------------------------------------------------------
integer,                                 intent(in)    :: isize, jsize, &
                                                          nlev_lsm
real, dimension(isize,jsize,nlev_lsm),   intent(in)    :: pfull
real, dimension(isize,jsize,nlev_lsm+1), intent(in)    :: phalf
type(donner_conv_type),                  intent(inout) :: Don_conv
real, dimension(isize,jsize,nlev_lsm),   intent(out)   :: dmeso_3d
real, dimension(isize,jsize,nlev_lsm+1), intent(out)   :: mhalf_3d
logical, dimension(isize, jsize),        intent(in)    :: exit_flag
character(len=*),                        intent(out)   :: ermesg
integer,                                 intent(out)   :: error
!----------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!
!   intent(inout) variables:
!     Don_conv
!
!   intent(out) variables:
!
!     dmeso_3d       detrainment rate from convective system 
!                    [ sec**(-1) ]
!     mhalf_3d       mass flux at model half-levels 
!                    [ kg / (m**2 sec) ]
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables:

      real      :: dnna, dnnb    ! weights for averaging full-level 
                                 ! mass fluxes to half levels
                                 ! [ dimensionless ]
      integer   :: i,j,k         ! do-loop indices

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!----------------------------------------------------------------------
!    calculate the detrainment rate from the convective system 
!    (dmeso_3d). dmeso_3d is (1/rho) dM/dz [ sec**(-1) ], where M is 
!    the detraining mass flux of the convective system (kg/(m**2 s)).
!----------------------------------------------------------------------
      do k=1,nlev_lsm
        do j=1,jsize
          do i=1,isize
            if ( .not. exit_flag(i,j) .and. &
                 Don_conv%xice(i,j,k) >= 1.0e-10) then
              dmeso_3d(i,j,k) = Don_conv%emes(i,j,k)/ &
                                Don_conv%xice(i,j,k)
            else
              dmeso_3d(i,j,k) = 0.
            end if
          end do
        end do
      end do

!---------------------------------------------------------------------
!    define the donner parameterization mass flux at model half-levels.
!    include both mesoscale and cell scale fluxes.
!---------------------------------------------------------------------
      do k=1,nlev_lsm-1
        do j=1,jsize
          do i=1,isize
            if (((Don_conv%uceml(i,j,k) <= 1.0e-10)  .and.   &
                (Don_conv%umeml(i,j,k) <= 1.0e-10) .and. &
                (Don_conv%dmeml(i,j,k) <= 1.0e-10)) .or. &
                  exit_flag(i,j)        ) then
              mhalf_3d(i,j,k)   = 0.
              mhalf_3d(i,j,k+1) = 0.
            else
              dnna = phalf(i,j,k+1) - pfull(i,j,k)
              dnnb = pfull(i,j,k+1) - phalf(i,j,k+1)
              mhalf_3d(i,j,k+1) =  &
                                 (dnnb*Don_conv%uceml(i,j,k) + &
                                  dnna*Don_conv%uceml(i,j,k+1) +  &
                                  dnnb*Don_conv%umeml(i,j,k) +   &
                                  dnna*Don_conv%umeml(i,j,k+1) +  &
                                  dnnb*Don_conv%dmeml(i,j,k) +   &
                                  dnna*Don_conv%dmeml(i,j,k+1))/ &
                                  (dnna + dnnb)
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------
!    define the donner parameterization mass fluxes at model top and
!    bottom to be 0.0.
!---------------------------------------------------------------------
      do j=1,jsize
        do i=1,isize
          mhalf_3d(i,j,nlev_lsm+1) = 0.
          mhalf_3d(i,j,1)      = 0.
        end do
      end do


!---------------------------------------------------------------------



end subroutine don_l_define_mass_flux_k 


!######################################################################

subroutine don_l_adjust_tiedtke_inputs_k   &
         (isize, jsize, nlev_lsm, Param, Col_diag, pfull, temp, &
          mixing_ratio, phalf, Don_conv, donner_humidity_factor, &
          donner_humidity_area, exit_flag, ermesg, error) 

!---------------------------------------------------------------------
!    subroutine adjust_tiedtke_inputs calculates the adjustments to 
!    the relative humidity used as a threshold for cloud formation in 
!    the Tiedtke stratiform cloud parameterization (u00), needed as a 
!    consequence of donner_deep_mod being active. see "Tiedtke u00 
!    adjustment" notes, 11/22/02.
!--------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_param_type, &
                             donner_column_diag_type
use sat_vapor_pres_k_mod, only: lookup_es_k

implicit none

!--------------------------------------------------------------------
integer,                     intent(in)     :: isize, jsize, nlev_lsm
type(donner_param_type),     intent(in)     :: Param
type(donner_column_diag_type),                                      &
                             intent(in)     :: Col_diag
logical, dimension(isize,jsize), intent(in) :: exit_flag   
real, dimension(isize,jsize,nlev_lsm),                              &
                             intent(in)     :: pfull, temp, mixing_ratio
real, dimension(isize,jsize,nlev_lsm+1),                            &
                             intent(in)     :: phalf
type(donner_conv_type),      intent(inout)  :: Don_conv
real, dimension(isize,jsize,nlev_lsm),                               &
                             intent(out)    :: donner_humidity_factor, &
                                               donner_humidity_area
character(len=*),            intent(out)    :: ermesg
integer,                     intent(out)    :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field on model half levels [ Pa ]
!     temp           temperature field at model full levels [ deg K ]
!     mixing_ratio   water vapor specific humidity at model full 
!                    levels [ kg(h2o) / kg(air) ]
!
!   intent(inout) variables:
!
!     Don_conv
!
!   intent(out) variables:
!
!     donner_humidity_ratio    
!                    ratio of large-scale specific humidity to the
!                    specific humidity in the environment outside
!                    of the convective system [ dimensionless ]
!     donner_humidity_area    
!                    fractional area of cell plus meso circulation
!                    associated with donner_deep_mod [ fraction ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:
 
     real                       :: acell, pzm, rfun
     integer                    :: i, j, k, n

!--------------------------------------------------------------------
!   local variables:
!
!     qrf           vapor specific humidity at model full levels 
!                   [ kg (h2o) /kg(air) ]
!     acell         fractional area of cells in grid box
!                   [ fraction ]
!     pzm           pressure level 200 hPa above top of mesoscale
!                   downdraft [ Pa ]
!     rfun          ratio of specific humidity to saturation specific
!                   humidity assumed in the region between the top of
!                   the mesoscale downdraft and the cloud base
!                   [ fraction ]
!     qsat          saturation specific humidity [ kg(h2o) / kg(air) ]
!     esat          saturation vapor pressure [ Pa ]
!     i,j,k,n       do-loop indices
!---------------------------------------------------------------------

      ermesg = ' ' ; error = 0

      do k=1,nlev_lsm           
        do j=1,jsize
          do i=1,isize

!---------------------------------------------------------------------
!    the output fields are given non-default values only in grid boxes
!    where donner_deep_mod has produced either cell or mesoscale cloud.
!---------------------------------------------------------------------
            if ( .not. exit_flag(i,j) .and. &
                 (Don_conv%cual(i,j,k)  > 0.0 .or. &
                  Don_conv%ampta1(i,j) > 0.0)) then

!-------------------------------------------------------------------
!    define the pressure at the base of the mesoscale updraft (pzm), 
!-------------------------------------------------------------------
              pzm = Don_conv%pzm_v(i,j) 

!--------------------------------------------------------------------
!    define the area of the convective cells in each grid box (acell).
!    if in the region of mesoscale updraft, the cell area is defined as
!    the difference between the total cloud area (Don_conv%cual) and 
!    the mesoscale area (Don_conv%ampta1). above and below this region 
!    it is defined as the total cloud area, since no mesoscale cloud is
!    present. be sure that the cell area is non-negative.
!--------------------------------------------------------------------
              if ((pfull(i,j,k) <= Don_conv%pzm_v(i,j)) .and.   &
                  (pfull(i,j,k) >= Don_conv%pztm_v(i,j))) then
                acell = Don_conv%cual(i,j,k) - Don_conv%ampta1(i,j)
              else
                acell = Don_conv%cual(i,j,k)
              endif
              acell = MAX (acell, 0.0)

!--------------------------------------------------------------------
!    define the fractional area of the grid box where the humidity
!    is affected by the deep convection in the mesoscale updraft region 
!    (donner_humidity_area). it is assumed to be the sum of the cell and
!    meso cloud areas. previously, below cloud base and above the top 
!    of the mesoscale updraft, there was no moisture-enhanced area due 
!    to deep convection.      
!      this is no longer the case after mods during AM2p9 -- in current
!      formulation there may be cell cloud above top of mesoscale
!      updraft (when cloud top exceeds plzb), and in grid box containing
!      cloud base level.
!--------------------------------------------------------------------
              donner_humidity_area(i,j,k) = Don_conv%cual(i,j,k)

!----------------------------------------------------------------------
!     between cloud base and the base of the mesoscale updraft, the
!     humidity in an area the size of the sum of the cell and meso 
!     cloud areas is assumed to be affected by the deep convection.
!----------------------------------------------------------------------
              if ((pfull(i,j,k) >  Don_conv%pzm_v(i,j)) .and.  &
                  (pfull(i,j,k) <= Don_conv%pb_v(i,j)) ) then
                donner_humidity_area(i,j,k) = Don_conv%cual(i,j,k) +   &
                                              Don_conv%ampta1(i,j)
              endif

!---------------------------------------------------------------------
!      limit the fractional area to be less than 1.0.
!---------------------------------------------------------------------
              donner_humidity_area(i,j,k) =     &
                                  MIN (donner_humidity_area(i,j,k), 1.0)

!---------------------------------------------------------------------
!    define the variable rfun, which serves as the ratio of specific
!    humidity to saturation specific humidity assumed to be present in
!    the part of the grid box with a mesoscale circulation. saturation 
!    is assumed in the cell portion and in the region from the top of 
!    the mesoscale downdraft to the top of the mesoscale updraft. from 
!    top of mesoscale downdraft to cloud base, the relative humidity 
!    goes from 100% to 70%, as a function of relative pressure distance
!    between the two end points.  
!---------------------------------------------------------------------
              if ((pfull(i,j,k) >= Don_conv%pztm_v(i,j)) .and. &
                  (pfull(i,j,k) <= Don_conv%pmd_v(i,j))) then
                rfun = 1.
              else if ((pfull(i,j,k) >= Don_conv%pmd_v(i,j)) .and.  &
                       (pfull(i,j,k) <= Don_conv%pb_v(i,j))) then
                rfun = 1. - 0.3*(Don_conv%pmd_v(i,j) - pfull(i,j,k))/ &
                                (Don_conv%pmd_v(i,j) -   &
                                                     Don_conv%pb_v(i,j))
              else if (pfull(i,j,k) < Don_conv%pztm_v(i,j)) then
                rfun = 0.
              else if (pfull(i,j,k) > Don_conv%pztm_v(i,j)) then
                rfun = 0.
              endif

              donner_humidity_factor(i,j,k) = Don_conv%ampta1(i,j)*rfun

!---------------------------------------------------------------------
!    if in diagnostics column, output the profiles of cell and mesoscale
!    cloud area, the pressure, temperature and specific humidity prof-
!    iles,  the saturation specific humidity, the large-scale and 
!    convective environment specific humidities.
!---------------------------------------------------------------------
              if (Col_diag%in_diagnostics_window) then
                do n=1, Col_diag%ncols_in_window
                  if (j == Col_diag%j_dc(n) .and.    &
                      i == Col_diag%i_dc(n)) then
                    if (pfull(i,j,k) .le. Don_conv%pb_v(i,j)) then
                      if (Don_conv%ampta1(i,j) .ne. 0.) then
                        write (Col_diag%unit_dc(n), '(a, i5, 2f20.10)') &
                           'k,cual,ampt= ', &
                            k, Don_conv%cual(i,j,k), Don_conv%ampta1(i,j)
                        write (Col_diag%unit_dc(n), '(a, i5, 4f20.10)') &
                                     'pb,prinp,trf, acell= ',  k, &
                                   Don_conv%pb_v(i,j), pfull(i,j,k),   &
                                    temp(i,j,k),acell
                        write (Col_diag%unit_dc(n), '(a, i5, f20.10)') &
                               'rfun,= ', k,  rfun
                      endif
                    endif
                  endif
                end do
              endif

!---------------------------------------------------------------------
!    if there is no cell or meso cloud in the grid box, set 
!    donner_humidity_ratio to 1.0 and donner_humidity_area to 0.0
!---------------------------------------------------------------------
            else
              donner_humidity_factor(i,j,k) = 0.0
              donner_humidity_area(i,j,k) = 0.0
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------


end subroutine don_l_adjust_tiedtke_inputs_k 



!#######################################################################

 subroutine don_l_strat_cloud_donner_tend_k    &
          (isize, jsize, nlev_lsm, Param, Col_diag, dmeso_3d, qlmeso, &
           qimeso, qlin, qiin, qain, mhalf_3d, phalf, dql, dqi, dqa,  &
           ermesg, error)

!--------------------------------------------------------------------
!    subroutine strat_cloud_donner is part of the linkage between
!    the deep convection parameterization of donner_deep_mod and the 
!    tiedtke/rotstayn cloud fraction/microphysics parameterization of 
!    strat_cloud_mod.
!-----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_column_diag_type

implicit none

!-----------------------------------------------------------------------
integer,                       intent(in)   :: isize, jsize, nlev_lsm
type(donner_param_type),       intent(in)   :: Param
type(donner_column_diag_type), intent(in)   :: Col_diag
real, dimension(isize,jsize,nlev_lsm),          &
                               intent(in)   :: dmeso_3d, qlmeso, qimeso,&
                                               qlin, qiin, qain
real, dimension(isize,jsize,nlev_lsm+1),         &
                               intent(in)   :: mhalf_3d, phalf
real, dimension(isize,jsize,nlev_lsm),           &
                               intent(out)  :: dql, dqi, dqa
character(len=*),              intent(out)  :: ermesg
integer,                       intent(out)  :: error
!---------------------------------------------------------------------
!   intent(in) variables:
!
!     dmeso_3d       mass detrainment rate from mesoscale region to 
!                    large-scale region [ sec**(-1) ]
!     qlmeso         mesoscale cloud liquid specific humidity   
!                    [ (kg (h2o) /kg air) ]
!     qimeso         mesoscale cloud ice specific humidity 
!                    [ (kg (h2o) /kg air) ]
!     mhalf_3d       total donner-related mass flux = 
!                    mesoscale_mass_flux + convective_mass_flux,
!                    on model interface levels [ kg / ((m**2) s) ]
!                    NOTE: mhalf_3d(:,:,1) and mhalf_3d(:,:,kdim+1) are 
!                    assumed to be 0.0 in this subroutine.
!     phalf          pressure at model half-levels [ Pa ]
!     qlin           large-scale cloud liquid specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qiin           large-scale cloud ice specific humidity
!                    [ kg(h2o) / kg(air) ]
!     qain           large-scale cloud fraction [ fraction ]
!
!   intent(out) variables:
!
!     dql            increment to large-scale cloud liquid field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqi            increment to large-scale cloud ice field from
!                    donner_deep_mod [ kg(h2o) / kg (air) ]
!     dqi            increment to large-scale cloud area field from
!                    donner_deep_mod [ fraction ]
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   local variables:
      real, dimension (isize, jsize,nlev_lsm) ::  mass  
      integer     :: k, n

!---------------------------------------------------------------------
!   local variables:
!
!         mass     mass of air per square meter in the given layer
!                  [ kg / m**2 ]
!         kdim     number of model layers
!         k,n      do-loop indices
!
!---------------------------------------------------------------------

      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    define the number of model layers (kdim) and the mass per unit area
!    contained in each layer (mass).
!---------------------------------------------------------------------
!     kdim = nlev_lsm
      do k=1,nlev_lsm
        mass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/Param%grav 
      end do

!---------------------------------------------------------------------
!    define the large scale cloud increments at level 1 to be 0.0.
!---------------------------------------------------------------------
      dql (:,:,1) = 0.
      dqi (:,:,1) = 0.
      dqa (:,:,1) = 0.
    
!---------------------------------------------------------------------
!    define the tendencies of cloud liquid, cloud ice and cloud area
!    due to the vertical mass flux associated with donner_deep con-
!    vection.
!---------------------------------------------------------------------
      do k=2,nlev_lsm

!---------------------------------------------------------------------
!    define the tendencies due to flux through the top interface of the
!    layer.
!---------------------------------------------------------------------
        dql(:,:,k) = mhalf_3d(:,:,k)*0.5*(qlin(:,:,k) + qlin(:,:,k-1))/ &
                                                           mass(:,:,k)
        dqi(:,:,k) = mhalf_3d(:,:,k)*0.5*(qiin(:,:,k) + qiin(:,:,k-1))/ &
                                                           mass(:,:,k)
        dqa(:,:,k) = mhalf_3d(:,:,k)*0.5*(qain(:,:,k) + qain(:,:,k-1))/ &
                                                           mass(:,:,k)

!---------------------------------------------------------------------
!    add the tendencies due to flux through the bottom interface of 
!    the layer.
!---------------------------------------------------------------------
        dql(:,:,k-1) = dql(:,:,k-1) - mhalf_3d(:,:,k)*0.5*     &
                            (qlin(:,:,k-1) + qlin(:,:,k))/mass(:,:,k-1)
        dqi(:,:,k-1) = dqi(:,:,k-1) - mhalf_3d(:,:,k)*0.5*      &
                            (qiin(:,:,k-1) + qiin(:,:,k))/mass(:,:,k-1)
        dqa(:,:,k-1) = dqa(:,:,k-1) - mhalf_3d(:,:,k)*0.5*       &
                            (qain(:,:,k-1) + qain(:,:,k))/mass(:,:,k-1)
      end do

!---------------------------------------------------------------------
!    add the effects of detrainment from the mesoscale region.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        dql (:,:,k) = dql (:,:,k) + dmeso_3d(:,:,k)*qlmeso(:,:,k)
        dqi (:,:,k) = dqi (:,:,k) + dmeso_3d(:,:,k)*qimeso(:,:,k)
        where (qlmeso(:,:,k) + qimeso(:,:,k) >= 1.e-10) 
          dqa (:,:,k) = dqa (:,:,k) + dmeso_3d(:,:,k)
        end where
      end do

!--------------------------------------------------------------------
!    if in diagnostics window, output values of the cloud variables and
!    the tendencies due to donner_deep convection.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          do k=1, nlev_lsm     
            if (dql(Col_diag%i_dc(n),Col_diag%j_dc(n), k) /= 0.0 .or. &
                dqi(Col_diag%i_dc(n),Col_diag%j_dc(n), k) /= 0.0 ) then 
              write (Col_diag%unit_dc(n), '(a, i5, 3f20.10)') &
                  'donner_deep,strat_cloud_donner_tend', k, &
                   qlin(Col_diag%i_dc(n),Col_diag%j_dc(n), k), &
                   qiin(Col_diag%i_dc(n),Col_diag%j_dc(n), k), &
                   qain(Col_diag%i_dc(n),Col_diag%j_dc(n), k)
              write (Col_diag%unit_dc(n), '(a, i5, 3f20.10)') &
                  'donner_deep,strat_cloud_donner_tend', k, &
                   dql(Col_diag%i_dc(n),Col_diag%j_dc(n), k), &
                   dqi(Col_diag%i_dc(n),Col_diag%j_dc(n), k), &
                   dqa(Col_diag%i_dc(n),Col_diag%j_dc(n), k)
            endif
          end do
        end do
      endif 

!-----------------------------------------------------------------------


end subroutine don_l_strat_cloud_donner_tend_k


!###################################################################



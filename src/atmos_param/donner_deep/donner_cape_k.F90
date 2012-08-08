
!VERSION NUMBER:
!  $Id: donner_cape_k.F90,v 19.0 2012/01/06 20:06:18 fms Exp $

!module donner_cape_inter_mod

!#include "donner_cape_interfaces.h"

!end module donner_cape_inter_mod

!####################################################################

subroutine don_c_def_conv_env_k          &
         (isize, jsize, nlev_lsm, nlev_hires, Nml, Param, Initialized, &
          Col_diag,    &
          temp, mixing_ratio, pfull, lag_cape_temp, lag_cape_vapor,    &
          lag_cape_press, current_displ, cbmf, Don_cape, Don_conv, ermesg, error)

!---------------------------------------------------------------------
!   subroutine don_c_def_conv_env_k manages the 
!   determination of the stability of moist ascent in each model column,
!   calculating various quantities which define the movement of the 
!   parcel in the given column. 
!---------------------------------------------------------------------

use donner_types_mod, only : donner_nml_type, donner_param_type, &
                             donner_column_diag_type, donner_cape_type,&
                             donner_initialized_type, &
                             donner_conv_type
implicit none

!----------------------------------------------------------------------
integer,                       intent(in)    :: isize, jsize, nlev_lsm, &
                                                nlev_hires
type(donner_nml_type),         intent(in)    :: Nml      
type(donner_param_type),       intent(in)    :: Param
type(donner_initialized_type),       intent(in)    :: Initialized
type(donner_column_diag_type), intent(in)    :: Col_diag
real,    dimension(isize,jsize,nlev_lsm),                    &
                               intent(in)    :: temp, mixing_ratio,  &
                                                pfull, lag_cape_temp, &
                                                lag_cape_vapor, &
                                                lag_cape_press
real,    dimension(isize,jsize),                             &
                               intent(in)    :: current_displ, cbmf 
type(donner_cape_type),        intent(inout) :: Don_cape
type(donner_conv_type),        intent(inout) :: Don_conv
character(len=*),              intent(out)   :: ermesg
integer,                       intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     lag_cape_temp  temperature field used in lag-time cape 
!                    calculation [ deg K ]
!     lag_cape_vapor vapor mixing ratio field used in lag-time
!                    cape calculation [ kg(h2o) / kg(dry air) ]
!     lag_cape_press model full-level pressure field used in 
!                    lag-time cape calculation  [ Pa ]
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     cbmf
!
!   intent(out) variables:
!
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape type derived type variable containing 
!                    diagnostics and intermediate results related to 
!                    the cape calculation associated with the donner 
!                    convection parameterization
!     Don_conv       donner_conv_type derived type variable containing 
!                    diagnostics and intermediate results describing 
!                    the nature of the convection produced by the 
!                    donner parameterization
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real,    dimension (isize, jsize, nlev_lsm) :: mid_cape_temp, &
                                                     mid_cape_vapor
      logical, dimension (isize, jsize)           :: no_convection
      integer                                     :: lowest_bl_index
      integer                                     :: i, j, k, n

!----------------------------------------------------------------------
!   local variables:
!     
!    mid_cape_temp       temperature field to be used in the cape 
!                        calculation using current time level values
!                        [ deg K ]
!    mid_cape_vapor      vapor mixing ratio field to be used in the cape
!                        calculation using current time level values
!                        [ kg(h2o) / kg (dry air) ]
!    no_convection       logical indicating columns in which cape calc-
!                        ulation may be skipped because low-level parcel
!                        displacement is downward 
!    lowest_bl_index     k index of topmost level within the surface 
!                        boundary layer (no local time variation of the
!                        temperature and moisture profiles is allowed 
!                        within the sbl) (index 1 at top of model) 
!    i, j, k, n          do-loop indices
!                        
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
       ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    define the vertical index of the highest level within the surface
!    boundary layer (lowest_bl_index).
!---------------------------------------------------------------------
      lowest_bl_index = nlev_lsm - Nml%model_levels_in_sfcbl + 1

!--------------------------------------------------------------------
!    if in diagnostics window, write message indicating lag-time cape
!    calculation is being done.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write (Col_diag%unit_dc(n), ' (//, a)')  &
              '               CAPE calculation for LAG time profile'
        end do
      endif
        
!--------------------------------------------------------------------
!    define a flag array which will be .true. when deep convection is
!    precluded in the column due to downward motion at the lowest model
!    level at the current time (no_convection).
!--------------------------------------------------------------------
      if (Nml%use_llift_criteria) then
        no_convection(:,:) = (current_displ(:,:) >= 0.0)
      else
        no_convection(:,:) = (current_displ(:,:) >  0.0)
      end if
      do j=1,jsize
        do i=1,isize
          if (Initialized%using_unified_closure .and.    &
                                             cbmf(i,j) == 0.) then
            no_convection(i,j) = .true.
          endif
        end do
      end do

!---------------------------------------------------------------------
!    call don_c_cape_calculation_driver_k with the lag-time-based 
!    profiles to define the lcl, the moist adiabat, convective inhibition
!    and cape for a parcel displaced upwards from the lowest model level.
!---------------------------------------------------------------------
      call don_c_cape_calculation_driver_k  &
           (isize, jsize, nlev_lsm, nlev_hires, Col_diag, Param, Nml, &
            lag_cape_temp, lag_cape_vapor, lag_cape_press,  &
            no_convection, Don_cape, ermesg, error)
 
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call don_c_cape_diagnostics_k to output the cape calculation 
!    diagnostics.
!---------------------------------------------------------------------
      call don_c_cape_diagnostics_k   &
           (isize, jsize, nlev_lsm, nlev_hires, Col_diag,  &
            Don_cape, no_convection, ermesg, error)
 
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return
  
!--------------------------------------------------------------------
!    save the values of cape and column integrated water vapor returned
!    from this calculation as the lag-time values of these quantities.
!    the convection parameterization requires the time-tendency of both
!    of these quantities.
!--------------------------------------------------------------------
      Don_cape%qint_lag (:,:) = Don_cape%qint(:,:)
      Don_cape%xcape_lag(:,:) = Don_cape%xcape(:,:)

!--------------------------------------------------------------------
!    if this is a diagnostics window, output a message indicating that 
!    the current-time cape calculation is beginning.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write (Col_diag%unit_dc(n), ' (//,  a)')  &
              '               CAPE calculation for CURRENT time profile'
        end do
      endif     

!--------------------------------------------------------------------
!    define the temperature and mixing ratio fields to be used in de-
!    termining the current-time cape sounding values at all levels 
!    above the surface boundary layer. they are the model values as 
!    input from the calling routine.
!--------------------------------------------------------------------
      do k=1,lowest_bl_index-1         
        mid_cape_temp(:,:,k)  = temp(:,:,k)
        mid_cape_vapor(:,:,k) = mixing_ratio(:,:,k)
      end do

!--------------------------------------------------------------------
!    in the surface boundary layer, use the temperature and water vapor
!    values used in the lag-time calculation of cape (the catendb clo-
!    sure). this prevents the production of cape time tendencies  
!    resulting from temporal noise from the surface boundary layer. 
!--------------------------------------------------------------------
      do k=lowest_bl_index, nlev_lsm
        mid_cape_temp(:,:,k) = lag_cape_temp(:,:,k)
        mid_cape_vapor(:,:,k) = lag_cape_vapor(:,:,k)
      end do

!---------------------------------------------------------------------
!    call don_c_cape_calculation_driver_k with the current-time-
!    based profiles to define the lcl, the moist adiabat, convective 
!    inhibition and cape for a parcel displaced upwards from the lowest 
!    model level.
!---------------------------------------------------------------------
      call don_c_cape_calculation_driver_k  &
           (isize, jsize, nlev_lsm, nlev_hires, Col_diag, Param, Nml, &
            mid_cape_temp, mid_cape_vapor, pfull,  &
            no_convection, Don_cape, ermesg, error)
 
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call don_c_cape_diagnostics_k to output the cape calculation 
!    diagnostics.
!---------------------------------------------------------------------
      call don_c_cape_diagnostics_k   &
           (isize, jsize, nlev_lsm, nlev_hires, Col_diag,  &
            Don_cape, no_convection, ermesg, error)
 
!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------


end subroutine don_c_def_conv_env_k



!#######################################################################

subroutine don_c_cape_calculation_driver_k  &
         (isize, jsize, nlev_lsm, nlev_hires, Col_diag, Param, Nml, &
          temperature, mixing_ratio, pfull, no_convection, Don_cape, &
          ermesg, error)

!--------------------------------------------------------------------
!    subroutine cape_calculation_driver defines high-resolution atmos-
!    pheric temperature and vapor mixing ratio profiles equivalent to 
!    those in the large-scale model, and determines the parameters 
!    defining parcel movement within this environment.
!--------------------------------------------------------------------

use donner_types_mod,only : donner_column_diag_type, donner_param_type,&
                            donner_cape_type, donner_nml_type

implicit none

!--------------------------------------------------------------------
integer,                               intent(in)    :: isize, jsize,  &
                                                        nlev_lsm, &
                                                        nlev_hires
type(donner_column_diag_type),         intent(in)    :: Col_diag
type(donner_param_type),               intent(in)    :: Param
type(donner_nml_type),                 intent(in)    :: Nml  
real, dimension(isize,jsize,nlev_lsm), intent(in)    :: temperature,  &
                                                        mixing_ratio, &
                                                        pfull
logical, dimension(isize,jsize),       intent(in)    :: no_convection
type(donner_cape_type),                intent(inout) ::  Don_cape
character(len=*),                      intent(out)   :: ermesg
integer,                               intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      temperature      temperature profile on large-scale model grid 
!                       that is to be used for the cape calculation 
!                       [ deg K ]
!      mixing_ratio     water vapor mixing ratio profile on the large-
!                       scale model grid that is to be used in the cape 
!                       calculation [ kg(h2o) / kg(dry air) ]
!      pfull            large-scale model full-level pressure profile 
!                       [ Pa ]
!      no_convection    logical indicating if convection is precluded 
!                       from column because of lowest-level downward 
!                       motion
!
!   intent(inout) variables:
!
!     Don_cape          donner_cape type derived type variable contain-
!                       ing diagnostics and intermediate results related
!                       to the cape calculation associated with the 
!                       donner convection parameterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

       logical :: debug_ijt     ! logical indicating if diagnostics are
                                ! desired in current column
       integer  :: diag_unit    ! unit number for diagnostic output for
                                ! current column
       integer  :: i, j, n      ! do-loop indices

      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    loop over columns in physics window.
!------------------------------------------------------------------
      do j=1,jsize
        do i=1,isize

!---------------------------------------------------------------------
!    if deep convection is not precluded in this column, continue the
!    calculation; if deep convection is precluded in this column, 
!    the variables contained in Don_cape will retain the values prev-
!    iously given them on initialization.
!---------------------------------------------------------------------
          if (.not. no_convection(i,j)) then

!---------------------------------------------------------------------
!    define a logical variable indicating if diagnostics are to be
!    produced for this column. if so, define the unit number for the
!    output file (diag_unit).
!---------------------------------------------------------------------
            debug_ijt = .false.
            if (Col_diag%in_diagnostics_window) then
               do n=1,Col_diag%ncols_in_window
                 if (j == Col_diag%j_dc(n) .and.     &
                     i == Col_diag%i_dc(n)) then
                   debug_ijt = .true.
                   diag_unit = Col_diag%unit_dc(n)
                   exit
                 endif
               end do
            endif

!--------------------------------------------------------------------
!    call generate_cape_sounding to produce a high-resolution atmos-
!    pheric sounding to be used to evaluate cape.
!--------------------------------------------------------------------
            call don_c_generate_cape_sounding_k &
                 (nlev_lsm, nlev_hires, temperature(i,j,:),   &
                  mixing_ratio(i,j,:), pfull(i,j,:),   &
                  Don_cape%model_t(i,j,:), Don_cape%model_r(i,j,:), &
                  Don_cape%model_p(i,j,:), Don_cape%cape_p(i,j,:), &
                  Don_cape%env_t(i,j,:), Don_cape%env_r(i,j,:), ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return

!--------------------------------------------------------------------
!    call displace_parcel to calculate the behavior of a parcel moving
!    upwards in the model column defined by Don_cape%cape_p, 
!    Don_cape%env_t and Don_cape%env_r.
!--------------------------------------------------------------------
            call don_c_displace_parcel_k   &
                 (nlev_hires, diag_unit, debug_ijt, Param,    &
                  Nml%do_freezing_for_cape, Nml%tfre_for_cape, &
                  Nml%dfre_for_cape, Nml%rmuz_for_cape, &
!      the following value of .true. corresponds to the dummy argument 
!      use_constant_rmuz. it currently must be .true. for this call to  
!      displace_parcel; it may be false when displace_parcel is called
!      for a closure calculation.
                  .true., &  ! (use_constant_rmuz)
!      the following value corresponds to dummy argument 
!      carry_condensate.  it currently must be set to .false. for cape
!      calculations used to determine the presence of convection; it 
!      may be set true in calls to displace_parcel used for closure
!      calculations.
                   .false., &  ! (carry_condensate)
                 Nml%closure_plume_condensate, & ! (is not used when 
                                                 ! carry_condensate
                                                 ! is .false.)
                  Don_cape%env_t(i,j,:), Don_cape%env_r(i,j,:), &
                  Don_cape%cape_p(i,j,:), .true.,  &
                  Don_cape%plfc(i,j), Don_cape%plzb(i,j),  &
                  Don_cape%plcl(i,j), Don_Cape%coin(i,j),   &
                  Don_cape%xcape(i,j), Don_cape%parcel_r(i,j,:), &
                  Don_cape%parcel_t(i,j,:), ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return

!---------------------------------------------------------------------
!   call integrate_vapor to produce the column integral of water vapor.
!--------------------------------------------------------------------
            call don_c_integrate_vapor_k  &
                 (nlev_hires, Param, Don_cape%env_r(i,j,:),  &
                  Don_cape%cape_p(i,j,:), Don_cape%qint(i,j), ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return
          endif
        end do
      end do

!--------------------------------------------------------------------



end subroutine don_c_cape_calculation_driver_k 


!###################################################################

subroutine don_c_displace_parcel_k   &
         (nlev_hires, diag_unit, debug_ijt, Param, do_freezing, &
          tfreezing, dfreezing, rmuz, use_constant_rmuz, &
          carry_condensate, condensate_carried, env_t, env_r,  &
          cape_p, coin_present, plfc, plzb, plcl, coin, xcape,       &
          parcel_r, parcel_t, ermesg, error)

!----------------------------------------------------------------------
!    displace_parcel moves a parcel upwards from the lowest model level
!    (istart) in an environment defined by env_t, env_r and cape_p, 
!    determining the critical transition levels during its ascent (plfc,
!    plzb, plcl), its temperature (parcel_t) and mixing ratio (parcel_r)
!    at each pressure level (cape_p) during its ascent, and if desired,
!    calculating the associated energy integrals for the parcel 
!    (xcape, coin).
!---------------------------------------------------------------------

use donner_types_mod, only:  donner_param_type 
implicit none

!---------------------------------------------------------------------
integer,                       intent(in)  :: nlev_hires
integer,                       intent(in)  :: diag_unit
logical,                       intent(in)  :: debug_ijt
type(donner_param_type),       intent(in)  :: Param
logical,                       intent(in)  :: do_freezing, &
                                              use_constant_rmuz, &
                                              carry_condensate
real,                          intent(in)  :: tfreezing, dfreezing, &
                                              rmuz, condensate_carried
real,   dimension(nlev_hires), intent(in)  :: env_t, env_r, cape_p
logical,                       intent(in)  :: coin_present
real,                          intent(out) :: plfc, plzb, plcl
real,                          intent(out) :: coin, xcape
real,   dimension(nlev_hires), intent(out) :: parcel_r, parcel_t   
character(len=*),              intent(out) :: ermesg
integer,                       intent(out) :: error

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      debug_ijt    logical indicating if diagnostics are desired for
!                   current column
!      diag_unit    i/o unit number for diagnostics file
!      env_t        envirionmental temperature field, index 1 nearest
!                   the surface [ deg K ]
!      env_r        envirionmental vapor mixing ratio field, index 1 
!                   nearest the surface [ kg(h2o) / kg(dry air) ]
!      cape_p       pressure profile, index 1 nearest the surface  
!                   [ Pa ]
!
!   intent(out) variables:
!
!      plfc         pressure at level of free convection [ Pa ]
!      plzb         pressure at level of zero buoyancy [ Pa ]
!      plcl         pressure at lifting condensation level [ Pa ]
!      parcel_r     vapor mixing ratio in parcel. it is set to the
!                   environmental value below level ISTART. index 1 
!                   is level nearest earth's surface.
!      parcel_t     temperature in parcel. it is set to the environ-
!                   mental value below level ISTART. index 1 is level 
!                   nearest earth's surface.
!
!   intent(out), optional:
!
!      coin         convective inhibition -- energy required to lift 
!                   parcel from level ISTART to the level of free 
!                   convection. if parcel becomes buoyant below lcl, 
!                   coin can be < 0.  [ Joules / kg ]
!      xcape        convective available potential energy -- energy 
!                   released as parcel moves from level of free 
!                   convection to level of zero buoyancy [ Joules / kg ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_hires)     ::     env_tv, parcel_tv, dtdp, &
                                             delt, sum_xcape, &
                                             rc, fact1, fact2, fact3
      integer ::  klcl, klfc, klzb
      real    ::  tlcl, rlcl
      logical ::  cape_exit
      integer  :: k

      ermesg = ' ' ; error = 0


!--------------------------------------------------------------------
!  define parcel departure point values. convert mixing ratio to 
!  specific humidity.
!--------------------------------------------------------------------
      parcel_t(1:Param%istart) = env_t(1:Param%istart)
      parcel_r(1:Param%istart) = env_r(1:Param%istart)

      call don_c_calculate_lcl_k   &
           (nlev_hires, Param,             &
             parcel_t(Param%istart),   cape_p,   &
                          env_r, env_t, parcel_r, parcel_t, plcl,  &
                          tlcl, rlcl, klcl, cape_exit, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      if (.not. cape_exit) then

!--------------------------------------------------------------------
!   if in debug mode, print out info on lcl in debug column.
!--------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, f19.10,i4, f20.14, e20.12)')  &
                           'in cape: plcl,klcl,tlcl,rlcl= ',   &
                             plcl    , klcl, tlcl, rlcl       
          write (diag_unit      , '(a, f19.10)')   &
                              'in cape: p(klcl)= ',  cape_p(klcl)
        endif
      else

!---------------------------------------------------------------------
!   if lcl not found, stop calculations in this column.
!---------------------------------------------------------------------
        plfc  = 0.0     
        plzb  = 0.0  
        coin  = 0.0
        xcape = 0.0
        return
      endif

!-------------------------------------------------------------------
!   calculate temperature along saturated adiabat, starting at p(klcl)
!   and a temperature tp to find the level of free convection and
!   the level of zero buoyancy. 
!--------------------------------------------------------------------
      call don_c_define_moist_adiabat_k  &
           (nlev_hires, klcl, Param, parcel_t(klcl), cape_p, env_r, &
            env_t, do_freezing, tfreezing, dfreezing, rmuz,  &
            use_constant_rmuz, carry_condensate, condensate_carried, &
            parcel_t, parcel_r, plfc, plzb, klfc, klzb,  &
            parcel_tv, env_tv, dtdp, rc, fact1, fact2, fact3, cape_exit,&
            ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return


      if (debug_ijt) then
        do k=klcl, klzb+1
          if (.not. (cape_exit)) then
            if  (k == klzb+1) then
            else
              write (diag_unit, '(a, i4, 2f20.14)')  &
                        'in cape: k,tv,tve= ',k,parcel_tv(k), env_tv(k)          
            endif
          else
            write (diag_unit, '(a, i4, 2f20.14)')  &
                        'in cape: k,tv,tve= ',k,parcel_tv(k), env_tv(k)          
          endif
          if (.not. cape_exit) then
            if (k == klzb ) then
              write (diag_unit, '(a, i4, 2f19.10)')  &
                                  'in cape: klzb,plzb,p(klzb)= ',  &
                                   klzb, plzb, cape_p(klzb)
            endif
          endif
          if (k /= klzb+1 .or. cape_exit) then
            if ( k == klzb .and. .not. cape_exit) then 
            else
              write (diag_unit, '(a, 3f17.10)')  &
                   'in cape: fact1,fact2,rc= ',fact1(k), fact2(k),rc(k)
              write (diag_unit, '(a, 2f17.10)') &
                           'in cape: fact1,fact3= ',fact1(k),fact3(k)
              write (diag_unit, '(a, f17.10)')  &
                            'in cape: dtdp= ',dtdp(k)
              write (diag_unit, '(a,  2f20.14)') &
                             'in cape: tc,t= ',parcel_t(k+1), env_t(k+1)
              write (diag_unit, '(a, f19.10, 2e20.12)') &
                   'in cape: p,r,rs= ',cape_p(k+1), env_r(k+1),   &
                                                         parcel_r(k+1)  
            endif
          endif
        end do
      endif
      if (cape_exit)  then
        coin  = 0.0
        xcape = 0.0
        return
      endif

!--------------------------------------------------------------------
!   if this was a call to calculate perturbed profile, bypass cape and
!   cin calculation, since only the tpc and rpc profiles are needed.
!--------------------------------------------------------------------
      if (coin_present) then

!-------------------------------------------------------------------
!   calculate convective inhibition.
!--------------------------------------------------------------------
        call don_c_calculate_cin_k   &
             (nlev_hires, diag_unit, debug_ijt, Param,  cape_p,  &
              parcel_r, env_r, parcel_t, env_t, plfc, cape_exit,   &
              coin, env_tv, parcel_tv, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

        if (.not. cape_exit) then
          if (debug_ijt) then
            do k=Param%istart,nlev_hires  
              if (parcel_tv(k) /= 0.0) then
                write (diag_unit, '(a, i4, 2f20.14)')  &
                            'in cape: k,tvc,tve= ', k,  &
                         parcel_tv(k), env_tv(k)
              endif
            end do
          endif
        else
          xcape = 0.
          return
        endif

!-------------------------------------------------------------------
!   if desired, print out lfc k index and pressure.
!-------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, f19.10)')  &
                    'in cape: klfc, p(klfc)= ', klfc, cape_p(klfc)
        endif

        call don_c_calculate_cape_k   &
             (nlev_hires, klfc, Param, plzb, cape_p, parcel_r, env_r,   &
              parcel_t, env_t, xcape, delt, sum_xcape, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

!---------------------------------------------------------------------
!   print out cape and cape contribution from this level.
!---------------------------------------------------------------------
        if (debug_ijt) then
          do k=1, size(parcel_t(:))
            if (delt(k) /= 0.0) then
              write (diag_unit, '(a,i4, 2f12.8)')  &
                    'in cape: k,delt,xcape= ',k, delt(k), sum_xcape(k)  
            endif
          end do
        endif

!--------------------------------------------------------------------
!  print out diagnostics (cape, cin, tot), if desired.
!--------------------------------------------------------------------
        if (debug_ijt) then
          if (coin /= 0.0 .or. &
              xcape /= 0.0) then
            write (diag_unit      , '(a, f12.6, a)')  &
                      'in cape: cin= ',coin                ,' J/kg'
            write (diag_unit      , '(a, f12.6, a)')  &
                      'in cape: xcape= ',xcape    ,  ' J/kg'
            write (diag_unit      , '(a, f12.6, a)')  &
                       'in cape: tot= ',xcape - coin     ,' J/kg'
          endif
        endif

!--------------------------------------------------------------------
!  check for error in cape calculation. stop execution if present.
!--------------------------------------------------------------------
        if (xcape        .lt. 0.) then            
          ermesg =  ' xcape error -- value < 0.0 '
          error  = 1
        endif
      endif  ! (present(coin))

!---------------------------------------------------------------------


end subroutine don_c_displace_parcel_k

!#####################################################################

subroutine don_c_define_moist_adiabat_k  &
         (nlev_hires, klcl, Param, starting_temp, press, env_r,  &
          env_t, do_freezing, tfreezing, dfreezing, rmuz_constant, &
          use_constant_rmuz, carry_condensate, condensate_carried, &
          parcel_t, parcel_r, plfc, plzb, klfc, &
          klzb, parcel_tv, env_tv, dtdp, rc, fact1, fact2, fact3, &
          cape_exit, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type 
use sat_vapor_pres_k_mod, only: compute_mrs_k

implicit none 

!----------------------------------------------------------------------
integer,                      intent(in)    :: nlev_hires, klcl
type(donner_param_type),      intent(in)    :: Param
real,                         intent(in)    :: starting_temp
real, dimension(nlev_hires),  intent(in)    :: press, env_r, env_t
logical,                      intent(in)    :: do_freezing,  &
                                               use_constant_rmuz, &
                                               carry_condensate
real,                         intent(in)    :: tfreezing, dfreezing, &
                                               rmuz_constant, &
                                               condensate_carried
real, dimension(nlev_hires),  intent(inout) :: parcel_t, parcel_r     
real,                         intent(out)   :: plfc, plzb 
integer,                      intent(out)   :: klfc, klzb
real, dimension(nlev_hires),  intent(out)   :: parcel_tv, env_tv, dtdp, &
                                               rc, fact1, fact2, fact3
logical,                      intent(out)   :: cape_exit
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error

      real     :: es_v_s, qe_v_s, rs_v_s, qs_v_s, pb, tp_s
      real,dimension(nlev_hires)     :: rmuz, z, fact7
      real     :: dz
      real     :: hlvls
      logical  :: capepos_s 
      integer  :: ieqv_s
      integer  :: k, nbad
      logical  :: not_all_frozen    
      real     :: cumulative_freezing, prev_cufr
      real     :: available_cd, r_at_cb, condensate_to_freeze
      logical  :: first_freezing_level

      ermesg = ' ' ; error = 0

      plfc = 0.0     
      plzb = 0.0  
      klfc = nlev_hires - 1
      klzb = nlev_hires - 1
      capepos_s = .false.
      cape_exit = .false.
      tp_s = starting_temp
      z = 0.

      first_freezing_level = .true.
      fact7 = 0.
      not_all_frozen = .true.
      prev_cufr = 0.
      if (carry_condensate .and. do_freezing) then
        do k=klcl,nlev_hires-1
          if (env_t(k) < tfreezing .and. not_all_frozen) then
            cumulative_freezing = (tfreezing-env_t(k))/dfreezing
            if (cumulative_freezing >= 1.0) then
              cumulative_freezing = 1.0         
              not_all_frozen = .false.
            endif
            fact7(k) = cumulative_freezing - prev_cufr  
            prev_cufr = cumulative_freezing
            if (.not. not_all_frozen) exit
          endif
        end do
      endif
          
!-------------------------------------------------------------------
!    calculate temperature along saturated adiabat, starting at p(klcl)
!    and a temperature tp to find the level of free convection and
!    the level of zero buoyancy. 
!--------------------------------------------------------------------
      do k=klcl,nlev_hires-1

!--------------------------------------------------------------------
!    if pressure has gone below the minimum at which deep convection 
!    is allowed, set flag to end calculation in this column.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    define saturation vapor pressure for the parcel.
!--------------------------------------------------------------------
        call compute_mrs_k (tp_s, press(k), Param%D622, Param%D608, &
                            rs_v_s, nbad,  esat = es_v_s)

!---------------------------------------------------------------------
!    save the cloud base r so that the condensate available when 
!    reaching the level where freezing begins may be determined. 
!---------------------------------------------------------------------
        if (k == klcl) then
          r_at_cb = rs_v_s
        endif
        if ( first_freezing_level) then
          available_cd = r_at_cb - rs_v_s
        endif
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (nbad /= 0) then
          ermesg = 'subroutine don_c_define_moist_adiabat_k: Temperatures out of range of esat table'
          error  = 1
          return
        endif

!--------------------------------------------------------------------
!    define the environmental and parcel virtual temperature and specific
!    humidity.
!--------------------------------------------------------------------
        qe_v_s = env_r(k)/(1. + env_r(k))
        env_tv(k) = env_t(k)*(1.+ Param%D608*qe_v_s   )
        qs_v_s = rs_v_s/(1. + rs_v_s)
        parcel_tv(k) = tp_s*(1. + Param%D608*qs_v_s   )

!--------------------------------------------------------------------
!    determine whether the parcel temperature is cooler or warmer than 
!    the environment.
!--------------------------------------------------------------------
        call don_u_numbers_are_equal_k  &
             (parcel_tv(k), env_tv(k), ermesg, error, ieqv_s)
   
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

!---------------------------------------------------------------------
!    integrate parcel upward, finding level of free convection and 
!    level of zero buoyancy.
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    determine if the level of free convection has been reached. 
!-------------------------------------------------------------------
        if ((ieqv_s >= 0) .and. (.not. capepos_s)) then
          capepos_s    = .true.
          plfc = press(k)  
          klfc = k
        endif

!-------------------------------------------------------------------
!   determine if the level of zero buoyancy has been reached.  if so,
!   set flag so that calculation will be ended in this column.
!-------------------------------------------------------------------
        if ((ieqv_s < 0) .and. (capepos_s)) then
          klzb = k
          plzb = (press(k) + press(k-1))/2.
          parcel_t(k+1:nlev_hires) = env_t(k+1:nlev_hires)
          parcel_r(k+1:nlev_hires) = env_r(k+1:nlev_hires)
          exit

!---------------------------------------------------------------------
!   if not, continue moving parcel up pseudo-adiabat to next cape-
!   calculation pressure level. define new parcel temperature and
!   mixing ratio at this level; if temperature is colder than allowed,
!   end integration.
!-------------------------------------------------------------------
        else   !  (cape is pos, parcel warmer than env)
          klzb    = 0
          if (do_freezing) then
            if (tp_s .gt. tfreezing) then
              hlvls = Param%hlv
            else if (tp_s .le. (tfreezing - dfreezing)) then
              hlvls = Param%hls
            else
              hlvls = (tfreezing - tp_s)*Param%hls + &
                      (tp_s - (tfreezing -dfreezing))*Param%hlv
              hlvls = hlvls/dfreezing 
            endif
          else
            hlvls = Param%hlv
          endif
          rc(k) = (1. - qs_v_s)*Param%rdgas + qs_v_s*Param%rvgas
          pb = 0.5*(press(k) + press(k+1))

!---------------------------------------------------------------------
!     define the entrainment coefficient (rmuz) [ m (-1) ] .
!     z is the height above cloud base [ m ] .
!     dz is the height increment for the current layer.
!---------------------------------------------------------------------
          if (use_constant_rmuz) then
            rmuz(k) = rmuz_constant
          else
            if (k == klcl)  then
              rmuz(klcl) = rmuz_constant
              z(klcl) = 0.  
            else
              dz = - alog((press(k)/press(k-1)))*Param%rdgas*  &
                    parcel_tv(k)/Param%grav
              z(k) = z(k-1) + dz
              rmuz(k) = rmuz_constant/( 1.0 + rmuz_constant*z(k))
            endif
          endif

          fact1(k) = Param%rdgas/Param%cp_air
          fact2(k) = parcel_tv(k) + (hlvls*rs_v_s/rc(k))
          fact1(k) = fact1(k)*fact2(k) +                  &
                     Param%rdgas*env_tv(k)*rmuz(k)*(tp_s-env_t(k)   &
                    + (hlvls*(rs_v_s-env_r(k))/Param%cp_air))/Param%grav
          fact3(k) = Param%d622*(hlvls**2)*es_v_s/    &
                     (Param%cp_air*pb*Param%rvgas*(parcel_tv(k)**2))
          fact3(k) = 1. + fact3(k)

!---------------------------------------------------------------------
!    calculate the term associated with the freezing of condensate 
!    carried along in the plume (fact7).  condensate_to_freeze is the
!    amount of condensate present at the level where freezing begins.
!    it is proportionately frozen over the specified range of freezing. 
!    term may be optionally included when cape is 
!    computed for cumulus closure determination.
!---------------------------------------------------------------------
          if (fact7(k) /= 0.0 .and. first_freezing_level) then
            condensate_to_freeze =    &
                                 MIN(available_cd, condensate_carried)
            first_freezing_level = .false.
          endif
          fact7(k) = Param%hlf*fact7(k)*condensate_to_freeze/  &
                    (Param%cp_air*alog(press(k+1)/press(k)))
          dtdp(k) = (fact1(k) + fact7(k))/fact3(k)
          tp_s    = tp_s + dtdp(k)*alog(press(k+1)/press(k))
          if (tp_s < Param%tmin)  then
            cape_exit = .true.
            parcel_t(k+1) = tp_s
            parcel_r(k+1) = rs_v_s
            parcel_t(k+2:nlev_hires) = env_t(k+2:nlev_hires)
            parcel_r(k+2:nlev_hires) = env_r(k+2:nlev_hires)
            exit  ! exit k loop
          else
            parcel_t(k+1) = tp_s   
            parcel_r(k+1) = rs_v_s   
          endif
        endif   !  (ieq < 0, capepos)
      end do   ! k loop


!-------------------------------------------------------------------


end subroutine don_c_define_moist_adiabat_k



!####################################################################

subroutine don_c_calculate_cin_k    &
         (nlev_hires, diag_unit, debug_ijt, Param, press, parcel_r,   &
          env_r, parcel_t, env_t, plfc, cape_exit, coin,   &
          env_tv, parcel_tv, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                     intent(in)  :: nlev_hires
integer,                     intent(in)  :: diag_unit
logical,                     intent(in)  :: debug_ijt            
type(donner_param_type),     intent(in)  :: Param
real, dimension(nlev_hires), intent(in)  :: press, parcel_r, env_r,   &
                                            parcel_t, env_t
real,                        intent(in)  :: plfc
logical,                     intent(out) :: cape_exit
real,                        intent(out) :: coin
real,dimension(nlev_hires),  intent(out) :: env_tv, parcel_tv
character(len=*),            intent(out) :: ermesg
integer,                     intent(out) :: error

      real      ::   rbc, rbe, qc, qe, tvc_v_s, tve_v_s, delt
      integer   ::   ieqv_s
      integer   ::   k

!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0
      coin = 0.
      cape_exit = .false.

      parcel_tv(1:Param%istart) = 0.
      env_tv(1:Param%istart) = 0.

!-------------------------------------------------------------------
!   calculate convective inhibition.
!--------------------------------------------------------------------
      do k=Param%istart,nlev_hires-1

!------------------------------------------------------------------
!   determine if sounding fails to produce a level of free convection.
!   if so, set flag to avoid cape calculation. If desired, print out
!   columns where lcl exists, but no lfc. No cin exists above the lfc,
!   and the lfc if it exists must be below UPPER_LIMIT_FOR_LFC.
!------------------------------------------------------------------
        if (press(k+1) <= Param%upper_limit_for_lfc) then
          cape_exit    = .true.
          exit   ! exit k loop
        endif

!--------------------------------------------------------------------
!    define the specific humidity and virtual temperature of the
!    parcel and environment.
!--------------------------------------------------------------------
        rbc = (parcel_r(k) + parcel_r(k+1))/2.
        rbe = (env_r(k) + env_r(k+1))/2.
        qc = rbc/(1. + rbc)
        qe = rbe/(1. + rbe)
        tvc_v_s  = (parcel_t(k) + parcel_t(k+1))/2.
        tve_v_s  = (env_t(k) + env_t(k+1))/2.
        tvc_v_s  = tvc_v_s*(1. + Param%d608*qc)
        tve_v_s  = tve_v_s*(1. + Param%d608*qe)

        env_tv(k) = tve_v_s
        parcel_tv(k) = tvc_v_s

!---------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
        call don_u_numbers_are_equal_k   &
             (tvc_v_s, tve_v_s, ermesg, error, ieqv_s )
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

!---------------------------------------------------------------------
!   add the contribution to cin from this pressure layer.
!---------------------------------------------------------------------
        if ((ieqv_s < 0) .or.      &
            (press(k) > plfc))  then
          delt  = Param%rdgas*(tvc_v_s - tve_v_s)*   &
                      alog(press(k)/press(k+1))
          coin = coin - delt   
        else

!------------------------------------------------------------------
!   determine if sounding fails to produce a level of free convection.
!   if so, set flag to avoid cape calculation. If desired, print out
!   columns where lcl exists, but no lfc.
!------------------------------------------------------------------
          parcel_tv(k+1:nlev_hires) = 0.0
          env_tv(k+1:nlev_hires) = 0.0
          exit
        endif
      end do  ! k loop

!--------------------------------------------------------------------


end subroutine don_c_calculate_cin_k



!######################################################################

subroutine don_c_calculate_cape_k  &
         (nlev_hires, klfc, Param,plzb, press, parcel_r, env_r, &
          parcel_t, env_t, xcape, delt, sum_xcape, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type

implicit none 

!----------------------------------------------------------------------
integer,                     intent(in)  :: nlev_hires, klfc
type(donner_param_type),     intent(in)  :: Param
real,                        intent(in)  :: plzb
real, dimension(nlev_hires), intent(in)  :: press, parcel_r, env_r,   &
                                            parcel_t, env_t
real,                        intent(out) :: xcape
real, dimension(nlev_hires), intent(out) :: delt, sum_xcape      
character(len=*),            intent(out) :: ermesg
integer,                     intent(out) :: error

      real    :: rbc, rbe, qc, qe, tvc_v_s, tve_v_s
      integer :: ieqv_s
      integer :: k

      ermesg = ' ' ; error = 0

      xcape = 0.
      delt(1:klfc) = 0.0
      sum_xcape(1:klfc) = 0.0

!--------------------------------------------------------------------
!  calculate convective available potential energy.
!--------------------------------------------------------------------
      do k=klfc,nlev_hires-1

!--------------------------------------------------------------------
!  define flag to indicate which columns are actively computing cape.
!-------------------------------------------------------------------
        if (press(k+1) > plzb) then

!--------------------------------------------------------------------
!  define virtual temperature and specific humidity of parcel and 
!  environment.
!-------------------------------------------------------------------
          rbc = (parcel_r(k)+parcel_r(k+1))/2.
          rbe = (env_r(k)+env_r(k+1))/2.
          qc = rbc/(1. + rbc)
          qe = rbe/(1. + rbe)
          tvc_v_s = (parcel_t(k)+parcel_t(k+1))/2.
          tve_v_s = (env_t(k)+env_t(k+1))/2.
          tvc_v_s = tvc_v_s*(1. + Param%d608*qc)
          tve_v_s = tve_v_s*(1. + Param%d608*qe)

!--------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
          call don_u_numbers_are_equal_k   &
               (tvc_v_s, tve_v_s, ermesg, error, ieqv_s)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return


!---------------------------------------------------------------------
!   add the contribution to column cape from this pressure layer.
!---------------------------------------------------------------------
          if (ieqv_s >= 0) then
            delt(k) = Param%rdgas*(tvc_v_s - tve_v_s)*    &
                      alog(press(k)/press(k+1))
            sum_xcape(k) = sum_xcape(k-1) + delt(k)   
          endif

!---------------------------------------------------------------------
!   print out cape and cape contribution from this level.
!---------------------------------------------------------------------
        else
          delt(k:nlev_hires) = 0.0
          sum_xcape(k:nlev_hires) = sum_xcape(k-1)
          xcape = sum_xcape(k-1)
          exit
        endif
      end do  ! end of k loop

!-----------------------------------------------------------------------


end subroutine don_c_calculate_cape_k



!######################################################################

subroutine don_c_generate_cape_sounding_k   &
         (nlev_lsm, nlev_hires, temp, mixing_ratio, pfull, model_t, &
          model_r, model_p, cape_p, env_t, env_r, ermesg, error)  

implicit none

!-------------------------------------------------------------------
!    subroutine generate_cape_sounding reverses the vertical index and
!    saves the input temperature (temp) and vapor mixing ratio 
!    (mixing_ratio) fields on the large-scale model pressure grid 
!    (pfull) as model_t, model_r and model_p. these output variables
!    have index 1 nearest the surface, while the input fields have index
!    1 nearest the upper boundary. after reversal, the fields are 
!    interpolated onto an enhanced vertical grid defined by cape_p, 
!    producing output fields of temperature (env_t) and mixing ratio 
!    (env_r) on the enhanced grid. these fields also have index 1 near-
!    est the surface.
!-------------------------------------------------------------------

integer,                     intent(in)  :: nlev_lsm, nlev_hires
real, dimension(nlev_lsm),   intent(in)  :: temp, mixing_ratio, pfull
real, dimension(nlev_lsm),   intent(out) :: model_t, model_r, model_p
real, dimension(nlev_hires), intent(out) :: cape_p, env_t, env_r
character(len=*),            intent(out) :: ermesg
integer,                     intent(out) :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      temperature      temperature profile on large-scale model grid 
!                       that is to be used for the cape calculation 
!                       [ deg K ]
!      mixing_ratio     water vapor mixing ratio profile on the large-
!                       scale model grid that is to be used in the cape 
!                       calculation [ kg(h2o) / kg(dry air) ]
!      pfull            large-scale model full-level pressure profile 
!                       [ Pa ]
!
!   intent(out) variables:
!
!      model_t          large-scale model temperature field with vert-
!                       ical index 1 nearest the surface [ deg K ]
!      model_r          large-scale model vapor mixing ratio field with
!                       vertical index 1 nearest the surface 
!                       [ kg(h2o) /kg(dry air) ]
!      model_p          large-scale model pressure field with vertical
!                       index 1 nearest the surface [ Pa ]
!      cape_p           high-resolution pressure profile used for cape 
!                       calculation. index 1 nearest surface.  [ Pa ]
!      env_t            high-resolution temperature profile used for 
!                       cape calculation. index 1 nearest surface. 
!                       [ deg K ]
!      env_r            high-resolution vapor mixing ratio profile used
!                       for cape calculation. index 1 nearest thea
!                       surface. [ kg(h2o) / kg(dry air) ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real     :: dp     ! pressure difference between levels in 
                         ! high-resolution model [ Pa ]
      integer  :: k      ! do-loop indices


      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    create arrays of moisture, temperature and pressure having vertical
!    index 1 being closest to the surface. require that the mixing ratio
!    be non-negative.
!--------------------------------------------------------------------
      do k=1,nlev_lsm
        model_r(nlev_lsm+1-k) = amax1 (mixing_ratio(k), 0.0e00)
        model_t(nlev_lsm+1-k) = temp (k)
        model_p(nlev_lsm+1-k) = pfull(k)
      end do

!-------------------------------------------------------------------
!   define the vertical resolution of the convection parameterization
!   grid. define the top level pressure in that grid to be zero.
!   interpolate to define the pressure levels of that grid.
!-------------------------------------------------------------------
      dp = (model_p(1) - model_p(nlev_lsm))/(nlev_hires-1)
      cape_p(nlev_hires) = 0.
      do k=1,nlev_hires-1
        cape_p(k) = model_p(1) - (k-1)*dp     
      end do

!--------------------------------------------------------------------
!   call map_lo_res_col_to_hi_res_col to interpolate the large-scale
!   model grid values of temperature and vapor mixing ratio to the vert-
!   ical grid to be used for the cape calculation. ensure that the 
!   vapor mixing ratio field is positive-definite.
!--------------------------------------------------------------------
      call don_u_lo1d_to_hi1d_k    &
           (nlev_lsm, nlev_hires, model_t, model_p, cape_p,   &
            env_t, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      call don_u_lo1d_to_hi1d_k    &
           (nlev_lsm, nlev_hires, model_r, model_p, cape_p,    &
            env_r, ermesg, error)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      do k=1,nlev_hires
         env_r(k) = MAX(env_r(k), 0.0)
      end do

!--------------------------------------------------------------------


end subroutine don_c_generate_cape_sounding_k

!#####################################################################

subroutine don_c_integrate_vapor_k   &
           (nlev_hires, Param, env_r, cape_p, qint, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
use donner_types_mod, only : donner_param_type

implicit none

!----------------------------------------------------------------------
integer,                      intent(in)  :: nlev_hires
type(donner_param_type),      intent(in)  :: Param
real, dimension(nlev_hires),  intent(in)  :: env_r, cape_p
real,                         intent(out) :: qint            
character(len=*),             intent(out) :: ermesg
integer,                      intent(out) :: error

      integer  ::  k
      real     ::  sum

      ermesg= ' ' ; error = 0

      sum  = env_r(1)*(cape_p(1) - cape_p(2))
      do k=2,nlev_hires-1
        sum = sum + env_r(k)*0.5*(cape_p(k-1) - cape_p(k+1))
      end do
      sum = sum + env_r(nlev_hires)*    &
                             (cape_p(nlev_hires-1) - cape_p(nlev_hires))
      qint = sum/Param%grav 

!---------------------------------------------------------------------


end subroutine don_c_integrate_vapor_k



!#####################################################################


subroutine don_c_cape_diagnostics_k   &
         (isize, jsize, nlev_lsm, nlev_hires, Col_diag, Don_cape,  &
          exit_flag, ermesg, error)

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

use donner_types_mod, only : donner_column_diag_type, donner_cape_type

implicit none

!----------------------------------------------------------------------
integer,                         intent(in)    :: isize, jsize,   &
                                                  nlev_lsm, nlev_hires
type(donner_column_diag_type),   intent(in)    :: Col_diag
type(donner_cape_type),          intent(inout) :: Don_cape
logical, dimension(isize,jsize), intent(in)    :: exit_flag
character(len=*),                intent(out)   :: ermesg
integer,                         intent(out)   :: error


      integer  :: idiag, jdiag, unitdiag
      integer  :: n, k

      ermesg= ' ' ; error = 0

      
!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          idiag = Col_diag%i_dc(n)
          jdiag = Col_diag%j_dc(n)
          unitdiag = Col_diag%unit_dc(n)
          if ( .not. exit_flag(idiag, jdiag)) then
            if (Don_cape%plfc(idiag,jdiag) /= 0.0 .and.  &
                Don_cape%plzb(idiag,jdiag) /= 0.0 ) then
              do k=1,nlev_lsm - Col_diag%kstart+1 
                write (unitdiag, '(a, i4, f20.14, e20.12, e15.5)') &
                       'calculate_cape input profiles:&
                        & k,temp, mixing ratio, pressure',  k,  &
                        Don_cape%model_t(idiag,jdiag,k), &
                        Don_cape%model_r(idiag,jdiag,k), &
                        Don_cape%model_p(idiag,jdiag,k)
               end do
               write (Col_diag%unit_dc(n), '(a, 2f19.10)')  &
                         'in donner_deep: plfc,plzb= ',  &
                          Don_cape%plfc(idiag,jdiag),  & 
                          Don_cape%plzb(idiag,jdiag)
            endif 
            write (unitdiag, '(a, 3f19.10)')  &
                     'in donner_deep: plcl,coin,xcape= ',   &
                     Don_cape%plcl(idiag,jdiag),  &
                     Don_cape%coin(idiag,jdiag),   &
                     Don_cape%xcape(idiag,jdiag)
            if (Don_cape%plfc(idiag,jdiag) /= 0.0 .and.  &
                Don_cape%plzb(idiag,jdiag) /= 0.0 ) then
              do k=1,nlev_hires              
                write (unitdiag, '(a, i4, f19.10)') &
                  'in donner_deep: k,cape_p= ',k,  &
                      Don_cape%cape_p(idiag,jdiag,k)
                write (unitdiag, '(a, i4, 2f20.14)')  &
                   'in donner_deep: k,tcape,tpca= ',k,   &
                    Don_cape%env_t(idiag,jdiag,k),   &
                    Don_cape%parcel_t(idiag,jdiag,k)
                write (unitdiag, '(a, i4, 2e20.12)')  &
                    'in donner_deep: k,rcape,rpca= ',k,   &
                    Don_cape%env_r(idiag,jdiag,k),    &
                    Don_cape%parcel_r(idiag,jdiag,k)
                if (Don_cape%cape_p(idiag,jdiag,k) <     &
                    Don_cape%plzb(idiag,jdiag))  exit
              end do
            endif
            write (unitdiag, '(a, f20.10)')  &
               'integrate_vapor: qint= ',    &
                   Don_cape%qint(idiag,jdiag)  
          endif 
        end do  ! (n loop)
      endif

!---------------------------------------------------------------------



end subroutine don_c_cape_diagnostics_k 


!####################################################################

subroutine don_c_calculate_lcl_k    &
         (nlev_hires, Param, starting_temp, press, env_r, env_t, &
          parcel_r, parcel_t, plcl, tlcl, rlcl, klcl, cape_exit, ermesg, error)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

use donner_types_mod, only :  donner_param_type
use sat_vapor_pres_k_mod, only: compute_qs_k

implicit none

!---------------------------------------------------------------------
integer,                     intent(in)    :: nlev_hires
type(donner_param_type),     intent(in)    :: Param
real,                        intent(in)    :: starting_temp
real, dimension(nlev_hires), intent(in)    :: press, env_r, env_t 
real, dimension(nlev_hires), intent(inout) :: parcel_r, parcel_t      
real,                        intent(out)   :: plcl, tlcl, rlcl
integer,                     intent(out)   :: klcl
logical,                     intent(out)   :: cape_exit
character(len=*),            intent(out)   :: ermesg
integer,                     intent(out)   :: error

      real     :: qs_v_s, dtdp_v_s, dt_v_s, cp_v_s, q_ve_s, tp_s
      integer  :: ieqv_s       
      integer  :: k, nbad
      
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      q_ve_s    = parcel_r(Param%istart)/( 1. + parcel_r(Param%istart))
      cp_v_s    = Param%cp_air*(1. +             &
                             ((Param%cp_vapor/Param%cp_air) - 1.)*q_ve_s)
      plcl = 0.0     
      klcl = nlev_hires - 1
      cape_exit = .false.
      tp_s = starting_temp

!--------------------------------------------------------------------
!  move the parcel upwards to find the lcl in the column.
!--------------------------------------------------------------------
      do k=Param%istart,nlev_hires  ! k loop to find lcl

!--------------------------------------------------------------------
!  if the temperature and pressure are still within limits, continue 
!  parcel movement. determine saturation specific humidity for parcels 
!  at this level.
!---------------------------------------------------------------------
        if (tp_s >= Param%tmin .and.      &
            press(k) >= Param%upper_limit_for_lcl) then
          call compute_qs_k (tp_s, press(k), Param%d622, Param%d608, &
                             qs_v_s, nbad, q = q_ve_s)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (nbad /= 0) then
            ermesg = 'subroutine don_c_calculate_lcl_k: Temperatures out of range of esat table'
            error  = 1
          endif

!--------------------------------------------------------------------
!  check if the parcel is now saturated.
!---------------------------------------------------------------------
          call don_u_numbers_are_equal_k  &
               (qs_v_s, q_ve_s, ermesg, error, ieqv_s)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!--------------------------------------------------------------------
!   if saturation is exact or if parcel is super-saturated at its 
!   starting level, save pressure, temp, mixing ratio and cloud base 
!   level. exit vertical loop.
!--------------------------------------------------------------------
          if ( (ieqv_s ==  0) .or. &
               (ieqv_s < 0 .and. k == Param%istart)) then
            plcl = press(k)  
            rlcl = env_r(k)
            tlcl = env_t(k)
            klcl = k
            exit
          endif

!--------------------------------------------------------------------
!   if parcel is super-saturated, define cloud-base pressure, temp 
!   and mixing ratio as the average value between the current level
!   and the next lower level, and this level as the cloud base level.
!   exit the column.
!--------------------------------------------------------------------
          if (ieqv_s < 0) then
            plcl = (press(k) + press(k-1))/2.
            tlcl = (env_t(k) + env_t(k-1))/2.
            rlcl = (env_r(k) + env_r(k-1))/2.
            klcl = k
            exit

!---------------------------------------------------------------------
!    if the parcel remains unsaturated at this level and the top of the 
!    model has not been reached, move parcel along dry adiabat to next 
!    pressure level. define temperature at this level; verify that it 
!    is warmer than tmin. save parcel temperature and mixing ratio at 
!    this next higher level.
!---------------------------------------------------------------------
          else  ! (ieqv_s < 0) 
            if (k < nlev_hires) then
              dtdp_v_s = Param%rdgas*tp_s/cp_v_s   
              dt_v_s = dtdp_v_s*alog( press(k+1)/press(k))
              tp_s = tp_s + dt_v_s    
              if (tp_s < Param%tmin)  then
                cape_exit = .true.
                parcel_t(k+1:nlev_hires) = env_t(k+1:nlev_hires)
                parcel_r(k+1:nlev_hires) = env_r(k+1:nlev_hires)
                exit
              else  
                parcel_t(k+1) = tp_s   
                parcel_r(k+1) = parcel_r(Param%istart)
              endif

!-------------------------------------------------------------------
!    if have reached top of model, set flag to stop integration in this
!    column.
!-------------------------------------------------------------------
            else
              cape_exit    = .true.
            endif
          endif ! (ieqv_s < 0)

!--------------------------------------------------------------------
!    if either parcel temperature or pressure is below cutoff values,
!    set remainder of parcel sounding to the environment and stop 
!    searching in this column.
!--------------------------------------------------------------------
        else
          parcel_t(k+1:nlev_hires) = env_t(k+1:nlev_hires)
          parcel_r(k+1:nlev_hires) = env_r(k+1:nlev_hires)
          cape_exit = .true.
        endif
      end do   ! k loop to find lcl

!------------------------------------------------------------------


end subroutine don_c_calculate_lcl_k



!#####################################################################



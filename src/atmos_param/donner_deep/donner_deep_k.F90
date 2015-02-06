!#VERSION NUMBER:
!  $Name: tikal $
!  $Id: donner_deep_k.F90,v 20.0 2013/12/13 23:17:18 fms Exp $

!module donner_deep_inter_mod

!#include "donner_deep_interfaces.h"

!end module donner_deep_inter_mod


!######################################################################

subroutine don_d_donner_deep_k   &
         (is, ie, js, je, isize, jsize, nlev_lsm, nlev_hires, ntr, me, &
          cloud_tracers_present, cbmf,  &
          dt, Param, Nml, temp, mixing_ratio, pfull, phalf,   &
          zfull, zhalf, omega, pblht, tkemiz, qstar, cush, coldT, qlin,&!miz
          qiin, qain, land, sfc_sh_flux, sfc_vapor_flux, tr_flux,  &
          tracers, cell_cld_frac, cell_liq_amt, cell_liq_size,  &
          cell_ice_amt, cell_ice_size, cell_droplet_number, &
          meso_cld_frac, meso_liq_amt,  &
          meso_liq_size, meso_ice_amt, meso_ice_size,    &
          meso_droplet_number, nsum, & 
          precip, delta_temp, delta_vapor, detf, uceml_inter, mtot,   &
          mfluxup, &
          donner_humidity_area, donner_humidity_factor, total_precip,  &
          temperature_forcing, moisture_forcing, parcel_rise, &
          delta_ql, delta_qi, delta_qa, qtrtnd, calc_conv_on_this_step, &
          mhalf_3d, &
          ermesg, error, Initialized, Col_diag, Don_rad, Don_conv, Don_cape, &
          Don_cem, Don_save, sd, Uw_p, ac, cp, ct, Don_budgets)
                        
!-------------------------------------------------------------------
!    subroutine don_d_donner_deep_k is the primary kernel sub-
!    routine of the donner deep convection parameterization. it receives
!    all input needed from donner_deep_mod and controls the generation
!    of output that is returned to donner_deep_mod, from which it is made
!    accessible to the rest of the model parameterizations, as needed.
!-------------------------------------------------------------------

use donner_types_mod, only : donner_initialized_type, donner_save_type,&
                             donner_rad_type, donner_nml_type, &
                             donner_param_type, donner_conv_type, &
                             donner_budgets_type, &
                             donner_column_diag_type, donner_cape_type, &
                             donner_cem_type
use  conv_utilities_k_mod,only : adicloud, sounding, uw_params
use  conv_plumes_k_mod,   only : cplume, ctend

implicit none

!--------------------------------------------------------------------
integer,                 intent(in)     :: is, ie, js, je, isize, jsize,&
                                           nlev_lsm, nlev_hires, ntr, me
logical,                 intent(in)     :: cloud_tracers_present
real, dimension(isize,jsize),    intent(inout)     :: cbmf
real,                    intent(in)     :: dt
type(donner_param_type), intent(in)     :: Param
type(donner_nml_type),   intent(inout)     :: Nml
real, dimension(isize,jsize,nlev_lsm),                                  &
                         intent(in)     :: temp, mixing_ratio, pfull,   zfull, &
                                           omega, qlin, qiin, qain, &
                                           cell_cld_frac,  cell_liq_amt,&
                                           cell_liq_size, cell_ice_amt, &
                                           cell_ice_size,  &
                                           cell_droplet_number, &
                                           meso_cld_frac,&
                                           meso_liq_amt, meso_liq_size, &
                                           meso_ice_amt, meso_ice_size,&
                                           meso_droplet_number
real,    dimension(isize,jsize,nlev_lsm+1),                            &
                         intent(in)     :: phalf, zhalf
real,    dimension(isize,jsize),                                      &
                         intent(in)     :: pblht, tkemiz, qstar, cush, land, &
                                           sfc_sh_flux, sfc_vapor_flux
logical, dimension(isize,jsize), intent(in) :: coldT    
real,    dimension(isize,jsize,ntr),                                 &
                         intent(in)     :: tr_flux 
real,    dimension(isize,jsize,nlev_lsm,ntr),                         &
                         intent(in)     :: tracers 
integer, dimension(isize,jsize),                                     &
                         intent(in)     :: nsum      
real,    dimension(isize,jsize),                                     & 
                         intent(out)    :: precip      
real, dimension(isize,jsize,nlev_lsm),                                 &
                         intent(out)    :: delta_temp, delta_vapor,&
                                           detf,  mtot, mfluxup, &
                                           donner_humidity_area,&
                                           donner_humidity_factor, &
                                           temperature_forcing,   &
                                           moisture_forcing, &
                                           delta_ql, delta_qi, &
                                           delta_qa
real, dimension(isize,jsize,nlev_lsm+1),                               &
                         intent(out)    :: uceml_inter
real, dimension(isize,jsize),                                         &
                         intent(out)    :: total_precip, parcel_rise
real,    dimension(isize,jsize,nlev_lsm,ntr),                        &
                         intent(out)    :: qtrtnd 
logical,                 intent(in)     :: calc_conv_on_this_step
real, dimension(isize,jsize,nlev_lsm+1),                               &
                         intent(out)    :: mhalf_3d
character(len=*),        intent(out)    :: ermesg
integer,                 intent(out)    :: error
type(donner_initialized_type),                            &
                         intent(inout)  :: Initialized
type(donner_column_diag_type),                               &
                         intent(inout)  :: Col_diag
type(donner_rad_type),   intent(inout)  :: Don_rad
type(donner_conv_type),  intent(inout)  :: Don_conv
type(donner_budgets_type),  intent(inout)  :: Don_budgets
type(donner_cape_type),  intent(inout)  :: Don_cape
type(donner_cem_type),   intent(inout)  :: Don_cem
type(donner_save_type),  intent(inout)  :: Don_save

type(sounding),          intent(inout)  :: sd
type(uw_params),          intent(inout)  :: Uw_p
type(adicloud),          intent(inout)  :: ac
type(cplume),            intent(inout)  :: cp
type(ctend),             intent(inout)  :: ct

!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     omega          model omega field at model full levels [ Pa / sec ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     cell_cld_frac  fractional coverage of convective cells in
!                    grid box [ dimensionless ]
!     cell_liq_amt   liquid water content of convective cells
!                    [ kg(h2o) / kg(air) ]
!     cell_liq_size  assumed effective size of cell liquid drops
!                    [ microns ]
!     cell_ice_amt   ice water content of cells
!                    [ kg(h2o) / kg(air) ]
!     cell_ice_size  generalized effective diameter for ice in
!                    convective cells [ microns ]
!     meso_cld_frac  fractional area of mesoscale clouds in grid
!                    box [ dimensionless ]
!     meso_liq_amt   liquid water content in mesoscale clouds
!                    [ kg(h2o) / kg(air) ]
!     meso_liq_size  assumed effective size of mesoscale drops
!                    [ microns ]
!     meso_ice_amt   ice water content of mesoscale elements
!                    [ kg(h2o) / kg(air) ]
!     meso_ice_size  generalized ice effective size for anvil ice
!                    [ microns ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux   sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     nsum           number of time levels over which the above variables
!                    have so far been summed
!
!   intent(out) variables:
!
!     precip         precipitation generated by deep convection
!                    [ kg(h2o) / m**2 ]
!     delta_temp     temperature increment due to deep convection 
!                    [ deg K ]
!     delta_vapor    water vapor mixing ratio increment due to deep 
!                    convection [ kg(h2o) / kg (dry air) ]
!     detf           detrained cell mass flux at model levels 
!                    [ (kg / (m**2 sec) ) ]
!     mtot           mass flux at model full levels, convective plus 
!                    mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     mfluxup        upward mass flux at model full levels, convective 
!                    plus mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing 
!                    vapor mixing ratio tendency due to donner
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     delta_ql       cloud water specific humidity increment due to 
!                    deep convection over the timestep
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qi       cloud ice specific humidity increment due to deep 
!                    convection over the timestep 
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qa       cloud area increment due to deep convection
!                    over the time step [ fraction ]
!     uceml_inter    upward cell mass flux at interface levels 
!                    [ (kg / (m**2 sec) ) ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!     parcel_rise    accumulated vertical displacement of a
!                    near-surface parcel as a result of the lowest
!                    model level omega field [ Pa ]
!     qtrtnd         tracer time tendencies due to deep convection
!                    during the time step
!                    [ kg(tracer) / (kg (dry air) sec) ]
!     calc_conv_on_this_step
!                    is this a step on which to calculate
!                    donner convection ?
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!   intent(inout) variables:
!
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Don_rad        donner_rad_type derived type variable used to hold 
!                    those fields needed to connect the donner deep 
!                    convection parameterization and the model radiation 
!                    package
!     Don_conv       donner_convection_type derived type variable
!                    containing diagnostics and intermediate results 
!                    describing the nature of the convection produced by
!                    the donner parameterization
!     Don_cape       donner_cape type derived type variable containing 
!                    diagnostics and intermediate results related to the
!                    cape calculation associated with the donner convec-
!                    tion parameterization
!     Don_save       donner_save_type derived type variable containing
!                    those variables which must be preserved across
!                    timesteps
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!     
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real,                                                      &
         dimension (isize, jsize, nlev_lsm) ::  lag_cape_temp,        &
                                                lag_cape_vapor,      &
                                                lag_cape_press, &
                                                dql, dqi, dqa
      real,                                                      &
         dimension (isize, jsize, nlev_lsm+1) ::  mhalf_3d_local  
      real,    dimension (isize, jsize)     ::  current_displ
      logical, dimension (isize, jsize)     ::  exit_flag
      integer                               ::  idiag, jdiag, unitdiag

      integer  :: i, j, k, n   

!--------------------------------------------------------------------
!   local variables:
!
!     lag_cape_temp        temperature field used in lag-time cape 
!                          calculation [ deg K ]
!     lag_cape_vapor       vapor mixing ratio field used in lag-time
!                          cape calculation [ kg(h2o) / kg(dry air) ]
!     lag_cape_press       model full-level pressure field used in 
!                          lag-time cape calculation 
!                          [ kg(h2o) / kg(dry air) ]
!     dql                  tendency of cloud liquid specific humidity
!                          due to donner convection 
!                          [ kg(h2o) / kg(moist air) / sec ]
!     dqi                  tendency of cloud ice specific humidity
!                          due to donner convection 
!                          [ kg(h2o) / kg(moist air) / sec ]
!     dqa                  tendency of large-scale cloud area
!                          due to donner convection 
!                          [ fraction / sec ]
!     current_displ        low-level parcel displacement to use in cape
!                          calculation on this step [ Pa ]
!     exit_flag            logical array indicating whether deep conv-
!                          ection exists in a column
!     idiag                physics window i index of current diagnostic
!                          column
!     jdiag                physics window j index of current diagnostic
!                          column
!     unitdiag             i/o unit assigned to current diagnostic
!                          column
!     i, j, k, n           do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!-------------------------------------------------------------------
!    verify that donner_deep_freq is an integral multiple of the model
!    physics time step.
!--------------------------------------------------------------------
      if (MOD (Nml%donner_deep_freq, int(dt)) /= 0) then
        ermesg = 'donner_deep_donner_deep_k: donner_deep timestep NOT &
            &an integral multiple of physics timestep'
        error = 1
        return
      endif

!---------------------------------------------------------------------
!    perform the following calculations only if this is a step upon
!    which donner convection is to be calculated.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then
 
!-------------------------------------------------------------------
!    call initialize_local_variables_k to allocate and initialize the
!    elements of the donner_conv, donner_cape and donner_rad derived type
!    variables.
!-------------------------------------------------------------------
        call don_d_init_loc_vars_k       &
             (isize, jsize, nlev_lsm, ntr, nlev_hires, cell_cld_frac, &
              cell_liq_amt, cell_liq_size, cell_ice_amt, cell_ice_size, &
              cell_droplet_number, &
              meso_cld_frac, meso_liq_amt, meso_liq_size, meso_ice_amt, &
              meso_ice_size, meso_droplet_number, nsum, Don_conv,   &
              Don_cape, Don_rad, Don_cem, Param, Don_budgets, Nml, &
              Initialized, sd, ac, cp, ct, ermesg, error) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return
      endif !(calc_conv_on_this_step)

!---------------------------------------------------------------------
!    allocate and initialize the components of the donner_budgets_type 
!    variable Don_budgets. definitions of these arrays are found in 
!    donner_types.h.
!---------------------------------------------------------------------
      allocate ( Don_budgets%liq_prcp     (isize, jsize, nlev_lsm) )
      allocate ( Don_budgets%frz_prcp     (isize, jsize, nlev_lsm) )
      Don_budgets%liq_prcp     = 0.
      Don_budgets%frz_prcp     = 0.
   if (Initialized%do_conservation_checks .or.   &
                                          Nml%do_budget_analysis) then
      allocate ( Don_budgets%lheat_precip (isize, jsize) )
      allocate ( Don_budgets%vert_motion  (isize, jsize) )
      Don_budgets%lheat_precip = 0.
      Don_budgets%vert_motion  = 0.
      allocate ( Don_budgets%water_budget (isize, jsize, nlev_lsm, &
                                        Don_budgets%n_water_budget) )
      allocate ( Don_budgets%enthalpy_budget (isize, jsize, nlev_lsm, &
                                     Don_budgets%n_enthalpy_budget) )
      allocate ( Don_budgets%precip_budget (isize, jsize, nlev_lsm, &
                                       Don_budgets%n_precip_paths, &
                                       Don_budgets%n_precip_types)  )
      Don_budgets%water_budget = 0.
      Don_budgets%enthalpy_budget = 0.
      Don_budgets%precip_budget = 0.
   endif

!---------------------------------------------------------------------
!    add the vertical displacement resulting from the current omega
!    field to the accumulated displacement of a parcel which originated
!    at the lowest model full level. prevent the parcel from moving 
!    below its starting point or going out the top of the atmosphere. 
!---------------------------------------------------------------------
      do j=1,jsize       
        do i=1,isize        
          parcel_rise(i,j) = Don_save%parcel_disp(i+is-1,j+js-1) +  &
                             omega(i,j,nlev_lsm)*dt
          parcel_rise(i,j) = MIN (0.0, parcel_rise(i,j))
          parcel_rise(i,j) = MAX (-phalf(i,j,nlev_lsm+1),    &
                                  parcel_rise(i,j))
        end do
      end do

!---------------------------------------------------------------------
!    if there are one or more diagnostic columns in the current physics
!    window, set a flag to so indicate. call don_d_column_input_fields
!    to print out the relevant input fields, location and control 
!    information for these diagnostics columns.   
!---------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        if (Col_diag%ncols_in_window > 0) then
          Col_diag%in_diagnostics_window = .true.
          call don_d_column_input_fields_k    &
               (isize, jsize, nlev_lsm, dt, calc_conv_on_this_step, &
                Col_diag, temp, mixing_ratio, pfull, omega, phalf,   &
                parcel_rise, ermesg, error) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!---------------------------------------------------------------------
!    if there are no diagnostic columns in the current physics
!    window, set a flag to so indicate. 
!---------------------------------------------------------------------
        else
          Col_diag%in_diagnostics_window = .false.
        endif

!---------------------------------------------------------------------
!    if column diagnostics are not desired in any model columns, set a 
!    flag to so indicate. 
!---------------------------------------------------------------------
      else
        Col_diag%in_diagnostics_window = .false.
      endif

!---------------------------------------------------------------------
!    perform the following calculations only if this is a step upon 
!    which donner convection is to be calculated.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then 

!---------------------------------------------------------------------
!    define the low-level displacement to be used on this step 
!    (current_displ). it is the current time-integrated value, unless 
!    the current lowest-level vertical velocity is downward, in which 
!    case the displacement to be used on the current step is set to 
!    zero, precluding deep convection on this step.
!---------------------------------------------------------------------
        do j=1,jsize       
          do i=1,isize        
            if (omega(i,j,nlev_lsm) > 0.)   then
              current_displ(i,j) = 0. 
            else
              current_displ(i,j) = parcel_rise(i,j)
            endif
          end do
        end do

!---------------------------------------------------------------------
!    define the temperature, vapor mixing ratio and pressure fields to
!    be used in calculating the lag values of cape so that a cape tend-
!    ency due to large-scale forcing may be computed.
!---------------------------------------------------------------------
        lag_cape_temp (:,:,:) = Don_save%lag_temp (is:ie,js:je,:)
        lag_cape_vapor(:,:,:) = Don_save%lag_vapor(is:ie,js:je,:)
        lag_cape_press(:,:,:) = Don_save%lag_press(is:ie,js:je,:)

!---------------------------------------------------------------------
!    call donner_convection_driver to calculate the effects of deep 
!    convection.  
!--------------------------------------------------------------------- 
        call don_d_convection_driver_k   &
             (isize, jsize, nlev_lsm, nlev_hires, ntr, me, &
              cloud_tracers_present, cbmf, dt, Nml, &
              Initialized, Param, Col_diag, temp, mixing_ratio,  &
              pfull, zfull, zhalf, pblht, tkemiz, qstar, cush, coldT,  &!miz
              qlin, qiin, qain, lag_cape_temp, lag_cape_vapor,    &
              lag_cape_press, phalf, current_displ, land, sfc_sh_flux, &
              sfc_vapor_flux, tr_flux, tracers, Don_cape, Don_conv, &
              Don_rad, Don_cem, temperature_forcing, moisture_forcing,  &
              total_precip, donner_humidity_factor, donner_humidity_area,&
              dql, dqi, dqa, mhalf_3d_local, exit_flag, ermesg, error, sd, Uw_p, ac, cp, ct, &
              Don_budgets)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) return

!--------------------------------------------------------------------- 
!    define the module variables used to preserve output fields across
!    physics timesteps, needed when the donner parameterization is not
!    executed on every physics step.
!
!    1) the variables defining the humidity disturbance caused by
!       donner convection. these variables are needed so that the large-
!       scale humidity may be properly adjusted in strat_cloud_mod.
!--------------------------------------------------------------------- 
        Don_save%humidity_area (is:ie,js:je,:) =    &
                                            donner_humidity_area (:,:,:)
        Don_save%humidity_factor(is:ie,js:je,:) =    &
                                           donner_humidity_factor(:,:,:)
 
!--------------------------------------------------------------------- 
!    2) the total precipitation produced by the donner parameter-
!       ization.
!--------------------------------------------------------------------- 
        Don_save%tprea1(is:ie,js:je) = total_precip(:,:)

!---------------------------------------------------------------------
!    3) the vapor and temperature forcing resulting from the donner
!       deep parameterization that will be output to the calling 
!       routine. NOTE: these values of cemetf and cememf have the terms
!       related to flux convergence of condensate and mesoscale 
!       detrainment removed when parameterization is run in a model 
!       using strat_cloud_mod, since in that case those terms will be 
!       calculated within that module.
!---------------------------------------------------------------------
        Don_save%cememf(is:ie,js:je,:) = moisture_forcing(:,:,:) 
        Don_save%cemetf(is:ie,js:je,:) = temperature_forcing(:,:,:)

!----------------------------------------------------------------------
!    4) the increments which must be applied to the strat_cloud vari-
!       ables as a result of donner convection. the returned tendencies 
!       are converted to increments. 
!----------------------------------------------------------------------
        if (cloud_tracers_present) then
        Don_save%dql_strat(is:ie,js:je,:) = dql(:,:,:)*dt
        Don_save%dqi_strat(is:ie,js:je,:) = dqi(:,:,:)*dt
        Don_save%dqa_strat(is:ie,js:je,:) = dqa(:,:,:)*dt
        endif

!--------------------------------------------------------------------
!    5) the net mass flux and detrained cell mass flux at model full 
!       levels that is associated with the donner deep parameterization. 
!       the net mass flux is needed as input to strat_cloud_mod, while
!       the detrained cell mass flux is needed by cu_mo_trans_mod. 
!--------------------------------------------------------------------
        do k=1,nlev_lsm
          do j=1,jsize
            do i=1,isize
              if ((Don_conv%uceml(i,j,k) <= 1.0e-10) .and.   &
                  (Don_conv%umeml(i,j,k) <= 1.0e-10) .and.   &
                  (Don_conv%dmeml(i,j,k) <= 1.0e-10) ) then
                Don_save%mass_flux(i+is-1,j+js-1,k) = 0.
              else
                Don_save%mass_flux(i+is-1,j+js-1,k) =   &
                                        Don_conv%uceml(i,j,k) + &
                                        Don_conv%umeml(i,j,k) + &
                                        Don_conv%dmeml(i,j,k)
              endif
              if ((Don_conv%uceml(i,j,k) <= 1.0e-10) .and.   &
                  (Don_conv%umeml(i,j,k) <= 1.0e-10) ) then
                Don_save%mflux_up(i+is-1,j+js-1,k) = 0.
              else
                Don_save%mflux_up(i+is-1,j+js-1,k) =   &
                                        Don_conv%uceml(i,j,k) + &
                                        Don_conv%umeml(i,j,k)
              endif
              if (Don_conv%detmfl(i,j,k) <= 1.0e-10) then
                Don_save%det_mass_flux(i+is-1,j+js-1,k) = 0.
              else
                Don_save%det_mass_flux(i+is-1,j+js-1,k) =      &
                                                  Don_conv%detmfl(i,j,k)
              endif
            end do
          end do
        end do

!--------------------------------------------------------------------
!    6) the upward mass flux at model interface levels that is 
!       associated with the convective cells present in the donner deep 
!       convction parameterization. this is needed by cu_mo_trans_mod.
!       define values at upper and lower boundary to be 0.0.
!--------------------------------------------------------------------
        Don_save%cell_up_mass_flux(:,:,1) = 0.
        Don_save%cell_up_mass_flux(:,:,nlev_lsm+1) = 0.
        do k=2,nlev_lsm
          do j=1,jsize
            do i=1,isize
              Don_save%cell_up_mass_flux(i+is-1,j+js-1,k) =  &
                                    0.5*(Don_Conv%uceml(i,j,k) + &
                                         Don_conv%uceml(i,j,k-1))
            end do
          end do
        end do
        
!--------------------------------------------------------------------
!    7) the tracer time tendencies resulting from the donner param-
!       eterization. 
!--------------------------------------------------------------------
        if (Initialized%do_donner_tracer) then
          Don_save%tracer_tends(is:ie,js:je,:,:) =      &
                   Don_conv%qtceme(:,:,:,:) + Don_conv%wetdept(:,:,:,:)
        endif
      else  ! (calc_conv_on_this_step)
        total_precip(:,:) = 0.0

!---------------------------------------------------------------------
!    end of if loop for code executed only on steps for which the donner
!    parameterization is requested.
!---------------------------------------------------------------------
      endif ! (calc_conv_on_this_step)

!--------------------------------------------------------------------- 
!    update the module variable containing the total low-level parcel 
!    displacement. this field is updated on every model physics step, 
!    regardless of whether donner convective tendencies are calculated 
!    or not.
!--------------------------------------------------------------------- 
      Don_save%parcel_disp(is:ie,js:je) = parcel_rise(:,:)

!----------------------------------------------------------------------
!    define the output fields to be passed back to the calling routine.
!    these fields are returned on every physics step, regardless of
!    whether or not the donner calculation is done, and so must be 
!    defined from module variables.
!
!    1) the temperature increment due to deep convection
!----------------------------------------------------------------------
      delta_temp(:,:,:) = Don_save%cemetf(is:ie,js:je,:)*dt

!----------------------------------------------------------------------
!    2) the moisture increment due to deep convection. if the moisture
!       tendency results in the production of a negative value, reduce 
!       the tendency to avoid producing the negative mixing ratio.
!----------------------------------------------------------------------
      do k=1,nlev_lsm
        do j=1,jsize      
          do i=1,isize       
            delta_vapor(i,j,k) = Don_save%cememf(i+is-1,j+js-1,k)*dt
!           if ((mixing_ratio(i,j,k) + delta_vapor(i,j,k)) < 0.) then
!             if (mixing_ratio(i,j,k) > 0.) then
!               delta_vapor (i,j,k) = -mixing_ratio(i,j,k)
!             else 
!               delta_vapor(i,j,k) = 0.0
!             endif
!           endif
          end do
        end do
      end do

!-------------------------------------------------------------------
!    3) the net mass flux, detrained cell mass flux and upward mass 
!       flux due to convective cells at interface levels resulting from
!       donner convection.
!-------------------------------------------------------------------
      mtot(:,:,:)        = Don_save%mass_flux(is:ie, js:je,:)
      mfluxup(:,:,:)     = Don_save%mflux_up(is:ie, js:je,:)
      detf(:,:,:)        = Don_save%det_mass_flux(is:ie,js:je,:)
      uceml_inter(:,:,:) = Don_save%cell_up_mass_flux(is:ie,js:je,:)
      mhalf_3d(:,:,:) = mhalf_3d_local

!-------------------------------------------------------------------
!    4) the increments of the large-scale cloud variables due to deep
!       convection and the variables describing the specific humidity
!       disturbance associated with donner convection.
!-------------------------------------------------------------------
      if (cloud_tracers_present) then
      delta_ql(:,:,:) = Don_save%dql_strat (is:ie, js:je,:)
      delta_qi(:,:,:) = Don_save%dqi_strat (is:ie, js:je,:)
      delta_qa(:,:,:) = Don_save%dqa_strat (is:ie, js:je,:)
      endif
      donner_humidity_area(:,:,:)  =             &
                                   Don_save%humidity_area(is:ie,js:je,:)
      donner_humidity_factor(:,:,:) =             &
                                 Don_save%humidity_factor(is:ie,js:je,:)

!----------------------------------------------------------------------
!    5) the precipitation accrued on the current timestep from deep
!       convection. 
!       note: precip    [mm/day] * 1/86400 [day/sec] * 1/1000 [ m/mm] * 
!                  1000 [kg(h2o)/m**3] * dt [sec] = kg/m2, as desired. 
!----------------------------------------------------------------------
      precip(:,:) = Don_save%tprea1(is:ie,js:je)*dt/Param%seconds_per_day

!--------------------------------------------------------------------
!    6) time tendencies of any tracers being transported by donner 
!       convection. if none have been defined, fill the output array
!       with zeroes.
!--------------------------------------------------------------------
      if (Initialized%do_donner_tracer) then
        qtrtnd(:,:,:,:) = Don_save%tracer_tends(is:ie,js:je,:,:)
      else
        qtrtnd(:,:,:,:) = 0.0                            
      endif

!---------------------------------------------------------------------
!    if this is a diagnostics window, output the increments to temper-
!    ature and vapor mixing ratio at levels where donner convection
!    has produced a modification.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          idiag = Col_diag%i_dc(n)
          jdiag = Col_diag%j_dc(n)
          unitdiag = Col_diag%unit_dc(n)
          do k=Col_diag%kstart,nlev_lsm
            if (delta_temp(idiag,jdiag,k) /= 0.0) then
              write (unitdiag, '(a, i4, f20.14, e20.12)') &
                   'in donner_deep: k,ttnd,qtnd',  k,  &
                   delta_temp(idiag,jdiag,k), delta_vapor(idiag,jdiag,k)
            endif
          end do
        end do
      endif

!---------------------------------------------------------------------
!    define the module variables containing the temperature, pressure 
!    and vapor fields that are to be used on the next time step to 
!    calculate a lag-time cape so that the time tendency of cape due
!    to large-scale forcing may be obtained. this field is updated on 
!    every model physics step, so that values are present in case the 
!    next step is a donner calculation step.
!--------------------------------------------------------------------- 
      Don_save%lag_temp (is:ie, js:je,:) = temp + delta_temp
      Don_save%lag_vapor(is:ie, js:je,:) = mixing_ratio + delta_vapor 
      Don_save%lag_press(is:ie, js:je,:) = pfull
  
!-------------------------------------------------------------------
!    perform the following calculations only if this is a step upon 
!    which donner convection is to be calculated.
!-------------------------------------------------------------------
      if (calc_conv_on_this_step) then

!---------------------------------------------------------------------
!    define the revised moisture tendency produced by donner convection
!    after it has been adjusted to prevent the generation of negative 
!    vapor mixing ratio.
!---------------------------------------------------------------------
        Don_conv%cememf_mod(:,:,:) = delta_vapor(:,:,:)/dt

!---------------------------------------------------------------------
!    if this is a diagnostics window, call donner_column_end_of_step
!    to output various diagnostic fields in the specified diagnostic
!    columns.
!---------------------------------------------------------------------
        if (Col_diag%in_diagnostics_window) then
          call don_d_column_end_of_step_k   &
               (isize, jsize, nlev_lsm, ntr, Col_diag, exit_flag,   &
                total_precip, parcel_rise, temperature_forcing,   &
                moisture_forcing, tracers, Don_cape,   &
                Don_conv, ermesg, error) 

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return
        endif
      endif !(calc_conv_on_this_step) 

!--------------------------------------------------------------------


end subroutine don_d_donner_deep_k




!####################################################################

subroutine don_d_init_loc_vars_k      &
         (isize, jsize, nlev_lsm, ntr, nlev_hires, cell_cld_frac,  &
          cell_liq_amt, cell_liq_size, cell_ice_amt, cell_ice_size,   &
          cell_droplet_number, &
          meso_cld_frac, meso_liq_amt, meso_liq_size, meso_ice_amt,   &
          meso_ice_size, meso_droplet_number, nsum, Don_conv,   &
          Don_cape, Don_rad, Don_cem, Param, Don_budgets, Nml, &
          Initialized, sd, ac, cp, ct, ermesg, error)

use donner_types_mod, only : donner_rad_type, donner_conv_type, &
                             donner_budgets_type, &
                             donner_initialized_type, &
                             donner_cape_type, donner_nml_type, &
                             donner_param_type, donner_cem_type
use  conv_utilities_k_mod,only : adicloud, sounding, ac_init_k, &
                                 sd_init_k
use  conv_plumes_k_mod,only    : cplume, ctend, cp_init_k, ct_init_k
implicit none

!--------------------------------------------------------------------
!   subroutine don_d_init_loc_vars_k allocates space 
!   for and initializes the array components of the donner_conv_type 
!   variable Don_conv, the donner_cape_type variable Don_cape, the 
!   donner_rad_type variable Don_rad, and the donner_cem_type 
!   variable Don_cem.
!--------------------------------------------------------------------

integer,                         intent(in)    :: isize, jsize,    &
                                                  nlev_lsm, ntr, &
                                                  nlev_hires
real,dimension(isize,jsize,nlev_lsm),                              &
                                 intent(in)    :: cell_cld_frac,  &
                                                  cell_liq_amt,   &
                                                  cell_liq_size, &
                                                  cell_ice_amt,   &
                                                  cell_ice_size, &
                                            cell_droplet_number, &
                                                  meso_cld_frac,    &
                                                  meso_liq_amt, &
                                                  meso_liq_size, &
                                                  meso_ice_amt,     &
                                                  meso_ice_size,  &
                                            meso_droplet_number 
integer, dimension(isize,jsize), intent(in)    :: nsum
type(donner_conv_type),          intent(inout) :: Don_conv
type(donner_cape_type),          intent(inout) :: Don_cape
type(donner_rad_type),           intent(inout) :: Don_rad
type(donner_cem_type),           intent(inout) :: Don_cem
type(donner_param_type),         intent(in)    :: Param
type(donner_budgets_type),       intent(inout) :: Don_budgets
type(donner_nml_type),           intent(in)    :: Nml
type(donner_initialized_type),   intent(in)    :: Initialized
type(sounding),               intent(inout) :: sd
type(adicloud),               intent(inout) :: ac
type(cplume),                 intent(inout) :: cp
type(ctend),                  intent(inout) :: ct
character(len=*),                intent(out)   :: ermesg
integer,                         intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize         size of x-dimension of physics window
!     jsize         size of y-dimension of physics window
!     nlev_lsm      number of layers in large-scale model
!     ntr           number of tracers to be transported by donner
!                   convection
!     nlev_hires    number of model layers in hi-res cloud model
!                   of the donner deep convection parameterization
!     cell_cld_frac fractional coverage of convective cells in
!                   grid box [ dimensionless ]
!     cell_liq_amt  liquid water content of convective cells
!                   [ kg(h2o) / kg(air) ]
!     cell_liq_size assumed effective size of cell liquid drops
!                   [ microns ]
!     cell_ice_amt  ice water content of cells
!                   [ kg(h2o) / kg(air) ]
!     cell_ice_size generalized effective diameter for ice in
!                   convective cells [ microns ]
!     meso_cld_frac fractional area of mesoscale clouds in grid
!                   box [ dimensionless ]
!     meso_liq_amt  liquid water content in mesoscale clouds
!                   [ kg(h2o) / kg(air) ]
!     meso_liq_size assumed effective size of mesoscale drops
!                   [ microns ]
!     meso_ice_amt  ice water content of mesoscale elements
!                   [ kg(h2o) / kg(air) ]
!     meso_ice_size generalized ice effective size for anvil ice
!                   [ microns ]
!     nsum          number of time levels over which the above variables
!                   have so far been summed
!
!   intent(inout) variables:
!
!     Don_conv     donner_conv_type derived type variable containing 
!                  diagnostics and intermediate results describing the 
!                  nature of the convection produced by the donner 
!                  parameterization
!     Don_cape     donner_cape type derived type variable containing 
!                  diagnostics and intermediate results related to the 
!                  cape calculation associated with the donner 
!                  convection parameterization
!     Don_rad      donner_rad_type derived type variable used to hold 
!                  those fields needed to connect the donner deep 
!                  convection parameterization and the model radiation 
!                  package
!     Don_cem      donner_cem_type derived type variable containing
!                  Donner cumulus ensemble member diagnostics
!     Param        donner_param_type variable containingthe parameters
!                  of the donner deep convection parameterization

!
!   intent(out) variables:
!
!     ermesg        character string containing any error message
!                   that is returned from a kernel subroutine
!
!---------------------------------------------------------------------
     integer :: ncem   ! Param%kpar, number of cumulus ensemble members

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg= ' ' ; error = 0

      call sd_init_k (nlev_lsm, ntr, sd)
      call ac_init_k (nlev_lsm, ac)
      call cp_init_k (nlev_lsm, ntr, cp)
      call ct_init_k (nlev_lsm, ntr, ct)

!---------------------------------------------------------------------
!    allocate the components of the donner_conv_type variable Don_conv.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_conv%cecon              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%ceefc              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cell_liquid_eff_diam     &
                                            (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cell_ice_geneff_diam     &
                                            (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cememf_mod         (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cemfc              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cmus               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%conv_temp_forcing  (isize, jsize, nlev_lsm) )
      allocate (Don_conv%conv_moist_forcing (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cual               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cuqi               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%cuql               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%detmfl             (isize, jsize, nlev_lsm) )
      allocate (Don_conv%dgeice             (isize, jsize, nlev_lsm) )
      allocate (Don_conv%dmeml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%ecds               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%eces               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%elt                (isize, jsize, nlev_lsm) )
      allocate (Don_conv%emds               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%emes               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%fre                (isize, jsize, nlev_lsm) )
      allocate (Don_conv%mrmes              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%tmes               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%uceml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%umeml              (isize, jsize, nlev_lsm) )
      allocate (Don_conv%wmms               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%wmps               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%xice               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%xliq               (isize, jsize, nlev_lsm) )
      allocate (Don_conv%a1                 (isize, jsize) )
      allocate (Don_conv%amax               (isize, jsize) )
      allocate (Don_conv%amos               (isize, jsize) )
      allocate (Don_conv%ampta1             (isize, jsize) )
      allocate (Don_conv%cell_precip        (isize, jsize) )
      allocate (Don_conv%dcape              (isize, jsize) )
      allocate (Don_conv%emdi_v             (isize, jsize) )
      allocate (Don_conv%meso_precip        (isize, jsize) )
      allocate (Don_conv%pb_v               (isize, jsize) )
      allocate (Don_conv%pmd_v              (isize, jsize) )
      allocate (Don_conv%przm               (isize, jsize) )
      allocate (Don_conv%prztm              (isize, jsize) )
      allocate (Don_conv%pzm_v              (isize, jsize) )
      allocate (Don_conv%pztm_v             (isize, jsize) )

      allocate (Don_conv%qtceme        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%qtmes1        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%qtren1        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%temptr        (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%wtp1          (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%wetdepc       (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%wetdepm       (isize, jsize, nlev_lsm, ntr) )
      allocate (Don_conv%wetdept       (isize, jsize, nlev_lsm, ntr) )

!---------------------------------------------------------------------
!    initialize the components of the donner_conv_type variable Don_conv.
!---------------------------------------------------------------------
      Don_conv%cecon                = 0.0
      Don_conv%ceefc                = 0.0
      Don_conv%cell_liquid_eff_diam = 0.0
      Don_conv%cell_ice_geneff_diam = 0.0
      Don_conv%cememf_mod           = 0.0
      Don_conv%cemfc                = 0.0
      Don_conv%cmus                 = 0.0
      Don_conv%conv_temp_forcing    = 0.0
      Don_conv%conv_moist_forcing   = 0.0
      Don_conv%cual                 = 0.0
      Don_conv%cuqi                 = 0.0
      Don_conv%cuql                 = 0.0
      Don_conv%detmfl               = 0.0
      Don_conv%dgeice               = 0.0
      Don_conv%dmeml                = 0.0
      Don_conv%ecds                 = 0.0
      Don_conv%eces                 = 0.0
      Don_conv%elt                  = 0.0
      Don_conv%emds                 = 0.0
      Don_conv%emes                 = 0.0
      Don_conv%fre                  = 0.0
      Don_conv%mrmes                = 0.0
      Don_conv%tmes                 = 0.0
      Don_conv%uceml                = 0.0
      Don_conv%umeml                = 0.0
      Don_conv%wmms                 = 0.0
      Don_conv%wmps                 = 0.0
      Don_conv%xice                 = 0.0
      Don_conv%xliq                 = 0.0
      Don_conv%a1                   = 0.0
      Don_conv%amax                 = 0.0
      Don_conv%amos                 = 0.0
      Don_conv%ampta1               = 0.0
      Don_conv%cell_precip          = 0.0
      Don_conv%dcape                = 0.0
      Don_conv%emdi_v               = 0.0
      Don_conv%meso_precip          = 0.0
      Don_conv%pb_v                 = 0.0
      Don_conv%pmd_v                = 0.0
      Don_conv%przm                 = 0.0
      Don_conv%prztm                = 0.0
      Don_conv%pzm_v                = 0.0
      Don_conv%pztm_v               = 0.0
      Don_conv%qtceme               = 0.0
      Don_conv%qtmes1               = 0.0
      Don_conv%qtren1               = 0.0
      Don_conv%temptr               = 0.0
      Don_conv%wtp1                 = 0.0
      Don_conv%wetdepc              = 0.0
      Don_conv%wetdepm              = 0.0
      Don_conv%wetdept              = 0.0

!---------------------------------------------------------------------
!    allocate the components of the donner_cape_type variable Don_cape.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_cape%coin       (isize, jsize) )
      allocate (Don_cape%plcl       (isize, jsize) )
      allocate (Don_cape%plfc       (isize, jsize) )
      allocate (Don_cape%plzb       (isize, jsize) )
      allocate (Don_cape%qint_lag   (isize, jsize) )
      allocate (Don_cape%qint       (isize, jsize) )
      allocate (Don_cape%xcape_lag  (isize, jsize) )
      allocate (Don_cape%xcape      (isize, jsize) )
      if (Nml%do_donner_cape) then
        allocate (Don_cape%cape_p     (isize, jsize, nlev_hires) )
        allocate (Don_cape%env_r      (isize, jsize, nlev_hires) )
        allocate (Don_cape%env_t      (isize, jsize, nlev_hires) )
        allocate (Don_cape%parcel_r   (isize, jsize, nlev_hires) )
        allocate (Don_cape%parcel_t   (isize, jsize, nlev_hires) )
      else
        allocate (Don_cape%cape_p     (isize, jsize, nlev_lsm) )!miz
        allocate (Don_cape%env_r      (isize, jsize, nlev_lsm) )!miz
        allocate (Don_cape%env_t      (isize, jsize, nlev_lsm) )!miz
        allocate (Don_cape%parcel_r   (isize, jsize, nlev_lsm) )!miz
        allocate (Don_cape%parcel_t   (isize, jsize, nlev_lsm) )!miz
      end if
      allocate (Don_cape%model_p    (isize, jsize, nlev_lsm) )
      allocate (Don_cape%model_r    (isize, jsize, nlev_lsm) )
      allocate (Don_cape%model_t    (isize, jsize, nlev_lsm) )

!---------------------------------------------------------------------
!    initialize the components of the donner_cape_type variable Don_cape.
!---------------------------------------------------------------------
      Don_cape%coin        = 0.0
      Don_cape%plcl        = 0.0
      Don_cape%plfc        = 0.0
      Don_cape%plzb        = 0.0   
      Don_cape%qint_lag    = 0.0
      Don_cape%qint        = 0.0
      Don_cape%xcape_lag   = 0.0 
      Don_cape%xcape       = 0.0 
      Don_cape%cape_p      = 0.0
      Don_cape%env_r       = 0.0
      Don_cape%env_t       = 0.0
      Don_cape%parcel_r    = 0.0
      Don_cape%parcel_t    = 0.0
      Don_cape%model_p     = 0.0
      Don_cape%model_r     = 0.0
      Don_cape%model_t     = 0.0

!---------------------------------------------------------------------
!    allocate the components of the donner_rad_type variable Don_rad.
!    definitions of these arrays are found in donner_types.h.
!---------------------------------------------------------------------
      allocate (Don_rad%cell_cloud_frac  (isize, jsize, nlev_lsm) )
      allocate (Don_rad%cell_ice_amt     (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_ice_size    (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_liquid_amt  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_liquid_size (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%cell_droplet_number (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_cloud_frac  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_ice_amt     (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_ice_size    (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_liquid_amt  (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_liquid_size (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%meso_droplet_number (isize, jsize, nlev_lsm ) )
      allocate (Don_rad%nsum             (isize, jsize) )

!---------------------------------------------------------------------
!    initialize the components of the donner_rad_type variable Don_rad 
!    using the input variables supplied.
!---------------------------------------------------------------------
      Don_rad%cell_cloud_frac  = cell_cld_frac
      Don_rad%cell_ice_amt     = cell_ice_amt
      Don_rad%cell_ice_size    = cell_ice_size
      Don_rad%cell_liquid_amt  = cell_liq_amt
      Don_rad%cell_liquid_size = cell_liq_size
      Don_rad%cell_droplet_number = cell_droplet_number
      Don_rad%meso_cloud_frac  = meso_cld_frac
      Don_rad%meso_ice_amt     = meso_ice_amt
      Don_rad%meso_ice_size    = meso_ice_size
      Don_rad%meso_liquid_amt  = meso_liq_amt
      Don_rad%meso_liquid_size = meso_liq_size
      Don_rad%meso_droplet_number = meso_droplet_number
      Don_rad%nsum             = nsum

  if (Nml%do_ensemble_diagnostics) then
!--------------------------------------------------------------------
!    allocate module variables for Donner cumulus ensemble member
!    diagnostics.  These are stored in the derived-type variable
!    Don_cem; see donner_types.h for description of these variables.
!    "ncem" is the number of cumulus ensemble members
!--------------------------------------------------------------------
      ncem = Param%kpar
      allocate ( Don_cem%pfull               (isize, jsize, nlev_lsm ) )
      allocate ( Don_cem%phalf               (isize, jsize, nlev_lsm+1 ) )
      allocate ( Don_cem%zfull               (isize, jsize, nlev_lsm ) )
      allocate ( Don_cem%zhalf               (isize, jsize, nlev_lsm+1 ) )
      allocate ( Don_cem%temp                (isize, jsize, nlev_lsm ) )
      allocate ( Don_cem%mixing_ratio        (isize, jsize, nlev_lsm ) )
      allocate ( Don_cem%cell_precip         (isize, jsize, ncem ) )
      allocate ( Don_cem%pb                  (isize, jsize, ncem ) )
      allocate ( Don_cem%ptma                (isize, jsize, ncem ) )
      allocate ( Don_cem%h1                  (isize, jsize, nlev_lsm, ncem ) )
      if (Nml%do_donner_plume) then
        allocate ( Don_cem%qlw                 (isize, jsize, nlev_hires, ncem ) )
        allocate ( Don_cem%cfracice            (isize, jsize, nlev_hires, ncem ) )
        allocate ( Don_cem%wv                  (isize, jsize, nlev_hires, ncem ) )
        allocate ( Don_cem%rcl                 (isize, jsize, nlev_hires, ncem ) )
      else
        allocate ( Don_cem%qlw                 (isize, jsize, nlev_lsm, ncem ) )
        allocate ( Don_cem%cfracice            (isize, jsize, nlev_lsm, ncem ) )
        allocate ( Don_cem%wv                  (isize, jsize, nlev_lsm, ncem ) )
        allocate ( Don_cem%rcl                 (isize, jsize, nlev_lsm, ncem ) )
      endif
      allocate ( Don_cem%a1                  (isize, jsize ) )
      allocate ( Don_cem%meso_precip         (isize, jsize ) )
      allocate ( Don_cem%cual                (isize, jsize, nlev_lsm ) )
      allocate ( Don_cem%temperature_forcing (isize, jsize, nlev_lsm ) )

!--------------------------------------------------------------------
!    initialize variables for Donner cumulus ensemble member
!    diagnostics.
!--------------------------------------------------------------------
      Don_cem%pfull               = 0.
      Don_cem%phalf               = 0.
      Don_cem%zfull               = 0.
      Don_cem%zhalf               = 0.
      Don_cem%temp                = 0.
      Don_cem%mixing_ratio        = 0.
      Don_cem%cell_precip         = 0.
      Don_cem%meso_precip         = 0.
      Don_cem%pb                  = 0.
      Don_cem%ptma                = 0.
      Don_cem%h1                  = 0.
      Don_cem%qlw                 = 0.
      Don_cem%cfracice            = 0.
      Don_cem%wv                  = 0.
      Don_cem%rcl                 = 0.
      Don_cem%a1                  = 0.
      Don_cem%cual                = 0.
      Don_cem%temperature_forcing = 0.

  endif ! (do_ensemble_diagnostics)

!----------------------------------------------------------------------


end subroutine don_d_init_loc_vars_k



!####################################################################

subroutine don_d_column_input_fields_k  &
         (isize, jsize, nlev_lsm, dt, calc_conv_on_this_step, Col_diag, &
          temp, mixing_ratio, pfull, omega, phalf, parcel_rise, ermesg, error)

use donner_types_mod, only : donner_column_diag_type     

implicit none

!---------------------------------------------------------------------
!    subroutine don_d_column_input_fields_k outputs the 
!    basic profile information for any diagnostic columns.
!---------------------------------------------------------------------

integer,                            intent(in)  :: isize, jsize, nlev_lsm
real,                               intent(in)  :: dt
logical,                            intent(in)  :: calc_conv_on_this_step
type(donner_column_diag_type),      intent(in)  :: Col_diag
real, dimension(isize,jsize,nlev_lsm),                            &
                                    intent(in)  :: temp, mixing_ratio, &
                                                   pfull, omega
real, dimension(isize,jsize,nlev_lsm+1),                            &
                                    intent(in)  :: phalf
real, dimension(isize,jsize),       intent(in)  :: parcel_rise              
character(len=*),                   intent(out) :: ermesg
integer,                            intent(out) :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize         size of x-dimension of physics window
!     jsize         size of y-dimension of physics window
!     nlev_lsm      number of layers in large-scale model
!     dt            physics time step [ sec ]
!     calc_conv_on_this_step
!                   logical indicating whether the deep convection
!                   calculation is to be done on this timestep
!     Col_diag      donner_column_diagtype variable containing the
!                   information defining the columns fro which diagnos-
!                   tics are desired.
!     temp          temperature field at model levels [ deg K ]
!     mixing_ratio  vapor mixing ratio field at model levels 
!                   [ kg(h20) / kg(dry air) ]
!     pfull         pressure field at full-levels 1:nlev_lsm    [ Pa ]
!     omega         model omega field at model full levels 
!                   [ Pa / sec ]
!     phalf         pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     parcel_rise   accumulated vertical displacement of a near-surface 
!                   parcel as a result of the lowest model level omega 
!                   field [ Pa ]
!
!   intent(out) variables:
!
!     ermesg        character string containing any error message
!                   that is returned from a kernel subroutine
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer :: idiag, jdiag, unitdiag
      integer :: n, k       

!----------------------------------------------------------------------
!   local variables:
!
!     idiag         physics window i index of current diagnostic column
!     jdiag         physics window j index of current diagnostic column
!     unitdiag      i/o unit assigned to current diagnostic column
!     n, k          do-loop indices
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    loop over the diagnostic columns in this physics window. output
!    the physics timestep and the window coordinates, and whether the 
!    convection calculation is to be done on this timestep.
!---------------------------------------------------------------------
      do n=1,Col_diag%ncols_in_window
        idiag = Col_diag%i_dc(n)
        jdiag = Col_diag%j_dc(n)
        unitdiag = Col_diag%unit_dc(n)
        write (unitdiag, '(a,f8.1, 2i4)')  &
                ' physics timestep, window i, window j= ',  &
                  dt, idiag, jdiag                         
        write (unitdiag,'(a,l4 )' )  ' conv_calc_on_this_step = ',   &
                  calc_conv_on_this_step

!----------------------------------------------------------------------
!    if the calculation is to be done on this timestep, and convection
!    in the column is not precluded by downward motion at the lowest 
!    level, output the column temperature and mixing ratio profiles
!    over the levels at which output has been requested.
!----------------------------------------------------------------------
        if (calc_conv_on_this_step) then
          if (omega(idiag,jdiag,nlev_lsm) < 0.0) then
            write (unitdiag, '(a)')  & 
                       '                      input profiles'
            write (unitdiag, '(a)')  &
                       '   k   press      temp           mixing ratio '
            write (unitdiag, '(a)')  &
                  '        hPa       deg K    g(h2o) / kg (dry air)  '
            do k=Col_diag%kstart,nlev_lsm
              write (unitdiag, '(i4, 2f10.4, 7x, 1pe13.5)')  &
                   k, 1.0E-02*pfull(idiag,jdiag,k), temp(idiag,jdiag,k),&
                   1.0e03*mixing_ratio(idiag,jdiag,k)
            end do
          endif 
        endif 

!---------------------------------------------------------------------
!    output the surface pressure, omega at the lowest level, and the 
!    accumulated parcel displacement.
!---------------------------------------------------------------------
        write (unitdiag,'(a,f13.4,1pe13.5)')  &
                  ' sfcprs (hPa),  omega_btm (Pa/sec)= ',   &
                  1.0E-02*phalf(idiag,jdiag,nlev_lsm+1),   &
                  omega(idiag,jdiag,nlev_lsm) 
        write (unitdiag,'(a,f13.6)')  ' omint (hPa)= ',   &
                  1.0E-02*parcel_rise(idiag,jdiag)
      end do

!---------------------------------------------------------------------


end subroutine don_d_column_input_fields_k 



!####################################################################

subroutine don_d_convection_driver_k    &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me,  &
          cloud_tracers_present, cbmf, dt, Nml,    &
          Initialized, Param, Col_diag, temp, mixing_ratio,  &
          pfull, zfull, zhalf, pblht, tkemiz, qstar, cush, coldT,  &!miz
          qlin, qiin, qain, lag_cape_temp, lag_cape_vapor,  &
          lag_cape_press, phalf, current_displ, land, sfc_sh_flux,  &
          sfc_vapor_flux, tr_flux, tracers, Don_cape, Don_conv, &
          Don_rad, Don_cem, temperature_forcing, moisture_forcing, &
          total_precip, &
          donner_humidity_factor, donner_humidity_area, dql, dqi, dqa, &
          mhalf_3d, &
          exit_flag, ermesg, error, sd, Uw_p, ac, cp, ct, Don_budgets)

use donner_types_mod, only : donner_initialized_type, donner_rad_type, &
                             donner_param_type, donner_conv_type, &
                             donner_nml_type, donner_budgets_type, &
                             donner_cem_type, &
                             donner_column_diag_type, donner_cape_type

use  conv_utilities_k_mod,only : adicloud, sounding, uw_params
use  conv_plumes_k_mod,only    : cplume, ctend
implicit none

!---------------------------------------------------------------------
!    subroutine don_d_convection_driver_k manages the cal-
!    culation of the effects of deep convection on atmospheric fields 
!    by calling routines to lift a parcel, determine if deep convection 
!    results and, if so, obtain the temperature and moisture forcing and 
!    precipitation produced, and the fields needed to assess the effects
!    of the deep convection on the radiative fluxes and heating and 
!    the large-scale cloud fields of the model.
!---------------------------------------------------------------------

integer,                         intent(in) :: isize, jsize, nlev_lsm,  &
                                               nlev_hires,    ntr, me
logical,                         intent(in) :: cloud_tracers_present
real, dimension(isize,jsize),    intent(inout) :: cbmf
real,                            intent(in) :: dt 
type(donner_nml_type),           intent(inout) :: Nml
type(donner_initialized_type),   intent(inout) :: Initialized
type(donner_param_type),         intent(in) :: Param
type(donner_column_diag_type),   intent(in) :: Col_diag
real,    dimension(isize,jsize,nlev_lsm),              &
                              intent(in)    :: temp, mixing_ratio,  &
                                               pfull, zfull, qlin, &
                                               qiin, qain,   &
                                               lag_cape_temp, &
                                               lag_cape_vapor, &
                                               lag_cape_press
real,    dimension(isize,jsize,nlev_lsm+1),                       &
                              intent(in)    ::  phalf, zhalf                  
real,    dimension(isize,jsize),                                     &
                              intent(in)    :: current_displ, pblht, &
                                               tkemiz, qstar, cush, land, &
                                               sfc_sh_flux,   &
                                               sfc_vapor_flux
logical, dimension(isize,jsize), intent(in) :: coldT
real,    dimension(isize,jsize,ntr),                             &
                              intent(in)    :: tr_flux        
real,    dimension(isize,jsize,nlev_lsm,ntr),                      &
                              intent(in)    :: tracers        
type(donner_cape_type),       intent(inout) :: Don_cape
type(donner_conv_type),       intent(inout) :: Don_conv
type(donner_budgets_type),       intent(inout) :: Don_budgets
type(donner_rad_type),        intent(inout) :: Don_rad
type(donner_cem_type),        intent(inout) :: Don_cem
real,    dimension(isize,jsize,nlev_lsm),                           &
                              intent(out)   :: temperature_forcing,&
                                               moisture_forcing
real,    dimension(isize,jsize),                                  &
                              intent(out)   :: total_precip
real,    dimension(isize,jsize,nlev_lsm),                              &
                              intent(out)   :: donner_humidity_factor, &
                                               donner_humidity_area, &
                                               dql, dqi, dqa
real,    dimension(isize,jsize,nlev_lsm+1),                          &
                              intent(out)   :: mhalf_3d  
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error
logical, dimension(isize,jsize),                                &
                              intent(out)   :: exit_flag
type(sounding),               intent(inout) :: sd
type(uw_params),               intent(inout) :: Uw_p
type(adicloud),               intent(inout) :: ac
type(cplume),                 intent(inout) :: cp
type(ctend),                  intent(inout) :: ct

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!                    tion parameterization
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     lag_cape_temp  temperature field used in lag-time cape 
!                    calculation [ deg K ]
!     lag_cape_vapor vapor mixing ratio field used in lag-time
!                    cape calculation [ kg(h2o) / kg(dry air) ]
!     lag_cape_press model full-level pressure field used in 
!                    lag-time cape calculation  [ Pa ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios of tracers transported by the
!                    donner deep convection parameterization
!                    [ kg(tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape type derived type variable containing 
!                    diagnostics and intermediate results related to the
!                    cape calculation associated with the donner convec-
!     Don_conv       donner_convection_type derived type variable
!                    containing diagnostics and intermediate results 
!                    describing the nature of the convection produced by
!                    the donner parameterization
!     Don_rad        donner_rad_type derived type variable used to hold 
!                    those fields needed to connect the donner deep 
!                    convection parameterization and the model radiation 
!                    package
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!
!   intent(out) variables:
!
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing  
!                    vapor mixing ratio tendency due to donner 
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     dql            tendency of cloud liquid specific humidity
!                    due to donner convection 
!                    [ kg(h2o) / kg(moist air) / sec ]
!     dqi            tendency of cloud ice specific humidity
!                    due to donner convection 
!                    [ kg(h2o) / kg(moist air) / sec ]
!     dqa            tendency of large-scale cloud area
!                    due to donner convection 
!                    [ fraction / sec ]
!     exit_flag      logical array indicating whether deep convection 
!                    exists in a column
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!----------------------------------------------------------------------
!   local variables:
 
     integer :: i, j, k, n
!    real, dimension(isize, jsize, nlev_lsm) :: alpha
     real    :: press0, press1, qlw_wd, qrw_wd, qrw_col_wd, delta_tracer, air_density_wd, tv_wd
 
!---------------------------------------------------------------------
!   local variables:
!
!      i, j, k, n       do-loop indices
!      press0, press1   pressure levels
!      qlw_wd           cloud ice (kg/kg)
!      qrw_wd           mesoscale precip per timestep (kg/m3)
!      qrw_col_wd       mesoscale precip per timestep (kg/kg)
!      delta_tracer     tracer change from mesoscale wet deposition (tracer units)
!      air_density_wd   air density (kg/m3)
!      tv_wd            virtual temperature (K)
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    call don_c_def_conv_env_k to determine how 
!    a parcel will behave with respect to deep convection in each model
!    column.
!---------------------------------------------------------------------
      if (Nml%do_donner_cape) then
        call don_c_def_conv_env_k          &
            (isize, jsize, nlev_lsm, nlev_hires, Nml, Param,  &
             Initialized, Col_diag, &
             temp, mixing_ratio, pfull, lag_cape_temp, lag_cape_vapor,  &
             lag_cape_press, current_displ, cbmf, Don_cape, Don_conv, ermesg, error)
      else
        call don_c_def_conv_env_miz   &
           (isize, jsize, nlev_lsm, ntr, dt, Nml, Param, Initialized, &
            Col_diag, tracers, pblht, tkemiz, qstar, cush, land, coldT,     &
           temp, mixing_ratio, pfull, phalf, zfull, zhalf,   &
           lag_cape_temp, lag_cape_vapor, current_displ, cbmf, Don_cape,  &
           Don_conv, sd, Uw_p, ac)
      endif

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call don_d_cupar to calculate the normalized deep convective 
!    forcing.
!---------------------------------------------------------------------
      call don_d_cupar_k     &
           (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Col_diag, &
            Param, Nml, Initialized, cbmf, current_displ, sfc_sh_flux, &
            sfc_vapor_flux, temp, mixing_ratio, pblht, tkemiz, qstar, &
            cush, land, coldT, pfull, phalf,&
            zfull, zhalf, sd, Uw_p, ac, cp, ct, tr_flux, tracers,    &
            Don_conv, Don_cape, Don_cem, temperature_forcing, &
            moisture_forcing, &
            total_precip, Don_budgets, ermesg, error, exit_flag)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    call define_donner_anvil_ice to define the ice content profile
!    (Don_conv%xice) and the pressures at top and bottom of mesoscale
!    anvil (Don_conv%prztm, Don_conv%przm).
!----------------------------------------------------------------------
      call don_m_define_anvil_ice_k   &
           (isize, jsize, nlev_lsm, Param, Col_diag, pfull, temp,    &
            exit_flag, Don_conv, ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!   Calculate wet deposition in mesoscale updraft here. Use:
!   Don_conv%meso_precip(isize,jsize): grid average (mm/day)
!   Don_conv%temptr(isize,jsize,nlev_lsm,ntr): tracer concentration in
!   mesoscale updraft (kg(tracer)/kg) This quantity is meaningful
!   only for pressure values p in the mesoscale updraft 
!   (pzm <= p <= pztm). It is defined on the pfull pressure surfaces.
!   pfull(isize,jsize,nlev_lsm): pressures (Pa) on full levels
!   phalf(isize,jsize,nlev_lsm+1): pressures (Pa) on half levels
!   Index convention for pressure levels:
!      -----  phalf(1)  (top of GCM grid)
!      -----  pfull(1)
!      -----  phalf(2)
!       ...
!      -----  pfull(k-1)
!      -----  phalf(k)
!      -----  pfull(k)
!      -----  phalf(k+1)
!       ...
!      -----  pfull(nlev_lsm)
!      -----  phalf(nlev_lsm+1)  (bottom of GCM grid)
!   Don_conv%pztm_v(isize,jsize): pressure at top of mesoscale updraft 
!                                 (Pa)
!   Don_conv%pzm_v(isize,jsize): pressure at base of mesoscale updraft 
!                                (Pa)
!   Don_conv%ampta1(isize,jsize): fractional area of mesoscale updraft
!   Don_conv%qtmes1(isize,jsize,nlev_lsm,ntr): grid-average tracer 
!                                              tendency
!   due to mesoscale updraft (kg(tracer)/(kg sec))
!   Don_conv%xice(isize,jsize,nlev_lsm): ice mass mixing ratio in
!   mesocale updfraft  (kg(ice)/kg)
!----------------------------------------------------------------------
      do j=1,jsize
      do i=1,isize
        if (Don_conv%meso_precip(i,j) > 1.e-20 .and. &
            Don_conv%ampta1(i,j) > 1.e-20) then
! convert precip from mm/day to kg/kg/timestep to kg/m3/timestep
           qrw_col_wd = Don_conv%meso_precip(i,j)*dt/Param%seconds_per_day * 1.e-3*Param%dens_h2o ! kg/m2/timestep
           qrw_col_wd = qrw_col_wd * Param%grav / ( Don_conv%pzm_v(i,j)-Don_conv%pztm_v(i,j) )    ! kg/kg/timestep
! convert precip from large-scale average to in-cloud rain amount
           qrw_col_wd = qrw_col_wd / Don_conv%ampta1(i,j)
           do k=1,nlev_lsm
             if (phalf(i,j,k) >= Don_conv%pzm_v(i,j)) then
               exit
             elseif (phalf(i,j,k+1) > Don_conv%pztm_v(i,j)) then
               press0 = MIN(phalf(i,j,k+1),Don_conv%pztm_v(i,j))
               press1 = MAX(phalf(i,j,k),Don_conv%pzm_v(i,j))
               qlw_wd = Don_conv%xice(i,j,k)
! convert precip from kg/kg/timestep to kg/m3/timestep
               tv_wd = temp(i,j,k) &
                       * ( 1 + Param%D608*mixing_ratio(i,j,k)/(1+mixing_ratio(i,j,k)) )
               air_density_wd = 0.5*(press0+press1)/(Param%rdgas*tv_wd)
               qrw_wd = qrw_col_wd * air_density_wd
               do n = 1,size(Don_conv%temptr,4)
                  if (Initialized%wetdep(n)%Lwetdep) then
                     call wet_deposition_0D( Initialized%wetdep(n)%Henry_constant, &
                                             Initialized%wetdep(n)%Henry_variable, &
                                             Initialized%wetdep(n)%frac_in_cloud, &
                                             Initialized%wetdep(n)%alpha_r, &
                                             Initialized%wetdep(n)%alpha_s, &
                                             temp(i,j,k), &
                                             press0, press1, air_density_wd, &
                                             qlw_wd, 0., qrw_wd, &
                                             Don_conv%temptr(i,j,k,n), &
                                             Initialized%wetdep(n)%Lgas, &
                                             Initialized%wetdep(n)%Laerosol, &
                                             Initialized%wetdep(n)%Lice, &
                                             delta_tracer )
                     Don_conv%wetdepm(i,j,k,n) = -Don_conv%ampta1(i,j)* delta_tracer/dt
                     Don_conv%wetdept(i,j,k,n) = Don_conv%wetdept(i,j,k,n) &
                                               + Don_conv%wetdepm(i,j,k,n)
                  end if
               end do
             end if
           end do
        end if
      end do
      end do

!---------------------------------------------------------------------
!    check for tracer realizability, and limit convective tendencies
!    if necessary.
!---------------------------------------------------------------------
      call don_d_check_trc_rlzbility( isize, jsize, nlev_lsm, ntr, dt, &
                                              tracers, Don_conv )

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call donner_rad_driver to define the cloud ice, cloud liquid and
!    cloud areas of the cell and mesoscale clouds associated with 
!    donner convection so as to make them available to the radiation
!    code.
!---------------------------------------------------------------------
      call don_r_donner_rad_driver_k   &
           (isize, jsize, nlev_lsm, Param, Col_diag, Initialized, &
            pfull, temp, land, exit_flag, Don_conv, Don_rad, Nml, ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call donner_lscloud_driver to provide the connection between 
!    donner convection and the large-scale cloud scheme. 
!---------------------------------------------------------------------
      if (Nml%do_donner_lscloud) then
      call don_l_lscloud_driver_k   &
           (isize, jsize, nlev_lsm, cloud_tracers_present, Param,  &
            Col_diag, pfull, temp, exit_flag,  &
            mixing_ratio, qlin, qiin, qain, phalf, Don_conv, &
            donner_humidity_factor, donner_humidity_area, dql, dqi,  &
            dqa, mhalf_3d, ermesg, error)
      else
       call don_l_lscloud_driver_miz   &
            (isize, jsize, nlev_lsm, cloud_tracers_present, Param,  &
             Col_diag, pfull, temp, exit_flag,  &
             mixing_ratio, qlin, qiin, qain, phalf, Don_conv, &
             donner_humidity_factor, donner_humidity_area, dql, dqi,  &
!            dqa,ermesg, error)
             dqa,mhalf_3d, ermesg, error)
      end if

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!!!  QUESTION 3:
!!     A BETTER COMMENT IS NEEDED HERE --
!!     WHY IS THIS REMOVAL NECESSARY ??
!    save the vapor and temperature forcing resulting from the donner
!    deep parameterization. remove from them the contributions due to  
!    vertical transport by the donner mass flux and the mesoscale flux.
!    also remove vapor and temperature tendencies
!    corresponding to these increments from the Donner cumulus
!    thermal forcing and moisture forcing, which included
!    them as evaporatation and/or sublimation in mulsub.
!    assumptions used in strat_cloud_donner_tend to relate detrainment
!    to net mass fluxes differ from those in mulsub, so the
!    increments here do not balance those in mulsub. the difference
!    remains as a phase change.
!    mulsub allowed ice and liquid from convective system to evaporate
!    and/or sublimate as part of thermal and moisture forcing terms
!    remove those tendencies here. different assumptions used to
!    calculate these increments/tendencies here and in mulsub, so
!    some residual phase change will generally remain.
!---------------------------------------------------------------------
      Don_conv%conv_temp_forcing(:,:,:)  = temperature_forcing(:,:,:)
      Don_conv%conv_moist_forcing(:,:,:) = moisture_forcing(:,:,:)

      if (cloud_tracers_present) then
        
        moisture_forcing(:,:,:) = moisture_forcing(:,:,:) - &
                       dql(:,:,:) - dqi(:,:,:)
        temperature_forcing(:,:,:) = temperature_forcing(:,:,:) +   &
                      (dql(:,:,:)*Param%hlv + dqi(:,:,:)*Param%hls)/  &
                                                          (Param%cp_air)
      endif
       
!---------------------------------------------------------------------


end subroutine don_d_convection_driver_k

!######################################################################

subroutine don_d_cupar_k     &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Col_diag, &
          Param, Nml, Initialized, cbmf, current_displ, sfc_sh_flux,   &
          sfc_vapor_flux,                                        &
          temp, mixing_ratio, pblht, tkemiz, qstar, cush, land, coldT, & !miz 
          pfull, phalf, zfull, zhalf,  &
          sd, Uw_p, ac,cp, ct, & !miz
          tr_flux, tracers, &
          Don_conv, Don_cape, Don_cem, temperature_forcing, &
          moisture_forcing,  &
          total_precip, Don_budgets, &
          ermesg, error, exit_flag)

!----------------------------------------------------------------------
!    subroutine cupar drives the parameterization for deep cumulus 
!    convection. it returns the temperature and moisture forcing assoc-
!    iated with deep convection, the total convective precipitation
!    and various diagnostics contained in Don_conv and Don_cape to the 
!    calling routine. it is based on (Donner, 1993, J.Atmos.Sci.).
!---------------------------------------------------------------------

use donner_types_mod, only : donner_initialized_type, donner_nml_type, &
                             donner_param_type, donner_conv_type, &
                             donner_budgets_type, donner_cem_type, &
                             donner_column_diag_type, donner_cape_type
use conv_utilities_k_mod, only : sounding, adicloud, uw_params
use  conv_plumes_k_mod,   only    : cplume, ctend
implicit none

!--------------------------------------------------------------------- 
integer,                           intent(in)    :: isize, jsize,    &
                                                    nlev_lsm,    &
                                                    nlev_hires,&
                                                    ntr, me
real,                              intent(in)    :: dt
type(donner_column_diag_type),     intent(in)    :: Col_diag
type(donner_param_type),           intent(in)    :: Param
type(donner_nml_type),             intent(in)    :: Nml
type(donner_initialized_type),     intent(inout) :: Initialized
type(sounding),                    intent(inout) :: sd
type(uw_params),                   intent(inout) :: Uw_p
type(adicloud),                    intent(inout) :: ac
type(cplume),                      intent(inout) :: cp
type(ctend),                       intent(inout) :: ct
real,    dimension(isize,jsize),   intent(in)    :: pblht, tkemiz, &
                                                    qstar, cush, land
logical, dimension(isize,jsize),   intent(in)    :: coldT

real,    dimension(isize,jsize),   intent(inout)    :: cbmf
real,    dimension(isize,jsize),   intent(in)    :: current_displ, &
                                                    sfc_sh_flux,  &
                                                    sfc_vapor_flux
real,    dimension(isize,jsize,nlev_lsm),                      &
                                   intent(in)    :: pfull, temp, &
                                                    mixing_ratio, zfull
real,    dimension(isize,jsize,nlev_lsm+1),                    &  
                                   intent(in)    :: phalf, zhalf
real,    dimension(isize,jsize,ntr),                           &
                                   intent(in)    :: tr_flux
real,    dimension(isize,jsize,nlev_lsm,ntr),               &
                                   intent(in)    :: tracers
type(donner_conv_type),            intent(inout) :: Don_conv
type(donner_budgets_type),            intent(inout) :: Don_budgets
type(donner_cape_type),            intent(inout) :: Don_cape
type(donner_cem_type),             intent(inout) :: Don_cem
real,    dimension(isize,jsize,nlev_lsm),                      &
                                   intent(out)   :: temperature_forcing,&
                                                    moisture_forcing
real,    dimension(isize,jsize),   intent(out)   :: total_precip
character(len=*),                  intent(out)   :: ermesg
integer,                           intent(out)   :: error
logical, dimension(isize,jsize),   intent(out)   :: exit_flag

!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!                    tion parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Initialized    donner_initialized_type variable containing
!                    variables which are defiuned during initialization.
!                    these values may be changed during model execution.
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     pfull          pressure field at model full levels [ Pa ]
!     temp           temperature field at model full levels [ deg K ]
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     tr_flux        flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer fields that are to be transported by donner
!                    convection [ kg (tracer) / kg (dry air) ]
!
!   intent(out) variables:
!    
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!   intent(inout) variables:
!
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:


      real, dimension (isize, jsize, nlev_lsm, ntr) :: xgcm_v
      real, dimension (isize, jsize, ntr)           :: sfc_tracer_flux
      integer                                       :: i, j, k, n 
      integer                                       :: kcb

!---------------------------------------------------------------------
!   local variables:
!
!      xgcm_v              tracer mixing ratio fields transported by 
!                          donner convection, index 1 nearest surface
!                          [ kg(tracer) / kg (dry air) ]
!      sfc_tracer_flux     tracer flux across the surface
!                          [ kg(tracer) / (m**2 sec) ]
!      i, j, k, n          do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    write a message to the output file for each diagnostic column in 
!    this window.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write (Col_diag%unit_dc(n), '(a, 2i4)')  &
           'cupar: entering cupar with i_dc, j_dc:',    &
                           Col_diag%i_dc(n), Col_diag%j_dc(n)
        end do
      endif

!----------------------------------------------------------------------
!    call donner_deep_check_for_deep_convection_k to determine if deep 
!    convection may at this time be precluded in any of the columns of 
!    this physics window. logical array exit_flag is returned, with a 
!    value of .false. if donner convection is still allowed, a value of 
!    .true. if deep convection is precluded in a particular coluumn.
!----------------------------------------------------------------------
      call don_d_check_for_deep_conv_k   &
           (isize, jsize, nlev_lsm, dt, Param, Nml, Col_diag, &
             Initialized, &
            current_displ, cbmf, Don_cape, Don_conv, exit_flag, ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!-------------------------------------------------------------------
!    for the tracers that are to be transported by donner_deep_mod,
!    define the tracer input fields that will be needed by the donner
!    cloud model.
!-------------------------------------------------------------------
      if (Initialized%do_donner_tracer) then

!-------------------------------------------------------------------
!    define the tracer fluxes across the surface.
!-------------------------------------------------------------------
        sfc_tracer_flux(:,:,:) = tr_flux(:,:,:)

!-------------------------------------------------------------------
!    define an inverted tracer profile (index 1 nearest ground) for use
!    in the cloud and convection routines.
!------------------------------------------------------------------
        do k=1,nlev_lsm
          xgcm_v(:,:,k,:) = tracers(:,:,nlev_lsm-k+1,:)
        end do

!--------------------------------------------------------------------
!    if tracers are not to be transported by donner_deep_mod, define
!    these tracer input fields to be 0.0.
!--------------------------------------------------------------------
      else
        xgcm_v = 0.
        sfc_tracer_flux = 0.0
      endif 

      if (Nml%do_ensemble_diagnostics) then
!--------------------------------------------------------------------
!    save "Don_cem" diagnostics
!--------------------------------------------------------------------
        Don_cem%pfull = pfull
        Don_cem%phalf = phalf
        Don_cem%zfull = zfull
        Don_cem%zhalf = zhalf
        Don_cem%temp = temp
        Don_cem%mixing_ratio = mixing_ratio
      endif
      
!---------------------------------------------------------------------
!    call subroutine mulsub to calculate normalized (in-cloud) cumulus 
!    forcings, one column at a time. the forcings are normalized by the 
!    cloud area at cloud base level a_1(p_b).
!---------------------------------------------------------------------
      call don_d_mulsub_k   &
           (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Param,   &
!++lwh
            Nml, Col_diag,   &
            Initialized,   &
            temp, mixing_ratio, pblht, tkemiz, qstar, cush, cbmf, land, coldT,  &
            phalf, pfull, zhalf, zfull,  &
            sd, Uw_p, ac, cp, ct, &
           sfc_vapor_flux, sfc_sh_flux, &
!--lwh
            sfc_tracer_flux, xgcm_v, Don_cape, Don_conv, Don_cem, exit_flag,  &
            total_precip, temperature_forcing, moisture_forcing, &
            Don_budgets, ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    if using cloud base mass flux calculated by uw_conv_mod in donner
!    deep parameterization closure, modify the cloud fractional area
!    to avoid exceeding this limit.
!---------------------------------------------------------------------
     if (Initialized%using_unified_closure) then
       
!---------------------------------------------------------------------
!    define the cbmf associated with the donner convection.
!---------------------------------------------------------------------
       do j=1,jsize
         do i=1,isize
           if ( .not. exit_flag(i,j)) then
!--------------------------------------------------------------------
!    if there is no cloud base mass flux, there will be no deep 
!    convection in this column.
!--------------------------------------------------------------------
!            if (cbmf(i,j) == 0.0) then
!              exit_flag(i,j) = .true.
!              Don_conv%a1(i,j) = 0.
!              cycle
!            endif

!---------------------------------------------------------------------
!    determine the cloud base level kcb.
!---------------------------------------------------------------------
             kcb = -6
             do k=1, nlev_lsm
               if (Don_conv%uceml(i,j,nlev_lsm-k+1) > 0.0) then
                 kcb = nlev_lsm -k + 1
                 exit
               endif
             end do
             if (kcb == -6) then
               ermesg = 'no cloud base level found'
               error = 1
               return
             endif

!--------------------------------------------------------------------
!    if the cloud base mass flux predicted by donner is more than is
!    available, reduce the donner cloud fraction so that only the
!    available mass flux is used by donner clouds. set the cbmf to be
!    returned and made available for uw shallow convection to 0.0.
!--------------------------------------------------------------------
             if (Don_conv%uceml(i,j,kcb)*Don_conv%a1(i,j) >=    &
                                                       cbmf(i,j)) then
               Don_conv%a1(i,j) = cbmf(i,j)/Don_conv%uceml(i,j,kcb)
               cbmf(i,j) = 0.0

!--------------------------------------------------------------------
!    if the cloud base mass flux predicted by donner is less than what
!    is available, reduce the available cbmf by the amount used by the
!    donner clouds. the remainder will be made available for uw shallow 
!    convection.
!--------------------------------------------------------------------
             else
               cbmf(i,j) = cbmf(i,j) -    &
                               Don_conv%uceml(i,j,kcb)*Don_conv%a1(i,j)
             endif
           endif
         end do
       end do
     endif

!---------------------------------------------------------------------
!    call remove_normalization to remove the normalization from the 
!    deep convection diagnostics and forcing terms by multiplying them 
!    by the fractional cloud area. the values thus become grid-box 
!    averages, rather than averages over the cloudy area, and so are 
!    ready to use in the large-scale model equations. all columns in
!    which exit_flag is .true. are given zero values for total_precip,
!    temperature_forcing and moisture_forcing.
!---------------------------------------------------------------------
      if (Nml%do_donner_plume) then
        call don_d_remove_normalization_k   &
              (isize, jsize, nlev_lsm, ntr, exit_flag, Don_conv, total_precip, &
              Initialized, &
              temperature_forcing, moisture_forcing, ermesg, error)
      else
        call don_d_remove_normalization_miz &
              (isize, jsize, nlev_lsm, ntr, exit_flag, Nml, Don_conv, total_precip, &
               Initialized, &
               temperature_forcing, moisture_forcing, ermesg, error)
      end if

      if (Nml%do_ensemble_diagnostics) then
!
!    save "Don_cem" diagnostics.
!
        Don_cem%a1 = Don_conv%a1
        Don_cem%cual = Don_conv%cual
        Don_cem%temperature_forcing = temperature_forcing
      endif

      do k=1,nlev_lsm
        Don_budgets%liq_prcp(:,:,k) =    &
                    Don_budgets%liq_prcp(:,:,k)*Don_conv%a1(:,:)
        Don_budgets%frz_prcp(:,:,k) =    &
                    Don_budgets%frz_prcp(:,:,k)*Don_conv%a1(:,:)
      end do
      if (Initialized%do_conservation_checks) then
        Don_budgets%vert_motion(:,:) =    &
                     Don_budgets%vert_motion(:,:)*Don_conv%a1(:,:)
        Don_budgets%lheat_precip(:,:) =   &
                    Don_budgets%lheat_precip(:,:)*Don_conv%a1(:,:)
      endif

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    in a diagnostics window, output various desired quantities.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          call don_d_output_cupar_diags_k    &
               (isize, jsize, nlev_lsm, Col_diag, n, exit_flag,   &
                total_precip, temperature_forcing, Don_conv, Don_cape, &
                ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return
        end do
      endif  ! (in_diagnostics_window)

!---------------------------------------------------------------------


end subroutine don_d_cupar_k


!#####################################################################


subroutine don_d_check_for_deep_conv_k   &
           (isize, jsize, nlev_lsm, dt, Param, Nml, Col_diag, &
             Initialized, &
            current_displ, cbmf, Don_cape, Don_conv, exit_flag, ermesg, error)

!---------------------------------------------------------------------
!    subroutine don_d_check_for_deep_conv_k tests for the 
!    sounding- and upward-motion-based criteria which will prevent deep 
!    convection from occurring in a column. if convection is precluded, 
!    the logical variable exit_flag is set to .true. and additional cal-
!    culations in that column will be skipped. if convection is deter-
!    mined to be possible, exit_flag is set to .false., and additional 
!    calculations in the column will be done.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_conv_type, &
                             donner_initialized_type, &
                             donner_column_diag_type, donner_cape_type,&
                             donner_nml_type

implicit none

!---------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize, nlev_lsm
real,                             intent(in)    :: dt
type(donner_param_type),          intent(in)    :: Param
type(donner_nml_type),            intent(in)    :: Nml
type(donner_column_diag_type),    intent(in)    :: Col_diag
type(donner_initialized_type),    intent(in)    :: Initialized
real,    dimension(isize,jsize),  intent(in)    :: current_displ, cbmf
type(donner_cape_type),           intent(inout) :: Don_cape
type(donner_conv_type),           intent(inout) :: Don_conv
logical, dimension (isize,jsize), intent(out)   :: exit_flag 
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!                    tion parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     current_displ  low-level parcel displacement to use in cape
!                    calculation on this step [ Pa ]
!     cbmf
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!
!   intent(out) variables:
!
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
 
      real, dimension (isize, jsize)  ::  pdeet1, pdeet2
      integer                         :: idiag, jdiag, unitdiag
      integer                         :: i, j, k, n 

!---------------------------------------------------------------------
!   local variables:
!
!     pdeet1              pressure depth between the level of free 
!                         convection and the level of zero buoyancy 
!                         [ Pa ]
!     pdeet2              pressure depth between the level of free 
!                         convection and the pressure at lowest 
!                         large-scale model grid level 
!                         [ Pa ]
!     idiag               physics window i index of current diagnostic
!                         column
!     jdiag               physics window j index of current diagnostic
!                         column
!     unitdiag            i/o unit assigned to current diagnostic
!                         column
!     i, j, k, n          do-loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    process each column in the physics window.
!---------------------------------------------------------------------
      do j=1,jsize     
        do i=1,isize     

!---------------------------------------------------------------------
!    define the time rates of change of column convective available 
!    potential energy (Don_conv%dcape_v). 
!---------------------------------------------------------------------
          Don_conv%dcape(i,j) = (Don_cape%xcape(i,j) -          &
                                 Don_cape%xcape_lag(i,j))/dt

!---------------------------------------------------------------------
!    define the pressure depth between the level of free convection
!    and the level of zero buoyancy (pdeet1) and the pressure depth 
!    between the level of free convection and the pressure at lowest 
!    large-scale model grid level (pdeet2).
!---------------------------------------------------------------------
          pdeet1(i,j) = Don_cape%plfc(i,j) - Don_cape%plzb(i,j)
          pdeet2(i,j) = Don_cape%plfc(i,j) - Don_cape%model_p(i,j,1)

!---------------------------------------------------------------------
!    check that all criteria for deep convection are satisfied; if so,
!    set exit_flag to be .false., if any of the criteria are not sat-
!    isfied, set exit_flag to .true. the criteria which can be evaluated
!    at this time are:
!       1) cape (Don_Cape%xcape) must be positive;
!       2) cape must be increasing with time (Don_conv%dcape > 0);
!       3) pressure depth between lfc and lzb (pdeet1) must be greater 
!          than pdeep_cv;
!       4) the time-integrated upward displacement of a parcel from the
!          lowest model level (current_displ) must be sufficient to 
!          allow the parcel to have reached the lfc;
!       5) convective inhibition must be less than cdeep_cv.
!---------------------------------------------------------------------
          if ((Don_cape%xcape(i,j) <= 0.)        .or.  &
              (Don_conv%dcape(i,j) <= 0. .and. Nml%do_dcape)   .or. &
              (pdeet1(i,j) < Param%pdeep_cv.and. Nml%use_pdeep_cv) .or.&
       (Initialized%using_unified_closure .and. cbmf(i,j) == 0.) .or. &
  ((pdeet2(i,j)<current_displ(i,j)) .and. Nml%use_llift_criteria) .or. &
              (Don_cape%coin(i,j) > Param%cdeep_cv) )   then
            exit_flag(i,j) = .true.
          else
            exit_flag(i,j) = .false.
          endif
        end do
      end do

!--------------------------------------------------------------------
!    if in diagnostics window, output info concerning the status of
!    deep convection in this column.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          idiag = Col_diag%i_dc(n)
          jdiag = Col_diag%j_dc(n)
          unitdiag = Col_diag%unit_dc(n)

!--------------------------------------------------------------------
!    for any diagnostic columns in the window for which deep convection
!    is possible, output the integrated upward displacement 
!    (current_displ), the value of cape (Don_Cape%xcape), the time rate
!    of change of cape (Don_conv%dcape) and the logical variable 
!    indicating whether deep convection is precluded in this column 
!    (exit_flag).
!--------------------------------------------------------------------
          if (.not. exit_flag(idiag,jdiag)) then
            write (unitdiag, '(a, f20.12, f20.12, e20.12, l4)')   &
                  'in cupar: omint,cape,dcape, exit_flag',    &
                         current_displ  (idiag,jdiag),   &
                         Don_Cape%xcape (idiag,jdiag),       &
                         Don_conv%dcape (idiag,jdiag),       &
                         exit_flag      (idiag,jdiag)  

!--------------------------------------------------------------------
!    output various thermodynamic parameters (cp_air, cp_vapor, d622,  &
!    rdgas, hlv, rvgas), various sounding levels and features (cape, 
!    convective inhibition, level of zero buoyancy, level of free 
!    convection and the model soundings (p, t, mixing ratio).
!--------------------------------------------------------------------
            write (unitdiag, '(a, 2f12.4)')   &
                   'in cupar: cpi,cpv= ',Param%cp_air, Param%cp_vapor
            write (unitdiag, '(a, 2f12.6, f12.2)')  &
                   'in cupar: rocp,rair,latvap= ',Param%d622, &
                                       Param%rdgas, Param%hlv   
            write (unitdiag, '(a, f12.7)') 'in cupar: rvap= ',Param%rvgas
            write (unitdiag, '(a, 2f14.7, f19.10)')  &
                    'in cupar: cape,cin,plzb= ',  &
                  Don_cape%xcape(idiag,jdiag), &
                  Don_cape%coin(idiag,jdiag), &
                  Don_cape%plzb(idiag,jdiag)
            write (unitdiag, '(a, f19.10)') 'in cupar: plfc= ', &
                  Don_cape%plfc(idiag,jdiag)
            do k=1,nlev_lsm-Col_diag%kstart+1
              write (unitdiag, '(a, i4, f19.10, f20.14, e20.12)') &
                                   'in cupar: k,pr,t,q= ',k,   &
                    Don_cape%model_p(idiag,jdiag,k),   &
                    Don_cape%model_t(idiag,jdiag,k),   &
                    Don_cape%model_r(idiag,jdiag,k)
            end do

!----------------------------------------------------------------------
!    if convection is precluded, output information indicating why.
!----------------------------------------------------------------------
          else
            write (unitdiag, '(a)')   &
               'in cupar: exit_flag is .true., no further calculations&
                              & in this column at this time'
            write (unitdiag, '(a)')   &
                'in cupar: reason(s) for no deep convection:'

!----------------------------------------------------------------------
!    case of no upward motion at lowest level:    
!----------------------------------------------------------------------
            if (current_displ(idiag,jdiag) == 0) then
              write (unitdiag, '(a)')   &
                    'no upward motion at lowest level'
            else 

!----------------------------------------------------------------------
!    case of non-positive cape:    
!----------------------------------------------------------------------
              if (Don_cape%xcape(idiag,jdiag) <= 0.) then      
                write (unitdiag, '(a, f20.12)')   &
                       'non-positive cape, cape = ', &
                           Don_Cape%xcape(idiag,jdiag)      
              endif

!----------------------------------------------------------------------
!    case of non-positive cape time tendency:    
!----------------------------------------------------------------------
              if (Don_conv%dcape(idiag,jdiag) <= 0.) then      
                write (unitdiag, '(a, f20.12)')   &
                      'non-positive cape time tendency, dcape = ', &
                            Don_conv%dcape(idiag,jdiag)      
              endif

              if (Don_cape%plfc(idiag,jdiag) == 0.0 .or.   &
                  Don_cape%plzb(idiag,jdiag) == 0.0) then

!----------------------------------------------------------------------
!    case of sounding not having a level of free convection for 
!    specified parcel:    
!----------------------------------------------------------------------
                if (Don_cape%plfc(idiag,jdiag) == 0.0 ) then
                  write (unitdiag, '(a)')   &
                    'lfc is not definable for parcel used in cape &
                        &calculation'
                endif

!----------------------------------------------------------------------
!    case of sounding not having a level of zero buoyancy for 
!    specified parcel:    
!----------------------------------------------------------------------
                if (Don_cape%plzb(idiag,jdiag) == 0.0) then
                  write (unitdiag, '(a)')   &
                    'lzb is not definable for parcel used in cape &
                        &calculation'
                endif
              else 


!----------------------------------------------------------------------
!    case of sounding not providing a deep enough layer of positive
!    buoyancy:    
!----------------------------------------------------------------------
                if (pdeet1(idiag,jdiag) < Param%pdeep_cv) then      
                  write (unitdiag, '(a, f20.12, a, f20.12,a)')   &
                       'depth of positive buoyancy too shallow, &
                         &plfc - plzb = ',    &
                           pdeet1(idiag,jdiag)*1.0e-02, ' hPa, &
                         & needed depth =', Param%pdeep_cv*1.0e-02, ' hPa'
                endif
                if (Don_cape%plfc(idiag,jdiag) ==  0.0) then 

!----------------------------------------------------------------------
!    case of parcel having insufficient displacement to reach the level
!    of free convection:
!----------------------------------------------------------------------
                else if        &
                  (pdeet2(idiag,jdiag) < current_displ(idiag,jdiag)) then
                  write (unitdiag, '(a, f20.12, a, f20.12, a)')   &
                      'parcel displacement insufficient to reach lfc, &
                       &displacement =',    &
                           current_displ(idiag,jdiag)*1.0e-02, ' hPa, &
                       &needed displacement = ',  &
                            pdeet2(idiag,jdiag)*1.0e-02, ' hPa'
                endif
              endif

!----------------------------------------------------------------------
!    case of sounding having too much convective inhibition:
!----------------------------------------------------------------------
              if (Don_cape%coin(idiag,jdiag) > Param%cdeep_cv) then      
                write (unitdiag, '(a, f20.12, a, f20.12)')   &
                       'convective inhibition too large, cin   = ', &
                             Don_cape%coin(idiag,jdiag), &
                            'max allowable =', Param%cdeep_cv
              endif
            endif
          endif  ! (not exit_flag)
        end do
      endif

!--------------------------------------------------------------------


end subroutine don_d_check_for_deep_conv_k



!#####################################################################

subroutine don_d_mulsub_k   &
         (isize, jsize, nlev_lsm, nlev_hires, ntr, me, dt, Param, Nml, &
!++lwh
          Col_diag, Initialized,   &
          temp, mixing_ratio, pblht, tkemiz, qstar, cush, cbmf, land, coldT, & !miz
          phalf, pfull, zhalf, zfull, sd,  &
           Uw_p, ac,         cp, ct, &
          sfc_vapor_flux, sfc_sh_flux, &
!--lwh
          sfc_tracer_flux, xgcm_v, Don_cape, Don_conv, Don_cem, exit_flag, &
          total_precip, temperature_forcing, moisture_forcing,  &
          Don_budgets, ermesg, error)

!--------------------------------------------------------------------
!    subroutine don_d_mulsub_k calculates the thermal and moisture
!    forcing produced by an ensemble of cumulus elements and any meso-
!    scale circulation which the ensemble induces, following Donner 
!    (1993, JAS). See also LJD notes, "Cu Closure A," 2/97. calculations
!    at and below this subroutine level are done a column at a time, in 
!    only those columns for which the possibility of deep convection has
!    not yet been ruled out.
!
!                L. Donner  GFDL 27 Apr 97
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_conv_type, &
                             donner_nml_type, donner_column_diag_type, &
                             donner_budgets_type, donner_cem_type, &
!++lwh
                             donner_cape_type, donner_initialized_type
!--lwh
use  conv_utilities_k_mod,only : pack_sd_lsm_k, extend_sd_k,  &
                                adicloud, sounding, uw_params
use  conv_plumes_k_mod,only    : cplume, ctend

implicit none

!---------------------------------------------------------------------
integer,                      intent(in)     ::  isize, jsize, nlev_lsm,&
                                                 nlev_hires, ntr, me
real,                         intent(in)     ::  dt
type(donner_param_type),      intent(in)     ::  Param
type(donner_nml_type),        intent(in)     ::  Nml  
type(donner_column_diag_type),                           &
                              intent(in)     ::  Col_diag
!++lwh
type(donner_initialized_type), intent(in)    :: Initialized
!--lwh
real,    dimension(isize,jsize,nlev_lsm+1),                    &
                              intent(in)     ::  phalf, zhalf
real,    dimension(isize,jsize,nlev_lsm),                      &
                              intent(in)     ::  pfull, zfull, temp, mixing_ratio
type(sounding),               intent(inout)  ::  sd
type(uw_params),               intent(inout)  ::  Uw_p
type(adicloud),               intent(inout)  ::  ac
type(cplume),                 intent(inout)  ::  cp
type(ctend),                  intent(inout)  ::  ct
real, dimension(isize,jsize), intent(in)     ::  pblht, tkemiz, qstar, &
                                                 cush, cbmf, land
logical, dimension(isize,jsize), intent(in)  ::  coldT
real,    dimension(isize,jsize),                                &
                              intent(in)     ::  sfc_vapor_flux,  &
                                                 sfc_sh_flux
real,    dimension(isize,jsize,ntr),                              &
                              intent(in)     ::  sfc_tracer_flux       
real,    dimension(isize,jsize,nlev_lsm,ntr),               &
                              intent(in)     ::  xgcm_v
type(donner_cape_type),       intent(inout)  ::  Don_cape
type(donner_conv_type),       intent(inout)  ::  Don_conv
type(donner_cem_type),        intent(inout)  ::  Don_cem
type(donner_budgets_type),       intent(inout)  ::  Don_budgets
logical, dimension(isize,jsize),                            &
                              intent(inout)  ::  exit_flag
real,    dimension(isize,jsize),                          &
                              intent(out)    ::  total_precip        
real,    dimension(isize,jsize,nlev_lsm),                        &
                              intent(out)    ::  temperature_forcing, &
                                                 moisture_forcing
character(len=*),             intent(out)    ::  ermesg
integer,                      intent(out)    ::  error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     dt             physics time step [ sec ]
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     phalf          pressure field at half-levels 1:nlev_lsm+1  [ Pa ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_tracer_flux 
!                    flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     xgcm_v         tracer fields that are to be transported by donner
!                    convection. index 1 nearest the ground.
!                    [ kg (tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     Don_cape       donner_cape_type derived type variable containing
!                    fields associated with the calculation of
!                    convective available potential energy (cape).
!     Don_cem        donner_cem_type derived type variable containing
!                    Donner cumulus ensemble member diagnostics
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!   intent(out) variables:
!    
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     ermesg         character string containing any error message
!                    that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!
!     On Output:
!     
!     ampt             mesoscale cloud fraction, normalized by a(1,p_b)
!     contot           ratio of convective to total precipitation
!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     cual(nlev)       cloud fraction, cells+meso, normalized by a(1,p_b)
!     cuq(nlev)        ice content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     cuqll(nlev)      liquid content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     ecds(nlev)       normalized convective downdraft evaporation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     eces(nlev)       normalzed convective-updraft evporation/sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     mrmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre_v           total normalized precipitation (mm/day)
!     detmfl(nlev)     normalized detrained mass flux from cell
!                      updrafts (kg/((m**2)*s)
!                      (index 1 at atmosphere bottom)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------



!!  UNITS
!!     ensmbl_anvil_cond  ! [mm / day ]
!!    ucemh  [kg /sec / m**2 ]
!!    detmfh [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]


!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        current_displ  integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb     pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!
!     On Input as Parameters:
!
!        kmax     number of vertical levels at Skyhi resolution
!        kpar     number of cumulus sub-ensembles
!        ncap     number of vertical levels in Cape.F resolution
!




!      disa_v              thermal forcing due to deep convection
!                          index 1 nearest surface, normalized by 
!                          cloud area  [ deg K / sec ]
!      dise_v              moisture forcing due to deep convection
!                          index 1 nearest surface, normalized by 
!                          cloud area  [ kg(h2o) / (kg(dry air) *sec ) ]

!----------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_lsm,Don_budgets%n_water_budget) ::  wat_budg
      real, dimension(nlev_lsm,Don_budgets%n_enthalpy_budget) :: &
                                                               ent_budg
      real, dimension(nlev_lsm,Don_budgets%n_precip_paths,   &
                               Don_budgets%n_precip_types)  :: &
                                                               prc_budg
      real, dimension(nlev_lsm)               ::        &
              ensmbl_cloud_area, cutotal, cmus_tot, cuq, cuql_v, disa, &
              disb, disd, disv, dise, dmeml, uceml, umeml, &
              ecds_liq, ecds_ice, eces_liq, eces_ice,  &
              disc_liq, disc_ice, dism_liq, dism_liq_frz, &
              dism_liq_remelt, dism_ice, dism_ice_melted, &
              disp_liq, disp_ice, disz, disz_remelt, disp_melted, &
              disze1, disze2, disze3,                               &
              emds_liq, emds_ice, emes_liq, emes_ice, &
              mrmes, mrmes_up, mrmes_dn, tmes, tmes_up, tmes_dn, &
              wmms, wmps, detmfl, meso_cloud_area, disf,            &
              disn, enctf, encmf, disg_liq, disg_ice, &
              enev, ensmbl_melt, ensmbl_melt_meso, anvil_precip_melt, &
              ensmbl_freeze, ensmbl_freeze_meso, temp_tend_melt,  &
              liq_prcp, frz_prcp

      real, dimension(isize, jsize, nlev_lsm) :: disa_v, dise_v
      real, dimension (nlev_lsm)              :: rlsm_miz, emsm_miz, &
                                                 cld_press_miz
      real, dimension (nlev_hires)            :: rlsm, emsm, cld_press
      real, dimension( nlev_lsm,ntr)          :: ensmbl_wetc
      real, dimension( nlev_lsm,ntr)          :: qtmes, qtren, wtp, &
                                                 temptr
      real, dimension (nlev_hires,ntr)        :: etsm
      real, dimension (nlev_lsm+1)            :: phalf_c               

      real, dimension (nlev_lsm)       :: model_tx, model_rx, model_px
      real, dimension (nlev_hires)     :: env_r, env_t, cape_p, &
                                          parcel_r, parcel_t
      real         ::  lofactor
      real         ::  ampta1, ensmbl_cond, pb, ensmbl_precip,  &
                       pt_ens, max_depletion_rate, dqls_v,  &
                       ensmbl_anvil_cond_liq, ensmbl_anvil_cond_ice, &
                       ensmbl_anvil_cond_liq_frz, &
                       qlsd_v, frz_frac, lprcp, vrt_mot
      logical      ::  meso_frz_intg_sum, melting_in_cloud, lmeso, &
                                      debug_ijt
      integer      ::  diag_unit
      integer      ::  kinv,  kcont
      integer      ::  i, j, k, n, kk


!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!----------------------------------------------------------------------
!    initialize the output arrays.
!----------------------------------------------------------------------
      temperature_forcing = 0.
      moisture_forcing    = 0.
      total_precip        = 0.

!---------------------------------------------------------------------
!    output a message to all diagnostic files indicating entry into
!    subroutine don_d_mulsub_k.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          write(Col_diag%unit_dc(n), '(a, 2i4)')    &
            'in mulsub: i_dc,j_dc= ', Col_diag%i_dc(n), Col_diag%j_dc(n)
        end do
      endif

!--------------------------------------------------------------------
!    LOOP OVER COLUMNS IN CURRENT PHYSICS WINDOW:
!--------------------------------------------------------------------
      do j=1,jsize               
        do i=1,isize                 

!--------------------------------------------------------------------
!    if it is already known that convection is not possible in this 
!    column, cycle to end of this loop and process the next column.
!--------------------------------------------------------------------
          if (exit_flag(i,j)) cycle

!--------------------------------------------------------------------
!    determine if column diagnostics are requested for this column.
!    define the output unit and set debug_ijt to .true. if it is.
!--------------------------------------------------------------------
          debug_ijt = .false.
          diag_unit = -99
          if (Col_diag%in_diagnostics_window ) then
            do n=1,Col_diag%ncols_in_window
              if (j == Col_diag%j_dc(n) .and.      &
                  i == Col_diag%i_dc(n)) then
                debug_ijt = .true.
                diag_unit = Col_diag%unit_dc(n)
                exit
              endif
            end do
          endif

!---------------------------------------------------------------------
!    define an inverted interface level pressure profile phalf_c 
!    (level 1 at the surface).
!---------------------------------------------------------------------
          do k=1,nlev_lsm+1
            phalf_c(k) = phalf(i,j,nlev_lsm+2-k)
          end do

          if (.not. Nml%do_donner_closure .or. &
              .not. Nml%do_donner_plume)  then
            call pack_sd_lsm_k (Nml%do_lands, land(i,j), coldT(i,j), &
                                dt, pfull(i,j,:), phalf(i,j,:),  &
                                zfull(i,j,:), zhalf(i,j,:), &
                                temp(i,j,:), mixing_ratio(i,j,:),  &
                                xgcm_v(i,j,:,:), sd)
            call extend_sd_k (sd, pblht(i,j), .false., Uw_p) 
          endif

!--------------------------------------------------------------------
!    define factor to modify entrainment coefficients if option is
!    activated.
!--------------------------------------------------------------------
   if (Nml%do_donner_closure .or. Nml%do_donner_plume) then
     if (Nml%lochoice == 0) then
       lofactor = 1. - land(i,j) *(1.0 - Nml%lofactor0)
     else if (Nml%lochoice == 1) then
       lofactor = Nml%pblht0/max (pblht(i,j), Nml%pblht0)
     else if (Nml%lochoice == 2) then
       lofactor = Nml%tke0/max (tkemiz(i,j), Nml%tke0)
     else if (Nml%lochoice == 3) then
       lofactor = Nml%tke0/max (tkemiz(i,j), Nml%tke0)
       lofactor = sqrt (lofactor)
     else
       lofactor = 1.0
     endif
   else
     lofactor = 1.0
   endif

!--------------------------------------------------------------------
!    call don_d_integ_cu_ensemble_k to determine the 
!    characteristics of the clouds in the cumulus ensemble defined in 
!    the current column.
!--------------------------------------------------------------------
          if (Nml%do_donner_plume) then
            call don_d_integ_cu_ensemble_k             &
                (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, &
                 lofactor, Param, Col_diag, Nml, Initialized,  &
                 Don_cape%model_t(i,j,:), Don_cape%model_r(i,j,:), &
                 Don_cape%model_p(i,j,:), phalf_c, xgcm_v(i,j,:,:), &
                 sfc_sh_flux(i,j), sfc_vapor_flux(i,j), &
                 sfc_tracer_flux(i,j,:), Don_cape%plzb(i,j), &
                 exit_flag(i,j),  ensmbl_precip, ensmbl_cond,       &
                 ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
                 ensmbl_anvil_cond_ice, pb,  &
                 pt_ens, ampta1, Don_conv%amax(i,j), emsm, rlsm,  &
                 cld_press, ensmbl_melt, ensmbl_melt_meso, &
                 ensmbl_freeze, ensmbl_freeze_meso, ensmbl_wetc, &
                 disb, disc_liq, disc_ice, dism_liq, dism_liq_frz, &
                 dism_liq_remelt, dism_ice, dism_ice_melted, &
                 disp_liq, disp_ice, disz, disz_remelt, disp_melted, &
                 disze1, disze2, disze3,                           &
                 disd, disv, disg_liq, disg_ice, enctf, encmf, enev,  &
                 ecds_liq, ecds_ice, eces_liq, eces_ice, &
                 ensmbl_cloud_area, cuq, cuql_v, detmfl, uceml, &
                 qtren, etsm, lmeso,                frz_frac,&
                 meso_frz_intg_sum, ermesg, error, melting_in_cloud, &
                 i, j, Don_cem)
          else
            call don_d_integ_cu_ensemble_miz             &
                (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, &
                 Param, Col_diag, Nml, Initialized, &
                 Don_cape%model_t(i,j, :), Don_cape%model_r(i,j,:), &
                 Don_cape%model_p(i,j,:), phalf_c, pblht(i,j),  &
                 tkemiz(i,j), qstar(i,j), cush(i,j), land(i,j), coldT(i,j), &
                 dt, sd, Uw_p, ac, cp, ct,  &
                 xgcm_v(i,j,:,:), sfc_sh_flux(i,j),                 &
                 sfc_vapor_flux(i,j), sfc_tracer_flux(i,j,:),   &
                 Don_cape%plzb(i,j), exit_flag(i,j),                &
                 ensmbl_precip, ensmbl_cond,                      &
                 ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
                 ensmbl_anvil_cond_ice, pb,  &
                 pt_ens, ampta1, Don_conv%amax(i,j), emsm_miz, &
                 rlsm_miz, cld_press_miz, ensmbl_melt,  &
                 ensmbl_melt_meso, ensmbl_freeze, ensmbl_freeze_meso, &
                 ensmbl_wetc, disb, disc_liq, disc_ice, dism_liq, &
                 dism_liq_frz, dism_liq_remelt, dism_ice, &
                 dism_ice_melted, disp_liq, disp_ice, disz, &
                 disz_remelt, disp_melted, disze1, disze2, disze3, &
                 disd, disv, disg_liq, disg_ice, &
                 enctf, encmf, enev, ecds_liq, ecds_ice,   &
                 eces_liq, eces_ice, ensmbl_cloud_area,  &
                 cuq, cuql_v, detmfl, uceml, qtren, etsm, lmeso, &
                 frz_frac, meso_frz_intg_sum, ermesg, error, melting_in_cloud, &
                 i, j, Don_cem)
          endif

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!--------------------------------------------------------------------
!    if the exit_flag was set within integrate_cumulus_ensemble (due to
!    an ensemble member either a) not reaching an acceptable level of 
!    free convection, b) not producing precipitation, c) having con-
!    densate evaporation within the cloud, or d) not having a net column
!    non-zero moisture forcing (the "moisture constraint") stop the 
!    calculations for this column -- deep convection is turned off here,
!    output fields will reflect the absence of the effects of deep 
!    convection in this column.
!--------------------------------------------------------------------
          if (exit_flag(i,j)) cycle

!--------------------------------------------------------------------
!    if mesoscale circulation is present, call subroutine meso_effects 
!    to obtain full ensemble output fields to be applied to large-scale
!    model fields.
!--------------------------------------------------------------------
          if (lmeso) then
            if (Nml%do_donner_plume) then
              call don_m_meso_effects_k  &
                 (me, nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, &
                  Param, Nml, Don_cape%model_p(i,j,:),   &
                  Don_cape%model_t(i,j,:), Don_cape%model_r(i,j,:),  &
                  phalf_c, rlsm, emsm, etsm, xgcm_v(i,j,:,:),   &
                  ensmbl_cond, ensmbl_precip, pb, Don_cape%plzb(i,j), &
                  pt_ens, ampta1,                     &
                  ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
                  ensmbl_anvil_cond_ice,  &
                  wtp, qtmes,  meso_frz_intg_sum,   &
                  anvil_precip_melt, meso_cloud_area, cmus_tot, dmeml, &
                  emds_liq, emds_ice, emes_liq, emes_ice, &
                  wmms, wmps, umeml, temptr, tmes,tmes_up,  &
                  tmes_dn,  mrmes, mrmes_up, mrmes_dn,  &
                  Don_conv%emdi_v(i,j), Don_conv%pmd_v(i,j),   &
                  Don_conv%pztm_v(i,j), Don_conv%pzm_v(i,j),    &
                  Don_conv%meso_precip(i,j), ermesg, error)
            else
              call don_m_meso_effects_miz  &
                 (me, nlev_lsm, nlev_hires, ntr, diag_unit, debug_ijt, &
                  Param, Nml, Don_cape%model_p(i,j,:),   &
                  Don_cape%model_t(i,j,:), Don_cape%model_r(i,j,:),  &
                  phalf_c, rlsm_miz, emsm_miz, etsm, xgcm_v(i,j,:,:),  &
                  ensmbl_cond, ensmbl_precip, pb, Don_cape%plzb(i,j), &
                  pt_ens, ampta1,                     &
                  ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
                  ensmbl_anvil_cond_ice,  &
                  wtp, qtmes, meso_frz_intg_sum,    &
                  anvil_precip_melt, meso_cloud_area, cmus_tot, dmeml, &
                  emds_liq, emds_ice, emes_liq, emes_ice, &
                  wmms, wmps, umeml, temptr, tmes, tmes_up, tmes_dn, &
                  mrmes, mrmes_up, mrmes_dn, Don_conv%emdi_v(i,j),    &
                  Don_conv%pmd_v(i,j), Don_conv%pztm_v(i,j), &
                  Don_conv%pzm_v(i,j), Don_conv%meso_precip(i,j), &
                  ermesg, error)
            endif


!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return

!--------------------------------------------------------------------
!    define cmus_tot   as the profile of total condensate source to the
!    large-scale flow from the mesoscale circulation; the sum of the 
!    water mass condensed in the mesoscale updraft plus the vapor 
!    transferred from cell to mesoscale and then condensed. 
!--------------------------------------------------------------------
          else
            qtmes = 0.
            wtp = 0.
            umeml = 0.
            dmeml = 0.
            cmus_tot = 0.
            tmes = 0.
            tmes_up = 0.
            tmes_dn = 0.
            wmms = 0.
            wmps = 0.
            mrmes = 0.
            mrmes_up = 0.
            mrmes_dn = 0.
            emds_liq = 0.
            emds_ice = 0.
            emes_liq = 0.
            emes_ice = 0.
            anvil_precip_melt = 0.
            meso_cloud_area = 0.
            meso_frz_intg_sum = .false.
          endif

        if (Nml%do_ensemble_diagnostics) then
           Don_cem%meso_precip = Don_conv%meso_precip
        endif

!---------------------------------------------------------------------
!    if in a diagnostics column, output the profiles of cell-scale 
!    tracer flux convergence (qtren). 
!---------------------------------------------------------------------
          if (debug_ijt) then
            do k=1,nlev_lsm
              do kcont=1,ntr  
                if (qtren(k,kcont) /= 0.00) then
                  write (diag_unit, '(a, 2i4, f19.10, e20.12)')  &
                  'in mulsub: jk, pr,qtren= ', k, kcont,              &
                            Don_cape%model_p(i,j,k), qtren(k,kcont)
                endif
              end do
            end do
          endif

!--------------------------------------------------------------------
!    if in diagnostics column, output the rate of condensate transfer 
!    from cells to anvil (ensmbl_anvil_cond), and the ratio of
!    convective precipitation to total precipitation (contotxx_v).
!--------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a,e20.12, a, e20.12)')  &
              'in mulsub: CATOT= ',ensmbl_anvil_cond_liq + &
               ensmbl_anvil_cond_liq_frz + ensmbl_anvil_cond_ice,   &
                               ' contot=',  &
                       ensmbl_precip/(ensmbl_precip +    &
                                              Don_conv%meso_precip(i,j))
          endif

!----------------------------------------------------------------------
!    call subroutine define_convective_forcing to combine the cell and
!    mesoscale contributions to the output fields and the time tendency
!    terms that will be returned to the large-scale model. it also
!    call subroutine output_diagnostic_profiles to print various 
!    output fields from the donner_deep parameterization in those 
!    columns for which diagnostics have been requested.
!----------------------------------------------------------------------
          call don_d_def_conv_forcing_k   &
              (nlev_lsm, diag_unit, debug_ijt, lmeso, Initialized,  &
               pb, Param, Nml, ensmbl_precip, &
               Don_conv%meso_precip(i,j), meso_cloud_area, &
               anvil_precip_melt, phalf_c, enev,  encmf, ensmbl_freeze,&
               ensmbl_freeze_meso, enctf, disg_liq, disg_ice, &
               ecds_liq, ecds_ice, eces_liq, eces_ice,  &
               emds_liq, emds_ice, emes_liq, emes_ice, mrmes, mrmes_up,&
               mrmes_dn, tmes, tmes_up, tmes_dn, wmps, &
               ensmbl_cloud_area, ensmbl_melt, ensmbl_melt_meso,&
               Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:),  &
               cmus_tot, wmms, disc_liq, disc_ice, dism_liq, &
               dism_liq_frz, dism_liq_remelt, dism_ice, &
               dism_ice_melted, meso_frz_intg_sum, &
               disp_liq, disp_ice, disb, disd, disv, total_precip(i,j),&
               disz, disz_remelt, disp_melted, disze1,disze2, disze3,&
                                         disf, disn, dise, disa, &
               cutotal, temp_tend_melt, lprcp, liq_prcp, frz_prcp, &
               vrt_mot, wat_budg, &
               Don_budgets%n_water_budget, ent_budg,   &
               Don_budgets%n_enthalpy_budget, prc_budg, &
               Don_budgets%n_precip_paths, Don_budgets%n_precip_types, &
               ermesg, error, melting_in_cloud)

          do k=1, nlev_lsm
            kk = nlev_lsm - k + 1
            Don_budgets%liq_prcp(i,j,k) = liq_prcp(kk)
            Don_budgets%frz_prcp(i,j,k) = frz_prcp(kk)
          end do
          if (Initialized%do_conservation_checks) then
            Don_budgets%lheat_precip(i,j) = lprcp
            Don_budgets%vert_motion(i,j)  = vrt_mot
          
          endif
          if (Initialized%do_conservation_checks .or.   &
                                          Nml%do_budget_analysis) then
            Don_budgets%water_budget(i,j,:,:) = wat_budg
            Don_budgets%enthalpy_budget(i,j,:,:) = ent_budg
            Don_budgets%precip_budget(i,j,:,:,:) = prc_budg
          endif

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!----------------------------------------------------------------------
!    call finalize_output_fields to convert to mks units and then store
!    profile arrays into the various components of the donner_conv type
!    derived-type variable Don_conv.
!----------------------------------------------------------------------
          call don_d_finalize_output_fields_k  &
               (nlev_lsm, ntr, i, j, Param, disb, disc_liq, disc_ice, &
                ensmbl_freeze, ensmbl_freeze_meso, &
                temp_tend_melt,  tmes, disd, cmus_tot, &
                ecds_liq, ecds_ice, eces_liq, eces_ice, emds_liq, &
                emds_ice, emes_liq, emes_ice, wmms, wmps, mrmes, &
                cutotal, dmeml, detmfl, temptr, uceml, umeml, cuq, &
                cuql_v, qtren, qtmes, wtp, ensmbl_wetc, Don_conv, &
                ermesg, error)

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
          if (error /= 0 ) return

!--------------------------------------------------------------------
!    store some additional output fields in the donner_conv type 
!    variable for later use.
!--------------------------------------------------------------------
          Don_conv%cell_precip(i,j) = ensmbl_precip
          Don_conv%pb_v(i,j) = pb
          Don_conv%ampta1(i,j) = ampta1
          dise_v(i,j,:) = dise(:)/(1.0E03*Param%SECONDS_PER_DAY)
          disa_v(i,j,:) = disa(:)/Param%SECONDS_PER_DAY

          do k=1,nlev_lsm
            kinv = nlev_lsm + 1 - k
            temperature_forcing(i,j,kinv) = disa(k)/ &
                                                Param%SECONDS_PER_DAY
            moisture_forcing(i,j,kinv) = dise(k)/    &
                                          (1.0e03*Param%SECONDS_PER_DAY)
          end do

!--------------------------------------------------------------------
!    for any diagnostic columns in the window in which deep convection
!    occurred, output the cloud anvil area (Don_conv%ampta1) and the
!    total precipitation produced (total_precip). also output the vert-
!    ical profile of total cloud fraction (Don_conv%cual).
!--------------------------------------------------------------------
          if (debug_ijt) then
            write  (diag_unit, '(a, 2e20.12)')   &
                  'in cupar:  ampt,tpre= ',  &
                            Don_conv%ampta1(i,j), total_precip(i,j)      
            do k=1,nlev_lsm-Col_diag%kstart+1    
              write (diag_unit, '(a, i4, e20.12)')  &
                   'in cupar: k,cual= ',k,  &
                                Don_conv%cual(i,j,nlev_lsm-k+1)
            end do
          endif

!---------------------------------------------------------------------
!    define the time rates of change of column-integrated water vapor
!    (dqls_v) and the time rate of change needed to deplete the column
!    water vapor in a single donner timestep (qlsd_v).
!---------------------------------------------------------------------
          dqls_v = (Don_cape%qint(i,j) - Don_cape%qint_lag(i,j))/dt
          qlsd_v = Don_cape%qint(i,j)/Nml%donner_deep_freq
          max_depletion_rate = dqls_v + qlsd_v

!--------------------------------------------------------------------
!    if in a diagnostic column, output these moisture tendency 
!    variables.
!--------------------------------------------------------------------
          if (debug_ijt) then
            write (diag_unit, '(a, 2e20.12)')   &
                  'in cupar: dqls,qlsd= ', dqls_v, qlsd_v     
          endif

!---------------------------------------------------------------------
!    call determine_cloud_area to define the cloud area of the convect-
!    ive clouds and so close the parameterization. note that exit_flag
!    may be set to .true. within determine_cloud_area, so that the 
!    if (exit_flag) loop must be closed after this call.
!---------------------------------------------------------------------
          if (.not. exit_flag(i,j)) then
            if (Nml%do_donner_closure) then
! this path used by donner_full parameterization:
              call don_d_determine_cloud_area_k  &
                (me, nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param,&
                 Initialized, &
                 Nml, lofactor, max_depletion_rate, Don_conv%dcape(i,j),   &
                 Don_conv%amax(i,j), dise_v(i,j,:), disa_v(i,j,:),    &
                 Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:), &
                 Don_cape%model_r(i,j,:), Don_cape%env_t(i,j,:), &
                 Don_cape%env_r(i,j,:), Don_cape%parcel_t(i,j,:), &
                 Don_cape%parcel_r(i,j,:), Don_cape%cape_p(i,j,:), &
                 exit_flag(i,j), Don_conv%amos(i,j), Don_conv%a1(i,j),&
                 ermesg, error)
            else  ! (do_donner_closure)
! these paths used by donner_lite parameterization:

              if (Nml%do_donner_cape) then
! if not do_donner_closure but do_donner_cape, then previous parcel cape
!  calculation will have used hires vertical grid.
                 call don_d_determine_cloud_area_miz  &
                (me, nlev_lsm, ntr, dt, nlev_hires, diag_unit,&
                 debug_ijt, Param, Initialized, Nml, xgcm_v(i,j,:,:), &
                 pfull(i,j,:), zfull(i,j,:), phalf(i,j,:),  &
                 zhalf(i,j,:), pblht(i,j), tkemiz(i,j), qstar(i,j), &
                 cush(i,j), cbmf(i,j), land(i,j),  coldT(i,j), sd, Uw_p, ac, &
                 max_depletion_rate, Don_cape%xcape(i,j), &
                 Don_conv%dcape(i,j),   &
                 Don_conv%amax(i,j), dise_v(i,j,:), disa_v(i,j,:),    &
                 Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:), &
                 Don_cape%model_r(i,j,:), Don_cape%env_t(i,j,:), &
                 Don_cape%env_r(i,j,:), Don_cape%parcel_t(i,j,:), &
                 Don_cape%parcel_r(i,j,:), Don_cape%cape_p(i,j,:), &
                 exit_flag(i,j), Don_conv%amos(i,j), Don_conv%a1(i,j),&
                 ermesg, error)
              else ! (do_donner_cape)
                if (Nml%do_hires_cape_for_closure) then
!  if not do_donner_cape (lo res cape calc for convection), but desire 
!  to use hires cape calc for closure:

!--------------------------------------------------------------------
!    call generate_cape_sounding to produce a high-resolution atmos-
!    pheric sounding to be used to evaluate cape.
!--------------------------------------------------------------------
                   call don_c_generate_cape_sounding_k &
                        (nlev_lsm, nlev_hires, temp(i,j,:),   &
                         mixing_ratio(i,j,:), pfull(i,j,:),   &
                         model_tx, model_rx, model_px, cape_p, &
                         env_t, env_r, ermesg,  error)

                   call don_d_determine_cloud_area_miz  &
                (me, nlev_lsm, ntr, dt, nlev_hires, diag_unit,&
                 debug_ijt, Param, Initialized, Nml, xgcm_v(i,j,:,:), &
                 pfull(i,j,:), zfull(i,j,:), phalf(i,j,:),  &
                 zhalf(i,j,:), pblht(i,j), tkemiz(i,j), qstar(i,j), &
                 cush(i,j), cbmf(i,j), land(i,j),  coldT(i,j), sd, Uw_p, ac, &
                 max_depletion_rate, Don_cape%xcape(i,j), &
                 Don_conv%dcape(i,j),   &
                 Don_conv%amax(i,j), dise_v(i,j,:), disa_v(i,j,:),    &
                 Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:), &
                 Don_cape%model_r(i,j,:), &

                 env_t, env_r, parcel_t, parcel_r, cape_p, & 

                 exit_flag(i,j), Don_conv%amos(i,j), Don_conv%a1(i,j),&
                 ermesg, error)

              else  ! (do_hires_cape)
!  lo res calc for cape in convection and in closure; standard 
!  donner_lite configuration
                 call don_d_determine_cloud_area_miz  &
!               (me, nlev_lsm, ntr, dt, nlev_hires, diag_unit,&
                (me, nlev_lsm, ntr, dt, nlev_lsm  , diag_unit,&
                 debug_ijt, Param, Initialized, Nml, xgcm_v(i,j,:,:), &
                 pfull(i,j,:), zfull(i,j,:), phalf(i,j,:),  &
                 zhalf(i,j,:), pblht(i,j), tkemiz(i,j), qstar(i,j), &
                 cush(i,j), cbmf(i,j), land(i,j),  coldT(i,j), sd, Uw_p, ac, &
                 max_depletion_rate, Don_cape%xcape(i,j), &
                 Don_conv%dcape(i,j),   &
                 Don_conv%amax(i,j), dise_v(i,j,:), disa_v(i,j,:),    &
                 Don_cape%model_p(i,j,:), Don_cape%model_t(i,j,:), &
                 Don_cape%model_r(i,j,:), Don_cape%env_t(i,j,:), &
                 Don_cape%env_r(i,j,:), Don_cape%parcel_t(i,j,:), &
                 Don_cape%parcel_r(i,j,:), Don_cape%cape_p(i,j,:), &
                 exit_flag(i,j), Don_conv%amos(i,j), Don_conv%a1(i,j),&
                 ermesg, error)
              endif
            endif
            endif

!             Don_budgets%liq_prcp(i,j,:) =    &
!                         Don_budgets%liq_prcp(i,j,:)*Don_conv%a1(i,j)
!             Don_budgets%frz_prcp(i,j,:) =    &
!                         Don_budgets%frz_prcp(i,j,:)*Don_conv%a1(i,j)
!           if (Initialized%do_conservation_checks) then
!             Don_budgets%vert_motion(i,j) =    &
!                          Don_budgets%vert_motion(i,j)*Don_conv%a1(i,j)
!             Don_budgets%lheat_precip(i,j) =   &
!                         Don_budgets%lheat_precip(i,j)*Don_conv%a1(i,j)
!           endif

!----------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!----------------------------------------------------------------------
            if (error /= 0 ) return
          endif

          if (.not.Nml%do_donner_lscloud) then 
             Don_conv%ecds(i,j,:) = Don_conv%ecds(i,j,:)*  &
                                        (1.0E03*Param%seconds_per_day)
             Don_conv%eces(i,j,:) = Don_conv%eces(i,j,:)*  &
                                        (1.0E03*Param%seconds_per_day)
            do k=1,nlev_lsm
              Don_conv%fre(i,j,nlev_lsm+1-k)=enev(k)
            end do
          endif

!--------------------------------------------------------------------
!    if the exit_flag was set within determine_cloud_area (due to
!    not having a net column non-zero moisture forcing (the "moisture 
!    constraint") set a flag to indicate that deep convection is turned
!    off; output fields will be made to reflect the absence of the 
!    effects of deep convection in this column. 
!--------------------------------------------------------------------
        end do
      end do

!--------------------------------------------------------------------



end subroutine don_d_mulsub_k



!######################################################################

subroutine don_d_integ_cu_ensemble_k             &
         (nlev_lsm, nlev_hires, ntr, me, diag_unit, debug_ijt, &
          lofactor, Param, Col_diag, Nml, Initialized, temp_c, &
          mixing_ratio_c, pfull_c, phalf_c,   &
          tracers_c, sfc_sh_flux_c, sfc_vapor_flux_c,   &
          sfc_tracer_flux_c, plzb_c, exit_flag_c, ensmbl_precip,    &
          ensmbl_cond, ensmbl_anvil_cond_liq,  &
          ensmbl_anvil_cond_liq_frz, &
          ensmbl_anvil_cond_ice, pb, pt_ens, ampta1, amax, &
          emsm, rlsm, cld_press, ensmbl_melt, ensmbl_melt_meso, &
          ensmbl_freeze, ensmbl_freeze_meso, ensmbl_wetc, &
          disb, disc_liq, disc_ice, dism_liq, dism_liq_frz, &
          dism_liq_remelt, dism_ice, dism_ice_melted, &
          disp_liq, disp_ice, disz, disz_remelt, disp_melted,        &
          disze1, disze2,disze3,                                      &
          disd, disv, disg_liq, disg_ice, enctf, encmf, enev,  &
          ecds_liq, ecds_ice, eces_liq, eces_ice, ensmbl_cloud_area,&
          cuq, cuql_v, detmfl, uceml, qtren, etsm, lmeso, &
          frz_frac, meso_frz_intg_sum,  ermesg, error, melting_in_cloud, &
          i, j, Don_cem)

!----------------------------------------------------------------------
!    subroutine integrate_cumulus_ensemble works on a single model 
!    column. all profile arrays used in this subroutine and below have 
!    index 1 nearest the surface. it first determines the lifting conden-
!    sation level (if one exists) of a parcel moving from the specified 
!    parcel_launch_level. if an lcl is found, subroutine 
!    donner_cloud_model_cloud_model is called to determine the behavior
!    of each of kpar cloud ensemble mem-
!    bers assumed present in the column (each ensemble member is ass-
!    umed to have a different entrainment rate). if all ensemble members
!    produce deep convection, the ensemble statistics are produced for 
!    use in the large-scale model; otherwise deep convection is not seen
!    in the large-scale model in this grid column. if the ensemble will 
!    support a mesoscale circulation, its impact on the large-scale model
!    fields is also determined. upon completion, the appropriate output 
!    fields needed by the large-scale model are returned to the calling 
!    routine.
!----------------------------------------------------------------------
use donner_types_mod, only : donner_param_type, &
                             donner_nml_type, donner_column_diag_type, &
                             donner_initialized_type, donner_cem_type

implicit none 

!----------------------------------------------------------------------
integer,                           intent(in)    :: nlev_lsm,    &
                                                    nlev_hires, ntr, &
                                                    me, diag_unit
logical,                           intent(in)    :: debug_ijt
type(donner_param_type),           intent(in)    :: Param
type(donner_column_diag_type),     intent(in)    :: Col_diag
type(donner_nml_type),             intent(in)    :: Nml   
type(donner_initialized_type),     intent(in)    :: Initialized
real,    dimension(nlev_lsm),      intent(in)    :: temp_c,   &
                                                    mixing_ratio_c,   &
                                                    pfull_c
real,    dimension(nlev_lsm+1),    intent(in)    :: phalf_c
real,    dimension(nlev_lsm,ntr),  intent(in)    :: tracers_c           
real,                              intent(in)    :: sfc_sh_flux_c,   &
                                                    sfc_vapor_flux_c 
real,    dimension(ntr),           intent(in)    :: sfc_tracer_flux_c 
real,                              intent(in)    :: plzb_c
real,                              intent(in)    :: lofactor
logical,                           intent(inout) :: exit_flag_c  
real,                              intent(out)   ::    &
                     ensmbl_precip, ensmbl_cond,                    &
                     ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
                     ensmbl_anvil_cond_ice, pb, pt_ens, ampta1, amax
real,    dimension(nlev_hires),    intent(out)   :: emsm, rlsm,  &
                                                    cld_press
real,    dimension(nlev_lsm),      intent(out)   :: ensmbl_melt,   &
                                                    ensmbl_melt_meso,&
                                                    ensmbl_freeze,&
                                                    ensmbl_freeze_meso,&
                                                    disb,       disd, &
                                                    disv, &
                                                    disc_liq, disc_ice,&
                                                    dism_liq, dism_ice,&
                                                    dism_ice_melted, &
                                                    dism_liq_frz, &
                                                    dism_liq_remelt, &
                                                    disp_liq, disp_ice,&
                                                    disp_melted, &
                                                    disz_remelt, &
                                                    disz, disze1,  &
                                                    disze2, disze3,&
                                                    enctf, encmf, &
                                                    disg_liq, disg_ice,&
                                                    enev,             &
                                                    ecds_liq, ecds_ice,&
                                                    eces_liq, eces_ice,&
                                                    ensmbl_cloud_area, &
                                                    cuq, cuql_v, &
                                                    detmfl, uceml
real,    dimension(nlev_lsm,ntr),  intent(out)   :: qtren, ensmbl_wetc
real,    dimension(nlev_hires,ntr),intent(out)   :: etsm
logical,                           intent(out)   :: lmeso       
real   ,                           intent(out)   :: frz_frac
logical,                           intent(out)   :: meso_frz_intg_sum 
character(len=*),                  intent(out)   :: ermesg
integer,                           intent(out)   :: error
logical ,                          intent(out)   :: melting_in_cloud
integer,                           intent(in)    :: i, j
type(donner_cem_type),             intent(inout) :: Don_cem

!---------------------------------------------------------------------
!   intent(in) variables:
! 
!     nlev_lsm       number of model layers in large-scale model
!     nlev_hires     number of model layers in hi-res cloud model
!                    of the donner deep convection parameterization
!     ntr            number of tracers to be transported by donner
!                    convection
!     me             local pe number
!     diag_unit      unit number for column diagnostics output, if 
!                    diagnostics are requested for the current column
!     debug_ijt      logical indicating whether current column requested
!                    column diagnostics
!     Param          donner_param_type variable containingthe parameters
!                    of the donner deep convection parameterization
!     Col_diag       donner_column_diagtype variable containing the
!                    information defining the columns fro which diagnos-
!                    tics are desired.
!     Nml            donner_nml_type variable containing the donner_nml
!                    variables that are needed outsied of donner_deep_mod
!     temp_c         temperature field at model full levels 
!                    index 1 nearest the surface [ deg K ]
!     mixing_ratio_c        vapor mixing ratio at model full levels 
!                    index 1 nearest the surface
!                    [ kg(h2o) / kg(dry air) ]
!     pfull_c         pressure field at large-scale model full levels 
!                    index 1 nearest the surface [ Pa ]
!     phalf_c        pressure field at large-scale model half-levels 
!                    index 1 nearest the surface [ Pa ]
!     tracers_c      tracer fields that are to be transported by donner
!                    convection.  index 1 nearest the surface 
!                    [ kg (tracer) / kg (dry air) ]
!     sfc_sh_flux_c  sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux_c water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     sfc_tracer_flux_c  
!                    flux across the surface of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     plzb_c         level of zero buoyancy for a parcel lifted from
!                    the parcel_launch_level.  [ Pa ]
!
!     cumulus ensemble member fields (see also donner_types.h):
!
!     --- single level ---
!
!     Don_cem_cell_precip 
!                    area weighted convective precipitation rate
!                    [ mm/day ]
!     Don_cem_pb     pressure at cloud base for ensemble (currently,
!                    all ensemble members have same base) [ Pa ]
!     Don_cem_ptma   pressure at cloud top for ensemble [ Pa ]
!
!     --- lo-res multi-level ---
! 
!     Don_cem_h1     condensation rate profile on lo-res grid
!                    for the current ensemble member
!                    [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!
!     --- hi-res multi-level ---
!
!     Don_cem_qlw    profile of cloud water for the current ensemble
!                    member [ kg(h2o) / kg(air) ]
!     Don_cem_cfracice
!                    fraction of condensate that is ice [ fraction ]
!     Don_cem_wv     vertical velocity profile [ m / s ]
!     Don_cem_rcl    cloud radius profile [ m ]
!
!   intent(inout) variables:
!
!     exit_flag_c    logical indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    current model column 
!
!   intent(out) variables:
!    
!     ensmbl_precip      sum of precipitation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_cond        sum of condensation rate over ensemble members,
!                        # 1 to the current, weighted by the area at 
!                        cloud base of each member
!                        [ mm / day ]
!     ensmbl_anvil_cond  sum of rate of transfer of condensate from cell 
!                        to anvil over ensemble members, # 1 to the c
!                        current, weighted by the area at cloud base of 
!                        each member [ mm / day ]
!     ensmbl_wetc        sum of wet-deposition rates from ensemble
!                        member #1 to current, weighted by the ratio
!                        of the member area to the area of member #1
!                        at cloud base
!                        [kg(tracer)/kg/sec]
!                        vertical index 1 at cloud base
!     pb                 pressure at cloud base for ensemble (all ensem-
!                        ble members have same base) [ Pa ]
!     pt_ens             pressure at cloud top for the ensemble (top 
!                        pressure of deepest ensemble member) [ Pa ]
!     ampta1             cloudtop anvil area (assumed to be five times
!                        larger than the sum of the cloud top areas of 
!                        the ensemble members, as in Leary and Houze 
!                        (1980).  [ fraction ]
!     amax               maximum allowable area of cloud base that is
!                        allowed; if cloud base area is larger than 
!                        amax, the cloud fractional area somewhere in
!                        the grid box would be greater than one, which 
!                        is non-physical.
!     emsm               vertical profile on the hi-res grid of vertical
!                        moisture flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's cont-
!                        ribution being weighted by its cloud area at 
!                        level k relative to the cloud base area of 
!                        ensemble member #1  
!                        [ kg (h2o) / ( kg(dry air) sec ) ]
!     rlsm               vertical profile on the hi-res grid of conden-
!                        sation rate, summed over ensemble members # 1 to
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!     cld_press          pressures at hi-res model levels [ Pa ]
!     ensmbl_melt        vertical profile on the lo-res grid of ice melt,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     ensmbl_freeze      vertical profile on the lo-res grid of freezing,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!     disg               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation
!                        associated with the evaporation of condensate
!                        in the convective downdraft and updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     enev               vertical profile on the lo-res grid of the      
!                        cloud-area-weighted profile of the potential
!                        cloud water evaporation, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at !
!                        level k relative to the cloud base area of 
!                        ensemble member #1.  this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ ( kg(h2o) ) / ( kg(dry air) sec ) ] 
!     enctf              vertical profile on the lo-res grid of the entr-
!                        opy forcing, consisting of the sum of the
!                        vertical entropy flux convergence and the latent
!                        heat release, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area
!                        of ensemble member #1
!                        [ deg K / day ]                        
!     encmf              vertical profile on the lo-res grid of the      
!                        moisture forcing, consisting of the sum of the
!                        vertical moisture flux convergence and the cond-
!                        ensation, summed over ensemble members # 1 to 
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ ( kg(h2o) ) / ( kg( dry air) day ) ] 
!     disb               vertical profile on the lo-res grid of the      
!                        temperature flux convergence, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area 
!                        of ensemble member #1.  
!                        [ deg K / day ] 
!     disc               vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation, 
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!     disd               vertical profile on the lo-res grid of the      
!                        vertical moisture flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        lo-res grid for the current ensemble member 
!                        [  g(h2o) / ( kg(dry air) day ) ]
!     ecds               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective downdraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     eces               vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!     ensmbl_cloud_area  total cloud area profile over all ensemble
!                        members on large_scale model grid [ fraction ]
!     cuq                ice water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     cuql_v             liquid water profile on large-scale model grid, 
!                        normalized by ensemble cloud area.
!     uceml              upward mass flux on large_scale model grid     
!                        [ kg (air) / (sec m**2) ]
!     detmfl             detrained mass flux on large-scale model grid
!                        normalized by ensemble cloud area
!                        [ kg (air) / (sec m**2) ]
!     etsm               vertical profile on the hi-res grid of vertical
!                        tracer flux convergence, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at i
!                        level k relative to the cloud base area of 
!                        ensemble member #1 
!                        [ kg (tracer) / ( kg(dry air) sec ) ]
!     qtren              vertical profile on the lo-res grid of the      
!                        vertical tracer flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!     lmeso              logical variable; if .false., then it has been
!                        determined that a mesoscale circulation cannot
!                        exist in the current column. final value not
!                        determined until all ensemble members have been
!                        integrated. 
!     ermesg             character string containing any error message
!                        that is returned from a kernel subroutine
!
!---------------------------------------------------------------------

!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (K/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     pb               pressure at base of cumulus updrafts (Pa)
!     pmd              pressure at top of mesoscale downdraft (Pa)
!     pztm             pressure at top of mesoscale updraft (Pa)
!     mrmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     qtmes(nlev,ncont)  tracer tendency due to mesoscale tracer-flux
!                        convergence (kg/kg/s) (normalized by a(1,p_b))
!                        index 1 at ground 
!     qtren_v          normalized tracer tendency due to cells...
!                      (lon,lat,vert,tracer index)
!                      Vertical index increases as height increases.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre_v           total normalized precipitation (mm/day)
!     detmfl(nlev)     detrained mass flux from cell updrafts
!                      (normalized by a(1,p_b))
!                      (index 1 near atmosphere bottom)
!                      (kg/((m**2)*s)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wtp_v            tracer redistributed by mesoscale processes
!                      (kg/kg/s) (normalized by a(1,p_b))
!                      vertical index increases with increasing height
!                      (lon,lat,vert,tracer index)
!--------------------------------------------------------------------


!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    detmfh [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    dpftr    [ (m**2 kg(tracer)) / (kg(air) sec) ] ==>   
!!                                          [ kg(tracer)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]




!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        current_displ  integrated low-level displacement (Pa)
!        cape_p   pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb_c   pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!

!----------------------------------------------------------------------
!   local variables:

      real,    dimension (nlev_hires)     ::                &
              efchr, emfhr, rcl, dpf, qlw, dfr, cfracice, &
              alp, cld_evap, flux, ucemh, cuql, cuqli, detmfh, tcc, wv

      real,    dimension (nlev_lsm)       ::           &
                  q1, cell_freeze, cell_melt,   &
              h1_liq, h1_ice, meso_melt, meso_freeze, h1_2, &
              evap_rate, ecd, ecd_liq, ecd_ice, &
              ece, ece_liq, ece_ice, sfcq, sfch
      real, dimension (nlev_lsm,ntr)  :: wetdepl

      real,    dimension (nlev_hires,ntr) :: etfhr, dpftr
      real,    dimension (nlev_lsm,ntr)   :: qtr
      real,    dimension (Param%kpar)     :: cuto, preto, ptma
      integer, dimension (Param%kpar)     :: ncca

      logical ::   lcl_reached                  
      integer ::   ncc_kou, ncc_ens
      integer ::   k,    kou
      integer ::   kk
      logical ::   meso_frz_intg                 
      real    ::   al, dp, mrb,  &
                   summel, ptt,   &
                   sbl, psmx, dint, cu, cell_precip,&
                   ca_liq, ca_ice, apt, &
                   tb, alpp,   &
                   pcsave, ensmbl_cld_top_area  
      real     ::  meso_frac, precip_frac, frz_frac_non_precip,  &
                   bak, meso_frz_frac, pmelt_lsm, precip_melt_frac, &
                   ecei_liq,  ci_liq_cond, ci_ice_cond
     

!----------------------------------------------------------------------
!   local variables:
!
!      ensmbl_cld_top_area  
!                       sum of the cloud top areas over ensemble members 
!                       # 1 to the current, normalized by the cloud base
!                       area of ensemble member # 1 [ dimensionless ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the character string which will contain any error mes-
!    sages returned through this subroutine.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    if in diagnostics column, output the large-scale model temperature,
!    vapor mixing ratio and full-level pressure profiles (index 1 near-
!    est the surface).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm-Col_diag%kstart+1
          write (diag_unit, '(a, i4, f20.14, e20.12, f19.10)')&
                'in mulsub: k,T,Q,P= ',k, temp_c(k),  &
                                      mixing_ratio_c(k), pfull_c(k)
        end do
      endif

!--------------------------------------------------------------------
!    call don_cm_lcl_k to calculate the temperature (tb), a
!    pressure (pb) and mixing ratio (mrb) at the lifting condensation 
!    level for a parcel starting from the parcel_launch_level. if a sat-
!    isfactory lcl is not reached for this parcel, the logical variable 
!    lcl_reached will be set to .false..
!--------------------------------------------------------------------
      call don_cm_lcl_k    &
           (Param, temp_c (Nml%parcel_launch_level),    &
            pfull_c       (Nml%parcel_launch_level),    &
            mixing_ratio_c(Nml%parcel_launch_level),   &
            tb, pb, mrb, lcl_reached, ermesg, error)     

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    if in diagnostics column and an lcl was defined, output the lcl 
!    temperature, pressure and mixing ratio. if an acceptble lcl was 
!    not reached, print a message.
!--------------------------------------------------------------------
      if (debug_ijt) then
        if (lcl_reached) then
          write (diag_unit, '(a, f20.14, f19.10, e20.12)') &
                                'in mulsub: tb,pb,qb= ',tb, pb, mrb  
        else
          write (diag_unit, '(a)') 'in mulsub: lcl not reached'
        endif
      endif

!--------------------------------------------------------------------
!    if an acceptable lcl was not reached, set exit_flag_c so that the
!    remaining computations for this column are bypassed, and return to
!    calling routine. 
!--------------------------------------------------------------------
      if (.not. lcl_reached) then
        exit_flag_c = .true.
        return
      endif
 
!---------------------------------------------------------------------
!    if calculations are continuing, initialize needed variables.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize variables which will accumulate scalar sums over all 
!    ensemble members.
!---------------------------------------------------------------------
      ensmbl_precip       = 0.
      ensmbl_cond         = 0.
      ensmbl_anvil_cond_liq   = 0.
      ensmbl_anvil_cond_liq_frz   = 0.
      ensmbl_anvil_cond_ice   = 0.
      ensmbl_cld_top_area = 0.

!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the cloud-model grid.
!---------------------------------------------------------------------
      do k=1,nlev_hires
        cuql(k)   = 0.
        cuqli(k)  = 0.
        ucemh(k)  = 0.
        detmfh(k) = 0.
        alp(k)    = 0.
        rlsm(k)   = 0.
        emsm(k)   = 0.
        etsm(k,:) = 0.
      end do

!---------------------------------------------------------------------
!    initialize the variables which will contain the sum over the 
!    ensemble members of the vertical profiles of various quantities 
!    on the large-scale model grid.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        ensmbl_freeze(k)    = 0.
        ensmbl_freeze_meso(k)    = 0.
        ensmbl_melt(k)    = 0.
        ensmbl_melt_meso(k)    = 0.
        disb(k)    = 0.
        disc_liq(k) = 0.
        disc_ice(k) = 0.
        dism_liq(k) = 0.
        dism_liq_frz(k) = 0.
        dism_liq_remelt(k) = 0.
        dism_ice(k) = 0.
        dism_ice_melted(k) = 0.
        disp_liq(k) = 0.
        disp_ice(k) = 0.
        disp_melted(k) = 0.
        disd(k)    = 0.
        disv(k)    = 0.
        disz(k) = 0.
        disz_remelt(k) = 0.
        disze1(k) = 0.
        disze2(k) = 0.
        disze3(k) = 0.
        ecds_liq(k)    = 0.
        ecds_ice(k)    = 0.
        eces_liq(k)    = 0.
        eces_ice(k)    = 0.
        enctf(k)   = 0.
        encmf(k)   = 0.
        disg_liq(k)    = 0.
        disg_ice(k)    = 0.
        enev(k)    = 0.
        qtren(k,:) = 0.
        ensmbl_wetc(k,:)  = 0.
      end do

      evap_rate = 0.

!--------------------------------------------------------------------
!    initialize a logical variable which will indicate whether a
!    mesoscale circulation is present in this column. this may be 
!    precluded via the nml variable allow_mesoscale_circulation. if any
!    ensemble members are unable to support a mesoscale circulation, 
!    lmeso will be set to .false. within the following loop over the kpar
!    ensemble members. if the first member of the ensemble (the most 
!    entraining) can, then it is likely (but not guaranteed) that the 
!    ensemble will be able to.
!--------------------------------------------------------------------
      if (Nml%allow_mesoscale_circulation) then
        lmeso = .true.
      else
        lmeso = .false.
      endif

!--------------------------------------------------------------------
!    define the array of cloud model pressure levels (cld_press).
!--------------------------------------------------------------------
      do k=1,nlev_hires
        cld_press(k) = pb + (k-1)*Param%dp_of_cloud_model
      end do

!--------------------------------------------------------------------
!    if this is the first ensemble member, initialize the variables
!    which are defined on this call and will be used by the other
!    ensemble members.
!--------------------------------------------------------------------
      pcsave = phalf_c(1)

!--------------------------------------------------------------------
!    loop over the KPAR members of the cumulus ensemble.
!--------------------------------------------------------------------
      meso_frz_intg_sum = .false.
      do kou=1,Param%kpar

!-------------------------------------------------------------------
!    define the appropriate entrainment factor (alpp) for this ensemble
!    member using values based on observations either obtained from
!    the GATE or KEP studies.
!-------------------------------------------------------------------
        if (trim(Nml%entrainment_constant_source) == 'gate') then
          alpp = Param%max_entrainment_constant_gate/  &
                           Param%ensemble_entrain_factors_gate(kou)
        else if (trim(Nml%entrainment_constant_source) == 'kep') then
          alpp = Param%max_entrainment_constant_kep/  &
                           Param%ensemble_entrain_factors_kep(kou)
        else
          ermesg = 'invalid entrainment_constant_source'
          error = 1
          return
        endif

        if (Nml%do_lands) then
          alpp = alpp*lofactor     
        endif

        if (debug_ijt) then
          write (diag_unit, '(a)')    &
                     'in mulsub: phalf, temp= :'
          do k=1,nlev_lsm 
          write (diag_unit, '(i4, 2f19.10)')    &
                      k, phalf_c(k), temp_c(k)
          end do
        endif

       pmelt_lsm = 2.0e05
       if( temp_c(1) >  Param%KELVIN ) then
       do k=1,nlev_lsm-1
        if ((temp_c(k) >= Param%KELVIN) .and.    &
           (temp_c(k+1) <= Param%KELVIN)) then
          pmelt_lsm = phalf_c(k+1)
          exit
        endif
       end do
       endif

       if (debug_ijt) then
         write (diag_unit, '(a, 2f19.10)')    &
         'before cm_cloud_model call pb,  pmelt_lsm    = ', &
                                    pb, pmelt_lsm
       endif
!--------------------------------------------------------------------
!    call cloud_model to obtain the in-cloud and environmental profiles
!    and fluxes and column integrals associated with this ensemble 
!    member.
!--------------------------------------------------------------------
        call don_cm_cloud_model_k   &
             (nlev_lsm, nlev_hires, ntr, kou, diag_unit, debug_ijt,   &
              Param, Col_diag, Initialized, tb, pb, alpp, cld_press, &
              temp_c, mixing_ratio_c, pfull_c, phalf_c, tracers_c, &
              pcsave,  exit_flag_c, wv, rcl, dpf, dpftr, qlw, dfr, flux, &
              ptma(kou), dint, cu, cell_precip, apt, cell_melt, &
              pmelt_lsm, summel, efchr, emfhr, cfracice, etfhr, &
              ncc_kou, tcc, ermesg, error)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return


!--------------------------------------------------------------------
!    if the cloud thickness is less than pdeep_mc, it will not
!    support a mesoscale circulation. set a logical flag to indicate
!    the absence of a mesoscale component for this column's cloud
!    ensemble.
!--------------------------------------------------------------------
        if (lmeso) then
          if ((pb - ptma(kou)) < Param%pdeep_mc)  then
            lmeso = .false.
          endif
        endif


        if (exit_flag_c) then
          if (lmeso) then
            cell_melt(:) = 0.0
          endif
           return
        endif

!--------------------------------------------------------------------
!    if calculations are continuing, 
!--------------------------------------------------------------------

!----------------------------------------------------------------------
!    if in diagnostics column, output the cloud base (pb) and cloud top
!    (ptma) pressures, and the mesoscale circulation logical variable
!    (lmeso).
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 2f19.10,1l4)')    &
         'in mulsub: PB,PT, lmeso= ', pb, ptma(kou), lmeso
        endif

!---------------------------------------------------------------------
!    define the cloud water from this ensemble member which must be 
!    evaporated if it turns out that there is no mesoscale circulation 
!    associated with the ensemble.
!---------------------------------------------------------------------
        cld_evap(:) = -dpf(:)*(1. - (cell_precip/cu))

!---------------------------------------------------------------------
!    define the pressure one cloud model level above cloud top (ptt).
!---------------------------------------------------------------------
        ptt = ptma(kou) + Param%dp_of_cloud_model

!----------------------------------------------------------------------
!    call define_lo_res_model_profiles to map profiles generated on the
!    cloud model grid to the vertical grid of the large-scale model for
!    this ensemble member.
!----------------------------------------------------------------------
        call don_d_def_lores_model_profs_k        &
             (nlev_lsm, nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt,&
              Nml, Param, pb, ptt, sfc_vapor_flux_c, sfc_sh_flux_c,  &
              sfc_tracer_flux_c, pfull_c, phalf_c, cld_press, tcc, dpf,&
              dpftr, dfr, cld_evap, qlw, emfhr, efchr, etfhr, &
              cell_freeze, evap_rate,     h1_liq, h1_ice, ci_liq_cond, &
              ci_ice_cond, h1_2, q1, qtr, wetdepl, ermesg, error)

        if (cu /= 0.0) then
          precip_frac = cell_precip/cu
        else 
          precip_frac = 0.
        endif

      if (Nml%do_ensemble_diagnostics) then
!----------------------------------------------------------------------
!    save "Don_cem" diagnostics for this ensemble member.
!----------------------------------------------------------------------
        Don_cem%cell_precip(i,j,kou) = cell_precip
        Don_cem%pb(i,j,kou) = pb
        Don_cem%ptma(i,j,kou) = ptma(kou)
! reverse index order
        do k=1,nlev_lsm
          Don_cem%h1(i,j,k,kou) = h1_liq(nlev_lsm-k + 1)  + &
                              h1_ice(nlev_lsm-k + 1)  
        end do
        do k=1,nlev_hires
          Don_cem%qlw(i,j,k,kou) = qlw(k)
          Don_cem%cfracice(i,j,k,kou) = cfracice(k)
          Don_cem%wv(i,j,k,kou) = wv(k)
          Don_cem%rcl(i,j,k,kou) = rcl(k)
        end do
     endif

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

!---------------------------------------------------------------------
!    if this member of the ensemble supports a mesoscale circulation,
!    call mesub to obtain various terms related to moving condensate
!    from the convective tower into the mesoscale anvil for this member.
!---------------------------------------------------------------------
        if (lmeso) then
          call don_cm_mesub_k     &
               (Nml, pfull_c, nlev_lsm, me, diag_unit, debug_ijt, Param, cu,   &
                ci_liq_cond, ci_ice_cond, pmelt_lsm, cell_precip, &
                dint, plzb_c, pb, ptma(kou), temp_c, phalf_c,     &
                ca_liq, ca_ice,  ecd, ecd_liq, ecd_ice, ecei_liq, &
                ece, ece_liq, ece_ice, meso_freeze, meso_melt, ermesg, error)
        else
          ca_liq = 0.
          ca_ice = 0.
          meso_freeze = 0.
          meso_melt   = 0.
        endif

       if (pmelt_lsm < pb) then
         melting_in_cloud = .true.
       else
         melting_in_cloud = .false.
       endif

         if (ci_ice_cond /= 0.0) then
           if (melting_in_cloud) then
             precip_melt_frac = summel/ci_ice_cond
           else
             precip_melt_frac = 0.
           endif  
         else
           precip_melt_frac = 0.
         endif

       if (debug_ijt) then
         write (diag_unit, '(a, 3e20.12)')  &
            'in mulsub: h1_ice intg, summel, precip_melt_frac', &
                       ci_ice_cond, summel, precip_melt_frac
       endif

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
       if (error /= 0 ) return

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_hires_k to add this member's 
!    contribution to those fields on the cloud model grid that are being
!    summed over all ensemble members.
!---------------------------------------------------------------------
       call don_d_add_to_ensmbl_sum_hires_k    &
             (nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, &
              Param%arat(kou), cfracice, rcl, flux, emfhr, dpf, &
              qlw, etfhr, cuql, cuqli, ucemh, alp, rlsm, emsm, detmfh, &
              etsm, ermesg, error)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

!----------------------------------------------------------------------
!    define the fraction of total condensate which is transferred to
!    the mesoscale circulation.
!----------------------------------------------------------------------
        if (cu /= 0.0) then
          meso_frac = (ca_liq + ca_ice)/cu
        else
          meso_frac = 0.
        endif

!----------------------------------------------------------------------
!    if there is a mesoscale circulation, define the fraction of total 
!    liquid condensate transferred to the mesoscale circulation that was
!    frozen (meso_frz_frac). define a logical indicating whether any 
!    such condensate exists (meso_frz_intg). 
!----------------------------------------------------------------------
        if (lmeso) then
          bak = 0.
          do kk=1,nlev_lsm
            dp = phalf_c(kk) - phalf_c(kk+1)
            bak = bak + meso_freeze(kk)*dp
          end do
          bak = bak/(Param%grav)
          bak = bak/(Param%seconds_per_day*1.0e3)
          if (debug_ijt) then
            write (diag_unit, '(a, 3e20.12)')  &
                  'in mulsub: column meso_freeze', bak
          endif
          if (bak > 0.0) then
             meso_frz_intg = .true.
          else
             meso_frz_intg = .false.
          endif
          if (ci_liq_cond /= 0.0) then
            meso_frz_frac = bak/ci_liq_cond
          else
            meso_frz_frac = 0.
          endif
        else
          meso_frz_intg = .false.
          meso_frz_frac = 0.        
        endif

!---------------------------------------------------------------------
!    if there has been liquid condensate, define the fraction of liquid
!    condensate which froze (frz_frac). define the fraction of
!    liquid condensate which froze but did not precipitate out and so
!    is available for evaporation and transfer to the mesoscale 
!    circulation (frz_frac_non_precip). 
!---------------------------------------------------------------------
        if (ci_liq_cond /= 0.0) then
          frz_frac = dint/ci_liq_cond    
          frz_frac_non_precip = frz_frac*(1.-precip_frac)   

!---------------------------------------------------------------------
!    deal with the case when the liquid condensate defined to be frozen
!    is more than the liquid condensate remaining after the appropriate
!    cell precipitation. in this case, limit the amount frozen to that 
!    which is still present in the atmosphere, and modify the 
!    cell_freeze profile and dint integral, and the frz_frac and 
!    frz_frac_non_precip ratios.
!--------------------------------------------------------------------- 
          if (.not. melting_in_cloud) then
            if (meso_frz_frac == 0. .and.  meso_frac > 0.) then
              if (meso_frac < frz_frac_non_precip) then
                do k=1,nlev_lsm
                  cell_freeze(k) = cell_freeze(k)*meso_frac/  &
                                    frz_frac_non_precip
                end do  
                dint = dint *meso_frac/frz_frac_non_precip
                frz_frac = meso_frac/(1.-precip_frac)
                frz_frac_non_precip = meso_frac
              endif
            endif
          endif
        else

!---------------------------------------------------------------------
!    if there is no liquid condensate in the column, then there is no
!    frozen liquid condensate in the column.
!---------------------------------------------------------------------
          frz_frac_non_precip = 0.
          frz_frac = 0.
        endif

        if (debug_ijt) then
          write (diag_unit, '(a, 3e20.12)')  &
                   'in mulsub pre anvil_cond_frz: h1_liq intg, dint,&
                       & frz_frac_non_precip           ', &
                   ci_liq_cond, dint, frz_frac_non_precip
          write (diag_unit, '(a, 1e20.12)')  &
                                 'in mulsub : frz_frac', &
                                    frz_frac

!----------------------------------------------------------------------
!    if there is a mesoscale circulation, define the fraction of total 
!    liquid condensate transferred to the mesoscale circulation which 
!    is frozen (meso_frz_frac). define a logical indicating whether any 
!    such condensate exists (meso_frz_intg). 
!----------------------------------------------------------------------
              write (diag_unit, '(a, i4, 2e20.12)')  &
           'in mulsub : kou,  meso_frz_frac, precip_melt_frac', &
                                 kou,  meso_frz_frac, precip_melt_frac
          endif

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_intgl_k to add this member's 
!    contribution to those integrals that are being summed over all 
!    ensemble members.
!---------------------------------------------------------------------
        call don_d_add_to_ensmbl_sum_intgl_k    &
             (diag_unit, debug_ijt, lmeso,                   &
                      Param%arat(kou),      &
              ca_liq, ca_ice, frz_frac_non_precip, meso_frac,  &
              cell_precip, cu, apt, ensmbl_precip, ensmbl_cond,   &
                                 ensmbl_anvil_cond_liq, &
              ensmbl_anvil_cond_liq_frz, meso_frz_intg, meso_frz_frac,&
              ensmbl_anvil_cond_ice, ensmbl_cld_top_area, ermesg, error)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3f19.10)')    &
                     'in mulsub: meso_frac, precip_frac,frz_frac_non_precip:', &  
                     kou, meso_frac, precip_frac, frz_frac_non_precip
          write (diag_unit, '(a, i4, 4f19.10)')    &
                     'in mulsub: cu, ca, cell_precip, dint   :', &  
                         kou, cu, ca_liq + ca_ice, cell_precip,   &
                         dint*Param%seconds_per_day
          write (diag_unit, '(a, 3f19.10)')    &
                     'in mulsub: pmelt_lsm, pb, summel   :', &  
                         pmelt_lsm, pb,  summel          
        endif

!---------------------------------------------------------------------
!    call don_d_add_to_ensmbl_sum_lores_k to add this member's 
!    contribution to those fields on the lrge-scale model grid that are 
!    being summed over all ensemble members.
!---------------------------------------------------------------------

        call don_d_add_to_ensmbl_sum_lores_k    &
             (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, &
                             frz_frac, Param, Nml,  &
              Param%arat(kou), dint, cell_freeze,         cell_melt, &
              wetdepl, temp_c,   &
              h1_2, ecd, ecd_liq, ecd_ice, ece, ece_liq, ece_ice, &
              evap_rate, q1,     h1_liq, h1_ice, pfull_c, meso_melt, &
              meso_freeze, phalf_c, qtr, ensmbl_melt, ensmbl_melt_meso,&
              ensmbl_freeze, ensmbl_freeze_meso, ensmbl_wetc, &
              meso_frac, precip_frac, frz_frac_non_precip, &
              disz, disz_remelt, disp_melted, disze1, disze2, disze3,  &
              disp_liq, disp_ice, enctf, encmf, enev, disg_liq,  &
              disg_ice, disb, disc_liq, disc_ice, dism_liq,  &
              dism_liq_frz, dism_liq_remelt, dism_ice, dism_ice_melted,&
              ecds_liq, ecds_ice, eces_liq, eces_ice, disd, disv, &
              qtren, ermesg, error, meso_frz_intg, melting_in_cloud, &
              precip_melt_frac, meso_frz_frac)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

!--------------------------------------------------------------------
!    save the cloud top (ptma) pressures, the total condensation (cuto),
!    total precpitation (preto) and cloud top index (ncca) from this !
!    ensemble member.
!--------------------------------------------------------------------
         if (meso_frz_intg) meso_frz_intg_sum = .true.
        cuto(kou)  = cu
        preto(kou) = cell_precip
        ncca(kou)  = ncc_kou
      end do   ! (kou loop over ensemble members)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! 31   CONTINUE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!--------------------------------------------------------------------
!    if calculations are continuing: 
!--------------------------------------------------------------------

!----------------------------------------------------------------------
!    define ensemble cloud top pressure (pt_ens) to be the cloud top of 
!    the most penetrative ensemble member. this is frequently, but not 
!    always, the ensemble member with the lowest entrainment rate. 
!    cloud base pressure (pb) is the same for all ensemble members. 
!    define the cloud top index(ncc_ens)  as the highest of any ensemble 
!    member.
!----------------------------------------------------------------------
      pt_ens  = MINVAL (ptma)
      ncc_ens = MAXVAL (ncca)

!----------------------------------------------------------------------
!    divide the ensemble mean ice and liquid condensate terms by the 
!    total cloud area to define the average cloud water and cloud ice 
!    concentrations within the cloudy area, as opposed to averaged over 
!    the entire grid box.
!----------------------------------------------------------------------
      do k=1,ncc_ens
        if (alp(k) > 0.) then
          cuql(k)  = cuql(k)/alp(k)
          cuqli(k) = cuqli(k)/alp(k)
        endif
      end do

!---------------------------------------------------------------------
!    define the cloudtop anvil area (ampta1), assumed to be five times 
!    larger than the sum of the cloud top areas of the ensemble members,
!    as in Leary and Houze (1980), 
!---------------------------------------------------------------------
      ampta1 = 5.*ensmbl_cld_top_area

!---------------------------------------------------------------------
!    if there is no precipitation production in this column, set the 
!    inverse of the max cloud area at any layer in the column to be 0.0.
!---------------------------------------------------------------------
      if (ensmbl_precip == 0.0) then
        amax      = 0.0
      else

!---------------------------------------------------------------------
!    if there is precip in the column, determine the maximum convective 
!    cell area at any level in the column (al). the total normalized 
!    cloud area in the column (cell area + mesoscale area) cannot be 
!    greater than 1.0. this constraint imposes a limit on the cloud area
!    at cloud base (amax). this limit will be imposed in subroutine
!    determine_cloud_area. see "a bounds notes" (7/6/97).
!---------------------------------------------------------------------
        al = MAXVAL (alp)
        amax = 1./(al + ampta1)
      endif

!---------------------------------------------------------------------
!    if in diagnostics column, output the total ensemble condensation,
!    (ensmbl_cond), precipitation (ensmbl_precip), and condensate 
!    transferred into the anvil (ensmbl_anvil_cond). also output 
!    surface pressure (phalf_c(1)), ensemble cloud base nd cloud top 
!    pressures (pb, pt_ens), the flag indicating if a mesoscale circul-
!    ation is present in the grid column (lmeso), and the cloud top anvil
!    area (ampta1).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12, a, e20.12)')  &
                      'in mulsub: CUTOT=', ensmbl_cond, ' PRETOT=', &
                                      ensmbl_precip
        write (diag_unit, '(a, 4e20.12)') &
               'in mulsub: CATOT, (sum, liq, frzliq,ice)=', &
              ensmbl_anvil_cond_liq  + ensmbl_anvil_cond_liq_frz  +  &
                                   ensmbl_anvil_cond_ice, &
                                     ensmbl_anvil_cond_liq, &
                                     ensmbl_anvil_cond_liq_frz, &
                                     ensmbl_anvil_cond_ice
        write (diag_unit, '(a, 3f19.10, 1l4)')  &
              'in mulsub: ps,pb,pt,lmeso= ',   &
                     phalf_c(1), pb, pt_ens, lmeso
        write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: ampt= ',ampta1     
      endif

!----------------------------------------------------------------------
!    define the pressure one level above cloud top (ptt).
!----------------------------------------------------------------------
      ptt = pt_ens + Param%dp_of_cloud_model

!--------------------------------------------------------------------
!    call define_ensemble_profiles to produce vertical profiles 
!    representing the ensemble-total cloud area (ensmbl_cloud_area), 
!    cloud liquid (cuql_v), cloud ice (cuq), mass flux(uceml) and
!    detrained mass flux (detmfl).
!--------------------------------------------------------------------
      call don_d_def_ensemble_profs_k    &
           (nlev_lsm, nlev_hires, ncc_ens, diag_unit, debug_ijt, ptt, &
            cld_press, alp, detmfh, ucemh, cuql, cuqli, phalf_c,  &
            ensmbl_cloud_area, cuql_v, cuq, detmfl, uceml, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, call error_mesg.
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    execute the following code if in a diagnostics column.
!---------------------------------------------------------------------
      if (debug_ijt) then

!----------------------------------------------------------------------
!    define the pressure at the large-scale model interface level at or 
!    just above cloud base (psmx).
!----------------------------------------------------------------------
        do k=1,nlev_lsm
          if ((phalf_c(k+1) <= pb) .and. (phalf_c(k) >= pb)) then
            psmx = phalf_c(k+1)
            exit
          endif
        end do

!----------------------------------------------------------------------
!    define the integrated boundary layer heating rate (sbl) due to the 
!    surface heat flux (sfcsf_v). it is defined in units of (deg K)/sec.
!    call don_u_map_hires_i_to_lores_c_k to distribute
!    this heating over the boundary layer.
!---------------------------------------------------------------------
       sbl = Param%grav*sfc_sh_flux_c/((phalf_c(1) - psmx)*Param%cp_air)
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
             'in cm_intgl_to_gcm_col: xav,p1,p2= ',sbl, phalf_c(1), psmx 
        call don_u_map_hires_i_to_lores_c_k   &
             (nlev_lsm, sbl, phalf_c(1), psmx, phalf_c, sfch, ermesg, error)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

        do k=1,size(sfch(:))
          if (sfch(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                            'in cm_intgl_to_gcm_col: k,x= ',k,sfch(k)
          endif
        end do

!----------------------------------------------------------------------
!    define the integrated boundary layer moistening rate (sbl) due to 
!    the surface moisture flux (sfcqf_v), which is defined in units of 
!    kg(h2o) per m**2 per sec. call 
!    don_u_map_hires_i_to_lores_c_k to distribute 
!    this moistening over the boundary layer.
!---------------------------------------------------------------------
        sbl = (sfc_vapor_flux_c*Param%grav)/(phalf_c(1) - psmx)
        write (diag_unit, '(a, e20.12, 2f19.10)')  &
             'in cm_intgl_to_gcm_col: xav,p1,p2= ',sbl, phalf_c(1), psmx 
        call don_u_map_hires_i_to_lores_c_k   &
             (nlev_lsm, sbl, phalf_c(1), psmx, phalf_c, sfcq, ermesg, error)

!---------------------------------------------------------------------
!    if an error message was returned from the kernel routine, return
!    to the calling program where it will be processed.
!---------------------------------------------------------------------
        if (error /= 0 ) return

        do k=1,size(sfcq(:))
          if (sfcq(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)') &
                            'in cm_intgl_to_gcm_col: k,x= ',k,sfcq(k)
          endif
        end do
      endif ! (debug_ijt)

!---------------------------------------------------------------------



end subroutine don_d_integ_cu_ensemble_k 

!#######################################################################

subroutine don_d_column_end_of_step_k  &
         (isize, jsize, nlev_lsm, ntr, Col_diag, exit_flag,   &
          total_precip, parcel_rise, temperature_forcing,&
          moisture_forcing, tracers, Don_cape, Don_conv, ermesg, error)       

!----------------------------------------------------------------------
!    subroutine don_d_column_end_of_step outputs the final values of
!    significant fields generated by donner_deep_mod in any columns
!    for which column diagnostics were requested, and in which deep
!    convection is present.
!----------------------------------------------------------------------

use donner_types_mod, only : donner_cape_type, donner_conv_type, &
                             donner_column_diag_type

implicit none

!----------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize,  &
                                                   nlev_lsm, ntr
type(donner_column_diag_type),    intent(in)    :: Col_diag
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
real,    dimension(isize,jsize),  intent(in)    :: total_precip,  &
                                                   parcel_rise
real,    dimension(isize,jsize,nlev_lsm),             &
                                  intent(in)    :: temperature_forcing, &
                                                   moisture_forcing
real,    dimension(isize,jsize,nlev_lsm,ntr),        &
                                  intent(in)    :: tracers
type(donner_cape_type),           intent(inout) :: Don_cape
type(donner_conv_type),           intent(inout) :: Don_conv
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     exit_flag      logical variable indicating whether deep convection
!                    is present or not in each column
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     parcel_rise    accumulated vertical displacement of a 
!                    near-surface parcel as a result of the lowest
!                    model level omega field [ Pa ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!
!   intent(inout) variables:
!
!     Don_cape       donner_cape type derived type variable containing
!                    diagnostics related to the cape calculation assoc-
!                    iated with the donner convection parameterization
!     Don_conv       donner_convection_type derived type variable con-
!                    taining diagnostics describing the nature of the 
!                    convection produced by the donner parameterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer :: k, n, kcont      ! do-loop indices
      integer :: i, j, unit

      ermesg= ' ' ; error = 0

!--------------------------------------------------------------------
!    determine if deep convection exists in any of the columns in the 
!    window for which column diagnostics were requested.
!--------------------------------------------------------------------
      do n=1,Col_diag%ncols_in_window

!--------------------------------------------------------------------
!    determine if deep convection exists in any of the columns in the 
!    window for which column diagnostics were requested. if deep
!    convection is present, output a multitude of values; if deep con-
!    vection is not present, cycle to check the next diagnostics 
!    column in the window.
!--------------------------------------------------------------------
        if (.not. exit_flag(Col_diag%i_dc(n), Col_diag%j_dc(n))) then
          i = Col_diag%i_dc(n)
          j = Col_diag%j_dc(n)
          unit = Col_diag%unit_dc(n)

!---------------------------------------------------------------------
!    output the pressures at lifting condensation level (plcl), at the
!    level of free convection (plfc), and at the level of zero buoyancy
!    (plzb).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
                'in donner_deep: plcl ', Don_cape%plcl(i,j)
          write (unit, '(a, e20.12)')  & 
                 'in donner_deep: plfc ', Don_cape%plfc(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: plzb ', Don_cape%plzb(i,j)

!---------------------------------------------------------------------
!    output the lag time value of cape (xcape_lag), the convective 
!    inhibition (coin), the time tendency of cape (dcape) and the lag
!    time column integrated water vapor (qint_lag).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: xcape ',   &
                              Don_cape%xcape_lag(i,j)
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: coin ', Don_cape%coin(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: dcape ', Don_conv%dcape(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: qint ',  Don_cape%qint_lag(i,j)

!---------------------------------------------------------------------
!    output the total cloud fractional area (a1), the maximum allowed
!    value for a1 (amax), the maximum cloud fractional area based on the
!    moisture constraint (amos), the total precipitation from the col-
!    umn (total_precip), the mesoscale cloud fractional area (ampta1),
!    the displacement of a parcel from its initial location due to 
!    accrued upward motion at the current time (parcel_rise), and the
!    convective precipitation rate (cell_precip).
!---------------------------------------------------------------------
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: a1   ', Don_conv%a1(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: amax ', Don_conv%amax(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: amos ', Don_conv%amos(i,j)
          write (unit, '(a, e20.12)')  & 
            'in donner_deep: tprea1 ', total_precip(i,j)
          write (unit, '(a, e20.12)')  & 
             'in donner_deep: ampta1 ', Don_conv%ampta1(i,j)
          write (unit, '(a, e20.12)')  & 
              'in donner_deep: omint', parcel_rise(i,j)
          write (unit, '(a, e20.12)')  & 
               'in donner_deep: rcoa1 ', Don_conv%cell_precip(i,j)

!---------------------------------------------------------------------
!    output various 3d fields between the specified highest index at
!    which diagnostics are to be output (kstart) and the nearest 
!    level to the surface (nlev_lsm), provided there has been some effect 
!    of deep convection at the level.
!---------------------------------------------------------------------
          do k=Col_diag%kstart,nlev_lsm
            if (temperature_forcing (i,j,k) == 0.0) cycle
            write (unit, '(a, i4)')'in donner_deep: k = ', k
            write (unit, '(a, e20.12)')  &
                 'in donner_deep: cemetf output to calling routine',  &
                     temperature_forcing(i,j,k)             
            write (unit, '(a, e20.12)')  &
                    'in donner_deep:TOTAL convective cemetf',  &
                     Don_conv%conv_temp_forcing(i,j,k)            
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: ceefc ',     &
                               Don_conv%ceefc(i,j,k)             
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: cecon ',  &
                                Don_conv%cecon(i,j,k)
            write (unit, '(a, e20.12)')  &
                     'in donner_deep: cemfc ',   &
                               Don_conv%cemfc(i,j,k)
            write (unit, '(a, e20.12)')  &
                 'in donner_deep: cememf output to calling routine',  &
                        moisture_forcing(i,j,k)               
            write (unit, '(a, e20.12)')  &
                     'in donner_deep: TOTAL convective cememf',  &
                        Don_conv%conv_moist_forcing (i,j,k)            
            write (unit, '(a, e20.12)')  &
                      'in donner_deep: cememf_mod',  &
                                Don_conv%cememf_mod(i,j,k)            
            write (unit, '(a, e20.12)')  &
                       'in donner_deep: cual  ',  &
                                  Don_conv%cual(i,j,k)              
            write (unit, '(a, e20.12)')  &
                       'in donner_deep: fre   ',   &
                                  Don_conv%fre(i,j,k)             
            write (unit, '(a, e20.12)')  &
                        'in donner_deep: elt   ',  &
                                     Don_conv%elt(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: cmus  ',    &
                                     Don_conv%cmus(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: ecds ',   &
                                    Don_conv%ecds(i,j,k)             
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: eces  ', &
                                     Don_conv%eces(i,j,k)            
            write (unit, '(a, e20.12)')  &
                         'in donner_deep: emds  ',  &
                                     Don_conv%emds(i,j,k)             
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: emes  ',  &
                                     Don_conv%emes(i,j,k)              
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: qmes  ',  &
                                    Don_conv%mrmes(i,j,k)            
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: wmps  ', &
                                      Don_conv%wmps(i,j,k)            
            write (unit, '(a, e20.12)')  &
                          'in donner_deep: wmms  ',  &
                                       Don_conv%wmms(i,j,k)            
            write (unit, '(a, e20.12)')  &
                           'in donner_deep: tmes  ',  &
                                        Don_conv%tmes(i,j,k)            
            write (unit, '(a, e20.12)')  &
                             'in donner_deep: dmeml ',   &
                                       Don_conv%dmeml(i,j,k)            
            write (unit, '(a, e20.12)')  &
                              'in donner_deep: uceml ',  &
                                       Don_conv%uceml(i,j,k)            
            write (unit, '(a, e20.12)')  &
                              'in donner_deep: detmfl ',  &
                                      Don_conv%detmfl(i,j,k)            
            write (unit, '(a, e20.12)')  &
                               'in donner_deep: umeml ',   &
                                      Don_conv%umeml(i,j,k)            

!---------------------------------------------------------------------
!    output various tracer-related fields for each tracer transported
!    by donner_deep_mod.
!---------------------------------------------------------------------
            do kcont=1,ntr     
              write (unit, '(a, e20.12)')  &
                              'in donner_deep: xgcm1 ',   &
                            tracers(i,j,k,kcont)                 
              write (unit, '(a, e20.12)')  &
                               'in donner_deep: qtren1 ',  &
                              Don_conv%qtren1(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                'in donner_deep: qtmes1 ',  &
                              Don_conv%qtmes1(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                'in donner_deep: temptr',  &
                              Don_conv%temptr(i,j,k,kcont)
              write (unit, '(a, e20.12)')  &
                                   'in donner_deep: qtceme ',   &
                               Don_conv%qtceme(i,j,k,kcont)             
              write (unit, '(a, e20.12)')  &
                                  'in donner_deep: wtp1 ',   &
                                Don_conv%wtp1(i,j,k,kcont)            
            end do
          end do  ! (k loop)
        endif
      end do    ! (n loop)

!--------------------------------------------------------------------


end subroutine don_d_column_end_of_step_k
 



!#####################################################################

subroutine don_d_convert_profile_k     &
         (name_hi, name_lo, n_lo, n_hi, ncc, profile_hi, press_hi, ptop,&
          include_set_value, include_sbl, include_conservation_factor, &
          set_value, sbl, conservation_factor, press_lo, diag_unit,  & 
          debug_ijt, profile_lo, ermesg, error)

!----------------------------------------------------------------------
!    subroutine don_d_convert_profile_k takes an input profile 
!    (profile_hi) associated with a character string name_hi on the 
!    hi-res model grid (press_hi) containing ncc_ens levels and extending
!    to a pressure level ptop and maps it to variable profile_lo assoc-
!    iated with character string name_lo on the lo-res model grid defined
!    by press_lo.
!    additonally, if desired, the integral of the profile on the lo-res
!    grid multiplied by conservation_factor may be set to set_value by 
!    modifying the profile below cloud base, or a specified sub-cloud 
!    source (sbl) may be added to the lo-res profile.
!    if column diagnostics are desired (debug_ijt), they are output to
!    diag_unit.
!-----------------------------------------------------------------------

implicit none

character(len=*),      intent(in)  :: name_hi, name_lo
integer,               intent(in)  :: n_lo, n_hi, ncc
real, dimension(n_hi), intent(in)  :: profile_hi, press_hi
real,                  intent(in)  :: ptop
logical,               intent(in)  :: include_set_value, include_sbl, &
                                      include_conservation_factor
real,                  intent(in)  :: set_value, sbl
real, dimension(n_lo), intent(in)  :: conservation_factor, press_lo
integer,               intent(in)  :: diag_unit
logical,               intent(in)  :: debug_ijt
real, dimension(n_lo), intent(out) :: profile_lo
character(len=*),      intent(out) :: ermesg
integer,               intent(out) :: error

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       name_hi       character string associated with input profile
!       name_lo       character string associated with output profile
!       n_lo          number of levels on lo-res grid
!       n_hi          number of levels on hi_res grid
!       ncc           number of layers in input profile that are affected
!                     by presence of cloud; it may be called with 
!                     ncc_kou for each ensemble member 
!                     (from define_lo_res_model_profiles or 
!                     add_to_ensemble_sum_hires), or with ncc_ens
!                     (from define_ensemble_profiles).
!       profile_hi    vertical profile on hi-res model grid
!       press_hi      full pressure levels of hi-res model [ Pa ]
!       ptop          pressure one level above cloud top  [ Pa ]
!       include_set_value
!                     it is desired to force the column integral to a
!                     specified value on the lo-res grid ?
!       include_sbl   it is desired to add a specified value to the
!                     profile in the layers below cloud base ?
!       include_conservation_factor
!                     the integrand which is to be set to set_value 
!                     includes a non-unity factor which multiplies the 
!                     profile ?
!       set_value     value desired for the integral of the
!                     output profile times conservation_factor      
!       sbl           value to be added to the profile in all layers
!                     below cloud base
!       conservation_factor
!                     the column integral of the product of the profile 
!                     and conservation_factor arrays is required to equal
!                     set_value
!       press_lo      interface pressure levels of lo-res model [ Pa ]
!       diag_unit     unit number for column diagnostics file
!       debug_ijt     column diagnostics are desired for this column ?
!
!   intent(out) variables:
!
!       profile_lo    vertical profile on lo-res model grid
!       ermesg        error message produced by any kernel routines
!                     called by this subroutine 
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
  
      real, dimension(n_lo) :: out       ! intermediate lo-res profile 
                                         ! after either setting column
                                         ! integral or adding boundary 
                                         ! layer source
      real, dimension(n_lo) :: conservation_factor_used                
                                         ! conservation_factor array 
                                         ! used in calculation; is array
                                         ! of 1.0 when 
                                         ! include_conservation_factor
                                         ! is .false.
      real                  :: intgl_hi  ! column integral of profile_hi 
      real                  :: intgl_lo  ! column integral of profile_lo
      integer               :: k         ! do-loop index

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!---------------------------------------------------------------------
!    if column diagnostics are desired, output a diagnostic message 
!    indicating the variable that is being processed.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)')  &
           'in mulsub: map_hi_res_col_to_lo_res_col: ' // trim(name_hi)
      endif

!----------------------------------------------------------------------
!   call don_u_map_hires_c_to_lores_c_k to map the 
!   profile from the hi-res model grid to the lo-res model grid.
!----------------------------------------------------------------------
      call don_u_map_hires_c_to_lores_c_k     &
          (n_lo, ncc+1, profile_hi(1:ncc+1), press_hi(1:ncc+1),  &
           ptop, press_lo, profile_lo, intgl_hi, intgl_lo, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
       if (error /= 0 ) then
        return
      endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals of the
!    profiles on both the hi- and lo-res grids.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 2e20.12)')  &
             'in mulsub: rintsum(' // trim(name_lo) // ' ) =',  &
                                              intgl_hi, intgl_lo

!---------------------------------------------------------------------
!    call don_u_compare_integrals_k to assess if the integrals
!    from the two grids are "equal", as they should be.
!---------------------------------------------------------------------
        call don_u_compare_integrals_k    &
                         (intgl_hi, intgl_lo, diag_unit, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) then
          return
        endif
      endif

!----------------------------------------------------------------------
!    if it is desired to set the integral value for the profile (i.e.,
!    set_value does not equal dummy_set_value), execute the following
!    code.
!----------------------------------------------------------------------
      if (include_set_value) then
        if (include_conservation_factor) then
          conservation_factor_used(:) = conservation_factor(:)
        else
          conservation_factor_used(:) = 1.0                   
        endif
        
!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrands at each
!    level on the lo-res grid.
!---------------------------------------------------------------------
        if (debug_ijt) then
          do k=1,n_lo                   
            if (profile_lo(k) /= 0.0) then
              write (diag_unit, '(a, i4, e20.12)') &
                 'in set_col_integral: k,phr,phr+= ', k, profile_lo(k)* &
                              conservation_factor_used(k)
            endif
          end do
        endif

!-----------------------------------------------------------------------
!    call don_u_set_column_integral_k to adjust the output
!    profile below cloud base so that the desired integral value is
!    obtained.
!-----------------------------------------------------------------------
        call don_u_set_column_integral_k    &
               (n_lo, profile_lo*conservation_factor_used, press_hi(1), &
                press_lo(1), set_value, press_lo, intgl_hi,     &
                intgl_lo, out, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) then
          return
        endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after the adjustment to the desired value.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, e20.12)')  &
                           'in set_col_integral: column(in)= ',intgl_hi
          write (diag_unit, '(a, e20.12)')  &
                          'in set_col_integral: column(out)= ',intgl_lo 
          do k=1,n_lo                 
            if (profile_lo(k)*conservation_factor_used(k) /= out(k)) then
              write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,qtr(in), qtr(out)= ', k,  &
                      profile_lo(k)*conservation_factor_used(k), out(k)
            endif
          end do
        endif

!---------------------------------------------------------------------
!    define the adjusted output profile by removing conservation_factor.
!---------------------------------------------------------------------
        profile_lo(:) = out(:)/conservation_factor_used(:)
      endif !(set_value /= dummy_set_value)

!----------------------------------------------------------------------
!    if a boundary layer source is to be added to the profile, execute
!    the following code.
!----------------------------------------------------------------------
      if (include_sbl .and. sbl /= 0.0) then

!----------------------------------------------------------------------
!    call don_u_apply_integral_source_k to apply the imposed 
!    subcloud source (sbl) to the input profile profile_out, resulting 
!    in the output profile out.  also returned are the column integrals
!    of the input profile (intgl_in) and the integral of the output
!    profile (intgl_out).
!    NOTE: in the original code, the subcloud source was not applied in 
!    the non-entropy case anywhere, and in the entropy case only to the 
!    model layer containing cloud base. 
!    I have MODIFIED THE CODE so that the value is APPLIED FROM SFC TO
!    TOP OF SPECIFIED REGION (CLOUD BASE) IS THIS CORRECT AND WHAT WAS
!    INTENDED ?
!----------------------------------------------------------------------
        call don_u_apply_integral_source_k     &
             (n_lo, profile_lo, press_hi(1), press_lo(1), sbl,  &
              press_lo, intgl_hi, intgl_lo, out, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
        if (error /= 0 ) then
          return
        endif

!---------------------------------------------------------------------
!    if column diagnostics are desired, output the integrals and 
!    profiles, both before and after adding the boundary layer source.
!---------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, 3e20.12)')  &
             'after apply_subcloud: column(in)= ',   &
                             intgl_hi, press_lo(1), press_hi(1)
          write (diag_unit, '(a, e20.12)')  &
                           'after apply_subcloud: column(out)= ',  &
                                                       intgl_lo 
          do k=1,n_lo                  
            if (profile_lo(k) /= out(k)) then
              write (diag_unit, '(a, i4, 2e20.12)') &
               'in set_col_integral: k,qtr(in), qtr(out)= ', k,  &
                                       profile_lo(k), out(k)
            endif
          end do
        endif

!----------------------------------------------------------------------
!    define the output profile on the lo-res model grid to be returned to
!    the calling routine.
!----------------------------------------------------------------------
        profile_lo(:) = out(:)
      endif  !(sbl /= 0.0)

!---------------------------------------------------------------------


end subroutine don_d_convert_profile_k



!#####################################################################

subroutine don_d_def_ensemble_profs_k    &
         (nlev_lsm, nlev_hires, ncc_ens, diag_unit, debug_ijt, ptt,  &
          cld_press, alp, detmfh, ucemh, cuql, cuqli, phalf_c,  &
          ensmbl_cloud_area, cuql_v, cuq, detmfl, uceml, ermesg, error)


!---------------------------------------------------------------------
!    subroutine don_d_def_ensemble_profs_k defines vertical 
!    profiles of cloud area, cloud ice, cloud liquid, vertical mass flux
!    and detrained vertical mass flux produced by the entire cumulus 
!    ensemble on the lo-res grid. 
!---------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------
integer,                        intent(in)   :: nlev_lsm, nlev_hires,&
                                                ncc_ens, diag_unit
logical,                        intent(in)   :: debug_ijt
real,                           intent(in)   :: ptt
real,    dimension(nlev_hires), intent(in)   :: cld_press, alp,  &
                                                detmfh, ucemh, cuql,  &
                                                cuqli
real,    dimension(nlev_lsm+1), intent(in)   :: phalf_c
real,    dimension(nlev_lsm),   intent(out)  :: ensmbl_cloud_area,  &
                                                cuql_v, cuq, detmfl, &
                                                uceml     
character(len=*),               intent(out)  :: ermesg
integer,                        intent(out)  :: error

!--------------------------------------------------------------------
!   intent(in) variables:
!
!        ptt        pressure one cloud model level above the ensemble
!                   cloud top [ Pa ]
!        nlev_lsm   number of levels in the low resolution grid
!        nlev_hires       number of levels in the high resolution grid
!        ncc_ens    cloud top index for the ensemble on the hi-res grid
!        cld_press  pressures at hi-res model levels [ Pa ]
!        alp        cloud area profile on hi-res model grid
!                   [ fraction ]
!        detmfh     detrained mass flux (layer above index level)
!                   (on cloud-model grid) (index 1 at cloud base)
!                   [ kg (air) / (sec m**2) ]
!        ucemh      upward mass flux on cloud model grid            
!                   [ kg (air) / (sec m**2) ]
!        cuql       ice water profile on cloud model grid; on input is
!                   normalized by total grid box area, on output is
!                   normalized by ensemble cloud area.
!        cuqli      liquid water profile on cloud model grid; on input 
!                   is normalized by total grid box area, on output is
!                   normalized by ensemble cloud area.
!        phalf_c    pressure at lo-res model half levels [ Pa ]
!        diag_unit  unit for column diagnostics output
!        debug_ijt  are column diagnostics desired in this column ?
!
!   intent(out) variables:
!
!        ensmbl_cloud_area  
!                   total cloud area profile over all ensemble members
!                   on large_scale model grid [ fraction ]
!        cuql_v     liquid water profile on large-scale model grid, 
!                   normalized by ensemble cloud area.
!        cuq        ice water profile on large-scale model grid, 
!                   normalized by ensemble cloud area.
!        detmfl     detrained mass flux on large-scale model grid
!                   normalized by ensemble cloud area
!                   index 1 near large-scale model base
!                   [ kg (air) / (sec m**2) ]
!        uceml      upward mass flux on large_scale model grid       
!                   [ kg (air) / (sec m**2) ]
!        ermesg     error message produced by any kernel routines
!                   called by this subroutine 
!
!--------------------------------------------------------------------

      real, dimension (nlev_lsm) :: conv_fact

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      conv_fact = 0.0

!----------------------------------------------------------------------
!    call don_d_convert_profile_k to map the ensemble-total cloud 
!    area profile from the cloud model grid (alp) to the large-scale 
!    model grid (ensmbl_cloud_area).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('alp', 'cual', nlev_lsm, nlev_hires, ncc_ens, alp, cld_press, &
          ptt, .false., .false., .false.,  0.0, 0.0, conv_fact, &
          phalf_c, diag_unit, debug_ijt, ensmbl_cloud_area, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total condensed ice 
!    profile from the cloud model grid (cuql) to the large-scale model 
!    grid (cuq).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('cuql', 'cuq', nlev_lsm, nlev_hires, ncc_ens, cuql, cld_press,&
          ptt, .false., .false., .false.,  0.0, 0.0, conv_fact, &
          phalf_c, diag_unit, debug_ijt, cuq, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total condensed liquid
!    profile from the cloud model grid (cuql) to the large-scale model 
!    grid (cuq).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('cuqli', 'cuql_v', nlev_lsm, nlev_hires, ncc_ens, cuqli, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,  &
          conv_fact, phalf_c, diag_unit, debug_ijt, cuql_v, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    call convert_profile to map the ensemble-total upward mass flux
!    profile from the cloud model grid (ucemh) to the large-scale model 
!    grid (uceml).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('ucemh', 'uceml', nlev_lsm, nlev_hires, ncc_ens, ucemh, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,   &
          conv_fact, phalf_c, diag_unit, debug_ijt, uceml, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    call convert_profile to map the ensemble-total detrained mass flux
!    profile from the cloud model grid (detmfh) to the large-scale model 
!    grid (detmfl).
!----------------------------------------------------------------------
      call don_d_convert_profile_k    &
         ('detmfh', 'detmfl', nlev_lsm, nlev_hires, ncc_ens, detmfh, &
          cld_press, ptt, .false., .false., .false.,  0.0, 0.0,  &
          conv_fact, phalf_c, diag_unit, debug_ijt, detmfl, ermesg, error)
 
!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------

end subroutine don_d_def_ensemble_profs_k
    

!#####################################################################

subroutine don_d_def_lores_model_profs_k               &
         (nlev_lsm, nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt,  &
          Nml, Param, pb, ptt, sfc_vapor_flux_c, sfc_sh_flux_c,   &
          sfc_tracer_flux_c, pfull_c, phalf_c, cld_press, tcc, dpf,  &
          dpftr, dfr, cld_evap, qlw, emfhr, efchr, etfhr, cell_freeze, &
          evap_rate, h1_liq, h1_ice,  ci_liq_cond, ci_ice_cond, h1_2, &
          q1, qtr, wetdepl, ermesg, error)
 
!---------------------------------------------------------------------
!    subroutine don_d_def_lores_model_profs_k maps vertical
!    profiles of various fields from the cloud-model grid to the large-
!    scale model grid. also, if desired, the sub-cloud base model levels
!    of the lo-res profiles may be modified so that the column integral 
!    equals a prescribed value (set_value), and / or a given value may
!    be assigned to the sub-cloud base levels.
!    this routine is called for each ensemble member individually, so 
!    that the input and output profiles are weighted by the cloud area 
!    of the ensemble member.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type

implicit none

!---------------------------------------------------------------------
integer,                            intent(in)    :: nlev_lsm,   &
                                                     nlev_hires,  &
                                                     ntr, ncc_kou, &
                                                     diag_unit
logical,                            intent(in)    :: debug_ijt
type(donner_param_type),            intent(in)    :: Param
type(donner_nml_type),              intent(in)    :: Nml
real,                               intent(in)    :: pb, ptt,  &
                                                     sfc_vapor_flux_c, &
                                                     sfc_sh_flux_c
real,    dimension(ntr),            intent(in)    :: sfc_tracer_flux_c
real,    dimension(nlev_lsm),       intent(in)    :: pfull_c
real,    dimension(nlev_lsm+1),     intent(in)    :: phalf_c
real,    dimension(nlev_hires),     intent(in)    :: cld_press, tcc,  &
                                                     dpf, &
                                                     dfr, cld_evap, qlw,&
                                                     emfhr, efchr
real,    dimension(nlev_hires,ntr), intent(in)    :: etfhr, dpftr
real,    dimension(nlev_lsm),       intent(out)   :: cell_freeze, &
                                                     evap_rate,      &
                                                     h1_liq, h1_ice, &
                                                     h1_2, q1
real,                               intent(out)   :: ci_liq_cond, &
                                                     ci_ice_cond
real,    dimension(nlev_lsm,ntr),   intent(out)   :: qtr, wetdepl
character(len=*),                   intent(out)   :: ermesg
integer,                            intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       nlev_lsm         number of levels on lo-res grid
!       nlev_hires             number of levels on hi_res grid
!       ntr              number of tracers being transported by the 
!                        donner deep convection parameterization
!       ncc_kou          number of layers in hi-res profile that are 
!                        affected by the presence of cloud
!       pb               pressure at cloud base [ Pa ]
!       ptt              pressure one model level above cloud top [ Pa ]
!       sfc_vapor_flux_c flux of water vapor from the surface into the 
!                        sub cloud layer [ kg(h2o) / (m**2 sec) ]
!       sfc_sh_flux_c    flux of sensible heat from the surface into the
!                        sub cloud layer [ W / m**2, or kg / m**3 ] 
!       sfc_tracer_flux_c  flux of tracer from the surface into the sub-
!                        cloud layer  [ kg(tracer) / (m**2 sec) ]
!       pfull_c           pressure at lo-res model full levels [ Pa ]
!       phalf_c         pressure at lo-res model half levels [ Pa ]
!       cld_press        pressure at hi-res model full levels [ Pa ]
!       dpf              condensation rate profile
!                        on hi-res grid  for the current ensemble member
!                        [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!      dpftr            wet-deposition rate profile
!                        on hi-res grid  for the current ensemble member        ,
!                        weighted by ratio of area to area at cloud base
!                        [ ( kg(tracer) ) / ( kg( dry air) sec ) ] 
!                        index 1 at physical base of cloud
!       dfr              profile of                     moisture tendency
!                        due to freezing in the convective updraft 
!                        on hi-res grid  for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       cld_evap                             profile of the potential
!                        cloud water evaporation for the curent ensemble
!                        member on th ehi-res grid. this amount of water!
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       qlw              profile of cloud water for the current ensemble
!                        member [ kg(h2o) / kg(air) ]
!       emfhr            vertical moisture flux convergence profile on 
!                        hi-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       efchr            vertical entropy flux convergence profile on
!                        hi-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       etfhr            vertical tracer flux convergence profile on
!                        hi-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       diag_unit        unit number of column diagnostics output file
!       debug_ijt        logical indicating whether diagnostics are 
!                        requested for this column 
!
!   intent(out) variables:
!
!       cell_freeze      profile of cloud-area-weighted moisture tendency
!                        due to freezing in the convective updraft 
!                        on lo-res grid for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       evap_rate        cloud-area-weighted profile of the potential
!                        cloud water evaporation for the current ensemble
!                        member on the lo-res grid. this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       h1                                   condensation rate profile
!                        on lo-res grid for the current ensemble member
!                        [ (      kg(h2o) ) / ( kg( dry air) sec ) ] 
!       h1_2             vertical entropy flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       q1               vertical moisture flux convergence profile on 
!                        lo-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       qtr              vertical tracer flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       wetdepl          wet-deposition rate on lo-res grid for the
!                        current ensemble member,
!                        weighted by ratio of area to area at cloud base
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!                        vertical index 1 at base of model
!       ermesg           error message produced by any kernel routines
!                        called by this subroutine 
!
!---------------------------------------------------------------------


!-----------------------------------------------------------------------
!   local variables:

      real, dimension (nlev_lsm)  ::   pi     ! inverse exner function
                                              ! used for setting column
                                              ! integral value (conserv-
                                              ! ation of theta)
      real, dimension (nlev_lsm)  ::   condensate ! liquid water profile on 
      real, dimension (nlev_lsm)  ::   conv_fact  
      real, dimension (nlev_hires) ::  dpftra
      real, dimension (nlev_hires) ::  dpf_warm, dpf_cold
      real, dimension (nlev_lsm)  ::   wetdepa
      integer                 ::   kcont, k,n ! do-loop indices
      real                    ::   sbl        ! value to be used for
                                              ! the profile at levels
                                              ! below cloud base
      real                    ::   set_value  ! desired column integral
                                              ! value 
      real  :: aak, dp

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

      conv_fact = 0.0
      dpf_cold = 0.
      dpf_warm = 0.

!--------------------------------------------------------------------
!    call don_d_convert_profile_k to map the tracer tendency due
!    to wet deposition from the cloud-model grid (dpftr) to the
!    large-scale model grid (wetdep) 
!--------------------------------------------------------------------
      do n=1,ntr
        dpftra(:) = dpftr(:,n)
        call don_d_convert_profile_k  &
             ('dpftra', 'wetdepa', nlev_lsm, nlev_hires, ncc_kou,   &
              dpftra(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,   &
              .false., .false., .false., 0.0, 0.0, conv_fact,   &
              phalf_c, diag_unit, debug_ijt, wetdepa, ermesg, error)
        wetdepl(:,n) = wetdepa(:)
      end do

!--------------------------------------------------------------------
!    call don_d_convert_profile_k to map the moisture tendency due
!    to freezing from the cloud model grid (dfr) to the large-scale model
!    grid (cell_freeze).
!--------------------------------------------------------------------
      call don_d_convert_profile_k   &
           ('DFR', 'frea', nlev_lsm, nlev_hires, ncc_kou,   &
            dfr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,   &
            .false., .false., .false., 0.0, 0.0, conv_fact,   &
            phalf_c, diag_unit, debug_ijt, cell_freeze, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    map the cloud condensate (qlw) from the cloud model to the large-
!    scale model (condensate). this field is only used for diagnostic 
!    purposes.
!----------------------------------------------------------------------
      if (debug_ijt) then
        call don_d_convert_profile_k    &
             ('QLW', 'evap', nlev_lsm, nlev_hires, ncc_kou,   &
              qlw(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, & 
              .false., .false., .false., 0.0, 0.0, conv_fact,&
              phalf_c, diag_unit, debug_ijt, condensate, ermesg, error)
      endif
      
!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    map the rate at which condensate which has not precipitated out
!    must evaporate from the cloud model grid (cld_evap) to the lo-res
!    model grid (evap_rate).
!---------------------------------------------------------------------
      call don_d_convert_profile_k    &
           ('QLW', 'evap_rate', nlev_lsm, nlev_hires, ncc_kou,   &
            cld_evap(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, &
            .false., .false., .false., 0.0, 0.0, conv_fact,&
            phalf_c,   diag_unit,  debug_ijt, evap_rate, ermesg, error)

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    if in diagnostics column, output profiles of the cloud evaporation
!    rate (cld_evap) and evaporation(qlw) on the hi-res model grid. 
!    cld_evap will be the actual evaporation rate if there turns out to 
!    be  no mesoscale circulation in the column, while qlw is the profile
!    of liquid water produced by the given ensemble member.
!----------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a)') 'in mulsub: P & CLD_EVAP'
        do k=1,ncc_kou-1
            write (diag_unit, '(a, i4, 2e20.12)')  &
                 'in mulsub: k, P & QLW', k, cld_press(k), cld_evap (k)
        end do
        write (diag_unit, '(a)') 'in mulsub: P & QLW'
        do k=1,ncc_kou-1
            write (diag_unit, '(a, i4, 2e20.12)')  &
                 'in mulsub: k, P & QLW', k, cld_press(k), qlw      (k)
        end do
      endif
      
      do k=1, nlev_hires
        if (tcc(k) > Param%tfre) then
          dpf_warm(k) = dpf(k)
        else
          dpf_cold(k) = dpf(k)
        endif
      end do

!----------------------------------------------------------------------
!    map the condensation rate from the cloud model (-dpf) to the 
!    large-scale model (h1). h1 is a term appropriate for use in the
!    water vapor equation; i.e., condensation is a loss term.
!----------------------------------------------------------------------
      call don_d_convert_profile_k                 &
           ('RLHR_warm', 'h1_liq', nlev_lsm, nlev_hires, ncc_kou,    &
            -dpf_warm(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,  &
            .false., .false., .false., 0.0, 0.0, conv_fact,&
            phalf_c, diag_unit, debug_ijt, h1_liq, ermesg, error)

      if (debug_ijt) then
          dp = cld_press(1) - cld_press(2)
        aak = -dpf_warm(1)*0.5*dp
        do k=2,ncc_kou 
          dp = cld_press(k) - cld_press(k+1)
          aak = aak - dpf_warm(k)*dp
        end do
        aak = aak/(Param%grav*1000.)
            write (diag_unit, '(a, e20.12, i6)')  &
                                 'in mulsub: dpf_warm intg, ncc_kou',  &
                                           aak, ncc_kou
      endif

        ci_liq_cond = 0.
        do k=1,nlev_lsm
          dp = phalf_c(k) - phalf_c(k+1)
          ci_liq_cond = ci_liq_cond + h1_liq(k)*dp
        end do
        ci_liq_cond = ci_liq_cond/(Param%grav      )
      if (debug_ijt) then
            write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: h1_liq intg', ci_liq_cond
      endif

!----------------------------------------------------------------------
!    map the condensation rate from the cloud model (-dpf) to the 
!    large-scale model (h1). h1 is a term appropriate for use in the
!    water vapor equation; i.e., condensation is a loss term.
!----------------------------------------------------------------------
      call don_d_convert_profile_k                 &
           ('RLHR_cold', 'h1_ice', nlev_lsm, nlev_hires, ncc_kou,    &
            -dpf_cold(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,  &
            .false., .false., .false., 0.0, 0.0, conv_fact,&
            phalf_c, diag_unit, debug_ijt, h1_ice, ermesg, error)

      if (debug_ijt) then
          dp = cld_press(1) - cld_press(2)
        aak = -dpf_cold(1)*0.5*dp
        do k=2,ncc_kou 
          dp = cld_press(k) - cld_press(k+1)
          aak = aak - dpf_cold(k)*dp
        end do
        aak = aak/(Param%grav*1000.)
            write (diag_unit, '(a, e20.12, i6)')  &
                                 'in mulsub: dpf_cold intg, ncc_kou',  &
                                           aak, ncc_kou
      endif
        ci_ice_cond = 0.
        do k=1,nlev_lsm
          dp = phalf_c(k) - phalf_c(k+1)
          ci_ice_cond = ci_ice_cond + h1_ice(k)*dp
        end do
        ci_ice_cond = ci_ice_cond/(Param%grav      )
      if (debug_ijt) then
            write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: h1_ice intg', ci_ice_cond
      endif

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    determine the vertical flux convergence of each tracer.
!---------------------------------------------------------------------
      do kcont=1,ntr  

!----------------------------------------------------------------------
!    calculate the imposed subcloud tracer-flux convergence (sbl) in 
!    units of kg(tracer) per kg(dry air) per sec. define set_value so
!    that the column integrated tracer flux convergence will be set to
!    zero.
!----------------------------------------------------------------------
        sbl = (sfc_tracer_flux_c(kcont)*Param%grav)/(phalf_c(1) - pb)
        set_value = 0.0

!----------------------------------------------------------------------
!    call convert_profile to map the vertical tracer flux convergence 
!    from the cloud model (etfhr) to the large-scale model grid (qtr). 
!    force the column integral of the flux convergence to be 0.0; then 
!    add the imposed sub-cloud convergence.
!----------------------------------------------------------------------
        call don_d_convert_profile_k     &
             ('qtrv', 'qtr', nlev_lsm, nlev_hires, ncc_kou,    &
              etfhr(1:ncc_kou+1,kcont), cld_press(1:ncc_kou+1), ptt,   & 
              .true., .true., .false., set_value, sbl, conv_fact, &
              phalf_c, diag_unit, debug_ijt, qtr(:,kcont),ermesg, error)
      end do

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    define the subcloud moisture flux convergence (sbl) in units of
!    kg(h2o) per kg(air) per sec. sfc_vapor_flux_c is the imposed bound-
!    ary layer moisture source in units of kg(h2o) per m**2 per sec. 
!----------------------------------------------------------------------
      sbl = (sfc_vapor_flux_c*Param%grav)/(phalf_c(1) - pb)
      set_value = 0.0

!----------------------------------------------------------------------
!    call don_d_convert_profile_k to map the vertical moisture 
!    flux convergence from the cloud model (emfhr) to the lo-res model 
!    grid (q1). force the column integral of the flux convergence to be 
!    0.0; then add the imposed sub-cloud convergence.
!----------------------------------------------------------------------
      call don_d_convert_profile_k       &
           ('EMFHR', 'q1', nlev_lsm, nlev_hires, ncc_kou,    &
            emfhr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt,   &
            .true., .true., .false., set_value, sbl, conv_fact, &
            phalf_c, diag_unit,  debug_ijt, q1, ermesg, error)
      if (debug_ijt) then
        aak = 0.
        do k=1,nlev_lsm
          dp = phalf_c(k) - phalf_c(k+1)
          aak = aak + q1(k)*dp
        end do
        aak = aak/(Param%grav*1000.)
            write (diag_unit, '(a, e20.12)')  &
                                 'in mulsub: q1 intg', aak
      endif
        

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------
!    calculate the subcloud entropy flux convergence (sbl) in units of
!    deg K per sec. sfc_sh_flux_c is the imposed boundary layer sensible
!    heat source in units of watts per square meter. define the inverse
!    exner function so that an integral constraint on theta may be
!    applied.
!----------------------------------------------------------------------
      sbl = Param%grav*sfc_sh_flux_c/((phalf_c(1) - pb)*Param%cp_air)
      set_value = 0.0
      do k=1,nlev_lsm               
        pi(k) = (1.0e05/pfull_c(k))**Param%kappa
      end do

!----------------------------------------------------------------------
!    map the temperature change due to vertical entropy flux convergence
!    from the cloud model (efchr) to the large-scale model grid (h1_2). 
!    force the column integral of the temperature change to be 0.0, thus
!    conserving enthalpy; then add any imposed sub-cloud convergence.
!----------------------------------------------------------------------
   if (Nml%frc_internal_enthalpy_conserv) then
      call don_d_convert_profile_k       &
           ('EFCHR', 'h1_2', nlev_lsm, nlev_hires, ncc_kou,   &
            efchr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, &
            .true., .true., .false., set_value, sbl, conv_fact, &
            phalf_c, diag_unit,  debug_ijt, h1_2, ermesg, error)
   else

!----------------------------------------------------------------------
!    map the temperature change due to vertical entropy flux convergence
!    from the cloud model (efchr) to the large-scale model grid (h1_2). 
!    force the column integral of the flux convergence times inverse 
!    exner function (i.e., theta) to be 0.0, thus conserving entropy;
!    then add any imposed sub-cloud convergence.
!----------------------------------------------------------------------

      call don_d_convert_profile_k       &
           ('EFCHR', 'h1_2', nlev_lsm, nlev_hires, ncc_kou,   &
            efchr(1:ncc_kou+1), cld_press(1:ncc_kou+1), ptt, &
            .true., .true., .true., set_value, sbl, pi, &
            phalf_c, diag_unit,  debug_ijt, h1_2, ermesg, error)
   endif

!---------------------------------------------------------------------
!    determine if an error message was returned from the kernel
!    routines. if so, return to calling routine.
!---------------------------------------------------------------------
      if (error /= 0 ) return

!---------------------------------------------------------------------
!    if in diagnostics column, output the profile of the entropy flux
!    convergence on the lo-res grid, at levels where it is non-zero.
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (debug_ijt) then
          if (h1_2(k) /= 0.0) then
            write (diag_unit, '(a, i4, e20.12)')  &
                                 'in mulsub: JK,H1= ', k, h1_2(k)
          endif
        endif
      end do

!----------------------------------------------------------------------



end subroutine don_d_def_lores_model_profs_k


!#####################################################################

subroutine don_d_add_to_ensmbl_sum_hires_k     &
         (nlev_hires, ntr, ncc_kou, diag_unit, debug_ijt, area_ratio,  &
          cfracice, rcl, flux, emfhr, dpf, qlw, etfhr, cuql, cuqli, &
          ucemh, alp, rlsm, emsm, detmfh, etsm, ermesg, error)

!-----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_hires_k adds the contrib-
!    utions from this ensemble member to various profiles on the hi-res
!    grid that are being summed over the ensemble.
!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------
integer,                         intent(in   )  :: nlev_hires, ntr, ncc_kou,&
                                                   diag_unit
logical,                         intent(in   )  :: debug_ijt
real,                            intent(in   )  :: area_ratio
real, dimension(nlev_hires),     intent(in   )  :: cfracice, rcl, flux, &
                                                   emfhr, dpf, qlw     
real, dimension(nlev_hires,ntr), intent(in   )  :: etfhr
real, dimension(nlev_hires),     intent(inout)  :: cuql, cuqli, ucemh, &
                                                   alp, rlsm, emsm,  &
                                                   detmfh
real, dimension(nlev_hires,ntr), intent(inout)  :: etsm         
character(len=*),                intent(  out)  :: ermesg
integer,                         intent(  out)  :: error
!---------------------------------------------------------------------
!   intent(in) variables:
!
!      nlev_hires            number of levels on hi_res grid
!      ntr             number of tracers being transported by the 
!                      donner deep convection parameterization
!      ncc_kou         number of layers in hi-res profile that are 
!                      affected by the presence of cloud
!      area_ratio      ratio of cloud base area of this ensemble member
!                      to that of ensemble member # 1. (ensemble member
!                      # 1 assumed to have largest cloud base area)
!                      [ dimensionless ]
!      cfracice        fraction of condensate that is ice [ fraction ]
!      rcl             profile of cloud radius for this ensemble member
!                      [ m ]
!      flux            upward mass flux profile in cloud for this
!                      ensemble member [ kg (air) / sec ]
!      emfhr           vertical moisture flux convergence for this
!                      ensemble member [ kg (h2o) / ( kg(dry air) sec ) ]
!      dpf             condensation rate profile on hi-res grid for the 
!                      current ensemble member
!                      [ kg(h2o) / ( kg( dry air) sec ) ] 
!      qlw             profile of cloud water for the current ensemble
!                      member [ kg(h2o) / kg(air) ]
!      etfhr           vertical tracer flux convergence profile on
!                      hi-res grid for the current ensemble member 
!                      [ kg(tracer) / ( kg(dry air) sec ) ]
!      debug_ijt       logical indicating whether diagnostics are 
!                      requested for this column 
!      diag_unit       unit number of column diagnostics output file
!
!   intent(inout) variables:
!
!      cuql            vertical profile on the hi-res grid of condensed 
!                      ice, summed over ensemble members # 1 to the cur-
!                      rent, each member's contribution being weighted by
!                      its cloud area at level k relative to the cloud 
!                      base area of ensemble member #1
!                      [ kg(h2o) / kg (dry air) ]
!      cuqli           vertical profile on the hi-res grid of condensed 
!                      liquid, summed over ensemble members # 1 to the 
!                      current, each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      [ kg(h2o) / kg (dry air) ]
!      ucemh           vertical profile on the hi-res grid of cell upward
!                      mass flux, summed over ensemble members # 1 to the
!                      current, each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      [ kg (dry air) / ( m**2 sec ) ]
!      alp             vertical profile on the hi-res grid of cloud area
!                      summed over ensemble members # 1 to the current, 
!                      each member's contribution being weighted
!                      by its cloud area at level k relative to the cloud
!                      base area of ensemble member #1
!                      as a result, the cloud area profile is expressed
!                      relative to the cloud base area of ensemble member
!                      # 1. [ dimensionless ]
!      rlsm            vertical profile on the hi-res grid of conden-
!                      sation rate, summed over ensemble members # 1 to
!                      the current, each member's contribution being 
!                      weighted by its cloud area at level k relative to 
!                      the cloud base area of ensemble member #1
!                      [ ( kg(h2o) ) / ( kg( dry air) sec ) ] 
!      emsm            vertical profile on the hi-res grid of vertical
!                      moisture flux convergence, summed over ensemble 
!                      members # 1 to the current, each member's contrib-
!                      ution being weighted by its cloud area at level k
!                      relative to the cloud base area of ensemble 
!                      member #1  [ kg (h2o) / ( kg(dry air) sec ) ]
!      detmfh          vertical profile on the hi-res grid of detrained
!                      mass flux in the layer above indexed level, summed
!                      over ensemble members # 1 to the current, each 
!                      member's contribution being weighted by its cloud 
!                      area at level k relative to the cloud base area of
!                      ensemble member #1 [ kg (dry air) / ( m**2 sec ) ]
!      etsm            vertical profile on the hi-res grid of vertical
!                      tracer flux convergence, summed over ensemble 
!                      members # 1 to the current, each member's contrib-
!                      ution being weighted by its cloud area at level k
!                      relative to the cloud base area of ensemble 
!                      member #1  [ kg (tracer) / ( kg(dry air) sec ) ]
!
!   intent(out) variables:
!
!      ermesg          error message produced by this subroutine or any
!                      kernel routines called by this subroutine 
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real       :: wt_factor   ! cloud area at level k for current 
                                ! ensemble member, relative to cloud
                                ! base area of ensemble member # 1
     integer     :: k, ktr      ! do-loop indices     

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!--------------------------------------------------------------------
!    add the contributions from this ensemble member to the arrays 
!    accumulating ensemble sums on the cloud model levels.
!--------------------------------------------------------------------
      do k=1,ncc_kou 

!----------------------------------------------------------------------
!    define the factor needed to normalize each ensemble member's con-
!    tribution by the cloud base area of ensemble member #1. wt_factor
!    is the cloud area at level k for ensemble member kou, relative to
!    the cloud area at cloud base (k=1) for ensemble member #1.
!-----------------------------------------------------------------------
        wt_factor = area_ratio*(rcl(k)/rcl(1))**2
        
!----------------------------------------------------------------------
!    add this ensemble member's appropriately weighted contribution to
!    the ensemble-total cloud area (alp), condensed ice (cuql), condensed
!    liquid (cuqli), cell upward mass flux (ucemh), cell detrained mass 
!    flux (detmfh), condensation rate (rlsm), vertical moisture flux 
!    convergence (emsm) and vertical tracer flux convergence (etsm). the
!    weighting factor area_ratio*(rcl(k)/rcl(1))**2 allows the contrib-
!    utions from each member to be added by normalizing each member's 
!    contribution by the cloud base area of ensemble member #1.
!    NOTE: several of the arrays already have some of the normalizing
!    factors already included and so here need only to be multiplied by 
!    a portion of wt_factor.
!----------------------------------------------------------------------
        alp(k)   = alp(k)   + wt_factor                      
        cuql(k)  = cuql(k)  + wt_factor*(cfracice(k)*qlw(k))
        cuqli(k) = cuqli(k) + wt_factor*((1.0 - cfracice(k))*qlw(k))
        ucemh(k) = ucemh(k) + area_ratio*flux(k)/(rcl(1)**2)
        if (k < ncc_kou) then
          if (flux(k+1) < flux(k)) then
            detmfh(k) = detmfh(k) + area_ratio*   &
                        ((flux(k)-flux(k+1))/(rcl(1)**2))
          endif
        endif
        rlsm(k)   = rlsm(k)   - area_ratio*dpf (k) 
        emsm(k)   = emsm(k)   + area_ratio*emfhr(k)

!----------------------------------------------------------------------
!    if in a diagnostics column, output the total cell upward mass flux 
!    (ucemh) and the cloud area at each level ( alp).
!----------------------------------------------------------------------
        if (debug_ijt) then
          write (diag_unit, '(a, i4, 2e20.12)')  &
                    'in mulsub: k,ucemh, alp= ',k,ucemh(k), alp(k)
        endif
      end do

!----------------------------------------------------------------------
!    add this ensemble member's appropriately weighted contribution to
!    the vertical tracer flux convergence (etsm). the weighting factor 
!    area_ratio allows the contributions from each member to be added by
!    normalizing each member's contribution by the cloud base area of 
!    ensemble member #1.
!----------------------------------------------------------------------
      do ktr=1,ntr
        do k=1,ncc_kou 
          etsm(k,ktr) = etsm(k,ktr) + area_ratio*etfhr(k,ktr)
        end do
      end do
      
!---------------------------------------------------------------------


end subroutine don_d_add_to_ensmbl_sum_hires_k   
 



!#####################################################################

subroutine don_d_add_to_ensmbl_sum_lores_k      &
         (nlev_lsm, ntr, diag_unit, debug_ijt, lmeso, &
                          frz_frac, Param, Nml,   &
          area_ratio, dint, cell_freeze,         cell_melt, wetdepl, &
          temp_c, h1_2, ecd, ecd_liq, ecd_ice, &
          ece, ece_liq, ece_ice, evap_rate, q1,     h1_liq, h1_ice, &
          pfull_c, meso_melt, meso_freeze, &
          phalf_c, qtr, ensmbl_melt, ensmbl_melt_meso, ensmbl_freeze,&
          ensmbl_freeze_meso, ensmbl_wetc,  &
          meso_frac, precip_frac, frz_frac_non_precip, disz, &
          disz_remelt, disp_melted, disze1, disze2, disze3,  &
          disp_liq, disp_ice, enctf, encmf, enev, disg_liq, disg_ice, &
          disb, disc_liq, disc_ice, dism_liq, dism_liq_frz, &
          dism_liq_remelt, dism_ice, dism_ice_melted, &
          ecds_liq, ecds_ice, eces_liq, eces_ice, disd, &
          disv, qtren, ermesg, error, meso_frz_intg, melting_in_cloud, &
          precip_melt_frac, meso_frz_frac)

!-----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_lores_k adds the contrib-
!    utions from this ensemble member to various profiles on the lo-res
!    grid that are being summed over the ensemble.
!-----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type

implicit none

!-----------------------------------------------------------------------
integer,                       intent(in   ) :: nlev_lsm, ntr, diag_unit
logical,                       intent(in   ) :: debug_ijt, lmeso
real,                          intent(in)   :: frz_frac
type(donner_param_type),       intent(in)    :: Param
type(donner_nml_type),       intent(in)    :: Nml  
real,                          intent(in   ) :: area_ratio, dint
real, dimension(nlev_lsm),     intent(in   ) :: cell_freeze, cell_melt, &
                                                temp_c, h1_2, ecd, ece, &
                                            ecd_liq, ecd_ice, ece_liq,&
                                              ece_ice, &
                                                evap_rate, q1,     &
                                                h1_liq, h1_ice, &
                                                pfull_c, meso_melt,   &
                                                meso_freeze
real, dimension(nlev_lsm+1),   intent(in   ) :: phalf_c  
real, dimension(nlev_lsm,ntr), intent(in   ) :: qtr, wetdepl
real,                          intent(in)    :: meso_frac, precip_frac
real,        intent(inout) ::                    frz_frac_non_precip
real, dimension(nlev_lsm),     intent(inout) :: ensmbl_melt,   &
                                                ensmbl_melt_meso, &
                                                ensmbl_freeze, enctf, &
                                                ensmbl_freeze_meso, &
                                                encmf, enev,       &
                                                disz_remelt, &
                                                disz,               &
                                                 disze1, disze2, &
                                                disze3,   &
                                                disg_liq, disg_ice, &
                                               disb,                   &
                                                ecds_liq, ecds_ice, &
                                                eces_liq, eces_ice, &
                                                disc_liq, disc_ice, &
                                                dism_liq, dism_ice, &
                                                dism_ice_melted, &
                                                dism_liq_frz, &
                                                dism_liq_remelt, &
                                                disp_liq, disp_ice, &
                                                disp_melted, &
                                                disd, disv
real, dimension(nlev_lsm,ntr), intent(inout) :: qtren, ensmbl_wetc
character(len=*),              intent(  out) :: ermesg
integer,                       intent(  out) :: error
logical,                       intent( in)   :: meso_frz_intg   
real,                          intent( in)   ::                &
                                                 meso_frz_frac
logical, intent(in) :: melting_in_cloud
real   , intent(in) ::            precip_melt_frac          

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     
!       nlev_lsm         number of levels on lo-res grid
!       ntr              number of tracers being transported by the 
!                        donner deep convection parameterization
!       area_ratio       ratio of cloud base area of this ensemble 
!                        member to that of ensemble member # 1. 
!                        (ensemble member # 1 assumed to have largest 
!                        cloud base area) [ dimensionless ]
!       dint             column sum of moisture tendency due to freezing
!                        in the convective updraft on hi-res grid for the
!                        current ensemble member
!!!!  CHECK ON THESE UNITS !!!!!
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!       cell_freeze      profile of cloud-area-weighted moisture tendency
!                        due to freezing in the convective updraft 
!                        on lo-res grid for the current ensemble member
!                        [ (      g(h2o) ) / ( kg(dry air) sec ) ] 
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       cell_melt        in-cloud melting of condensate associated with
!                        convective cells. made up of two parts, 1) that
!                        due to the freezing of liquid carried upwards
!                        in the cell updraft, 2) that due to the melting 
!                        of condensed ice that precipitates out. if meso
!                        circulation is present, this component is zero;
!                        melting will be determined in subroutine mesub.
!                        [ (      g(h2o) ) / ( kg(dry air) day ) ] 
!!   CHECK UNITS HERE !!!!
!!!!!!!!======>>>>>>>    NOTE UNITS OF g(h2o). (Verify and change.)
!       temp_c           temperature at model levels [ deg K ]
!       h1_2             vertical entropy flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ deg K / sec ]                        
!       ecd              profile of condensate evaporated in convective
!                        downdraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       ece              profile of condensate evaporated in convective
!                        updraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       evap_rate        cloud-area-weighted profile of the potential
!                        cloud water evaporation for the current ensemble
!                        member on the lo-res grid. this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       phalf_c          pressure at lo-res model interface levels [ Pa ]
!       q1               vertical moisture flux convergence profile on 
!                        lo-res grid for the current ensemble member 
!                        [ kg(h2o) / ( kg(dry air) sec ) ]
!       h1                                   condensation rate profile
!                        on lo-res grid for the current ensemble member
!                        [ (      kg(h2o) ) / ( kg( dry air) sec ) ] 
!       pfull_c          pressure on lo-res model full levels [ Pa ]
!       meso_melt        profile of condensate which is melted in meso-
!                        scale downdraft on large-scale model grid
!                        [ g(h2o) / kg(air) / day ]
!       meso_freeze      profile of condensate which is frozen upon 
!                        entering the anvil on the large-scale grid
!                        [ g(h2o) / kg(air) / day ]
!       qtr              vertical tracer flux convergence profile on
!                        lo-res grid for the current ensemble member 
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!       lmeso            a mesoscale circulation exists in the current
!                        grid box ?
!       wetdepl          wet deposition for current ensemble member,
!                        weighted by ratio of area to area at cloud base
!                        [ kg(tracer)/ ( kg sec) ]
!                        vertical index 1 at base of model
!       debug_ijt        logical indicating whether diagnostics are 
!                        requested for this column 
!       diag_unit        unit number of column diagnostics output file
!
!   intent(inout) variables:
!
!       ensmbl_melt      vertical profile on the lo-res grid of ice melt,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!       ensmbl_freeze    vertical profile on the lo-res grid of freezing,
!                        both from the cells and any mesoscale circul-
!                        ation, summed over ensemble members # 1 to the 
!                        current, each member's contribution being 
!                        weighted by its cloud area at level k relative !
!                        to the cloud base area of ensemble member #1
!                        [ kg(h2o) / kg (dry air) ]
!       ensmbl_wetc      vertical profile on the lo-res grid of wet
!                        deposition from cells, summed over ensemble
!                        members #1 to current, each members contributio        n
!                        weighted by the ratio of its area 
!                        to the area of ensemble member #1 at cloud base
!                        [ kg(tracer) / kg s ]
!                        vertical index 1 at model base
!       enctf            vertical profile on the lo-res grid of the entr-
!                        opy forcing, consisting of the sum of the
!                        vertical entropy flux convergence and the latent
!                        heat release, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area
!                        of ensemble member #1
!                        [ deg K / day ]                        
!       encmf            vertical profile on the lo-res grid of the      
!                        moisture forcing, consisting of the sum of the
!                        vertical moisture flux convergence and the cond-
!                        ensation, summed over ensemble members # 1 to 
!                        the current, each member's contribution being 
!                        weighted by its cloud area at level k relative 
!                        to the cloud base area of ensemble member #1
!                        [ (      kg(h2o) ) / ( kg( dry air) day ) ] 
!       enev             vertical profile on the lo-res grid of the      
!                        cloud-area-weighted profile of the potential
!                        cloud water evaporation, summed over ensemble 
!                        members # 1 to the current, each member's con-
!                        tribution being weighted by its cloud area at !
!                        level k relative to the cloud base area of 
!                        ensemble member #1.  this amount of water
!                        must be evaporated if it turns out that there is
!                        no mesoscale circulation generated in the 
!                        column.
!                        [ (      kg(h2o) ) / ( kg(dry air) sec ) ] 
!       disg             vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation
!                        associated with the evaporation of condensate
!                        in the convective downdraft and updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!       disb             vertical profile on the lo-res grid of the      
!                        temperature flux convergence, summed over 
!                        ensemble members # 1 to the current, each mem-
!                        ber's contribution being weighted by its cloud 
!                        area at level k relative to the cloud base area 
!                        of ensemble member #1.  
!                        [ deg K / day ] 
!       disc             vertical profile on the lo-res grid of the      
!                        latent heat term in the temperature equation, 
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ deg K / day ] 
!       ecds             vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective downdraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!       eces             vertical profile on the lo-res grid of the      
!                        condensate evaporated in convective updraft,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ g(h2o) / kg(air) / day ]
!       disd             vertical profile on the lo-res grid of the      
!                        vertical moisture flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        lo-res grid for the current ensemble member 
!                        [  g(h2o) / ( kg(dry air) day ) ]
!       qtren            vertical profile on the lo-res grid of the      
!                        vertical tracer flux convergence,
!                        summed over ensemble members # 1 to the current,
!                        each member's contribution being weighted by its
!                        cloud area at level k relative to the cloud base
!                        area of ensemble member #1.  
!                        [ kg(tracer) / ( kg(dry air) sec ) ]
!
!    intent(out) variables:
!
!       ermesg           error message produced by any kernel routines
!                        called by this subroutine 
!
!---------------------------------------------------------------------- 

!--------------------------------------------------------------------
!   local variables:


      real, dimension (nlev_lsm) :: rlh  
                                     !  condensation term in temperature
                                     !  equation on lo-res grid for cur-
                                     !  rent ensemble member 
                                     !  [ deg K / day ]
      real, dimension (nlev_lsm) :: cmf 
                                     !  forcing term for moisture 
                                     !  equation on lo-res grid for
                                     !  current ensemble member; sum 
                                     !  of vertical flux convergence 
                                     !  and condensation terms 
                                     !  [ g(h2o) / ( kg(air) day ) ]

     real     :: convrat   !  latent heat factor, appropriate for the 
                           !  temperature at a given model level 
                           !  ( = L / cp ) [ deg K ]
     real     :: qtrsum    !  sum of tracer flux convergence over all 
                           !  tracers and all levels and all ensemble 
                           !  members up to the current. used as a diag-
                           !  nostic; sum should be 0.0, if no boundary 
                           !  source term.
                           !  [ kg(tracer) / ( kg(dry air) sec ) ]
     integer  :: k, kcont  !  do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!--------------------------------------------------------------------
!    sum up various cloud-base-area weighted contributions to vertical
!    profiles on the large-scale grid that are being summed over the 
!    ensemble.
!--------------------------------------------------------------------
      do k=1,nlev_lsm       

!---------------------------------------------------------------------
!    define the moisture forcing term (sum of condensation h1 and 
!    vertical flux convergence q1) on the large-scale grid. convert to
!    units of g(h20) per kg(air) per day, requiring multiplication by
!    1.0e3 g(h2o) /kg(h2o) times SECONDS_PER_DAY. add this member's 
!    contribution to the sum over the ensemble (encmf). 
!----------------------------------------------------------------------
        cmf(k) = (-(h1_liq(k) + h1_ice(k)) + q1(k))*  &
                                     (Param%SECONDS_PER_DAY*1.0e03)
        encmf(k) = encmf(k) + area_ratio*cmf(k)

!----------------------------------------------------------------------
!   define the cell precipitation forms.
!   disz : frozen liquid condensate
!   disz_remelt: frozen liquid condensate which then melts
!   disp_liq: liquid condensate
!   disp_ice: frozen condensate
!   disp_melted: frozen condensate which melts
!----------------------------------------------------------------------
        if (.not. melting_in_cloud) then
          disz(k) = disz(k) + area_ratio*h1_liq(k)* &
                    Param%seconds_per_day*     frz_frac*precip_frac
        else  ! (not melting in cloud)
          disz_remelt(k) = disz_remelt(k) + area_ratio*h1_liq(k)*  &
                           Param%seconds_per_day*frz_frac*precip_frac
        endif  ! (not melting in cloud)

        disp_liq(k) = disp_liq(k) + area_ratio*h1_liq(k)*   &
                      Param%seconds_per_day*(1.0-frz_frac)* &
                      precip_frac              
        disp_ice(k) = disp_ice(k) + area_ratio*h1_ice(k)*  &
                      Param%seconds_per_day*           &
                      (1.0-precip_melt_frac)*precip_frac
        disp_melted(k) = disp_melted(k) + area_ratio*h1_ice(k)*  &
                         Param%seconds_per_day*           &
                         (precip_melt_frac)*precip_frac

!---------------------------------------------------------------------
!    define the heating rate due to liquid and ice condensation.
!---------------------------------------------------------------------
        disc_liq(k) = disc_liq(k) + area_ratio*h1_liq(k)* &
                      Param%seconds_per_day*Param%hlv/Param%cp_air
        disc_ice(k) = disc_ice(k) + area_ratio*h1_ice(k)* &
                      Param%seconds_per_day*Param%hls/Param%cp_air

!---------------------------------------------------------------------
!    define the  heating rates associated with the evaporation of
!    non-precipitating condensate. these terms are used when a meso-
!    scale circulation is not present; the _chgd variables are used
!    when the mesoscale circulation is present for the initial ensemble
!    member, but is not sustainable by one of the remaining members.
!---------------------------------------------------------------------
        disze1(k) = disze1(k) + area_ratio*h1_liq(k)*  &
                    Param%seconds_per_day*Param%hls*  &
                    frz_frac           *(1.-precip_frac)/Param%cp_air
        disze2(k) = disze2(k) + area_ratio*h1_liq(k)*  &
                    Param%seconds_per_day*Param%hlv*  &
                 (1.0-frz_frac           )*(1.-precip_frac)/Param%cp_air
        disze3(k) = disze3(k) - area_ratio*h1_ice(k)* &
                    Param%seconds_per_day*Param%hls*   &
                    (1.-precip_frac)/Param%cp_air

!----------------------------------------------------------------------
!    define the components of precipitation from the cell condensate 
!    that was transferred to the mesoscale circulation.
!    dism_liq: precip from liquid condensate that does not freeze
!              occurs when no melting and no freezing in cloud OR 
!              is the liquid condensate which does not freeze when 
!              there is melting in the cloud
!    dism_liq_frz: precip from liquid condensate that freezes
!              occurs when there is no melting but is freezing
!    dism_ice:  precip from frozen condensate
!               occurs when there is no melting
!    dism_ice_melted: precip from frozen condensate that melts
!               occurs when there is melting
!    dism_liq_remelt: precip from liquid condensate that freezes and 
!                     then melts
!               occurs when there is both melting and freezing 
!----------------------------------------------------------------------
        if (.not. melting_in_cloud) then
          if (.not. (meso_frz_intg) .and. frz_frac_non_precip == 0.) then
            dism_liq(k) = dism_liq(k) + area_ratio*h1_liq(k)*  &
                          meso_frac*Param%seconds_per_day                          
          else
            dism_liq_frz(k) = dism_liq_frz(k) + area_ratio*h1_liq(k)* &
                          meso_frac*Param%seconds_per_day
          endif
          dism_ice(k) = dism_ice(k) + area_ratio*h1_ice(k)*  &
                        meso_frac*Param%seconds_per_day
        else   ! (not melting in cloud)
          dism_liq(k) = dism_liq(k) + area_ratio*h1_liq(k)* &
                meso_frac*(1.-frz_frac_non_precip)*Param%seconds_per_day
          dism_ice_melted(k) = dism_ice_melted(k) +  &
                                    area_ratio*h1_ice(k)*  &
                        meso_frac*Param%seconds_per_day
          dism_liq_remelt(k) = dism_liq_remelt(k) +      &
                               area_ratio*h1_liq(k)*  &
                 meso_frac*frz_frac_non_precip*Param%seconds_per_day
        endif

        if (debug_ijt) then
          write (diag_unit, '(a, i4, 3e20.12)')  &
                              'in mulsub: precip profiles', &
                               k, disp_liq(k)*Param%hlv/Param%cp_air, &
                             Param%hls/Param%cp_air*disp_ice(k), &
                             Param%hls/Param%cp_air*disz(k)
          write (diag_unit, '(a, i4, 2e20.12)')  &
                           'in mulsub: remelt, melt precip profiles', &
                           k, Param%hlv/Param%cp_air*disz_remelt(k), &
                                  Param%hlv/Param%cp_air*disp_melted(k)
          write (diag_unit, '(a, i4, 3e20.12)')  &
                              'in mulsub: evap   profiles', &
                               k, disze1(k), disze2(k),  -disze3(k)     
          write (diag_unit, '(a, i4, 2e20.12)')  &
                              'in mulsub: cd     profiles', &
                               k, disc_liq(k), disc_ice(k)
        endif

!--------------------------------------------------------------------
!    add this member's weighted contribution to the ensemble's temper-
!    ature flux convergence (disb), the ensemble's water vapor flux 
!    convergence (disd) and the ensemble's entropy flux convergence 
!    (enctf). convert the rates to units of per day, and for disd from
!    kg(h2o) per kg(air) to g(h2o) per kg(air).
!--------------------------------------------------------------------
        disb(k) = disb(k) + area_ratio*(h1_2(k)*Param%SECONDS_PER_DAY)
        disd(k) = disd(k) + area_ratio*(q1(k)*  &
                             (Param%SECONDS_PER_DAY*1.0e3))
        disv(k) = disv(k) + area_ratio*((h1_liq(k) + h1_ice(k))* &
                             (Param%SECONDS_PER_DAY*1.0e3))
!   change enctf to reflect need for both ice and liq cd in layer of
!   tfre
        enctf(k) = enctf(k) + area_ratio*    &
                       (h1_2(k)*Param%SECONDS_PER_DAY + &
                       (h1_liq(k)*Param%hlv + h1_ice(k)*Param%hls)*  &
                           Param%seconds_per_day/ Param%cp_air )

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add this member's contribution
!    to the mesoscale condensate's evaporation associated with convect-
!    ive downdrafts (ecds) and that associated with evaporation into 
!    the environment (eces). if there has been no freezing associated
!    with the mesoscale condensate, define the condensation term for
!    the temperature equation using the latent heat of vaporization
!    (disg). if there has been freezing, then the convective downdraft 
!    heating uses the latent heat of vaporization, whereas the entrain-
!    ment evaporation is of ice and so uses the latent heat of 
!    sublimation.
!--------------------------------------------------------------------
        if (lmeso) then
          ecds_liq(k) = ecds_liq(k) + area_ratio*ecd_liq(k)
          ecds_ice(k) = ecds_ice(k) + area_ratio*ecd_ice(k)
          eces_liq(k) = eces_liq(k) + area_ratio*ece_liq(k)
          eces_ice(k) = eces_ice(k) + area_ratio*ece_ice(k)

!---------------------------------------------------------------------
!    evaporation of the frozen condensate removes hls.
!---------------------------------------------------------------------
          disg_ice(k) = disg_ice(k) - area_ratio* &
                          ((ecd_ice(k) + ece_ice(k))*  &
                                Param%hls/(Param%CP_AIR*1000.))

!---------------------------------------------------------------------
!    if there has been cell freezing and either melting or meso 
!    freezing, then the liquid condensate evaporated in the environment 
!    would have previously frozen and so removes hls. the liquid con-
!    densate evaporated in the downdraft would not have been frozen and
!    so removes hlv.
!---------------------------------------------------------------------
          if (dint /= 0. .and. &       
                (melting_in_cloud .or. meso_frz_intg)) then
            disg_liq(k) = disg_liq(k) - area_ratio*  &
                    ((ecd_liq(k)*Param%HLV + ece_liq(k)*Param%hls)/   &
                                           (Param%CP_AIR*1000.))

          else

!---------------------------------------------------------------------
!    if there has been no freezing in the cells or freezing in the 
!    cells only and no melting, then the evaporation of liquid conden-
!    sate in both environment and downdraft removes hlv.
!---------------------------------------------------------------------
            disg_liq(k) = disg_liq(k) - area_ratio*  &
                          ((ecd_liq(k) + ece_liq(k))*  &
                                Param%hlv/(Param%CP_AIR*1000.))
          endif
        endif    ! (lmeso)

!---------------------------------------------------------------------
!    add this member's cloud water evaporation rate to the sum over 
!    the ensemble (enev).
!---------------------------------------------------------------------
        enev(k) = enev(k) + area_ratio*evap_rate(k)

!--------------------------------------------------------------------
!    if a mesoscale circulation exists, add the appropriately-weighted
!    anvil freezing and melting terms to the arrays accumulating their 
!    sums over the ensemble (ensmbl_melt, ensmbl_freeze). if in a diag-
!    nostic column, output the anvil (meso_freeze) and ensemble-sum
!    (ensmbl_freeze) freezing profiles.
!--------------------------------------------------------------------
        if (lmeso) then
          ensmbl_melt_meso(k) = ensmbl_melt_meso(k) - area_ratio*  &
                           meso_melt(k)
          ensmbl_freeze_meso(k) = ensmbl_freeze_meso(k) + area_ratio* &
                                  meso_freeze(k)
          if (debug_ijt) then
            if (meso_freeze(k) /= 0.0) then
              write (diag_unit, '(a, i4, 2e20.12)')  &
                              'in mulsub: jk,fres,fre= ',   &
                               k, ensmbl_freeze_meso(k), meso_freeze(k)
            endif
          endif
        endif

!--------------------------------------------------------------------
!    add the appropriately-weighted convective cell freezing and 
!    melting terms to the arrays accumulating vertical profiles of 
!    total cloud melting (ensmbl_melt) and freezing (ensmbl_freeze) 
!    over the entire ensemble.  if in diagnostic column, output the 
!    convective cell (cell_freeze) and accumulated (ensmbl_freeze) 
!    freezing profiles.
!--------------------------------------------------------------------
        ensmbl_freeze(k) = ensmbl_freeze(k) + area_ratio*cell_freeze(k)
        ensmbl_melt(k) = ensmbl_melt(k) - area_ratio*cell_melt(k)
        if (debug_ijt) then
          if (cell_freeze(k) /= 0.0) then
            write (diag_unit, '(a, i4, 2e20.12)')  &
                     'in mulsub: jk,fres,frea= ',    &
                                     k, ensmbl_freeze(k), cell_freeze(k)
          endif
        endif
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, initialize a variable to sum the 
!    pressure-weighted tracer flux convergence, summed over all tracers
!    and all levels, for this ensemble member. upon completion of the 
!    loop, qtrsum should be 0.0, if there are no imposed tracer sources 
!    or sinks.
!---------------------------------------------------------------------
      if (debug_ijt) then
        qtrsum = 0.
        do k=1,nlev_lsm       
          do kcont=1,ntr  
            write (diag_unit, '(a,  i4, 2e20.12)')  &
                      'in mulsub: jk,    qtr,qtren=         ', &
                              k,    qtr(k,kcont),qtren(k,kcont)
            qtrsum = qtrsum + qtr(k,kcont)*(phalf_c(k) - phalf_c(k+1))
            write (diag_unit, '(a,  i4, e20.12)')  &
                         'in mulsub: jk,    qtrsum= ', k,    qtrsum
          end do
        end do
      endif 

!--------------------------------------------------------------------
!   add this ensemble member's appropriately-weighted contributions 
!   to the tracer flux convergence profiles being summed over the 
!   ensemble. 
!--------------------------------------------------------------------
      do k=1,nlev_lsm        
        ensmbl_wetc(k,:) = ensmbl_wetc(k,:) + area_ratio*wetdepl(k,:)
        qtren(k,:) = qtren(k,:) + area_ratio*qtr(k,:)
      end do

!----------------------------------------------------------------------
!    if in diagnostics column, output the profiles of the amount of
!    cloud  water evaporated (if lmeso is true) or the cloud water evap-
!    oration rate (if lmeso is false),  the array ctf (forcing for 
!    entropy eqn) and cmf (forcing for moisture equation) for this 
!    ensemble member.
!----------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm        
          if (h1_liq(k) /= 0.0  .or. h1_ice(k) /= 0.0 .or.  &
              h1_2(k) /= 0.0) then
            if (temp_c(k) >= Param%tfre) then
              convrat = Param%HLV/Param%CP_AIR
              rlh(k) = (h1_liq(k) + h1_ice(k))*Param%SECONDS_PER_DAY* &
                                                                convrat
            else
              convrat = Param%HLS/Param%CP_AIR
              rlh(k) = (h1_liq(k) + h1_ice(k))*Param%SECONDS_PER_DAY*  &
                                                                convrat
            endif
            write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & ctf', k, pfull_c(    k),  &
                                               h1_2(k)*86400. + rlh(k)
          endif
        end do
        do k=1,nlev_lsm        
          if (cmf(k) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10, e20.12)') &
                       'in mulsub: k, p & cmf', k, pfull_c(k), cmf(k)
          endif
        end do
      endif

!---------------------------------------------------------------------


end subroutine don_d_add_to_ensmbl_sum_lores_k 




!#####################################################################

subroutine don_d_add_to_ensmbl_sum_intgl_k        &
         (diag_unit, debug_ijt, lmeso,  area_ratio,  ca_liq, ca_ice,&
          frz_frac_non_precip, meso_frac, cell_precip, cu, &
          apt, ensmbl_precip, ensmbl_cond,                     &  
          ensmbl_anvil_cond_liq, ensmbl_anvil_cond_liq_frz, &
          meso_frz_intg, meso_frz_frac, ensmbl_anvil_cond_ice, &
          ensmbl_cld_top_area, ermesg, error)

!----------------------------------------------------------------------
!    subroutine don_d_add_to_ensmbl_sum_intgl_k adds the contrib-
!    utions from this ensemble member to various global integrals.
!----------------------------------------------------------------------

implicit none

!----------------------------------------------------------------------
integer,          intent(in   ) :: diag_unit
logical,          intent(in   ) :: debug_ijt, lmeso
real,             intent(in   ) :: area_ratio,     ca_liq, ca_ice, &
                                   frz_frac_non_precip, meso_frac, &
                                   cell_precip, cu, apt
real,             intent(inout) :: ensmbl_precip, ensmbl_cond, &
                                   ensmbl_cld_top_area, &
                                   ensmbl_anvil_cond_liq,  &
                                   ensmbl_anvil_cond_liq_frz,  &
                                   ensmbl_anvil_cond_ice
real,             intent(in   ) :: meso_frz_frac
logical,          intent(in   ) :: meso_frz_intg               
character(len=*), intent(  out) :: ermesg
integer,          intent(  out) :: error 
!----------------------------------------------------------------------
!   intent(in) variables:
!
!      area_ratio       ratio of cloud base area of this ensemble member
!                       to that of ensemble member # 1. (ensemble member
!                       # 1 assumed to have largest cloud base area)
!                       [ dimensionless ]
!      ca               rate of transfer of condensate from cell to 
!                       anvil for this ensemble member 
!                       [ mm / day ]
!      cell_precip      precipitation rate for this ensemble member
!                       [ mm / day ]
!      cu               condensation rate for this ensemble member
!                       [ mm / day ]
!      apt              ratio of cloud top area to cloud base area
!                       for this ensemble member [ dimensionless ]
!      lmeso            logical indicating if mesoscale circulation 
!                       is present
!      debug_ijt        logical indicating whether diagnostics are 
!                       requested for this column 
!      diag_unit        unit number of column diagnostics output file
! 
!   intent(inout) variables:
!
!      ensmbl_precip    sum of precipitation rate over ensemble members,
!                       # 1 to the current, weighted by the area at 
!                       cloud base of each member
!                       [ mm / day ]
!      ensmbl_cond      sum of condensation rate over ensemble members,
!                       # 1 to the current, weighted by the area at 
!                       cloud base of each member
!                       [ mm / day ]
!      ensmbl_anvil_cond 
!                       sum of rate of transfer of condensate from cell 
!                       to anvil over ensemble members, # 1 to the c
!                       current, weighted by the area at cloud base of 
!                       each member [ mm / day ]
!      ensmbl_cld_top_area  
!                       sum of the cloud top areas over ensemble members 
!                       # 1 to the current, normalized by the cloud base
!                       area of ensemble member # 1 [ dimensionless ]
!
!   intent(out) variables:
!
!      ermesg           error message produced by this subroutine or any
!                       kernel routines called by this subroutine 
!
!---------------------------------------------------------------------

      real :: local_frz_frac

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!--------------------------------------------------------------------
!    if a mesoscale circulation is present, add this member's cloud-
!    base_area-weighted contribution of condensate transferred to the 
!    anvil (ensmbl_anvil_cond) and cloud top cloud fraction 
!    (ensmbl_cld_top_area) to the arrays accumulating the ensemble sums.
!--------------------------------------------------------------------
      if (lmeso) then
        if (meso_frac /= 0.0) then
          local_frz_frac = (frz_frac_non_precip + meso_frz_frac)/  &
                                                              meso_frac
        else
          local_frz_frac = 0.0
        endif
        ensmbl_anvil_cond_liq   = ensmbl_anvil_cond_liq   +   &
                                  area_ratio*ca_liq*(1.-local_frz_frac)
        ensmbl_anvil_cond_liq_frz   = ensmbl_anvil_cond_liq_frz   +   &
                                        area_ratio*ca_liq*local_frz_frac
        ensmbl_anvil_cond_ice   = ensmbl_anvil_cond_ice   +  &
                                                   area_ratio*ca_ice
        ensmbl_cld_top_area = ensmbl_cld_top_area + area_ratio*apt
      endif

!--------------------------------------------------------------------
!    add this ensemble member's weighted contribution to the total 
!    precipitation (ensmbl_precip) and condensation (ensmbl_cond). 
!--------------------------------------------------------------------
      ensmbl_precip = ensmbl_precip + area_ratio*cell_precip
      ensmbl_cond   = ensmbl_cond   + area_ratio*cu

!---------------------------------------------------------------------



end subroutine don_d_add_to_ensmbl_sum_intgl_k 

!#####################################################################


!######################################################################

subroutine don_d_output_diag_profs_k    &
         (nlev_lsm, diag_unit, pfull_c, disc_liq, disc_ice, disb,  &
          disd, disn, encmf, temp_tend_freeze, temp_tend_freeze_meso, &
          temp_tend_melt, cmus_tot, emds_liq, &
          emds_ice,  emes_liq, emes_ice, wmms, &
          wmps, tmes, mrmes, eces_liq, eces_ice, ecds_liq, ecds_ice, &
          disa, dise, disg_2liq, disg_2ice, disf, &
          ermesg, error)

!---------------------------------------------------------------------
!    subroutine output_diagnostic_profiles prints out vertical profiles
!    of various donner_deep variables in those columns for which 
!    diagnostics have been requested.
!---------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------
integer,                      intent(in)   :: nlev_lsm, diag_unit
real,    dimension(nlev_lsm), intent(in)   :: pfull_c, disc_liq, &
                                              disc_ice, disb, disd,&
                                              disn, encmf,  &
                                              temp_tend_freeze,    &
                                              temp_tend_freeze_meso,  &
                                              temp_tend_melt, cmus_tot,&
                                                          wmms, wmps, &
                                              emds_liq, emds_ice, &
                                              emes_liq, emes_ice, &
                                              tmes, mrmes,             &
                                              ecds_liq, ecds_ice, &
                                              eces_liq, eces_ice, &
                                              disa, dise, disg_2liq,  &
                                              disg_2ice, disf 
character(len=*),              intent(out) :: ermesg
integer,                       intent(out) :: error

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     diag_unit
!     pfull_c
!     disc
!     disb
!     disd
!     disn
!     encmf
!     temp_tend_freeze
!     temp_tend_melt
!     cmus_tot
!     emds
!     emes
!     wmms
!     wmps
!     tmes
!     mrmes
!     eces
!     ecds
!     disa
!     dise
!     disg_2
!     disf
!
!----------------------------------------------------------------------

!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]



!-------------------------------------------------------------------
!   local variables:

      integer  ::  k      ! do-loop index
     
!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!  disc: cloud ensemble cell condensation heating rate [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ( (disc_liq(k)+ disc_ice(k))  /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, 2e20.12)')  &
              'in mulsub: k, P & liq/ice =',  k, pfull_c(k),disc_liq(k), &
                                                        disc_ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  disb  : cloud ensemble cell vertical entropy flux convergence 
!          [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disb(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
              'in mulsub: k, P & EFC =',  k, pfull_c(k),disb(k)
        endif
      end do

!---------------------------------------------------------------------
!  disd: cloud ensemble cell vertical moisture flux convergence 
!        [ kg(h2o)/ (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disd(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
               'in mulsub: k, P & EMF =',  k, pfull_c(k),disd(k)
        endif
      end do

!---------------------------------------------------------------------
!  disn : cloud ensemble cell thermal forcing 
!         [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disn(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & cell thermal forcing =',    &
                        k, pfull_c(k),disn(k)
        endif
      end do

!---------------------------------------------------------------------
!  encmf : cloud ensemble cell moisture forcing 
!          [ kg(h2o) / (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (encmf(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & cell moisture forcing =',    &
                        k, pfull_c(k),encmf(k)
        endif
      end do

!---------------------------------------------------------------------
!  temp_tend_freeze  : cloud ensemble temperature tendency due to   
!                      freezing [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((temp_tend_freeze(k) +   &
                          temp_tend_freeze_meso(k))/= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up freeze        =',    &
           k, pfull_c(k),temp_tend_freeze(k) + temp_tend_freeze_meso(k)
        endif
      end do

!---------------------------------------------------------------------
!  temp_tend_melt  : cloud ensemble plus mesoscale temperature 
!                    tendency due to melting [ deg K / sec ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (temp_tend_melt(k)  /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down melt        =',    &
                       k, pfull_c(k),temp_tend_melt(k) 
        endif
      end do

!---------------------------------------------------------------------
!  cmus_tot : water mass condensed in mesoscale updraft
!             [ kg(h2o) / (kg(air) sec) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (cmus_tot(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                      'in mulsub: k, P & meso up con           =',    &
                         k, pfull_c(k),cmus_tot  (    k)
        endif
      end do

!---------------------------------------------------------------------
!  emds : evaporation in mesoscale downdrafts.
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((emds_liq(k) + emds_ice(k)) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                     'in mulsub: k, P & meso down evap        =',    &
                       k, pfull_c(k), emds_liq(k) + emds_ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  emes : evaporation in mesoscale updrafts.
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((emes_liq(k) + emes_ice(k)) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                        'in mulsub: k, P & meso up evap        =',    &
                          k, pfull_c(k),emes_liq(k) + emes_ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  wmms : vapor transferred from cell to mesoscale circulation and 
!         then condensed [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (wmms(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & meso cell con       =',    &
                       k, pfull_c(k),wmms(k)
        endif
      end do

!---------------------------------------------------------------------
!  wmps : vapor transferred from cell to mesoscale circulation  
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (wmps(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso vap redist     =',    &
                          k, pfull_c(k),wmps(k)
        endif
      end do

!---------------------------------------------------------------------
!  tmes   : mesoscale temperature flux convergence                
!           [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (tmes(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                       'in mulsub: k, P & meso efc            =',    &
                        k, pfull_c(k),tmes(k)
        endif
      end do

!---------------------------------------------------------------------
!  mrmes   : mesoscale moisture flux convergence                
!            [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
!! WAS ORIGINALLY    :    if (tmes  (    k) /= 0.00 ) then
        if (mrmes(k) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)') &
                        'in mulsub: k, P & meso mfc            =',    &
                         k, pfull_c(k),mrmes(k)
        endif
      end do

!---------------------------------------------------------------------
!  eces   : sublimation in mesoscale updrafts               
!           [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((eces_liq(k) + eces_ice(k)) /= 0.00 ) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & up con evap         =',    &
                          k, pfull_c(k), eces_liq(k) + eces_ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  ecds_v : sublimation in mesoscale downdrafts             
!           [ kg(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((ecds_liq(k) + ecds_ice(k)) /= 0.00 ) then
           write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, P & down con evap         =',    &
                       k, pfull_c(k), ecds_liq(k) + ecds_ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  disa   : total temperature tendency due to deep convection
!           [ deg K / day ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disa(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                           'in mulsub: k, p & ens thermal forc', &
                              k, pfull_c(k), disa(k)
         endif
       end do

!---------------------------------------------------------------------
!  dise : total moisture tendency due to deep convection
!         [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (dise(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                          'in mulsub: k, p & ens moisture forc', &
                           k, pfull_c(k), dise(k)
        endif
      end do

!---------------------------------------------------------------------
!  disg2_v : total moisture tendency due to remaining cell condensate
!            [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if ((disg_2liq(k) + disg_2ice(k)) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                         'in mulsub: k, p & thermal modifications', &
                          k, pfull_c(k), disg_2liq(k) + disg_2ice(k)+ &
                         temp_tend_freeze(k) + temp_tend_melt(k)
          write (diag_unit, '(a, i4, f19.10, 2e20.12)')  &
            'in mulsub: k, p & thermal modifications -- liq, ice', &
                          k, pfull_c(k), disg_2liq(k),  disg_2ice(k)
        endif
      end do

!---------------------------------------------------------------------
!  disf_v : total moisture tendency due to evaporation in cell updrafts
!           [ g(h2o) / (kg(air) day) ]
!---------------------------------------------------------------------
      do k=1,nlev_lsm
        if (disf(k) /= 0.0) then
          write (diag_unit, '(a, i4, f19.10, e20.12)')  &
                      'in mulsub: k, p & moisture modifications', &
                       k, pfull_c(k), disf(k)
        endif
      end do

!---------------------------------------------------------------------



end subroutine don_d_output_diag_profs_k


!#####################################################################

subroutine don_d_def_conv_forcing_k  &
         (nlev_lsm, diag_unit, debug_ijt, lmeso, Initialized,   &
          pb, Param, Nml, ensmbl_precip, meso_precip, meso_cloud_area, &
          anvil_precip_melt, phalf_c, enev, encmf, ensmbl_freeze, &
          ensmbl_freeze_meso, enctf, disg_liq, disg_ice, ecds_liq, &
          ecds_ice, eces_liq, eces_ice, emds_liq, emds_ice, emes_liq,&
          emes_ice, mrmes, mrmes_up, mrmes_dn, tmes, tmes_up, tmes_dn, &
          wmps, ensmbl_cloud_area, ensmbl_melt, ensmbl_melt_meso, &
          pfull_c, temp_c, cmus_tot, wmms, disc_liq, disc_ice, &
          dism_liq, dism_liq_frz, dism_liq_remelt, dism_ice, &
          dism_ice_melted, meso_frz_intg_sum, disp_liq, disp_ice, &
          disb, disd, disv, total_precip_c, disz, disz_remelt, &
          disp_melted, disze1, disze2, disze3,              &
          disf, disn, dise, disa, cutotal, temp_tend_melt,&
          lprcp, liq_prcp, frz_prcp, vrt_mot, water_budget,  &
          n_water_budget, &
          enthalpy_budget, n_enthalpy_budget, precip_budget, &
          n_precip_paths, n_precip_types, ermesg, error, melting_in_cloud)

!---------------------------------------------------------------------
!    subroutine define_convective_forcing produces the effects of
!    the donner_deep parameterization on the large-scale flow, defining
!    the time tendency terms and integral quantities resulting from the
!    parameterization.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_budgets_type,  &
                             donner_initialized_type

implicit none

!---------------------------------------------------------------------
integer,                      intent(in)  :: nlev_lsm, diag_unit
logical,                      intent(in)  :: debug_ijt, lmeso 
type(donner_param_type),      intent(in)  :: Param
type(donner_initialized_type), intent(in)  :: Initialized
type(donner_nml_type),        intent(in)  :: Nml    
real,                         intent(in)  :: ensmbl_precip, meso_precip
real,                         intent(in)  :: pb                        
logical,                      intent(in)  :: meso_frz_intg_sum         
real,    dimension(nlev_lsm), intent(inout)  :: ensmbl_melt,      &
                                                ensmbl_melt_meso, &
                                                anvil_precip_melt, &
                                                ensmbl_freeze, &
                                                ensmbl_freeze_meso
real,    dimension(nlev_lsm + 1), intent(in)  :: phalf_c
real,    dimension(nlev_lsm), intent(in)  :: meso_cloud_area,  &
                                             enev, encmf, &
                                             disz_remelt, disz,       &
                                             disze1, disze2, disze3, &
                                             enctf, ecds_liq, ecds_ice,&
                                             eces_liq, eces_ice, &
                                             disg_liq, disg_ice, &
                                             mrmes,  &
                                             emds_liq, emds_ice, &
                                             emes_liq, emes_ice, &
                                             mrmes_up, mrmes_dn, tmes, &
                                             tmes_up, tmes_dn, &
                                             wmps, ensmbl_cloud_area,  &
                                             pfull_c,   &
                                             temp_c, cmus_tot, wmms,   &
                                             disb, disd, disv, &
                                             disc_liq, disc_ice, &
                                             dism_liq, dism_ice, &
                                             dism_ice_melted, &
                                             dism_liq_frz, &
                                             dism_liq_remelt, &
                                             disp_liq, disp_ice, &
                                             disp_melted
real,                         intent(out) :: total_precip_c
real,    dimension(nlev_lsm), intent(out) :: disf,         &
                                             disn, dise,  &
                                             disa, cutotal,   &
                                             temp_tend_melt, &
                                             liq_prcp, frz_prcp
real,                         intent(out) :: lprcp, vrt_mot
integer, intent(in) :: n_water_budget, n_enthalpy_budget, &
                       n_precip_paths, n_precip_types
real, dimension(nlev_lsm,n_water_budget), intent(out) :: &
                                                           water_budget
real, dimension(nlev_lsm,n_enthalpy_budget), intent(out) ::&
                                                       enthalpy_budget
real, dimension(nlev_lsm,n_precip_paths, n_precip_types),  &
                                              intent(out) ::&
                                                       precip_budget
character(len=*),             intent(out) :: ermesg
integer,                      intent(out) :: error
logical,                      intent( in) :: melting_in_cloud

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        diag_unit         i/o unit for column diagnostics output
!        ensmbl_precip
!        meso_precip       
!        lmeso              a mesoscale circulation is present in this
!                           column ?     
!        debug_ijt          column diagnostics are desired in this 
!                           column ?
!        meso_cloud_area
!        anvil_precip_melt
!        phalf_c            pressure at large-scale model half levels 
!                           [ Pa ]
!        enev
!        encmf
!        ensmbl_freeze
!        enctf
!        disg
!        ecds
!        eces
!        emds
!        emes
!        mrmes
!        tmes
!        wmps
!        ensmbl_cloud_area
!        ensmbl_melt         
!        pfull_c
!        temp_c
!        cmus_tot
!
!    intent(out) variables:
!
!        total_precip_c
!        disf
!        disg_2
!        disn
!        dise
!        disa
!        cutotal
!        temp_tend_melt
!        temp_tend_freeze
!
!---------------------------------------------------------------------



!!  UNITS
!!    ucemh  [kg /sec / m**2 ]
!!    conint [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    precip [ kg / sec ] ===> [ kg / sec / m**2 ]
!!    q1     [ kg(h2o) / kg(air) / sec ]
!!    h1     [ kg(h2o) / kg(air) / sec ]
!!    cmf    [ g(h2o) / kg(air) /day ]
!!    rlh    [ kg(h2o) / kg(air) / day ]  * [ L / Cp ] = [ deg K / day ]
!!    h1_2   [ deg K / sec ]
!!    efc    [ deg K / day ]
!!    efchr  [ deg K / sec ]
!!    ehfh   [ kg(air) (deg K) / (sec**3 m)
!!    ctf    [ deg K / day ]
!!    disb_v [ deg K / day ]
!!    disc_v [ deg K / day ] 
!!    disn   [ deg K / day ] 
!!    ecd    [ g(h2o) / kg(air) / day ]
!!    ece    [ g(h2o) / kg(air) / day ]
!!    ecds_v [ g(h2o) / kg(air) / day ]
!!    eces_v [ g(h2o) / kg(air) / day ]
!!    enctf  [ deg K / day ]
!!    encmf  [ g(h2o) / kg(air) /day ]
!!    pf     [ (m**2 kg(h2o)) / (kg(air) sec) ]
!!    dpf    [ (m**2 kg(h2o)) / (kg(air) sec) ] ==>   
!!                                          [ kg(h2o)) / (kg(air) sec) ]
!!    qlw2   [ kg(h2o)) / (kg(air) sec) ]
!!    qlw    [ kg(h2o)) / kg(air) ]
!!    evap   [ kg(h2o)) / kg(air) ]
!!    evap_rate [ kg(h2o)) / (kg(air) sec) ]
!!    disg   [ deg K / day ]


!--------------------------------------------------------------------
!   local variables:

      real,    dimension(nlev_lsm) :: disl_liq, &
                                      disl_ice, disl_ice_melted, &
                                      disl_liq_depo, disl_liq_cd, &
                                      disl_ice_depo, disl_ice_cd, &
                                      disga_liq_up, disga_liq_dn, &
                                      disga_ice_up, disga_ice_dn, &
                                      disg_2liq, disg_2ice, &
                                      cell_evap, meso_cd, meso_depo, &
                                      meso_evap
      real    ::   vsuma, vsumb, vsumc, vsumd, vsumd1, vsumd2, vsume, &
                   vsumf, vsumg, vsumg1, vsumg2, vsumh, vsumi, vsumi1,&
                   vsumi2
      real    ::   tsuma, tsumb, tsumiup, tsumidn 
      real    ::   tsumj, tsumj1, tsumk1, tsumk2, tsumk3
      real    ::   tsummliq, tsummice, tsummfliq
      real    ::   tsumaliq, tsumcliq, tsumdliq, tsumeliq,   &
                   tsumfliq, tsumg1liq, tsumg2liq
      real    ::   tsumaice, tsumcice, tsumdice, tsumeice,  &
                   tsumfice, tsumg1ice, tsumg2ice
      real    ::   liq_ice
      real    ::   esumc, sumf, summ, sumn
      real    ::   dp
      real    ::                   x7, x8
      real    ::       x8a, x5a, x5b, x5c, x5d, x6a, x6b, x6c, x6d
      real    ::      v5,     v6
      integer ::  k

!--------------------------------------------------------------------
!   local variables:
!
!         esum
!         esuma
!         esumc
!         sumf
!         summ
!         sumqme
!         sumg
!         sumn
!         sumelt
!         sumfre
!         summes
!         disl
!         disga
!         nlev            number of layers in large-scale model
!         k               do-loop index
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = '  ' ; error = 0
      lprcp = 0.
      vrt_mot = 0.
      liq_ice = 0.

!--------------------------------------------------------------------
!    define the total precipitation (total_precip_c) from the parameter-
!    ization as the sum of the convective (ensmbl_precip) and mesoscale
!    (meso_precip) precipitation. 
!--------------------------------------------------------------------
      total_precip_c = ensmbl_precip + meso_precip    

!----------------------------------------------------------------------
!    add the mesoscale cloud area to the cell-ensemble cloud area to
!    obtain the total cloud area profile.    
!----------------------------------------------------------------------
      do k=1,nlev_lsm
        cutotal (k) = ensmbl_cloud_area(k) + meso_cloud_area(k)
      end do

!---------------------------------------------------------------------
!    if in a diagnostics column, output the profiles of ensemble-total 
!    cloud area (ensmbl_cloud_area) and mesoscale cloud area 
!    (meso_cloud_area), total cloud area (cu_total), deposition in 
!    mesoscale updrafts (cmus), evaporation in mesoscale downdrafts 
!    (emds), evaporation from mesoscale updrafts (emes), water vapor 
!    supplied to mesoscale circulation (wmps), melted anvil precip-
!    itation ( anvil_precip_melt), mesoscale temperature flux
!    convergence (tmes) and mesoscale vapor-flux convergence (mrmes).
!---------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_lsm
          write (diag_unit, '(a, i4, f19.10, 3e20.12)') &
                 'in mulsub: jk, pr,cual,cuml, cutot= ', k, pfull_c(k),&
                  ensmbl_cloud_area (k), meso_cloud_area(k), cutotal(k)
          write (diag_unit, '(a, i4, 3e20.12)')  &
                     'in mulsub: jk,cmu,emd,eme= ', k, cmus_tot(k), &
                   emds_liq(k) + emds_ice(k), emes_liq(k) + emes_ice(k)
          write (diag_unit, '(a, i4, 2e20.12)') &
                     'in mulsub: jk,wmm,wmp,elt= ', k,           &
                      wmps(k),  -anvil_precip_melt(k)
          write (diag_unit, '(a, i4, f20.14, e20.12)')  &
                      'in mulsub: jk,tmes,qmes= ', k, tmes(k), mrmes(k)
        end do
      endif

!---------------------------------------------------------------------
!    define terms which will appear in the large-scale model equations.
!---------------------------------------------------------------------
      do k=1,nlev_lsm

!----------------------------------------------------------------------
!    combine several of the moisture tendency terms associated with the
!    donner_deep parameterization (disf). if a mesoscale circulation is
!    present, the following terms are included : 1) transfer of vapor 
!    from mesoscale to large-scale flow ( the sum of the water mass 
!    condensed in the mesoscale updraft plus the vapor transferred from
!    cell to mesoscale and then condensed -- cmus_tot), 2) evaporation 
!    in cumulus downdrafts (ecds), 3) evaporation from cumulus updrafts 
!    (eces), 4)  vapor transferred from cells to mesoscale (wmps), 5) 
!    evaporation from mesoscale updrafts (emes), 6) evaporation from 
!    mesoscale downdrafts (emds), and 7) mesoscale moisture-flux 
!    convergence (mrmes). 
!----------------------------------------------------------------------
        if (lmeso) then
          cell_evap(k) = (ecds_liq(k) + ecds_ice(k)) +  &
                         (eces_liq(k) + eces_ice(k))
          meso_cd(k) =  wmms(k)
          meso_depo(k) = -cmus_tot(k) - wmms(k)
          meso_evap(k) = emes_liq(k) + emds_liq(k) + &
                         emes_ice(k) + emds_ice(k)

!----------------------------------------------------------------------
!    if a mesoscale circulation is not present, disf is simply the 
!    moisture tendency associated with the evaporation of the condensed
!    cloud water that did not precipitate out (enev). convert to units 
!    of g(h2o) per kg(air) per day.
!----------------------------------------------------------------------
        else
          cell_evap(k) = enev(k)*(1.0E03*Param%seconds_per_day)
          meso_cd(k) = 0.
          meso_depo(k) = 0.
          meso_evap(k) = 0.
        endif

        disf(k) = meso_depo(k) + meso_cd(k) + cell_evap(k) +    &
                  wmps(k) + meso_evap(k)  + mrmes_up(k) + mrmes_dn(k)

!---------------------------------------------------------------------
!    define the sum of disf and the term containing the tendency due 
!    to cell-scale vertical moisture-flux convergence and associated
!    condensation (encmf), and store in array dise.
!---------------------------------------------------------------------
        if (Nml%do_donner_lscloud) then
          dise(k) = encmf(k) + disf(k)
        else
          dise(k) = encmf(k)
        endif

!----------------------------------------------------------------------
!    define the temperature tendencies associated with the freezing
!    of updraft liquid (temp_tend_freeze) and the melting of ice
!    falling from the anvil (temp_tend_melt). combine several of the 
!    temperature tendencies associated with the cell component of the
!    donner_deep parameterization (disn). disn is composed of 1) a term
!    combining the vertical flux convergence of temperature and cloud 
!    condensation (enctf), 2) evaporation of liquid in the cell updraft
!    and downdraft (disg), 3) the freezing of updraft liquid 
!    (temp_tend_freeze) and 4) the melting of ice (temp_tend_melt). 
!    separately define the temperature tendency resulting from the 
!    latent heat release associated with sublimation occurring in the 
!    mesoscale updraft and downdraft (disga).
!----------------------------------------------------------------------
        ensmbl_freeze (k) = ensmbl_freeze(k)*Param%hlf/     &
                               (Param%cp_air*1000.)
        ensmbl_freeze_meso (k) = ensmbl_freeze_meso(k)*Param%hlf/  &
                               (Param%cp_air*1000.)
        anvil_precip_melt(k) =  (anvil_precip_melt(k))*  &
                                          Param%hlf/(Param%cp_air*1000.)
        ensmbl_melt(k) =  ensmbl_melt(k)*  &
                                          Param%hlf/(Param%cp_air*1000.)
        ensmbl_melt_meso(k) =  ensmbl_melt_meso(k)* &
                                          Param%hlf/(Param%cp_air*1000.)
        if (lmeso) then
          disga_ice_up(k) =  -(emes_ice(k) )*  &
                                    Param%hls/(Param%cp_air*1000.)
          disga_ice_dn(k) =  -( emds_ice(k))*  &
                                    Param%hls/(Param%cp_air*1000.)
          if (melting_in_cloud) then
            disl_ice_depo(k) =     &
                           -meso_depo(k)*Param%hls/(Param%cp_air*1000.)
            disl_ice_cd(k) =       &
                          -( meso_cd(k))*Param%hls/(Param%cp_air*1000.)
            disl_liq_depo(k) = 0.
            disl_liq_cd(k) = 0.
            disga_liq_up(k) =  -(emes_liq(k) )*  &
                                    Param%hls/(Param%cp_air*1000.)
            disga_liq_dn(k) =  -(emds_liq(k))*  &
                                    Param%hls/(Param%cp_air*1000.)
          else
            disga_liq_up(k) =  -emes_liq(k) * &
                                    Param%hlv/(Param%cp_air*1000.)
            disga_liq_dn(k) =  -emds_liq(k) * &
                                     Param%hlv/(Param%cp_air*1000.)
            if (.not. meso_frz_intg_sum      ) then
              disl_liq_depo(k) = -(meso_depo(k))*    &
                                     Param%hlv/(Param%cp_air*1000.)
              disl_liq_cd(k) = -( meso_cd(k))*    &
                                     Param%hlv/(Param%cp_air*1000.)
              disl_ice_depo(k) = 0.
              disl_ice_cd(k) = 0.
            else
! if no melting but freezing, then hls carried out
              disl_ice_depo(k) = -meso_depo(k)*  &
                                      Param%hls/(Param%cp_air*1000.)
              disl_ice_cd(k) = -( meso_cd(k))*     &
                                      Param%hls/(Param%cp_air*1000.)
              disl_liq_depo(k) = 0.
              disl_liq_cd(k) = 0.
            endif
          endif
        else ! (lmeso)
          disl_liq_depo(k) = 0.
          disl_liq_cd(k) = 0.
          disl_ice_depo(k) = 0.
          disl_ice_cd(k) = 0.
          disga_liq_up(k) = 0.
          disga_liq_dn(k) = 0.
          disga_ice_up(k) = 0.
          disga_ice_dn(k) = 0.
        endif

!--------------------------------------------------------------------
!    if in a diagnostics column, output the profile of temperature 
!    change associated with evaporation in the mesoscale circulation 
!    (disga).
!--------------------------------------------------------------------
        if (debug_ijt) then
          if ( (disga_liq_up(k) + disga_liq_dn(k) +   &
                disga_ice_up(k) + disga_ice_dn(k)) /= 0.0) then
            write (diag_unit, '(a, i4, f19.10,  e20.12)')  &
                    'in mulsub: jk,pr,disga= ', k, pfull_c(k),  &
                      -(disga_liq_up(k) + disga_liq_dn(k) +   &
                      disga_ice_up(k) + disga_ice_dn(k))
            write (diag_unit, '(a, i4, f19.10,  2e20.12)')  &
                    'in mulsub: jk,pr,disgal, disgai= ', k,   &
            pfull_c(k), -(disga_liq_up(k) + disga_liq_dn(k)),  &
                    -(disga_ice_up(k) + disga_ice_dn(k))
          endif
        endif

!--------------------------------------------------------------------
!    if a mesoscale circulation is present, define the heating terms
!    equivalent to the disf array (disg_2). included in disg_2 are 
!    terms associated with 1) transfer of vapor from mesoscale to 
!    large-scale flow (disl), 2) evaporation in cumulus updrafts and 
!    downdrafts (disg), 3) freezing of liquid in the updraft 
!    (temp_tend_freeze), 4) melting of ice (temp_tend_melt), 5) evap-
!    oration in the mesoscale circulation (disga), and 6) mesoscale 
!    temperature flux convergence (tmes).
!--------------------------------------------------------------------
        if (lmeso) then
          disg_2liq(k) = disl_liq_depo(k)  + disl_liq_cd(k) +    &
                         disg_liq(k) + disga_liq_up(k) + disga_liq_dn(k)
          disg_2ice(k) = disl_ice_depo(k) + disl_ice_cd(k) +   &
                         disg_ice(k) + disga_ice_up(k) +disga_ice_dn(k)

!---------------------------------------------------------------------
!    define the sum of disg_2 and the term containing the tendency due 
!    to cell-scale vertical temperature-flux convergence and associated
!    condensation (enctf), and store in array disa.
!---------------------------------------------------------------------
          disn(k) = enctf(k) + disg_liq(k) + disg_ice(k) +  &
                     ensmbl_freeze_meso(k) + &
                    ensmbl_freeze(k) + ensmbl_melt(k)  + &
                    ensmbl_melt_meso(k)                         
          if (Nml%do_donner_lscloud) then
            disa(k) = enctf(k) + (disl_liq_depo(k)  + disl_liq_cd(k) + &
                    disl_ice_depo(k) + disl_ice_cd(k)) +   &
                    disg_liq(k) + disga_liq_up(k) + disga_liq_dn(k) + &
                    disg_ice(k) + disga_ice_up(k) +disga_ice_dn(k) + &
                    ensmbl_freeze_meso(k) + &
                    ensmbl_freeze(k) + ensmbl_melt(k) +  &
                    ensmbl_melt_meso(k) + &
                    anvil_precip_melt(k) + tmes_up(k) + tmes_dn(k)
          else
            disa(k) = enctf(k)
          endif
        else ! (lmeso)
          disg_2liq(k) = -disze2(k)         
          disg_2ice(k) = -disze1(k)

!---------------------------------------------------------------------
!    define the sum of disg_2 and the term containing the tendency due 
!    to cell-scale vertical temperature-flux convergence and associated
!    condensation (enctf), and store in array disa.
!---------------------------------------------------------------------
          if (Nml%do_donner_lscloud) then
            disa(k) = enctf(k) + disg_2liq(k) + disg_2ice(k) + &
                     disze3(k) + ensmbl_freeze(k) + ensmbl_melt(k) 
          else
            disa(k) = enctf(k)
          endif

          disn(k) = disa(k)
        endif

        temp_tend_melt(k) = ensmbl_melt(k) +   &
                              ensmbl_melt_meso(k) + &
                                               anvil_precip_melt(k)

!--------------------------------------------------------------------
!    define precip types.
!--------------------------------------------------------------------
        if (lmeso) then
          if (melting_in_cloud) then 
            disl_ice_melted(k) = -(meso_depo(k) + meso_cd(k))/(1000.)
            disl_liq(k) = 0.
            disl_ice(k) = 0.
          else
            if (.not. meso_frz_intg_sum      ) then
              disl_liq(k) = -(meso_depo(k) + meso_cd(k))/(1000.)
              disl_ice(k) = 0.
              disl_ice_melted(k) = 0.
            else
! if no melting but freezing, then hls carried out
              disl_ice    (k) = -(meso_depo(k) + meso_cd(k))/(1000.)
              disl_ice_melted(k) = 0.
              disl_liq(k) = 0.
            endif
          endif
        else ! (lmeso)
          disl_liq(k) = 0.
          disl_ice(k) = 0.
          disl_ice_melted(k) = 0.
        endif
        if (lmeso) then
            liq_prcp(k) = Param%anvil_precip_efficiency*(dism_liq(k) + &
                          disl_liq(k) +    &
                          dism_liq_remelt(k) + dism_ice_melted(k) + &
                          disl_ice_melted(k)) + disp_liq(k) + &
                          disp_melted(k) + disz_remelt(k)
            frz_prcp(k) = Param%anvil_precip_efficiency*  &
                          (dism_liq_frz(k) + dism_ice(k) +   &
                          disl_ice(k)) + disp_ice(k) + disz(k)
          else
            liq_prcp(k) = disp_liq(k) + disp_melted(k) + disz_remelt(k)
            frz_prcp(k) = disp_ice(k) + disz(k)
          endif
     end do

      if (Nml%do_budget_analysis .or.    &
                             Initialized%do_conservation_checks) then
        tsumb = 0.
        tsumiup = 0.
        tsumidn = 0.
        do k=1,nlev_lsm
          dp = Param%cp_air*(phalf_c(k) - phalf_c(k+1))/Param%grav
          tsumb = tsumb + disb(k)*dp
          tsumiup  = tsumiup  + tmes_up(k)*dp
          tsumidn  = tsumidn  + tmes_dn(k)*dp
        end do

        if (lmeso) then
          vrt_mot = tsumiup + tsumidn + tsumb
        else
          vrt_mot = tsumb
        endif

        x5a = 0.
        x5b = 0.
        x5c = 0.
        x5d = 0.
        x6a = 0.
        x6b = 0.
        x6c = 0.
        x6d = 0.
        x7 = 0.
        x8 = 0.
        x8a = 0.
        v5 = 0.
        v6 = 0.

        do k=1,nlev_lsm
          dp = (phalf_c(k) - phalf_c(k+1))/Param%grav
          v5 = v5 + disz (k)*dp
          v6 = v6 + disz_remelt (k)*dp
          if (lmeso) then
          x5a = x5a + (dism_liq(k)                     )*dp
          x5b = x5b + (disl_liq(k) )*dp
          x5c = x5c + (dism_liq_frz(k) )*dp
          x5d = x5d + (dism_liq_remelt(k) )*dp
          x6a = x6a + (dism_ice(k)                     )*dp
          x6b = x6b + (disl_ice(k))*dp
          x6c = x6c + (dism_ice_melted(k) )*dp
          x6d = x6d + (disl_ice_melted(k))*dp
          endif
          x7 = x7 + disp_liq(k)*dp
          x8 = x8 + disp_ice(k)*dp
          x8a = x8a + disp_melted(k)*dp
        end do
        x5a = x5a*Param%anvil_precip_efficiency
        x5b = x5b*Param%anvil_precip_efficiency
        x5c = x5c*Param%anvil_precip_efficiency
        x5d = x5d*Param%anvil_precip_efficiency
        x6a = x6a*Param%anvil_precip_efficiency
        x6b = x6b*Param%anvil_precip_efficiency
        x6c = x6c*Param%anvil_precip_efficiency
        x6d = x6d*Param%anvil_precip_efficiency

      if (debug_ijt) then
         
      do k=1, nlev_lsm 
         write (diag_unit, '(1(a, i4, 2e20.12))') &
          
             'total precip-- k, liq, frz', k, liq_prcp(k), frz_prcp(k) 
      end do
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  &
            'PRECIPITATION SOURCES  --  UNITS OF  mm / day'
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a, e20.12))')  &
           'TOTAL PRECIPITATION: ', x5a + x5b + x5c + x5d + x6a + &
                                     x6b + x6c + x6d + x7 +  &
                                     x8 + v5 + x8a + v6
         write (diag_unit, '(1(a, e20.12))')  &
           'MESO  PRECIPITATION: ', x5a + x5b + x5c + x5d + x6a + &
                                     x6b + x6c + x6d 
         write (diag_unit, '(1(a, e20.12))')  &
           'CELL  PRECIPITATION: ', x7 + x8 + x8a + v5 + v6
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a, e20.12))')  &
           'cell liquid condensate precipitated out as liquid', x7 
         write (diag_unit, '(1(a, e20.12))')  &
          'cell liquid condensate precipitated out as frozen liquid', v5
         write (diag_unit, '(1(a, e20.12))')  &
           'cell liquid condensate which froze, remelted and  &
               &precipitated out as liquid', v6 
         write (diag_unit, '(1(a, e20.12))')  &
           'cell liquid condensate transferred to the mesoscale &
                &circulation and precipitated out as liquid', x5a
         write (diag_unit, '(1(a, e20.12))')  &
           'cell liquid condensate transferred to the mesoscale &
                  &circulation, then frozen and precipitated out &
                  &as frozen liquid', x5c
         write (diag_unit, '(1(a, e20.12))')  &
           'cell liquid condensate transferred to the mesoscale &
              &circulation, frozen, remelted and precipitated &
                                       &out as liquid', x5d
         write (diag_unit, '(1(a, e20.12))')  &
           'mesoscale liquid condensate precipitated out as liquid', x5b
         write (diag_unit, '(1(a, e20.12))')  &
           'cell ice    condensate precipitated out as ice   ', x8 
         write (diag_unit, '(1(a, e20.12))')  &
           'cell ice    condensate which melted and  precipitated &
                             &out as liquid', x8a
         write (diag_unit, '(1(a, e20.12))')  &
           'cell ice transferred to mesoscale and precipitated &
                        &out as ice   ', x6a
         write (diag_unit, '(1(a, e20.12))')  &
           'cell ice transferred to mesoscale which melted and &
                      & precipitated out as liquid', x6c
         write (diag_unit, '(1(a, e20.12))')  &
           'mesoscale ice condensate  precipitated out as ice   ', x6b
         write (diag_unit, '(1(a, e20.12))')  &
           'mesoscale ice condensate which melted and  &
                    &precipitated out as liquid', x6d
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
       
   endif

        v5 = v5*Param%hls   
        v6 = v6*Param%hlv    
        x5a = x5a*Param%hlv
        x5b = x5b*Param%hlv
        x5c = x5c*Param%hls
        x5d = x5d*Param%hlv
        x6a = x6a*Param%hls
        x6b = x6b*Param%hls
        x6c = x6c*Param%hlv
        x6d = x6d*Param%hlv
        x7 = x7*Param%hlv   
        x8 = x8*Param%hls   
        x8a = x8a*Param%hlv      
        lprcp = x5a + x5b + x5c + x5d + x6a + x6b + x6c + x6d + x7 +  &
                x8 + v5 + x8a + v6

      if (debug_ijt) then
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  &
            'LATENT HEAT REMOVED BY THE PRECIPITATION SOURCES  --&
                           &  UNITS OF  J / (m**2 day)'
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a, e20.12))')  &
           'TOTAL PRECIPITATION: ', x5a + x5b + x5c + x5d + x6a + &
                                     x6b + x6c + x6d + x7 +  &
                                     x8 + v5 + x8a + v6
         write (diag_unit, '(1(a, e20.12))')  &
           'MESO  PRECIPITATION: ', x5a + x5b + x5c + x5d + x6a + &
                                     x6b + x6c + x6d 
         write (diag_unit, '(1(a, e20.12))')  &
           'CELL  PRECIPITATION: ', x7 + x8 + x8a + v5 + v6
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate precipitated &
                             &out as liquid', x7 
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate precipitated &
                             &out as frozen liquid', v5 
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate which froze, &
                       &remelted and precipitated out as liquid', v6 
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate transferred to &
                &the mesoscale circulation and precipitated out&
                & as liquid', x5a
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate transferred to &
             &the mesoscale circulation, then frozen and precipitated &
               &out as frozen liquid', x5c
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell liquid condensate transferred &
               &to the mesoscale circulation, frozen, remelted &
                 &and precipitated out as liquid', x5d
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by mesoscale liquid condensate &
                          &precipitated out as liquid', x5b
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell ice condensate precipitated &
                     &out as ice   ', x8 
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell ice    condensate which melted &
                  &and  precipitated out as liquid', x8a
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell ice transferred to mesoscale and &
                      &precipitated out as ice   ', x6a
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by cell ice transferred to mesoscale &
                &which melted and  precipitated out as liquid', x6c
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by mesoscale ice condensate  precipitated &
                       &out as ice   ', x6b
         write (diag_unit, '(1(a, e20.12))')  &
           'heat removed by mesoscale ice condensate which melted &
                        &and  precipitated out as liquid', x6d
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
       
   endif
      endif

      if (Initialized%do_conservation_checks .or.   &
                                         Nml%do_budget_analysis) then
! units of water budget: g(h20) / kg(air) / day
        do k=1,nlev_lsm
          water_budget(k,1) = dise(nlev_lsm-k+1)
          water_budget(k,2) = encmf(nlev_lsm-k+1)
          water_budget(k,3) = meso_depo(nlev_lsm-k+1)
          water_budget(k,4) = meso_cd(nlev_lsm-k+1)
          water_budget(k,5) = cell_evap(nlev_lsm-k+1)
          water_budget(k,6) = wmps (nlev_lsm-k+1)
          water_budget(k,7) = meso_evap(nlev_lsm-k+1)
          water_budget(k,8) = mrmes_up(nlev_lsm-k+1)
          water_budget(k,9) = mrmes_dn(nlev_lsm-k+1)

!  enthalpy_budget terms are in units of deg K / day.
          if (lmeso) then
            enthalpy_budget(k,1) = disa(nlev_lsm-k+1)
            enthalpy_budget(k,2) = enctf(nlev_lsm-k+1)
            enthalpy_budget(k,3) = disl_liq_depo(nlev_lsm-k+1)
            enthalpy_budget(k,4) = disl_liq_cd(nlev_lsm-k+1)
            enthalpy_budget(k,5) = disg_liq(nlev_lsm-k+1)
            enthalpy_budget(k,6) = disga_liq_up(nlev_lsm-k+1)
            enthalpy_budget(k,7) = disga_liq_dn(nlev_lsm-k+1)
            enthalpy_budget(k,8) = disl_ice_depo(nlev_lsm-k+1)
            enthalpy_budget(k,9) = disl_ice_cd  (nlev_lsm-k+1)
            enthalpy_budget(k,10) = disg_ice(nlev_lsm-k+1)
            enthalpy_budget(k,11) = disga_ice_up(nlev_lsm-k+1)
            enthalpy_budget(k,12) = disga_ice_dn(nlev_lsm-k+1)
            enthalpy_budget(k,13) = ensmbl_freeze_meso(nlev_lsm-k+1)
            enthalpy_budget(k,14) = ensmbl_freeze(nlev_lsm-k+1)
            enthalpy_budget(k,15) = ensmbl_melt(nlev_lsm-k+1)
            enthalpy_budget(k,16) = ensmbl_melt_meso(nlev_lsm-k+1)
            enthalpy_budget(k,17) = anvil_precip_melt(nlev_lsm-k+1)
            enthalpy_budget(k,18) = tmes_up(nlev_lsm-k+1)
            enthalpy_budget(k,19) = tmes_dn(nlev_lsm-k+1)
          else
            enthalpy_budget(k,1) = disa(nlev_lsm-k+1)
            enthalpy_budget(k,2) = enctf(nlev_lsm-k+1)
            enthalpy_budget(k,5) = disg_2liq(nlev_lsm-k+1) +  &
                                               disg_2ice(nlev_lsm-k+1)
            enthalpy_budget(k,10) = disze3   (nlev_lsm-k+1)
            enthalpy_budget(k,14) = ensmbl_freeze(nlev_lsm-k+1)
            enthalpy_budget(k,15) = ensmbl_melt(nlev_lsm-k+1)
          endif
        end do

!--------------------------------------------------------------------
!    compute the column integrals of the various tendency terms for 
!    the vapor equation. 
!   vsuma  : total vapor tendency from donner_deep parameterization
!    sumf  : total vapor tendency less the vertical flux convergence and
!            condensation      
!    summ  : vapor tendency due to vertical flux convergence and
!            condensation 
!    sumqme: mesoscale moisture flux convergence
!--------------------------------------------------------------------

        sumf   = 0.
        summ   = 0.
        vsuma  = 0.
        vsumb  = 0.
        vsumc  = 0.
        vsumd  = 0.
        vsumd1 = 0.
        vsumd2 = 0.
        vsume  = 0.
        vsumf  = 0.
        vsumg  = 0.
        vsumg1 = 0.
        vsumg2 = 0.
        vsumh  = 0.
        vsumi  = 0.
        vsumi1 = 0.
        vsumi2 = 0.
         
        do k=1,nlev_lsm
          dp = (phalf_c(k) - phalf_c(k+1))/Param%grav
          sumf   = sumf   + disf(k)*dp
          summ   = summ   + encmf(k)*dp
          vsuma = vsuma + dise(k)*dp
          vsumb = vsumb + disd(k)*dp
          vsumc = vsumc - disv(k)*dp
          vsumd = vsumd + cell_evap(k)*dp
          vsumd1 = vsumd1 + (ecds_liq(k) + ecds_ice(k))*dp
          vsumd2 = vsumd2 + (eces_liq(k) + eces_ice(k))*dp
          vsume  = vsume + meso_cd(k)*dp
          vsumf = vsumf + meso_depo(k)*dp
          vsumg = vsumg + meso_evap(k)*dp
          vsumg1 = vsumg1 + (emes_liq(k) + emes_ice(k))*dp
          vsumg2 = vsumg2 + (emds_liq(k) + emds_ice(k))*dp
          vsumh = vsumh + wmps(k)*dp
          vsumi1 = vsumi1 + mrmes_up(k)*dp
          vsumi2 = vsumi2 + mrmes_dn(k)*dp
        end do
!---------------------------------------------------------------------
!    convert the moisture terms to units of mm(h2o) per day.
!---------------------------------------------------------------------
        sumf   = sumf/(1000.)
        summ   = summ/(1000.)
        vsuma   = vsuma/(1000.)
        vsumb   = vsumb/(1000.)
        vsumc   = vsumc/(1000.)
        vsumd   = vsumd/(1000.)
        vsumd1  = vsumd1/(1000.)
        vsumd2  = vsumd2/(1000.)
        vsume   = vsume/(1000.)
        vsumf   = vsumf/(1000.)
        vsumg   = vsumg/(1000.)
        vsumg1  = vsumg1/(1000.)
        vsumg2  = vsumg2/(1000.)
        vsumh   = vsumh/(1000.)
        vsumi1  = vsumi1/(1000.)
        vsumi2  = vsumi2/(1000.)
        vsumi = vsumi1 + vsumi2


! units for precip_budget: (kg(h2o) / kg (air) / day
       do k=1,nlev_lsm
         precip_budget(k,1,1) = disp_liq(nlev_lsm-k+1)
         precip_budget(k,2,1) = disz(nlev_lsm-k+1)
         precip_budget(k,3,1) = disz_remelt(nlev_lsm-k+1)
         precip_budget(k,4,1) = disp_ice(nlev_lsm-k+1)           
         precip_budget(k,5,1) = disp_melted(nlev_lsm-k+1)            

         precip_budget(k,1,2) = dism_liq(nlev_lsm-k+1)
         precip_budget(k,2,2) = dism_liq_frz(nlev_lsm-k+1)
         precip_budget(k,3,2) = dism_liq_remelt(nlev_lsm-k+1)
         precip_budget(k,4,2) = dism_ice(nlev_lsm-k+1)
         precip_budget(k,5,2) = dism_ice_melted(nlev_lsm-k+1)

         precip_budget(k,1,3) = disl_liq(nlev_lsm-k+1)
         precip_budget(k,2,3) = 0.0  
         precip_budget(k,3,3) = 0.0
         precip_budget(k,4,3) = disl_ice(nlev_lsm-k+1)
         precip_budget(k,5,3) = disl_ice_melted(nlev_lsm-k+1)
       end do
         precip_budget(:,:,2) = precip_budget(:,:,2)*   &
                                Param%anvil_precip_efficiency
         precip_budget(:,:,3) = precip_budget(:,:,3)*   &
                                Param%anvil_precip_efficiency
!--------------------------------------------------------------------
!    output the various column integrals.
!--------------------------------------------------------------------
       if (debug_ijt) then
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(a, e20.12, a)') &
              'in mulsub: CELL MOISTURE FORCING= ', summ, ' MM/DAY'
         write (diag_unit, '(a, e20.12, a)') &
          'in mulsub: TOTAL TENDENCY LESS CELL MOISTURE FORCING= ', &
                                                       sumf, ' MM/DAY'
         write (diag_unit, '(a, e20.12, a)') &
          'in mulsub: TOTAL CELL MOISTURE TENDENCY= ', &
                             summ + vsumd, ' MM/DAY'
         if (lmeso) then
           write (diag_unit, '(a, e20.12, a)') &
            'in mulsub: TOTAL MESO MOISTURE TENDENCY= ', &
                             vsuma - summ - vsumd, ' MM/DAY'
         endif
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  &
            'COLUMN INTEGRAL VAPOR BUDGET TERMS --  &
                                           &UNITS OF  mm / day'
         write (diag_unit, '(1(a))')   '***************************'
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(a, e20.12)') &
             'TOTAL VAPOR TENDENCY=', vsuma 
         write (diag_unit, '(a, e20.12)') &
              '    DYNAMICAL TENDENCY=', vsumb + vsumh + vsumi
         write (diag_unit, '(a, e20.12)') &
             '         VAPOR CONVERGENCE IN CELLS=', vsumb 
         if (lmeso) then
           write (diag_unit, '(a, e20.12)') &
             '         TRANSFER FROM CELL UPDRAFTS TO MESO=', vsumh 
           write (diag_unit, '(a, e20.12)') &
             '         TRANSFER BY MESOSCALE EDDY FLUXES= ', vsumi
           write (diag_unit, '(a, e20.12)') &
             '              TRANSFER BY UPWARD EDDIES= ', vsumi1
           write (diag_unit, '(a, e20.12)') &
             '              TRANSFER BY DOWNWARD EDDIES= ', vsumi2
         endif
         write (diag_unit, '(a, e20.12)') &
             '    CONDENSATION IN CELLS=', vsumc 
         write (diag_unit, '(a, e20.12)') &
             '    EVAPORATION IN CELLS= ', vsumd

         if (lmeso) then
           write (diag_unit, '(a, e20.12)') &
             '         CELL DOWNDRAFT EVAP=', vsumd1  
           write (diag_unit, '(a, e20.12)') &
             '         CELL UPDRAFT EVAP=', vsumd2  
           write (diag_unit, '(a, e20.12)') &
             '    MESOSCALE CONDENSATION  =', vsume
           write (diag_unit, '(a, e20.12)') &
             '    MESOSCALE DEPOSITION=', vsumf 
           write (diag_unit, '(a, e20.12)') &
             '    MESOSCALE EVAPORATION=', vsumg  
           write (diag_unit, '(a, e20.12)') &
             '         MESO EVAPORATION IN UPDRAFTS =', vsumg1
           write (diag_unit, '(a, e20.12)') &
             '         MESO EVAPORATION IN DOWNDRAFTS=', vsumg2  
         endif

         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a))')  
         write (diag_unit, '(1(a,3l4         ))')  &
                  'in mulsub:  lmeso, anvil melts?, meso freezes?', &
                  lmeso, melting_in_cloud, meso_frz_intg_sum  


       endif ! (debug_ijt)

!--------------------------------------------------------------------
!    compute the column integrals of the various tendency terms for 
!    the temperature equation. 
!    tsuma : total temperature tendency from donner_deep parameter-
!            ization
!    sumg  : total temperature tendency less the vertical flux conver-
!            gence and condensation
!    esumc : temperature tendency due to vertical flux convergence and
!            condensation  
!    summes: temperature tendency due to mesoscale temperature flux
!            convergence
!    sumelt: temperature tendency due to melting within column
!    sumfre: temperature tendency due to freezing within the column
!    sumn  : temperature tendency associated with the cell component of
!            the donner_deep parameterization
!--------------------------------------------------------------------

       esumc  = 0.
       sumn = 0.
       tsumj = 0.
       tsumj1 = 0.
       tsumk1 = 0.
       tsumk2 = 0.
       tsumk3 = 0.
       tsumcliq = 0.
       tsumdliq = 0.
       tsumeliq = 0.
       tsumfliq = 0.
       tsumg1liq = 0.
       tsumg2liq = 0.
       tsumcice = 0.
       tsumdice = 0.
       tsumeice = 0.
       tsumfice = 0.
       tsumg1ice = 0.
       tsumg2ice = 0.
       tsummliq = 0.
       tsummfliq = 0.
       tsummice = 0.

       do k=1,nlev_lsm
         dp = Param%cp_air*(phalf_c(k) - phalf_c(k+1))/Param%grav
         esumc  = esumc  + enctf(k)*dp        
         sumn   = sumn   + disn(k)*dp

         tsumcliq = tsumcliq + disc_liq(k)*dp
         tsumdliq = tsumdliq + disg_liq(k)*dp
         tsumeliq = tsumeliq + disl_liq_cd(k)*dp
         tsumfliq = tsumfliq + disl_liq_depo(k)*dp
         tsumg1liq = tsumg1liq + disga_liq_up(k)*dp
         tsumg2liq = tsumg2liq + disga_liq_dn(k)*dp
         tsumcice = tsumcice + disc_ice(k)*dp
         tsumdice = tsumdice + disg_ice(k)*dp
         tsumeice = tsumeice + disl_ice_cd(k)*dp
         tsumfice = tsumfice + disl_ice_depo(k)*dp
         tsumg1ice = tsumg1ice + disga_ice_up(k)*dp
         tsumg2ice = tsumg2ice + disga_ice_dn(k)*dp
         if (.not. lmeso) then
           tsummliq = tsummliq + disg_2liq(k)*dp
           tsummfliq = tsummfliq + disg_2ice(k)*dp
           tsummice = tsummice +  disze3(k)*dp
         endif
          
         tsumj    = tsumj    + ensmbl_freeze(k)*dp
         tsumj1   = tsumj1   + ensmbl_freeze_meso(k)*dp
         tsumk1   = tsumk1   + ensmbl_melt(k)*dp
         tsumk3   = tsumk3   + ensmbl_melt_meso(k)*dp
         tsumk2   = tsumk2 +  anvil_precip_melt(k) *dp
       end do

       if (lmeso) then
         tsumaliq = tsumcliq + tsumdliq + tsumeliq + tsumfliq +  &
                    tsumg1liq + tsumg2liq               
         tsumaice = tsumcice + tsumeice + tsumfice + tsumdice +  &
                    tsumg1ice + tsumg2ice
         liq_ice = tsumj + tsumj1 + tsumk1 + tsumk2 + tsumk3
       else
         tsumaliq = tsumcliq + tsummliq 
         tsumaice = tsumcice + tsummice  + tsummfliq
         liq_ice = tsumj + tsumk1
       endif

       tsuma = tsumaliq + tsumaice + vrt_mot + liq_ice     


       if (debug_ijt) then
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: CELL TEMPERATURE FORCING= ', esumc,   &
                                               ' Joules / (m**2 * DAY)'
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: CELL TENDENCY LESS FORCING ', sumn - esumc,   &
                                               ' Joules / (m**2 * DAY)'
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: TOTAL CELL TEMPERATURE TENDENCY ', sumn,   &
                                               ' Joules / (m**2 * DAY)'
        if (lmeso) then
        write (diag_unit, '(a, e20.12, a)') &
             'in mulsub: TOTAL MESO TEMPERATURE TENDENCY ', &
                                         tsuma - sumn,  &
                                               ' Joules / (m**2 * DAY)'
        endif
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(1(a))')   '***************************'
        write (diag_unit, '(1(a))')  &
            'COLUMN INTEGRAL TEMPERATURE BUDGET TERMS --  &
                                &UNITS OF (Joules / (m**2 * day)'
        write (diag_unit, '(1(a))')   '***************************'
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(1(a,e20.12))')  &
            'TOTAL TEMPERATURE TENDENCY=',tsuma

        if (lmeso) then
        write (diag_unit, '(1(a,e20.12))')  &
            '    DYNAMICAL TENDENCY=', vrt_mot
        write (diag_unit, '(1(a,e20.12))')  &
            '       CELL ENTROPY CONVERGENCE=', tsumb 
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESO ENTROPY CONVERGENCE=', tsumiup + tsumidn
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO UP ENTROPY CONVERGENCE=', tsumiup
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO DOWN ENTROPY CONVERGENCE=', tsumidn
        write (diag_unit, '(1(a,e20.12))')  &
            '    LATENT HEATING: VAPOR/LIQUID=', tsumaliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       CONDENSATION IN CELLS=', tsumcliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       EVAPORATION IN CELLS=', tsumdliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESOSCALE CONDENSATION=', tsumeliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESOSCALE DEPOSITION=', tsumfliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESO EVAPORATION=', tsumg1liq + tsumg2liq
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO EVAPORATION UPDRAFT=', tsumg1liq
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO EVAPORATION DOWNDRAFT=', tsumg2liq
        write (diag_unit, '(1(a,e20.12))')  &
            '    LATENT HEATING: VAPOR/ICE=',tsumaice
        write (diag_unit, '(1(a,e20.12))')  &
            '       CELL CONDENSATION=', tsumcice
        write (diag_unit, '(1(a,e20.12))')  &
            '       CELL UPDRAFT   EVAPORATION=', tsumdice
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESOSCALE CONDENSATION=', tsumeice
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESOSCALE DEPOSITION=', tsumfice
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESOSCALE EVAPORATION=', tsumg1ice + tsumg2ice
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO EVAPORATION IN UPDRAFTS=', tsumg1ice
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO EVAPORATION IN DOWNDRAFTS=', tsumg2ice
        write (diag_unit, '(1(a,e20.12,a))')  &
            '    LATENT HEATING: LIQUID/ICE=', liq_ice
        write (diag_unit, '(1(a,e20.12))')  &
            '       CELL FREEZING',   tsumj      
        write (diag_unit, '(1(a,e20.12))')  &
            '       MESO FREEZING',   tsumj1      
        write (diag_unit, '(1(a,e20.12))')  &
            '       TOTAL MELTING',   tsumk1 + tsumk2 + tsumk3      
        write (diag_unit, '(1(a,e20.12))')  &
            '             CELL MELTING',   tsumk1     
        write (diag_unit, '(1(a,e20.12))')  &
            '             MESO MELTING (FOR CONSRV OF ICE)',   tsumk3 
        write (diag_unit, '(1(a,e20.12))')  &
            '             ANVIL PRECIP MELTING',   tsumk2      
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(1(a))')  



        else ! (lmeso)
        write (diag_unit, '(1(a,e20.12))')  &
            '    DYNAMICAL TENDENCY=', vrt_mot
        write (diag_unit, '(1(a,e20.12))')  &
            '       CONVERGENCE FROM CELL-SCALE MOTIONS=', tsumb 
        write (diag_unit, '(1(a,e20.12))')  &
            '    LATENT HEATING: VAPOR/LIQUID=', tsumaliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       LIQUID CONDENSATION IN CELLS=', tsumcliq
        write (diag_unit, '(1(a,e20.12))')  &
            '       LIQUID EVAPORATION IN CELLS=', tsummliq
        write (diag_unit, '(1(a,e20.12))')  &
            '    LATENT HEATING: VAPOR/ICE=',tsumaice
        write (diag_unit, '(1(a,e20.12))')  &
            '       ICE CONDENSATION IN CELLS=', tsumcice
        write (diag_unit, '(1(a,e20.12))')  &
            '       ICE EVAPORATION IN CELLS=', tsummice
        write (diag_unit, '(1(a,e20.12))')  &
            '       FROZEN LIQUID EVAPORATION IN CELLS=', tsummfliq
        write (diag_unit, '(1(a,e20.12,a))')  &
            '    LATENT HEATING: LIQUID/ICE=', liq_ice
        write (diag_unit, '(1(a,e20.12))')  &
            '       FREEZING IN CELLS',   tsumj      
        write (diag_unit, '(1(a,e20.12))')  &
            '       MELTING IN CELLS',   tsumk1     
        write (diag_unit, '(1(a))')  
        write (diag_unit, '(1(a))')  
        endif   ! (lmeso)

      endif ! (debug_ijt)

   endif ! (do_budget_analysis)

!---------------------------------------------------------------------
!    call subroutine output_diagnostic_profiles to print various 
!    output fields from the donner_deep parameterization in those 
!    columns for which diagnostics have been requested.
!---------------------------------------------------------------------
      if (debug_ijt) then
         call don_d_output_diag_profs_k    &
              (nlev_lsm, diag_unit, pfull_c,  disc_liq, disc_ice,  &
               disb, disd, disn,  &
              encmf, ensmbl_freeze, ensmbl_freeze_meso, &
              temp_tend_melt,  cmus_tot,  &
               emds_liq, emds_ice, &
               emes_liq, emes_ice, wmms, wmps, tmes, mrmes,  &
               eces_liq, eces_ice, ecds_liq, ecds_ice, disa, &
               dise, disg_2liq, disg_2ice, disf, ermesg, error)
      endif

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!----------------------------------------------------------------------


end subroutine don_d_def_conv_forcing_k



!#####################################################################

subroutine don_d_finalize_output_fields_k   &
         (nlev_lsm, ntr, i, j, Param, disb, disc_liq, disc_ice,  &
          ensmbl_freeze, ensmbl_freeze_meso, &
        temp_tend_melt,  tmes, disd, cmus_tot, ecds_liq, ecds_ice, &
         eces_liq, eces_ice, emds_liq, emds_ice, emes_liq, emes_ice, &
          wmms, wmps, mrmes, cutotal, dmeml, detmfl, temptr, uceml, &
          umeml, cuq, cuql_v, qtren, qtmes, wtp, ensmbl_wetc,   &
          Don_conv, ermesg, error)

!----------------------------------------------------------------------
!    subroutine finalize_output_fields stores output variables from 
!    columns with active deep convection into the appropriate elements 
!    of the donner_conv_type variable Don_conv.
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_conv_type 

implicit none

!---------------------------------------------------------------------
integer,                          intent(in)    :: nlev_lsm, ntr
integer,                          intent(in)    :: i, j
type(donner_param_type),          intent(in)    :: Param
real,    dimension(nlev_lsm),     intent(in)    :: disb, disc_liq, disc_ice,   &
                                                   ensmbl_freeze, &
                                                   ensmbl_freeze_meso, &
                                                   temp_tend_melt,  &
                                                   tmes, disd, cmus_tot,&
                                                   emds_liq, emds_ice, &
                                                   ecds_liq, ecds_ice, &
                                                   eces_liq, eces_ice, &
                                                         wmms, wmps,  &
                                                   emes_liq, emes_ice, &
                                                   mrmes, cutotal, &
                                                   dmeml, detmfl, uceml,&
                                                   umeml, cuq, cuql_v
real,    dimension(nlev_lsm,ntr), intent(in)    :: qtren, qtmes, wtp, &
                                                   temptr, ensmbl_wetc
type(donner_conv_type),           intent(inout) :: Don_conv
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error
!---------------------------------------------------------------------
!   intent(in) variables:
!
!       i, j           i, j indices of the current grid column
!       wmms
!       wmps
!       mrmes
!       emds
!       emes
!       ecds
!       eces
!       disd
!       cmus_tot
!       disb
!       disc
!       temp_tend_melt
!       temp_tend_freeze
!       tmes
!       cutotal
!       cuq
!       cuql_v
!       dmeml
!       detmfl
!       uceml
!       umeml
!       exit_flag
!       total_precip
!       meso_precip
!       qtren
!       qtmes
!       wtp
!
!   intent(inout) variables:
!
!       Don_conv       donner_convection_type derived type variable 
!                      containing fields produced by the donner_deep
!                      convection mod 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: k, kinv    ! do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    if deep convection occurred in this column, save various output
!    fields. if it did not, then these components of the Don_conv
!    derived-type variable will retain their default values.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    save the following arrays as elements of the donner_conv type 
!    variable Don_Conv. make sure all arrays are in mks units, requiring
!    conversion of some arrays from per day to per second, g(h2o) to 
!    kg(h2o) and / or mm to m. reverse the vertical index, making these
!    profile arrays compatible with the large-scale model grid ( index 1
!    nearest upper boundary) rather than the cloud model grid (index 1 
!    nearest sfc).
!---------------------------------------------------------------------
      do k=1,nlev_lsm            
        kinv = nlev_lsm + 1 - k
        Don_conv%ceefc (i,j,kinv)   = disb(k)/Param%seconds_per_day
        Don_conv%cecon (i,j,kinv)   = (disc_liq(k) + disc_ice(k))/Param%seconds_per_day
        Don_conv%tmes  (i,j,kinv)   = tmes(k)/Param%seconds_per_day
        Don_conv%fre   (i,j,kinv)   = (ensmbl_freeze(k) +  &
                                        ensmbl_freeze_meso(k))/  &
                                                 Param%seconds_per_day
        Don_conv%elt   (i,j,kinv)   = temp_tend_melt(k) / &
                                                   Param%seconds_per_day
        Don_conv%cmus  (i,j,kinv)   = cmus_tot(k)/       &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%ecds  (i,j,kinv)   = (ecds_liq(k) + ecds_ice(k))/    &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%eces  (i,j,kinv)   = (eces_liq(k) + eces_ice(k))/    &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%emds  (i,j,kinv)   = (emds_liq(k) + emds_ice(k))/    &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%emes  (i,j,kinv)   = (emes_liq(k) + emes_ice(k))/    &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%mrmes  (i,j,kinv)   = mrmes(k)/          &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%wmps  (i,j,kinv)   = wmps(k)/             &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%wmms  (i,j,kinv)   = wmms(k)/                 &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%cemfc (i,j,kinv)   = disd(k)/                 &
                                           (1.0E03*Param%seconds_per_day)
        Don_conv%cual  (i,j,kinv)   = cutotal(k)
        Don_conv%dmeml (i,j,kinv)   = dmeml(k)
        Don_conv%uceml (i,j,kinv)   = uceml(k)
        if (detmfl(k) <= 1.0e-10) then
          Don_conv%detmfl(i,j,kinv) = 0.
        else
          Don_conv%detmfl(i,j,kinv)   = detmfl(k)
        endif
        Don_conv%umeml (i,j,kinv)   = umeml(k)
        Don_conv%cuqi  (i,j,kinv)   = cuq(k)
        Don_conv%cuql  (i,j,kinv)   = cuql_v(k)
        Don_conv%qtren1(i,j,kinv,:) = qtren(k,:)
        Don_conv%qtmes1(i,j,kinv,:) = qtmes(k,:)
        Don_conv%temptr(i,j,kinv,:) = temptr(k,:)
        Don_conv%wtp1  (i,j,kinv,:) = wtp(k,:)
        Don_conv%wetdepc(i,j,kinv,:)= ensmbl_wetc(k,:)
      end do
        

!--------------------------------------------------------------------


end subroutine don_d_finalize_output_fields_k 



!#####################################################################

!#####################################################################

subroutine don_d_determine_cloud_area_k            &
         (me, nlev_lsm, nlev_hires, diag_unit, debug_ijt, Param,  &
          Initialized, Nml, lofactor, &
          max_depletion_rate, dcape, amax, dise_v, disa_v, pfull_c,  &
          temp_c, mixing_ratio_c, env_t, env_r, parcel_t, parcel_r, &
          cape_p, exit_flag, amos, a1, ermesg, error)

!---------------------------------------------------------------------
!    subroutine determine_cloud_area defines the convective cloud area
!    and so closes the donner_deep parameterization. The arrays 
!    Don_conv%a1 and Don_conv%amos are output by this routine.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_initialized_type
use conv_utilities_k_mod, only : sounding, adicloud 

implicit none

!-----------------------------------------------------------------------
integer,                      intent(in)    :: me, nlev_lsm,     &
                                               nlev_hires, diag_unit
logical,                      intent(in)    :: debug_ijt
type(donner_param_type),      intent(in)    :: Param
type(donner_initialized_type),intent(in)    :: Initialized
type(donner_nml_type),        intent(in)    :: Nml      
real,                         intent(in)    :: max_depletion_rate,   &
                                               lofactor, &
                                               dcape, amax
real, dimension(nlev_lsm),    intent(in)    :: dise_v, disa_v, &
                                               pfull_c, temp_c,  &
                                               mixing_ratio_c 
real, dimension(nlev_hires),  intent(in)    :: env_t, env_r, parcel_t,  &
                                               parcel_r, cape_p
logical,                      intent(inout) :: exit_flag
real,                         intent(out)   :: amos, a1
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       diag_unit          unit number for column diagnostics file
!       debug_ijt          column_diagnostics are requested in 
!                          current column  ?
!       max_depletion_rate rate of moisture depletion due to convection
!                          that would result in a column without vapor
!                          [ kg(h2o) / ( kg(air) sec ) ]      
!       dcape              time tendency of cape
!       amax
!       dise_v
!       disa_v
!       pfull_c
!       temp_c 
!       mixing_ratio_c
!       env_t
!       env_r
!       parcel_t
!       parcel_r
!       cape_p
!
!   intent(inout) variables:
!
!       exit_flag
!
!   intent(out) variables:
!
!       amos
!       a1
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
 
      real, dimension (nlev_lsm)         :: a1_vk              
      real, dimension(nlev_hires)        :: qli0_v, qli1_v, qt_v,  &
                                            qr_v, rl_v, ri_v
      real                               :: qtest, tfint, disbar
      integer                            :: k
!----------------------------------------------------------------------
!   local variables:
!
!         a1_vk
!         qli0      normalized component of cumulus condensate forcing
!         qli1      un-normalized component of condensate forcing
!         qt_v      temperature tendency due to deep convection on
!                   cape grid [ deg K / sec ]
!         qr_v      vapor mixing ratio tendency due to deep convection
!                   on cape grid [ kg(h2o) / ( kg(air) sec ]
!         rl_v      large-scale liquid mixing ratio
!         ri_v      large-scale ice mixing ratio 
!         qtest
!         tfint     column integral of moisture time tendency due to
!                   convection  [ mm / sec , or  kg / (m**2 sec ) ]
!         disbar    water vapor time tendency due to deep convection at 
!                   large-scale model interface levels
!                   [ kg(h2o) / ( kg(air) sec ) ]
!         nlev      number of layers in large-scale model
!         k         do-loop index

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    call map_lo_res_col_to_hi_res_col to interpolate moisture and
!    temperature forcings from large-scale model grid (dise_v, disa_v)
!    to the vertical grid used in the cape calculation (qr_v, qt_v). 
!--------------------------------------------------------------------
      call don_u_lo1d_to_hi1d_k   &
            (nlev_lsm, nlev_hires, disa_v, pfull_c, cape_p, qt_v, ermesg, error)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

      call don_u_lo1d_to_hi1d_k   &
            (nlev_lsm, nlev_hires, dise_v, pfull_c, cape_p, qr_v, ermesg, error)


!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    if in a diagnostic column, output the temperature and moisture 
!    forcings on both the cape grid (qt_v, qr_v) and the large-scale
!    model grid (disa_v, dise_v).
!--------------------------------------------------------------------
      if (debug_ijt) then
        do k=1,nlev_hires
          if (qr_v(k) /= 0.0 .or. qt_v(k) /= 0.0) then 
            write (diag_unit, '(a, i4, e20.12, f20.14)')  &
                     'in cupar: k,qr,qt= ',k, qr_v(k), qt_v(k)
          endif
        end do
        do k=1,nlev_lsm
          if (dise_v(k) /= 0.0 .or. disa_v(k) /= 0.0) then 
            write (diag_unit, '(a, i4, 2e20.12)')  &
                    'in cupar: k,dise,disa= ',k, dise_v(k), disa_v(k)
          endif
        end do
      endif

!--------------------------------------------------------------------
!   define condensate variables on the cape grid (qli0, qli1, rl_v, 
!   ri_v). these variables are not used in the current version of the
!   cumulus closure scheme implemented in subroutine cumulus_closure, 
!   so they are given values of 0.0.
!--------------------------------------------------------------------
      do k=1,nlev_hires
        qli0_v(k) = 0.
        qli1_v(k) = 0.
        rl_v(k)   = 0.
        ri_v(k)   = 0.
      end do

!--------------------------------------------------------------------
!    call subroutine cumulus_closure to determine cloud base cloud
!    fraction and so close the deep-cumulus parameterization.
!--------------------------------------------------------------------
      call cu_clo_cumulus_closure_k   &
           (nlev_hires, diag_unit, debug_ijt, Param, Initialized, &
            Nml, lofactor, dcape, &
            cape_p, qli0_v, qli1_v, qr_v, qt_v, env_r, ri_v, &
            rl_v, parcel_r, env_t, parcel_t, a1, ermesg, error)     

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine. 
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    calculate the vertical integral of normalized moisture forcing 
!    in the column (tfint) in units of kg (h2o) per m**2 per second, or
!    mm (h2o) per second.
!-------------------------------------------------------------------
      tfint = 0.0
      do k=2,nlev_lsm
        disbar = 0.5*(dise_v(k-1) + dise_v(k))
        tfint = tfint - disbar*(pfull_c(k-1) - pfull_c(k))
      end do
      tfint = tfint/Param%grav

!--------------------------------------------------------------------
!    restrict the cloud-base area fraction produced by subroutine
!    cumulus_closure to be no larger than the cloud base area that 
!    results in total grid box coverage at some higher level (amax). 
!--------------------------------------------------------------------
      a1 = MIN (amax, a1)

!---------------------------------------------------------------------
!    set the cloud-base area fraction to be 0.0 if there is no net
!    column integral of moisture forcing in the column. this is 
!    referred to as the moisture constraint. see "Moisture Constraint",
!    8/8/97. set the exit_flag to .true., turning off convection in
!    this column, output a message, and return to calling subprogram.
!---------------------------------------------------------------------
      if (tfint == 0.) then      
        a1 = 0.
        exit_flag      = .true.
        if (debug_ijt) then
          write (diag_unit, '(a)')  &
                 'convection turned off in column because of moist&
                  &ure constraint; cloud area being set to 0.0'
        endif
        return
      endif

!---------------------------------------------------------------------
!    if in a diagnostic column, output the column integral of the 
!    moisture forcing (tfint) and the fractional cloud area (a1) after
!    assuring that moisture forcing is present in the column.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)')  &
                      'in cupar: tfint= ',tfint       
        write (diag_unit, '(a, e20.12)')  &
                      'in cupar: a1_v = ',a1       
      endif

!---------------------------------------------------------------------
!    restrict cloud fractional area by the moisture constraint. this
!    requirement limits the cloud area so that the moisture tendency 
!    due to the deep convection (tfint - which occurs only within the 
!    cloud fractional area) will not remove more vapor from the column 
!    than is available. here amos is the cloud area over which applic-
!    ation of the convective moisture tendency will result in total
!    vapor depletion in the column.
!---------------------------------------------------------------------
      amos = max_depletion_rate/tfint     
      if (a1 > amos)  then    
        a1 = max(amos, 0.)
      endif 

!---------------------------------------------------------------------
!    for any diagnostic columns in the window in which deep convection
!    was possible, output the column integral of the moisture forcing 
!    (tfint), the max cloud area allowed by the moisture constraint 
!    (amos) and the fractional cloud area after applying the moisture
!    constraint (a1).
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, 3e20.12)')  &
                   'in cupar: tfint,amos,a1= ',  &
                                       tfint, amos, a1  
      endif

!---------------------------------------------------------------------
!    verify that the current value of a1 will not produce negative
!    value of vapor mixing ratio at any level in the column when the
!    convective moisture tendency is applied. determine the large-scale
!    model mixing ratio for the current value of a1 (qtest). if qtest
!    is negative at any level for this value of a1, reset the value 
!    of a1, so that no negative mixing ratios will be produced.
!--------------------------------------------------------------------
      do k=1,nlev_lsm
        qtest = mixing_ratio_c(k) + a1*Nml%donner_deep_freq*dise_v(k)
        if (qtest < 0.) then
          a1_vk(k) = -mixing_ratio_c(k)/(dise_v(k)*Nml%donner_deep_freq)
        else
          a1_vk(k) = a1     
        endif
      end do

!--------------------------------------------------------------------
!    define the a1 for the column as the smallest of those defined
!    in the column. 
!--------------------------------------------------------------------
      a1 = MINVAL (a1_vk)

!---------------------------------------------------------------------
!    if in a diagnostic column, output the final value of a1, after 
!    all necessary constraints have been applied.
!---------------------------------------------------------------------
      if (debug_ijt) then
        write (diag_unit, '(a, e20.12)') 'in cupar: a1= ',a1        
      endif


!--------------------------------------------------------------------


end subroutine don_d_determine_cloud_area_k 





!####################################################################

!######################################################################



!######################################################################

subroutine don_d_remove_normalization_k   &
      (isize, jsize, nlev_lsm, ntr, exit_flag, Don_conv, total_precip, &
       Initialized, &
       temperature_forcing, moisture_forcing, ermesg, error)

!---------------------------------------------------------------------
!    subroutine remove_normalization removes the normalization by the
!    cloud base fractional area from the various convective diagnostics
!    and output fields so that they are ready fro use in the large-scale
!    model equations.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_nml_type, &
                              donner_initialized_type, DET_MASS_FLUX, &
                              MASS_FLUX, CELL_UPWARD_MASS_FLUX, &
                              TEMP_FORCING, MOIST_FORCING, PRECIP, &
                              FREEZING, RADON_TEND

implicit none 

!---------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize, nlev_lsm, ntr
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
type(donner_conv_type),           intent(inout) :: Don_conv
real   , dimension(isize,jsize),  intent(inout) :: total_precip
type(donner_initialized_type),    intent(inout) :: Initialized
real   , dimension(isize,jsize,nlev_lsm),                 &
                                  intent(inout) :: temperature_forcing, &
                                                   moisture_forcing
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error
!----------------------------------------------------------------------
!   intent(in) variables:
!
!     exit_flag      logical array indicating whether donner convection
!                    is not active (.true.) or is active (.false.) in
!                    each model column 
!
!   intent(inout) variables:
!    
!     Don_conv       donner_convection_type derived type variable 
!                    containing fields produced by the donner_deep
!                    convection mod 
!     total_precip   precipitation generated by deep convection
!                    [ kg / m**2 ]
!     moisture_forcing
!                    time tendency of vapor mixing ratio due to deep 
!                    convection [ kg(h2o) / kg(dry air) / sec ]
!     temperature_forcing
!                    time tendency of temperature due to deep 
!                    convection [ deg K / sec ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension(nlev_lsm) :: variable
      integer :: i, j, k, n    ! do-loop indices

!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    remove normalization from the cumulus diagnostics and forcing terms
!    by multiplying them by the fractional cloud base area. these values
!    thus become grid-box averages, rather than averages over the cloudy
!    area, and so are appropriate to use in the large-scale model
!    equations. 
!---------------------------------------------------------------------
      do j=1,jsize                          
        do i=1,isize

!---------------------------------------------------------------------
!    if deep convection is present in the column, denormalize the 
!    convective fields.
!---------------------------------------------------------------------
          if (.not. exit_flag(i,j)) then
            if (Initialized%monitor_output) then
              do n=1, size(Initialized%Don_monitor, 1)
                select case (Initialized%Don_monitor(n)%index)
                  case (DET_MASS_FLUX)
                     variable(:) = Don_conv%detmfl(i,j,:)*   &
                                                       Don_conv%a1(i,j)
                     call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (MASS_FLUX)
                     variable(:) =   &
                      (Don_conv%umeml(i,j,:) + Don_conv%dmeml(i,j,:) + &
                        Don_conv%uceml(i,j,:))*Don_conv%a1(i,j)
                     call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (CELL_UPWARD_MASS_FLUX)
                    variable(:) = Don_conv%uceml(i,j,:)*Don_conv%a1(i,j)
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (TEMP_FORCING)
                    variable(:) =   &
                           temperature_forcing(i,j,:)*Don_conv%a1(i,j)
                     call don_u_process_monitor_k (variable, i, j,  &
                             nlev_lsm, Initialized%Don_monitor(n))
                  case (MOIST_FORCING)
                     variable(:) =   &
                              moisture_forcing(i,j,:)*Don_conv%a1(i,j)
                     call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (PRECIP)
                    variable(:) = total_precip(i,j)*Don_conv%a1(i,j)
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
                  case (FREEZING)
                    variable(:) = Don_conv%fre(i,j,:)*Don_conv%a1(i,j)
                    call don_u_process_monitor_k (variable, i, j,  &
                                 nlev_lsm, Initialized%Don_monitor(n))
                end select
             end do
          endif
         total_precip(i,j) =  total_precip(i,j)*Don_conv%a1(i,j)
            Don_conv%ampta1(i,j) =  Don_conv%ampta1(i,j)*Don_conv%a1(i,j)
            Don_conv%cell_precip(i,j) =              &
                             Don_conv%cell_precip (i,j)*Don_conv%a1(i,j)
            Don_conv%meso_precip(i,j) =              &
                             Don_conv%meso_precip (i,j)*Don_conv%a1(i,j)
            Don_conv%emdi_v(i,j) = Don_conv%emdi_v(i,j)*Don_conv%a1(i,j)
            do k=1,nlev_lsm                           
               Don_conv%wetdepc(i,j,k,:) = &
                              Don_conv%wetdepc(i,j,k,:)*Don_conv%a1(i,j)
               Don_conv%wetdept(i,j,k,:) = &
                               Don_conv%wetdepc(i,j,k,:)
              temperature_forcing(i,j,k) =   &
                             temperature_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ceefc(i,j,k) =   &
                                  Don_conv%ceefc(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cecon(i,j,k) =        &
                                  Don_conv%cecon(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cemfc(i,j,k) =      &
                                  Don_conv%cemfc(i,j,k)*Don_conv%a1(i,j)
              moisture_forcing(i,j,k) =      &
                                moisture_forcing(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cual (i,j,k) =       &
                                   Don_conv%cual(i,j,k)*Don_conv%a1(i,j)
              Don_conv%fre(i,j,k) = Don_conv%fre(i,j,k)*Don_conv%a1(i,j)
              Don_conv%elt(i,j,k) = Don_conv%elt(i,j,k)*Don_conv%a1(i,j)
              Don_conv%cmus(i,j,k) =      &
                                   Don_conv%cmus(i,j,k)*Don_conv%a1(i,j)
              Don_conv%ecds(i,j,k) =      &
                                   Don_conv%ecds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%eces(i,j,k) =      &
                                   Don_conv%eces(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emds(i,j,k) =       &
                                   Don_conv%emds(i,j,k)*Don_conv%a1(i,j)
              Don_conv%emes(i,j,k) =       &
                                   Don_conv%emes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%mrmes(i,j,k) =       &
                                   Don_conv%mrmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmps(i,j,k) =       &
                                   Don_conv%wmps(i,j,k)*Don_conv%a1(i,j)
              Don_conv%wmms(i,j,k) =      &
                                   Don_conv%wmms(i,j,k)*Don_conv%a1(i,j)
              Don_conv%tmes(i,j,k) =      &
                                   Don_conv%tmes(i,j,k)*Don_conv%a1(i,j)
              Don_conv%dmeml(i,j,k) =      &
                                  Don_conv%dmeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%uceml(i,j,k) =      &
                                  Don_conv%uceml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%detmfl(i,j,k) =      &
                                  Don_conv%detmfl(i,j,k)*Don_conv%a1(i,j)
              Don_conv%umeml(i,j,k) =      &
                                  Don_conv%umeml(i,j,k)*Don_conv%a1(i,j)
              Don_conv%qtren1(i,j,k,:) =     &
                               Don_conv%qtren1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtmes1(i,j,k,:) =     &
                               Don_conv%qtmes1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%wtp1(i,j,k,:) =       &
                                 Don_conv%wtp1(i,j,k,:)*Don_conv%a1(i,j)
              Don_conv%qtceme(i,j,k,:) =   &
                     Don_conv%qtmes1(i,j,k,:) + Don_conv%qtren1(i,j,k,:)
            end do
        if (Initialized%monitor_output) then
              do n=1, size(Initialized%Don_monitor, 1)
                select case (Initialized%Don_monitor(n)%index)
                  case (RADON_TEND)
                    variable(:) = Don_conv%qtceme   &
                         (i,j,:,Initialized%Don_monitor(n)%tracer_index)
                    call don_u_process_monitor_k (variable, i, j,  &
                                   nlev_lsm, Initialized%Don_monitor(n))
 
                 end select
               end do
            endif

!---------------------------------------------------------------------
!    if deep convection is not present in the column, define the output
!    fields appropriately.
!---------------------------------------------------------------------
          else
            total_precip(i,j) = 0.
            do k=1,nlev_lsm
              temperature_forcing(i,j,k) = 0.
              moisture_forcing(i,j,k) = 0.
            end do
          endif

        end do
      end do

!---------------------------------------------------------------------


end subroutine don_d_remove_normalization_k



!######################################################################

subroutine don_d_output_cupar_diags_k    &
         (isize, jsize, nlev_lsm, Col_diag, n, exit_flag, &
          total_precip, temperature_forcing, Don_conv, Don_cape, ermesg, error)

!----------------------------------------------------------------------
!----------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_cape_type, &
                             donner_column_diag_type

implicit none

!----------------------------------------------------------------------
integer,                          intent(in)    :: isize, jsize,  &
                                                   nlev_lsm
type(donner_column_diag_type),    intent(in)    :: Col_diag
integer,                          intent(in)    :: n
logical, dimension(isize,jsize),  intent(in)    :: exit_flag
real, dimension (isize,jsize),    intent(in)    :: total_precip
real, dimension (isize,jsize,nlev_lsm),                      &
                                  intent(in)    :: temperature_forcing
type(donner_conv_type),           intent(inout) :: Don_conv
type(donner_cape_type),           intent(inout) :: Don_cape
character(len=*),                 intent(out)   :: ermesg
integer,                          intent(out)   :: error

      integer  :: idiag, jdiag, unitdiag
      integer  :: i,j,k


!-----------------------------------------------------------------------
!    initialize the error message character string.
!-----------------------------------------------------------------------
      ermesg = ' ' ; error = 0
      idiag = Col_diag%i_dc(n)
      jdiag = Col_diag%j_dc(n)
      unitdiag = Col_diag%unit_dc(n)

!---------------------------------------------------------------------
!    find any other columns in the current physics window in which deep
!    convection has produced a precipitation rate of over 1 mm/day. 
!    output the precipitation rate and cloud areas for each of these 
!    columns.
!---------------------------------------------------------------------
      do j=1,jsize
        do i=1,isize       
          if (.not. exit_flag(i,j) ) then
            if (Don_conv%cell_precip(i,j) > 1.) then
              write (unitdiag, '(a)')  &
                     ' the following columns in the current physics&
                         & window contain deep convection producing&
                         & rainfall rates of over 1mm per day '
              write (unitdiag, '(a, 2i4, 3e20.12)')  &
                        'in cupar: i,j, precip rate, cloud area, &
                                      &anvil area = ,',  &
                               i, j, Don_conv%cell_precip(i,j),    &
                               Don_conv%a1(i,j), Don_conv%ampta1(i,j)
            endif
          endif
        end do
      end do

!---------------------------------------------------------------------
!    if in a diagnostic window, output convection-related upper tropo-
!    spheric heating rates if there is convective precipitation in any
!    of the diagnostic columns.
!---------------------------------------------------------------------
      if (total_precip(idiag,jdiag) /= 0.) then
        do k=Col_diag%kstart,nlev_lsm
          if ((Don_cape%model_p(idiag,jdiag,k) > 100.e02) .and.&
              (Don_cape%model_p(idiag,jdiag,k) < 500.e02)) then 
            if (temperature_forcing(idiag,jdiag,nlev_lsm-k+1) /= 0.) then
              write (unitdiag, '(a, 3i4, f20.14)')    &
                     'in cupar: j_dc,i_dc,k,t= ',  &
                               jdiag, idiag, k,    &
                                    Don_cape%model_t(idiag,jdiag,k)
              write (unitdiag, '(a, e20.12, i4, 2e20.12)')&
                     'in cupar: tprea1,k,pr,cemetf= ',  &
                            total_precip(idiag,jdiag), k,    &
                            Don_cape%model_p(idiag,jdiag,k),   &
                       temperature_forcing(idiag,jdiag,nlev_lsm-k+1 )
            endif
          endif
        end do
      endif

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of convective and total 
!    precipitation and cloud areas,
!----------------------------------------------------------------------
      if (.not. exit_flag(idiag,jdiag) ) then
        write (unitdiag, '(a, 2e20.12)')  &
                     'in cupar: contot,tpre=', &
                 Don_conv%cell_precip(idiag,jdiag) /  &
                                           (total_precip(idiag,jdiag)),&
                        total_precip(idiag,jdiag)
        write (unitdiag, '(a, 2e20.12)') 'in cupar: a1,ampt =',  &
                         Don_conv%a1 (idiag,jdiag), &
                         Don_conv%ampta1(idiag,jdiag)
        write (unitdiag, '(a, e20.12)')  'in cupar: amax= ', &
                          Don_conv%amax(idiag,jdiag)

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of mesoscale and 
!    cell-scale mass fluxes.
!----------------------------------------------------------------------
        do k=Col_diag%kstart,nlev_lsm
          write (unitdiag, '(a, i4, f19.10, 3e20.12)')  &
                 'in cupar: k,pr,uceml,dmeml,umeml= ',  &
                     k,  Don_cape%model_p(idiag,jdiag,nlev_lsm-k+1),  &
                         Don_conv%uceml(idiag,jdiag,k), &
                         Don_conv%dmeml(idiag,jdiag,k),  &
                         Don_conv%umeml(idiag,jdiag,k)
        end do

!----------------------------------------------------------------------
!    if in a diagnostic window, output values of cloud liquid (cuql).
!    at any levels at which heating associated with the donner deep
!    convection is greater than 0.002 deg K / sec, output heating rate,
!    cloud area, cape, cape tendency, and cloud area.
!----------------------------------------------------------------------
        do k=Col_diag%kstart,nlev_lsm
          write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k,cuql', &
                              k,Don_conv%cuql (idiag,jdiag    ,k)
          if (ABS(temperature_forcing(idiag,jdiag,k)) > 0.002) then
            write (unitdiag, '(a, i4, e20.12)')  &
                             'in donner_deep: k, cemetf= ',k,   &
                                 temperature_forcing(idiag,jdiag,k)
            write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k, cual= ',k,    &
                                    Don_conv%cual(idiag,jdiag,k )
            write (unitdiag, '(a, i4, e20.12)')  &
                            'in donner_deep: k, xcape= ',k,    &
                                  Don_cape%xcape_lag(idiag,jdiag) 
            write (unitdiag, '(a, i4, e20.12)')   &
                              'in donner_deep: k, dcape = ',k,    &
                                  Don_conv%dcape(idiag,jdiag)
            write (unitdiag, '(a, i4, e20.12)')  &
                              'in donner_deep: k,a1    = ',k,    &
                                  Don_conv%a1 (idiag,jdiag)
            write (unitdiag, '(a, i4, e20.12)')   &
                              'in donner_deep: k, amax  = ',k,   &
                                  Don_conv%amax(idiag,jdiag)
          endif
        end do
      endif   ! (not exit_flag)

!--------------------------------------------------------------------



end subroutine don_d_output_cupar_diags_k



!####################################################################

subroutine don_d_dealloc_loc_vars_k   &
         (Don_conv, Don_cape, Don_rad, Don_cem, Don_budgets, Nml,   &
          Initialized, sd, ac, cp, ct, ermesg, error)

!----------------------------------------------------------------------
!    subroutine don_d_dealloc_loc_vars_k deallocates the
!    local variables found in subroutine donner_deep of donner_deep_mod.
!    these are limited to the pointer components of the donner_conv_type,
!    donner_cape_type, donner_rad_type and donner_cem_type arrays 
!    resident there.
!----------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type, donner_cape_type, &
                             donner_rad_type, donner_budgets_type, &
                             donner_cem_type, &
                             donner_nml_type, donner_initialized_type
use  conv_utilities_k_mod,only : adicloud, sounding, ac_end_k, &
                                 sd_end_k
use  conv_plumes_k_mod,only    : cplume, ctend, cp_end_k, ct_end_k

implicit none

!----------------------------------------------------------------------
type(donner_conv_type),         intent(inout) :: Don_conv
type(donner_cape_type),         intent(inout) :: Don_cape
type(donner_rad_type),          intent(inout) :: Don_rad 
type(donner_cem_type),          intent(inout) :: Don_cem
type(donner_budgets_type),      intent(inout) :: Don_budgets
type(donner_nml_type),          intent(inout) :: Nml         
type(donner_initialized_type),  intent(inout) :: Initialized 
type(sounding),                 intent(inout) ::  sd
type(adicloud),                 intent(inout) ::  ac
type(cplume),                   intent(inout) ::  cp
type(ctend),                    intent(inout) ::  ct
character(len=*),               intent(out)   :: ermesg
integer,                        intent(out)   :: error
!----------------------------------------------------------------------
!   intent(inout) variables:
!
!     Don_conv             donner_convection_type derived type variable
!                          containing diagnostics and intermediate
!                          results describing the nature of the convec-
!                          tion produced by the donner parameterization
!     Don_cape             donner_cape type derived type variable con-
!                          taining diagnostics and intermediate results
!                          related to the cape calculation associated
!                          with the donner convection parameterization
!     Don_rad              donner_rad_type derived type variable used
!                          to hold those fields needed to connect the
!                          donner deep convection parameterization and
!                          the model radiation package
!     Don_cem              donner_cem_type derived type variable 
!                          containing Donner cumulus ensemble member 
!                          diagnostics
!
!  intent(out) variables:
! 
!     ermesg               character string containing any error message
!                          to be returned to calling routine
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the error message string.
!---------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    deallocate the components of the donner lite derived types.  
!------------------------------------------------------------------
      call sd_end_k(sd)
      call ac_end_k(ac)
      call cp_end_k(cp)
      call ct_end_k(ct)
!----------------------------------------------------------------------
!    deallocate the components of the donner_conv_type variable.
!----------------------------------------------------------------------
      deallocate (Don_conv%conv_temp_forcing    )
      deallocate (Don_conv%conv_moist_forcing   )
      deallocate (Don_conv%ceefc                )
      deallocate (Don_conv%cecon                )
      deallocate (Don_conv%cemfc                )
      deallocate (Don_conv%cememf_mod           )
      deallocate (Don_conv%cual                 )
      deallocate (Don_conv%fre                  )
      deallocate (Don_conv%elt                  )
      deallocate (Don_conv%cmus                 )
      deallocate (Don_conv%ecds                 )
      deallocate (Don_conv%eces                 )
      deallocate (Don_conv%emds                 )
      deallocate (Don_conv%emes                 )
      deallocate (Don_conv%mrmes                )
      deallocate (Don_conv%wmps                 )
      deallocate (Don_conv%wmms                 )
      deallocate (Don_conv%tmes                 )
      deallocate (Don_conv%dmeml                )
      deallocate (Don_conv%uceml                )
      deallocate (Don_conv%detmfl               )
      deallocate (Don_conv%umeml                )
      deallocate (Don_conv%xice                 ) 
      deallocate (Don_conv%xliq                 )
      deallocate (Don_conv%qtren1               )
      deallocate (Don_conv%qtceme               )
      deallocate (Don_conv%qtmes1               )
      deallocate (Don_conv%temptr               )
      deallocate (Don_conv%wtp1                 )
      deallocate (Don_conv%wetdepc              )
      deallocate (Don_conv%wetdepm              )
      deallocate (Don_conv%wetdept              )
      deallocate (Don_conv%dgeice               )
      deallocate (Don_conv%cuqi                 )
      deallocate (Don_conv%cuql                 )
      deallocate (Don_conv%cell_liquid_eff_diam )
      deallocate (Don_conv%cell_ice_geneff_diam )
      deallocate (Don_conv%dcape                )  
      deallocate (Don_conv%a1                   )
      deallocate (Don_conv%amax                 )
      deallocate (Don_conv%amos                 )
      deallocate (Don_conv%ampta1               )
      deallocate (Don_conv%cell_precip          )
      deallocate (Don_conv%meso_precip          )
      deallocate (Don_conv%emdi_v               )
      deallocate (Don_conv%prztm                )
      deallocate (Don_conv%przm                 )
      deallocate (Don_conv%pb_v                 )
      deallocate (Don_conv%pmd_v                )
      deallocate (Don_conv%pztm_v               )
      deallocate (Don_conv%pzm_v                )

!----------------------------------------------------------------------
!    deallocate the components of the donner_cape_type variable.
!----------------------------------------------------------------------
      deallocate (Don_cape%coin       )
      deallocate (Don_cape%plcl       )
      deallocate (Don_cape%plfc       )
      deallocate (Don_cape%plzb       )
      deallocate (Don_cape%xcape      )
      deallocate (Don_cape%xcape_lag  )
      deallocate (Don_cape%parcel_r   )
      deallocate (Don_cape%parcel_t   )
      deallocate (Don_cape%cape_p     )
      deallocate (Don_cape%env_r      )
      deallocate (Don_cape%env_t      )
      deallocate (Don_cape%model_p    )
      deallocate (Don_cape%model_r    )
      deallocate (Don_cape%model_t    )
      deallocate (Don_cape%qint       )     
      deallocate (Don_cape%qint_lag   ) 

!----------------------------------------------------------------------
!    deallocate the components of the donner_rad_type variable.
!----------------------------------------------------------------------
      deallocate (Don_rad%cell_cloud_frac  )
      deallocate (Don_rad%cell_liquid_amt  )
      deallocate (Don_rad%cell_liquid_size )
      deallocate (Don_rad%cell_ice_amt     )
      deallocate (Don_rad%cell_ice_size    )
      deallocate (Don_rad%cell_droplet_number )
      deallocate (Don_rad%meso_cloud_frac  )
      deallocate (Don_rad%meso_liquid_amt  )
      deallocate (Don_rad%meso_liquid_size )
      deallocate (Don_rad%meso_ice_amt     )
      deallocate (Don_rad%meso_ice_size    )
      deallocate (Don_rad%meso_droplet_number )
      deallocate (Don_rad%nsum             )        

   if (Nml%do_ensemble_diagnostics) then
!--------------------------------------------------------------------
!    deallocate the components of the donner_cem_type variable.
!--------------------------------------------------------------------
      deallocate (Don_cem%pfull       )
      deallocate (Don_cem%phalf       )
      deallocate (Don_cem%zfull       )
      deallocate (Don_cem%zhalf       )
      deallocate (Don_cem%temp        )
      deallocate (Don_cem%mixing_ratio )
      deallocate (Don_cem%cell_precip )
      deallocate (Don_cem%meso_precip )
      deallocate (Don_cem%pb          )
      deallocate (Don_cem%ptma        )
      deallocate (Don_cem%h1          )
      deallocate (Don_cem%qlw         )
      deallocate (Don_cem%cfracice    )
      deallocate (Don_cem%wv          )
      deallocate (Don_cem%rcl         )
      deallocate (Don_cem%a1          )
      deallocate (Don_cem%cual        )
      deallocate (Don_cem%temperature_forcing )
   endif

!----------------------------------------------------------------------
!    deallocate the components of the donner_budgets_type variable.
!----------------------------------------------------------------------
      deallocate (Don_budgets%liq_prcp    )
      deallocate (Don_budgets%frz_prcp    )
   if (Initialized%do_conservation_checks .or.    &
                                           Nml%do_budget_analysis) then
      deallocate (Don_budgets%lheat_precip)
      deallocate (Don_budgets%vert_motion )
      deallocate (Don_budgets%water_budget)
      deallocate (Don_budgets%enthalpy_budget)
      deallocate (Don_budgets%precip_budget)
   endif

!----------------------------------------------------------------------


end subroutine don_d_dealloc_loc_vars_k 


!######################################################################

!++lwh
subroutine don_d_check_trc_rlzbility( isize, jsize, nlev_lsm, ntr, dt, &
                                             tracers, Don_conv )
!---------------------------------------------------------------------
!  Check for tracer realizability. If convective tendencies would
!  produce negative tracer mixing ratios, scale down tracer tendency
!  terms uniformly for this tracer throughout convective column. This is
!  equivalent to limiting the cell areas.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_conv_type

!---------------------------------------------------------------------
!  Dummy arguments
!---------------------------------------------------------------------
integer,                 intent(in)     :: isize, jsize, nlev_lsm, ntr
real,                    intent(in)     :: dt 
real, dimension(isize,jsize,nlev_lsm,ntr), &
                         intent(in)     :: tracers        
type(donner_conv_type),  intent(inout)  :: Don_conv

!---------------------------------------------------------------------
!   intent(in) variables:
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     isize          x-direction size of the current physics window
!     jsize          y-direction size of the current physics window
!     nlev_lsm       number of model layers in large-scale model
!     dt             physics time step [ sec ]
!
!   intent(inout) variables:
!     Don_conv       donner_convection_type derived type variable
!                    containing diagnostics and intermediate results 
!                    describing the nature of the convection produced by
!                    the donner parameterization
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  Local variables
!---------------------------------------------------------------------

   integer :: i,j,n,k
   real, dimension(nlev_lsm) :: tracer0, trtend, trtendw, tracer1,  &
                                tracer1w
   real :: ratio, tracer_max, tracer_min

!---------------------------------------------------------------------
!   local variables:
!
!     tracers        tracer mixing ratios of tracers transported by the
!                    donner deep convection parameterization
!                    [ tracer units, e.g., kg(tracer) / kg (dry air) ]
!     tracer0        column tracer mixing ratios before convection
!     trtend         column tracer mixing ratio tendencies due to convective transport [ (tracer units) / s ]
!     trtendw        column tracer mixing ratio tendencies due to convective transport + wetdep [ (tracer units) / s ]
!     tracer1        column tracer mixing ratios after convective transport only
!     tracer1w       column tracer mixing ratios after convective transport + wet deposition
!     i, j, k, n     do-loop indices
!     ratio          ratio by which tracer convective tendencies need to 
!                    be reduced to permit realizability (i.e., to prevent
!                    negative tracer mixing ratios)
!
!---------------------------------------------------------------------

   do n = 1,ntr
   do i = 1,isize
   do j = 1,jsize
      
      tracer0(:)  = tracers(i,j,:,n)
      trtend(:)   = Don_conv%qtceme(i,j,:,n)
      trtendw(:)  = trtend(:) + Don_conv%wetdept(i,j,:,n)
      tracer1(:)  = tracer0 + dt * trtend(:)
      tracer1w(:) = tracer0 + dt * trtendw(:)
 
      tracer_min = 1.e20
      tracer_max = -1.e20

      do k = 1,nlev_lsm
         if (trtend(k) /= 0.) then
            tracer_max = max(tracer0(k),tracer_max)
            tracer_min = min(tracer0(k),tracer_min)
         end if
      end do
       
      ratio = 1.
      do k = 1,nlev_lsm
         if (tracer0(k) > 0. .and. tracer1w(k)<0.) then
            ratio = MIN( ratio,tracer0(k)/(-trtendw(k)*dt) )
         end if
         if (tracer1(k)<tracer_min .and. trtend(k) /= 0.0 ) then
           ratio = MIN( ratio,(tracer0(k)-tracer_min)/(-trtend(k)*dt) )
         end if
         if (tracer1(k)>tracer_max  .and. trtend(k) /= 0.0 ) then
            ratio = MIN( ratio,(tracer_max-tracer0(k))/(trtend(k)*dt) )
         end if
      end do
      ratio = MAX(0.,MIN(1.,ratio))
      if (ratio /= 1.) then
         Don_conv%qtceme(i,j,:,n)  = Don_conv%qtceme(i,j,:,n)  * ratio
         Don_conv%qtren1(i,j,:,n)  = Don_conv%qtren1(i,j,:,n)  * ratio
         Don_conv%qtmes1(i,j,:,n)  = Don_conv%qtmes1(i,j,:,n)  * ratio
         Don_conv%wtp1(i,j,:,n)    = Don_conv%wtp1(i,j,:,n)    * ratio
         Don_conv%wetdepc(i,j,:,n) = Don_conv%wetdepc(i,j,:,n) * ratio
         Don_conv%wetdepm(i,j,:,n) = Don_conv%wetdepm(i,j,:,n) * ratio
         Don_conv%wetdept(i,j,:,n) = Don_conv%wetdept(i,j,:,n) * ratio
      end if
   end do
   end do
   end do


end subroutine don_d_check_trc_rlzbility
!--lwh


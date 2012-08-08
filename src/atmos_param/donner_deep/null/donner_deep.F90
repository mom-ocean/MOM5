                       module donner_deep_mod

use time_manager_mod,       only: time_type

use fms_mod,                only: error_mesg, FATAL, WARNING

implicit none
private

character(len=128)  :: version =  '$Id: donner_deep.F90,v 19.0 2012/01/06 20:09:14 fms Exp $'
character(len=128)  :: tagname =  '$Name: siena_201207 $'

public :: donner_deep_init, donner_deep, donner_deep_end, &
          donner_deep_restart, donner_deep_time_vary, donner_deep_endts

                          contains

!#####################################################################
subroutine donner_deep_init (lonb, latb, pref, axes, secs, days, &
                             tracers_in_donner, do_conservation_checks,&
                             using_unified_closure, using_fms_code)
real,            dimension(:,:), intent(in)   :: lonb, latb
real,            dimension(:),   intent(in)   :: pref
integer,         dimension(4),   intent(in)   :: axes
integer,                         intent(in)   :: secs, days
logical,         dimension(:),   intent(in)   :: tracers_in_donner
logical,                         intent(in)   :: do_conservation_checks
logical,                         intent(in)   :: using_unified_closure
logical,                         intent(in), optional :: &
                                                 using_fms_code

call error_mesg('donner_deep_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep_init

!###################################################################

subroutine donner_deep (is, ie, js, je, dt, temp, mixing_ratio, pfull, &
                        phalf, zfull, zhalf, omega, pblht, tkemiz, &
                        qstar, cush, coldT, land, sfc_sh_flux,  &
                        sfc_vapor_flux, tr_flux, tracers, secs, days, &
                        cbmf, cell_cld_frac,  &
                        cell_liq_amt, cell_liq_size, cell_ice_amt,   &
                        cell_ice_size, cell_droplet_number, &
                        meso_cld_frac, meso_liq_amt, &
                        meso_liq_size, meso_ice_amt, meso_ice_size,  &
                        meso_droplet_number, &
                        nsum, precip, delta_temp, delta_vapor, detf, &
                        uceml_inter, mtot, mfluxup, mhalf_3d,&
                        donner_humidity_area,    &
                        donner_humidity_factor, qtrtnd, donner_wetdep,&
                        lheat_precip, vert_motion,        &
                        total_precip, liquid_precip, frozen_precip, &
                        frz_meso, liq_meso, frz_cell, liq_cell, &
                        qlin, qiin, qain,              &      ! optional
                        delta_ql, delta_qi, delta_qa)         ! optional
                        
!-------------------------------------------------------------------
!    donner_deep is the prognostic driver subroutine of donner_deep_mod.
!    it takes as input the temperature (temp), vapor mixing ratio 
!    (mixing_ratio), pressure at full and half-levels (pfull, phalf),
!    vertical velocity at full levels (omega), the large scale cloud 
!    variables (qlin, qiin, qain), the land fraction (land),  the heat 
!    (sfc_sh_flux) , moisture (sfc_vapor_flux) and tracer (tr_flux) 
!    fluxes across the surface that are to be seen by this parameter-
!    ization, the tracers to be transported by the donner convection
!    parameterization (tracers), and the current time (as time_type 
!    variable Time). the routine returns the precipitation (precip),
!    increments to the temperature (delta_temp) and mixing ratio 
!    (delta_vapor), the detrained mass flux (detf), upward cell mass 
!    flux at interface levels  (uceml_inter) and total mass flux at full
!    levels (mtot), two arrays needed to connect the donner convection 
!    and strat cloud parameterizations (donner_humidity_area, 
!    donner_humidity_ratio), increments to the cloudwater (delta_ql), 
!    cloudice (delta_qi) and cloud area (delta_qa) fields and tendencies
!    for those tracers that are to be transported by the donner convect-
!    ion parameterization (qtrtnd). there are an additional eleven arrays
!    defining the donner scheme cloud characteristics needed by the rad-
!    iation package, which are passed in and updated on donner calcul-
!    ation steps.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
real,                         intent(in)    :: dt
real, dimension(:,:,:),       intent(in)    :: temp, mixing_ratio, &
                                               pfull, phalf, zfull, zhalf, omega
real, dimension(:,:),         intent(in)    :: pblht, tkemiz, qstar,cush
real, dimension(:,:),         intent(in)    :: land
logical, dimension(:,:),      intent(in)    :: coldT
real, dimension(:,:),         intent(in)    :: sfc_sh_flux, &
                                               sfc_vapor_flux
real, dimension(:,:,:),       intent(in)    :: tr_flux 
real, dimension(:,:,:,:),     intent(in)    :: tracers 
integer,                      intent(in)    :: secs, days
real, dimension(:,:),         intent(inout) :: cbmf              
real, dimension(:,:,:),       intent(inout) :: cell_cld_frac,  &
                                               cell_liq_amt,  &
                                               cell_liq_size, &
                                               cell_ice_amt,  &
                                               cell_ice_size, &
                                           cell_droplet_number, &
                                               meso_cld_frac,  &
                                               meso_liq_amt, &
                                               meso_liq_size, &
                                               meso_ice_amt,   &
                                               meso_ice_size, &
                                           meso_droplet_number
integer, dimension(:,:),      intent(inout) :: nsum
real, dimension(:,:),         intent(out)   :: precip, &
                                               lheat_precip, &
                                               vert_motion, &
                                               total_precip
real, dimension(:,:,:),       intent(out)   :: delta_temp, delta_vapor,&
                                               detf, uceml_inter, &
                                               mtot, mfluxup, &
                                               mhalf_3d, &
                                               donner_humidity_area,&
                                               donner_humidity_factor, &
                                               liquid_precip, &
                                               frozen_precip, frz_meso,&
                                            liq_meso, frz_cell, liq_cell
real, dimension(:,:,:,:),     intent(out)   :: qtrtnd 
real, dimension(:,:,:),       intent(out)   :: donner_wetdep
real, dimension(:,:,:),       intent(in),                &
                                   optional :: qlin, qiin, qain
real, dimension(:,:,:),       intent(out),               &
                                   optional :: delta_ql, delta_qi, &
                                               delta_qa


call error_mesg('donner_deep', &
      'This module is not supported as part of the public release', FATAL)

end subroutine donner_deep

!####################################################################

subroutine donner_deep_end

call error_mesg('donner_deep', &
      'This module is not supported as part of the public release', WARNING)

end subroutine donner_deep_end

!######################################################################

!#######################################################################
! <SUBROUTINE NAME="donner_deep_restart">
!
! <DESCRIPTION>
!  Dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine donner_deep_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

end subroutine donner_deep_restart
! </SUBROUTINE> NAME="donner_deep_restart"

!######################################################################
subroutine donner_deep_time_vary (dt)  
                                  
real, intent(in) :: dt

end subroutine donner_deep_time_vary

!######################################################################
subroutine donner_deep_endts

end subroutine donner_deep_endts

end module donner_deep_mod


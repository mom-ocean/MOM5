module auscom_ice_parameters_mod

! Set parameters for ice formation/melting

implicit none

integer :: il_out
integer :: dt_cpl
logical :: pop_icediag = .true.
logical :: do_ice = .true. !determine if ice formation should be calculated!    
logical :: do_ice_once = .true. 
           !if .t. call ice_formation only once per coupling interval
           !   .f. call ice_formation every time step!   (21/07/2008)
logical :: fixmeltT = .true.
logical :: limit_srfstress = .false.  
           !.t. in case surface stress from ice is too big, upsetting the model!
           !should not be an issue after ice-ocn stress error fixed (20090311) 
logical :: use_ioaice = .true. 
           !use aice from cice to determine if ice exists
real    :: aice_cutoff = .15
           !if ioaice >= aice_cutoff then ice exists
real    :: tlthk0 = 10.0  !top layer thickness (m) (in 'IC', akward, 20100423)
integer :: kmxice = 1     !max layer that allows for ice formation	
real    :: Tmelt = -1.0   !ice melting point (C)
real    :: Mstress = 2.0  !maximum surface stress allowed into ocn.
                          !(only applys when limit_srfstress = .t.)
real    :: icemlt_factor = 1.0  !potential ice melting reduction factor
!real    :: frazil_factor = 0.5  !ocn/ice 'time-stepping' associated factor
   ! frazil_factor is associated with the difference between ocean model 
   ! and ice model time-stepping: for mom4, two-level frog-leap is used 
   ! and frazil heat flux is calculated and accumulated with 
   ! frazil_factor = 1, which is supposed to be used for a ice model with 
   ! the same two-level time-stepping scheme such as SIS. 
   ! but cice uses forward time-stepping, which means we need 'correct' 
   ! the frazil energy by multiplying 0.5 before sending to cice.
real    :: frazil_factor = 1.0  !CH: mom4 and cice use same 2-level time-steeping!!!
real    :: sign_stflx = 1.0    !if saltflux from cice needs change sign.
logical :: iceform_adj_salt = .false.   
   ! option for adjusting salinity when frazil is calculated in mom4 using POP 
   ! approach (see comment in ice_formation_new routine)
!20110802: low background_diffusivity in the equatorial zone (defined by lat_low_bgdiff) 
!real    :: lat_low_bgdiff = 10    !degrees S/N (beyond this zone set to the 'background_diffusivity') 
!real    :: bg_diff_eq = 1.0e-6    !lowest bg_diff at equator
!

!maximum depth/level
integer :: ksmax = 5        !deepest level of the Red Sea/Gulf Bay 
!
integer :: sfix_hours = 12  !do s mixing every sfix_hours.
!
logical :: redsea_gulfbay_sfix = .true.
logical :: do_sfix_now = .true.
logical :: chk_i2o_fields = .false.
logical :: chk_o2i_fields = .false.
! How often to dump the coupling fields if either of the above options are .true.
! The unit of time is seconds. By default fields are dumped every timestep.
integer :: chk_fields_period =  1
! The time in seconds after which the field dumps should begin.
integer :: chk_fields_start_time = 0

namelist /auscom_ice_nml/  dt_cpl, &
                   tlthk0,                              & !23/04/2010
                   pop_icediag,                         &
                   do_ice_once,                         & !21/07/2008
                   kmxice,                              & !12/03/2008
                   fixmeltT,                            & !24/04/2008  
                   Tmelt,                               & !24/04/2008
                   limit_srfstress,                     & !20090319
                   Mstress,                             & !20090319
                   use_ioaice,                          & !20090718
                   aice_cutoff,                         & !20090718
                   icemlt_factor,                       & !NO longer needed!
                   frazil_factor,                       & !16/07/2008	     
                   iceform_adj_salt,                    & !20100410
                   sign_stflx,                          & !20100410
                   ksmax,                               &
                   sfix_hours,                          &
                   chk_i2o_fields,                      &
                   chk_o2i_fields,                      &
                   chk_fields_period,                   &
                   chk_fields_start_time
integer :: int_sec

!===========================================================================
end module auscom_ice_parameters_mod



                    module mcm_moist_processes_mod

! MS: Cleaned up moist_processes_op_mod from PJK and renamed it
!     mcm_moist_processes_mod. Deleted unused options and set the default values
!     of the namelist options to agree with supersource values.

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             rel humidity cloud scheme 
!             Diagnostic cloud scheme 
!
!-----------------------------------------------------------------------

use             fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                               write_version_number, stdlog, close_file, &
                               open_namelist_file, file_exist, check_nml_error

use     lscale_cond_mod, only: lscale_cond, lscale_cond_init

use     mcm_mca_lsc_mod, only: mcm_mca_lsc

use  sat_vapor_pres_mod, only: lookup_es

use    time_manager_mod, only: time_type

use    diag_manager_mod, only: register_diag_field, send_data

use         dry_adj_mod, only: dry_adj, dry_adj_init

use       rh_clouds_mod, only: rh_clouds_init, rh_clouds_end, rh_clouds_sum

use      diag_cloud_mod, only: diag_cloud_init, diag_cloud_end, diag_cloud_sum

use   diag_integral_mod, only: diag_integral_field_init, sum_diag_integral_field

use       constants_mod, only: grav, rdgas, rvgas

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   mcm_moist_processes, mcm_moist_processes_init, mcm_moist_processes_end

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

   logical           :: module_is_initialized=.false.
   

   real, parameter :: d622 = rdgas/rvgas
   real, parameter :: d378 = 1.-d622

!--------------------- version number ----------------------------------
   character(len=128) :: version = '$Id: mcm_moist_processes.F90,v 10.0 2003/10/24 22:00:58 fms Exp $'
   character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------

   logical :: do_mca=.false., do_lsc=.false.,  &
              do_dryadj=.true.,               &
              do_rh_clouds=.true., do_diag_clouds=.false.

!------ tracer mapping for stratiform cloud scheme ------

   integer :: nql=0, nqi=0, nqa=0

!  MS mod (use the supersource value for pdepth)
!! real :: pdepth = 150.e2
   real :: pdepth = 4250.
   real :: tfreeze = 273.16

!---------------- namelist variable definitions ------------------------
!
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=true];
!              disabled for mcm_moist_processes.
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true];
!              disabled for mcm_moist_processes.
! do_rh_clouds = switch to turn on/off simple relative humidity cloud scheme
!                [logical, default: do_rh_clouds=true ]
! do_diag_clouds = switch to turn on/off (Gordon's) diagnostic cloud scheme
!                [logical, default: do_diag_clouds=false ]
!  do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=true ]
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!   notes: 1) do_mca  option is disabled for mcm_moist_processes.
!          2) do_lsc  option is disabled for mcm_moist_processes.
!             setting do_mca=.true. or do_lsc=.true. in the namelist will
!             cause error stop.
!          3) pdepth and tfreeze are used to determine liquid vs. solid
!             precipitation for mca, lsc, and ras schemes.
!
!-----------------------------------------------------------------------

namelist /mcm_moist_processes_nml/ do_mca, do_lsc, do_dryadj, pdepth, &
                                   tfreeze, do_rh_clouds, do_diag_clouds

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_tdt_dadj, id_rh

character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine mcm_moist_processes (is, ie, js, je, Time, dt,            &
                                phalf, pfull, t, q, lprec, fprec) ! PJK Mod

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!   out:  lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
! 
!-----------------------------------------------------------------------
integer,         intent(in)              :: is,ie,js,je
type(time_type), intent(in)              :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull
   real, intent(inout),dimension(:,:,:)  :: t, q                         ! PJK
   real, intent(out), dimension(:,:)     :: lprec, fprec

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
logical,dimension(size(t,1),size(t,2))         :: coldT

real, dimension(size(t,1),size(t,2),size(t,3)) :: RH
real, dimension(size(t,1),size(t,2))           :: rain, precip
integer unit

integer :: i, j, k, ix, jx, kx
real    :: dtinv
logical :: used
!-----------------------------------------------------------------------

! The following local quantitities are used exclusively for diagnostic clouds
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CONVPRC  Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IX x JX)

real, dimension(size(t,1),size(t,2),size(t,3)) :: lgscldelq,cnvcntq
real, dimension(size(t,1),size(t,2)) :: convprc

!-----------------------------------------------------------------------

      if (.not.module_is_initialized) call error_mesg ('mcm_moist_processes',  &
                     'mcm_moist_processes_init has not been called.', FATAL)

!-----------------------------------------------------------------------
!-------- input array size and position in global storage --------------

      ix=size(t,1); jx=size(t,2); kx=size(t,3)

!-----------------------------------------------------------------------

      dtinv=1./dt
      lprec=0.0; fprec=0.0; precip=0.0

      tin (:,:,:)=t(:,:,:)
      qin (:,:,:)=q(:,:,:)

!----------------- mean temp in lower atmosphere -----------------------
!----------------- used for determination of rain vs. snow -------------
!----------------- use input temp ? ------------------------------------

   call tempavg (pdepth, phalf, t, tsnow)

   where (tsnow <= tfreeze)
          coldT=.TRUE.
   elsewhere
          coldT=.FALSE.
   endwhere
      
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- dry adjustment scheme -------------------------------

if (do_dryadj) then

         call dry_adj (tin, pfull, phalf, ttnd)
         tin=tin+ttnd
         ttnd=ttnd*dtinv

!------------- diagnostics for dt/dt_dry_adj ---------------------------
     if ( id_tdt_dadj > 0 ) then
        used = send_data ( id_tdt_dadj, ttnd, Time, is, js, 1)
     endif
! ----------------------------------------------------------------------
end if

!-----------------------------------------------------------------------
!***********************************************************************
!----------------- moist convective adjustment -------------------------

!----------------- supersource moist convective adjustment and large-scale condensation

!-----------------------------------------------------------------------

!  TK mod:
!  call mcm_mca_lsc(tin,qin,phalf,ttnd,qtnd,rain,snow)
   call mcm_mca_lsc(tin,qin,phalf,ttnd,qtnd,rain,snow, RH)


   !pass RH to rh_clouds_sum
   call rh_clouds_sum (is, js, RH)
   

!------- update input values and compute tendency -------
                    
      tin=tin+ttnd;    qin=qin+qtnd
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv
      
!------- add on tendency ----------
!    tdt=tdt+ttnd; qdt=qdt+qtnd


!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow
!------- compute rh clouds if desired ------
! Deleted; see TK NOTE above. RH calculation now done in mcm_mca_lsc
!    if (do_rh_clouds) then
!          
!          !calculate relative humidity
!          call rh_calc(pfull,tin,qin,RH,mask)
!          
!          !pass RH to rh_clouds_sum
!          call rh_clouds_sum (is, js, RH)
!           
!    end if

!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR CONVECTIVE SCHEME ---------------------
!-----------------------------------------------------------------------
!------- diagnostics for tdt_conv -------
      if ( id_tdt_conv > 0 ) then
        used = send_data ( id_tdt_conv, ttnd, Time, is, js, 1)
      endif
!------- diagnostics for qdt_conv -------
      if ( id_qdt_conv > 0 ) then
        used = send_data ( id_qdt_conv, qtnd, Time, is, js, 1)
      endif
!------- diagnostics for precip -------
      if ( id_prec_conv > 0 ) then
        used = send_data ( id_prec_conv, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow -------
      if ( id_snow_conv > 0 ) then
        used = send_data ( id_snow_conv, snow, Time, is, js )
      endif

!-----------------------------------------------------------------------
                      if (do_diag_clouds) then
!  capture convective precip and convective spec hum changes (and calculate, 
!  convective spec hum counter) which are needed as predictors 
!  for Gordon's diagnostic clouds
      where (qtnd(:,:,:) < 0.0)
            cnvcntq (:,:,:) = 1.0
      else where
            cnvcntq (:,:,:) = 0.0
      end where
      convprc = precip
                      endif

!-----------------------------------------------------------------------
!***********************************************************************
!--------------------- GENERAL DIAGNOSTICS -----------------------------
!-----------------------------------------------------------------------

!------- diagnostics for total precip -------
   if ( id_precip > 0 ) then
        used = send_data ( id_precip, precip, Time, is, js )
   endif


!-----------------------------------------------------------------------
!------- diagnostics for relative humidity -------

   if ( id_rh > 0 ) then
      if (.not.(do_rh_clouds.or.do_diag_clouds)) &
          call rh_calc (pfull,tin,qin,RH)
      used = send_data ( id_rh, RH*100., Time, is, js, 1 )
   endif

!-----------------------------------------------------------------------
!---- accumulate global integral of precipiation (mm/day) -----

call sum_diag_integral_field ('prec', precip*86400., is, js)

!-----------------------------------------------------------------------
  t = tin                          ! PJK
  q = qin                          ! PJK
end subroutine mcm_moist_processes

!#######################################################################

subroutine mcm_moist_processes_init ( id, jd, kd, axes, Time )

!-----------------------------------------------------------------------
integer,         intent(in) :: id, jd, kd, axes(4)
type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer  unit,io,ierr,nt
!-----------------------------------------------------------------------

       if ( module_is_initialized ) return

       if ( file_exist('input.nml')) then

         unit = open_namelist_file()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=mcm_moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'mcm_moist_processes_nml')
         enddo
  10     call close_file (unit)

!------------------- dummy checks --------------------------------------

         if ( do_lsc ) call error_mesg   &
                   ('mcm_moist_processes_init ',  &
                    'do_lsc is not an option in mcm_moist_processes', FATAL)
         if ( do_mca ) call error_mesg   &
                   ('mcm_moist_processes_init ',  &
                    'do_mca is not an option in mcm_moist_processes', FATAL)

      endif

!--------- write namelist ------------------

      call write_version_number(version, tagname)
      if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=mcm_moist_processes_nml)

!------------ initialize various schemes ----------

      if (do_dryadj) call dry_adj_init ()
      if (do_rh_clouds)  call rh_clouds_init (id,jd,kd)     ! PJK: added
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')

!----- initialize quantities for diagnostics output -----

   call diag_field_init ( axes, Time )


   module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine mcm_moist_processes_init

!#######################################################################

subroutine mcm_moist_processes_end

integer  unit
!-----------------------------------------------------------------------

      if ( .not.module_is_initialized ) return

!----------------close various schemes-----------------

      if (do_rh_clouds)   call rh_clouds_end
      if (do_diag_clouds) call diag_cloud_end

!-----------------------------------------------------------------------

      module_is_initialized = .false.

end subroutine mcm_moist_processes_end

!#######################################################################

      subroutine tempavg (pdepth,phalf,temp,tsnow)

!-----------------------------------------------------------------------
!
!    computes a mean atmospheric temperature for the bottom
!    "pdepth" pascals of the atmosphere.
!
!   input:  pdepth     atmospheric layer in pa.
!           phalf      pressure at model layer interfaces
!           temp       temperature at model layers
!
!   output:  tsnow     mean model temperature in the lowest
!                      "pdepth" pascals of the atmosphere
!
!-----------------------------------------------------------------------
      real, intent(in)  :: pdepth
      real, intent(in) , dimension(:,:,:) :: phalf,temp
      real, intent(out), dimension(:,:)   :: tsnow
!-----------------------------------------------------------------------
 real, dimension(size(temp,1),size(temp,2)) :: prsum, done, pdel, pdep
 real  sumdone
 integer  k
!-----------------------------------------------------------------------

      tsnow=0.0; prsum=0.0; done=1.0; pdep=pdepth

      do k=size(temp,3),1,-1

         pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*done(:,:)

         where ((prsum(:,:)+pdel(:,:))  >  pdep(:,:))
            pdel(:,:)=pdepth-prsum(:,:)
            done(:,:)=0.0
            pdep(:,:)=0.0
         endwhere

         tsnow(:,:)=tsnow(:,:)+pdel(:,:)*temp(:,:,k)
         prsum(:,:)=prsum(:,:)+pdel(:,:)

         sumdone=sum(done(:,:))
         if (sumdone < 1.e-4) exit

      enddo

         tsnow(:,:)=tsnow(:,:)/prsum(:,:)

!-----------------------------------------------------------------------

      end subroutine tempavg

!#######################################################################

      subroutine rh_calc(pfull,T,qv,RH)

        IMPLICIT NONE


        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH

        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat
        
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = Rdgas/RVgas = d622
!
!       and where 1- epsilon = d378
!
!       Note that RH does not have its proper value
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        call lookup_es(t(:,:,:),esat(:,:,:))
        
        !calculate denominator in qsat formula
        rh(:,:,:) = pfull(:,:,:)-d378*esat(:,:,:)
     
        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        rh(:,:,:) = max(rh(:,:,:),esat(:,:,:)) 
        
        !calculate rh
        rh(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/rh(:,:,:))
      

END SUBROUTINE rh_calc



!#######################################################################

subroutine diag_field_init ( axes, Time )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  integer, dimension(3) :: half = (/1,2,4/)

!------------ initializes diagnostic fields in this module -------------

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from Manabe Climate Model conv adj',     'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Mixing ratio tendency from Manabe Climate Model conv adj',   'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from Manabe Climate Model conv adj',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from Manabe Climate Model conv adj',       'kg/m2/s' )

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s' )

   id_tdt_dadj = register_diag_field ( mod_name, &
     'tdt_dadj', axes(1:3), Time, &
   'Temperature tendency from dry conv adj',       'deg_K/s',  &
                        missing_value=missing_value               )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  & 
                        missing_value=missing_value               )

!-----------------------------------------------------------------------

end subroutine diag_field_init

!#######################################################################

                 end module mcm_moist_processes_mod


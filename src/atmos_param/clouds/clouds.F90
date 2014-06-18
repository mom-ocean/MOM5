
                      module clouds_mod

!=======================================================================
!
!            determines cloud properties necessary for 
!                    fels-schwartzkopf radiation
!
!=======================================================================

use    cloud_rad_mod, only:  cloud_rad_init, cloud_summary
use  cloud_zonal_mod, only:  cloud_zonal
use    cloud_obs_mod, only:  cloud_obs, cloud_obs_init
use time_manager_mod, only:  time_type
use          mpp_mod, only:  input_nml_file
use          fms_mod, only:  error_mesg, FATAL, file_exist,   &
                             check_nml_error, open_namelist_file,      &
                             mpp_pe, mpp_root_pe, close_file, &
                             write_version_number, stdlog
use    rh_clouds_mod, only:  do_rh_clouds, rh_clouds, rh_clouds_avg
use  strat_cloud_mod, only:  do_strat_cloud, strat_cloud_avg
use   diag_cloud_mod, only:  do_diag_cloud, diag_cloud_driver, &
                             diag_cloud_avg
use diag_manager_mod, only:  register_diag_field, send_data

implicit none
private

!------------------- public interfaces ---------------------------------

public   clouds, clouds_init, clouds_end

!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------
 character(len=128) :: version = '$Id: clouds.F90,v 19.0 2012/01/06 20:02:48 fms Exp $'
 character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   note:  the fels-schwarzkopf radiation code permits bi-spectral
!          cloud reflectivity associated with cloud cdwtr droplets:
!            -->  visible band - (cao3sw and cuvrf);
!            -->  near infra-red band - (cah2sw and cirrf).
!          the f-s code assumes that all gaseous absorption by
!          cdwtr vapor occurs in the near infra-red band.
!          thus, original code contains cbsw and cirab.
!          we shall include cbo3sw and cuvab and let cbsw = cbh2sw.
!          however, these spectral absorptivities will be set to zero.

      real, dimension(3) :: cao3sw = (/ 0.210, 0.450, 0.590 /)
      real, dimension(3) :: cah2sw = (/ 0.210, 0.450, 0.590 /)
      real, dimension(3) :: cbsw   = (/ 0.005, 0.020, 0.035 /)


!-----------------------------------------------------------------------

      logical :: module_is_initialized =.false.

      integer :: id_tot_cld_amt, id_high_cld_amt, id_mid_cld_amt, &
                 id_low_cld_amt, id_cld_amt, id_em_cld,  &
                 id_alb_uv_cld, id_alb_nir_cld,          &
                 id_abs_uv_cld, id_abs_nir_cld

      character(len=6), parameter :: mod_name = 'clouds'

      real :: missing_value = -999.

!-----------------------------------------------------------------------
!------------------------- namelist ------------------------------------

      logical :: do_zonal_clouds = .false.
      logical :: do_obs_clouds   = .false.
      logical :: do_no_clouds    = .false.
      logical :: do_isccp_cloud_diags = .false.
      
      namelist /clouds_nml/ do_zonal_clouds,  &
                            do_obs_clouds,    &
                            do_no_clouds,     &
                            do_isccp_cloud_diags

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine clouds  (is, js, clear_sky, Time, Time_diag, lat, &
                    land, tsfc, pfull, phalf, t, q, cosz,    &
                    nclds, ktopsw, kbtmsw, ktoplw, kbtmlw,   &
                    cldamt, cuvrf, cirrf, cirab, emcld, mask, kbot)

!-----------------------------------------------------------------------
        integer, intent(in)                    :: is, js
        logical, intent(in)                    :: clear_sky
type(time_type), intent(in)                    :: Time, Time_diag

   real, intent(in), dimension(:,:)    :: lat
   real, intent(in), dimension(:,:)    :: land,tsfc
   real, intent(in), dimension(:,:,:)  :: pfull,phalf,t,q
   real, intent(in), dimension(:,:)    :: cosz
integer, intent(out), dimension(:,:)   :: nclds
integer, intent(out), dimension(:,:,:) :: ktopsw,kbtmsw,ktoplw,kbtmlw
   real, intent(out), dimension(:,:,:) :: cldamt,cuvrf,cirrf,cirab,emcld
   real, intent(in),  dimension(:,:,:),optional :: mask
integer, intent(in),  dimension(:,:),  optional :: kbot
!-----------------------------------------------------------------------
   real,dimension(size(cirab,1),size(cirab,2),size(cirab,3)) :: cuvab
integer,dimension(size(ktoplw,1),size(ktoplw,2),size(ktoplw,3)) ::  &
                       ktop, kbtm
!      TCA_CA   array for total clouds diagnositc
!      HML_CA   array for high, middle and low clouds diagnostics
!               as per 3rd index values of 1,2 and 3 respectively.
   real,dimension(size(pfull,1),size(pfull,2))   :: tca
   real,dimension(size(pfull,1),size(pfull,2),3) :: hml_ca
 
!      pflux    array for the flux pressure levels for isccp cloud
!               diagnostics calculations
   real,dimension(size(phalf,1),size(phalf,2),size(phalf,3))  :: pflux
   

   real,dimension(size(pfull,1),size(pfull,2),size(pfull,3)) ::  &
                                                     ql,qi,cf,rh,cloud
   real,dimension(size(phalf,1),size(phalf,2),size(phalf,3)) :: phaf
   real :: rad2deg
integer :: i,j,k,kb,kx,kp1,n,ierr
logical :: used
!-----------------------------------------------------------------------

! The following local quantitities are used exclusively for diagnostic clouds
!      TEMP     Averaged temperature (Deg K) at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      QMIX     Averaged specific humidity at full model levels 
!                   (dimensioned IDIM x JDIM x kx)
!      OMEGA  Averaged pressure vertical velocity at full model levels
!                   (dimensioned IDIM x JDIM x kx)
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CONVPRC  Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IX x JX)
!      Note:    RH is also a predictor but since it is used in the rh cloud
!               scheme as well it has already been declared.
!      PSFC     Surface pressure
!               (dimensioned IX x JX)
real, dimension(size(t,1),size(t,2),size(t,3)) :: temp,qmix,omega
real, dimension(size(t,1),size(t,2),size(t,3)) :: lgscldelq,cnvcntq
real, dimension(size(t,1),size(t,2)) :: convprc,psfc

!-----------------------------------------------------------------------

  ierr = 1

  kx  = size(ktopsw,3)-1
  kp1 = kx+1
    
  if (kx /= size(pfull,3)) call error_mesg ('clouds in clouds_mod', &
                       'input arrays have the incorrect size.',FATAL)

!-----------------------------------------------------------------------
!----------- default clouds values ----------

      call default_clouds (nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                           cldamt,cuvrf,cirrf,cuvab,cirab,emcld)

!-------------- no clouds ----------------------------------------------

      if (clear_sky .or. do_no_clouds) then
          if (present(kbot)) call step_mtn_clouds (kx,kbot,          &
                                 nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                                 cldamt,cuvrf,cirrf,cirab,emcld)
          return
      endif

!-----------------------------------------------------------------------
!--------------- determine rh cloud properties -----------------
 if ( do_rh_clouds() ) then
!-----------------------------------------------------------------------
!---- compute rh_clouds -----

     call rh_clouds_avg (is, js, rh, ierr)

     if (ierr == 0) then
         rad2deg = 90./acos(0.0)
         call rh_clouds(rh,pfull,phalf(:,:,kx+1),cosz,lat*rad2deg,&
                        nclds,ktop(:,:,2:kp1),kbtm(:,:,2:kp1),  &
                        cldamt(:,:,2:kp1),cuvrf(:,:,2:kp1),  &
                        cirrf(:,:,2:kp1),cuvab(:,:,2:kp1),  &
                        cirab(:,:,2:kp1),emcld(:,:,2:kp1))
     endif     

!-----------------------------------------------------------------------
 endif
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!--------------- determine prognostic cloud properties -----------------
 if ( do_strat_cloud() ) then
!-----------------------------------------------------------------------
      
     call strat_cloud_avg (is, js, ql, qi, cf, ierr)

     if (ierr == 0) then
         call cloud_summary (is,js,land,ql,qi,cf,q,pfull, &
                             phalf,t,cosz,tsfc,&
                             nclds,ktop(:,:,2:kp1),kbtm(:,:,2:kp1),  &
                             cldamt(:,:,2:kp1), &
                             Time=Time_diag, &
                             r_uv=cuvrf(:,:,2:kp1), &
                             r_nir=cirrf(:,:,2:kp1), &
                             ab_uv=cuvab(:,:,2:kp1), &
                             ab_nir=cirab(:,:,2:kp1), &
                             em_lw=emcld(:,:,2:kp1))
     endif     

!-----------------------------------------------------------------------
 endif     
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!--------------- determine diagnostic cloud properties -----------------
 if ( do_diag_cloud() ) then
!-----------------------------------------------------------------------
      
     call diag_cloud_avg (is, js, temp,qmix,rh,omega, &
                          lgscldelq,cnvcntq,convprc,ierr)

     psfc(:,:) = phalf(:,:,kp1)
     if (ierr == 0) then
         call diag_cloud_driver (is,js, &
                    temp,qmix,rh,omega,lgscldelq,cnvcntq,convprc, &
                    pfull,phalf,psfc,cosz,lat,Time, &
                    nclds,ktop(:,:,2:kp1),kbtm(:,:,2:kp1), &
                    cldamt(:,:,2:kp1),cuvrf(:,:,2:kp1),  &
                    cirrf(:,:,2:kp1),cuvab(:,:,2:kp1),  &
                    cirab(:,:,2:kp1),emcld(:,:,2:kp1) )
!                    cirab(:,:,2:kp1),emcld(:,:,2:kp1) ,kbot)
     endif     

!-----------------------------------------------------------------------
 endif     
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  ----------- zonal or observed clouds ??? --------------
!  (also do if avg cloud properties could not be returned)
!-----------------------------------------------------------------------
 if ( (do_zonal_clouds .or. do_obs_clouds) .and. ierr /= 0  ) then
!-----------------------------------------------------------------------

!    ---- constrain phalf to:  0. <= phalf <= 101325. ----
      if (present(kbot)) then
         do k=1,kx+1; do j=1,size(phalf,2); do i=1,size(phalf,1)
            kb=kbot(i,j)
            phaf(i,j,k)=101325.*phalf(i,j,k)/phalf(i,j,kb+1)
         enddo; enddo; enddo
      else
         do k=1,kx+1
            phaf(:,:,k)=101325.*phalf(:,:,k)/phalf(:,:,kx+1)
         enddo
      endif

!   ---- get three cloud levels (store in k=2,4) -----

      call cloud_zonal (Time, lat, phaf, nclds,  &
                        ktopsw(:,:,2:4), kbtmsw(:,:,2:4), &
                        ktoplw(:,:,2:4), kbtmlw(:,:,2:4), &
                        cldamt(:,:,2:4), cuvrf(:,:,2:4),  &
                        cirrf(:,:,2:4), cirab(:,:,2:4),  &
                        emcld(:,:,2:4))

      if (do_obs_clouds) call cloud_obs (is, js, Time, cldamt(:,:,2:4))

!-----------------------------------------------------------------------
 endif
!-----------------------------------------------------------------------


!----- store longwave and shortwave indices -----

   if ( (do_rh_clouds() .or. do_strat_cloud() .or. do_diag_cloud()) &
        .and. ierr == 0 ) then

         do j=1,size(ktoplw,2)
         do i=1,size(ktoplw,1)
            if (nclds(i,j) > 0) then
               n=nclds(i,j)
               ktoplw(i,j,2:n+1)=ktop(i,j,2:n+1)
               kbtmlw(i,j,2:n+1)=kbtm(i,j,2:n+1)
               ktopsw(i,j,2:n+1)=ktop(i,j,2:n+1)
               kbtmsw(i,j,2:n+1)=kbtm(i,j,2:n+1)+1
            endif
         enddo
         enddo
   endif


!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!---- total cloud diagnostic ----
      if ( id_tot_cld_amt > 0 ) then
         call compute_tca_random ( nclds, cldamt(:,:,2:kp1), tca )
         used = send_data ( id_tot_cld_amt, tca, Time_diag, is, js )
      endif

!---- high,mid,low cloud diagnostics ----
      if ( id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
            id_low_cld_amt > 0 ) then

         if (do_isccp_cloud_diags) then
!           Do alternative (isccp) high,mid,low cloud methodology.            
!           Use methodology adopted from cloudrad_package:

!           Construct the pflux array, which is the flux pressure
!           level midway between the pressure levels in the model
!           where radiation is actually computed.  Add surface
!           pressure as the final vertical level of the 3-D array.
!           The first element of the array is set to zero.
!           The array which uses these is based on cgs units, 
!           so also multiply by 10 to convert from Pa to dynes/cm2.

            if (present(kbot)) then
               call error_mesg ('clouds in clouds_mod', &
                  'compute_isccp_clds not set up/tested on case where kbot is present.',FATAL)
            end if
     
            pflux(:,:,1) = 0.
            do k = 2, kx
               pflux(:,:,k) = 0.5 * (pfull(:,:,k-1)+pfull(:,:,k)) * 10.0
            end do
            pflux(:,:,kp1) = phalf(:,:,kp1) * 10.0  
       
            call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             cldamt(:,:,2:kp1), cloud )
                             
            call compute_isccp_clds ( pflux, cloud, hml_ca)               

         else
!           Use original methodology contained in this module:            
            call compute_hml_ca_random ( nclds, cldamt(:,:,2:kp1), &
                    kbtmlw(:,:,2:kp1), pfull, phalf, hml_ca, kbot )
         end if           
                    
         if ( id_high_cld_amt > 0 ) used = send_data &
                     ( id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js )
         if ( id_mid_cld_amt > 0 ) used = send_data &
                      ( id_mid_cld_amt, hml_ca(:,:,2), Time_diag, is, js )
         if ( id_low_cld_amt > 0 ) used = send_data &
                      ( id_low_cld_amt, hml_ca(:,:,3), Time_diag, is, js )
      endif

!------- cloud amount -------------------------
      if ( id_cld_amt > 0 ) then
         call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             cldamt(:,:,2:kp1), cloud )
         used = send_data ( id_cld_amt, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- cloud emissivity ---------------------------------------
      if ( id_em_cld > 0 ) then
         call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             emcld(:,:,2:kp1), cloud )
         used = send_data ( id_em_cld, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- ultra-violet reflected by cloud -----------------------------
      if ( id_alb_uv_cld > 0 ) then
         call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             cuvrf(:,:,2:kp1), cloud )
         used = send_data ( id_alb_uv_cld, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- infra-red reflected by cloud -----------------------------
      if ( id_alb_nir_cld > 0 ) then
         call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             cirrf(:,:,2:kp1), cloud )
         used = send_data ( id_alb_nir_cld, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- ultra-violet absorbed by cloud (not implemented)------------
!     if ( id_abs_uv_cld > 0 ) then
!        call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
!                            cuvab(:,:,2:kp1), cloud )
!        used = send_data ( id_abs_uv_cld, cloud, Time_diag, is, js, 1, &
!                           rmask=mask )
!     endif

!------- infra-red absorbed by cloud -----------------------------
      if ( id_abs_nir_cld > 0 ) then
         call expand_cloud (nclds, ktoplw(:,:,2:kp1),kbtmlw(:,:,2:kp1), &
                             cirab(:,:,2:kp1), cloud )
         used = send_data ( id_abs_nir_cld, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!--------------END OF DIAGNOSTICS --------------------------------------
!-----------------------------------------------------------------------
!---- step mountain clouds in underground levels ----

     if (present(kbot)) call step_mtn_clouds (kx,kbot,              &
                                nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                                cldamt,cuvrf,cirrf,cirab,emcld)

!-----------------------------------------------------------------------

      end subroutine clouds

!#######################################################################

 subroutine compute_tca_random ( nclds, cldamt, tca )

   integer, intent(in)  :: nclds (:,:)
   real,    intent(in)  :: cldamt(:,:,:)
   real,    intent(out) :: tca   (:,:)

   integer :: k, max_cld

!---- find maximum number of clouds -----

    max_cld = maxval(nclds)

!---- compute total cloud amount assuming that -----
!       independent clouds overlap randomly

    tca = 1.0

    do k = 1, max_cld
       tca(:,:) = tca(:,:) * (1. - cldamt(:,:,k))
    enddo

    tca = (1. - tca) * 100.


 end subroutine compute_tca_random

!#######################################################################

 subroutine compute_hml_ca_random ( nclds, cldamt, kbtm, pfull, phalf, &
                                    hml_ca, kbot )

   integer, intent(in)  :: nclds (:,:), kbtm(:,:,:)
   real,    intent(in)  :: cldamt(:,:,:), pfull(:,:,:), phalf(:,:,:)
   real,    intent(out) :: hml_ca   (:,:,:)
   integer, intent(in), optional :: kbot(:,:)

! hml_ca - array for high, middle and low clouds as per 3rd index values
!          of 1,2 and 3 respectively.

! ---------------------------------------------------------------------

   integer,  parameter :: mid_btm = 7.0e4,high_btm = 4.0e4

! local array
real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: pfull_norm
! pfull_norm is normalized pressure at full model levels

   integer :: i, j, k, kx, kb

   kx   = size(pfull,3)

!  calculate normalized presure to allow full range of high middle low 
!  clouds independent of orography

   if (present(kbot)) then
         do k=1,kx+1; do j=1,size(phalf,2); do i=1,size(phalf,1)
            kb=kbot(i,j)
            pfull_norm(i,j,k)=101325.*phalf(i,j,k)/phalf(i,j,kb+1)
         enddo; enddo; enddo
   else
         do k=1,kx
            pfull_norm(:,:,k)=101325.*pfull(:,:,k)/phalf(:,:,kx+1)
         enddo
   endif

!---- compute high, middle and low cloud amounts assuming that -----
!       independent clouds overlap randomly

    hml_ca = 1.0

    do j=1,size(nclds,2)
    do i=1,size(nclds,1)
!   if (nclds(i,j) == 0) cycle

      do k = 1, nclds(i,j)
         if (pfull_norm(i,j,kbtm(i,j,k)) .le. high_btm) then
           hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - cldamt(i,j,k))
         else if ( (pfull_norm(i,j,kbtm(i,j,k)) .le. mid_btm) .and.  &
                   (pfull_norm(i,j,kbtm(i,j,k)) .gt. high_btm) ) then
           hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - cldamt(i,j,k))
         else 
           hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - cldamt(i,j,k))
         endif
      enddo

    enddo
    enddo

    hml_ca = (1. - hml_ca) * 100.


 end subroutine compute_hml_ca_random

!#######################################################################

subroutine compute_isccp_clds ( pflux, cloud, hml_ca)

real,  dimension(:,:,:),   intent(in)  :: pflux, cloud
real,  dimension(:,:,:),   intent(out) :: hml_ca

!
!   define arrays giving the fractional cloudiness for clouds with
!   tops within the ISCCP definitions of high (10-440 hPa), middle
!   (440-680 hPa) and low (680-1000 hPa).
 
!    note that at this point pflux is in cgs units. change this later.

!    This routine is copied from cloudrad_package_mod, where
!    it is private and thus non-accessible from here directly.
!    Modified to work within this routine...
 
! ---------------------------------------------------------------------

   real,  parameter :: mid_btm = 6.8e5,high_btm = 4.4e5
  
! local array

integer :: i, j, k

 
!---- compute high, middle and low cloud amounts assuming that -----
!       independent clouds overlap randomly

    hml_ca = 1.0
    
 
   do j=1, size(cloud,2)
    do i=1, size(cloud,1)
      do k = 1, size(cloud,3)
        if (pflux(i,j,k)  <=  high_btm) then
          hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - cloud(i,j,k))
        else if ( (pflux(i,j,k) >  high_btm) .and.  &
           (pflux(i,j,k) <=  mid_btm) ) then
         hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - cloud(i,j,k))
       else  if ( pflux(i,j,k) > mid_btm ) then
         hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - cloud(i,j,k))
       endif
    enddo
  enddo
  enddo

    hml_ca = 1. - hml_ca
    hml_ca = 100. * hml_ca
  
end subroutine compute_isccp_clds
!#######################################################################

      subroutine default_clouds (nclds,ktopsw,kbtmsw,ktop,kbtm,  &
                                 cldamt,cuvrf,cirrf,cuvab,cirab,emcld)

!----------------------------------------------------------------------
   integer, intent(inout), dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktopsw,kbtmsw,ktop,kbtm
      real, intent(inout), dimension(:,:,:) :: cldamt,cuvrf,cirrf,  &
                                               cuvab,cirab,emcld
!----------------------------------------------------------------------
      integer  kp1

      kp1=size(ktopsw,3)

      nclds(:,:)=0

      cldamt=0.0; emcld =1.0
      cuvrf =0.0; cirrf =0.0; cuvab =0.0; cirab =0.0
      ktopsw=kp1; kbtmsw=kp1
!     ktop  =kp1; kbtm  =kp1
      ktop  =kp1-1; kbtm  =kp1-1
!     ---- reset top properties ----
      ktopsw(:,:,1)=1
      kbtmsw(:,:,1)=0
      ktop  (:,:,1)=1
      kbtm  (:,:,1)=0

!----------------------------------------------------------------------

      end subroutine default_clouds

!#######################################################################

      subroutine step_mtn_clouds (kmax,kbot,nclds,  &
                                  ktopsw,kbtmsw,ktop,kbtm,  &
                                  cldamt,cuvrf,cirrf,cirab,emcld)

!-----------------------------------------------------------------------
   integer, intent(in)                      :: kmax
   integer, intent(in),    dimension(:,:)   :: kbot
   integer, intent(inout), dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktopsw,kbtmsw,ktop,kbtm
      real, intent(inout), dimension(:,:,:) :: cldamt,cuvrf,cirrf,  &
                                                      cirab,emcld
!-----------------------------------------------------------------------
   integer  i,j,n,kp1

   kp1=kmax+1

         do j=1,size(kbot,2)
         do i=1,size(kbot,1)
            if (kbot(i,j) < kmax) then
               nclds(i,j)=nclds(i,j)+1
               n=nclds(i,j)+1
               ktopsw(i,j,n)=kbot(i,j)+1
               kbtmsw(i,j,n)=kp1
               ktop  (i,j,n)=kbot(i,j)+1
               kbtm  (i,j,n)=kmax
               cldamt(i,j,n)=1.0
               cuvrf (i,j,n)=0.0
               cirrf (i,j,n)=0.0
               cirab (i,j,n)=0.0
               emcld (i,j,n)=1.0
            endif
         enddo
         enddo

!-----------------------------------------------------------------------

      end subroutine step_mtn_clouds

!#######################################################################

      subroutine remove_cloud_overlap (nclds,ktop,kbtm)

!-----------------------------------------------------------------------
   integer, intent(in),    dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktop,kbtm
!-----------------------------------------------------------------------
!                   removes sw cloud overlap
!
!    nclds = number of clouds
!    ktop  = sw indice for cloud top
!    kbtm  = sw indice for cloud bottom
!
!-----------------------------------------------------------------------
   integer  i,j,kc

         do j=1,size(nclds,2)
         do i=1,size(nclds,1)
            if (nclds(i,j) <= 1) cycle

            do kc=2,nclds(i,j)
               if (ktop(i,j,kc+1) >= kbtm(i,j,kc)) cycle
!    case 1: thin or thick upper cloud, thick lower cloud
                  if (ktop(i,j,kc+1) <  kbtm(i,j,kc+1)) then
                     ktop(i,j,kc+1)=ktop(i,j,kc+1)+1
                  else
!    case 2: thick upper cloud, thin lower cloud
                     kbtm(i,j,kc)=kbtm(i,j,kc)-1
                  endif
            enddo

         enddo
         enddo

!-----------------------------------------------------------------------

      end subroutine remove_cloud_overlap

!#######################################################################

      subroutine expand_cloud ( nclds, ktop, kbtm, cloud_in, cloud_out )

      integer, intent(in)  :: nclds(:,:), ktop(:,:,:), kbtm(:,:,:)
      real,    intent(in)  :: cloud_in (:,:,:)
      real,    intent(out) :: cloud_out(:,:,:)

      integer :: i, j, n

         cloud_out = 0.0
         do j=1,size(nclds,2)
         do i=1,size(nclds,1)
            do n=1,nclds(i,j)
              cloud_out(i,j,ktop(i,j,n):kbtm(i,j,n)) = cloud_in(i,j,n)
            enddo
         enddo
         enddo

      end subroutine expand_cloud

!#######################################################################

      subroutine clouds_init ( lonb, latb, axes, Time )

!-----------------------------------------------------------------------
           real, intent(in), dimension(:,:) :: lonb, latb
        integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)               :: Time
!-----------------------------------------------------------------------
      integer  unit,io,ierr, logunit

!-------------- read namelist --------------

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=clouds_nml, iostat=io)
      ierr = check_nml_error(io,"clouds_nml")
#else
      if ( file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=clouds_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'clouds_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

!      ----- write namelist -----

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           logunit = stdlog()
           write (logunit, nml=clouds_nml)
      endif


      if (do_obs_clouds) call cloud_obs_init (lonb, latb)

!------------ initialize diagnostic fields -----------------------------

      call diag_field_init ( Time, axes )

!-----------------------------------------------------------------------
       if (do_strat_cloud() ) then
         call cloud_rad_init
       endif

      module_is_initialized =.true.

!-----------------------------------------------------------------------

      end subroutine clouds_init

!#######################################################################

      subroutine clouds_end

!-----------------------------------------------------------------------

      module_is_initialized =.false.

!-----------------------------------------------------------------------

      end subroutine clouds_end

!#######################################################################

   subroutine diag_field_init ( Time, axes )

     type(time_type), intent(in) :: Time
     integer        , intent(in) :: axes(4)

!-----------------------------------------------------------------------

    id_tot_cld_amt = &
    register_diag_field ( mod_name, 'tot_cld_amt', axes(1:2), Time, &
                         'total cloud amount', 'percent'            )

    id_high_cld_amt = &
    register_diag_field ( mod_name, 'high_cld_amt', axes(1:2), Time, &
                         'high cloud amount', 'percent'            )

    id_mid_cld_amt = &
    register_diag_field ( mod_name, 'mid_cld_amt', axes(1:2), Time, &
                         'mid cloud amount', 'percent'            )

    id_low_cld_amt = &
    register_diag_field ( mod_name, 'low_cld_amt', axes(1:2), Time, &
                         'low cloud amount', 'percent'            )

    id_cld_amt = &
    register_diag_field ( mod_name, 'cld_amt', axes(1:3), Time, &
                         'cloud amount', 'percent',             &
                         missing_value=missing_value            )

    id_em_cld = &
    register_diag_field ( mod_name, 'em_cld', axes(1:3), Time, &
                         'cloud emissivity', 'none',           &
                          missing_value=missing_value          )

    id_alb_uv_cld = &
    register_diag_field ( mod_name, 'alb_uv_cld', axes(1:3), Time, &
                         'UV reflected by cloud', 'none',          &
                          missing_value=missing_value              )

    id_alb_nir_cld = &
    register_diag_field ( mod_name, 'alb_nir_cld', axes(1:3), Time, &
                         'IR reflected by cloud', 'none',           &
                          missing_value=missing_value               )

!   --- do not output this field ---
!   id_abs_uv_cld = &
!   register_diag_field ( mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                        'UV absorbed by cloud', 'none',           &
!                         missing_value=missing_value              )

    id_abs_nir_cld = &
    register_diag_field ( mod_name, 'abs_nir_cld', axes(1:3), Time, &
                         'IR absorbed by cloud', 'none',            &
                          missing_value=missing_value               )

!-----------------------------------------------------------------------

   end subroutine diag_field_init

!#######################################################################

end module clouds_mod

